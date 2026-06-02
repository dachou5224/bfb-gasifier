"""Reactor orchestration for thesis-aligned NR solver."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import logging
from time import perf_counter
from typing import List

import numpy as np

logger = logging.getLogger(__name__)

from src.core.cell import (
    Cell,
    CellGeometry,
    SolidProps,
    S_CHAR,
    S_VM,
    S_MOISTURE,
    S_ASH,
    N_SOLID_COMP,
)
from src.core.connectivity import (
    apply_all_nr_boundary_data as _apply_all_nr_boundary_data_exec,
    apply_bottom_recycle as _apply_bottom_recycle_exec,
    apply_local_nr_boundary_data as _apply_local_nr_boundary_data_exec,
    propagate_upstream as _propagate_upstream_exec,
    propagated_solid_stream as _propagated_solid_stream,
    recycled_solid_stream as _recycled_solid_stream,
    set_bottom_cell_feeds as _set_bottom_cell_feeds_exec,
)
from src.core.connectivity_graph import affected_nr_residual_cells as _affected_nr_cells_graph
from src.core.freeboard_bridge import (
    refresh_explicit_freeboard_transport_from_closure as _refresh_fb_transport_bridge,
    sync_freeboard_cells_from_closure as _sync_fb_cells_bridge,
)
from src.core.side_block_bridge import initialize_explicit_side_block_states as _init_side_blocks_bridge
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    TarFuelType,
    configure_tar_components_by_fuel,
)
from src.core.feed_inlet import compute_gas_feeds_mol_s
from src.core.legacy_gs_support import (
    _capture_cells_state,
    _clamp_soft_rollback_monotonic,
    _gs_profile_metrics,
    _positive_rebound_penalty,
    _restore_cells_state,
)


def _compound_heat_loss_fraction(total_frac: float, weight: float) -> float:
    """Map a total heat-loss fraction onto a weighted axial segment."""
    total = float(np.clip(total_frac, 0.0, 1.0))
    seg_weight = float(max(weight, 0.0))
    if total <= 0.0 or seg_weight <= 0.0:
        return 0.0
    return float(1.0 - math.pow(max(1.0 - total, 0.0), seg_weight))


def _resolve_axial_heat_loss_distribution(cfg: "ReactorConfig") -> tuple[list[float], float]:
    """Distribute reactor-level heat loss across bed cells and freeboard.

    `heat_loss_frac` is interpreted as the total loss fraction over the full
    `bed + freeboard` height. `freeboard_heat_loss_frac` remains available as an
    additional freeboard-only increment for legacy sensitivity studies.
    """
    bed_height = float(max(cfg.H_bed, 0.0))
    freeboard_height = float(max(cfg.H_freeboard, 0.0))
    total_height = bed_height + freeboard_height

    if total_height <= 0.0:
        return [0.0] * int(max(cfg.n_cells, 0)), float(np.clip(cfg.freeboard_heat_loss_frac, 0.0, 1.0))

    n_bed = int(max(cfg.n_cells, 0))
    bed_weight_each = (bed_height / total_height) / max(n_bed, 1)
    bed_losses = [_compound_heat_loss_fraction(cfg.heat_loss_frac, bed_weight_each) for _ in range(n_bed)]

    freeboard_loss = _compound_heat_loss_fraction(cfg.heat_loss_frac, freeboard_height / total_height)
    freeboard_extra = float(np.clip(cfg.freeboard_heat_loss_frac, 0.0, 1.0))
    if freeboard_extra > 0.0:
        freeboard_loss = float(1.0 - (1.0 - freeboard_loss) * (1.0 - freeboard_extra))

    return bed_losses, freeboard_loss


def _cell_solid_outflow_component(cell: Cell, comp_idx: int, *, direction: str = "total") -> float:
    if direction == "up":
        rates = cell._solid_upflow_rates()
    elif direction == "down":
        rates = cell._solid_downflow_rates()
    else:
        rates = cell._solid_outflow_rates()
    return float(np.sum(np.maximum(rates[:, comp_idx], 0.0)))


from src.solvers.cell_solver import solve_cell, evaluate_cell_state  # Re-export for test monkeypatch compatibility

def _resolve_nr_init_strategy(
    init_strategy: str | None,
    gs_warmup_steps: int | None = None,
    allow_legacy_gs: bool = False,
) -> str:
    """Resolve global-NR initialization strategy (NR-only policy).

    ``gs_warmup_steps`` / ``allow_legacy_gs`` are kept for backward compatibility
    and ignored under NR-only policy.
    """
    if init_strategy is None:
        return "vorabrechnung"
    strategy = str(init_strategy).strip().lower()
    aliases = {
        "paper_vorab": "vorabrechnung",
        "paper_vorabrechnung": "vorabrechnung",
        "vorab": "vorabrechnung",
    }
    strategy = aliases.get(strategy, strategy)
    if strategy != "vorabrechnung":
        raise ValueError(
            f"Unsupported nr_init_strategy={init_strategy!r}; "
            "expected 'vorabrechnung' (NR-only policy)"
        )
    return strategy

@dataclass
class ReactorConfig:
    """反应器配置。"""
    n_cells: int = 10             # [-]    轴向 cell 数
    H_bed: float = 5.0            # [m]    流化床高度
    H_freeboard: float = 9.5      # [m]    自由板区高度
    n_freeboard_cells: int = 0    # [-]    自由板离散段数（0=不求解自由板）
    D_bed: float = 0.6            # [m]    床层直径

    P: float = 2_500_000.0        # [Pa]   操作压力
    T_inlet: float = 300.0        # [K]    进气温度

    fuel_type: TarFuelType = "coal"

    # 固体属性
    rho_s: float = 1400.0         # [kg/m³]
    d_p: float = 0.00225          # [m]  文献范围 1.5–3.0 mm，取中值
    d_p_min: float | None = None  # [m]  可选粒径离散下限
    d_p_max: float | None = None  # [m]  可选粒径离散上限
    phi_s: float = 0.86
    eps_mf: float = 0.45
    n_age_classes: int = 10

    # 进料
    fuel_feed: float = 0.938      # [kg/s] (3377 kg/h)

    # 燃料分析
    moisture_wt: float = 16.9     # [wt%]
    ash_dry_wt: float = 11.41     # [wt%]
    VM_daf: float = 53.42         # [wt%]
    C_dry: float = 61.5           # [wt%]
    H_dry: float = 4.1            # [wt%]
    O_dry: float = 21.8           # [wt%]
    S_dry: float = 0.0            # [wt%]
    HHV_dry: float = 22.0         # [MJ/kg]
    u0_target: float | None = None
    hydrodynamics_u_d_closure: str = "wein_1992_eq312"
    hydrodynamics_bubble_diameter_model: str = "hilligardt_ode"
    hydrodynamics_psi_b_strategy: str = "wein_1992"
    hydrodynamics_lambda_strategy: str = "hamel_280"
    hydrodynamics_xi_strategy: str = "hamel_regime"
    hydrodynamics_bubble_velocity_strategy: str = "heinbockel_eq343"
    hydrodynamics_bubble_ode_strategy: str = "heinbockel_eq341"
    heat_loss_frac: float = 0.1
    freeboard_heat_loss_frac: float = 0.0  # Optional extra freeboard-only increment
    freeboard_enabled_reactions: tuple[str, ...] = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R10", "R11", "R12")
    freeboard_closure_enable_char_hetero: bool = False
    freeboard_explicit_enable_char_hetero: bool = True
    freeboard_use_shrunk_bed_top_d_p: bool = True
    freeboard_trajectory_model: str = "analytical_wirsum"
    freeboard_trajectory_coeff_model: str = "stable_mixed_drag_split"
    freeboard_u_gb_scale: float = 1.0
    freeboard_beta_a_scale: float = 1.0
    freeboard_velocity_sigma: float = 0.60
    freeboard_velocity_bins: int = 5
    freeboard_cyclone_capture_char_frac: float = 0.90
    freeboard_cyclone_capture_ash_frac: float = 0.95
    freeboard_secondary_injection_xi: float | None = None
    freeboard_secondary_O2_mol_s: float = 0.0
    freeboard_secondary_N2_mol_s: float = 0.0
    freeboard_secondary_H2O_mol_s: float = 0.0
    freeboard_secondary_T_K: float = 293.15
    freeboard_secondary_local_refine: int = 1
    freeboard_secondary_injection_mode: str = "lumped"
    nitrogen_fraction: float = 0.0
    use_gibbs_minor: bool = True
    nr_jacobian_lag_steps: int = 1
    thesis_mode: bool = False
    # Deprecated no-op (kept only for backward compatibility): GS path removed.
    allow_legacy_gs: bool = False
    # Hamel A1-style major-species Gibbs seed in Vorabrechnung x0 generation.
    vorab_major_gibbs_x0: bool = True
    # Major-gibbs init solver mode: augmented | hamel_reduced | shadow_compare
    major_gibbs_solver_mode: str = "augmented"
    # Hamel strict Vorabrechnung: drying/DAEM sources are budgeted once and fixed.
    thesis_vorab_sources_single_shot: bool = True

    # 气化剂
    ER: float | None = None
    primary_agent: str = "air_steam"
    steam_to_o2_molar: float = 0.8
    o2_steam_n2_frac_of_o2: float = 0.01

    O2_feed: float = 0.5
    H2O_feed: float = 0.3
    N2_feed: float = 0.02

    # 底部主气化剂分配到悬浮相的比例（其余进入气泡相）
    gas_inlet_dense_frac: float = 0.7

    recirculation_frac: float = 0.1
    recycle_gas: bool = True
    enable_r12: bool = True
    solid_lower_inlet_frac: float = 0.0
    top_solid_inlet_frac: float = 0.0
    allow_reactive_solid_propagation: bool = False
    reactive_solid_cutoff_xi: float = 0.35
    r4_scale: float = 1.0
    r5_scale: float = 1.0
    r6_scale: float = 1.0
    r7_scale: float = 1.0

    def __post_init__(self) -> None:
        if bool(self.thesis_mode):
            # Hamel 口径：fresh-feed 的 VM/moisture 可在下部若干 cells 继续释放，
            # 但上部与 recycle 仍由传播函数/回收函数保持隔离。
            self.allow_reactive_solid_propagation = True
            self.reactive_solid_cutoff_xi = float(np.clip(self.reactive_solid_cutoff_xi, 0.0, 1.0))
            # Thesis 主线默认启用 reduced solver（shadow compare 仅显式配置时启用）。
            if str(self.major_gibbs_solver_mode).strip().lower() == "augmented":
                self.major_gibbs_solver_mode = "hamel_reduced"
            # thesis 口径下自由板反应集必须包含 R8（与文档一致）。
            if "R8" not in self.freeboard_enabled_reactions:
                self.freeboard_enabled_reactions = tuple((*self.freeboard_enabled_reactions, "R8"))
        if self.ER is not None:
            o2, h2o, n2 = compute_gas_feeds_mol_s(
                fuel_feed_kg_s=self.fuel_feed,
                moisture_wt=self.moisture_wt,
                C_dry=self.C_dry,
                H_dry=self.H_dry,
                O_dry=self.O_dry,
                ER=float(self.ER),
                primary_agent=self.primary_agent,
                S_dry=self.S_dry,
                steam_to_o2_molar=self.steam_to_o2_molar,
                o2_steam_n2_frac_of_o2=self.o2_steam_n2_frac_of_o2,
            )
            self.O2_feed = o2
            self.H2O_feed = h2o
            self.N2_feed = n2


class Reactor:
    """BFB 反应器：多 cell 串联求解。"""

    def __init__(self, config: ReactorConfig) -> None:
        self.config = config
        if bool(self.config.thesis_mode) and ("R8" not in self.config.freeboard_enabled_reactions):
            # 兼容调用方在 ReactorConfig 构造后才打开 thesis_mode 的场景。
            self.config.freeboard_enabled_reactions = tuple((*self.config.freeboard_enabled_reactions, "R8"))
        if bool(self.config.thesis_mode):
            # Keep thesis defaults active even when thesis_mode is toggled post-construction.
            self.config.allow_reactive_solid_propagation = True
            self.config.reactive_solid_cutoff_xi = float(np.clip(self.config.reactive_solid_cutoff_xi, 0.0, 1.0))
            if str(self.config.major_gibbs_solver_mode).strip().lower() == "augmented":
                self.config.major_gibbs_solver_mode = "hamel_reduced"
        self.cells: List[Cell] = []
        self.freeboard_cells: List[Cell] = []
        self._last_explicit_freeboard_closure: dict | None = None
        self.side_cells: List[Cell] = []
        self.cyclone_cell: Cell | None = None
        self.return_leg_cell: Cell | None = None
        self._thesis_fixed_vorab_sources_ready: bool = False
        self._build_cells()
        self._apply_heat_loss_distribution()
        configure_tar_components_by_fuel(config.fuel_type)

    def _build_cells(self) -> None:
        cfg = self.config
        dh = cfg.H_bed / cfg.n_cells
        nk = max(1, int(cfg.n_age_classes))
        if nk > 1 and cfg.d_p_min is not None and cfg.d_p_max is not None:
            d_min = float(min(cfg.d_p_min, cfg.d_p_max))
            d_max = float(max(cfg.d_p_min, cfg.d_p_max))
            d_p_classes = np.linspace(d_min, d_max, nk, dtype=np.float64)
        else:
            d_p_classes = np.full(nk, float(cfg.d_p), dtype=np.float64)
        mass_fractions = np.ones(nk, dtype=np.float64) / float(nk)

        solid = SolidProps(
            rho_s=cfg.rho_s,
            d_p=cfg.d_p,
            phi_s=cfg.phi_s,
            eps_mf=cfg.eps_mf,
            n_size_classes=nk,
            d_p_classes=d_p_classes,
            mass_fractions=mass_fractions,
            moisture_wt=cfg.moisture_wt,
            ash_dry_wt=cfg.ash_dry_wt,
            VM_daf=cfg.VM_daf,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            O_dry=cfg.O_dry,
            HHV_dry_MJ_kg=cfg.HHV_dry,
            nitrogen_fraction=cfg.nitrogen_fraction,
        )

        for i in range(cfg.n_cells):
            geo = CellGeometry(D_bed=cfg.D_bed, dh=dh, h_center=(i + 0.5) * dh)
            cell = Cell(geo=geo, solid=solid, fuel_type=cfg.fuel_type)
            cell.P = cfg.P
            cell.u0_target = cfg.u0_target
            cell.use_gibbs_minor = bool(cfg.use_gibbs_minor)
            cell.u_d_closure = str(cfg.hydrodynamics_u_d_closure)
            cell.bubble_diameter_model = str(cfg.hydrodynamics_bubble_diameter_model)
            cell.psi_b_strategy = str(cfg.hydrodynamics_psi_b_strategy)
            cell.lambda_strategy = str(cfg.hydrodynamics_lambda_strategy)
            cell.xi_strategy = str(cfg.hydrodynamics_xi_strategy)
            cell.bubble_velocity_strategy = str(cfg.hydrodynamics_bubble_velocity_strategy)
            cell.bubble_ode_strategy = str(cfg.hydrodynamics_bubble_ode_strategy)
            if bool(cfg.thesis_mode):
                # Thesis mode enforces the Hamel-aligned hydrodynamics closure chain.
                cell.u_d_closure = "wein_1992_eq312"
                cell.bubble_diameter_model = "hilligardt_ode"
                cell.psi_b_strategy = "wein_1992"
                cell.lambda_strategy = "hamel_280"
                cell.xi_strategy = "hamel_regime"
                cell.bubble_velocity_strategy = "heinbockel_eq343"
                cell.bubble_ode_strategy = "heinbockel_eq341"
                cell.solid_state_model = "holdup_transport"
            cell.enable_r12 = bool(cfg.enable_r12)
            cell.r4_scale = float(max(cfg.r4_scale, 0.0))
            cell.r5_scale = float(max(cfg.r5_scale, 0.0))
            cell.r6_scale = float(max(cfg.r6_scale, 0.0))
            cell.r7_scale = float(max(cfg.r7_scale, 0.0))
            self.cells.append(cell)
        if self._use_explicit_side_block_cells():
            side_geo = CellGeometry(D_bed=cfg.D_bed, dh=dh, h_center=cfg.H_bed + 0.5 * dh)
            cyclone = Cell(geo=side_geo, solid=solid, fuel_type=cfg.fuel_type)
            cyclone.P = cfg.P
            cyclone.cell_type = "cyclone"
            cyclone.solid_state_model = "holdup_transport"
            cyclone.use_gibbs_minor = bool(cfg.use_gibbs_minor)
            cyclone.enable_r12 = bool(cfg.enable_r12)
            cyclone.r4_scale = float(max(cfg.r4_scale, 0.0))
            cyclone.r5_scale = float(max(cfg.r5_scale, 0.0))
            cyclone.r6_scale = float(max(cfg.r6_scale, 0.0))
            cyclone.r7_scale = float(max(cfg.r7_scale, 0.0))

            return_leg_geo = CellGeometry(D_bed=cfg.D_bed, dh=dh, h_center=cfg.H_bed + 1.5 * dh)
            return_leg = Cell(geo=return_leg_geo, solid=solid, fuel_type=cfg.fuel_type)
            return_leg.P = cfg.P
            return_leg.cell_type = "return_leg"
            return_leg.solid_state_model = "holdup_transport"
            return_leg.use_gibbs_minor = bool(cfg.use_gibbs_minor)
            return_leg.enable_r12 = bool(cfg.enable_r12)
            return_leg.r4_scale = float(max(cfg.r4_scale, 0.0))
            return_leg.r5_scale = float(max(cfg.r5_scale, 0.0))
            return_leg.r6_scale = float(max(cfg.r6_scale, 0.0))
            return_leg.r7_scale = float(max(cfg.r7_scale, 0.0))

            self.cyclone_cell = cyclone
            self.return_leg_cell = return_leg
            self.side_cells = [cyclone, return_leg]
        if self._use_explicit_freeboard_cells():
            dh_fb = float(cfg.H_freeboard) / max(int(cfg.n_freeboard_cells), 1)
            for i in range(int(cfg.n_freeboard_cells)):
                geo_fb = CellGeometry(
                    D_bed=cfg.D_bed,
                    dh=dh_fb,
                    h_center=float(cfg.H_bed) + (i + 0.5) * dh_fb,
                )
                fb_cell = Cell(geo=geo_fb, solid=solid, fuel_type=cfg.fuel_type)
                fb_cell.P = cfg.P
                fb_cell.cell_type = "freeboard"
                fb_cell.use_gibbs_minor = bool(cfg.use_gibbs_minor)
                fb_cell.enable_r12 = bool(cfg.enable_r12)
                fb_cell.r4_scale = float(max(cfg.r4_scale, 0.0))
                fb_cell.r5_scale = float(max(cfg.r5_scale, 0.0))
                fb_cell.r6_scale = float(max(cfg.r6_scale, 0.0))
                fb_cell.r7_scale = float(max(cfg.r7_scale, 0.0))
                fb_cell.freeboard_u_gb_scale = float(max(cfg.freeboard_u_gb_scale, 1e-6))
                fb_cell.freeboard_explicit_char_hetero_enabled = bool(cfg.freeboard_explicit_enable_char_hetero)
                fb_cell.solid_state_model = "freeboard_closure"
                self.freeboard_cells.append(fb_cell)

    def _use_explicit_side_block_cells(self) -> bool:
        return bool(self.config.thesis_mode)

    def _use_explicit_freeboard_cells(self) -> bool:
        return bool(self.config.thesis_mode) and float(self.config.H_freeboard) > 0.0 and int(self.config.n_freeboard_cells) > 0

    def _use_explicit_freeboard_solver_graph(self) -> bool:
        return self._use_explicit_freeboard_cells()

    def _solver_cells_for_nr(self) -> List[Cell]:
        if self._use_explicit_side_block_cells():
            if self._use_explicit_freeboard_solver_graph():
                return self.cells + self.freeboard_cells + self.side_cells
            if not self._use_explicit_freeboard_cells():
                return self.cells + self.side_cells
        return self.cells

    def _use_side_blocks_in_nr_boundary_path(self) -> bool:
        return self._use_explicit_side_block_cells() and (
            not self._use_explicit_freeboard_cells() or self._use_explicit_freeboard_solver_graph()
        )

    def _default_nr_jacobian_strategy(self) -> str:
        """Default Jacobian strategy selected from the active solver graph topology."""
        if self._use_side_blocks_in_nr_boundary_path():
            return "band_plus_side_elements_structured"
        return "block_tridiag_structured"

    def _initialize_explicit_side_block_states(self) -> None:
        """委托 ``side_block_bridge``：主链顶端状态 -> side-block 初始化。"""
        _init_side_blocks_bridge(self)

    def _affected_nr_residual_cells(self, changed_cell_idx: int) -> tuple[int, ...]:
        return _affected_nr_cells_graph(
            changed_cell_idx=changed_cell_idx,
            n_bed=len(self.cells),
            n_freeboard=len(self.freeboard_cells),
            explicit_freeboard_graph=self._use_explicit_freeboard_solver_graph(),
            side_blocks_in_boundary_path=self._use_side_blocks_in_nr_boundary_path(),
            has_side_blocks=bool(
                self._use_explicit_side_block_cells()
                and self.cyclone_cell is not None
                and self.return_leg_cell is not None
            ),
        )

    def _sync_freeboard_cells_from_closure(self, fb: dict) -> None:
        """委托 ``freeboard_bridge``：closure → freeboard_cells 状态。"""
        _sync_fb_cells_bridge(self, fb)

    def _refresh_explicit_freeboard_transport_from_closure(self) -> dict | None:
        """委托 ``freeboard_bridge``：床顶 → simulate_freeboard → 同步显式自由板链。"""
        return _refresh_fb_transport_bridge(self)

    def _apply_heat_loss_distribution(self) -> None:
        """Refresh per-bed-cell heat-loss fractions from reactor-level settings."""
        bed_losses, _ = _resolve_axial_heat_loss_distribution(self.config)
        for cell, loss_frac in zip(self.cells, bed_losses):
            cell.heat_loss_frac = float(loss_frac)

    def _compute_carbon_conversion(self) -> float:
        from src.solvers.result_builder import compute_carbon_conversion

        return compute_carbon_conversion(self, _cell_solid_outflow_component)

    def _build_exit_summary(self) -> dict:
        from src.solvers.result_builder import build_exit_summary

        return build_exit_summary(
            self,
            resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
            cell_solid_outflow_component_fn=_cell_solid_outflow_component,
        )

    def _set_bottom_cell_feeds(self) -> None:
        _set_bottom_cell_feeds_exec(self.cells, self.config)

    def _apply_bottom_recycle(self, relax: float | None = None) -> None:
        if self._use_side_blocks_in_nr_boundary_path() and self.return_leg_cell is not None:
            _apply_all_nr_boundary_data_exec(
                self.cells,
                self.config,
                cyclone_cell=self.cyclone_cell,
                return_leg_cell=self.return_leg_cell,
            )
            return
        _apply_bottom_recycle_exec(self.cells, self.config, relax)

    def _propagate_upstream(self, i: int) -> None:
        _propagate_upstream_exec(self.cells, self.config, i)

    def _apply_all_bc_for_nr(self) -> None:
        _apply_all_nr_boundary_data_exec(
            self.cells,
            self.config,
            freeboard_cells=self.freeboard_cells if self._use_explicit_freeboard_solver_graph() else None,
            cyclone_cell=self.cyclone_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            return_leg_cell=self.return_leg_cell if self._use_side_blocks_in_nr_boundary_path() else None,
        )

    def _apply_local_bc_for_nr(self, changed_cell_idx: int) -> None:
        _apply_local_nr_boundary_data_exec(
            self.cells,
            self.config,
            changed_cell_idx,
            freeboard_cells=self.freeboard_cells if self._use_explicit_freeboard_solver_graph() else None,
            cyclone_cell=self.cyclone_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            return_leg_cell=self.return_leg_cell if self._use_side_blocks_in_nr_boundary_path() else None,
        )

    def _refresh_vorabrechnung_sources_for_nr(self, force: bool = False) -> None:
        """Outer-loop Vorabrechnung refresh for thesis-aligned global NR."""
        from src.solvers.vorabrechnung import refresh_vorabrechnung_for_cells

        # Explicit freeboard cells are part of the NR unknown vector in thesis mode.
        # Re-syncing them from closure during outer refresh overwrites NR-updated states
        # and can trigger large residual jumps between outer iterations.
        if not self._use_explicit_freeboard_solver_graph():
            self._refresh_explicit_freeboard_transport_from_closure()
        # Side-block states are seeded once during init/precalc.
        # Re-seeding each outer loop overwrites NR-updated side states and breaks Check2.
        self._apply_all_bc_for_nr()
        strict_single_shot = bool(self.config.thesis_mode) and bool(self.config.thesis_vorab_sources_single_shot)
        refresh_vorabrechnung_for_cells(
            self._solver_cells_for_nr(),
            force=force and (not strict_single_shot),
            refresh_sources=(not strict_single_shot),
        )

    def _initialize_nr_hydrodynamics_freeze_cache(self) -> None:
        """Initialize frozen hydrodynamics cache from the current initialized state.

        This makes hydrodynamics explicitly part of initialization (x0 stage),
        and ensures inner NR always starts from a prepared frozen snapshot.
        """
        for cell in self._solver_cells_for_nr():
            cell.calc_hydrodynamics()
            cell._snapshot_vorabrechnung_hydrodynamics()

    def _initialize_thesis_single_shot_vorab_sources(self) -> None:
        """Initialize fixed drying/DAEM source caches once using bed-mean temperature."""
        if not (bool(self.config.thesis_mode) and bool(self.config.thesis_vorab_sources_single_shot)):
            return
        if self._thesis_fixed_vorab_sources_ready:
            return
        from src.solvers.vorabrechnung import initialize_fixed_vorabrechnung_sources_for_cells

        T_bed_avg = float(np.mean([float(c.T) for c in self.cells])) if self.cells else float(self.config.T_inlet)
        initialize_fixed_vorabrechnung_sources_for_cells(
            self.cells,
            T_reference=T_bed_avg,
        )
        self._thesis_fixed_vorab_sources_ready = True

    def _vorabrechnung_snapshot_signature(self) -> str:
        """Build a lightweight deterministic signature for current Vorabrechnung state."""
        payload: list[float] = []
        vm_valid_count = 0
        for cell in self._solver_cells_for_nr():
            payload.extend(
                [
                    float(cell.T),
                    float(cell.u_mf),
                    float(cell.u_d),
                    float(cell.u_b),
                    float(cell.d_b),
                    float(cell.eps_b),
                    float(cell.eps_d_voidage),
                    float(cell.K_bd),
                ]
            )
            if bool(getattr(cell, "_vm_cache_valid", False)):
                vm_valid_count += 1
                payload.append(float(np.sum(getattr(cell, "_vm_gas_source_cache", 0.0))))
                payload.append(float(np.sum(getattr(cell, "_vm_solid_sink_cache", 0.0))))
            else:
                payload.extend([0.0, 0.0])
        arr = np.asarray(payload, dtype=np.float64)
        if arr.size == 0:
            return "vorab:v1:empty"
        s1 = float(np.sum(arr))
        s2 = float(np.dot(arr, arr))
        return f"vorab:v1:n={arr.size};vm={vm_valid_count};sum={s1:.8e};sumsq={s2:.8e}"

    def _evaluate_current_gs_state(self) -> tuple[float, float]:
        """评估当前全堆 GS 状态的残差（不更新状态）。"""
        self._set_bottom_cell_feeds()
        self._apply_all_bc_for_nr()
        
        res_list = []
        rms_list = []
        for cell in self.cells:
            metrics = evaluate_cell_state(cell)
            res_list.append(metrics["residual"])
            rms_list.append(metrics["rms_scaled"])
            
        return float(np.max(res_list)), float(np.mean(rms_list))

    def solve(
        self,
        max_global_iter: int = 20,
        tol_global: float = 1e-4,
        solver: str = "global_nr",
        verbose: bool = False,
        nr_init_strategy: str | None = None,
        nr_jacobian_strategy: str | None = None,
        nr_linear_solver_backend: str | None = None,
        nr_jacobian_lag_steps: int | None = None,
    ) -> dict:
        """稳态求解入口。

        与 Hamel (1999) 论文一致，仅支持全域联立 ``solver="global_nr"``。
        """
        self._apply_heat_loss_distribution()
        solver = str(solver).strip().lower()
        if solver != "global_nr":
            raise ValueError(
                "Only solver='global_nr' is supported (NR-only policy aligned with Hamel)."
            )
        return self._solve_global_nr(
            max_iter=max_global_iter,
            tol=tol_global,
            verbose=verbose,
            init_strategy=nr_init_strategy,
            jacobian_strategy=nr_jacobian_strategy,
            linear_solver_backend=nr_linear_solver_backend,
            jacobian_lag_steps=nr_jacobian_lag_steps,
        )

    def _solve_gauss_seidel(
        self,
        max_global_iter: int,
        tol_global: float,
        verbose: bool = False,
    ) -> dict:
        """Legacy Gauss-Seidel solver path with compatibility hooks."""
        if not bool(self.config.allow_legacy_gs):
            raise RuntimeError(
                "Legacy Gauss-Seidel path is explicitly disabled in ReactorConfig. "
                "Use allow_legacy_gs=True to enable for debugging."
            )
        from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0

        # Warm-start reference profile used by bottom-cell temperature cap logic.
        T_ref = estimate_axial_T_profile(
            n_cells=self.config.n_cells,
            T_inlet=self.config.T_inlet,
            O2_feed=self.config.O2_feed,
            H2O_feed=self.config.H2O_feed,
            N2_feed=self.config.N2_feed,
            fuel_feed_kg_s=self.config.fuel_feed,
            C_dry=self.config.C_dry,
            H_dry=self.config.H_dry,
            moisture_wt=self.config.moisture_wt,
            P=self.config.P,
        )
        generate_initial_x0(
            cells=self.cells,
            O2_feed=self.config.O2_feed,
            H2O_feed=self.config.H2O_feed,
            N2_feed=self.config.N2_feed,
            fuel_feed_kg_s=self.config.fuel_feed,
            C_dry=self.config.C_dry,
            H_dry=self.config.H_dry,
            O_dry=self.config.O_dry,
            moisture_wt=self.config.moisture_wt,
            ash_dry_wt=self.config.ash_dry_wt,
            VM_daf=self.config.VM_daf,
            T_profile=T_ref,
            fuel_type=self.config.fuel_type,
            major_gibbs_solver_mode=str(self.config.major_gibbs_solver_mode),
        )

        history: list[dict] = []
        best_state = _capture_cells_state(self.cells)
        _, best_rms = self._evaluate_current_gs_state()
        best_profile = _gs_profile_metrics(self.cells, o2_feed=self.config.O2_feed)
        best_score = float(best_rms + best_profile["penalty"])
        best_iter = 0
        restored_best_state = False
        residual_gs = 0.0
        rms_scaled_gs = best_rms

        for g_iter in range(max_global_iter):
            self._set_bottom_cell_feeds()
            cell_traces: list[dict] = []
            for i, cell in enumerate(self.cells):
                self._propagate_upstream(i)
                # Fast stiff step (used in GS legacy loop).
                fast_res = solve_cell(
                    cell,
                    stiff_stabilization=True,
                    skip_homotopy=True,
                    skip_multistart=True,
                )
                full_caps: list[float | None] = []
                if not bool(fast_res.get("physically_converged", False)):
                    candidate_caps: list[float | None] = [float(T_ref[i] + 75.0)] if i == 0 else [None]
                    for cap in candidate_caps:
                        if cap in full_caps:
                            continue
                        full_caps.append(cap)
                        full_res = solve_cell(
                            cell,
                            stiff_stabilization=True,
                            skip_homotopy=False,
                            skip_multistart=False,
                            temperature_cap_K=cap,
                        )
                        if bool(full_res.get("physically_converged", False)):
                            break

                _clamp_soft_rollback_monotonic(cell)
                cell_traces.append(
                    {
                        "cell_index": i,
                        "full_temperature_caps_K": [float(c) for c in full_caps if c is not None],
                        "full_temperature_cap_K": next((float(c) for c in full_caps if c is not None), None),
                    }
                )

            # Optional upper-pair corrective sweep for high-RMS top cells.
            cell_metrics = [evaluate_cell_state(c) for c in self.cells]
            upper_pair = sorted(
                range(len(self.cells)),
                key=lambda idx: float(cell_metrics[idx].get("rms_scaled", 0.0)),
                reverse=True,
            )[:2]
            upper_before = max((float(cell_metrics[idx].get("rms_scaled", 0.0)) for idx in upper_pair), default=0.0)
            upper_correction = {"attempted": False, "before_rms": float(upper_before)}
            if upper_before >= 0.2 and len(upper_pair) == 2:
                upper_correction["attempted"] = True
                for idx in sorted(upper_pair):
                    solve_cell(
                        self.cells[idx],
                        stiff_stabilization=True,
                        skip_homotopy=True,
                        skip_multistart=True,
                    )
                    _clamp_soft_rollback_monotonic(self.cells[idx])

            residual_gs, rms_scaled_gs = self._evaluate_current_gs_state()
            profile = _gs_profile_metrics(
                self.cells,
                o2_feed=self.config.O2_feed,
                T_ref=np.asarray(T_ref, dtype=np.float64),
                m_char_in=float(np.sum(np.maximum(self.cells[0].m_solid_in[:, S_CHAR] + self.cells[0].m_solid_zu[:, S_CHAR], 0.0))),
            )
            score = float(rms_scaled_gs + profile["penalty"])
            if score < best_score:
                best_score = score
                best_iter = g_iter + 1
                best_state = _capture_cells_state(self.cells)
                best_profile = profile
            elif g_iter >= 1:
                history.append(
                    {
                        "iter": g_iter + 1,
                        "residual_gs": float(residual_gs),
                        "rms_scaled_gs": float(rms_scaled_gs),
                        "profile_metrics": profile,
                        "cells": cell_traces,
                        "upper_pair_correction": upper_correction,
                    }
                )
                _restore_cells_state(self.cells, best_state)
                restored_best_state = True
                break

            history.append(
                {
                    "iter": g_iter + 1,
                    "residual_gs": float(residual_gs),
                    "rms_scaled_gs": float(rms_scaled_gs),
                    "profile_metrics": profile,
                    "cells": cell_traces,
                    "upper_pair_correction": upper_correction,
                }
            )
            if float(rms_scaled_gs) <= float(tol_global):
                break

        out = self._build_exit_summary()
        out.update(
            {
                "converged": bool(rms_scaled_gs <= tol_global),
                "n_iter": len(history),
                "best_iter": int(best_iter if best_iter > 0 else 1),
                "restored_best_state": bool(restored_best_state),
                "residual_gs": float(residual_gs),
                "rms_scaled_gs": float(rms_scaled_gs),
                "history": history,
                "best_profile_metrics": best_profile,
            }
        )
        return out

    def _solve_global_nr(
        self,
        max_iter: int,
        tol: float,
        verbose: bool = False,
        init_strategy: str | None = None,
        gs_warmup_steps: int | None = None,
        jacobian_strategy: str | None = None,
        linear_solver_backend: str | None = None,
        jacobian_lag_steps: int | None = None,
    ) -> dict:
        from src.solvers.convergence import InnerConvergence
        from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
        from src.workflow.steps.nr_inner_step import build_global_nr_inner_solve_fn
        from src.workflow.steps.outer_abgleich_step import run_outer_abgleich_for_global_nr

        cfg = self.config
        if gs_warmup_steps is not None and int(max(gs_warmup_steps, 0)) > 0:
            raise ValueError("gs_warmup_steps is not supported under NR-only policy.")
        self._thesis_fixed_vorab_sources_ready = False
        solver_cells = self._solver_cells_for_nr()
        solve_started = perf_counter()

        precalc = run_init_and_precalc_for_global_nr(
            self,
            init_strategy=init_strategy,
            gs_warmup_steps=gs_warmup_steps,
        )
        resolved_init_strategy = precalc.resolved_init_strategy
        nr_init_s_total = precalc.nr_init_s_total
        nr_vorabrechnung_s = precalc.nr_vorabrechnung_s
        if bool(cfg.thesis_mode) and bool(cfg.thesis_vorab_sources_single_shot):
            # Thesis strict path keeps drying/DAEM fixed after init; allowing too many
            # outer slots over-fragments inner-NR budget and hurts Check1 progress.
            outer_max = 2
        else:
            outer_max = max(2, max_iter // 4)
        tol_rms = float(np.clip(0.01 * max(tol, 1.0), 0.005, 0.02))
        total_inner_budget = max(int(max_iter), 1)
        jacobian_mode = (
            self._default_nr_jacobian_strategy()
            if jacobian_strategy is None
            else str(jacobian_strategy)
        )
        jacobian_lag = int(
            max(
                1,
                cfg.nr_jacobian_lag_steps if jacobian_lag_steps is None else jacobian_lag_steps,
            )
        )
        # Check1：与历史 ``norm_F < tol_rms`` 等价（rtol=0 时阈值为 atol）
        inner_check1 = InnerConvergence(atol=float(tol_rms), rtol=0.0)
        _solve_inner = build_global_nr_inner_solve_fn(
            self,
            solver_cells,
            cfg,
            tol_rms,
            inner_check1,
            verbose,
        )

        outer_result = run_outer_abgleich_for_global_nr(
            cells=solver_cells,
            outer_max=outer_max,
            total_inner_budget=total_inner_budget,
            tol=float(tol),
            jacobian_mode=jacobian_mode,
            jacobian_lag=jacobian_lag,
            tol_rms=tol_rms,
            refresh_fn=self._refresh_vorabrechnung_sources_for_nr,
            snapshot_signature_fn=self._vorabrechnung_snapshot_signature,
            solve_inner_fn=_solve_inner,
        )

        last_nr_result = outer_result.last_nr_result
        outer_converged = outer_result.converged_outer
        outer_iters = outer_result.outer_iters
        outer_history = outer_result.outer_history
        all_nr_history = outer_result.all_nr_history
        all_nr_lambdas = outer_result.all_nr_lambdas
        all_nr_line_search_trials = outer_result.all_nr_line_search_trials
        all_nr_clip_history = outer_result.all_nr_clip_history
        agg_nr_timing = outer_result.agg_nr_timing
        agg_nr_counts = outer_result.agg_nr_counts
        outer_refresh_s_total = outer_result.outer_refresh_s_total
        inner_solve_s_total = outer_result.inner_solve_s_total
        used_inner_budget = outer_result.used_inner_budget
        return self._finalize_global_nr_result(
            dict(
                last_nr_result,
                converged_outer=outer_converged,
                converged_inner_nr=bool(last_nr_result.get("converged")),
                converged=outer_converged,
                norm_history=all_nr_history,
                accepted_lambda_history=all_nr_lambdas,
                line_search_trial_counts=all_nr_line_search_trials,
                clip_history=all_nr_clip_history,
                n_iter=len(all_nr_history),
                nr_init_strategy=resolved_init_strategy,
                nr_init_s_total=nr_init_s_total,
                nr_vorabrechnung_s=nr_vorabrechnung_s,
                nr_vorabrechnung_policy=(
                    "single_shot_sources_outer_refresh_hydrodynamics_fixed_inner"
                    if (bool(cfg.thesis_mode) and bool(cfg.thesis_vorab_sources_single_shot))
                    else "outer_refresh_fixed_inner_sources"
                ),
                nr_gs_warmup_steps=0,
                nr_outer_max=outer_max,
                nr_outer_iters=outer_iters,
                nr_inner_tol_rms=tol_rms,
                nr_jacobian_strategy=jacobian_mode,
                nr_jacobian_strategy_requested=("auto" if jacobian_strategy is None else str(jacobian_strategy)),
                nr_linear_solver_backend=(
                    str(linear_solver_backend)
                    if linear_solver_backend is not None
                    else (
                        "structured_direct"
                        if str(jacobian_mode) in {"block_tridiag_structured", "band_plus_side_elements_structured"}
                        else "sparse_direct_fallback"
                    )
                ),
                nr_jacobian_lag_steps=jacobian_lag,
                nr_gs_warmup_s=0.0,
                nr_outer_refresh_s_total=outer_refresh_s_total,
                nr_inner_solve_s_total=inner_solve_s_total,
                nr_outer_history=outer_history,
                nr_vorabrechnung_signatures=[str(item.get("vorabrechnung_signature", "")) for item in outer_history],
                nr_inner_budget_total=total_inner_budget,
                nr_inner_budget_used=used_inner_budget,
                nr_total_s=perf_counter() - solve_started,
                timing=agg_nr_timing or last_nr_result.get("timing"),
                counts=agg_nr_counts or last_nr_result.get("counts"),
            )
        )

    def _finalize_global_nr_result(self, nr_result: dict) -> dict:
        from src.solvers.result_builder import finalize_global_nr_result

        return finalize_global_nr_result(
            self,
            nr_result,
            resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
            cell_solid_outflow_component_fn=_cell_solid_outflow_component,
        )
