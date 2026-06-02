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
    align_bottom_primary_gas_state_to_inlet_split as _align_bottom_primary_gas_state_exec,
    apply_all_nr_boundary_data as _apply_all_nr_boundary_data_exec,
    apply_bottom_recycle as _apply_bottom_recycle_exec,
    apply_local_nr_boundary_data as _apply_local_nr_boundary_data_exec,
    propagate_upstream as _propagate_upstream_exec,
    propagated_solid_stream as _propagated_solid_stream,
    recycled_solid_stream as _recycled_solid_stream,
    seed_bed_holdup_chain_from_transport as _seed_bed_holdup_chain_from_transport_exec,
    set_bottom_cell_feeds as _set_bottom_cell_feeds_exec,
    snapshot_bottom_gas_inlet_split_from_vorabrechnung as _snapshot_bottom_gas_split_exec,
)
from src.core.connectivity_graph import affected_nr_residual_cells as _affected_nr_cells_graph
from src.core.freeboard_bridge import (
    refresh_explicit_freeboard_transport_from_closure as _refresh_fb_transport_bridge,
    sync_freeboard_hydrodynamic_bridge_from_bed_top as _sync_fb_hydro_bridge,
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
) -> str:
    """Resolve global-NR initialization strategy (NR-only policy)."""
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
    bed_dh_profile: tuple[float, ...] | None = None  # [m] Optional per-cell bed heights; sum must equal H_bed
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
    freeboard_age_quadrature_bins: int = 1
    freeboard_age_quadrature_max_age: float = 0.98
    freeboard_cyclone_capture_char_frac: float = 0.90
    freeboard_cyclone_capture_ash_frac: float = 0.95
    freeboard_secondary_injection_xi: float | None = None
    freeboard_secondary_O2_mol_s: float = 0.0
    freeboard_secondary_N2_mol_s: float = 0.0
    freeboard_secondary_H2O_mol_s: float = 0.0
    freeboard_secondary_T_K: float = 293.15
    freeboard_secondary_local_refine: int = 1
    freeboard_secondary_injection_mode: str = "lumped"
    explicit_freeboard_solver_graph_enabled: bool = True
    nitrogen_fraction: float = 0.0
    use_gibbs_minor: bool = True
    nr_jacobian_lag_steps: int = 1
    nr_damping_halvings_thesis: int = 14
    nr_lambda_min_thesis: float = 1.0 / 4096.0
    nr_prefer_full_step_thesis: bool = True
    nr_step_model_thesis: str = "newton"  # newton | lm | ptc | equil_newton
    nr_lm_mu0_thesis: float = 1.0e-4
    nr_lm_mu_growth_thesis: float = 10.0
    nr_ptc_alpha0_thesis: float = 1.0
    nr_ptc_alpha_growth_thesis: float = 2.0
    nr_equil_iters_thesis: int = 3
    nr_equil_scale_clip_thesis: float = 1.0e3
    nr_line_search_max_trials_thesis: int = 0  # 0 => use n_damp_halvings
    nr_nonmonotone_enabled_thesis: bool = False
    nr_nonmonotone_window_thesis: int = 5
    nr_nonmonotone_relax_thesis: float = 1.02
    nr_strict_check1_before_refresh_thesis: bool = True
    nr_reaction_continuation_stages_thesis: tuple[float, ...] = (1.0,)
    nr_reaction_continuation_height_m_thesis: float = 0.0
    nr_reaction_continuation_inner_iter_cap_thesis: int = 0
    nr_temperature_fence_enabled_thesis: bool = True
    nr_temperature_fence_lower_margin_K_thesis: float = 150.0
    nr_temperature_fence_upper_margin_K_thesis: float = 350.0
    thesis_mode: bool = False
    explicit_side_blocks_enabled: bool = True
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

    # 底部主气化剂分配到悬浮相的比例（其余进入气泡相）。
    # ``gas_inlet_split_strategy="fixed"`` 时作为直接比例；Hamel-aligned
    # ``precalc_hydrodynamic_flux`` 时仅作为 Vorabrechnung 尚未可用的 fallback。
    gas_inlet_dense_frac: float = 0.7
    gas_inlet_split_strategy: str = "fixed"

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
        if cfg.n_cells <= 0:
            raise ValueError(f"n_cells must be positive, got {cfg.n_cells}")
        dh_profile_cfg = cfg.bed_dh_profile
        if dh_profile_cfg is None:
            bed_dh = np.full(int(cfg.n_cells), float(cfg.H_bed) / float(cfg.n_cells), dtype=np.float64)
        else:
            bed_dh = np.asarray(tuple(float(v) for v in dh_profile_cfg), dtype=np.float64)
            if bed_dh.size != int(cfg.n_cells):
                raise ValueError(
                    f"bed_dh_profile length ({bed_dh.size}) must equal n_cells ({int(cfg.n_cells)})"
                )
            if np.any(bed_dh <= 0.0):
                raise ValueError("bed_dh_profile entries must be strictly positive")
            if not np.isclose(float(np.sum(bed_dh)), float(cfg.H_bed), rtol=1e-8, atol=1e-10):
                raise ValueError(
                    f"sum(bed_dh_profile)={float(np.sum(bed_dh)):.12g} must equal H_bed={float(cfg.H_bed):.12g}"
                )
        z0 = np.concatenate(([0.0], np.cumsum(bed_dh)[:-1]))
        zc = z0 + 0.5 * bed_dh
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
            geo = CellGeometry(D_bed=cfg.D_bed, dh=float(bed_dh[i]), h_center=float(zc[i]))
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
            cell._n_vorab_cells = cfg.n_cells  # legacy uniform-mesh tau fallback
            cell._vorab_bed_height = float(cfg.H_bed)
            self.cells.append(cell)
        top_dh = float(bed_dh[-1])
        if self._use_explicit_side_block_cells():
            side_geo = CellGeometry(D_bed=cfg.D_bed, dh=top_dh, h_center=cfg.H_bed + 0.5 * top_dh)
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

            return_leg_geo = CellGeometry(D_bed=cfg.D_bed, dh=top_dh, h_center=cfg.H_bed + 1.5 * top_dh)
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
        return bool(self.config.thesis_mode) and bool(self.config.explicit_side_blocks_enabled)

    def _use_explicit_freeboard_cells(self) -> bool:
        return bool(self.config.thesis_mode) and float(self.config.H_freeboard) > 0.0 and int(self.config.n_freeboard_cells) > 0

    def _use_explicit_freeboard_solver_graph(self) -> bool:
        return self._use_explicit_freeboard_cells() and bool(self.config.explicit_freeboard_solver_graph_enabled)

    def _solver_cells_for_nr(self) -> List[Cell]:
        if self._use_explicit_side_block_cells():
            if self._use_explicit_freeboard_solver_graph():
                return self.cells + self.freeboard_cells + self.side_cells
            return self.cells + self.side_cells
        return self.cells

    def _resolve_nr_reaction_continuation_stages(self) -> tuple[float, ...]:
        """Return monotone reaction-source continuation stages for thesis NR.

        The stages are a path-control device only.  The final stage is always
        forced to 1.0 so the reported residuals remain the full physical model.
        """
        if not bool(self.config.thesis_mode):
            return (1.0,)
        raw = getattr(self.config, "nr_reaction_continuation_stages_thesis", (1.0,))
        try:
            stages = tuple(float(x) for x in raw)
        except TypeError:
            stages = (float(raw),)
        stages = tuple(float(np.clip(x, 0.0, 1.0)) for x in stages if np.isfinite(float(x)))
        if not stages:
            stages = (1.0,)
        if stages[-1] < 1.0:
            stages = stages + (1.0,)
        return stages

    def _nr_reaction_continuation_target_cells(self, solver_cells: List[Cell]) -> list[Cell]:
        height = float(max(getattr(self.config, "nr_reaction_continuation_height_m_thesis", 0.0), 0.0))
        if height <= 0.0:
            return []
        return [
            c
            for c in solver_cells
            if str(getattr(c, "cell_type", "bed")) == "bed" and float(c.geo.h_center) <= height
        ]

    def _use_side_blocks_in_nr_boundary_path(self) -> bool:
        # Side blocks participate in the NR boundary path whenever they are part of
        # the active thesis solver graph, including the "closure-owned freeboard"
        # mode where freeboard cells are excluded from NR unknowns.
        return self._use_explicit_side_block_cells()

    def _default_nr_jacobian_strategy(self) -> str:
        """Default Jacobian strategy selected from the active solver graph topology.

        Recycle gas (top-cell → bottom-cell coupling) creates a non-tridiagonal
        Jacobian structure.  The structured path declares that recycle closure as a
        side block, avoiding the much slower dense finite-difference fallback.
        """
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

    def _sync_freeboard_cells_from_closure(self, fb: dict, *, preserve_gas_state: bool = False) -> None:
        """委托 ``freeboard_bridge``：closure → freeboard_cells 状态。"""
        _sync_fb_cells_bridge(self, fb, preserve_gas_state=preserve_gas_state)

    def _refresh_explicit_freeboard_transport_from_closure(self, *, preserve_gas_state: bool = False) -> dict | None:
        """委托 ``freeboard_bridge``：床顶 → simulate_freeboard → 同步显式自由板链。"""
        return _refresh_fb_transport_bridge(self, preserve_gas_state=preserve_gas_state)

    def _sync_freeboard_hydrodynamic_bridge_from_bed_top(self) -> None:
        """委托 ``freeboard_bridge``：只同步床顶 ghost-bubble 边界。"""
        _sync_fb_hydro_bridge(self)

    def _seed_initialized_holdup_chain_for_nr(self) -> None:
        """Seed bed holdup states once from initialized thesis transport inflows.

        After x0 switched from total-bed-holdup semantics to active-solid inflow
        semantics, upper bed cells no longer receive a resident solid inventory
        automatically. A one-time bottom→top sweep restores the intended Eq. 2.6
        transport chain before explicit freeboard/side-block initialization.
        """
        top_above_cell = self.freeboard_cells[0] if self.freeboard_cells else None
        _seed_bed_holdup_chain_from_transport_exec(self.cells, top_above_cell=top_above_cell)

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

    def _snapshot_bottom_gas_inlet_split_from_vorabrechnung(self) -> None:
        _snapshot_bottom_gas_split_exec(self.cells, self.config)
        _align_bottom_primary_gas_state_exec(self.cells, self.config)

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
        strict_single_shot = bool(self.config.thesis_mode) and bool(self.config.thesis_vorab_sources_single_shot)
        refresh_sources = not strict_single_shot
        force_sources = force and (not strict_single_shot)
        self._apply_all_bc_for_nr()
        if self._use_explicit_freeboard_solver_graph():
            # Refresh bed hydrodynamics first, then update the Hamel freeboard
            # ghost-bubble bridge from the current bed-surface state before
            # snapshotting freeboard hydrodynamics.  This keeps explicit
            # freeboard gas/energy NR state intact.
            refresh_vorabrechnung_for_cells(
                self.cells,
                force=force_sources,
                refresh_sources=refresh_sources,
            )
            self._snapshot_bottom_gas_inlet_split_from_vorabrechnung()
            self._sync_freeboard_hydrodynamic_bridge_from_bed_top()
            trailing_cells = self.freeboard_cells + (self.side_cells if self._use_explicit_side_block_cells() else [])
            refresh_vorabrechnung_for_cells(
                trailing_cells,
                force=force_sources,
                refresh_sources=refresh_sources,
            )
            self._apply_all_bc_for_nr()
            return
        refresh_vorabrechnung_for_cells(
            self._solver_cells_for_nr(),
            force=force_sources,
            refresh_sources=refresh_sources,
        )
        self._snapshot_bottom_gas_inlet_split_from_vorabrechnung()
        self._apply_all_bc_for_nr()

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
        self._apply_all_bc_for_nr()
        self._seed_initialized_holdup_chain_for_nr()
        self._apply_all_bc_for_nr()

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
        explicit_outer_cap = getattr(cfg, "nr_outer_iter_max", None)
        explicit_inner_cap = getattr(cfg, "nr_inner_iter_max_thesis", None)
        short_probe_mode = int(max_iter) <= 4
        if explicit_outer_cap is not None:
            outer_max = max(1, int(explicit_outer_cap))
        elif bool(cfg.thesis_mode) and bool(cfg.thesis_vorab_sources_single_shot) and not short_probe_mode:
            # Thesis single-shot: drying/DAEM sources fixed after init, but
            # hydrodynamics (K_bd, eps_b, u0, …) must be refreshed after every
            # accepted Newton step because frozen hydro becomes inconsistent with
            # the evolving gas/solid/temperature state.  Hamel's outer Abgleich
            # loop is conceptually "one hydro refresh → one Newton step → repeat",
            # so we allow up to max_iter outer slots (inner_iter_cap per outer will
            # naturally be capped at 1–2 before line-search failure triggers the
            # next outer refresh).
            outer_max = max(10, max_iter)
        else:
            outer_max = max(2, max_iter // 4)
        tol_rms = float(np.clip(0.01 * max(tol, 1.0), 0.005, 0.02))
        # thesis single-shot 模式：每个 outer 允许多次 NR 尝试（默认 5× 内层 budget），
        # 确保 Tikhonov / GD 兜底有足够的迭代机会在外层刷新前收敛。
        if explicit_inner_cap is not None:
            total_inner_budget = max(int(explicit_inner_cap), 1)
        elif bool(cfg.thesis_mode) and bool(cfg.thesis_vorab_sources_single_shot) and not short_probe_mode:
            total_inner_budget = max(int(max_iter) * 5, 20)
        else:
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

        continuation_stages = self._resolve_nr_reaction_continuation_stages()
        continuation_target = self._nr_reaction_continuation_target_cells(solver_cells)
        continuation_active = len(continuation_stages) > 1 and len(continuation_target) > 0
        continuation_warmup_cap = int(max(getattr(cfg, "nr_reaction_continuation_inner_iter_cap_thesis", 0), 0))

        def _set_reaction_multiplier(multiplier: float) -> None:
            for c in solver_cells:
                c.nr_reaction_rate_multiplier = 1.0
            for c in continuation_target:
                c.nr_reaction_rate_multiplier = float(multiplier)

        def _run_outer(stage_idx: int, stage_multiplier: float, inner_budget: int):
            _set_reaction_multiplier(stage_multiplier)
            result = run_outer_abgleich_for_global_nr(
                cells=solver_cells,
                outer_max=outer_max,
                total_inner_budget=int(max(inner_budget, 1)),
                tol=float(tol),
                jacobian_mode=jacobian_mode,
                jacobian_lag=jacobian_lag,
                tol_rms=tol_rms,
                refresh_fn=self._refresh_vorabrechnung_sources_for_nr,
                snapshot_signature_fn=self._vorabrechnung_snapshot_signature,
                solve_inner_fn=_solve_inner,
                strict_check1_before_refresh=(
                    bool(cfg.thesis_mode) and bool(getattr(cfg, "nr_strict_check1_before_refresh_thesis", True))
                ),
            )
            for item in result.outer_history:
                item["reaction_continuation_stage"] = int(stage_idx)
                item["reaction_rate_multiplier"] = float(stage_multiplier)
            return result

        outer_results = []
        try:
            if continuation_active:
                for stage_idx, stage_multiplier in enumerate(continuation_stages[:-1], start=1):
                    stage_budget = continuation_warmup_cap if continuation_warmup_cap > 0 else max(1, total_inner_budget // 2)
                    outer_results.append(_run_outer(stage_idx, stage_multiplier, stage_budget))
            outer_results.append(_run_outer(len(continuation_stages), continuation_stages[-1], total_inner_budget))
        finally:
            for c in solver_cells:
                c.nr_reaction_rate_multiplier = 1.0

        outer_result = outer_results[-1]

        all_nr_history = []
        all_nr_lambdas = []
        all_nr_line_search_trials = []
        all_nr_clip_history = []
        outer_history = []
        agg_nr_timing: dict[str, float] = {}
        agg_nr_counts: dict[str, int] = {}
        outer_refresh_s_total = 0.0
        inner_solve_s_total = 0.0
        used_inner_budget = 0
        outer_iters = 0
        for stage_idx, stage_result in enumerate(outer_results, start=1):
            all_nr_history.extend(stage_result.all_nr_history)
            all_nr_lambdas.extend(stage_result.all_nr_lambdas)
            all_nr_line_search_trials.extend(stage_result.all_nr_line_search_trials)
            for item in stage_result.all_nr_clip_history:
                entry = dict(item)
                entry["reaction_continuation_stage"] = int(stage_idx)
                entry["reaction_rate_multiplier"] = float(continuation_stages[min(stage_idx - 1, len(continuation_stages) - 1)])
                all_nr_clip_history.append(entry)
            outer_history.extend(stage_result.outer_history)
            for key, val in stage_result.agg_nr_timing.items():
                agg_nr_timing[key] = float(agg_nr_timing.get(key, 0.0)) + float(val)
            for key, val in stage_result.agg_nr_counts.items():
                if key == "jacobian_nnz_last":
                    agg_nr_counts[key] = int(val)
                else:
                    agg_nr_counts[key] = int(agg_nr_counts.get(key, 0)) + int(val)
            outer_refresh_s_total += float(stage_result.outer_refresh_s_total)
            inner_solve_s_total += float(stage_result.inner_solve_s_total)
            used_inner_budget += int(stage_result.used_inner_budget)
            outer_iters += int(stage_result.outer_iters)

        last_nr_result = outer_result.last_nr_result
        outer_converged = outer_result.converged_outer
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
                    str(last_nr_result.get("linear_solver_backend"))
                    if last_nr_result.get("linear_solver_backend") is not None
                    else (
                    str(linear_solver_backend)
                    if linear_solver_backend is not None
                    else (
                        "structured_direct"
                        if str(jacobian_mode) in {"block_tridiag_structured", "band_plus_side_elements_structured"}
                        else "sparse_direct_fallback"
                    )
                    )
                ),
                nr_jacobian_lag_steps=jacobian_lag,
                nr_reaction_continuation_active=bool(continuation_active),
                nr_reaction_continuation_stages=list(continuation_stages),
                nr_reaction_continuation_target_cells=len(continuation_target),
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
