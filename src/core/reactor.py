"""Reactor orchestration for thesis-aligned NR solver."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import logging
from time import perf_counter
from typing import List, TYPE_CHECKING

if TYPE_CHECKING:
    from src.core.connectivity_graph import ConnectivityGraph

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
from src.core.hydrodynamics_defaults import (
    HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL,
    HAMEL_HYDRO_BUBBLE_ODE_STRATEGY,
    HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY,
    HAMEL_HYDRO_LAMBDA_STRATEGY,
    HAMEL_HYDRO_PSI_B_STRATEGY,
    HAMEL_HYDRO_U_D_CLOSURE,
    HAMEL_HYDRO_V_B,
    HAMEL_HYDRO_XI_STRATEGY,
    apply_hamel_hydrodynamics_to_cell,
)
from src.core.hamel_parity_defaults import (
    HAMEL_BED_SOLID_STATE_MODEL,
    HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE,
    HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL,
    HAMEL_SIZE_MIGRATION_INVENTORY_MODE,
    PREPROJ_GAS_REACTION_RATE_MULTIPLIER,
    PREPROJ_PHASE_SPLIT_REACTION_RATE_MULTIPLIER,
    PREPROJ_SOLID_HOLDUP_REACTION_RATE_MULTIPLIER,
    resolve_bed_solid_state_model,
    resolve_feed_solid_mass_fractions,
    resolve_size_migration_inventory_mode,
)
from src.core.energy_balance_defaults import (
    DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    DEFAULT_SOLID_ENTHALPY_MODE,
    resolve_reactor_Q_wall_W_total,
)
from src.core.connectivity import (
    apply_all_nr_boundary_data as _apply_all_nr_boundary_data_exec,
    apply_bottom_recycle as _apply_bottom_recycle_exec,
    apply_local_nr_boundary_data as _apply_local_nr_boundary_data_exec,
    propagate_upstream as _propagate_upstream_exec,
    propagated_solid_stream as _propagated_solid_stream,
    recycled_solid_stream as _recycled_solid_stream,
    seed_bed_holdup_chain_from_transport as _seed_bed_holdup_chain_from_transport_exec,
    reconcile_upper_bed_solid_holdup_for_init as _reconcile_upper_bed_solid_holdup_for_init_exec,
    set_bottom_cell_feeds as _set_bottom_cell_feeds_exec,
    snapshot_bottom_gas_inlet_split_from_vorabrechnung as _snapshot_bottom_gas_split_exec,
)
from src.core.connectivity_graph import affected_nr_residual_cells as _affected_nr_cells_graph
from src.core.connectivity_graph import build_solver_cell_registry as _build_solver_cell_registry
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


def _resolve_axial_explicit_heat_loss_watts(cfg: "ReactorConfig") -> tuple[list[float], list[float], float, float]:
    """按床高/自由板高度比例分配 reactor 级 Q_W / Q_WÜ [W]。

    Ref: Hamel (1999) Eq. (2.7) 中 Q_W、Q_WÜ 为显式边界热流项。
    """
    from src.core.energy_balance_defaults import distribute_axial_explicit_heat_loss_watts

    return distribute_axial_explicit_heat_loss_watts(
        H_bed=float(cfg.H_bed),
        H_freeboard=float(cfg.H_freeboard),
        n_bed_cells=int(cfg.n_cells),
        Q_wall_W_total=resolve_reactor_Q_wall_W_total(cfg),
        Q_heat_exchanger_W_total=float(getattr(cfg, "Q_heat_exchanger_W_total", 0.0)),
    )


def resolve_freeboard_heat_loss_bundle(cfg: "ReactorConfig") -> "FreeboardHeatLossBundle":
    """解析自由板轨迹热损 bundle（轴向比例 + 显式 Q）。"""
    from src.core.energy_balance_defaults import FreeboardHeatLossBundle

    _, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
    _, _, fb_q_wall, fb_q_hex = _resolve_axial_explicit_heat_loss_watts(cfg)
    return FreeboardHeatLossBundle(
        heat_loss_mode=str(getattr(cfg, "energy_balance_heat_loss_mode", DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE)),
        heat_loss_frac=float(freeboard_loss),
        Q_wall_W=float(fb_q_wall),
        Q_heat_exchanger_W=float(fb_q_hex),
    )


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
    hydrodynamics_u_d_closure: str = HAMEL_HYDRO_U_D_CLOSURE
    hydrodynamics_bubble_diameter_model: str = HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL
    hydrodynamics_psi_b_strategy: str = HAMEL_HYDRO_PSI_B_STRATEGY
    hydrodynamics_lambda_strategy: str = HAMEL_HYDRO_LAMBDA_STRATEGY
    hydrodynamics_xi_strategy: str = HAMEL_HYDRO_XI_STRATEGY
    hydrodynamics_bubble_velocity_strategy: str = HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY
    hydrodynamics_bubble_ode_strategy: str = HAMEL_HYDRO_BUBBLE_ODE_STRATEGY
    hydrodynamics_v_b: float = HAMEL_HYDRO_V_B  # [-] Hamel Eq. (3.43) 单泡速度权重 v_b
    heat_loss_frac: float = 0.1
    freeboard_heat_loss_frac: float = 0.0  # Optional extra freeboard-only increment
    energy_balance_heat_loss_mode: str = DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE
    solid_enthalpy_mode: str = DEFAULT_SOLID_ENTHALPY_MODE
    Q_wall_W_total: float = 0.0  # [W] Eq. (2.7) Q_W；显式模式时按轴向高度分配至各 cell
    Q_heat_exchanger_W_total: float = 0.0  # [W] Eq. (2.7) Q_WÜ
    freeboard_enabled_reactions: tuple[str, ...] = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R10", "R11", "R12")
    freeboard_closure_enable_char_hetero: bool = False
    freeboard_explicit_enable_char_hetero: bool = True
    freeboard_use_shrunk_bed_top_d_p: bool = True
    freeboard_trajectory_model: str = "analytical_wirsum"
    bed_solid_state_model: str = HAMEL_BED_SOLID_STATE_MODEL
    feed_solid_size_distribution_mode: str = HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE
    nr_init_bottom_zone_solid_energy_preprojection_thesis: bool = False
    nr_init_bottom_primary_exchange_preprojection_thesis: bool = False
    nr_init_bottom_major_local_balance_preprojection_thesis: bool = False
    nr_init_bottom_major_phase_split_preprojection_thesis: bool = False
    nr_init_bed1_phase_split_preprojection_thesis: bool = False
    nr_init_bottom_joint_preprojection_thesis: bool = False
    preproj_solid_holdup_reaction_rate_multiplier_thesis: float = PREPROJ_SOLID_HOLDUP_REACTION_RATE_MULTIPLIER
    preproj_gas_reaction_rate_multiplier_thesis: float = PREPROJ_GAS_REACTION_RATE_MULTIPLIER
    preproj_phase_split_reaction_rate_multiplier_thesis: float = PREPROJ_PHASE_SPLIT_REACTION_RATE_MULTIPLIER
    size_migration_inventory_mode: str = HAMEL_SIZE_MIGRATION_INVENTORY_MODE
    freeboard_trajectory_coeff_model: str = HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL
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
    # Default on (2026-07-07): hybrid_* structured Jacobian; set False to force pure FD overlay off.
    nr_hybrid_jacobian_thesis: bool = True
    nr_damping_halvings_thesis: int = 14
    nr_lambda_init_thesis: float = 1.0  # [-] Hamel Eq. 2.9 first damped trial when prefer_full_step off
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
    # Check1 continuation: if energy RMS rises by this factor, re-allow Vorab refresh next outer.
    # 0 disables (default). Enabling (e.g. 1.25) cuts mid-trajectory maxE but raises final split.
    nr_check1_energy_regress_refresh_ratio_thesis: float = 0.0
    # While Check1 continuation skips outer Vorab refresh, optional tighter inner |ΔT| cap.
    # 0 disables (default). Non-zero values (100–150 K) raise final split in LU seed850 probes.
    nr_inner_abs_dT_budget_check1_continuation_K_thesis: float = 0.0
    # 可选：覆盖 ``clip(0.01*tol_global)`` 的内层 Check1 RMS 阈值（如 VTT 侧块回路用 0.02）。
    nr_inner_tol_rms_thesis: float | None = None
    nr_reaction_continuation_stages_thesis: tuple[float, ...] = (1.0,)
    nr_reaction_continuation_height_m_thesis: float = 0.0
    nr_reaction_continuation_inner_iter_cap_thesis: int = 0
    nr_temperature_fence_enabled_thesis: bool = True
    nr_temperature_fence_lower_margin_K_thesis: float = 150.0
    nr_temperature_fence_upper_margin_K_thesis: float = 350.0
    # 仅底格覆盖相对下沿 [K]；None = 用全局 lower_margin。handoff §5.58。
    nr_temperature_fence_bed0_lower_margin_K_thesis: float | None = None
    # 底格以上床格覆盖相对下沿 [K]；None = 用全局 lower_margin。勿与 upper_margin（上沿）混淆。handoff §5.64。
    nr_temperature_fence_upper_bed_lower_margin_K_thesis: float | None = None
    # 底格 NR 见到的 Startwert T [K]；None = 保持 Vorab 剖面。围栏仍以 T_est 为中心。handoff §5.59。
    vorab_bed0_nr_startwert_T_K_thesis: float | None = None
    # NR 温度下界 [K]：relative 模式下与 ``t_ref - lower_margin`` 取较大者；absolute 模式忽略相对 margin
    nr_temperature_fence_absolute_min_K_thesis: float = 300.0
    nr_temperature_fence_lower_mode_thesis: str = "relative"
    nr_enforce_gas_phase_split_convergence_thesis: bool = True
    # Wirsum (1997/98) Eq. 2.23：线搜索仅当残差欧氏/RMS 范数下降才接受。
    # split/energy 的 max-merit 为工程桥（默认关）；需要时 opt-in。
    # Ref: docs/wirsum_1997_nr_solver_notes.md；handoff §6.6
    nr_line_search_gas_phase_split_merit_thesis: bool = False
    nr_line_search_energy_merit_thesis: bool = False
    # Hydro 语义（P1）：Hamel 外—内结构要求 Inner 在固定映射 F(x; h) 上求解。
    # False（默认 / outer_fixed）：仅 outer refresh 与 Inner 入口 prepare 更新 h；
    #   接受步不刷 hydro，避免「trial 用冻结 h、接受后刷 h 却混用旧 F」的混合语义。
    # True（accept_refresh，工程 opt-in）：每接受步刷 h，且必须重算 F 作下一步 RHS（P0）。
    # Ref: Hamel (1999) Bild 2.2 / §2.1 Vorabrechnung vs Zellenmodell
    nr_refresh_hydrodynamics_on_accepted_step_thesis: bool = False
    # 线搜索试探快氧化投影（工程桥，非 Wirsum；默认关。Vorab x₀ 快氧化仍开）
    nr_line_search_fast_oxidation_projection_thesis: bool = False
    # 线搜索按相快氧化（O₂/H₂ 只在同一相内计量消除）。不是合计 fastox，默认关。
    nr_line_search_per_phase_fastox_thesis: bool = False
    # 按相 fastox 作用的 bed 格数（0=全部 bed）。仅底格会留下上段 R12 悬崖。
    nr_line_search_per_phase_fastox_bed_limit_thesis: int = 0
    # 禁止 Newton 步在同一相内同时增加 O₂ 与合成气（分相切向 clip）。工程桥，默认关。
    nr_clip_same_phase_oxidizer_fuel_step_thesis: bool = False
    # 允许 dN_O2>0 进入已有合成气的相（按相 fastox 烧掉 overlap，剩余 O₂ 点亮 R1）。
    # 仍禁止合成气进入已有 O₂ 的相。工程桥，默认关；需同时开按相 fastox。
    nr_clip_allow_oxidizer_into_fuel_thesis: bool = False
    # 拒绝合成气伪下降（工程桥，非 Wirsum；默认关）
    nr_line_search_syngas_collapse_guard_thesis: bool = False
    nr_line_search_syngas_collapse_ratio_thesis: float = 0.05
    # 固相 FTB 用相对地板，避免 ~1e-6 kg 库存把全局 α 压到 1e-3 以下（§5.6）。
    nr_clip_solid_ftb_relative_floor_thesis: bool = False
    # LS 试探合成气下限（0=关闭）。>0 时在 fastox 前后把 CO/H₂/CH₄ 托到该比例×当前点。
    nr_line_search_syngas_floor_ratio_thesis: float = 0.0
    # Refresh-aware 联合接受（工程桥，默认关）：energy 主导且 energy↓、但 pre-refresh
    # max-merit 未降时，仅当「试探 + 水力刷新」后 merit/split 均不劣于当前才接受。
    # 裸 energy-descent 已证伪（phase2_extent_raw_energy_desc_accept_ab_o8）。
    nr_line_search_refresh_aware_accept_thesis: bool = False
    # Energy 主导时冻结气相 holdup Newton 分量（只动 T/固相），抑制 split 爬升。
    nr_line_search_split_constrained_step_thesis: bool = False
    nr_bottom_joint_preprojection_gas_relax_thesis: float = 0.35
    nr_bottom_joint_bed1_energy_guard_ratio_thesis: float = 1.25
    nr_bottom_joint_neighbor_energy_weight_thesis: float = 2.0
    nr_bottom_joint_energy_aware_gate_thesis: bool = True
    nr_bottom_joint_global_rms_gate_thesis: bool = True
    nr_bottom_joint_split_aware_gate_thesis: bool = True
    nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis: bool = True
    nr_bottom_major_phase_split_preprojection_thesis: bool = True
    nr_bottom_zone_solid_energy_preprojection_thesis: bool = True
    nr_bottom_zone_preprojection_max_bed_index_thesis: int = 2
    nr_bottom_zone_align_downstream_temperature_thesis: bool = False
    nr_bottom_zone_bed1_energy_temperature_preprojection_thesis: bool = True
    nr_bottom_zone_downstream_temperature_blend_thesis: float = 0.75
    nr_outer_preinner_phase_split_preprojection_thesis: bool = True
    nr_outer_preinner_phase_split_max_bed_index_thesis: int = 0
    nr_outer_preinner_phase_split_upper_bed_min_rms_thesis: float = 0.06
    # bed0 split LSQ 可降局部 split 却引爆 NR max_abs（bubble/H2O）；默认守卫回滚
    nr_outer_preinner_split_global_residual_guard_thesis: bool = True
    nr_outer_preinner_split_global_rms_factor_thesis: float = 1.25
    nr_outer_preinner_split_global_max_abs_factor_thesis: float = 1.25
    nr_outer_preinner_split_global_max_abs_delta_thesis: float = 0.5
    # True：bed0 也只做 major 相份额（冻各物种总量），禁止产物 Nd/Nb 独立到 80 凭空造 CO
    nr_outer_preinner_bed0_phase_only_majors_thesis: bool = True
    nr_outer_preinner_bed1_phase_split_thesis: bool = False
    # After bed0 T2, bed1 may become the sole hotspot (H₂); follow-up T2 only then.
    # Always-on bed1 T2 raises energy (falsified). Mid-solve follow-up also hurts E
    # unless gated to late-floor (bed0 cold + global S already low).
    nr_outer_preinner_bed1_hotspot_followup_thesis: bool = True
    nr_outer_preinner_bed1_hotspot_min_rms_thesis: float = 0.05
    nr_outer_preinner_bed1_hotspot_max_bed0_rms_thesis: float = 0.01
    nr_outer_preinner_bed1_hotspot_max_global_rms_thesis: float = 0.04
    # 上段密相 CO/H2：把正的 total residual 吸入 N_d（关 res_sum）。默认关；非 bed0
    # phase-split LSQ 冻物种总量，无法吸收库存（handoff §5.4）。
    nr_outer_preinner_upper_dense_syngas_absorb_thesis: bool = False
    nr_outer_preinner_upper_dense_absorb_min_bed_index_thesis: int = 6
    nr_outer_preinner_upper_dense_absorb_species_thesis: tuple[str, ...] = ("CO", "H2")
    # 全量吸收会抬 split（平台 oneshot）；软吸收 α∈(0,1] 可选。
    nr_outer_preinner_upper_dense_absorb_fraction_thesis: float = 1.0
    # "dense" | "phase_share"；phase_share 按当前相份额分给两相（handoff §5.11）
    nr_outer_preinner_upper_dense_absorb_mode_thesis: str = "dense"
    # True：吸收在全局门禁/local2 之后跑（phase_share 推荐）
    nr_outer_preinner_upper_dense_absorb_post_gate_thesis: bool = False
    # True：跳过吸收后的局部 split/merit 回滚，仅靠全局 preinner merit 门禁。
    nr_outer_preinner_upper_dense_absorb_skip_local_gate_thesis: bool = False
    # 上段 CO 局部 2×2 Newton（密相+气泡 holdup，FD J）；联合库存与 N_ex（§5.7）。
    nr_outer_preinner_upper_co_local2x2_thesis: bool = False
    nr_outer_preinner_upper_co_local2x2_min_bed_index_thesis: int = 6
    nr_outer_preinner_upper_co_local2x2_fraction_thesis: float = 1.0
    # True：local2 在全局门禁/fallback 之后单独跑，避免与 split LSQ 同捆被回滚（handoff §5.8）
    nr_outer_preinner_upper_co_local2x2_post_gate_thesis: bool = True
    # True：床顶 y_CO 达阈值后冻结 local2/phase_share（同一次 solve 内；handoff §5.13）
    nr_outer_freeze_upper_co_bridge_after_y_co_thesis: bool = False
    nr_outer_freeze_upper_co_bridge_y_co_threshold_thesis: float = 0.14
    # >0：累计接受 local2 达 N 次后冻结 upper 桥（比 y_CO 阈值更早；0=关）
    nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis: int = 0
    # True：首轮 outer 结束后若床顶 y_CO≥阈值，嵌套一次
    # ``_solve_global_nr(skip_precalc=True)``（不重 Vorab/init；§5.14；默认关）
    nr_outer_hot_restart_after_y_co_thesis: bool = False
    nr_outer_hot_restart_y_co_threshold_thesis: float = 0.14
    # 续算 max_iter 提示：>0 时 cont_max_iter≈bonus/5；0=复用本轮 max_iter
    nr_outer_hot_restart_bonus_inner_budget_thesis: int = 0
    # 保留：outer_loop 中途门禁复位用（生产路径不读；该路径已证伪）
    nr_outer_hot_restart_bonus_outer_slots_thesis: int = 0
    # True：无入口支撑时将微量 O₂ holdup 粘滞为 0，避免裸动力学速率面跳跃（§5.13；默认关）
    nr_snap_trace_o2_holdup_thesis: bool = False
    nr_snap_trace_o2_atol_mol_s_thesis: float = 1.0e-12  # [mol/s]
    # 仅 bed 序号 ≥ 该值启用（2=专治 bed2+；0=所有 bed）
    nr_snap_trace_o2_min_bed_index_thesis: int = 0
    # True：速率评价时把微量 O₂ 当 0（不动 holdup；§5.15；默认关）
    nr_trace_o2_rate_gate_thesis: bool = False
    nr_trace_o2_rate_gate_atol_mol_s_thesis: float = 1.0e-12  # [mol/s]
    nr_trace_o2_rate_gate_min_bed_index_thesis: int = 2
    # True：仅当床顶 y_CO≥阈值后才武装 rate-gate（先抬 CO，再压 bed2 速率面）
    nr_trace_o2_rate_gate_after_y_co_thesis: bool = False
    nr_trace_o2_rate_gate_y_co_threshold_thesis: float = 0.14
    nr_outer_preinner_solid_energy_preprojection_thesis: bool = True
    nr_outer_preinner_solid_energy_max_bed_index_thesis: int = 2
    nr_outer_preinner_on_check1_continuation_thesis: bool = True
    # Bild 2.2: Check1 continuation stays in NR; holdup remap is Node-A style (default off).
    nr_outer_preinner_holdup_on_check1_continuation_thesis: bool = False
    nr_outer_preinner_global_merit_gate_thesis: bool = True
    nr_outer_preinner_merit_driven_thesis: bool = True
    nr_outer_preinner_split_tier_min_rms_thesis: float = 0.03
    # Also run T2 when a single species/cell dominates (global RMS can sit just
    # below tau_split while max|split| remains large — bed0 CO₂ N_ex gap).
    nr_outer_preinner_split_tier_min_max_abs_thesis: float = 0.05
    nr_outer_preinner_energy_tier_min_rms_thesis: float = 0.05
    nr_outer_preinner_solid_tier_min_rms_thesis: float = 0.01
    nr_outer_preinner_reconcile_holdup_before_zone_thesis: bool = True
    nr_outer_preinner_split_shock_rms_thesis: float = 0.15
    nr_outer_preinner_max_bed_energy_scaled_thesis: float = 0.08
    nr_outer_preinner_holdup_flux_mismatch_min_thesis: float = 0.15
    nr_outer_preinner_holdup_stream_imbalance_min_thesis: float = 0.15
    nr_outer_preinner_split_shock_bed0_improvement_thesis: float = 0.65
    nr_outer_preinner_split_plateau_bed0_improvement_thesis: float = 0.85
    nr_outer_postinner_holdup_reconcile_thesis: bool = False
    # Post-inner 自由板焓投影会经返料耦合拉偏床层（出口 CO 回归）；默认关。
    # 求解结束后的一次性投影见 ``nr_finalize_freeboard_energy_T_projection_thesis``。
    nr_outer_postinner_freeboard_energy_T_projection_thesis: bool = False
    # Finalize：单次轴向 T→Eq.2.7，清除 fb 能量鬼峰而不进入外层 Newton 反馈。
    nr_finalize_freeboard_energy_T_projection_thesis: bool = True
    # True：床顶 bed 的 Eq.2.7 固体焓出口改用入口透传（类 freeboard_closure；§5.19；默认关）
    nr_bed_top_solid_energy_passthrough_thesis: bool = False
    # True：仅当床顶 y_CO≥阈值后才武装床顶焓透传（§5.20；须同时开 passthrough）
    nr_bed_top_solid_energy_passthrough_after_y_co_thesis: bool = False
    nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis: float = 0.14
    # True：preinner post-gate 对床顶 bed 跑 solid holdup LSQ（§5.21；默认关）
    nr_outer_preinner_bed_top_solid_holdup_thesis: bool = False
    nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis: bool = False
    nr_outer_preinner_bed_top_solid_holdup_y_co_threshold_thesis: float = 0.14
    # 自床顶向下计入的 bed 数（1=仅 bed9；2=bed8+bed9；§5.23）
    nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis: int = 1
    # post-gate 强制对指定 bed 跑 phase-split LSQ（§5.24；默认空；勿默认开）
    nr_outer_preinner_targeted_phase_split_beds_thesis: tuple[int, ...] = ()
    nr_outer_preinner_targeted_phase_split_after_y_co_thesis: bool = False
    nr_outer_preinner_targeted_phase_split_y_co_threshold_thesis: float = 0.14
    # post-gate 指定 bed 物种总量闭合（含负 res_sum 放出；§5.25；默认空）
    nr_outer_preinner_targeted_species_sum_beds_thesis: tuple[int, ...] = ()
    nr_outer_preinner_targeted_species_sum_species_thesis: tuple[str, ...] = ("N2",)
    nr_outer_preinner_targeted_species_sum_fraction_thesis: float = 1.0
    nr_outer_preinner_targeted_species_sum_mode_thesis: str = "phase_share"
    nr_outer_preinner_targeted_species_sum_bidirectional_thesis: bool = True
    nr_outer_preinner_targeted_species_sum_after_y_co_thesis: bool = False
    nr_outer_preinner_targeted_species_sum_y_co_threshold_thesis: float = 0.14
    nr_inner_record_group_history_thesis: bool = False
    # [K] 每步 |ΔT| 帽。600 允许一刀砸进围栏；CORE 覆写 40（handoff §5.53）。
    nr_inner_t_step_default_cap_K_thesis: float = 600.0
    nr_inner_split_dominant_t_step_cap_K_thesis: float = 150.0
    nr_inner_split_dominant_merit_ratio_thesis: float = 0.85
    # Absolute |ΔT| budget per inner solve (0 = disabled). Caps multi-iter accumulation.
    nr_inner_abs_dT_budget_K_thesis: float = 0.0
    # F3: keep 150 K T-cap across outers after split-dominant engage until release.
    nr_inner_t_cap_split_hysteresis_thesis: bool = True
    nr_inner_split_dominant_t_cap_release_ratio_thesis: float = 0.70
    nr_outer_merit_early_stop_thesis: bool = True
    nr_outer_merit_early_stop_min_iters_thesis: int = 8
    nr_outer_merit_early_stop_regression_rel_thesis: float = 0.25
    nr_outer_merit_early_stop_regression_patience_thesis: int = 1
    nr_outer_merit_early_stop_plateau_rel_thesis: float = 0.002
    # Solid-holdup LSQ in bottom-zone / pre-inner bridges (Phase2 perf defaults).
    nr_bottom_zone_solid_holdup_lsq_max_nfev_thesis: int = 12
    nr_bottom_zone_solid_holdup_skip_rms_scaled_thesis: float = 0.0
    nr_bottom_zone_solid_repass_thesis: bool = True
    nr_outer_preinner_solid_holdup_skip_rms_scaled_thesis: float = 0.01
    nr_outer_preinner_solid_repass_thesis: bool = False
    nr_bottom_zone_bed1_energy_lsq_max_nfev_thesis: int = 12
    # C-section: NH3/TAR1/TAR2 enter global NR unknowns (H2S stays Gibbs-only).
    nr_include_pyrolysis_gas_species_in_unknowns_thesis: bool = False
    # Bottom x0: seed NH3/TAR dense holdup from fixed Vorabrechnung pyrolysis sources.
    nr_seed_pyrolysis_gas_dense_holdup_x0_thesis: bool | None = None
    nr_pyrolysis_dense_holdup_seed_scale_thesis: float = 1.0
    vorab_bottom_temperature_cap_K_thesis: float = 1180.0
    vorab_bottom_temperature_floor_K_thesis: float = 1050.0
    # HTW Vorab：沿床焓步进 T 剖面（取代全局绝热 clip 作为 bed 初温；不读 validation JSON）
    vorab_bed_temperature_march_thesis: bool = False
    vorab_bed_temperature_march_heat_loss_frac_thesis: float = 0.05
    vorab_bed_temperature_march_dH_comb_J_per_mol_O2_thesis: float = 393_500.0
    vorab_bed_temperature_march_drying_enthalpy_J_per_kg_thesis: float = 2.26e6
    vorab_bed_temperature_march_pyrolysis_enthalpy_J_per_kg_daf_vm_thesis: float = 1.0e6
    vorab_bed_temperature_march_Cp_gas_J_per_mol_K_thesis: float = 32.0
    # HTW 多段主气化剂 N_zu 分级进料（锥段喷嘴占位；与 holdup 轴向 profile 互补）
    vorab_staged_primary_gas_injection_thesis: bool = False
    vorab_staged_primary_gas_bed_cells_thesis: int = 3
    vorab_staged_primary_gas_weight_decay_thesis: float = 1.5
    # 可选：显式 bed0..bed(n-1) 主气化剂份额 [-]，非 None 时覆盖 decay 指数权重
    vorab_staged_primary_gas_layer_weights_thesis: tuple[float, ...] | None = None
    vorab_staged_primary_gas_dense_frac_thesis: float = 0.70
    # 床底脱挥发分区：DAEM 总 VM 沿前 N 个 bed cell 按高度展开
    # （Hamel 1999 Kap.2.1：Freisetzung … in Abhängigkeit der Reaktorhöhe → 各 cell 份额）。
    # 0 = legacy 瞬时全在 bed0（工程偏离）。
    vorab_vm_devolatilization_zone_cells_thesis: int = 0
    # linear=均分；front_loaded=ζ^0.5 前重（利于 bed0 点火，仍跨 zone 格）
    vorab_vm_devolatilization_zone_profile_thesis: str = "linear"
    # 沿床焓步进时 bed 初温可低于 Gibbs 可靠区间；仅用于 major-Gibbs x0 求解温度下限
    vorab_major_gibbs_min_temperature_K_thesis: float = 750.0
    # 沿床焓步进初温偏低时，NR 温度围栏上界需锚定到可操作床温而非 t_ref+350K
    vorab_nr_temperature_fence_upper_anchor_K_thesis: float = 1400.0
    # NR 初温 seed（与 Vorab 冷预算 t_budget_K 解耦；不读 validation JSON）
    vorab_nr_temperature_seed_bottom_K_thesis: float = 750.0
    vorab_nr_temperature_seed_top_K_thesis: float = 1100.0
    vorab_nr_temperature_seed_progress_power_thesis: float = 0.65
    vorab_nr_temperature_seed_o2_progress_weight_thesis: float = 0.35
    vorab_nr_temperature_seed_bed0_cold_anchor_thesis: bool = True
    # cold-anchor 时 bed0 = T_bottom + soft*(T_top-T_bottom)*O₂/ξ进度，soft∈[0,1]
    vorab_nr_temperature_seed_bed0_soft_progress_thesis: float = 0.0
    vorab_nr_temperature_seed_mode_thesis: str = "adiabatic_anchor"
    vorab_nr_temperature_seed_adiabatic_blend_thesis: float = 0.85
    vorab_nr_temperature_fence_lower_margin_K_thesis: float = 80.0
    # Init 末段：bed1+ char/ash holdup 对齐 ``m_solid_auf_in / K_auf`` 输运支撑（Eq.2.6–2.7）
    nr_init_align_upper_bed_solid_holdup_thesis: bool = True
    # True：holdup support/align 用 in+R（含负反应汇；§5.22；默认关）
    nr_holdup_support_include_negative_R_thesis: bool = False
    nr_init_align_upper_bed_holdup_allow_reduce_thesis: bool = False
    thesis_mode: bool = False
    explicit_side_blocks_enabled: bool = True
    # B-tier Verbindungsmatrix: ``from_config`` builds graph from segment counts;
    # set ``connectivity_graph_path`` to load JSON instead.
    connectivity_preset: str = "from_config"
    connectivity_graph_path: str | None = None
    n_cyclone_stages_thesis: int = 1
    n_return_legs_thesis: int = 1
    return_leg_bed_targets_thesis: tuple[int, ...] | None = None
    # Hamel A1-style major-species Gibbs seed in Vorabrechnung x0 generation (legacy per-cell extension).
    vorab_major_gibbs_x0: bool = False
    # x₀ holdup mode: transport_from_vorab | major_gibbs_per_cell | legacy_heuristic
    vorab_x0_holdup_mode_thesis: str = "transport_from_vorab"
    # transport x₀ 叠加 o2_remaining+热解 H₂ 后的快氧化 Startwert（R12→R5→R6 + 残氧 R1）；同构裸动力学必需
    vorab_transport_x0_fast_oxidation_closure_thesis: bool = True
    # 气相快氧化后，残氧按 R1 接到炭库存（Startwert；关掉则留给底格 solve_cell）
    vorab_transport_x0_char_oxidation_o2_closure_thesis: bool = True
    # 底格 Eq.2.4/2.5 总量闭合 + 两相分相 Startwert（O₂ 气泡 / 合成气悬浮）。不改 R(x)。
    vorab_bed0_eq24_two_phase_startwert_thesis: bool = False
    # 分相后从气泡挪 stoich+seed 的 O₂ 进悬浮相并按相 fastox，留下点亮 R1 的残氧 [mol/s]。0=关。
    vorab_bed0_r1_oxidizer_seed_mol_s_thesis: float = 0.0
    # 底格单格 least_squares（opt-in）；冷根/重叠变差回滚。默认关：裸 R12 会卡住。
    vorab_init_solve_bottom_cell_thesis: bool = False
    vorab_init_solve_bottom_cell_T_floor_K_thesis: float = 800.0
    # Hamel §4.3 Gibbs for incremental pyrolysis gas composition in Vorab budget/x₀
    vorab_pyrolysis_gibbs_composition_thesis: bool = True
    vorab_pyrolysis_gibbs_solver_mode_thesis: str = "hamel_reduced"
    # Major-gibbs init solver mode: augmented | hamel_reduced | shadow_compare
    major_gibbs_solver_mode: str = "augmented"
    # Hamel strict Vorabrechnung: drying/DAEM sources are budgeted once and fixed.
    thesis_vorab_sources_single_shot: bool = True
    # A-tier Startwertwahl: unified Vorabrechnung cell budgets (tau_bed, cumulative O2/VM/dry).
    vorab_a_tier_startwert_budget_thesis: bool = True
    # Init 阶段 3–5：bed 格 ``u0_target`` 取自 Vorab 气量（主气 + 累积热解），再刷新 K/BC。
    vorab_init_bed_vorab_gas_u0_target_thesis: bool = True
    # P0 Vorab→cell mapping (bed): opt-in; default off until full zu/ab stream mapping lands.
    vorab_a_tier_cell_mapping_thesis: bool = False
    # A-tier bed zu/ab stream holdup mapping (bed0 primary anchored to N_zu; default on with budget).
    vorab_a_tier_bed_stream_mapping_thesis: bool = True
    # Upstream bed cells: inflow-chain holdup + local phase-split closure (experimental; default off).
    vorab_a_tier_bed_stream_upstream_closure_thesis: bool = False
    # Axial primary holdup profile: O2 from Vorab budget, H2O/N2 from upstream inflow (opt-in).
    vorab_a_tier_bed_axial_primary_profile_thesis: bool = False
    vorab_a_tier_bed_axial_primary_profile_cells_thesis: int = 0
    vorab_a_tier_bed_axial_primary_blend_thesis: float = 0.35
    # True：无 O2 补气格 H2O 也按 inflow 透传（§5.28）。默认 False——与 N2 透传叠开时
    # 存在 CO≈0.15 / 0.21 双吸引子；opt-in 且须配 LSQ2+延迟焓透传（§5.29）。
    vorab_a_tier_post_staged_h2o_passthrough_thesis: bool = False
    # True：无 O2 补气格 CO/CO2/H2/CH4/TAR 在 fastox 后再贴到 inflow（§5.63）。
    # 默认 False。与 H2O 透传同一时序洞，须单独 A/B，勿叠开验收。
    vorab_a_tier_post_staged_syngas_passthrough_thesis: bool = False
    # bed0–bed1 major phase-split after product/source Vorab closure (default off; enable when solid chain stable).
    vorab_a_tier_bed01_phase_split_after_vorab_closure_thesis: bool = False
    vorab_a_tier_bed01_phase_split_max_bed_index_thesis: int = 0
    # bed1 major phase-split after zu/ab stream mapping (recycle / mid-bed closure; VTT opt-in).
    vorab_a_tier_bed1_phase_split_after_stream_mapping_thesis: bool = False
    # Align CO/CO2/H2/CH4 holdup with fixed Vorab gas sources on bottom bed cells.
    vorab_a_tier_product_holdup_source_closure_thesis: bool = True
    vorab_a_tier_product_holdup_source_closure_bed_cells_thesis: int = 0
    # Init Check2-lite: Vorab budget vs cell x0 alignment gate + corrective bridge.
    vorab_a_tier_init_abgleich_bridge_thesis: bool = True
    vorab_a_tier_init_abgleich_max_passes_thesis: int = 2
    vorab_a_tier_init_abgleich_primary_rel_tol_thesis: float = 0.05
    vorab_a_tier_init_abgleich_temperature_K_thesis: float = 30.0
    vorab_a_tier_init_abgleich_product_rel_tol_thesis: float = 0.40
    vorab_a_tier_init_abgleich_o2_profile_rel_tol_thesis: float = 0.85
    vorab_a_tier_init_abgleich_solid_rel_tol_thesis: float = 0.25
    vorab_a_tier_init_abgleich_temperature_blend_thesis: float = 0.25
    vorab_a_tier_init_abgleich_temperature_blend_march_thesis: float = 0.0
    # Outer Check2 fail: lightweight Vorab ↔ cell bridge before next refresh (Hamel Node-B).
    vorab_a_tier_outer_abgleich_bridge_thesis: bool = True
    vorab_a_tier_outer_abgleich_max_passes_thesis: int = 1
    # Bed reactive solid (VM/moisture) holdup from Vorab dry/VM budget after transport chain seed.
    vorab_a_tier_bed_solid_reactive_profile_thesis: bool = True
    # Outer Check2 (beyond ΔT): bed major phase-split RMS gate (thesis A-tier).
    nr_outer_check2_rel_tol: float = 1e-4
    nr_outer_check2_abs_tol_T_K: float = 1.0  # [K] Check2 ΔT 绝对容差（Vorab↔cell 对齐）
    nr_outer_check2_gas_phase_split_gate_thesis: bool = True
    nr_outer_check2_max_gas_phase_split_rms_thesis: float = 0.05
    nr_outer_check2_bed0_primary_gate_thesis: bool = False
    nr_outer_check2_max_bed0_primary_rel_drift_thesis: float = 0.10
    nr_bottom_joint_preprojection_when_cell_mapped_thesis: bool = False
    nr_bottom_primary_exchange_preprojection_when_cell_mapped_thesis: bool = False
    nr_bottom_phase_split_preprojection_when_cell_mapped_thesis: bool = False

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
    # 若设置，仅 O₂ 主进料走该密相份额；H₂O/N₂ 仍用 ``gas_inlet_dense_frac``。
    # 贫氧床底应将 O₂ 优先送入密相（R1 异相燃烧），减少气泡相 R5 CO 氧化。
    gas_inlet_o2_dense_frac: float | None = None
    gas_inlet_split_strategy: str = "fixed"

    recirculation_frac: float = 0.1
    recycle_gas: bool = True
    enable_r12: bool = True
    solid_lower_inlet_frac: float = 0.0
    top_solid_inlet_frac: float = 0.0
    allow_reactive_solid_propagation: bool = False
    reactive_solid_cutoff_xi: float = 0.35
    r4_scale: float = 1.0
    r2_scale: float = 1.0  # 异相水煤气 (R2) 强度倍率；Char 还原标定用
    r5_scale: float = 1.0
    # Hamel 同构默认关：床层代数源项燃料预封顶 + 比例 O₂/H₂O 限幅为工程路径（opt-in）。
    extent_limiters_enabled_thesis: bool = False
    # 自由板轴向轨迹（显式 dt 步进）的库存界：与床层代数限幅解耦；默认开，防焓步进炸温。
    freeboard_extent_dt_bounds_thesis: bool = True
    # 仅 bed0 气泡相 R5 额外倍率（1.0=不变）；与 ``cell.r5_bubble_scale`` 相乘。
    r5_bubble_scale_bed0: float = 1.0
    r6_scale: float = 1.0
    r7_scale: float = 1.0
    r8_scale: float = 1.0  # WGSR (R8) 强度倍率；CO/CO2 标定用
    # 灵敏度：Gibbs 热解后将 split_frac×CH4 重映射为 2 H2 + CO（不改 §4.3 求解）
    pyrolysis_ch4_to_h2_co_split_frac: float = 0.0
    # 灵敏度：削减 split_frac×CO 并注入 2×trim H2（不增加 CO；C 隐式转 char）
    pyrolysis_co_reduce_h2_inject_frac: float = 0.0

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
            self.nr_include_pyrolysis_gas_species_in_unknowns_thesis = True
            if self.nr_seed_pyrolysis_gas_dense_holdup_x0_thesis is None:
                self.nr_seed_pyrolysis_gas_dense_holdup_x0_thesis = True
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
            self.config.nr_include_pyrolysis_gas_species_in_unknowns_thesis = True
            if self.config.nr_seed_pyrolysis_gas_dense_holdup_x0_thesis is None:
                self.config.nr_seed_pyrolysis_gas_dense_holdup_x0_thesis = True
        self.cells: List[Cell] = []
        self.freeboard_cells: List[Cell] = []
        self._last_explicit_freeboard_closure: dict | None = None
        self.side_cells: List[Cell] = []
        self.cyclone_cells: List[Cell] = []
        self.return_leg_cells: List[Cell] = []
        self.side_block_cells_by_id: dict[str, Cell] = {}
        self.cyclone_cell: Cell | None = None
        self.return_leg_cell: Cell | None = None
        self.connectivity_graph: ConnectivityGraph | None = None
        self._thesis_fixed_vorab_sources_ready: bool = False
        self._build_cells()
        self._apply_heat_loss_distribution()
        from src.core.connectivity_graph import build_connectivity_from_reactor_config

        self.connectivity_graph = build_connectivity_from_reactor_config(self.config, reactor=self)
        self._build_side_block_cells_from_graph()
        configure_tar_components_by_fuel(config.fuel_type)
        from src.solvers.nr_indexing import install_main_nr_gas_indexing_on_cells

        install_main_nr_gas_indexing_on_cells(self._solver_cells_for_nr(), self.config)
        from src.solvers.nr_indexing import install_nr_jacobian_metadata_on_cells

        install_nr_jacobian_metadata_on_cells(self._solver_cells_for_nr(), self.config)

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
        mass_fractions = resolve_feed_solid_mass_fractions(
            n_size_classes=nk,
            mode=str(cfg.feed_solid_size_distribution_mode),
        )

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

        bed_solid_model = resolve_bed_solid_state_model(
            configured=str(cfg.bed_solid_state_model),
            thesis_mode=bool(cfg.thesis_mode),
        )
        size_migration_mode = resolve_size_migration_inventory_mode(
            configured=str(cfg.size_migration_inventory_mode),
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
                # Thesis mode enforces Hamel hydrodynamics even if config was overridden.
                apply_hamel_hydrodynamics_to_cell(cell)
            cell.solid_state_model = bed_solid_model
            cell.size_migration_inventory_mode = size_migration_mode
            cell.hydrodynamics_v_b = float(cfg.hydrodynamics_v_b)
            cell.enable_r12 = bool(cfg.enable_r12)
            cell.extent_limiters_enabled = bool(
                getattr(cfg, "extent_limiters_enabled_thesis", False)
            )
            cell.r4_scale = float(max(cfg.r4_scale, 0.0))
            cell.r2_scale = float(max(cfg.r2_scale, 0.0))
            cell.r5_scale = float(max(cfg.r5_scale, 0.0))
            cell.r5_bubble_scale = float(max(getattr(cfg, "r5_bubble_scale_bed0", 1.0), 0.0)) if i == 0 else 1.0
            cell.pyrolysis_ch4_to_h2_co_split_frac = float(
                np.clip(getattr(cfg, "pyrolysis_ch4_to_h2_co_split_frac", 0.0), 0.0, 1.0)
            )
            cell.pyrolysis_co_reduce_h2_inject_frac = float(
                np.clip(getattr(cfg, "pyrolysis_co_reduce_h2_inject_frac", 0.0), 0.0, 1.0)
            )
            cell.r6_scale = float(max(cfg.r6_scale, 0.0))
            cell.r7_scale = float(max(cfg.r7_scale, 0.0))
            cell.r8_scale = float(max(cfg.r8_scale, 0.0))
            cell._n_vorab_cells = cfg.n_cells  # legacy uniform-mesh tau fallback
            cell._vorab_bed_height = float(cfg.H_bed)
            self.cells.append(cell)
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
                fb_cell.extent_limiters_enabled = bool(
                    getattr(cfg, "extent_limiters_enabled_thesis", False)
                )
                fb_cell.r4_scale = float(max(cfg.r4_scale, 0.0))
                fb_cell.r2_scale = float(max(cfg.r2_scale, 0.0))
                fb_cell.r5_scale = float(max(cfg.r5_scale, 0.0))
                fb_cell.r6_scale = float(max(cfg.r6_scale, 0.0))
                fb_cell.r7_scale = float(max(cfg.r7_scale, 0.0))
                fb_cell.r8_scale = float(max(cfg.r8_scale, 0.0))
                fb_cell.freeboard_u_gb_scale = float(max(cfg.freeboard_u_gb_scale, 1e-6))
                fb_cell.freeboard_explicit_char_hetero_enabled = bool(cfg.freeboard_explicit_enable_char_hetero)
                fb_cell.solid_state_model = "freeboard_closure"
                self.freeboard_cells.append(fb_cell)

    def _build_side_block_cells_from_graph(self) -> None:
        """Instantiate cyclone/return-leg NR cells declared in the frozen connectivity graph."""
        if not self._use_explicit_side_block_cells():
            return
        graph = self.connectivity_graph
        if graph is None:
            return
        from src.core.connectivity_graph import ordered_solver_block_ids

        cfg = self.config
        top_dh = float(self.cells[-1].geo.dh) if self.cells else float(cfg.H_bed) / max(int(cfg.n_cells), 1)
        solid = self.cells[0].solid if self.cells else SolidProps(
            rho_s=cfg.rho_s,
            eps_mf=cfg.eps_mf,
            n_size_classes=max(1, int(cfg.n_age_classes)),
        )

        def _make_side_cell(block_id: str, cell_type: str, stage_index: int) -> Cell:
            h_offset = 0.5 + 0.5 * float(stage_index)
            geo = CellGeometry(
                D_bed=cfg.D_bed,
                dh=top_dh,
                h_center=float(cfg.H_bed) + h_offset * top_dh,
            )
            cell = Cell(geo=geo, solid=solid, fuel_type=cfg.fuel_type)
            cell.P = cfg.P
            cell.cell_type = cell_type
            cell.solid_state_model = "holdup_transport"
            cell.use_gibbs_minor = bool(cfg.use_gibbs_minor)
            cell.enable_r12 = bool(cfg.enable_r12)
            cell.extent_limiters_enabled = bool(
                getattr(cfg, "extent_limiters_enabled_thesis", False)
            )
            cell.r4_scale = float(max(cfg.r4_scale, 0.0))
            cell.r2_scale = float(max(cfg.r2_scale, 0.0))
            cell.r5_scale = float(max(cfg.r5_scale, 0.0))
            cell.r6_scale = float(max(cfg.r6_scale, 0.0))
            cell.r7_scale = float(max(cfg.r7_scale, 0.0))
            cell.r8_scale = float(max(cfg.r8_scale, 0.0))
            cell._side_block_graph_id = block_id
            return cell

        cyclone_ids = ordered_solver_block_ids(graph, "cyclone")
        return_ids = ordered_solver_block_ids(graph, "return_leg")
        self.cyclone_cells = [
            _make_side_cell(block_id, "cyclone", i) for i, block_id in enumerate(cyclone_ids)
        ]
        self.return_leg_cells = [
            _make_side_cell(block_id, "return_leg", i) for i, block_id in enumerate(return_ids)
        ]
        self.side_block_cells_by_id = {
            block_id: cell
            for block_id, cell in zip(cyclone_ids, self.cyclone_cells)
        }
        self.side_block_cells_by_id.update(
            {block_id: cell for block_id, cell in zip(return_ids, self.return_leg_cells)}
        )
        side_blocks = [
            block
            for block in graph.blocks
            if block.kind in {"cyclone", "return_leg"} and block.solver_cell_index is not None
        ]
        side_blocks.sort(key=lambda block: int(block.solver_cell_index))
        self.side_cells = [self.side_block_cells_by_id[block.block_id] for block in side_blocks]
        self.cyclone_cell = self.cyclone_cells[0] if self.cyclone_cells else None
        self.return_leg_cell = self.return_leg_cells[0] if self.return_leg_cells else None

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
            base = "band_plus_side_elements_structured"
        else:
            base = "block_tridiag_structured"
        if bool(getattr(self.config, "nr_hybrid_jacobian_thesis", False)):
            from src.solvers.global_nr_solver import resolve_hybrid_jacobian_strategy

            return resolve_hybrid_jacobian_strategy(base)
        return base

    def _nr_jacobian_uses_structured_direct_linear_solver(self, jacobian_mode: str) -> bool:
        from src.solvers.global_nr_solver import structured_jacobian_uses_direct_linear_solver

        return structured_jacobian_uses_direct_linear_solver(jacobian_mode)

    def _initialize_explicit_side_block_states(self) -> None:
        """委托 ``side_block_bridge``：主链顶端状态 -> side-block 初始化。"""
        _init_side_blocks_bridge(self)

    def _solver_cell_registry(self) -> dict:
        return _build_solver_cell_registry(self)

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
            connectivity_graph=self.connectivity_graph,
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

        Upper-bed ``K_auf`` reconciliation runs separately via
        ``_reconcile_upper_bed_solid_holdup_for_nr`` after abgleich / preprojection
        so mid-bridge pyrolysis Gibbs seeding is not disturbed.
        """
        top_above_cell = self.freeboard_cells[0] if self.freeboard_cells else None
        _seed_bed_holdup_chain_from_transport_exec(self.cells, top_above_cell=top_above_cell)

    def _reconcile_upper_bed_solid_holdup_for_nr(
        self,
        *,
        allow_reduce: bool | None = None,
    ) -> bool:
        """Raise bed1+ char/ash holdup toward current transport inflow support."""
        if not bool(getattr(self.config, "nr_init_align_upper_bed_solid_holdup_thesis", True)):
            return False
        top_above_cell = self.freeboard_cells[0] if self.freeboard_cells else None
        return _reconcile_upper_bed_solid_holdup_for_init_exec(
            self.cells,
            self.config,
            top_above_cell=top_above_cell,
            allow_reduce=allow_reduce,
        )

    def _apply_heat_loss_distribution(self) -> None:
        """Refresh per-cell energy-balance heat-loss settings from reactor config."""
        cfg = self.config
        bed_losses, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
        bed_q_wall, bed_q_hex, fb_q_wall, fb_q_hex = _resolve_axial_explicit_heat_loss_watts(cfg)
        heat_loss_mode = str(getattr(cfg, "energy_balance_heat_loss_mode", DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE))
        solid_enthalpy_mode = str(getattr(cfg, "solid_enthalpy_mode", DEFAULT_SOLID_ENTHALPY_MODE))
        for cell, loss_frac, q_wall, q_hex in zip(self.cells, bed_losses, bed_q_wall, bed_q_hex):
            cell.heat_loss_frac = float(loss_frac)
            cell.heat_loss_mode = heat_loss_mode
            cell.solid_enthalpy_mode = solid_enthalpy_mode
            cell.Q_wall_W = float(q_wall)
            cell.Q_heat_exchanger_W = float(q_hex)
        n_fb = len(self.freeboard_cells)
        # Match bed compounding: split the freeboard height share across NR cells so
        # ``inlet_fraction`` mode matches ``simulate_freeboard`` heat_loss_frac_total.
        total_height = float(max(cfg.H_bed, 0.0) + max(cfg.H_freeboard, 0.0))
        fb_weight_each = (
            (float(max(cfg.H_freeboard, 0.0)) / total_height) / max(n_fb, 1) if total_height > 0.0 else 0.0
        )
        fb_loss_each = (
            _compound_heat_loss_fraction(cfg.heat_loss_frac, fb_weight_each) if n_fb > 0 else 0.0
        )
        if float(getattr(cfg, "freeboard_heat_loss_frac", 0.0) or 0.0) > 0.0 and n_fb > 0:
            # Keep total freeboard_loss (incl. extra) while preferring per-cell compound weights.
            fb_loss_each = float(freeboard_loss) / float(n_fb) if n_fb > 1 else float(freeboard_loss)
        fb_q_wall_each = float(fb_q_wall) / max(n_fb, 1)
        fb_q_hex_each = float(fb_q_hex) / max(n_fb, 1)
        for cell in self.freeboard_cells:
            cell.heat_loss_frac = float(fb_loss_each)
            cell.heat_loss_mode = heat_loss_mode
            cell.solid_enthalpy_mode = solid_enthalpy_mode
            cell.Q_wall_W = fb_q_wall_each
            cell.Q_heat_exchanger_W = fb_q_hex_each
        for cell in self.side_cells:
            cell.heat_loss_frac = 0.0
            cell.heat_loss_mode = heat_loss_mode
            cell.solid_enthalpy_mode = solid_enthalpy_mode
            cell.Q_wall_W = 0.0
            cell.Q_heat_exchanger_W = 0.0

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
        """Freeze bottom inlet dense fraction after Vorab hydro refresh.

        Hamel Bild 2.2: outer Vorab refresh updates hydro/zu-ab boundaries only;
        NR inventories (N_d/N_b) stay as accepted unknowns until the next inner solve.
        ``align_bottom_primary_gas_state_to_inlet_split`` is init/mapping-only
        (see ``apply_vorabrechnung_bed_cell_mapping`` / init_precalc_step).
        """
        _snapshot_bottom_gas_split_exec(self.cells, self.config)

    def _apply_bottom_recycle(self, relax: float | None = None) -> None:
        if self._use_side_blocks_in_nr_boundary_path() and self.return_leg_cell is not None:
            _apply_all_nr_boundary_data_exec(
                self.cells,
                self.config,
                cyclone_cell=self.cyclone_cell,
                return_leg_cell=self.return_leg_cell,
                connectivity_graph=self.connectivity_graph,
                cell_registry=self._solver_cell_registry(),
            )
            return
        _apply_bottom_recycle_exec(self.cells, self.config, relax)

    def _propagate_upstream(self, i: int) -> None:
        _propagate_upstream_exec(
            self.cells,
            self.config,
            i,
            connectivity_graph=self.connectivity_graph,
            cell_registry=self._solver_cell_registry(),
        )

    def _apply_all_bc_for_nr(self) -> None:
        registry = self._solver_cell_registry()
        _apply_all_nr_boundary_data_exec(
            self.cells,
            self.config,
            freeboard_cells=self.freeboard_cells if self._use_explicit_freeboard_solver_graph() else None,
            cyclone_cell=self.cyclone_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            return_leg_cell=self.return_leg_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            connectivity_graph=self.connectivity_graph,
            cell_registry=registry,
        )

    def _apply_local_bc_for_nr(self, changed_cell_idx: int) -> None:
        registry = self._solver_cell_registry()
        _apply_local_nr_boundary_data_exec(
            self.cells,
            self.config,
            changed_cell_idx,
            freeboard_cells=self.freeboard_cells if self._use_explicit_freeboard_solver_graph() else None,
            cyclone_cell=self.cyclone_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            return_leg_cell=self.return_leg_cell if self._use_side_blocks_in_nr_boundary_path() else None,
            connectivity_graph=self.connectivity_graph,
            cell_registry=registry,
        )

    def _refresh_vorabrechnung_sources_for_nr(self, force: bool = False) -> None:
        """Outer-loop Vorabrechnung refresh for thesis-aligned global NR."""
        from src.solvers.vorabrechnung import refresh_vorabrechnung_for_cells
        from src.solvers.vorabrechnung.vorab_policy import (
            thesis_vorab_sources_single_shot_enabled,
            vorab_outer_may_refresh_drying_pyro_sources,
        )

        # Explicit freeboard cells are part of the NR unknown vector in thesis mode.
        # Re-syncing them from closure during outer refresh overwrites NR-updated states
        # and can trigger large residual jumps between outer iterations.
        if not self._use_explicit_freeboard_solver_graph():
            self._refresh_explicit_freeboard_transport_from_closure()
        # Side-block states are seeded once during init/precalc.
        # Re-seeding each outer loop overwrites NR-updated side states and breaks Check2.
        strict_single_shot = thesis_vorab_sources_single_shot_enabled(self.config)
        refresh_sources = vorab_outer_may_refresh_drying_pyro_sources(self.config)
        force_sources = bool(force) and refresh_sources
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
        from src.solvers.vorabrechnung import refresh_cells_hydrodynamics_for_frozen_inner

        refresh_cells_hydrodynamics_for_frozen_inner(
            self._solver_cells_for_nr(),
            cfg=self.config,
        )

    def _set_outer_abgleich_context(self, *, refresh_skipped: bool = False) -> None:
        """Publish Check1-continuation context for inner |ΔT| budget selection."""
        self._nr_outer_refresh_skipped = bool(refresh_skipped)

    def _initialize_thesis_single_shot_vorab_sources(self) -> None:
        """Initialize fixed drying/DAEM source caches once using bed-mean temperature."""
        if not (bool(self.config.thesis_mode) and bool(self.config.thesis_vorab_sources_single_shot)):
            return
        if self._thesis_fixed_vorab_sources_ready:
            return
        from src.solvers.vorabrechnung import _bed_cell_indices, initialize_fixed_vorabrechnung_sources_for_cells
        from src.solvers.vorabrechnung.vorab_policy import freeze_vorab_drying_pyro_sources_for_cells

        T_bed_avg = float(np.mean([float(c.T) for c in self.cells])) if self.cells else float(self.config.T_inlet)
        n_zone = int(getattr(self.config, "vorab_vm_devolatilization_zone_cells_thesis", 0))
        zone_indices = _bed_cell_indices(self.cells)[:n_zone] if n_zone > 0 else []
        initialize_fixed_vorabrechnung_sources_for_cells(
            self.cells,
            T_reference=T_bed_avg,
            per_cell_temperature=n_zone > 0,
            vm_zone_cell_indices=zone_indices,
        )
        freeze_vorab_drying_pyro_sources_for_cells(self.cells)
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
        skip_precalc: bool = False,
    ) -> dict:
        from src.solvers.convergence import InnerConvergence, OuterMeritEarlyStop
        from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
        from src.workflow.steps.nr_inner_step import build_global_nr_inner_solve_fn
        from src.workflow.steps.outer_abgleich_step import run_outer_abgleich_for_global_nr
        from src.solvers.vorabrechnung.vorab_policy import vorab_refresh_policy_label

        cfg = self.config
        if gs_warmup_steps is not None and int(max(gs_warmup_steps, 0)) > 0:
            raise ValueError("gs_warmup_steps is not supported under NR-only policy.")
        from src.solvers.vorabrechnung.vorab_policy import clear_vorab_drying_pyro_source_freeze_for_cells

        self._thesis_fixed_vorab_sources_ready = False
        clear_vorab_drying_pyro_source_freeze_for_cells(self._solver_cells_for_nr())
        # 每次 solve 清掉 CO 带宽冻结运行时旗标（§5.13）
        if hasattr(cfg, "_upper_co_bridge_frozen_runtime"):
            try:
                delattr(cfg, "_upper_co_bridge_frozen_runtime")
            except Exception:
                setattr(cfg, "_upper_co_bridge_frozen_runtime", False)
        if hasattr(cfg, "_upper_co_local2_accept_count_runtime"):
            try:
                delattr(cfg, "_upper_co_local2_accept_count_runtime")
            except Exception:
                setattr(cfg, "_upper_co_local2_accept_count_runtime", 0)
        solver_cells = self._solver_cells_for_nr()
        from src.solvers.global_nr_solver import (
            install_bed_top_solid_energy_passthrough_on_cells,
            install_trace_o2_rate_gate_on_cells,
            install_trace_o2_snap_on_cells,
        )

        install_trace_o2_snap_on_cells(solver_cells, cfg)
        install_trace_o2_rate_gate_on_cells(solver_cells, cfg)
        install_bed_top_solid_energy_passthrough_on_cells(self.cells, cfg)
        solve_started = perf_counter()

        if skip_precalc:
            resolved_init_strategy = _resolve_nr_init_strategy(init_strategy)
            nr_init_s_total = 0.0
            nr_vorabrechnung_s = 0.0
        else:
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
        tol_rms = float(
            getattr(cfg, "nr_inner_tol_rms_thesis", None)
            if getattr(cfg, "nr_inner_tol_rms_thesis", None) is not None
            else np.clip(0.01 * max(tol, 1.0), 0.005, 0.02)
        )
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

            def _pre_inner_bridge(*, check1_continuation: bool = False) -> dict[str, bool | float]:
                from src.core.connectivity import run_outer_preinner_start_value_bridge
                from src.solvers.global_nr_solver import (
                    arm_bed_top_solid_energy_passthrough_if_y_co,
                    arm_trace_o2_rate_gate_if_y_co,
                )

                # §5.15：CO 带宽后武装微量 O₂ 速率门控（先抬 CO，再压 bed2 速率面）
                gate_arm = arm_trace_o2_rate_gate_if_y_co(solver_cells, cfg)
                # §5.20：CO 带宽后武装床顶固体焓透传
                pt_arm = arm_bed_top_solid_energy_passthrough_if_y_co(self.cells, cfg)
                diag = run_outer_preinner_start_value_bridge(
                    solver_cells,
                    cfg,
                    apply_bc_fn=self._apply_all_bc_for_nr,
                    apply_local_bc_fn=self._apply_local_bc_for_nr,
                    top_above_cell=self.freeboard_cells[0] if self.freeboard_cells else None,
                    check1_continuation=bool(check1_continuation),
                )
                diag = dict(diag or {})
                diag["trace_o2_rate_gate_armed"] = bool(gate_arm.get("armed", False))
                diag["trace_o2_rate_gate_newly_armed"] = bool(gate_arm.get("newly_armed", False))
                diag["trace_o2_rate_gate_y_co_top"] = float(gate_arm.get("y_co_top", float("nan")))
                diag["bed_top_energy_passthrough_armed"] = bool(pt_arm.get("armed", False))
                diag["bed_top_energy_passthrough_newly_armed"] = bool(pt_arm.get("newly_armed", False))
                diag["bed_top_energy_passthrough_y_co_top"] = float(pt_arm.get("y_co_top", float("nan")))
                self._outer_preinner_last_diag = diag
                return diag

            def _post_inner_holdup() -> bool:
                changed = False
                if bool(getattr(cfg, "nr_outer_postinner_holdup_reconcile_thesis", True)):
                    if self._reconcile_upper_bed_solid_holdup_for_nr():
                        changed = True
                if bool(cfg.thesis_mode) and bool(
                    getattr(cfg, "nr_outer_postinner_freeboard_energy_T_projection_thesis", True)
                ):
                    from src.core.freeboard_bridge import (
                        project_explicit_freeboard_temperatures_to_energy_closure,
                    )

                    fb_diag = project_explicit_freeboard_temperatures_to_energy_closure(self)
                    self._freeboard_energy_T_projection_diag = fb_diag
                    if int(fb_diag.get("n_projected", 0)) > 0:
                        changed = True
                if changed:
                    self._apply_all_bc_for_nr()
                return changed

            pre_inner_enabled = bool(
                getattr(cfg, "nr_outer_preinner_phase_split_preprojection_thesis", True)
                or getattr(cfg, "nr_outer_preinner_solid_energy_preprojection_thesis", True)
            )
            from src.solvers.vorab_cell_abgleich import (
                measure_outer_abgleich_gas_metrics,
                outer_convergence_from_config,
                run_outer_vorab_cell_abgleich_bridge,
            )

            outer_checker = outer_convergence_from_config(cfg)
            frozen_budgets = getattr(self, "_vorab_a_tier_cell_budgets", None)

            def _outer_gas_metrics_fn() -> dict[str, float]:
                return measure_outer_abgleich_gas_metrics(
                    solver_cells,
                    cfg,
                    apply_bc_fn=self._apply_all_bc_for_nr,
                )

            def _outer_check2_bridge_fn() -> dict[str, bool | float]:
                if not frozen_budgets:
                    return {"enabled": False, "aligned": False, "bridge_passes": 0.0}
                return run_outer_vorab_cell_abgleich_bridge(
                    solver_cells,
                    frozen_budgets,
                    cfg,
                    apply_bc_fn=self._apply_all_bc_for_nr,
                    reseed_holdup_chain_fn=self._seed_initialized_holdup_chain_for_nr,
                )

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
                pre_inner_fn=_pre_inner_bridge if pre_inner_enabled else None,
                post_inner_holdup_fn=_post_inner_holdup,
                pre_inner_on_continuation=bool(
                    getattr(cfg, "nr_outer_preinner_on_check1_continuation_thesis", True)
                ),
                outer_convergence=outer_checker,
                gas_metrics_fn=_outer_gas_metrics_fn,
                check2_bridge_fn=_outer_check2_bridge_fn,
                strict_check1_before_refresh=(
                    bool(cfg.thesis_mode) and bool(getattr(cfg, "nr_strict_check1_before_refresh_thesis", True))
                ),
                check1_energy_regress_refresh_ratio=(
                    float(getattr(cfg, "nr_check1_energy_regress_refresh_ratio_thesis", 0.0))
                    if bool(cfg.thesis_mode)
                    else 0.0
                ),
                outer_context_fn=self._set_outer_abgleich_context,
                merit_early_stop=OuterMeritEarlyStop(
                    enabled=bool(cfg.thesis_mode)
                    and bool(getattr(cfg, "nr_outer_merit_early_stop_thesis", True)),
                    min_outer_iters=int(
                        getattr(cfg, "nr_outer_merit_early_stop_min_iters_thesis", 8)
                    ),
                    regression_rel_tol=float(
                        getattr(cfg, "nr_outer_merit_early_stop_regression_rel_thesis", 0.25)
                    ),
                    regression_patience=int(
                        getattr(cfg, "nr_outer_merit_early_stop_regression_patience_thesis", 1)
                    ),
                    plateau_rel_tol=float(
                        getattr(cfg, "nr_outer_merit_early_stop_plateau_rel_thesis", 0.002)
                    ),
                    use_gas_phase_split_merit=bool(
                        getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)
                    ),
                    use_energy_merit=bool(
                        getattr(cfg, "nr_line_search_energy_merit_thesis", False)
                    ),
                ),
            )
            for item in result.outer_history:
                item["reaction_continuation_stage"] = int(stage_idx)
                item["reaction_rate_multiplier"] = float(stage_multiplier)
            return result

        def _y_co_top_for_hot_restart() -> float:
            from src.core.species import GAS_SPECIES_INDEX

            bed = [c for c in solver_cells if str(getattr(c, "cell_type", "bed")) == "bed"]
            if not bed:
                return float("nan")
            top = bed[-1]
            j_co = int(GAS_SPECIES_INDEX["CO"])
            n_co = float(top.N_d[j_co] + top.N_b[j_co])
            n_gas = float(sum(top.N_d) + sum(top.N_b))
            return float(n_co / max(n_gas, 1.0e-30))

        outer_results = []
        # skip_precalc 续算本身不再套娃热重启
        hot_restart_enabled = bool(
            getattr(cfg, "nr_outer_hot_restart_after_y_co_thesis", False)
        ) and (not bool(skip_precalc))
        hot_restart_meta: dict[str, object] = {
            "enabled": bool(hot_restart_enabled),
            "triggered": False,
            "y_co_top": float("nan"),
            "threshold": float(getattr(cfg, "nr_outer_hot_restart_y_co_threshold_thesis", 0.14)),
            "mode": "nested_skip_precalc",
        }
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
        preinner_s_total = 0.0
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
            preinner_s_total += float(stage_result.preinner_s_total)
            used_inner_budget += int(stage_result.used_inner_budget)
            outer_iters += int(stage_result.outer_iters)

        last_nr_result = outer_result.last_nr_result
        outer_converged = outer_result.converged_outer

        # §5.14：首轮结束后若 CO 进带宽，嵌套一次 skip_precalc 续算（完整重建 inner 闭包）。
        # 同闭包二次 `_run_outer` / 中途门禁复位均证伪（bed2 不愈）。
        if bool(hot_restart_enabled) and (not bool(outer_converged)):
            y_co_now = float(_y_co_top_for_hot_restart())
            thr = float(hot_restart_meta["threshold"])
            hot_restart_meta["y_co_top"] = y_co_now
            if np.isfinite(y_co_now) and y_co_now >= thr:
                bonus_inner = int(
                    max(getattr(cfg, "nr_outer_hot_restart_bonus_inner_budget_thesis", 0) or 0, 0)
                )
                cont_max_iter = int(max_iter)
                if bonus_inner > 0:
                    cont_max_iter = int(max(cont_max_iter, (bonus_inner + 4) // 5))
                hot_restart_meta["triggered"] = True
                hot_restart_meta["cont_max_iter"] = int(cont_max_iter)
                hot_restart_meta["y_co_top_at_trigger"] = float(y_co_now)
                # 对齐显式两段 oracle：首轮结束后先做 finalize 自由板 T 投影再续算
                if bool(cfg.thesis_mode) and bool(
                    getattr(cfg, "nr_finalize_freeboard_energy_T_projection_thesis", True)
                ):
                    from src.core.freeboard_bridge import (
                        project_explicit_freeboard_temperatures_to_energy_closure,
                    )

                    self._freeboard_energy_T_projection_diag = (
                        project_explicit_freeboard_temperatures_to_energy_closure(self)
                    )
                    hot_restart_meta["pre_cont_freeboard_T_projection"] = True
                stage1_history = list(outer_history)
                stage1_lambdas = list(all_nr_lambdas)
                stage1_trials = list(all_nr_line_search_trials)
                stage1_clips = list(all_nr_clip_history)
                stage1_norms = list(all_nr_history)
                stage1_outer_iters = int(outer_iters)
                stage1_inner_used = int(used_inner_budget)
                stage2 = self._solve_global_nr(
                    max_iter=int(cont_max_iter),
                    tol=float(tol),
                    verbose=bool(verbose),
                    init_strategy=resolved_init_strategy,
                    jacobian_strategy=jacobian_strategy,
                    linear_solver_backend=linear_solver_backend,
                    jacobian_lag_steps=jacobian_lag_steps,
                    skip_precalc=True,
                )
                assert isinstance(stage2, dict)
                h2 = list(stage2.get("nr_outer_history") or [])
                for item in h2:
                    item["outer_hot_restarted"] = True
                    item["outer_hot_restart_stage"] = 2
                    item["y_co_top_at_hot_restart_trigger"] = float(y_co_now)
                stage2["nr_outer_history"] = stage1_history + h2
                stage2["nr_accepted_lambda_history"] = stage1_lambdas + list(
                    stage2.get("nr_accepted_lambda_history") or []
                )
                stage2["nr_line_search_trial_counts"] = stage1_trials + list(
                    stage2.get("nr_line_search_trial_counts") or []
                )
                stage2["nr_clip_history"] = stage1_clips + list(
                    stage2.get("nr_clip_history") or []
                )
                stage2["norm_history"] = stage1_norms + list(stage2.get("norm_history") or [])
                stage2["nr_outer_iters"] = int(stage1_outer_iters) + int(
                    stage2.get("nr_outer_iters") or 0
                )
                stage2["nr_inner_budget_used"] = int(stage1_inner_used) + int(
                    stage2.get("nr_inner_budget_used") or 0
                )
                stage2["n_iter"] = len(list(stage2.get("norm_history") or []))
                stage2["nr_outer_hot_restart"] = dict(hot_restart_meta)
                stage2["nr_init_s_total"] = float(nr_init_s_total)
                stage2["nr_vorabrechnung_s"] = float(nr_vorabrechnung_s)
                stage2["nr_total_s"] = float(perf_counter() - solve_started)
                return stage2

        from src.solvers.global_nr_solver import compute_live_nr_scaled_metrics

        ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
        ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
        ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
        live_metrics = compute_live_nr_scaled_metrics(
            solver_cells,
            self._apply_all_bc_for_nr,
            ref_gas_mol_s=ref_gas,
            ref_solid_kg_s=ref_solid,
            ref_energy_W=ref_energy,
            use_gas_phase_split_merit=bool(
                getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)
                and getattr(cfg, "nr_enforce_gas_phase_split_convergence_thesis", False)
            ),
            use_energy_merit=bool(
                getattr(cfg, "nr_line_search_energy_merit_thesis", False)
                and getattr(cfg, "nr_enforce_gas_phase_split_convergence_thesis", False)
            ),
        )
        last_nr_result = dict(last_nr_result, **live_metrics)
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
                nr_vorabrechnung_policy=vorab_refresh_policy_label(cfg),
                nr_gs_warmup_steps=0,
                nr_outer_max=outer_max,
                nr_outer_iters=outer_iters,
                nr_merit_early_stop=bool(outer_result.merit_early_stop),
                nr_merit_early_stop_reason=str(outer_result.merit_early_stop_reason),
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
                        if self._nr_jacobian_uses_structured_direct_linear_solver(jacobian_mode)
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
                nr_preinner_s_total=preinner_s_total,
                nr_outer_history=outer_history,
                nr_vorabrechnung_signatures=[str(item.get("vorabrechnung_signature", "")) for item in outer_history],
                nr_inner_budget_total=total_inner_budget,
                nr_inner_budget_used=used_inner_budget,
                nr_outer_hot_restart=dict(hot_restart_meta),
                nr_total_s=perf_counter() - solve_started,
                timing=agg_nr_timing or last_nr_result.get("timing"),
                counts=agg_nr_counts or last_nr_result.get("counts"),
            )
        )

    def _finalize_global_nr_result(self, nr_result: dict) -> dict:
        from src.solvers.result_builder import finalize_global_nr_result

        cfg = self.config
        if bool(cfg.thesis_mode) and bool(
            getattr(cfg, "nr_finalize_freeboard_energy_T_projection_thesis", True)
        ):
            from src.core.freeboard_bridge import (
                project_explicit_freeboard_temperatures_to_energy_closure,
            )

            # One-shot after NR: clean freeboard T/energy ghosts for profiles without
            # feeding colder freeboard enthalpy back into outer Newton / recycle.
            self._freeboard_energy_T_projection_diag = (
                project_explicit_freeboard_temperatures_to_energy_closure(self)
            )

        return finalize_global_nr_result(
            self,
            nr_result,
            resolve_axial_heat_loss_distribution_fn=_resolve_axial_heat_loss_distribution,
            cell_solid_outflow_component_fn=_cell_solid_outflow_component,
        )
