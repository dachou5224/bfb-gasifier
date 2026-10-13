"""Phase 2 freeboard-aware 与 VTT 回路 NR 验收 ``ReactorConfig`` builder。"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from src.core.constants import Rg
from src.core.reactor import ReactorConfig

from tests.validation.case_lu_loader import load_case_LU
from tests.validation.config_mutators import apply_case_secondary_stream
from tests.validation.config_phase1 import build_phase1_htw_lu_global_nr_reactor_config

VTT_MINIMAL_LOOP_SOLVE_KWARGS: dict[str, object] = {
    "max_iter": 40,
    "tol": 1.0,
    "init_strategy": "vorabrechnung",
}

# R12（thesis R6）开启时 staged-march 底温种子；关 R12 时若仍为此指纹则回退冷锚。
_R12_NR_TEMPERATURE_SEED_BOTTOM_K: float = 1150.0
_R12_OFF_NR_TEMPERATURE_SEED_BOTTOM_K: float = 900.0


def apply_enable_r12_temperature_seed_policy(cfg: ReactorConfig) -> None:
    """使 bed0 NR 温度种子与 ``enable_r12`` 一致。

    冷锚 900 K + R12 会把 bed0 能量残差钉在 preinner/NR 无法离开的平台
    （outer RMS≈0.037 平坦）。开启 R12 时去掉 cold_anchor，并将底温种子抬到
    ≥1150 K（不降低已显式抬高的种子，如 audit seed1250）。

    Source: data/phase2_r12_seed_window_o8.json; Hamel Eq.5.38–5.39
    """
    if not bool(getattr(cfg, "vorab_staged_primary_gas_injection_thesis", False)):
        return
    bottom = float(getattr(cfg, "vorab_nr_temperature_seed_bottom_K_thesis", 900.0))
    cold = bool(getattr(cfg, "vorab_nr_temperature_seed_bed0_cold_anchor_thesis", False))
    if bool(getattr(cfg, "enable_r12", False)):
        cfg.vorab_nr_temperature_seed_bed0_cold_anchor_thesis = False
        if bottom < _R12_NR_TEMPERATURE_SEED_BOTTOM_K:
            cfg.vorab_nr_temperature_seed_bottom_K_thesis = _R12_NR_TEMPERATURE_SEED_BOTTOM_K
        return
    # 关 R12：仅当仍是本策略写入的指纹时，恢复 staged-march 冷锚基线
    if (not cold) and abs(bottom - _R12_NR_TEMPERATURE_SEED_BOTTOM_K) < 1e-6:
        cfg.vorab_nr_temperature_seed_bed0_cold_anchor_thesis = True
        cfg.vorab_nr_temperature_seed_bottom_K_thesis = _R12_OFF_NR_TEMPERATURE_SEED_BOTTOM_K


def build_phase2_htw_lu_freeboard_reactor_config(
    case: dict | None = None,
    n_freeboard_cells: int = 8,
    *,
    enable_vorab_march_staged: bool = True,
    enable_vorab_cold_budget_march: bool = False,
) -> ReactorConfig:
    """Phase 2 LU freeboard-aware 开发口径。"""
    if case is None:
        case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.explicit_side_blocks_enabled = True
    cfg.H_freeboard = max(float(case["reactor_height_m"]) - float(case["bed_height_m"]), 0.0)
    cfg.n_freeboard_cells = int(max(n_freeboard_cells, 0))
    cfg.D_bed = float(case.get("bed_diameter_m", case["D_freeboard"]))
    cfg.freeboard_heat_loss_frac = 0.0
    cfg.freeboard_enabled_reactions = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R12")
    cfg.freeboard_trajectory_model = "analytical_wirsum"
    cfg.freeboard_u_gb_scale = 1.0
    cfg.freeboard_beta_a_scale = 1.0
    cfg.freeboard_velocity_sigma = 0.60
    cfg.freeboard_velocity_bins = 5
    cfg.freeboard_age_quadrature_bins = 8
    cfg.freeboard_age_quadrature_max_age = 0.98
    cfg.freeboard_cyclone_capture_char_frac = 0.90
    cfg.freeboard_cyclone_capture_ash_frac = 0.95
    cfg.freeboard_secondary_injection_xi = case.get("secondary_injection_xi")
    # Jacobian FD 复用：lag=6 W2 ~134 s vs lag=4 ~710 s；W1 lag4/6/7 均 ~53 s 收敛（data/w2_jacobian_lag_sweep.json, w1_jacobian_lag_sweep.json）
    cfg.nr_jacobian_lag_steps = 6
    cfg.nr_line_search_max_trials_thesis = 12
    # Check1：Phase2 HTW 在 tol=0.01 时 W2 需 ~74 outer；0.03 约 52 outer，出口组分/温度偏差可接受（见 data/w2_check1_tol_sweep.json）
    cfg.nr_inner_tol_rms_thesis = 0.03
    # Pre-inner：关 solid-energy 预投影可显著减 wall time（W2 ~668→507 s），出口仍 converged（见 data/w2_walltime_options_sweep.json）
    cfg.nr_outer_preinner_solid_energy_preprojection_thesis = False
    cfg.freeboard_trajectory_coeff_model = "exact_hamel"
    cfg.explicit_freeboard_solver_graph_enabled = True
    # Hamel Kap.2.1：Trocknung/Pyrolyse 按反应器高度分配到各 cell（非 bed0 瞬时释放）。
    # n_zone=3 与底部分级进气格数同量级；profile=linear 展开 DAEM 总量。
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 3
    cfg.vorab_vm_devolatilization_zone_profile_thesis = "linear"
    # Hamel 同构默认不缩放；历史 0.25 为防 breakthrough 的工程压低，非文献值
    cfg.r5_scale = 1.0
    # Hamel 同构：原始动力学直通 F(x)；工程限幅仅 opt-in
    cfg.extent_limiters_enabled_thesis = False
    # P1：Inner 固定 h（outer refresh + Inner 入口 prepare）；勿在接受步刷 hydro
    cfg.nr_refresh_hydrodynamics_on_accepted_step_thesis = False
    cfg.gas_inlet_split_strategy = "fixed"
    cfg.gas_inlet_dense_frac = 0.50
    if enable_vorab_march_staged:
        cfg.vorab_staged_primary_gas_injection_thesis = True
        cfg.vorab_staged_primary_gas_bed_cells_thesis = 3
        cfg.vorab_a_tier_bed_axial_primary_profile_thesis = True
        # post-staged H2O 透传默认关（§5.29：与 LSQ2+Ept 叠开存在 CO≈0.15/0.21
        # 双吸引子）。要冲 rms 时 opt-in + OMP=1 + PYTHONHASHSEED=0。
        cfg.vorab_nr_temperature_seed_mode_thesis = "budget_progress"
        # 900 K: lowers init/final energy floor vs 850 while keeping split ~0.047
        # (850→900 A/B: final E 0.155→0.093, maxE 1.06→0.64; Tin stays LU 293 K).
        cfg.vorab_nr_temperature_seed_bottom_K_thesis = 900.0
        cfg.vorab_nr_temperature_seed_top_K_thesis = 1120.0
        cfg.vorab_nr_temperature_seed_bed0_cold_anchor_thesis = True
        # Cap multi-iter |ΔT|; 350 K keeps warmer bed0 seed from outer-level ~600 K jumps
        # while still allowing the O1 heating that lowers split below the ~0.086 floor.
        cfg.nr_inner_abs_dT_budget_K_thesis = 350.0
        cfg.vorab_a_tier_init_abgleich_temperature_blend_thesis = 0.0
        cfg.vorab_staged_primary_gas_dense_frac_thesis = 0.50
        cfg.vorab_x0_holdup_mode_thesis = "transport_from_vorab"
        cfg.vorab_major_gibbs_x0 = False
        cfg.vorab_pyrolysis_gibbs_composition_thesis = True
        cfg.vorab_a_tier_product_holdup_source_closure_bed_cells_thesis = 0
    if enable_vorab_cold_budget_march:
        cfg.vorab_bed_temperature_march_thesis = True
    cfg.nr_init_bottom_major_phase_split_preprojection_thesis = True
    cfg.nr_init_bed1_phase_split_preprojection_thesis = True
    # 分级进气末端（默认 bed0–2）之后 bed3 N₂/H₂O 曾因 axial ``max(inflow,current)``
    # 保留过填 x0（handoff §5.26–5.28；vorab_cell_mapping 对无 O2 格透传 H2O/N2）。
    # outer preinner T2 覆盖到 bed2（max=2）。bed1 须显式 include：同构裸动力学
    # + 快氧化 x₀ 下 A/B 显示 bed1 T2 把 o8 split 0.070→0.037
    # （data/phase2_extent_raw_split_preinner_ab.json）。max3 抬能量残差；勿把
    # upper_tau 降到 0.05（会过早选中 bed2，回归变差）。
    if bool(getattr(cfg, "vorab_staged_primary_gas_injection_thesis", False)):
        # max=2：覆盖分级进气区；max9 虽降 split 但抬 energy（o8 E≈0.068）。
        # bed1 显式 include：同构路径 A/B o8 split 0.070→0.037
        # （data/phase2_extent_raw_split_preinner_ab.json / bed1_preinner_o8.json）。
        cfg.nr_outer_preinner_phase_split_max_bed_index_thesis = 2
        cfg.nr_outer_preinner_bed1_phase_split_thesis = True
    cfg = apply_case_secondary_stream(cfg, case)
    apply_enable_r12_temperature_seed_policy(cfg)
    return cfg


def build_vtt_minimal_loop_nr_acceptance_config(
    case: dict | None = None,
    *,
    n_bed: int = 3,
    n_freeboard: int = 1,
) -> ReactorConfig:
    """B 档 VTT 单回路 NR 验收口径：3 bed + 1 fb + cyclone + return。"""
    if case is None:
        case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.thesis_mode = True
    cfg.explicit_side_blocks_enabled = True
    cfg.n_cells = int(n_bed)
    cfg.n_freeboard_cells = int(n_freeboard)
    cfg.H_freeboard = max(float(case["reactor_height_m"]) - float(case["bed_height_m"]), 0.0)
    cfg.explicit_freeboard_solver_graph_enabled = True
    cfg.recycle_gas = True
    cfg.recirculation_frac = 0.1
    cfg.freeboard_cyclone_capture_char_frac = 0.90
    cfg.freeboard_cyclone_capture_ash_frac = 0.95
    cfg.vorab_major_gibbs_x0 = False
    cfg.major_gibbs_solver_mode = "hamel_reduced"
    cfg.nr_jacobian_lag_steps = 2
    cfg.nr_inner_tol_rms_thesis = 0.02
    cfg.nr_outer_preinner_bed1_phase_split_thesis = True
    cfg.vorab_a_tier_bed1_phase_split_after_stream_mapping_thesis = True
    graph_path = Path(__file__).resolve().parents[2] / "data" / "connectivity" / "vtt_minimal_loop.json"
    cfg.connectivity_graph_path = str(graph_path)
    # Hybrid is default (``nr_hybrid_jacobian_thesis=True``); set False for pure structured FD.
    return cfg


def build_phase2_htw_lu_local_refined_damped_config(
    case: dict | None = None,
    *,
    n_freeboard_cells: int = 8,
    reactive_zone_height_m: float = 1.0,
    tau_cell_target_s: float = 0.20,
    max_refine_cells: int = 10,
    thesis_damp_halvings: int = 18,
    thesis_lambda_min: float = 1.0 / 16384.0,
) -> ReactorConfig:
    """Phase 2: 床层入口局部加密 + thesis 阻尼联动配置。"""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case=case, n_freeboard_cells=n_freeboard_cells)
    n_base = int(cfg.n_cells)
    dh_base = float(cfg.H_bed) / float(max(n_base, 1))

    A_bed = 0.25 * np.pi * float(cfg.D_bed) ** 2
    n_tot_in = float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed)
    u_eff = max(
        (n_tot_in * float(Rg) * float(cfg.T_inlet) / max(float(cfg.P), 1e-9)) / max(A_bed, 1e-12),
        0.05,
    )
    dh_tau = max(float(tau_cell_target_s) * u_eff, 0.05)
    dh_fine = min(0.5 * dh_base, dh_tau)

    h_reactive = float(np.clip(reactive_zone_height_m, 0.0, float(cfg.H_bed) - 1e-6))
    n_fine = int(max(1, math.ceil(h_reactive / max(dh_fine, 1e-9))))
    n_fine = int(min(n_fine, max(int(max_refine_cells), 1)))
    h_fine = n_fine * dh_fine
    if h_fine >= float(cfg.H_bed):
        n_fine = max(1, int(float(cfg.H_bed) / max(dh_fine, 1e-9)) - 1)
        h_fine = n_fine * dh_fine

    h_rem = float(cfg.H_bed) - h_fine
    n_upper = int(max(1, math.ceil(h_rem / max(dh_base, 1e-9))))
    dh_upper = h_rem / float(n_upper)

    bed_profile = tuple([float(dh_fine)] * n_fine + [float(dh_upper)] * n_upper)
    cfg.bed_dh_profile = bed_profile
    cfg.n_cells = len(bed_profile)
    cfg.n_freeboard_cells = max(1, int(round(float(n_freeboard_cells) * cfg.n_cells / max(n_base, 1))))

    cfg.nr_damping_halvings_thesis = int(max(thesis_damp_halvings, 1))
    cfg.nr_lambda_min_thesis = float(max(thesis_lambda_min, 1e-12))
    cfg.nr_prefer_full_step_thesis = True
    return cfg
