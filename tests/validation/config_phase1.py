"""Phase 1 HTW LU 基线与 global NR 开发口径 ``ReactorConfig`` builder。"""

from __future__ import annotations

import math

from src.core.reactor import ReactorConfig

from tests.validation.case_lu_loader import load_case_LU


def build_phase1_htw_lu_reactor_config(case: dict | None = None) -> ReactorConfig:
    """构建 Phase 1（Table 2 LU）与 ``test_table2_LU`` slow 用例一致的 ``ReactorConfig``。"""
    if case is None:
        case = load_case_LU()
    return ReactorConfig(
        n_age_classes=int(case.get("n_size_classes", 10)),
        n_cells=10,
        H_bed=case["H_bed"],
        H_freeboard=0.0,
        D_bed=case["bed_diameter_m"],
        P=case["P"],
        T_inlet=case["T_inlet"],
        fuel_type=str(case.get("fuel_type", "coal")),
        rho_s=1200.0,
        d_p=0.5 * sum(case.get("particle_diameter_mm_range", (1.5, 3.0))) * 1.0e-3,
        d_p_min=float(case.get("particle_diameter_mm_range", (1.5, 3.0))[0]) * 1.0e-3,
        d_p_max=float(case.get("particle_diameter_mm_range", (1.5, 3.0))[1]) * 1.0e-3,
        phi_s=0.75,
        eps_mf=0.45,
        fuel_feed=case["fuel_feed"] / 3600.0,
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
        moisture_wt=case["moisture_wt"],
        ash_dry_wt=case["ash_dry_wt"],
        VM_daf=case["VM_daf"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        nitrogen_fraction=case.get("N_dry", 0.0),
        use_gibbs_minor=False,
        gas_inlet_dense_frac=0.30,
        recirculation_frac=0.1 if case.get("recirculation", True) else 0.0,
        recycle_gas=bool(case.get("recirculation", True)),
        # Hamel Eq.5.38–5.39（thesis R6 H₂ 氧化）；仓库标签 R12。同构默认开启。
        enable_r12=True,
        top_solid_inlet_frac=0.50,
        solid_lower_inlet_frac=1.0,
        r4_scale=0.50,
        r5_scale=0.75,
        r7_scale=2.50,
    )


def build_phase1_htw_lu_global_nr_reactor_config(case: dict | None = None) -> ReactorConfig:
    """Phase 1 LU shared global NR 开发口径。"""
    cfg = build_phase1_htw_lu_reactor_config(case)
    cfg.gas_inlet_dense_frac = 0.20
    cfg.gas_inlet_split_strategy = "precalc_hydrodynamic_flux"
    cfg.hydrodynamics_u_d_closure = "wein_1992_eq312"
    cfg.hydrodynamics_bubble_diameter_model = "hilligardt_ode"
    cfg.hydrodynamics_psi_b_strategy = "wein_1992"
    cfg.hydrodynamics_lambda_strategy = "hamel_280"
    cfg.hydrodynamics_xi_strategy = "hamel_regime"
    cfg.hydrodynamics_bubble_velocity_strategy = "heinbockel_eq343"
    cfg.hydrodynamics_bubble_ode_strategy = "heinbockel_eq341"
    cfg.use_gibbs_minor = True
    cfg.thesis_mode = True
    cfg.explicit_side_blocks_enabled = False
    return cfg


def build_phase1_htw_lu_refined_config(
    case: dict | None = None,
    n_fine: int = 3,
    dh_fine: float = 0.15,
) -> ReactorConfig:
    """Phase 1 LU 配置 + 底部轴向加密。"""
    if case is None:
        case = load_case_LU()

    H_bed = float(case["H_bed"])
    H_remaining = H_bed - n_fine * dh_fine
    if H_remaining <= 0:
        raise ValueError(
            f"n_fine={n_fine} × dh_fine={dh_fine} = {n_fine*dh_fine:.2f} m "
            f"already exceeds H_bed={H_bed} m"
        )
    n_upper = max(1, math.ceil(H_remaining / 0.6))
    dh_upper = H_remaining / n_upper
    n_cells = n_fine + n_upper

    return ReactorConfig(
        n_age_classes=int(case.get("n_size_classes", 10)),
        n_cells=n_cells,
        H_bed=H_bed,
        H_freeboard=0.0,
        D_bed=case["bed_diameter_m"],
        P=case["P"],
        T_inlet=case["T_inlet"],
        fuel_type=str(case.get("fuel_type", "coal")),
        rho_s=1200.0,
        d_p=0.5 * sum(case.get("particle_diameter_mm_range", (1.5, 3.0))) * 1.0e-3,
        d_p_min=float(case.get("particle_diameter_mm_range", (1.5, 3.0))[0]) * 1.0e-3,
        d_p_max=float(case.get("particle_diameter_mm_range", (1.5, 3.0))[1]) * 1.0e-3,
        phi_s=0.75,
        eps_mf=0.45,
        fuel_feed=case["fuel_feed"] / 3600.0,
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
        moisture_wt=case["moisture_wt"],
        ash_dry_wt=case["ash_dry_wt"],
        VM_daf=case["VM_daf"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        nitrogen_fraction=case.get("N_dry", 0.0),
        gas_inlet_dense_frac=0.40,
        recirculation_frac=0.1 if case.get("recirculation", True) else 0.0,
        recycle_gas=bool(case.get("recirculation", True)),
        # Hamel Eq.5.38–5.39（thesis R6 H₂ 氧化）；仓库标签 R12。同构默认开启。
        enable_r12=True,
        top_solid_inlet_frac=0.50,
        r4_scale=0.75,
        r5_scale=0.50,
    )
