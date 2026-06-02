"""Table 2 LU 工况 — 全局 NR 验证脚本。

使用 ``global_nr`` 求解器运行 LU 工况，并对比 validation_cases.json。
"""

from __future__ import annotations

import os
import sys
import time
import numpy as np
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

try:
    import pytest
except ImportError:
    class _PytestShim:
        class mark:
            @staticmethod
            def slow(f): return f
    pytest = _PytestShim()

from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from src.core.reactor import Reactor
from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0

from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    apply_case_secondary_stream,
    apply_secondary_air_repartition,
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    extract_case_inlet_topology_hints,
    load_case_LU,
    load_validation_case_node,
    json_numeric_or_none
)

def run_case_LU_global_nr(
    *,
    max_global_iter: int | None = None,
    tol_global: float | None = None,
    recirculation_frac: float = 0.8, 
    n_cells: int = 5,
    solver: str | None = None,
    heat_loss_frac: float = 0.1,
    enable_r12: bool = True,
    u0_target: float | None = None,
) -> dict:
    case = load_case_LU()

    # 与 tests/test_table2_LU.py 使用同一 Phase-1 基线，避免 GS/NR 参数漂移
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.n_cells = n_cells
    cfg.recirculation_frac = recirculation_frac
    cfg.heat_loss_frac = heat_loss_frac
    cfg.enable_r12 = enable_r12
    if u0_target is not None:
        cfg.u0_target = u0_target

    print("=" * 60)
    print(f"BFB Gasifier - Table 2 LU (Solver: {solver})")
    print("=" * 60)
    print(f"燃料: {case['fuel']} @ {case['fuel_feed']} kg/h")
    print(f"压力: {case['P']/1e6:.1f} MPa | ER: {case['ER']}")
    print(f"床层高度: {case['bed_height_m']} m | 循环率: {recirculation_frac}")
    print(f"计算进料: O2={cfg.O2_feed:.2f}, H2O={cfg.H2O_feed:.2f}, N2={cfg.N2_feed:.2f} mol/s")
    print()

    reactor = Reactor(cfg)
    
    T_prof_init = estimate_axial_T_profile(
        n_cells=cfg.n_cells,
        T_inlet=cfg.T_inlet,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        moisture_wt=case["moisture_wt"],
        P=case["P"],
    )
    generate_initial_x0(
        cells=reactor.cells,
        O2_feed=cfg.O2_feed,
        H2O_feed=cfg.H2O_feed,
        N2_feed=cfg.N2_feed,
        fuel_feed_kg_s=cfg.fuel_feed,
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        moisture_wt=case["moisture_wt"],
        ash_dry_wt=case["ash_dry_wt"],
        VM_daf=case["VM_daf"],
        T_profile=T_prof_init,
        fuel_type="coal",
    )

    t0 = time.perf_counter()
    solve_kwargs = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    if max_global_iter is not None:
        solve_kwargs["max_global_iter"] = max_global_iter
    if tol_global is not None:
        solve_kwargs["tol_global"] = tol_global
    if solver is not None:
        solve_kwargs["solver"] = solver
    result = reactor.solve(
        **solve_kwargs,
        verbose=True,
    )
    elapsed = time.perf_counter() - t0

    T_exit = result["T_profile"][-1]
    carbon_conv = result["carbon_conv"]
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outs = raw["outputs"]
    T_json = float(outs["exit_temperature_K"])
    y_dry = result["exit_gas_dry"]
    exit_gas = result["exit_gas"]
    dry_ref = outs.get("exit_gas_dry_mol_frac", {})

    print(f"\n wall time: {elapsed:.1f} s | 迭代次数: {result['n_iter']}")
    print(f" 最终收敛: {result.get('converged_fully', result.get('converged'))}")
    print(f" NR init: {result.get('nr_init_strategy')}")
    print(f" NR warmup: {result.get('nr_gs_warmup_steps')}")
    print(f" NR Jacobian: {result.get('nr_jacobian_strategy')}")
    print(f" 碳转化率: {carbon_conv*100.0:.1f}% (JSON: {outs.get('carbon_conversion_pct', 'N/A')}%)")
    print(f" 出口温度: {T_exit:.1f} K (JSON: {T_json} K)")
    
    print("\n湿基组成 (exit_gas):")
    for sp in GAS_SPECIES:
        y = exit_gas.get(sp, 0.0)
        print(f"  {sp:6s}: {y:.4f}")

    print("\n[DEBUG] 顶格原始摩尔流 [mol/s]:")
    top = reactor.cells[-1]
    for i, sp in enumerate(GAS_SPECIES):
        print(f"  {sp:6s}: N_d={top.N_d[i]:.4e}, N_b={top.N_b[i]:.4e}, Sum={top.N_d[i]+top.N_b[i]:.4e}")

    print("\n干基组成 (模拟 vs JSON):")
    for sp in ("CO", "CO2", "H2", "CH4"):
        t = json_numeric_or_none(dry_ref.get(sp))
        s = y_dry.get(sp, 0.0)
        err = abs(s - t) / t if t else 0.0
        print(f"  {sp:4s}: sim={s:.4f} json={t if t else 0.0:.4f} err={err:.2%}")
    
    return result


@pytest.mark.slow
def test_phase1_lu_global_nr_shared_policy_reports_explicit_init_strategy():
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    result = Reactor(cfg).solve(**PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)

    assert result["nr_init_strategy"] == PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS["nr_init_strategy"]
    assert result["nr_gs_warmup_steps"] == 0
    assert result["nr_jacobian_strategy"] == PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS["nr_jacobian_strategy"]
    assert result["nr_outer_iters"] >= 1
    assert result["nr_outer_max"] >= result["nr_outer_iters"]
    assert result["nr_init_s_total"] >= result["nr_vorabrechnung_s"] >= 0.0
    assert result["nr_timing"]["jacobian_build_s"] >= 0.0
    assert result["nr_counts"]["global_residual_calls"] >= 2
    assert 500.0 < float(result["T_profile"][-1]) < 2000.0
    assert 0.0 <= float(result["carbon_conv"]) <= 1.0
    if result.get("converged_outer", result.get("converged", False)):
        assert float(result["carbon_conv"]) > 0.0


def test_phase1_lu_global_nr_shared_config_uses_separate_dense_fraction():
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    assert cfg.gas_inlet_dense_frac == pytest.approx(0.20)


def test_phase1_lu_global_nr_shared_config_uses_hamel_wein_hydrodynamics_chain():
    cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    assert cfg.thesis_mode is True
    assert cfg.hydrodynamics_u_d_closure == "wein_1992_eq312"
    assert cfg.hydrodynamics_bubble_diameter_model == "hilligardt_ode"
    assert cfg.hydrodynamics_psi_b_strategy == "wein_1992"
    assert cfg.hydrodynamics_lambda_strategy == "hamel_280"
    assert cfg.hydrodynamics_xi_strategy == "hamel_regime"
    assert cfg.hydrodynamics_bubble_velocity_strategy == "heinbockel_eq343"
    assert cfg.hydrodynamics_bubble_ode_strategy == "heinbockel_eq341"
    assert cfg.use_gibbs_minor is True


def test_load_case_lu_exposes_secondary_and_recycle_hints_from_json():
    case = load_case_LU()
    assert case["recirculation"] is True
    assert case["cyclone_present"] is True
    assert case["recycle_topology"] == "cyclone return to bed"
    assert case["bed_diameter_m"] == pytest.approx(0.6)
    assert case["freeboard_diameter_m"] == pytest.approx(0.6)
    assert case["total_air_Nm3_h"] == pytest.approx(5536.9)
    assert case["secondary_air_Nm3_h"] == pytest.approx(30.1)
    assert case["primary_air_Nm3_h"] == pytest.approx(5536.9 - 30.1)
    assert case["secondary_injection_xi"] == pytest.approx(0.60)
    assert case["secondary_injection_agent"] == "Air"


def test_extract_case_inlet_topology_hints_converts_numeric_secondary_air():
    node = load_validation_case_node("CASE_VTT_PRESSURIZED_PEAT_13")
    hints = extract_case_inlet_topology_hints(node)
    assert hints["secondary_injection_xi"] is None
    assert hints["secondary_air_Nm3_h"] == pytest.approx(41.6)
    assert hints["secondary_air_O2_mol_s"] == pytest.approx(0.21 * 41.6 * (1000.0 / 22.414) / 3600.0)
    assert hints["secondary_air_N2_mol_s"] == pytest.approx(0.79 * 41.6 * (1000.0 / 22.414) / 3600.0)
    assert hints["recirculation"] is True
    assert hints["cyclone_present"] is True
    assert hints["freeboard_diameter_m"] == pytest.approx(0.25)


def test_extract_case_inlet_topology_hints_supports_wesseling_secondary_o2():
    node = load_validation_case_node("CASE_HTW_WESSELING_2")
    hints = extract_case_inlet_topology_hints(node)
    o2_mol_s = 37.5 * (1000.0 / 22.414) / 3600.0
    assert hints["total_O2_Nm3_h"] == pytest.approx(1921.6)
    assert hints["secondary_O2_Nm3_h"] == pytest.approx(37.5)
    assert hints["primary_O2_Nm3_h"] == pytest.approx(1921.6 - 37.5)
    assert hints["secondary_injection_agent"] == "O2"
    assert hints["secondary_agent_Nm3_h"] == pytest.approx(37.5)
    assert hints["secondary_agent_O2_mol_s"] == pytest.approx(o2_mol_s)
    assert hints["secondary_agent_N2_mol_s"] == pytest.approx(0.0)


def test_extract_case_inlet_topology_hints_resolves_wsv_shared_reactor_block():
    node = load_validation_case_node("CASE_WSV400_B")
    hints = extract_case_inlet_topology_hints(node)
    assert hints["bed_diameter_m"] == pytest.approx(0.4)
    assert hints["freeboard_diameter_m"] == pytest.approx(0.6)
    assert hints["reactor_height_m"] == pytest.approx(3.0)
    assert hints["fuel_inlet_mm_above_nozzle"] == pytest.approx(200.0)
    assert hints["recirculation"] is False


def test_apply_secondary_air_repartition_preserves_total_air_budget():
    cfg = build_phase2_htw_lu_freeboard_reactor_config(load_case_LU())
    o2_total = cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s
    n2_total = cfg.N2_feed + cfg.freeboard_secondary_N2_mol_s
    apply_secondary_air_repartition(cfg, secondary_air_frac=0.15, injection_xi=0.60)
    assert cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s == pytest.approx(o2_total)
    assert cfg.N2_feed + cfg.freeboard_secondary_N2_mol_s == pytest.approx(n2_total)
    assert cfg.freeboard_secondary_H2O_mol_s == pytest.approx(0.0)
    assert cfg.freeboard_secondary_injection_xi == pytest.approx(0.60)


def test_apply_case_secondary_stream_uses_explicit_case_budget():
    case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    o2_total = cfg.O2_feed
    n2_total = cfg.N2_feed
    apply_case_secondary_stream(cfg, case)
    assert cfg.freeboard_secondary_O2_mol_s == pytest.approx(case["secondary_agent_O2_mol_s"])
    assert cfg.freeboard_secondary_N2_mol_s == pytest.approx(case["secondary_agent_N2_mol_s"])
    assert cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s == pytest.approx(o2_total)
    assert cfg.N2_feed + cfg.freeboard_secondary_N2_mol_s == pytest.approx(n2_total)
    assert cfg.freeboard_secondary_injection_xi == pytest.approx(case["secondary_injection_xi"])


def test_phase2_lu_freeboard_config_restores_reactor_height_segment():
    case = load_case_LU()
    cfg = build_phase2_htw_lu_freeboard_reactor_config(case)
    assert cfg.H_freeboard == pytest.approx(case["reactor_height_m"] - case["bed_height_m"])
    assert cfg.D_bed == pytest.approx(case["bed_diameter_m"])
    assert cfg.n_freeboard_cells > 0
    assert cfg.freeboard_secondary_injection_xi == pytest.approx(case["secondary_injection_xi"])
    assert cfg.freeboard_enabled_reactions == ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R12")
    assert cfg.freeboard_trajectory_model == "analytical_wirsum"
    assert cfg.freeboard_u_gb_scale == pytest.approx(1.0)
    assert cfg.freeboard_beta_a_scale == pytest.approx(1.0)
    assert cfg.freeboard_velocity_sigma == pytest.approx(0.60)
    assert cfg.freeboard_velocity_bins == 5
    assert cfg.freeboard_cyclone_capture_char_frac == pytest.approx(0.90)
    assert cfg.freeboard_cyclone_capture_ash_frac == pytest.approx(0.95)
    assert cfg.freeboard_secondary_O2_mol_s == pytest.approx(case["secondary_agent_O2_mol_s"])
    assert cfg.freeboard_secondary_N2_mol_s == pytest.approx(case["secondary_agent_N2_mol_s"])

if __name__ == "__main__":
    # 使用 global_nr 验证
    run_case_LU_global_nr()
