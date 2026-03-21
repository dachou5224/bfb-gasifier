"""LU 工况对比验证：Gauss-Seidel vs 全局 NR，对照 validation_cases.json。

运行：pytest tests/test_lu_gs_vs_global_nr.py -m slow --tb=short
或：python tests/test_lu_gs_vs_global_nr.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

# 允许 `python tests/test_lu_gs_vs_global_nr.py` 直接运行（与 test_table2_LU 一致）
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from src.core.reactor import Reactor, ReactorConfig

_REPO_ROOT = Path(__file__).resolve().parent.parent
_VALIDATION_JSON = _REPO_ROOT / "data" / "validation_cases.json"


@pytest.mark.slow
def test_lu_gauss_seidel_vs_global_nr_integration():
    data = json.loads(_VALIDATION_JSON.read_text(encoding="utf-8"))
    case = data["CASE_HTW_WESSELING_1"]
    outs = case["outputs"]
    T_ref = float(outs["exit_temperature_K"])
    XC_ref = float(outs["carbon_conversion_pct"])
    CO_ref = float(outs["exit_gas_dry_mol_frac"]["CO"])

    cfg = ReactorConfig(
        n_cells=10,
        H_bed=float(case["inputs"]["reactor"]["bed_height_m"]),
        D_bed=float(case["inputs"]["reactor"]["diameter_m"]),
        P=float(case["inputs"]["operating_conditions"]["pressure_MPa"]) * 1e6,
        T_inlet=float(case["inputs"]["operating_conditions"]["T_inlet_K"]),
        fuel_type="coal",
        rho_s=1000.0,
        d_p=2.25e-3,
        phi_s=0.86,
        eps_mf=0.45,
        fuel_feed=float(case["inputs"]["fuel"]["feed_rate_kg_h"]) / 3600.0,
        ER=float(case["inputs"]["operating_conditions"]["ER"]),
        primary_agent="air_steam",
        S_dry=float(case["inputs"]["fuel"]["ultimate_analysis_dry_wt_pct"]["S"]),
        moisture_wt=float(
            case["inputs"]["fuel"]["proximate_analysis"]["moisture_wt_pct"]
        ),
        ash_dry_wt=float(
            case["inputs"]["fuel"]["proximate_analysis"]["ash_dry_wt_pct"]
        ),
        VM_daf=float(case["inputs"]["fuel"]["proximate_analysis"]["VM_daf_pct"]),
        C_dry=float(case["inputs"]["fuel"]["ultimate_analysis_dry_wt_pct"]["C"]),
        H_dry=float(case["inputs"]["fuel"]["ultimate_analysis_dry_wt_pct"]["H"]),
        O_dry=float(case["inputs"]["fuel"]["ultimate_analysis_dry_wt_pct"]["O"]),
        recirculation_frac=0.1,
    )

    print("\n" + "=" * 60)
    print("Running Gauss-Seidel (baseline)...")
    r_gs = Reactor(cfg)
    out_gs = r_gs.solve(max_global_iter=15, tol_global=5.0, solver="gauss_seidel")
    T_gs = out_gs["T_profile"][-1]
    XC_gs = out_gs["carbon_conv"] * 100
    y_CO_gs = out_gs["exit_gas"].get("CO", 0.0)
    y_H2O_gs = out_gs["exit_gas"].get("H2O", 0.0)
    CO_dry_gs = y_CO_gs / max(1 - y_H2O_gs, 0.01)
    print(f"  T_exit = {T_gs:.1f} K  (target {T_ref:.1f} K)")
    print(f"  X_C    = {XC_gs:.1f}%  (target {XC_ref:.1f}%)")
    print(f"  CO dry = {CO_dry_gs:.4f}  (target {CO_ref:.4f})")

    print("\nRunning Global NR...")
    r_nr = Reactor(cfg)
    out_nr = r_nr.solve(max_global_iter=25, tol_global=1.0, solver="global_nr")
    T_nr = out_nr["T_profile"][-1]
    XC_nr = out_nr["carbon_conv"] * 100
    y_CO_nr = out_nr["exit_gas"].get("CO", 0.0)
    y_H2O_nr = out_nr["exit_gas"].get("H2O", 0.0)
    CO_dry_nr = y_CO_nr / max(1 - y_H2O_nr, 0.01)
    print(f"  T_exit = {T_nr:.1f} K  (target {T_ref:.1f} K)")
    print(f"  X_C    = {XC_nr:.1f}%  (target {XC_ref:.1f}%)")
    print(f"  CO dry = {CO_dry_nr:.4f}  (target {CO_ref:.4f})")
    print(
        f"  NR 外迭代={out_nr.get('converged_outer', out_nr['converged'])}  "
        f"内层={out_nr.get('converged_inner_nr')}  "
        f"iters={out_nr['n_iter']}"
    )
    if out_nr.get("rms_scaled_final") is not None:
        print(f"  RMS(‖F̂‖) 末次内层: {out_nr['rms_scaled_final']:.4e}")
    print(f"  Norm history: {[f'{v:.2e}' for v in out_nr['norm_history']]}")

    # 物理合理性
    assert 500 < T_nr < 2000, f"T_exit 非物理: {T_nr:.1f} K"
    assert 0 < XC_nr < 100, f"碳转化率非物理: {XC_nr:.1f}%"
    assert 0 < CO_dry_nr < 0.6, f"CO 非物理: {CO_dry_nr:.4f}"

    # 多外迭代拼接的 norm_history 不要求单调

    sanity_script = _REPO_ROOT / "tests" / "sanity_checks.py"
    result = subprocess.run(
        [sys.executable, str(sanity_script)],
        cwd=str(_REPO_ROOT),
        capture_output=True,
        text=True,
    )
    print("\n" + result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    assert result.returncode == 0, "sanity_checks 退出码非 0"
    assert "PASSED" in result.stdout, "sanity_checks 未全部通过"


if __name__ == "__main__":
    test_lu_gauss_seidel_vs_global_nr_integration()
    print("\n" + "=" * 60)
    print("STEP 4 PASS — integration test complete")
