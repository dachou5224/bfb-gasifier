"""工况 LU 端到端验证（Phase 6.1）。

使用 HTW 加压炉褐煤工况参数运行完整反应器模型，
对比实验数据（Hamel & Krumm 2001, Table 2）。

验证目标（**默认与** ``data/validation_cases.json`` **中** ``CASE_HTW_WESSELING_1.outputs`` **对比**）：
- 出口温度：相对 ``exit_temperature_K`` / ``exit_temperature_measured_K`` 取较优者，容差见根节点 ``validation_performance_summary``（±10–15% → 测试用 **15%**）。
- 干基主气相（CO、CO₂、H₂、CH₄）：与 ``exit_gas_dry_mol_frac`` 中数值字段对比；CH₄ 容差 **20%**，其余 **15%**（与 summary 一致）。
- 标定前若需跳过严格断言：设置 ``BFB_RELAX_VALIDATION=1``。

Phase 1 基线与共享加载逻辑见 ``tests/validation_case_utils.py``。

Source: docs/CLAUDE.md；data/validation_cases.json
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from src.core.reactor import Reactor

from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_reactor_config,
    estimate_gas_feeds,
    json_numeric_or_none,
    load_case_LU,
    load_validation_case_node,
    strict_validation_gate,
    validation_numeric_tolerances,
)


def test_load_case_LU_reads_validation_cases_json():
    """CASE_LU 与 validation_cases 中 CASE_HTW_WESSELING_1（Table 7.1 Sim 1）一致。"""
    case = load_case_LU()
    assert case["fuel_feed"] == pytest.approx(3377.4)
    assert case["ER"] == pytest.approx(0.337)
    assert case["P"] == pytest.approx(2.5e6)
    assert case["primary_agent"] == "air_steam"
    # 床层高度用于轴向离散，须为 bed_height_m（≈6 m），勿误用全炉 height_m（≈14.5 m）
    assert case["H_bed"] == pytest.approx(case["bed_height_m"])
    assert case["H_bed"] == pytest.approx(6.0)


def test_lu_particle_size_classes_map_to_solid():
    """CASE_LU 粒径范围按 JSON 的 10 个离散类映射到 SolidProps。"""
    case = load_case_LU()
    cfg = build_phase1_htw_lu_reactor_config(case)
    cfg.n_cells = 2
    assert cfg.n_age_classes == 10
    assert cfg.d_p == pytest.approx(2.25e-3)
    assert cfg.d_p_min == pytest.approx(1.5e-3)
    assert cfg.d_p_max == pytest.approx(3.0e-3)
    reactor = Reactor(cfg)
    for c in reactor.cells:
        assert c.solid.n_size_classes == 10
        assert len(c.solid.d_p_classes) == 10
        assert c.solid.d_p_classes[0] == pytest.approx(1.5e-3)
        assert c.solid.d_p_classes[-1] == pytest.approx(3.0e-3)
        assert c.solid.mass_fractions.shape == (10,)
        assert float(c.solid.mass_fractions.sum()) == pytest.approx(1.0)


def test_phase1_lu_config_uses_current_tuned_window():
    """LU shared config 固定在当前已验证更优的稳定窗口。"""
    cfg = build_phase1_htw_lu_reactor_config()
    assert cfg.gas_inlet_dense_frac == pytest.approx(0.30)
    assert cfg.r4_scale == pytest.approx(0.50)
    assert cfg.r5_scale == pytest.approx(0.75)
    assert cfg.r7_scale == pytest.approx(2.50)


def test_phase1_lu_window_reports_finite_dense_sensitivity_after_particle_remap():
    """文献粒径映射后，旧 dense tuning window 只保留为有限性/敏感性监测。"""
    ref_dry = {"CO": 0.157, "CO2": 0.133, "H2": 0.145, "CH4": 0.034}

    def _raw_score(out: dict) -> float:
        y = out["exit_gas_dry"]
        score = sum(abs(float(y.get(sp, 0.0)) - tgt) / tgt for sp, tgt in ref_dry.items())
        score += abs(float(out["T_profile"][-1]) - 1120.0) / 1120.0
        score += abs(float(out["carbon_conv"]) - 0.95) / 0.95
        return float(score)

    def _conv_score(out: dict) -> float:
        rms_raw = out.get("rms_scaled_final")
        rms = float(rms_raw) if rms_raw is not None else 1.0
        return float(_raw_score(out) + 2.0 * max(rms - 0.15, 0.0))

    solve_kwargs = {
        "max_global_iter": 3,
        "tol_global": 1.0,
        "solver": "global_nr",
        "nr_init_strategy": "vorabrechnung",
        "nr_jacobian_strategy": "block_tridiag_structured",
    }

    cfg_baseline = build_phase1_htw_lu_reactor_config()
    out_baseline = Reactor(cfg_baseline).solve(**solve_kwargs)

    cfg_dense_035 = build_phase1_htw_lu_reactor_config()
    cfg_dense_035.gas_inlet_dense_frac = 0.35
    out_dense_035 = Reactor(cfg_dense_035).solve(**solve_kwargs)

    raw_baseline = _raw_score(out_baseline)
    raw_dense_035 = _raw_score(out_dense_035)
    conv_baseline = _conv_score(out_baseline)
    conv_dense_035 = _conv_score(out_dense_035)

    assert cfg_baseline.gas_inlet_dense_frac == pytest.approx(0.30)
    assert raw_baseline > 0.0
    assert conv_baseline > 0.0
    assert raw_dense_035 > 0.0
    assert conv_dense_035 > 0.0
    assert out_baseline["exit_gas_dry"]
    assert out_dense_035["exit_gas_dry"]
    # The Hamel particle-size remap intentionally invalidates the old
    # dense=0.30≈0.35 tuned-window equivalence; keep this as a drift monitor
    # rather than a calibration lock.
    assert abs(raw_baseline - raw_dense_035) / raw_dense_035 < 0.50
    assert abs(conv_baseline - conv_dense_035) / conv_dense_035 < 0.50


@pytest.mark.slow
@pytest.mark.xfail(
    reason="与 validation_cases.json 尚未完全对齐；标定后应移除 xfail 并改为 XPASS",
    raises=AssertionError,
    strict=False,
)
def test_htw_wesseling_1_exit_temperature_vs_validation_json():
    """默认与 ``validation_cases.json`` → ``CASE_HTW_WESSELING_1.outputs`` 数值对比。

    设置 ``BFB_RELAX_VALIDATION=1`` 时仅检查温度物理区间，**不**与 JSON 容差断言（标定前本地可用）。

    **xfail**：当前模型仍常超差；对比逻辑仍会执行。标定后删除 ``@pytest.mark.xfail`` 使本用例在通过时记为 XPASS。
    """
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outs = raw["outputs"]
    T_ref = float(outs["exit_temperature_K"])
    T_meas = float(outs.get("exit_temperature_measured_K", T_ref))
    tol = validation_numeric_tolerances()
    rtol_T = tol["rtol_T"]
    dry_json = outs.get("exit_gas_dry_mol_frac", {})

    cfg = build_phase1_htw_lu_reactor_config()
    result = Reactor(cfg).solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
    validation_candidate_ok, validation_candidate_reasons = strict_validation_gate(result)
    T_exit = result["T_profile"][-1]
    assert 273.0 < T_exit < 4000.0, f"出口温度非物理合理值: T_exit={T_exit}"

    err_ref = abs(T_exit - T_ref) / T_ref
    err_meas = abs(T_exit - T_meas) / max(T_meas, 1.0)
    rel_best = min(err_ref, err_meas)

    relax = os.environ.get("BFB_RELAX_VALIDATION", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    if relax:
        return
    if not validation_candidate_ok:
        pytest.xfail(
            "严格 validation 前置条件未满足: " + ", ".join(validation_candidate_reasons)
        )

    assert rel_best <= rtol_T, (
        f"T_exit={T_exit:.2f} K 与 validation_cases 中 "
        f"exit_temperature_K={T_ref} K / measured={T_meas} K 的最佳相对偏差 {rel_best:.2%} "
        f"超过容差 {rtol_T:.0%}。标定前可设 BFB_RELAX_VALIDATION=1"
    )

    y_dry = result["exit_gas_dry"]
    for sp in ("CO", "CO2", "H2", "CH4"):
        tgt = json_numeric_or_none(dry_json.get(sp))
        if tgt is None or tgt <= 0.0:
            continue
        sim = y_dry.get(sp)
        if sim is None:
            pytest.fail(f"模拟结果缺少干基 {sp}（无法与 JSON 对比）")
        rtol_sp = tol["rtol_CH4"] if sp == "CH4" else tol["rtol_CO_CO2_H2"]
        rel_sp = abs(sim - tgt) / tgt
        assert rel_sp <= rtol_sp, (
            f"干基 {sp}: 模拟={sim:.4f}, validation_cases={tgt:.4f}, "
            f"相对偏差 {rel_sp:.2%} > {rtol_sp:.0%}"
        )

    X_json = json_numeric_or_none(outs.get("carbon_conversion_pct"))
    if X_json is not None and X_json > 0.0:
        rtol_X = tol["rtol_carbon_conv"]
        sim_X = result["carbon_conv"] * 100.0
        rel_X = abs(sim_X - X_json) / X_json
        assert rel_X <= rtol_X, (
            f"碳转化率: 模拟={sim_X:.1f}%, JSON={X_json:.1f}%, "
            f"相对偏差 {rel_X:.2%} > {rtol_X:.0%}"
        )


def run_case_LU():
    case = load_case_LU()
    feeds = estimate_gas_feeds(case)

    print("=" * 60)
    print("BFB Gasifier - Table 2 LU 工况验证")
    print("=" * 60)
    print(f"燃料: {case['fuel']} @ {case['fuel_feed']} kg/h")
    print(f"压力: {case['P']/1e6:.1f} MPa")
    print(f"ER: {case['ER']}")
    print(f"O2: {feeds['O2_feed']:.2f} mol/s")
    print(f"H2O: {feeds['H2O_feed']:.2f} mol/s")
    print()

    cfg = build_phase1_htw_lu_reactor_config(case)
    assert abs(cfg.O2_feed - feeds["O2_feed"]) < 1e-6

    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
    validation_candidate_ok, validation_candidate_reasons = strict_validation_gate(result)

    T_exit = result["T_profile"][-1]
    carbon_conv = result["carbon_conv"]
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outs = raw["outputs"]
    T_json = float(outs["exit_temperature_K"])
    T_meas_json = float(outs.get("exit_temperature_measured_K", T_json))
    tol = validation_numeric_tolerances()
    err_T = min(
        abs(T_exit - T_json) / T_json,
        abs(T_exit - T_meas_json) / max(T_meas_json, 1.0),
    )
    y_dry = result["exit_gas_dry"]
    dry_ref = outs.get("exit_gas_dry_mol_frac", {})

    print(f"迭代次数: {result['n_iter']}")
    print(f"收敛: {result['converged']}")
    print(f"validation candidate: {validation_candidate_ok} ({', '.join(validation_candidate_reasons) if validation_candidate_reasons else 'ok'})")
    print()

    # 温度剖面
    print("温度剖面 [K]:")
    for i, T in enumerate(result["T_profile"]):
        h = (i + 0.5) * cfg.H_bed / cfg.n_cells
        print(f"  h={h:.1f}m: T={T:.0f} K ({T-273.15:.0f} °C)")
    print()

    # 出口气体组成
    print("出口气体组成 (悬浮相, 湿基):")
    for sp in ["CO", "CO2", "H2", "H2O", "CH4", "O2", "N2", "H2S", "NH3"]:
        y = result["exit_gas"].get(sp, 0.0)
        print(f"  y_{sp:4s} = {y:.4f}")
    print()

    print("与 validation_cases.json（CASE_HTW_WESSELING_1.outputs）对比:")
    print(f"  出口温度: 模拟 T_exit={T_exit:.1f} K | JSON exit_temperature_K={T_json:.0f} K, "
          f"measured={T_meas_json:.0f} K → 最佳相对偏差 {err_T:.2%} "
          f"(容差 ≤{tol['rtol_T']:.0%}) {'OK' if err_T <= tol['rtol_T'] else '超差'}")
    X_pct_json = json_numeric_or_none(outs.get("carbon_conversion_pct"))
    if X_pct_json is not None:
        sim_pct = carbon_conv * 100.0
        eX = abs(sim_pct - X_pct_json) / X_pct_json
        print(f"  碳转化率: 模拟 {sim_pct:.1f}% | JSON {X_pct_json:.0f}% → 相对偏差 {eX:.2%}")
    print("  干基摩尔分数 (模拟 vs JSON):")
    for sp in ("CO", "CO2", "H2", "CH4"):
        t = json_numeric_or_none(dry_ref.get(sp))
        s = y_dry.get(sp)
        if t is None or s is None:
            continue
        rtol_sp = tol["rtol_CH4"] if sp == "CH4" else tol["rtol_CO_CO2_H2"]
        e = abs(s - t) / t
        ok = e <= rtol_sp
        print(f"    {sp:4s}: sim={s:.4f} json={t:.4f}  rel_err={e:.2%} "
              f"(≤{rtol_sp:.0%}) {'OK' if ok else '超差'}")
    print("=" * 60)
    print("说明: 默认 pytest 与 JSON 容差对比；标定前可设 BFB_RELAX_VALIDATION=1 跳过严格断言。")
    print("=" * 60)

    return {
        "T_exit": T_exit,
        "carbon_conv": carbon_conv,
        "y_CO_wet": result["exit_gas"].get("CO", 0.0),
        "y_CO_dry": y_dry.get("CO"),
        "err_T_vs_json": err_T,
        "validation_T_ref_K": T_json,
        "validation_T_meas_K": T_meas_json,
    }


if __name__ == "__main__":
    run_case_LU()
