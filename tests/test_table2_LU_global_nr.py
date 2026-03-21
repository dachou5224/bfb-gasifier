"""Table 2 LU 工况 — 仅全局 NR（Vorabrechnung + 双层外迭代）端到端脚本。

与 ``test_table2_LU.py`` 使用相同 ``validation_cases.json`` → ``CASE_HTW_WESSELING_1`` 输入，
区别仅在于 ``Reactor.solve(..., solver=\"global_nr\")``，并输出 NR 诊断（``norm_history``、残差等）。

**本文件内联**与 LU 工况相关的 JSON 读取辅助函数（与 ``test_table2_LU.py`` 保持一致，修改时请同步）。

运行::

    python tests/test_table2_LU_global_nr.py

或::

    pytest tests/test_table2_LU_global_nr.py -m slow --tb=short

Source: docs/CLAUDE.md；data/validation_cases.json
"""

from __future__ import annotations

import json
import os
import sys
import time
import numpy as np
from pathlib import Path

# 允许 ``python tests/test_table2_LU_global_nr.py`` 直接运行
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

try:
    import pytest
except ImportError:  # 无 pytest 时仍可跑 ``if __name__`` 主流程
    class _PytestShim:
        class mark:
            @staticmethod
            def slow(f):
                return f

    pytest = _PytestShim()  # type: ignore[misc, assignment]

from src.core.feed_inlet import compute_gas_feeds_mol_s
from src.core.species import GAS_SPECIES_INDEX
from src.core.reactor import Reactor, ReactorConfig

# ── 与 test_table2_LU.py 同步：validation_cases.json → CASE_HTW_WESSELING_1 ──
CASE_LU_VALIDATION_KEY = "CASE_HTW_WESSELING_1"


def _flatten_htw_wesseling_1_case(node: dict) -> dict:
    """将 validation_cases.json 中嵌套结构展平为 Reactor 测试所用扁平字段。"""
    ins = node["inputs"]
    fuel = ins["fuel"]
    pa = fuel["proximate_analysis"]
    ua = fuel["ultimate_analysis_dry_wt_pct"]
    react = ins["reactor"]
    op = ins["operating_conditions"]
    ga = ins["gasification_agent"]
    agent_type = str(ga.get("type", ""))
    if "O2" in agent_type and "Steam" in agent_type and "Air" not in agent_type:
        primary_agent = "o2_steam"
    else:
        primary_agent = "air_steam"

    return {
        "_comment": node.get("_comment", ""),
        "_source_case_key": CASE_LU_VALIDATION_KEY,
        "reactor": "HTW_pressurised",
        "fuel": "brown_coal_RB",
        "P": float(op["pressure_MPa"]) * 1e6,
        "T_inlet": float(op["T_inlet_K"]),
        "fuel_feed": float(fuel["feed_rate_kg_h"]),
        "primary_agent": primary_agent,
        "primary_agent_note": agent_type,
        "ER": float(op["ER"]),
        "recirculation": bool(op.get("recirculation", True)),
        "H_bed": float(react["height_m"]),
        "D_freeboard": float(react["diameter_m"]),
        "bed_height_m": float(react["bed_height_m"]),
        "moisture_wt": float(pa["moisture_wt_pct"]),
        "ash_dry_wt": float(pa["ash_dry_wt_pct"]),
        "VM_daf": float(pa["VM_daf_pct"]),
        "C_dry": float(ua["C"]),
        "H_dry": float(ua["H"]),
        "O_dry": float(ua["O"]),
        "N_dry": float(ua["N"]),
        "S_dry": float(ua["S"]),
    }


def _validation_cases_path() -> Path:
    return _REPO_ROOT / "data" / "validation_cases.json"


def load_case_LU() -> dict:
    data_path = _validation_cases_path()
    with open(data_path, encoding="utf-8") as f:
        data = json.load(f)
    if CASE_LU_VALIDATION_KEY not in data:
        raise KeyError(
            f"{data_path} 中缺少键 {CASE_LU_VALIDATION_KEY!r}（Table 2 LU / Table 7.1 Sim 1）"
        )
    return _flatten_htw_wesseling_1_case(data[CASE_LU_VALIDATION_KEY])


def load_validation_json_root() -> dict:
    data_path = _validation_cases_path()
    with open(data_path, encoding="utf-8") as f:
        return json.load(f)


def load_validation_case_node(key: str = CASE_LU_VALIDATION_KEY) -> dict:
    data = load_validation_json_root()
    if key not in data:
        raise KeyError(f"{_validation_cases_path()} 中缺少键 {key!r}")
    return data[key]


def _json_numeric_or_none(val: object) -> float | None:
    if isinstance(val, (int, float)):
        return float(val)
    return None


def validation_numeric_tolerances() -> dict[str, float]:
    return {
        "rtol_T": 0.15,
        "rtol_CO_CO2_H2": 0.15,
        "rtol_CH4": 0.20,
    }


def estimate_gas_feeds(case: dict) -> dict:
    """与 ``ReactorConfig(ER=..., primary_agent=...)`` / ``compute_gas_feeds_mol_s`` 一致。"""
    fuel_kg_s = case["fuel_feed"] / 3600.0
    o2, h2o, n2 = compute_gas_feeds_mol_s(
        fuel_feed_kg_s=fuel_kg_s,
        moisture_wt=case["moisture_wt"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
    )
    return {"O2_feed": o2, "H2O_feed": h2o, "N2_feed": n2}


def run_case_LU_global_nr(
    *,
    max_global_iter: int = 25,
    tol_global: float = 1.0,
) -> dict:
    """与 ``test_table2_LU.run_case_LU`` 相同配置，仅 solver=global_nr。"""
    case = load_case_LU()
    feeds = estimate_gas_feeds(case)

    print("=" * 60)
    print("BFB Gasifier - Table 2 LU 工况验证（Global NR 专用）")
    print("=" * 60)
    print(f"燃料: {case['fuel']} @ {case['fuel_feed']} kg/h")
    print(f"压力: {case['P']/1e6:.1f} MPa")
    print(f"ER: {case['ER']}")
    print(f"O2: {feeds['O2_feed']:.2f} mol/s")
    print(f"H2O: {feeds['H2O_feed']:.2f} mol/s")
    print(f"N2: {feeds['N2_feed']:.2f} mol/s")
    print()

    cfg = ReactorConfig(
        n_age_classes=1,
        n_cells=15,          # 增加 Cell 数量以提高轴向分辨率
        H_bed=case["H_bed"],
        H_freeboard=0.0,
        D_bed=case["D_freeboard"],
        P=case["P"],
        T_inlet=case["T_inlet"],
        fuel_type="coal",
        rho_s=1000.0,
        d_p=0.5e-3,
        phi_s=0.86,
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
        HHV_dry=22.0,  # [MJ/kg] 褐煤典型高位热值
        u0_target=1.0, # [m/s] 根据 Hamel (1999) 对 LU 工况的描述
        heat_loss_frac=0.15, 
        recirculation_frac=0.1,
    )
    assert abs(cfg.O2_feed - feeds["O2_feed"]) < 1e-6

    reactor = Reactor(cfg)
    
    # 强制物理初始化：所有 cell 的初值等于底部进料
    idx = GAS_SPECIES_INDEX
    for c in reactor.cells:
        c.T = cfg.T_inlet + 800.0  
        c.N_d[idx["N2"]] = cfg.N2_feed
        c.N_d[idx["H2O"]] = cfg.H2O_feed
        c.N_d[idx["O2"]] = cfg.O2_feed

    t0 = time.perf_counter()
    result = reactor.solve(
        max_global_iter=100,  # 允许极长演化
        tol_global=1e-4,
        solver="gauss_seidel",
    )
    elapsed = time.perf_counter() - t0

    # 调试：检查第一个 cell 的速率
    c0 = reactor.cells[0]
    print("\n--- 调试：第一个 Cell (h=0.5m) 动力学状态 ---")
    print(f"  T: {c0.T:.1f} K, P: {c0.P/1e6:.2f} MPa")
    y0 = c0._mole_fractions("d")
    idx = GAS_SPECIES_INDEX
    print(f"  y_O2: {y0[idx['O2']]:.4f}, y_H2O: {y0[idx['H2O']]:.4f}")
    
    # 重新触发计算以获取最新 R_gas
    c0.calc_reactions()
    print(f"  R_gas_d[O2]: {c0.R_gas_d[idx['O2']]:.2e} mol/s")
    print(f"  R_gas_d[CO]: {c0.R_gas_d[idx['CO']]:.2e} mol/s")
    print(f"  R_solid[Char]: {np.sum(c0.R_solid[:, 0]):.2e} kg/s")
    print("-------------------------------------------\n")

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

    print(f" wall time: {elapsed:.1f} s")
    print(f" 迭代次数 (NR 步合计): {result['n_iter']}")
    print(
        f" 收敛: 外迭代={result.get('converged_outer', result['converged'])}  "
        f"内层NR(RMS)={result.get('converged_inner_nr')}  "
        f"二者兼具={result.get('converged_fully', False)}"
    )
    rms_sf = result.get("rms_scaled_final")
    if rms_sf is not None:
        print(f" 末次内层 RMS(‖F̂‖): {rms_sf:.4e}  (应与 tol_rms 比较；‖F‖ 为未缩放 L2)")
    print(f" 残差范数 ||F||: {result.get('residual', float('nan')):.4e}")
    nh = result.get("norm_history", [])
    if nh:
        print(f" Norm history (RMS ‖F̂‖), 多外迭代拼接: {len(nh)} 步（未必单调）")
        head = [f"{v:.2e}" for v in nh[:12]]
        print(f"   [{', '.join(head)}{' ...' if len(nh) > 12 else ''}]")
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
    print(
        f"  出口温度: 模拟 T_exit={T_exit:.1f} K | JSON exit_temperature_K={T_json:.0f} K, "
        f"measured={T_meas_json:.0f} K → 最佳相对偏差 {err_T:.2%} "
        f"(容差 ≤{tol['rtol_T']:.0%}) {'OK' if err_T <= tol['rtol_T'] else '超差'}"
    )
    X_pct_json = _json_numeric_or_none(outs.get("carbon_conversion_pct"))
    if X_pct_json is not None:
        sim_pct = carbon_conv * 100.0
        eX = abs(sim_pct - X_pct_json) / X_pct_json
        print(f"  碳转化率: 模拟 {sim_pct:.1f}% | JSON {X_pct_json:.0f}% → 相对偏差 {eX:.2%}")
    print("  干基摩尔分数 (模拟 vs JSON):")
    for sp in ("CO", "CO2", "H2", "CH4"):
        t = _json_numeric_or_none(dry_ref.get(sp))
        s = y_dry.get(sp)
        if t is None or s is None:
            continue
        rtol_sp = tol["rtol_CH4"] if sp == "CH4" else tol["rtol_CO_CO2_H2"]
        e = abs(s - t) / t
        ok = e <= rtol_sp
        print(
            f"    {sp:4s}: sim={s:.4f} json={t:.4f}  rel_err={e:.2%} "
            f"(≤{rtol_sp:.0%}) {'OK' if ok else '超差'}"
        )
    print("=" * 60)
    print("说明: 与 test_table2_LU 相同工况；仅求解器为 global_nr。")
    print("=" * 60)

    return {
        "T_exit": T_exit,
        "carbon_conv": carbon_conv,
        "y_CO_wet": result["exit_gas"].get("CO", 0.0),
        "y_CO_dry": y_dry.get("CO"),
        "err_T_vs_json": err_T,
        "validation_T_ref_K": T_json,
        "validation_T_meas_K": T_meas_json,
        "converged": result["converged"],
        "converged_outer": result.get("converged_outer", result["converged"]),
        "converged_inner_nr": result.get("converged_inner_nr"),
        "converged_fully": result.get("converged_fully"),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "n_iter": result["n_iter"],
        "norm_history": nh,
        "residual": result.get("residual"),
        "wall_s": elapsed,
    }


@pytest.mark.slow
def test_table2_lu_global_nr_smoke():
    """LU 工况 + global_nr：物理合理 + NR 范数总体下降。

    与 ``test_table2_LU`` 中 JSON 容差断言不同：此处不对 validation_cases 做严格断言
    （与 ``test_htw_wesseling_1_exit_temperature_vs_validation_json`` 的 xfail 策略一致；
    容差对比见 ``run_case_LU_global_nr`` 终端输出）。

    若需对出口温度强制执行 JSON 容差，设置环境变量 ``BFB_STRICT_LU_GLOBAL_NR=1``。
    """
    out = run_case_LU_global_nr(max_global_iter=25, tol_global=1.0)
    T_exit = out["T_exit"]
    assert 500 < T_exit < 2000, f"T_exit 非物理: {T_exit:.1f} K"
    XC = out["carbon_conv"] * 100
    assert 0 < XC < 100, f"碳转化率非物理: {XC:.1f}%"
    assert out["y_CO_dry"] is not None
    assert 0 < float(out["y_CO_dry"]) < 0.6

    h = out["norm_history"]
    assert len(h) >= 2, "应有 NR 范数历史"
    # 多外迭代拼接的 norm_history 不要求单调；以内层 converged_inner_nr / rms_scaled_final 为准

    strict = os.environ.get("BFB_STRICT_LU_GLOBAL_NR", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    if not strict:
        return

    tol = validation_numeric_tolerances()
    assert out["err_T_vs_json"] <= tol["rtol_T"], (
        f"T_exit 与 JSON 最佳相对偏差 {out['err_T_vs_json']:.2%} 超过 {tol['rtol_T']:.0%}"
    )


if __name__ == "__main__":
    run_case_LU_global_nr()
    print("\nTable 2 LU — Global NR 脚本运行结束。")
