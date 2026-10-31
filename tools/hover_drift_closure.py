"""悬停闭合的研究（历史工具，固定在 T26j/T26k 时期）：闭合 A 的悬停平衡、温度敏感性与漂移闭合估算。

**本工具不在当前模型上运行。** 它假定的两个前提都已撤销：

- 捕获点原在 ADR 0049 第 3 步的自由空间焦炭块解出之后（灰块之前）。ADR 0052 删除了整个
  预启动，当前初值构造的第一次预计算之前只有 ADR 0045 第 3 项第 1 步的零存量；把零存量当作
  该工具的 A 平衡对照点，会把 stock_a、ratio_a 与 others 算成零起点的量并写回
  `reports/research-hover-drift/`，覆盖已有的非零研究数据。
- 工具沿用 λ = α(1 + q_k) 求单位消耗速率的分母。ADR 0051 起悬停状态的粒级迁移改用 α_s·μ，
  该关系在一般情形下不再成立。

复现请 checkout T26k 合并前的提交（`cf5db31` 之前，ADR 0049 与 ADR 0044 仍在时）。
若要在当前模型上重做这项研究，须另开研究卡，按 ADR 0051 的 μ 与 ADR 0052 的状态重写全部
计算、标签与输出路径，不能只移动捕获钩子。

以下历史代码保留作版本对照。
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.residual as residual
from hamel.freeboard.motion import terminal_settling_velocity
from hamel.solver.vector import pack_state, unknown_slots

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-hover-drift"
_CLASS = 2  # 粒级 3（0,25–0,5 mm）
_TEMPERATURES = (1073.0, 1123.0, 1173.0, 1223.0, 1273.0)
_SLOPES = (1.0, 1.2, 1.4)


def _capture() -> dict:
    """在当前模型上明确拒绝运行。见模块说明。

    审查（PR #139 P2）指出：ADR 0052 删除预启动后，把捕获点挪到固相牛顿之前取到的是
    ADR 0045 的零起点，不是本工具假设的“自由空间焦炭块已解出”的状态。
    """
    raise RuntimeError(
        "hover_drift_closure 固定在 T26j/T26k 时期（ADR 0049 预启动、ADR 0044 迁移）；"
        "ADR 0052 删除了预启动、ADR 0051 改用 α_s·μ，当前模型上不能运行。"
        "复现请 checkout cf5db31 之前的提交；在当前模型上重做须另开研究卡并重写全部计算与输出。"
    )



def _hover_balance(state, frozen) -> list[dict]:
    """粒级 3 各悬停焦炭行：存量 M、总出流系数 λ = −∂r/∂M、流入 I = r + λM。"""
    slots = [
        slot for slot in unknown_slots(frozen)
        if slot.kind == "char" and slot.zone == "freeboard" and "hover" in slot.name
        and slot.name.endswith(f".k{_CLASS + 1}")
    ]
    x = np.asarray(pack_state(state, frozen), dtype=float)
    base = _values(x, frozen)
    rows = []
    for slot in slots:
        stock = float(x[slot.index])
        step = 1.0e-6 * max(stock, 1.0)
        trial = x.copy()
        trial[slot.index] += step
        rate = -float(_values(trial, frozen)[slot.index] - base[slot.index]) / step
        rows.append({"name": slot.name, "stock_kg": stock, "lambda_per_s": rate, "inflow_kg_s": float(base[slot.index]) + rate * stock})
    return rows


def _lambda_at(state, frozen, temperature: float, name: str) -> float:
    hot = replace(state, freeboard_temperature_k=np.full(np.shape(state.freeboard_temperature_k), temperature))
    slot = next(item for item in unknown_slots(frozen) if item.name == name)
    x = np.asarray(pack_state(hot, frozen), dtype=float)
    base = _values(x, frozen)[slot.index]
    trial = x.copy()
    trial[slot.index] += 1.0e-6
    return -float(_values(trial, frozen)[slot.index] - base) / 1.0e-6


def _values(x, frozen) -> np.ndarray:
    """原残差在各变量上的值。历史代码读 ``initial._values``，该函数随 ADR 0052 删除。"""
    return np.asarray(residual.evaluate_residual(np.asarray(x, dtype=float), frozen).values, dtype=float)


def _drift(inflow, alpha, widths, slope):
    """漂移闭合的三种离散。界面 j 的窗口 W_j = ln(ū_g,j / ū_g,j+1)，按几何格自下而上，g_j = 3W_j/n。

    缩粒时 d ln d/dt = −α/3（m ∝ d³，焦炭密度恒定），d ln v_t/dt = −n·α/3。
    精确链：入射点在窗口内按 ln v_t 均匀；一阶迎风：κ_j = α/g_j；
    上沿入射（ADR 0051 的出口因子，不含 ρ_k 映射）：κ_j = α/(e^{g_j} − 1)。
    最上一个稳定界面的出流进入上方格的自由状态。
    """
    exact = np.zeros(widths.size)
    upwind = np.zeros(widths.size)
    top = np.zeros(widths.size)
    carried = 0.0
    flux = 0.0
    passed = 0.0
    for j, width in enumerate(widths):
        g = 3.0 * width / slope
        kept = (1.0 - np.exp(-g)) / g
        exact[j] = inflow[j] / alpha * (1.0 - kept) + carried * (1.0 - np.exp(-g)) / alpha
        carried = inflow[j] * kept + carried * np.exp(-g)
        kappa = slope * alpha / (3.0 * width)
        upwind[j] = (inflow[j] + flux) / (alpha + kappa)
        flux = kappa * upwind[j]
        exit_rate = alpha / np.expm1(g)
        top[j] = (inflow[j] + passed) / (alpha + exit_rate)
        passed = exit_rate * top[j]
    return exact, carried, upwind, flux, top, passed


def main() -> None:
    """历史入口。在当前模型上 `_capture()` 直接拒绝，见模块说明。"""
    data = _capture()
    _OUT.mkdir(parents=True, exist_ok=True)
    state, frozen = data["state"], data["frozen"]
    check = residual._Context.domain
    residual._Context.domain = lambda self: None  # 只作诊断：求系数时暂停有效域检查
    try:
        rows = _hover_balance(state, frozen)
        bottom = rows[-1]["name"]
        lambdas = {temperature: _lambda_at(state, frozen, temperature, bottom) for temperature in _TEMPERATURES}
    finally:
        residual._Context.domain = check

    diameters = np.asarray(frozen.diameter_m, dtype=float)
    masses = frozen.rho_p_kg_m3 * np.pi * diameters ** 3 / 6.0
    q = masses[_CLASS - 1] / (masses[_CLASS] - masses[_CLASS - 1])
    alpha = rows[-1]["lambda_per_s"] / (1.0 + q)
    # 衡算下标 0 为顶格；几何格自下而上，界面 j 在几何格 j 与 j+1 之间、悬停记在几何格 j。
    u_g = data["u_g"]
    count = u_g.size - 1  # 最上一格之上是出口，没有稳定界面
    geometric = list(reversed(rows))[:count]
    inflow = np.array([row["inflow_kg_s"] for row in geometric])
    stock_a = np.array([row["stock_kg"] for row in geometric])
    widths = np.log(u_g[:count] / u_g[1:count + 1])
    capacity = np.asarray(frozen.v_freeboard_m3, dtype=float)[::-1] * (1.0 - frozen.epsilon_mf) * frozen.rho_p_kg_m3
    solids = (
        np.asarray(state.freeboard_char_free_kg).sum(axis=1) + np.asarray(state.freeboard_char_hover_kg).sum(axis=1)
    )[::-1]
    others = solids.copy()
    others[:count] -= stock_a
    ratio_a = solids / capacity

    rho, mu, rho_p = data["rho_g"], data["mu_g"], data["rho_p"]
    settling = []
    for diameter in (0.25e-3, 0.30e-3, 0.35e-3, 0.40e-3, 0.45e-3, 0.50e-3):
        values = [terminal_settling_velocity(rho[1], mu[1], diameter, rho_p, phi) for phi in (0.6, 0.8, 1.0)]
        nudged = terminal_settling_velocity(rho[1], mu[1], 1.01 * diameter, rho_p, 0.8)
        settling.append({"d_mm": diameter * 1.0e3, "v_t": values, "n": float(np.log(nudged / values[1]) / np.log(1.01))})

    drift = []
    for slope in _SLOPES:
        exact, out_exact, upwind, out_upwind, top, out_top = _drift(inflow, alpha, widths, slope)
        drift.append({
            "n": slope,
            "exact_kg": exact.tolist(), "exact_out_kg_s": out_exact,
            "upwind_kg": upwind.tolist(), "upwind_out_kg_s": out_upwind,
            "top_kg": top.tolist(), "top_out_kg_s": out_top,
            "ratio_top": ((others + np.r_[top, np.zeros(solids.size - count)]) / capacity).tolist(),
            "ratio_exact": ((others + np.r_[exact, np.zeros(solids.size - count)]) / capacity).tolist(),
            "ratio_upwind": ((others + np.r_[upwind, np.zeros(solids.size - count)]) / capacity).tolist(),
        })
    payload = {
        "case": "vtt_pressurized 工况 1，T26j 预启动焦炭块解出后（灰块之前）；历史工具，当前模型不适用",
        "hover_rows": rows,
        "q_k": float(q),
        "alpha_per_s": float(alpha),
        "u_g_m_s": u_g.tolist(),
        "window_widths": widths.tolist(),
        "capacity_kg": capacity.tolist(),
        "ratio_a": ratio_a.tolist(),
        "lambda_by_temperature": {str(int(key)): value for key, value in lambdas.items()},
        "settling": settling,
        "drift": drift,
        "freeboard_temperature_k": np.asarray(state.freeboard_temperature_k, dtype=float).tolist(),
        "bed_char_kg": float(np.asarray(state.bed_char_kg).sum()),
    }
    (_OUT / "study.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps({key: payload[key] for key in ("alpha_per_s", "q_k", "ratio_a", "lambda_by_temperature")}, indent=1))
    for item in drift:
        print(item["n"], np.round(item["exact_kg"], 2), round(sum(item["exact_kg"]), 2), np.round(item["ratio_exact"], 3))
        print(item["n"], np.round(item["upwind_kg"], 2), round(sum(item["upwind_kg"]), 2), np.round(item["ratio_upwind"], 3))
        print(item["n"], np.round(item["top_kg"], 2), round(sum(item["top_kg"]), 3), np.round(item["ratio_top"], 3))


if __name__ == "__main__":
    sys.exit(main())
