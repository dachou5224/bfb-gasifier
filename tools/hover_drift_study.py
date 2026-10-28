"""悬停颗粒随转化的漂移：研究脚本，只报告，不进入求解。

取工况 1 初值点第 1 轮预计算的自由空间气速剖面 u_gb(h) 与各格气体物性，
跟踪一个悬停粒级颗粒在三种转化方式下终端速度 v_t 与平衡高度 h* 随焦炭转化率 X 的变化：

- SP：缩粒（Hamel 模拟中采用，S. 85）。组成不变，直径按剩余质量缩小，密度不变。
- SD：缩密度（Hamel S. 77–78 对气化的物理图像）。直径不变，焦炭从内部消耗，密度下降。
- AR：保灰收缩（目标所有者提出的图像）。只按焦炭体积收缩，灰留下，密度升高。

焦炭与灰的密度 Hamel 未给 VTT 的数值，这里按区间扫描，只作探索。
同时给出 Reh 图口径的 Ar、Re、Y = 0.75·Re²/Ar（与 reh_app 相同的定义）。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard.motion import terminal_settling_velocity
from hamel.freeboard.transport import ghost_bubble_velocity
from hamel.solver.precalc import precalculate

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-hover-density"
G = 9.81


def conversion_state(mode: str, x: float, d0: float, rho0: float, f_c0: float, rho_ash: float) -> tuple[float, float]:
    """返回 (直径, 颗粒密度)。x 为颗粒焦炭的转化率，f_c0 为初始焦炭质量分数。"""
    remaining = 1.0 - f_c0 * x
    if mode == "SP":
        # 组成不变：焦炭消耗时同体积的灰一起释放（ADR 0031、0044），颗粒质量按 1 − x 缩小。
        return d0 * (1.0 - x) ** (1.0 / 3.0), rho0
    if mode == "SD":
        return d0, rho0 * remaining
    if mode == "AR":
        # 1/ρ0 = f_c0/ρ_C + (1 − f_c0)/ρ_A，由 ρ0 与 ρ_A 反求焦炭表观密度。
        inverse_char = (1.0 / rho0 - (1.0 - f_c0) / rho_ash) / f_c0
        if inverse_char <= 0.0:
            raise ValueError("给定的 ρ_A 与 f_C0 使焦炭密度无解")
        rho_char = 1.0 / inverse_char
        char_volume = f_c0 / rho_char * (1.0 - x)
        ash_volume = (1.0 - f_c0) / rho_ash
        volume_ratio = (char_volume + ash_volume) * rho0
        return d0 * volume_ratio ** (1.0 / 3.0), rho0 * remaining / volume_ratio
    raise ValueError(mode)


def _velocity_profile(column):
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)

    def velocity(height):
        cell = int(min(max(np.searchsorted(edges, height, side="right") - 1, 0), edges.size - 2))
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, column["bubble_m_s"], column["bubble_diameter_m"]), cell
    return edges, velocity


def stable_root(edges, velocity, rho_g, mu_g, diameter, rho_p, phi_s, samples=4000):
    """自下而上第一个 u_g − v_t 由正变负的点（稳定平衡）。没有则返回 None。"""
    heights = np.linspace(edges[0], edges[-1], samples)
    previous = None
    for height in heights:
        gas, cell = velocity(float(height))
        settle = terminal_settling_velocity(float(rho_g[cell]), float(mu_g[cell]), diameter, rho_p, phi_s)
        residual = gas - settle
        if previous is not None and previous[1] > 0.0 >= residual:
            h0, r0 = previous
            return h0 + (height - h0) * r0 / (r0 - residual)
        previous = (float(height), residual)
    return None


def reh_point(velocity, diameter, rho_p, rho_g, mu_g):
    nu = mu_g / rho_g
    re = velocity * diameter / nu
    ar = G * diameter**3 * (rho_p - rho_g) / (nu**2 * rho_g)
    return {"Re": re, "Ar": ar, "Y": 0.75 * re**2 / ar}


def study(column, klass=2, f_c0_values=(0.9, 0.5, 0.1), rho_ash_values=(1500.0, 2000.0, 2500.0), alpha=7.0e-4):
    edges, velocity = _velocity_profile(column)
    rho_g = np.asarray(column["rho_g"], dtype=float)
    mu_g = np.asarray(column["mu_g"], dtype=float)
    rho0 = float(column["rho_p"])
    phi_s = float(column["phi_s"])
    d0 = float(np.asarray(column["diameter_m"])[klass])
    h0 = stable_root(edges, velocity, rho_g, mu_g, d0, rho0, phi_s)
    home = None if h0 is None else int(min(np.searchsorted(edges, h0, side="right") - 1, edges.size - 2))
    profile = [{"h_m": float(h), "u_g_m_s": float(velocity(float(h))[0])} for h in np.linspace(edges[0], edges[2], 25)]
    cases = []
    grid = np.linspace(0.0, 0.99, 100)
    for f_c0 in f_c0_values:
        for mode in ("SP", "SD", "AR"):
            for rho_ash in (rho_ash_values if mode == "AR" else (None,)):
                track = []
                exit_x = None
                exit_kind = None
                for x in grid:
                    try:
                        diameter, rho_p = conversion_state(mode, float(x), d0, rho0, f_c0, rho_ash or 2000.0)
                    except ValueError as exc:
                        exit_kind = f"无解：{exc}"
                        break
                    root = stable_root(edges, velocity, rho_g, mu_g, diameter, rho_p, phi_s)
                    cell_gas = velocity(root if root is not None else edges[0])
                    settle = terminal_settling_velocity(float(rho_g[cell_gas[1]]), float(mu_g[cell_gas[1]]), diameter, rho_p, phi_s)
                    reh = reh_point(cell_gas[0], diameter, rho_p, float(rho_g[cell_gas[1]]), float(mu_g[cell_gas[1]]))
                    track.append({"X": float(x), "d_m": diameter, "rho_p": rho_p, "v_t_m_s": settle, "h_star_m": root, **reh})
                    if exit_x is None and home is not None:
                        if root is None:
                            surface_gas = velocity(float(edges[0]))[0]
                            exit_kind = "落回床面" if settle >= surface_gas else "无平衡点，被带出"
                            exit_x = float(x)
                        elif root >= edges[home + 1]:
                            exit_kind = "平衡点上移出原悬停格"
                            exit_x = float(x)
                        elif root < edges[home]:
                            exit_kind = "平衡点下移出原悬停格"
                            exit_x = float(x)
                residence = None if exit_x is None else -math.log(max(1.0 - exit_x, 1e-12)) / alpha
                cases.append({
                    "mode": mode, "f_c0": f_c0, "rho_ash": rho_ash, "exit_X": exit_x, "exit": exit_kind,
                    "residence_s_at_alpha": residence, "track": track,
                })
    return {"d0_m": d0, "rho0": rho0, "phi_s": phi_s, "h_star0_m": h0, "home_cell": home, "edges_m": [float(e) for e in edges],
            "alpha_per_s": alpha, "profile": profile, "cases": cases}


def main():
    prepared = precalculate(_CASE, "1", None, adiabatic=False, iteration=1, solver=SolverConfig())
    column = prepared.reported["trajectory"]
    result = study(column)
    _OUT.mkdir(parents=True, exist_ok=True)
    for case in result["cases"]:
        case["track"] = case["track"][::3]
    (_OUT / "drift.json").write_text(json.dumps(result, ensure_ascii=False, separators=(",", ":"), default=float), encoding="utf-8")
    sys.stdout.write(json.dumps({k: result[k] for k in ("d0_m", "rho0", "h_star0_m", "home_cell")}, default=float) + "\n")
    for case in result["cases"]:
        sys.stdout.write(f"{case['mode']} f_C0={case['f_c0']} ρ_A={case['rho_ash']}: {case['exit']} at X={case['exit_X']}, t≈{case['residence_s_at_alpha']}\n")


if __name__ == "__main__":
    main()
