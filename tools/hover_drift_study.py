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

from hamel.freeboard.motion import terminal_settling_velocity
from hamel.freeboard.transport import ghost_bubble_velocity

_ROOT = Path(__file__).resolve().parents[1]
_RECORD = _ROOT / "reports" / "t26c-diagnosis" / "baseline" / "record.json"
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


def _profile(column):
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)

    def gas(cell, height):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, column["bubble_m_s"], column["bubble_diameter_m"])

    def at(height):
        cell = int(min(max(np.searchsorted(edges, height, side="right") - 1, 0), edges.size - 2))
        return gas(cell, height), cell
    return edges, gas, at


def roots(edges, gas, rho_g, mu_g, diameter, rho_p, phi_s) -> list:
    """连续剖面 u_gb(h) 上的稳定根（du_g/dh < 0）。研究脚本自己的根搜索。

    ADR 0047 已把模型的悬停判定改为格界面（那里的 stable_equilibria 随之改写）；
    本脚本研究的是连续剖面，故在这里保留自己的二分，不调用模型。
    """
    count = len(edges) - 1
    settling = [terminal_settling_velocity(float(rho_g[i]), float(mu_g[i]), diameter, rho_p, phi_s) for i in range(count)]
    found = []
    for face in range(1, count):
        height = float(edges[face])
        if (float(gas(face - 1, height)) - settling[face - 1]) * (float(gas(face, height)) - settling[face]) < 0.0:
            found.append(height)
    for cell in range(count):
        left, right = float(edges[cell]), float(edges[cell + 1])
        previous_h = left
        previous_r = float(gas(cell, left)) - settling[cell]
        for step in range(1, 33):
            height = left + (right - left) * step / 32.0
            residual = float(gas(cell, height)) - settling[cell]
            if previous_r * residual < 0.0:
                lo, hi, r_lo = previous_h, height, previous_r
                for _ in range(60):
                    middle = 0.5 * (lo + hi)
                    r_mid = float(gas(cell, middle)) - settling[cell]
                    if r_lo * r_mid <= 0.0:
                        hi = middle
                    else:
                        lo, r_lo = middle, r_mid
                root = 0.5 * (lo + hi)
                if left < root < right:
                    found.append(root)
            previous_h, previous_r = height, residual
    return sorted(found)


def _cell_of(edges, height):
    return int(min(max(np.searchsorted(edges, height, side="right") - 1, 0), edges.size - 2))


def reh_point(velocity, diameter, rho_p, rho_g, mu_g):
    nu = mu_g / rho_g
    re = velocity * diameter / nu
    ar = G * diameter**3 * (rho_p - rho_g) / (nu**2 * rho_g)
    return {"Re": re, "Ar": ar, "Y": 0.75 * re**2 / ar}


def _fate_after_loss(edges, at, rho_g, mu_g, diameter, rho_p, phi_s, previous, current_roots):
    """原分支消失后颗粒的去向：在原高度上比较 u_g 与 v_t，向下落到下一个根或床面，向上升到下一个根或出口。"""
    gas, cell = at(previous)
    settle = terminal_settling_velocity(float(rho_g[cell]), float(mu_g[cell]), diameter, rho_p, phi_s)
    if settle > gas:
        lower = [h for h in current_roots if h < previous]
        return ("下落到下方平衡点", max(lower)) if lower else ("落回床面", None)
    upper = [h for h in current_roots if h > previous]
    return ("上升到上方平衡点", min(upper)) if upper else ("被带出自由空间", None)


def track_branch(edges, at, gas, rho_g, mu_g, d0, rho0, phi_s, start, mode, f_c0, rho_ash, grid, jump=0.005, min_dx=1e-6):
    """沿同一条平衡分支跟踪。

    相邻两点的最近根相差超过 jump 时，把 X 步长减半重算；步长已小于 min_dx 仍跳动（或无根）才判为分支消失。
    这样平衡点连续快速移动时不会被误判为换分支。
    """
    home = _cell_of(edges, start)
    state = {"previous": start, "exit_x": None, "exit_kind": None, "loss": None}
    track = []

    def evaluate(x):
        diameter, rho_p = conversion_state(mode, float(x), d0, rho0, f_c0, rho_ash or 2000.0)
        return diameter, rho_p, roots(edges, gas, rho_g, mu_g, diameter, rho_p, phi_s)

    def advance(x_from, x_to):
        diameter, rho_p, present = evaluate(x_to)
        previous = state["previous"]
        nearest = min(present, key=lambda h: abs(h - previous)) if present else None
        if nearest is None or abs(nearest - previous) > jump:
            if x_to - x_from > min_dx:
                middle = 0.5 * (x_from + x_to)
                return advance(x_from, middle) and advance(middle, x_to)
            state["loss"] = {"X": float(x_to), "from_m": previous,
                             "to": _fate_after_loss(edges, at, rho_g, mu_g, diameter, rho_p, phi_s, previous, present)}
            if state["exit_x"] is None:
                state["exit_x"], state["exit_kind"] = float(x_to), "分支消失：" + state["loss"]["to"][0]
            return False
        if state["exit_x"] is None and _cell_of(edges, nearest) != home:
            state["exit_x"] = float(x_to)
            state["exit_kind"] = "平衡点上移出原悬停格" if nearest > start else "平衡点下移出原悬停格"
        state["previous"] = nearest
        return True

    x_prev = float(grid[0])
    for x in grid:
        try:
            if x > x_prev and not advance(x_prev, float(x)):
                break
            diameter, rho_p, present = evaluate(x)
        except ValueError as exc:
            state["exit_kind"] = f"无解：{exc}"
            break
        here = state["previous"]
        u_here, cell_here = at(here)
        settle = terminal_settling_velocity(float(rho_g[cell_here]), float(mu_g[cell_here]), diameter, rho_p, phi_s)
        reh = reh_point(u_here, diameter, rho_p, float(rho_g[cell_here]), float(mu_g[cell_here]))
        track.append({"X": float(x), "d_m": diameter, "rho_p": rho_p, "v_t_m_s": settle, "h_star_m": here, "roots_m": present, **reh})
        x_prev = float(x)
    return {"start_m": start, "home_cell": home, "exit_X": state["exit_x"], "exit": state["exit_kind"],
            "branch_lost": state["loss"], "track": track}


def study(column, klass=2, f_c0_values=(0.9, 0.5, 0.1), rho_ash_values=(1500.0, 2000.0, 2500.0), alpha=7.0e-4):
    edges, gas, at = _profile(column)
    rho_g = np.asarray(column["rho_g"], dtype=float)
    mu_g = np.asarray(column["mu_g"], dtype=float)
    rho0 = float(column["rho_p"])
    phi_s = float(column["phi_s"])
    d0 = float(np.asarray(column["diameter_m"])[klass])
    starts = roots(edges, gas, rho_g, mu_g, d0, rho0, phi_s)
    grid = np.linspace(0.0, 0.99, 100)
    cases = []
    for start in starts:
        for f_c0 in f_c0_values:
            for mode in ("SP", "SD", "AR"):
                for rho_ash in (rho_ash_values if mode == "AR" else (None,)):
                    branch = track_branch(edges, at, gas, rho_g, mu_g, d0, rho0, phi_s, start, mode, f_c0, rho_ash, grid)
                    x = branch["exit_X"]
                    branch.update(mode=mode, f_c0=f_c0, rho_ash=rho_ash,
                                  residence_s_at_alpha=None if x is None else -math.log(max(1.0 - x, 1e-12)) / alpha)
                    cases.append(branch)
    return {"d0_m": d0, "rho0": rho0, "phi_s": phi_s, "initial_roots_m": starts, "edges_m": [float(e) for e in edges],
            "alpha_per_s": alpha, "cases": cases}


def main():
    with _RECORD.open(encoding="utf-8") as handle:
        record = json.load(handle)
    column = record["diagnostics"]["rounds"][0]["info"]["sensitivity"]["trajectory"]
    result = study(column)
    for case in result["cases"]:
        case["track"] = case["track"][::3]
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "drift.json").write_text(json.dumps(result, ensure_ascii=False, separators=(",", ":"), default=float), encoding="utf-8")
    sys.stdout.write(json.dumps({k: result[k] for k in ("d0_m", "rho0", "initial_roots_m")}, default=float) + "\n")
    for case in result["cases"]:
        lost = case["branch_lost"]
        sys.stdout.write(f"起点 {case['start_m']:.6f} {case['mode']} f_C0={case['f_c0']} ρ_A={case['rho_ash']}: {case['exit']} at X={case['exit_X']}, t≈{case['residence_s_at_alpha']}, 分支消失={None if lost is None else (lost['X'], lost['to'])}\n")


if __name__ == "__main__":
    main()
