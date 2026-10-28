"""逐格常系数的颗粒轨迹：闭式推导的数值核验，以及逐格常气速下的“夹住”界面。研究脚本，只报告。

式 (3.88) du/dt = A·u² + B·u + C，Δ = 4AC − B²（Hamel S. 46–47）。在系数为常数的区间内：

- Δ > 0：θ₀ = arctan((2A·u₀ + B)/√Δ)，
  u(τ) = [√Δ·tan(θ₀ + √Δ·τ/2) − B]/(2A)（与式 3.93 相同），
  h(τ) = h₀ − (1/A)·ln[cos(θ₀ + √Δ·τ/2)/cos θ₀] − B·τ/(2A)（与式 3.95 相同）。
- Δ < 0：s = √(−Δ)，Ψ = 2A·u + B，由 dτ = 2·dΨ/(Ψ² − s²) 得 E(τ) = [(Ψ₀ + s)/(Ψ₀ − s)]·exp(−s·τ)，
  Ψ(τ) = s·(E + 1)/(E − 1)，u = (Ψ − B)/(2A)，
  h(τ) = h₀ + [s·τ − B·τ − 2·ln|(E₀ − exp(s·τ))/(E₀ − 1)|]/(2A)，E₀ = E(0)。
- h(u)，两种情形：h = h₀ + (1/(2A))·ln(Q/Q₀) − (B/(2A))·I(u₀→u)，Q = A·u² + B·u + C，
  I 为 ∫du/Q：Δ > 0 时 (2/√Δ)·[arctan(Ψ/√Δ) − arctan(Ψ₀/√Δ)]，
  Δ < 0 时 (1/s)·ln[((Ψ − s)(Ψ₀ + s))/((Ψ + s)(Ψ₀ − s))]。

与印刷式比较：式 (3.93)–(3.95) 正确；式 (3.96) 的 ln 项前应为 +1/(2A)；式 (3.98) 在 τ = 0 不返回 u₀；
式 (3.103) 的 Ψ、Ψ₀ 定义对调且对数的分子分母配对不对（与 ADR 0016 第 11 条 T17 的核验一致）。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "research-cellwise-trajectory"


def closed_u(a, b, c, u0, tau):
    delta = 4.0 * a * c - b * b
    if delta > 0.0:
        root = math.sqrt(delta)
        theta = math.atan((2.0 * a * u0 + b) / root)
        return (root * math.tan(theta + root * tau / 2.0) - b) / (2.0 * a)
    s = math.sqrt(-delta)
    psi0 = 2.0 * a * u0 + b
    e = (psi0 + s) / (psi0 - s) * math.exp(-s * tau)
    return (s * (e + 1.0) / (e - 1.0) - b) / (2.0 * a)


def closed_h(a, b, c, u0, tau, h0=0.0):
    delta = 4.0 * a * c - b * b
    if delta > 0.0:
        root = math.sqrt(delta)
        theta = math.atan((2.0 * a * u0 + b) / root)
        return h0 - math.log(math.cos(theta + root * tau / 2.0) / math.cos(theta)) / a - b * tau / (2.0 * a)
    s = math.sqrt(-delta)
    psi0 = 2.0 * a * u0 + b
    e0 = (psi0 + s) / (psi0 - s)
    return h0 + (s * tau - b * tau - 2.0 * math.log(abs((e0 - math.exp(s * tau)) / (e0 - 1.0)))) / (2.0 * a)


def closed_h_of_u(a, b, c, u0, u, h0=0.0):
    delta = 4.0 * a * c - b * b
    q, q0 = a * u * u + b * u + c, a * u0 * u0 + b * u0 + c
    psi, psi0 = 2.0 * a * u + b, 2.0 * a * u0 + b
    if delta > 0.0:
        root = math.sqrt(delta)
        integral = 2.0 / root * (math.atan(psi / root) - math.atan(psi0 / root))
    else:
        s = math.sqrt(-delta)
        integral = math.log(abs(((psi - s) * (psi0 + s)) / ((psi + s) * (psi0 - s)))) / s
    return h0 + math.log(abs(q / q0)) / (2.0 * a) - b * integral / (2.0 * a)


def printed_u_398(a, b, c, u0, tau):
    """式 (3.98)（含 3.99）按印刷式，用于说明 τ = 0 不返回 u₀。"""
    s = math.sqrt(-(4.0 * a * c - b * b))
    eta = math.log(abs((s + (2.0 * a * u0 + b)) / (s - (2.0 * a * u0 + b)))) / s
    e = math.exp(eta * s - tau * s)
    return ((s - b) * e - s - b) / (2.0 * a * e - 2.0 * a)


def _rk4(a, b, c, u0, tau, steps=20000):
    dt = tau / steps
    u, h = u0, 0.0
    f = lambda v: a * v * v + b * v + c
    for _ in range(steps):
        k1 = f(u); k2 = f(u + dt * k1 / 2); k3 = f(u + dt * k2 / 2); k4 = f(u + dt * k3)
        l1 = u; l2 = u + dt * k1 / 2; l3 = u + dt * k2 / 2; l4 = u + dt * k3
        u += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        h += dt * (l1 + 2 * l2 + 2 * l3 + l4) / 6
    return u, h


def verify() -> list:
    """取几组物理量级的系数，比较闭式与 RK4。"""
    rows = []
    cases = [
        ("Δ>0，颗粒快于气体上行", -2.0, -30.0, -150.0, 6.0, 0.05),
        ("Δ<0，颗粒慢于气体（被带起）", -2.0, -40.0, 50.0, 0.5, 0.2),
        ("Δ<0，颗粒下落", 2.0, -40.0, -10.0, -0.5, 0.2),
    ]
    for name, a, b, c, u0, tau in cases:
        u_num, h_num = _rk4(a, b, c, u0, tau)
        u_cf, h_cf = closed_u(a, b, c, u0, tau), closed_h(a, b, c, u0, tau)
        row = {
            "case": name, "A": a, "B": b, "C": c, "u0": u0, "tau": tau, "delta": 4 * a * c - b * b,
            "u_rk4": u_num, "u_closed": u_cf, "h_rk4": h_num, "h_closed": h_cf,
            "h_of_u_closed": closed_h_of_u(a, b, c, u0, u_cf),
            "closed_u_at_0": closed_u(a, b, c, u0, 0.0),
        }
        if row["delta"] < 0.0:
            row["printed_398_at_0"] = printed_u_398(a, b, c, u0, 0.0)
        rows.append(row)
    return rows


def traps(column, mode="mean") -> list:
    """逐格常气速下，每个粒级的“夹住”界面：下格 u_g > v_t > 上格 u_g。"""
    from hamel.freeboard.motion import terminal_settling_velocity
    from hamel.freeboard.transport import ghost_bubble_velocity
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)
    rho, mu = np.asarray(column["rho_g"], dtype=float), np.asarray(column["mu_g"], dtype=float)

    def gas(cell, h):
        return ghost_bubble_velocity(cell, h, edges, diameters, flows, column["bubble_m_s"], column["bubble_diameter_m"])
    cells = edges.size - 1
    if mode == "mean":
        u_cell = [float(np.mean([gas(i, h) for h in np.linspace(edges[i], edges[i + 1], 201)[:-1] + 0.5 * (edges[i + 1] - edges[i]) / 200])) for i in range(cells)]
    else:
        u_cell = [gas(i, 0.5 * (edges[i] + edges[i + 1])) for i in range(cells)]
    rows = []
    for klass, d in enumerate(np.asarray(column["diameter_m"], dtype=float)):
        vt = [terminal_settling_velocity(float(rho[i]), float(mu[i]), float(d), float(column["rho_p"]), float(column["phi_s"])) for i in range(cells)]
        caught = [i for i in range(cells - 1) if u_cell[i] > vt[i] and u_cell[i + 1] < vt[i + 1]]
        rows.append({"k": klass + 1, "d_mm": float(d) * 1e3, "v_t": vt, "trap_interfaces": [float(edges[i + 1]) for i in caught]})
    return [{"u_cell": u_cell, "mode": mode}, rows]


def main():
    _OUT.mkdir(parents=True, exist_ok=True)
    out = {"verify": verify()}
    if len(sys.argv) > 1 and sys.argv[1] == "traps":
        from hamel.cases.config import SolverConfig
        from hamel.solver.precalc import precalculate
        prepared = precalculate(_ROOT / "data" / "cases" / "vtt_pressurized.json", "1", None, adiabatic=False, iteration=1, solver=SolverConfig())
        column = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in prepared.reported["trajectory"].items()}
        out["column"] = column
        out["traps_mean"] = traps(column, "mean")
        out["traps_mid"] = traps(column, "mid")
    (_OUT / "cellwise.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    for row in out["verify"]:
        sys.stdout.write(json.dumps(row, ensure_ascii=False) + "\n")
    for key in ("traps_mean", "traps_mid"):
        if key in out:
            info, rows = out[key]
            sys.stdout.write(f"{key} u_cell={[round(x, 4) for x in info['u_cell']]}\n")
            for r in rows:
                sys.stdout.write(f"  k={r['k']} d={r['d_mm']:.4f} mm v_t(底格)={r['v_t'][0]:.4f} 夹住界面={r['trap_interfaces']}\n")


if __name__ == "__main__":
    main()
