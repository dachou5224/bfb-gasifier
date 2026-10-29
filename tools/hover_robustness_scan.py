"""悬停的稳健性扫描：研究脚本，只报告，不进入求解。

在工况 1 初值点第 1 轮保存的自由空间剖面（reports/research-cellwise-trajectory/cellwise.json 的 column）上，
按 ADR 0047 的逐格常气速判定界面夹住（下格 ū_g > v_t > 上格 ū_g），
扫描 ρ_P、球形度、格内气速取法、d_b,WS（鬼泡衰减长度）与二次风入口所在界面，
给出各粒级被夹住的质量份额（粒级内按 ADR 0023 质量在 ln d 上均匀；球形度为研究假设：[0,3; 1,0] 截断正态、
截断后均值 0,8，即 ADR 0048 提议中的取法，尚未批准）与夹住界面；另给单一代表粒径、单一 φ_S 时的全部夹住界面。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import truncnorm

from hamel.freeboard.motion import terminal_settling_velocity
from hamel.freeboard.transport import ghost_bubble_velocity

_ROOT = Path(__file__).resolve().parents[1]
_COLUMN = _ROOT / "reports" / "research-cellwise-trajectory" / "cellwise.json"
_OUT = _ROOT / "reports" / "research-hover-robustness"
BOUNDS_MM = [(0.0291, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 0.71), (0.71, 1.0), (1.0, 2.0), (2.0, 3.15)]


def phi_nodes(sigma, n=8):
    if sigma is None:
        return np.array([0.8])
    loc = brentq(lambda m: truncnorm.mean((0.3 - m) / sigma, (1.0 - m) / sigma, loc=m, scale=sigma) - 0.8, 0.5, 1.2)
    a, b = (0.3 - loc) / sigma, (1.0 - loc) / sigma
    return truncnorm.ppf((np.arange(n) + 0.5) / n, a, b, loc=loc, scale=sigma)


def cell_velocity(column, rule, decay_scale=1.0, sa_shift=0):
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float).copy()
    if sa_shift:
        # 二次风入口移到相邻界面：把第 1 格与第 2 格之间的流量台阶整体移动 sa_shift 格（近似，只改体积流）。
        step = flows[1] - flows[0]
        base = flows[0]
        flows = np.array([base if i < 1 + sa_shift else base + step for i in range(flows.size)])
    d_b = float(column["bubble_diameter_m"]) * decay_scale
    gas = lambda i, h: ghost_bubble_velocity(i, h, edges, diameters, flows, column["bubble_m_s"], d_b)
    out = []
    for i in range(edges.size - 1):
        lo, hi = edges[i], edges[i + 1]
        if rule == "mean":
            hs = lo + (np.arange(200) + 0.5) / 200 * (hi - lo)
            out.append(float(np.mean([gas(i, h) for h in hs])))
        elif rule == "mid":
            out.append(gas(i, 0.5 * (lo + hi)))
        elif rule == "lower":
            out.append(gas(i, lo))
        else:
            out.append(gas(i, hi - 1e-9))
    return edges, np.array(out)


def trapped_share(column, u_cell, edges, rho_p, phis, klass, n_d=24):
    a, b = BOUNDS_MM[klass]
    ds = np.exp(np.log(a) + (np.arange(n_d) + 0.5) / n_d * (np.log(b) - np.log(a))) * 1e-3
    rho, mu = column["rho_g"], column["mu_g"]
    hit, where = 0, {}
    for d in ds:
        for phi in phis:
            for i in range(len(u_cell) - 1):
                below = terminal_settling_velocity(rho[i], mu[i], d, rho_p, phi)
                above = terminal_settling_velocity(rho[i + 1], mu[i + 1], d, rho_p, phi)
                if u_cell[i] > below and u_cell[i + 1] < above:
                    hit += 1
                    key = round(float(edges[i + 1]), 4)
                    where[key] = where.get(key, 0) + 1
                    break
    total = len(ds) * len(phis)
    return hit / total, {k: v / total for k, v in sorted(where.items())}


def representative_traps(column, u_cell, edges, rho_p, phi, klass):
    """单一代表粒径（ADR 0023）、单一 φ_S 时被夹住的全部界面。"""
    d = float(np.asarray(column["diameter_m"])[klass])
    rho, mu = column["rho_g"], column["mu_g"]
    return [round(float(edges[i + 1]), 4) for i in range(len(u_cell) - 1)
            if u_cell[i] > terminal_settling_velocity(rho[i], mu[i], d, rho_p, phi)
            and u_cell[i + 1] < terminal_settling_velocity(rho[i + 1], mu[i + 1], d, rho_p, phi)]


def scan():
    column = json.loads(_COLUMN.read_text(encoding="utf-8"))["column"]
    rho0 = float(column["rho_p"])
    cases = [("基准（高度平均，φ 截断正态 σ 0,1）", {})]
    for f in (0.75, 0.9, 0.98, 1.02, 1.1, 1.13):
        cases.append((f"ρ_P × {f}", {"rho": f}))
    cases.append(("φ_S 单一 0,78（代表值扰动）", {"sigma": None, "phi": 0.78}))
    for s in (None, 0.05, 0.15):
        cases.append(("φ_S 单一 0,8" if s is None else f"φ σ = {s}", {"sigma": s}))
    for r in ("mid", "lower", "upper"):
        cases.append((f"格内取 {r}", {"rule": r}))
    for k in (0.5, 2.0):
        cases.append((f"d_b,WS × {k}（衰减长度）", {"decay": k}))
    for sh in (-1, 1):
        cases.append((f"二次风入口移 {sh:+d} 格", {"sa": sh}))
    rows = []
    for name, opt in cases:
        edges, u_cell = cell_velocity(column, opt.get("rule", "mean"), opt.get("decay", 1.0), opt.get("sa", 0))
        phis = phi_nodes(opt.get("sigma", 0.1)) if "phi" not in opt else np.array([opt["phi"]])
        rho_p = rho0 * opt.get("rho", 1.0)
        per = []
        for k in range(len(BOUNDS_MM)):
            share, where = trapped_share(column, u_cell, edges, rho_p, phis, k)
            per.append({"k": k + 1, "share": share, "interfaces": where})
        rep_phi = opt.get("phi", 0.8)
        rep = representative_traps(column, u_cell, edges, rho_p, rep_phi, 2)
        rows.append({"case": name, "u_cell": [float(x) for x in u_cell], "classes": per, "class3_representative_traps": rep})
    return rows


def write(rows):
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "scan.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = ["# 悬停稳健性扫描（研究，只报告）", "",
             "剖面：工况 1 初值点第 1 轮（`reports/research-cellwise-trajectory/cellwise.json`）。判据：ADR 0047 的逐格界面夹住。",
             "份额：粒级内质量在 ln d 上均匀（ADR 0023，每级 24 点）× 球形度（**研究假设**：[0,3; 1,0] 截断正态、截断后均值 0,8、σ 0,1，8 点；即 ADR 0048 提议中的取法，尚未批准），满足夹住条件的潜在质量份额。",
             "只计夹住条件，不积分起始速度与真实捕获轨迹，不是模型悬停存量或捕获通量的验证；各扰动是冻结剖面上的单因素研究。",
             "二次风入口移动只把体积流台阶移到相邻界面，是近似。初值点剖面受第 1 轮冻结量扭曲，数值只说明敏感性。", "",
             "| 情形 | 粒级 1 | 2 | 3 | 4 | 5 | 6 | 7 | 夹住界面（m） | 粒级 3 单一代表值的夹住界面（m） |", "|---|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        shares = " | ".join(f"{c['share']:.2f}" for c in row["classes"])
        where = sorted({k for c in row["classes"] for k in c["interfaces"]})
        rep = ", ".join(str(w) for w in row["class3_representative_traps"]) or "无"
        lines.append(f"| {row['case']} | {shares} | {', '.join(str(w) for w in where) or '—'} | {rep} |")
    (_OUT / "scan.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    rows = scan()
    write(rows)
    sys.stdout.write((_OUT / "scan.md").read_text(encoding="utf-8"))
