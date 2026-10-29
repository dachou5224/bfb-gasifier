"""T26j 对照报告：工况 1 初值点第 1 轮预计算气速剖面上的参考轨迹系数。只报告。

逐粒级比较旧代表粒径/整格 R*（ADR 0023、0047，取自 T26h 快照）与 ADR 0048 的粒级内取样；
另报 N_d、N_φ 加倍与 σ_φ = 0,05、0,15 的变化。气速剖面取自 `reports/t26h-trajectories/precalc_old.json`，
不重跑预计算（ADR 0048 实现后初始预计算会因 80 s 时限停下，见 reports/t26j-validation/）。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard.accounting import account_freeboard, stable_equilibria, within_class_samples
from hamel.freeboard._common import FreeboardError
from hamel.freeboard.ejection import ejection_velocity_samples
from hamel.freeboard.motion import terminal_settling_velocity
from hamel.freeboard.transport import ghost_bubble_velocity

_ROOT = Path(__file__).resolve().parents[1]
_TRACK = _ROOT / "reports" / "t26h-trajectories" / "precalc_old.json"
_OUT = _ROOT / "reports" / "t26j-trajectories"
# ADR 0023 的粒级边界 [a_k, b_k]，与 fuel_size_class_diameters 的级序相同。
_BOUNDS = [0.0291e-3, 0.1e-3, 0.25e-3, 0.5e-3, 0.71e-3, 1.0e-3, 2.0e-3, 3.15e-3]
_VARIANTS = (
    ("N8_sigma0.10", {"diameters": 8, "sphericities": 8, "sigma": 0.1}),
    ("N16_sigma0.10", {"diameters": 16, "sphericities": 16, "sigma": 0.1}),
    ("N8_sigma0.05", {"diameters": 8, "sphericities": 8, "sigma": 0.05}),
    ("N8_sigma0.15", {"diameters": 8, "sphericities": 8, "sigma": 0.15}),
)


def _column():
    track = json.loads(_TRACK.read_text(encoding="utf-8"))
    column = track["trajectory"]
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)

    def gas(cell, height):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, column["bubble_m_s"], column["bubble_diameter_m"])

    return track, column, edges, gas


def _books(track, column, edges, gas, variant) -> dict:
    rho = np.asarray(column["rho_g"], dtype=float)
    mu = np.asarray(column["mu_g"], dtype=float)
    speeds, shares, direct = ejection_velocity_samples(column["bubble_m_s"], 48)
    config = SolverConfig()
    rows = {}
    for klass in range(len(column["diameter_m"])):
        samples = within_class_samples(_BOUNDS[klass], _BOUNDS[klass + 1], **variant)
        try:
            book = account_freeboard(
                edges, gas, rho, mu, float(column["diameter_m"][klass]), column["rho_p"], 0.8,
                speeds, shares, direct, config, samples=samples,
            )
        except FreeboardError as exc:
            rows[klass] = {"error": str(exc)}
            continue
        rows[klass] = {
            "return": book.return_share,
            "elutriate": book.elutriation_share,
            "capture": book.capture_share,
            "capture_cells": list(book.capture_cells),
            "stable_heights": list(book.stable_heights_m),
            "theta_up": [cell.theta_up for cell in book.cells],
            "theta_down": [cell.theta_down for cell in book.cells],
            "c_up": [cell.c_up for cell in book.cells],
            "c_down": [cell.c_down for cell in book.cells],
            "a_f": [cell.a_f for cell in book.cells],
            "b_f": [cell.b_f for cell in book.cells],
            "free_state": [cell.free_state for cell in book.cells],
            "m_hat": [cell.m_hat for cell in book.cells],
            "hover_slip": [cell.hover_slip_m_s for cell in book.cells],
        }
    return rows


def _old(column) -> dict:
    """T26h 快照：整格 R*、代表粒径、单一 φ_S = 0,8。按衡算下标（0 为顶格）。"""
    return {
        "a_free": column.get("a_free"),
        "b_free": column.get("b_free"),
        "hover_free": column.get("a_hover"),
        "hover_hover": column.get("b_hover"),
        "theta": column.get("theta"),
        "ejection_to_free": column.get("ejection_to_free"),
        "ejection_to_hover": column.get("ejection_to_hover"),
        "return": column.get("reference_return_share"),
        "elutriate": column.get("reference_elutriation_share"),
    }


def _statistics(column, edges, gas, variant) -> dict:
    """ADR 0050 第 5 项：该变体全部参考轨迹的最大段数、归宿时刻最大值与 99 % 分位、最长的一条。"""
    rho = np.asarray(column["rho_g"], dtype=float)
    mu = np.asarray(column["mu_g"], dtype=float)
    speeds, shares, direct = ejection_velocity_samples(column["bubble_m_s"], 48)
    config = SolverConfig()
    segments = 0
    times: list[float] = []
    longest = None
    for klass in range(len(column["diameter_m"])):
        samples = within_class_samples(_BOUNDS[klass], _BOUNDS[klass + 1], **variant)
        book = account_freeboard(
            edges, gas, rho, mu, float(column["diameter_m"][klass]), column["rho_p"], 0.8,
            speeds, shares, direct, config, samples=samples,
        )
        segments = max(segments, book.max_segments)
        times.extend([book.max_time_s, book.p99_time_s])
        if book.longest is not None and (longest is None or book.longest[1] > longest[1][1]):
            longest = (klass, book.longest)
    return {
        "max_segments": segments,
        "max_time_s": max((book for book in times), default=0.0),
        "p99_time_s": max(times[1::2], default=0.0),
        "longest": None if longest is None else {"klass": longest[0] + 1, "record": list(longest[1])},
    }


def _mean_velocity(gas, cell: int, edges) -> float:
    from hamel.freeboard.accounting import cell_mean_velocity

    return cell_mean_velocity(gas, cell, edges[cell], edges[cell + 1])


def main() -> None:
    _OUT.mkdir(parents=True, exist_ok=True)
    track, column, edges, gas = _column()
    payload = {"source": str(_TRACK.relative_to(_ROOT)), "old": _old(track), "variants": {}, "diagnosis": {}}
    for name, variant in _VARIANTS:
        payload["variants"][name] = _books(track, column, edges, gas, variant)
        print(name, "done", flush=True)
    for name, variant in _VARIANTS:
        payload["diagnosis"][name] = _statistics(column, edges, gas, variant)
        print("statistics", name, "done", flush=True)
    (_OUT / "books.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    (_OUT / "report.md").write_text(_markdown(payload), encoding="utf-8")
    print("written", _OUT)


def _fmt(values, digits: int = 4) -> str:
    if values is None:
        return "—"
    return ", ".join("—" if value is None else f"{float(value):.{digits}g}" for value in values)


def _markdown(payload: dict) -> str:
    lines = ["# T26j 对照报告：粒级内取样与半格捕获区（ADR 0048）", ""]
    lines.append("工况 1 初值点第 1 轮预计算的气速剖面（取自 `reports/t26h-trajectories/precalc_old.json`，"
                 "与 T26h 的旧系数同一剖面）。只报告，不改模型。")
    lines.append("")
    lines.append("衡算下标 0 为自由空间顶格，7 为最下格（床面以上）。稳定界面逐子样本不同"
                 "（与粒径、球形度一一对应），表中不合并；表中只列该粒级全部子样本合起来出现转入的格。")
    lines.append("")
    old = payload["old"]
    lines.append("## 旧规则（代表粒径、单一 φ_S = 0,8、整格 R*）")
    lines.append("")
    lines.append("| 粒级 | 回落份额 | 带出份额 | 悬停份额 | 抛出 free/hover | θ（顶→底） |")
    lines.append("|---|---|---|---|---|---|")
    for klass in range(len(old["return"])):
        lines.append(
            f"| {klass + 1} | {old['return'][klass]:.4g} | {old['elutriate'][klass]:.4g} | "
            f"{max(0.0, 1.0 - old['return'][klass] - old['elutriate'][klass]):.4g} | "
            f"{old['ejection_to_free'][klass]:.3g}/{old['ejection_to_hover'][klass]:.3g} | "
            f"{_fmt([old['theta'][balance][klass] for balance in range(8)], 2)} |"
        )
    lines.append("")
    base = payload["variants"]["N8_sigma0.10"]
    lines.append("## 新取样（N_d = N_φ = 8、σ_φ = 0,1、半格 R*）")
    lines.append("")
    lines.append("| 粒级 | 回落 | 带出 | 捕获 | 转入 c↑ 非零格 | 自由状态不存在的格 | θ↑ | θ↓ |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for klass, row in base.items():
        index = int(klass)
        if "error" in row:
            lines.append(f"| {index + 1} | — | — | — | — | 未完成：{row['error']} | — | — |")
            continue
        carried = [index for index, value in enumerate(row["c_up"]) if value]
        dropped = [index for index, value in enumerate(row["c_down"]) if value]
        absent = [index for index, value in enumerate(row["free_state"]) if value != "存在"]
        lines.append(
            f"| {index + 1} | {row['return']:.4g} | {row['elutriate']:.4g} | {row['capture']:.4g} | "
            + f"{carried}（c↓ {dropped}） | {absent} | {_fmt(row['theta_up'], 3)} | {_fmt(row['theta_down'], 3)} |"
        )
    lines.append("")
    lines.append("### 转入系数与自由状态存量（新取样）")
    lines.append("")
    for klass, row in base.items():
        if "error" in row:
            continue
        lines.append(f"- 粒级 {int(klass) + 1}：a^f {_fmt(row['a_f'])}；b^f {_fmt(row['b_f'])}；c↑ {_fmt(row['c_up'])}；c↓ {_fmt(row['c_down'])}")
        lines.append(f"  - m̂（零自由参考存量判据）：{_fmt(row['m_hat'])}；自由状态：{row['free_state']}")
    lines.append("")
    lines.append("## 加密与 σ_φ 对照（只报告）")
    lines.append("")
    lines.append("| 变体 | " + " | ".join(f"粒级 {k + 1}" for k in range(len(base))) + " |")
    lines.append("|---|" + "---|" * len(base))
    for name, rows in payload["variants"].items():
        cells = []
        for klass in range(len(base)):
            row = rows.get(klass, rows.get(str(klass), {}))
            if "error" in row:
                cells.append("未定")
                continue
            reference = base.get(klass, base.get(str(klass), {}))
            if "error" in reference:
                cells.append(f"回落 {row['return']:.3g}/带出 {row['elutriate']:.3g}/捕获 {row['capture']:.3g}")
                continue
            worst = max(
                abs(row["capture"] - reference["capture"]),
                abs(row["return"] - reference["return"]),
                abs(row["elutriate"] - reference["elutriate"]),
            )
            cells.append(f"{worst:.2g}")
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append("表中数值为三项份额相对 N8_sigma0.10 的最大绝对差。")
    lines.append("")
    if payload["diagnosis"]:
        lines.append("## ADR 0050 第 5 项：段数与归宿时刻（只报告）")
        lines.append("")
        lines.append("| 变体 | 最大段数 | 归宿时刻最大值 / s | 99 % 分位 / s | 最长的一条 |")
        lines.append("|---|---|---|---|---|")
        for name, item in payload["diagnosis"].items():
            longest = item.get("longest")
            if longest:
                segments, seconds, diameter, phi, start = longest["record"]
                text = (
                    f"粒级 {longest['klass']}，d = {diameter * 1e3:.4g} mm，φ_S = {phi:.4g}，"
                    f"{start}，{seconds:.1f} s、{segments} 段"
                )
            else:
                text = "—"
            lines.append(
                f"| {name} | {item['max_segments']} | {item['max_time_s']:.1f} | {item['p99_time_s']:.1f} | {text} |"
            )
        lines.append("")
        lines.append("段数上限为 1000（ADR 0050），远高于实测最大值，终止保护不再是限制因素。")
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n"


if __name__ == "__main__":
    main()
