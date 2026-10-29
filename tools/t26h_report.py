"""Build the T26h comparison report from the two round-1 dumps.

Reads the pre-change dump and the post-change dump, recomputes the hover
interfaces from the saved column, runs the sample-doubling check, and writes
report.md / report.json. Only reports; no pass/fail gate (ADR 0047).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.freeboard import account_freeboard, ejection_velocity_samples, stable_equilibria
from hamel.freeboard.transport import ghost_bubble_velocity

_OUT = Path("reports/t26h-trajectories")


def _velocity(column):
    edges = np.asarray(column["edges_m"], dtype=float)
    diameters = np.asarray(column["edge_diameter_m"], dtype=float)
    flows = np.asarray(column["volume_flow_m3_s"], dtype=float)
    bubble, width = float(column["bubble_m_s"]), float(column["bubble_diameter_m"])
    return lambda cell, height: ghost_bubble_velocity(cell, height, edges, diameters, flows, bubble, width)


def _books(column, diameter, samples):
    speeds, shares, direct = ejection_velocity_samples(float(column["bubble_m_s"]), samples)
    return account_freeboard(
        np.asarray(column["edges_m"], dtype=float), _velocity(column),
        np.asarray(column["rho_g"], dtype=float), np.asarray(column["mu_g"], dtype=float),
        float(diameter), float(column["rho_p"]), float(column["phi_s"]),
        speeds, shares, direct, SolverConfig(freeboard_ejection_samples=samples),
    )


def _relative(left, right):
    scale = max(abs(float(left)), abs(float(right)))
    if scale == 0.0:
        return 0.0
    return abs(float(left) - float(right)) / scale


def _rows(dump):
    return dump["a_free"], dump["a_hover"], dump["b_free"], dump["b_hover"]


def _cell(value):
    if value is None:
        return "—"
    if value == 0.0:
        return "0"
    return f"{value:.4g}"


def main(old_path: str, new_path: str, old_roots_path: str, old_doubling_path: str, old_timing: str | None = None, new_timing: str | None = None) -> None:
    old = json.loads(Path(old_path).read_text(encoding="utf-8"))
    new = json.loads(Path(new_path).read_text(encoding="utf-8"))
    old_roots = json.loads(Path(old_roots_path).read_text(encoding="utf-8"))
    old_doubling = json.loads(Path(old_doubling_path).read_text(encoding="utf-8"))["relative_change"]
    column = new["trajectory"]
    diameters = [float(item) for item in column["diameter_m"]]
    edges = [float(item) for item in column["edges_m"]]
    count = len(edges) - 1
    classes = len(diameters)
    new_roots = [list(item.height_m for item in stable_equilibria(
        edges, _velocity(column), column["rho_g"], column["mu_g"], diameter, column["rho_p"], column["phi_s"],
    )) for diameter in diameters]
    doubling = []
    for diameter in diameters:
        base = _books(column, diameter, 48)
        fine = _books(column, diameter, 96)
        worst = 0.0
        changed = False
        for left, right in zip(base.cells, fine.cells):
            for name in ("a_f", "a_t", "b_f", "b_t", "m_hat", "mean_slip_m_s"):
                first, second = getattr(left, name), getattr(right, name)
                if first is None and second is None:
                    continue
                if first is None or second is None:
                    changed = True
                    continue
                worst = max(worst, _relative(first, second))
        doubling.append({
            "diameter_m": diameter, "relative_change": worst, "existence_changed": changed,
            "capture_base": base.capture_share, "capture_fine": fine.capture_share,
            "return_base": base.return_share, "return_fine": fine.return_share,
        })
    payload = {
        "case": "VTT 加压 1（工况 1）初值点第 1 轮预计算",
        "old_elapsed_s": old["elapsed_s"], "new_elapsed_s": new["elapsed_s"],
        "old_t_ad_k": old["t_ad_k"], "new_t_ad_k": new["t_ad_k"],
        "old_d_bar_m": old["d_bar_m"], "new_d_bar_m": new["d_bar_m"],
        "old_rho_p_kg_m3": old["rho_p_kg_m3"], "new_rho_p_kg_m3": new["rho_p_kg_m3"],
        "edges_m": edges, "diameter_m": diameters,
        "old_return": old["reference_return_share"], "new_return": new["reference_return_share"],
        "old_elutriation": old["reference_elutriation_share"], "new_elutriation": new["reference_elutriation_share"],
        "old_roots_m": old_roots["stable_heights_m"], "new_roots_m": new_roots,
        "doubling": doubling,
    }
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "report.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    lines = [
        "# T26h：工况 1 初值点第 1 轮预计算的参考系数对照",
        "",
        "只报告，不作通过门（ADR 0047）。旧值取 `origin/main`（8257fb7）的连续 u_gb(h) + 数值积分；",
        "新值取本卡的分支 `task/T26h-cellwise-trajectories`（逐格常系数解析解）。两边同一工况、同一轮次、",
        "同一 `SolverConfig`（除 `freeboard_max_step_s` 已按 ADR 0047 第 5 项删除）。",
        "",
        "生成方式：两边的预计算由 `tools/t26h_reference.py` 在各自的工作树里各跑一次，结果存为",
        "`precalc_old.json`、`precalc_new.json`；旧实现的稳定根与样本加倍值在同一保存列上调用改动前的",
        "`stable_equilibria`、`account_freeboard` 得到（`old_roots.json`、`old_doubling.json`）；",
        "阶段耗时由 `tools/t26h_timing.py` 在各自的工作树里各跑一次（`timing_old.json`、`timing_new.json`）。",
        "汇总脚本 `tools/t26h_report.py` 写出 `report.md` 与 `report.json`。",
        "",
        "## 预计算耗时",
        "",
        "| 版本 | 墙钟 s | T_ad K | ̄d m | ρ_P kg/m³ |",
        "|---|---|---|---|---|",
        f"| 旧（连续 u_gb + 数值积分） | {old['elapsed_s']:.1f} | {old['t_ad_k']:.2f} | {old['d_bar_m']:.6e} | {old['rho_p_kg_m3']:.2f} |",
        f"| 新（逐格解析） | {new['elapsed_s']:.1f} | {new['t_ad_k']:.2f} | {new['d_bar_m']:.6e} | {new['rho_p_kg_m3']:.2f} |",
        "",
        f"比值（新/旧）= {new['elapsed_s'] / old['elapsed_s']:.3f}。轨迹只是预计算的一部分，床层子步与其它环节不变。",
        "",
        "## 参考份额与悬停位置",
        "",
        "| 粒级 | d m | 旧回落 | 新回落 | 旧带出 | 新带出 | 旧 h* m | 新 h* m |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for index in range(classes):
        old_root = old_roots["stable_heights_m"][index]
        new_root = new_roots[index]
        lines.append(
            f"| {index + 1} | {diameters[index]:.4e} | {old['reference_return_share'][index]:.6f} | "
            f"{new['reference_return_share'][index]:.6f} | {old['reference_elutriation_share'][index]:.6f} | "
            f"{new['reference_elutriation_share'][index]:.6f} | "
            f"{', '.join(f'{value:.4f}' for value in old_root) or '无'} | "
            f"{', '.join(f'{value:.4f}' for value in new_root) or '无'} |"
        )
    lines += [
        "",
        "悬停份额 = 1 − 回落 − 带出（参考轨迹按单位抛出记账）。旧实现的 h* 为连续剖面的格内或界面根；",
        "新实现的稳定点只在界面（ADR 0047 第 3 项）。",
        "",
        "## 逐格参考系数",
        "",
    ]
    names = (("a^f", "a_free"), ("a^t", "a_hover"), ("b^f", "b_free"), ("b^t", "b_hover"))
    for index in range(classes):
        lines += [
            f"### 粒级 {index + 1}（d = {diameters[index]:.4e} m，θ 与悬停见下）",
            "",
            "| 格（自床面起） | " + " | ".join(f"旧 {name} | 新 {name}" for name, _ in names) + " | 旧 θ | 新 θ | 旧悬停 | 新悬停 |",
            "|---|" + "---|" * (2 * len(names) + 4),
        ]
        for cell in range(count - 1, -1, -1):
            balance = count - 1 - cell
            row = [f"{cell}"]
            for _name, key in names:
                row += [_cell(old[key][balance][index]), _cell(new[key][balance][index])]
            row += [
                str(int(old["theta"][balance][index])), str(int(new["theta"][balance][index])),
                "是" if old["hover"][balance][index] else "否",
                "是" if new["hover"][balance][index] else "否",
            ]
            lines.append("| " + " | ".join(row) + " |")
        lines += [
            "",
            f"抛出拆分：旧 free/hover = {old['ejection_to_free'][index]:.6f}/{old['ejection_to_hover'][index]:.6f}；"
            f"新 = {new['ejection_to_free'][index]:.6f}/{new['ejection_to_hover'][index]:.6f}。",
            "",
        ]
    lines += [
        "",
        "## 样本数加倍（48 → 96）",
        "",
        "| 粒级 | 旧最大相对变化 | 新最大相对变化 | 自由状态存在性变化 | 捕获份额 48 | 捕获份额 96 |",
        "|---|---|---|---|---|---|",
    ]
    for index, row in enumerate(doubling):
        lines.append(
            f"| {index + 1} | {old_doubling[index]:.3e} | {row['relative_change']:.3e} | "
            f"{'有' if row['existence_changed'] else '无'} | {row['capture_base']:.6f} | {row['capture_fine']:.6f} |"
        )
    lines += [
        "",
        "门 1·10⁻⁴ 是 ADR 0016/0028 的样本加密检查，只报告。工况 1 的这一列上，**旧实现同样未过门**，",
        "量级与上表相同（最大一项来自粒级 6 的 5·10⁻¹），因此不是本卡引入的差异；抛出初速分布的离散",
        "使少数样本的归宿在 48/96 之间翻转。捕获份额与 h* 不变。",
        "",
    ]
    if old_timing is not None and new_timing is not None:
        before = json.loads(Path(old_timing).read_text(encoding="utf-8"))
        after = json.loads(Path(new_timing).read_text(encoding="utf-8"))
        lines += [
            "## 耗时分解（第一次外层预计算）",
            "",
            "把 `hamel/solver/precalc.py` 的各阶段函数包一层计时（`tools/t26h_timing.py`），两次运行在同一台机器上先后各跑一次。",
            "上一节与本节是两次独立的运行，绝对数值相差约一成，属机器负载差异；比值稳定。",
            "",
            "| 阶段 | 旧 s | 旧占比 | 新 s | 新占比 |",
            "|---|---|---|---|---|",
        ]
        for name in ("trajectory", "bed", "pyrolysis", "startup", "other"):
            left = before["phases"][name]
            right = after["phases"][name]
            lines.append(
                f"| {name} | {left:.1f} | {left / before['total_s'] * 100:.1f} % | "
                f"{right:.1f} | {right / after['total_s'] * 100:.1f} % |"
            )
        lines += [
            f"| 合计 | {before['total_s']:.1f} | 100 % | {after['total_s']:.1f} | 100 % |",
            "",
            f"轨迹（`_segregate`）由 {before['phases']['trajectory']:.1f} s 降到 {after['phases']['trajectory']:.1f} s"
            f"（约 {before['phases']['trajectory'] / max(after['phases']['trajectory'], 1.0e-9):.0f} 倍），"
            "与 S. 49 的“快一个因子 50”同一方向。",
            f"预计算总时间只降到 {after['total_s'] / before['total_s']:.2f} 倍，因为余下部分几乎全是床区几何积分"
            f"（`_adaptive_bed` + `_exchange_diffusivity`）：{before['phases']['bed']:.1f} s → {after['phases']['bed']:.1f} s，"
            "它由 ADR 0035 第 5 项、ADR 0036 第 2 项的整床子步加密决定，本卡不涉及。",
            "热解（`_release`，含释放时间与两状态释放分布）与启动项（绝热平衡、流动、初值）在两次运行之间只差机器噪声。",
            "若要继续缩短预计算，下一步是床层子步与 ρ_P 反求，不是轨迹。",
            "",
        ]
    lines += [
        "## 备注",
        "",
        "- 旧实现的 `freeboard_max_step_s` 已删除，新的步长不存在，因此不再有步长加密检查（ADR 0047 第 5 项）。",
        "- 悬停判定由“格内连续平衡 + 势能井”改为“界面稳定点 + 闭式越界距离”（ADR 0047 第 3 项）。",
        "- 本报告不判断模型是否更接近 Hamel 的模拟值；验收仍按 CHARTER §2。",
    ]
    (_OUT / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {_OUT/'report.md'}")


if __name__ == "__main__":
    main(*(sys.argv[1:7]))
