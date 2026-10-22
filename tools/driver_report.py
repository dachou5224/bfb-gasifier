"""把一次运行表写成 Markdown 和 JSON。数字来自调用方，这里不求解。"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.cases.judge import Measurement, Premises, judge
from hamel.cases.targets import acceptance_targets

_CASE_ID = "CASE_VTT_PRESSURIZED_1"
_DIGITIZED = "data/digitized"
_FROZEN = (
    ("F1", "式 (3.43) 的 −1 为误印；α_b = u_b·ε_d/u_d", "ADR 0011 第 7、10 条；ADR 0022 A1"),
    ("F2", "R8 用式 (5.45)，a_R8 = 0.02", "ADR 0013 第 3 项"),
    ("F3", "Serio 第一对", "ADR 0013 第 4 项"),
    ("F4", "料斗吹扫氮计入，298.15 K", "ADR 0013 第 5 项；ADR 0014 C2"),
    ("F5", "8×0.15 m + 8×0.2875 m", "ADR 0014 B"),
    ("F6", "床高 1.2 m；床料 3.963 kg；ρ_P 外层反求", "ADR 0014 A2–A4"),
    ("F7", "壁面散热 21.7 kW，单一 k_W，T_amb = 298.15 K", "ADR 0014 X1；ADR 0029"),
    ("F8", "旋风纯分离；连接管不计反应与传热", "ADR 0014 X2、X2a"),
    ("F9", "第一块布风板，n_0 = 139.5，h_A = 0", "ADR 0019"),
    ("F10", "ε_mf = 0.46，φ_S = 0.8", "ADR 0024"),
    ("F11", "粒级代表径为 Sauter 平均，最细档下限 0.0291 mm", "ADR 0023"),
    ("F12", "ε_P = 0.85", "ADR 0017 第 17 项"),
    ("F13", "D_g 用 Hamel 式 (5.20)，cell 气体温度", "ADR 0021 第 4 项"),
    ("F14", "灰的比热 950 J/(kg·K)", "ADR 0018 第 6 项"),
    ("F15", "ρ_C = ρ_P", "ADR 0017 第 19 项"),
    ("F16", "自由空间焦油反应只用 Serio", "ADR 0022 C3"),
    ("F17", "悬停闭合为反应限制的滞留", "ADR 0026"),
    ("F18", "热解时间用平均释放时间", "ADR 0015 第 10 项；ADR 0027"),
    ("F19", "基准初值扣除 21.7 kW；只有 S1 不扣除", "ADR 0017 第 9 项；ADR 0029"),
    ("F20", "内层 1e-6；外层 (a)(b)；元素 0.1 %；能量 0.5 %", "ADR 0012"),
    ("F21", "复合颗粒 shrinking particle", "ADR 0031"),
    ("F22", "燃料用同粒级焦炭的系数", "ADR 0032"),
    ("F23", "平均泡间距用式 (3.57)", "ADR 0028"),
    ("F24", "床区 u_r = u_d/ε_d；自由状态为轨迹加权滑移；悬停为 v_t", "ADR 0035 第 4 项"),
    ("F25", "未收敛时非 B 类用 B 类系数延拓；收敛后按真实类别复核", "ADR 0035 第 2 项"),
)


def write_report(report_dir, records, *, title: str, toy: bool) -> dict:
    payload = _document(records, title, toy)
    folder = Path(report_dir)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "driver_report.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (folder / "driver_report.md").write_text(_markdown(payload), encoding="utf-8")
    return payload


def convert(value):
    """数组先展开成同形状的列表，再把 NumPy 标量收成 Python 标量。"""
    if isinstance(value, dict):
        return {str(key): convert(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return convert(value.tolist())
    if isinstance(value, (list, tuple)):
        return [convert(item) for item in value]
    if isinstance(value, np.generic):
        return convert(value.item())
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    return str(value)


def _document(records, title, toy) -> dict:
    root = Path(__file__).resolve().parents[1]
    targets = acceptance_targets(_CASE_ID, root / _DIGITIZED)
    judged = [_judge(record, targets) for record in records]
    return {
        "title": title, "toy": toy, "timeout_s": SolverConfig().wall_timeout_s,
        "frozen": [{"id": item[0], "choice": item[1], "source": item[2]} for item in _FROZEN],
        "numerical": _numerical(), "runs": judged, "differences": _differences(judged),
    }


def _judge(record, targets) -> dict:
    measured = {}
    for target in targets:
        if not target.present:
            continue
        value = (record.get("measurements") or {}).get(target.quantity)
        if value is not None:
            measured[target.quantity] = Measurement(float(value), target.unit, target.basis)
    premises = Premises(
        residuals_met=record.get("status") == "CONVERGED",
        balance_closed=bool((record.get("balance") or {}).get("ok")),
        state_finite_nonnegative=bool(record.get("nonnegative")),
    )
    case = judge(_CASE_ID, measured, targets, premises)
    body = dict(record)
    body["acceptance"] = {
        "status": case.status, "premises_met": case.premises_met,
        "quantities": [
            {"quantity": item.quantity, "status": item.status, "value": item.value, "distance": item.distance_to_boundary}
            for item in case.quantities
        ],
    }
    return convert(body)


def _differences(runs) -> list:
    baseline = next((item for item in runs if item["name"] == "baseline"), None)
    base = None if baseline is None else (baseline.get("measurements") or {})
    rows = []
    for item in runs:
        if item["name"] == "baseline":
            continue
        current = item.get("measurements") or {}
        delta = {}
        if base and current and item.get("status") == "CONVERGED" and baseline.get("status") == "CONVERGED":
            for key, value in current.items():
                if key in base and value is not None and base[key] is not None:
                    delta[key] = float(value) - float(base[key])
        rows.append({"name": item["name"], "status": item["status"], "delta": delta, "reason": item.get("reason", "")})
    return rows


def _numerical() -> dict:
    config = SolverConfig()
    names = (
        "pyrolysis_radial_intervals", "pyrolysis_time_rtol", "freeboard_ejection_samples",
        "freeboard_max_step_s", "freeboard_delta_h_m", "freeboard_delta_u_m_s",
        "particle_tolerance_w", "particle_max_iterations", "coupling_tolerance_k",
        "coupling_max_iterations", "bed_integration_intervals_per_cell",
        "bed_integration_max_intervals_per_cell",
    )
    noted = {name: getattr(config, name) for name in names}
    noted["note"] = "逐根 δ 与床层子步 n 来自该次运行的预计算。玩具工况没有这些量。"
    return noted


def _cell(value):
    if isinstance(value, float):
        return f"{value:.6g}"
    return value


def _markdown(payload) -> str:
    lines = [f"# {payload['title']}", ""]
    if payload["toy"]:
        lines += [
            "这是玩具工况，用来核对驱动、收敛门和报告。数字不是第 1 关的解。",
            "S1、S2、S6、S7、S8 已能通过现有接口表达。样例不启动工况进程，表中不列这几次的装置结果。",
            "",
        ]
    lines += [f"单次运行墙钟上限 {payload['timeout_s']} s（ADR 0037）。", "", "## 验收", ""]
    for run in payload["runs"]:
        acceptance = run.get("acceptance") or {}
        lines += [f"### {run['name']}：{run['status']}", "", f"改动：{run.get('change', '')}。来源：{run.get('source', '')}。"]
        if run.get("reason"):
            lines.append(f"说明：{run['reason']}。")
        premise = "满足" if acceptance.get("premises_met") else "不满足"
        lines += [
            f"验收状态：{acceptance.get('status', '未判定')}。前提：{premise}。", "",
            "| 量 | 状态 | 值 | 到边界的距离 |", "|---|---|---|---|",
        ]
        for item in acceptance.get("quantities") or []:
            lines.append(
                f"| {item['quantity']} | {item['status']} | {_cell(item['value'])} | {_cell(item['distance'])} |"
            )
        balance = run.get("balance") or {}
        lines += ["", f"外部平衡：{balance.get('ratios')}。非负：{run.get('nonnegative')}。残差比：{run.get('ratios')}。", ""]
    lines += [
        "## 诊断", "",
        "每次运行的 JSON `diagnostics` 写入 T25b 公开结果里已有的量：Geldart 记录与式 (8) 距离、T_ad、T_init、y_ad、ρ_P、排灰与灰分量、抛出与回落带出系数、悬停标记与 h* 所在的逐根 δ、ū_r/v_t、产率越界与延拓、床层子步 n、k_W。",
        "气泡相总流、d_b/D、与 KURKELA 表 5 的气速对照、自由空间元素与质量闭合、数值积分和闭式的比对、悬停焦炭消耗占比，不在这份公开结果上。本卡不改 T25b 的接口，这些量留空，不另算。",
        "S9–S12 未接线，因此 ρ_P 可行区间没有对应的运行。",
        "", "## 敏感性差值", "", "| 运行 | 状态 | 差值 |", "|---|---|---|",
    ]
    for row in payload["differences"]:
        lines.append(f"| {row['name']} | {row['status']} | {row['delta'] or row['reason']} |")
    lines += ["", "## 数值加密", ""]
    for key, value in payload["numerical"].items():
        lines.append(f"- {key}：{value}")
    lines += ["", "## 冻结选择", "", "| 编号 | 选择 | 来源 |", "|---|---|---|"]
    for item in payload["frozen"]:
        lines.append(f"| {item['id']} | {item['choice']} | {item['source']} |")
    lines.append("")
    return "\n".join(lines)
