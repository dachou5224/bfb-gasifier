"""把一次运行表写成 Markdown 和 JSON。数字来自调用方，这里不求解。"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.cases.judge import Measurement, Premises, judge
from hamel.cases.targets import acceptance_targets, half_width
from tools.run_prep import _ACCEPTANCE, closed_form_report, constructed_books, refinement_report, toy_applicability_rows

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
    ("F26", "T̄_FRB 按自由空间格子高度加权，含锥段", "ADR 0040 第 1 项"),
    ("F27", "X_C 为出口气体（含焦油）碳／燃料碳", "ADR 0040 第 1 项"),
    ("F28", "床区流体力学气体物性用摩尔流量加权平均；D_g 逐格", "ADR 0040 第 2 项"),
)
_KINDS = {
    "T_FRB_K": "temperature", "X_C": "carbon_conversion", "CH4": "ch4",
    "CO": "major_gas", "CO2": "major_gas", "H2": "major_gas", "H2O": "major_gas",
}
_SHARED = "各验收量共享壁面散热这一校准输入。"
_TOY_EXTRA = None


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
        "shared_input": _SHARED,
        "numerical": _numerical(), "runs": judged, "differences": _differences(judged, targets),
        "closed_form": _recorded(records, "closed_form", toy, "closed_form"),
        "refinement": _recorded(records, "refinement", toy, "refinement"),
    }


def _recorded(records, key, toy: bool, extra: str):
    for record in records:
        value = (record.get("diagnostics") or {}).get(key)
        if value:
            return convert(value)
    if not toy:
        return None
    return _toy_extra()[extra]


def _toy_extra() -> dict:
    global _TOY_EXTRA
    if _TOY_EXTRA is None:
        config = SolverConfig()
        example = constructed_books()
        _TOY_EXTRA = {
            "closed_form": closed_form_report(config, toy_applicability_rows()),
            "refinement": refinement_report(config, particle=(example["state"], example["frozen"], 1.0, 1.0, 1.0)),
            "diagnostics": example["diagnostics"],
            "measurements": example["measurements"],
        }
    return _TOY_EXTRA


def publish_toy_report(report_dir) -> dict:
    """写出 reports/toy。基准和 S1 用构造数字，S2–S16 只列接线。"""
    from hamel.solver.campaign import run_specs
    from hamel.solver.sensitivity import pack_reported

    extra = _toy_extra()
    profile = SimpleNamespace(
        epsilon_b=np.array([0.2, 0.1]), u_b=np.array([1.5, 1.2]), d_b=np.array([0.2, 0.05]),
    )
    flow = SimpleNamespace(
        u0=np.array([0.8, 0.7]), free_mol_s=np.array([[0.02, 0.01]]), free_temperature=np.array([1000.0]),
        pressure=101300.0,
    )
    released = SimpleNamespace(
        mass_gap_kg_s=np.array([1.0e-6, -2.0e-6, 0.0]),
        element_gap_kg_s={"C": np.array([1.0e-7, 0.0, 3.0e-7])},
    )
    case = SimpleNamespace(bed_diameter_m=0.15, free_diameter_m=0.25, n_free=2)
    diagnostics = dict(extra["diagnostics"])
    diagnostics["sensitivity"] = pack_reported(profile, flow, case, released)
    measurements = dict(extra["measurements"])
    shifted = dict(measurements)
    shifted["T_FRB_K"] = 1115.0
    records = []
    for spec in run_specs():
        blank = {
            "name": spec.name, "change": spec.change, "source": spec.source,
            "rounds": 0, "iterations": 0, "ratios": {}, "diagnostics": {},
        }
        if spec.name == "baseline":
            records.append({
                **blank, "status": "CONVERGED", "reason": "构造例，不是第 1 关的解",
                "balance": {"ok": True, "ratios": {"C": 0.0, "H": 0.0, "O": 0.0, "energy": 0.0}},
                "nonnegative": True, "measurements": measurements, "diagnostics": diagnostics,
            })
        elif spec.name == "S1":
            records.append({
                **blank, "status": "CONVERGED", "reason": "构造例，k_W = 0 的并列数字，不是装置结果",
                "balance": {"ok": True, "ratios": {"C": 0.0, "H": 0.0, "O": 0.0, "energy": 0.0}},
                "nonnegative": True, "measurements": shifted, "diagnostics": {},
            })
        else:
            records.append({
                **blank, "status": "WIRED", "reason": "已接线。样例不启动工况进程，不列装置结果。",
                "balance": None, "nonnegative": None, "measurements": None,
            })
    return write_report(report_dir, records, title="玩具工况", toy=True)


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


def _over_half(quantity: str, delta: float, targets) -> bool:
    kind = _KINDS.get(quantity)
    target = next((item for item in targets if item.quantity == quantity and item.present and item.sim is not None), None)
    if kind is None or target is None:
        return False
    return abs(float(delta)) > 0.5 * half_width(kind, float(target.sim))


def _differences(runs, targets) -> list:
    baseline = next((item for item in runs if item["name"] == "baseline"), None)
    base = None if baseline is None else (baseline.get("measurements") or {})
    rows = []
    for item in runs:
        if item["name"] == "baseline":
            continue
        current = item.get("measurements") or {}
        delta = {}
        marked = []
        if base and current and item.get("status") == "CONVERGED" and baseline.get("status") == "CONVERGED":
            for key, value in current.items():
                if key in base and value is not None and base[key] is not None:
                    gap = float(value) - float(base[key])
                    delta[key] = gap
                    if _over_half(key, gap, targets):
                        marked.append(key)
        rows.append({
            "name": item["name"], "status": item["status"], "delta": delta,
            "over_half_tolerance": marked, "reason": item.get("reason", ""),
        })
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


def _sensitivity_example(runs) -> dict | None:
    for run in runs:
        sensitivity = (run.get("diagnostics") or {}).get("sensitivity")
        if sensitivity:
            return sensitivity
    return None


def _cell(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return value


def _run(runs, name):
    return next((item for item in runs if item.get("name") == name), None)


def _measured(run, quantity):
    if run is None:
        return None
    return (run.get("measurements") or {}).get(quantity)


def _acceptance_pair(runs) -> list:
    baseline, adiabatic = _run(runs, "baseline"), _run(runs, "S1")
    lines = [
        "基准与 S1 的验收量：", "",
        "| 量 | 基准 | S1 |", "|---|---|---|",
    ]
    for name in _ACCEPTANCE:
        lines.append(f"| {name} | {_cell(_measured(baseline, name))} | {_cell(_measured(adiabatic, name))} |")
    lines.append("")
    return lines


def _stocks(runs) -> list:
    lines = []
    found = False
    for run in runs:
        stocks = (run.get("diagnostics") or {}).get("stocks")
        if not stocks:
            continue
        found = True
        lines.append(f"### {run['name']}")
        lines.append("")
        if stocks.get("note"):
            lines.append(stocks["note"])
            lines.append("")
        deviation = (run.get("diagnostics") or {}).get("bed_gas_temperature_deviation_k")
        if deviation is not None:
            lines.append(f"床区气体温度与平均值的最大偏差：{_cell(float(deviation))} K。")
            lines.append("")
        lines += ["| 粒级 | 参考回落份额 | 参考带出份额 | 实际回落 kg/s | 实际带出 kg/s |", "|---|---|---|---|---|"]
        for row in stocks.get("classes") or []:
            lines.append(
                f"| {row['class']} | {_cell(row['reference_return_share'])} | {_cell(row['reference_elutriation_share'])} | {_cell(row['return_kg_s'])} | {_cell(row['elutriation_kg_s'])} |"
            )
        lines += ["", "| 格 | 固体 kg | 反应消耗 kg/s |", "|---|---|---|"]
        for row in stocks.get("cells") or []:
            lines.append(f"| {row['cell']} | {_cell(row['solid_kg'])} | {_cell(row['char_consumption_kg_s'])} |")
        lines += ["", "| 格 | 粒级 | 悬停 kg | 占该格固体 |", "|---|---|---|---|"]
        for row in stocks.get("hover") or []:
            lines.append(f"| {row['cell']} | {row['class']} | {_cell(row['hover_kg'])} | {_cell(row['hover_fraction'])} |")
        lines.append("")
    if not found:
        lines += ["这次运行没有收敛存量。", ""]
    return lines


def _stream_table(closure) -> list:
    lines = ["| 流 | 方向 | 质量 kg/s | C | H | O |", "|---|---|---|---|---|---|"]
    for row in closure.get("streams") or []:
        lines.append(
            f"| {row['name']} | {row['direction']} | {_cell(row['mass'])} | {_cell(row['C'])} | {_cell(row['H'])} | {_cell(row['O'])} |"
        )
    lines += ["", "| 量 | 入 | 出 | 差值 | 相对差 |", "|---|---|---|---|---|"]
    for key, row in (closure.get("totals") or {}).items():
        lines.append(f"| {key} | {_cell(row['in'])} | {_cell(row['out'])} | {_cell(row['difference'])} | {_cell(row['relative'])} |")
    lines.append("")
    lines.append(f"气固反应消耗的焦炭（不进入上表）：{_cell(closure.get('reaction_char_kg_s'))} kg/s。")
    lines.append("")
    return lines


def _closure_and_form(payload) -> list:
    lines = ["## 全自由空间闭合", ""]
    if payload["toy"]:
        lines += ["下表是构造状态的逐股记账。反应项没有加进差值。构造状态不是收敛解，差值不必为零。", ""]
    found = False
    for run in payload["runs"]:
        closure = (run.get("diagnostics") or {}).get("freeboard_closure")
        if not closure:
            continue
        found = True
        lines += [f"### {run['name']}", ""]
        lines.extend(_stream_table(closure))
    if not found:
        lines += ["这次运行没有闭合记账。", ""]
    form = payload.get("closed_form")
    lines += ["## 数值积分与闭式", ""]
    if not form:
        lines += ["这次报告没有做闭式比对。", ""]
        return lines
    benchmark = form["benchmark"]
    lines += [
        f"固定构造基准：Δ = {_cell(benchmark['delta'])}，步长 {_cell(benchmark['step_s'])} s，样本数设置 {benchmark['samples']}。",
        f"最大相对偏差 {_cell(benchmark['max_relative'])}（速度 {_cell(benchmark['velocity_relative'])}，高度 {_cell(benchmark['height_relative'])}）。",
        "",
        "| 粒径 m | 分支 | Δ | 较快支 Δ | 较慢支 Δ | 状态 | 截止 | 速度相对偏差 | 高度相对偏差 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in form["classes"]:
        lines.append(
            f"| {_cell(row['diameter_m'])} | {row['branch']} | {_cell(row['delta'])} | {_cell(row['delta_faster'])} | {_cell(row['delta_slower'])} | {row['status']} | {row['stopped']} | {_cell(row['max_velocity_relative'])} | {_cell(row['max_height_relative'])} |"
        )
    lines.append("")
    return lines


def _refinement_lines(payload) -> list:
    lines = []
    refinement = payload.get("refinement")
    if refinement:
        trajectory = refinement["trajectory"]
        lines.append(
            f"- 轨迹样本加倍、步长减半：a、b 相对变化 {_cell(trajectory['relative_change'])}，门 {trajectory['gate']}，{'满足' if trajectory['met'] else '未满足'}。"
        )
        pyrolysis = refinement["pyrolysis"]
        if pyrolysis.get("diameter_m") is not None:
            lines.append(f"- 热解粒径 {_cell(pyrolysis['diameter_m'])} m，温度 {_cell(pyrolysis.get('temperature_k'))} K。")
        for row in pyrolysis.get("rows") or []:
            if "name" in row:
                lines.append(
                    f"- 热解 {row['name']}：相对变化 {_cell(row['relative_change'])}，门 {row['gate']}，{'满足' if row['met'] else '未满足'}。"
                )
            else:
                lines.append(
                    f"- 热解粒径 {_cell(row.get('diameter_m'))} m：τ_pyr = {_cell(row.get('tau_pyr_s'))} s，温度 {_cell(row.get('temperature_k'))} K。"
                )
        particle = refinement.get("particle")
        if particle:
            lines.append(
                f"- 颗粒温度（构造状态，不是收敛解）：最大变化 {_cell(particle['temperature_change_k'])} K，重算残差相对范数 {_cell(particle['residual_relative_norm'])}。温度门 0.1 K，残差门 1e-6，{'满足' if particle['met'] else '未满足'}。构造状态的残差不要求过门。"
            )
    for run in payload["runs"]:
        info = ((run.get("diagnostics") or {}).get("final") or {})
        if "bed_intervals" not in info and "substep_history" not in info and "delta_notes" not in info:
            continue
        lines.append(f"- {run['name']} 的 ADR 0036：床层子步 n = {info.get('bed_intervals')}，加倍记录 {info.get('substep_history')}，逐根 δ {info.get('delta_notes')}。")
    if not any("ADR 0036" in line for line in lines):
        lines.append("- ADR 0036 的逐根 δ 与床层子步：这次运行没有这两项记录。")
    return lines


def _markdown(payload) -> str:
    lines = [f"# {payload['title']}", ""]
    if payload["toy"]:
        lines += [
            "这是玩具工况，用来核对驱动、收敛门和报告。数字不是第 1 关的解。",
            "S1–S16 都已接线。样例不启动工况进程。基准与 S1 的数字来自构造例。S2–S16 不列装置结果。",
            "",
        ]
    lines += [
        f"单次运行墙钟上限 {payload['timeout_s']} s（ADR 0037）。", "",
        "## 验收", "", payload["shared_input"], "",
    ]
    lines.extend(_acceptance_pair(payload["runs"]))
    for run in payload["runs"]:
        acceptance = run.get("acceptance") or {}
        if payload["toy"] and run.get("status") == "WIRED":
            lines += [
                f"### {run['name']}：已接线", "",
                f"改动：{run.get('change', '')}。来源：{run.get('source', '')}。",
                "样例不启动工况进程，不列装置结果。", "",
            ]
            continue
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
        log = run.get("stderr_log")
        if run.get("status") == "ERROR" and log:
            lines.append(f"stderr 日志：{log}。")
            lines.append("")
        used = (run.get("diagnostics") or {}).get("used")
        if used:
            lines += ["实际使用的冻结值：", ""]
            for key, value in used.items():
                lines.append(f"- {key}：{_cell(value)}")
            lines.append("")
    lines += [
        "## 存量与流量", "",
        "参考回落与带出份额来自无反应轨迹。实际流量是收敛存量乘输运系数。悬停占比按 ADR 0026。",
        "",
    ]
    lines.extend(_stocks(payload["runs"]))
    lines += [
        "## 诊断", "",
        "每次运行的 JSON `diagnostics` 写入预计算已经算出的量，以及求解后由这些量整理的诊断。",
        "气泡相总流之比、d_b/D 与 d_b>D 的格数、与 KURKELA 表 5 的气速对照、悬停焦炭消耗占比，都取自这些中间量。",
        "床区各格气体温度与摩尔流量加权平均温度的最大偏差写在 `bed_gas_temperature_deviation_k`。",
        "全自由空间闭合按 ADR 0040 第 5 项逐股记账，不加反应项。闭式比对只在报告工具里做，不设通过门。",
        "",
    ]
    example = _sensitivity_example(payload["runs"])
    if example:
        lines += ["构造例不启动工况进程。已有中间量仍由打包函数给出。", ""]
        for key, value in example.items():
            lines.append(f"- {key}：{_cell(value) if not isinstance(value, (dict, list)) else value}")
        lines.append("")
    lines.extend(_closure_and_form(payload))
    lines += ["## 敏感性差值", "", "| 运行 | 状态 | 差值 | 超过容差一半 |", "|---|---|---|---|"]
    for row in payload["differences"]:
        marked = "、".join(row.get("over_half_tolerance") or []) or ""
        lines.append(f"| {row['name']} | {row['status']} | {row['delta'] or row['reason']} | {marked} |")
    lines += ["", "## 数值加密", ""]
    for key, value in payload["numerical"].items():
        lines.append(f"- {key}：{value}")
    lines.extend(_refinement_lines(payload))
    lines += ["", "## 冻结选择", "", "| 编号 | 选择 | 来源 |", "|---|---|---|"]
    for item in payload["frozen"]:
        lines.append(f"| {item['id']} | {item['choice']} | {item['source']} |")
    lines.append("")
    return "\n".join(lines)
