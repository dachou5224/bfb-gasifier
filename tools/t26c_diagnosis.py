"""T26c 失败点上的五项诊断。只重算指定比较，不改牛顿步、容差或冻结量。"""
from __future__ import annotations

import json
import math
import pickle
from dataclasses import replace
from pathlib import Path

import numpy as np

from hamel.kinetics.common import KineticsInputError
from hamel.properties.data import PropertyError
from hamel.solver.driver import _feasible, _infeasible_clause
from hamel.solver.frozen import AssemblyError, FreeboardDomainError
from hamel.solver.gates import ResidualBlocks
from hamel.solver.newton import _SQRT_EPS, _at, normalized_merit
from hamel.solver.residual import evaluate_residual

_CAUGHT = (AssemblyError, FreeboardDomainError, KineticsInputError, PropertyError)
_AGREE = 1.0e-2
_H = tuple(10.0 ** (-power) for power in range(1, 13))
_MULTIPLIERS = (1.0, 1.0e2, 1.0e4)
_ROOT = Path(__file__).resolve().parents[1]


def segment_root(start: float, end: float):
    """start + h·(end − start) = 0 在 (0, 1] 上的 h。start 为 0 时返回 0。"""
    if start == 0.0:
        return 0.0
    delta = end - start
    if delta == 0.0:
        return None
    step = -start / delta
    if 0.0 < step <= 1.0:
        return float(step)
    return None


def diagnose(record_path, destination) -> dict:
    record = json.loads(Path(record_path).read_text(encoding="utf-8"))
    folder = _resolve(record["diagnostics"]["failure_npz"]).parent
    with np.load(folder / "failure.npz") as saved:
        payload = {key: np.array(saved[key]) for key in saved.files}
    with (folder / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    blocks = ResidualBlocks(
        payload["gas"], payload["solid"], payload["energy"],
        float(payload["n_dot_ref"]), float(payload["m_dot_ref"]), float(payload["q_dot_ref"]),
    )
    names = [str(item) for item in payload["names"].tolist()]
    x, residual, matrix, delta = payload["x"], payload["F"], payload["J"], payload["dx"]
    scales = payload["scales"]
    fresh = _evaluate(frozen, x)
    direction = _direction(frozen, blocks, x, residual, matrix, delta)
    noise = _noise(frozen, blocks, x, residual, matrix, scales)
    columns = _columns(frozen, blocks, x, residual, matrix, scales, names, record)
    switches = _switches(frozen, x, delta)
    trials = _trials(record)
    precomputed = _precomputed(record, frozen)
    body = {
        "failure_npz": _relative(folder / "failure.npz"),
        "frozen_path": _relative(folder / "frozen.pkl"),
        "repeated_residual_norm": None if fresh is None else normalized_merit(fresh - residual, blocks),
        "repeated_residual_reason": None if fresh is not None else "同一冻结量下重算 F(x) 没有返回",
        "direction": direction,
        "noise": noise,
        "columns": columns,
        "switches": switches,
        "trials": trials,
        "precomputed": precomputed,
        "hypotheses": _hypotheses(direction, noise, switches, trials, record),
    }
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "diagnosis.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    (destination / "diagnosis.md").write_text(_markdown(record, body), encoding="utf-8")
    return body


def _evaluate(frozen, values):
    try:
        return evaluate_residual(np.asarray(values, dtype=float), frozen).values
    except _CAUGHT:
        return None


def _probe(frozen, point):
    def evaluate(trial):
        return evaluate_residual(np.asarray(trial, dtype=float), frozen).values

    try:
        values, status, reason = _at(np.asarray(point, dtype=float), evaluate, lambda trial: _feasible(trial, frozen), lambda trial: _infeasible_clause(trial, frozen))
    except _CAUGHT as exc:
        return None, "reject", f"{type(exc).__name__}: {exc}"
    return values, status, reason


def _direction(frozen, blocks, x, residual, matrix, delta) -> dict:
    predicted = matrix @ delta
    merit0 = normalized_merit(residual, blocks)
    scaled_f = residual.astype(float).copy()
    _scale_into(scaled_f, blocks)
    scaled_step = predicted.astype(float).copy()
    _scale_into(scaled_step, blocks)
    denom = float(np.linalg.norm(scaled_f))
    slope = None if denom == 0.0 else float(np.dot(scaled_f, scaled_step) / denom)
    rows = []
    for step in _H:
        values, status, reason = _probe(frozen, x + step * delta)
        row = {"h": step, "status": status, "reason": reason, "merit": None, "relative": None}
        if status == "ok":
            row["merit"] = normalized_merit(values, blocks)
            quotient = (values - residual) / step
            norm = float(np.linalg.norm(predicted))
            row["relative"] = None if norm == 0.0 else float(np.linalg.norm(predicted - quotient) / norm)
        rows.append(row)
    return {"merit": merit0, "slope": slope, "rows": rows}


def _noise(frozen, blocks, x, residual, matrix, scales) -> dict:
    tight = replace(
        frozen,
        particle_tolerance_w=frozen.particle_tolerance_w / 10.0,
        coupling_tolerance_k=frozen.coupling_tolerance_k / 10.0,
        particle_max_iterations=frozen.particle_max_iterations * 2,
        coupling_max_iterations=frozen.coupling_max_iterations * 2,
    )
    values, status, reason = _probe(tight, x)
    steps = _steps(x, scales)
    column_norms = [normalized_merit(matrix[:, index] * steps[index], blocks) for index in range(x.size)]
    finite = [item for item in column_norms if math.isfinite(item)]
    return {
        "status": status,
        "reason": reason,
        "norm": None if status != "ok" else normalized_merit(values - residual, blocks),
        "tolerances": {
            "particle_tolerance_w": tight.particle_tolerance_w,
            "coupling_tolerance_k": tight.coupling_tolerance_k,
            "particle_max_iterations": tight.particle_max_iterations,
            "coupling_max_iterations": tight.coupling_max_iterations,
        },
        "column_norm_min": None if not finite else float(min(finite)),
        "column_norm_median": None if not finite else float(np.median(finite)),
        "column_norm_max": None if not finite else float(max(finite)),
    }


def _columns(frozen, blocks, x, residual, matrix, scales, names, record) -> list:
    iteration = _failed_iteration(record)
    chosen = [item["index"] for item in iteration["largest"]]
    nominal = _steps(x, scales)
    rows = []
    for index in chosen:
        computed = {}
        for factor in _MULTIPLIERS:
            step = nominal[index] * factor
            point = np.array(x, dtype=float, copy=True)
            point[index] += step
            values, status, reason = _probe(frozen, point)
            column = None if status != "ok" else (values - residual) / step
            computed[factor] = {"status": status, "reason": reason, "column": column}
        base = computed[1.0]["column"]
        saved = matrix[:, index]
        saved_relative = None
        if base is not None:
            denom = float(np.linalg.norm(saved))
            saved_relative = None if denom == 0.0 else float(np.linalg.norm(saved - base) / denom)
        relatives = {}
        for factor in _MULTIPLIERS:
            column = computed[factor]["column"]
            if base is None or column is None:
                relatives[str(factor)] = None
                continue
            denom = float(np.linalg.norm(base))
            relatives[str(factor)] = 0.0 if denom == 0.0 else float(np.linalg.norm(column - base) / denom)
        rows.append({
            "index": int(index), "name": names[index] if index < len(names) else str(index),
            "saved_relative": saved_relative,
            "relative": relatives,
            "status": {str(factor): computed[factor]["status"] for factor in _MULTIPLIERS},
            "reason": {str(factor): computed[factor]["reason"] for factor in _MULTIPLIERS},
        })
    return rows


def _switches(frozen, x, delta) -> dict:
    end = np.asarray(x, dtype=float) + np.asarray(delta, dtype=float)
    state0 = _safe_unpack(x, frozen)
    state1 = _safe_unpack(end, frozen)
    char_rows, bounds = [], []
    if state0 is not None and state1 is not None:
        char_rows.extend(_char_switches("床区", state0.bed_char_kg, state1.bed_char_kg, None))
        char_rows.extend(_char_switches("自由", state0.freeboard_char_free_kg, state1.freeboard_char_free_kg, ~frozen.free_absent))
        char_rows.extend(_char_switches("悬停", state0.freeboard_char_hover_kg, state1.freeboard_char_hover_kg, frozen.hover))
        bounds.extend(_sign_switches("床区焦炭", state0.bed_char_kg, state1.bed_char_kg, "m < 0"))
        bounds.extend(_sign_switches("自由焦炭", state0.freeboard_char_free_kg, state1.freeboard_char_free_kg, "m < 0"))
        bounds.extend(_sign_switches("悬停焦炭", state0.freeboard_char_hover_kg, state1.freeboard_char_hover_kg, "m < 0"))
        bounds.extend(_sign_switches("自由灰", state0.freeboard_ash_free_kg, state1.freeboard_ash_free_kg, "m < 0"))
        bounds.extend(_sign_switches("悬停灰", state0.freeboard_ash_hover_kg, state1.freeboard_ash_hover_kg, "m < 0"))
        bounds.extend(_sign_switches("床区气相", state0.bed_gas_mol_s, state1.bed_gas_mol_s, "Ṅ < 0"))
        bounds.extend(_sign_switches("自由空间气相", state0.freeboard_gas_mol_s, state1.freeboard_gas_mol_s, "Ṅ < 0"))
        bounds.extend(_temperature_switches(state0, state1))
        bounds.extend(_inert_switches(frozen, state0, state1))
    feasible0, clause0 = _endpoint_feasible(x, frozen)
    feasible1, clause1 = _endpoint_feasible(end, frozen)
    return {
        "x_feasible": feasible0, "x_clause": clause0, "x_plus_dx_feasible": feasible1, "x_plus_dx_clause": clause1,
        "char": char_rows, "bounds": bounds,
    }


def _trials(record) -> dict:
    iteration = _failed_iteration(record)
    counts = {}
    ok = 0
    for trial in iteration["trials"]:
        if trial["status"] == "ok":
            ok += 1
            key = "ok，归一残差没有严格下降"
        else:
            key = trial["reason"] or trial["status"]
        counts[key] = counts.get(key, 0) + 1
    return {"iteration": iteration["iteration"], "merit": iteration["merit"], "ok": ok, "counts": counts, "n": len(iteration["trials"])}


def _precomputed(record, frozen) -> dict:
    info = next(item["info"] for item in record["diagnostics"]["rounds"] if item["inner"] == "DAMPING_FAILED")
    sensitivity = info.get("sensitivity") or {}
    pyrolysis = sensitivity.get("pyrolysis") or {}
    return {
        "hover": _pairs(frozen.hover),
        "free_absent": _pairs(frozen.free_absent),
        "theta_one": _pairs(np.asarray(frozen.theta) == 1.0),
        "yield_notes": info.get("yield_notes"),
        "pyrolysis_temperature_k": pyrolysis.get("temperature_k"),
        "pyrolysis_pressure_pa": pyrolysis.get("pressure_pa"),
        "geldart": info.get("geldart"),
        "hydro_group": info.get("hydro_group"),
        "geldart_record": info.get("geldart_record"),
        "bottom_residual_kg_s": info.get("bottom_residual_kg_s"),
    }


def _hypotheses(direction, noise, switches, trials, record) -> list:
    iteration = _failed_iteration(record)
    merit = float(iteration["merit"])
    largest = iteration["largest"][0]
    share = 0.0 if merit == 0.0 else (float(largest["value"]) ** 2) / (merit ** 2)
    noise_norm = noise["norm"]
    typical = noise["column_norm_median"]
    if noise_norm is None or typical is None:
        noise_verdict = "不能区分"
    elif noise_norm >= typical:
        noise_verdict = "支持"
    elif noise_norm < 0.1 * typical:
        noise_verdict = "排除"
    else:
        noise_verdict = "不能区分"
    char_hit = any(item["h"] is not None or item["at_0"] <= 0.0 for item in switches["char"])
    bound_hit = any(item["h"] is not None or item["at_0"] < 0.0 for item in switches["bounds"])
    smooth = "支持" if char_hit or bound_hit or switches["x_feasible"] is False else "排除"
    reasons = trials["counts"]
    feasible_reject = all(_feasible_reason(key) for key in reasons) and trials["ok"] == 0 and trials["n"] > 0
    if trials["ok"] > 0:
        feasible_verdict = "排除"
    elif feasible_reject:
        feasible_verdict = "支持"
    else:
        feasible_verdict = "不能区分"
    return [
        {"id": "noise", "verdict": noise_verdict, "noise": noise_norm, "typical": typical},
        {"id": "nonsmooth", "verdict": smooth, "char": len(switches["char"]), "bounds": len(switches["bounds"])},
        {"id": "feasible", "verdict": feasible_verdict, "ok": trials["ok"], "n": trials["n"]},
        {"id": "dominant", "verdict": "支持" if share > 0.5 else "排除", "share": share, "name": largest["name"], "block": largest["block"]},
        {"id": "direction", "slope": direction["slope"], "decreased": [row["h"] for row in direction["rows"] if row["merit"] is not None and row["merit"] < direction["merit"]],
         "agreed": [row["h"] for row in direction["rows"] if row["relative"] is not None and row["relative"] < _AGREE]},
    ]


def _markdown(record, body) -> str:
    direction = body["direction"]
    lines = [
        "# T26c 基准失败点诊断",
        "",
        f"运行状态：{record['status']}。{record['reason']}。外层 {record['rounds']} 轮，内层迭代 {record['iterations']} 次。",
        f"失败点文件：`{body['failure_npz']}`。冻结量：`{body['frozen_path']}`。",
        "",
        "床底守恒残差不在内层 F 里。记录中的 bottom 是该轮预计算的 |残差| / ṁ_ref，内层迭代期间不变。",
        "",
        "## 6. 方向",
        "",
        f"失败点归一范数 g(0) = {direction['merit']:.6e}。保存的雅可比给出的线性化斜率 g'(0) = {_num(direction['slope'])}。",
        "",
        "| h | 状态 | g(h) | ‖JΔx − [F(x+hΔx)−F(x)]/h‖ / ‖JΔx‖ | 原因 |",
        "|---|---|---|---|---|",
    ]
    for row in direction["rows"]:
        lines.append(f"| {row['h']:.0e} | {row['status']} | {_num(row['merit'])} | {_num(row['relative'])} | {row['reason']} |")
    decreased = [row["h"] for row in direction["rows"] if row["merit"] is not None and row["merit"] < direction["merit"]]
    agreed = [row["h"] for row in direction["rows"] if row["relative"] is not None and row["relative"] < _AGREE]
    lines.extend([
        "",
        f"g(h) < g(0) 的步长：{_list(decreased)}。相对差小于 {_AGREE:.0e} 的步长：{_list(agreed)}。",
        _direction_sentence(direction, decreased, agreed),
        "",
        "## 7. 内层残差噪声",
        "",
    ])
    noise = body["noise"]
    lines.append(
        f"加密评价状态 {noise['status']}。‖F_tight − F‖ 的归一值 = {_num(noise['norm'])}。"
        f"各列 ‖F(x+s_j e_j)−F(x)‖ 的归一值：最小 {_num(noise['column_norm_min'])}，"
        f"中位 {_num(noise['column_norm_median'])}，最大 {_num(noise['column_norm_max'])}。"
    )
    if noise["reason"]:
        lines.append(f"加密评价的返回说明：{noise['reason']}。")
    lines.extend(["", "同一份冻结量、同一容差下重算 F(x)，与失败时保存的 F 之差的归一范数 = "
                  f"{_num(body['repeated_residual_norm'])}。", "", "## 8. 差分步长", ""])
    lines.append("每一列是与该残差分量同下标的未知量。残差与未知量按同一变量名表对齐，这里不展开该方程的其余依赖列。")
    lines.append("")
    lines.append("| 分量 | 1× 与保存列的相对差 | 1× | 10²× | 10⁴× |")
    lines.append("|---|---|---|---|---|")
    for row in body["columns"]:
        relative = row["relative"]
        lines.append(
            f"| {row['name']} | {_num(row['saved_relative'])} | {_num(relative.get('1.0'))} | "
            f"{_num(relative.get('100.0'))} | {_num(relative.get('10000.0'))} |"
        )
    char = body["switches"]["char"]
    bounds = body["switches"]["bounds"]
    lines.extend([
        "",
        "## 9. 不光滑点",
        "",
        f"x 可行：{body['switches']['x_feasible']}。{body['switches']['x_clause'] or ''}",
        f"x+Δx 可行：{body['switches']['x_plus_dx_feasible']}。{body['switches']['x_plus_dx_clause'] or ''}",
        f"m_C ≤ 0 的位置 {len(char)} 个；可行性边界 {len(bounds)} 个。",
        "",
    ])
    lines.extend(_switch_table("m_C", char))
    lines.extend(_switch_table("可行性", bounds))
    pre = body["precomputed"]
    lines.extend([
        "",
        "预计算量在这一轮内层中不变。悬停存在的 (衡算格, 粒级)：" + _pairs_text(pre["hover"]) + "。",
        "自由状态不存在的 (衡算格, 粒级)：" + _pairs_text(pre["free_absent"]) + "。",
        f"产率式记录：{pre['yield_notes']}。代入温度 {pre['pyrolysis_temperature_k']} K，压力 {pre['pyrolysis_pressure_pa']} Pa。",
        f"Geldart 记录 {pre['geldart_record']}。真实类别 {pre['geldart']}，流体力学调用类别 {pre['hydro_group']}。",
        f"床底守恒残差 {pre['bottom_residual_kg_s']} kg/s。",
        "",
        "## 10. 试探点",
        "",
        f"失败迭代 {body['trials']['iteration']}，归一范数 {body['trials']['merit']:.6e}，试探点 {body['trials']['n']} 个，其中 ok {body['trials']['ok']} 个。",
        "",
    ])
    for reason, count in body["trials"]["counts"].items():
        lines.append(f"- {count}：{reason}")
    lines.extend(["", "## 假说", ""])
    text = {
        "noise": "内层迭代容差带来的残差噪声大于差分步长引起的变化",
        "nonsmooth": "失败点或 x 到 x+Δx 的线段上有 m_C ≤ 0 切换或可行性边界",
        "feasible": "失败迭代的试探点普遍违反可行条件",
        "dominant": "归一范数由一个分量主导（该分量的平方超过范数平方的一半）",
    }
    for item in body["hypotheses"]:
        if item["id"] == "direction":
            lines.append(
                f"- 方向：g'(0) = {_num(item['slope'])}。g(h) < g(0) 的 h 为 {_list(item['decreased'])}。"
                f"相对差 < {_AGREE:.0e} 的 h 为 {_list(item['agreed'])}。"
            )
            continue
        lines.append(f"- {text[item['id']]}：{item['verdict']}。")
    lines.append("")
    lines.append("比较规则：噪声归一值不低于列差分增量的中位数则记为支持，低于该中位数的十分之一则记为排除，介于两者之间记为不能区分。")
    lines.append("")
    return "\n".join(lines)


def _switch_table(title, rows) -> list:
    if not rows:
        return [f"{title}：无。", ""]
    lines = [f"### {title}", "", "| 位置 | h | 起点 | 终点 |", "|---|---|---|---|"]
    for item in rows:
        lines.append(f"| {item['where']} | {_num(item['h'])} | {item['at_0']:.6e} | {item['at_1']:.6e} |")
    lines.append("")
    return lines


def _char_switches(zone, start, end, mask) -> list:
    rows = []
    for cell in range(start.shape[0]):
        for klass in range(1, start.shape[1]):
            if mask is not None and not bool(mask[cell, klass]):
                continue
            left, right = float(start[cell, klass]), float(end[cell, klass])
            root = segment_root(left, right)
            if left <= 0.0 or right <= 0.0 or root not in (None, 0.0):
                rows.append({"where": f"{zone} 衡算格 {cell} k={klass + 1}", "h": None if root in (None, 0.0) and left != 0.0 else root, "at_0": left, "at_1": right})
    return rows


def _sign_switches(zone, start, end, label) -> list:
    rows = []
    array0, array1 = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    for index, (left, right) in enumerate(zip(array0.ravel(), array1.ravel())):
        root = segment_root(float(left), float(right))
        crossed = root not in (None, 0.0)
        if float(left) < 0.0 or float(right) < 0.0 or crossed:
            rows.append({"where": f"{label} {zone} 平坦下标 {index}", "h": root if crossed else None, "at_0": float(left), "at_1": float(right)})
    return rows


def _temperature_switches(state0, state1) -> list:
    rows = []
    for zone, left, right in (
        ("床区", state0.bed_temperature_k, state1.bed_temperature_k),
        ("自由空间", state0.freeboard_temperature_k, state1.freeboard_temperature_k),
    ):
        for index, (start, end) in enumerate(zip(np.asarray(left, dtype=float), np.asarray(right, dtype=float))):
            root = segment_root(float(start), float(end))
            if float(start) <= 0.0 or float(end) <= 0.0 or root not in (None, 0.0):
                rows.append({"where": f"T ≤ 0 {zone} 格 {index}", "h": root, "at_0": float(start), "at_1": float(end)})
    return rows


def _inert_switches(frozen, state0, state1) -> list:
    rows = []
    fuel = np.asarray(frozen.fuel_bed_kg, dtype=float).sum(axis=1)
    mass = np.asarray(frozen.geometric_mass_kg, dtype=float)
    left = mass - fuel - state0.bed_char_kg.sum(axis=1)
    right = mass - fuel - state1.bed_char_kg.sum(axis=1)
    for index, (start, end) in enumerate(zip(left, right)):
        root = segment_root(float(start), float(end))
        if float(start) < 0.0 or float(end) < 0.0 or root not in (None, 0.0):
            rows.append({"where": f"m_inert < 0 床区衡算格 {index}", "h": root if root not in (None, 0.0) else None, "at_0": float(start), "at_1": float(end)})
    return rows


def _endpoint_feasible(values, frozen):
    try:
        ok = _feasible(values, frozen)
    except _CAUGHT as exc:
        return False, f"{type(exc).__name__}: {exc}"
    if ok:
        return True, ""
    return False, _infeasible_clause(values, frozen)


def _safe_unpack(values, frozen):
    from hamel.solver.vector import unpack_state
    try:
        return unpack_state(np.asarray(values, dtype=float), frozen)
    except AssemblyError:
        return None


def _failed_iteration(record) -> dict:
    for item in reversed(record["diagnostics"]["rounds"]):
        newton = item.get("newton") or []
        if newton and item["inner"] == "DAMPING_FAILED":
            return newton[-1]
    raise KeyError("运行记录里没有阻尼失败的迭代")


def _steps(x, scales) -> np.ndarray:
    return np.array([_SQRT_EPS * max(abs(float(x[index])), float(scales[index])) for index in range(x.size)])


def _scale_into(scaled, blocks) -> None:
    from hamel.solver.newton import _divide
    _divide(scaled, blocks.gas, blocks.n_dot_ref)
    _divide(scaled, blocks.solid, blocks.m_dot_ref)
    _divide(scaled, blocks.energy, blocks.q_dot_ref)


def _resolve(path) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = _ROOT / candidate
    return candidate


def _relative(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(_ROOT))
    except ValueError:
        return str(path)


def _pairs(mask) -> list:
    array = np.asarray(mask, dtype=bool)
    return [[int(cell), int(klass)] for cell, klass in zip(*np.nonzero(array))]


def _pairs_text(pairs) -> str:
    if not pairs:
        return "无"
    return "、".join(f"({cell}, k={klass + 1})" for cell, klass in pairs)


def _feasible_reason(reason: str) -> bool:
    return reason.startswith(("Ṅ < 0", "m < 0", "m_inert < 0", "T ≤ 0", "T 非有限", "x 非有限", "feasible 不满足", "可行条件"))


def _direction_sentence(direction, decreased, agreed) -> str:
    slope = direction["slope"]
    merit = direction["merit"]
    if slope is not None and slope < 0 and _num(slope) == _num(-merit):
        linear = "线性化斜率是负的，与 −g(0) 相同到表中的位数，所以 Δx 是线性模型的下降方向。"
    elif slope is not None and slope < 0:
        linear = "线性化斜率为负，所以 Δx 是线性模型的下降方向。"
    else:
        linear = "线性化斜率不是负的，线性模型不把 Δx 标成下降方向。"
    observed = "测得的 g(h) 在全部可行 h 上都没有低于 g(0)。" if not decreased else "有可行 h 使 g(h) 低于 g(0)。"
    match = "差分与 J 在所测的 h 上都没有达到相对差 1e-2。" if not agreed else "有 h 使差分与 J 的相对差小于 1e-2。"
    return linear + observed + match


def _num(value) -> str:
    if value is None:
        return "—"
    return f"{float(value):.6e}"


def _list(values) -> str:
    if not values:
        return "无"
    return "、".join(f"{item:.0e}" for item in values)
