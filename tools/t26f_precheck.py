"""固体子系统可解性预检。固定气体与温度，只解固相行。只报告。

求解手段：ADR 0043。冻结量是第 1 轮预计算。不改残差，也不提出修正。
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.campaign import _process
from hamel.solver.driver import _blocks, _feasible, _references, _scales
from hamel.solver.gates import ResidualBlocks, block_ratios, inner_met
from hamel.solver.newton import SolveTimeout, damped_newton, normalized_merit
from hamel.solver.precalc import precalculate
from hamel.solver.residual import _cell_volume_m3, evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots, unpack_state

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "t26f-precheck"
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_FAILURE = _ROOT / "reports" / "t26c-diagnosis" / "baseline"
_LABELS = ("initial", "failure")
_SCOPE = "子系统收敛只说明在给定状态下找到满足既有门的固相近似根；不收敛或停止不能证明有效域内无正解。"


def _hover_terms(m_c, consumption, migration, residual, generation, migration_in) -> dict:
    """由现有悬停焦炭行的各项作除法。m_C ≤ 0 时不计算比速率。

    ``frozen_rate_estimate_kg`` 把当前比消耗、比迁移速率冻结后作一次除法，只是局部估算。
    现行闭合下比速率随 m_C、m_A 变化，这个数不是原非线性平衡的稳态根，
    不能据此判断实际稳态存量或有效域。
    """
    inflow = float(residual - generation + consumption + migration - migration_in)
    source = inflow + float(generation) + float(migration_in)
    specific = None if m_c <= 0.0 else float(consumption) / float(m_c)
    specific_migration = None if m_c <= 0.0 else float(migration) / float(m_c)
    removal = None if specific is None else specific + specific_migration
    estimate = None if removal in (None, 0.0) else source / removal
    return {
        "inflow_kg_s": inflow, "generation_kg_s": float(generation), "migration_in_kg_s": float(migration_in),
        "specific_consumption_per_s": specific, "specific_migration_per_s": specific_migration,
        "frozen_rate_estimate_kg": estimate, "source_kg_s": source,
    }


def _volume_parts(solid_kg, rho_p, cell_volume, epsilon_mf) -> dict:
    """固体体积与 (1−ε_mf)·V_cell 之比。与有效域使用同一格体积。"""
    limit = float(cell_volume) * (1.0 - float(epsilon_mf))
    solid_volume = float(solid_kg) / float(rho_p)
    ratio = None if limit == 0.0 else solid_volume / limit
    return {"solid_kg": float(solid_kg), "solid_m3": solid_volume, "limit_m3": limit, "ratio": ratio}


def solve(label: str) -> None:
    """解一个给定状态，把完整记录写到 JSON，并向 stdout 打一行摘要。"""
    started = time.perf_counter()
    deadline = time.monotonic() + SolverConfig().wall_timeout_s
    try:
        frozen, full_x = _load(label)
        body = _solve_loaded(label, frozen, full_x, lambda: time.monotonic() >= deadline)
    except SolveTimeout:
        body = _stopped(label, "TIMEOUT", "运行超过墙钟上限")
    body["wall_clock_s"] = time.perf_counter() - started
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / f"{label}.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    kept = {key: body.get(key) for key in ("label", "status", "reason", "iterations", "solid_ratio", "wall_clock_s")}
    sys.stdout.write(json.dumps(kept, ensure_ascii=False))


def _solve_loaded(label, frozen, full_x, deadline) -> dict:
    config = SolverConfig()
    n_ref, m_ref, q_ref = _references(assemble_case(_CASE, "1").physical)
    m_tot = float(frozen.geometric_mass_kg.sum())
    slots = unknown_slots(frozen)
    solid_slots = [slot for slot in slots if slot.kind not in ("gas_outflow", "temperature")]
    solid_index = np.array([slot.index for slot in solid_slots], dtype=int)
    scales = _scales(frozen, n_ref, m_tot)[solid_index]
    empty = np.array([], dtype=int)
    blocks = ResidualBlocks(empty, np.arange(solid_index.size), empty, n_ref, m_ref, q_ref)
    base = np.array(full_x, dtype=float, copy=True)
    calls = {"n": 0}

    def evaluate(solid):
        calls["n"] += 1
        if calls["n"] % 20 == 0:
            print(f"{label} 残差评价 {calls['n']}", file=sys.stderr, flush=True)
            (_OUT / f"{label}.progress").write_text(f"{calls['n']}\n", encoding="utf-8")
        trial = base.copy()
        trial[solid_index] = solid
        return evaluate_residual(trial, frozen).values[solid_index]

    def feasible(solid):
        trial = base.copy()
        trial[solid_index] = solid
        return _feasible(trial, frozen)

    print(f"{label} 固体未知量 {solid_index.size}", file=sys.stderr, flush=True)
    _OUT.mkdir(parents=True, exist_ok=True)
    result = damped_newton(
        base[solid_index], evaluate, scales,
        lambda residual: inner_met(residual, blocks, config.residual_tolerance),
        feasible, lambda residual: normalized_merit(residual, blocks),
        config, deadline, blocks=blocks,
    )
    solved = base.copy()
    solved[solid_index] = result.x
    assembly, assembly_error = _assembly(solved, frozen)
    body = _record(label, frozen, solved, slots, assembly, n_ref, m_ref, q_ref)
    body.update(
        status=result.status, reason=result.reason, iterations=int(result.iterations),
        merit=None if result.merit is None else float(result.merit), residual_calls=calls["n"],
        assembly_error=assembly_error,
    )
    if result.residual is not None and result.residual.size:
        body["subsystem_solid_ratio"] = float(np.max(np.abs(result.residual)) / m_ref)
        order = np.argsort(-np.abs(result.residual))[:8]
        body["largest_solid"] = [{
            "name": solid_slots[int(item)].name, "kg_s": float(result.residual[item]),
        } for item in order]
    return body


def _assembly(full_x, frozen):
    try:
        return evaluate_residual(full_x, frozen), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"


def _record(label, frozen, full_x, slots, assembly, n_ref, m_ref, q_ref) -> dict:
    state = unpack_state(full_x, frozen)
    ratios = None
    if assembly is not None:
        ratios = block_ratios(assembly.values, _blocks(frozen, n_ref, m_ref, q_ref))
    inventories = _inventories(state, frozen)
    positive = [row for row in inventories if row["m_c_kg"] > 0.0 and row["ash_over_char"] is not None]
    volumes = _volumes(state, frozen)
    return {
        "label": label, "scope": _SCOPE, "unknowns": int(full_x.size),
        "block_ratios": ratios, "solid_ratio": None if ratios is None else ratios["solid"],
        "min_char": _minimum(inventories), "min_positive_char": _minimum(positive),
        "inventories": inventories, "volumes": volumes,
        "volume_over_limit": [row for row in volumes if row["ratio"] is not None and row["ratio"] > 1.0],
        "interface": _interface(state, frozen), "hover": _hover(state, frozen, slots, assembly),
    }


def _inventories(state, frozen) -> list:
    rows = []
    inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    for cell in range(frozen.n_bed):
        total = float(state.bed_char_kg[cell].sum())
        for klass in range(frozen.n_class):
            m_c = float(state.bed_char_kg[cell, klass])
            m_a = None if total <= 0.0 else float(inert[cell]) * m_c / total
            rows.append(_pair("bed", cell, frozen.n_bed, klass, "bed", m_c, m_a, "床区灰按该格焦炭质量比例分到粒级"))
    for cell in range(frozen.n_freeboard):
        for klass in range(frozen.n_class):
            if not frozen.free_absent[cell, klass]:
                rows.append(_pair(
                    "freeboard", cell, frozen.n_freeboard, klass, "free",
                    float(state.freeboard_char_free_kg[cell, klass]), float(state.freeboard_ash_free_kg[cell, klass]),
                    "该状态的灰未知量",
                ))
            if frozen.hover[cell, klass]:
                rows.append(_pair(
                    "freeboard", cell, frozen.n_freeboard, klass, "hover",
                    float(state.freeboard_char_hover_kg[cell, klass]), float(state.freeboard_ash_hover_kg[cell, klass]),
                    "该状态的灰未知量",
                ))
    return rows


def _pair(zone, cell, count, klass, phase, m_c, m_a, ash_rule) -> dict:
    ratio = None if m_a is None or m_c <= 0.0 else float(m_a) / float(m_c)
    return {
        "zone": zone, "balance_index": int(cell), "geometric_index": int(count - 1 - cell),
        "k": int(klass + 1), "state": phase, "m_c_kg": float(m_c),
        "m_a_kg": None if m_a is None else float(m_a), "ash_over_char": ratio, "ash_rule": ash_rule,
    }


def _minimum(rows):
    if not rows:
        return None
    lowest = min(row["m_c_kg"] for row in rows)
    return {"m_c_kg": lowest, "rows": [row for row in rows if row["m_c_kg"] == lowest]}


def _volumes(state, frozen) -> list:
    rows = []
    for cell in range(frozen.n_freeboard):
        solid = float(
            state.freeboard_char_free_kg[cell].sum() + state.freeboard_char_hover_kg[cell].sum()
            + state.freeboard_ash_free_kg[cell].sum() + state.freeboard_ash_hover_kg[cell].sum()
        )
        parts = _volume_parts(solid, frozen.rho_p_kg_m3, _cell_volume_m3("freeboard", cell, frozen), frozen.epsilon_mf)
        parts.update(balance_index=int(cell), geometric_index=int(frozen.n_freeboard - 1 - cell))
        rows.append(parts)
    return rows


def _interface(state, frozen) -> dict:
    """床区各界面与床面的上、下行流。衡算下标 0 是床顶格。

    格 i 与其下方格 i + 1 之间：上行取下方格的 Ṁ_auf[i+1] = K_auf[i+1]·M[i+1]，
    下行取上方格的 Ṁ_ab[i]（ADR 0009 第 2 项；与 ``solid_mass_flow_residual`` 的入流下标一致）。
    床面：上行取床顶格的 Ṁ_auf[0]（床面抛出），下行取自由空间最下格的回落。
    另列第 1 轮冻结的 Ṁ_ab,0。焦炭按各源格的焦炭质量比例分。
    """
    mass = np.asarray(frozen.geometric_mass_kg, dtype=float)
    k_auf = np.asarray(frozen.k_auf_per_s, dtype=float)
    m_ab = np.asarray(frozen.m_dot_ab_kg_s, dtype=float)
    char = np.asarray(state.bed_char_kg, dtype=float)
    bottom = frozen.n_freeboard - 1
    returning = (frozen.b_free_per_s[bottom] + frozen.b_hover_per_s[bottom]) * state.freeboard_char_free_kg[bottom]
    rows = []
    for upper in range(frozen.n_bed - 1):
        lower = upper + 1
        up = float(k_auf[lower] * mass[lower])
        down = float(m_ab[upper])
        char_up = k_auf[lower] * char[lower]
        char_down = None if mass[upper] == 0.0 else m_ab[upper] * char[upper] / mass[upper]
        rows.append({
            "face": f"床区 {upper}/{lower}", "upper_balance_index": upper, "lower_balance_index": lower,
            "up_kg_s": up, "down_kg_s": down, "net_up_kg_s": up - down,
            "char_up_kg_s": [float(item) for item in char_up],
            "char_down_kg_s": None if char_down is None else [float(item) for item in char_down],
        })
    surface_up = float(k_auf[0] * mass[0])
    char_surface_up = k_auf[0] * char[0]
    rows.append({
        "face": "床面", "upper_balance_index": None, "lower_balance_index": 0,
        "up_kg_s": surface_up, "down_kg_s": None, "net_up_kg_s": None,
        "char_up_kg_s": [float(item) for item in char_surface_up],
        "char_down_kg_s": [float(item) for item in returning],
        "char_net_up_kg_s": float(char_surface_up.sum() - returning.sum()),
    })
    classes = []
    for klass in range(frozen.n_class):
        up = float(char_surface_up[klass])
        classes.append({
            "k": int(klass + 1), "char_up_kg_s": up, "return_char_kg_s": float(returning[klass]),
            "char_net_up_kg_s": up - float(returning[klass]),
            "ejection_to_free_kg_s": up * float(frozen.ejection_to_free[klass]),
            "ejection_to_hover_kg_s": up * float(frozen.ejection_to_hover[klass]),
        })
    return {
        "m_dot_ab_0_kg_s": float(m_ab[0]), "m_dot_auf_0_kg_s": surface_up,
        "faces": rows, "surface_by_class": classes,
        "ejection_to_free_kg_s": float(sum(row["ejection_to_free_kg_s"] for row in classes)),
        "ejection_to_hover_kg_s": float(sum(row["ejection_to_hover_kg_s"] for row in classes)),
        "return_char_kg_s": float(returning.sum()),
    }


def _hover(state, frozen, slots, assembly) -> list:
    if assembly is None:
        return []
    residual_of = {slot.index: float(assembly.values[slot.index]) for slot in slots if slot.kind == "char" and slot.phase == "hover"}
    slot_of = {(slot.balance_index, slot.size_class): slot for slot in slots if slot.kind == "char" and slot.phase == "hover"}
    rows = []
    for cell in range(frozen.n_freeboard):
        for klass in range(frozen.n_class):
            if not frozen.hover[cell, klass]:
                continue
            slot = slot_of[(cell, klass)]
            m_c = float(state.freeboard_char_hover_kg[cell, klass])
            migration_in = 0.0
            if klass + 1 < frozen.n_class:
                migration_in = float(frozen.theta[cell, klass]) * float(
                    assembly.char_migration_free_kg_s[cell, klass + 1] + assembly.char_migration_hover_kg_s[cell, klass + 1]
                )
            terms = _hover_terms(
                m_c, float(assembly.char_consumption_hover_kg_s[cell, klass]),
                float(assembly.char_migration_hover_kg_s[cell, klass]), residual_of[slot.index],
                float(frozen.char_gen_hover_kg_s[cell, klass]), migration_in,
            )
            terms.update(
                balance_index=int(cell), geometric_index=int(frozen.n_freeboard - 1 - cell),
                k=int(klass + 1), name=slot.name, m_c_kg=m_c,
                m_a_kg=float(state.freeboard_ash_hover_kg[cell, klass]),
                consumption_kg_s=float(assembly.char_consumption_hover_kg_s[cell, klass]),
                migration_kg_s=float(assembly.char_migration_hover_kg_s[cell, klass]),
                residual_kg_s=residual_of[slot.index],
            )
            rows.append(terms)
    return rows


def _load(label: str):
    if label == "initial":
        print("初值点第 1 轮预计算", file=sys.stderr, flush=True)
        _OUT.mkdir(parents=True, exist_ok=True)
        (_OUT / "initial.progress").write_text("precalc\n", encoding="utf-8")
        prepared = precalculate(_CASE, "1", None, adiabatic=False, iteration=1, solver=SolverConfig())
        print("预计算结束", file=sys.stderr, flush=True)
        return prepared.frozen, pack_state(prepared.state, prepared.frozen)
    if label != "failure":
        raise SystemExit("状态只能是 initial 或 failure")
    with np.load(_FAILURE / "failure.npz", allow_pickle=True) as saved:
        payload = {key: np.array(saved[key]) for key in saved.files}
    with (_FAILURE / "frozen.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    names = [slot.name for slot in unknown_slots(frozen)]
    if names != [str(item) for item in payload["names"].tolist()]:
        raise SystemExit("失败点变量名与冻结量不一致")
    return frozen, np.asarray(payload["x"], dtype=float)


def _stopped(label, status, reason) -> dict:
    return {"label": label, "status": status, "reason": reason, "scope": _SCOPE, "iterations": None}


def supervise(label: str, timeout: float) -> dict:
    """子进程加硬超时。超时记 TIMEOUT，不换策略。"""
    folder = _OUT / label
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed = _process([sys.executable, "-m", "tools.t26f_precheck", "solve", label], timeout)
    if completed is None:
        body = _stopped(label, "TIMEOUT", "墙钟上限内没有返回")
    elif completed.returncode != 0:
        (folder / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
        body = _stopped(label, "ERROR", f"子进程退出码 {completed.returncode}")
        body["stderr_log"] = str((folder / "stderr.log").relative_to(_OUT))
    else:
        body = json.loads((_OUT / f"{label}.json").read_text(encoding="utf-8"))
        if completed.stderr:
            (folder / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    body["supervised_wall_clock_s"] = time.perf_counter() - started
    (folder / "record.json").write_text(json.dumps(body, ensure_ascii=False, indent=2), encoding="utf-8")
    return body


def write_report(records: dict) -> None:
    lines = [
        "# T26f 固体子系统预检", "", _SCOPE, "",
        "气体流量、组成与各格温度固定为给定状态。未知量只含固相存量：床区焦炭，自由空间的自由、悬停焦炭与灰。",
        "残差取 `evaluate_residual` 的固相行。阻尼牛顿、差分步长与行列缩放用 ADR 0043。收敛门是 ADR 0012 第 1 项的固相块。",
        "冻结量是第 1 轮预计算。初值点按 ADR 0017 第 9 项现算；失败点用 `reports/t26c-diagnosis/baseline/` 保存的状态和冻结量。",
        "",
    ]
    for label, title in (("initial", "工况 1 初值点"), ("failure", "T26c 基准失败点")):
        lines += [f"## {title}", ""]
        lines.extend(_state_lines(records.get(label)))
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "precheck.md").write_text("\n".join(lines).rstrip("\n") + "\n", encoding="utf-8")
    (_OUT / "precheck.json").write_text(json.dumps({"scope": _SCOPE, "states": records}, ensure_ascii=False, indent=2), encoding="utf-8")


def _state_lines(body) -> list:
    if not body:
        return ["没有记录。", ""]
    lines = [
        f"状态 {body.get('status')}。原因：{body.get('reason') or '无'}。迭代 {body.get('iterations')}。墙钟 {_fmt(body.get('supervised_wall_clock_s', body.get('wall_clock_s')))} s。",
        f"固相 ‖r‖∞/ṁ_ref = {_fmt(body.get('solid_ratio'))}。分块：{_fmt(body.get('block_ratios'))}。",
    ]
    if body.get("assembly_error"):
        lines.append(f"终态残差装配停止：{body['assembly_error']}。存量与体积仍按该向量计算；比消耗速率、入流与稳态存量空缺。")
    minimum = body.get("min_char")
    if minimum:
        lines += ["", f"最小 m_C = {_fmt(minimum['m_c_kg'])} kg。", "", _inventory_table(minimum["rows"])]
    positive = body.get("min_positive_char")
    if positive and (minimum is None or positive["m_c_kg"] != minimum["m_c_kg"]):
        lines += ["", f"最小正 m_C = {_fmt(positive['m_c_kg'])} kg。", "", _inventory_table(positive["rows"])]
    volumes = body.get("volumes") or []
    over = body.get("volume_over_limit") or []
    lines += ["", f"自由与悬停固体体积之和相对 (1−ε_mf)·V_cell。超过 1 的格子：{len(over)}。逐项存量在同目录 JSON 的 inventories。"]
    if volumes:
        lines += ["", "| 衡算下标 | 几何下标 | 固体 kg | 固体 m³ | 上限 m³ | 比值 |", "|---|---|---|---|---|---|"]
        lines += [f"| {row['balance_index']} | {row['geometric_index']} | {_fmt(row['solid_kg'])} | {_fmt(row['solid_m3'])} | {_fmt(row['limit_m3'])} | {_fmt(row['ratio'])} |" for row in volumes]
    face = body.get("interface")
    if face:
        lines += [
            "", "床区界面的冻结整体流。格 i/i+1 之间上行取下方格的 Ṁ_auf[i+1]，下行取上方格的 Ṁ_ab[i]；衡算下标 0 是床顶格。",
            f"第 1 轮冻结的 Ṁ_ab,0 = {_fmt(face['m_dot_ab_0_kg_s'])} kg/s（床顶格向下进入格 1 的流）。",
            "", "| 界面 | 上行 kg/s | 下行 kg/s | 净上行 kg/s |", "|---|---|---|---|",
        ]
        lines += [
            f"| {row['face']} | {_fmt(row['up_kg_s'])} | {_fmt(row['down_kg_s'])} | {_fmt(row['net_up_kg_s'])} |"
            for row in face["faces"]
        ]
        lines += [
            "", "床面：上行为床顶格的抛出 Ṁ_auf,0，下行为自由空间最下格的焦炭回落；床面不以 Ṁ_ab,0 代替回落。",
            f"焦炭抛向自由状态 {_fmt(face['ejection_to_free_kg_s'])} kg/s，抛向悬停状态 {_fmt(face['ejection_to_hover_kg_s'])} kg/s，回落 {_fmt(face['return_char_kg_s'])} kg/s。",
            "", "| k | 床面焦炭上行 | 回落 | 净上行 | 抛向自由 | 抛向悬停 |", "|---|---|---|---|---|---|",
        ]
        lines += [
            f"| {row['k']} | {_fmt(row['char_up_kg_s'])} | {_fmt(row['return_char_kg_s'])} | {_fmt(row['char_net_up_kg_s'])} | {_fmt(row['ejection_to_free_kg_s'])} | {_fmt(row['ejection_to_hover_kg_s'])} |"
            for row in face["surface_by_class"]
        ]
    hover = body.get("hover") or []
    focus = [row for row in hover if row["k"] == 3]
    lines += ["", f"悬停状态 {len(hover)} 个，其中 k = 3 的 {len(focus)} 个。比消耗速率 = 消耗 / m_C。冻结比速率估算 = (入流 + 冻结生成 + θ·更大粒级迁入) / (比消耗速率 + 比迁移速率)，m_C ≤ 0 或分母为 0 时不计算。",
              "这只是把当前比速率冻结后的局部估算：现行闭合下比速率随 m_C、m_A 变化，该估算不满足原非线性平衡，不是稳态存量，不能据此判断实际稳态存量或有效域。"]
    if focus:
        lines += ["", _hover_table(focus)]
    if hover and len(hover) != len(focus):
        lines.append("其余悬停状态在 JSON 的 hover。")
    largest = body.get("largest_solid") or []
    if largest and body.get("status") != "CONVERGED":
        lines += ["", "固相残差绝对值最大的几项：", "", "| 名称 | kg/s |", "|---|---|"]
        lines += [f"| {row['name']} | {_fmt(row['kg_s'])} |" for row in largest]
    lines += ["", _SCOPE, ""]
    return lines


def _inventory_table(rows) -> str:
    lines = ["| 区 | 衡算下标 | k | 状态 | m_C kg | m_A kg | m_A/m_C |", "|---|---|---|---|---|---|---|"]
    lines += [
        f"| {row['zone']} | {row['balance_index']} | {row['k']} | {row['state']} | {_fmt(row['m_c_kg'])} | {_fmt(row['m_a_kg'])} | {_fmt(row['ash_over_char'])} |"
        for row in rows
    ]
    return "\n".join(lines)


def _hover_table(rows) -> str:
    lines = ["| 衡算下标 | k | m_C kg | 比消耗 1/s | 入流 kg/s | 冻结比速率估算 kg | 残差 kg/s |", "|---|---|---|---|---|---|---|"]
    lines += [
        f"| {row['balance_index']} | {row['k']} | {_fmt(row['m_c_kg'])} | {_fmt(row['specific_consumption_per_s'])} | {_fmt(row['inflow_kg_s'])} | {_fmt(row['frozen_rate_estimate_kg'])} | {_fmt(row['residual_kg_s'])} |"
        for row in rows
    ]
    return "\n".join(lines)


def _fmt(value) -> str:
    if isinstance(value, dict):
        return "{" + ", ".join(f"{key}={_fmt(item)}" for key, item in value.items()) + "}"
    if value is None:
        return "—"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return "—"
    return f"{number:.8e}"


def main(argv) -> None:
    command = argv[1] if len(argv) > 1 else "run"
    if command == "solve":
        solve(argv[2])
        return
    if command != "run":
        raise SystemExit("用法：run | solve initial | solve failure")
    timeout = float(SolverConfig().wall_timeout_s)
    records = {label: supervise(label, timeout) for label in _LABELS}
    write_report(records)


if __name__ == "__main__":
    main(sys.argv)
