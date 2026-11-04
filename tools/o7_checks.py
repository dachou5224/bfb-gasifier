"""O7（灰的闭合嫌疑）的低成本检查 K1–K4（M5 审计；目标所有者 2026-11-04 决定）。只报告，不改 hamel/。

输入：`tools/solver_path_diagnostic.py outer` 写出的各轮解与冻结量（只在本地）：
    <dir>/outer-x-round{k}.npy、<dir>/outer-frozen-round{k}.pkl；初值缓存 <cache>/x0.npy。
只用质量流与存量，不读取、不比较任何出口验收量。

- k4：每轮（第 k 轮冻结量 + 第 k 轮内层解）的床区焦炭与粒级份额、床灰、自由空间顶格各粒级灰、
  带出灰 E = Σ a^f_顶·m_灰,自由,顶、灰生成 G = Σ ash_gen（床区 + 自由 + 悬停）、余额 D = G − E。
  A5 的排灰灰分量 = 燃料带入灰 − 带出灰 − 随未转化燃料离开的灰；T26n 后后两项之和恒等于
  燃料带入灰 − Σ ash_gen，故同一状态上 A5 余额 = G − E（kg/s）。
- k3：同一状态上床区灰的边界收支差 = Σ ash_gen_bed + 回落灰 − 床面抛出灰 − 冻结的排灰灰分量。
- k1：固定第 k 轮的冻结量、气相与温度，把床区焦炭按比例缩放到给定总量 C（床灰 = M − C − 燃料，
  按 ADR 0017 第 19 项分粒级），只解自由空间固相（焦炭与灰，自由与悬停），算 D(C)，求 D = 0 的 C*。
- ledger：逐粒级的床灰 A_k、带出灰、带出系数 e_k = 带出灰 / A_k 与热解释放到床区的灰 S_k（闭合复核用）。
- k2：固定第 k 轮的冻结量，从不同起点跑内层（gaslog,o2pow,relax），比较收敛的床区焦炭。
  --drain-zero：把冻结的排灰与排焦炭取 0（反事实设定，不是模型写法）。在耦合问题中它不是
  床区焦炭的严格上界；解出的状态只看它自己的 A5 余额，不与 k1 的 C*（只适用于 k1 的扫描族）比较。
--corellaX：第 1 轮状态的最底床格高于 900 °C，评价残差须同诊断一样外推 Corella（仅诊断）。
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.newton as newton
from hamel.solver.gates import ResidualBlocks, inner_met
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state, unknown_slots, unpack_state
from tools.solver_path_diagnostic import Problem, run_inner

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "research-solver-path" / "o7"


def _load(folder: Path, k: int):
    with (folder / f"outer-frozen-round{k}.pkl").open("rb") as handle:
        frozen = pickle.load(handle)
    return frozen, np.load(folder / f"outer-x-round{k}.npy")


def balance(frozen, x) -> dict:
    """K4 与 K3：同一状态上的灰收支（kg/s）与存量（kg）。"""
    state = unpack_state(x, frozen)
    char = np.asarray(state.bed_char_kg, dtype=float)
    masses = np.asarray(frozen.geometric_mass_kg, dtype=float)
    inert = masses - char.sum(axis=1) - np.asarray(frozen.fuel_bed_kg).sum(axis=1)
    share = np.where(char.sum(axis=1, keepdims=True) > 0, char / char.sum(axis=1, keepdims=True), 0.0)
    bed_ash_by_class = inert[:, None] * share
    top = np.asarray(state.freeboard_ash_free_kg[0], dtype=float)
    carried = float(np.dot(frozen.a_free_per_s[0], top))
    generated = float(np.sum(frozen.ash_gen_bed_kg_s) + np.sum(frozen.ash_gen_free_kg_s) + np.sum(frozen.ash_gen_hover_kg_s))
    ejected = np.asarray(frozen.k_auf_per_s[0] * bed_ash_by_class[0], dtype=float)
    returning = float(np.dot(frozen.b_free_per_s[-1], state.freeboard_ash_free_kg[-1]))
    return {
        "bed_char_kg": float(char.sum()), "bed_char_class_share": (char.sum(axis=0) / char.sum()).tolist(),
        "bed_ash_kg": float(inert.sum()), "bed_ash_by_cell_kg": inert.tolist(),
        "geometric_mass_kg": float(masses.sum()), "fuel_bed_kg": float(np.sum(frozen.fuel_bed_kg)),
        "freeboard_top_free_ash_kg": top.tolist(), "a_free_top_per_s": np.asarray(frozen.a_free_per_s[0]).tolist(),
        "carried_ash_kg_s": carried, "ash_generated_kg_s": generated, "a5_balance_kg_s": generated - carried,
        "frozen_drain_ash_kg_s": float(frozen.drain_ash_kg_s),
        "bed_ash_boundary_gap_kg_s": float(np.sum(frozen.ash_gen_bed_kg_s)) + returning - float(ejected.sum()) - float(frozen.drain_ash_kg_s),
        "ejected_ash_kg_s": float(ejected.sum()), "returning_ash_kg_s": returning,
        "bed_temperature_c": (np.asarray(state.bed_temperature_k) - 273.15).tolist(),
    }


def ledger(frozen, x) -> dict:
    """逐粒级灰账：床灰按现行闭合（随焦炭比例）分粒级时的 A_k、带出灰与 e_k，以及床区灰供给 S_k。"""
    state = unpack_state(x, frozen)
    char = np.asarray(state.bed_char_kg, dtype=float)
    inert = np.asarray(frozen.geometric_mass_kg) - char.sum(axis=1) - np.asarray(frozen.fuel_bed_kg).sum(axis=1)
    by_class = (inert[:, None] * char / char.sum(axis=1, keepdims=True)).sum(axis=0)
    carried = np.asarray(frozen.a_free_per_s[0]) * np.asarray(state.freeboard_ash_free_kg[0])
    with np.errstate(invalid="ignore", divide="ignore"):
        e = np.where(by_class > 0, carried / by_class, 0.0)
    return {
        "diameter_mm": (np.asarray(frozen.diameter_m) * 1e3).tolist(),
        "bed_ash_by_class_kg": by_class.tolist(), "carried_by_class_kg_s": carried.tolist(),
        "elutriation_per_s": e.tolist(),
        "ash_supply_bed_by_class_kg_s": np.asarray(frozen.ash_gen_bed_kg_s).sum(axis=0).tolist(),
        "ash_supply_freeboard_kg_s": float(np.sum(frozen.ash_gen_free_kg_s) + np.sum(frozen.ash_gen_hover_kg_s)),
    }


def response(frozen, x, totals) -> list:
    """K1：床区焦炭缩放到各总量 C，只解自由空间固相，返回 D(C)。"""
    problem = Problem(frozen)
    slots = unknown_slots(frozen)
    rows = np.array([s.index for s in slots if s.kind in ("char", "ash") and s.name.startswith("freeboard")], dtype=int)
    blocks = ResidualBlocks(np.array([], dtype=int), np.arange(rows.size), np.array([], dtype=int),
                            problem.n_ref, problem.m_ref, problem.q_ref)
    scales = np.full(rows.size, problem.m_tot)
    out = []
    base_state = unpack_state(x, frozen)
    base_total = float(base_state.bed_char_kg.sum())
    for total in totals:
        state = unpack_state(x, frozen)
        state.bed_char_kg = np.asarray(base_state.bed_char_kg, dtype=float) * (total / base_total)
        start = np.asarray(pack_state(state, frozen), dtype=float)

        def embed(y, start=start):
            z = start.copy()
            z[rows] = y
            return z

        result = newton.damped_newton(
            start[rows], lambda y: evaluate_residual(embed(y), frozen).values[rows], scales,
            lambda r: inner_met(r, blocks, problem.config.residual_tolerance),
            lambda y: bool(np.isfinite(y).all()), lambda r: newton.normalized_merit(r, blocks),
            problem.config, None, None, blocks, None, False,
        )
        entry = {"bed_char_kg": total, "status": result.status, "iterations": result.iterations,
                 "solid_merit": result.merit}
        if result.status == "CONVERGED":
            entry.update(balance(frozen, embed(result.x)))
            entry["min_freeboard_stock_kg"] = float(np.min(result.x))
        out.append(entry)
    return out


def zero_crossing(points) -> float | None:
    """D(C) 的线性插值零点；没有变号时返回 None。"""
    good = [(p["bed_char_kg"], p["a5_balance_kg_s"]) for p in points if p.get("status") == "CONVERGED"]
    for (c0, d0), (c1, d1) in zip(good, good[1:]):
        if d0 == 0.0:
            return c0
        if d0 * d1 < 0.0:
            return c0 + (c1 - c0) * d0 / (d0 - d1)
    return None


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("check", choices=("k4", "k1", "k2", "ledger"))
    parser.add_argument("--dir", required=True, help="outer 写出的目录")
    parser.add_argument("--cache", default=None, help="build 的缓存目录（k2 的初值起点）")
    parser.add_argument("--round", type=int, default=2)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--totals", default="0.5,1,1.5,2,2.5,3,3.5")
    parser.add_argument("--limit-s", type=float, default=5400.0)
    parser.add_argument("--out", default=str(_OUT))
    parser.add_argument("--drain-zero", action="store_true")
    parser.add_argument("--start-char", default="", help="k2 的附加起点：床区焦炭总量（kg），逗号分隔")
    parser.add_argument("--corellaX", action="store_true")
    args = parser.parse_args(argv)
    if args.corellaX:
        import hamel.kinetics.homogeneous as homogeneous
        homogeneous._CORELLA_HIGH_K = float("inf")
    folder, out = Path(args.dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    if args.check == "ledger":
        summary = {"rounds": [dict(round=k, **ledger(*_load(folder, k))) for k in range(1, args.rounds + 1)]}
    elif args.check == "k4":
        summary = {"rounds": [dict(round=k, **balance(*_load(folder, k))) for k in range(1, args.rounds + 1)]}
    elif args.check == "k1":
        frozen, x = _load(folder, args.round)
        points = response(frozen, x, [float(v) for v in args.totals.split(",")])
        summary = {"round": args.round, "points": points, "c_star_kg": zero_crossing(points)}
    else:
        frozen, x = _load(folder, args.round)
        if args.drain_zero:
            frozen = replace(frozen, drain_ash_kg_s=0.0, drain_char_kg_s=np.zeros_like(frozen.drain_char_kg_s))
        problem = Problem(frozen)
        starts = {"round_solution": x}
        size = x.size
        if args.cache and np.load(Path(args.cache) / "x0.npy").size == size:
            starts["initial_x0"] = np.load(Path(args.cache) / "x0.npy")
        # 同一轮解上把床区焦炭按比例缩放到给定总量作起点（不同轮的未知量个数可能不同）。
        base = unpack_state(x, frozen)
        for total in [float(v) for v in args.start_char.split(",") if v]:
            state = unpack_state(x, frozen)
            state.bed_char_kg = np.asarray(base.bed_char_kg) * (total / float(base.bed_char_kg.sum()))
            starts[f"char_{total:g}kg"] = np.asarray(pack_state(state, frozen), dtype=float)
        summary = {"round": args.round, "drain_zero": args.drain_zero, "corellaX": args.corellaX, "runs": {}}
        for name, start in starts.items():
            local: dict = {}
            run = run_inner(problem, start, {"gaslog", "o2pow", "relax"}, args.limit_s, 150, local)
            stage = run["stages"][-1]
            entry = {"status": stage["status"], "reason": stage["reason"], "iterations": len(stage["iterations"]),
                     "final_merit": stage["final_merit"], "final_check_ok": (stage.get("final_check") or {}).get("ok"),
                     "start": balance(frozen, start)["bed_char_kg"]}
            if "x_final" in local:
                entry.update(balance(frozen, local["x_final"]))
            summary["runs"][name] = entry
    summary["elapsed_s"] = time.monotonic() - started
    suffix = ("-drain0" if args.check == "k2" and args.drain_zero else "")
    name = f"{args.check}.json" if args.check in ("k4", "ledger") else f"{args.check}-round{args.round}{suffix}.json"
    summary["corellaX"] = args.corellaX
    (out / name).write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    print(name, round(summary["elapsed_s"]))


if __name__ == "__main__":
    main()
