"""限时求解路径诊断（2026-11）。只报告，不改 ``hamel/``。

求解手段（气体对数变量、O₂ 幂变量、三层试探规则）已并入正式实现（ADR 0054），本工具不再
自带一套牛顿：``inner`` 与 ``outer`` 直接调用 ``driver.solve_aligned`` 的同一组
``newton.build_transform``／``newton.build_rule``。只保留诊断专用的两项：

- ``corellaX``：**仅诊断、属模型有效域放宽**：Corella 高温段外推到 900 °C 以上（ADR 0010 第 9 条
  未批准外推，不进入 ``hamel/``）。
- ``o2mix``：只改起点，把进风格 O₂ 出流取进风 O₂ 的 ``--o2-share``。

子命令
- build：跑预计算与初值构造，把交给内层的起点 x0 与冻结量存进 --cache。
- inner：读 --cache，跑第 1 轮内层，逐次记录牛顿迭代、阻尼试探与障碍普查。
- outer：完整 solve_case（外层 --rounds 轮），逐轮记录状态与收敛解检查。
- check：对已保存的 x 做第 3 层检查。

所有结果写在 --out（默认 reports/research-solver-path/diagnostic/）。既有运行记录不重生成。
"""

from __future__ import annotations

import argparse
import json
import pickle
import signal
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
import hamel.solver.newton as newton
import hamel.solver.precalc as precalc
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.properties.data import GAS_SPECIES
from hamel.solver.gates import block_ratios, inner_met
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import jacobian_layout, pack_state, unknown_slots, unpack_state

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-solver-path" / "diagnostic"
_O2 = GAS_SPECIES.index("O2")


class _Timeout(BaseException):
    """继承 BaseException：任何 except Exception 都不会吞掉全局期限（CHARTER §7）。"""


class _Built(BaseException):
    """初值构造已结束：跳出 solve_case。"""


def _corella_extrapolate() -> None:
    """仅诊断：Corella 高温段参数外推到 900 °C 以上。ADR 0010 第 9 条未批准。"""
    import hamel.kinetics.homogeneous as homogeneous
    homogeneous._CORELLA_HIGH_K = float("inf")


def _alarm(limit_s):
    def expire(signum, frame):
        raise _Timeout()

    previous = signal.signal(signal.SIGALRM, expire)
    signal.alarm(int(limit_s) + 1)

    def restore():
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    return restore


def build(cache: Path, loops: int, limit_s: float) -> dict:
    """预计算 + 初值构造，存 x0 与冻结量。"""
    captured: dict = {}
    original = driver.construct

    def construct(*args, **kwargs):
        kwargs["max_loops"] = loops
        out = original(*args, **kwargs)
        captured.update(result=out[0], loops=out[1], note=out[2])
        raise _Built()

    driver.construct = construct
    restore_alarm = _alarm(limit_s)
    started = time.monotonic()
    status = ["BUILT", ""]
    try:
        driver.solve_case(str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=1))
        status = ["NOT_BUILT", "solve_case 未进入初值构造"]
    except _Built:
        pass
    except _Timeout:
        status = ["TIMEOUT", ""]
    except Exception as exc:  # noqa: BLE001 - 诊断要记录任何停止
        status = [type(exc).__name__, str(exc)[:800]]
    finally:
        restore_alarm()
        driver.construct = original
    summary = {"status": status, "elapsed_s": time.monotonic() - started}
    if "result" in captured:
        result = captured["result"]
        cache.mkdir(parents=True, exist_ok=True)
        with (cache / "frozen.pkl").open("wb") as handle:
            pickle.dump(result.frozen, handle, protocol=pickle.HIGHEST_PROTOCOL)
        np.save(cache / "x0.npy", pack_state(result.state, result.frozen))
        summary["loops"] = captured["loops"]
        summary["note"] = captured["note"]
        summary["drain_ash_kg_s"] = float(result.frozen.drain_ash_kg_s)
    return summary


class Problem:
    """第 1 轮内层的原问题：与 driver.solve_aligned 的 prepare 同一组块、尺度与规则。"""

    def __init__(self, frozen, config=None):
        assembled = assemble_case(str(_CASE), "1")
        self.config = assembled.solver if config is None else config
        self.frozen = frozen
        self.n_ref, self.m_ref, self.q_ref = driver._references(assembled.physical)
        self.blocks = driver._blocks(frozen, self.n_ref, self.m_ref, self.q_ref)
        self.m_tot = float(frozen.geometric_mass_kg.sum())
        self.slots = unknown_slots(frozen)
        self.names = [slot.name for slot in self.slots]
        self.layout = jacobian_layout(self.slots, frozen)
        self.scales = driver._scales(frozen, self.n_ref, self.m_tot)
        self.o2 = np.array([i for i, name in enumerate(self.names) if name.endswith(".O2")], dtype=int)
        self.evaluations = 0

    def transform(self, x0):
        return newton.build_transform(self.slots, x0, self.n_ref)

    def rule(self, transform):
        """正式实现的三层规则。求解手段：ADR 0054 第 2 项。"""
        return newton.build_rule(
            self.slots, self.frozen, transform,
            tau=float(self.config.stock_roundoff_relative) * self.m_tot,
            tau_n=float(self.config.stock_roundoff_relative) * self.n_ref,
        )

    def residual(self, x):
        self.evaluations += 1
        return evaluate_residual(x, self.frozen).values

    def merit(self, values):
        return newton.normalized_merit(values, self.blocks)

    def converged(self, values):
        return inner_met(values, self.blocks, self.config.residual_tolerance)

    def ratios(self, values):
        return block_ratios(values, self.blocks)


def o2_mix(problem: Problem, x0, share: float):
    """进风格 O₂ 出流取进风 O₂ × share；其余不变。只改起点。"""
    frozen = problem.frozen
    x = np.array(x0, dtype=float)
    state = unpack_state(x, frozen)
    primary = float(np.asarray(frozen.primary_feed_mol_s)[_O2])
    secondary = float(np.asarray(frozen.secondary_feed_mol_s)[_O2])
    changed = {}
    bed = np.asarray(state.bed_gas_mol_s, dtype=float)
    total = bed[-1].sum(axis=1)
    for phase in range(2):
        target = share * primary * total[phase] / total.sum()
        if bed[-1, phase, _O2] < target:
            changed[f"bed[-1,{phase}]"] = (float(bed[-1, phase, _O2]), target)
            bed[-1, phase, _O2] = target
    state.bed_gas_mol_s = bed
    free = np.asarray(state.freeboard_gas_mol_s, dtype=float)
    cell = int(frozen.secondary_balance_index)
    target = share * secondary
    if free[cell, _O2] < target:
        changed[f"free[{cell}]"] = (float(free[cell, _O2]), target)
        free[cell, _O2] = target
    state.freeboard_gas_mol_s = free
    return pack_state(state, frozen), changed


def projected_damp(y, delta, merit, evaluate, feasible, merit_of, limit, deadline, trials=None, explain=None, transform=newton.IDENTITY):
    """o2proj（仅诊断）：试探点上把 O₂ 幂变量 s 投影到 s ≥ 0，其余与 newton._damp 相同。

    s = 0 即 Ṅ_O2 = 0，残差按原式评价（R5 在 c_O2 = 0 处有定义）；只改试探点，不改残差。
    不是取绝对值折返。ADR 0054 审查记录：本条只作诊断，不进入正式实现。
    """
    for index in range(limit + 1):
        newton._clock(deadline)
        point = y + delta / (2.0 ** index)
        point[transform.power] = np.maximum(point[transform.power], 0.0)
        values, status, reason = newton._at(point, evaluate, feasible, transform, explain)
        trial_merit = merit_of(values) if status == "ok" else None
        if trials is not None:
            trials.append({"i": index, "step": float(2.0 ** (-index)), "status": status, "reason": reason,
                           "merit": None if trial_merit is None else float(trial_merit), "accepted": False})
        if status != "ok":
            continue
        if trial_merit < merit:
            if trials is not None:
                trials[-1]["accepted"] = True
            return point, values, trial_merit
    return None


def census(problem: Problem, x, dx) -> dict:
    """沿 Δx 越过 0 下界的分量按组计数与各组最大可行 α。"""
    out = {}
    lower = np.ones(x.size, dtype=bool)
    lower[[slot.index for slot in problem.slots if slot.kind == "temperature"]] = False
    crossing = lower & (x + dx < 0.0)
    for index in np.flatnonzero(crossing):
        name = problem.names[index]
        group = "O2" if name.endswith(".O2") else "gas" if ".gas" in name else "solid"
        alpha = float(x[index] / -dx[index]) if dx[index] < 0 else np.inf
        entry = out.setdefault(group, {"count": 0, "alpha": np.inf, "worst": ""})
        entry["count"] += 1
        if alpha < entry["alpha"]:
            entry["alpha"], entry["worst"] = alpha, name
    return out


def run_inner(problem: Problem, x0, variant: set, limit_s: float, nu_max: int) -> dict:
    """按正式实现的变换与三层规则跑内层；fixedT 时先在温度固定的子系统上收敛，再放开。"""
    started = time.monotonic()
    deadline = lambda: time.monotonic() - started >= limit_s
    config = replace(problem.config, nu_max=nu_max)
    stages = []
    x = np.asarray(x0, dtype=float)
    if "fixedT" in variant:
        free = np.array([slot.index for slot in problem.slots if slot.kind != "temperature"], dtype=int)
        stage = _sub_newton(problem, x, free, config, deadline, "fixedT")
        stages.append(stage)
        x = stage.pop("x") if stage["status"] == "CONVERGED" else x
        if stages[-1]["status"] != "CONVERGED":
            return {"stages": stages, "elapsed_s": time.monotonic() - started, "x": x}
    transform = problem.transform(x)
    rule = problem.rule(transform)
    log = {"iterations": [], "blocks": problem.blocks, "names": problem.names, "explain": rule.explain}
    original = newton._damp
    censuses: list = []
    # o2proj 只用于诊断（外层联立中痕量 O₂ 使幂变量越零，ADR 0054 审查记录）。
    inner_damp = projected_damp if "o2proj" in variant else original

    def damp(y, delta, *args, **kwargs):
        censuses.append(census(problem, transform.to_x(y), transform.to_x(y + delta) - transform.to_x(y)))
        return inner_damp(y, delta, *args, **kwargs)

    newton._damp = damp
    try:
        result = newton.damped_newton(
            x, problem.residual, transform.scaled(problem.scales), problem.converged, rule.trial,
            problem.merit, config, deadline, log, problem.blocks, problem.layout, False,
            transform, rule.final,
        )
    finally:
        newton._damp = original
    summary = _summary("full", result, log, problem)
    summary["final_check"] = rule.final(np.asarray(result.x, dtype=float))
    for entry, found in zip(summary["iterations"], censuses):
        entry["census"] = found
    stages.append(summary)
    return {"stages": stages, "elapsed_s": time.monotonic() - started, "x": result.x}


def _sub_newton(problem, x, free, config, deadline, label):
    """固定 free 以外的分量（温度），对其余行做同一阻尼牛顿。差分逐列。"""
    base = np.asarray(x, dtype=float).copy()

    def embed(y):
        z = base.copy()
        z[free] = y
        return z

    blocks = _sub_blocks(problem, free)
    transform = newton.build_transform(problem.slots, x, problem.n_ref)
    rule = newton.build_rule(
        problem.slots, problem.frozen, transform,
        tau=float(config.stock_roundoff_relative) * problem.m_tot,
        tau_n=float(config.stock_roundoff_relative) * problem.n_ref, unknowns=free,
    )
    log = {"iterations": [], "blocks": blocks, "names": [problem.names[i] for i in free]}
    result = newton.damped_newton(
        base[free], lambda y: problem.residual(embed(y))[free], transform.scaled(problem.scales)[free],
        lambda r: inner_met(r, blocks, config.residual_tolerance), rule.trial,
        lambda r: newton.normalized_merit(r, blocks), config, deadline, log, blocks, None, False,
    )
    summary = _summary(label, result, log, problem)
    summary["x"] = embed(result.x)
    return summary


def _sub_blocks(problem, rows):
    from hamel.solver.gates import ResidualBlocks
    return ResidualBlocks(
        np.flatnonzero(np.isin(rows, problem.blocks.gas)),
        np.flatnonzero(np.isin(rows, problem.blocks.solid)),
        np.flatnonzero(np.isin(rows, problem.blocks.energy)),
        problem.n_ref, problem.m_ref, problem.q_ref,
    )


def _summary(label, result, log, problem) -> dict:
    iterations = []
    for entry in log.get("iterations", []):
        trials = entry.get("trials") or []
        accepted = next((item for item in trials if item.get("accepted")), None)
        reasons: dict = {}
        for trial in trials:
            if trial["status"] != "ok":
                key = str(trial.get("reason", ""))[:60]
                reasons[key] = reasons.get(key, 0) + 1
        iterations.append({
            "iteration": entry["iteration"], "merit": entry["merit"], "blocks": entry["blocks"],
            "step_norm": entry["step"]["norm"], "step_max": entry["step"]["max"],
            "jacobian_ok": entry.get("jacobian_ok"),
            "final_check": entry.get("final_check"),
            "accepted_i": None if accepted is None else accepted["i"],
            "rejections": reasons, "largest": entry["largest"][:5],
        })
    return {"label": label, "status": result.status, "reason": result.reason,
            "iterations": iterations, "final_merit": result.merit,
            "final_ratios": None if result.residual is None or result.residual.size != len(problem.names) else problem.ratios(result.residual)}


def outer(rounds: int, limit_s: float, out: Path, extrapolate: bool) -> dict:
    """完整 solve_case（外层 rounds 轮）。求解手段已并入正式实现，只按需外推 Corella。"""
    if extrapolate:
        _corella_extrapolate()
    rounds_log: list = []
    original = driver.precalculate

    def precalculate(*args, **kwargs):
        result = original(*args, **kwargs)
        frozen = result.frozen
        rounds_log.append({"round": len(rounds_log) + 1, "drain_ash_kg_s": float(frozen.drain_ash_kg_s)})
        return result

    driver.precalculate = precalculate
    restore_alarm = _alarm(limit_s * (rounds + 1) + 7200)
    started = time.monotonic()
    try:
        result = driver.solve_case(str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=rounds))
        status = [result.status, result.reason]
        extra = {"rounds_driver": result.rounds, "nonnegative": result.nonnegative, "balance": result.balance,
                 "measurements": result.measurements, "final_check": result.diagnostics.get("final_check")}
        for row, record in zip(result.diagnostics.get("rounds") or [], rounds_log):
            record.update({"inner": row.get("inner"), "reason": row.get("reason")})
    except (_Timeout, newton.SolveTimeout):
        status, extra = ["TIMEOUT", ""], {}
    except Exception as exc:  # noqa: BLE001 - 诊断要记录任何停止
        status, extra = [type(exc).__name__, str(exc)[:800]], {}
    finally:
        restore_alarm()
        driver.precalculate = original
    return {"status": status, "elapsed_s": time.monotonic() - started, "rounds": rounds_log, **extra}


def _final_check(problem, x) -> dict:
    transform = problem.transform(x)
    rule = problem.rule(transform)
    return rule.final(np.asarray(x, dtype=float))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("--cache", required=True)
    b.add_argument("--loops", type=int, default=6)
    b.add_argument("--limit-s", type=float, default=7200.0)
    i = sub.add_parser("inner")
    i.add_argument("--cache", required=True)
    i.add_argument("--variant", default="")
    i.add_argument("--o2-share", type=float, default=0.5)
    i.add_argument("--nu-max", type=int, default=200)
    i.add_argument("--limit-s", type=float, default=7200.0)
    i.add_argument("--tag", default=None)
    c = sub.add_parser("check")
    c.add_argument("--frozen", required=True)
    c.add_argument("--x", required=True)
    c.add_argument("--tag", required=True)
    o = sub.add_parser("outer")
    o.add_argument("--rounds", type=int, default=1)
    o.add_argument("--limit-s", type=float, default=5400.0)
    o.add_argument("--corellax", action="store_true")
    for item in (b, i, o, c):
        item.add_argument("--out", default=str(_OUT))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.command == "check":
        with Path(args.frozen).open("rb") as handle:
            problem = Problem(pickle.load(handle))
        summary = _final_check(problem, np.load(args.x))
        name = f"check-{args.tag}.json"
    elif args.command == "outer":
        summary = outer(args.rounds, args.limit_s, out, args.corellax)
        name = "outer.json"
    elif args.command == "build":
        summary = build(Path(args.cache), args.loops, args.limit_s)
        name = "build.json"
    else:
        cache = Path(args.cache)
        with (cache / "frozen.pkl").open("rb") as handle:
            problem = Problem(pickle.load(handle))
        x0 = np.load(cache / "x0.npy")
        variant = set(filter(None, args.variant.split(",")))
        if "corellaX" in variant:
            _corella_extrapolate()
        changed = {}
        if "o2mix" in variant:
            x0, changed = o2_mix(problem, x0, args.o2_share)
        restore = _alarm(args.limit_s + 60)
        try:
            summary = run_inner(problem, x0, variant, args.limit_s, args.nu_max)
        except (_Timeout, newton.SolveTimeout):
            summary = {"stages": [], "status": "TIMEOUT"}
        finally:
            restore()
        summary.update(variant=sorted(variant), o2_share=args.o2_share, start_changes=changed)
        tag = args.tag or (args.variant.replace(",", "+") or "base")
        name = f"inner-{tag}.json"
        if summary.get("x") is not None:
            np.save(out / f"x-{tag}.npy", summary["x"])
    (out / name).write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=_plain), encoding="utf-8")
    print(name, summary.get("status") or [stage.get("status") for stage in summary.get("stages", [])])


def _plain(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return str(value)


if __name__ == "__main__":
    main()
