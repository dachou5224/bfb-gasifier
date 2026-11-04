"""限时求解路径诊断（2026-11 目标所有者决定第 4 项）。只报告，不改 hamel/。

判据：从冻结模型出发，基准外层第 1 轮内层 CONVERGED；外层顺带看 2–3 轮的趋势。

子命令
- build：跑预计算与初值构造（ADR 0052），把交给内层的起点 x0 与冻结量存进 --cache。
  与正式写法只差一处（候选 2，三层试探规则）：预计算中 ADR 0014 A5 的余额为负时不停止，
  取同式负值继续交替，A5 的非负只在收敛解上看（tools/alternation_diagnostic.signed_drain）。
- inner：读 --cache，按 --variant 跑第 1 轮内层，逐次记录牛顿迭代、阻尼试探与障碍普查。
- outer：完整 solve_case（外层 --rounds 轮），内层按 --variant 替换起点或求解手段。

候选手段（--variant，可用逗号组合）
- base：原样（ADR 0043 差分步长、ADR 0017 第 7 项可行条件）。
- o2mix：进风格 O₂ 出流取进风 O₂ 的份额 --o2-share，只改起点。
- fixedT：先固定全部温度解气相与固相，再放开温度（Hamel S. 15 的分阶段）。
- o2step：O₂ 列的差分步长改为 √ε·max(|x|, 1e-12·ṅ_ref)。
- o2log：O₂ 出流取 ln Ṅ。
- o2pow：O₂ 出流取 Ṅ^n，n 为该相 R5 的 O₂ 级数（0,5／0,3）。
- gaslog：其余（起点为正的）气体出流取 ln Ṅ。
- o2inlet：进风格的 O₂ 取幂变量，其余格的 O₂ 取对数（外层联立中试过：痕量 O₂ 下溢，矩阵奇异）。
- o2proj：试探点上把 O₂ 幂变量 s 投影到 0（不折返）。
- relax：三层试探规则（试探点只要求气体 ≥ 0、T > 0；非负在收敛解上检查）。
- corellaX：**仅诊断、属模型有效域放宽**：Corella 高温段外推到 900 °C 以上。
  outer 默认不带；用到时须在结论中注明（M5 审计）。

所有结果写在 --out（默认 reports/research-solver-path/diagnostic/）。
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
from tools.alternation_diagnostic import signed_drain

_ROOT = Path(__file__).resolve().parents[1]
_CASE = _ROOT / "data" / "cases" / "vtt_pressurized.json"
_OUT = _ROOT / "reports" / "research-solver-path" / "diagnostic"
_O2 = GAS_SPECIES.index("O2")


class _Timeout(BaseException):
    """继承 BaseException：任何 except Exception 都不会吞掉全局期限（CHARTER §7）。"""


class _Built(BaseException):
    """初值构造已结束：跳出 solve_case。"""


def _patch_a5():
    """预计算的 A5 余额为负时取负值继续（候选 2）。返回还原函数。"""
    original = precalc._drain

    def drain(*args):
        try:
            return original(*args)
        except precalc.PrecalcError as exc:
            if "带出灰" not in str(exc):
                raise
            return signed_drain(*args)

    precalc._drain = drain
    return lambda: setattr(precalc, "_drain", original)


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

    restore_a5 = _patch_a5()
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
    except Exception as exc:
        status = [type(exc).__name__, str(exc)[:800]]
    finally:
        restore_alarm()
        restore_a5()
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


# ---------------------------------------------------------------- 内层

def final_violations(x, names, stock_mask, gas_mask, inert, tau, tau_n) -> dict:
    """第 3 层：存量 ≥ −τ、气体出流 ≥ −τ_n、各床格惰性物 ≥ 0。

    分开报告舍入容差内的负值（不算违反）与真正的违反；ok 只看后者。
    """
    x = np.asarray(x, dtype=float)
    stock_round = np.flatnonzero(stock_mask & (x < 0.0) & (x >= -tau))
    stock_bad = np.flatnonzero(stock_mask & (x < -tau))
    gas_bad = np.flatnonzero(gas_mask & (x < -tau_n))
    inert_bad = np.flatnonzero(np.asarray(inert) < 0.0)
    return {
        "ok": bool(stock_bad.size == 0 and gas_bad.size == 0 and inert_bad.size == 0),
        "tau_kg": tau, "tau_n_mol_s": tau_n,
        "stocks_within_tau": {names[i]: float(x[i]) for i in stock_round},
        "stocks_below_tau": {names[i]: float(x[i]) for i in stock_bad[np.argsort(x[stock_bad])][:8]},
        "gas_below_tau_n": {names[i]: float(x[i]) for i in gas_bad[:8]},
        "inert_kg": np.asarray(inert, dtype=float).tolist(),
        "inert_negative_cells": inert_bad.tolist(),
    }


class Problem:
    """第 1 轮内层的原问题：与 driver.solve_aligned 的 prepare 同一组块、尺度与可行条件。"""

    def __init__(self, frozen):
        assembled = assemble_case(str(_CASE), "1")
        self.config = assembled.solver
        self.frozen = frozen
        self.n_ref, self.m_ref, self.q_ref = driver._references(assembled.physical)
        self.blocks = driver._blocks(frozen, self.n_ref, self.m_ref, self.q_ref)
        self.m_tot = float(frozen.geometric_mass_kg.sum())
        self.slots = unknown_slots(frozen)
        self.names = [slot.name for slot in self.slots]
        self.layout = jacobian_layout(self.slots, frozen)
        self.scales = driver._scales(frozen, self.n_ref, self.m_tot)
        self.o2 = np.array([i for i, name in enumerate(self.names) if name.endswith(".O2")], dtype=int)
        self.temperature = np.array([s.index for s in self.slots if s.kind == "temperature"], dtype=int)
        self.evaluations = 0

    def residual(self, x):
        self.evaluations += 1
        return evaluate_residual(x, self.frozen).values

    relaxed = False

    def feasible(self, x):
        if not self.relaxed:
            return driver._feasible(x, self.frozen)
        # 候选 2（三层规则）：试探点只要求气体出流 ≥ 0、温度 > 0（速率式的定义域）；
        # 存量与床格惰性物的非负只在收敛解上检查（final_check）。
        state = unpack_state(x, self.frozen)
        gas = np.concatenate([state.bed_gas_mol_s.ravel(), state.freeboard_gas_mol_s.ravel()])
        temperature = np.concatenate([state.bed_temperature_k, state.freeboard_temperature_k])
        return bool(np.all(gas >= 0.0) and np.all(temperature > 0.0))

    def final_check(self, x) -> dict:
        """收敛解上的第 3 层检查（#145 审查 P2）。惰性物不在 x 中，逐格由 M − 焦炭 − 燃料算出。"""
        state = unpack_state(x, self.frozen)
        inert = self.frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - self.frozen.fuel_bed_kg.sum(axis=1)
        lower = np.ones(x.size, dtype=bool)
        lower[self.temperature] = False
        stocks = np.zeros(x.size, dtype=bool)
        stocks[[s.index for s in self.slots if s.kind in ("char", "ash")]] = True
        return final_violations(
            x, self.names, stocks, lower & ~stocks, np.asarray(inert, dtype=float),
            float(self.config.stock_roundoff_relative) * self.m_tot,
            float(self.config.stock_roundoff_relative) * self.n_ref,
        )

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
    # 床区最底格（衡算下标 -1）：两相各按一次风分配。
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


def round_handoff(problem: Problem, x0):
    """障碍 O1：初值构造的收敛解允许存量 ≥ −τ（ADR 0052 修订），内层可行条件要求 ≥ 0。

    交接处把 [−τ, 0) 的存量记为 0（τ = stock_roundoff_relative·M_tot），返回被改的项。
    """
    x = np.array(x0, dtype=float)
    tau = float(problem.config.stock_roundoff_relative) * problem.m_tot
    lower = np.ones(x.size, dtype=bool)
    lower[problem.temperature] = False
    hit = np.flatnonzero(lower & (x < 0.0) & (x >= -tau))
    rounded = {problem.names[i]: float(x[i]) for i in hit}
    x[hit] = 0.0
    return x, rounded


def census(problem: Problem, x, dx) -> dict:
    """沿 Δx 越过 0 下界的分量按组计数与各组最大可行 α。"""
    out = {}
    lower = np.ones(x.size, dtype=bool)
    lower[problem.temperature] = False
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


def _inlet_o2(problem: Problem) -> set:
    """进风格的 O₂ 出流下标：床区最底格两相（一次风）与二次风所在的自由空间格。"""
    frozen = problem.frozen
    bottom = f"bed[b={frozen.n_bed - 1},"
    secondary = f"freeboard[b={int(frozen.secondary_balance_index)},"
    return {int(i) for i in problem.o2 if problem.names[i].startswith((bottom, secondary))}


def transforms(problem: Problem, variant: set, x, scales):
    """候选 4：只换牛顿所用的自变量，原方程不变（CHARTER §5 的求解手段，先例 ADR 0042）。

    - o2log：O₂ 出流取 z = ln Ṅ。
    - o2pow：O₂ 出流取 s = Ṅ^n，n 为该相 R5 的 O₂ 级数（悬浮相 0,5；气泡相与自由空间 0,3，
      Hamel Tab. 5.3, S. 91），R5 在 s 上近似线性。
    - gaslog：其余气体出流（起点为正者）取 z = ln Ṅ；恒为 0 的分量（Ar）不变。
    负 s 映到负 Ṅ，由原可行条件拒绝。差分与列缩放的尺度：对数列取 1，幂列取 ṅ_ref^n。
    """
    log_idx, pow_idx, order = [], [], []
    inlet = _inlet_o2(problem)
    for index in problem.o2:
        if "o2inlet" in variant and index not in inlet:
            # 非进风格的 O₂ 只来自下游残余，真解为痕量：取对数变量（与 gaslog 同）。
            log_idx.append(index)
        elif "o2pow" in variant or "o2inlet" in variant:
            pow_idx.append(index)
            order.append(0.5 if ".suspension." in problem.names[index] else 0.3)
        elif "o2log" in variant:
            log_idx.append(index)
    if "gaslog" in variant:
        o2 = set(problem.o2.tolist())
        log_idx += [i for i in problem.blocks.gas.tolist() if i not in o2 and x[i] > 0.0]
    log_idx, pow_idx, order = np.array(log_idx, dtype=int), np.array(pow_idx, dtype=int), np.array(order)

    def to_x(y):
        v = np.array(y, dtype=float)
        with np.errstate(over="ignore"):
            v[log_idx] = np.exp(v[log_idx])
        v[pow_idx] = np.sign(v[pow_idx]) * np.abs(v[pow_idx]) ** (1.0 / order)
        return v

    def from_x(v):
        y = np.array(v, dtype=float)
        y[log_idx] = np.log(y[log_idx])
        y[pow_idx] = np.sign(y[pow_idx]) * np.abs(y[pow_idx]) ** order
        return y

    scales = np.array(scales, dtype=float)
    scales[log_idx] = 1.0
    scales[pow_idx] = problem.n_ref ** order
    return to_x, from_x, scales


def _projected_damp(x, delta, bound, merit, evaluate, feasible, merit_of, limit, deadline, trials=None, explain=None):
    """o2proj：与 newton._damp 相同的减半与严格下降，只是试探点上把 O₂ 幂变量 s 投影到 s ≥ 0。

    s = 0 即 Ṅ_O2 = 0，残差按原式评价（R5 在 c_O2 = 0 有定义）；投影不改残差，只改试探点。
    不是取绝对值折返：负 s 只截到边界。
    """
    for index in range(limit + 1):
        newton._clock(deadline)
        point = x + delta / (2.0 ** index)
        point[bound] = np.maximum(point[bound], 0.0)
        values, status, reason = newton._at(point, evaluate, feasible, explain)
        trial_merit = merit_of(values) if status == "ok" else None
        if trials is not None:
            trials.append({"i": index, "step": float(2.0 ** (-index)), "status": status, "reason": reason,
                           "merit": None if trial_merit is None else float(trial_merit), "accepted": False})
        if status == "domain":
            raise newton.FreeboardDomainError("试探点落在自由空间有效域外")
        if status != "ok":
            continue
        if trial_merit < merit:
            if trials is not None:
                trials[-1]["accepted"] = True
            return point, values, trial_merit
    return None


def run_inner(problem: Problem, x0, variant: set, limit_s: float, nu_max: int, journal: dict) -> dict:
    """候选组合下跑内层。fixedT 时先在温度固定的子系统上收敛，再放开。"""
    started = time.monotonic()
    deadline = lambda: time.monotonic() - started >= limit_s
    config = replace(problem.config, nu_max=nu_max)
    if "o2step" in variant:
        scales = problem.scales.copy()
        scales[problem.o2] = 1.0e-12 * problem.n_ref
    else:
        scales = problem.scales
    stages = []
    problem.relaxed = "relax" in variant
    x = np.asarray(x0, dtype=float)
    if "fixedT" in variant:
        free = np.setdiff1d(np.arange(x.size), problem.temperature)
        stage = _sub_newton(problem, x, free, scales, config, deadline, "fixedT")
        stages.append(stage)
        x = stage.pop("x")
        if stage["status"] != "CONVERGED":
            return {"stages": stages, "elapsed_s": time.monotonic() - started, "evaluations": problem.evaluations}
    to_x, from_x, scales = transforms(problem, variant, x, scales)
    x = from_x(x)
    log = {"iterations": [], "blocks": problem.blocks, "names": problem.names,
           "explain": lambda trial: driver._infeasible_clause(to_x(trial), problem.frozen)}
    original = newton._damp
    censuses = []
    pow_idx = problem.o2 if ("o2pow" in variant or "o2inlet" in variant) else np.array([], dtype=int)

    def damp(x, delta, *args, **kwargs):
        censuses.append(census(problem, to_x(x), to_x(x + delta) - to_x(x)))
        if "o2proj" not in variant:
            return original(x, delta, *args, **kwargs)
        return _projected_damp(x, delta, pow_idx, *args, **kwargs)

    newton._damp = damp
    try:
        result = newton.damped_newton(
            x, lambda y: problem.residual(to_x(y)), scales, problem.converged,
            lambda y: bool(np.isfinite(to_x(y)).all()) and problem.feasible(to_x(y)), problem.merit,
            config, deadline, log, problem.blocks, problem.layout, False,
        )
    finally:
        newton._damp = original
    result.x = to_x(result.x)
    summary = _summary("full", result, log, problem)
    summary["final_check"] = problem.final_check(result.x)
    for entry, found in zip(summary["iterations"], censuses):
        entry["census"] = found
    stages.append(summary)
    journal["x_final"] = result.x
    return {"stages": stages, "elapsed_s": time.monotonic() - started, "evaluations": problem.evaluations}


def _sub_newton(problem, x, free, scales, config, deadline, label):
    """固定 free 以外的分量（温度），对其余行做同一阻尼牛顿。差分逐列（布局只对全问题有效）。"""
    base = np.asarray(x, dtype=float).copy()
    rows = free

    def embed(y):
        z = base.copy()
        z[rows] = y
        return z

    blocks = _sub_blocks(problem, rows)
    log = {"iterations": [], "blocks": blocks, "names": [problem.names[i] for i in rows]}
    result = newton.damped_newton(
        base[rows], lambda y: problem.residual(embed(y))[rows], scales[rows],
        lambda r: inner_met(r, blocks, config.residual_tolerance),
        lambda y: problem.feasible(embed(y)), lambda r: newton.normalized_merit(r, blocks),
        config, deadline, log, blocks, None, False,
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
        accepted = next((t for t in trials if t.get("accepted")), None)
        reasons: dict = {}
        for trial in trials:
            if trial["status"] != "ok":
                key = str(trial.get("reason", ""))[:60]
                reasons[key] = reasons.get(key, 0) + 1
        iterations.append({
            "iteration": entry["iteration"], "merit": entry["merit"], "blocks": entry["blocks"],
            "step_norm": entry["step"]["norm"], "step_max": entry["step"]["max"],
            "accepted_i": None if accepted is None else accepted["i"],
            "accepted_merit": None if accepted is None else accepted["merit"],
            "rejections": reasons, "largest": entry["largest"][:5],
        })
    return {"label": label, "status": result.status, "reason": result.reason,
            "iterations": iterations, "final_merit": result.merit,
            "final_ratios": None if result.residual is None or result.residual.size != len(problem.names) else problem.ratios(result.residual)}


def outer(variant: set, rounds: int, limit_s: float, out: Path) -> dict:
    """完整 solve_case（外层 rounds 轮）：初值构造同 build，内层换成 variant 的组合。只报告。"""
    frozen_box: dict = {}
    original_pre, original_newton = driver.precalculate, driver.damped_newton
    rounds_log: list = []

    def pre(*args, **kwargs):
        result = original_pre(*args, **kwargs)
        frozen_box["frozen"] = result.frozen
        return result

    def inner(x0, evaluate, scales, converged, feasible, merit_of, config, deadline=None, journal=None, blocks=None, layout=None, verify=False):
        problem = Problem(frozen_box["frozen"])
        start, rounded = round_handoff(problem, x0)
        local: dict = {}
        summary = run_inner(problem, start, variant, limit_s, int(config.nu_max), local)
        stage = summary["stages"][-1]
        x = local.get("x_final", start)
        state = unpack_state(x, problem.frozen)
        rounds_log.append({
            "round": len(rounds_log) + 1, "status": stage["status"], "reason": stage["reason"],
            "iterations": len(stage["iterations"]), "final_ratios": stage["final_ratios"],
            "final_check": stage.get("final_check"), "handoff_rounded": rounded,
            "elapsed_s": summary["elapsed_s"],
            "bed_temperature_c": (state.bed_temperature_k - 273.15).tolist(),
            "freeboard_temperature_c": (state.freeboard_temperature_k - 273.15).tolist(),
            "drain_ash_kg_s": float(problem.frozen.drain_ash_kg_s),
            "m_dot_ab_kg_s": np.asarray(problem.frozen.m_dot_ab_kg_s).tolist(),
        })
        np.save(out / f"outer-x-round{len(rounds_log)}.npy", x)
        with (out / f"outer-frozen-round{len(rounds_log)}.pkl").open("wb") as handle:
            pickle.dump(problem.frozen, handle, protocol=pickle.HIGHEST_PROTOCOL)
        (out / "outer.json").write_text(json.dumps({"rounds": rounds_log}, ensure_ascii=False, indent=1, default=_plain), encoding="utf-8")
        values = None if x is None else problem.residual(x)
        status = stage["status"]
        if status == "CONVERGED" and not (stage.get("final_check") or {}).get("ok", False):
            status = "FINAL_CHECK_FAILED"
            rounds_log[-1]["status"] = status
        return newton.NewtonResult(status, stage["reason"], x, values, len(stage["iterations"]), stage["final_merit"])

    restore_a5 = _patch_a5()
    driver.precalculate, driver.damped_newton = pre, inner
    restore_alarm = _alarm(limit_s * (rounds + 1) + 7200)
    started = time.monotonic()
    try:
        result = driver.solve_case(str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=rounds))
        status = [result.status, result.reason]
        extra = {"rounds_driver": result.rounds, "nonnegative": result.nonnegative, "balance": result.balance,
                 "measurements": result.measurements}
    except (_Timeout, newton.SolveTimeout):
        status, extra = ["TIMEOUT", ""], {}
    except Exception as exc:
        status, extra = [type(exc).__name__, str(exc)[:800]], {}
    finally:
        restore_alarm()
        restore_a5()
        driver.precalculate, driver.damped_newton = original_pre, original_newton
    return {"status": status, "elapsed_s": time.monotonic() - started, "rounds": rounds_log, **extra}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("--cache", required=True)
    b.add_argument("--loops", type=int, default=6)
    b.add_argument("--limit-s", type=float, default=7200.0)
    i = sub.add_parser("inner")
    i.add_argument("--cache", required=True)
    i.add_argument("--variant", default="base")
    i.add_argument("--o2-share", type=float, default=0.5)
    i.add_argument("--nu-max", type=int, default=60)
    i.add_argument("--limit-s", type=float, default=7200.0)
    i.add_argument("--tag", default=None)
    c = sub.add_parser("check")
    c.add_argument("--frozen", required=True)
    c.add_argument("--x", required=True)
    c.add_argument("--tag", required=True)
    o = sub.add_parser("outer")
    o.add_argument("--variant", default="gaslog,o2pow,relax")
    o.add_argument("--rounds", type=int, default=3)
    o.add_argument("--limit-s", type=float, default=5400.0)
    for p in (b, i, o, c):
        p.add_argument("--out", default=str(_OUT))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.command == "check":
        with Path(args.frozen).open("rb") as handle:
            problem = Problem(pickle.load(handle))
        summary = problem.final_check(np.load(args.x))
        name = f"check-{args.tag}.json"
    elif args.command == "outer":
        variant = set(args.variant.split(","))
        if "corellaX" in variant:
            import hamel.kinetics.homogeneous as homogeneous
            homogeneous._CORELLA_HIGH_K = float("inf")
        summary = outer(variant, args.rounds, args.limit_s, out)
        name = "outer.json"
    elif args.command == "build":
        summary = build(Path(args.cache), args.loops, args.limit_s)
        name = "build.json"
    else:
        cache = Path(args.cache)
        with (cache / "frozen.pkl").open("rb") as handle:
            frozen = pickle.load(handle)
        problem = Problem(frozen)
        x0 = np.load(cache / "x0.npy")
        variant = set(args.variant.split(","))
        if "corellaX" in variant:
            # 仅诊断：Corella 高温段参数外推到 900 °C 以上（ADR 0010 第 9 条未批准，不进入正式路径）。
            import hamel.kinetics.homogeneous as homogeneous
            homogeneous._CORELLA_HIGH_K = float("inf")
        x0, handoff = round_handoff(problem, x0)
        changed = {"handoff_rounded": handoff}
        if "o2mix" in variant:
            x0, changed = o2_mix(problem, x0, args.o2_share)
        journal: dict = {}
        restore = _alarm(args.limit_s + 60)
        try:
            summary = run_inner(problem, x0, variant, args.limit_s, args.nu_max, journal)
        except (_Timeout, newton.SolveTimeout):
            summary = {"stages": [], "status": "TIMEOUT"}
        finally:
            restore()
        summary.update(variant=sorted(variant), o2_share=args.o2_share, start_changes=changed)
        tag = args.tag or args.variant.replace(",", "+")
        name = f"inner-{tag}.json"
        if "x_final" in journal:
            np.save(out / f"x-{tag}.npy", journal["x_final"])
    (out / name).write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=_plain), encoding="utf-8")
    print(name, json.dumps(summary.get("status") or summary.get("ok") if "ok" in summary else summary.get("status") or [s.get("status") for s in summary.get("stages", [])], ensure_ascii=False))


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
