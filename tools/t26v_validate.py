"""T26v：主线 Corella 外推的基准外层、S20 截断与 T26t D.3 的逐位对照。

不在 ``tools`` 里打 Corella 补丁：基准用默认配置，S20 只改 ``corella_high_temperature``。不记录出口组成、T̄_FRB、X_C。
``python -m tools.t26v_validate baseline|s20|compare [--limit-s …]``
"""
from __future__ import annotations
import argparse, dataclasses, json, pickle, signal, time
from dataclasses import fields
from pathlib import Path
import numpy as np
import hamel.solver.driver as driver, hamel.solver.initial as initial
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.kinetics.common import KineticsInputError
from hamel.solver.gates import _abs_max
from hamel.solver.newton import SolveTimeout
from hamel.solver.vector import unpack_state

_ROOT = Path(__file__).resolve().parents[1]
_SCRATCH, _CASE = _ROOT / ".scratch/t26v", _ROOT / "data/cases/vtt_pressurized.json"
_OLD = _ROOT / ".scratch/t26t/final_args.pkl"
_ROUNDS, _LIMIT, _OLD_ITERATIONS = 26, 8.0 * 3600.0, 181  # T26t D.3：第 17 轮、内层 181 次


class _Timeout(BaseException):
    """硬墙钟。继承 BaseException：任何 ``except Exception`` 都不会吞掉它（CHARTER §7）。"""


def _write(name, payload):
    _SCRATCH.mkdir(parents=True, exist_ok=True)
    (_SCRATCH / name).write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=str), encoding="utf-8")


def run(label, truncate, limit_s=_LIMIT):
    """一次外层运行：记录预计算入参、试探点最底床格温度与 Corella 拒绝次数。"""
    assembled = assemble_case(str(_CASE), "1")
    calls, monitors, trial = [], [], {"bottom": -np.inf, "reject_k": -np.inf, "rejections": 0}
    real_pre, real_eval, real_mon = driver.precalculate, driver.evaluate_residual, driver._monitor

    def evaluate(values, frozen):
        try:
            bottom = float(unpack_state(values, frozen).bed_temperature_k[-1])
        except Exception:  # noqa: BLE001 - 不可解包的试探点交给正式路径报错
            return real_eval(values, frozen)
        trial["bottom"] = max(trial["bottom"], bottom)
        try:
            return real_eval(values, frozen)
        except KineticsInputError as exc:
            if "Corella" not in str(exc):
                raise
            trial["rejections"] += 1
            trial["reject_k"] = max(trial["reject_k"], bottom)
            raise

    def precalculate(*args, **kwargs):
        calls.append((args, dict(kwargs)))
        return real_pre(*args, **kwargs)

    def monitor(state, frozen, assembly=None):
        value = real_mon(state, frozen, assembly)
        monitors.append((value, float(np.sum(frozen.geometric_mass_kg))))
        return value

    real_initial_eval, initial.evaluate_residual = initial.evaluate_residual, evaluate
    driver.precalculate, driver.evaluate_residual, driver._monitor = precalculate, evaluate, monitor
    config = SolverConfig(outer_round_limit=_ROUNDS, corella_high_temperature="truncate" if truncate else "extrapolate")
    started, result, status, reason = time.monotonic(), None, "EXCEPTION", ""

    def expire(signum, frame):
        raise _Timeout()

    # 硬墙钟：初值构造、预计算与单次残差／雅可比都不轮询 deadline，只能靠信号中断（CHARTER §7）。
    # 用 itimer 取亚秒精度（alarm 只有整秒），并留 50 ms 余量避免与收尾清理竞争。
    previous = signal.signal(signal.SIGALRM, expire)
    signal.setitimer(signal.ITIMER_REAL, max(float(limit_s), 0.0) + 0.05)
    try:
        try:
            result = driver.solve_case(str(_CASE), "1", config=config, deadline=lambda: time.monotonic() - started >= limit_s)
            status, reason = result.status, result.reason
        except (SolveTimeout, _Timeout):
            status, reason = "TIMEOUT", f"运行超过墙钟上限 {limit_s:.0f} s"
        except Exception as exc:  # noqa: BLE001 - 求解异常记进报告，不当作 TIMEOUT
            reason = str(exc)[:800] or type(exc).__name__
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        driver.precalculate, driver.evaluate_residual, driver._monitor = real_pre, real_eval, real_mon
        initial.evaluate_residual = real_initial_eval
    summary = {"label": label, "status": status, "reason": reason, "elapsed_s": time.monotonic() - started,
               "trial_bottom_max_k": trial["bottom"], "trial_rejections": trial["rejections"],
               "trial_reject_bottom_max_k": trial["reject_k"]}
    if result is None:
        _write(f"{label}.json", summary)
        return summary
    diagnostics, (n_ref, m_ref) = result.diagnostics or {}, driver._references(assembled.physical)[:2]
    trace = []
    for index, item in enumerate(diagnostics.get("rounds") or []):
        row = {"round": item.get("round"), "inner": item.get("inner"), "reason": item.get("reason"),
               "iterations": len(item.get("newton") or []), "switched": bool(item.get("switched"))}
        if 0 < index < len(monitors):
            now, m_tot, was = monitors[index][0], monitors[index][1], monitors[index - 1][0]
            row["gate"] = {"dT": _abs_max(now.temperature_k, was.temperature_k),
                           "solid": _abs_max(now.solid_flow_kg_s, was.solid_flow_kg_s) / m_ref,
                           "inventory": _abs_max(now.inventory_kg, was.inventory_kg) / m_tot,
                           "gas": _abs_max(now.gas_flow_mol_s, was.gas_flow_mol_s) / n_ref}
        trace.append(row)
    summary.update({
        "rounds": result.rounds, "iterations": result.iterations, "ratios": result.ratios,
        "balance": result.balance, "nonnegative": result.nonnegative, "trace": trace,
        "initial": diagnostics.get("initial"), "consistent": diagnostics.get("consistent"),
        "final_check": diagnostics.get("final_check"), "consistent_ash_relative": diagnostics.get("consistent_ash_relative"),
        "bottom_relative": diagnostics.get("bottom_relative"), "negative_floor": diagnostics.get("negative_floor"),
        "bed_temperature_peak": diagnostics.get("bed_temperature_peak"),
        "final_bed_k": monitors[-1][0].temperature_k[len(((diagnostics.get("final") or {}).get("a_free_per_s"))):].tolist()
        if monitors and (diagnostics.get("final") or {}).get("a_free_per_s") is not None else None,
    })
    first = next(((args, kwargs) for args, kwargs in calls
                  if kwargs.get("iteration") == result.rounds and kwargs.get("consistency")), None)
    if first is not None:
        with (_SCRATCH / f"final_args_{label}.pkl").open("wb") as handle:
            pickle.dump({"args": first[0], "kwargs": first[1]}, handle, protocol=pickle.HIGHEST_PROTOCOL)
        summary["final_args_saved"] = sum(1 for args, kw in calls
                                          if kw.get("iteration") == result.rounds and kw.get("consistency"))
    _write(f"{label}.json", summary)
    return summary


def _solver(saved):
    """旧进程 pickle 的 SolverConfig 没有 ADR 0059 的新字段：补入默认外推后再比较。"""
    names = {item.name for item in fields(SolverConfig)}
    data = {key: value for key, value in vars(saved).items() if key in names}
    data.setdefault("corella_high_temperature", "extrapolate")
    return SolverConfig(**data)


def _differences(left, right, path, out):
    """数组用 ``np.array_equal``，其它值用 ``==``；dataclass 与字典逐层进入。"""
    if dataclasses.is_dataclass(left) and not isinstance(left, type):
        if type(left) is not type(right):
            return out.append(f"{path}: 类型不同")
        for item in dataclasses.fields(left):
            _differences(getattr(left, item.name), getattr(right, item.name), f"{path}.{item.name}", out)
    elif isinstance(left, dict) and isinstance(right, dict):
        if set(left) != set(right):
            return out.append(f"{path}: 键不同")
        for key in left:
            _differences(left[key], right[key], f"{path}.{key}", out)
    elif isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            return out.append(f"{path}: 长度不同")
        for index, (first, second) in enumerate(zip(left, right)):
            _differences(first, second, f"{path}[{index}]", out)
    elif isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        if not np.array_equal(np.asarray(left), np.asarray(right)):
            out.append(path)
    elif left != right:
        out.append(f"{path}: {left!r} / {right!r}")


def compare():
    with _OLD.open("rb") as handle:
        old = pickle.load(handle)
    with (_SCRATCH / "final_args_baseline.pkl").open("rb") as handle:
        new = pickle.load(handle)
    old_kwargs, new_kwargs = dict(old["kwargs"]), dict(new["kwargs"])
    backfill = "corella_high_temperature" not in vars(old_kwargs["solver"])
    old_kwargs["solver"] = _solver(old_kwargs["solver"])
    differences = []
    _differences(list(old["args"]), list(new["args"]), "args", differences)
    _differences(old_kwargs, new_kwargs, "kwargs", differences)
    baseline = json.loads((_SCRATCH / "baseline.json").read_text(encoding="utf-8"))
    s20 = json.loads((_SCRATCH / "s20.json").read_text(encoding="utf-8"))
    report = {"solver_backfill": backfill, "equal": not differences, "differences": differences,
              "old_rounds": old_kwargs.get("iteration"), "old_iterations": _OLD_ITERATIONS,
              "new_rounds": baseline.get("rounds"), "new_iterations": baseline.get("iterations"),
              "baseline_status": baseline.get("status"), "s20_status": s20.get("status")}
    if baseline.get("final_bed_k") and s20.get("final_bed_k"):
        left, right = np.asarray(baseline["final_bed_k"]), np.asarray(s20["final_bed_k"])
        report["s20_bed_delta_k"] = (right - left).tolist() if right.shape == left.shape else None
    _write("compare.json", report)
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("baseline", "s20", "compare"))
    parser.add_argument("--limit-s", type=float, default=_LIMIT)
    args = parser.parse_args(argv)
    if args.mode == "compare":
        print(json.dumps(compare(), ensure_ascii=False, default=str))
        return 0
    summary = run(args.mode, args.mode == "s20", args.limit_s)
    print(json.dumps({key: summary.get(key) for key in
                      ("label", "status", "reason", "rounds", "iterations", "elapsed_s", "trial_rejections")},
                     ensure_ascii=False, default=str))
    # CHARTER §7：超时记 TIMEOUT 并以非零状态退出；求解异常也非零。
    return {"TIMEOUT": 2}.get(summary.get("status"), 0 if summary.get("status") != "EXCEPTION" else 3)


if __name__ == "__main__":
    raise SystemExit(main())
