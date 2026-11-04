"""O7 的外层不动点检查：欠松弛外层（M5 审计 codex 第 1 节；目标所有者 2026-11-04 决定）。只报告，不改 hamel/。

问题：是否存在外层不动点 x* = 内层(预计算(x*))，并且在 x* 上 A5 排灰灰分量 ≥ 0？
K1–K2（`tools/o7_checks.py`）只在外层第 1、2 轮的两组冻结系数上做了检查。本工具让冻结系数随外层更新。

写法（诊断，与正式路径的差别逐条列出）：
1. 初值构造同 `solver_path_diagnostic.py build`（A5 中间余额取负值继续）。
2. 第 k 轮：用状态 s_k 做预计算（与 driver 相同的输入：s_k、床区反应汇、上一次的排灰、
   燃料回落与释放分布），得冻结量 q_k；内层从**上一轮的内层解** s̃_{k−1} 出发（第 1 轮从初值构造的状态），
   用 `gaslog,o2pow,relax[,corellaX]` 求解，得 s̃_k。
3. 欠松弛只作用于预计算的输入：s_{k+1} = s_k + ω·(s̃_k − s_k)（ReactorState 的全部数组逐元素；
   反应汇同样欠松弛）。内层起点不取混合态：气相线性混合会在进风格造出 H₂ 与 O₂ 并存的状态，
   R6 残差达 10⁵ 量级（第一次试跑 ω = 0,25、0,5 的第 2 轮都在 6 次后雅可比不可解）。
4. 预计算中 A5 余额为负时取负值进入冻结的排灰（同诊断）。传给下一次预计算的“上一次排灰”
   （决定出口系数 k_exit 与向下流递推的外排）为负时取 0，**这是路径辅助**：
   负的出口系数使燃料释放的输运系数为负，预计算按定义域停止（第一次外层运行的第 3 轮）。
   这一替换只在中间轮起作用；不动点若满足 A5 ≥ 0，该替换不起作用。每轮记录它是否起作用。
5. 停止：相邻两轮温度变化 < 1 K、床区焦炭与 A5 余额的相对变化 < 1·10⁻³，或做满 --rounds。

不动点只在满足以下各条时才算：内层 CONVERGED、收敛解第 3 层检查通过、A5 ≥ 0、
第 4 条的替换未起作用；corellaX 下还须报告最高床温是否 ≤ 900 °C。
只用质量流与存量判 O7，不读取出口验收量。
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

import hamel.solver.driver as driver
import hamel.solver.precalc as precalc
from hamel.cases.config import SolverConfig
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import ReactorState, pack_state, unpack_state
from tools.o7_checks import balance
from tools.solver_path_diagnostic import _CASE, Problem, _Built, _patch_a5, run_inner

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "research-solver-path" / "outer-coupled"
_FIELDS = (
    "bed_gas_mol_s", "freeboard_gas_mol_s", "bed_temperature_k", "freeboard_temperature_k", "bed_char_kg",
    "freeboard_char_free_kg", "freeboard_char_hover_kg", "freeboard_ash_free_kg", "freeboard_ash_hover_kg",
)


def construct() -> object:
    """初值构造，返回最后一次预计算的 PrecalcResult。"""
    captured: dict = {}
    original = driver.construct

    def capture(*args, **kwargs):
        out = original(*args, **kwargs)
        captured["result"] = out[0]
        raise _Built()

    driver.construct = capture
    try:
        driver.solve_case(str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=1))
    except _Built:
        pass
    finally:
        driver.construct = original
    return captured["result"]


def blend(old: ReactorState, new: ReactorState, omega: float) -> ReactorState:
    values = {name: np.asarray(getattr(old, name), dtype=float) + omega * (
        np.asarray(getattr(new, name), dtype=float) - np.asarray(getattr(old, name), dtype=float)) for name in _FIELDS}
    return ReactorState(**values)


def run(rounds: int, omega: float, variant: set, limit_s: float, out: Path) -> dict:
    restore = _patch_a5()
    log: dict = {"omega": omega, "variant": sorted(variant), "rounds": []}
    started = time.monotonic()
    try:
        first = construct()
        log["construct_s"] = time.monotonic() - started
        frozen0 = first.frozen
        state = first.state
        reaction, _ = driver.solid_feedback(state, frozen0)
        drain = float(frozen0.drain_ash_kg_s + frozen0.drain_char_kg_s.sum())
        fuel_return = float(frozen0.fuel_down_free_kg_s[-1].sum())
        release = first.solid_context.released
        previous = None
        warm = state
        for k in range(1, rounds + 1):
            substituted = drain < 0.0
            if k == 1:
                # 第 1 轮与 driver 相同：直接用初值构造最后一次预计算的冻结量。
                result = first
            else:
                result = precalc.precalculate(
                    str(_CASE), "1", state, reaction_solid_kg_s=reaction,
                    previous_drain_kg_s=max(drain, 0.0), previous_fuel_return_kg_s=fuel_return,
                    previous_release=release, iteration=k, solver=SolverConfig(),
                )
            frozen = result.frozen
            problem = Problem(frozen)
            local: dict = {}
            inner = run_inner(problem, pack_state(warm, frozen), variant, limit_s, 150, local)
            stage = inner["stages"][-1]
            x = local.get("x_final")
            fallback = None
            if stage["status"] == "DAMPING_FAILED" and "o2proj" not in variant and x is not None:
                # 只在阻尼失败时，从停下的点改用 O₂ 投影阻尼续解（痕量 O₂ 的幂变量越零）。
                again: dict = {}
                retry = run_inner(problem, x, variant | {"o2proj"}, limit_s, 150, again)
                fallback = {"from_merit": stage["final_merit"], "status": retry["stages"][-1]["status"],
                            "iterations": len(retry["stages"][-1]["iterations"])}
                stage, x = retry["stages"][-1], again.get("x_final", x)
            entry = {
                "round": k, "inner": stage["status"], "inner_reason": stage["reason"],
                "iterations": len(stage["iterations"]), "final_merit": stage["final_merit"],
                "final_check_ok": (stage.get("final_check") or {}).get("ok"),
                "drain_substituted": substituted, "previous_drain_kg_s": drain,
                "frozen_drain_ash_kg_s": float(frozen.drain_ash_kg_s), "inner_s": inner["elapsed_s"],
                "o2proj_fallback": fallback,
                "newton_tail": [{key: it.get(key) for key in ("iteration", "merit", "accepted_i", "rejections", "census", "largest")}
                                for it in stage["iterations"][-4:]],
            }
            if x is not None:
                np.save(out / f"x-round{k}.npy", x)
            if x is None or stage["status"] != "CONVERGED":
                if x is not None:
                    entry.update({key: value for key, value in balance(frozen, x).items() if key in (
                        "bed_char_kg", "bed_ash_kg", "carried_ash_kg_s", "ash_generated_kg_s", "a5_balance_kg_s",
                        "bed_temperature_c")})
                log["rounds"].append(entry)
                log["stop"] = "内层未收敛"
                break
            solved = unpack_state(x, frozen)
            entry.update({key: value for key, value in balance(frozen, x).items() if key in (
                "bed_char_kg", "bed_ash_kg", "carried_ash_kg_s", "ash_generated_kg_s", "a5_balance_kg_s",
                "bed_temperature_c")})
            entry["freeboard_temperature_c"] = (solved.freeboard_temperature_k - 273.15).tolist()
            entry["max_bed_temperature_c"] = float(np.max(solved.bed_temperature_k) - 273.15)
            log["rounds"].append(entry)
            (out / "coupled.json").write_text(json.dumps(log, ensure_ascii=False, indent=1), encoding="utf-8")
            assembly = evaluate_residual(x, frozen)
            new_reaction = -assembly.char_consumption_bed_kg_s.sum(axis=1)
            next_state = blend(state, solved, omega)
            warm = solved
            reaction = np.asarray(reaction) + omega * (np.asarray(new_reaction) - np.asarray(reaction))
            drain = float(frozen.drain_ash_kg_s + frozen.drain_char_kg_s.sum())
            fuel_return = float(frozen.fuel_down_free_kg_s[-1].sum())
            release = result.solid_context.released
            if previous is not None:
                temperature = np.concatenate([solved.bed_temperature_k, solved.freeboard_temperature_k])
                change_t = float(np.max(np.abs(temperature - previous["t"])))
                change_c = abs(entry["bed_char_kg"] - previous["c"]) / max(entry["bed_char_kg"], 1e-12)
                change_d = abs(entry["a5_balance_kg_s"] - previous["d"]) / max(entry["ash_generated_kg_s"], 1e-12)
                entry["change"] = {"temperature_k": change_t, "char_relative": change_c, "a5_relative_to_g": change_d}
                if change_t < 1.0 and change_c < 1e-3 and change_d < 1e-3:
                    log["stop"] = "相邻两轮变化小于判据"
                    break
            previous = {"t": np.concatenate([solved.bed_temperature_k, solved.freeboard_temperature_k]),
                        "c": entry["bed_char_kg"], "d": entry["a5_balance_kg_s"]}
            state = next_state
        else:
            log["stop"] = "做满轮数"
    except Exception as exc:
        log["stop"] = f"{type(exc).__name__}: {str(exc)[:600]}"
    finally:
        restore()
    log["elapsed_s"] = time.monotonic() - started
    last = log["rounds"][-1] if log["rounds"] else {}
    log["fixed_point_admissible"] = bool(
        log.get("stop") == "相邻两轮变化小于判据" and last.get("inner") == "CONVERGED" and last.get("final_check_ok")
        and last.get("a5_balance_kg_s", -1.0) >= 0.0 and not last.get("drain_substituted")
    )
    return log


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rounds", type=int, default=12)
    parser.add_argument("--omega", type=float, default=0.5)
    parser.add_argument("--variant", default="gaslog,o2pow,relax,corellaX")
    parser.add_argument("--limit-s", type=float, default=3600.0)
    parser.add_argument("--out", default=str(_OUT))
    args = parser.parse_args(argv)
    variant = set(args.variant.split(","))
    if "corellaX" in variant:
        import hamel.kinetics.homogeneous as homogeneous
        homogeneous._CORELLA_HIGH_K = float("inf")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    log = run(args.rounds, args.omega, variant, args.limit_s, out)
    (out / "coupled.json").write_text(json.dumps(log, ensure_ascii=False, indent=1), encoding="utf-8")
    print(log.get("stop"), log["fixed_point_admissible"], round(log["elapsed_s"]))


if __name__ == "__main__":
    main()
