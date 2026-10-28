"""预计算耗时的分期统计（只报告，不进入求解）。

把第一次外层预计算的各阶段函数包一层计时，给出轨迹、床区剖面、热解与其余部分
的墙钟占比。ADR 0047 只改轨迹，因此其余几项应当基本不变。
"""
from __future__ import annotations

import json
import sys
import time

from hamel.cases.config import SolverConfig
from hamel.solver import precalc

# 阶段名 → precalc 里的函数名。同一阶段的函数可以嵌套，计时按最外层计入。
PHASES = {
    "trajectory": ("_segregate",),
    "bed": ("_adaptive_bed", "_exchange_diffusivity"),
    "pyrolysis": ("_release",),
    "startup": ("_equilibrium", "_flows", "_bed_gas", "_feed_share", "_initial_inventory", "_drain"),
}


def main(out: str) -> None:
    totals: dict[str, float] = {name: 0.0 for name in PHASES}
    depth: dict[str, int] = {name: 0 for name in PHASES}

    for phase, names in PHASES.items():
        for name in names:
            original = getattr(precalc, name)

            def wrapped(*args, __original=original, __phase=phase, **kwargs):
                # 阶段内嵌套（例如 _adaptive_bed 内的 _density）只在最外层计时一次。
                if depth[__phase]:
                    return __original(*args, **kwargs)
                depth[__phase] += 1
                started = time.perf_counter()
                try:
                    return __original(*args, **kwargs)
                finally:
                    totals[__phase] += time.perf_counter() - started
                    depth[__phase] -= 1

            setattr(precalc, name, wrapped)

    started = time.perf_counter()
    prepared = precalc.precalculate(
        "data/cases/vtt_pressurized.json", "1", None, adiabatic=False, iteration=1, solver=SolverConfig()
    )
    elapsed = time.perf_counter() - started
    attributed = sum(totals.values())
    payload = {
        "total_s": elapsed,
        "phases": {**dict(totals), "other": elapsed - attributed},
        "t_ad_k": float(prepared.t_ad_k),
        "rho_p_kg_m3": float(prepared.frozen.rho_p_kg_m3),
    }
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=1)
    for name, value in payload["phases"].items():
        print(f"{name:12s} {value:8.1f} s  {value / elapsed * 100:5.1f} %")
    print(f"{'total':12s} {elapsed:8.1f} s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "timing.json")
