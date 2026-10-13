#!/usr/bin/env python3
"""定位 bed2 R_H2O 翻转发生在 outer 的哪一段窗口（handoff §5.13 下一刀）。

窗口：refresh 前/后、preinner 前/后、inner NR 后、postinner 后。
在 outer abgleich 入口注入钩子（避免 import 绑定导致 monkeypatch 失效）。

用法：
  python3 scripts/audit_bed2_flip_newton_window.py --max-global-iter 4
"""

from __future__ import annotations

import argparse
from typing import Any, Callable

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.species import GAS_SPECIES_INDEX
from src.workflow.steps import outer_abgleich_step as oas

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
    "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
    "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.3,
    "nr_outer_preinner_upper_dense_absorb_mode_thesis": "phase_share",
    "nr_outer_preinner_upper_dense_absorb_post_gate_thesis": True,
}

MAJORS = ("H2O", "H2", "O2", "CO", "CO2", "N2", "CH4")


def _bed2_R(cells: list[Any]) -> dict[str, float]:
    c = cells[2]
    c.residuals()
    out: dict[str, float] = {"T": float(c.T)}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    return out


def _bed_vec(cell: Any) -> dict[str, float]:
    out: dict[str, float] = {"T": float(cell.T)}
    for sp in MAJORS:
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"{sp}_tot"] = float(cell.N_d[j] + cell.N_b[j])
        out[f"{sp}_in"] = float(cell.N_d_in[j] + cell.N_b_in[j])
    return out


def _all_beds(cells: list[Any], n: int = 10) -> dict[str, dict[str, float]]:
    return {f"bed{i}": _bed_vec(cells[i]) for i in range(min(n, len(cells)))}


def _top_deltas(
    a: dict[str, dict[str, float]], b: dict[str, dict[str, float]], k: int = 20
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for bed, vb in b.items():
        va = a.get(bed) or {}
        for key, val in vb.items():
            old = float(va.get(key, 0.0))
            new = float(val)
            d = new - old
            if abs(d) < 1e-12:
                continue
            rows.append(
                {
                    "key": f"{bed}.{key}",
                    "old": old,
                    "new": new,
                    "d": d,
                    "score": abs(d) if key != "T" else abs(d) * 0.01,
                }
            )
    rows.sort(key=lambda r: float(r["score"]), reverse=True)
    return rows[:k]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=4)
    args = ap.parse_args()

    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=dict(ISOMORPH))
    run_phase2_init(reactor)
    cells = reactor.cells

    events: list[dict[str, Any]] = []
    snap_ref: dict[str, Any] = {"beds": None, "label": ""}
    flip_info: dict[str, Any] | None = None
    ctr = {"refresh": 0, "preinner": 0, "inner": 0, "post": 0}

    def _note(tag: str, i: int, before: dict[str, float], after: dict[str, float], **extra: Any) -> None:
        nonlocal flip_info
        events.append(
            {
                "tag": tag,
                "i": i,
                "R_before": before,
                "R_after": after,
                **extra,
            }
        )
        print(
            f"[{tag}#{i}] R_H2O {before['R_H2O']:+.3f}→{after['R_H2O']:+.3f}"
            + (f" {extra}" if extra else "")
        )
        # 首次跨越：负/小 → 大正
        if (
            flip_info is None
            and before["R_H2O"] < 5.0
            and after["R_H2O"] > 20.0
        ):
            beds_now = _all_beds(cells)
            flip_info = {
                "window": f"{tag}#{i}",
                "R_before": before,
                "R_after": after,
                "top_deltas_vs_ref": (
                    None
                    if snap_ref["beds"] is None
                    else _top_deltas(snap_ref["beds"], beds_now)
                ),
                "ref_label": snap_ref["label"],
            }
            print(f"  *** FLIP at {tag}#{i} (ref={snap_ref['label']}) ***")
            if flip_info["top_deltas_vs_ref"]:
                for row in flip_info["top_deltas_vs_ref"][:12]:
                    print(
                        f"    {row['key']}: {row['old']:+.4g}→{row['new']:+.4g} "
                        f"(d={row['d']:+.4g})"
                    )

    orig = oas.run_outer_abgleich_for_global_nr

    def wrapped_abgleich(**kwargs: Any) -> Any:
        cells_kw = kwargs["cells"]
        refresh_fn = kwargs["refresh_fn"]
        pre_inner_fn = kwargs.get("pre_inner_fn")
        post_fn = kwargs.get("post_inner_holdup_fn")
        solve_inner_fn = kwargs["solve_inner_fn"]

        def wrap_refresh(force: bool) -> None:
            ctr["refresh"] += 1
            i = ctr["refresh"]
            before = _bed2_R(cells_kw)
            refresh_fn(force)
            after = _bed2_R(cells_kw)
            _note("refresh", i, before, after)

        def wrap_preinner(**pkw: Any) -> dict[str, bool | float] | None:
            ctr["preinner"] += 1
            i = ctr["preinner"]
            before = _bed2_R(cells_kw)
            assert pre_inner_fn is not None
            diag = pre_inner_fn(**pkw) if pkw else pre_inner_fn()
            after = _bed2_R(cells_kw)
            _note(
                "preinner",
                i,
                before,
                after,
                local2=bool((diag or {}).get("upper_co_local2x2_accepted", False)),
                absorb=bool((diag or {}).get("upper_dense_absorb_accepted", False)),
            )
            # 以未翻转的 preinner 出口为差分参考
            if after["R_H2O"] < 5.0:
                snap_ref["beds"] = _all_beds(cells_kw)
                snap_ref["label"] = f"after_preinner#{i}"
            return diag

        def wrap_inner(
            max_iter: int, lambda_seed: float, jac_mode: str, jac_lag: int
        ) -> dict:
            ctr["inner"] += 1
            i = ctr["inner"]
            before = _bed2_R(cells_kw)
            beds_before = _all_beds(cells_kw)
            result = solve_inner_fn(max_iter, lambda_seed, jac_mode, jac_lag)
            after = _bed2_R(cells_kw)
            extra: dict[str, Any] = {
                "n_newton": int(result.get("n_newton_iters_attempted", result.get("n_iter", -1))),
                "rms": float(result.get("rms_scaled_final", float("nan"))),
            }
            if before["R_H2O"] < 5.0 and after["R_H2O"] > 20.0:
                extra["top_deltas"] = _top_deltas(beds_before, _all_beds(cells_kw))
            _note("inner", i, before, after, **extra)
            return result

        def wrap_post() -> bool:
            ctr["post"] += 1
            i = ctr["post"]
            before = _bed2_R(cells_kw)
            assert post_fn is not None
            out = bool(post_fn())
            after = _bed2_R(cells_kw)
            _note("postinner", i, before, after)
            return out

        kwargs = dict(kwargs)
        kwargs["refresh_fn"] = wrap_refresh
        kwargs["solve_inner_fn"] = wrap_inner
        if pre_inner_fn is not None:
            kwargs["pre_inner_fn"] = wrap_preinner
        if post_fn is not None:
            kwargs["post_inner_holdup_fn"] = wrap_post
        return orig(**kwargs)

    oas.run_outer_abgleich_for_global_nr = wrapped_abgleich  # type: ignore[assignment]
    try:
        result, wall = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter))
        )
    finally:
        oas.run_outer_abgleich_for_global_nr = orig  # type: ignore[assignment]

    # reactor 内是 from ... import run_outer_abgleich_for_global_nr —— 需同步补丁 reactor 绑定
    # 上面若未生效，events 会很少；下面再试一次经 reactor 模块补丁重跑太贵。
    # 改：补丁 reactor 模块里的同名引用。
    exit_co = float((result.get("bed_exit_gas") or {}).get("CO", float("nan")))
    out = {
        "note": "bed2 flip window via outer_abgleich_step hook",
        "exit_CO": exit_co,
        "final_R": _bed2_R(cells),
        "wall_s": float(wall),
        "n_events": len(events),
        "events": events,
        "flip_info": flip_info,
    }
    path = data_path("phase2_bed2_flip_newton_window_o4.json")
    write_json(path, out)
    print(
        f"wrote {path} n_events={len(events)} CO={exit_co:.4f} "
        f"flip={None if flip_info is None else flip_info.get('window')}"
    )


if __name__ == "__main__":
    main()
