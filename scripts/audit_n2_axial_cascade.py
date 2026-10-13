#!/usr/bin/env python3
"""全床 N2 轴向 cascade 通盘审计（handoff §5.26）。

假设：局部放出 bed k 的过填 holdup → 出流下降 → bed k+1 入口下降而 holdup
未动 → gap 上移。验证全床剖面 + 链式放出 A/B。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_n2_axial_cascade.py
"""

from __future__ import annotations

from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from audit_bed3_n2_on_lsq_n2_ab import LSQ2, STACK
from audit_bed9_holdup_closure_ab import _metrics
from audit_rate_gate_rms_stack_ab import ISOMORPH, PS

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import oneshot_bed_dense_species_sum_absorb
from src.core.species import GAS_SPECIES_INDEX, N_GAS


def _n2_row(cell: Any, idx: int) -> dict[str, float]:
    cell.residuals()
    j = int(GAS_SPECIES_INDEX["N2"])
    n_d = float(cell.N_d[j])
    n_b = float(cell.N_b[j])
    gas_d = float(np.sum(cell.N_d))
    gas_b = float(np.sum(cell.N_b))
    in_d = float(cell.N_zu_d[j] + cell.N_rez_d[j] + cell.N_d_in[j])
    in_b = float(cell.N_zu_b[j] + cell.N_rez_b[j] + cell.N_b_in[j])
    r_d = float(cell.R_gas_d[j])
    r_b = float(cell.R_gas_b[j])
    n_ex = float(cell.N_ex[j])
    n_hold = n_d + n_b
    n_in = in_d + in_b
    res_sum = n_in + r_d + r_b - n_hold
    # 两相出流 = holdup（稳态 CSTR 式未知量即出口流量）
    n_out = n_hold
    return {
        "bed": float(idx),
        "T": float(cell.T),
        "N_hold": n_hold,
        "N_in": n_in,
        "N_out": n_out,
        "gap": n_hold - n_in,  # = -res_sum if R=0
        "res_sum": res_sum,
        "R": r_d + r_b,
        "N_ex": n_ex,
        "dy": n_d / max(gas_d, 1e-30) - n_b / max(gas_b, 1e-30),
        "y_d": n_d / max(gas_d, 1e-30),
        "y_b": n_b / max(gas_b, 1e-30),
        "in_d": in_d,
        "in_b": in_b,
        "N_d": n_d,
        "N_b": n_b,
    }


def _profile(reactor: Any) -> list[dict[str, float]]:
    reactor._apply_all_bc_for_nr()
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    rows = [_n2_row(c, i) for i, c in enumerate(beds)]
    # 轴向接口：本格 out vs 下一格 in（总量）
    for i in range(len(rows) - 1):
        rows[i]["out_minus_next_in"] = float(rows[i]["N_out"] - rows[i + 1]["N_in"])
    rows[-1]["out_minus_next_in"] = float("nan")
    return rows


def _print_profile(label: str, rows: list[dict[str, float]]) -> None:
    print(f"--- {label} ---", flush=True)
    print(
        f"{'bed':>3} {'N_in':>8} {'N_hold':>8} {'gap':>8} {'dy':>8} "
        f"{'N_ex':>8} {'Δout→in+1':>10}",
        flush=True,
    )
    for r in rows:
        print(
            f"{int(r['bed']):3d} {r['N_in']:8.3f} {r['N_hold']:8.3f} {r['gap']:+8.3f} "
            f"{r['dy']:+8.4f} {r['N_ex']:+8.3f} "
            f"{r['out_minus_next_in']:+10.3f}"
            if np.isfinite(r["out_minus_next_in"])
            else f"{int(r['bed']):3d} {r['N_in']:8.3f} {r['N_hold']:8.3f} {r['gap']:+8.3f} "
            f"{r['dy']:+8.4f} {r['N_ex']:+8.3f} {'n/a':>10}",
            flush=True,
        )
    gaps = [r["gap"] for r in rows]
    print(
        f"  gap: max@{int(np.argmax(np.abs(gaps)))}={gaps[int(np.argmax(np.abs(gaps)))]:+.3f} "
        f"sum={sum(gaps):+.3f} n_overfill={sum(1 for g in gaps if g > 0.5)}",
        flush=True,
    )


def _posthoc_sequential_bleed(reactor: Any, beds: tuple[int, ...]) -> dict[str, Any]:
    """自下而上逐格双向放出，记录每步 gap 剖面峰值。"""
    snaps = []
    for c in reactor.cells:
        snaps.append((c.N_d.copy(), c.N_b.copy()))

    def _restore() -> None:
        for c, (nd, nb) in zip(reactor.cells, snaps):
            c.N_d[:] = nd
            c.N_b[:] = nb
            c._thermo_cache_valid = False
        reactor._apply_all_bc_for_nr()

    traj: list[dict[str, Any]] = []
    _restore()
    base = _profile(reactor)
    traj.append({"step": "base", "profile": base, "peak_bed": int(np.argmax([abs(r["gap"]) for r in base]))})

    _restore()
    for b in beds:
        cell = reactor.cells[b]
        cell.residuals()
        added = oneshot_bed_dense_species_sum_absorb(
            cell, species=("N2",), mode="phase_share", bidirectional=True, fraction=1.0
        )
        reactor._apply_all_bc_for_nr()
        prof = _profile(reactor)
        gaps = [r["gap"] for r in prof]
        peak = int(np.argmax([abs(g) for g in gaps]))
        traj.append(
            {
                "step": f"bleed_bed{b}",
                "added_N2": float(added.get("N2", 0.0)),
                "profile": prof,
                "peak_bed": peak,
                "peak_gap": float(gaps[peak]),
                "gap_sum": float(sum(gaps)),
            }
        )
        print(
            f"  [seq bleed bed{b}] added={added.get('N2', 0):+.2f} "
            f"peak@bed{peak} gap={gaps[peak]:+.3f} gap_sum={sum(gaps):+.3f}",
            flush=True,
        )
    _restore()
    return {"trajectory": traj}


def _sum_cfg(beds: tuple[int, ...]) -> dict[str, object]:
    return {
        "nr_outer_preinner_targeted_species_sum_beds_thesis": beds,
        "nr_outer_preinner_targeted_species_sum_species_thesis": ("N2",),
        "nr_outer_preinner_targeted_species_sum_fraction_thesis": 1.0,
        "nr_outer_preinner_targeted_species_sum_mode_thesis": "phase_share",
        "nr_outer_preinner_targeted_species_sum_bidirectional_thesis": True,
        "nr_outer_preinner_targeted_species_sum_after_y_co_thesis": True,
        "nr_outer_preinner_targeted_species_sum_y_co_threshold_thesis": 0.14,
    }


def _run(label: str, overrides: dict[str, object], max_iter: int = 8) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    met = _metrics(reactor, result)
    prof = _profile(reactor)
    _print_profile(label, prof)
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"E={met['groups']['energy_rms']:.4f} wall={wall_s:.0f}s | "
        + ", ".join(
            f"c{t['cell']}:{t.get('species') or t.get('kind')}={t['F_hat']:+.3f}"
            for t in met["topk"][:3]
        ),
        flush=True,
    )
    return {
        "label": label,
        "wall_s": float(wall_s),
        "metrics": met,
        "profile": prof,
        "overrides": dict(overrides),
        "_reactor": reactor,
    }


def main() -> None:
    arms: list[dict[str, Any]] = []
    seq: dict[str, Any] = {}

    base = _run("lsq_n2_base", LSQ2)
    reactor = base.pop("_reactor")
    print("sequential bottom-up bleed on LSQ n=2 end-state:", flush=True)
    seq = _posthoc_sequential_bleed(reactor, tuple(range(10)))
    arms.append(base)

    for label, beds in (
        ("sum_3_5", (3, 4, 5)),
        ("sum_3_9", (3, 4, 5, 6, 7, 8, 9)),
        ("sum_1_9", tuple(range(1, 10))),
    ):
        arm = _run(label, {**LSQ2, **_sum_cfg(beds)})
        arm.pop("_reactor", None)
        arms.append(arm)

    out = data_path("phase2_n2_axial_cascade_o8.json")
    write_json(
        out,
        {
            "note": "full-bed N2 axial cascade + chained bleed A/B; §5.26",
            "arms": arms,
            "posthoc_sequential_bleed": seq,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
