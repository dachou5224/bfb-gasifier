#!/usr/bin/env python3
"""床内 O₂ 轴向：气泡/密相拆分 + 源(N_zu)/汇(R12/R5/R6/R1) 对照 Hamel「第二格耗尽」。

Usage
-----
    python3 scripts/audit_bed_o2_axial_phase_plateau_lu.py
    python3 scripts/audit_bed_o2_axial_phase_plateau_lu.py --max-global-iter 8 --seed-T0 1150
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import numpy as np

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from _lu_harness import (  # noqa: E402
    add_json_output_arg,
    add_max_global_iter_arg,
    build_phase2_reactor,
    phase2_solve_kwargs,
    timed_phase2_solve,
    write_json,
)
import src.workflow.steps.init_precalc_step as ips  # noqa: E402
from audit_r1_o2_phase_bottleneck_lu import _bundle_for_cell  # noqa: E402
from src.core.species import GAS_SPECIES_INDEX as idx  # noqa: E402


def _run_solve(*, profile: str, seed_T0: float | None, max_iter: int) -> tuple[Any, dict, float]:
    _orig = ips._install_temperature_fence_from_vorabrechnung

    def _install(reactor: Any, T_est: Any) -> None:
        _orig(reactor, T_est)
        if seed_T0 is None:
            return
        beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
        if not beds:
            return
        b0 = beds[0]
        b0._nr_temperature_min_K = float(max(700.0, seed_T0 - 200.0))
        b0._nr_temperature_max_K = float(min(2500.0, max(seed_T0 + 400.0, 1400.0)))
        b0.T = float(seed_T0)
        b0._vorab_temperature_estimate_K = float(seed_T0)

    overrides: dict[str, object] = {
        "vorab_vm_devolatilization_zone_cells_thesis": 3,
        "vorab_vm_devolatilization_zone_profile_thesis": profile,
        "enable_r12": True,
    }
    if seed_T0 is not None:
        overrides["vorab_nr_temperature_seed_bed0_cold_anchor_thesis"] = False
        overrides["vorab_nr_temperature_seed_bottom_K_thesis"] = float(seed_T0)

    ips._install_temperature_fence_from_vorabrechnung = _install
    try:
        reactor = build_phase2_reactor(
            ladder_level="p2",
            builder_kwargs={"enable_vorab_march_staged": True},
            cfg_overrides=overrides,
        )
        result, wall = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=max_iter)
        )
        return reactor, result, wall
    finally:
        ips._install_temperature_fence_from_vorabrechnung = _orig


def _audit_beds(reactor: Any) -> list[dict[str, float | int]]:
    j = idx["O2"]
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    rows: list[dict[str, float | int]] = []
    for i, cell in enumerate(beds):
        cell.calc_hydrodynamics()
        cell.calc_exchange()
        bundle, thermo, _areas, _dg = _bundle_for_cell(cell)
        y_b = np.asarray(thermo["y_b"], dtype=np.float64)
        y_d = np.asarray(thermo["y_d"], dtype=np.float64)
        ymix = float(cell._mole_fractions("combined")[j])

        def g(name: str, default: float = 0.0) -> float:
            return float(getattr(bundle, name, default) or default)

        r12b, r12d = g("extent_r12_bubble"), g("extent_r12_dense")
        r12 = r12b + r12d
        r5b, r5d = g("extent_r5_bubble"), g("extent_r5_dense")
        r6b, r6d = g("extent_r6_bubble"), g("extent_r6_dense")
        r1o2 = g("extent_r1_o2")
        sink = 0.5 * r12 + 0.5 * (r5b + r5d) + 1.5 * (r6b + r6d) + r1o2
        nb, nd = float(cell.N_b[j]), float(cell.N_d[j])
        zu_b, zu_d = float(cell.N_zu_b[j]), float(cell.N_zu_d[j])
        rows.append(
            {
                "i": i,
                "T_K": float(cell.T),
                "y_O2_b": float(y_b[j]),
                "y_O2_d": float(y_d[j]),
                "y_O2_mix": ymix,
                "N_b_O2": nb,
                "N_d_O2": nd,
                "N_tot_O2": nb + nd,
                "N_zu_b_O2": zu_b,
                "N_zu_d_O2": zu_d,
                "N_zu_O2": zu_b + zu_d,
                "N_ex_O2": float(cell.N_ex[j]),
                "K_bd": float(cell.K_bd),
                "eps_b": float(cell.eps_b),
                "R12": r12,
                "R12b": r12b,
                "R12d": r12d,
                "R5b": r5b,
                "R5d": r5d,
                "R6": r6b + r6d,
                "R1_O2": r1o2,
                "O2_sink_est": sink,
                "lf_b": g("limit_factor_o2_bubble", 1.0),
                "lf_d": g("limit_factor_o2_dense", 1.0),
                "supply_b": g("o2_supply_bubble"),
                "supply_d": g("o2_supply_dense"),
                "slack_b": g("o2_slack_bubble"),
                "slack_d": g("o2_slack_dense"),
            }
        )
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    add_max_global_iter_arg(p, default=8)
    add_json_output_arg(p, "phase2_bed_o2_axial_phase_plateau_o8.json")
    p.add_argument("--profile", default="front_loaded", choices=("front_loaded", "linear"))
    p.add_argument("--seed-T0", type=float, default=1150.0)
    args = p.parse_args()
    seed = None if float(args.seed_T0) < 0 else float(args.seed_T0)

    print(
        f"=== O2 axial phase/source/sink  profile={args.profile} "
        f"seed={seed} iter={args.max_global_iter} ===",
        flush=True,
    )
    reactor, result, wall = _run_solve(
        profile=str(args.profile), seed_T0=seed, max_iter=int(args.max_global_iter)
    )
    rows = _audit_beds(reactor)

    print(
        f"wall={wall:.1f}s rms={result.get('rms_scaled_final')} "
        f"conv={result.get('converged_fully')}",
        flush=True,
    )
    print(
        f"{'i':>2} {'T':>5} {'yOb':>7} {'yOd':>7} {'ymix':>7} "
        f"{'Nb':>6} {'Nd':>6} {'zu':>6} {'Nex':>7} | "
        f"{'R12':>5} {'R5b':>5} {'R5d':>5} {'sink':>6} {'lfd':>5}",
        flush=True,
    )
    for row in rows:
        print(
            f"{row['i']:2d} {row['T_K']:5.0f} {row['y_O2_b']:7.4f} {row['y_O2_d']:7.4f} "
            f"{row['y_O2_mix']:7.4f} {row['N_b_O2']:6.2f} {row['N_d_O2']:6.2f} "
            f"{row['N_zu_O2']:6.2f} {row['N_ex_O2']:7.2f} | "
            f"{row['R12']:5.2f} {row['R5b']:5.2f} {row['R5d']:5.2f} "
            f"{row['O2_sink_est']:6.2f} {row['lf_d']:5.2f}",
            flush=True,
        )

    ymix = [float(r["y_O2_mix"]) for r in rows]
    ntot = [float(r["N_tot_O2"]) for r in rows]
    print("\n--- vs Hamel (O2 verzehrt in 2. Zelle) ---", flush=True)
    print(
        f"ymix bed1/bed2 = {ymix[1]:.5f} / {ymix[2]:.5f}  "
        f"(target ~0 by cell 2)",
        flush=True,
    )
    if len(ntot) > 3 and ntot[3] > ntot[2] * 1.2:
        print(
            f"REBOUND N_tot_O2 bed2={ntot[2]:.4f} -> bed3={ntot[3]:.4f}",
            flush=True,
        )
    payload = {
        "note": "Bed O2 axial phase/source/sink vs Hamel 2nd-cell depletion",
        "profile": args.profile,
        "seed_T0": seed,
        "max_iter": int(args.max_global_iter),
        "wall_s": wall,
        "rms_scaled_final": result.get("rms_scaled_final"),
        "converged_fully": result.get("converged_fully"),
        "O2_feed": float(reactor.config.O2_feed),
        "sec_O2": float(getattr(reactor.config, "freeboard_secondary_O2_mol_s", 0.0) or 0.0),
        "staged_primary": bool(
            getattr(reactor.config, "vorab_staged_primary_gas_injection_thesis", False)
        ),
        "r5_scale": float(getattr(reactor.config, "r5_scale", 1.0) or 1.0),
        "exit_dry": {k: float(v) for k, v in (result.get("exit_gas_dry") or {}).items()},
        "rows": rows,
        "hamel_quote": "Der Sauerstoff ist bereits in der zweiten Zelle verzehrt (Bild 7.8)",
    }
    out = write_json(args.output, payload)
    print(f"\nWrote {out}", flush=True)


if __name__ == "__main__":
    main()
