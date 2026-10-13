#!/usr/bin/env python3
"""Phase2 LU：Vorab march / staged primary / layer weights 对比（合并入口）。

Usage
-----
    PYTHONUNBUFFERED=1 python3 scripts/compare_vorab_march_lu.py --scenario staged
    PYTHONUNBUFFERED=1 python3 scripts/compare_vorab_march_lu.py --scenario weights
    PYTHONUNBUFFERED=1 python3 scripts/compare_vorab_march_lu.py --scenario both
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from typing import Literal

from _lu_harness import (
    REF_CO,
    REF_CO2,
    add_json_output_arg,
    add_max_global_iter_arg,
    apply_baseline_kinetics,
    bed_cell_indices,
    build_phase2_config,
    co_co2_score,
    ensure_repo_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from src.core.cell_thermo_helpers import gas_mole_fractions
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx

ensure_repo_path()

ScenarioKind = Literal["staged", "weights", "both"]
DRY_SPECIES = ("CO", "CO2", "H2", "CH4", "N2", "O2")


@dataclass
class StagedScenarioResult:
    label: str
    vorab_bed_temperature_march: bool
    vorab_staged_primary_gas: bool
    bed0_T_init_K: float
    bed0_T_marched_budget_K: float
    bed0_T_solved_K: float
    bed_T_init_K: list[float]
    bed_T_solved_K: list[float]
    staged_bed_cells_with_O2_zu: int
    converged_fully: bool
    rms_scaled_final: float
    co_exit: float
    co2_exit: float
    bed_co: float
    bed_co2: float
    score_co_co2: float
    wall_time_s: float


@dataclass
class WeightsScenarioSpec:
    label: str
    layer_weights: tuple[float, ...] | None
    enable_march: bool


@dataclass
class BedRow:
    bed_i: int
    xi: float
    T_K: float
    O2_zu_pct: float
    y_dry_CO: float
    y_dry_CO2: float
    CO_share_pct: float
    staged_injection: bool


@dataclass
class WeightsScenarioResult:
    label: str
    layer_weights: list[float] | None
    vorab_bed_temperature_march: bool
    bed0_T_init_K: float
    bed0_T_marched_budget_K: float
    bed0_T_solved_K: float
    bed1_T_solved_K: float
    bed2_T_solved_K: float
    converged_fully: bool
    rms_scaled_final: float
    co_exit: float
    co2_exit: float
    score_co_co2: float
    bed_rows: list[BedRow]
    wall_time_s: float


def _count_staged_o2_zu(reactor: Reactor) -> int:
    n = 0
    for cell in reactor.cells:
        if float(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]]) > 1.0e-12:
            n += 1
    return n


def _dry_y(N_b, N_d) -> dict[str, float]:
    y = gas_mole_fractions(N_b, N_d, "mixed")
    total = float(sum(y[idx[s]] for s in DRY_SPECIES))
    if total < 1e-15:
        return {s: 0.0 for s in DRY_SPECIES}
    return {s: float(y[idx[s]] / total) for s in DRY_SPECIES}


def run_staged_scenario(*, label: str, enable_march_staged: bool, max_global_iter: int) -> StagedScenarioResult:
    cfg = build_phase2_config(enable_vorab_march_staged=enable_march_staged)
    apply_baseline_kinetics(cfg)

    reactor = Reactor(cfg)
    init = run_phase2_init(reactor)
    bed_idx = bed_cell_indices(reactor)
    bed_T_init = [float(init.T_est[i]) for i in bed_idx]
    marched = getattr(reactor, "_vorab_marched_t_budget_K", None)
    bed0_marched = float(marched[bed_idx[0]]) if marched is not None and bed_idx else float("nan")

    out, wall_s = timed_phase2_solve(
        reactor,
        **phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False),
    )

    bed_T_solved = [float(reactor.cells[i].T) for i in bed_idx]
    y = out.get("exit_gas_dry") or {}
    bed = out.get("bed_exit_gas_dry") or {}
    co = float(y.get("CO", 0.0))
    co2 = float(y.get("CO2", 0.0))

    return StagedScenarioResult(
        label=label,
        vorab_bed_temperature_march=bool(cfg.vorab_bed_temperature_march_thesis),
        vorab_staged_primary_gas=bool(cfg.vorab_staged_primary_gas_injection_thesis),
        bed0_T_init_K=float(bed_T_init[0]) if bed_T_init else float("nan"),
        bed0_T_marched_budget_K=bed0_marched,
        bed0_T_solved_K=float(reactor.cells[bed_idx[0]].T) if bed_idx else float("nan"),
        bed_T_init_K=bed_T_init,
        bed_T_solved_K=bed_T_solved,
        staged_bed_cells_with_O2_zu=_count_staged_o2_zu(reactor),
        converged_fully=bool(out.get("converged_fully", False)),
        rms_scaled_final=float(out.get("rms_scaled_final", float("nan"))),
        co_exit=co,
        co2_exit=co2,
        bed_co=float(bed.get("CO", 0.0)),
        bed_co2=float(bed.get("CO2", 0.0)),
        score_co_co2=co_co2_score(co, co2),
        wall_time_s=float(wall_s),
    )


def run_weights_scenario(spec: WeightsScenarioSpec, *, max_global_iter: int) -> WeightsScenarioResult:
    cfg = build_phase2_config(
        enable_vorab_march_staged=True,
        enable_vorab_cold_budget_march=bool(spec.enable_march),
    )
    apply_baseline_kinetics(cfg)
    cfg.vorab_staged_primary_gas_layer_weights_thesis = spec.layer_weights

    reactor = Reactor(cfg)
    init = run_phase2_init(reactor)
    bed_idx = bed_cell_indices(reactor)
    marched = getattr(reactor, "_vorab_marched_t_budget_K", None)
    bed0_marched = float(marched[bed_idx[0]]) if marched is not None and bed_idx else float("nan")
    bed0_init = float(init.T_est[bed_idx[0]]) if bed_idx else float("nan")

    out, wall_s = timed_phase2_solve(
        reactor,
        **phase2_solve_kwargs(max_global_iter=int(max_global_iter), verbose=False),
    )

    o2_feed = float(cfg.O2_feed)
    h_bed = float(cfg.H_bed)
    rows: list[BedRow] = []
    for bi in bed_idx[:6]:
        cell = reactor.cells[bi]
        y = _dry_y(cell.N_b, cell.N_d)
        co_sum = y["CO"] + y["CO2"]
        o2_zu = float(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]])
        rows.append(
            BedRow(
                bed_i=int(bi),
                xi=float(cell.geo.h_center) / h_bed if h_bed > 0 else 0.0,
                T_K=float(cell.T),
                O2_zu_pct=100.0 * o2_zu / o2_feed if o2_feed > 0 else 0.0,
                y_dry_CO=y["CO"],
                y_dry_CO2=y["CO2"],
                CO_share_pct=100.0 * y["CO"] / co_sum if co_sum > 1e-12 else 0.0,
                staged_injection=o2_zu > 1e-9,
            )
        )

    exit_y = out.get("exit_gas_dry") or {}
    co = float(exit_y.get("CO", 0.0))
    co2 = float(exit_y.get("CO2", 0.0))
    w_list = None if spec.layer_weights is None else [float(x) for x in spec.layer_weights]

    return WeightsScenarioResult(
        label=spec.label,
        layer_weights=w_list,
        vorab_bed_temperature_march=bool(cfg.vorab_bed_temperature_march_thesis),
        bed0_T_init_K=bed0_init,
        bed0_T_marched_budget_K=bed0_marched,
        bed0_T_solved_K=float(reactor.cells[bed_idx[0]].T) if bed_idx else float("nan"),
        bed1_T_solved_K=float(reactor.cells[bed_idx[1]].T) if len(bed_idx) > 1 else float("nan"),
        bed2_T_solved_K=float(reactor.cells[bed_idx[2]].T) if len(bed_idx) > 2 else float("nan"),
        converged_fully=bool(out.get("converged_fully", False)),
        rms_scaled_final=float(out.get("rms_scaled_final", float("nan"))),
        co_exit=co,
        co2_exit=co2,
        score_co_co2=co_co2_score(co, co2),
        bed_rows=rows,
        wall_time_s=float(wall_s),
    )


def _run_staged_block(*, max_global_iter: int) -> tuple[list[StagedScenarioResult], dict]:
    scenarios = [
        ("legacy_clip_T", False),
        ("march_staged", True),
    ]
    rows: list[StagedScenarioResult] = []
    print(
        f"Phase2 LU Vorab staged compare | r2=r4=r5=r6=r7=r8=1.0 | max_global_iter={max_global_iter}",
        flush=True,
    )
    print(f"ref CO={REF_CO} CO2={REF_CO2}", flush=True)
    print("=" * 100, flush=True)

    for label, enabled in scenarios:
        row = run_staged_scenario(
            label=label,
            enable_march_staged=enabled,
            max_global_iter=int(max_global_iter),
        )
        rows.append(row)
        print(
            f"{row.label:16s} | init bed0={row.bed0_T_init_K:.1f}K "
            f"marched={row.bed0_T_marched_budget_K:.1f}K solved={row.bed0_T_solved_K:.1f}K "
            f"| O2_zu_cells={row.staged_bed_cells_with_O2_zu} "
            f"| conv={row.converged_fully} rms={row.rms_scaled_final:.4g} "
            f"| CO={row.co_exit:.4f} CO2={row.co2_exit:.4f} score={row.score_co_co2:.3f} "
            f"| {row.wall_time_s:.0f}s",
            flush=True,
        )

    best = min(rows, key=lambda r: r.score_co_co2)
    print("=" * 100, flush=True)
    print(f"best staged score: {best.label} (CO={best.co_exit:.4f} CO2={best.co2_exit:.4f})", flush=True)
    payload = {
        "note": "baseline kinetics r2=r4=r5=r6=r7=r8=1.0; vorab_major_gibbs_x0=default",
        "ref_co": REF_CO,
        "ref_co2": REF_CO2,
        "max_global_iter": int(max_global_iter),
        "best_by_co_co2_score": best.label,
        "rows": [asdict(r) for r in rows],
    }
    return rows, payload


def _run_weights_block(*, max_global_iter: int) -> tuple[list[WeightsScenarioResult], dict]:
    scenarios = [
        WeightsScenarioSpec("baseline_decay", None, False),
        WeightsScenarioSpec("weights_50_30_20", (0.50, 0.30, 0.20), False),
        WeightsScenarioSpec("march_decay", None, True),
        WeightsScenarioSpec("march_50_30_20", (0.50, 0.30, 0.20), True),
    ]

    print(
        f"Phase2 LU staged×march weights | r2–r8=1.0 | max_global_iter={max_global_iter}",
        flush=True,
    )
    print(f"ref CO={REF_CO} CO2={REF_CO2}", flush=True)
    print("=" * 110, flush=True)

    rows: list[WeightsScenarioResult] = []
    for spec in scenarios:
        row = run_weights_scenario(spec, max_global_iter=int(max_global_iter))
        rows.append(row)
        print(
            f"{row.label:18s} march={row.vorab_bed_temperature_march} "
            f"w={row.layer_weights or 'decay'} "
            f"| init={row.bed0_T_init_K:.0f}K budget={row.bed0_T_marched_budget_K:.0f}K "
            f"solved0/1/2={row.bed0_T_solved_K:.0f}/{row.bed1_T_solved_K:.0f}/{row.bed2_T_solved_K:.0f}K "
            f"| conv={row.converged_fully} rms={row.rms_scaled_final:.4g} "
            f"CO={row.co_exit:.4f} CO2={row.co2_exit:.4f} score={row.score_co_co2:.3f} "
            f"| {row.wall_time_s:.0f}s",
            flush=True,
        )
        for br in row.bed_rows[:3]:
            print(
                f"  bed{br.bed_i} xi={br.xi:.2f} T={br.T_K:.0f}K O2in={br.O2_zu_pct:5.1f}% "
                f"yCO={br.y_dry_CO:.3f} yCO2={br.y_dry_CO2:.3f} CO%={br.CO_share_pct:.0f}",
                flush=True,
            )

    converged = [r for r in rows if r.converged_fully]
    best_score = min(rows, key=lambda r: r.score_co_co2)
    best_conv = min(converged, key=lambda r: r.score_co_co2) if converged else None

    print("=" * 110, flush=True)
    print(f"best weights score (any): {best_score.label} (score={best_score.score_co_co2:.3f})", flush=True)
    if best_conv:
        print(
            f"best weights converged: {best_conv.label} (score={best_conv.score_co_co2:.3f})",
            flush=True,
        )
    else:
        print("best weights converged: none", flush=True)

    payload = {
        "note": "decay default vs 50/30/20 × cold budget march; dual-track NR seed when march on",
        "ref_co": REF_CO,
        "ref_co2": REF_CO2,
        "max_global_iter": int(max_global_iter),
        "best_by_score": best_score.label,
        "best_converged_by_score": None if best_conv is None else best_conv.label,
        "rows": [asdict(r) for r in rows],
    }
    return rows, payload


def default_output_for_scenario(scenario: ScenarioKind) -> str:
    if scenario == "staged":
        return "phase2_vorab_march_staged_compare_baseline_kinetics.json"
    if scenario == "weights":
        return "phase2_staged_march_weights_compare.json"
    return "phase2_vorab_march_compare.json"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Compare Vorab march / staged / layer weights on Phase2 LU")
    parser.add_argument(
        "--scenario",
        choices=("staged", "weights", "both"),
        default="both",
        help="staged=march vs clip; weights=decay vs 50/30/20 × march; both=run both blocks",
    )
    add_max_global_iter_arg(parser, default=25)
    add_json_output_arg(parser, default_output_for_scenario("both"))
    args = parser.parse_args(argv)

    if args.output == default_output_for_scenario("both") and args.scenario != "both":
        args.output = default_output_for_scenario(args.scenario)

    payload: dict = {
        "scenario": args.scenario,
        "max_global_iter": int(args.max_global_iter),
        "ref_co": REF_CO,
        "ref_co2": REF_CO2,
    }

    if args.scenario in {"staged", "both"}:
        _, staged_payload = _run_staged_block(max_global_iter=int(args.max_global_iter))
        payload["staged"] = staged_payload

    if args.scenario in {"weights", "both"}:
        if args.scenario == "both":
            print("", flush=True)
        _, weights_payload = _run_weights_block(max_global_iter=int(args.max_global_iter))
        payload["weights"] = weights_payload

    if args.scenario == "staged":
        payload = payload["staged"]
    elif args.scenario == "weights":
        payload = payload["weights"]

    write_json(args.output, payload, ensure_ascii=False)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
