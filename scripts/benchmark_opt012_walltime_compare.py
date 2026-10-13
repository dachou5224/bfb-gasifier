#!/usr/bin/env python3
"""OPT-012 walltime A/B：partial residual + local BC vs legacy full residual + full BC。

同机、同配置，仅 shim 预投影 helper，不改物理接受准则。
"""

from __future__ import annotations

import argparse
import json
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, Iterator

import numpy as np

import _lu_harness  # noqa: F401

from _lu_harness import (
    LADDER_P2,
    apply_config_overrides,
    build_phase2_config,
    phase2_solve_kwargs,
    repo_root,
    run_phase2_init,
    write_json,
)
from src.core.connectivity import (
    run_bottom_zone_solid_energy_preprojection,
    run_outer_preinner_start_value_bridge,
)
from src.core.reactor import Reactor


@dataclass(frozen=True)
class TimingRow:
    label: str
    init_s: float
    bottom_zone_s: float
    preinner_bridge_s: float
    init_plus_preinner_s: float


def _build_reactor() -> Reactor:
    cfg = build_phase2_config()
    apply_config_overrides(cfg, LADDER_P2)
    return Reactor(cfg)


@contextmanager
def _legacy_preproj_shim() -> Iterator[None]:
    import src.core.connectivity as conn

    saved_vector = conn._preproj_cell_residual_vector
    saved_energy = conn._preproj_cell_energy_residual_abs
    saved_trial_bc = conn._preproj_refresh_trial_bc

    def _legacy_vector(cell: Any, *groups: str) -> np.ndarray:
        del groups
        return np.asarray(cell.residuals(), dtype=np.float64)

    def _legacy_energy_abs(cell: Any) -> float:
        return abs(float(np.asarray(cell.residuals(), dtype=np.float64)[-1]))

    def _legacy_trial_bc(
        cell_index: int,
        *,
        apply_bc_fn,
        apply_local_bc_fn=None,
        neighbor_indices: tuple[int, ...] = (),
    ) -> None:
        del cell_index, apply_local_bc_fn, neighbor_indices
        if apply_bc_fn is not None:
            apply_bc_fn()

    conn._preproj_cell_residual_vector = _legacy_vector
    conn._preproj_cell_energy_residual_abs = _legacy_energy_abs
    conn._preproj_refresh_trial_bc = _legacy_trial_bc
    try:
        yield
    finally:
        conn._preproj_cell_residual_vector = saved_vector
        conn._preproj_cell_energy_residual_abs = saved_energy
        conn._preproj_refresh_trial_bc = saved_trial_bc


def _time_paths(*, label: str, legacy: bool) -> TimingRow:
    ctx = _legacy_preproj_shim() if legacy else _nullcontext()
    with ctx:
        reactor = _build_reactor()
        sk = phase2_solve_kwargs()
        init_strategy = str(sk.get("nr_init_strategy", "vorabrechnung"))

        t0 = perf_counter()
        run_phase2_init(reactor, init_strategy=init_strategy)
        init_s = perf_counter() - t0

        cells = reactor._solver_cells_for_nr()
        apply = reactor._apply_all_bc_for_nr
        apply_local = reactor._apply_local_bc_for_nr

        t0 = perf_counter()
        run_bottom_zone_solid_energy_preprojection(
            cells,
            reactor.config,
            apply_bc_fn=apply,
            apply_local_bc_fn=apply_local,
        )
        bottom_zone_s = perf_counter() - t0

        t0 = perf_counter()
        run_outer_preinner_start_value_bridge(
            cells,
            reactor.config,
            apply_bc_fn=apply,
            apply_local_bc_fn=apply_local,
        )
        preinner_bridge_s = perf_counter() - t0

    return TimingRow(
        label=label,
        init_s=init_s,
        bottom_zone_s=bottom_zone_s,
        preinner_bridge_s=preinner_bridge_s,
        init_plus_preinner_s=init_s + bottom_zone_s + preinner_bridge_s,
    )


@contextmanager
def _nullcontext() -> Iterator[None]:
    yield


def _speedup(before: float, after: float) -> float | None:
    if before <= 0.0 or after <= 0.0:
        return None
    return before / after


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Each mode (current/legacy) is timed this many times; report median",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root() / "data" / "opt012_walltime_compare_latest.json",
    )
    args = parser.parse_args()
    repeats = max(int(args.repeats), 1)

    current_rows: list[TimingRow] = []
    legacy_rows: list[TimingRow] = []
    # 交替顺序，减轻单次冷启动偏差
    for i in range(repeats):
        if i % 2 == 0:
            print(f"[{i + 1}/{repeats}] current ...", flush=True)
            current_rows.append(_time_paths(label="opt012_current", legacy=False))
            print(f"[{i + 1}/{repeats}] legacy ...", flush=True)
            legacy_rows.append(_time_paths(label="legacy_shim", legacy=True))
        else:
            print(f"[{i + 1}/{repeats}] legacy ...", flush=True)
            legacy_rows.append(_time_paths(label="legacy_shim", legacy=True))
            print(f"[{i + 1}/{repeats}] current ...", flush=True)
            current_rows.append(_time_paths(label="opt012_current", legacy=False))

    def _median(rows: list[TimingRow], field: str) -> float:
        vals = sorted(float(getattr(r, field)) for r in rows)
        mid = len(vals) // 2
        if len(vals) % 2:
            return vals[mid]
        return 0.5 * (vals[mid - 1] + vals[mid])

    current = TimingRow(
        label="opt012_current_median",
        init_s=_median(current_rows, "init_s"),
        bottom_zone_s=_median(current_rows, "bottom_zone_s"),
        preinner_bridge_s=_median(current_rows, "preinner_bridge_s"),
        init_plus_preinner_s=_median(current_rows, "init_plus_preinner_s"),
    )
    legacy = TimingRow(
        label="legacy_shim_median",
        init_s=_median(legacy_rows, "init_s"),
        bottom_zone_s=_median(legacy_rows, "bottom_zone_s"),
        preinner_bridge_s=_median(legacy_rows, "preinner_bridge_s"),
        init_plus_preinner_s=_median(legacy_rows, "init_plus_preinner_s"),
    )

    payload = {
        "repeats": repeats,
        "note": "legacy_shim restores full cell.residuals() and apply_bc_fn() inside preprojection trial loops",
        "current_runs": [r.__dict__ for r in current_rows],
        "legacy_runs": [r.__dict__ for r in legacy_rows],
        "current_median": current.__dict__,
        "legacy_median": legacy.__dict__,
        "speedup_current_vs_legacy": {
            "init_s": _speedup(legacy.init_s, current.init_s),
            "bottom_zone_s": _speedup(legacy.bottom_zone_s, current.bottom_zone_s),
            "preinner_bridge_s": _speedup(legacy.preinner_bridge_s, current.preinner_bridge_s),
            "init_plus_preinner_s": _speedup(legacy.init_plus_preinner_s, current.init_plus_preinner_s),
        },
        "delta_s": {
            "init_s": round(legacy.init_s - current.init_s, 3),
            "bottom_zone_s": round(legacy.bottom_zone_s - current.bottom_zone_s, 3),
            "preinner_bridge_s": round(legacy.preinner_bridge_s - current.preinner_bridge_s, 3),
            "init_plus_preinner_s": round(legacy.init_plus_preinner_s - current.init_plus_preinner_s, 3),
        },
        "pct_saved": {
            "init_s": round(100.0 * (legacy.init_s - current.init_s) / legacy.init_s, 1)
            if legacy.init_s > 0
            else None,
            "bottom_zone_s": round(100.0 * (legacy.bottom_zone_s - current.bottom_zone_s) / legacy.bottom_zone_s, 1)
            if legacy.bottom_zone_s > 0
            else None,
            "preinner_bridge_s": round(
                100.0 * (legacy.preinner_bridge_s - current.preinner_bridge_s) / legacy.preinner_bridge_s,
                1,
            )
            if legacy.preinner_bridge_s > 0
            else None,
            "init_plus_preinner_s": round(
                100.0 * (legacy.init_plus_preinner_s - current.init_plus_preinner_s) / legacy.init_plus_preinner_s,
                1,
            )
            if legacy.init_plus_preinner_s > 0
            else None,
        },
    }

    write_json(args.output, payload)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
