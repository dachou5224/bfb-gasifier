"""LU stable-window local calibration with convergence-aware scoring.

Purpose
-------
Search a small neighborhood around the current LU baseline while penalizing
high GS residual/rms.  This avoids selecting parameter sets that look closer
to validation targets but are materially less self-consistent.

Usage
-----
    python3 scripts/calibrate_lu_local_window.py
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config


REF_DRY = {"CO": 0.157, "CO2": 0.133, "H2": 0.145, "CH4": 0.034}
REF_T_EXIT = 1120.0
REF_XC = 0.95


@dataclass
class Row:
    dense: float
    r4: float
    r5: float
    r7: float
    n_iter: int
    t_exit: float
    xc: float
    rms: float
    co: float
    co2: float
    h2: float
    ch4: float
    score_raw: float
    score_conv: float


def _score_raw(out: dict) -> float:
    y = out["exit_gas_dry"]
    score = sum(abs(float(y.get(sp, 0.0)) - tgt) / tgt for sp, tgt in REF_DRY.items())
    score += abs(float(out["T_profile"][-1]) - REF_T_EXIT) / REF_T_EXIT
    score += abs(float(out["carbon_conv"]) - REF_XC) / REF_XC
    return float(score)


def _score_conv(out: dict) -> float:
    """Validation score + explicit GS convergence penalty."""
    raw = _score_raw(out)
    rms = float(out.get("rms_scaled_gs", 0.0))
    # Below ~0.15 treat as acceptable stable branch; above that, penalize quickly.
    penalty = 2.0 * max(rms - 0.15, 0.0)
    return float(raw + penalty)


def run_one(dense: float, r4: float, r5: float, r7: float) -> Row:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.gas_inlet_dense_frac = dense
    cfg.r4_scale = r4
    cfg.r5_scale = r5
    cfg.r7_scale = r7

    out = Reactor(cfg).solve(max_global_iter=3, tol_global=1.0, solver="gauss_seidel")
    y = out["exit_gas_dry"]
    return Row(
        dense=dense,
        r4=r4,
        r5=r5,
        r7=r7,
        n_iter=int(out["n_iter"]),
        t_exit=float(out["T_profile"][-1]),
        xc=float(out["carbon_conv"]),
        rms=float(out["rms_scaled_gs"]),
        co=float(y.get("CO", 0.0)),
        co2=float(y.get("CO2", 0.0)),
        h2=float(y.get("H2", 0.0)),
        ch4=float(y.get("CH4", 0.0)),
        score_raw=_score_raw(out),
        score_conv=_score_conv(out),
    )


def main() -> int:
    dense_grid = [0.30, 0.35]
    r4_grid = [0.50]
    r5_grid = [0.75]
    r7_grid = [2.00, 2.50, 3.00]

    rows: list[Row] = []
    print("=" * 132)
    print("LU local calibration window (convergence-aware)")
    print("=" * 132)

    n_total = len(dense_grid) * len(r4_grid) * len(r5_grid) * len(r7_grid)
    n = 0
    for dense in dense_grid:
        for r4 in r4_grid:
            for r5 in r5_grid:
                for r7 in r7_grid:
                    n += 1
                    row = run_one(dense, r4, r5, r7)
                    rows.append(row)
                    print(
                        f"[{n:02d}/{n_total}] dense={dense:.2f} r4={r4:.2f} r5={r5:.2f} r7={r7:.2f} | "
                        f"T={row.t_exit:.1f} Xc={row.xc:.3f} rms={row.rms:.3f} "
                        f"CO={row.co:.3f} CO2={row.co2:.3f} H2={row.h2:.3f} CH4={row.ch4:.3f} | "
                        f"raw={row.score_raw:.3f} conv={row.score_conv:.3f}"
                    )

    rows_raw = sorted(rows, key=lambda r: r.score_raw)
    rows_conv = sorted(rows, key=lambda r: r.score_conv)

    print("\n" + "=" * 132)
    print("Top By Raw Validation Score")
    print("=" * 132)
    for i, r in enumerate(rows_raw[:5], start=1):
        print(
            f"#{i} dense={r.dense:.2f} r4={r.r4:.2f} r5={r.r5:.2f} r7={r.r7:.2f} | "
            f"T={r.t_exit:.1f} Xc={r.xc:.3f} rms={r.rms:.3f} "
            f"CO={r.co:.3f} CO2={r.co2:.3f} H2={r.h2:.3f} CH4={r.ch4:.3f} | "
            f"raw={r.score_raw:.3f} conv={r.score_conv:.3f}"
        )

    print("\n" + "=" * 132)
    print("Top By Convergence-Aware Score")
    print("=" * 132)
    for i, r in enumerate(rows_conv[:5], start=1):
        print(
            f"#{i} dense={r.dense:.2f} r4={r.r4:.2f} r5={r.r5:.2f} r7={r.r7:.2f} | "
            f"T={r.t_exit:.1f} Xc={r.xc:.3f} rms={r.rms:.3f} "
            f"CO={r.co:.3f} CO2={r.co2:.3f} H2={r.h2:.3f} CH4={r.ch4:.3f} | "
            f"raw={r.score_raw:.3f} conv={r.score_conv:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
