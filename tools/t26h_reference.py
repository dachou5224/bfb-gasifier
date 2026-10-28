"""Dump the round-1 reference coefficients of case 1 for the T26h comparison report.

Only reads and prints; does not change any model quantity.
"""
from __future__ import annotations

import json
import time

import numpy as np

from hamel.cases.config import SolverConfig
from hamel.solver.precalc import precalculate


def main(out: str = "t26h-baseline.json") -> None:
    started = time.perf_counter()
    prepared = precalculate("data/cases/vtt_pressurized.json", "1", None, adiabatic=False, iteration=1, solver=SolverConfig())
    elapsed = time.perf_counter() - started
    frozen = prepared.frozen
    payload = {
        "elapsed_s": elapsed,
        "t_ad_k": float(prepared.t_ad_k),
        "t_init_k": float(prepared.t_init_k),
        "d_bar_m": float(prepared.d_bar_m),
        "d_sv_m": float(prepared.d_sv_m),
        "rho_p_kg_m3": float(frozen.rho_p_kg_m3),
        "class_diameter_m": [float(value) for value in frozen.diameter_m],
        "edges_m": [float(value) for value in prepared.reported["trajectory"]["edges_m"]],
        "a_free": np.asarray(frozen.a_free_per_s, dtype=float).tolist(),
        "a_hover": np.asarray(frozen.a_hover_per_s, dtype=float).tolist(),
        "b_free": np.asarray(frozen.b_free_per_s, dtype=float).tolist(),
        "b_hover": np.asarray(frozen.b_hover_per_s, dtype=float).tolist(),
        "hover": np.asarray(frozen.hover, dtype=bool).tolist(),
        "theta": np.asarray(frozen.theta, dtype=float).tolist(),
        "ejection_to_free": np.asarray(frozen.ejection_to_free, dtype=float).tolist(),
        "ejection_to_hover": np.asarray(frozen.ejection_to_hover, dtype=float).tolist(),
        "reference_return_share": np.asarray(prepared.reported["reference_return_share"], dtype=float).tolist(),
        "reference_elutriation_share": np.asarray(prepared.reported["reference_elutriation_share"], dtype=float).tolist(),
        "trajectory": {
            key: (value.tolist() if hasattr(value, "tolist") else value)
            for key, value in prepared.reported["trajectory"].items()
        },
    }
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=1)
    print(f"elapsed {elapsed:.1f} s -> {out}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else "t26h-baseline.json")
