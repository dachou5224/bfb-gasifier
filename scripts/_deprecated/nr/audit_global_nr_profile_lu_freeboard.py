#!/usr/bin/env python3
"""LU 工况 freeboard-aware global NR profile 审计。"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    json_numeric_or_none,
    load_validation_case_node,
)


PROFILE_KEYS = {
    "T_K": {"kind": "temperature", "model_key": "T"},
    "CO_mol_wet": {"kind": "species", "species": "CO"},
    "CO2_mol_wet": {"kind": "species", "species": "CO2"},
    "H2_mol_wet": {"kind": "species", "species": "H2"},
    "CH4_mol_wet": {"kind": "species", "species": "CH4"},
    "O2_mol_wet": {"kind": "species", "species": "O2"},
    "H2O_mol_wet": {"kind": "species", "species": "H2O"},
}


def _interp_profile(x_query: np.ndarray, x_data: np.ndarray, y_data: np.ndarray) -> np.ndarray:
    return np.interp(x_query, x_data, y_data, left=y_data[0], right=y_data[-1])


def _available_reference_points(axial_profiles: dict[str, Any], key: str) -> tuple[np.ndarray, np.ndarray]:
    pairs = [
        (float(x), float(y))
        for x, y in zip(axial_profiles["xi"], axial_profiles[key])
        if y is not None
    ]
    return np.array([x for x, _ in pairs], dtype=float), np.array([y for _, y in pairs], dtype=float)


def _rel_err(sim: float, ref: float) -> float:
    return abs(sim - ref) / max(abs(ref), 1e-12)


def _build_model_profiles() -> tuple[dict[str, Any], dict[str, np.ndarray], dict[str, Any]]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    xi_bed = np.array(result["bed_axial_xi"], dtype=float)
    xi_fb = np.array(result["freeboard_axial_xi"], dtype=float)
    xi_model = np.concatenate([xi_bed, xi_fb, np.array([1.0], dtype=float)])

    wet_profiles: dict[str, list[float]] = {sp: [] for sp in GAS_SPECIES}
    for cell in reactor.cells:
        y = cell._mole_fractions("combined")
        for j, sp in enumerate(GAS_SPECIES):
            wet_profiles[sp].append(float(y[j]))
    for sp in GAS_SPECIES:
        wet_profiles[sp].extend([float(v) for v in result["freeboard_gas_profiles_wet"].get(sp, [])])
        wet_profiles[sp].append(float(result["exit_gas"].get(sp, 0.0)))

    model_profiles: dict[str, np.ndarray] = {
        "xi": xi_model,
        "T": np.concatenate(
            [
                np.array(result["bed_T_profile"], dtype=float),
                np.array(result["freeboard_T_profile"], dtype=float),
                np.array([float(result["reactor_exit_T"])]),
            ]
        ),
    }
    for sp in GAS_SPECIES:
        model_profiles[sp] = np.array(wet_profiles[sp], dtype=float)

    return result, model_profiles, monitor


def run_audit() -> dict[str, Any]:
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outputs = raw["outputs"]
    axial_profiles = outputs["axial_profiles"]
    result, model, monitor = _build_model_profiles()

    profile_summary: dict[str, Any] = {}
    for ref_key, meta in PROFILE_KEYS.items():
        xi_ref, y_ref = _available_reference_points(axial_profiles, ref_key)
        y_model = _interp_profile(
            xi_ref,
            model["xi"],
            model["T"] if meta["kind"] == "temperature" else model[meta["species"]],
        )
        err = y_model - y_ref
        profile_summary[ref_key] = {
            "n_points": int(len(xi_ref)),
            "mae": float(np.mean(np.abs(err))),
            "rmse": float(math.sqrt(float(np.mean(err**2)))),
            "max_abs_err": float(np.max(np.abs(err))),
        }

    exit_species = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        ref = json_numeric_or_none(outputs.get("exit_gas_dry_mol_frac", {}).get(sp))
        sim = float(result["exit_gas_dry"].get(sp, 0.0))
        if ref is None:
            continue
        exit_species[sp] = {
            "ref": ref,
            "model": sim,
            "rel_err": _rel_err(sim, ref),
        }

    t_exit_ref = float(outputs["exit_temperature_K"])
    t_exit_model = float(result["reactor_exit_T"])
    return {
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
        "nr_total_s": float(result.get("nr_total_s", float("nan"))),
        "profile_summary": profile_summary,
        "exit_temperature": {
            "ref_K": t_exit_ref,
            "model_K": t_exit_model,
            "rel_err": _rel_err(t_exit_model, t_exit_ref),
        },
        "exit_species_dry": exit_species,
        "bed_exit_T": float(result["bed_T_profile"][-1]),
        "reactor_exit_T": float(result["reactor_exit_T"]),
        "freeboard_active": bool(result["freeboard_active"]),
        "n_bed": len(result["bed_T_profile"]),
        "n_freeboard": len(result["freeboard_T_profile"]),
        "nr_monitor": monitor,
    }


def main() -> int:
    out = run_audit()
    print("=" * 120)
    print("LU freeboard-aware global NR profile audit")
    print("=" * 120)
    print_nr_monitor(out["nr_monitor"], prefix="NR monitor")
    print(
        f"freeboard_active={out['freeboard_active']} n_bed={out['n_bed']} n_freeboard={out['n_freeboard']} "
        f"converged={out['converged']} n_iter={out['n_iter']} rms={out['rms_scaled_final']:.3e} nr_total_s={out['nr_total_s']:.2f}s"
    )
    print(
        f"bed_exit_T={out['bed_exit_T']:.1f}K reactor_exit_T={out['reactor_exit_T']:.1f}K "
        f"ref_exit_T={out['exit_temperature']['ref_K']:.1f}K rel_err={out['exit_temperature']['rel_err']:.2%}"
    )
    print("exit gas dry:")
    for sp, row in out["exit_species_dry"].items():
        print(f"  {sp:4s} model={row['model']:.4f} ref={row['ref']:.4f} rel_err={row['rel_err']:.2%}")
    print("axial profile summary:")
    for key in ("T_K", "CO_mol_wet", "CO2_mol_wet", "H2_mol_wet", "CH4_mol_wet", "O2_mol_wet", "H2O_mol_wet"):
        row = out["profile_summary"][key]
        unit = "K" if key == "T_K" else "mol/mol"
        print(
            f"  {key:12s} n={row['n_points']:>2d} "
            f"MAE={row['mae']:.4f} {unit} RMSE={row['rmse']:.4f} {unit} MAX={row['max_abs_err']:.4f} {unit}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
