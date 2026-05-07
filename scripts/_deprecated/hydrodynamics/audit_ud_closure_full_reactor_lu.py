#!/usr/bin/env python3
"""LU 工况 `u_d` closure 对整炉结果的影响审计。"""

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
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
    json_numeric_or_none,
    load_validation_case_node,
)


_CLOSURES = ("current", "backsolve_visible_epsb")
_PROFILE_KEYS = {
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
    pairs = [(float(x), float(y)) for x, y in zip(axial_profiles["xi"], axial_profiles[key]) if y is not None]
    return np.array([x for x, _ in pairs], dtype=float), np.array([y for _, y in pairs], dtype=float)


def _rel_err(sim: float, ref: float) -> float:
    return abs(sim - ref) / max(abs(ref), 1e-12)


def _build_model_profiles(closure: str) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.hydrodynamics_u_d_closure = closure
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
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

    merged = dict(result)
    merged["_nr_monitor"] = monitor
    return merged, model_profiles


def _summarize_closure(closure: str, outputs: dict[str, Any]) -> dict[str, Any]:
    axial_profiles = outputs["axial_profiles"]
    result, model = _build_model_profiles(closure)

    profile_summary: dict[str, Any] = {}
    for ref_key, meta in _PROFILE_KEYS.items():
        xi_ref, y_ref = _available_reference_points(axial_profiles, ref_key)
        y_model = _interp_profile(
            xi_ref,
            model["xi"],
            model["T"] if meta["kind"] == "temperature" else model[meta["species"]],
        )
        err = y_model - y_ref
        profile_summary[ref_key] = {
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
        exit_species[sp] = {"ref": ref, "model": sim, "rel_err": _rel_err(sim, ref)}

    return {
        "closure": closure,
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
        "wall_time_s": float((result.get("_nr_monitor") or {}).get("wall_time_s", float("nan"))),
        "nr_total_s": float(result.get("nr_total_s", float("nan"))),
        "x0_ok": bool(((result.get("_nr_monitor") or {}).get("x0_sanity") or {}).get("ok", False)),
        "bed_exit_T": float(result["bed_T_profile"][-1]),
        "reactor_exit_T": float(result["reactor_exit_T"]),
        "exit_species_dry": exit_species,
        "profile_summary": profile_summary,
        "freeboard_active": bool(result["freeboard_active"]),
    }


def main() -> int:
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outputs = raw["outputs"]
    summaries = [_summarize_closure(closure, outputs) for closure in _CLOSURES]

    print("=" * 172)
    print("LU full-reactor u_d-closure audit")
    print("=" * 172)
    print(
        f"{'closure':<24} {'Texit[K]':>9} {'Tbed[K]':>9} {'n_iter':>7} {'rms':>10} {'wall[s]':>8} {'x0_ok':>7} "
        f"{'CO':>8} {'CO2':>8} {'H2':>8} {'CH4':>8} {'MAE_T':>10} {'MAE_CO':>10} {'MAE_H2O':>10}"
    )
    print("-" * 172)
    for row in summaries:
        print(
            f"{row['closure']:<24} {row['reactor_exit_T']:>9.1f} {row['bed_exit_T']:>9.1f} "
            f"{row['n_iter']:>7d} {row['rms_scaled_final']:>10.3e} {row['wall_time_s']:>8.2f} {str(row['x0_ok']):>7} "
            f"{row['exit_species_dry']['CO']['model']:>8.4f} {row['exit_species_dry']['CO2']['model']:>8.4f} "
            f"{row['exit_species_dry']['H2']['model']:>8.4f} {row['exit_species_dry']['CH4']['model']:>8.4f} "
            f"{row['profile_summary']['T_K']['mae']:>10.2f} {row['profile_summary']['CO_mol_wet']['mae']:>10.4f} "
            f"{row['profile_summary']['H2O_mol_wet']['mae']:>10.4f}"
        )

    print("-" * 172)
    print("exit dry-gas relative errors:")
    for row in summaries:
        print(f"  {row['closure']}:")
        for sp in ("CO", "CO2", "H2", "CH4"):
            ent = row["exit_species_dry"][sp]
            print(f"    {sp:4s} model={ent['model']:.4f} ref={ent['ref']:.4f} rel_err={ent['rel_err']:.2%}")
    print("-" * 172)
    print("readout:")
    print("  - 若 `backsolve_visible_epsb` 让 `MAE_T` / `CO/CO2/H2/H2O` 同时改善，则可以把它继续当成 bed-hydrodynamics audit branch。")
    print("  - 若它只改善 hydrodynamics 但把整炉 profile 明显推坏，则说明 `u_d` 之外还有 `d_b / lambda_b / K_bd` 的协同偏差。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
