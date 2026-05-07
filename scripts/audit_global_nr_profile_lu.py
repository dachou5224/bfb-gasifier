#!/usr/bin/env python3
"""LU 工况 global NR profile 审计（freeboard-complete Jacobian 口径）。

固定 shared ``global_nr`` 开发口径（Vorabrechnung + block-tridiag Jacobian），
对照 ``validation_cases.json`` 中 ``CASE_HTW_WESSELING_1.outputs.axial_profiles`` 的：

- 轴向温度 `T_K`
- 湿基主气相 `CO/CO2/H2/CH4/O2/H2O`
- 出口干基主气相 `CO/CO2/H2/CH4`

目标不是继承 GS 时代的 branch judgement，而是把当前 NR 路径的 profile / exit 偏差
集中量化，便于后续针对 syngas 组成和出口结果做修正。
"""

from __future__ import annotations

import json
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
    PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS,
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


def _build_model_profiles() -> tuple[dict[str, Any], dict[str, np.ndarray], dict[str, Any], str]:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    solve_kwargs = dict(PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS)
    result, monitor = solve_with_nr_monitor(
        reactor,
        solve_kwargs,
        check_x0=True,
    )

    xi_model = np.array(result.get("axial_xi_reactor", []), dtype=float)
    T_model = np.array(result.get("T_profile", []), dtype=float)
    if xi_model.size != T_model.size or xi_model.size == 0:
        xi_model = np.array([float(c.geo.h_center / cfg.H_bed) for c in reactor.cells], dtype=float)
        T_model = np.array([float(c.T) for c in reactor.cells], dtype=float)

    bed_wet: dict[str, list[float]] = {sp: [] for sp in GAS_SPECIES}
    for cell in reactor.cells:
        y = cell._mole_fractions("combined")
        for j, sp in enumerate(GAS_SPECIES):
            bed_wet[sp].append(float(y[j]))

    freeboard_wet = result.get("freeboard_gas_profiles_wet", {})
    model_profiles: dict[str, np.ndarray] = {"xi": xi_model, "T": T_model}
    for sp in GAS_SPECIES:
        fb_vals = freeboard_wet.get(sp, [])
        vals = list(bed_wet[sp]) + [float(v) for v in fb_vals]
        if len(vals) != len(T_model):
            # Fallback for partial profile payloads: use available axial cells + reactor exit point.
            vals = list(bed_wet[sp]) + [float(result["exit_gas"].get(sp, 0.0))]
        model_profiles[sp] = np.array(vals, dtype=float)

    return result, model_profiles, monitor, str(reactor.config.major_gibbs_solver_mode)


def _available_reference_points(axial_profiles: dict[str, Any], key: str) -> tuple[np.ndarray, np.ndarray]:
    xi_raw = axial_profiles["xi"]
    y_raw = axial_profiles[key]
    pairs = [(float(x), float(y)) for x, y in zip(xi_raw, y_raw) if y is not None]
    xi = np.array([x for x, _ in pairs], dtype=float)
    y = np.array([y for _, y in pairs], dtype=float)
    return xi, y


def _rel_err(sim: float, ref: float) -> float:
    return abs(sim - ref) / max(abs(ref), 1e-12)


def run_audit() -> dict[str, Any]:
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outputs = raw["outputs"]
    axial_profiles = outputs["axial_profiles"]
    result, model, monitor, major_gibbs_solver_mode = _build_model_profiles()

    profile_summary: dict[str, Any] = {}
    for ref_key, meta in PROFILE_KEYS.items():
        xi_ref, y_ref = _available_reference_points(axial_profiles, ref_key)
        if meta["kind"] == "temperature":
            y_model = _interp_profile(xi_ref, model["xi"], model["T"])
        else:
            y_model = _interp_profile(xi_ref, model["xi"], model[meta["species"]])

        err = y_model - y_ref
        profile_summary[ref_key] = {
            "n_points": int(len(xi_ref)),
            "mae": float(np.mean(np.abs(err))),
            "rmse": float(math.sqrt(float(np.mean(err**2)))),
            "max_abs_err": float(np.max(np.abs(err))),
            "rows": [
                {
                    "xi": float(x),
                    "ref": float(y_r),
                    "model": float(y_m),
                    "error": float(e),
                }
                for x, y_r, y_m, e in zip(xi_ref, y_ref, y_model, err)
            ],
        }

    dry_ref = outputs.get("exit_gas_dry_mol_frac", {})
    exit_species: dict[str, Any] = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        ref = json_numeric_or_none(dry_ref.get(sp))
        sim = float(result["exit_gas_dry"].get(sp, 0.0))
        if ref is None:
            continue
        exit_species[sp] = {
            "ref": ref,
            "model": sim,
            "abs_err": abs(sim - ref),
            "rel_err": _rel_err(sim, ref),
        }

    t_exit_ref = float(outputs["exit_temperature_K"])
    t_exit_model = float(result["T_profile"][-1])

    return {
        "case_key": CASE_LU_VALIDATION_KEY,
        "solver": result.get("nr_init_strategy"),
        "major_gibbs_solver_mode": major_gibbs_solver_mode,
        "solve_kwargs": dict(PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS),
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "nr_total_s": result.get("nr_total_s"),
        "nr_outer_iters": result.get("nr_outer_iters"),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "T_profile_model_K": [float(v) for v in model["T"]],
        "xi_model": [float(v) for v in model["xi"]],
        "profile_summary": profile_summary,
        "exit_temperature": {
            "ref_K": t_exit_ref,
            "model_K": t_exit_model,
            "abs_err_K": abs(t_exit_model - t_exit_ref),
            "rel_err": _rel_err(t_exit_model, t_exit_ref),
        },
        "exit_species_dry": exit_species,
        "exit_gas_dry_model": {
            sp: float(result["exit_gas_dry"].get(sp, 0.0))
            for sp in ("CO", "CO2", "H2", "CH4")
        },
        "profile_observations": axial_profiles.get("_observations", {}),
        "nr_monitor": monitor,
    }


def _strict_failures(out: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    if str(out.get("major_gibbs_solver_mode", "")).strip().lower() != "hamel_reduced":
        failures.append(
            "major_gibbs_solver_mode 非 hamel_reduced（当前 gate 要求 reduced 主线）"
        )
    if not bool(out.get("converged", False)):
        failures.append("global NR 未收敛（converged=False）")
    if not bool(out.get("nr_monitor", {}).get("x0_sanity", {}).get("ok", False)):
        x0 = out.get("nr_monitor", {}).get("x0_sanity", {})
        failures.append(
            "x0 sanity 未通过: "
            f"finite={x0.get('finite_ok')} nonneg={x0.get('nonneg_ok')} "
            f"molar_scale={x0.get('molar_scale_ok')} ub={x0.get('ub_ok')} kbd={x0.get('kbd_ok')}"
        )
    return failures


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Audit LU global NR profile with freeboard-complete Jacobian")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return non-zero when convergence or x0 sanity gate fails",
    )
    args = parser.parse_args()

    try:
        out = run_audit()
    except Exception as exc:
        print("=" * 108)
        print("LU global NR profile audit")
        print("=" * 108)
        print(f"precalc/solve failed before profile summary: {exc}")
        print("-" * 108)
        print("STRICT FAILURES:")
        print("  - global NR 初始化/求解阶段失败")
        print(f"  - {exc}")
        return 1
    print("=" * 108)
    print("LU global NR profile audit")
    print("=" * 108)
    print_nr_monitor(out["nr_monitor"], prefix="NR monitor")
    print(
        f"init={out['solve_kwargs'].get('nr_init_strategy')} "
        f"major_gibbs={out.get('major_gibbs_solver_mode')} "
        f"jacobian={out['solve_kwargs'].get('nr_jacobian_strategy')} "
        f"converged={out['converged']} outer_iters={out['nr_outer_iters']} "
        f"n_iter={out['n_iter']} rms={out['rms_scaled_final']:.3e} "
        f"nr_total_s={float(out['nr_total_s']):.2f}s"
    )
    t_exit = out["exit_temperature"]
    print(
        f"exit temperature: model={t_exit['model_K']:.1f}K ref={t_exit['ref_K']:.1f}K "
        f"abs_err={t_exit['abs_err_K']:.1f}K rel_err={t_exit['rel_err']:.2%}"
    )
    print("exit gas dry:")
    for sp, row in out["exit_species_dry"].items():
        print(
            f"  {sp:4s} model={row['model']:.4f} ref={row['ref']:.4f} "
            f"abs_err={row['abs_err']:.4f} rel_err={row['rel_err']:.2%}"
        )
    print("axial profile summary:")
    for key in ("T_K", "CO_mol_wet", "CO2_mol_wet", "H2_mol_wet", "CH4_mol_wet", "O2_mol_wet", "H2O_mol_wet"):
        row = out["profile_summary"][key]
        unit = "K" if key == "T_K" else "mol/mol"
        print(
            f"  {key:12s} n={row['n_points']:>2d} "
            f"MAE={row['mae']:.4f} {unit} "
            f"RMSE={row['rmse']:.4f} {unit} "
            f"MAX={row['max_abs_err']:.4f} {unit}"
        )
    print("-" * 108)
    print("worst profile points:")
    for key in ("T_K", "CO_mol_wet", "CO2_mol_wet", "H2_mol_wet", "CH4_mol_wet", "O2_mol_wet", "H2O_mol_wet"):
        worst = max(out["profile_summary"][key]["rows"], key=lambda row: abs(row["error"]))
        print(
            f"  {key:12s} xi={worst['xi']:.2f} model={worst['model']:.4f} "
            f"ref={worst['ref']:.4f} err={worst['error']:+.4f}"
        )
    if args.strict:
        strict_failures = _strict_failures(out)
        if strict_failures:
            print("-" * 108)
            print("STRICT FAILURES:")
            for item in strict_failures:
                print(f"  - {item}")
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
