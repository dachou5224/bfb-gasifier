#!/usr/bin/env python3
"""LU 工况下 VM 气源分相灵敏度审计。

目的：
- 固定 shared ``global_nr`` 开发口径（Vorabrechnung + block-tridiag Jacobian）
- 只改变热解气源 ``gas_src_vm`` 在 dense / bubble 两相的分配
- 观察对出口 syngas、轴向温度/组分 profile，以及 NR 迭代步数的影响

说明：
- 当前主模型实现把热解气完全写入 ``R_gas_d``。
- 这里通过临时 monkeypatch ``src.core.cell.build_reaction_sources`` 做
  dev-only sensitivity，不改默认实现。
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import src.core.cell as cell_mod
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
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
    xi = np.array([x for x, _ in pairs], dtype=float)
    y = np.array([y for _, y in pairs], dtype=float)
    return xi, y


def _rel_err(sim: float, ref: float) -> float:
    return abs(sim - ref) / max(abs(ref), 1e-12)


def _make_phase_split_wrapper(vm_to_bubble_frac: float):
    original = cell_mod.build_reaction_sources

    def wrapped(**kwargs):
        gas = np.asarray(kwargs["gas_src_vm"], dtype=np.float64)
        gas_to_bubble = gas * float(vm_to_bubble_frac)
        gas_to_dense = gas - gas_to_bubble
        kwargs = dict(kwargs)
        kwargs["gas_src_vm"] = gas_to_dense
        out = original(**kwargs)
        out.R_gas_b += gas_to_bubble
        return out

    return original, wrapped


def _solve_fraction(vm_to_bubble_frac: float) -> dict[str, Any]:
    original, wrapped = _make_phase_split_wrapper(vm_to_bubble_frac)
    cell_mod.build_reaction_sources = wrapped
    try:
        cfg = build_phase1_htw_lu_global_nr_reactor_config()
        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
            check_x0=True,
        )
    finally:
        cell_mod.build_reaction_sources = original

    xi_cells = np.array([float(c.geo.h_center / cfg.H_bed) for c in reactor.cells], dtype=float)
    xi_model = np.concatenate([xi_cells, np.array([1.0], dtype=float)])
    wet_profiles: dict[str, list[float]] = {sp: [] for sp in GAS_SPECIES}
    for cell in reactor.cells:
        y = cell._mole_fractions("combined")
        for j, sp in enumerate(GAS_SPECIES):
            wet_profiles[sp].append(float(y[j]))
    for sp in GAS_SPECIES:
        wet_profiles[sp].append(float(result["exit_gas"].get(sp, 0.0)))

    model_profiles: dict[str, np.ndarray] = {
        "xi": xi_model,
        "T": np.concatenate(
            [np.array(result["T_profile"], dtype=float), np.array([float(result["T_profile"][-1])])]
        ),
    }
    for sp in GAS_SPECIES:
        model_profiles[sp] = np.array(wet_profiles[sp], dtype=float)

    cell0 = reactor.cells[0]
    cell0_y = cell0._mole_fractions("combined")
    out = {
        "vm_to_bubble_frac": float(vm_to_bubble_frac),
        "Texit": float(result["T_profile"][-1]),
        "Tpeak": float(max(result["T_profile"])),
        "Tpeak_cell": int(np.argmax(np.array(result["T_profile"], dtype=float))),
        "carbon_conv": float(result["carbon_conv"]),
        "rms_scaled_final": float(result.get("rms_scaled_final", float("nan"))),
        "n_iter": int(result.get("n_iter", 0)),
        "nr_outer_iters": int(result.get("nr_outer_iters", 0)),
        "wall_time_s": float(monitor.get("wall_time_s", float("nan"))),
        "nr_total_s": float(result.get("nr_total_s", float("nan"))),
        "x0_ok": bool((monitor.get("x0_sanity") or {}).get("ok", False)),
        "exit_gas_dry": {
            sp: float(result["exit_gas_dry"].get(sp, 0.0))
            for sp in ("CO", "CO2", "H2", "CH4")
        },
        "cell0": {
            "y_O2": float(cell0_y[GAS_SPECIES_INDEX["O2"]]),
            "y_CH4": float(cell0_y[GAS_SPECIES_INDEX["CH4"]]),
            "Nd_CH4": float(cell0.N_d[GAS_SPECIES_INDEX["CH4"]]),
            "Nb_CH4": float(cell0.N_b[GAS_SPECIES_INDEX["CH4"]]),
            "Nd_O2": float(cell0.N_d[GAS_SPECIES_INDEX["O2"]]),
            "Nb_O2": float(cell0.N_b[GAS_SPECIES_INDEX["O2"]]),
        },
        "model_profiles": model_profiles,
    }
    return out


def run_audit() -> dict[str, Any]:
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outputs = raw["outputs"]
    axial_profiles = outputs["axial_profiles"]
    fractions = (0.0, 0.5, 1.0)
    rows = [_solve_fraction(frac) for frac in fractions]

    for row in rows:
        profile_summary: dict[str, Any] = {}
        for ref_key, meta in PROFILE_KEYS.items():
            xi_ref, y_ref = _available_reference_points(axial_profiles, ref_key)
            if meta["kind"] == "temperature":
                y_model = _interp_profile(xi_ref, row["model_profiles"]["xi"], row["model_profiles"]["T"])
            else:
                y_model = _interp_profile(
                    xi_ref,
                    row["model_profiles"]["xi"],
                    row["model_profiles"][meta["species"]],
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
            sim = float(row["exit_gas_dry"][sp])
            if ref is None:
                continue
            exit_species[sp] = {
                "ref": ref,
                "model": sim,
                "rel_err": _rel_err(sim, ref),
            }
        row["profile_summary"] = profile_summary
        row["exit_species_summary"] = exit_species
        row.pop("model_profiles", None)

    return {"fractions": rows}


def main() -> int:
    out = run_audit()
    print("=" * 150)
    print("LU VM phase-split audit (shared global NR)")
    print("=" * 150)
    print(
        f"{'vm->bubble':>10} {'Texit[K]':>10} {'Tpeak[K]':>10} {'peak@':>6} "
        f"{'COdry':>8} {'CO2dry':>8} {'H2dry':>8} {'CH4dry':>8} "
        f"{'rms':>10} {'n_iter':>8} {'outer':>7} {'wall[s]':>9} {'nr[s]':>8} {'x0_ok':>7} {'cell0 yO2':>10}"
    )
    print("-" * 150)
    for row in out["fractions"]:
        print(
            f"{row['vm_to_bubble_frac']:>10.2f} {row['Texit']:>10.1f} {row['Tpeak']:>10.1f} {row['Tpeak_cell']:>6d} "
            f"{row['exit_gas_dry']['CO']:>8.4f} {row['exit_gas_dry']['CO2']:>8.4f} {row['exit_gas_dry']['H2']:>8.4f} {row['exit_gas_dry']['CH4']:>8.4f} "
            f"{row['rms_scaled_final']:>10.3e} {row['n_iter']:>8d} {row['nr_outer_iters']:>7d} {row['wall_time_s']:>9.2f} "
            f"{row['nr_total_s']:>8.2f} {str(row['x0_ok']):>7} {row['cell0']['y_O2']:>10.4f}"
        )
    print("-" * 150)
    print("profile MAE:")
    for row in out["fractions"]:
        s = row["profile_summary"]
        print(
            f"  vm->bubble={row['vm_to_bubble_frac']:.2f} "
            f"T={s['T_K']['mae']:.1f}K CO={s['CO_mol_wet']['mae']:.4f} CO2={s['CO2_mol_wet']['mae']:.4f} "
            f"H2={s['H2_mol_wet']['mae']:.4f} CH4={s['CH4_mol_wet']['mae']:.4f} O2={s['O2_mol_wet']['mae']:.4f} H2O={s['H2O_mol_wet']['mae']:.4f}"
        )
    print("-" * 150)
    print("exit rel err vs LU ref:")
    for row in out["fractions"]:
        s = row["exit_species_summary"]
        print(
            f"  vm->bubble={row['vm_to_bubble_frac']:.2f} "
            f"CO={s['CO']['rel_err']:.2%} CO2={s['CO2']['rel_err']:.2%} "
            f"H2={s['H2']['rel_err']:.2%} CH4={s['CH4']['rel_err']:.2%}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
