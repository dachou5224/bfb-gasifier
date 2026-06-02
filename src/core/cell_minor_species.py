"""Minor-species / Gibbs helper logic for ``Cell``."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from src.core.elemental_ledger import vm_element_molar_rates
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, get_atom_count
from src.thermodynamics import minor_species as minor_species_mod

if TYPE_CHECKING:
    from src.core.cell import Cell


def element_moles_in_gas(cell: "Cell", element: str) -> float:
    """Total element molar flow in dense-phase tracked gas species [mol/s]."""
    total = 0.0
    for sp in GAS_SPECIES:
        ni = max(float(cell.N_d[GAS_SPECIES_INDEX[sp]]), 0.0)
        total += ni * float(get_atom_count(sp, element))
    return total


def tracked_minor_element_moles_in_gas(cell: "Cell", element: str) -> float:
    """Tracked minor-element inventory inside the 11-species state vector [mol/s]."""
    tracked_species = {
        "S": ("H2S",),
        "N": ("NH3",),
    }.get(element.upper(), ())
    if not tracked_species:
        return 0.0

    total = 0.0
    el = element.upper()
    for sp in tracked_species:
        ni = max(float(cell.N_d[GAS_SPECIES_INDEX[sp]]), 0.0)
        total += ni * float(get_atom_count(sp, el))
    return total


def get_element_release(cell: "Cell", element: str) -> float:
    """Element amount available to the Gibbs minor-species subsystem [mol/s]."""
    from src.core.cell import S_VM

    el = element.upper()
    if el in {"S", "N"}:
        n = tracked_minor_element_moles_in_gas(cell, el)
    else:
        n = element_moles_in_gas(cell, el)
    m_vm = float(np.sum(cell.m_solid_zu[:, S_VM] + cell.m_solid_in[:, S_VM]))
    vm_atoms = vm_element_molar_rates(
        m_vm=m_vm,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
    )
    if el == "S":
        n += max(vm_atoms["S"], 0.0)
    elif el == "N":
        n += max(vm_atoms["N"], 0.0)
    return max(n, 0.0)


def calc_minor_species_gibbs(
    cell: "Cell",
    *,
    use_sulfur: bool = True,
    use_nitrogen: bool = True,
    k_relax: float = 1.0,
) -> dict[str, float]:
    """Gibbs minor-species relaxation sources [mol/s], written to dense phase."""
    idx = GAS_SPECIES_INDEX
    T = cell.T
    P = cell.P

    s_at = 0.0
    n_at = 0.0
    h_for_minor = c_for_minor = o_for_minor = 0.0
    n_cur_h2s = 0.0
    n_cur_nh3 = 0.0
    if use_sulfur:
        s_at = cell.get_element_release("S")
        if s_at > 1e-20:
            h_for_minor = max(cell._element_moles_in_gas("H"), 1e-6)
            o_for_minor = max(cell._element_moles_in_gas("O"), 1e-6)
            c_for_minor = max(cell._element_moles_in_gas("C"), 1e-6)
            n_cur_h2s = max(float(cell.N_d[idx["H2S"]]), 0.0)
    if use_nitrogen:
        n_at = cell.get_element_release("N")
        if n_at > 1e-20:
            if h_for_minor <= 0.0:
                h_for_minor = max(cell._element_moles_in_gas("H"), 1e-6)
            if c_for_minor <= 0.0:
                c_for_minor = max(cell._element_moles_in_gas("C"), 1e-6)
            if o_for_minor <= 0.0:
                o_for_minor = max(cell._element_moles_in_gas("O"), 1e-6)
            n_cur_nh3 = max(float(cell.N_d[idx["NH3"]]), 0.0)

    cache_sig = (
        bool(use_sulfur),
        bool(use_nitrogen),
        float(k_relax),
        float(T),
        float(P),
        float(s_at),
        float(n_at),
        float(h_for_minor),
        float(c_for_minor),
        float(o_for_minor),
        float(n_cur_h2s),
        float(n_cur_nh3),
    )
    if getattr(cell, "_minor_gibbs_cache_valid", False) and cache_sig == getattr(cell, "_minor_gibbs_cache_signature", None):
        return dict(getattr(cell, "_minor_gibbs_cache_out", {}))

    out: dict[str, float] = {}

    if use_sulfur:
        if s_at > 1e-20:
            elems_s = {
                "S": s_at,
                "H": h_for_minor,
                "O": o_for_minor,
                "C": c_for_minor,
            }
            ws_s = cell._minor_gibbs_warm_start.get("S", {})
            try:
                dist_s, diag_s = minor_species_mod.solve_sulfur_distribution(
                    T,
                    P,
                    elems_s,
                    lambda0=np.asarray(ws_s.get("lambda", []), dtype=np.float64).tolist(),
                    ln_N0=float(ws_s.get("ln_N", 0.0)),
                    return_diag=True,
                )
                cell._minor_gibbs_warm_start["S"] = {
                    "lambda": np.asarray(diag_s.get("lambda", np.zeros(4, dtype=np.float64)), dtype=np.float64),
                    "ln_N": float(diag_s.get("ln_N", 0.0)),
                }
            except TypeError:
                dist_s = minor_species_mod.solve_sulfur_distribution(T, P, elems_s)
            for sp in ("H2S",):
                if sp in dist_s and sp in idx:
                    n_eq = min(
                        max(float(dist_s[sp]), 0.0),
                        s_at / max(float(get_atom_count(sp, "S")), 1.0),
                    )
                    n_cur = max(float(cell.N_d[idx[sp]]), 0.0)
                    out[sp] = out.get(sp, 0.0) + k_relax * (n_eq - n_cur)

    if use_nitrogen:
        if n_at > 1e-20:
            elems_n = {
                "N": n_at,
                "H": h_for_minor,
                "C": c_for_minor,
                "O": o_for_minor,
            }
            ws_n = cell._minor_gibbs_warm_start.get("N", {})
            try:
                dist_n, diag_n = minor_species_mod.solve_nitrogen_distribution(
                    T,
                    P,
                    elems_n,
                    lambda0=np.asarray(ws_n.get("lambda", []), dtype=np.float64).tolist(),
                    ln_N0=float(ws_n.get("ln_N", 0.0)),
                    return_diag=True,
                )
                cell._minor_gibbs_warm_start["N"] = {
                    "lambda": np.asarray(diag_n.get("lambda", np.zeros(4, dtype=np.float64)), dtype=np.float64),
                    "ln_N": float(diag_n.get("ln_N", 0.0)),
                }
            except TypeError:
                dist_n = minor_species_mod.solve_nitrogen_distribution(T, P, elems_n)
            for sp in ("NH3",):
                if sp in dist_n and sp in idx:
                    n_eq = min(
                        max(float(dist_n[sp]), 0.0),
                        n_at / max(float(get_atom_count(sp, "N")), 1.0),
                    )
                    n_cur = max(float(cell.N_d[idx[sp]]), 0.0)
                    out[sp] = out.get(sp, 0.0) + k_relax * (n_eq - n_cur)

    cell._minor_gibbs_cache_valid = True
    cell._minor_gibbs_cache_signature = cache_sig
    cell._minor_gibbs_cache_out = dict(out)
    return out
