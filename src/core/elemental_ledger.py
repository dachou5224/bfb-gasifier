"""Central elemental ledger helpers for gas/solid source audits."""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import (
    GAS_SPECIES,
    MOLECULAR_WEIGHT,
    TAR_SURROGATE_FORMULA,
    configure_tar_components_by_fuel,
    get_atom_count,
    get_tar_component_mapping,
)

TRACKED_ELEMENTS: tuple[str, ...] = ("C", "H", "O", "N", "S")
ATOMIC_MASS_KG_PER_MOL: dict[str, float] = {
    "C": 12.011e-3,
    "H": 1.00794e-3,
    "O": 15.999e-3,
    "N": 14.0067e-3,
    "S": 32.065e-3,
}
H2O_MOLAR_MASS_KG_PER_MOL: float = MOLECULAR_WEIGHT["H2O"] / 1000.0


def build_gas_species_registry(fuel_type: str) -> dict[str, dict[str, Any]]:
    """Return element counts and molar masses for all gas species."""

    configure_tar_components_by_fuel(fuel_type)
    tar_mapping = get_tar_component_mapping(fuel_type)
    registry: dict[str, dict[str, Any]] = {}

    for species in GAS_SPECIES:
        if species in ("TAR1", "TAR2"):
            surrogate = tar_mapping[species]
            c_atoms, h_atoms = TAR_SURROGATE_FORMULA[surrogate]
            atoms = {"C": c_atoms, "H": h_atoms, "O": 0, "N": 0, "S": 0}
        else:
            atoms = {element: int(get_atom_count(species, element)) for element in TRACKED_ELEMENTS}
        registry[species] = {
            "molar_mass_kg_per_mol": float(MOLECULAR_WEIGHT[species] / 1000.0),
            "atoms": atoms,
        }

    return registry


def gas_element_molar_rates(
    gas_rates: npt.ArrayLike,
    *,
    fuel_type: str,
) -> dict[str, float]:
    """Project gas species rates [mol/s] into elemental rates [mol-atom/s]."""

    registry = build_gas_species_registry(fuel_type)
    gas_arr = np.asarray(gas_rates, dtype=np.float64)
    assert gas_arr.shape == (len(GAS_SPECIES),), f"gas_rates shape mismatch: {gas_arr.shape}"

    out = {element: 0.0 for element in TRACKED_ELEMENTS}
    for i, species in enumerate(GAS_SPECIES):
        atoms = registry[species]["atoms"]
        for element in TRACKED_ELEMENTS:
            out[element] += float(gas_arr[i]) * float(atoms[element])
    return out


def vm_element_molar_rates(
    *,
    m_vm: float,
    ash_dry_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    nitrogen_fraction: float,
    sulfur_fraction: float,
    sulfur_volatile_frac: float,
) -> dict[str, float]:
    """Convert VM mass flow [kg/s] into elemental molar rates [mol-atom/s]."""

    to_daf = 1.0 / max(1.0 - float(ash_dry_wt) / 100.0, 1e-9)
    m_vm_val = float(m_vm)
    return {
        "C": m_vm_val * (float(C_dry) / 100.0) * to_daf / ATOMIC_MASS_KG_PER_MOL["C"],
        "H": m_vm_val * (float(H_dry) / 100.0) * to_daf / ATOMIC_MASS_KG_PER_MOL["H"],
        "O": m_vm_val * (float(O_dry) / 100.0) * to_daf / ATOMIC_MASS_KG_PER_MOL["O"],
        "N": m_vm_val * (float(nitrogen_fraction) / 100.0) * to_daf / ATOMIC_MASS_KG_PER_MOL["N"],
        "S": m_vm_val * float(sulfur_fraction) * float(sulfur_volatile_frac) * to_daf / ATOMIC_MASS_KG_PER_MOL["S"],
    }


def solid_element_molar_rates(
    solid_rates: npt.ArrayLike,
    *,
    ash_dry_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    nitrogen_fraction: float,
    sulfur_fraction: float,
    sulfur_volatile_frac: float,
    char_index: int,
    vm_index: int,
    moisture_index: int,
) -> dict[str, float]:
    """Project solid component rates [kg/s] into elemental rates [mol-atom/s]."""

    solid_arr = np.asarray(solid_rates, dtype=np.float64)
    assert solid_arr.ndim == 2, f"solid_rates must be 2D, got {solid_arr.shape}"

    out = {element: 0.0 for element in TRACKED_ELEMENTS}
    char_mass = float(np.sum(solid_arr[:, char_index]))
    out["C"] += char_mass / ATOMIC_MASS_KG_PER_MOL["C"]

    vm_mass = float(np.sum(solid_arr[:, vm_index]))
    vm_rates = vm_element_molar_rates(
        m_vm=vm_mass,
        ash_dry_wt=ash_dry_wt,
        C_dry=C_dry,
        H_dry=H_dry,
        O_dry=O_dry,
        nitrogen_fraction=nitrogen_fraction,
        sulfur_fraction=sulfur_fraction,
        sulfur_volatile_frac=sulfur_volatile_frac,
    )
    for element, value in vm_rates.items():
        out[element] += float(value)

    moisture_mass = float(np.sum(solid_arr[:, moisture_index]))
    n_h2o = moisture_mass / H2O_MOLAR_MASS_KG_PER_MOL
    out["H"] += 2.0 * n_h2o
    out["O"] += 1.0 * n_h2o
    return out


def combined_element_molar_rates(
    *,
    gas_rates: npt.ArrayLike,
    solid_rates: npt.ArrayLike,
    fuel_type: str,
    ash_dry_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    nitrogen_fraction: float,
    sulfur_fraction: float,
    sulfur_volatile_frac: float,
    char_index: int,
    vm_index: int,
    moisture_index: int,
) -> dict[str, float]:
    """Combine gas and solid source terms into net elemental rates."""

    gas_part = gas_element_molar_rates(gas_rates, fuel_type=fuel_type)
    solid_part = solid_element_molar_rates(
        solid_rates,
        ash_dry_wt=ash_dry_wt,
        C_dry=C_dry,
        H_dry=H_dry,
        O_dry=O_dry,
        nitrogen_fraction=nitrogen_fraction,
        sulfur_fraction=sulfur_fraction,
        sulfur_volatile_frac=sulfur_volatile_frac,
        char_index=char_index,
        vm_index=vm_index,
        moisture_index=moisture_index,
    )
    return {
        element: float(gas_part[element] + solid_part[element])
        for element in TRACKED_ELEMENTS
    }


def max_absolute_element_rate(element_rates: dict[str, float]) -> float:
    """Largest absolute elemental imbalance [mol-atom/s]."""

    return max(abs(float(v)) for v in element_rates.values())
