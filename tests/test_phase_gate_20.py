from __future__ import annotations

from src.core.elemental_ledger import TRACKED_ELEMENTS, build_gas_species_registry
from src.core.phase_gate_checks import gate_20_elemental_closure


def test_gas_species_registry_covers_tar_atoms_for_both_fuels():
    for fuel_type in ("coal", "biomass"):
        registry = build_gas_species_registry(fuel_type)
        for tar_label in ("TAR1", "TAR2"):
            atoms = registry[tar_label]["atoms"]
            assert atoms["C"] > 0
            assert atoms["H"] > 0
            assert set(atoms) == set(TRACKED_ELEMENTS)


def test_gate_20_elemental_closure_smoke():
    result = gate_20_elemental_closure()
    assert result.passed, result.hard_failures
    assert result.metrics["pyrolysis"]["max_abs_closure"] <= 1e-8
    assert result.metrics["reaction"]["max_abs_closure"] <= 1e-7
    assert result.metrics["reaction"]["reaction_numbering"]["authority"] == "hamel_1999_thesis"
    assert "R11" in result.metrics["reaction"]["net_molar_gas_source_by_thesis"]
