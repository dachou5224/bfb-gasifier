from __future__ import annotations

import pytest

from src.core.phase_gate_checks import (
    build_phase1_structural_contract,
    gate_10_structural_parity,
    validate_phase1_structural_contract,
)


def test_phase1_structural_contract_has_unique_expected_final_owners():
    contract = build_phase1_structural_contract()
    failures = validate_phase1_structural_contract(contract)
    assert not failures, failures


@pytest.mark.slow
def test_gate_10_structural_parity_smoke():
    result = gate_10_structural_parity()
    assert result.passed, result.hard_failures
    assert "structural_contract" in result.metrics
    assert result.metrics["sanity"]["passed"]
