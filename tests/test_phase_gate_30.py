from __future__ import annotations

import pytest

from src.core.phase_gate_checks import gate_30_phase_partition_closure


def test_gate_30_phase_partition_closure_smoke():
    result = gate_30_phase_partition_closure()
    assert result.passed, result.hard_failures
    ledger = result.metrics["phase_partition_ledger"]
    assert abs(ledger["N_ex_sum"]) <= 1e-12
    assert ledger["o2_slack_bubble"] >= 0.0
    assert ledger["o2_shortfall_dense"] >= 0.0
    assert ledger["o2_transfer_potential_bd"] >= 0.0
    assert ledger["o2_demand_bubble_limited"] <= ledger["o2_demand_bubble_raw"] + 1e-12
    assert ledger["o2_demand_dense_limited"] <= ledger["o2_demand_dense_raw"] + 1e-12
    assert ledger["h2o_demand_limited"] <= ledger["h2o_demand_raw"] + 1e-12
    assert ledger["reaction_numbering"]["authority"] == "hamel_1999_thesis"
    assert "R11" in ledger["net_molar_gas_source_by_thesis"]
    breakdown_sum = (
        ledger["net_molar_gas_source_vm"]
        + ledger["net_molar_gas_source_r5"]
        + ledger["net_molar_gas_source_r6"]
        + ledger["net_molar_gas_source_r7"]
        + ledger["net_molar_gas_source_r8"]
        + ledger["net_molar_gas_source_r9"]
        + ledger["net_molar_gas_source_r10"]
        + ledger["net_molar_gas_source_r11"]
        + ledger["net_molar_gas_source_r12"]
        + ledger["net_molar_gas_source_char"]
    )
    assert ledger["net_molar_gas_source_total"] == pytest.approx(breakdown_sum)
