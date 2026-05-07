from __future__ import annotations

from tests.validation_case_utils import strict_validation_gate


def test_strict_validation_gate_accepts_fully_converged_low_rms_result():
    ok, reasons = strict_validation_gate(
        {
            "converged_fully": True,
            "rms_scaled_final": 0.05,
        }
    )
    assert ok is True
    assert reasons == []


def test_strict_validation_gate_rejects_unconverged_result_even_if_kpis_exist():
    ok, reasons = strict_validation_gate(
        {
            "converged": False,
            "converged_outer": False,
            "converged_fully": False,
            "rms_scaled_final": 0.05,
            "T_profile": [1100.0],
            "carbon_conv": 0.9,
        }
    )
    assert ok is False
    assert "solver_not_converged" in reasons


def test_strict_validation_gate_rejects_high_rms_even_when_marked_converged():
    ok, reasons = strict_validation_gate(
        {
            "converged_fully": True,
            "rms_scaled_final": 0.25,
        },
        rms_max=0.15,
    )
    assert ok is False
    assert "rms_scaled_final>0.150" in reasons
