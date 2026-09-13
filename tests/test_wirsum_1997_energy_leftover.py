"""holdup Startwert 之后，bed0 −5 MW 可加总；最大桶是 VM，不是再给 VM 加 K。"""

from __future__ import annotations

import inspect

from src.solvers.wirsum_1997_energy_leftover import (
    attribute_eq27_leftover,
    build_energy_leftover_payload,
)


def test_leftover_buckets_sum_to_residual() -> None:
    row = {
        "cell": 0,
        "residual_W": -10.0,
        "H_in_parts_W": {
            "H_gas_axial_in": 0.0,
            "H_gas_zu": -3.0,
            "H_gas_rez": 0.0,
            "H_solid_rez": 0.0,
            "H_solid_axial_in": 0.0,
            "H_solid_auf_in": 0.0,
            "H_solid_ab_in": 0.0,
            "H_solid_zu": -8.0,
        },
        "H_zu_by_comp_W": {"char": -5.0, "vm": -3.0, "moisture": 0.0, "ash": 0.0},
        "H_solid_out_by_comp_W": {"char": -5.0, "vm": 0.0, "moisture": 0.0, "ash": 0.0},
        "H_gas_out_W": -1.0,
        "Q_loss_W": 3.0,
        "holdup_zero": False,
        "vm_in_not_in_energy_out": True,
    }
    # E = H_in - H_out - Q = (-3-8) - (-1-5) - 3 = -11 + 6 - 3 = -8
    row["residual_W"] = -8.0
    got = attribute_eq27_leftover(row)
    assert got["sum_matches_residual"] is True
    assert got["buckets_W"]["vm_H_zu_minus_H_sout"] == -3.0
    assert got["buckets_W"]["minus_Q"] == -3.0
    assert got["buckets_W"]["gas_H_in_minus_H_out"] == -2.0


def test_lu_leftover_after_holdup_is_mostly_vm_then_gas_and_q() -> None:
    payload = build_energy_leftover_payload()
    bed0 = payload["per_cell"][0]
    assert bed0["sum_matches_residual"] is True
    assert bed0["holdup_zero"] is False
    assert abs(float(bed0["residual_W"] or 0.0)) > 1.0e6
    assert bed0["dominant_bucket"] == "vm_H_zu_minus_H_sout"
    assert float(bed0["buckets_W"]["vm_H_zu_minus_H_sout"] or 0.0) < -1.0e6
    assert float(bed0["abs_share"]["vm_H_zu_minus_H_sout"] or 0.0) > 0.4
    assert abs(float(bed0["buckets_W"]["gas_H_in_minus_H_out"] or 0.0)) > 1.0e5
    assert float(bed0["buckets_W"]["minus_Q"] or 0.0) < 0.0


def test_leftover_does_not_add_k_vm_or_experiment_co() -> None:
    import src.solvers.wirsum_1997_energy_leftover as mod

    source = inspect.getsource(mod)
    assert "K_vm" in source
    assert "vm_in_eq27_outlet" in source
    dumped = str(build_energy_leftover_payload())
    assert "experiment" not in dumped.lower()
    assert "y_CO_target" not in dumped
