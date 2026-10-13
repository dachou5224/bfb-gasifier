"""正的微量氧离平衡很远时，只在联合步失去下降方向后做一次标量初值。"""

import json
from pathlib import Path

import numpy as np
import pytest

from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import (
    bracket_positive_scalar,
    oxygen_rows_needing_restart,
    solve_named_controlled,
)


def test_positive_mismatch_finds_one_root() -> None:
    found = bracket_positive_scalar(1.0e-30, lambda level: 1.0e-22 - float(level))
    assert found["status"] == "root"
    assert found["level"] == pytest.approx(1.0e-22, rel=1.0e-6)


def test_nonpositive_oxygen_is_not_lifted() -> None:
    found = bracket_positive_scalar(0.0, lambda level: 1.0)
    assert found["status"] == "nonpositive"
    assert "level" not in found


def test_missing_bracket_does_not_invent_a_floor() -> None:
    found = bracket_positive_scalar(1.0e-20, lambda level: 1.0)
    assert found["status"] == "no_bracket"
    assert "level" not in found


def test_multiple_crossings_are_declined() -> None:
    def row(level: float) -> float:
        value = float(np.log10(level))
        if value < -12.0:
            return -1.0
        if value < 4.0:
            return 1.0
        return -1.0

    found = bracket_positive_scalar(1.0e-8, row, max_samples=9)
    assert found["status"] == "multiple_roots"
    assert "level" not in found


def test_balanced_oxygen_does_not_request_a_restart() -> None:
    residual = np.zeros(50)
    residual[0] = 2.0e-6
    residual[5] = 1.0e-12
    residual[21] = 1.0e-12
    assert oxygen_rows_needing_restart(residual, residual.copy()) == []


def test_saved_trace_oxygen_stall_passes_after_one_restart() -> None:
    saved = json.loads(
        Path("docs/audit_reports/equation_loop_2026-10-08/advance_0.952.json").read_text()
    )
    x0 = np.asarray(saved["trial_x"], dtype=float)
    out = solve_named_controlled(
        q_on_wesseling1_path(float(saved["s"])),
        x0,
        deadline_s=30.0,
        temperature_radius_scale=4.0,
    )
    assert out["passed"] is True
    assert out["stats"]["oxygen_restart"]["status"] == "accepted"
    q = outer_P_named(out["x"], q_on_wesseling1_path(float(saved["s"])))
    residual = coupled_residual_named(out["x"], q)
    gap = external_imbalance(out["x"], q)["gap_mol_s"]
    assert max(abs(float(gap[key])) for key in ("C", "H", "O")) < 1.0e-4
    assert float(np.linalg.norm(residual[5:6])) < 1.0e-6
    assert not np.array_equal(out["x"], x0)
