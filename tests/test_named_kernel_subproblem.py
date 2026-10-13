"""局部子问题的无量纲保护和切线右端。不跑参数延拓。"""

import numpy as np
import pytest

from src.solvers.named_kernel_sources import named_kernel_source_files
from src.solvers.wirsum_1997_named_f import (
    apply_tar_reaction_transport_rows,
    linear_tar_row_map,
    original_residual_from_tar_rows,
)
from src.solvers.wirsum_1997_named_f_e8 import (
    _gate_progress_step,
    material_protection_margin,
)


def test_material_protection_rejects_a_block_ten_times_over_the_gate() -> None:
    assert material_protection_margin(1.0e-6) == pytest.approx(0.0)
    assert material_protection_margin(1.0e-5) == pytest.approx(-99.0)
    assert material_protection_margin(1.0e-5) < -1.0e-8
    old_margin = 1.0e-12 - (1.0e-5) ** 2
    assert old_margin >= -1.0e-8


def test_protected_block_stays_inside_the_gate_when_the_quadratic_wants_to_leave() -> None:
    n = 50
    matrix = np.zeros((n, 1), dtype=np.float64)
    matrix[0, 0] = -1.0e-2
    matrix[14, 0] = 1.0
    entrance = np.zeros(n, dtype=np.float64)
    entrance[0] = 1.0e-7
    entrance[14] = 1.0
    quality: dict = {}
    step, predicted, status = _gate_progress_step(
        matrix,
        entrance,
        np.ones(n, dtype=np.float64),
        np.array([-2.0]),
        np.array([2.0]),
        quality,
    )
    assert status in (1, 2)
    assert abs(float(predicted[0])) <= 1.0e-6 * (1.0 + 1.0e-4)
    assert abs(float(step[0]) + 1.0) > 0.5
    assert float(quality["max_nd_violation"]) <= 1.0e-8


def test_expired_deadline_does_not_run_the_inner_optimizer() -> None:
    import time

    n = 50
    matrix = np.zeros((n, 1), dtype=np.float64)
    matrix[0, 0] = -1.0e-2
    matrix[14, 0] = 1.0
    entrance = np.zeros(n, dtype=np.float64)
    entrance[0] = 1.0e-7
    entrance[14] = 1.0
    quality: dict = {}
    t0 = time.monotonic()
    _step, _predicted, status = _gate_progress_step(
        matrix,
        entrance,
        np.ones(n, dtype=np.float64),
        np.array([-2.0]),
        np.array([2.0]),
        quality,
        deadline_at=time.monotonic() - 1.0,
    )
    assert time.monotonic() - t0 < 0.5
    assert quality["deadline_hit"] is True
    assert quality["solved"] is False
    assert status in (0, 2)


def test_linear_tar_map_is_homogeneous_and_inverts() -> None:
    vector = np.linspace(-0.2, 0.3, 50)
    assert np.all(linear_tar_row_map(np.zeros(50)) == 0.0)
    other = np.linspace(0.1, -0.4, 50)
    mapped = linear_tar_row_map(vector + other)
    assert mapped == pytest.approx(linear_tar_row_map(vector) + linear_tar_row_map(other))
    assert original_residual_from_tar_rows(linear_tar_row_map(vector)) == pytest.approx(vector)
    assert linear_tar_row_map(original_residual_from_tar_rows(vector)) == pytest.approx(vector)


def test_state_assembly_is_not_the_linear_map_of_an_increment() -> None:
    from dataclasses import replace

    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    base = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    q = outer_P_named(base, nq)
    increment = np.zeros(50)
    increment[40] = 1.0e-6
    increment[41] = -2.0e-6
    state_rows = apply_tar_reaction_transport_rows(base, q, increment)
    zero_rows = apply_tar_reaction_transport_rows(base, q, np.zeros(50))
    linear_rows = linear_tar_row_map(increment)
    assert state_rows[41] == pytest.approx(zero_rows[41])
    assert linear_rows[41] != pytest.approx(zero_rows[41])


def test_source_fingerprint_includes_execution_dependencies() -> None:
    files = named_kernel_source_files()
    assert "src/solvers/wirsum_1997_named_f_p2.py" in files
    assert "src/physics/wirsum_1997_freeboard.py" in files
    assert "src/thermal/devolatilization.py" in files
    assert "data/validation_cases.json" in files
