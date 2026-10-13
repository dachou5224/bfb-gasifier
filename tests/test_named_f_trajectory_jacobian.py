"""校正器里影响飞行时间的列用完整 P 差商，不再冻住 tau。"""

import json
from pathlib import Path

import numpy as np
import pytest

import src.solvers.wirsum_1997_named_f_e8 as e8
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    original_residual_from_tar_rows,
    outer_P_named,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_SOLID = (14, 15, 37, 38, 39, 46, 47)


def _capture(x, q, *, coupled: bool):
    cap: dict = {}
    original = e8._gate_progress_step

    def grab(matrix, entrance_row, scale, lo, hi, quality=None, **kwargs):
        cap["matrix"] = np.asarray(matrix, dtype=float).copy()
        cap["scale"] = np.asarray(scale, dtype=float).copy()
        return (
            np.zeros(matrix.shape[1]),
            original_residual_from_tar_rows(entrance_row * scale),
            0,
        )

    e8._gate_progress_step = grab
    try:
        e8._e8_log_scaled_correct(
            q, x, stats={}, temperature_radius_scale=4.0, deadline_s=30.0,
            couple_trajectory_jacobian=coupled,
        )
    finally:
        e8._gate_progress_step = original
    z0 = e8._seed_downstream_o2(q, np.asarray(x, dtype=float).copy())
    log_ids = [int(i) for i in e8._e8_positive_log_ids(z0) if i not in range(40, 46)]
    return cap, z0, log_ids


def _column(cap, log_ids, state_index: int) -> np.ndarray:
    j = log_ids.index(int(state_index))
    return original_residual_from_tar_rows(cap["matrix"][:, j] * cap["scale"])


def _full_fd(z0, nq, index: int, h: float) -> np.ndarray:
    xp = np.asarray(z0, dtype=float).copy()
    xm = xp.copy()
    xp[index] *= np.exp(h)
    xm[index] *= np.exp(-h)
    fp = coupled_residual_named(xp, outer_P_named(xp, nq))
    fm = coupled_residual_named(xm, outer_P_named(xm, nq))
    return (fp - fm) / (2.0 * h)


def _load(path: str):
    cert = json.loads(Path(path).read_text())
    x = np.asarray(cert["x"], dtype=float)
    nq = q_on_wesseling1_path(float(cert["s"]))
    q = outer_P_named(x, nq)
    return x, nq, q


@pytest.fixture(scope="module")
def root():
    return _load("docs/audit_reports/equation_loop_2026-10-08/certified.json")


def test_return_root_temperature_column_flips_without_the_trajectory() -> None:
    x, nq, q = _load("docs/audit_reports/s0698_review_2026-10-10/root_s0698.json")
    coupled, z0, ids = _capture(x, q, coupled=True)
    frozen, _z, ids_f = _capture(x, q, coupled=False)
    assert ids == ids_f
    full = _full_fd(z0, nq, 49, 1.0e-5)
    full_fine = _full_fd(z0, nq, 49, 1.0e-6)
    assert full[47] > 0.0
    assert full_fine[47] > 0.0
    assert abs(full[47] - full_fine[47]) / abs(full[47]) < 0.01
    got = _column(coupled, ids, 49)
    old = _column(frozen, ids_f, 49)
    assert old[47] < 0.0
    assert got[47] == pytest.approx(full[47], rel=0.02, abs=1.0e-4)
    solid = np.asarray(_SOLID, dtype=int)
    assert float(np.max(np.abs(got[solid] - full[solid]))) < 1.0e-4


def test_freeboard_temperature_column_matches_the_full_difference(root) -> None:
    x, nq, q = root
    coupled, z0, ids = _capture(x, q, coupled=True)
    full = _full_fd(z0, nq, 49, 1.0e-5)
    got = _column(coupled, ids, 49)
    assert got[47] == pytest.approx(full[47], rel=0.02, abs=1.0e-4)
    solid = np.asarray(_SOLID, dtype=int)
    assert float(np.max(np.abs(got[solid] - full[solid]))) < 1.0e-4


@pytest.mark.parametrize("index", (48, 16, 19))
def test_trajectory_columns_match_full_difference_on_the_solid_block(root, index: int) -> None:
    x, nq, q = root
    coupled, z0, ids = _capture(x, q, coupled=True)
    full = _full_fd(z0, nq, index, 1.0e-5 if index >= 48 else 1.0e-4)
    got = _column(coupled, ids, index)
    solid = np.asarray(_SOLID, dtype=int)
    assert float(np.max(np.abs(got[solid] - full[solid]))) < 2.0e-4
    assert got[47] == pytest.approx(float(full[47]), rel=0.05, abs=2.0e-4)
