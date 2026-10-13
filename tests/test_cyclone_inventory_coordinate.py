"""旋风近零库存在入流出现后要用有尺度的线性坐标，不能留在对数列里。"""

import json
from pathlib import Path

import numpy as np

import src.solvers.wirsum_1997_named_f_e8 as e8
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    original_residual_from_tar_rows,
    outer_P_named,
)
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_L = 1.7628e-5  # [kg/s] s=0.700 失败点的反应后旋风入流
_K_AUF = 0.25  # [1/s]
_K_AB = 0.75  # [1/s]


def test_equilibrium_mass_matches_the_two_by_two_block() -> None:
    mass = e8.cyclone_equilibrium_mass_kg(_L, _K_AUF, _K_AB)
    assert abs(float(mass[0]) - 2.1696e-5) / 2.1696e-5 < 1.0e-4
    assert abs(float(mass[1]) - 5.4240e-6) / 5.4240e-6 < 1.0e-4
    matrix = np.array([
        [_K_AUF + _K_AB, -_K_AB],
        [-_K_AUF, _K_AUF + _K_AB],
    ])
    solved = np.linalg.solve(matrix, np.array([_L, 0.0]))
    assert np.allclose(mass, solved)


def test_near_zero_and_exact_zero_both_leave_the_log_coordinate() -> None:
    equilibrium = e8.cyclone_equilibrium_mass_kg(_L, _K_AUF, _K_AB)
    near = np.array([2.2241e-25, 1.8615e-27])
    exact = np.zeros(2)
    for mass in (near, exact):
        layout = e8.cyclone_reactivation_layout(
            mass, equilibrium,
            k_auf=_K_AUF, k_ab=_K_AB, inflow_kg_s=_L, log_span=0.5,
        )
        assert list(layout["index"]) == [37, 38]
        assert np.all(layout["lo"] >= -1.0e-12)
        assert np.all(layout["hi"] >= 1.0)
        reached = layout["base"] + layout["scale"] * ((equilibrium - np.maximum(mass, 0.0)) / layout["scale"])
        assert np.allclose(reached, equilibrium)


def test_inventory_already_near_its_inflow_stays_logarithmic() -> None:
    equilibrium = e8.cyclone_equilibrium_mass_kg(_L, _K_AUF, _K_AB)
    layout = e8.cyclone_reactivation_layout(
        equilibrium * np.exp(0.2), equilibrium,
        k_auf=_K_AUF, k_ab=_K_AB, inflow_kg_s=_L, log_span=0.5,
    )
    assert layout["index"].size == 0


def test_dust_below_the_solid_gate_is_not_rewritten() -> None:
    dust = np.array([2.2241e-25, 1.8615e-27])
    layout = e8.cyclone_reactivation_layout(
        dust, dust,
        k_auf=_K_AUF, k_ab=_K_AB, inflow_kg_s=1.0e-25, log_span=0.5,
    )
    assert layout["index"].size == 0


def test_vanished_inflow_can_return_a_visible_inventory_to_zero() -> None:
    layout = e8.cyclone_reactivation_layout(
        np.array([1.0e-2, 1.0e-20]), np.zeros(2),
        k_auf=_K_AUF, k_ab=_K_AB, inflow_kg_s=0.0, log_span=0.5,
    )
    assert list(layout["index"]) == [37]
    assert abs(float(layout["base"][0] + layout["scale"][0] * layout["lo"][0])) <= 1.0e-15


def _capture(x, q):
    cap: dict = {}
    stats: dict = {}
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
            q, x, stats=stats, log_span=0.5, temperature_radius_scale=4.0,
            couple_trajectory_jacobian=False, deadline_s=60.0,
        )
    finally:
        e8._gate_progress_step = original
    return cap, stats


def test_certified_root_keeps_logarithmic_cyclone_mass() -> None:
    cert = json.loads(
        Path("docs/audit_reports/equation_loop_2026-10-08/certified.json").read_text()
    )
    x = np.asarray(cert["x"], dtype=float)
    q = outer_P_named(x, q_on_wesseling1_path(float(cert["s"])))
    layout = e8._cyclone_reactivation_for_state(x, q, 0.5)
    assert layout["index"].size == 0


def test_saved_trial_cyclone_column_keeps_recycle_and_energy() -> None:
    trial = json.loads(
        Path("docs/audit_reports/equation_loop_2026-10-08/advance_0.700.json").read_text()
    )
    x = np.asarray(trial["trial_x"], dtype=float)
    nq = q_on_wesseling1_path(0.7)
    q = outer_P_named(x, nq)
    cap, stats = _capture(x, q)
    assert stats["cyclone_reactivation_ids"] == [37, 38]
    columns = []
    for col, index, unit in zip(
        stats["cyclone_reactivation_columns"],
        stats["cyclone_reactivation_ids"],
        stats["cyclone_reactivation_scale_kg"],
    ):
        got = original_residual_from_tar_rows(cap["matrix"][:, int(col)] * cap["scale"])
        h = 1.0e-4
        xp = x.copy()
        xp[int(index)] = float(x[int(index)]) + float(unit) * h
        fp = coupled_residual_named(xp, q)
        fm = coupled_residual_named(x, q)
        columns.append(got)
        assert np.allclose(got, (fp - fm) / h, rtol=1.0e-4, atol=1.0e-8)
    recycle = _K_AB * float(stats["cyclone_reactivation_scale_kg"][0])
    assert abs(float(columns[0][14] + columns[0][15]) - recycle) / recycle < 1.0e-3
    assert abs(float(columns[1][49])) > 0.1
