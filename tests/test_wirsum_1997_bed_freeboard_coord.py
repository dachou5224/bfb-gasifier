"""阶段 F 后半：Lapple 旋风 + 内层气/T 与外层固体协调。"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from src.physics.wirsum_1997_cyclone import (
    CYCLONE_MODEL_CONSTANT,
    CYCLONE_MODEL_LAPPLE,
    CycloneGeometry,
    cyclone_eta,
    grade_efficiency_lapple,
    lapple_cut_size,
)
from src.solvers.wirsum_1997_bed_freeboard_coord import (
    OUTER_CRITERIA,
    build_stage_f2_payload,
    coordinate_once,
    inner_temperature,
    outer_criterion_coverage_coord,
)
from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry
from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED


def test_lapple_cut_is_half_efficiency() -> None:
    d50 = lapple_cut_size(
        mu=4.0e-5, inlet_width=0.011, n_turns=5.0, v_i=15.0, rho_p=1400.0, rho_g=0.8
    )
    assert d50 > 0.0
    assert grade_efficiency_lapple(d50, d50) == pytest.approx(0.5)
    assert grade_efficiency_lapple(10.0 * d50, d50) > 0.99
    assert grade_efficiency_lapple(0.1 * d50, d50) < 0.02


def test_lapple_eta_is_not_constant_0p95() -> None:
    geo = CycloneGeometry()
    eta, d50, name = cyclone_eta(
        mode=CYCLONE_MODEL_LAPPLE,
        d_p=2.0e-4,
        mu=4.0e-5,
        rho_p=1400.0,
        rho_g=0.80,
        N_tot=20.0,
        T=1100.0,
        P=2.5e6,
        geo=geo,
    )
    assert name == CYCLONE_MODEL_LAPPLE
    assert math.isfinite(d50)
    assert eta != pytest.approx(0.95, abs=1e-6)


def test_inner_T_drops_with_heat_loss() -> None:
    T0, r0 = inner_temperature(T_in=1100.0, N_tot=20.0, Q_wall=0.0)
    T1, r1 = inner_temperature(T_in=1100.0, N_tot=20.0, Q_wall=6.0e4)
    assert T0 == pytest.approx(1100.0)
    assert r0 == pytest.approx(0.0, abs=1e-9)
    assert T1 < T0
    assert r1 == pytest.approx(0.0, abs=1e-8)


def test_coordinate_inner_does_not_recompute_eject() -> None:
    from src.solvers.wirsum_1997_bed_freeboard_coord import inner_solve_gas_T

    src = inspect.getsource(inner_solve_gas_T)
    assert "eject_bed_top" not in src
    assert "simulate_freeboard" not in src


def test_repeat_outer_is_history_free() -> None:
    payload = build_stage_f2_payload()
    assert payload["repeat_outer_gap"]["eta"] == pytest.approx(0.0)
    assert payload["repeat_outer_gap"]["recycle_kg_s"] == pytest.approx(0.0)
    assert payload["inner_Q0"]["residual_gas_max"] == pytest.approx(0.0, abs=1e-12)
    assert payload["inner_Q0"]["residual_T"] == pytest.approx(0.0, abs=1e-8)
    assert payload["inner_Q0"]["N_fb"] == pytest.approx(
        payload["inner_Q0"]["N_d"] + payload["inner_Q0"]["N_b"]
    )


def test_heat_loss_lowers_eta_at_old_cut() -> None:
    payload = build_stage_f2_payload()
    assert payload["inner_Qwall"]["T"] < payload["inner_Q0"]["T"]
    assert payload["inner_Qwall"]["d50_m"] > payload["outer_Q0"]["d50_m"]
    assert payload["inner_Qwall"]["eta_at_Q0_cut"] < 0.5
    assert payload["reaction"]["model"] == "named_absent"
    assert payload["reaction"]["gap_mol_s"] == 0.0


def test_coverage_marks_temperature_not_furnace_energy() -> None:
    cov = outer_criterion_coverage_coord()
    assert set(cov) == set(OUTER_CRITERIA)
    assert cov["temperature"] == "covered"
    assert cov["energy"] == "not_covered"
    assert cov["pyrolysis_release"] == "not_covered"
    assert cov["element_balance"] == "not_covered"


def test_payload_skips_experiment_co() -> None:
    dumped = str(build_stage_f2_payload())
    assert "y_CO_target" not in dumped
    assert "experiment" not in dumped.lower()
    assert "fastox" not in dumped


def test_fine_leave_recycle_uses_lapple_not_constant() -> None:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    N_d_vec = np.array([0.2, 0.8])
    N_b_vec = np.array([0.6, 0.2])
    common = dict(
        geo=geo,
        geo_ex=geo_ex,
        cyclone_geo=CycloneGeometry(),
        N_d_in=1.0,
        N_b_in=0.0,
        N_tot_d=10.0,
        N_tot_b=10.0,
        N_d_vec=N_d_vec,
        N_b_vec=N_b_vec,
        m_char=0.32,
        m_ash=0.08,
        recycle=True,
        Q_wall=0.0,
        freeboard_model=FREEBOARD_MODEL_HAMEL_EXTENDED,
    )
    lap = coordinate_once(**common, cyclone_mode=CYCLONE_MODEL_LAPPLE)
    cst = coordinate_once(**common, cyclone_mode=CYCLONE_MODEL_CONSTANT)
    assert lap["streams"].fate == "reaches_next_cell"
    assert float(np.sum(lap["streams"].m_leave)) > 0.0
    assert lap["eta"] != pytest.approx(cst["eta"])
    np.testing.assert_allclose(lap["streams"].m_recycle, lap["eta"] * lap["streams"].m_leave)
