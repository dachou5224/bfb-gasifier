"""P2：有限水力 P(x;b)。不改生产默认，不更新旋风 K，不宣称炉温。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import (
    default_x0_named,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q, jacobian_named_e8, residual_named_e8
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q
from src.solvers.wirsum_1997_named_f_r11 import freeze_e6_r10_r11_q
from src.solvers.wirsum_1997_named_f_p2 import (
    eq372_launch_nodes,
    finite_hydro_from_x,
    gas_velocity_and_density,
)


def test_reference_inputs_immutable_and_p_has_no_history() -> None:
    nq = freeze_e8_q(Q_W=0.0)
    assert nq.reference is not None
    x = default_x0_named(nq)
    x[-1] = 800.0
    q1 = outer_P_named(x, nq)
    q2 = outer_P_named(x, q1)
    assert q1.reference is not None and q2.reference is not None
    assert q1.reference.rho_g_ref == pytest.approx(nq.reference.rho_g_ref)
    assert q2.reference.rho_g_ref == pytest.approx(nq.reference.rho_g_ref)
    np.testing.assert_allclose(q1.energy.species.geo_fb.rho_g, q2.energy.species.geo_fb.rho_g)
    np.testing.assert_allclose(q1.energy.species.geo_fb.u_g, q2.energy.species.geo_fb.u_g)
    np.testing.assert_allclose(q1.energy.species.geo_ex.K_bd, q2.energy.species.geo_ex.K_bd)
    assert q1.energy.species.k == nq.energy.species.k
    assert q1.energy.species.k.model == "frozen_context"
    x[-1] = 1420.0
    q3 = outer_P_named(x, q1)
    q4 = outer_P_named(x, nq)
    np.testing.assert_allclose(q3.energy.species.geo_fb.rho_g, q4.energy.species.geo_fb.rho_g)
    np.testing.assert_allclose(q3.energy.species.geo_fb.u_g, q4.energy.species.geo_fb.u_g)


def test_total_gas_flow_changes_u_g() -> None:
    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[-1] = 1200.0
    q0 = outer_P_named(x, nq)
    x2 = x.copy()
    x2[0:14] *= 1.5
    q1 = outer_P_named(x2, nq)
    assert q1.energy.species.geo_fb.u_g > q0.energy.species.geo_fb.u_g
    b = nq.reference
    assert b is not None
    u0, rho0, n0, _mw0 = gas_velocity_and_density(
        x[0:7], x[7:14], T=1200.0, P=b.P, area=b.area
    )
    u1, _rho1, n1, _mw1 = gas_velocity_and_density(
        x2[0:7], x2[7:14], T=1200.0, P=b.P, area=b.area
    )
    assert n1 == pytest.approx(1.5 * n0)
    assert u1 == pytest.approx(1.5 * u0)
    assert q0.energy.species.geo_fb.u_g == pytest.approx(u0)
    assert q1.energy.species.geo_fb.u_g == pytest.approx(u1)
    assert abs(rho1 := q1.energy.species.geo_fb.rho_g) > 0.0
    _ = rho1
    assert q0.hydro is not None and q1.hydro is not None
    assert q1.hydro.K_bd != pytest.approx(float(nq.energy.species.geo_ex.K_bd), abs=0.0)


def test_eq372_nodes_normalized_and_refine() -> None:
    u, w = eq372_launch_nodes(1.2, 0.6, 11)
    assert u.shape == w.shape
    assert float(np.min(u)) >= 0.0
    assert float(np.sum(w)) == pytest.approx(1.0)
    u2, w2 = eq372_launch_nodes(1.2, 0.6, 41)
    mean1 = float(np.dot(u, w))
    mean2 = float(np.dot(u2, w2))
    assert abs(mean2 - mean1) < 0.05
    u1, w1 = eq372_launch_nodes(1.2, 0.6, 1)
    assert u1[0] == pytest.approx(1.2)
    assert w1[0] == pytest.approx(1.0)


def test_launch_distribution_refinement_converges() -> None:
    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[-1] = 1200.0
    b = nq.reference
    assert b is not None
    rows = []
    for n_bins in (5, 11, 21):
        snap, _geo, _kbd = finite_hydro_from_x(x, b, n_bins=n_bins, dt=2.0e-3, t_max=5.0)
        rows.append(snap)
        np.testing.assert_allclose(
            snap.w_leave + snap.w_return + snap.w_undecided, 1.0, atol=1.0e-12
        )
    d_leave = np.linalg.norm(rows[2].w_leave - rows[1].w_leave)
    d_leave_coarse = np.linalg.norm(rows[1].w_leave - rows[0].w_leave)
    assert d_leave <= d_leave_coarse + 0.05
    d_tau = float(np.linalg.norm(rows[2].tau_s - rows[1].tau_s))
    d_tau_c = float(np.linalg.norm(rows[1].tau_s - rows[0].tau_s))
    assert d_tau <= d_tau_c + 0.5
    snap_dt, _, _ = finite_hydro_from_x(x, b, n_bins=11, dt=1.0e-3, t_max=5.0)
    assert np.linalg.norm(snap_dt.w_leave - rows[1].w_leave) < 0.15
    snap_t, _, _ = finite_hydro_from_x(x, b, n_bins=11, dt=2.0e-3, t_max=8.0)
    assert float(np.sum(snap_t.w_undecided)) <= float(np.sum(rows[1].w_undecided)) + 1.0e-12


def test_undecided_mass_is_not_leave_in_distribution() -> None:
    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[-1] = 800.0
    q = outer_P_named(x, nq)
    assert q.hydro is not None
    assert float(q.hydro.w_undecided[1] + q.energy.species.w_return[1]) == pytest.approx(
        1.0 - float(q.energy.species.w_leave[1]), abs=1.0e-12
    )
    assert float(q.energy.species.w_leave[1]) < 0.5


def test_f_and_jacobian_use_same_q() -> None:
    nq8 = freeze_e8_q(Q_W=0.0)
    nq6 = freeze_e6_q(Q_W=0.0)
    nq_r10 = freeze_e6_r10_q(Q_W=0.0)
    nq_r11 = freeze_e6_r10_r11_q(Q_W=0.0)
    nq7 = freeze_e7_q(Q_W=0.0)
    x = default_x0_named(nq8)
    x[-1] = 1200.0
    qP = outer_P_named(x, nq8)
    F_h = residual_named_e8(x, nq8)
    from src.solvers.wirsum_1997_named_f import coupled_residual_named

    F_same = coupled_residual_named(x, qP)
    np.testing.assert_allclose(F_h, F_same)
    J = jacobian_named_e8(x, nq6, nq_r10, nq_r11, nq7, nq8)
    i = 15  # M_coarse：避开 R12 大数行
    eps = max(1.0e-10, 1.0e-8 * max(abs(float(x[i])), 0.05))
    xp = x.copy()
    xp[i] = float(x[i]) + eps
    fd = (coupled_residual_named(xp, qP) - F_same) / eps
    solid = np.array([14, 15, 37, 38, 39], dtype=int)
    np.testing.assert_allclose(J[solid, i], fd[solid], rtol=0.05, atol=1.0e-6)
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.vorab_outer_source_policy != "hamel_outer_refresh"
