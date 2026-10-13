"""P0：R7/R8 灰基准、逐反应元素零和、E8 全非 T 列、Jv、线性步证据。不求 E7 链。"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from src.core.constants import Rg
from src.core.species import configure_tar_components_by_fuel, gibbs_molar, get_atom_count
from src.kinetics.gas_reactions import rate_R7, rate_R8, wgsr_equilibrium_constant
from src.kinetics.tar_reactions import (
    TAR_STOICH_R10,
    TAR_STOICH_R11,
    get_lumped_tar_stoichiometry,
    rate_R11_suspension,
)
from src.solvers.wirsum_1997_bed_freeboard import split_eject_by_fate
from src.solvers.wirsum_1997_named_f import (
    R10_FUEL_TYPE,
    coupled_residual_named,
    default_x0_named,
    homog_sources,
    outer_P_named,
    r10_nu_major,
    r11_nu_major,
)
from src.solvers.wirsum_1997_named_f_continuation import bound_feasible_direction
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import (
    _E8_STEP_EVIDENCE,
    _e8_newton_step,
    apply_e8_certified_rollback,
    freeze_e8_q,
    jacobian_named_e8,
    residual_named_e8,
)
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q
from src.solvers.wirsum_1997_named_f_r11 import freeze_e6_r10_r11_q
from src.solvers.wirsum_1997_named_f_r12 import NU_R12_MAJOR
from src.solvers.wirsum_1997_named_f_r12_r5 import NU_R5_MAJOR
from src.solvers.wirsum_1997_named_f_r12_r5_r8 import NU_R8_MAJOR
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import NU_R6_MAJOR
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6_r7 import NU_R7_MAJOR
from src.solvers.wirsum_1997_named_f_r234 import NU_R2_MAJOR, NU_R3_MAJOR, NU_R4_MAJOR
from src.solvers.wirsum_1997_stiffness_full_gas import IDX, MAJOR, N_MAJ

_SURR_ATOMS: dict[str, dict[str, float]] = {
    "C6H6": {"C": 6.0, "H": 6.0, "O": 0.0},
    "C10H8": {"C": 10.0, "H": 8.0, "O": 0.0},
    "C16H34": {"C": 16.0, "H": 34.0, "O": 0.0},
}


def _atoms_major(nu: np.ndarray, extra: dict[str, float] | None = None) -> dict[str, float]:
    out = {"C": 0.0, "H": 0.0, "O": 0.0}
    for name, coeff in zip(MAJOR, np.asarray(nu, dtype=np.float64)):
        for el in out:
            out[el] += float(coeff) * float(get_atom_count(name, el))
    if extra:
        for sp, coeff in extra.items():
            if sp == "C_solid":
                out["C"] += float(coeff)
                continue
            for el in out:
                out[el] += float(coeff) * float(get_atom_count(sp, el))
    return out


def _assert_zero(atoms: dict[str, float], tag: str) -> None:
    for el, val in atoms.items():
        assert abs(val) < 1.0e-12, f"{tag} {el}={val}"


def _e8_pack():
    nq8 = freeze_e8_q(Q_W=0.0)
    return nq8, freeze_e6_q(), freeze_e6_r10_q(), freeze_e6_r10_r11_q(), freeze_e7_q()


def test_r7_zero_without_ash_and_element_neutral() -> None:
    T = 1200.0
    C_ch4 = 0.05
    assert rate_R7(T, C_ch4, 0.2, 0.0, 0.0, rho_ash_kg_m3=0.0) == 0.0
    r = rate_R7(T, C_ch4, 0.2, 0.0, 0.0, rho_ash_kg_m3=75.0)
    assert r > 0.0
    np.testing.assert_allclose(2.0 * r, rate_R7(T, C_ch4, 0.2, 0.0, 0.0, rho_ash_kg_m3=150.0))
    _assert_zero(_atoms_major(NU_R7_MAJOR), "R7")


def test_r7_absent_from_bubble_and_ash_free_fb() -> None:
    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[IDX["CH4"]] = 0.2
    x[7 + IDX["CH4"]] = 0.2
    qT = nq.energy.species
    sd, sb = homog_sources(x[:39], qT, rate_scale=1.0, rho_ash=0.0, enabled=frozenset({"R7"}))
    np.testing.assert_allclose(sd, np.zeros(N_MAJ))
    np.testing.assert_allclose(sb, np.zeros(N_MAJ))
    sd75, sb75 = homog_sources(x[:39], qT, rate_scale=1.0, rho_ash=75.0, enabled=frozenset({"R7"}))
    assert float(sd75[IDX["CH4"]]) < 0.0
    np.testing.assert_allclose(sb75, np.zeros(N_MAJ))


def test_catalytic_rates_zero_without_ash() -> None:
    T = 1200.0
    P = 2.5e6
    assert rate_R8(T, P, 0.2, 0.2, 0.1, 0.1, rho_ash_kg_m3=0.0) == 0.0
    assert rate_R8(T, P, 0.2, 0.2, 0.1, 0.1, rho_ash_kg_m3=75.0) != 0.0
    assert rate_R11_suspension(T, 0.01, 0.0) == 0.0
    assert rate_R11_suspension(T, 0.01, 75.0) > 0.0
    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[IDX["CO"]] = 0.2
    x[IDX["H2O"]] = 0.2
    qT = nq.energy.species
    sd, sb = homog_sources(x[:39], qT, rate_scale=1.0, rho_ash=0.0, enabled=frozenset({"R8"}))
    np.testing.assert_allclose(sd, np.zeros(N_MAJ))
    np.testing.assert_allclose(sb, np.zeros(N_MAJ))


def test_r1_r12_element_zero_sum() -> None:
    for alpha in (0.5, 0.7, 1.0):
        nu = np.zeros(N_MAJ, dtype=np.float64)
        nu[IDX["O2"]] = -alpha
        nu[IDX["CO"]] = 2.0 * (1.0 - alpha)
        nu[IDX["CO2"]] = 2.0 * alpha - 1.0
        _assert_zero(_atoms_major(nu, {"C_solid": -1.0}), f"R1 alpha={alpha}")
    _assert_zero(_atoms_major(NU_R2_MAJOR, {"C_solid": -1.0}), "R2")
    _assert_zero(_atoms_major(NU_R3_MAJOR, {"C_solid": -1.0}), "R3")
    _assert_zero(_atoms_major(NU_R4_MAJOR, {"C_solid": -1.0}), "R4")
    _assert_zero(_atoms_major(NU_R5_MAJOR), "R5")
    _assert_zero(_atoms_major(NU_R6_MAJOR), "R6")
    _assert_zero(_atoms_major(NU_R7_MAJOR), "R7")
    _assert_zero(_atoms_major(NU_R8_MAJOR), "R8")
    _assert_zero(_atoms_major(NU_R12_MAJOR), "R12")


def test_r10_r11_element_zero_sum() -> None:
    configure_tar_components_by_fuel("coal")
    for surr, sto in TAR_STOICH_R10.items():
        atoms = {el: -n for el, n in _SURR_ATOMS[surr].items()}
        for sp, nu in sto.items():
            for el in atoms:
                atoms[el] += float(nu) * float(get_atom_count(sp, el))
        _assert_zero(atoms, f"R10 {surr}")
    for surr, sto in TAR_STOICH_R11.items():
        atoms = {el: -n for el, n in _SURR_ATOMS[surr].items()}
        for sp, nu in sto.items():
            for el in atoms:
                atoms[el] += float(nu) * float(get_atom_count(sp, el))
        _assert_zero(atoms, f"R11 {surr}")
    lumped10 = get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)
    extra10 = {k: v for k, v in lumped10.items() if k not in IDX}
    _assert_zero(_atoms_major(r10_nu_major(), extra10), "R10 lumped")
    lumped11 = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
    extra11 = {k: v for k, v in lumped11.items() if k not in IDX}
    _assert_zero(_atoms_major(r11_nu_major(), extra11), "R11 lumped")


def test_r8_keq_listed_not_ocr_and_not_fitted_to_1100k() -> None:
    """Eq. 5.47：采用式单列；不改系数、不用 1100 K 反推。"""
    T = 1000.0
    k_code = float(wgsr_equilibrium_constant(T))
    k_ocr = float(np.exp(-36.893 + 4.019 / T))
    k_fit_1100 = float(np.exp(-4019.0 / 1100.0 + 4019.0 / T))
    dG = (
        gibbs_molar("CO2", T) + gibbs_molar("H2", T)
        - gibbs_molar("CO", T) - gibbs_molar("H2O", T)
    )
    k_5_43 = float(np.exp(-dG / (Rg * T)))
    assert abs(k_code - float(np.exp(-3.6893 + 4019.0 / T))) < 1.0e-12
    assert k_ocr < 1.0e-10
    assert k_code > 0.1
    assert abs(k_code - k_fit_1100) > 1.0e-3
    assert 0.1 < k_code / k_5_43 < 10.0


def test_outer_p_rho_g_is_idempotent() -> None:
    nq8 = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq8)
    x[-1] = 800.0
    q1 = outer_P_named(x, nq8)
    q2 = outer_P_named(x, q1)
    np.testing.assert_allclose(q1.energy.species.geo_fb.rho_g, q2.energy.species.geo_fb.rho_g)
    x[-1] = 1420.0
    q1 = outer_P_named(x, nq8)
    q2 = outer_P_named(x, q1)
    np.testing.assert_allclose(q1.energy.species.geo_fb.rho_g, q2.energy.species.geo_fb.rho_g)


def test_undecided_fate_is_not_leave() -> None:
    m = np.array([0.4, 0.6], dtype=np.float64)
    ret, leave = split_eject_by_fate(m, "undecided")
    np.testing.assert_allclose(ret, 0.0)
    np.testing.assert_allclose(leave, 0.0)
    ret_r, leave_r = split_eject_by_fate(m, "returns")
    np.testing.assert_allclose(ret_r, m)
    np.testing.assert_allclose(leave_r, 0.0)
    ret_l, leave_l = split_eject_by_fate(m, "reaches_next_cell")
    np.testing.assert_allclose(ret_l, 0.0)
    np.testing.assert_allclose(leave_l, m)


def test_e8_non_t_jacobian_matches_frozen_p_fd() -> None:
    nq8, nq6, nq_r10, nq_r11, nq7 = _e8_pack()
    x0 = default_x0_named(nq8)
    sample = (0, IDX["CO"], 14, 15, 37, 38, 39, 40, 41)
    fd_rel = 1.0e-8
    fd_abs = 1.0e-10
    for T in (800.0, 1420.0):
        z = x0.copy()
        z[-1] = T
        nqP = outer_P_named(z, nq8)
        J = jacobian_named_e8(z, nq6, nq_r10, nq_r11, nq7, nq8)
        scale = np.maximum(np.abs(z), 0.05)
        f0 = coupled_residual_named(z, nqP, rate_scale=1.0)
        f_init = coupled_residual_named(z, nq8, rate_scale=1.0)
        for i in sample:
            eps = max(fd_abs, fd_rel * max(abs(float(z[i])), abs(float(scale[i]))))
            zp = z.copy()
            zp[i] = float(z[i]) + eps
            fd = (coupled_residual_named(zp, nqP, rate_scale=1.0) - f0) / eps
            np.testing.assert_allclose(
                J[:, i], fd, rtol=2.0e-6, atol=2.0e-8,
                err_msg=f"T={T} col={i} w_leave={nqP.energy.species.w_leave.tolist()}",
            )
        wP = np.asarray(nqP.energy.species.w_leave, dtype=np.float64)
        w0 = np.asarray(nq8.energy.species.w_leave, dtype=np.float64)
        if np.linalg.norm(wP - w0) > 1.0e-12:
            i = 15
            eps = max(fd_abs, fd_rel * max(abs(float(z[i])), abs(float(scale[i]))))
            zp = z.copy()
            zp[i] = float(z[i]) + eps
            fd_wrong = (coupled_residual_named(zp, nq8, rate_scale=1.0) - f_init) / eps
            assert np.linalg.norm(J[:, i] - fd_wrong) > 1.0e-8


def test_e8_jv_and_feasible_direction_fd_on_both_branches() -> None:
    nq8, nq6, nq_r10, nq_r11, nq7 = _e8_pack()
    x0 = default_x0_named(nq8)
    solid = np.array([14, 15, 37, 38, 39], dtype=int)
    for T in (800.0, 1420.0):
        z = x0.copy()
        z[-1] = T
        nqP = outer_P_named(z, nq8)
        J = jacobian_named_e8(z, nq6, nq_r10, nq_r11, nq7, nq8)
        v = np.zeros_like(z)
        v[solid] = np.array([1.0, 0.5, 0.2, 0.2, 0.1])
        v = v / np.linalg.norm(v)
        eps = 1.0e-7
        f0 = coupled_residual_named(z, nqP, rate_scale=1.0)
        fd = (coupled_residual_named(z + eps * v, nqP, rate_scale=1.0) - f0) / eps
        jv = J @ v
        np.testing.assert_allclose(jv[solid], fd[solid], rtol=0.05, atol=1.0e-6)
        dx = np.zeros_like(z)
        dx[solid] = -0.1 * np.maximum(z[solid], 0.05)
        p, _changed = bound_feasible_direction(z, dx)
        p[-1] = 0.0
        nrm = float(np.linalg.norm(p))
        assert nrm > 0.0
        p_hat = p / nrm
        fd_p = (
            coupled_residual_named(z + eps * p_hat, nqP, rate_scale=1.0) - f0
        ) / eps
        jp = J @ p_hat
        np.testing.assert_allclose(jp[solid], fd_p[solid], rtol=0.05, atol=1.0e-6)


def test_e8_t_column_fd_sensitivity_all_rows() -> None:
    nq8, nq6, nq_r10, nq_r11, nq7 = _e8_pack()
    z = default_x0_named(nq8)
    z[-1] = 1420.0
    J = jacobian_named_e8(z, nq6, nq_r10, nq_r11, nq7, nq8)
    f0 = residual_named_e8(z, nq8)
    fd1 = None
    rel = {}
    for hT in (0.1, 0.5, 1.0):
        zp = z.copy()
        zp[-1] = float(z[-1]) + hT
        fd = (residual_named_e8(zp, nq8) - f0) / hT
        assert np.all(np.isfinite(fd))
        if hT == 1.0:
            np.testing.assert_allclose(J[:, -1], fd, rtol=1.0e-9, atol=1.0e-9)
            fd1 = fd
        else:
            rel[hT] = float(np.linalg.norm(fd - J[:, -1]) / max(np.linalg.norm(J[:, -1]), 1.0e-12))
    assert fd1 is not None
    assert rel[0.1] >= 0.0
    assert rel[0.5] >= 0.0


def test_e8_linear_step_records_four_evidence_blocks() -> None:
    nq8 = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq8)
    n = int(x.size)
    J = np.eye(n, dtype=np.float64)
    F = np.ones(n, dtype=np.float64)
    _e8_newton_step(J, F, x)
    ev = dict(_E8_STEP_EVIDENCE)
    assert set(ev) >= {"raw_J", "modified_matrix", "feasible_direction", "rollback"}
    assert "raw_rank" in ev["raw_J"]
    assert "used_trace_diag" in ev["modified_matrix"]
    assert "p" in ev["feasible_direction"]
    assert ev["rollback"]["returned_state_source"] == "pending"
    last = SimpleNamespace(x=x + 1.0, F=F + 1.0)
    cert = {"endpoint_certified": False}
    cert0 = {"endpoint_certified": True}
    last2, cert2, src = apply_e8_certified_rollback(last, x, F, cert, cert0)
    assert src == "e7_certified_root"
    np.testing.assert_allclose(last2.x, x)
    np.testing.assert_allclose(last2.F, F)
    assert cert2 is cert0
    last3, _c3, src3 = apply_e8_certified_rollback(
        SimpleNamespace(x=x + 1.0, F=F + 1.0),
        x, F,
        {"endpoint_certified": True},
        {"endpoint_certified": True},
    )
    assert src3 == "newton"
    np.testing.assert_allclose(last3.x, x + 1.0)
