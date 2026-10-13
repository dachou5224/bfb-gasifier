"""阶段 E：最小联立阶梯。延拓不是物理模型；未实现 rung 必须显式失败。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.solvers.wirsum_1997_stiffness_ladder import (
    ExchangeGeometry,
    assert_rung_implemented,
    build_stage_e_payload,
    classify_newton_stop,
    default_r12_inlets,
    directional_derivative_error,
    exact_tracer_exchange,
    implemented_rungs,
    r12_exchange_jacobian,
    r12_exchange_residual,
    tracer_exchange_jacobian,
    tracer_exchange_residual,
)
from src.solvers.wirsum_1997_strict import WirsumStrictSolveResult


def test_later_rungs_are_implemented_and_named() -> None:
    for name in ("full_gas", "char", "energy"):
        assert name in implemented_rungs()
        assert_rung_implemented(name)


def test_module_does_not_restore_fastox() -> None:
    import src.solvers.wirsum_1997_stiffness_ladder as mod

    src = inspect.getsource(mod)
    assert "fastox" in src
    assert "snap_trace_o2" not in src
    assert "physical_convergence_on_continuation" in src


def test_small_step_is_not_physical_convergence() -> None:
    fake = WirsumStrictSolveResult(
        x=np.zeros(2),
        F=np.ones(2),
        converged=True,
        iterations=8,
        failure_reason="",
    )
    report = classify_newton_stop(fake, model_is_full_rate=True)
    assert report["failure_reason"] == "small_step"
    assert report["physical_convergence"] is False


def test_line_search_is_not_physical_convergence() -> None:
    fake = WirsumStrictSolveResult(
        x=np.zeros(2),
        F=np.ones(2),
        converged=False,
        iterations=3,
        failure_reason="line_search",
    )
    report = classify_newton_stop(fake, model_is_full_rate=True)
    assert report["physical_convergence"] is False
    assert report["failure_reason"] == "line_search"


def test_continuation_scale_is_not_physical_even_if_newton_converges() -> None:
    fake = WirsumStrictSolveResult(
        x=np.zeros(2),
        F=np.zeros(2),
        converged=True,
        iterations=1,
        failure_reason="",
    )
    report = classify_newton_stop(fake, model_is_full_rate=False)
    assert report["converged_newton"] is True
    assert report["physical_convergence"] is False


def test_transport_and_exchange_payload_linear_errors() -> None:
    payload = build_stage_e_payload()
    tr = payload["transport"]
    assert tr["converged_newton"] is True
    assert tr["physical_convergence"] is True
    assert tr["backward_error_at_root"] < 1.0e-10
    assert tr["newton_step_backward_error_at_x0"] < 1.0e-10
    assert tr["jacobian_rank"] == 3

    ex = payload["exchange"]
    assert ex["raw_l2"]["converged_newton"] is True
    assert ex["scaled_rms"]["converged_newton"] is True
    assert ex["max_abs_raw_vs_exact"] < 1.0e-12
    assert ex["F_xx_max_abs"] == pytest.approx(0.0, abs=0.0)
    assert ex["dir_deriv_eps"]["1e-3"] < 1.0e-10
    assert ex["dir_deriv_eps"]["1e-5"] < 1.0e-8


def test_tracer_jacobian_column_permutation_is_consistent() -> None:
    geo = ExchangeGeometry()
    x = np.array([0.4, 0.6], dtype=np.float64)

    def res(z: np.ndarray) -> np.ndarray:
        return tracer_exchange_residual(
            z, geo=geo, N_d_in=1.0, N_b_in=0.0, N_tot_d=10.0, N_tot_b=10.0
        )

    def jac(z: np.ndarray) -> np.ndarray:
        return tracer_exchange_jacobian(z, geo=geo, N_tot_d=10.0, N_tot_b=10.0)

    perm = np.array([1, 0])
    P = np.eye(2)[perm]
    J = jac(x)
    y = x[perm]

    def res_perm(z: np.ndarray) -> np.ndarray:
        return res(z[perm])[perm]

    J_perm_fd = np.zeros((2, 2), dtype=np.float64)
    eps = 1.0e-8
    f0 = res_perm(y)
    for k in range(2):
        z = y.copy()
        z[k] += eps
        J_perm_fd[:, k] = (res_perm(z) - f0) / eps
    np.testing.assert_allclose(J_perm_fd, P @ J @ P.T, rtol=1e-5, atol=1e-7)


def test_r12_residual_uses_production_rate() -> None:
    from src.kinetics.gas_reactions import rate_R12
    from src.solvers.wirsum_1997_stiffness_ladder import SPECIES, _Y_DUMMY, _phase_conc, unpack_r12_state

    geo = ExchangeGeometry()
    N_in_d, N_in_b = default_r12_inlets()
    x = np.array([0.25, 1.0e-4, 0.05, 10.0, 0.05, 0.10, 0.02, 10.0], dtype=np.float64)
    N_d, N_b = unpack_r12_state(x)
    C_d = _phase_conc(N_d, geo)
    F = r12_exchange_residual(x, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=1.0)
    r = rate_R12(geo.T, float(C_d[0]), float(C_d[1]), geo.P, _Y_DUMMY, None)
    N_ex = geo.K_bd * geo.V_b * (_phase_conc(N_b, geo) - C_d)
    expected_h2_d = N_in_d[0] - N_d[0] + N_ex[0] - r * geo.V_d
    assert F[0] == pytest.approx(expected_h2_d, rel=1e-12)
    assert SPECIES[0] == "H2"


def test_r12_jacobian_directional_error_shrinks() -> None:
    geo = ExchangeGeometry()
    N_in_d, N_in_b = default_r12_inlets()
    x = np.array([0.25, 1.0e-4, 0.05, 10.0, 0.05, 0.10, 0.02, 10.0], dtype=np.float64)

    def res(z: np.ndarray) -> np.ndarray:
        return r12_exchange_residual(z, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=1.0)

    def jac(z: np.ndarray) -> np.ndarray:
        return r12_exchange_jacobian(z, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=1.0)

    e3 = directional_derivative_error(res, jac, x, np.ones(8), 1e-3)
    e5 = directional_derivative_error(res, jac, x, np.ones(8), 1e-5)
    assert e5 < e3


def test_r12_continuation_is_named_and_full_rate_is_distinct() -> None:
    payload = build_stage_e_payload()
    cont = payload["r12"]["continuation"]
    assert cont["scale_0"]["model_is_full_rate"] is False
    assert cont["scale_0"]["physical_convergence"] is False
    assert cont["scale_1e+00"]["model_is_full_rate"] is True
    if not cont["scale_1e+00"]["converged_newton"]:
        assert cont["scale_1e+00"]["failure_reason"] in {
            "line_search",
            "singular_jacobian",
            "max_iter",
        }
        assert cont["scale_1e+00"]["physical_convergence"] is False
    uni = payload["r12"]["uniqueness_at_scale_1e-8"]
    if uni["max_abs_x_spread"] is not None:
        assert uni["max_abs_x_spread"] < 1.0e-6
    full_u = payload["r12"]["uniqueness_at_full_rate"]
    if full_u["from_continuation"]["physical_convergence"]:
        assert full_u["from_perturbation"]["physical_convergence"] is True
        assert full_u["max_abs_x_vs_continuation"] < 1.0e-10
        assert full_u["direct_separated_at_full_rate"]["physical_convergence"] is False
        assert full_u["direct_separated_at_full_rate"]["failure_reason"] == "line_search"


def test_payload_does_not_touch_experiment_co() -> None:
    dumped = str(build_stage_e_payload())
    assert "y_CO_target" not in dumped
    assert "experiment" not in dumped.lower()


def test_r6_is_partial_oxidation_not_complete() -> None:
    from src.solvers.wirsum_1997_stiffness_full_gas import (
        IDX,
        NU_R6,
        _gas_source,
    )
    from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry
    from src.core.species import GAS_SPECIES_INDEX, N_GAS

    geo = ExchangeGeometry(T=1400.0)
    C = np.zeros(7, dtype=np.float64)
    C[IDX["CH4"]] = 0.2
    C[IDX["O2"]] = 0.4
    C[IDX["N2"]] = 8.0
    y = np.zeros(N_GAS, dtype=np.float64)
    ntot = float(np.sum(C))
    y[GAS_SPECIES_INDEX["CH4"]] = C[IDX["CH4"]] / ntot
    y[GAS_SPECIES_INDEX["O2"]] = C[IDX["O2"]] / ntot
    y[GAS_SPECIES_INDEX["N2"]] = C[IDX["N2"]] / ntot
    src_d = _gas_source(C, y, geo=geo, V=1.0, phase="dense", rate_scale=1.0, rho_ash=0.0)
    src_b = _gas_source(C, y, geo=geo, V=1.0, phase="bubble", rate_scale=1.0, rho_ash=75.0)
    np.testing.assert_allclose(NU_R6[IDX["CO"]], 1.0)
    np.testing.assert_allclose(NU_R6[IDX["CO2"]], 0.0)
    np.testing.assert_allclose(NU_R6[IDX["O2"]], -1.5)
    np.testing.assert_allclose(NU_R6[IDX["H2O"]], 2.0)
    # 仅 CH4/O2：R6 应产 CO 不产 CO2；气泡 R8 灰催化为 0
    assert src_d[IDX["CO"]] > 0.0
    assert src_d[IDX["CO2"]] == pytest.approx(0.0, abs=1e-14)
    assert src_d[IDX["O2"]] / src_d[IDX["CH4"]] == pytest.approx(1.5, rel=1e-8)
    assert src_d[IDX["H2O"]] / (-src_d[IDX["CH4"]]) == pytest.approx(2.0, rel=1e-8)
    np.testing.assert_allclose(src_b, src_d)


def test_full_gas_elements_close_when_newton_converges() -> None:
    from src.solvers.wirsum_1997_stiffness_full_gas import build_full_gas_char_energy_section

    later = build_full_gas_char_energy_section()
    gas = later["full_gas"]
    assert gas["reactions"] == ["R5", "R6", "R7", "R8", "R12"]
    assert "R9" in gas["not"]
    assert gas["continuation"]["scale_0"]["physical_convergence"] is False
    assert "last_successful_scale" in gas["continuation"]
    if gas["continuation"]["scale_1e+00"]["physical_convergence"]:
        for el, gap in gas["element_gap_mol_s"].items():
            assert gap < 1.0e-8, el
    else:
        assert gas["continuation"]["scale_1e+00"]["failure_reason"] in {
            "line_search",
            "singular_jacobian",
            "max_iter",
            "small_step",
        }
        assert gas["continuation"]["scale_1e+00"]["physical_convergence"] is False
    char = later["char"]
    assert char["continuation"]["scale_0"]["physical_convergence"] is False
    energy = later["energy"]
    assert energy["solid_sensible"] == "omitted_named_not_hamel_cp_integral"
    if energy["continuation"]["scale_1e+00"]["physical_convergence"]:
        assert energy["T_K"] is not None
        assert 300.0 < float(energy["T_K"]) < 2500.0
    else:
        assert energy["continuation"]["scale_1e+00"]["physical_convergence"] is False


def test_energy_kernel_explicit_scales_are_not_physical_convergence() -> None:
    from src.solvers.wirsum_1997_stiffness_full_gas import (
        ENERGY_F_SCALE,
        ENERGY_X_SCALE,
        N_ENERGY,
        default_full_gas_inlets,
        energy_jacobian_quality_at,
    )

    N_d, N_b = default_full_gas_inlets()
    M = 0.010 / 0.080  # [kg]
    T = 1100.0  # [K]
    x = np.concatenate([N_d, N_b, np.array([M, T], dtype=np.float64)])
    inert = energy_jacobian_quality_at(x, rate_scale=0.0)
    full = energy_jacobian_quality_at(x, rate_scale=1.0)
    assert inert["n_unknowns"] == full["n_unknowns"] == N_ENERGY == 16
    assert full["T_retained"] is True
    assert full["fictitious_Q_W"] is False
    assert full["physical_convergence"] is False
    assert full["merit_kind"] == "raw_l2"
    assert full["linear_solve_quality_checked"] is True
    assert ENERGY_X_SCALE.shape == ENERGY_F_SCALE.shape == (16,)
    # 无反应时数值满秩；完整速率下原始秩判断可失败，均衡后也不等于已物理收敛
    assert inert["raw_rank"] == 16
    assert inert["equilibrated_rank"] == 16
    assert full["raw_rank"] <= 16
    assert full["raw_condition"] > inert["raw_condition"]
    dE_an = float(full["analytic_dE_dT_W_per_K"])
    dE_fd = float(full["fd_dE_dT_W_per_K"])
    assert dE_an < 0.0
    assert abs(dE_fd - dE_an) / max(abs(dE_an), 1.0) < 0.05
    assert full["energy_row_M_W_per_kg"] == pytest.approx(0.0, abs=1.0e-8)


def test_full_gas_module_skips_fastox() -> None:
    import src.solvers.wirsum_1997_stiffness_full_gas as mod

    src = inspect.getsource(mod)
    assert "fastox" in src
    assert "snap_trace_o2" not in src


def test_exact_tracer_matches_residual_zero() -> None:
    geo = ExchangeGeometry()
    exact = exact_tracer_exchange(
        geo=geo, N_d_in=1.0, N_b_in=0.0, N_tot_d=10.0, N_tot_b=10.0
    )
    F = tracer_exchange_residual(
        exact, geo=geo, N_d_in=1.0, N_b_in=0.0, N_tot_d=10.0, N_tot_b=10.0
    )
    np.testing.assert_allclose(F, 0.0, atol=1e-14)
