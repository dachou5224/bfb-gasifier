from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from src.core.cell import Cell, S_ASH, S_CHAR, SolidProps
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import configure_tar_components_by_fuel
from src.core.species import N_GAS
from src.solvers.cell_solver import _build_residual_scales, evaluate_cell_state
from src.solvers.vorabrechnung import refresh_cell_vorabrechnung


def test_cell_residuals_pass_rate_multiplier_through_reaction_stage(monkeypatch):
    cell = Cell()
    order: list[str] = []
    seen: list[float] = []

    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: order.append("hyd"))
    monkeypatch.setattr(cell, "calc_exchange", lambda: order.append("exchange"))

    def _calc_reactions(rate_multiplier: float = 1.0):
        order.append("reactions")
        seen.append(rate_multiplier)

    monkeypatch.setattr(cell, "calc_reactions", _calc_reactions)
    monkeypatch.setattr(cell, "calc_gas_balance", lambda: order.append("gas") or np.zeros(2 * N_GAS))
    monkeypatch.setattr(cell, "calc_solid_balance", lambda: order.append("solid") or np.zeros_like(cell.m_solid))
    monkeypatch.setattr(cell, "calc_energy_balance", lambda: order.append("energy") or 0.0)

    cell.residuals(rate_multiplier=0.25)

    assert seen == [0.25]
    assert order == ["hyd", "exchange", "reactions", "gas", "solid", "energy"]


def test_drying_pyrolysis_source_overwrites_solid_sink_instead_of_accumulating(monkeypatch):
    cell = Cell()
    cell.R_solid[:, :] = 5.0
    solid_sink = np.ones_like(cell.R_solid)
    gas_source = np.full(N_GAS, 2.0)

    monkeypatch.setattr(
        "src.core.cell.calc_drying_pyrolysis_sources",
        lambda **_: SimpleNamespace(gas_source=gas_source, solid_sink=solid_sink),
    )

    out = cell._calc_drying_pyrolysis_gas_source(0.1)

    np.testing.assert_allclose(out, gas_source)
    np.testing.assert_allclose(cell.R_solid, solid_sink)


def test_calc_reactions_uses_frozen_vorabrechnung_cache_when_valid(monkeypatch):
    cell = Cell()
    cell.u_mf = 1.0
    cell._vm_cache_valid = True
    cell._vm_gas_source_cache[:] = 3.0
    cell._vm_solid_sink_cache[:, :] = 4.0

    captured: dict[str, np.ndarray] = {}

    monkeypatch.setattr(
        cell,
        "_calc_drying_pyrolysis_gas_source",
        lambda tau: (_ for _ in ()).throw(AssertionError("dynamic source recompute should stay outside inner NR")),
    )
    monkeypatch.setattr(cell, "_calc_char_surface_area_per_class", lambda: np.zeros(cell.R_solid.shape[0]))
    monkeypatch.setattr(cell, "_compute_char_conversion", lambda: 0.0)
    monkeypatch.setattr(cell, "_catalyst_bulk_density", lambda: 0.0)
    monkeypatch.setattr(cell, "calc_minor_species_gibbs", lambda: None)

    def _fake_build_reaction_sources(**kwargs):
        captured["gas_src_vm"] = np.array(kwargs["gas_src_vm"], copy=True)
        captured["solid_sink_vm"] = np.array(kwargs["solid_sink_vm"], copy=True)
        return SimpleNamespace(
            R_gas_b=np.zeros(N_GAS),
            R_gas_d=np.zeros(N_GAS),
            R_solid=np.zeros_like(cell.R_solid),
        )

    monkeypatch.setattr("src.core.cell.build_reaction_sources", _fake_build_reaction_sources)

    cell.calc_reactions()

    np.testing.assert_allclose(captured["gas_src_vm"], cell._vm_gas_source_cache)
    np.testing.assert_allclose(captured["solid_sink_vm"], cell._vm_solid_sink_cache)


@pytest.mark.parametrize("cell_type", ["cyclone", "return_leg"])
def test_calc_reactions_is_disabled_for_side_block_cells(cell_type: str):
    cell = Cell()
    cell.cell_type = cell_type
    cell.R_gas_b.fill(1.0)
    cell.R_gas_d.fill(-2.0)
    cell.R_solid.fill(3.0)

    cell.calc_reactions()

    np.testing.assert_allclose(cell.R_gas_b, 0.0)
    np.testing.assert_allclose(cell.R_gas_d, 0.0)
    np.testing.assert_allclose(cell.R_solid, 0.0)


def test_cell_hydrodynamics_reuses_local_state_cache(monkeypatch):
    cell = Cell()
    cell.T = 1000.0
    cell.P = 2.5e6
    cell.N_d[0] = 1.0

    calls = {"n": 0}

    def _fake_hydro(**kwargs):
        calls["n"] += 1
        return SimpleNamespace(
            u_mf=0.1,
            u_b=0.2,
            d_b=0.3,
            eps_b=0.4,
            eps_d=0.6,
            eps_d_voidage=0.5,
            K_bd=1.2,
            V_b=0.7,
            V_d=0.8,
            u0=0.9,
            u_d=0.11,
            n_rz=2.7,
        )

    monkeypatch.setattr("src.core.cell.calc_cell_hydrodynamics", _fake_hydro)

    cell.calc_hydrodynamics()
    cell.calc_hydrodynamics()
    assert calls["n"] == 1

    # Inlet changes do not alter the local hydrodynamics state.
    cell.N_d_in[0] = 5.0
    cell.calc_hydrodynamics()
    assert calls["n"] == 1

    # Local primitive-state changes must invalidate the cache.
    cell.N_d[0] = 2.0
    cell.calc_hydrodynamics()
    assert calls["n"] == 2


def test_cell_residuals_use_frozen_vorabrechnung_hydrodynamics_in_inner_nr(monkeypatch):
    cell = Cell()
    cell.u_mf = 1.1
    cell.u_b = 2.2
    cell.d_b = 0.3
    cell.eps_b = 0.4
    cell.eps_d = 0.6
    cell.eps_d_voidage = 0.5
    cell.K_bd = 3.3
    cell.V_b = 0.7
    cell.V_d = 0.8
    cell.u0 = 0.9
    cell.u_d = 0.12
    cell.n_rz = 4.4
    cell.K_solid_auf[:, S_CHAR] = 5.5
    cell.K_solid_ab[:, S_CHAR] = 6.6
    cell.compute_vorabrechnung(0.1)
    cell.enable_inner_nr_vorabrechnung_freeze()

    # Disturb live hydrodynamics to confirm residuals restores the frozen snapshot.
    cell.u_mf = -1.0
    cell.K_bd = -1.0
    cell.K_solid_auf[:, S_CHAR] = -1.0
    cell.K_solid_ab[:, S_CHAR] = -1.0

    monkeypatch.setattr(
        cell,
        "calc_hydrodynamics",
        lambda: (_ for _ in ()).throw(AssertionError("inner NR should use frozen Vorabrechnung hydrodynamics")),
    )
    monkeypatch.setattr(cell, "calc_exchange", lambda: None)
    monkeypatch.setattr(cell, "calc_reactions", lambda rate_multiplier=1.0: None)
    monkeypatch.setattr(cell, "calc_gas_balance", lambda: np.zeros(2 * N_GAS))
    monkeypatch.setattr(cell, "calc_solid_balance", lambda: np.zeros_like(cell.m_solid))
    monkeypatch.setattr(cell, "calc_energy_balance", lambda: 0.0)

    cell.residuals()

    assert cell.u_mf == 1.1
    assert cell.K_bd == 3.3
    np.testing.assert_allclose(cell.K_solid_auf[:, S_CHAR], 5.5)
    np.testing.assert_allclose(cell.K_solid_ab[:, S_CHAR], 6.6)


def test_cell_thermo_cache_reuses_local_bundle_for_exchange_and_reactions(monkeypatch):
    cell = Cell()
    cell.T = 1000.0
    cell.P = 2.5e6
    cell.N_d[0] = 1.0
    cell.N_b[1] = 2.0
    cell.K_bd = 1.0
    cell.V_b = 0.5
    cell.u_mf = 1.0
    cell._vm_cache_valid = True

    calls = {"n": 0}

    def _fake_diffusivity(T, P):
        calls["n"] += 1
        return 1.23e-4

    monkeypatch.setattr("src.core.cell.gas_diffusivity_correlation", _fake_diffusivity)
    monkeypatch.setattr(cell, "_calc_char_surface_area_per_class", lambda: np.zeros(cell.R_solid.shape[0]))
    monkeypatch.setattr(cell, "_compute_char_conversion", lambda: 0.0)
    monkeypatch.setattr(cell, "_catalyst_bulk_density", lambda: 0.0)
    monkeypatch.setattr(cell, "calc_minor_species_gibbs", lambda: None)
    monkeypatch.setattr(
        "src.core.cell.build_reaction_sources",
        lambda **kwargs: SimpleNamespace(
            R_gas_b=np.zeros(N_GAS),
            R_gas_d=np.zeros(N_GAS),
            R_solid=np.zeros_like(cell.R_solid),
        ),
    )

    cell.calc_exchange()
    cell.calc_reactions()
    assert calls["n"] == 1

    # Inlet-only changes do not alter the local thermo state.
    cell.N_d_in[0] = 5.0
    cell.calc_exchange()
    cell.calc_reactions()
    assert calls["n"] == 1

    # Local gas-state changes must invalidate the thermo cache.
    cell.N_d[0] = 2.0
    cell.calc_exchange()
    cell.calc_reactions()
    assert calls["n"] == 2


def test_char_surface_area_uses_holdup_state_distribution_under_holdup_transport():
    cell = Cell(
        solid=SolidProps(
            n_size_classes=2,
            d_p_classes=np.array([0.001, 0.002], dtype=float),
            mass_fractions=np.array([0.5, 0.5], dtype=float),
        )
    )
    cell.solid_state_model = "holdup_transport"
    cell.V_d = 0.1
    cell.solid.rho_s = 1000.0
    cell.solid.eps_mf = 0.5
    cell.m_solid[:, :] = 0.0
    cell.m_solid[:, S_CHAR] = np.array([2.0, 1.0], dtype=float)
    cell.m_solid[:, S_ASH] = np.array([1.0, 0.0], dtype=float)

    areas = cell._calc_char_surface_area_per_class()

    m_bed = cell.solid.rho_s * (1.0 - cell.solid.eps_mf) * cell.V_d
    expected_inventory = m_bed * np.array([2.0, 1.0], dtype=float) / 4.0
    expected_areas = expected_inventory * 6.0 / (cell.solid.rho_s * cell.solid.d_p_classes)
    np.testing.assert_allclose(areas, expected_areas)


def test_solid_outflow_rates_use_holdup_times_frozen_transport_coefficients():
    cell = Cell()
    cell.solid_state_model = "holdup_transport"
    cell.m_solid[0, :] = np.array([2.0, 1.0, 0.5, 4.0], dtype=float)
    cell.K_solid_auf.fill(0.2)
    cell.K_solid_ab.fill(0.3)

    np.testing.assert_allclose(cell._solid_upflow_rates(), 0.2 * cell.m_solid)
    np.testing.assert_allclose(cell._solid_downflow_rates(), 0.3 * cell.m_solid)
    np.testing.assert_allclose(cell._solid_outflow_rates(), 0.5 * cell.m_solid)


def test_cell_residual_scales_include_frozen_transport_inflows_for_holdup_mode():
    configure_tar_components_by_fuel("coal")
    cell = Cell()
    cell.solid_state_model = "holdup_transport"
    cell.m_solid_zu[0, S_CHAR] = 1.0
    cell.m_solid_rez[0, S_ASH] = 2.0
    cell.m_solid_in[0, S_CHAR] = 3.0
    cell.m_solid_auf_in[0, S_CHAR] = 4.0
    cell.m_solid_ab_in[0, S_ASH] = 5.0

    scales = _build_residual_scales(cell)

    assert scales[2 * N_GAS] == 15.0


def test_evaluate_cell_state_uses_current_residual_vector(monkeypatch):
    cell = Cell()
    residual = np.linspace(1.0, float(cell._work_res.size), cell._work_res.size)

    monkeypatch.setattr(cell, "residuals", lambda rate_multiplier=1.0: residual.copy())

    metrics = evaluate_cell_state(cell, scales=np.ones_like(residual))

    assert metrics["residual"] == float(np.max(np.abs(residual)))
    assert metrics["rms_scaled"] == float(np.sqrt(np.mean(residual**2)))


def test_reactor_reports_globally_reevaluated_gs_metrics():
    reactor = Reactor(ReactorConfig(n_cells=2, allow_legacy_gs=True))
    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=1, tol_global=1e-4, solver="gauss_seidel")


def test_reactor_runs_upper_pair_corrective_sweep_for_high_upper_pair_rms():
    reactor = Reactor(ReactorConfig(n_cells=3, allow_legacy_gs=True))
    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=1, tol_global=1e-4, solver="gauss_seidel")


def test_reactor_caps_bottom_full_stiff_solve_to_reference_window():
    reactor = Reactor(ReactorConfig(n_cells=2, allow_legacy_gs=True))
    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=1, tol_global=1e-4, solver="gauss_seidel")


def test_reactor_deduplicates_bottom_full_temperature_caps_when_default_is_wide():
    reactor = Reactor(ReactorConfig(n_cells=1, allow_legacy_gs=True))
    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=1, tol_global=1e-4, solver="gauss_seidel")


def test_reactor_stops_after_first_clear_degradation_from_best_iter():
    reactor = Reactor(ReactorConfig(n_cells=1, allow_legacy_gs=True))
    with pytest.raises(ValueError, match="Only solver='global_nr'"):
        reactor.solve(max_global_iter=5, tol_global=1e-4, solver="gauss_seidel")


def test_refresh_cell_vorabrechnung_runs_hydrodynamics_then_cache_update(monkeypatch):
    cell = Cell()
    order: list[str] = []

    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: order.append("hyd"))
    monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: order.append("vorab"))

    refresh_cell_vorabrechnung(cell, force=False)

    assert order == ["hyd", "vorab"]


def test_refresh_cell_vorabrechnung_force_invalidates_before_recompute(monkeypatch):
    cell = Cell()
    order: list[str] = []

    monkeypatch.setattr(cell, "calc_hydrodynamics", lambda: order.append("hyd"))
    monkeypatch.setattr(cell, "invalidate_vorabrechnung_cache", lambda: order.append("invalidate"))
    monkeypatch.setattr(cell, "compute_vorabrechnung", lambda tau: order.append("vorab"))

    refresh_cell_vorabrechnung(cell, force=True)

    assert order == ["hyd", "invalidate", "vorab"]
