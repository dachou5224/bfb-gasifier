from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_balances import calc_energy_balance_H_in, calc_energy_balance_residual, calc_solid_enthalpy_flow
from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
    SOLID_ENTHALPY_MODE_HAMEL_EQ28,
    SOLID_ENTHALPY_MODE_LEGACY_CP_MIX,
)
from src.core.species import configure_tar_components_by_fuel, cp_ash, cp_char, cp_sand, GAS_SPECIES_INDEX


def _matched_energy_kwargs(*, heat_loss_frac: float = 0.0) -> dict:
    gas = np.zeros(11, dtype=np.float64)
    gas[4] = 2.0
    solid = np.array([[0.5, 0.0, 0.0, 0.0]], dtype=np.float64)
    return dict(
        N_b_in=gas * 0.3,
        N_d_in=gas * 0.7,
        T_in_gas=1100.0,
        N_zu_b=np.zeros_like(gas),
        N_zu_d=np.zeros_like(gas),
        T_zu_gas=1100.0,
        N_rez_b=np.zeros_like(gas),
        N_rez_d=np.zeros_like(gas),
        T_rez_gas=1100.0,
        m_solid_rez=np.zeros_like(solid),
        T_rez_solid=1100.0,
        m_solid_in=np.zeros_like(solid),
        T_in_solid=1100.0,
        m_solid_zu=solid,
        T_zu_solid=1100.0,
        N_b=gas * 0.3,
        N_d=gas * 0.7,
        m_solid=solid,
        T=1100.0,
        heat_loss_frac=heat_loss_frac,
        ash_dry_wt=10.0,
        VM_daf=50.0,
        h_f_dry=-5.0e5,
        h_cache={},
    )


def test_solid_enthalpy_eq28_uses_component_cp_not_mix() -> None:
    m = np.array([[1.0, 0.0, 0.0, 0.0]], dtype=np.float64)
    common = dict(
        ash_dry_wt=10.0,
        VM_daf=50.0,
        h_f_dry=-5.0e5,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
    )
    h_eq28 = calc_solid_enthalpy_flow(
        m,
        1200.0,
        solid_enthalpy_mode=SOLID_ENTHALPY_MODE_HAMEL_EQ28,
        **common,
    )
    h_legacy = calc_solid_enthalpy_flow(
        m,
        1200.0,
        solid_enthalpy_mode=SOLID_ENTHALPY_MODE_LEGACY_CP_MIX,
        **common,
    )
    assert h_eq28 != h_legacy


def test_energy_balance_hamel_eq27_explicit_heat_loss() -> None:
    configure_tar_components_by_fuel("coal")
    kw = _matched_energy_kwargs()
    res0 = calc_energy_balance_residual(**kw, heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT)
    assert res0 == 0.0
    q_loss = 1.0e5
    res = calc_energy_balance_residual(
        **kw,
        heat_loss_mode=HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
        Q_wall_W=q_loss,
    )
    assert res == pytest.approx(-q_loss)


def test_inlet_fraction_scales_linearly_with_heat_loss_frac() -> None:
    configure_tar_components_by_fuel("coal")
    gas = np.zeros(11, dtype=np.float64)
    gas[GAS_SPECIES_INDEX["H2"]] = 20.0
    solid = np.zeros((1, 4), dtype=np.float64)
    base = dict(
        N_b_in=gas,
        N_d_in=np.zeros_like(gas),
        T_in_gas=1100.0,
        N_zu_b=np.zeros_like(gas),
        N_zu_d=np.zeros_like(gas),
        T_zu_gas=1100.0,
        N_rez_b=np.zeros_like(gas),
        N_rez_d=np.zeros_like(gas),
        T_rez_gas=1100.0,
        m_solid_rez=solid,
        T_rez_solid=1100.0,
        m_solid_in=solid,
        T_in_solid=1100.0,
        m_solid_zu=solid,
        T_zu_solid=1100.0,
        N_b=gas * 0.5,
        N_d=gas * 0.5,
        m_solid=solid,
        T=1100.0,
        ash_dry_wt=10.0,
        VM_daf=50.0,
        h_f_dry=-5.0e5,
        h_cache={},
    )
    H_in = calc_energy_balance_H_in(
        N_b_in=base["N_b_in"],
        N_d_in=base["N_d_in"],
        T_in_gas=base["T_in_gas"],
        N_zu_b=base["N_zu_b"],
        N_zu_d=base["N_zu_d"],
        T_zu_gas=base["T_zu_gas"],
        N_rez_b=base["N_rez_b"],
        N_rez_d=base["N_rez_d"],
        T_rez_gas=base["T_rez_gas"],
        m_solid_rez=base["m_solid_rez"],
        T_rez_solid=base["T_rez_solid"],
        m_solid_in=base["m_solid_in"],
        T_in_solid=base["T_in_solid"],
        m_solid_zu=base["m_solid_zu"],
        T_zu_solid=base["T_zu_solid"],
        ash_dry_wt=base["ash_dry_wt"],
        VM_daf=base["VM_daf"],
        h_f_dry=base["h_f_dry"],
        h_cache=base["h_cache"],
    )
    assert H_in > 0.0
    r0 = calc_energy_balance_residual(**{**base, "heat_loss_frac": 0.0}, heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION)
    r1 = calc_energy_balance_residual(**{**base, "heat_loss_frac": 0.1}, heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION)
    r2 = calc_energy_balance_residual(**{**base, "heat_loss_frac": 0.2}, heat_loss_mode=HEAT_LOSS_MODE_INLET_FRACTION)
    assert (r1 - r0) / (r2 - r0) == pytest.approx(0.5)


def test_project_freeboard_temperatures_clears_energy_ghost() -> None:
    """轴向焓投影把自由板过热鬼峰拉回 Eq.2.7 能量根，并刷新下游 T_in。"""
    from src.core.cell import Cell
    from src.core.freeboard_bridge import project_explicit_freeboard_temperatures_to_energy_closure
    from src.core.energy_balance_defaults import HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT

    configure_tar_components_by_fuel("coal")

    class _FakeReactor:
        def __init__(self, cells, freeboard_cells):
            self.cells = cells
            self.freeboard_cells = freeboard_cells
            self.connectivity_graph = None

        def _solver_cell_registry(self):
            return {}

    bed = Cell()
    bed.cell_type = "bed"
    bed.T = 1000.0
    gas = np.zeros(11, dtype=np.float64)
    gas[GAS_SPECIES_INDEX["N2"]] = 20.0
    gas[GAS_SPECIES_INDEX["CO"]] = 5.0
    gas[GAS_SPECIES_INDEX["H2O"]] = 5.0
    bed.N_d[:] = gas
    bed.N_b[:] = 0.0

    fb0 = Cell()
    fb0.cell_type = "freeboard"
    fb0.solid_state_model = "freeboard_closure"
    fb0.heat_loss_mode = HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
    fb0.Q_wall_W = 1.5e5
    fb0.N_d_in[:] = gas
    fb0.N_d[:] = gas
    fb0.N_b_in[:] = 0.0
    fb0.N_b[:] = 0.0
    fb0.T_in_gas = 1000.0
    fb0.T = 1000.0
    # Force an energy ghost above the enthalpy root.
    fb0.T = 1300.0
    assert abs(fb0.calc_energy_balance()) > 1.0e4

    fb1 = Cell()
    fb1.cell_type = "freeboard"
    fb1.solid_state_model = "freeboard_closure"
    fb1.heat_loss_mode = HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
    fb1.Q_wall_W = 1.5e5
    fb1.N_d[:] = gas
    fb1.N_b[:] = 0.0
    fb1.T = 1200.0

    reactor = _FakeReactor([bed], [fb0, fb1])
    diag = project_explicit_freeboard_temperatures_to_energy_closure(reactor)
    assert diag["n_projected"] >= 1
    assert abs(fb0.calc_energy_balance()) < 1.0e3
    assert fb0.T < 1300.0
    assert fb1.T_in_gas == pytest.approx(fb0.T, abs=1.0)


def test_freeboard_zero_ash_holdup_does_not_use_dense_r8_fallback() -> None:
    """自由板上段灰持量为 0 时 ρ_ash=0，禁止密相灰密度回退放大 R8。"""
    from src.core.cell import Cell, S_ASH

    cell = Cell()
    cell.cell_type = "freeboard"
    cell.solid_state_model = "freeboard_closure"
    cell.V_d = 0.3  # [m³]
    cell.m_solid[:, S_ASH] = 0.0
    assert cell._ash_bulk_density() == pytest.approx(0.0)

    cell.m_solid[0, S_ASH] = 0.03  # [kg]
    assert cell._ash_bulk_density() == pytest.approx(0.03 / 0.3)

    bed = Cell()
    bed.cell_type = "bed"
    bed.V_d = 0.3
    bed.m_solid[:, S_ASH] = 0.0
    assert bed._ash_bulk_density() > 1.0  # legacy dense fallback still active


def test_freeboard_closure_energy_uses_solid_inlet_passthrough() -> None:
    """freeboard_closure 能量出口用入口固体透传，避免 m_in≠K·m 造成 300 K 假根。"""
    from src.core.cell import Cell
    from src.core.energy_balance_defaults import HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT

    configure_tar_components_by_fuel("coal")
    cell = Cell()
    cell.cell_type = "freeboard"
    cell.solid_state_model = "freeboard_closure"
    cell.heat_loss_mode = HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
    cell.heat_loss_frac = 0.0
    cell.Q_wall_W = 0.0
    cell.Q_heat_exchanger_W = 0.0
    cell.T = 951.4
    cell.T_in_gas = 951.4
    cell.T_in_solid = 951.4
    cell.T_zu_gas = 951.4
    cell.T_zu_solid = 951.4
    cell.T_rez_gas = 951.4
    cell.T_rez_solid = 951.4

    gas = np.zeros(11, dtype=np.float64)
    gas[GAS_SPECIES_INDEX["N2"]] = 20.0
    gas[GAS_SPECIES_INDEX["H2O"]] = 5.0
    cell.N_d_in[:] = gas
    cell.N_d[:] = gas
    cell.N_b_in[:] = 0.0
    cell.N_b[:] = 0.0

    # Inlet solid mass larger than K·m holdup outflow (trajectory mismatch).
    cell.m_solid_in[:] = 0.0
    cell.m_solid_in[0, 0] = 0.40  # char
    cell.m_solid[:] = 0.0
    cell.m_solid[0, 0] = 0.05
    cell.K_solid_auf[:] = 0.0
    cell.K_solid_ab[:] = 0.0
    cell.K_solid_ab[0, 0] = 3.6  # downflow ≈ 0.18 kg/s ≠ 0.40

    outlet = cell._solid_energy_outlet_rates()
    transport = cell._solid_outflow_rates()
    np.testing.assert_allclose(outlet, cell.m_solid_in)
    assert float(np.sum(outlet)) == pytest.approx(0.40)
    assert float(np.sum(transport)) == pytest.approx(0.18)

    # Matched T / gas: residual should be ~0 (no structural solid enthalpy hole).
    assert cell.calc_energy_balance() == pytest.approx(0.0, abs=1.0e3)


def test_resolve_inlet_fraction_calibration_falls_back_to_abs_bridge(monkeypatch) -> None:
    from src.core.energy_balance_defaults import resolve_inlet_fraction_calibration_W

    def _fake_summarize(_reactor, *, use_abs_H_in: bool = False) -> dict[str, float]:
        if not use_abs_H_in:
            return {"bed_W": 0.0, "freeboard_W": 0.0, "total_W": 0.0}
        return {"bed_W": 100.0, "freeboard_W": 20.0, "total_W": 120.0}

    monkeypatch.setattr(
        "src.core.energy_balance_defaults.summarize_reactor_inlet_fraction_heat_loss_W",
        _fake_summarize,
    )
    implied, mode = resolve_inlet_fraction_calibration_W(object())
    assert mode == "abs_H_in_audit_bridge"
    assert implied["total_W"] == pytest.approx(120.0)
