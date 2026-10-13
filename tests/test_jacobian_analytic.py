from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import Cell
from src.core.cell_balances import S_CHAR
from src.solvers.jacobian_analytic import (
    JACOBIAN_DERIVATIVE_SOURCES,
    analytic_packed_derivative,
    calc_side_temperature_closure_cp_dot,
    inventory_self_derivative_unscaled,
    side_temperature_reference_K,
)
from src.solvers.nr_indexing import has_temperature_var, n_gas_var, temperature_offset


def test_jacobian_derivative_sources_cover_side_and_gas_groups() -> None:
    groups = {entry.residual_group for entry in JACOBIAN_DERIVATIVE_SOURCES}
    assert "energy_side_closure" in groups
    assert "gas_dense" in groups
    assert "gas_axial_inflow" in groups
    assert "side_recycle_coupling" in groups


def test_axial_gas_inflow_coupling_is_plus_one_on_immediate_downstream() -> None:
    from src.solvers.jacobian_analytic import inventory_inflow_coupling_derivative

    row = Cell()
    row.cell_type = "bed"
    col = Cell()
    col.cell_type = "bed"
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=0, col_local=0, row_block=1, col_block=0
    ) == 1.0
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=0, col_local=0, row_block=0, col_block=0
    ) is None
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=0, col_local=1, row_block=1, col_block=0
    ) is None


def test_legacy_solid_band_inflow_uses_lower_and_upper_fractions() -> None:
    from src.solvers.jacobian_analytic import inventory_inflow_coupling_derivative
    from src.solvers.nr_indexing import n_gas_var

    row = Cell()
    row.cell_type = "bed"
    row.solid_state_model = "legacy_stream"
    row._nr_solid_lower_inlet_frac = 0.35
    col = Cell()
    col.cell_type = "bed"
    col.solid_state_model = "legacy_stream"
    sol_off = n_gas_var(row)
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=sol_off, col_local=sol_off, row_block=1, col_block=0
    ) == pytest.approx(0.35)
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=sol_off, col_local=sol_off, row_block=0, col_block=1
    ) == pytest.approx(0.65)


def test_holdup_solid_band_inflow_uses_k_auf_and_k_ab() -> None:
    from src.solvers.jacobian_analytic import inventory_inflow_coupling_derivative

    row = Cell()
    row.cell_type = "bed"
    row.solid_state_model = "holdup_transport"
    col = Cell()
    col.cell_type = "bed"
    col.solid_state_model = "holdup_transport"
    col.K_solid_auf[0, S_CHAR] = 0.25
    col.K_solid_ab[0, S_CHAR] = 0.15
    col.m_solid[0, S_CHAR] = 1.0
    sol_off = n_gas_var(row)
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=sol_off, col_local=sol_off, row_block=1, col_block=0
    ) == pytest.approx(0.25)
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=sol_off, col_local=sol_off, row_block=0, col_block=1
    ) == pytest.approx(0.15)
    col.m_solid[0, S_CHAR] = 0.0
    assert inventory_inflow_coupling_derivative(
        row, col, row_local=sol_off, col_local=sol_off, row_block=1, col_block=0
    ) == pytest.approx(0.0)


def test_holdup_transport_hybrid_jacobian_matches_fd_on_lu_bed_chain() -> None:
    from tests.global_nr_fixtures import run_holdup_generate_initial_x0, setup_holdup_transport_reactor
    from src.solvers.global_nr_solver import (
        build_equation_scales,
        build_jacobian_fd,
        global_residual,
        pack_reactor,
    )

    def _assemble(strategy: str):
        reactor, cfg, t_profile = setup_holdup_transport_reactor(n_cells=3)
        run_holdup_generate_initial_x0(reactor, cfg, t_profile)
        reactor._apply_all_bc_for_nr()
        cells = reactor.cells
        x = pack_reactor(cells)
        f0 = global_residual(x, cells, reactor._apply_all_bc_for_nr)
        eq_scale = build_equation_scales(
            cells,
            ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
            ref_solid_kg_s=cfg.fuel_feed,
            ref_energy_W=cfg.fuel_feed * 20e6,
        )
        return build_jacobian_fd(
            x,
            f0,
            cells,
            reactor._apply_all_bc_for_nr,
            eq_scale,
            strategy=strategy,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
            affected_residual_cells_fn=reactor._affected_nr_residual_cells,
        )[0]

    J_fd = _assemble("block_tridiag_structured")
    J_hybrid = _assemble("hybrid_block_tridiag_structured")
    diff = (J_fd - J_hybrid).toarray()
    assert float(np.max(np.abs(diff))) < 1e-8


def test_jacobian_transport_k_stays_frozen_during_structured_column_perturbation() -> None:
    from tests.global_nr_fixtures import run_holdup_generate_initial_x0, setup_holdup_transport_reactor
    from src.solvers.global_nr_solver import pack_reactor, unpack_reactor
    from src.solvers.jacobian_transport_freeze import (
        enter_jacobian_transport_k_freeze,
        exit_jacobian_transport_k_freeze,
        restore_jacobian_frozen_transport_k,
    )

    reactor, cfg, t_profile = setup_holdup_transport_reactor(n_cells=3)
    run_holdup_generate_initial_x0(reactor, cfg, t_profile)
    reactor._apply_all_bc_for_nr()
    cells = reactor.cells
    x = pack_reactor(cells)
    k_ref = float(cells[1].K_solid_auf[0, S_CHAR])

    enter_jacobian_transport_k_freeze(cells)
    try:
        unpack_reactor(x, cells)
        reactor._apply_all_bc_for_nr()
        assert float(cells[1].K_solid_auf[0, S_CHAR]) == pytest.approx(k_ref)

        x_p = x.copy()
        x_p[0] += 1e-5
        unpack_reactor(x_p, cells)
        reactor._apply_local_bc_for_nr(0)
        restore_jacobian_frozen_transport_k(cells)
        assert float(cells[1].K_solid_auf[0, S_CHAR]) == pytest.approx(k_ref)
    finally:
        exit_jacobian_transport_k_freeze(cells)


def test_side_block_temperature_analytic_diagonal_matches_finite_difference() -> None:
    from src.solvers.nr_indexing import has_temperature_var

    cell = Cell()
    cell.cell_type = "return_leg"
    cell.solid_state_model = "holdup_transport"
    cell.m_solid.fill(0.5)
    cell.m_solid_ab_in.fill(0.2)
    cell.T = 1100.0
    cell.T_solid_ab_in = 1150.0
    cell.T_in_solid = 900.0
    assert has_temperature_var(cell)

    t_off = temperature_offset(cell)
    raw = analytic_packed_derivative(cell, row_local=t_off, col_local=t_off)
    assert raw is not None
    assert raw < 0.0

    h = 1e-5
    f0 = cell._calc_side_temperature_closure_residual()
    cell.T += h
    f1 = cell._calc_side_temperature_closure_residual()
    fd = (f1 - f0) / h
    assert raw == pytest.approx(fd, rel=1e-4, abs=1e-2)


def test_calc_side_temperature_closure_cp_dot_is_positive_floor() -> None:
    cell = Cell()
    cell.cell_type = "return_leg"
    cell.m_solid_ab_in.fill(1.0)
    cell.T_solid_ab_in = 1000.0
    assert calc_side_temperature_closure_cp_dot(cell) >= 1.0e3


def test_side_temperature_reference_prefers_solid_ab_inflow() -> None:
    cell = Cell()
    cell.cell_type = "return_leg"
    cell.m_solid_ab_in.fill(1.0)
    cell.T_solid_ab_in = 1200.0
    cell.T_in_solid = 800.0
    assert side_temperature_reference_K(cell) == pytest.approx(1200.0)


def test_bed_gas_dense_inventory_diagonal_is_minus_one() -> None:
    cell = Cell()
    cell.cell_type = "bed"
    assert inventory_self_derivative_unscaled(cell, row_local=0, col_local=0) == -1.0
    assert inventory_self_derivative_unscaled(cell, row_local=1, col_local=0) is None


def test_holdup_solid_diagonal_uses_transport_coefficients() -> None:
    from src.solvers.jacobian_analytic import _holdup_solid_self_derivative
    from src.solvers.nr_indexing import solid_flat_indices_for_nr

    cell = Cell()
    cell.cell_type = "bed"
    cell.solid_state_model = "holdup_transport"
    cell.K_solid_auf[0, S_CHAR] = 0.2
    cell.K_solid_ab[0, S_CHAR] = 0.3
    flat = solid_flat_indices_for_nr(cell)
    raw = _holdup_solid_self_derivative(cell, int(flat[0]))
    assert raw == pytest.approx(-0.5)
    sol_off = n_gas_var(cell)
    assert inventory_self_derivative_unscaled(cell, row_local=sol_off, col_local=sol_off) == pytest.approx(-0.5)


def test_holdup_anchor_diagonal_is_minus_one_without_transport() -> None:
    from src.solvers.jacobian_analytic import _holdup_solid_self_derivative
    from src.solvers.nr_indexing import solid_flat_indices_for_nr

    cell = Cell()
    cell.cell_type = "bed"
    cell.solid_state_model = "holdup_transport"
    flat = solid_flat_indices_for_nr(cell)
    assert _holdup_solid_self_derivative(cell, int(flat[0])) == -1.0
    cell.m_solid_in[0, S_CHAR] = 0.05
    assert _holdup_solid_self_derivative(cell, int(flat[0])) == -1.0


def test_holdup_k0_anchor_none_diagonal_is_minus_K() -> None:
    from src.solvers.jacobian_analytic import _holdup_solid_self_derivative
    from src.solvers.nr_indexing import solid_flat_indices_for_nr

    cell = Cell()
    cell.cell_type = "bed"
    cell.solid_state_model = "holdup_transport"
    cell.holdup_k0_anchor = "none"
    flat = solid_flat_indices_for_nr(cell)
    assert _holdup_solid_self_derivative(cell, int(flat[0])) == pytest.approx(0.0)
    cell.m_solid_in[0, S_CHAR] = 0.05
    assert _holdup_solid_self_derivative(cell, int(flat[0])) == pytest.approx(0.0)
    cell.K_solid_auf[0, S_CHAR] = 0.2
    cell.K_solid_ab[0, S_CHAR] = 0.3
    assert _holdup_solid_self_derivative(cell, int(flat[0])) == pytest.approx(-0.5)


def test_hybrid_overlay_matches_full_fd_on_inventory_column() -> None:
    from tests.global_nr_fixtures import build_initialized_lu_reactor
    from src.solvers.global_nr_solver import (
        build_equation_scales,
        build_jacobian_fd,
        global_residual,
        pack_reactor,
    )

    def _assemble(strategy: str):
        reactor = build_initialized_lu_reactor(n_cells=3)
        cells = reactor.cells
        cfg = reactor.config
        x = pack_reactor(cells)
        f0 = global_residual(x, cells, reactor._apply_all_bc_for_nr)
        eq_scale = build_equation_scales(
            cells,
            ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
            ref_solid_kg_s=cfg.fuel_feed,
            ref_energy_W=cfg.fuel_feed * 20e6,
        )
        return build_jacobian_fd(
            x,
            f0,
            cells,
            reactor._apply_all_bc_for_nr,
            eq_scale,
            strategy=strategy,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
            affected_residual_cells_fn=reactor._affected_nr_residual_cells,
        )[0]

    # 两次组装不能共用同一 reactor：FD 扰动会改写 K_solid 等非 NR 状态，unpack 不会恢复。
    J_fd = _assemble("block_tridiag_structured")
    J_hybrid = _assemble("hybrid_block_tridiag_structured")
    diff = (J_fd - J_hybrid).toarray()
    assert float(np.max(np.abs(diff))) < 1e-8


def test_bed_cell_temperature_diagonal_uses_outflow_cp_dot() -> None:
    cell = Cell()
    cell.cell_type = "bed"
    cell.T = 1100.0
    cell.N_d.fill(0.2)
    cell.N_b.fill(0.1)
    cell.m_solid.fill(0.05)
    t_off = temperature_offset(cell)
    raw = analytic_packed_derivative(cell, row_local=t_off, col_local=t_off)
    assert raw is not None
    assert raw < 0.0

    h = 1e-5
    t_save = float(cell.T)
    f0 = cell.calc_energy_balance()
    cell.T = t_save + h
    f1 = cell.calc_energy_balance()
    cell.T = t_save
    fd = (f1 - f0) / h
    assert raw == pytest.approx(fd, rel=1e-3, abs=1.0)
