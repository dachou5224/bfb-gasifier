"""P0.3：给定 packed x，global_residual 必须是纯函数。Q9 无投影解包。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import S_CHAR
from src.core.reactor import Reactor
from src.solvers.global_nr_solver import (
    global_residual,
    global_residual_strict_domain_no_projection,
    pack_reactor,
)
from src.solvers.nr_indexing import (
    CLOSED_DOMAIN_T_OPEN_LOWER_K,
    DomainViolation,
    UNPACK_MODE_PRODUCTION_CLIP,
    UNPACK_MODE_STRICT_DOMAIN_NO_PROJECTION,
    cell_packed_fields_equal,
    max_feasible_alpha_closed_domain,
    packed_cell_in_closed_domain,
    packed_cell_safe_to_evaluate_F,
    pack_cell,
    packed_reactor_in_closed_domain,
    snapshot_cell_packed_fields,
    unpack_cell,
    unpack_cell_strict_domain_no_projection,
    unpack_reactor,
    unpack_reactor_strict_domain_no_projection,
)
from src.solvers.wirsum_1997_residual_projection_audit import (
    DOWNSTREAM_CLIPS,
    strict_path_still_projects_negative_state,
    unpack_projection_removed_from_strict_path,
)
from tests.global_nr_fixtures import call_generate_initial_x0_phase1, phase1_axial_t_profile
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config


def _two_cell_reactor() -> Reactor:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.n_cells = 2
    reactor = Reactor(cfg)
    call_generate_initial_x0_phase1(reactor, cfg, phase1_axial_t_profile(cfg))
    return reactor


def _cell_physics_snapshot(cells) -> dict[str, np.ndarray]:
    return {
        "N_d": np.stack([np.array(c.N_d, dtype=np.float64, copy=True) for c in cells]),
        "N_b": np.stack([np.array(c.N_b, dtype=np.float64, copy=True) for c in cells]),
        "T": np.array([float(c.T) for c in cells], dtype=np.float64),
        "m_solid_char": np.stack(
            [np.array(c.m_solid[:, S_CHAR], dtype=np.float64, copy=True) for c in cells]
        ),
        "eps_b": np.array([float(c.eps_b) for c in cells], dtype=np.float64),
        "K_bd": np.array([float(c.K_bd) for c in cells], dtype=np.float64),
    }


def test_global_residual_does_not_mutate_packed_x() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x0 = pack_reactor(cells)
    x_in = x0.copy()
    _ = global_residual(x_in, cells, reactor._apply_all_bc_for_nr)
    np.testing.assert_array_equal(x_in, x0)


def test_global_residual_is_repeatable_on_same_x() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x0 = pack_reactor(cells)
    f1 = global_residual(x0.copy(), cells, reactor._apply_all_bc_for_nr)
    snap1 = _cell_physics_snapshot(cells)
    f2 = global_residual(x0.copy(), cells, reactor._apply_all_bc_for_nr)
    snap2 = _cell_physics_snapshot(cells)
    np.testing.assert_allclose(f1, f2, rtol=0.0, atol=0.0)
    for key in snap1:
        np.testing.assert_allclose(snap1[key], snap2[key], rtol=0.0, atol=0.0)


def test_legal_pack_unpack_identity_on_strict_path() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x0 = pack_reactor(cells)
    ok, why = packed_reactor_in_closed_domain(x0, cells)
    assert ok, why
    unpack_reactor_strict_domain_no_projection(x0.copy(), cells)
    np.testing.assert_allclose(pack_reactor(cells), x0, rtol=0.0, atol=0.0)


def test_zero_species_is_legal_and_roundtrips() -> None:
    reactor = _two_cell_reactor()
    cell = reactor.cells[0]
    x = pack_cell(cell)
    x[0] = 0.0  # [mol/s] 物理允许的零组分
    ok, why = packed_cell_in_closed_domain(x, cell)
    assert ok, why
    ok_f, why_f = packed_cell_safe_to_evaluate_F(x, cell)
    assert ok_f, why_f
    unpack_cell_strict_domain_no_projection(x, cell)
    np.testing.assert_allclose(pack_cell(cell), x, rtol=0.0, atol=0.0)


def test_illegal_negative_flow_does_not_rewrite_cell() -> None:
    reactor = _two_cell_reactor()
    cell = reactor.cells[0]
    before = snapshot_cell_packed_fields(cell)
    x = pack_cell(cell)
    x[0] = -1.0e-6  # [mol/s]
    with pytest.raises(DomainViolation):
        unpack_cell_strict_domain_no_projection(x, cell)
    assert cell_packed_fields_equal(snapshot_cell_packed_fields(cell), before)


def test_illegal_reactor_unpack_rewrites_no_cell() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    snaps = [snapshot_cell_packed_fields(c) for c in cells]
    x = pack_reactor(cells)
    x[0] = -1.0  # [mol/s]
    with pytest.raises(DomainViolation):
        unpack_reactor_strict_domain_no_projection(x, cells)
    for cell, snap in zip(cells, snaps):
        assert cell_packed_fields_equal(snapshot_cell_packed_fields(cell), snap)


def test_production_unpack_still_clips_negative_and_cold_T() -> None:
    """Q8：生产默认仍是 clip，不是 thesis 无投影。"""
    assert UNPACK_MODE_PRODUCTION_CLIP != UNPACK_MODE_STRICT_DOMAIN_NO_PROJECTION
    reactor = _two_cell_reactor()
    cell = reactor.cells[0]
    x = pack_cell(cell)
    x[0] = -0.5  # [mol/s]
    x[-1] = 200.0  # [K] 低于生产 fence
    unpack_cell(x, cell)
    packed = pack_cell(cell)
    assert packed[0] == 0.0
    assert packed[-1] == pytest.approx(300.0)
    # 同一非法向量在 strict 路径必须拒绝
    with pytest.raises(DomainViolation):
        unpack_cell_strict_domain_no_projection(x, cell)


def test_T_between_zero_and_production_fence_roundtrips_on_strict_path() -> None:
    reactor = _two_cell_reactor()
    cell = reactor.cells[0]
    x = pack_cell(cell)
    x[-1] = 200.0  # [K]
    assert x[-1] > CLOSED_DOMAIN_T_OPEN_LOWER_K
    unpack_cell_strict_domain_no_projection(x, cell)
    assert pack_cell(cell)[-1] == pytest.approx(200.0)


def test_strict_residual_Fx_Fy_Fx() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x = pack_reactor(cells)
    y = x.copy()
    y[0] = max(0.5 * float(y[0]), 0.0)  # [mol/s]
    y[-1] = float(y[-1]) + 5.0  # [K]
    apply_bc = reactor._apply_all_bc_for_nr
    fx = global_residual_strict_domain_no_projection(x.copy(), cells, apply_bc)
    fy = global_residual_strict_domain_no_projection(y.copy(), cells, apply_bc)
    fx2 = global_residual_strict_domain_no_projection(x.copy(), cells, apply_bc)
    np.testing.assert_allclose(fx, fx2, rtol=0.0, atol=0.0)
    assert not np.array_equal(fx, fy)


def test_production_residual_Fx_Fy_Fx() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x = pack_reactor(cells)
    y = x * 1.01
    y = np.maximum(y, 0.0)
    apply_bc = reactor._apply_all_bc_for_nr
    fx = global_residual(x.copy(), cells, apply_bc)
    _ = global_residual(y.copy(), cells, apply_bc)
    fx2 = global_residual(x.copy(), cells, apply_bc)
    np.testing.assert_allclose(fx, fx2, rtol=0.0, atol=0.0)


def test_fractional_boundary_step_is_named_and_allows_zero() -> None:
    reactor = _two_cell_reactor()
    cells = reactor.cells
    x = pack_reactor(cells)
    dx = np.zeros_like(x)
    dx[0] = -2.0 * float(x[0])  # 全步会变负
    alpha, policy = max_feasible_alpha_closed_domain(x, dx, cells)
    assert policy == "fractional_boundary_to_closed_domain"
    assert 0.0 < alpha <= 0.5 + 1.0e-12
    x_trial = x + alpha * dx
    unpack_reactor_strict_domain_no_projection(x_trial, cells)
    dx0 = np.zeros_like(x)
    dx0[0] = -1.0
    x_zero = x.copy()
    x_zero[0] = 0.0
    alpha0, reason = max_feasible_alpha_closed_domain(x_zero, dx0, cells)
    assert alpha0 == 0.0
    assert "active_set" in reason


def test_downstream_clip_audit_does_not_claim_rates_are_clean() -> None:
    assert unpack_projection_removed_from_strict_path()
    assert strict_path_still_projects_negative_state()
    kinds = {rec["kind"] for rec in DOWNSTREAM_CLIPS}
    assert "state_projection" in kinds
    assert "maps_negative_concentration_to_zero" in kinds
    unpack_recs = [rec for rec in DOWNSTREAM_CLIPS if rec["symbol"] == "unpack_cell"]
    assert len(unpack_recs) == 1
    assert unpack_recs[0]["on_strict_unpack_path"] is False
    assert global_residual.__kwdefaults__["unpack_fn"] is unpack_reactor
