"""P0 seam：NR 残差映射完整性（refresh / fastox / 隐藏状态回滚）。

验证求解器契约：接受步水力刷新后的 F 必须成为下一步 RHS；
线搜索 fastox 投影不得留下无法随 packed x 回滚的隐藏状态。
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp

from src.core.cell import Cell
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.global_nr_iteration import (
    _bed_species_holdup_mol,
    _bed_syngas_holdup_mol,
    _floor_syngas_holdups_in_state,
    _restore_excluded_gas_holdups,
    _snapshot_excluded_gas_holdups,
)
from src.solvers.global_nr_solver import (
    _clip_dx,
    diagnose_newton_step_consistency,
    pack_reactor,
    solve_global_nr,
    unpack_reactor,
)
from src.solvers.nr_indexing import install_main_nr_gas_indexing_on_cells
from src.solvers.vorabrechnung import (
    hydrodynamics_freeze_fingerprint,
    live_hydrodynamics_fingerprint,
    prepare_inner_nr_hydrodynamics_freeze,
    refresh_cells_hydrodynamics_for_frozen_inner,
)
from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def test_phase2_default_uses_outer_fixed_hydro_semantics():
    """P1: production Phase2 isomorphic path keeps h fixed during Inner."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    assert cfg.nr_refresh_hydrodynamics_on_accepted_step_thesis is False


def test_bed_syngas_holdup_sums_co_h2_ch4_only():
    cell = Cell()
    idx = GAS_SPECIES_INDEX
    cell.cell_type = "bed"
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["CO"]] = 1.0
    cell.N_b[idx["H2"]] = 2.0
    cell.N_d[idx["CH4"]] = 3.0
    cell.N_b[idx["N2"]] = 100.0
    assert _bed_syngas_holdup_mol([cell]) == pytest.approx(6.0)


def test_floor_syngas_holdups_keeps_co_fraction():
    """LS packed-state floor restores CO/H₂/CH₄ before fastox/merit."""
    cell = Cell()
    idx = GAS_SPECIES_INDEX
    cell.cell_type = "bed"
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["CO"]] = 10.0
    cell.N_b[idx["H2"]] = 4.0
    install_main_nr_gas_indexing_on_cells([cell], None)
    x_base = pack_reactor([cell])
    x_trial = x_base.copy()
    x_trial[:] = 0.0  # wipe everything
    n_hit = _floor_syngas_holdups_in_state(x_trial, x_base, [cell], ratio=0.05)
    assert n_hit >= 1
    unpack_reactor(x_trial, [cell])
    assert float(cell.N_d[idx["CO"]]) == pytest.approx(0.5)
    assert float(cell.N_b[idx["H2"]]) == pytest.approx(0.2)


def test_clip_dx_global_alpha_keeps_significant_co():
    """Global holdup α keeps significant CO above 20% of binding inventory."""
    from src.solvers.nr_indexing import n_gas_var

    cell = Cell()
    idx = GAS_SPECIES_INDEX
    cell.cell_type = "bed"
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["CO"]] = 10.0
    cell.N_b[idx["CO"]] = 5.0
    install_main_nr_gas_indexing_on_cells([cell], None)
    x = pack_reactor([cell])
    dx = np.zeros_like(x)
    dx[: n_gas_var(cell)] = -1.0e3
    clipped = _clip_dx(dx, [cell], ref_gas_mol_s=100.0, ref_solid_kg_s=1.0)
    unpack_reactor(np.maximum(x + clipped, 0.0), [cell])
    assert float(cell.N_b[idx["CO"]]) >= 1.0 - 1.0e-9
    assert float(cell.N_d[idx["CO"]]) >= 2.0 - 1.0e-9


def test_syngas_collapse_guard_rejects_zero_co_accept(monkeypatch):
    """LS must not accept a fastox-repaired trial that wiped bed CO (H₂ may remain)."""
    cell = Cell()
    idx = GAS_SPECIES_INDEX
    cell.cell_type = "bed"
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["CO"]] = 5.0
    cell.N_b[idx["H2"]] = 5.0
    x0 = pack_reactor([cell])
    n_total = len(x0)

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (
            sp.eye(n_total, format="csr"),
            {"residual_cell_calls": 0, "nnz": n_total, "zero_rows": 0, "zero_cols": 0},
        ),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )

    base_co = _bed_species_holdup_mol([cell], "CO")
    base_h2 = _bed_species_holdup_mol([cell], "H2")

    def _fake_residual(x, cells, apply_bc_fn):
        unpack_reactor(np.asarray(x, dtype=np.float64), cells)
        # Pathology: CO wiped, H2 kept — total syngas may still look OK.
        co = _bed_species_holdup_mol(cells, "CO")
        if abs(co - base_co) > 1.0e-9:
            for c in cells:
                c.N_d[idx["CO"]] = 0.0
                c.N_b[idx["CO"]] = 0.0
            return np.full(n_total, 0.01)
        return np.full(n_total, 1.0)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _fake_residual)
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._solve_linear_step",
        lambda J, rhs, **kwargs: (
            np.full(n_total, -10.0),
            "test",
            False,
            {"band_lu_s": 0.0, "side_update_s": 0.0, "schur_size": 0, "fallback_used": False},
        ),
    )

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=2,
        tol_rms=0.1,
        lambda_init=1.0,
        prefer_full_step=True,
        line_search_syngas_collapse_guard=True,
        line_search_syngas_collapse_ratio=0.05,
        n_damp_halvings=4,
    )
    assert int(result["counts"].get("line_search_syngas_collapse_rejects", 0)) >= 1
    assert _bed_species_holdup_mol([cell], "CO") == pytest.approx(base_co)
    assert _bed_species_holdup_mol([cell], "H2") == pytest.approx(base_h2)


def test_outer_fixed_semantics_skips_accept_step_hydro_refresh(monkeypatch):
    """Hamel outer_fixed: accepted Newton steps must not refresh hydro mid-Inner."""
    cell = Cell()
    x0 = pack_reactor([cell])
    n_total = len(x0)
    refresh_calls: list[int] = []

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (
            sp.eye(n_total, format="csr"),
            {"residual_cell_calls": 0, "nnz": n_total},
        ),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )
    monkeypatch.setattr(
        "src.solvers.vorabrechnung.refresh_cells_hydrodynamics_for_frozen_inner",
        lambda cells, *, cfg=None: refresh_calls.append(len(cells)),
    )

    calls = {"n": 0}

    def _fake_global_residual(x, cells, apply_bc_fn):
        calls["n"] += 1
        if calls["n"] == 1:
            return np.full(n_total, 2.0)
        return np.full(n_total, 0.05)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _fake_global_residual)

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=2,
        tol_rms=0.1,
        lambda_init=1.0,
        prefer_full_step=True,
        refresh_hydrodynamics_on_accepted_step=False,
    )

    assert refresh_calls == []
    assert result["counts"].get("hydrodynamics_refresh_on_accept", 0) in (0, None)
    assert result["counts"].get("hydrodynamics_refresh_residual_recompute", 0) in (0, None)
    assert result["refresh_hydrodynamics_on_accepted_step"] is False
    assert float(result["rms_scaled_final"]) == pytest.approx(0.05)


def test_accepted_step_hydro_refresh_recomputes_residual_for_next_rhs(monkeypatch):
    """Opt-in accept_refresh: Accept → refresh hydro → F recomputed for next RHS/J."""
    cell = Cell()
    x0 = pack_reactor([cell])
    n_total = len(x0)
    refresh_calls: list[int] = []
    residual_after_refresh: list[np.ndarray] = []

    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_equation_scales",
        lambda *args, **kwargs: np.ones(n_total),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver.build_jacobian_fd",
        lambda *args, **kwargs: (
            sp.eye(n_total, format="csr"),
            {"residual_cell_calls": 0, "nnz": n_total},
        ),
    )
    monkeypatch.setattr(
        "src.solvers.global_nr_solver._clip_dx",
        lambda dx, cells, ref_gas_mol_s, ref_solid_kg_s, **kwargs: dx,
    )
    monkeypatch.setattr(
        "src.solvers.vorabrechnung.refresh_cells_hydrodynamics_for_frozen_inner",
        lambda cells, *, cfg=None: refresh_calls.append(len(cells)),
    )

    calls = {"n": 0}

    def _fake_global_residual(x, cells, apply_bc_fn):
        calls["n"] += 1
        if calls["n"] == 1:
            return np.full(n_total, 2.0)
        if refresh_calls:
            F = np.full(n_total, 0.05)
            residual_after_refresh.append(F.copy())
            return F
        return np.full(n_total, 0.5)

    monkeypatch.setattr("src.solvers.global_nr_solver.global_residual", _fake_global_residual)

    result = solve_global_nr(
        cells=[cell],
        apply_bc_fn=lambda: None,
        max_iter=2,
        tol_rms=0.1,
        lambda_init=1.0,
        prefer_full_step=True,
        refresh_hydrodynamics_on_accepted_step=True,
    )

    assert refresh_calls == [1]
    assert result["counts"].get("hydrodynamics_refresh_residual_recompute", 0) == 1
    assert residual_after_refresh, "post-refresh residual must be evaluated"
    assert float(result["rms_scaled_final"]) == pytest.approx(0.05)
    hist = result.get("clip_history") or []
    assert hist, "clip_history should record accepted step"
    assert float(hist[0]["merit_pre_refresh"]) == pytest.approx(0.5)
    assert float(hist[0]["merit_post_refresh"]) == pytest.approx(0.05)


@pytest.mark.parametrize("phase_update", ["full_hydro", "local_participating"])
def test_fastox_projection_is_idempotent_on_coexistence(phase_update: str):
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 0.5
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 1.0
    cell.N_d[idx["N2"]] = 10.0
    cell.N_b[idx["N2"]] = 3.0

    project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update=phase_update
    )
    n_d_1 = cell.N_d.copy()
    n_b_1 = cell.N_b.copy()
    project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update=phase_update
    )
    assert np.allclose(cell.N_d, n_d_1, atol=1e-12, rtol=0.0)
    assert np.allclose(cell.N_b, n_b_1, atol=1e-12, rtol=0.0)


def test_fastox_projection_is_identity_when_no_fuel_o2_coexistence():
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    # Feasible: fuel present, no O2 → projection must be a no-op.
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["CO"]] = 0.5
    cell.N_d[idx["N2"]] = 8.0
    cell.N_b[idx["N2"]] = 2.0
    n_d_0 = cell.N_d.copy()
    n_b_0 = cell.N_b.copy()

    diag = project_bed_holdup_fast_oxidation_closure(reactor.cells, cfg, bed_cell_limit=1)
    assert float(diag["cells_touched"]) == 0.0
    assert np.allclose(cell.N_d, n_d_0, atol=1e-15, rtol=0.0)
    assert np.allclose(cell.N_b, n_b_0, atol=1e-15, rtol=0.0)


def test_fastox_full_hydro_redistributes_inert_n2_when_oxidation_fires():
    """x0 full_hydro mode may redistribute the full gas vector by bubble share."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 1.0
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 2.0
    cell.N_d[idx["N2"]] = 9.0
    cell.N_b[idx["N2"]] = 1.0
    n2_dense_before = float(cell.N_d[idx["N2"]])
    n2_total = float(cell.N_d[idx["N2"]] + cell.N_b[idx["N2"]])

    project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update="full_hydro"
    )
    n2_dense_after = float(cell.N_d[idx["N2"]])
    n2_total_after = float(cell.N_d[idx["N2"]] + cell.N_b[idx["N2"]])
    assert n2_total_after == pytest.approx(n2_total, abs=1e-12)
    assert abs(n2_dense_after - n2_dense_before) > 1e-9


def test_fastox_local_participating_preserves_inert_n2_phase_split():
    """P2 LS projection must not touch non-participating N2 phase holdups."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 1.0
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 2.0
    cell.N_d[idx["N2"]] = 9.0
    cell.N_b[idx["N2"]] = 1.0
    n2_d0 = float(cell.N_d[idx["N2"]])
    n2_b0 = float(cell.N_b[idx["N2"]])

    diag = project_bed_holdup_fast_oxidation_closure(
        reactor.cells, cfg, bed_cell_limit=1, phase_update="local_participating"
    )
    assert float(diag["cells_touched"]) >= 1.0
    assert float(cell.N_d[idx["N2"]]) == pytest.approx(n2_d0, abs=1e-15)
    assert float(cell.N_b[idx["N2"]]) == pytest.approx(n2_b0, abs=1e-15)
    assert float(cell.N_d[idx["H2"]] + cell.N_b[idx["H2"]]) < 1e-12


def test_excluded_gas_snapshot_rolls_back_h2s_after_rejected_full_hydro_fastox():
    """full_hydro x0-style P can move H2S; rejected trial must restore it."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.nr_include_pyrolysis_gas_species_in_unknowns_thesis = True
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    install_main_nr_gas_indexing_on_cells(reactor.cells, cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX

    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 1.0
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 2.0
    cell.N_d[idx["H2S"]] = 0.8
    cell.N_b[idx["H2S"]] = 0.2
    cell.N_d[idx["N2"]] = 5.0
    cell.N_b[idx["N2"]] = 5.0

    x_base = pack_reactor([cell])
    snap = _snapshot_excluded_gas_holdups([cell])
    h2s_d0 = float(cell.N_d[idx["H2S"]])
    h2s_b0 = float(cell.N_b[idx["H2S"]])

    x_trial = x_base.copy()
    _restore_excluded_gas_holdups([cell], snap)
    unpack_reactor(x_trial, [cell])
    project_bed_holdup_fast_oxidation_closure(
        [cell],
        SimpleNamespace(
            vorab_transport_x0_fast_oxidation_closure_thesis=True,
            fuel_type="coal",
        ),
        bed_cell_limit=1,
        phase_update="full_hydro",
    )
    assert abs(float(cell.N_d[idx["H2S"]]) - h2s_d0) > 1e-12 or abs(
        float(cell.N_b[idx["H2S"]]) - h2s_b0
    ) > 1e-12

    _restore_excluded_gas_holdups([cell], snap)
    unpack_reactor(x_base, [cell])
    assert float(cell.N_d[idx["H2S"]]) == pytest.approx(h2s_d0, abs=1e-15)
    assert float(cell.N_b[idx["H2S"]]) == pytest.approx(h2s_b0, abs=1e-15)


def test_local_fastox_leaves_h2s_untouched_without_rollback():
    """P2 local LS projection must not mutate NR-excluded H2S."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.nr_include_pyrolysis_gas_species_in_unknowns_thesis = True
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    install_main_nr_gas_indexing_on_cells(reactor.cells, cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 1.0
    cell.N_d[idx["O2"]] = 2.0
    cell.N_b[idx["O2"]] = 2.0
    cell.N_d[idx["H2S"]] = 0.8
    cell.N_b[idx["H2S"]] = 0.2
    h2s_d0 = float(cell.N_d[idx["H2S"]])
    h2s_b0 = float(cell.N_b[idx["H2S"]])

    project_bed_holdup_fast_oxidation_closure(
        [cell],
        SimpleNamespace(
            vorab_transport_x0_fast_oxidation_closure_thesis=True,
            fuel_type="coal",
        ),
        bed_cell_limit=1,
        phase_update="local_participating",
    )
    assert float(cell.N_d[idx["H2S"]]) == pytest.approx(h2s_d0, abs=1e-15)
    assert float(cell.N_b[idx["H2S"]]) == pytest.approx(h2s_b0, abs=1e-15)


def test_residual_base_trial_base_reproducible_with_fastox_rollback():
    """base → fastox trial → restore+base must reproduce the base residual vector."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 0
    cfg.nr_include_pyrolysis_gas_species_in_unknowns_thesis = True
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cells = list(reactor.cells)
    install_main_nr_gas_indexing_on_cells(cells, cfg)
    idx = GAS_SPECIES_INDEX

    for cell in cells:
        cell.calc_hydrodynamics()
        cell._snapshot_vorabrechnung_hydrodynamics()
        cell.enable_inner_nr_vorabrechnung_freeze()
        cell.N_d[idx["H2S"]] = 0.3
        cell.N_b[idx["H2S"]] = 0.1

    from src.solvers.global_nr_solver import global_residual

    x_base = pack_reactor(cells)
    F_base_1 = global_residual(x_base, cells, lambda: None).copy()
    snap = _snapshot_excluded_gas_holdups(cells)

    # Trial: perturb packed state and apply fastox.
    x_trial = x_base.copy()
    x_trial[: min(4, x_trial.size)] *= 1.05
    _restore_excluded_gas_holdups(cells, snap)
    unpack_reactor(x_trial, cells)
    project_bed_holdup_fast_oxidation_closure(
        cells,
        SimpleNamespace(
            vorab_transport_x0_fast_oxidation_closure_thesis=True,
            fuel_type="coal",
        ),
        phase_update="local_participating",
    )
    _ = global_residual(pack_reactor(cells), cells, lambda: None)

    # Rollback to base mapping.
    _restore_excluded_gas_holdups(cells, snap)
    F_base_2 = global_residual(x_base, cells, lambda: None)
    assert np.allclose(F_base_2, F_base_1, atol=1e-10, rtol=1e-10)


def test_prepare_inner_freeze_realigns_after_preinner_like_state_change():
    """Pre-Inner-like N/T change leaves stale freeze; Inner prepare must resnapshot."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 0
    reactor = Reactor(cfg)
    cells = list(reactor.cells)
    refresh_cells_hydrodynamics_for_frozen_inner(cells, cfg=cfg)
    freeze_before = hydrodynamics_freeze_fingerprint(cells)

    # Simulate Pre-Inner mutation after outer refresh snapshot.
    for cell in cells:
        cell.T = float(cell.T) + 120.0
        cell.N_d *= 1.15
        cell.N_b *= 0.90
        cell._hydro_cache_valid = False
    for cell in cells:
        cell.calc_hydrodynamics()
    live_after = live_hydrodynamics_fingerprint(cells)
    freeze_stale = hydrodynamics_freeze_fingerprint(cells)
    assert abs(float(freeze_stale["sum"]) - float(live_after["sum"])) > 1e-6
    assert abs(float(freeze_stale["sum"]) - float(freeze_before["sum"])) < 1e-12

    prep = prepare_inner_nr_hydrodynamics_freeze(cells, cfg=cfg)
    freeze_after = hydrodynamics_freeze_fingerprint(cells)
    live_aligned = live_hydrodynamics_fingerprint(cells)
    assert float(prep["freeze_matches_live"]) == 1.0
    assert abs(float(freeze_after["sum"]) - float(live_aligned["sum"])) <= 1e-9
    assert abs(float(freeze_after["sum"]) - float(freeze_stale["sum"])) > 1e-6


def test_diagnose_newton_step_consistency_identity_and_clip_defect():
    n = 4
    J = sp.eye(n, format="csr")
    F_hat = np.array([1.0, -2.0, 0.5, 0.25], dtype=np.float64)
    dx = -F_hat.copy()
    good = diagnose_newton_step_consistency(J, F_hat, dx, dx)
    assert float(good["linear_defect_rel"]) == pytest.approx(0.0, abs=1e-14)
    assert float(good["clip_linear_defect_rel"]) == pytest.approx(0.0, abs=1e-14)

    dx_clipped = dx.copy()
    dx_clipped[0] = 0.0  # hard clip destroys descent for first equation
    bad = diagnose_newton_step_consistency(J, F_hat, dx, dx_clipped)
    assert float(bad["linear_defect_rel"]) == pytest.approx(0.0, abs=1e-14)
    assert float(bad["clip_linear_defect_rel"]) > 0.1
    assert float(bad["clipped_components"]) == 1.0


def test_local_fastox_displacement_is_o_lambda_from_feasible_base():
    """P2: from feasible x, ||P(x+λd)-x|| scales as O(λ) (no inert redistribution)."""
    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_transport_x0_fast_oxidation_closure_thesis = True
    reactor = Reactor(cfg)
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_d[:] = 0.0
    cell.N_b[:] = 0.0
    # Feasible: fuel only, no O2.
    cell.N_d[idx["H2"]] = 1.0
    cell.N_b[idx["H2"]] = 0.5
    cell.N_d[idx["N2"]] = 10.0
    cell.N_b[idx["N2"]] = 2.0
    x0 = pack_reactor([cell])
    # Direction that introduces O2 into the trial (creates fuel–O2 overlap ~ λ).
    dx = np.zeros_like(x0)
    # Packed layout: dense major gases then bubble; O2 is among majors.
    # Use a dense O2 slot if present; otherwise first DOF as fallback.
    from src.solvers.nr_indexing import install_main_nr_gas_indexing_on_cells, main_nr_gas_idx_for_cell

    install_main_nr_gas_indexing_on_cells([cell], cfg)
    x0 = pack_reactor([cell])
    dx = np.zeros_like(x0)
    gas_idx = [int(i) for i in main_nr_gas_idx_for_cell(cell)]
    o2_pos = gas_idx.index(int(idx["O2"]))
    dx[o2_pos] = 1.0  # add O2 in dense block

    disps = []
    for lam in (1.0e-3, 1.0e-4, 1.0e-5):
        unpack_reactor(x0 + lam * dx, [cell])
        n2_d = float(cell.N_d[idx["N2"]])
        n2_b = float(cell.N_b[idx["N2"]])
        project_bed_holdup_fast_oxidation_closure(
            [cell],
            SimpleNamespace(
                vorab_transport_x0_fast_oxidation_closure_thesis=True,
                fuel_type="coal",
            ),
            bed_cell_limit=1,
            phase_update="local_participating",
        )
        assert float(cell.N_d[idx["N2"]]) == pytest.approx(n2_d, abs=1e-15)
        assert float(cell.N_b[idx["N2"]]) == pytest.approx(n2_b, abs=1e-15)
        disps.append(float(np.linalg.norm(pack_reactor([cell]) - x0)))

    # Displacement must shrink roughly with λ (allow generous factor for stoich scaling).
    assert disps[0] > disps[1] > 0.0
    assert disps[1] > disps[2] > 0.0
    assert disps[0] / disps[2] == pytest.approx(100.0, rel=0.5)
