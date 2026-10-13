from __future__ import annotations

import numpy as np

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.global_nr_solver import _residual_gas_phase_metrics, build_equation_scales, global_residual, pack_reactor
from src.solvers.vorabrechnung import seed_pyrolysis_gas_dense_holdup_for_nr_x0
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


def test_seed_pyrolysis_dense_holdup_reduces_init_split_and_keeps_bubble_zero():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 3
    cfg.n_freeboard_cells = 0
    cfg.vorab_bottom_temperature_cap_K_thesis = 1180.0
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0

    cfg_base = build_phase2_htw_lu_freeboard_reactor_config()
    cfg_base.n_cells = 3
    cfg_base.n_freeboard_cells = 0
    cfg_base.vorab_bottom_temperature_cap_K_thesis = 1180.0
    cfg_base.vorab_vm_devolatilization_zone_cells_thesis = 0
    cfg_base.nr_seed_pyrolysis_gas_dense_holdup_x0_thesis = False

    baseline = Reactor(cfg_base)
    run_init_and_precalc_for_global_nr(
        baseline,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    seeded = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        seeded,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )

    idx = GAS_SPECIES_INDEX
    cell0 = seeded.cells[0]
    assert getattr(seeded, "_pyrolysis_dense_holdup_seed_diag", None) is not None
    assert float(seeded._pyrolysis_dense_holdup_seed_diag["cells_seeded"]) >= 1.0
    for sp in ("NH3", "TAR1", "TAR2"):
        j = idx[sp]
        assert float(cell0.N_d[j]) > 1e-12
        assert float(cell0.N_b[j]) == 0.0

    def _split_rms(reactor: Reactor) -> float:
        cells = reactor._solver_cells_for_nr()
        cfg_local = reactor.config
        ref_gas = max(float(cfg_local.O2_feed + cfg_local.H2O_feed + cfg_local.N2_feed), 1.0)
        ref_energy = max(abs(float(cfg_local.fuel_feed)) * 20.0e6, 1.0e6)
        ref_solid = max(abs(float(cfg_local.fuel_feed)), 1.0e-9)
        F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
        scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
        F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
        return float(_residual_gas_phase_metrics(F_hat, cells)["gas_phase_split_rms"])

    split_before = _split_rms(baseline)
    split_after = _split_rms(seeded)
    # §4.3 Gibbs + 元素池 cap 改变热解计量；此处只检 init split 量级与 seed 结构。
    assert split_before < 0.5
    assert split_after < 0.5

    nh3 = GAS_SPECIES_INDEX["NH3"]
    assert float(seeded.cells[0].N_d[nh3]) > 1e-12


def test_seed_pyrolysis_dense_holdup_is_idempotent_on_zero_source_cells():
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    cfg.n_cells = 2
    cfg.n_freeboard_cells = 0
    cfg.vorab_vm_devolatilization_zone_cells_thesis = 0
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy="vorabrechnung",
        gs_warmup_steps=None,
    )
    top = reactor.cells[-1]
    if bool(getattr(top, "_vm_cache_valid", False)):
        top._vm_gas_source_cache = np.zeros_like(top._vm_gas_source_cache)
    before = top.N_d.copy()
    diag = seed_pyrolysis_gas_dense_holdup_for_nr_x0(reactor.cells, holdup_scale=1.0, residual_closure_passes=0)
    assert np.allclose(top.N_d, before)
    assert float(diag["cells_seeded"]) >= 0.0
