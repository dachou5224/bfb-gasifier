"""F_E1/F_E2 捆包之外：按反应逐个接入。不改生产默认，不改具名阶段编号。"""

from __future__ import annotations

import numpy as np

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import (
    NAMED_PROBE_RXNS,
    coupled_residual_named,
    default_enabled_rxns,
    default_x0_named,
    freeze_stage_q,
    homog_sources,
    probe_reaction_ladder,
    r234_dense_source,
    stage_for_enabled,
)
from src.solvers.wirsum_1997_stiffness_full_gas import HOMOG_RXN_NAMES, IDX


def test_named_stage_bundles_unchanged() -> None:
    assert default_enabled_rxns(0) == frozenset()
    assert default_enabled_rxns(1) == frozenset({"R2", "R3", "R4"})
    assert default_enabled_rxns(2) == frozenset({"R2", "R3", "R4", *HOMOG_RXN_NAMES})
    nq1 = freeze_stage_q(1)
    nq2 = freeze_stage_q(2)
    assert nq1.enabled_rxns == default_enabled_rxns(1)
    assert nq2.enabled_rxns == default_enabled_rxns(2)
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True


def test_isolated_r2_matches_masked_source_not_full_e1() -> None:
    nq_r2 = freeze_stage_q(1, enabled_rxns=frozenset({"R2"}))
    nq_e1 = freeze_stage_q(1)
    x = default_x0_named(nq_e1)
    x[IDX["H2O"]] = 0.20
    x[IDX["CO2"]] = 0.15
    F_r2 = coupled_residual_named(x, nq_r2)
    F_e1 = coupled_residual_named(x, nq_e1)
    assert not np.allclose(F_r2, F_e1)
    src_r2, _ = r234_dense_source(x[:-1], nq_r2.energy.species, rate_scale=1.0, enabled=frozenset({"R2"}))
    src_all, _ = r234_dense_source(x[:-1], nq_e1.energy.species, rate_scale=1.0)
    assert abs(float(src_r2[IDX["H2O"]])) > 0.0
    assert abs(float(src_all[IDX["CO2"]])) >= abs(float(src_r2[IDX["CO2"]])) - 1.0e-15


def test_isolated_r12_does_not_turn_on_r5() -> None:
    nq = freeze_stage_q(2, enabled_rxns=frozenset({"R12"}))
    x = default_x0_named(nq)
    sd12, sb12 = homog_sources(
        x[:39], nq.energy.species, rate_scale=1.0, rho_ash=75.0, enabled=frozenset({"R12"}),
    )
    sd5, _ = homog_sources(
        x[:39], nq.energy.species, rate_scale=1.0, rho_ash=75.0, enabled=frozenset({"R5"}),
    )
    assert stage_for_enabled(frozenset({"R12"})) == 2
    F12 = coupled_residual_named(x, nq)
    F5 = coupled_residual_named(x, freeze_stage_q(2, enabled_rxns=frozenset({"R5"})))
    assert not np.allclose(F12, F5)
    assert float(np.linalg.norm(sd12)) > 0.0 or float(np.linalg.norm(sb12)) > 0.0
    assert not np.allclose(sd12, sd5)


def test_probe_ladder_contract_no_furnace_claim() -> None:
    out = probe_reaction_ladder(mode="isolated", max_iter=1)
    assert out["not_a_furnace_convergence_claim"] is True
    assert out["physical_convergence"] is False
    assert out["baseline"]["label"] == "F_E0_R1_only"
    names = [row["label"] for row in out["isolated"]]
    assert names == [f"only_{n}" for n in NAMED_PROBE_RXNS]
    assert len(out["isolated"]) == 8
