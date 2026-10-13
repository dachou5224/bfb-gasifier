"""F_E1：在 F_E0 上接入 R2/R3/R4。不改生产默认。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.core.species import MW_C_KG_PER_MOL
from src.solvers.wirsum_1997_coupled_f_energy import coupled_residual as residual_e0
from src.solvers.wirsum_1997_coupled_f_energy import freeze_energy_q, pack, default_x0
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    freeze_stage_q,
    r234_dense_source,
    scan_E_named,
)
from src.solvers.wirsum_1997_named_f_registry import f_name, n_unknowns, omitted_named
from src.solvers.wirsum_1997_stiffness_full_gas import IDX


def _zero_reactants(x: np.ndarray) -> np.ndarray:
    y = x.copy()
    for base in (0, 7):
        y[base + IDX["CO2"]] = 0.0
        y[base + IDX["H2"]] = 0.0
        y[base + IDX["H2O"]] = 0.0
        y[base + IDX["O2"]] = 0.0
    return y


def test_e1_name_and_dimension() -> None:
    nq = freeze_stage_q(1)
    assert f_name(1) == "F_E1"
    assert n_unknowns(1) == 40
    assert "R2_R12" not in omitted_named(1)
    assert omitted_named(1)["R5_R8_R12"] == "omitted, named"
    assert omitted_named(1)["energy_wall_heat"] == "omitted, named"
    x = default_x0(nq.energy)
    F = coupled_residual_named(x, nq)
    assert F.size == 40


def test_e1_degenerates_to_e0_without_reactants() -> None:
    nq0 = freeze_stage_q(0)
    nq1 = freeze_stage_q(1)
    x = _zero_reactants(default_x0(nq0.energy))
    F0 = residual_e0(x, nq0.energy)
    F1 = coupled_residual_named(x, nq1)
    np.testing.assert_allclose(F0, F1, atol=1e-12)
    sp = x[:-1]
    src, cons = r234_dense_source(sp, nq1.energy.species, rate_scale=1.0)
    np.testing.assert_allclose(src, 0.0, atol=1e-15)
    np.testing.assert_allclose(cons, 0.0, atol=1e-15)


def test_e1_r2_consumes_char_and_makes_co_h2() -> None:
    nq = freeze_stage_q(1)
    q = nq.energy.species
    x = default_x0(nq.energy)
    x = _zero_reactants(x)
    x[IDX["H2O"]] = 0.20  # [mol/s] 仅密相蒸汽
    src, cons = r234_dense_source(x[:-1], q, rate_scale=1.0)
    assert float(np.sum(cons)) > 0.0
    assert src[IDX["CO"]] > 0.0
    assert src[IDX["H2"]] > 0.0
    assert src[IDX["H2O"]] < 0.0
    e2 = float(np.sum(cons)) / MW_C_KG_PER_MOL
    assert src[IDX["CO"]] == pytest.approx(e2, rel=1e-4)
    assert src[IDX["H2"]] == pytest.approx(e2, rel=1e-4)
    assert src[IDX["H2O"]] == pytest.approx(-e2, rel=1e-4)


def test_e1_E_of_T_not_a_furnace_target() -> None:
    nq = freeze_stage_q(1)
    out = scan_E_named(nq, np.array([500.0, 1100.0]))
    assert out["not_a_furnace_convergence_claim"] is True
    assert out["f_name"] == "F_E1"
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
