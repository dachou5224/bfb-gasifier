"""阶段 F3：离开炭与内层气源汇共用 ``rate_R1``。"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from src.core.species import MW_C_KG_PER_MOL
from src.solvers.wirsum_1997_reaction_coord import (
    K_FB_DEFAULT,
    K_FB_MODEL,
    char_r1_extent,
    default_fb_inlets,
    gas_source_from_extent,
    leaving_char_from_outer,
    r1_fb_residual,
    solve_r1_fb,
    build_stage_f3_payload,
)
from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry
from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry


def test_extent_is_the_only_rate_entry() -> None:
    import src.solvers.wirsum_1997_reaction_coord as mod

    src = inspect.getsource(mod)
    assert "from src.kinetics.char_reactions import rate_R1" in src
    assert "A0 * np.exp" not in src
    assert "np.exp(-" not in src
    assert "from src.core.freeboard" not in src


def test_same_inputs_give_zero_r_gap() -> None:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    N = default_fb_inlets()
    M = 0.2
    a = char_r1_extent(T=geo_ex.T, P=geo_ex.P, N=N, M_kg=M, d_p=geo.d_p, rho_s=geo.rho_p)
    b = char_r1_extent(T=geo_ex.T, P=geo_ex.P, N=N, M_kg=M, d_p=geo.d_p, rho_s=geo.rho_p)
    assert a.e1 == pytest.approx(b.e1, rel=0.0, abs=0.0)
    src = gas_source_from_extent(a)
    assert src[0] == pytest.approx(-a.alpha * a.e1)


def test_residual_uses_shared_extent() -> None:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    N_in = default_fb_inlets()
    m_leave = 0.15
    M = m_leave / K_FB_DEFAULT
    x = np.concatenate([N_in, np.array([M])])
    F = r1_fb_residual(
        x, T=geo_ex.T, P=geo_ex.P, N_in=N_in, m_leave_char=m_leave,
        K_fb=K_FB_DEFAULT, d_p=geo.d_p, rho_s=geo.rho_p,
    )
    ext = char_r1_extent(
        T=geo_ex.T, P=geo_ex.P, N=N_in, M_kg=M, d_p=geo.d_p, rho_s=geo.rho_p
    )
    src = gas_source_from_extent(ext)
    np.testing.assert_allclose(F[:4], src)
    assert F[4] == pytest.approx(-MW_C_KG_PER_MOL * ext.e1, rel=1e-9)


def test_outer_leave_is_positive_for_fines() -> None:
    assert leaving_char_from_outer(recycle=False) > 0.0


def test_solve_coordinates_r_and_carbon() -> None:
    geo = BedFreeboardGeometry(d_p=2.0e-4)
    geo_ex = ExchangeGeometry()
    m_leave = leaving_char_from_outer(recycle=False)
    out = solve_r1_fb(
        T=geo_ex.T,
        P=geo_ex.P,
        N_in=default_fb_inlets(),
        m_leave_char=m_leave,
        K_fb=K_FB_DEFAULT,
        d_p=geo.d_p,
        rho_s=geo.rho_p,
    )
    assert out["r_gap_mol_s"] == 0.0
    assert out["K_fb_model"] == K_FB_MODEL
    if out["physical_convergence"]:
        assert abs(out["carbon_gas_minus_e1"]) < 1.0e-8
        assert out["residual_l2"] < 1.0e-6
    else:
        assert out["failure_reason"] in {
            "small_step",
            "line_search",
            "singular_jacobian",
            "max_iter",
        }
        assert out["physical_convergence"] is False


def test_payload_names_r1_only_and_skips_experiment() -> None:
    payload = build_stage_f3_payload()
    assert payload["solved"] is False
    assert payload["reactions"] == ["R1"]
    assert "R12" in payload["not"]
    dumped = str(payload)
    assert "y_CO_target" not in dumped
    assert "experiment" not in dumped.lower()
    assert payload["root"]["r_gap_mol_s"] == 0.0
    assert payload["inlet_minus_root_e1"] > 1.0
    assert payload["newton"]["physical_convergence"] is True
