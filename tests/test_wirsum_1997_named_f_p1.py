"""P1：Hegermann 燃料 waf 碳基准、DAEM 剩余物出口、焦油稀组分/EOS、独立账。

验收按 docs/audit_reports/named_f_e8_review_2026-09-16.md §5 P1：
多个 T、τ、进料与灰含量下，独立于内部残差 F 的 C/H/O/灰及总焓账；
覆盖明显未完全热解的状态。不改生产默认，不开放水力 P。
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from src.core.species import MW_C_KG_PER_MOL
from src.solvers.wirsum_1997_coupled_f_energy import enthalpy_buckets_W, freeze_energy_q
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import freeze_q
from src.solvers.wirsum_1997_named_f import (
    default_x0_named,
    freeze_stage_q,
    outer_P_named,
    r10_c_tar,
    tar_concentration_diagnostics,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_p1_ledger import mapping_ledger, named_stream_ledger
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    CARBON_BASIS_FUEL_WAF,
    REMAINDER_DESTINATION,
    TAR_CONC_DILUTE_SHARED_SLOT,
    dry_mass_from_cho_atoms,
    example_lu_like_feed,
    map_one_state,
)
from src.thermal.hegermann_yield import (
    T_WS_MAX_K,
    T_WS_MIN_K,
    hegermann_applicability,
    xc_koks_hegermann,
    xc_teer_hegermann,
)

_P_PA: float = 2.5e6  # [Pa]
_ATOM: float = 1.0e-9
_H_ID: float = 1.0e-6  # [W]


def _scale_feed(feed, *, scale: float = 1.0, ash_dry_wt: float | None = None):
    ash = float(feed.ash_dry_wt) if ash_dry_wt is None else float(ash_dry_wt)
    out = replace(
        feed,
        nC_fuel=float(feed.nC_fuel) * scale,
        nH_fuel=float(feed.nH_fuel) * scale,
        nO_fuel=float(feed.nO_fuel) * scale,
        nN_fuel=float(feed.nN_fuel) * scale,
        nS_fuel=float(feed.nS_fuel) * scale,
        nC_fixed=float(feed.nC_fixed) * scale,
        nC_vm=float(feed.nC_vm) * scale,
        m_moist_kg_s=float(feed.m_moist_kg_s) * scale,
        ash_dry_wt=ash,
    )
    return replace(out, m_dry_kg_s=dry_mass_from_cho_atoms(out))


def test_hegermann_coke_not_proximate_fixed_and_varies_with_T() -> None:
    feed = example_lu_like_feed()
    rows = []
    for T in (1073.0, 1173.0, 1273.0):
        led = mapping_ledger(feed, P_Pa=_P_PA, T_K=T, tau_s=1.0, X_end=1.0)
        assert led["carbon_basis"] == CARBON_BASIS_FUEL_WAF
        assert led["hegermann_coke_basis_gap"] < _ATOM
        assert led["x_C_koks_recovered"] == pytest.approx(
            xc_koks_hegermann(_P_PA, T), abs=1e-12
        )
        assert led["x_C_teer_recovered"] == pytest.approx(
            xc_teer_hegermann(_P_PA, T), abs=1e-9
        )
        assert abs(led["hegermann_koks_minus_fixed"]) > 1.0e-6
        rows.append(led["nC_koks"])
    assert rows[0] != pytest.approx(rows[1], abs=1e-12)
    assert rows[1] != pytest.approx(rows[2], abs=1e-12)
    q = freeze_q()
    assert q.carbon_basis == CARBON_BASIS_FUEL_WAF
    assert q.m_feed[1] == pytest.approx(q.nC_koks * MW_C_KG_PER_MOL)
    assert abs(q.nC_koks - q.nC_fixed) > 1.0e-9


def test_hegermann_range_declared_and_clipped() -> None:
    appl = hegermann_applicability(1.0e4, 800.0)
    assert appl["fuel"] == "Rheinische Braunkohle"
    assert appl["fuel_type"] == "coal"
    assert appl["T_clipped"] is True
    assert appl["P_clipped"] is True
    assert float(appl["T_used_K"]) == T_WS_MIN_K
    appl_ok = hegermann_applicability(_P_PA, 1200.0)
    assert appl_ok["T_clipped"] is False
    assert appl_ok["P_clipped"] is False
    assert T_WS_MIN_K < 1200.0 < T_WS_MAX_K
    with pytest.raises(ValueError, match="Rheinische"):
        map_one_state(
            example_lu_like_feed(),
            P_Pa=_P_PA,
            T_K=1200.0,
            X_vm=np.array([0.0, 1.0]),
            tau_s=np.array([1.0]),
            occupancy_source="named_p1",
            fuel_type="biomass",
        )


def test_mapping_elemental_ledger_multiple_T_tau_feed_ash() -> None:
    base = example_lu_like_feed()
    cases: list[tuple[float, float, float, float | None]] = []
    for T in (800.0, 1100.0, 1200.0, 1420.0):
        for tau in (0.1, 1.0):
            cases.append((T, tau, 1.0, None))
    for scale in (0.5, 2.0):
        cases.append((1100.0, 1.0, scale, None))
    for ash in (5.0, 11.41, 25.0):
        cases.append((1100.0, 1.0, 1.0, ash))
    cases.append((800.0, 0.1, 2.0, 25.0))
    incomplete = 0
    for T, tau, scale, ash in cases:
        feed = _scale_feed(base, scale=scale, ash_dry_wt=ash)
        led = mapping_ledger(feed, P_Pa=_P_PA, T_K=T, tau_s=tau)
        assert led["not_from_F"] is True
        assert led["atom_C_gap"] < _ATOM
        assert led["atom_H_gap"] < _ATOM
        assert led["atom_O_gap"] < _ATOM
        assert led["remainder_destination"] == REMAINDER_DESTINATION
        assert np.isfinite(led["H_map_gap_W"])
        H_parts = (
            led["H_prod_gas_tar_W"]
            + led["H_koks_W"]
            + led["H_unreleased_W"]
            + led["H_ash_W"]
        )
        assert np.isfinite(H_parts)
        if led["X_end"] < 0.99:
            incomplete += 1
            assert led["nC_unreleased"] > 1.0e-8
            assert led["H_unreleased_W"] != pytest.approx(0.0, abs=1e-12)
            assert led["m_unreleased_kg_s"] > 0.0
        if abs(led["X_end"] - 1.0) < 1.0e-8:
            assert led["nC_unreleased"] == pytest.approx(0.0, abs=_ATOM)
            assert led["x_C_koks_recovered"] == pytest.approx(
                led["x_C_koks_fuel_waf"], abs=1e-12
            )
            assert led["x_C_teer_recovered"] == pytest.approx(
                led["x_C_teer_fuel_waf"], abs=1e-9
            )
        ash_frac = float(feed.ash_dry_wt) / 100.0
        assert led["m_ash_kg_s"] == pytest.approx(
            float(feed.m_dry_kg_s) * ash_frac, rel=1e-12
        )
        if T < T_WS_MIN_K:
            assert led["hegermann"]["T_clipped"] is True
    assert incomplete >= 3


def test_daem_remainder_enters_named_q_and_energy() -> None:
    nq = freeze_stage_q(4)
    x = default_x0_named(nq)
    x[-1] = 800.0
    nq_p = outer_P_named(x, nq)
    assert nq_p.energy.species.nC_unreleased > 1.0e-4
    assert nq_p.energy.species.m_unreleased_kg_s > 0.0
    assert nq_p.energy.species.remainder_destination == REMAINDER_DESTINATION
    assert nq_p.energy.species.nC_koks == pytest.approx(
        xc_koks_hegermann(float(nq.energy.species.geo_ex.P), 800.0)
        * nq.energy.feed.nC_fuel,
        abs=1e-12,
    )
    streams = named_stream_ledger(x, nq_p)
    assert streams["nC_unreleased"] == pytest.approx(
        nq_p.energy.species.nC_unreleased, abs=0.0
    )
    assert streams["H_enthalpy_gap_W"] == pytest.approx(
        abs(streams["energy_row_W"]), abs=_H_ID
    )
    x_hot = x.copy()
    x_hot[-1] = 1420.0
    nq_hot = outer_P_named(x_hot, nq)
    assert nq_hot.energy.species.nC_unreleased < 1.0e-5
    assert nq_hot.energy.species.nC_koks != pytest.approx(
        nq_p.energy.species.nC_koks, abs=1e-12
    )


def test_tar_dilute_shared_slot_vs_strict_eos() -> None:
    nq = freeze_stage_q(6)
    x = default_x0_named(nq)
    u = unpack_named(x, 6)
    T = float(u["T"])
    P = float(nq.energy.species.geo_ex.P)
    n_tar = float(u["N_tar1"]) + float(u["N_tar2"])
    diag = tar_concentration_diagnostics(u["N_d"], u["N_b"], n_tar, T, P)
    assert diag["model_used_for_R10_R11"] == TAR_CONC_DILUTE_SHARED_SLOT
    assert diag["strict_eos_total_C_over_P_RT"] == pytest.approx(1.0)
    assert diag["dilute_total_C_over_P_RT"] == pytest.approx(
        1.0 + diag["y_tar_dilute"], abs=1e-15
    )
    assert diag["eos_relative_excess"] == pytest.approx(diag["y_tar_dilute"], abs=1e-15)
    assert diag["phase_allocation"] == "instantaneous_well_mixed_equal_C_dense_and_bubble"
    assert diag["not_freeboard_tar"] is True
    C_tar = r10_c_tar(u["N_d"], u["N_b"], n_tar, T, P)
    assert C_tar == pytest.approx(diag["C_tar_dilute_mol_m3"], abs=1e-15)
    streams = named_stream_ledger(x, nq)
    assert streams["tar_eos"]["model_used_for_R10_R11"] == TAR_CONC_DILUTE_SHARED_SLOT


def test_named_enthalpy_identity_includes_remainder() -> None:
    qe = freeze_energy_q()
    from src.solvers.wirsum_1997_coupled_f_energy import default_x0

    b = enthalpy_buckets_W(default_x0(qe), qe)
    assert "H_unreleased_vm_W" in b
    assert b["sum_W"] == pytest.approx(b["energy_W"], abs=_H_ID)
    nq = freeze_stage_q(5)
    x = default_x0_named(nq)
    x[-1] = 800.0
    nq_p = outer_P_named(x, nq)
    streams = named_stream_ledger(x, nq_p)
    assert nq_p.energy.species.m_unreleased_kg_s > 0.0
    assert streams["H_enthalpy_gap_W"] == pytest.approx(
        abs(streams["energy_row_W"]), abs=_H_ID
    )
    assert streams["carbon_basis"] == CARBON_BASIS_FUEL_WAF
