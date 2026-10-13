"""阶段 C：无反应等温双向输运。15% 跳变与混合分母是工程行为，不是 Eq.(2-6)。"""

from __future__ import annotations

import numpy as np
import pytest

from src.core.cell import S_ASH, S_CHAR
from src.core.connectivity.solid_transport import (
    SOLID_TRANSPORT_INVENTORY_ENGINEERING,
    SOLID_TRANSPORT_INVENTORY_SPECIES_OWN,
    cell_resident_holdup_for_transport_closure,
    cell_total_solid_holdup,
    update_bed_solid_transport_coefficients,
    update_bed_solid_transport_inflows,
)
from src.core.reactor import Reactor, ReactorConfig
from src.core.solid_transport_reduced import (
    BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING,
    MIGRATION_MODEL_COHORT,
    MIGRATION_MODEL_UNVERIFIED_F_MIG_3,
    cohort_split_rates,
    engineering_closure_jump_kg,
    engineering_closure_mass_kg,
    mixed_denominator_gap,
    radius_retreat_m_per_s,
    reconstruct_k_down_species_own,
    single_particle_crossing,
    species_own_identity,
    streams_from_species_km,
    two_bin_mass_residual,
    unverified_f_mig_3,
)


def test_engineering_closure_jumps_at_15_percent() -> None:
    geo = 187.958
    jump = engineering_closure_jump_kg(geo)
    assert jump["closure_below_threshold_kg"] == pytest.approx(
        (BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING - 1.0e-6) * geo
    )
    assert jump["closure_at_threshold_kg"] == pytest.approx(geo)
    assert jump["jump_ratio"] == pytest.approx(1.0 / BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING, rel=1e-4)


def test_upper_cells_never_use_active_inventory_in_engineering_mode() -> None:
    geo, active = 191.261, 6.113
    assert engineering_closure_mass_kg(geo, active, is_bottom=False) == pytest.approx(geo)
    assert engineering_closure_mass_kg(geo, active, is_bottom=True) == pytest.approx(active)


def test_species_own_km_matches_reconstructed_downflow() -> None:
    M = np.array([[4.08, 1.13], [4.79, 1.32], [4.63, 1.28]], dtype=np.float64)
    K_up = np.array([0.079, 0.079, 0.032], dtype=np.float64)
    feed = np.array([[0.322, 0.089], [0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    ident = species_own_identity(M, K_up, feed)
    assert ident["km_vs_down_gap_kg_s"] == pytest.approx(0.0, abs=1e-12)
    assert ident["bottom_down_kg_s"] == pytest.approx(0.0)
    assert ident["max_k_down_bottom"] == pytest.approx(0.0)
    assert ident["max_abs_residual_upper_kg_s"] == pytest.approx(0.0, abs=1e-12)
    assert ident["max_up_edge_gap_kg_s"] == pytest.approx(0.0)
    assert ident["max_down_edge_gap_kg_s"] == pytest.approx(0.0)


def test_engineering_mixed_denominator_does_not_recover_component_flow() -> None:
    M = np.array([[4.08, 1.13], [4.79, 1.32], [4.63, 1.28]], dtype=np.float64)
    M_closure = np.array([5.204, 191.261, 184.872], dtype=np.float64)
    K_up = np.array([0.079, 0.079, 0.032], dtype=np.float64)
    feed = np.array([[0.322, 0.089], [0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    gap = mixed_denominator_gap(M, M_closure, K_up, feed)
    assert gap["max_abs_gap_kg_s"] > 0.1


def test_species_own_streams_have_no_bottom_penetration() -> None:
    M = np.ones((3, 2), dtype=np.float64)
    K_up = np.array([0.1, 0.1, 0.1], dtype=np.float64)
    K_down = reconstruct_k_down_species_own(
        M, K_up, np.array([[1.0, 0.2], [0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    )
    streams = streams_from_species_km(M, K_up, K_down)
    assert float(np.sum(streams.m_down[0])) == pytest.approx(0.0)


def _two_bed(*, inventory_mode: str, no_penetration: bool = False) -> Reactor:
    reactor = Reactor(
        ReactorConfig(
            n_cells=2,
            H_bed=2.0,
            D_bed=1.0,
            solid_transport_inventory_mode=inventory_mode,
            solid_transport_bottom_no_penetration=no_penetration,
        )
    )
    for cell in reactor.cells:
        cell.solid_state_model = "holdup_transport"
        cell.eps_b = 0.2
        cell.u_b = 0.5
        cell.eps_d_voidage = 0.5
        cell.m_solid.fill(0.0)
        cell.m_solid_zu.fill(0.0)
        cell.m_solid_rez.fill(0.0)
        cell.m_solid_in.fill(0.0)
        cell.R_solid.fill(0.0)
    return reactor


def test_production_default_still_uses_15_percent_bottom_switch() -> None:
    reactor = _two_bed(inventory_mode=SOLID_TRANSPORT_INVENTORY_ENGINEERING)
    lower, upper = reactor.cells
    lower.m_solid[0, S_CHAR] = 0.8
    lower.m_solid[0, S_ASH] = 0.1
    geo = cell_total_solid_holdup(lower)
    resident = cell_resident_holdup_for_transport_closure(lower, is_bottom_bed_cell=True)
    assert resident == pytest.approx(0.9)
    assert geo > resident * 10.0
    upper_geo = cell_resident_holdup_for_transport_closure(upper, is_bottom_bed_cell=False)
    assert upper_geo == pytest.approx(cell_total_solid_holdup(upper))


def test_species_own_mode_uses_char_ash_on_every_cell_and_km_identity() -> None:
    reactor = _two_bed(
        inventory_mode=SOLID_TRANSPORT_INVENTORY_SPECIES_OWN,
        no_penetration=True,
    )
    lower, upper = reactor.cells
    lower.m_solid[0, S_CHAR] = 4.0
    lower.m_solid[0, S_ASH] = 1.0
    upper.m_solid[0, S_CHAR] = 5.0
    upper.m_solid[0, S_ASH] = 1.2
    lower.m_solid_zu[0, S_CHAR] = 0.3
    update_bed_solid_transport_coefficients(reactor.cells, reactor.config)
    update_bed_solid_transport_inflows(reactor.cells)

    assert cell_resident_holdup_for_transport_closure(
        lower, is_bottom_bed_cell=True, inventory_mode=SOLID_TRANSPORT_INVENTORY_SPECIES_OWN
    ) == pytest.approx(5.0)
    assert cell_resident_holdup_for_transport_closure(
        upper, is_bottom_bed_cell=False, inventory_mode=SOLID_TRANSPORT_INVENTORY_SPECIES_OWN
    ) == pytest.approx(6.2)

    assert float(np.sum(lower.K_solid_ab[:, S_CHAR])) == pytest.approx(0.0)
    applied_down = upper.K_solid_ab * np.maximum(upper.m_solid, 0.0)
    np.testing.assert_allclose(lower.m_solid_ab_in, applied_down)
    applied_up = lower.K_solid_auf * np.maximum(lower.m_solid, 0.0)
    np.testing.assert_allclose(upper.m_solid_auf_in, applied_up)


def test_strict_config_selects_species_own_and_no_bottom_penetration() -> None:
    reactor = Reactor(ReactorConfig(n_cells=2, nr_solver_mode="wirsum_1997_strict"))
    assert reactor.config.solid_transport_inventory_mode == SOLID_TRANSPORT_INVENTORY_SPECIES_OWN
    assert reactor.config.solid_transport_bottom_no_penetration is True
    assert reactor.config.holdup_k0_anchor == "none"


def test_single_particle_crossing_time_and_mass_ratio() -> None:
    d0, d1, u = 1.0e-3, 2.0e-4, 1.0e-6
    path = single_particle_crossing(d0, d1, u)
    assert path["t_cross_s"] == pytest.approx((d0 - d1) / (2.0 * u))
    assert path["m_remaining_over_m0"] == pytest.approx((d1 / d0) ** 3)
    mid = 5.0e-4
    t_sum = (
        single_particle_crossing(d0, mid, u)["t_cross_s"]
        + single_particle_crossing(mid, d1, u)["t_cross_s"]
    )
    assert t_sum == pytest.approx(path["t_cross_s"])
    rem = (
        single_particle_crossing(d0, mid, u)["m_remaining_over_m0"]
        * single_particle_crossing(mid, d1, u)["m_remaining_over_m0"]
    )
    assert rem == pytest.approx(path["m_remaining_over_m0"])


def test_cohort_splits_consume_and_migrate() -> None:
    M, d0, d1, u = 2.0, 1.0e-3, 2.5e-4, 5.0e-7
    split = cohort_split_rates(M, d_from=d0, d_to=d1, u_radius=u)
    t = split["t_cross_s"]
    assert split["m_consume_kg_s"] + split["m_migrate_kg_s"] == pytest.approx(M / t)
    assert split["m_migrate_kg_s"] / (M / t) == pytest.approx((d1 / d0) ** 3)
    f3 = unverified_f_mig_3(u, d0, d1) * M
    assert f3 != pytest.approx(split["m_migrate_kg_s"], rel=1e-3)


def test_two_bin_ash_preserved_and_feed_closes() -> None:
    d = np.array([2.0e-4, 1.0e-3], dtype=np.float64)
    M = np.array([0.4, 1.6], dtype=np.float64)
    A = np.array([0.1, 0.3], dtype=np.float64)
    rho = 1400.0
    R = np.array([0.01, 0.05], dtype=np.float64)
    u = radius_retreat_m_per_s(float(R[1]), M_kg=float(M[1]), rho=rho, d_p=float(d[1]))
    split = cohort_split_rates(float(M[1]), d_from=float(d[1]), d_to=float(d[0]), u_radius=u)
    feed_c = np.array([0.01 - split["m_migrate_kg_s"], split["m_consume_kg_s"] + split["m_migrate_kg_s"]])
    feed_a = np.array([0.0, 0.02], dtype=np.float64)
    out = two_bin_mass_residual(
        M, A, d_p=d, rho=rho, feed_char=feed_c, feed_ash=feed_a, R_char=R, model=MIGRATION_MODEL_COHORT
    )
    np.testing.assert_allclose(out["F_char_kg_s"], 0.0, atol=1e-12)
    assert out["ash_gap_kg_s"] == pytest.approx(0.0, abs=1e-12)
    assert out["not_production_coeff_3"] is True
    alt = two_bin_mass_residual(
        M, A, d_p=d, rho=rho, feed_char=feed_c, feed_ash=feed_a, R_char=R,
        model=MIGRATION_MODEL_UNVERIFIED_F_MIG_3,
    )
    assert alt["m_migrate_char_kg_s"] != pytest.approx(out["m_migrate_char_kg_s"], rel=1e-4)
