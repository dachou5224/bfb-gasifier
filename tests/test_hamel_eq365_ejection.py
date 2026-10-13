"""式 (3.65) 的分段和式 (3.67) 的总抛射。数值锚在印刷页 40。"""

import math

import pytest

from src.physics.wirsum_1997_freeboard import (
    WAKE_F_W,
    WAKE_ZETA_W,
    calc_roof_mass_flux_eq359,
    calc_wake_mass_flux_engineering,
    f0_hamel_367,
    zeta_w_f_w_hamel_365,
)


def test_eq365_has_only_the_printed_breakpoint() -> None:
    assert zeta_w_f_w_hamel_365(0.0) == pytest.approx(0.1)
    assert zeta_w_f_w_hamel_365(9.999) == pytest.approx(0.1)
    at_10 = 0.5 * (1.0 - math.exp(9.5 - 10.0))
    assert zeta_w_f_w_hamel_365(10.0) == pytest.approx(at_10)
    assert zeta_w_f_w_hamel_365(10.0) == pytest.approx(0.1967346701436833)
    high = zeta_w_f_w_hamel_365(40.0)
    assert 0.1 < high < 0.5
    assert zeta_w_f_w_hamel_365(200.0) == pytest.approx(0.5)


def test_eq367_matches_the_printed_combination() -> None:
    common = dict(rho_p=1400.0, eps_d=0.45, eps_b=0.13, u_b=0.8)
    d_bar = 2.5e-4
    d_b = 0.21
    zeta = zeta_w_f_w_hamel_365(12.0)
    s_sum = 0.3
    pref = common["eps_b"] * common["u_b"] * common["rho_p"] * (1.0 - common["eps_d"])
    roof_factor = 3.0 * d_bar / d_b
    printed = pref * (roof_factor + (zeta - roof_factor) * s_sum)
    assert f0_hamel_367(
        d_bar=d_bar, d_b=d_b, zeta_f=zeta, S=s_sum, **common,
    ) == pytest.approx(printed)
    roof = calc_roof_mass_flux_eq359(d_p=d_bar, d_b=d_b, **common)
    wake = calc_wake_mass_flux_engineering(**common, zeta_w=zeta, f_w=1.0)
    assert printed == pytest.approx((1.0 - s_sum) * roof + s_sum * wake)


def test_probability_chain_keeps_the_printed_recurrence_and_open_base() -> None:
    from src.solvers.wesseling_case_contract import ejection_and_wesseling_open_items

    book = ejection_and_wesseling_open_items()
    assert book["eq_3_63_rewritten_as_geometric_power"] is False
    assert "P''_{k-1}" in book["eq_3_63_literal"]
    assert book["P_double_prime_1"] is None
    assert book["approach_angle_deg"]["selected"] is None
    assert book["live_hamel_extended_uses_eq_367"] is False
    assert book["steam_status"] == "missing"
    assert book["fuel_size_applied_to_bins"] is False
    assert book["secondary_air_placed"] is False


def test_pyrolysis_temperature_column_is_a_finite_coupling() -> None:
    from dataclasses import replace

    import numpy as np

    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_e8 import (
        _pyrolysis_temperature_column,
        freeze_e8_q,
    )

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    q = outer_P_named(x, nq)
    column = _pyrolysis_temperature_column(q, x)
    assert column.shape == x.shape
    assert np.all(np.isfinite(column))
    assert float(np.linalg.norm(column)) > 0.0


def test_tar_sum_fraction_moves_a_trace_component_without_changing_the_total() -> None:
    from src.solvers.wirsum_1997_named_f import (
        inventories_from_sum_fraction,
        tar_consumption_fractions,
        tar_sum_and_fraction,
    )

    a1, a2 = tar_consumption_fractions()
    total, frac = tar_sum_and_fraction(5.45e-12, 5.32e-6, a1, a2)
    n1, n2 = inventories_from_sum_fraction(total, frac, a1, a2)
    assert n1 + n2 == pytest.approx(total)
    assert n1 == pytest.approx(5.45e-12, rel=1.0e-6, abs=1.0e-18)
    shifted = inventories_from_sum_fraction(total, frac + 0.2, a1, a2)
    assert shifted[0] + shifted[1] == pytest.approx(total)
    assert shifted[0] > 1.0e-7
    assert min(shifted) >= 0.0
    assert inventories_from_sum_fraction(0.0, 0.4, a1, a2) == (0.0, 0.0)
    backward = tar_sum_and_fraction(5.32e-6, 5.45e-12, a1, a2)
    other = inventories_from_sum_fraction(backward[0], backward[1] - 0.2, a1, a2)
    assert other[0] + other[1] == pytest.approx(backward[0])
    assert other[1] > 1.0e-7
    assert min(other) >= 0.0


def test_tar_row_transform_keeps_the_original_zero_and_drops_reaction() -> None:
    from src.solvers.wirsum_1997_named_f import tar_consumption_fractions, tar_transport_ratio

    a1, a2 = tar_consumption_fractions()
    assert a1 + a2 == pytest.approx(1.0)
    transport = tar_transport_ratio(0.4, -0.2, 0.1, 0.05, a1=a1, a2=a2)
    xi = 3.5
    f1 = (0.4 - 0.1) + a1 * xi
    f2 = (-0.2 - 0.05) + a2 * xi
    assert a2 * f1 - a1 * f2 == pytest.approx(transport)
    total = f1 + f2
    # 逆变换回到原残差
    back1 = a1 * total + transport
    back2 = a2 * total - transport
    assert back1 == pytest.approx(f1)
    assert back2 == pytest.approx(f2)
    assert tar_transport_ratio(0.0, 0.0, 0.0, 0.0, a1=a1, a2=a2) == 0.0


def test_undeveloped_bubble_tar_uses_linear_coordinates() -> None:
    import numpy as np

    from src.solvers.wirsum_1997_named_f_e8 import (
        _e8_positive_log_ids,
        _linear_undeveloped_bubble_tar_ids,
    )

    x = np.zeros(50)
    x[42] = 1.0e-18
    x[43] = 1.0e-18
    f = np.zeros(50)
    f[42] = -1.3e-5
    f[43] = -2.0e-5
    lin = _linear_undeveloped_bubble_tar_ids(x, f)
    assert set(int(i) for i in lin) == {42, 43}
    kept = _e8_positive_log_ids(x)
    kept = kept[~np.isin(kept, lin)]
    assert 42 not in kept and 43 not in kept
    x[21] = 1.0e-16
    f[21] = 1.0e-6
    assert 21 not in set(int(i) for i in _linear_undeveloped_bubble_tar_ids(x, f))


def test_bubble_h2_zero_enters_the_joint_problem_without_a_seed() -> None:
    from dataclasses import replace

    from src.solvers.wirsum_1997_named_f import (
        coupled_residual_named,
        default_x0_named,
        outer_P_named,
    )
    from src.solvers.wirsum_1997_named_f_e8 import (
        _linear_bubble_zero_ids,
        freeze_e8_q,
    )

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = default_x0_named(nq)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q, rate_scale=1.0)
    assert 9 in _linear_bubble_zero_ids(x, f)
    assert float(x[9]) == 0.0


def test_extended_wake_product_is_only_the_low_speed_branch() -> None:
    product = WAKE_ZETA_W * WAKE_F_W
    assert product == pytest.approx(0.1)
    assert product == pytest.approx(zeta_w_f_w_hamel_365(2.0))
    assert product != pytest.approx(zeta_w_f_w_hamel_365(20.0))
