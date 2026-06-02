from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_kinetics import build_reaction_sources


def test_o2_limiter_is_phase_specific() -> None:
    n_gas = 11
    idx_co = 0
    idx_o2 = 5

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_b[idx_co] = 8.0
    c_b[idx_o2] = 5.0
    y_b[idx_co] = 0.45
    y_b[idx_o2] = 0.55

    c_d[idx_co] = 0.5
    c_d[idx_o2] = 0.3
    y_d[idx_co] = 0.4
    y_d[idx_o2] = 0.6

    bundle = build_reaction_sources(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.15,
        V_d=0.15,
        C_b=c_b,
        C_d=c_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )

    assert bundle.limit_factor_o2_bubble < 1e-8
    assert bundle.limit_factor_o2_dense > 0.0
    assert abs(bundle.R_gas_b[idx_co]) < 1e-5
    assert bundle.R_gas_d[idx_co] < 0.0


def test_o2_exchange_counts_as_dense_phase_supply() -> None:
    n_gas = 11
    idx_co = 0
    idx_o2 = 5

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_d[idx_co] = 0.5
    c_d[idx_o2] = 0.3
    y_d[idx_co] = 0.4
    y_d[idx_o2] = 0.6

    common_kwargs = dict(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.15,
        V_d=0.15,
        C_b=c_b,
        C_d=c_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.zeros(n_gas, dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )

    no_exchange = build_reaction_sources(
        **common_kwargs,
        N_ex=np.zeros(n_gas, dtype=np.float64),
    )
    with_exchange = build_reaction_sources(
        **common_kwargs,
        N_ex=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
    )

    assert with_exchange.limit_factor_o2_dense > no_exchange.limit_factor_o2_dense
    assert with_exchange.R_gas_d[idx_co] < no_exchange.R_gas_d[idx_co]


def test_o2_diagnostics_capture_stranded_bubble_slack_vs_dense_shortfall() -> None:
    n_gas = 11
    idx_co = 0
    idx_o2 = 5

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_b[idx_o2] = 1.0
    y_b[idx_o2] = 1.0
    c_d[idx_co] = 8.0
    c_d[idx_o2] = 5.0
    y_d[idx_co] = 0.45
    y_d[idx_o2] = 0.55

    bundle = build_reaction_sources(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.15,
        V_d=0.15,
        C_b=c_b,
        C_d=c_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )

    assert bundle.o2_slack_bubble > 0.0
    assert bundle.o2_shortfall_dense > 0.0
    assert bundle.o2_transfer_potential_bd == min(bundle.o2_slack_bubble, bundle.o2_shortfall_dense)


def test_net_molar_source_breakdown_matches_gas_source_sum_for_forward_r7_source() -> None:
    n_gas = 11
    idx_ch4 = 1
    idx_h2o = 2
    idx_co = 0
    idx_h2 = 4

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_d[idx_ch4] = 0.2
    c_d[idx_h2o] = 0.2
    c_d[idx_co] = 8.0
    c_d[idx_h2] = 12.0
    y_d[idx_ch4] = 0.01
    y_d[idx_h2o] = 0.09
    y_d[idx_co] = 0.35
    y_d[idx_h2] = 0.55

    bundle = build_reaction_sources(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.0,
        V_d=0.15,
        C_b=c_b,
        C_d=c_d,
        y_b=y_b,
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.array([0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )

    assert bundle.net_molar_gas_source_r7 > 0.0
    assert bundle.net_molar_gas_source_r5 == 0.0
    assert bundle.net_molar_gas_source_r6 == 0.0
    assert bundle.net_molar_gas_source_r10 == 0.0
    assert bundle.net_molar_gas_source_r11 == 0.0
    assert bundle.net_molar_gas_source_char == 0.0
    assert bundle.net_molar_gas_source_total == pytest.approx(
        float(np.sum(bundle.R_gas_b) + np.sum(bundle.R_gas_d))
    )
    assert bundle.net_molar_gas_source_total == pytest.approx(bundle.net_molar_gas_source_r7)
