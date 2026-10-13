from __future__ import annotations

import numpy as np
import pytest

from src.core.cell_kinetics import build_reaction_sources
from src.core.species import GAS_SPECIES_INDEX


def test_pyrolysis_gas_source_enters_dense_phase_not_bubble_phase() -> None:
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    gas_src_vm = np.zeros(n_gas, dtype=np.float64)
    gas_src_vm[idx["CO"]] = 2.5
    gas_src_vm[idx["CH4"]] = 0.75
    gas_src_vm[idx["H2O"]] = 1.2

    bundle = build_reaction_sources(
        T=1100.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.1,
        V_d=0.1,
        C_b=np.zeros(n_gas, dtype=np.float64),
        C_d=np.zeros(n_gas, dtype=np.float64),
        y_b=np.zeros(n_gas, dtype=np.float64),
        y_d=np.zeros(n_gas, dtype=np.float64),
        gas_src_vm=gas_src_vm,
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.0,
        rho_cat=0.0,
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.zeros(n_gas, dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )

    assert bundle.R_gas_d[idx["CO"]] == pytest.approx(gas_src_vm[idx["CO"]])
    assert bundle.R_gas_d[idx["CH4"]] == pytest.approx(gas_src_vm[idx["CH4"]])
    assert bundle.R_gas_d[idx["H2O"]] == pytest.approx(gas_src_vm[idx["H2O"]])
    assert np.max(np.abs(bundle.R_gas_b)) == pytest.approx(0.0)


def test_o2_limiter_is_phase_specific() -> None:
    """工程限幅：气泡 R5 燃料充足、O₂ 近零 → lf_b≪1；密相 R5 悬相速率小 → lf_d≈1。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    idx_co = idx["CO"]
    idx_h2o = idx["H2O"]
    idx_o2 = idx["O2"]

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_b[idx_co] = 8.0
    c_b[idx_h2o] = 1.0
    c_b[idx_o2] = 5.0
    y_b[idx_co] = 0.45
    y_b[idx_h2o] = 0.10
    y_b[idx_o2] = 0.55

    c_d[idx_co] = 0.5
    c_d[idx_h2o] = 0.5
    c_d[idx_o2] = 0.3
    y_d[idx_co] = 0.4
    y_d[idx_h2o] = 0.2
    y_d[idx_o2] = 0.6

    # 燃料预封顶：气泡须有 CO 进料，否则 r5b→0、lf_b 恒为 1；无 O₂ 进料 → budget≈floor
    n_zu_b = np.zeros(n_gas, dtype=np.float64)
    n_zu_b[idx_co] = 20.0
    n_zu_d = np.zeros(n_gas, dtype=np.float64)
    n_zu_d[idx_o2] = 0.2
    n_zu_d[idx_co] = 2.0

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
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=n_zu_d,
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=n_zu_b,
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
        extent_limiters_enabled=True,
    )

    # O₂ 供给 floor=1e-6 → lf_b ≈ 1e-7（燃料预封顶后 dem_b≈9.95）
    assert bundle.limit_factor_o2_bubble < 1e-6
    assert bundle.limit_factor_o2_dense > 0.99
    assert abs(bundle.R_gas_b[idx_co]) < 1e-5
    assert bundle.R_gas_d[idx_co] < 0.0


def test_o2_exchange_counts_as_dense_phase_supply() -> None:
    """密相用 R6（CH4）造真实 O₂ 需求；N_ex(O₂) 抬高 budget → lf_d 上升。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    idx_ch4 = idx["CH4"]
    idx_o2 = idx["O2"]

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_d[idx_ch4] = 5.0
    c_d[idx_o2] = 5.0
    y_d[idx_ch4] = 0.5
    y_d[idx_o2] = 0.5

    n_zu_d = np.zeros(n_gas, dtype=np.float64)
    n_zu_d[idx_ch4] = 10.0
    n_zu_d[idx_o2] = 0.05

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
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=n_zu_d,
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
        extent_limiters_enabled=True,
    )

    no_exchange = build_reaction_sources(
        **common_kwargs,
        N_ex=np.zeros(n_gas, dtype=np.float64),
    )
    n_ex = np.zeros(n_gas, dtype=np.float64)
    n_ex[idx_o2] = 0.5
    with_exchange = build_reaction_sources(
        **common_kwargs,
        N_ex=n_ex,
    )

    assert no_exchange.limit_factor_o2_dense < 0.1
    assert with_exchange.limit_factor_o2_dense > no_exchange.limit_factor_o2_dense
    # 更多 O₂ → R6 消耗更多 CH4
    assert with_exchange.R_gas_d[idx_ch4] < no_exchange.R_gas_d[idx_ch4]


def test_o2_diagnostics_capture_stranded_bubble_slack_vs_dense_shortfall() -> None:
    """气泡仅 O₂（slack）；密相 CH4 燃料充足、O₂ 不足（shortfall）。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    idx_ch4 = idx["CH4"]
    idx_o2 = idx["O2"]

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_b[idx_o2] = 1.0
    y_b[idx_o2] = 1.0
    c_d[idx_ch4] = 5.0
    c_d[idx_o2] = 5.0
    y_d[idx_ch4] = 0.5
    y_d[idx_o2] = 0.5

    n_zu_d = np.zeros(n_gas, dtype=np.float64)
    n_zu_d[idx_ch4] = 10.0
    n_zu_d[idx_o2] = 0.05
    n_zu_b = np.zeros(n_gas, dtype=np.float64)
    n_zu_b[idx_o2] = 0.2
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
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=n_zu_d,
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=n_zu_b,
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
        extent_limiters_enabled=True,
    )

    assert bundle.o2_slack_bubble > 0.0
    assert bundle.o2_shortfall_dense > 0.0
    assert bundle.o2_transfer_potential_bd == min(bundle.o2_slack_bubble, bundle.o2_shortfall_dense)


def test_fuel_precap_blocks_phantom_o2_demand_from_high_tar_concentration() -> None:
    """高 C_TAR、零 tar 进料：燃料预封顶后 O₂ 需求不得由 raw R10 主导。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    idx_tar1 = idx["TAR1"]
    idx_o2 = idx["O2"]

    c_d = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)
    c_d[idx_tar1] = 2.0
    c_d[idx_o2] = 5.0
    y_d[idx_tar1] = 0.3
    y_d[idx_o2] = 0.7

    n_zu_d = np.zeros(n_gas, dtype=np.float64)
    n_zu_d[idx_o2] = 1.0  # 有 O₂、无 tar 进料

    bundle = build_reaction_sources(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.15,
        V_d=0.15,
        C_b=np.zeros(n_gas, dtype=np.float64),
        C_d=c_d,
        y_b=np.zeros(n_gas, dtype=np.float64),
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=n_zu_d,
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
        extent_limiters_enabled=True,
    )

    # tar_supply floor=1e-6 → R10 O₂ 贡献 ≤ ν·0.995·1e-6，不得出现 ~10⁶ 伪需求
    assert bundle.o2_demand_dense_raw < 1.0
    assert bundle.limit_factor_o2_dense > 0.99


def test_hamel_raw_kinetics_passthrough_no_extent_limiters() -> None:
    """Hamel 同构默认：无燃料预封顶/比例 lf；写入 extent ≈ 动力学；R12 ≫ R5。"""
    from src.kinetics.gas_reactions import rate_R5_bubble, rate_R12

    n_gas = 11
    idx = GAS_SPECIES_INDEX
    c_b = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    for sp, v in (("H2", 1.0), ("CO", 1.0), ("O2", 1.0), ("H2O", 1.0)):
        c_b[idx[sp]] = v
        y_b[idx[sp]] = 0.25
    n_zu_b = np.zeros(n_gas, dtype=np.float64)
    n_zu_b[idx["H2"]] = 0.01
    n_zu_b[idx["CO"]] = 0.01
    n_zu_b[idx["O2"]] = 0.01
    T, P, V_b = 1200.0, 2.5e6, 0.2
    kin_r12 = rate_R12(T, c_b[idx["H2"]], c_b[idx["O2"]], P, y_b, idx) * V_b
    kin_r5 = rate_R5_bubble(
        T, c_b[idx["CO"]], c_b[idx["O2"]], c_b[idx["H2O"]], P, y_b, idx
    ) * V_b
    assert kin_r12 > 1e3 * kin_r5

    bundle = build_reaction_sources(
        T=T,
        P=P,
        fuel_type="coal",
        V_b=V_b,
        V_d=0.2,
        C_b=c_b,
        C_d=np.zeros(n_gas, dtype=np.float64),
        y_b=y_b,
        y_d=np.zeros(n_gas, dtype=np.float64),
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.0,
        rho_cat=0.0,
        rho_ash_kg_m3=0.0,
        enable_r12=True,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=0.0,
        r2_scale=0.0,
        r5_scale=1.0,
        r6_scale=0.0,
        r7_scale=0.0,
        r8_scale=0.0,
        rate_multiplier=1.0,
        N_zu_d=np.zeros(n_gas, dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=n_zu_b,
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
        extent_limiters_enabled=False,
    )
    assert bundle.limit_factor_o2_bubble == pytest.approx(1.0)
    assert bundle.extent_r12_bubble == pytest.approx(kin_r12, rel=1e-9)
    assert bundle.extent_r5_bubble == pytest.approx(kin_r5, rel=1e-9)
    # 若误开燃料预封顶，R12 会被压到 ~0.01 量级
    assert bundle.extent_r12_bubble > 1.0


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
        rho_ash_kg_m3=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=1.0,
        r8_scale=1.0,
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


def test_char_kinetics_r2_surface_rate_independent_of_particle_diameter() -> None:
    """Hamel Eq. (5.24)–(5.25)：R2 面速率不含 d_p；等面积时各粒径类炭消耗相同。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX

    c_b = np.zeros(n_gas, dtype=np.float64)
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)

    c_d[idx["H2O"]] = 5.0
    y_d[idx["H2O"]] = 1.0

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
        solid_sink_vm=np.zeros((2, 4), dtype=np.float64),
        areas=np.ones(2, dtype=np.float64),
        solid_d_p=2.0e-3,
        solid_d_p_classes=np.array([1.0e-3, 3.0e-3], dtype=np.float64),
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=0.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.array([0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(2, 4),
        char_index=0,
    )

    assert bundle.extent_r2 > 0.0
    assert bundle.R_solid[0, 0] == pytest.approx(bundle.R_solid[1, 0])

    c_d[idx["O2"]] = 2.0
    bundle_r1 = build_reaction_sources(
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
        solid_sink_vm=np.zeros((2, 4), dtype=np.float64),
        areas=np.ones(2, dtype=np.float64),
        solid_d_p=2.0e-3,
        solid_d_p_classes=np.array([1.0e-3, 3.0e-3], dtype=np.float64),
        D_g=5.0e-5,
        char_conversion=0.5,
        rho_cat=0.0,
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=1.0,
        r5_scale=1.0,
        r6_scale=1.0,
        r7_scale=0.0,
        r8_scale=1.0,
        rate_multiplier=1.0,
        N_zu_d=np.array([0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(2, 4),
        char_index=0,
    )
    assert bundle_r1.R_solid[0, 0] != pytest.approx(bundle_r1.R_solid[1, 0])
    assert abs(bundle_r1.R_solid[0, 0]) > abs(bundle_r1.R_solid[1, 0])


def test_r4_uses_char_mass_not_area_hamel_eq533() -> None:
    """Boudouard extent 按 Eq.5.33 质量基；面积变化不改变 R4。"""
    from src.core.constants import Rg
    from src.kinetics.char_reactions import extent_R4

    n_gas = 11
    idx = GAS_SPECIES_INDEX
    c_d = np.zeros(n_gas, dtype=np.float64)
    y_d = np.zeros(n_gas, dtype=np.float64)
    T = 1200.0
    P = 2.5e6
    y_d[idx["CO2"]] = 0.12
    y_d[idx["CO"]] = 0.05
    y_d[idx["N2"]] = 0.83
    c_d[idx["CO2"]] = y_d[idx["CO2"]] * P / (Rg * T)
    c_d[idx["CO"]] = y_d[idx["CO"]] * P / (Rg * T)
    m_char = np.array([1.0, 2.0], dtype=np.float64)  # [kg]
    p_co2 = float(c_d[idx["CO2"]] * Rg * T)
    p_co = float(c_d[idx["CO"]] * Rg * T)
    expected = extent_R4(T, p_co2, p_co, float(np.sum(m_char)))

    common = dict(
        T=T,
        P=P,
        fuel_type="coal",
        V_b=0.0,
        V_d=0.15,
        C_b=np.zeros(n_gas, dtype=np.float64),
        C_d=c_d,
        y_b=np.zeros(n_gas, dtype=np.float64),
        y_d=y_d,
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((2, 4), dtype=np.float64),
        solid_d_p=2.0e-3,
        D_g=5.0e-5,
        char_conversion=0.3,
        rho_cat=0.0,
        rho_ash_kg_m3=75.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=1.0,
        r2_scale=0.0,
        r5_scale=0.0,
        r6_scale=0.0,
        r7_scale=0.0,
        r8_scale=0.0,
        rate_multiplier=1.0,
        N_zu_d=np.zeros(n_gas, dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=np.zeros(n_gas, dtype=np.float64),
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(2, 4),
        char_index=0,
        m_char_classes=m_char,
    )
    b_small_area = build_reaction_sources(areas=np.array([0.1, 0.1], dtype=np.float64), **common)
    b_large_area = build_reaction_sources(areas=np.array([10.0, 10.0], dtype=np.float64), **common)
    assert b_small_area.extent_r4 == pytest.approx(expected, rel=1e-9)
    assert b_large_area.extent_r4 == pytest.approx(expected, rel=1e-9)
    assert b_small_area.extent_r4 == pytest.approx(b_large_area.extent_r4)


def test_r5_hamel_co_basis_stoichiometry() -> None:
    """Hamel §5.2.1: CO + 1/2 O2 → CO2；extent = mol CO/s（非 2CO+O2）。"""
    n_gas = 11
    idx = GAS_SPECIES_INDEX
    c_b = np.zeros(n_gas, dtype=np.float64)
    y_b = np.zeros(n_gas, dtype=np.float64)
    c_b[idx["CO"]] = 2.0
    c_b[idx["O2"]] = 2.0
    c_b[idx["H2O"]] = 1.0
    y_b[idx["CO"]] = 0.4
    y_b[idx["O2"]] = 0.4
    y_b[idx["H2O"]] = 0.2
    n_zu_b = np.zeros(n_gas, dtype=np.float64)
    n_zu_b[idx["CO"]] = 50.0
    n_zu_b[idx["O2"]] = 50.0
    n_zu_b[idx["H2O"]] = 10.0

    bundle = build_reaction_sources(
        T=1200.0,
        P=2.5e6,
        fuel_type="coal",
        V_b=0.2,
        V_d=0.2,
        C_b=c_b,
        C_d=np.zeros(n_gas, dtype=np.float64),
        y_b=y_b,
        y_d=np.zeros(n_gas, dtype=np.float64),
        gas_src_vm=np.zeros(n_gas, dtype=np.float64),
        solid_sink_vm=np.zeros((1, 4), dtype=np.float64),
        areas=np.zeros(1, dtype=np.float64),
        solid_d_p=1.0e-3,
        D_g=5.0e-5,
        char_conversion=0.0,
        rho_cat=0.0,
        rho_ash_kg_m3=0.0,
        enable_r12=False,
        use_gibbs_minor=False,
        gibbs_minor_sources=None,
        r4_scale=0.0,
        r2_scale=0.0,
        r5_scale=1.0,
        r6_scale=0.0,
        r7_scale=0.0,
        r8_scale=0.0,
        rate_multiplier=1.0,
        N_zu_d=np.zeros(n_gas, dtype=np.float64),
        N_d_in=np.zeros(n_gas, dtype=np.float64),
        N_zu_b=n_zu_b,
        N_b_in=np.zeros(n_gas, dtype=np.float64),
        N_rez_d=np.zeros(n_gas, dtype=np.float64),
        N_rez_b=np.zeros(n_gas, dtype=np.float64),
        N_ex=np.zeros(n_gas, dtype=np.float64),
        solid_shape=(1, 4),
        char_index=0,
    )
    e5 = float(bundle.extent_r5_bubble)
    assert e5 > 0.0
    # ΔCO = −e5, ΔO2 = −0.5 e5, ΔCO2 = +e5
    assert bundle.R_gas_b[idx["CO"]] == pytest.approx(-e5, rel=1e-9)
    assert bundle.R_gas_b[idx["O2"]] == pytest.approx(-0.5 * e5, rel=1e-9)
    assert bundle.R_gas_b[idx["CO2"]] == pytest.approx(e5, rel=1e-9)
    assert bundle.o2_demand_bubble_limited == pytest.approx(0.5 * e5, rel=1e-9)
