"""F_E8：K_eject 进 P(x)。不改旋风 K，不宣称炉温，不开放生产 refresh。"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from src.core.reactor import ReactorConfig
from src.solvers.wirsum_1997_named_f import GEO_FB_RHO_G_T_REF_K, default_x0_named, outer_P_named
from src.solvers.wirsum_1997_named_f_e8 import (
    e8_source_snapshot,
    freeze_e8_q,
    residual_named_e8,
    solve_named_e8_from_e7,
)
from src.solvers.wirsum_1997_named_f_endpoint_io import (
    SNAPSHOT_PATH,
    load_named_e8_snapshot,
    recertify_named_e7,
    recertify_named_e8,
)
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_registry import n_unknowns, omitted_named


def _snapshot() -> dict:
    if not SNAPSHOT_PATH.is_file():
        pytest.skip(
            f"missing {SNAPSHOT_PATH.name}; run data/named_e8_recompute_2026-09-17.py"
        )
    return load_named_e8_snapshot()


def _e7_from_snapshot() -> dict:
    snap = _snapshot()
    x = np.asarray(snap["e7"]["x"], dtype=np.float64)
    cert = recertify_named_e7(x)
    if not cert["endpoint_certified"]:
        pytest.fail("saved E7 x failed live recertify; rerun named_e8_recompute")
    return {
        "x": x,
        "T_K": float(x[-1]),
        "physical_convergence": True,
        "endpoint_certified": True,
    }


def test_e8_outer_p_updates_rho_g_not_cyclone_k() -> None:
    nq8 = freeze_e8_q(Q_W=0.0)
    assert n_unknowns(8) == 45
    assert "K_eject_in_this_P" not in omitted_named(8)
    assert omitted_named(8)["cyclone_K_update_P"] == "omitted, named"
    assert omitted_named(8)["production_hamel_outer_refresh"] == "omitted, named"
    x = default_x0_named(nq8)
    x[-1] = 1500.0
    x_copy = x.copy()
    q2 = outer_P_named(x, nq8)
    np.testing.assert_allclose(x, x_copy)
    assert q2.energy.species.k == nq8.energy.species.k
    assert q2.energy.species.geo_fb.rho_g != nq8.energy.species.geo_fb.rho_g
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.vorab_outer_source_policy != "hamel_outer_refresh"


def test_e8_coarse_fate_can_flip_below_t_ref() -> None:
    nq8 = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq8)
    x_hi = x.copy()
    x_hi[-1] = float(GEO_FB_RHO_G_T_REF_K) + 300.0
    x_lo = x.copy()
    x_lo[-1] = 800.0
    q_hi = outer_P_named(x_hi, nq8)
    q_lo = outer_P_named(x_lo, nq8)
    assert float(q_hi.energy.species.w_leave[0]) >= float(q_hi.energy.species.w_leave[1])
    assert float(q_hi.energy.species.w_leave[1] + q_hi.energy.species.w_return[1]) <= 1.0 + 1.0e-12
    assert q_hi.energy.species.k == nq8.energy.species.k
    # 800 K 粗档未判定不是逃逸。
    assert q_lo.hydro is not None
    assert float(q_lo.energy.species.w_leave[1]) < 0.5
    np.testing.assert_allclose(
        q_lo.energy.species.w_leave + q_lo.energy.species.w_return + q_lo.hydro.w_undecided,
        1.0,
        atol=1.0e-12,
    )


def test_e8_at_e7_root_hydro_updates_and_keeps_cyclone_k() -> None:
    e7 = _e7_from_snapshot()
    snap = e8_source_snapshot(e7["x"], freeze_e7_q(Q_W=0.0), freeze_e8_q(Q_W=0.0))
    assert snap["cyclone_k_unchanged"] is True
    assert snap["rho_g8"] != snap["rho_g7"]
    nq8 = freeze_e8_q(Q_W=0.0)
    q8 = outer_P_named(np.asarray(e7["x"], dtype=np.float64), nq8)
    assert q8.hydro is not None
    np.testing.assert_allclose(
        q8.energy.species.w_leave + q8.energy.species.w_return + q8.hydro.w_undecided,
        1.0,
        atol=1.0e-12,
    )


def test_e8_saved_endpoint_recertifies_live_not_cached_boolean() -> None:
    snap = _snapshot()
    assert "endpoint_certified" not in snap
    x = np.asarray(snap["e8"]["x"], dtype=np.float64)
    stored = bool(snap["e8"]["endpoint_certified_at_write"])
    cert = recertify_named_e8(x)
    assert cert["endpoint_certified"] is True
    assert stored is True
    assert cert["gas_mol_s_l2"] < 1.0e-6
    assert cert["energy_W_abs"] < 1.0
    np.testing.assert_allclose(cert["F"], snap["e8"]["F"], rtol=1.0e-9, atol=1.0e-9)
    q = snap["e8"]["q_P"]
    assert q["cyclone_k"]["model"] == "frozen_context"
    live_q = e8_source_snapshot(x, freeze_e7_q(Q_W=0.0), freeze_e8_q(Q_W=0.0))
    np.testing.assert_allclose(q["K_eject"], live_q["K_eject8"])
    np.testing.assert_allclose(q["w_leave"], live_q["w_leave8"])


def test_e8_solve_does_not_change_rate_or_claim_furnace() -> None:
    e7 = _e7_from_snapshot()
    sol = solve_named_e8_from_e7(Q_W=0.0, max_iter=20, e7=e7)
    assert sol["rate_formula_changed"] is False
    assert sol["f_name"] == "F_E8"
    assert sol["n_unknowns"] == 45
    assert sol["K_eject_in_this_P"] is True
    assert sol["finite_hydro_P"] is True
    assert sol["cyclone_K_update_P"] is False
    assert sol["freeboard_reaction"] is True
    assert sol["Q_W_is_fitted_to_1100K"] is False
    assert sol["not_a_furnace_convergence_claim"] is True
    assert sol["production_hamel_outer_refresh"] is False
    assert sol["x_nonneg_projection"] is False
    assert sol["coupling_mode"] == "composite_P"
    assert sol["endpoint_cert"]["clip_applied"] is False
    assert sol["omitted_named"]["cyclone_K_update_P"] == "omitted, named"
    live = recertify_named_e8(sol["x"])
    assert live["endpoint_certified"] is True
    assert live["gas_mol_s_l2"] < 1.0e-6
    assert live["energy_W_abs"] < 1.0
    assert sol["equilibrated_rank"] == 45
    assert abs(float(sol["T_K"]) - 1100.0) > 50.0
    assert float(sol["w_leave"][0]) + float(sol["w_leave"][1]) <= 2.0
    f8 = residual_named_e8(sol["x"], freeze_e8_q(Q_W=0.0))
    assert float(np.linalg.norm(f8[0:14])) < 1.0e-6
    cfg = ReactorConfig()
    assert cfg.thesis_vorab_sources_single_shot is True
    assert cfg.vorab_outer_source_policy != "hamel_outer_refresh"
    ev = sol["linear_step_last"]
    assert set(ev) >= {"raw_J", "modified_matrix", "feasible_direction", "rollback"}
    assert ev["rollback"]["returned_state_source"] in ("newton", "e7_certified_root")


def test_block_gate_accepts_energy_above_mixed_l2() -> None:
    """尺度 0.5 的合格态能量约 5e-5 W，整体无权 L2 过不了 1e-6，分块门要通过。"""
    from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
    from src.solvers.wirsum_1997_named_f_e8 import as_e8_state

    f = np.zeros(45, dtype=np.float64)
    f[-1] = 5.3814e-5
    blocks = named_endpoint_blocks(f)
    assert blocks["passed"] is True
    assert float(np.linalg.norm(f)) > 1.0e-6
    f[0] = 2.0e-6
    assert named_endpoint_blocks(f)["passed"] is False
    z7 = np.linspace(0.1, 1.4, 43)
    z8 = as_e8_state(z7)
    assert z8.size == 45
    assert z8[42] == 0.0 and z8[43] == 0.0
    assert z8[-1] == pytest.approx(z7[-1])


def test_e11_freeboard_char_holdup_matches_leaving_flux() -> None:
    """无 R2–R4 时，自由空间炭库存等于抛出流乘停留时间，旋风入口不另扣消耗。"""
    from src.solvers.wirsum_1997_named_f import (
        coupled_residual_named,
        fb_char_balance,
        outer_P_named,
    )
    from src.solvers.wirsum_1997_named_f_registry import f_name

    x10 = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e10.npz")["x"]
    x11 = np.concatenate([x10[:-2], [0.0, 0.0], x10[-2:]])
    assert n_unknowns(11) == 50
    assert f_name(11) == "F_E11"
    nq = replace(
        freeze_e8_q(Q_W=0.0),
        stage=11,
        enabled_rxns=frozenset({"R5", "R6", "R8", "R12", "R10", "R11"}),
    )
    q = outer_P_named(x11, nq)
    resid, _src, cons = fb_char_balance(x11, q, rate_scale=1.0)
    assert float(np.sum(cons)) == pytest.approx(0.0, abs=1.0e-12)
    tau = np.asarray(q.hydro.tau_s, dtype=np.float64)
    m_eject = resid.copy()
    assert float(m_eject[1]) > 0.0
    x11[46] = float(m_eject[0]) * float(tau[0])
    x11[47] = float(m_eject[1]) * float(tau[1])
    assert x11[47] > 0.0
    f = coupled_residual_named(x11, outer_P_named(x11, nq))
    assert f[46] == pytest.approx(0.0, abs=1.0e-8)
    assert f[47] == pytest.approx(0.0, abs=1.0e-8)
    x_rx = x11.copy()
    x_rx[19] = 0.2
    x_rx[46] = 0.05
    x_rx[47] = 0.05
    nq_r = replace(nq, enabled_rxns=nq.enabled_rxns | frozenset({"R2", "R3", "R4"}))
    _r0, _s0, c0 = fb_char_balance(x11, outer_P_named(x11, nq_r), rate_scale=1.0)
    _r1, _s1, c1 = fb_char_balance(x_rx, outer_P_named(x_rx, nq_r), rate_scale=1.0)
    assert float(np.sum(c1)) > float(np.sum(c0))


def test_e10_freeboard_tar_closes_transport_and_oxidizes() -> None:
    """自由空间焦油等于床层两相之和时，无反应则输运残差为零；有 R10 时消耗自由空间氧。"""
    from src.solvers.wirsum_1997_named_f import coupled_residual_named
    from src.solvers.wirsum_1997_named_f_registry import f_name

    snap = _snapshot()
    x8 = np.asarray(snap["e8"]["x"], dtype=np.float64)
    inlet1 = float(x8[40] + x8[42])
    inlet2 = float(x8[41] + x8[43])
    T = float(x8[-1])
    x10 = np.concatenate([x8[:-1], [inlet1, inlet2, T, T]])
    assert x10.size == 48
    assert n_unknowns(10) == 48
    assert n_unknowns(9) == 46
    assert f_name(10) == "F_E10"
    nq = replace(
        freeze_e8_q(Q_W=0.0),
        stage=10,
        enabled_rxns=frozenset({"R5", "R6", "R8", "R12"}),
    )
    f = coupled_residual_named(x10, nq)
    assert f[44] == pytest.approx(0.0, abs=1.0e-12)
    assert f[45] == pytest.approx(0.0, abs=1.0e-12)
    nq_r = replace(nq, enabled_rxns=frozenset({"R5", "R6", "R8", "R12", "R10", "R11"}))
    f_r = coupled_residual_named(x10, nq_r)
    assert abs(float(f_r[44])) > 1.0e-8
    assert abs(float(f_r[45])) > 1.0e-8
    assert abs(float(f_r[16] - f[16])) > 1.0e-8


def test_e9_two_temperatures_keep_e8_species_when_zones_match() -> None:
    """E9 增加自由空间温度。两温相等时物种行与 E8 相同，两区能量之和为原总能量。"""
    from src.solvers.wirsum_1997_named_f import coupled_residual_named, energy_row_named
    from src.solvers.wirsum_1997_named_f_registry import f_name

    snap = _snapshot()
    x8 = np.asarray(snap["e8"]["x"], dtype=np.float64)
    x9 = np.concatenate([x8, x8[-1:]])
    nq8 = freeze_e8_q(Q_W=0.0)
    nq9 = replace(nq8, stage=9)
    assert n_unknowns(8) == 45
    assert n_unknowns(9) == 46
    assert f_name(9) == "F_E9"
    f8 = coupled_residual_named(x8, nq8)
    f9 = coupled_residual_named(x9, nq9)
    np.testing.assert_allclose(f9[:-2], f8[:-1], rtol=0.0, atol=1.0e-8)
    assert f9[-2] + f9[-1] == pytest.approx(float(f8[-1]), abs=1.0e-6)
    assert f9[-2] + f9[-1] == pytest.approx(float(energy_row_named(x9, nq9)), abs=1.0e-6)
    x9[ -1] = float(x8[-1]) + 40.0
    f_hot = coupled_residual_named(x9, nq9)
    assert f_hot[-1] != pytest.approx(f9[-1])
    assert abs(float(f_hot[16] - f9[16])) + abs(float(f_hot[21] - f9[21])) > 0.0


def test_bed_surface_bubble_ignores_freeboard_temperature() -> None:
    """床面气泡速度和总抛射能力不随自由空间温度变。床温改变时床面速度可以变。"""
    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_p2 import bed_solids_ejection_split

    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e9.npz")["x"].copy()
    nq = replace(freeze_e8_q(Q_W=0.0), stage=9)
    base = outer_P_named(x, nq)
    cap = bed_solids_ejection_split(base.energy.species.geo_fb, 0.017, 0.206, 1.0)
    hot_fb = x.copy()
    hot_fb[-1] = float(hot_fb[-1]) + 80.0
    fb = outer_P_named(hot_fb, nq)
    cap_fb = bed_solids_ejection_split(fb.energy.species.geo_fb, 0.017, 0.206, 1.0)
    assert fb.hydro.u_b == pytest.approx(base.hydro.u_b, abs=1.0e-12)
    assert cap_fb["total_capacity_kg_s"] == pytest.approx(cap["total_capacity_kg_s"])
    assert fb.reference.d_b == pytest.approx(0.21)
    assert fb.reference.eps_b == pytest.approx(0.13)
    hot_bed = x.copy()
    hot_bed[-2] = float(hot_bed[-2]) + 80.0
    bed = outer_P_named(hot_bed, nq)
    assert bed.hydro.u_b != pytest.approx(base.hydro.u_b, abs=1.0e-4)


def test_freeboard_bundle_passes_bed_top_flow_without_second_ejection() -> None:
    """自由空间入口透传床顶质量流。改变气泡速度不再生成第二份抛射。"""
    from src.core.freeboard.transport import _build_size_class_bundles

    d_p = np.array([2.0e-4, 3.0e-4])
    m_char = np.array([0.0011, 0.0022])
    m_ash = np.zeros(2)
    totals = []
    for factor in (0.5, 1.0, 2.0):
        _d, _c, flows, total = _build_size_class_bundles(
            d_p_classes=d_p,
            m_char_classes=m_char,
            m_ash_classes=m_ash,
            rho_p=1400.0,
            eps_b=0.13,
            eps_d_void=0.45,
            u_b=1.1 * factor,
            d_b=0.21,
            area=0.126,
        )
        np.testing.assert_allclose(flows, m_char)
        totals.append(total)
    assert totals[0] == pytest.approx(totals[1]) == pytest.approx(totals[2]) == pytest.approx(0.0033)
    doubled = m_char * 2.0
    _d, _c, flows, total = _build_size_class_bundles(
        d_p_classes=d_p,
        m_char_classes=doubled,
        m_ash_classes=m_ash,
        rho_p=1400.0,
        eps_b=0.13,
        eps_d_void=0.45,
        u_b=1.1,
        d_b=0.21,
        area=0.126,
    )
    np.testing.assert_allclose(flows, doubled)
    assert total == pytest.approx(0.0066)


def test_bed_face_return_enters_once_and_internal_edge_cancels() -> None:
    """上面几格的向下流不进床。床顶与第一格的内部边流相消。"""
    from src.core.freeboard.transport import bed_face_class_flows_kg_s

    down = np.array([[0.20, 0.05], [0.04, 0.01]], dtype=np.float64)
    face = bed_face_class_flows_kg_s(down)
    np.testing.assert_allclose(face, [0.20, 0.05])
    up = np.array([0.30, 0.10], dtype=np.float64)
    f_bed = -up + face
    f_fb = up - face
    np.testing.assert_allclose(f_bed + f_fb, 0.0, atol=0.0)


def test_hamel_367_limits_are_not_the_sum_of_roof_and_wake() -> None:
    """无聚并只剩 roof，S=1 只剩 wake。各粒级 roof 相加不等于平均粒径的 roof。"""
    from src.physics.wirsum_1997_freeboard import (
        calc_roof_mass_flux_eq359,
        calc_wake_mass_flux_engineering,
        f0_hamel_367,
    )

    common = dict(rho_p=1400.0, eps_d=0.45, eps_b=0.13, u_b=1.1)
    roof = calc_roof_mass_flux_eq359(d_p=2.5e-4, d_b=0.21, **common)
    wake = calc_wake_mass_flux_engineering(zeta_w=0.1, f_w=1.0, **common)
    assert f0_hamel_367(d_bar=2.5e-4, d_b=0.21, zeta_f=0.1, S=0.0, **common) == pytest.approx(roof)
    assert f0_hamel_367(d_bar=2.5e-4, d_b=0.21, zeta_f=0.1, S=1.0, **common) == pytest.approx(wake)
    assert f0_hamel_367(d_bar=2.5e-4, d_b=0.21, zeta_f=0.1, S=0.0, **common) != pytest.approx(roof + wake)
    fine = calc_roof_mass_flux_eq359(d_p=2.0e-4, d_b=0.21, **common)
    coarse = calc_roof_mass_flux_eq359(d_p=3.0e-4, d_b=0.21, **common)
    mean = calc_roof_mass_flux_eq359(d_p=2.5e-4, d_b=0.21, **common)
    assert fine + coarse != pytest.approx(mean)


def test_geometric_bed_gap_is_not_assigned_to_char() -> None:
    """几何床料差额不是炭。炭通量随炭库存变，总抛射能力不跟着变成炭。"""
    from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
    from src.solvers.wirsum_1997_named_f_p2 import bed_solids_ejection_split

    geo = BedFreeboardGeometry()
    base = bed_solids_ejection_split(geo, 0.017, 0.206, 1.0)
    more = bed_solids_ejection_split(geo, 0.034, 0.206, 1.0)
    assert base["gap_is_explicit_species"] is False
    assert base["geometric_gap_kg"] > 50.0
    assert base["inert_kg_s"] == pytest.approx(0.0)
    assert base["char_kg_s"] + base["ash_kg_s"] < base["total_capacity_kg_s"]
    assert more["total_capacity_kg_s"] == pytest.approx(base["total_capacity_kg_s"])
    assert more["char_kg_s"] == pytest.approx(2.0 * base["char_kg_s"])


def test_bed_face_ejection_is_shared_and_independent_of_cell_height() -> None:
    """床面抛出与第一格厚度无关。全回落时顶部出流为零，格内库存仍可为正。"""
    from src.solvers.wirsum_1997_named_f import outer_P_named
    from src.solvers.wirsum_1997_named_f_p2 import bed_face_char_interface

    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e11.npz")["x"]
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    geo = q.energy.species.geo_fb
    m_char = np.array([float(x[14]), float(x[15])], dtype=np.float64)
    d_p = np.array([2.0e-4, 3.0e-4], dtype=np.float64)
    full = bed_face_char_interface(geo, m_char, d_p, H_bed=1.0, H_cell=2.0)
    thin = bed_face_char_interface(geo, m_char, d_p, H_bed=1.0, H_cell=0.5)
    np.testing.assert_allclose(full["m_eject"], thin["m_eject"], atol=1.0e-15)
    assert float(full["m_eject"][1]) > 0.0
    assert float(full["m_cross"][1]) == pytest.approx(0.0, abs=1.0e-12)
    assert float(full["m_return"][1]) == pytest.approx(float(full["m_eject"][1]), rel=1.0e-6)
    assert float(full["m_cell"][1]) > 0.0
    rho = float(geo.rho_p)
    eps_d = float(geo.eps_d)
    eps_b = float(geo.eps_b)
    area = float(geo.area)
    m_geom = rho * (1.0 - eps_d) * (1.0 - eps_b) * area * 1.0
    assert float(full["m_eject"][1]) < 0.01 * m_geom


def test_first_cell_top_is_not_the_furnace_exit() -> None:
    """第一格变薄不改变抛出流。上面还有格子时，越格流量不算出炉。"""
    from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
    from src.solvers.wirsum_1997_named_f_p2 import bed_face_char_interface

    geo = BedFreeboardGeometry(u_b=4.0, u_g=0.05, d_p=2.0e-4)
    m_char = np.array([0.01], dtype=np.float64)
    d_p = np.array([2.0e-4], dtype=np.float64)
    kw = dict(n_bins=5, dt=2.0e-3, t_max=2.0)
    tall = bed_face_char_interface(geo, m_char, d_p, H_bed=1.0, H_cell=2.0, is_roof=True, **kw)
    short = bed_face_char_interface(geo, m_char, d_p, H_bed=1.0, H_cell=0.02, is_roof=False, **kw)
    np.testing.assert_allclose(tall["m_eject"], short["m_eject"], atol=1.0e-15)
    assert float(short["m_cross"][0]) > float(tall["m_cross"][0])
    assert float(short["m_furnace_exit"][0]) == pytest.approx(0.0)
    assert float(short["m_to_next"][0]) == pytest.approx(float(short["m_cross"][0]))
    assert float(short["m_return"][0]) < float(short["m_eject"][0])
    assert float(tall["m_furnace_exit"][0]) == pytest.approx(float(tall["m_cross"][0]))
    assert float(tall["m_to_next"][0]) == pytest.approx(0.0)


def test_constant_ug_return_time_matches_independent_ode() -> None:
    """恒定表观气速下，解析步与独立加速度积分的回落时间相差小于 1%。"""
    from scipy.integrate import solve_ivp

    from src.physics.wirsum_1997_freeboard import (
        integrate_particle_fate_with_time,
        particle_acceleration_cw,
        u_gb_hamel_374,
    )
    from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
    from src.solvers.wirsum_1997_named_f_p2 import integrate_launch_distribution

    kw = dict(u_g=0.16, rho_g=0.4, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.4)
    fate, t_ev, _z, _u = integrate_particle_fate_with_time(
        0.4, cell_height=2.0, dt=2.0e-3, t_max=2.0, **kw,
    )
    assert fate == "returns"

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        return [float(y[1]), particle_acceleration_cw(float(y[1]), **kw)]

    def back(_t: float, y: np.ndarray) -> float:
        return float(y[0])

    back.terminal = True
    back.direction = -1.0
    sol = solve_ivp(
        rhs, (0.0, 2.0), np.array([0.0, 0.4]), events=back,
        rtol=1.0e-6, atol=1.0e-8, dense_output=False,
    )
    assert sol.t_events[0].size == 1
    t_ode = float(sol.t_events[0][0])
    assert abs(t_ev - t_ode) / t_ode < 0.01

    geo = BedFreeboardGeometry(u_g=0.16, u_b=1.13, d_b=0.15, d_p=3.0e-4, rho_g=0.4, rho_p=1400.0, mu_g=4.0e-5)
    coarse = integrate_launch_distribution(geo, n_bins=11, dt=2.0e-3, t_max=2.0, gas_profile="superficial")
    fine = integrate_launch_distribution(geo, n_bins=81, dt=2.0e-3, t_max=2.0, gas_profile="superficial")
    assert coarse["w_return"] == pytest.approx(1.0, abs=1.0e-8)
    assert fine["w_return"] == pytest.approx(1.0, abs=1.0e-8)
    assert abs(coarse["tau_s"] - fine["tau_s"]) / fine["tau_s"] < 0.01
    bubble = integrate_launch_distribution(geo, n_bins=11, dt=2.0e-3, t_max=0.2, gas_profile="hamel_374")
    assert bubble["w_undecided"] > 0.9
    assert bubble["tau_s"] < 0.01 * bubble["occupancy_s"]
    assert u_gb_hamel_374(0.0, 0.16, 1.13, 0.15) == pytest.approx(1.13)
    assert u_gb_hamel_374(20.0, 0.16, 1.13, 0.15) == pytest.approx(0.16, abs=1.0e-3)


def test_upper_cell_return_reaches_the_bed_only_through_the_lower_cell() -> None:
    """两格串联与一整段的回落时间一致。上一格的回落不直接算进床。"""
    from src.physics.wirsum_1997_freeboard import (
        integrate_particle_fate_with_time,
        integrate_particle_through_cells,
    )

    kw = dict(u_g=0.16, rho_g=0.4, rho_p=1400.0, d_p=3.0e-4, mu_g=4.0e-5, z=0.4)
    one, t_one, _z, _u = integrate_particle_fate_with_time(
        0.4, cell_height=2.0, dt=2.0e-3, t_max=2.0, **kw,
    )
    chained, t_chain, cell = integrate_particle_through_cells(
        0.4, [1.0, 1.0], dt=2.0e-3, t_max=2.0, **kw,
    )
    assert one == "returns"
    assert chained == "returns"
    assert cell == 0
    assert abs(t_chain - t_one) / t_one < 0.01
    crossed, t_up, cell_up = integrate_particle_through_cells(
        4.0, [0.02, 0.02], dt=2.0e-3, t_max=2.0, **kw,
    )
    assert crossed == "reaches_next_cell"
    assert cell_up == 1
    assert t_up > 0.0


def test_undecided_weight_is_not_a_bed_or_furnace_exit() -> None:
    """截断未判定的份额不进回床，也不进出炉。它的时间不进库存。"""
    from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
    from src.solvers.wirsum_1997_named_f_p2 import integrate_launch_distribution

    geo = BedFreeboardGeometry(
        u_g=0.16, u_b=1.13, d_b=0.15, d_p=3.0e-4, rho_g=0.4, rho_p=1400.0, mu_g=4.0e-5,
    )
    out = integrate_launch_distribution(
        geo, n_bins=11, dt=2.0e-3, t_max=0.05, gas_profile="superficial", n_cells=2,
    )
    assert out["w_undecided"] > 0.5
    assert out["w_return"] + out["w_leave"] + out["w_undecided"] == pytest.approx(1.0)
    assert out["tau_s"] < 0.05 * out["w_undecided"]


def test_p3_e11_saved_state_passes_block_gates_live() -> None:
    """P3 的 E11 端点按当前残差重算通过分块门，不读文件里的布尔值。"""
    from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
    from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
    from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    blocks = named_endpoint_blocks(coupled_residual_named(x, q, rate_scale=1.0))
    assert blocks["passed"] is True
    assert float(q.hydro.w_return[1]) == pytest.approx(1.0)
    assert float(q.hydro.w_leave[1]) == pytest.approx(0.0)


def test_bed_and_freeboard_temperatures_drive_different_hydro() -> None:
    """床温改变床层气速。只改自由空间温度不改变床层气速，但可以改变轨迹时间。"""
    from src.solvers.wirsum_1997_named_f import outer_P_named

    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e9.npz")["x"].copy()
    nq = replace(freeze_e8_q(Q_W=0.0), stage=9)
    base = outer_P_named(x, nq).hydro
    hot_fb = x.copy()
    hot_fb[-1] = float(hot_fb[-1]) + 80.0
    fb = outer_P_named(hot_fb, nq).hydro
    hot_bed = x.copy()
    hot_bed[-2] = float(hot_bed[-2]) + 80.0
    bed = outer_P_named(hot_bed, nq).hydro
    assert fb.u_g == pytest.approx(base.u_g, rel=0.0, abs=1.0e-12)
    assert bed.u_g != pytest.approx(base.u_g, rel=0.0, abs=1.0e-6)
    assert not np.allclose(fb.tau_s, base.tau_s, rtol=0.0, atol=1.0e-6)


def test_char_cross_enthalpy_is_zero_when_all_return_and_isothermal() -> None:
    """全部回落、无反应、两区同温时，跨床面炭焓流为零。有反应时回落流扣掉同一份消耗。"""
    from src.solvers.wirsum_1997_named_f import (
        char_face_enthalpy_w,
        char_face_flows_kg_s,
        edges_from_x,
        fb_char_balance,
        outer_P_named,
        path_char_consumption_kg_s,
    )

    x10 = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e10.npz")["x"]
    x = np.concatenate([x10[:-2], [0.05, 0.05], x10[-2:]])
    nq = replace(
        freeze_e8_q(Q_W=0.0),
        stage=11,
        enabled_rxns=frozenset({"R5", "R6", "R8", "R12", "R10", "R11"}),
    )
    q = outer_P_named(x, nq)
    T = float(x[-2])
    assert char_face_enthalpy_w(x, q, T, T) == pytest.approx(0.0, abs=1.0e-6)
    nq_r = replace(nq, enabled_rxns=nq.enabled_rxns | frozenset({"R2", "R3", "R4"}))
    q_r = outer_P_named(x, nq_r)
    _eject, m_ret, m_leave = char_face_flows_kg_s(x, q_r)
    _resid, _src, cons = fb_char_balance(x, q_r, rate_scale=1.0)
    on_return, on_leave = path_char_consumption_kg_s(
        cons, q_r.hydro.tau_return_s, q_r.hydro.tau_leave_s,
    )
    ed = edges_from_x(x[:39], q_r.energy.species)
    np.testing.assert_allclose(m_ret, np.asarray(ed["m_return"]) - on_return, atol=1.0e-12)
    np.testing.assert_allclose(m_leave, np.asarray(ed["m_leave"]) - on_leave, atol=1.0e-12)
    assert float(np.sum(on_return + on_leave)) == pytest.approx(float(np.sum(cons)))
    assert float(np.sum(cons)) > 0.0


def test_mixed_fate_consumption_follows_flight_time_once() -> None:
    """混合命运时消耗按飞行时间拆。床、离开和气体里的碳加起来只有一份。"""
    from src.solvers.wirsum_1997_named_f import path_char_consumption_kg_s

    cons = np.array([0.0, 0.04])
    tau_return = np.array([0.0, 0.10])
    tau_leave = np.array([0.0, 0.30])
    w_return = np.array([0.0, 0.80])
    on_return, on_leave = path_char_consumption_kg_s(cons, tau_return, tau_leave)
    assert float(on_return[1]) == pytest.approx(0.01)
    assert float(on_leave[1]) == pytest.approx(0.03)
    assert float(on_return[1]) != pytest.approx(float(cons[1] * w_return[1]))
    m_eject = 0.20
    m_undecided = 0.02
    m_return = w_return[1] * m_eject - float(on_return[1])
    m_leave = 0.10 * m_eject - float(on_leave[1])
    assert m_return + m_leave + m_undecided + float(cons[1]) == pytest.approx(m_eject)


def test_freeboard_char_holdup_uses_lumped_closure_not_path_integral() -> None:
    """库存零点是集总式。同样的 kτ 下，它和轨迹剩余质量积分不是同一个数。"""
    m_eject = 3.0e-4  # [kg/s]
    tau = 0.25  # [s]
    k = 0.20  # [1/s]
    M = m_eject * tau / (1.0 + k * tau)
    consume = k * M
    assert m_eject - M / tau - consume == pytest.approx(0.0, abs=1.0e-15)
    M_path = m_eject / k * (1.0 - np.exp(-k * tau))
    assert M != pytest.approx(M_path, rel=1.0e-3)
    assert abs(M - M_path) / M < 0.05


def test_zone_energies_sum_to_single_row_and_share_the_interface() -> None:
    """两区温度相同时，床层加自由空间等于原来的总能量。交界面焓流等大反号。"""
    from src.solvers.wirsum_1997_named_f import energy_row_named, energy_zones_named

    snap = _snapshot()
    x = np.asarray(snap["e8"]["x"], dtype=np.float64)
    nq = freeze_e8_q(Q_W=0.0)
    T = float(x[-1])
    e_bed, e_fb = energy_zones_named(x, nq, T, T)
    assert e_bed + e_fb == pytest.approx(float(energy_row_named(x, nq)), rel=0.0, abs=1.0e-6)
    e_bed_hi, e_fb_hi = energy_zones_named(x, nq, T + 50.0, T)
    e_bed_lo, e_fb_lo = energy_zones_named(x, nq, T - 50.0, T)
    assert (e_bed_hi + e_fb_hi) - (e_bed_lo + e_fb_lo) == pytest.approx(0.0, abs=1.0e-4)
    assert e_bed_hi != pytest.approx(e_bed)
    assert e_fb_hi == pytest.approx(-e_bed_hi + (e_bed + e_fb), abs=1.0e-4)


def test_trace_freeboard_o2_stays_in_log_ids() -> None:
    """低于 1e-12 mol/s 的自由空间氧仍进入对数坐标。"""
    from src.solvers.wirsum_1997_named_f_e8 import _e8_positive_log_ids

    z = np.zeros(45, dtype=np.float64)
    z[21] = 1.0e-30
    z[-1] = 1400.0
    ids = _e8_positive_log_ids(z)
    assert 21 in set(ids.tolist())
    assert 28 not in set(ids.tolist())


def test_e8_tar_is_two_phase_with_exchange() -> None:
    """Hamel §2.2 / §5.2：焦油两相衡算，交换等大反号；R10 用各相浓度。"""
    from src.solvers.wirsum_1997_named_f import (
        coupled_residual_named,
        r10_major_and_tar_source,
        tar_sources_by_phase,
        unpack_named,
    )
    from src.solvers.wirsum_1997_stiffness_full_gas import IDX

    nq = freeze_e8_q(Q_W=0.0)
    x = default_x0_named(nq)
    x[0:7] = 0.3
    x[7:14] = 0.3
    x[IDX["O2"]] = 0.05
    x[7 + IDX["O2"]] = 0.05
    x[40:44] = 0.0
    x[42] = 1.0e-3
    x[-1] = 1200.0
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q)
    s1d, _s2d, s1b, _s2b = tar_sources_by_phase(x, q, rate_scale=1.0)
    feed1 = float(q.energy.n_tar.get("TAR1", 0.0))
    u = unpack_named(x, 8)
    assert f[40] + f[42] == pytest.approx(
        float(u["N_tar1"]) + float(u["N_tar1_b"]) - feed1 - s1d - s1b
    )
    assert f[42] != pytest.approx(f[40])

    x_d = x.copy()
    x_d[40] = 1.0e-3
    x_d[42] = 0.0
    _sd_b, sb_b, _, _ = r10_major_and_tar_source(x, q, rate_scale=1.0)
    _sd_d, sb_d, _, _ = r10_major_and_tar_source(x_d, q, rate_scale=1.0)
    assert abs(float(sb_b[IDX["O2"]]) - float(sb_d[IDX["O2"]])) > 1.0e-12
    assert n_unknowns(7) == 43
    assert n_unknowns(6) == 43


def test_cold_start_opens_bubble_gas_and_freeboard_tar() -> None:
    """冷初值里整组为零、但进料为正的气泡气体和自由空间焦油会被打开。"""
    from src.solvers.wirsum_1997_named_f import default_x0_named
    from src.solvers.wirsum_1997_named_f_control import reactivate_zero_unknowns

    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = default_x0_named(nq)
    q = outer_P_named(x, nq)
    opened_x, opened = reactivate_zero_unknowns(q, x)
    assert 7 in opened and 44 in opened
    assert float(opened_x[7]) > 0.0
    assert float(opened_x[44]) > 0.0
    assert float(x[7]) == 0.0


def test_wp1_three_cases_match_the_saved_acceptance() -> None:
    """温度 +50 K 过门。减半和冷启动失败时留下分块倍数和调用次数。"""
    import json
    from pathlib import Path

    rows = {
        row["name"]: row
        for row in json.loads(
            Path("docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp1_three_cases.json").read_text(
                encoding="utf-8"
            )
        )["cases"]
    }
    hot = rows["T_plus_50K"]
    assert hot["passed"] is True
    assert hot["reason"] == "block_gates"
    assert hot["blocks"]["gas_l2"] < 1.0e-6
    half = rows["positive_half"]
    assert half["passed"] is False
    assert half["stats"]["F_calls"] > 0.0
    assert half["gate_ratio"]["gas"] > 1.0
    cold = rows["cold_start"]
    assert 7 in cold["opened_zeros"] and 44 in cold["opened_zeros"]
    assert cold["T_fb"] > 250.0
    assert cold["passed"] is False
    assert cold["blocks"]["gas_l2"] > 1.0e-6


def test_stage11_residual_calls_char_balance_once() -> None:
    """一次复合残差只装配一份自由空间炭平衡。结果与基线 F 一致。"""
    import json
    from pathlib import Path

    from src.solvers.wirsum_1997_named_f import coupled_residual_named, fb_char_balance

    pkg = json.loads(Path("data/named_p3_baseline_2026-09-26.json").read_text(encoding="utf-8"))
    state = next(row for row in pkg["states"] if row["name"] == "P3_E11")
    x = np.asarray(state["x"], dtype=np.float64)
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    calls = {"n": 0}
    real = fb_char_balance

    def counted(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    import src.solvers.wirsum_1997_named_f as named

    named.fb_char_balance = counted
    try:
        f = coupled_residual_named(x, q, rate_scale=1.0)
    finally:
        named.fb_char_balance = real
    assert calls["n"] == 1
    np.testing.assert_allclose(f, np.asarray(state["F"], dtype=np.float64), rtol=0.0, atol=1.0e-8)


def test_trace_freeboard_oxygen_column_stays_finite() -> None:
    """自由空间氧的中心差分列保持有限。没有改成解析雅可比。"""
    import json
    from pathlib import Path

    from src.solvers.wirsum_1997_named_f import coupled_residual_named

    pkg = json.loads(Path("data/named_p3_baseline_2026-09-26.json").read_text(encoding="utf-8"))
    state = next(row for row in pkg["states"] if row["name"] == "P3_E11")
    x = np.asarray(state["x"], dtype=np.float64)
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    h = 1.0e-4
    xp = x.copy()
    xm = x.copy()
    xp[21] = x[21] * np.exp(h)
    xm[21] = x[21] * np.exp(-h)
    column = (
        coupled_residual_named(xp, q, rate_scale=1.0)
        - coupled_residual_named(xm, q, rate_scale=1.0)
    ) / (2.0 * h)
    assert np.all(np.isfinite(column))
    assert abs(float(column[21])) > 0.0


def test_block_merit_does_not_prefer_a_large_energy_row() -> None:
    """能量行很大时，原单位范数可以更小，分块尺度不能因此接受它。"""
    from src.solvers.wirsum_1997_named_f_e8 import residual_block_merit

    gas_bad = np.zeros(50)
    gas_bad[0] = 0.02
    energy_bad = np.zeros(50)
    energy_bad[-1] = 0.5
    assert float(np.linalg.norm(gas_bad)) < float(np.linalg.norm(energy_bad))
    assert residual_block_merit(energy_bad) < residual_block_merit(gas_bad)


def test_named_p3_baseline_recomputes_e8_and_p3() -> None:
    """基线包里的两份状态按当前代码重算仍过门，且没有覆盖 E8 快照。"""
    import json
    from pathlib import Path

    from src.solvers.wirsum_1997_named_f import coupled_residual_named
    from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
    from src.solvers.wirsum_1997_named_f_endpoint_io import SNAPSHOT_PATH

    pkg = json.loads(Path("data/named_p3_baseline_2026-09-26.json").read_text(encoding="utf-8"))
    assert pkg["e8_snapshot_overwritten"] is False
    assert pkg["e8_snapshot_path"].endswith(SNAPSHOT_PATH.name)
    for state in pkg["states"]:
        x = np.asarray(state["x"], dtype=np.float64)
        nq = replace(freeze_e8_q(Q_W=0.0), stage=int(state["stage"]))
        f = coupled_residual_named(x, outer_P_named(x, nq), rate_scale=1.0)
        blocks = named_endpoint_blocks(f)
        assert blocks["passed"] is True
        np.testing.assert_allclose(f, np.asarray(state["F"], dtype=np.float64), rtol=0.0, atol=1.0e-9)


def test_p4_readonly_same_frozen_state_explains_deltas_and_does_not_migrate() -> None:
    """同一冻结状态下对照传质、抛射、均相源和床面焓。不切生产入口，不跑全炉。"""
    from types import SimpleNamespace

    from src.core.connectivity.solid_transport import bed_upflow_coefficient
    from src.core.constants import n_b
    from src.kinetics.gas_reactions import rate_R5_bubble
    from src.physics.mass_transfer import calc_kbd, calc_u_br
    from src.solvers.wirsum_1997_named_f import (
        char_face_enthalpy_w,
        fb_gas_source,
        q_at_T,
    )
    from src.solvers.wirsum_1997_stiffness_full_gas import IDX

    x = np.load("docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    b = q.reference
    hydro = q.hydro
    u_d = float(b.u_mf) / max(float(b.eps_mf), 1.0e-12)
    u_br = calc_u_br(u_d, float(b.P), n_b_factor=n_b)
    k_shared = calc_kbd(
        u_br, float(b.d_b), float(hydro.D_g), float(b.eps_mf), float(hydro.u_b),
        model="hamel_350",
    )
    assert float(hydro.K_bd) == pytest.approx(k_shared, rel=0.0, abs=0.0)

    prod_cell = SimpleNamespace(
        eps_b=float(b.eps_b),
        u_b=float(hydro.u_b),
        geo=SimpleNamespace(dh=float(b.H_bed)),
    )
    k_wake = bed_upflow_coefficient(prod_cell, is_top_bed_cell=True)
    k_named = float(q.energy.species.K_eject[1])
    assert k_wake > 0.0 and k_named > 0.0
    assert k_named == pytest.approx(k_wake, rel=0.10)
    assert k_named != pytest.approx(k_wake, rel=0.01)

    q_fb = q_at_T(q.energy.species, float(x[-1]))
    src = fb_gas_source(
        x[:39], q_fb, rate_scale=1.0, enabled=frozenset({"R5"}),
    )
    n_fb = np.maximum(x[16:23], 0.0)
    n_tot = max(float(np.sum(n_fb)), 1.0e-30)
    y = n_fb / n_tot
    c = y * (float(b.P) / (8.314 * float(x[-1])))
    from src.core.species import GAS_SPECIES_INDEX
    v_fb = float(b.area) * float(b.H_freeboard)
    r5 = rate_R5_bubble(
        float(x[-1]), float(c[IDX["CO"]]), float(c[IDX["O2"]]), float(c[IDX["H2O"]]),
        float(b.P), y, GAS_SPECIES_INDEX,
    )
    assert float(src[IDX["CO"]]) == pytest.approx(-r5 * v_fb, rel=1.0e-12, abs=1.0e-18)
    h_cross = char_face_enthalpy_w(x, q, float(x[-2]), float(x[-1]))
    h_same = char_face_enthalpy_w(x, q, float(x[-1]), float(x[-1]))
    assert h_cross != pytest.approx(h_same, abs=1.0)
