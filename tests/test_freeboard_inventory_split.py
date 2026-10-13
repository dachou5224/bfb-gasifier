"""长停留时间下集总库存与路径分摊。不改速率，不推进同伦。"""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
from src.solvers.wirsum_1997_named_f import (
    char_face_flows_kg_s,
    fb_char_balance,
    outer_P_named,
    path_char_consumption_kg_s,
)
from src.solvers.wirsum_1997_named_f_p2 import integrate_launch_distribution
from src.solvers.wesseling_two_zone import q_on_wesseling1_path


def _path_flows(
    m_eject: float,
    k: float,
    w_return: float,
    t_return: float,
    w_leave: float,
    t_leave: float,
) -> tuple[float, float]:
    """平衡库存上的反应后路径流 [kg/s]。时间是单颗粒飞行时间。"""
    tau_return = float(w_return) * float(t_return)
    tau_leave = float(w_leave) * float(t_leave)
    tau = tau_return + tau_leave
    consume = float(m_eject) * k * tau / (1.0 + k * tau)
    on_return, on_leave = path_char_consumption_kg_s(
        np.array([consume]), np.array([tau_return]), np.array([tau_leave]),
    )
    return (
        float(m_eject) * float(w_return) - float(on_return[0]),
        float(m_eject) * float(w_leave) - float(on_leave[0]),
    )


def test_single_path_equilibrium_keeps_the_open_flow_positive() -> None:
    back, _leave = _path_flows(1.0, k=0.45, w_return=1.0, t_return=56.0, w_leave=0.0, t_leave=0.0)
    _back, leave = _path_flows(1.0, k=0.45, w_return=0.0, t_return=0.0, w_leave=1.0, t_leave=7144.0)
    assert back == pytest.approx(1.0 / (1.0 + 0.45 * 56.0))
    assert back > 0.0
    assert leave == pytest.approx(1.0 / (1.0 + 0.45 * 7144.0))
    assert leave > 0.0


def test_equilibrium_time_split_can_still_make_one_path_negative() -> None:
    """回落占时间 10/10.1、质量只占 0.20。平衡消耗仍把回落流打成负的。"""
    back, leave = _path_flows(
        1.0, k=1.0, w_return=0.20, t_return=50.0, w_leave=0.80, t_leave=0.125,
    )
    assert back < 0.0
    assert back + leave == pytest.approx(1.0 / 11.1)


def test_uniform_field_near_the_threshold_does_not_mix_both_paths() -> None:
    """0.2 mm 在这段自由空间里，终速变号时不会同时有两条带时间的路径。"""
    geo = BedFreeboardGeometry(
        u_g=0.455076, u_b=1.352955, d_b=0.21, d_p=2.0e-4,
        rho_g=5.51522, rho_p=1400.0, mu_g=4.0e-5, H_freeboard=4.2688,
    )
    for shift in (0.0, 0.0009, 0.002):
        dist = integrate_launch_distribution(
            replace(geo, u_g=geo.u_g + shift),
            n_cells=1, resolve_events=True, t_max=5.0,
        )
        assert min(dist["w_return"], dist["w_leave"]) < 1.0e-3
        if dist["w_return"] < 1.0e-3:
            assert dist["tau_return_s"] == pytest.approx(0.0, abs=1.0e-12)


def test_certified_coarse_inventory_matches_the_lumped_return() -> None:
    path = Path("docs/audit_reports/s0698_review_2026-10-10/root_s0698.json")
    cert = json.loads(path.read_text())
    x = np.asarray(cert["x"], dtype=float)
    q = outer_P_named(x, q_on_wesseling1_path(float(cert["s"])))
    resid, src, consume = fb_char_balance(x, q, rate_scale=1.0)
    m_eject, m_return, m_leave = char_face_flows_kg_s(
        x, q, rate_scale=1.0, char_pack=(resid, src, consume),
    )
    tau = float(q.hydro.tau_s[1])
    m_coarse = float(x[47])
    k = float(consume[1]) / m_coarse
    m_formula = float(m_eject[1]) * tau / (1.0 + k * tau)
    assert q.hydro.w_leave[1] == pytest.approx(0.0)
    assert q.hydro.w_return[1] == pytest.approx(1.0)
    assert m_coarse == pytest.approx(m_formula, rel=0.002)
    assert float(m_return[1]) > 0.0
    assert float(m_leave[1]) == pytest.approx(0.0)
    assert abs(float(m_return[0])) < 1.0e-12
