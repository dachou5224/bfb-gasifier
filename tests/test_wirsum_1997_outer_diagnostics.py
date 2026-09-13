"""P4：Vorabrechnung outer map 场增量与收缩比。不依赖 Reactor / LU。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from src.solvers.wirsum_1997_outer_diagnostics import (
    HYDRO_FIELDS,
    OuterIterationRecorder,
    apply_linear_outer_map,
    build_outer_map_ab_payload,
    cell_solid_transport_K_s,
    contractivity_ratio,
    field_l2_deltas,
    make_vorab_fields,
    snapshot_vorab_fields,
)


def _star() -> object:
    return make_vorab_fields(
        eps_b=[0.20, 0.18, 0.15],
        K_bd=[4.0, 3.5, 3.0],
        K_s=[0.40, 0.30, 0.20],
        u_b=[0.80, 0.85, 0.90],
        T=[1100.0, 1080.0, 1050.0],
    )


def _offset() -> object:
    return make_vorab_fields(
        eps_b=[0.35, 0.30, 0.28],
        K_bd=[1.0, 1.2, 1.5],
        K_s=[0.10, 0.10, 0.10],
        u_b=[0.40, 0.45, 0.50],
        T=[900.0, 920.0, 940.0],
    )


def test_field_l2_deltas_match_hand_calc() -> None:
    a = make_vorab_fields(
        eps_b=[0.2, 0.3],
        K_bd=[1.0, 2.0],
        K_s=[0.1, 0.2],
        u_b=[0.5, 0.5],
        T=[1000.0, 1100.0],
    )
    b = make_vorab_fields(
        eps_b=[0.2, 0.4],
        K_bd=[1.0, 2.0],
        K_s=[0.1, 0.2],
        u_b=[0.5, 0.5],
        T=[1000.0, 1100.0],
    )
    delta = field_l2_deltas(a, b)
    assert delta["eps_b"] == pytest.approx(0.1)
    assert delta["K_bd"] == pytest.approx(0.0)
    assert delta["T"] == pytest.approx(0.0)


def test_contractivity_ratio_undefined_when_prev_is_zero() -> None:
    assert contractivity_ratio(0.2, 0.0) is None
    assert contractivity_ratio(0.2, 0.4) == pytest.approx(0.5)


def test_linear_map_contractive_ratios_are_exactly_rho() -> None:
    h_star = _star()
    h0 = _offset()
    h1 = apply_linear_outer_map(h0, h_star, 0.4)
    h2 = apply_linear_outer_map(h1, h_star, 0.4)
    d1 = field_l2_deltas(h0, h1)
    d2 = field_l2_deltas(h1, h2)
    for name in HYDRO_FIELDS:
        assert contractivity_ratio(d2[name], d1[name]) == pytest.approx(0.4)


def test_expanding_map_is_flagged_even_if_inner_residual_drops() -> None:
    payload = build_outer_map_ab_payload()
    shrink = payload["contractive_rho_0p4"]
    grow = payload["expanding_rho_1p3"]
    assert shrink["verdict"]["non_contractive"] is False
    assert grow["verdict"]["non_contractive"] is True
    assert shrink["verdict"]["mean_hydro_rho_tail"] == pytest.approx(0.4)
    assert grow["verdict"]["mean_hydro_rho_tail"] == pytest.approx(1.3)

    r_shrink = [step["inner_residual"] for step in shrink["steps"]]
    r_grow = [step["inner_residual"] for step in grow["steps"]]
    assert r_shrink[-1] < r_shrink[0]
    assert r_grow[-1] < r_grow[0]
    assert r_shrink[-1] == pytest.approx(r_grow[-1])


def test_recorder_writes_history_json(tmp_path) -> None:
    rec = OuterIterationRecorder()
    h_star = _star()
    h = _offset()
    rec.record(h, 1.0)
    rec.record(apply_linear_outer_map(h, h_star, 0.5), 0.5)
    path = rec.write_json(tmp_path / "outer_iteration_history.json", label="unit")
    text = path.read_text(encoding="utf-8")
    assert "eps_b" in text
    assert "K_bd" in text
    assert "inner_residual" in text
    payload = rec.payload(label="unit")
    assert payload["steps"][0]["delta"] is None
    assert payload["steps"][1]["delta"]["eps_b"] > 0.0


def test_snapshot_from_cells_reads_solid_K() -> None:
    @dataclass
    class _Cell:
        T: float
        eps_b: float
        K_bd: float
        u_b: float
        K_solid_auf: np.ndarray
        K_solid_ab: np.ndarray

    cells = [
        _Cell(
            T=1000.0,
            eps_b=0.2,
            K_bd=3.0,
            u_b=0.8,
            K_solid_auf=np.array([[0.10, 0.10]], dtype=np.float64),
            K_solid_ab=np.array([[0.20, 0.20]], dtype=np.float64),
        )
    ]
    assert cell_solid_transport_K_s(cells[0]) == pytest.approx(0.30)
    snap = snapshot_vorab_fields(cells)
    assert snap.K_s[0] == pytest.approx(0.30)
    assert snap.T[0] == pytest.approx(1000.0)
