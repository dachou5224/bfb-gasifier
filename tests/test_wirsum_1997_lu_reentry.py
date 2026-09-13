"""P5：LU 重诊清单。不用实验 CO 做 numerical bridge 判据。"""

from __future__ import annotations

import json

from src.solvers.wirsum_1997_lu_reentry import (
    CHECKLIST_ORDER,
    FORBIDDEN_BRIDGE_KEYS,
    audit_composition_observation,
    audit_lu_reentry,
    build_p5_lu_reentry_payload,
)


def test_checklist_order_is_fixed() -> None:
    assert CHECKLIST_ORDER == (
        "elemental_closure",
        "mass_closure",
        "energy_closure",
        "inner_nr",
        "outer_convergence",
        "freeboard_closure",
        "composition_observation",
        "temperature_observation",
    )


def test_composition_observation_has_no_experiment_bridge_keys() -> None:
    class _Cell:
        T = 1100.0

        def _mole_fractions(self, _kind: str):
            return [0.1] * 20

    payload = audit_composition_observation([_Cell()], p0_exit=None)
    assert payload["forbidden"] == list(FORBIDDEN_BRIDGE_KEYS)
    assert payload["role"] == "observation_vs_p0_baseline"
    for key in FORBIDDEN_BRIDGE_KEYS:
        assert key not in payload["exit_gas_dry"]
        assert payload.get(key) is None


def test_reduced_3cell_x0_follows_checklist_and_opens_closures() -> None:
    payload = build_p5_lu_reentry_payload(include_jacobian=False)
    assert payload["solved"] is False
    assert payload["checklist_order"] == list(CHECKLIST_ORDER)
    assert tuple(payload["items"]) == CHECKLIST_ORDER

    elem = payload["items"]["elemental_closure"]
    assert elem["open"] is True
    assert elem["dominant_element"] in {"H", "O"}
    carbon = [row["rates_mol_atom_s"]["C"] for row in elem["per_cell"]]
    assert max(abs(v) for v in carbon) < 1e-9

    mass = payload["items"]["mass_closure"]
    energy = payload["items"]["energy_closure"]
    assert mass["open"] is True
    assert energy["open"] is True
    assert mass["max_gas_residual_l2"] > 1.0

    inner = payload["items"]["inner_nr"]
    assert inner["solved"] is False
    assert inner["residual_l2"] > 0.0

    outer = payload["items"]["outer_convergence"]
    assert outer["status"] == "single_frame_only"
    assert outer["contractivity"] is None

    fb = payload["items"]["freeboard_closure"]
    assert fb["resolved_name"] == "engineering_freeboard"
    assert fb["wirsum_1997_exact_wired"] is False

    comp = payload["items"]["composition_observation"]
    temp = payload["items"]["temperature_observation"]
    assert comp["role"] == "observation_vs_p0_baseline"
    assert temp["role"] == "observation_vs_p0_baseline"
    for key in FORBIDDEN_BRIDGE_KEYS:
        assert key not in comp
        assert key not in temp
    assert "elemental_closure" in payload["open_items"]


def test_audit_does_not_call_solve(monkeypatch) -> None:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)

    def _boom(*_args, **_kwargs):
        raise AssertionError("P5 must not call reactor.solve")

    monkeypatch.setattr(reactor, "solve", _boom)
    payload = audit_lu_reentry(reactor, include_jacobian=False)
    assert payload["solved"] is False
