"""T26a：口径、逐格 D_g、缺字段、闭合记账和报告标记。"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from hamel.cases.config import SolverConfig
from hamel.hydrodynamics.bed import BedConditions, BedProfile, local_bed_state
from hamel.properties.data import GAS_SPECIES, atomic_weight, atoms
from hamel.properties.transport import interphase_diffusivity_m2_s
from hamel.solver.driver import _measurements
from hamel.solver.precalc import SEARCH, PrecalcError, _exchange_diffusivity, _load
from hamel.solver.residual import evaluate_residual
from hamel.solver.vector import pack_state
from tests.test_assembly import make_frozen, make_state
from tools.driver_report import write_report
from tools.run_prep import applicability_row, freeboard_closure, refine_particle

_CASE = Path(__file__).resolve().parents[1] / "data" / "cases" / "vtt_pressurized.json"


def _flow():
    flow = np.zeros((1, len(GAS_SPECIES)))
    index = {name: position for position, name in enumerate(GAS_SPECIES)}
    flow[0, index["CO"]] = 2.0
    flow[0, index["C6H6"]] = 0.5
    flow[0, index["N2"]] = 8.0
    return flow


def _physical():
    return SimpleNamespace(
        fuel_mass_flow=1.0, moisture_mass_fraction=0.0,
        carbon_mass_fraction_wf=0.6, hydrogen_mass_fraction_wf=0.04, oxygen_mass_fraction_wf=0.2,
    )


def _state(temperatures):
    return SimpleNamespace(freeboard_gas_mol_s=_flow(), freeboard_temperature_k=np.asarray(temperatures, dtype=float))


def test_freeboard_temperature_is_height_weighted_and_x_c_ignores_solids():
    temperatures = np.array([1000.0, 1100.0, 1200.0])
    equal = _measurements(_state(temperatures), _physical(), np.full(3, 0.2875))
    assert equal["T_FRB_K"] == float(np.mean(temperatures))
    heights = np.array([0.1, 0.2, 0.4])
    weighted = _measurements(_state(temperatures), _physical(), heights)
    assert weighted["T_FRB_K"] == float(np.dot(temperatures, heights) / np.sum(heights))
    carbon = 2.0 + 0.5 * 6.0
    fuel = 0.6 / (atomic_weight("C") * 1.0e-3)
    assert equal["X_C"] == 100.0 * carbon / fuel
    assert weighted["X_C"] == equal["X_C"]


def test_bed_exchange_uses_each_cell_temperature():
    pressure = 101300.0
    temperatures = np.array([900.0, 1300.0])
    diffusivity = interphase_diffusivity_m2_s(float(np.mean(temperatures)), pressure)
    conditions = BedConditions(1.0e-3, 0.46, 1600.0, 1.2, 2.8e-5, 2.8e-5 / 1.2, pressure, diffusivity, 0.15, "technical", "B", True)
    profile = BedProfile(
        0.2, np.full(2, 0.2), np.full(2, 1.2), np.full(2, 0.04), np.full(2, 0.7),
        np.full(2, 1.0), np.full(2, 0.4), np.full(2, 0.01),
    )
    flow = SimpleNamespace(bed_temperature=temperatures, pressure=pressure, u0=np.full(2, 0.8))
    updated = _exchange_diffusivity(SimpleNamespace(solver=SolverConfig()), flow, profile, conditions)
    assert updated.k_bd[0] != updated.k_bd[1]
    for index, temperature in enumerate(temperatures):
        local = local_bed_state(
            0.8, 0.04, 0.2, conditions.__class__(**{**conditions.__dict__, "diffusivity": interphase_diffusivity_m2_s(float(temperature), pressure)}),
            SEARCH,
        )
        assert updated.k_bd[index] == local.k_bd
    assert np.array_equal(updated.d_b, profile.d_b)


def test_missing_voidage_or_sphericity_raises(tmp_path):
    doc = json.loads(_CASE.read_text(encoding="utf-8"))
    _delete(doc, "bed_sphericity")
    path = tmp_path / "missing.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(PrecalcError, match="bed_sphericity"):
        _load(path, "1")
    doc = json.loads(_CASE.read_text(encoding="utf-8"))
    _delete(doc, "bed_minimum_fluidization_voidage")
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(PrecalcError, match="bed_minimum_fluidization_voidage"):
        _load(path, "1")


def test_freeboard_closure_counts_reaction_and_redirect_once():
    fractions = {"C": 0.5, "H": 0.05, "O": 0.3}
    for closure in ("B", "C"):
        hover = np.zeros((2, 2))
        hover[1, 0] = 0.4
        frozen = make_frozen(
            hover_closure=closure, a_hover_per_s=hover,
            ejection_to_free=np.array([0.75, 0.75]), ejection_to_hover=np.array([0.25, 0.25]),
            fuel_up_bed_kg_s=np.full((2, 2), 1.0e-4), fuel_up_free_kg_s=np.full((2, 2), 2.0e-5),
            fuel_down_free_kg_s=np.full((2, 2), 1.0e-5),
        )
        state = make_state(frozen)
        assembly = evaluate_residual(pack_state(state, frozen), frozen)
        reaction = float(np.sum(assembly.char_consumption_free_kg_s) + np.sum(assembly.char_consumption_hover_kg_s))
        assert reaction > 0.0
        book = freeboard_closure(state, frozen, assembly, fractions)
        assert book["reaction_char_kg_s"] == pytest.approx(reaction)
        assert book["ejection_hover_char_kg_s"] > 0.0
        assert book["redirect_char_kg_s"] + 1.0e-15 >= book["ejection_hover_char_kg_s"]
        carbon_in, carbon_out = _hand_carbon(state, frozen, assembly, fractions)
        assert book["totals"]["C"]["in"] == pytest.approx(carbon_in)
        assert book["totals"]["C"]["out"] == pytest.approx(carbon_out)
        assert book["totals"]["C"]["difference"] == pytest.approx(carbon_in - carbon_out)
        assert book["reaction_char_kg_s"] not in (book["totals"]["C"]["in"], book["totals"]["C"]["out"])
        outgoing = [item for item in book["streams"] if item["name"].startswith(closure)]
        assert len(outgoing) == 1 and outgoing[0]["direction"] == "出" and outgoing[0]["C"] > 0.0
    plain = make_frozen(hover_closure="A")
    quiet = freeboard_closure(make_state(plain), plain, evaluate_residual(pack_state(make_state(plain), plain), plain), fractions)
    assert quiet["redirect_char_kg_s"] == 0.0


def test_closed_form_marks_nonpositive_delta_and_compares_the_other():
    config = SolverConfig()
    positive = applicability_row(8.0e-3, 1.1, 1.8e-5, 900.0, 0.8, 1.5, 1.0, 0.4, config.freeboard_max_step_s, 0.05)
    negative = applicability_row(1.5e-4, 1.2, 2.0e-5, 2000.0, 1.0, 0.5, 1.0, 0.5, config.freeboard_max_step_s, 0.05)
    assert positive["delta"] > 0.0 and positive["status"] == "compared"
    assert positive["max_height_relative"] is not None
    assert negative["delta"] <= 0.0 and negative["status"] == "N/A"
    assert negative["max_velocity_relative"] is None


def test_half_tolerance_is_marked_and_error_keeps_stderr(tmp_path, monkeypatch):
    baseline = {
        "name": "baseline", "change": "不改", "source": "ADR 0030", "status": "CONVERGED",
        "measurements": {"T_FRB_K": 1100.0, "X_C": 50.0, "CO": 0.1, "CO2": 0.1, "H2": 0.1, "CH4": 0.01, "H2O": 0.1},
        "balance": {"ok": True, "ratios": {}}, "nonnegative": True, "ratios": {}, "diagnostics": {},
    }
    shifted = dict(baseline)
    shifted["name"] = "S1"
    shifted["measurements"] = dict(baseline["measurements"])
    shifted["measurements"]["T_FRB_K"] = 1115.0
    shifted["measurements"]["CO"] = 0.1
    payload = write_report(tmp_path, [baseline, shifted], title="标记", toy=False)
    row = payload["differences"][0]
    assert "T_FRB_K" in row["over_half_tolerance"]
    assert "CO" not in row["over_half_tolerance"]
    from hamel.solver.campaign import run_campaign

    def fail(_args, _timeout):
        return SimpleNamespace(returncode=2, stdout="", stderr="完整失败\n第二行\n")

    monkeypatch.setattr("hamel.solver.campaign._process", fail)
    records = run_campaign(_CASE, tmp_path / "run", timeout_s=1.0)
    assert records[0]["status"] == "ERROR"
    log = Path(records[0]["stderr_log"])
    assert log.read_text(encoding="utf-8") == "完整失败\n第二行\n"
    text = (tmp_path / "run" / "driver_report.md").read_text(encoding="utf-8")
    assert "ERROR" in text and str(log) in text


def test_committed_toy_report_lists_the_run_prep_items():
    text = (Path(__file__).resolve().parents[1] / "reports" / "toy" / "driver_report.md").read_text(encoding="utf-8")
    for item in (
        "各验收量共享壁面散热这一校准输入", "F26", "F28", "全自由空间闭合",
        "气固反应消耗的焦炭", "数值积分与闭式", "超过容差一半", "悬停占比来自闭合 A",
    ):
        assert item in text


def test_particle_refinement_reports_the_temperature_gate():
    frozen = make_frozen()
    state = make_state(frozen)
    report = refine_particle(state, frozen, 1.0, 1.0, 1.0)
    assert report["temperature_change_k"] >= 0.0
    assert "residual_relative_norm" in report
    assert report["temperature_gate_k"] == 0.1
    assert report["residual_gate"] == 1.0e-6


def _hand_carbon(state, frozen, assembly, fractions) -> tuple[float, float]:
    scale = atomic_weight("C") * 1.0e-3

    def gas(flow) -> float:
        total = 0.0
        for index, name in enumerate(GAS_SPECIES):
            total += float(flow[index]) * atoms(name).get("C", 0) * scale
        return total

    up = frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]
    down = frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]
    incoming = (
        gas(state.bed_gas_mol_s[0, 0] + state.bed_gas_mol_s[0, 1])
        + gas(frozen.secondary_feed_mol_s)
        + float(np.sum(frozen.k_auf_per_s[0] * state.bed_char_kg[0]))
        + float(np.sum(frozen.fuel_up_bed_kg_s[0])) * fractions["C"]
    )
    outgoing = (
        gas(state.freeboard_gas_mol_s[0])
        + float(np.sum(up * state.freeboard_char_free_kg[0]))
        + float(np.sum(frozen.fuel_up_free_kg_s[0])) * fractions["C"]
        + float(np.sum(down * state.freeboard_char_free_kg[-1]))
        + float(np.sum(frozen.fuel_down_free_kg_s[-1])) * fractions["C"]
        + float(assembly.redirect_char_kg_s)
        + float(assembly.redirect_fuel_kg_s) * fractions["C"]
    )
    return incoming, outgoing


def _delete(node, name) -> None:
    if isinstance(node, dict):
        node.pop(name, None)
        for item in list(node.values()):
            _delete(item, name)
    elif isinstance(node, list):
        for item in node:
            _delete(item, name)
