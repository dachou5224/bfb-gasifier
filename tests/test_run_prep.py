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
from tools.run_prep import (
    _CLOSED_FORM_WINDOW_S, _book_change, _fuel_packet, applicability_row, closed_form_from_state,
    freeboard_closure, fuel_parts, refine_particle, refine_trajectory, refinement_from_case, toy_particle,
    toy_trajectory,
)

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


def _as_received(frozen, assembly, fractions):
    redirect = assembly.redirect_fuel_by_class_kg_s
    if redirect is None:
        redirect = np.zeros(frozen.n_class)
    return {
        "received": fractions, "dried": dict(fractions),
        "bed_w": np.array(frozen.fuel_up_bed_kg_s[0], dtype=float), "bed_p": np.zeros(frozen.n_class),
        "up_w": np.array(frozen.fuel_up_free_kg_s[0], dtype=float), "up_p": np.zeros(frozen.n_class),
        "down_w": np.array(frozen.fuel_down_free_kg_s[-1], dtype=float), "down_p": np.zeros(frozen.n_class),
        "redirect_w": np.array(redirect, dtype=float), "redirect_p": np.zeros(frozen.n_class),
    }


def test_freeboard_closure_counts_reaction_and_redirect_once():
    fractions = {"C": 0.5, "H": 0.05, "O": 0.3}
    for closure in ("B", "C"):
        # 进入悬停状态的流由格内转入系数给出。ADR 0048 第 2 项。
        c_up = np.zeros((2, 2))
        c_up[1, 0] = 0.4
        frozen = make_frozen(
            hover_closure=closure, c_up_per_s=c_up,
            fuel_up_bed_kg_s=np.full((2, 2), 1.0e-4), fuel_up_free_kg_s=np.full((2, 2), 2.0e-5),
            fuel_down_free_kg_s=np.full((2, 2), 1.0e-5),
        )
        state = make_state(frozen)
        assembly = evaluate_residual(pack_state(state, frozen), frozen)
        reaction = float(np.sum(assembly.char_consumption_free_kg_s) + np.sum(assembly.char_consumption_hover_kg_s))
        assert reaction > 0.0
        book = freeboard_closure(state, frozen, assembly, _as_received(frozen, assembly, fractions))
        assert book["reaction_char_kg_s"] == pytest.approx(reaction)
        assert book["hover_inflow_char_kg_s"] > 0.0
        assert book["redirect_char_kg_s"] + 1.0e-15 >= book["hover_inflow_char_kg_s"]
        carbon_in, carbon_out = _hand_carbon(state, frozen, assembly, fractions)
        assert book["totals"]["C"]["in"] == pytest.approx(carbon_in)
        assert book["totals"]["C"]["out"] == pytest.approx(carbon_out)
        assert book["totals"]["C"]["difference"] == pytest.approx(carbon_in - carbon_out)
        assert book["reaction_char_kg_s"] not in (book["totals"]["C"]["in"], book["totals"]["C"]["out"])
        outgoing = [item for item in book["streams"] if item["name"].startswith(closure)]
        assert len(outgoing) == 1 and outgoing[0]["direction"] == "出" and outgoing[0]["C"] > 0.0
    plain = make_frozen(hover_closure="A")
    quiet_state = make_state(plain)
    quiet_assembly = evaluate_residual(pack_state(quiet_state, plain), plain)
    quiet = freeboard_closure(quiet_state, plain, quiet_assembly, _as_received(plain, quiet_assembly, fractions))
    assert quiet["redirect_char_kg_s"] == 0.0


def test_closed_form_marks_nonpositive_delta_and_compares_the_other():
    config = SolverConfig()
    positive = applicability_row(8.0e-3, 1.1, 1.8e-5, 900.0, 0.8, 1.5, 1.0, 0.4, _CLOSED_FORM_WINDOW_S)
    negative = applicability_row(1.5e-4, 1.2, 2.0e-5, 2000.0, 1.0, 0.5, 1.0, 0.5, _CLOSED_FORM_WINDOW_S)
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


def test_dried_fuel_keeps_its_own_composition():
    received = {"C": 0.48, "H": 0.048, "O": 0.16}
    dried = {"C": 0.6, "H": 0.05, "O": 0.2}
    assert _fuel_packet(0.0, 0.8, received, dried)["C"] == pytest.approx(0.48)
    assert _fuel_packet(0.8, 0.0, received, dried)["C"] == pytest.approx(0.384)
    upward = np.zeros((2, 2))
    upward[0, 0] = 1.0
    frozen = make_frozen(hover_closure="A", a_free_per_s=upward)
    reported = {
        "moisture_mass_fraction": 0.2,
        "undried_free_kg": np.zeros((2, 2)),
        "dried_free_kg": np.array([[1.0, 0.0], [0.0, 0.0]]),
        "undried_bed_kg": np.zeros((2, 2)),
        "dried_bed_kg": np.zeros((2, 2)),
    }
    physical = SimpleNamespace(
        moisture_mass_fraction=0.2, carbon_mass_fraction_wf=0.6,
        hydrogen_mass_fraction_wf=0.05, oxygen_mass_fraction_wf=0.2,
    )
    parts = fuel_parts(frozen, reported, physical)
    assert parts["up_w"][0] == pytest.approx(0.0)
    assert parts["up_p"][0] == pytest.approx(0.8)
    packet = _fuel_packet(parts["up_w"], parts["up_p"], parts["received"], parts["dried"])
    assert packet["C"] == pytest.approx(0.48)


def test_redirect_is_added_to_the_actual_return_or_elutriation():
    from tools.run_prep import stock_report

    c_up = np.zeros((2, 2))
    c_up[1, 0] = 0.4
    shared = dict(
        c_up_per_s=c_up,
        fuel_up_bed_kg_s=np.full((2, 2), 1.0e-4), fuel_up_free_kg_s=np.full((2, 2), 2.0e-5),
        fuel_down_free_kg_s=np.full((2, 2), 1.0e-5),
    )
    for closure, field in (("B", "return_kg_s"), ("C", "elutriation_kg_s")):
        frozen = make_frozen(hover_closure=closure, **shared)
        state = make_state(frozen)
        assembly = evaluate_residual(pack_state(state, frozen), frozen)
        assert assembly.redirect_char_by_class_kg_s is not None
        assert float(np.sum(assembly.redirect_char_by_class_kg_s)) == pytest.approx(assembly.redirect_char_kg_s)
        assert float(np.sum(assembly.redirect_fuel_by_class_kg_s)) == pytest.approx(assembly.redirect_fuel_kg_s)
        stocks = stock_report(state, frozen, assembly, {})
        from tools.run_prep import _mobile_ash

        klass = 0
        upward = float(frozen.a_free_per_s[0, klass])
        downward = float(frozen.b_free_per_s[-1, klass])
        top = float(state.freeboard_char_free_kg[0, klass] + _mobile_ash(state, frozen, 0)[klass])
        bottom = float(state.freeboard_char_free_kg[-1, klass] + _mobile_ash(state, frozen, -1)[klass])
        ordinary_return = downward * bottom + float(frozen.fuel_down_free_kg_s[-1, klass])
        ordinary_out = upward * top + float(frozen.fuel_up_free_kg_s[0, klass])
        extra = float(assembly.redirect_char_by_class_kg_s[klass] + assembly.redirect_ash_by_class_kg_s[klass] + assembly.redirect_fuel_by_class_kg_s[klass])
        assert extra > 0.0
        row = stocks["classes"][klass]
        if closure == "B":
            assert row["return_kg_s"] == pytest.approx(ordinary_return + extra)
            assert row["elutriation_kg_s"] == pytest.approx(ordinary_out)
        else:
            assert row[field] == pytest.approx(ordinary_out + extra)
            assert row["return_kg_s"] == pytest.approx(ordinary_return)


def test_unconverted_fuel_is_part_of_the_cell_solid():
    from tools.run_prep import stock_report

    present = np.zeros((2, 2), dtype=bool)
    present[0, 1] = True
    theta = np.zeros((2, 2))
    theta[0, 1] = 1.0
    quiet = np.zeros((2, 2))
    frozen = make_frozen(
        hover_closure="A", hover=present, free_absent=present, theta=theta,
        a_free_per_s=quiet, a_hover_per_s=quiet, b_free_per_s=quiet, b_hover_per_s=quiet,
    )
    state = make_state(frozen)
    assembly = evaluate_residual(pack_state(state, frozen), frozen)
    undried = np.zeros((2, 2))
    undried[0, 1] = 0.002
    reported = {
        "moisture_mass_fraction": 0.2, "undried_hover_kg": undried, "dried_hover_kg": np.zeros((2, 2)),
        "undried_free_kg": np.zeros((2, 2)), "dried_free_kg": np.zeros((2, 2)),
    }
    stocks = stock_report(state, frozen, assembly, reported)
    pair = next(item for item in stocks["hover"] if item["cell"] == 0 and item["class"] == 2)
    char_ash = float(state.freeboard_char_hover_kg[0, 1] + state.freeboard_ash_hover_kg[0, 1])
    assert pair["hover_kg"] == pytest.approx(char_ash + 0.002)
    assert pair["cell_solid_kg"] >= pair["hover_kg"]
    assert pair["hover_fraction"] == pytest.approx(pair["hover_kg"] / pair["cell_solid_kg"])


def test_refinement_uses_the_supplied_case_and_the_report_prints_it(tmp_path):
    from dataclasses import replace

    config = replace(SolverConfig(), freeboard_ejection_samples=4)
    base_gas = float(toy_trajectory()["gas_m_s"])

    def case(gas, temperature):
        return {
            "trajectory": {**toy_trajectory(), "gas_m_s": gas, "diameter_m": [2.0e-4]},
            "pyrolysis": {**toy_particle(), "temperature_k": temperature, "diameter_m": [2.0e-4]},
        }

    slow = refinement_from_case(config, case(base_gas, 1000.0), None, None, None)
    fast = refinement_from_case(config, case(base_gas + 8.0, 1200.0), None, None, None)
    from hamel.properties.thermo import density_kg_m3

    frozen = make_frozen()
    state = make_state(frozen)
    state.freeboard_temperature_k = np.array([500.0, 1500.0])
    reported = {"phi_s": 0.8, "trajectory": {"bubble_m_s": 1.0, "phi_s": 0.8}}
    short = replace(config, freeboard_max_segments=2)
    formed = closed_form_from_state(state, frozen, reported, short)
    flow = state.freeboard_gas_mol_s[-1]
    total = float(np.sum(flow))
    fractions = {name: float(flow[index]) / total for index, name in enumerate(GAS_SPECIES)}
    bottom = density_kg_m3(1500.0, float(frozen.freeboard_pressure_pa[-1]), fractions)
    top = density_kg_m3(500.0, float(frozen.freeboard_pressure_pa[0]), fractions)
    assert formed["classes"][0]["rho_g"] == pytest.approx(bottom)
    assert formed["classes"][0]["rho_g"] != pytest.approx(top)
    assert slow["trajectory"]["rows"][0]["a_f"] != fast["trajectory"]["rows"][0]["a_f"]
    assert slow["pyrolysis"]["rows"][0]["tau_pyr_s"] != fast["pyrolysis"]["rows"][0]["tau_pyr_s"]
    record = {
        "name": "基准", "change": "不改", "source": "ADR 0030", "status": "CONVERGED",
        "balance": {"ok": True, "ratios": {}}, "nonnegative": True, "ratios": {}, "reason": "构造",
        "measurements": {"T_FRB_K": 1100.0, "X_C": 50.0, "CO": 0.1, "CO2": 0.1, "H2": 0.1, "CH4": 0.01, "H2O": 0.1},
        "diagnostics": {
            "closed_form": {
                "benchmark": {"delta": 1.0, "step_s": 0.02, "samples": 4, "max_relative": 1.25e-9, "velocity_relative": 1.25e-9, "height_relative": 1.25e-9},
                "classes": [],
            },
            "refinement": {
                **fast,
                "pyrolysis": {**fast["pyrolysis"], "met": False},
                "particle": {"temperature_change_k": 1.2, "residual_relative_norm": 3.0, "met": False},
            },
        },
    }
    write_report(tmp_path, [record], title="真实报告", toy=False)
    rendered = (tmp_path / "driver_report.md").read_text(encoding="utf-8")
    assert "1.25e-09" in rendered or "1.25e-9" in rendered
    for name in ("radial", "time", "activation"):
        assert name in rendered
    assert "未满足" in rendered
    assert "构造状态" not in rendered


def test_missing_free_state_stays_in_the_refinement():
    from dataclasses import replace

    from tests.test_freeboard_accounting import _DENSITY, _DIAMETER, _GAS, _MU, _PHI, _column, _cellwise

    edges, rho, mu = _column(3, 1.5, _GAS, _MU)
    column = {
        "edges_m": edges, "gas": _cellwise(1.0, 1.0, 0.2), "rho_g": rho, "mu_g": mu,
        "diameter_m": _DIAMETER, "rho_p": _DENSITY, "phi_s": _PHI, "bubble_m_s": 0.05,
    }
    config = replace(SolverConfig(), freeboard_ejection_samples=4)
    result = refine_trajectory(config, column)
    # 半格捕获区下最上一格的自由状态不存在（ADR 0048 第 2 项），其余格有自由状态。
    assert result["a_f"][-1] is None
    assert result["existence_changed"] is False
    absent = SimpleNamespace(a_f=None, b_f=None, c_up=None, c_down=None)
    present = SimpleNamespace(a_f=0.2, b_f=None, c_up=None, c_down=None)
    same, unchanged = _book_change(SimpleNamespace(cells=[absent]), SimpleNamespace(cells=[absent]))
    gap, changed = _book_change(SimpleNamespace(cells=[absent]), SimpleNamespace(cells=[present]))
    assert same == 0.0 and unchanged is False
    assert changed is True and gap == 0.0


def _delete(node, name) -> None:
    if isinstance(node, dict):
        node.pop(name, None)
        for item in list(node.values()):
            _delete(item, name)
    elif isinstance(node, list):
        for item in node:
            _delete(item, name)
