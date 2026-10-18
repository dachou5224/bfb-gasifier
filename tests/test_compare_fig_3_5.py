# -*- coding: utf-8 -*-
"""ADR 0013 的判定。图上的点是合成的，数字化文件尚未进入 main。"""
from __future__ import annotations

import dataclasses
import json

import pytest

from hamel.hydrodynamics.bed import (
    BedConditions,
    minimum_fluidization_velocity,
    state_at_height,
    throughflow_velocity,
    suspension_superficial_velocity,
)
from tests.test_hydrodynamics import _search
from tools.compare_fig_3_5 import (
    bubble_diameter,
    compare,
    decide,
    main,
    normalize_diameter_figure,
    normalize_throughflow_figure,
    throughflow,
)


def _row(name, error, prime, solvable=True):
    return {"name": name, "E": error, "E_prime": prime, "solvable": solvable}


def test_unique_primary_candidate_is_selected():
    rows = [
        _row("adopted", 0.4, 2.0),
        _row("other", 2.2, 0.1),
        _row("far", 4.0, 0.1),
    ]
    result = decide(rows)
    assert result["status"] == "selected"
    assert result["selected"] == "adopted"


def test_secondary_figure_breaks_a_primary_tie():
    rows = [
        _row("left", 1.0, 0.2),
        _row("right", 1.5, 2.5),
    ]
    result = decide(rows)
    assert result["selected"] == "left"


def test_tied_secondary_stays_unresolved():
    rows = [_row("left", 1.0, 0.4), _row("right", 1.2, 0.8)]
    result = decide(rows)
    assert result["status"] == "unresolved"
    assert result["selected"] is None
    assert set(result["primary"]) == {"left", "right"}


def test_no_qualified_combination_stays_unresolved():
    result = decide([_row("a", 3.1, 0.1), _row("b", 4.0, 0.1)])
    assert result["status"] == "unresolved"
    assert result["selected"] is None


def test_unsolvable_combination_blocks_selection():
    result = decide([_row("good", 0.1, 0.1), _row("bad", None, None, solvable=False)])
    assert result["status"] == "unresolved"
    assert result["selected"] is None


def test_boundary_of_the_band_is_included():
    result = decide([_row("best", 1.0, 0.4), _row("edge", 2.0, 1.0), _row("out", 2.01, 0.1)])
    assert result["status"] == "unresolved"
    assert set(result["primary"]) == {"best", "edge"}


def test_synthetic_figure_selects_the_adopted_reading():
    model = _synthetic_model()
    diameter_points = []
    through_points = []
    for excess in (0.25, 0.4):
        for pressure in (1.0e6, 1.5e6):
            rho_g, mu_g, nu_g = model["properties_at"](pressure)
            u_mf = minimum_fluidization_velocity(
                model["epsilon_mf"], nu_g, model["phi_s"], model["d_p"], model["rho_p"], rho_g
            )
            conditions = BedConditions(
                model["d_p"],
                model["epsilon_mf"],
                model["rho_p"],
                rho_g,
                mu_g,
                nu_g,
                pressure,
                model["diffusivity"],
                model["vessel_diameter"],
                "technical",
                "B",
                True,
            )
            state = state_at_height(
                model["height_m"],
                u_mf + excess,
                u_mf,
                model["area_m2"],
                model["n_0"],
                model["h_a"],
                model["n_steps"],
                conditions,
                model["search"],
            )
            tool_diameter = bubble_diameter(
                model["height_m"], excess, pressure, "interstitial", "power", model
            )
            assert tool_diameter == pytest.approx(state.d_b, rel=1e-9)
            diameter_points.append(
                {
                    "excess_velocity_m_s": excess,
                    "pressure_Pa": pressure,
                    "bubble_diameter_m": state.d_b,
                }
            )
            suspension = suspension_superficial_velocity(u_mf, u_mf + excess, "technical", "B", True)
            through_points.append(
                {
                    "excess_velocity_m_s": excess,
                    "pressure_Pa": pressure,
                    "throughflow_velocity_m_s": throughflow_velocity(suspension, pressure),
                }
            )
            assert throughflow(excess, pressure, model) == pytest.approx(through_points[-1]["throughflow_velocity_m_s"])
    decision = compare(
        {
            "height_m": model["height_m"],
            "est_error": {"bubble_diameter_m": 1.0e-4},
            "points": diameter_points,
        },
        {"est_error": {"throughflow_velocity_m_s": 1.0e-3}, "points": through_points},
        model,
    )
    adopted = next(row for row in decision["rows"] if row["name"] == "interstitial/power")
    others = [row for row in decision["rows"] if row["name"] != "interstitial/power"]
    assert adopted["E"] == pytest.approx(0.0, abs=1e-6)
    assert all(row["E"] > adopted["E"] for row in others)
    assert decision["status"] == "selected"
    assert decision["selected"] == "interstitial/power"


def test_digitized_curves_drop_unseparated_points():
    from tools.compare_fig_3_5 import load_figure

    diameter = normalize_diameter_figure(load_figure("data/digitized/fig_3_5.json"))
    through = normalize_throughflow_figure(load_figure("data/digitized/fig_3_3.json"))
    assert diameter["height_m"] == pytest.approx(0.55)
    assert len(diameter["points"]) == 30
    assert diameter["est_error"]["bubble_diameter_m"] == pytest.approx(0.005)
    assert len(through["points"]) == 22
    assert through["est_error"]["throughflow_velocity_m_s"] == pytest.approx(0.02)


def test_dict_search_is_rejected_before_the_six_readings():
    model = _synthetic_model()
    model["search"] = dataclasses.asdict(model["search"])
    with pytest.raises(ValueError, match="HydroSearch"):
        compare(
            {"height_m": 0.55, "est_error": {"bubble_diameter_m": 1.0e-4}, "points": []},
            {"est_error": {"throughflow_velocity_m_s": 1.0e-3}, "points": []},
            model,
        )


def test_model_json_binds_search_and_does_not_mark_every_reading_unsolvable(tmp_path, capsys):
    model = _synthetic_model()
    pressure = 1.0e6
    excess = 0.25
    rho_g, mu_g, nu_g = model["properties_at"](pressure)
    u_mf = minimum_fluidization_velocity(
        model["epsilon_mf"], nu_g, model["phi_s"], model["d_p"], model["rho_p"], rho_g
    )
    conditions = BedConditions(
        model["d_p"],
        model["epsilon_mf"],
        model["rho_p"],
        rho_g,
        mu_g,
        nu_g,
        pressure,
        model["diffusivity"],
        model["vessel_diameter"],
        "technical",
        "B",
        True,
    )
    state = state_at_height(
        model["height_m"],
        u_mf + excess,
        u_mf,
        model["area_m2"],
        model["n_0"],
        model["h_a"],
        model["n_steps"],
        conditions,
        model["search"],
    )
    suspension = suspension_superficial_velocity(u_mf, u_mf + excess, "technical", "B", True)
    diameter_path = tmp_path / "fig_3_5.json"
    through_path = tmp_path / "fig_3_3.json"
    model_path = tmp_path / "model.json"
    diameter_path.write_text(
        json.dumps(
            {
                "height_m": model["height_m"],
                "est_error": {"bubble_diameter_m": 1.0e-4},
                "points": [
                    {
                        "excess_velocity_m_s": excess,
                        "pressure_Pa": pressure,
                        "bubble_diameter_m": state.d_b,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    through_path.write_text(
        json.dumps(
            {
                "est_error": {"throughflow_velocity_m_s": 1.0e-3},
                "points": [
                    {
                        "excess_velocity_m_s": excess,
                        "pressure_Pa": pressure,
                        "throughflow_velocity_m_s": throughflow_velocity(suspension, pressure),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    payload = {
        name: model[name]
        for name in (
            "d_p",
            "rho_p",
            "epsilon_mf",
            "phi_s",
            "diffusivity",
            "vessel_diameter",
            "area_m2",
            "n_0",
            "h_a",
            "height_m",
            "n_steps",
        )
    }
    payload["search"] = dataclasses.asdict(model["search"])
    payload["pressures"] = [
        {"pressure_Pa": pressure, "rho_g": rho_g, "mu_g": mu_g, "nu_g": nu_g}
    ]
    model_path.write_text(json.dumps(payload), encoding="utf-8")
    assert main(
        ["--fig-3-5", str(diameter_path), "--fig-3-3", str(through_path), "--model", str(model_path)]
    ) == 0
    reported = json.loads(capsys.readouterr().out)
    assert reported["selected"] == "interstitial/power"
    assert all(row["solvable"] for row in reported["rows"])


def test_incomplete_search_json_stops_before_comparing(tmp_path):
    model = _synthetic_model()
    payload = {"pressures": [], "search": dataclasses.asdict(model["search"])}
    del payload["search"]["velocity_floor"]
    path = tmp_path / "model.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SystemExit, match="尚未比较"):
        main(["--fig-3-5", "unused.json", "--fig-3-3", "unused.json", "--model", str(path)])


def _synthetic_model():
    pressures = {
        1.0e6: (1.2 * 1.0e6 / 101_300.0, 2.8e-5),
        1.5e6: (1.2 * 1.5e6 / 101_300.0, 2.8e-5),
    }

    def properties_at(pressure):
        rho_g, mu_g = pressures[float(pressure)]
        return rho_g, mu_g, mu_g / rho_g

    return {
        "d_p": 0.7e-3,
        "rho_p": 2600.0,
        "epsilon_mf": 0.46,
        "phi_s": 0.8,
        "diffusivity": 2.0e-5,
        "vessel_diameter": 0.15,
        "area_m2": 0.02,
        "n_0": 80,
        "h_a": 0.0,
        "height_m": 0.55,
        "n_steps": 24,
        "search": _search(),
        "properties_at": properties_at,
    }
