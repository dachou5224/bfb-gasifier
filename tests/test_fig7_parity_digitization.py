"""Tests for Fig 7.10 / 7.13 / 7.17 parity digitization merge."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES_PATH = ROOT / "data" / "validation_cases.json"
PARITY_PATH = ROOT / "data" / "hamel_fig7_parity_digitized.json"


def _load() -> dict:
    return json.loads(CASES_PATH.read_text(encoding="utf-8"))


def test_parity_ssot_exists():
    data = json.loads(PARITY_PATH.read_text(encoding="utf-8"))
    assert "fig_7_10_vtt_pressurized" in data
    assert "13" in data["fig_7_10_vtt_pressurized"]["cases"]
    assert "F" in data["fig_7_17_wsv400"]["cases"]


def test_vtt_peat_13_outputs_filled():
    out = _load()["CASE_VTT_PRESSURIZED_PEAT_13"]["outputs"]
    assert out["exit_temperature_K"] == 1175
    assert out["carbon_conversion_pct"] == 97.0
    assert out["exit_gas_dry_mol_frac"]["CO"] == pytest.approx(0.113, abs=0.001)
    assert out["parity_plot_data"]["available"] is True


def test_vtt_sawdust_14_outputs_filled():
    out = _load()["CASE_VTT_PRESSURIZED_SAWDUST_14"]["outputs"]
    assert out["exit_temperature_K"] == 1045
    dry = out["exit_gas_dry_mol_frac"]
    assert dry["CO"] == pytest.approx(0.174, abs=0.001)
    assert dry["CH4"] == pytest.approx(0.058, abs=0.001)


def test_wsv_b_f_outputs_from_fig717():
    cases = _load()
    b = cases["CASE_WSV400_B"]["outputs"]
    f = cases["CASE_WSV400_F"]["outputs"]
    assert b["exit_temperature_K"] == 1010
    assert b["exit_gas_dry_mol_frac"]["CO"] == 0.155
    assert f["exit_gas_dry_mol_frac"]["H2"] == 0.15
    assert b["parity_plot_data"]["source"].startswith("Fig. 7.17")


def test_existing_axial_profiles_not_overwritten():
    """Peat6 exit T from Fig 7.9 axial profile must remain."""
    out = _load()["CASE_VTT_PRESSURIZED_PEAT_6"]["outputs"]
    assert out["exit_temperature_K"] == 1050
    assert out["axial_profiles"]["available"] is True
    assert out["parity_plot_data"]["T_FRB_K"]["meas"] == 1055


def test_fig713_ssot_exists():
    data = json.loads(PARITY_PATH.read_text(encoding="utf-8"))
    assert "fig_7_13_vtt_atmospheric" in data
    assert "E" in data["fig_7_13_vtt_atmospheric"]["cases"]


def test_vtt_atm_peat_e_outputs_filled():
    out = _load()["CASE_VTT_ATMOSPHERIC_PEAT_E"]["outputs"]
    assert out["exit_temperature_K"] == 1100
    assert out["carbon_conversion_pct"] == 91.0
    assert out["exit_gas_wet_mol_frac"]["CO"] == 0.112
    assert out["exit_gas_dry_mol_frac"]["CH4"] == pytest.approx(0.052, abs=0.002)
    assert out["parity_plot_data"]["source"].startswith("Fig. 7.13")
    assert out["parity_plot_data"]["T_FRB_K"]["meas"] == 1100


def test_vtt_atm_peat_d_keeps_fig712_axial():
    """Case D must keep Fig 7.12 axial outlet, not overwritten by parity."""
    out = _load()["CASE_VTT_ATMOSPHERIC_PEAT_D"]["outputs"]
    assert out["exit_temperature_K"] == 920
    assert out["axial_profiles"]["available"] is True
    assert "parity_plot_data" not in out or out.get("parity_plot_data", {}).get("available") is not True


def test_wsv_a_keeps_fig716_exit_t():
    out = _load()["CASE_WSV400_A"]["outputs"]
    assert out["exit_temperature_K"] == 980
    assert out["parity_plot_data"]["T_FRB_K"]["meas"] == 1135
