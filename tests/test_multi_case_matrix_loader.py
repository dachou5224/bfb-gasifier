"""Multi-case 矩阵 loader 契约（不求解）。"""

from __future__ import annotations

import pytest

from tests.validation.case_lu_loader import (
    CASE_LU_VALIDATION_KEY,
    MULTI_CASE_MATRIX_KEYS,
    load_case_LU,
    load_validation_case_flat,
)


@pytest.mark.parametrize("case_key", MULTI_CASE_MATRIX_KEYS)
def test_load_validation_case_flat_has_reactor_fields(case_key: str) -> None:
    case = load_validation_case_flat(case_key)
    assert case["_source_case_key"] == case_key
    assert case["H_bed"] > 0.0
    assert case["P"] > 0.0
    assert case["ER"] > 0.0
    assert "primary_agent" in case
    assert "fuel_type" in case


def test_lu_flat_unchanged_secondary_hints() -> None:
    case = load_case_LU()
    assert case["_source_case_key"] == CASE_LU_VALIDATION_KEY
    assert case["primary_agent"] == "air_steam"
    assert case["secondary_injection_xi"] == pytest.approx(0.60)
    assert case["total_air_Nm3_h"] == pytest.approx(5536.9)


def test_wesseling_2_flat_o2_steam() -> None:
    case = load_validation_case_flat("CASE_HTW_WESSELING_2")
    assert case["primary_agent"] == "o2_steam"
    assert case["fuel_feed"] == pytest.approx(6858.9)
    assert case["secondary_injection_agent"] == "O2"


def test_vtt_peat_6_flat_biomass_and_fallbacks() -> None:
    case = load_validation_case_flat("CASE_VTT_PRESSURIZED_PEAT_6")
    assert case["fuel_type"] == "biomass"
    assert case["recirculation"] is False
    assert case["moisture_wt"] == pytest.approx(14.8)
    assert case["H_bed"] == pytest.approx(1.2)
    assert case["bed_diameter_m"] == pytest.approx(0.15)
    assert case["secondary_injection_xi"] == pytest.approx(0.6)
