"""Master table JSON/MD 与 validation_cases.json 一致性契约。"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES_PATH = ROOT / "data" / "validation_cases.json"
LITERATURE_PATH = ROOT / "data" / "validation_plant_literature.json"
MASTER_JSON = ROOT / "data" / "validation_cases_master_table.json"
MASTER_MD = ROOT / "docs" / "validation_cases_master_table.md"
BUILD_SCRIPT = ROOT / "scripts" / "build_validation_cases_master_table.py"

EXPECTED_CASE_IDS: tuple[str, ...] = (
    "CASE_HTW_WESSELING_1",
    "CASE_HTW_WESSELING_2",
    "CASE_HTW_BERRENRATH",
    "CASE_VTT_PRESSURIZED_COAL_1",
    "CASE_VTT_PRESSURIZED_PEAT_6",
    "CASE_VTT_PRESSURIZED_PEAT_13",
    "CASE_VTT_PRESSURIZED_SAWDUST_14",
    "CASE_VTT_ATMOSPHERIC_PEAT_D",
    "CASE_VTT_ATMOSPHERIC_PEAT_E",
    "CASE_WSV400_A",
    "CASE_WSV400_B",
    "CASE_WSV400_C",
    "CASE_WSV400_D",
    "CASE_WSV400_E",
    "CASE_WSV400_F",
)


def _load_master() -> dict:
    return json.loads(MASTER_JSON.read_text(encoding="utf-8"))


def test_master_table_has_15_cases():
    data = _load_master()
    assert data["_meta"]["case_count"] == 15
    assert data["_meta"]["plant_group_count"] == 4
    ids = [c["case_id"] for c in data["cases"]]
    assert ids == list(EXPECTED_CASE_IDS)


def test_four_validation_plant_groups():
    data = _load_master()
    plants = data["validation_plants"]
    assert set(plants) == {
        "PLANT_HTW_WESSELING",
        "PLANT_HTW_BERRENRATH",
        "PLANT_VTT",
        "PLANT_WSV400",
    }
    assert plants["PLANT_HTW_WESSELING"]["case_count"] == 2
    assert plants["PLANT_HTW_BERRENRATH"]["case_count"] == 1
    assert plants["PLANT_VTT"]["case_count"] == 7
    assert plants["PLANT_WSV400"]["case_count"] == 6


@pytest.mark.parametrize(
    ("case_id", "plant_group_id"),
    [
        ("CASE_HTW_WESSELING_1", "PLANT_HTW_WESSELING"),
        ("CASE_HTW_WESSELING_2", "PLANT_HTW_WESSELING"),
        ("CASE_HTW_BERRENRATH", "PLANT_HTW_BERRENRATH"),
        ("CASE_VTT_PRESSURIZED_COAL_1", "PLANT_VTT"),
        ("CASE_VTT_ATMOSPHERIC_PEAT_E", "PLANT_VTT"),
        ("CASE_WSV400_F", "PLANT_WSV400"),
    ],
)
def test_case_plant_group_mapping(case_id: str, plant_group_id: str):
    rec = next(c for c in _load_master()["cases"] if c["case_id"] == case_id)
    assert rec["plant_group"]["plant_group_id"] == plant_group_id


def test_master_md_exists_and_lists_all_cases():
    text = MASTER_MD.read_text(encoding="utf-8")
    assert "Validation Cases 主表" in text
    assert "## 7. 原始文献" in text
    for case_id in EXPECTED_CASE_IDS:
        assert case_id in text


def test_literature_ssot_covers_all_plants():
    lit = json.loads(LITERATURE_PATH.read_text(encoding="utf-8"))
    plant_map = lit["plant_literature_map"]
    refs = lit["refs"]
    assert set(plant_map) == {
        "PLANT_HTW_WESSELING",
        "PLANT_HTW_BERRENRATH",
        "PLANT_VTT",
        "PLANT_WSV400",
    }
    for plant_id, mapping in plant_map.items():
        for role in ("primary_data", "secondary_compilation"):
            for ref_id in mapping[role]:
                assert ref_id in refs, f"{plant_id}.{role} unknown ref {ref_id}"


def test_master_json_includes_literature_with_key_dois():
    data = _load_master()
    lit = data["validation_literature"]
    assert lit["refs"]["kurkela_stahlberg_1992_part1"]["doi"] == "10.1016/0378-3820(92)90038-R"
    assert lit["refs"]["leppalahti_kurkela_1991_fuel"]["doi"] == "10.1016/0016-2361(91)90026-7"
    assert data["_meta"]["literature_ref_count"] >= 10


@pytest.mark.parametrize("case_id", EXPECTED_CASE_IDS)
def test_each_case_has_core_fields(case_id: str):
    rec = next(c for c in _load_master()["cases"] if c["case_id"] == case_id)
    assert rec["plant"]
    assert rec["fuel"]["type"]
    assert rec["reactor"]["bed_diameter_m"]
    assert rec["inlets"]["pressure_MPa"]
    assert rec["inlets"]["ER_lambda_O2"]
    assert rec["outlets"]["exit_temperature_K"]


def test_lu1_exit_gas_matches_ssot():
    master = next(c for c in _load_master()["cases"] if c["case_id"] == "CASE_HTW_WESSELING_1")
    ssot = json.loads(CASES_PATH.read_text(encoding="utf-8"))["CASE_HTW_WESSELING_1"]["outputs"]
    dry = ssot["exit_gas_dry_mol_frac"]
    assert master["outlets"]["gas_basis"] == "dry"
    assert float(master["outlets"]["CO"]) == pytest.approx(float(dry["CO"]), rel=1e-3)
    assert master["outlets"]["exit_temperature_K"] == str(ssot["exit_temperature_K"])


def test_vtt_atm_e_wet_basis():
    rec = next(c for c in _load_master()["cases"] if c["case_id"] == "CASE_VTT_ATMOSPHERIC_PEAT_E")
    assert rec["outlets"]["gas_basis"] == "wet"
    assert rec["outlets"]["H2O"] == "0.142"
    assert rec["outlets"]["carbon_conversion_pct"] == "91"


def test_build_script_regenerates_json():
    before = MASTER_JSON.read_text(encoding="utf-8")
    subprocess.run(
        [sys.executable, str(BUILD_SCRIPT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    after = MASTER_JSON.read_text(encoding="utf-8")
    assert json.loads(before) == json.loads(after)
