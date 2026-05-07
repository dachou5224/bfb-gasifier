from __future__ import annotations

import json
from pathlib import Path


def _load_validation_cases() -> dict:
    root = Path(__file__).resolve().parent.parent
    with open(root / "data" / "validation_cases.json", encoding="utf-8") as f:
        return json.load(f)


def test_literature_regression_baseline_case_exists_and_is_stable() -> None:
    """文献回归基线：关键工况及核心参考量必须固定存在。"""
    data = _load_validation_cases()
    assert "CASE_HTW_WESSELING_1" in data
    case = data["CASE_HTW_WESSELING_1"]

    inp = case["inputs"]
    out = case["outputs"]
    oc = inp["operating_conditions"]
    prox = inp["fuel"]["proximate_analysis"]

    # 基线字段稳定性（防止静默漂移/误改）
    assert abs(float(oc["pressure_MPa"]) - 2.5) < 1e-12
    assert abs(float(oc["ER"]) - 0.337) < 1e-12
    assert abs(float(prox["moisture_wt_pct"]) - 16.9) < 1e-12

    # 回归目标字段必须存在（具体误差带由端到端测试约束）
    assert out.get("exit_temperature_K") is not None
    assert out.get("carbon_conversion_pct") is not None
    assert "exit_gas_dry_mol_frac" in out and "CO" in out["exit_gas_dry_mol_frac"]
