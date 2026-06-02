from __future__ import annotations

from pathlib import Path

import pytest

from src.core.phase_gate_checks import gate_00_baseline_frozen, load_phase_gate_baselines


def test_phase_gate_baseline_json_contains_gate_00_schema():
    data = load_phase_gate_baselines()
    gate = data["gate_00"]
    assert gate["schema_version"] == 1
    assert "cell_snapshot" in gate
    assert "reactor_snapshot" in gate
    assert "tolerances" in gate
    repo_root = Path(__file__).resolve().parent.parent
    assert (repo_root / "data" / "phase_gate_baselines.json").exists()


@pytest.mark.slow
def test_gate_00_baseline_frozen_smoke():
    result = gate_00_baseline_frozen()
    assert result.metrics["sanity"]["passed"], result.metrics["sanity"]["stdout"]
    assert result.passed, result.hard_failures
    assert result.artifacts["baseline_json"].endswith("data/phase_gate_baselines.json")
    assert "current_reactor_snapshot" in result.metrics
    assert "current_cell_snapshot" in result.metrics
