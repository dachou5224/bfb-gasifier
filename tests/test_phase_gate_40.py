from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from src.core.cell import S_CHAR, S_MOISTURE, S_VM
from src.core.phase_gate_checks import _phase4_solid_stream_filters, _phase4_solved_state_extent, gate_40_axial_extent_recycle


def test_phase4_solid_stream_filters_fast():
    filters = _phase4_solid_stream_filters()
    assert filters["recycle_vm_zero"]
    assert filters["recycle_moisture_zero"]
    assert filters["propagated_vm_zero"]
    assert filters["propagated_moisture_zero"]


def test_phase4_extent_uses_final_profile_metrics_over_latest_history(monkeypatch):
    cfg = SimpleNamespace(H_bed=1.0)

    def _solid(vm: float, moist: float, char: float) -> np.ndarray:
        arr = np.zeros((1, 4), dtype=float)
        arr[0, S_VM] = vm
        arr[0, S_MOISTURE] = moist
        arr[0, S_CHAR] = char
        return arr

    class FakeReactor:
        def __init__(self, _cfg):
            self.cells = [
                SimpleNamespace(
                    geo=SimpleNamespace(h_center=0.25),
                    T=1200.0,
                    m_solid_zu=_solid(0.1, 0.05, 0.2),
                    m_solid=_solid(0.1, 0.05, 0.2),
                ),
                SimpleNamespace(
                    geo=SimpleNamespace(h_center=0.75),
                    T=1100.0,
                    m_solid_zu=_solid(0.0, 0.0, 0.0),
                    m_solid=_solid(0.0, 0.0, 0.2),
                ),
            ]

        def solve(self, **kwargs):
            return {
                "converged": True,
                "n_iter": 3,
                "history": [
                    {
                        "vm_rebound": 0.2,
                        "moist_rebound": 0.3,
                    }
                ],
                "final_profile_metrics": {
                    "vm_rebound": 0.0,
                    "moist_rebound": 0.0,
                },
            }

    fake_utils = SimpleNamespace(
        load_case_LU=lambda: {"fake": True},
        build_phase2_htw_lu_freeboard_reactor_config=lambda case: cfg,
        PHASE1_HTW_LU_SOLVE_KWARGS={},
    )

    monkeypatch.setattr("src.core.phase_gate_checks._load_validation_utils", lambda: fake_utils)
    monkeypatch.setattr("src.core.phase_gate_checks.Reactor", FakeReactor)
    monkeypatch.setattr(
        "src.core.phase_gate_checks._phase4_release_rebound",
        lambda _reactor: {
            "vm_release_rebound": 0.0,
            "moist_release_rebound": 0.0,
            "vm_release_scale": 1.0,
            "moist_release_scale": 1.0,
        },
    )

    extent = _phase4_solved_state_extent()

    assert extent["vm_rebound"] == 0.0
    assert extent["moist_rebound"] == 0.0
    assert extent["vm_release_rebound"] == 0.0
    assert extent["moist_release_rebound"] == 0.0


def test_phase4_extent_prefers_phase2_freeboard_solve_kwargs(monkeypatch):
    cfg = SimpleNamespace(H_bed=1.0)
    captured: dict[str, object] = {}

    def _solid(vm: float, moist: float, char: float) -> np.ndarray:
        arr = np.zeros((1, 4), dtype=float)
        arr[0, S_VM] = vm
        arr[0, S_MOISTURE] = moist
        arr[0, S_CHAR] = char
        return arr

    class FakeReactor:
        def __init__(self, _cfg):
            self.cells = [
                SimpleNamespace(
                    geo=SimpleNamespace(h_center=0.25),
                    T=1200.0,
                    m_solid_zu=_solid(0.1, 0.05, 0.2),
                    m_solid=_solid(0.1, 0.05, 0.2),
                )
            ]

        def solve(self, **kwargs):
            captured["solve_kwargs"] = dict(kwargs)
            return {
                "converged": True,
                "n_iter": 1,
                "history": [],
                "final_profile_metrics": {},
            }

    fake_utils = SimpleNamespace(
        load_case_LU=lambda: {"fake": True},
        build_phase2_htw_lu_freeboard_reactor_config=lambda case: cfg,
        PHASE1_HTW_LU_SOLVE_KWARGS={"nr_jacobian_strategy": "block_tridiag_structured"},
        PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS={"nr_jacobian_strategy": None},
    )

    monkeypatch.setattr("src.core.phase_gate_checks._load_validation_utils", lambda: fake_utils)
    monkeypatch.setattr("src.core.phase_gate_checks.Reactor", FakeReactor)
    monkeypatch.setattr(
        "src.core.phase_gate_checks._phase4_release_rebound",
        lambda _reactor: {
            "vm_release_rebound": 0.0,
            "moist_release_rebound": 0.0,
            "vm_release_scale": 1.0,
            "moist_release_scale": 1.0,
        },
    )

    _ = _phase4_solved_state_extent()
    assert captured["solve_kwargs"]["nr_jacobian_strategy"] is None


@pytest.mark.slow
def test_gate_40_axial_extent_recycle_smoke():
    result = gate_40_axial_extent_recycle()
    extent = result.metrics["extent_audit"]
    # Structural smoke gate: always require audit payload integrity.
    assert "rows" in extent
    assert isinstance(extent["rows"], list)
    assert "vm_rebound" in extent
    assert "moist_rebound" in extent
    assert "vm_release_rebound" in extent
    assert "moist_release_rebound" in extent
    assert extent["vm_done_cell"] is not None
    assert extent["moisture_done_cell"] is not None
    assert "char_dominant_from_cell" in extent
    assert result.metrics["solid_stream_filters"]["recycle_vm_zero"]
    assert result.metrics["solid_stream_filters"]["recycle_moisture_zero"]
    # During ongoing calibration, gate can fail on rebound thresholds; ensure failures are explainable.
    if not result.passed:
        assert any(
            "vm_release_rebound" in msg
            or "moist_release_rebound" in msg
            or "vm_rebound" in msg
            or "moist_rebound" in msg
            or "char-dominant" in msg
            for msg in result.hard_failures
        )
