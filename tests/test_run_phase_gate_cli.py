"""run_phase_gate CLI 与 harness 脚本 import 烟雾测试。"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


def test_run_phase_gate_registry_covers_all_ids() -> None:
    from run_phase_gate import GATE_ORDER, GATE_REGISTRY

    assert set(GATE_REGISTRY) == set(GATE_ORDER)


def test_run_gate_unknown_id_raises() -> None:
    from run_phase_gate import run_gate

    with pytest.raises(ValueError, match="Unknown gate"):
        run_gate("99")


def test_run_gate_json_output_smoke(monkeypatch: pytest.MonkeyPatch) -> None:
    """不跑真实 gate 物理检查，只验证 JSON 路径与 exit code。"""
    import run_phase_gate as rpg

    class _FakeGateResult:
        passed = True

        def to_dict(self) -> dict:
            return {
                "gate_id": "10",
                "phase": "smoke",
                "passed": True,
                "hard_failures": [],
                "warnings": [],
                "metrics": {},
            }

    monkeypatch.setitem(rpg.GATE_REGISTRY, "10", (lambda: _FakeGateResult(), rpg._print_gate_10))
    code = rpg.run_gate("10", json_only=True)
    assert code == 0


@pytest.mark.parametrize(
    "script_name",
    [
        "run_phase_gate.py",
        "audit_zone_profile_lu.py",
        "audit_flowchart_mapping.py",
        "audit_cell_conservation.py",
        "audit_freeboard_direction_lu.py",
        "validation_performance_report.py",
        "compare_vorab_march_lu.py",
        "benchmark_phase2_convergence_ladder.py",
    ],
)
def test_migrated_scripts_import_with_scripts_on_path(script_name: str) -> None:
    path = _SCRIPTS / script_name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
