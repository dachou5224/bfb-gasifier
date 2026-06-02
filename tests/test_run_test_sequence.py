from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_run_test_sequence_module():
    repo_root = Path(__file__).resolve().parent.parent
    mod_path = repo_root / "scripts" / "run_test_sequence.py"
    spec = importlib.util.spec_from_file_location("run_test_sequence", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_stage_order_includes_module_audits_and_convergence_ladder():
    mod = _load_run_test_sequence_module()
    stage_ids = [stage.id for stage in mod.STAGES]
    assert "15_module_audits" in stage_ids
    assert "35_convergence_ladder" in stage_ids
    assert stage_ids.index("15_module_audits") < stage_ids.index("20_coupled_cell")
    assert stage_ids.index("35_convergence_ladder") < stage_ids.index("40_phase_gates")


def test_missing_promotion_prereqs_ignores_non_whole_model_stage():
    mod = _load_run_test_sequence_module()
    missing = mod._missing_promotion_prereqs("30_solver_policy", {"00_sanity"})
    assert missing == []


def test_missing_promotion_prereqs_reports_all_unpassed_gates():
    mod = _load_run_test_sequence_module()
    missing = mod._missing_promotion_prereqs(
        "50_whole_model",
        {"00_sanity", "10_independent_modules", "15_module_audits"},
    )
    assert missing == [
        "20_coupled_cell",
        "30_solver_policy",
        "35_convergence_ladder",
        "40_phase_gates",
    ]


def test_missing_promotion_prereqs_empty_when_all_prereqs_passed():
    mod = _load_run_test_sequence_module()
    missing = mod._missing_promotion_prereqs(
        "50_whole_model",
        set(mod.WHOLE_MODEL_PROMOTION_PREREQS),
    )
    assert missing == []


def test_convergence_and_whole_model_use_strict_freeboard_nr_audit():
    mod = _load_run_test_sequence_module()
    stage_map = {stage.id: stage for stage in mod.STAGES}
    assert "python3 scripts/audit_global_nr_profile_lu.py --strict" in stage_map["35_convergence_ladder"].commands
    assert "python3 scripts/audit_global_nr_profile_lu.py --strict" in stage_map["50_whole_model"].commands
    assert all("audit_phase1_htw_lu.py" not in cmd for cmd in stage_map["50_whole_model"].commands)
