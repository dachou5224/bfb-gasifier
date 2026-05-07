from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_run_module_audits_module():
    repo_root = Path(__file__).resolve().parent.parent
    mod_path = repo_root / "scripts" / "run_module_audits.py"
    spec = importlib.util.spec_from_file_location("run_module_audits", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_module_audit_matrix_has_unique_ids_and_commands():
    mod = _load_run_module_audits_module()
    module_ids = [item.module_id for item in mod.MODULE_AUDITS]
    assert len(module_ids) == len(set(module_ids))
    for item in mod.MODULE_AUDITS:
        assert item.tests
        assert item.audits
        assert item.pass_rule
        assert item.triage_hint


def test_pick_modules_selects_subset_and_validates_unknown():
    mod = _load_run_module_audits_module()
    picked = mod._pick_modules(["physics", "workflow"])
    assert [item.module_id for item in picked] == ["physics", "workflow"]

    try:
        mod._pick_modules(["unknown_module"])
    except ValueError as exc:
        assert "unknown_module" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError for unknown module id")
