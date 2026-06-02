from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_audit_module():
    repo_root = Path(__file__).resolve().parent.parent
    mod_path = repo_root / "scripts" / "audit_global_nr_profile_lu.py"
    spec = importlib.util.spec_from_file_location("audit_global_nr_profile_lu", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_strict_failures_requires_hamel_reduced_mode():
    mod = _load_audit_module()
    failures = mod._strict_failures(
        {
            "major_gibbs_solver_mode": "shadow_compare",
            "converged": True,
            "nr_monitor": {"x0_sanity": {"ok": True}},
        }
    )
    assert any("major_gibbs_solver_mode" in msg for msg in failures)


def test_strict_failures_empty_when_core_checks_pass():
    mod = _load_audit_module()
    failures = mod._strict_failures(
        {
            "major_gibbs_solver_mode": "hamel_reduced",
            "converged": True,
            "nr_monitor": {"x0_sanity": {"ok": True}},
        }
    )
    assert failures == []
