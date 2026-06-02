from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_audit_module():
    repo_root = Path(__file__).resolve().parent.parent
    mod_path = repo_root / "scripts" / "audit_phase1_htw_lu.py"
    spec = importlib.util.spec_from_file_location("audit_phase1_htw_lu", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_overall_validation_pass_requires_convergence_gate():
    mod = _load_audit_module()
    ok = mod._overall_validation_pass(
        validation_candidate_ok=False,
        pass_T=True,
        pass_species=True,
        carbon_ok=True,
        relax=False,
    )
    assert ok is False


def test_overall_validation_pass_accepts_all_checks_when_candidate_is_valid():
    mod = _load_audit_module()
    ok = mod._overall_validation_pass(
        validation_candidate_ok=True,
        pass_T=True,
        pass_species=True,
        carbon_ok=True,
        relax=False,
    )
    assert ok is True


def test_overall_validation_pass_relax_mode_overrides_failures():
    mod = _load_audit_module()
    ok = mod._overall_validation_pass(
        validation_candidate_ok=False,
        pass_T=False,
        pass_species=False,
        carbon_ok=False,
        relax=True,
    )
    assert ok is True
