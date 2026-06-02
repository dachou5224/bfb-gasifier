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


def test_char_mass_ledger_reasons_gate_large_residual():
    mod = _load_audit_module()
    reasons = mod._char_mass_ledger_reasons(
        {
            "fresh_char_kg_s": 1.0,
            "residual_char_kg_s": 0.05,
            "size_migration_char_kg_s": 0.0,
        }
    )
    assert reasons == ["char_mass_residual>2.0%_fresh"]


def test_char_mass_ledger_reasons_accept_conservative_small_residual():
    mod = _load_audit_module()
    reasons = mod._char_mass_ledger_reasons(
        {
            "fresh_char_kg_s": 1.0,
            "residual_char_kg_s": 0.001,
            "size_migration_char_kg_s": 0.0,
        }
    )
    assert reasons == []


def test_char_mass_ledger_reasons_prefers_bed_transport_summary_when_available():
    mod = _load_audit_module()
    reasons = mod._char_mass_ledger_reasons(
        {
            "fresh_char_kg_s": 1.0,
            "residual_char_kg_s": 0.5,
            "size_migration_char_kg_s": 0.0,
        },
        bed_transport_summary={
            "max_cell_char_residual_norm_local": 0.05,
            "max_abs_auf_in_gap_vs_below_kg_s": 0.0,
            "max_abs_ab_in_gap_vs_above_kg_s": 0.0,
        },
    )
    assert reasons == []


def test_char_mass_ledger_reasons_flags_large_bed_transport_residual():
    mod = _load_audit_module()
    reasons = mod._char_mass_ledger_reasons(
        {
            "fresh_char_kg_s": 1.0,
            "residual_char_kg_s": 0.001,
            "size_migration_char_kg_s": 0.0,
        },
        bed_transport_summary={
            "max_cell_char_residual_norm_local": 0.12,
            "max_abs_auf_in_gap_vs_below_kg_s": 0.0,
            "max_abs_ab_in_gap_vs_above_kg_s": 0.0,
        },
    )
    assert reasons == ["bed_char_cell_residual>10.0%_local_transport"]
