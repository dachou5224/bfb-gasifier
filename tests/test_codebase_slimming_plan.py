"""Post-cleanup validation for the codebase slimming plan."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tests.slimming_plan_utils import (
    REPO_ROOT,
    collect_gate_pytest_references,
    collect_module_audit_script_paths,
    grep_repo,
    list_root_scripts,
    parse_scripts_index_active,
    pytest_node_exists,
)

INDEX_ACTIVE_SCRIPTS = {
    "audit_char_mass_conservation_lu.py",
    "audit_flowchart_mapping.py",
    "audit_hamel_consistency.py",
}

MOVED_TO_DEPRECATED = {
    "audit_bottom_bed_energy_balance_lu.py": "scripts/_deprecated/thermal/audit_bottom_bed_energy_balance_lu.py",
    "audit_phase2_freeboard_coeff_compare.py": "scripts/_deprecated/freeboard/audit_phase2_freeboard_coeff_compare.py",
    "audit_psi_b_strategy_initialized_lu.py": "scripts/_deprecated/hydrodynamics/audit_psi_b_strategy_initialized_lu.py",
    "benchmark_structured_jacobian.py": "scripts/_deprecated/nr/benchmark_structured_jacobian.py",
    "audit_gs_history_lu.py": "scripts/_deprecated/nr/audit_gs_history_lu.py",
    "audit_solver_parity_lu.py": "scripts/_deprecated/nr/audit_solver_parity_lu.py",
}

REMOVED_RUNTIME_MODULES = (
    "src/core/legacy_gs_support.py",
    "src/solvers/gauss_seidel_reactor_solver.py",
    "src/solvers/cell_outer_loop.py",
)

STALE_GATE_PYTEST_NODES = {
    "tests/test_global_nr_solver.py::test_thesis_mode_forces_global_nr_even_when_solver_arg_is_gauss_seidel",
    "tests/test_global_nr_solver.py::test_gauss_seidel_solver_allowed_with_legacy_flag",
    "tests/test_global_nr_solver.py::test_non_thesis_mode_still_disables_gs_warmup_when_legacy_gs_is_off",
    "tests/test_global_nr_solver.py::test_global_nr_stops_when_outer_is_matched_and_inner_stalls",
}


@pytest.mark.slimming_audit
class TestPostCleanupGateIntegrity:
    def test_all_gate_pytest_nodes_exist(self) -> None:
        missing = [
            node
            for _, _, node in collect_gate_pytest_references()
            if not pytest_node_exists(node)
        ]
        assert missing == []

    def test_stale_gate_nodes_are_not_referenced(self) -> None:
        referenced = {node for _, _, node in collect_gate_pytest_references()}
        assert referenced.isdisjoint(STALE_GATE_PYTEST_NODES)


@pytest.mark.slimming_audit
class TestPostCleanupScriptLayout:
    def test_moved_scripts_not_in_scripts_root(self) -> None:
        root = list_root_scripts()
        assert set(MOVED_TO_DEPRECATED).isdisjoint(root)

    def test_moved_scripts_exist_under_deprecated(self) -> None:
        for rel in MOVED_TO_DEPRECATED.values():
            assert (REPO_ROOT / rel).is_file(), rel

    def test_scripts_index_covers_gate_active_scripts(self) -> None:
        indexed = parse_scripts_index_active()
        assert INDEX_ACTIVE_SCRIPTS.issubset(indexed)

    def test_module_audit_scripts_exist_on_disk(self) -> None:
        missing = [
            rel
            for rel in collect_module_audit_script_paths()
            if not (REPO_ROOT / rel).is_file()
        ]
        assert missing == []


@pytest.mark.slimming_audit
class TestPostCleanupSourceLayout:
    def test_removed_runtime_modules_are_gone(self) -> None:
        present = [rel for rel in REMOVED_RUNTIME_MODULES if (REPO_ROOT / rel).is_file()]
        assert present == []

    def test_reactor_has_no_gauss_seidel_solver_method(self) -> None:
        source = (REPO_ROOT / "src/core/reactor.py").read_text(encoding="utf-8")
        assert "_solve_gauss_seidel" not in source
        assert "allow_legacy_gs" not in source

    def test_phase_gate_deep_audit_script_points_to_deprecated_path(self) -> None:
        from src.core import phase_gate_checks

        source = Path(phase_gate_checks.__file__).read_text(encoding="utf-8")
        assert "scripts/_deprecated/thermal/audit_drying_pyrolysis_extent_lu.py" in source

    def test_no_active_python_imports_deprecated_scripts_dir(self) -> None:
        hits = [
            path
            for path in grep_repo("scripts/_deprecated")
            if not path.startswith("scripts/_deprecated/")
            and path
            not in {
                "scripts/SCRIPTS_INDEX.md",
                "scripts/_deprecated/README.md",
                "tests/test_codebase_slimming_plan.py",
                "tests/slimming_plan_utils.py",
                "scripts/validate_slimming_plan.py",
                "src/core/phase_gate_checks.py",
                "tests/validation_case_utils.py",
            }
        ]
        assert hits == []


@pytest.mark.slimming_audit
class TestActiveCoreHealth:
    def test_reactor_rejects_gauss_seidel_solver(self) -> None:
        from src.core.reactor import Reactor
        from tests.validation_case_utils import build_phase1_htw_lu_reactor_config

        reactor = Reactor(build_phase1_htw_lu_reactor_config())
        with pytest.raises(ValueError, match="global_nr"):
            reactor.solve(solver="gauss_seidel", max_global_iter=1)

    def test_sanity_checks_pass(self) -> None:
        proc = subprocess.run(
            [sys.executable, "tests/sanity_checks.py"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
