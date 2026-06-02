"""Shared helpers for codebase slimming-plan validation tests."""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_INDEX = REPO_ROOT / "scripts" / "SCRIPTS_INDEX.md"
PYTEST_NODE_RE = re.compile(
    r"(tests/test_[\w.]+\.py(?:::[\w]+)?)",
)


def load_run_test_sequence_module():
    mod_path = REPO_ROOT / "scripts" / "run_test_sequence.py"
    spec = importlib.util.spec_from_file_location("run_test_sequence", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {mod_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_run_module_audits_module():
    mod_path = REPO_ROOT / "scripts" / "run_module_audits.py"
    spec = importlib.util.spec_from_file_location("run_module_audits", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {mod_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def parse_scripts_index_active() -> set[str]:
    text = SCRIPTS_INDEX.read_text(encoding="utf-8")
    active_block = text.split("## Active Scripts", 1)[1].split("## Deprecated", 1)[0]
    names: set[str] = set()
    for line in active_block.splitlines():
        line = line.strip()
        if not line.startswith("- `"):
            continue
        inner = line[3:]
        if "`" not in inner:
            continue
        names.add(inner.split("`", 1)[0])
    return names


def list_root_scripts() -> set[str]:
    return {path.name for path in (REPO_ROOT / "scripts").glob("*.py")}


def extract_pytest_nodes_from_command(cmd: str) -> list[str]:
    return PYTEST_NODE_RE.findall(cmd)


def pytest_node_exists(node: str) -> bool:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", node],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return False
    combined = proc.stdout + proc.stderr
    # ``--collect-only`` 仍会在摘要行打印 "no tests ran"，不能据此判失败。
    if "::test_" in combined or "::Test" in combined:
        return True
    return bool(re.search(r"\d+\s+tests?\s+collected", combined))


def collect_gate_pytest_references() -> list[tuple[str, str, str]]:
    refs: list[tuple[str, str, str]] = []
    seq = load_run_test_sequence_module()
    for stage in seq.STAGES:
        for cmd in stage.commands:
            if "pytest" not in cmd:
                continue
            for node in extract_pytest_nodes_from_command(cmd):
                refs.append(("run_test_sequence", stage.id, node))

    audits = load_run_module_audits_module()
    for mod in audits.MODULE_AUDITS:
        for cmd in mod.tests:
            if "pytest" not in cmd:
                continue
            for node in extract_pytest_nodes_from_command(cmd):
                refs.append(("run_module_audits", mod.module_id, node))
    return refs


def collect_module_audit_script_paths() -> list[str]:
    audits = load_run_module_audits_module()
    paths: list[str] = []
    for mod in audits.MODULE_AUDITS:
        for cmd in mod.audits:
            token = cmd.split()[1] if len(cmd.split()) > 1 else ""
            if token.startswith("scripts/"):
                paths.append(token)
    return paths


def grep_repo(pattern: str, *, globs: tuple[str, ...] = ("*.py",)) -> list[str]:
    hits: list[str] = []
    skip_prefixes = (".git/", ".venv/", "node_modules/", "archive/")
    for glob in globs:
        for path in REPO_ROOT.rglob(glob):
            rel = path.relative_to(REPO_ROOT).as_posix()
            if any(rel.startswith(prefix) for prefix in skip_prefixes):
                continue
            if "/__pycache__/" in rel:
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            if pattern in text:
                hits.append(rel)
    return sorted(set(hits))
