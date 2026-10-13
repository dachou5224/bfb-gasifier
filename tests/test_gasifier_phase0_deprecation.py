"""OPT-011：gasifier Phase0 GlobalState 路径弃用契约。"""

from __future__ import annotations

import importlib
import warnings

import pytest


def test_correlations_import_without_deprecation_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        importlib.import_module("src.gasifier.models.correlations.reaction_kinetics")


@pytest.mark.parametrize(
    "import_path,attr",
    [
        ("src.gasifier.models", "ConnectivityModel"),
        ("src.gasifier.domain", "GlobalState"),
    ],
)
def test_phase0_public_import_emits_deprecation(import_path: str, attr: str) -> None:
    mod = importlib.import_module(import_path)
    with pytest.warns(DeprecationWarning, match="Phase0"):
        getattr(mod, attr)


def test_gasifier_workflow_solvers_remain_canonical_reexports():
    from src.gasifier.solvers import solve_global_nr
    from src.gasifier.workflow import run_global_nr_workflow

    assert solve_global_nr.__module__.startswith("src.solvers")
    assert run_global_nr_workflow.__module__.startswith("src.workflow")
