"""阶段 A：两条初始化路径的上下文契约快照。

路径
----
- ``reduced_3cell_x0``：``build_initialized_lu_reactor``，只 ``generate_initial_x0``
- ``phase2_3bed_2fb_init``：``build_initialized_thesis_freeboard_reactor``，走生产 ``init_precalc``

只初始化，不跑 Newton。数字是路径指纹，不是物性或收敛证明。

Source: docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 A
"""

from __future__ import annotations

import subprocess
from typing import Any

from src.core.cell import N_SOLID_COMP, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.species import GAS_SPECIES
from src.physics.wirsum_1997_freeboard import resolve_freeboard_model
from src.solvers.nr_indexing import (
    MAJOR_NR_GAS_SPECIES,
    PYROLYSIS_NR_GAS_SPECIES,
    main_nr_excluded_species,
    n_var,
    solid_component_indices_for_nr,
)
from src.solvers.vorabrechnung.vorab_policy import (
    resolve_vorab_outer_source_policy,
    thesis_vorab_sources_single_shot_enabled,
    vorab_outer_may_refresh_drying_pyro_sources,
    vorab_refresh_policy_label,
    vorab_tau_model_label,
)
from src.solvers.wirsum_1997_vorab_context import freeboard_trajectory_contract

CONTRACT_KNOBS: tuple[str, ...] = (
    "thesis_mode",
    "thesis_vorab_sources_single_shot",
    "vorab_outer_source_policy",
    "nr_solver_mode",
    "holdup_k0_anchor",
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis",
    "bed_solid_state_model",
    "energy_balance_heat_loss_mode",
    "solid_enthalpy_mode",
    "freeboard_trajectory_model",
    "explicit_freeboard_solver_graph_enabled",
    "n_cells",
    "n_freeboard_cells",
)

SOLID_COMPONENTS: tuple[str, ...] = ("char", "vm", "moisture", "ash")
_SOLID_IDX = (S_CHAR, S_VM, S_MOISTURE, S_ASH)


def _git_head() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    if out.returncode != 0:
        return None
    return str(out.stdout).strip() or None


def _knobs(cfg: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in CONTRACT_KNOBS:
        out[key] = getattr(cfg, key, None)
    return out


def snapshot_reactor_contract(reactor: Any, *, path: str, init_kind: str) -> dict[str, Any]:
    """从已初始化 reactor 抽出上下文契约（不评残差正确性）。"""
    cfg = reactor.config
    cells = list(reactor.cells)
    fb = list(getattr(reactor, "freeboard_cells", []) or [])
    all_cells = cells + fb
    excluded = main_nr_excluded_species(cfg)
    n_unknowns = int(sum(n_var(cell) for cell in all_cells))
    vm_cache = [bool(getattr(cell, "_vm_cache_valid", False)) for cell in cells]
    sources_frozen = [
        bool(getattr(cell, "_vorab_drying_pyro_sources_frozen", False)) for cell in cells
    ]
    bed_solid_nr = [
        [SOLID_COMPONENTS[int(i)] for i in solid_component_indices_for_nr(cell)]
        for cell in cells
    ]
    fb_solid_nr = [
        [SOLID_COMPONENTS[int(i)] for i in solid_component_indices_for_nr(cell)]
        for cell in fb
    ]
    return {
        "path": path,
        "init_kind": init_kind,
        "not_a_correctness_claim": True,
        "units": {
            "gas_holdup": "mol/s",
            "solid_holdup": "kg",
            "temperature": "K",
            "pressure": "Pa",
        },
        "knobs": _knobs(cfg),
        "n_bed_cells": int(len(cells)),
        "n_freeboard_cells": int(len(fb)),
        "n_unknowns": n_unknowns,
        "gas_species": list(GAS_SPECIES),
        "solid_components": list(SOLID_COMPONENTS),
        "solid_column_index": {name: int(idx) for name, idx in zip(SOLID_COMPONENTS, _SOLID_IDX)},
        "n_solid_comp": int(N_SOLID_COMP),
        "nr_gas_major": list(MAJOR_NR_GAS_SPECIES),
        "nr_gas_pyrolysis": list(PYROLYSIS_NR_GAS_SPECIES),
        "nr_gas_excluded": list(excluded),
        "nr_solid_bed": bed_solid_nr,
        "nr_solid_freeboard": fb_solid_nr,
        "vm_cache_valid": vm_cache,
        "drying_pyro_sources_frozen": sources_frozen,
        "thesis_vorab_sources_single_shot": thesis_vorab_sources_single_shot_enabled(cfg),
        "vorab_outer_source_policy": resolve_vorab_outer_source_policy(cfg),
        "vorab_outer_source_policy_configured": getattr(cfg, "vorab_outer_source_policy", "auto"),
        "outer_may_refresh_drying_pyro": vorab_outer_may_refresh_drying_pyro_sources(cfg),
        "vorab_refresh_policy": vorab_refresh_policy_label(cfg),
        "vorab_tau_model": vorab_tau_model_label(),
        "freeboard_trajectory": freeboard_trajectory_contract(reactor),
        "explicit_freeboard_solver_graph": bool(
            getattr(reactor, "_use_explicit_freeboard_solver_graph", lambda: False)()
        ),
        "freeboard_model_resolved": resolve_freeboard_model(
            str(getattr(cfg, "freeboard_trajectory_model", "analytical_wirsum"))
        ).name,
        "note": (
            "reduced 与 Phase2 生产 init 不是同一对象；"
            "single-shot 源项冻结是工程策略，不是已核对的 Hamel 外层。"
        ),
    }


def build_reduced_3cell_contract() -> dict[str, Any]:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    return snapshot_reactor_contract(
        reactor,
        path="reduced_3cell_x0",
        init_kind="generate_initial_x0_only",
    )


def build_phase2_init_contract() -> dict[str, Any]:
    from tests.global_nr_fixtures import build_initialized_thesis_freeboard_reactor

    reactor = build_initialized_thesis_freeboard_reactor(n_bed=3, n_freeboard=2)
    return snapshot_reactor_contract(
        reactor,
        path="phase2_3bed_2fb_init",
        init_kind="run_init_and_precalc_for_global_nr",
    )


def build_stage_a_baseline_payload() -> dict[str, Any]:
    return {
        "source": "Wirsum alignment stage A baseline contracts",
        "solved": False,
        "git_head": _git_head(),
        "paths": {
            "reduced_3cell_x0": build_reduced_3cell_contract(),
            "phase2_3bed_2fb_init": build_phase2_init_contract(),
        },
    }
