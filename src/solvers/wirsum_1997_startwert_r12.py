"""同相 R12 Startwert：x0 上消掉同一相里的 H₂ ∩ O₂。

Hamel §2.1 的 Startwertwahl，不是 F(x) 限幅，也不是工程 fastox。
本核只做 thesis R6 / 仓库 R12 计量 ``H2 + 0.5 O2 → H2O``，按相独立。
不做 R5 / R6 / R10 / R1，也不改 residual。

与已有工程路径的命名边界：

- ``project_bed_holdup_fast_oxidation_closure``：R12→R5→R6→R10 + 可选 R1
- ``nr_line_search_fast_oxidation_projection_thesis``：Newton 线搜索投影
- 本核：``same_phase_r12_startwert``，只处理同相 H₂/O₂

Source: Hamel (1999) §2.1 Startwertwahl; Eq. 5.38–5.39
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import GAS_SPECIES_INDEX
from src.solvers.wirsum_1997_elemental_decomp import decompose_reactor_elemental_sources

_IDX_H2 = GAS_SPECIES_INDEX["H2"]
_IDX_O2 = GAS_SPECIES_INDEX["O2"]
_IDX_H2O = GAS_SPECIES_INDEX["H2O"]
_IDX_CO = GAS_SPECIES_INDEX["CO"]
_IDX_CH4 = GAS_SPECIES_INDEX["CH4"]


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def same_phase_r12_overlap_mol_s(N: npt.ArrayLike) -> float:
    """同相可发生的 R12 进度 ξ = min(N_H2, 2 N_O2) [mol/s]。"""
    n = np.asarray(N, dtype=np.float64).reshape(-1)
    n_h2 = max(float(n[_IDX_H2]), 0.0)  # [mol/s]
    n_o2 = max(float(n[_IDX_O2]), 0.0)  # [mol/s]
    return float(min(n_h2, 2.0 * n_o2))


def phase_fuel_oxygen_snapshot(N: npt.ArrayLike) -> dict[str, float]:
    """单相 holdup 上的燃料/氧重叠（只观测，不投影）。

    ξ12：H2 + 0.5 O2 → H2O，进度按消耗的 H2 [mol/s]。
    ξ5 / ξ6 只报告，本核不烧掉。
    """
    n = np.asarray(N, dtype=np.float64).reshape(-1)
    n_h2 = max(float(n[_IDX_H2]), 0.0)  # [mol/s]
    n_o2 = max(float(n[_IDX_O2]), 0.0)  # [mol/s]
    n_h2o = float(n[_IDX_H2O])  # [mol/s]
    n_co = max(float(n[_IDX_CO]), 0.0)  # [mol/s]
    n_ch4 = max(float(n[_IDX_CH4]), 0.0)  # [mol/s]
    return {
        "N_H2": n_h2,
        "N_O2": n_o2,
        "N_H2O": n_h2o,
        "N_CO": n_co,
        "N_CH4": n_ch4,
        "xi12_possible_mol_s": float(min(n_h2, 2.0 * n_o2)),
        "xi5_possible_mol_s": float(min(n_co, 2.0 * n_o2)),
        "xi6_possible_mol_s": float(min(n_ch4, n_o2 / 1.5)) if n_o2 > 0.0 else 0.0,
    }


def apply_same_phase_r12_startwert(N: npt.NDArray[np.float64]) -> dict[str, float]:
    """就地按 R12 计量消掉同相 H₂ ∩ O₂。只改 holdup，不改 F(x)。

    ξ = min(N_H2, 2 N_O2) [mol/s]：
    N_H2 ← N_H2 − ξ，N_O2 ← N_O2 − 0.5 ξ，N_H2O ← N_H2O + ξ。
    """
    n = np.asarray(N, dtype=np.float64)
    before = phase_fuel_oxygen_snapshot(n)
    xi = float(before["xi12_possible_mol_s"])  # [mol/s]
    n[_IDX_H2] = max(float(n[_IDX_H2]) - xi, 0.0)
    n[_IDX_O2] = max(float(n[_IDX_O2]) - 0.5 * xi, 0.0)
    n[_IDX_H2O] = float(n[_IDX_H2O]) + xi
    N[:] = n
    after = phase_fuel_oxygen_snapshot(n)
    return {
        "xi12_mol_s": xi,
        "N_H2_before": before["N_H2"],
        "N_O2_before": before["N_O2"],
        "N_H2_after": after["N_H2"],
        "N_O2_after": after["N_O2"],
        "overlap_after_mol_s": after["xi12_possible_mol_s"],
    }


def apply_same_phase_r12_startwert_to_cells(cells: list[Any]) -> dict[str, Any]:
    """对每格 dense / bubble holdup 独立做同相 R12 Startwert。"""
    per_cell: list[dict[str, Any]] = []
    xi_total = 0.0  # [mol/s]
    for i, cell in enumerate(cells):
        n_d = np.asarray(cell.N_d, dtype=np.float64)
        n_b = np.asarray(cell.N_b, dtype=np.float64)
        dense = apply_same_phase_r12_startwert(n_d)
        bubble = apply_same_phase_r12_startwert(n_b)
        cell.N_d[:] = n_d
        cell.N_b[:] = n_b
        if hasattr(cell, "_thermo_cache_valid"):
            cell._thermo_cache_valid = False
        xi_cell = float(dense["xi12_mol_s"]) + float(bubble["xi12_mol_s"])
        xi_total += xi_cell
        per_cell.append({"cell": i, "dense": dense, "bubble": bubble, "xi12_mol_s": xi_cell})
    return {
        "name": "same_phase_r12_startwert",
        "units": "mol/s",
        "xi12_total_mol_s": _finite(xi_total),
        "cells_touched": int(sum(1 for row in per_cell if float(row["xi12_mol_s"]) > 0.0)),
        "per_cell": per_cell,
    }


def _cell_overlap_row(cell: Any, index: int) -> dict[str, Any]:
    dense = phase_fuel_oxygen_snapshot(cell.N_d)
    bubble = phase_fuel_oxygen_snapshot(cell.N_b)
    return {
        "cell": index,
        "dense": dense,
        "bubble": bubble,
        "xi12_possible_mol_s": float(dense["xi12_possible_mol_s"] + bubble["xi12_possible_mol_s"]),
        "xi5_possible_mol_s": float(dense["xi5_possible_mol_s"] + bubble["xi5_possible_mol_s"]),
        "xi6_possible_mol_s": float(dense["xi6_possible_mol_s"] + bubble["xi6_possible_mol_s"]),
    }


def _max_field(rows: list[dict[str, Any]], key: str) -> float:
    return max(float(row[key]) for row in rows) if rows else 0.0


def snapshot_lu_r12_startwert_state(reactor: Any) -> dict[str, Any]:
    """x0 上的同相重叠 + 源项 R12 分解 + packed ‖F‖₂。不跑 Newton。"""
    from src.solvers.global_nr_solver import global_residual, pack_reactor

    reactor._apply_all_bc_for_nr()
    cells = list(reactor.cells)
    overlap = [_cell_overlap_row(cell, i) for i, cell in enumerate(cells)]
    for cell in cells:
        cell.residuals()
    decomp = decompose_reactor_elemental_sources(cells)
    nr_cells = reactor._solver_cells_for_nr()
    x = pack_reactor(nr_cells)
    F = global_residual(x, nr_cells, reactor._apply_all_bc_for_nr)
    return {
        "overlap": overlap,
        "max_xi12_possible_mol_s": _finite(_max_field(overlap, "xi12_possible_mol_s")),
        "max_xi5_possible_mol_s": _finite(_max_field(overlap, "xi5_possible_mol_s")),
        "max_xi6_possible_mol_s": _finite(_max_field(overlap, "xi6_possible_mol_s")),
        "max_r12_extent_mol_s": decomp["max_r12_extent_mol_s"],
        "max_HO_after_r12": decomp["max_HO_after_r12"],
        "residual_l2": _finite(float(np.linalg.norm(F, ord=2))),
        "residual_linf": _finite(float(np.linalg.norm(F, ord=np.inf))),
        "decomp": decomp,
    }


def build_startwert_r12_ab_payload() -> dict[str, Any]:
    """reduced 3-cell LU：关 / 开同相 R12 Startwert 的 A/B。不跑 150-iter NR。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    before = snapshot_lu_r12_startwert_state(reactor)
    applied = apply_same_phase_r12_startwert_to_cells(list(reactor.cells))
    after = snapshot_lu_r12_startwert_state(reactor)
    return {
        "source": "Wirsum alignment same-phase R12 Startwert A/B",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "same_phase_r12_startwert",
        "not": [
            "project_bed_holdup_fast_oxidation_closure",
            "nr_line_search_fast_oxidation_projection_thesis",
            "residual_clip",
        ],
        "units": {"holdup": "mol/s", "r12_extent": "mol/s", "residual": "mixed SI packed F"},
        "applied": applied,
        "before": before,
        "after": after,
        "note": (
            "只消同相 H2∩O2。CO/CH4 与 O2 的重叠留下给后续诊断。"
            "裸 R12 源项应从 10¹³ 落到热解/其它均相量级。"
            "不用实验 CO 做判据。"
        ),
    }
