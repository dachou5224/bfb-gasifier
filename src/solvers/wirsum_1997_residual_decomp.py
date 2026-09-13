"""R12 Startwert 之后拆 packed ‖F‖₂：气 / 固 / 能，以及 R5 对消。

同相 R12 把 ‖F‖₂ 从 10¹³ 降到 ~10⁷ 后，linf ≈ l2，说明只剩一条大方程。
本核只切片与分解，不再投影、不加 residual clip。

R5：CO + 0.5 O2 → CO2（Hamel Eq.5.36）。
与 ``same_phase_r12_startwert`` / 工程 fastox 分开命名。

Source: Hamel (1999) Eq. 2.4–2.7, 5.36–5.37; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, N_GAS
from src.solvers.nr_indexing import (
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_solid_var,
    n_var,
)
from src.solvers.wirsum_1997_elemental_decomp import extract_r12_extent
from src.solvers.wirsum_1997_startwert_r12 import apply_same_phase_r12_startwert_to_cells

_IDX_CO = GAS_SPECIES_INDEX["CO"]
_IDX_CO2 = GAS_SPECIES_INDEX["CO2"]
_IDX_O2 = GAS_SPECIES_INDEX["O2"]


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def extract_r5_extent(gas: npt.ArrayLike) -> dict[str, float]:
    """从气相源项抽出最接近 R5 化学计量的进度 ξ [mol/s]。

    理想：``ΔCO = −ξ``，``ΔO2 = −0.5 ξ``，``ΔCO2 = +ξ``。
    用 ``ξ = −ΔCO``。
    """
    g = np.asarray(gas, dtype=np.float64).reshape(-1)
    d_co = float(g[_IDX_CO])
    d_o2 = float(g[_IDX_O2])
    d_co2 = float(g[_IDX_CO2])
    xi = -d_co
    return {
        "extent_mol_s": xi,
        "d_CO": d_co,
        "d_O2": d_o2,
        "d_CO2": d_co2,
        "o2_over_co2": (d_o2 / d_co2) if d_co2 != 0.0 else float("nan"),
        "co_over_co2": (d_co / d_co2) if d_co2 != 0.0 else float("nan"),
        "stoich_residual_o2": d_o2 + 0.5 * xi,
        "stoich_residual_co2": d_co2 - xi,
    }


def _gas_linf(F_gas: npt.NDArray[np.float64]) -> dict[str, Any]:
    flat = np.asarray(F_gas, dtype=np.float64).reshape(-1)
    if flat.size == 0:
        return {"species": None, "phase": None, "value": None}
    n = min(int(flat.size), 2 * N_GAS)
    i = int(np.argmax(np.abs(flat[:n])))
    phase = "dense" if i < N_GAS else "bubble"
    sp = GAS_SPECIES[i % N_GAS]
    return {"species": sp, "phase": phase, "value": _finite(float(flat[i]))}


def decompose_cell_residual_families(cell: Any, index: int) -> dict[str, Any]:
    """一格的气 / 固 / 能残差切片。假定 ``residuals()`` 已跑过。"""
    F_gas = np.asarray(cell.residuals(groups=("gas",)), dtype=np.float64)
    F_solid = np.asarray(cell.residuals(groups=("solid",)), dtype=np.float64)
    F_e = np.asarray(cell.residuals(groups=("energy",)), dtype=np.float64)
    energy_W = float(F_e.reshape(-1)[0]) if F_e.size else 0.0
    gas = np.asarray(cell.R_gas_b, dtype=np.float64) + np.asarray(cell.R_gas_d, dtype=np.float64)
    r5 = extract_r5_extent(gas)
    r12 = extract_r12_extent(gas)
    gas_linf = _gas_linf(F_gas)
    xi5 = abs(float(r5["extent_mol_s"]))
    return {
        "cell": index,
        "gas_residual_l2": _finite(float(np.linalg.norm(F_gas))),
        "gas_residual_linf": _finite(float(np.linalg.norm(F_gas, ord=np.inf))),
        "gas_linf": gas_linf,
        "solid_residual_l2": _finite(float(np.linalg.norm(F_solid))),
        "energy_residual_W": _finite(energy_W),
        "r5": {k: _finite(v) for k, v in r5.items()},
        "r12": {k: _finite(v) for k, v in r12.items()},
        "r5_looks_like_pair": bool(
            xi5 > 1.0
            and abs(float(r5["co_over_co2"]) + 1.0) < 1e-3
            and abs(float(r5["o2_over_co2"]) + 0.5) < 1e-3
        ),
    }


def _packed_entry_labels(cell: Any, cell_index: int) -> list[str]:
    labels: list[str] = []
    gas_idx = main_nr_gas_idx_for_cell(cell)
    n_gas = int(n_gas_var(cell))
    if n_gas > 0:
        names = [GAS_SPECIES[int(j)] for j in gas_idx]
        if n_gas == 2 * len(names):
            labels.extend(f"c{cell_index}.N_d.{sp}" for sp in names)
            labels.extend(f"c{cell_index}.N_b.{sp}" for sp in names)
        else:
            labels.extend(f"c{cell_index}.N.{sp}" for sp in names)
    n_sol = int(n_solid_var(cell))
    labels.extend(f"c{cell_index}.solid[{k}]" for k in range(n_sol))
    if has_temperature_var(cell):
        labels.append(f"c{cell_index}.energy")
    expected = int(n_var(cell))
    if len(labels) != expected:
        raise AssertionError(
            f"packed residual labels {len(labels)} != n_var {expected} for cell {cell_index}"
        )
    return labels


def decompose_packed_residual(reactor: Any) -> dict[str, Any]:
    """packed F 的方程族切片 + |F| 最大项。不跑 Newton。"""
    from src.solvers.global_nr_solver import global_residual, pack_reactor

    reactor._apply_all_bc_for_nr()
    cells = list(reactor._solver_cells_for_nr())
    for cell in cells:
        cell.residuals()
    x = pack_reactor(cells)
    F = global_residual(x, cells, reactor._apply_all_bc_for_nr)
    labels: list[str] = []
    for i, cell in enumerate(cells):
        labels.extend(_packed_entry_labels(cell, i))
    if len(labels) != int(F.size):
        raise AssertionError(f"label count {len(labels)} != packed F size {F.size}")
    imax = int(np.argmax(np.abs(F)))
    families = [decompose_cell_residual_families(cell, i) for i, cell in enumerate(cells)]
    gas_l2 = float(np.sqrt(sum(float(row["gas_residual_l2"] or 0.0) ** 2 for row in families)))
    solid_l2 = float(np.sqrt(sum(float(row["solid_residual_l2"] or 0.0) ** 2 for row in families)))
    energy_l2 = float(np.sqrt(sum(float(row["energy_residual_W"] or 0.0) ** 2 for row in families)))
    parts = {"gas": gas_l2, "solid": solid_l2, "energy": energy_l2}
    dominant = max(parts, key=parts.get)
    linf_label = labels[imax]
    if linf_label.endswith(".energy"):
        linf_family = "energy"
    elif ".solid[" in linf_label:
        linf_family = "solid"
    else:
        linf_family = "gas"
    return {
        "residual_l2": _finite(float(np.linalg.norm(F, ord=2))),
        "residual_linf": _finite(float(np.linalg.norm(F, ord=np.inf))),
        "linf_index": imax,
        "linf_label": linf_label,
        "linf_family": linf_family,
        "linf_over_l2": _finite(abs(float(F[imax])) / float(np.linalg.norm(F, ord=2))),
        "family_l2": {k: _finite(v) for k, v in parts.items()},
        "dominant_family": dominant,
        "per_cell": families,
        "max_r5_extent_mol_s": _finite(
            max(abs(float(row["r5"]["extent_mol_s"] or 0.0)) for row in families)
        ),
        "max_r12_extent_mol_s": _finite(
            max(abs(float(row["r12"]["extent_mol_s"] or 0.0)) for row in families)
        ),
    }


def build_residual_decomp_payload() -> dict[str, Any]:
    """reduced 3-cell：R12 Startwert 前 / 后的 packed 残差族分解。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    before = decompose_packed_residual(reactor)
    apply_same_phase_r12_startwert_to_cells(list(reactor.cells))
    after = decompose_packed_residual(reactor)
    return {
        "source": "Wirsum alignment residual family decomposition after R12 Startwert",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "packed_residual_family_decomp",
        "not": [
            "same_phase_r5_startwert",
            "project_bed_holdup_fast_oxidation_closure",
            "residual_clip",
        ],
        "units": {
            "gas": "mol/s",
            "solid": "kg/s",
            "energy": "W",
            "r5_extent": "mol/s",
        },
        "before": before,
        "after": after,
        "note": (
            "R12 Startwert 之后看哪一族撑住 ‖F‖₂。"
            "能量残差单位是 W，与气相 mol/s 不可直接比大小；"
            "raw ‖F‖₂ 被 W 主导不证明质量已闭合。"
            "不用实验 CO。"
        ),
    }
