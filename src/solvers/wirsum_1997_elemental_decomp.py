"""P5 续：把源项元素口拆成热解 / R12 对消 / 其余。

LU x0 上 H/O 开口（~10⁻² mol-atom/s）不是热解 An=b 失败。
裸 thesis R6/仓库 R12（Hamel Eq.5.38）在 2.5 MPa、~1300 K 给出 10¹³ mol/s
量级；H₂ / O₂ / H₂O 按 ``H2 + 0.5 O2 → H2O`` 对消后，账本只剩 float64 残差。

本核只分解，不加 limiter / fastox / residual clip。

Source: Hamel (1999) Eq. 5.38–5.39（thesis R6 = 仓库 R12）
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.cell import S_CHAR, S_MOISTURE, S_VM
from src.core.elemental_ledger import (
    TRACKED_ELEMENTS,
    combined_element_molar_rates,
    gas_element_molar_rates,
    max_absolute_element_rate,
    solid_element_molar_rates,
)
from src.core.species import GAS_SPECIES_INDEX

_IDX_H2 = GAS_SPECIES_INDEX["H2"]
_IDX_O2 = GAS_SPECIES_INDEX["O2"]
_IDX_H2O = GAS_SPECIES_INDEX["H2O"]


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _as_rates(values: dict[str, float]) -> dict[str, float | None]:
    return {el: _finite(values[el]) for el in TRACKED_ELEMENTS}


def cell_ledger_kwargs(cell: Any) -> dict[str, Any]:
    return dict(
        fuel_type=cell.fuel_type,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )


def ledger_gas_solid(
    gas: npt.ArrayLike,
    solid: npt.ArrayLike,
    cell: Any,
) -> dict[str, float]:
    return combined_element_molar_rates(
        gas_rates=np.asarray(gas, dtype=np.float64),
        solid_rates=np.asarray(solid, dtype=np.float64),
        **cell_ledger_kwargs(cell),
    )


def extract_r12_extent(gas: npt.ArrayLike) -> dict[str, float]:
    """从气相源项抽出最接近 R12 化学计量的进度 ξ [mol/s]。

    理想：``ΔH2 = −ξ``，``ΔO2 = −0.5 ξ``，``ΔH2O = +ξ``。
    用 ``ξ = −ΔH2``，再报与 O2 / H2O 的相对偏离。
    """
    g = np.asarray(gas, dtype=np.float64).reshape(-1)
    d_h2 = float(g[_IDX_H2])
    d_o2 = float(g[_IDX_O2])
    d_h2o = float(g[_IDX_H2O])
    xi = -d_h2
    return {
        "extent_mol_s": xi,
        "d_H2": d_h2,
        "d_O2": d_o2,
        "d_H2O": d_h2o,
        "o2_over_h2o": (d_o2 / d_h2o) if d_h2o != 0.0 else float("nan"),
        "h2_over_h2o": (d_h2 / d_h2o) if d_h2o != 0.0 else float("nan"),
        "stoich_residual_o2": d_o2 + 0.5 * xi,
        "stoich_residual_h2o": d_h2o - xi,
    }


def subtract_r12_pair(gas: npt.ArrayLike, extent: float) -> npt.NDArray[np.float64]:
    """从气相源项去掉 ``ξ (H2 + 0.5 O2 → H2O)``。"""
    out = np.asarray(gas, dtype=np.float64).copy()
    xi = float(extent)
    out[_IDX_H2] += xi
    out[_IDX_O2] += 0.5 * xi
    out[_IDX_H2O] -= xi
    return out


def decompose_cell_elemental_sources(cell: Any) -> dict[str, Any]:
    """分解一格的气+固源项元素账。不改 cell，假定 residuals() 已跑过。"""
    gas = np.asarray(cell.R_gas_b, dtype=np.float64) + np.asarray(cell.R_gas_d, dtype=np.float64)
    solid = np.asarray(cell.R_solid, dtype=np.float64)
    pyro_g = np.asarray(getattr(cell, "_vm_gas_source_cache", np.zeros_like(gas)), dtype=np.float64)
    pyro_s = np.asarray(getattr(cell, "_vm_solid_sink_cache", np.zeros_like(solid)), dtype=np.float64)
    if pyro_g.shape != gas.shape:
        pyro_g = np.zeros_like(gas)
    if pyro_s.shape != solid.shape:
        pyro_s = np.zeros_like(solid)

    total = ledger_gas_solid(gas, solid, cell)
    pyro = ledger_gas_solid(pyro_g, pyro_s, cell)
    gas_atoms = gas_element_molar_rates(gas, fuel_type=cell.fuel_type)
    solid_atoms = solid_element_molar_rates(solid, **{
        k: v for k, v in cell_ledger_kwargs(cell).items() if k != "fuel_type"
    })
    r12 = extract_r12_extent(gas)
    gas_wo_r12 = subtract_r12_pair(gas, r12["extent_mol_s"])
    after_r12 = ledger_gas_solid(gas_wo_r12, solid, cell)
    rest_gas = gas - pyro_g
    rest_solid = solid - pyro_s
    rest = ledger_gas_solid(rest_gas, rest_solid, cell)

    xi = abs(float(r12["extent_mol_s"]))
    leftover = max_absolute_element_rate(after_r12)
    solid_h = abs(float(solid_atoms["H"]))
    return {
        "total": _as_rates(total),
        "gas_only": _as_rates(gas_atoms),
        "solid_only": _as_rates(solid_atoms),
        "pyrolysis_cache": _as_rates(pyro),
        "pyrolysis_cache_used": bool(np.any(np.abs(pyro_g) > 0.0) or np.any(np.abs(pyro_s) > 0.0)),
        "rest_minus_cache": _as_rates(rest),
        "r12": {k: _finite(v) for k, v in r12.items()},
        "after_removing_r12": _as_rates(after_r12),
        "leftover_over_r12_extent": _finite(leftover / xi) if xi > 0.0 else None,
        "gas_solid_H_mismatch_over_solid_H": (
            _finite(abs(float(total["H"])) / solid_h) if solid_h > 1e-12 else None
        ),
        "has_solid_vm_or_moisture": bool(
            abs(float(np.sum(solid[:, S_VM]))) > 1e-16
            or abs(float(np.sum(solid[:, S_MOISTURE]))) > 1e-16
        ),
        "r12_looks_like_pair": bool(
            xi > 1.0
            and abs(float(r12["h2_over_h2o"]) + 1.0) < 1e-6
            and abs(float(r12["o2_over_h2o"]) + 0.5) < 1e-6
        ),
    }


def decompose_reactor_elemental_sources(cells: list[Any]) -> dict[str, Any]:
    per_cell = [decompose_cell_elemental_sources(c) for c in cells]
    max_xi = max(abs(float(row["r12"]["extent_mol_s"] or 0.0)) for row in per_cell)
    max_leftover = max(
        abs(float(row["after_removing_r12"][el] or 0.0))
        for row in per_cell
        for el in ("H", "O")
    )
    return {
        "source": "Wirsum alignment P5 elemental decomposition",
        "units": {"ledger": "mol-atom/s", "r12_extent": "mol/s"},
        "per_cell": per_cell,
        "max_r12_extent_mol_s": _finite(max_xi),
        "max_HO_after_r12": _finite(max_leftover),
        "note": (
            "R12（Eq.5.38）元素中性，只解释 10¹³ 量级质量残差。"
            "bed0 的 H/O 净口是气相原子流与固相 VM/水分账的差分；"
            "无固相的格子是其它均相反应在大数对消后的残差。"
            "都不是实验 CO，也不先加 limiter。"
        ),
    }


def build_p5_elemental_decomp_payload() -> dict[str, Any]:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    return decompose_reactor_elemental_sources(list(reactor.cells))
