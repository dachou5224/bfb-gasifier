"""Eq. 2.7 能量项分解：R12 Startwert 之后 bed0 的 −10 MW。

只切片 ``H_in − H_out − Q``，不改 residual，不加 T fence / energy clip。
固相出口用 cell 现用的 ``_solid_energy_outlet_rates``（与 ``calc_energy_balance`` 同一载体）。

Source: Hamel (1999) Eq. 2.7–2.8; docs/hamel_submodels/02_cell_balances_and_exchange.md
"""

from __future__ import annotations

from typing import Any

import numpy as np

from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_energy_balance_H_in, calc_gas_enthalpy_flow
from src.core.energy_balance_defaults import (
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
)
from src.solvers.wirsum_1997_startwert_r12 import apply_same_phase_r12_startwert_to_cells

_SOLID_NAMES = ("char", "vm", "moisture", "ash")
_SOLID_IDX = (S_CHAR, S_VM, S_MOISTURE, S_ASH)


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _h_gas(cell: Any, N: Any, T: float) -> float:
    return float(calc_gas_enthalpy_flow(np.asarray(N, dtype=np.float64), float(T), cell._h_cache))


def _h_solid(cell: Any, m: Any, T: float) -> float:
    return float(cell._calc_solid_enthalpy_flow(np.asarray(m, dtype=np.float64), float(T)))


def _q_loss_W(cell: Any, H_in_W: float) -> tuple[float, str]:
    mode = str(getattr(cell, "heat_loss_mode", HEAT_LOSS_MODE_INLET_FRACTION)).strip().lower()
    if mode == HEAT_LOSS_MODE_INLET_FRACTION:
        return float(cell.heat_loss_frac) * max(float(H_in_W), 0.0), mode
    if mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT:
        return float(cell.Q_wall_W) + float(cell.Q_heat_exchanger_W), mode
    raise ValueError(f"Unsupported heat_loss_mode={mode!r}")


def _component_mass(m: Any) -> dict[str, float]:
    arr = np.asarray(m, dtype=np.float64)
    totals = np.sum(np.maximum(arr, 0.0), axis=0)
    return {name: float(totals[idx]) for name, idx in zip(_SOLID_NAMES, _SOLID_IDX)}


def _component_only(m: Any, comp: int) -> np.ndarray:
    out = np.zeros_like(np.asarray(m, dtype=np.float64))
    out[:, comp] = np.maximum(np.asarray(m, dtype=np.float64)[:, comp], 0.0)
    return out


def decompose_cell_eq27(cell: Any, index: int) -> dict[str, Any]:
    """一格 Eq. 2.7 的焓流分项 [W]。不改 cell。"""
    H_gas_axial_in = _h_gas(cell, cell.N_b_in + cell.N_d_in, cell.T_in_gas)
    H_gas_zu = _h_gas(cell, cell.N_zu_b + cell.N_zu_d, cell.T_zu_gas)
    H_gas_rez = _h_gas(cell, cell.N_rez_b + cell.N_rez_d, cell.T_rez_gas)
    H_solid_rez = _h_solid(cell, cell.m_solid_rez, cell.T_rez_solid)
    H_solid_in = _h_solid(cell, cell.m_solid_in, cell.T_in_solid)
    H_solid_auf = _h_solid(cell, cell.m_solid_auf_in, cell.T_solid_auf_in)
    H_solid_ab = _h_solid(cell, cell.m_solid_ab_in, cell.T_solid_ab_in)
    H_solid_zu = _h_solid(cell, cell.m_solid_zu, cell.T_zu_solid)
    H_in_parts = {
        "H_gas_axial_in": H_gas_axial_in,
        "H_gas_zu": H_gas_zu,
        "H_gas_rez": H_gas_rez,
        "H_solid_rez": H_solid_rez,
        "H_solid_axial_in": H_solid_in,
        "H_solid_auf_in": H_solid_auf,
        "H_solid_ab_in": H_solid_ab,
        "H_solid_zu": H_solid_zu,
    }
    H_in = float(sum(H_in_parts.values()))
    H_in_api = float(
        calc_energy_balance_H_in(
            N_b_in=cell.N_b_in,
            N_d_in=cell.N_d_in,
            T_in_gas=cell.T_in_gas,
            N_zu_b=cell.N_zu_b,
            N_zu_d=cell.N_zu_d,
            T_zu_gas=cell.T_zu_gas,
            N_rez_b=cell.N_rez_b,
            N_rez_d=cell.N_rez_d,
            T_rez_gas=cell.T_rez_gas,
            m_solid_rez=cell.m_solid_rez,
            T_rez_solid=cell.T_rez_solid,
            m_solid_in=cell.m_solid_in,
            T_in_solid=cell.T_in_solid,
            m_solid_zu=cell.m_solid_zu,
            T_zu_solid=cell.T_zu_solid,
            h_cache=cell._h_cache,
            m_solid_auf_in=cell.m_solid_auf_in,
            T_solid_auf_in=cell.T_solid_auf_in,
            m_solid_ab_in=cell.m_solid_ab_in,
            T_solid_ab_in=cell.T_solid_ab_in,
            ash_dry_wt=cell.solid.ash_dry_wt,
            VM_daf=cell.solid.VM_daf,
            h_f_dry=cell.solid.h_f_dry,
            solid_enthalpy_mode=cell.solid_enthalpy_mode,
        )
    )
    gas_out = np.asarray(cell.N_b, dtype=np.float64) + np.asarray(cell.N_d, dtype=np.float64)
    solid_out = np.asarray(cell._solid_energy_outlet_rates(), dtype=np.float64)
    H_gas_out = _h_gas(cell, gas_out, cell.T)
    H_solid_out = _h_solid(cell, solid_out, cell.T)
    H_out = H_gas_out + H_solid_out
    q_loss, loss_mode = _q_loss_W(cell, H_in)
    residual = H_in - H_out - q_loss
    live = float(cell.calc_energy_balance())

    m_in = (
        np.maximum(cell.m_solid_zu, 0.0)
        + np.maximum(cell.m_solid_in, 0.0)
        + np.maximum(cell.m_solid_rez, 0.0)
        + np.maximum(cell.m_solid_auf_in, 0.0)
        + np.maximum(cell.m_solid_ab_in, 0.0)
    )
    m_in_comp = _component_mass(m_in)
    m_out_comp = _component_mass(solid_out)
    H_zu_by_comp = {
        name: _h_solid(cell, _component_only(cell.m_solid_zu, idx), cell.T_zu_solid)
        for name, idx in zip(_SOLID_NAMES, _SOLID_IDX)
    }
    H_sout_by_comp = {
        name: _h_solid(cell, _component_only(solid_out, idx), cell.T)
        for name, idx in zip(_SOLID_NAMES, _SOLID_IDX)
    }
    holdup = _component_mass(cell.m_solid)
    k_auf = np.asarray(getattr(cell, "K_solid_auf", np.zeros_like(cell.m_solid)), dtype=np.float64)
    k_ab = np.asarray(getattr(cell, "K_solid_ab", np.zeros_like(cell.m_solid)), dtype=np.float64)
    k_auf_comp = {
        name: float(np.max(k_auf[:, idx])) if k_auf.ndim == 2 else float(np.max(k_auf))
        for name, idx in zip(_SOLID_NAMES, _SOLID_IDX)
    }
    k_ab_comp = {
        name: float(np.max(k_ab[:, idx])) if k_ab.ndim == 2 else float(np.max(k_ab))
        for name, idx in zip(_SOLID_NAMES, _SOLID_IDX)
    }
    dominant_in = max(H_in_parts, key=lambda k: abs(H_in_parts[k]))
    return {
        "cell": index,
        "T_K": _finite(cell.T),
        "T_in_gas_K": _finite(cell.T_in_gas),
        "T_zu_solid_K": _finite(cell.T_zu_solid),
        "heat_loss_mode": loss_mode,
        "H_in_parts_W": {k: _finite(v) for k, v in H_in_parts.items()},
        "H_in_W": _finite(H_in),
        "H_in_api_W": _finite(H_in_api),
        "H_in_parts_match_api": bool(abs(H_in - H_in_api) <= 1.0),
        "H_gas_out_W": _finite(H_gas_out),
        "H_solid_out_W": _finite(H_solid_out),
        "H_out_W": _finite(H_out),
        "Q_loss_W": _finite(q_loss),
        "residual_W": _finite(residual),
        "live_energy_residual_W": _finite(live),
        "decomp_matches_live": bool(abs(residual - live) <= 1.0),
        "dominant_H_in_term": dominant_in,
        "H_zu_by_comp_W": {k: _finite(v) for k, v in H_zu_by_comp.items()},
        "H_solid_out_by_comp_W": {k: _finite(v) for k, v in H_sout_by_comp.items()},
        "m_in_kg_s": {k: _finite(v) for k, v in m_in_comp.items()},
        "m_energy_out_kg_s": {k: _finite(v) for k, v in m_out_comp.items()},
        "dm_in_minus_energy_out_kg_s": {
            k: _finite(m_in_comp[k] - m_out_comp[k]) for k in _SOLID_NAMES
        },
        "vm_in_not_in_energy_out": bool(
            m_in_comp["vm"] > 1e-6 and m_out_comp["vm"] < 0.1 * m_in_comp["vm"]
        ),
        "solid_state_model": str(getattr(cell, "solid_state_model", "")),
        "m_holdup_kg": {k: _finite(v) for k, v in holdup.items()},
        "K_auf_max_1_per_s": {k: _finite(v) for k, v in k_auf_comp.items()},
        "K_ab_max_1_per_s": {k: _finite(v) for k, v in k_ab_comp.items()},
        "holdup_zero": bool(sum(holdup.values()) <= 1e-12),
        "H_solid_out_zero": bool(abs(H_solid_out) <= 1.0),
    }


def decompose_reactor_eq27(cells: list[Any]) -> dict[str, Any]:
    per_cell = [decompose_cell_eq27(cell, i) for i, cell in enumerate(cells)]
    bed0 = per_cell[0]
    return {
        "per_cell": per_cell,
        "bed0_residual_W": bed0["residual_W"],
        "bed0_dominant_H_in_term": bed0["dominant_H_in_term"],
        "bed0_vm_in_not_in_energy_out": bed0["vm_in_not_in_energy_out"],
        "bed0_holdup_zero": bed0["holdup_zero"],
        "bed0_H_solid_out_zero": bed0["H_solid_out_zero"],
        "max_abs_residual_W": _finite(
            max(abs(float(row["residual_W"] or 0.0)) for row in per_cell)
        ),
    }


def build_eq27_energy_decomp_payload() -> dict[str, Any]:
    """reduced 3-cell：R12 Startwert 前 / 后的 Eq. 2.7 分项。不跑 Newton。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    before = decompose_reactor_eq27(list(reactor.cells))
    apply_same_phase_r12_startwert_to_cells(list(reactor.cells))
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    after = decompose_reactor_eq27(list(reactor.cells))
    return {
        "source": "Wirsum alignment Eq.2.7 energy term decomposition",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "eq27_energy_term_decomp",
        "not": [
            "temperature_fence",
            "energy_clip",
            "vm_in_eq27_outlet",
            "same_phase_r5_startwert",
        ],
        "units": {"enthalpy": "W", "solid_mass": "kg/s", "temperature": "K"},
        "before": before,
        "after": after,
        "note": (
            "bed0 的 −10 MW 在 R12 Startwert 前已在。"
            "holdup_transport 下 H_aus,solid = K·m；reduced x0 的 m_solid=0，故 H_sout=0。"
            "H_zu 仍带入燃料焓（char+VM）。K_VM=0 是第二层；先有库存为零。"
            "禁止把 VM 补进 Eq.2.7 出口，禁止用 T clip 抹掉残差。"
            "不用实验 CO。"
        ),
    }
