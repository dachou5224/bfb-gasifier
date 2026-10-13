"""Connectivity submodule — see package ``__init__`` for public API."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, cp_ash, cp_char, cp_sand
from src.core.connectivity_graph import (
    ConnectivityBlock,
    ConnectivityEdge,
    ConnectivityGraph,
    bed_top_freeboard_ab_source_block_id,
    block_id_bed_index,
    build_connectivity_graph_from_reactor_params,
    cyclone_inlet_source_block_id,
    axial_upstream_block_id,
    freeboard_chain_inlet_source_block_id,
    infer_stream_class,
    ordered_freeboard_block_ids,
    resolve_bed_neighbor_block_ids,
    resolve_registry_cell,
    solid_ab_freeboard_downstream_block_id,
    solid_ab_source_block_id,
    solid_auf_source_block_id,
)
from src.solvers.nr_indexing import (
    MAJOR_NR_GAS_SPECIES,
    nr_dense_only_pyrolysis_support_species,
    nr_projection_gas_species,
)

_F_W = 0.25
_ZETA_W = 0.40

from src.core.connectivity.solid_transport import recycled_solid_stream_from_cell


def _bottom_hydrodynamic_gas_inlet_dense_fraction(cell: Cell) -> float | None:
    """Return dense-phase superficial flux share from current hydrodynamics."""
    u0 = float(getattr(cell, "u0", 0.0))
    u_d = float(getattr(cell, "u_d", 0.0))
    eps_b = float(getattr(cell, "eps_b", 0.0))
    dense_flux = max(u_d, 0.0) * max(1.0 - eps_b, 0.0)
    if not np.isfinite(u0) or not np.isfinite(dense_flux) or u0 <= 1e-12:
        return None
    return float(np.clip(dense_flux / u0, 0.0, 1.0))

def snapshot_bottom_gas_inlet_split_from_vorabrechnung(cells: list[Cell], cfg: Any) -> None:
    """Freeze the bottom gas split after Vorabrechnung hydrodynamics has been refreshed."""
    if not cells:
        return
    cell = cells[0]
    strategy = str(getattr(cfg, "gas_inlet_split_strategy", "fixed")).strip().lower()
    if strategy not in {
        "precalc_hydrodynamic_flux",
        "hydrodynamic_flux",
        "vorabrechnung_hydrodynamic_flux",
    }:
        cell._vorab_bottom_gas_inlet_dense_frac = np.nan
        return
    frac = _bottom_hydrodynamic_gas_inlet_dense_fraction(cell)
    if frac is not None:
        cell._vorab_bottom_gas_inlet_dense_frac = float(frac)

def align_bottom_primary_gas_state_to_inlet_split(cells: list[Cell], cfg: Any) -> None:
    """Align bottom-cell primary gas inventory split with the frozen inlet split.

    Init / Vorab→cell x₀ mapping only (Hamel Kap.2 p.16). Preserves each primary
    species' total molar flow and only changes the dense/bubble partition when
    building the initial cell state from Vorabrechnung budgets.

    Must not run on outer Abgleich refresh: NR unknowns N_d/N_b must stay consistent
    with ex_bd closure at the current linearization point.
    """
    if not cells:
        return
    cell = cells[0]
    dense_frac = resolve_bottom_gas_inlet_dense_fraction(cell, cfg)
    o2_dense_frac = resolve_bottom_o2_dense_fraction(cell, cfg)
    idx = GAS_SPECIES_INDEX
    for sp in ("H2O", "N2"):
        j = idx[sp]
        total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
        cell.N_d[j] = max(total * dense_frac, 1e-12)
        cell.N_b[j] = max(total * (1.0 - dense_frac), 1e-12)
    j = idx["O2"]
    total_o2 = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
    cell.N_d[j] = max(total_o2 * o2_dense_frac, 1e-12)
    cell.N_b[j] = max(total_o2 * (1.0 - o2_dense_frac), 1e-12)

def resolve_bottom_gas_inlet_dense_fraction(cell: Cell, cfg: Any) -> float:
    """Resolve bottom gas split from fixed config or Vorabrechnung hydrodynamic flux.

    Ref: Hamel (1999) Kapitel 2.1 Vorabrechnung and Kapitel 3.1.2 Eq.3.7-3.10.
    The suspension-phase share follows the dense-phase superficial flux
    ``u_d * (1 - eps_b)`` divided by the total superficial velocity ``u0``.
    """
    fallback = float(np.clip(getattr(cfg, "gas_inlet_dense_frac", 0.0), 0.0, 1.0))
    strategy = str(getattr(cfg, "gas_inlet_split_strategy", "fixed")).strip().lower()
    if strategy in {"", "fixed", "legacy_fixed"}:
        return fallback
    if strategy not in {
        "precalc_hydrodynamic_flux",
        "hydrodynamic_flux",
        "vorabrechnung_hydrodynamic_flux",
    }:
        raise ValueError(
            "Unsupported gas_inlet_split_strategy="
            f"{getattr(cfg, 'gas_inlet_split_strategy', None)!r}; "
            "expected 'fixed' or 'precalc_hydrodynamic_flux'."
        )

    frozen = float(getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", np.nan))
    if np.isfinite(frozen):
        return float(np.clip(frozen, 0.0, 1.0))
    live = _bottom_hydrodynamic_gas_inlet_dense_fraction(cell)
    if live is not None:
        return live
    return fallback

def resolve_bottom_o2_dense_fraction(cell: Cell, cfg: Any) -> float:
    """Resolve bottom O₂ split; falls back to general gas split when unset."""
    o2_specific = getattr(cfg, "gas_inlet_o2_dense_frac", None)
    if o2_specific is not None:
        return float(np.clip(float(o2_specific), 0.05, 0.95))
    return resolve_bottom_gas_inlet_dense_fraction(cell, cfg)

def _bed_cell_indices_from_cells(cells: list[Cell]) -> list[int]:
    return [i for i, cell in enumerate(cells) if str(getattr(cell, "cell_type", "bed")) == "bed"]

def _staged_primary_gas_layer_weights(n_layers: int, *, decay: float = 1.5) -> np.ndarray:
    """锥段多喷嘴占位：低层 bed 格权重更大 [-]。"""
    if n_layers <= 0:
        return np.array([], dtype=np.float64)
    if n_layers == 1:
        return np.array([1.0], dtype=np.float64)
    raw = np.exp(-float(decay) * np.arange(n_layers, dtype=np.float64))
    return raw / max(float(np.sum(raw)), 1e-12)

def _resolve_staged_primary_gas_layer_weights(cfg: Any, n_layers: int) -> np.ndarray:
    """解析分级进料层权重：显式 tuple 优先，否则指数衰减。"""
    custom = getattr(cfg, "vorab_staged_primary_gas_layer_weights_thesis", None)
    if custom is not None:
        raw = np.asarray(tuple(custom), dtype=np.float64).reshape(-1)
        if raw.size < n_layers:
            raise ValueError(
                f"vorab_staged_primary_gas_layer_weights_thesis 长度 {raw.size} "
                f"< staged bed cells {n_layers}"
            )
        weights = raw[:n_layers]
        if np.any(weights < 0.0):
            raise ValueError("vorab_staged_primary_gas_layer_weights_thesis 不得含负值")
        total = float(np.sum(weights))
        if total <= 0.0:
            raise ValueError("vorab_staged_primary_gas_layer_weights_thesis 之和须 > 0")
        return weights / total
    decay = float(getattr(cfg, "vorab_staged_primary_gas_weight_decay_thesis", 1.5))
    return _staged_primary_gas_layer_weights(n_layers, decay=decay)

def _clear_bed_primary_gas_zu(cells: list[Cell], bed_indices: list[int]) -> None:
    idx = GAS_SPECIES_INDEX
    for bi in bed_indices:
        cell = cells[bi]
        for sp in ("O2", "H2O", "N2"):
            j = int(idx[sp])
            cell.N_zu_d[j] = 0.0
            cell.N_zu_b[j] = 0.0

def _assign_primary_gas_zu_fraction(
    cell: Cell,
    *,
    O2_mol_s: float,
    H2O_mol_s: float,
    N2_mol_s: float,
    dense_frac: float,
    o2_dense_frac: float | None = None,
) -> None:
    idx = GAS_SPECIES_INDEX
    bubble_frac = 1.0 - float(np.clip(dense_frac, 0.05, 0.95))
    dense_frac = 1.0 - bubble_frac
    o2_dense = float(np.clip(o2_dense_frac if o2_dense_frac is not None else dense_frac, 0.05, 0.95))
    o2_bubble = 1.0 - o2_dense
    for sp, total, d_frac, b_frac in (
        ("O2", O2_mol_s, o2_dense, o2_bubble),
        ("H2O", H2O_mol_s, dense_frac, bubble_frac),
        ("N2", N2_mol_s, dense_frac, bubble_frac),
    ):
        j = int(idx[sp])
        cell.N_zu_d[j] += float(total) * d_frac
        cell.N_zu_b[j] += float(total) * b_frac

def distribute_staged_primary_gas_to_bed_cells(cells: list[Cell], cfg: Any) -> dict[str, float]:
    """沿前 n 个床格分配主气化剂 ``N_zu``（HTW 多段喷嘴占位）。

    bed0 密/气泡份额仍走 ``resolve_bottom_gas_inlet_dense_fraction``；
    上层格使用 ``vorab_staged_primary_gas_dense_frac_thesis``。

    Source: Hamel (1999) §7.1.1 HTW 多段气化剂进床；thesis 工程占位
    """
    bed_indices = _bed_cell_indices_from_cells(cells)
    if not bed_indices:
        return {"staged_primary_gas_bed_cells": 0.0}

    n_stage = int(getattr(cfg, "vorab_staged_primary_gas_bed_cells_thesis", 3))
    n_stage = int(np.clip(n_stage, 1, len(bed_indices)))
    weights = _resolve_staged_primary_gas_layer_weights(cfg, n_stage)

    _clear_bed_primary_gas_zu(cells, bed_indices)

    O2_total = float(cfg.O2_feed)
    H2O_total = float(cfg.H2O_feed)
    N2_total = float(cfg.N2_feed)
    staged_dense = float(getattr(cfg, "vorab_staged_primary_gas_dense_frac_thesis", 0.70))

    for k in range(n_stage):
        bi = bed_indices[k]
        cell = cells[bi]
        w = float(weights[k])
        if k == 0:
            dense_frac = resolve_bottom_gas_inlet_dense_fraction(cell, cfg)
            o2_dense_frac = resolve_bottom_o2_dense_fraction(cell, cfg)
        else:
            dense_frac = float(np.clip(staged_dense, 0.05, 0.95))
            o2_dense_frac = dense_frac
        _assign_primary_gas_zu_fraction(
            cell,
            O2_mol_s=O2_total * w,
            H2O_mol_s=H2O_total * w,
            N2_mol_s=N2_total * w,
            dense_frac=dense_frac,
            o2_dense_frac=o2_dense_frac,
        )

    return {
        "staged_primary_gas_bed_cells": float(n_stage),
        "staged_primary_gas_layer_weights": [float(w) for w in weights],
    }

def set_bottom_cell_feeds(cells: list[Cell], cfg: Any) -> None:
    """Apply primary gas and fresh solid feeds to the bottom bed cell."""
    if not cells:
        return
    bed_indices = _bed_cell_indices_from_cells(cells)
    cell = cells[bed_indices[0]] if bed_indices else cells[0]

    use_staged = bool(getattr(cfg, "vorab_staged_primary_gas_injection_thesis", False)) and len(
        bed_indices
    ) >= 2
    if use_staged:
        distribute_staged_primary_gas_to_bed_cells(cells, cfg)
    else:
        _clear_bed_primary_gas_zu(cells, bed_indices or [0])
        dense_frac = resolve_bottom_gas_inlet_dense_fraction(cell, cfg)
        o2_dense_frac = resolve_bottom_o2_dense_fraction(cell, cfg)
        _assign_primary_gas_zu_fraction(
            cell,
            O2_mol_s=float(cfg.O2_feed),
            H2O_mol_s=float(cfg.H2O_feed),
            N2_mol_s=float(cfg.N2_feed),
            dense_frac=dense_frac,
            o2_dense_frac=o2_dense_frac,
        )

    cell.m_solid_zu.fill(0.0)
    w_m = float(cfg.moisture_wt) / 100.0
    w_ash_dry = float(cfg.ash_dry_wt) / 100.0
    w_vm_dry = (float(cfg.VM_daf) / 100.0) * (1.0 - w_ash_dry)
    w_char_dry = 1.0 - w_ash_dry - w_vm_dry
    nk = int(cfg.n_age_classes)
    feed_frac = np.asarray(cell.solid.mass_fractions, dtype=np.float64)
    feed_frac = feed_frac / max(float(np.sum(feed_frac)), 1.0e-12)
    char_total = float(cfg.fuel_feed) * (1.0 - w_m) * w_char_dry
    vm_total = float(cfg.fuel_feed) * (1.0 - w_m) * w_vm_dry
    moist_total = float(cfg.fuel_feed) * w_m
    ash_total = float(cfg.fuel_feed) * (1.0 - w_m) * w_ash_dry
    cell.m_solid_zu[:, S_CHAR] = char_total * feed_frac
    cell.m_solid_zu[:, S_VM] = vm_total * feed_frac
    cell.m_solid_zu[:, S_MOISTURE] = moist_total * feed_frac
    cell.m_solid_zu[:, S_ASH] = ash_total * feed_frac

    cell.T_zu_gas = float(cfg.T_inlet)
    cell.T_zu_solid = 293.15
    cell.T_in_gas = float(cfg.T_inlet)
    cell.T_in_solid = float(cfg.T_inlet)

def apply_bottom_recycle(cells: list[Cell], cfg: Any, relax: float | None = None) -> None:
    """Map the top-cell recycle stream back onto the bottom-cell boundary."""
    if not cells:
        return
    bot = cells[0]
    if float(cfg.recirculation_frac) <= 0.0:
        bot.N_rez_d.fill(0.0)
        bot.N_rez_b.fill(0.0)
        bot.m_solid_rez.fill(0.0)
        bot.T_rez_gas = bot.T_zu_gas
        bot.T_rez_solid = bot.T_zu_solid
        return

    top = cells[-1]
    frac = float(cfg.recirculation_frac)
    if bool(cfg.recycle_gas):
        new_N_rez_d = np.maximum(top.N_d, 0.0) * frac
        new_N_rez_b = np.maximum(top.N_b, 0.0) * frac
    else:
        new_N_rez_d = np.zeros_like(bot.N_rez_d)
        new_N_rez_b = np.zeros_like(bot.N_rez_b)
    new_m_rez = recycled_solid_stream_from_cell(top, frac)
    new_T_rez = float(top.T)

    if relax is None:
        bot.N_rez_d[:] = new_N_rez_d
        bot.N_rez_b[:] = new_N_rez_b
        bot.m_solid_rez[:] = new_m_rez
        bot.T_rez_gas = new_T_rez
        bot.T_rez_solid = new_T_rez
        return

    omega = float(np.clip(relax, 0.0, 1.0))
    bot.N_rez_d[:] = omega * new_N_rez_d + (1.0 - omega) * bot.N_rez_d
    bot.N_rez_b[:] = omega * new_N_rez_b + (1.0 - omega) * bot.N_rez_b
    bot.m_solid_rez[:] = omega * new_m_rez + (1.0 - omega) * bot.m_solid_rez
    bot.T_rez_gas = omega * new_T_rez + (1.0 - omega) * bot.T_rez_gas
    bot.T_rez_solid = omega * new_T_rez + (1.0 - omega) * bot.T_rez_solid

def apply_bottom_recycle_from_return_leg(
    bed_cells: list[Cell],
    return_leg_cell: Cell,
    cfg: Any,
    relax: float | None = None,
    *,
    target_bed_index: int = 0,
) -> None:
    """Map an explicit return-leg solid stream onto a target bed cell."""
    if not bed_cells:
        return
    bed_idx = int(np.clip(target_bed_index, 0, len(bed_cells) - 1))
    bot = bed_cells[bed_idx]
    new_N_rez_d = np.zeros_like(bot.N_rez_d)
    new_N_rez_b = np.zeros_like(bot.N_rez_b)
    new_m_rez = np.maximum(return_leg_cell._solid_downflow_rates(), 0.0)
    new_T_rez = float(return_leg_cell.T)

    if relax is None:
        bot.N_rez_d[:] = new_N_rez_d
        bot.N_rez_b[:] = new_N_rez_b
        bot.m_solid_rez[:] = new_m_rez
        bot.T_rez_gas = bot.T_zu_gas
        bot.T_rez_solid = new_T_rez
        return

    omega = float(np.clip(relax, 0.0, 1.0))
    bot.N_rez_d[:] = omega * new_N_rez_d + (1.0 - omega) * bot.N_rez_d
    bot.N_rez_b[:] = omega * new_N_rez_b + (1.0 - omega) * bot.N_rez_b
    bot.m_solid_rez[:] = omega * new_m_rez + (1.0 - omega) * bot.m_solid_rez
    bot.T_rez_gas = bot.T_zu_gas
    bot.T_rez_solid = omega * new_T_rez + (1.0 - omega) * bot.T_rez_solid

