"""Thesis-mode connectivity and routing for Hamel-style reactor blocks."""

from __future__ import annotations

from typing import Any
from dataclasses import asdict, dataclass

import numpy as np

from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX

_F_W = 0.25
_ZETA_W = 0.40


@dataclass(frozen=True)
class ConnectivityBlock:
    block_id: str
    kind: str
    axial_index: int | None = None
    xi_reactor: float | None = None
    solver_coupling: str = "explicit"


@dataclass(frozen=True)
class ConnectivityEdge:
    from_block: str
    to_block: str
    kind: str
    direction: str
    solver_coupling: str
    phase_scope: str | None = None
    mechanism: str | None = None


def recycled_solid_stream(m_solid_top: np.ndarray, frac: float) -> np.ndarray:
    """Recycle stream for the current thesis-aligned bed model."""
    m_recycle = np.zeros_like(m_solid_top)
    m_recycle[:, S_CHAR] = np.maximum(m_solid_top[:, S_CHAR], 0.0) * float(frac)
    m_recycle[:, S_ASH] = np.maximum(m_solid_top[:, S_ASH], 0.0) * float(frac)
    return m_recycle


def propagated_solid_stream(
    m_solid: np.ndarray,
    frac: float = 1.0,
    *,
    include_reactive: bool = False,
) -> np.ndarray:
    """Axial solid propagation for the current thesis-aligned bed model."""
    out = np.zeros_like(m_solid)
    out[:, S_CHAR] = np.maximum(m_solid[:, S_CHAR], 0.0) * float(frac)
    out[:, S_ASH] = np.maximum(m_solid[:, S_ASH], 0.0) * float(frac)
    if include_reactive:
        out[:, S_VM] = np.maximum(m_solid[:, S_VM], 0.0) * float(frac)
        out[:, S_MOISTURE] = np.maximum(m_solid[:, S_MOISTURE], 0.0) * float(frac)
    return out


def propagated_solid_stream_from_cell(
    cell: Cell,
    frac: float = 1.0,
    *,
    include_reactive: bool = False,
) -> np.ndarray:
    """Project the cell's outgoing solid transport stream onto routed components."""
    upflow = cell._solid_upflow_rates()
    routed = propagated_solid_stream(upflow, frac=frac, include_reactive=False)
    if include_reactive:
        if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
            # Keep units consistent ([kg/s]): propagate fresh-feed reactive inlet support
            # rather than holdup inventory ([kg]).
            reactive_source = np.maximum(cell.m_solid_zu + cell.m_solid_in + cell.m_solid_rez, 0.0)
        else:
            reactive_source = np.maximum(upflow, 0.0)
        routed[:, S_VM] = reactive_source[:, S_VM] * float(frac)
        routed[:, S_MOISTURE] = reactive_source[:, S_MOISTURE] * float(frac)
    return routed


def recycled_solid_stream_from_cell(cell: Cell, frac: float) -> np.ndarray:
    """Recycle only the routed fraction of the cell's outgoing solid stream."""
    return recycled_solid_stream(cell._solid_upflow_rates(), frac)


def project_zero_source_solid_components(cell: Cell) -> None:
    """Project VM/moisture state onto the physically allowed inlet support."""
    vm_in = float(
        np.sum(
            np.maximum(
                cell.m_solid_zu[:, S_VM]
                + cell.m_solid_in[:, S_VM]
                + cell.m_solid_rez[:, S_VM],
                0.0,
            )
        )
    )
    if vm_in <= 1e-12:
        cell.m_solid[:, S_VM] = 0.0

    moist_in = float(
        np.sum(
            np.maximum(
                cell.m_solid_zu[:, S_MOISTURE]
                + cell.m_solid_in[:, S_MOISTURE]
                + cell.m_solid_rez[:, S_MOISTURE],
                0.0,
            )
        )
    )
    if moist_in <= 1e-12:
        cell.m_solid[:, S_MOISTURE] = 0.0


def cell_total_solid_holdup(cell: Cell) -> float:
    """Estimate total solid holdup [kg] from local hydrodynamics and geometry.

    Ref: Hamel (1999) Eq. 2.6 coupling interpretation; user-provided bed discretization
    for ``M_solid = A * dh * (1-eps_b) * rho_p * (1-eps_d)``.
    """
    area = float(np.pi * cell.geo.D_bed**2 / 4.0)
    return float(
        max(area, 0.0)
        * max(float(cell.geo.dh), 0.0)
        * max(1.0 - float(cell.eps_b), 0.0)
        * max(float(cell.solid.rho_s), 0.0)
        * max(1.0 - float(cell.eps_d_voidage), 0.0)
    )


def bed_upflow_coefficient(cell: Cell, *, is_top_bed_cell: bool) -> float:
    """Hamel-style wake carrying coefficient ``K_auf`` for a bed cell [1/s]."""
    denom = max((1.0 - float(cell.eps_b)) * float(cell.geo.dh), 1e-12)
    k_auf = _F_W * max(float(cell.eps_b), 0.0) * max(float(cell.u_b), 0.0) / denom
    if is_top_bed_cell:
        k_auf *= _ZETA_W
    return float(max(k_auf, 0.0))


def update_bed_solid_transport_coefficients(
    bed_cells: list[Cell],
    cfg: Any,
    *,
    top_downflow_total: float = 0.0,
) -> None:
    """Update frozen bed ``K_auf`` / ``K_ab`` transport coefficients.

    The current implementation follows the thesis-aligned bed formulas supplied by the
    user:
    - ``K_auf = f_w * eps_b * u_b / ((1-eps_b) * dh)``
    - top bed cell additionally applies the surface ejection factor ``zeta_w``
    - ``K_ab`` is reconstructed from top-down continuity using total solid holdup.

    This helper is intentionally scoped to bed cells only. Freeboard trajectory-based
    coefficients will be added separately.
    """
    if not bed_cells:
        return

    holdups = np.array([cell_total_solid_holdup(cell) for cell in bed_cells], dtype=np.float64)
    k_auf_vals = np.array(
        [bed_upflow_coefficient(cell, is_top_bed_cell=(i == len(bed_cells) - 1)) for i, cell in enumerate(bed_cells)],
        dtype=np.float64,
    )
    m_auf_totals = k_auf_vals * holdups

    m_ab_from_above = float(max(top_downflow_total, 0.0))
    for i in range(len(bed_cells) - 1, -1, -1):
        cell = bed_cells[i]
        up_in = float(m_auf_totals[i - 1]) if i > 0 else 0.0
        # Keep K_ab reconstruction consistent with the solid Eq.2.6 residual terms:
        # include side-feed/recycle and any explicit solid inlet source.
        ext_in = float(np.sum(np.maximum(cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in, 0.0)))
        reaction_total = float(np.sum(np.asarray(cell.R_solid, dtype=np.float64)))
        m_ab_total = max(up_in - float(m_auf_totals[i]) + m_ab_from_above + ext_in + reaction_total, 0.0)
        k_ab = m_ab_total / max(float(holdups[i]), 1e-12) if holdups[i] > 1e-12 else 0.0

        cell.K_solid_auf.fill(0.0)
        cell.K_solid_ab.fill(0.0)
        # Thesis holdup transport: only char/ash are transported across cells.
        cell.K_solid_auf[:, S_CHAR] = float(k_auf_vals[i])
        cell.K_solid_auf[:, S_ASH] = float(k_auf_vals[i])
        cell.K_solid_ab[:, S_CHAR] = float(k_ab)
        cell.K_solid_ab[:, S_ASH] = float(k_ab)
        m_ab_from_above = m_ab_total


def update_bed_solid_transport_inflows(bed_cells: list[Cell]) -> None:
    """Map neighboring bed-cell holdup onto thesis-style axial solid inflow terms."""
    update_bed_solid_transport_inflows_from_above(bed_cells)


def update_bed_solid_transport_inflows_from_above(
    bed_cells: list[Cell],
    top_above_cell: Cell | None = None,
) -> None:
    """Map neighboring bed/freeboard holdup onto thesis-style bed inflow terms."""
    if not bed_cells:
        return
    for cell in bed_cells:
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)

    for i, cell in enumerate(bed_cells):
        if i > 0:
            below = bed_cells[i - 1]
            cell.m_solid_auf_in[:, :] = below.K_solid_auf * np.maximum(below.m_solid, 0.0)
        if i + 1 < len(bed_cells):
            above = bed_cells[i + 1]
            cell.m_solid_ab_in[:, :] = above.K_solid_ab * np.maximum(above.m_solid, 0.0)
        elif top_above_cell is not None:
            cell.m_solid_ab_in[:, :] = np.maximum(top_above_cell._solid_downflow_rates(), 0.0)
        # Guardrail: axial transport channels carry only char/ash in thesis holdup mode.
        cell.m_solid_auf_in[:, S_VM] = 0.0
        cell.m_solid_auf_in[:, S_MOISTURE] = 0.0
        cell.m_solid_ab_in[:, S_VM] = 0.0
        cell.m_solid_ab_in[:, S_MOISTURE] = 0.0


def propagate_explicit_freeboard_chain(
    bed_cells: list[Cell],
    freeboard_cells: list[Cell],
) -> None:
    """Propagate axial gas/solid inlets from bed top through explicit freeboard cells."""
    if not bed_cells or not freeboard_cells:
        return

    prev = bed_cells[-1]
    for cell in freeboard_cells:
        _set_explicit_freeboard_inlet_from_prev(prev, cell)
        prev = cell


def _set_explicit_freeboard_inlet_from_prev(prev: Cell, cell: Cell) -> None:
    """Refresh one explicit freeboard cell inlet from its upstream source cell."""
    cell.N_b_in[:] = np.maximum(prev.N_b, 0.0)
    cell.N_d_in[:] = np.maximum(prev.N_d, 0.0)
    cell.T_in_gas = float(prev.T)
    cell.m_solid_in.fill(0.0)
    if str(getattr(cell, "solid_state_model", "legacy_stream")) != "freeboard_closure":
        cell.m_solid_in[:] = propagated_solid_stream_from_cell(prev, include_reactive=False)
    active_mask = getattr(cell, "_freeboard_active_char_ash_mask", None)
    if isinstance(active_mask, np.ndarray) and active_mask.shape[0] == cell.solid.n_size_classes:
        mask = np.maximum(active_mask, 0.0).reshape(-1)
        cell.m_solid_in[:, S_CHAR] *= mask
        cell.m_solid_in[:, S_ASH] *= mask
    if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        vm_col = np.maximum(cell.m_solid_in[:, S_VM], 0.0)
        moist_col = np.maximum(cell.m_solid_in[:, S_MOISTURE], 0.0)
        cell.m_solid_in.fill(0.0)
        cell.m_solid_in[:, S_VM] = vm_col
        cell.m_solid_in[:, S_MOISTURE] = moist_col
    cell.T_in_solid = float(prev.T)


def update_freeboard_solid_transport_inflows(
    freeboard_cells: list[Cell],
    *,
    bottom_below_cell: Cell,
) -> None:
    """Project explicit freeboard K_auf/K_ab into neighboring inflow terms.

    Hamel-aligned ``freeboard_closure`` cells keep closure hold-up and closure
    outflow coefficients, but do not participate in explicit Eulerian solid
    inflow chaining inside the NR boundary path.
    """
    if not freeboard_cells:
        return
    for i, cell in enumerate(freeboard_cells):
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
        if i + 1 < len(freeboard_cells):
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
        # Guardrail: freeboard axial transport carries only char/ash.
        cell.m_solid_auf_in[:, S_VM] = 0.0
        cell.m_solid_auf_in[:, S_MOISTURE] = 0.0
        cell.m_solid_ab_in[:, S_VM] = 0.0
        cell.m_solid_ab_in[:, S_MOISTURE] = 0.0
        _apply_freeboard_active_class_mask_to_transport_inflows(cell)


def _refresh_freeboard_solid_transport_inflows_for_indices(
    freeboard_cells: list[Cell],
    *,
    bottom_below_cell: Cell,
    indices: set[int],
) -> None:
    """Refresh only selected explicit freeboard transport inflow rows."""
    if not freeboard_cells:
        return
    n = len(freeboard_cells)
    for i in sorted(int(idx) for idx in indices if 0 <= int(idx) < n):
        cell = freeboard_cells[i]
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
        if i + 1 < n:
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
        cell.m_solid_auf_in[:, S_VM] = 0.0
        cell.m_solid_auf_in[:, S_MOISTURE] = 0.0
        cell.m_solid_ab_in[:, S_VM] = 0.0
        cell.m_solid_ab_in[:, S_MOISTURE] = 0.0
        _apply_freeboard_active_class_mask_to_transport_inflows(cell)


def _apply_freeboard_active_class_mask_to_transport_inflows(cell: Cell) -> None:
    """Mask char/ash transport inflows by closure-active freeboard size classes."""
    active_mask = getattr(cell, "_freeboard_active_char_ash_mask", None)
    if not (isinstance(active_mask, np.ndarray) and active_mask.shape[0] == cell.solid.n_size_classes):
        return
    mask = np.maximum(active_mask, 0.0).reshape(-1)
    cell.m_solid_auf_in[:, S_CHAR] *= mask
    cell.m_solid_auf_in[:, S_ASH] *= mask
    cell.m_solid_ab_in[:, S_CHAR] *= mask
    cell.m_solid_ab_in[:, S_ASH] *= mask


def _sideblock_capture_efficiency(cfg: Any, comp_idx: int) -> float:
    """Cyclone capture efficiency for routed solid components."""
    if comp_idx == S_CHAR:
        return float(np.clip(getattr(cfg, "freeboard_cyclone_capture_char_frac", 0.0), 0.0, 1.0))
    if comp_idx == S_ASH:
        return float(np.clip(getattr(cfg, "freeboard_cyclone_capture_ash_frac", 0.0), 0.0, 1.0))
    return 0.0


def _estimate_cyclone_residence_time_s(cyclone_cell: Cell) -> float:
    """Estimate cyclone residence time from gas holdup volume and volumetric throughput."""
    gas_total = float(np.sum(np.maximum(cyclone_cell.N_b + cyclone_cell.N_d, 0.0)))
    if gas_total <= 1e-12:
        gas_total = float(np.sum(np.maximum(cyclone_cell.N_b_in + cyclone_cell.N_d_in, 0.0)))
    T_ref = float(cyclone_cell.T if gas_total > 1e-12 and cyclone_cell.T > 0.0 else cyclone_cell.T_in_gas)
    if gas_total <= 1e-12 or T_ref <= 0.0:
        return 1.0
    gas_vol_flow = gas_total * Rg * T_ref / max(float(cyclone_cell.P), 1e-12)
    V_cyc = max(float(cyclone_cell.V_d), float(np.pi * cyclone_cell.geo.D_bed**2 * cyclone_cell.geo.dh / 4.0), 1e-12)
    dt = V_cyc / max(gas_vol_flow, 1e-12)
    return float(np.clip(dt, 0.1, 1.0))


def update_side_block_solid_transport_coefficients(
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    cfg: Any,
) -> None:
    """Update thesis-style side-block transport coefficients ``K = 1/dt``.

    Ref: Hamel (1999) Eq. 2.6 coupling interpretation from user-provided thesis notes.
    Cyclone uses a short gas-throughput residence time; return leg uses dense-leg holdup
    capacity divided by the current cyclone underflow throughput.
    """
    K_cyc = 1.0 / max(_estimate_cyclone_residence_time_s(cyclone_cell), 1e-12)
    cyclone_cell.K_solid_auf.fill(0.0)
    cyclone_cell.K_solid_ab.fill(0.0)
    for comp_idx in (S_CHAR, S_ASH):
        eta = _sideblock_capture_efficiency(cfg, comp_idx)
        cyclone_cell.K_solid_auf[:, comp_idx] = (1.0 - eta) * K_cyc
        cyclone_cell.K_solid_ab[:, comp_idx] = eta * K_cyc

    leg_capacity = float(
        max(return_leg_cell.V_d, np.pi * return_leg_cell.geo.D_bed**2 * return_leg_cell.geo.dh / 4.0)
        * max(1.0 - float(return_leg_cell.solid.eps_mf), 0.0)
        * max(float(return_leg_cell.solid.rho_s), 0.0)
    )
    leg_throughput = float(np.sum(np.maximum(cyclone_cell._solid_downflow_rates(), 0.0)))
    K_leg = leg_throughput / max(leg_capacity, 1e-12) if leg_capacity > 1e-12 else 0.0
    return_leg_cell.K_solid_auf.fill(0.0)
    return_leg_cell.K_solid_ab.fill(0.0)
    for comp_idx in (S_CHAR, S_ASH):
        return_leg_cell.K_solid_ab[:, comp_idx] = K_leg


def seed_side_block_holdup_from_inflows(cell: Cell) -> None:
    """Seed a side-block holdup state from its current thesis inflows and frozen ``K``."""
    if str(cell.solid_state_model) not in {"holdup_transport", "freeboard_closure"}:
        return
    if float(np.sum(np.maximum(cell.m_solid, 0.0))) > 1e-12:
        return
    inflow = np.maximum(cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in + cell.m_solid_auf_in + cell.m_solid_ab_in, 0.0)
    K_total = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    seeded = np.zeros_like(cell.m_solid)
    mask = K_total > 1e-12
    seeded[mask] = inflow[mask] / K_total[mask]
    cell.m_solid[:, :] = seeded


def apply_frozen_bed_solid_transport_coefficients(bed_cells: list[Cell]) -> bool:
    """Replay frozen Vorabrechnung hydrodynamics/transport coefficients if available."""
    if not bed_cells:
        return False
    if not all(bool(cell._freeze_vorabrechnung_inner_nr) and bool(cell._vorab_solid_transport_cache_valid) for cell in bed_cells):
        return False
    for cell in bed_cells:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    return True


def apply_frozen_side_block_solid_transport_coefficients(
    cyclone_cell: Cell,
    return_leg_cell: Cell,
) -> bool:
    """Replay frozen side-block transport coefficients inside inner NR."""
    side_cells = (cyclone_cell, return_leg_cell)
    if not all(bool(cell._freeze_vorabrechnung_inner_nr) and bool(cell._vorab_solid_transport_cache_valid) for cell in side_cells):
        return False
    for cell in side_cells:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    return True


def set_bottom_cell_feeds(cells: list[Cell], cfg: Any) -> None:
    """Apply primary gas and fresh solid feeds to the bottom bed cell."""
    if not cells:
        return
    cell = cells[0]
    idx = GAS_SPECIES_INDEX
    cell.N_zu_d.fill(0.0)
    cell.N_zu_b.fill(0.0)
    cell.m_solid_zu.fill(0.0)

    dense_frac = float(np.clip(cfg.gas_inlet_dense_frac, 0.0, 1.0))
    bubble_frac = 1.0 - dense_frac
    cell.N_zu_d[idx["O2"]] = float(cfg.O2_feed) * dense_frac
    cell.N_zu_d[idx["H2O"]] = float(cfg.H2O_feed) * dense_frac
    cell.N_zu_d[idx["N2"]] = float(cfg.N2_feed) * dense_frac
    cell.N_zu_b[idx["O2"]] = float(cfg.O2_feed) * bubble_frac
    cell.N_zu_b[idx["H2O"]] = float(cfg.H2O_feed) * bubble_frac
    cell.N_zu_b[idx["N2"]] = float(cfg.N2_feed) * bubble_frac

    w_m = float(cfg.moisture_wt) / 100.0
    w_ash_dry = float(cfg.ash_dry_wt) / 100.0
    w_vm_dry = (float(cfg.VM_daf) / 100.0) * (1.0 - w_ash_dry)
    w_char_dry = 1.0 - w_ash_dry - w_vm_dry
    nk = int(cfg.n_age_classes)
    cell.m_solid_zu[:, S_CHAR] = float(cfg.fuel_feed) * (1.0 - w_m) * w_char_dry / nk
    cell.m_solid_zu[:, S_VM] = float(cfg.fuel_feed) * (1.0 - w_m) * w_vm_dry / nk
    cell.m_solid_zu[:, S_MOISTURE] = float(cfg.fuel_feed) * w_m / nk
    cell.m_solid_zu[:, S_ASH] = float(cfg.fuel_feed) * (1.0 - w_m) * w_ash_dry / nk

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
) -> None:
    """Map the explicit return-leg solid stream back onto the bottom bed cell."""
    if not bed_cells:
        return
    bot = bed_cells[0]
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


def propagate_upstream(cells: list[Cell], cfg: Any, i: int) -> None:
    """Propagate axial gas and suspension-phase solid coupling into cell ``i``."""
    if i == 0:
        return
    prev, curr = cells[i - 1], cells[i]
    curr.N_b_in[:] = prev.N_b
    curr.N_d_in[:] = prev.N_d
    curr.T_in_gas = prev.T
    xi_curr = float(curr.geo.h_center / max(float(cfg.H_bed), 1e-12))
    reactive_inlet_allowed = bool(cfg.allow_reactive_solid_propagation) and (
        xi_curr <= float(np.clip(cfg.reactive_solid_cutoff_xi, 0.0, 1.0))
    )
    if i < len(cells) - 1:
        above = cells[i + 1]
        lower_frac = float(np.clip(cfg.solid_lower_inlet_frac, 0.0, 1.0))
        upper_frac = 1.0 - lower_frac
        above_solid = propagated_solid_stream_from_cell(above, include_reactive=reactive_inlet_allowed)
        lower_solid = propagated_solid_stream_from_cell(prev, include_reactive=reactive_inlet_allowed)
        above_empty = float(np.sum(above_solid)) <= 1e-12
        above_reactive_empty = float(
            np.sum(np.maximum(above_solid[:, S_VM] + above_solid[:, S_MOISTURE], 0.0))
        ) <= 1e-12
        if above_empty:
            above_solid = lower_solid
            above_T = prev.T
        elif reactive_inlet_allowed and above_reactive_empty:
            # Keep counter-current char/ash routing from above, but pull reactive fresh-feed
            # support from lower cell so lower-zone VM/moisture release can continue.
            above_solid[:, S_VM] = lower_solid[:, S_VM]
            above_solid[:, S_MOISTURE] = lower_solid[:, S_MOISTURE]
            above_T = above.T
        else:
            above_T = above.T
        curr.m_solid_in[:] = (
            upper_frac * above_solid
            + lower_frac * lower_solid
        )
        curr.T_in_solid = upper_frac * above_T + lower_frac * prev.T
    else:
        curr.m_solid_in[:] = propagated_solid_stream_from_cell(
            prev,
            frac=float(np.clip(cfg.top_solid_inlet_frac, 0.0, 1.0)),
            include_reactive=reactive_inlet_allowed,
        )
        curr.T_in_solid = prev.T

    # Thesis holdup transport path uses m_solid_auf_in/m_solid_ab_in for axial transport
    # terms. Keep m_solid_in as a dedicated reactive-feed carrier only (VM/moisture) to
    # avoid char/ash double-counting in Eq.2.6 when both paths are active.
    if str(getattr(curr, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        vm_col = np.maximum(curr.m_solid_in[:, S_VM], 0.0)
        moist_col = np.maximum(curr.m_solid_in[:, S_MOISTURE], 0.0)
        curr.m_solid_in.fill(0.0)
        curr.m_solid_in[:, S_VM] = vm_col
        curr.m_solid_in[:, S_MOISTURE] = moist_col


def route_auxiliary_side_blocks(
    source_cells: list[Cell],
    cyclone_cell: Cell,
    return_leg_cell: Cell,
    cfg: Any,
) -> None:
    """Route the current top source outlet into cyclone, then cyclone underflow into return leg."""
    if not source_cells:
        return
    top = source_cells[-1]
    cyclone_cell.N_b_in[:] = top.N_b
    cyclone_cell.N_d_in[:] = top.N_d
    cyclone_cell.T_in_gas = float(top.T)
    cyclone_cell.m_solid_in.fill(0.0)
    cyclone_cell.m_solid_auf_in[:, :] = propagated_solid_stream_from_cell(top, include_reactive=False)
    cyclone_cell.m_solid_ab_in.fill(0.0)
    cyclone_cell.T_in_solid = float(top.T)
    cyclone_cell.N_zu_b.fill(0.0)
    cyclone_cell.N_zu_d.fill(0.0)
    cyclone_cell.N_rez_b.fill(0.0)
    cyclone_cell.N_rez_d.fill(0.0)
    cyclone_cell.m_solid_zu.fill(0.0)
    cyclone_cell.m_solid_rez.fill(0.0)
    cyclone_cell.T_zu_gas = float(top.T)
    cyclone_cell.T_rez_gas = float(top.T)
    cyclone_cell.T_zu_solid = float(top.T)
    cyclone_cell.T_rez_solid = float(top.T)

    return_leg_cell.N_b_in.fill(0.0)
    return_leg_cell.N_d_in.fill(0.0)
    return_leg_cell.T_in_gas = float(cyclone_cell.T)
    return_leg_cell.m_solid_in.fill(0.0)
    return_leg_cell.m_solid_auf_in.fill(0.0)
    return_leg_cell.m_solid_ab_in.fill(0.0)
    return_leg_cell.T_in_solid = float(cyclone_cell.T)
    return_leg_cell.N_zu_b.fill(0.0)
    return_leg_cell.N_zu_d.fill(0.0)
    return_leg_cell.N_rez_b.fill(0.0)
    return_leg_cell.N_rez_d.fill(0.0)
    return_leg_cell.m_solid_zu.fill(0.0)
    return_leg_cell.m_solid_rez.fill(0.0)
    return_leg_cell.T_zu_gas = float(cyclone_cell.T)
    return_leg_cell.T_rez_gas = float(cyclone_cell.T)
    return_leg_cell.T_zu_solid = float(cyclone_cell.T)
    return_leg_cell.T_rez_solid = float(cyclone_cell.T)


def _refresh_bed_transport_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell] | None = None,
) -> None:
    """Refresh bed-side transport coefficients and axial inflow terms."""
    top_downflow_total = 0.0
    if freeboard_cells:
        top_downflow_total = float(np.sum(np.maximum(freeboard_cells[0]._solid_downflow_rates(), 0.0)))
    if not apply_frozen_bed_solid_transport_coefficients(cells):
        update_bed_solid_transport_coefficients(cells, cfg, top_downflow_total=top_downflow_total)
    update_bed_solid_transport_inflows_from_above(cells, freeboard_cells[0] if freeboard_cells else None)


def _refresh_side_block_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell],
    cyclone_cell: Cell,
    return_leg_cell: Cell,
) -> None:
    """Refresh explicit cyclone/return-leg routing and bottom recycle."""
    route_auxiliary_side_blocks(freeboard_cells if freeboard_cells else cells, cyclone_cell, return_leg_cell, cfg)
    if float(np.sum(np.maximum(cyclone_cell.N_d + cyclone_cell.N_b, 0.0))) <= 1e-12:
        cyclone_cell.N_d[:] = np.maximum(cyclone_cell.N_d_in, 0.0)
        cyclone_cell.N_b[:] = np.maximum(cyclone_cell.N_b_in, 0.0)
        cyclone_cell.T = float(cyclone_cell.T_in_gas)
    side_transport_frozen = apply_frozen_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell)
    if not side_transport_frozen:
        update_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell, cfg)
    seed_side_block_holdup_from_inflows(cyclone_cell)
    # Do not overwrite the side-block temperature state during NR boundary updates.
    # Cyclone / return-leg T are packed in the Newton vector and must remain free
    # to move under their own energy residuals; forcing T := T_in_* here creates
    # dead temperature DOFs and Jacobian zero columns.
    if float(np.sum(np.maximum(cyclone_cell.m_solid, 0.0))) <= 1e-12 and float(
        np.sum(np.maximum(cyclone_cell.N_d + cyclone_cell.N_b, 0.0))
    ) <= 1e-12:
        cyclone_cell.T = float(cyclone_cell.T_in_solid)
    return_leg_cell.m_solid_ab_in[:, :] = np.maximum(cyclone_cell._solid_downflow_rates(), 0.0)
    return_leg_cell.T_in_solid = float(cyclone_cell.T)
    if not side_transport_frozen:
        update_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell, cfg)
    seed_side_block_holdup_from_inflows(return_leg_cell)
    apply_bottom_recycle_from_return_leg(cells, return_leg_cell, cfg, relax=None)
    project_zero_source_solid_components(cyclone_cell)
    project_zero_source_solid_components(return_leg_cell)


def apply_all_nr_boundary_data(
    cells: list[Cell],
    cfg: Any,
    *,
    freeboard_cells: list[Cell] | None = None,
    cyclone_cell: Cell | None = None,
    return_leg_cell: Cell | None = None,
) -> None:
    """Apply the full thesis-mode axial routing / recycle boundary update."""
    if not cells:
        return
    freeboard_cells = freeboard_cells or []
    set_bottom_cell_feeds(cells, cfg)
    cells[-1].m_solid_in[:] = 0.0
    for i in range(len(cells)):
        propagate_upstream(cells, cfg, i)
    if freeboard_cells:
        propagate_explicit_freeboard_chain(cells, freeboard_cells)
    _refresh_bed_transport_boundary_data(cells, cfg, freeboard_cells=freeboard_cells)
    if freeboard_cells:
        update_freeboard_solid_transport_inflows(freeboard_cells, bottom_below_cell=cells[-1])
    if cyclone_cell is not None and return_leg_cell is not None:
        _refresh_side_block_boundary_data(
            cells,
            cfg,
            freeboard_cells=freeboard_cells,
            cyclone_cell=cyclone_cell,
            return_leg_cell=return_leg_cell,
        )
    else:
        apply_bottom_recycle(cells, cfg, relax=None)
    for cell in cells:
        project_zero_source_solid_components(cell)
    for cell in freeboard_cells:
        project_zero_source_solid_components(cell)


def apply_local_nr_boundary_data(
    cells: list[Cell],
    cfg: Any,
    changed_cell_idx: int,
    *,
    freeboard_cells: list[Cell] | None = None,
    cyclone_cell: Cell | None = None,
    return_leg_cell: Cell | None = None,
) -> None:
    """Apply only the BC updates touched by a single perturbed cell."""
    n_cells = len(cells)
    if n_cells == 0:
        return
    if freeboard_cells or cyclone_cell is not None or return_leg_cell is not None:
        freeboard_cells = freeboard_cells or []
        idx = int(changed_cell_idx)
        n_bed = len(cells)
        n_fb = len(freeboard_cells)
        touched_bed: set[int] = set()
        touched_fb: set[int] = set()
        touched_side = False

        set_bottom_cell_feeds(cells, cfg)
        cells[-1].m_solid_in[:] = 0.0

        if idx < n_bed:
            touched_bed.add(idx)
            if idx > 0:
                touched_bed.add(idx - 1)
            if idx + 1 < n_bed:
                propagate_upstream(cells, cfg, idx + 1)
                touched_bed.add(idx + 1)
            if idx == n_bed - 1:
                if n_fb > 0:
                    _set_explicit_freeboard_inlet_from_prev(cells[-1], freeboard_cells[0])
                    _refresh_freeboard_solid_transport_inflows_for_indices(
                        freeboard_cells,
                        bottom_below_cell=cells[-1],
                        indices={0},
                    )
                    touched_fb.add(0)
                elif cyclone_cell is not None and return_leg_cell is not None:
                    touched_side = True
            _refresh_bed_transport_boundary_data(cells, cfg, freeboard_cells=freeboard_cells)
        elif idx < n_bed + n_fb:
            fb_idx = idx - n_bed
            touched_fb.add(fb_idx)
            if fb_idx > 0:
                _refresh_freeboard_solid_transport_inflows_for_indices(
                    freeboard_cells,
                    bottom_below_cell=cells[-1],
                    indices={fb_idx - 1},
                )
                touched_fb.add(fb_idx - 1)
            if fb_idx + 1 < n_fb:
                _set_explicit_freeboard_inlet_from_prev(freeboard_cells[fb_idx], freeboard_cells[fb_idx + 1])
                _refresh_freeboard_solid_transport_inflows_for_indices(
                    freeboard_cells,
                    bottom_below_cell=cells[-1],
                    indices={fb_idx + 1},
                )
                touched_fb.add(fb_idx + 1)
            if fb_idx == 0:
                _refresh_bed_transport_boundary_data(cells, cfg, freeboard_cells=freeboard_cells)
                touched_bed.update(range(n_bed))
            if fb_idx == n_fb - 1 and cyclone_cell is not None and return_leg_cell is not None:
                touched_side = True
        elif cyclone_cell is not None and idx == n_bed + n_fb:
            touched_side = True
        elif return_leg_cell is not None and idx == n_bed + n_fb + 1:
            apply_bottom_recycle_from_return_leg(cells, return_leg_cell, cfg, relax=None)
            touched_bed.add(0)
            if n_bed > 1:
                touched_bed.add(1)
            project_zero_source_solid_components(return_leg_cell)

        if n_fb > 0 and not touched_fb and idx < n_bed and idx != n_bed - 1:
            # Keep freeboard transport coefficients/inflows coherent only when explicit
            # freeboard cells are already active in the solver graph.
            _refresh_freeboard_solid_transport_inflows_for_indices(
                freeboard_cells,
                bottom_below_cell=cells[-1],
                indices=set(),
            )

        if touched_side and cyclone_cell is not None and return_leg_cell is not None:
            _refresh_side_block_boundary_data(
                cells,
                cfg,
                freeboard_cells=freeboard_cells,
                cyclone_cell=cyclone_cell,
                return_leg_cell=return_leg_cell,
            )
            touched_bed.add(0)
            if n_bed > 1:
                touched_bed.add(1)
            if n_fb > 0:
                touched_fb.add(n_fb - 1)

        for bed_idx in sorted(touched_bed):
            if 0 <= bed_idx < n_bed:
                project_zero_source_solid_components(cells[bed_idx])
        for fb_idx in sorted(touched_fb):
            if 0 <= fb_idx < n_fb:
                project_zero_source_solid_components(freeboard_cells[fb_idx])
        return

    set_bottom_cell_feeds(cells, cfg)
    cells[-1].m_solid_in[:] = 0.0
    apply_bottom_recycle(cells, cfg, relax=None)
    project_zero_source_solid_components(cells[0])

    touched = {int(changed_cell_idx)}
    if int(changed_cell_idx) - 1 > 0:
        touched.add(int(changed_cell_idx) - 1)
    if int(changed_cell_idx) + 1 < n_cells:
        touched.add(int(changed_cell_idx) + 1)
    if int(changed_cell_idx) == n_cells - 1:
        touched.add(0)

    for idx in sorted(touched):
        if idx > 0:
            propagate_upstream(cells, cfg, idx)
        project_zero_source_solid_components(cells[idx])
    if not apply_frozen_bed_solid_transport_coefficients(cells):
        update_bed_solid_transport_coefficients(cells, cfg)
    update_bed_solid_transport_inflows(cells)


def build_thesis_connectivity_topology(
    *,
    n_bed_cells: int,
    n_freeboard_cells: int,
    H_bed: float,
    H_freeboard: float,
    recycle_gas: bool,
    recirculation_frac: float,
    cyclone_capture_char_frac: float,
    cyclone_capture_ash_frac: float,
    secondary_injection_applied: bool,
    secondary_injection_segment: int | None,
    explicit_side_blocks: bool = False,
) -> dict:
    """Build the minimal Hamel-style connectivity graph for thesis mode."""
    total_height = max(float(H_bed) + max(float(H_freeboard), 0.0), 1e-12)
    freeboard_active = int(max(n_freeboard_cells, 0)) > 0 and float(H_freeboard) > 0.0

    blocks: list[ConnectivityBlock] = [
        ConnectivityBlock("primary_gas_feed", "side_feed", solver_coupling="boundary_condition"),
        ConnectivityBlock("fuel_feed", "side_feed", solver_coupling="boundary_condition"),
    ]
    edges: list[ConnectivityEdge] = [
        ConnectivityEdge("primary_gas_feed", "bed_cell_0", "side_feed", "into_reactor", "boundary_condition"),
        ConnectivityEdge("fuel_feed", "bed_cell_0", "side_feed", "into_reactor", "boundary_condition"),
    ]

    for i in range(int(max(n_bed_cells, 0))):
        xi = ((i + 0.5) * float(H_bed) / max(float(n_bed_cells), 1.0)) / total_height
        blocks.append(
            ConnectivityBlock(
                block_id=f"bed_cell_{i}",
                kind="bed_cell",
                axial_index=i,
                xi_reactor=float(xi),
                solver_coupling="explicit_state",
            )
        )
        if i > 0:
            edges.append(
                ConnectivityEdge(
                    f"bed_cell_{i - 1}",
                    f"bed_cell_{i}",
                    "axial_upflow",
                    "upward",
                    "explicit_state",
                )
            )
            edges.append(
                ConnectivityEdge(
                    f"bed_cell_{i - 1}",
                    f"bed_cell_{i}",
                    "wake_solid_upflow",
                    "upward",
                    "explicit_state",
                    phase_scope="suspension_phase_only",
                    mechanism="bubble_wake_carrying",
                )
            )
            edges.append(
                ConnectivityEdge(
                    f"bed_cell_{i}",
                    f"bed_cell_{i - 1}",
                    "internal_zirkulation",
                    "downward",
                    "explicit_state",
                    phase_scope="suspension_phase_only",
                    mechanism="mass_continuity_backmixing",
                )
            )

    top_source_block = f"bed_cell_{max(int(n_bed_cells) - 1, 0)}"
    if freeboard_active:
        interface_xi = float(H_bed / total_height)
        blocks.append(
            ConnectivityBlock(
                "bed_freeboard_interface",
                "interface_block",
                axial_index=int(max(n_bed_cells, 0)),
                xi_reactor=interface_xi,
                solver_coupling="explicit_handoff",
            )
        )
        edges.append(
            ConnectivityEdge(
                top_source_block,
                "bed_freeboard_interface",
                "bubble_rupture_handoff",
                "upward",
                "explicit_handoff",
            )
        )
        for i in range(int(max(n_freeboard_cells, 0))):
            xi = (float(H_bed) + (i + 0.5) * float(H_freeboard) / max(float(n_freeboard_cells), 1.0)) / total_height
            blocks.append(
                ConnectivityBlock(
                    block_id=f"freeboard_cell_{i}",
                    kind="freeboard_cell",
                    axial_index=i,
                    xi_reactor=float(xi),
                    solver_coupling="explicit_state",
                )
            )
            if i == 0:
                edges.append(
                    ConnectivityEdge(
                        "bed_freeboard_interface",
                        "freeboard_cell_0",
                        "freeboard_entry",
                        "upward",
                        "explicit_state",
                    )
                )
            else:
                edges.append(
                    ConnectivityEdge(
                        f"freeboard_cell_{i - 1}",
                        f"freeboard_cell_{i}",
                        "axial_upflow",
                        "upward",
                        "explicit_state",
                    )
                )
        top_source_block = f"freeboard_cell_{max(int(n_freeboard_cells) - 1, 0)}"

    if secondary_injection_applied and secondary_injection_segment is not None:
        target_block = (
            f"freeboard_cell_{int(secondary_injection_segment)}"
            if freeboard_active
            else top_source_block
        )
        blocks.append(
            ConnectivityBlock(
                "secondary_inlet_block",
                "side_feed",
                axial_index=int(secondary_injection_segment),
                solver_coupling="explicit_side_feed",
            )
        )
        edges.append(
            ConnectivityEdge(
                "secondary_inlet_block",
                target_block,
                "side_feed",
                "into_reactor",
                "explicit_side_feed",
            )
        )

    blocks.extend(
        [
            ConnectivityBlock(
                "cyclone_block",
                "cyclone",
                solver_coupling="explicit_state" if explicit_side_blocks else "post_freeboard_separation",
            ),
            ConnectivityBlock(
                "return_leg_block",
                "return_leg",
                solver_coupling="explicit_state" if explicit_side_blocks else "boundary_condition",
            ),
            ConnectivityBlock("system_exit", "system_exit", solver_coupling="postprocess"),
        ]
    )
    edges.extend(
        [
            ConnectivityEdge(
                top_source_block,
                "cyclone_block",
                "cyclone_inlet",
                "upward",
                "explicit_state" if explicit_side_blocks else "post_freeboard_separation",
            ),
            ConnectivityEdge("cyclone_block", "system_exit", "gas_exit", "out_of_system", "postprocess"),
        ]
    )

    recycle_active = bool(recycle_gas) or float(recirculation_frac) > 0.0 or float(cyclone_capture_char_frac) > 0.0 or float(cyclone_capture_ash_frac) > 0.0
    if recycle_active:
        edges.append(
            ConnectivityEdge(
                "cyclone_block",
                "return_leg_block",
                "solid_separation",
                "downward",
                "explicit_state" if explicit_side_blocks else "post_freeboard_separation",
            )
        )
        edges.append(
            ConnectivityEdge(
                "return_leg_block",
                "bed_cell_0",
                "external_zirkulation",
                "downward_to_bottom",
                "boundary_condition",
            )
        )

    return {
        "mode": "hamel_minimal_connectivity",
        "blocks": [asdict(block) for block in blocks],
        "edges": [asdict(edge) for edge in edges],
        "counts": {
            "n_blocks": len(blocks),
            "n_edges": len(edges),
            "n_bed_blocks": int(max(n_bed_cells, 0)),
            "n_freeboard_blocks": int(max(n_freeboard_cells, 0)) if freeboard_active else 0,
        },
        "freeboard_active": bool(freeboard_active),
        "recycle_active": bool(recycle_active),
        "secondary_injection_active": bool(secondary_injection_applied and secondary_injection_segment is not None),
    }
