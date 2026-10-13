"""Hamel-aligned analytical Jacobian partials for global NR.

Ref: Hamel (1999) §2.3 block-tridiagonal + Nebenelemente; Eq. (2-4)–(2-7).

Hybrid assembly:
- **Stencil** (zu/rez/in − inventory, no R/ex): ``project_stencil_balance_residual``
- **Analytic inventory diagonal**: ``inventory_self_derivative_unscaled``
- **Reaction / exchange / migration**: retained from full FD via
  ``J = J_analytic_diag + (J_fd − J_stencil_fd)`` on matching rows
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.core.cell import Cell, N_SOLID_COMP
from src.core.cell_balances import (
    S_ASH,
    S_CHAR,
    HOLDUP_K0_ANCHOR_NONE,
    calc_gas_balance_residual,
    calc_gas_enthalpy_flow,
    calc_solid_balance_residual,
    calc_solid_enthalpy_flow,
    holdup_transport_char_ash_mask,
)
from src.core.species import cp_ash, cp_char, cp_sand
from src.solvers.nr_indexing import (
    cell_kind,
    gas_mode_for_nr,
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_solid_var,
    n_var,
    solid_flat_indices_for_nr,
    temperature_offset,
)

__all__ = [
    "JACOBIAN_DERIVATIVE_SOURCES",
    "calc_bed_energy_outflow_cp_dot",
    "calc_side_temperature_closure_cp_dot",
    "side_temperature_reference_K",
    "analytic_packed_derivative",
    "apply_analytic_overlay_to_band_column",
    "apply_analytic_overlay_to_column",
    "inventory_inflow_coupling_derivative",
    "inventory_self_derivative_unscaled",
    "project_stencil_balance_residual",
]


@dataclass(frozen=True)
class JacobianDerivativeSource:
    residual_group: str
    hamel_ref: str
    state_variable: str
    analytic_status: str
    notes: str


JACOBIAN_DERIVATIVE_SOURCES: tuple[JacobianDerivativeSource, ...] = (
    JacobianDerivativeSource(
        "gas_dense",
        "Eq. (2-4)",
        "N_d,j",
        "analytic −1 inventory + FD(R, ex_bd)",
        "Hybrid: ∂(−N_d)/∂N_d analytic; reactions/exchange from full−stencil FD",
    ),
    JacobianDerivativeSource(
        "gas_bubble",
        "Eq. (2-5)",
        "N_b,j",
        "analytic −1 inventory + FD(R, ex_bd)",
        "Hybrid split on bubble inventory diagonal",
    ),
    JacobianDerivativeSource(
        "solid_holdup",
        "Eq. (2-6)",
        "m_jk",
        "analytic −(K_auf+K_ab) or anchor −1 + FD(R, mig)",
        "bed holdup with Jacobian K freeze; side-block solids stay FD",
    ),
    JacobianDerivativeSource(
        "energy_bed",
        "Eq. (2-7)",
        "T",
        "analytic −∂H_out/∂T + FD(hydro/ex coupling)",
        "Hybrid: outlet enthalpy cp_dot; H_in/q_loss fixed w.r.t. local T",
    ),
    JacobianDerivativeSource(
        "energy_side_closure",
        "Eq. (2-7) side Nebenelement",
        "T",
        "analytic",
        "cp_dot·(T_ref − T); ∂/∂T = −cp_dot",
    ),
    JacobianDerivativeSource(
        "gas_axial_inflow",
        "Eq. (2-4)/(2-5) N_d_in/N_b_in",
        "N upstream",
        "analytic +1 band + FD(R, ex)",
        "Immediate downstream row w.r.t. upstream packed gas DOF (bed/freeboard chain)",
    ),
    JacobianDerivativeSource(
        "solid_axial_inflow",
        "Eq. (2-6) m_in",
        "m upstream/above",
        "analytic lower/upper_frac (legacy) or K_auf/K_ab (holdup) band + FD",
        "legacy/holdup band overlay; holdup requires Jacobian K freeze during column FD",
    ),
    JacobianDerivativeSource(
        "side_recycle_coupling",
        "Nebenelemente rez",
        "N, m, T upstream",
        "FD on declared side blocks",
        "Graph-driven side_block_pairs; band LU + Woodbury solve",
    ),
)


def _project_balance_parts(
    cell: Cell,
    *,
    gas_parts: np.ndarray | None,
    solid_parts: np.ndarray | None,
    energy: float | None,
) -> np.ndarray:
    gas_mode = gas_mode_for_nr(cell)
    gas_idx = main_nr_gas_idx_for_cell(cell)
    parts: list[np.ndarray] = []
    if gas_parts is not None and gas_mode == "two_phase":
        half = gas_parts.size // 2
        parts.extend([gas_parts[:half], gas_parts[half:]])
    elif gas_parts is not None and gas_mode == "single":
        parts.append(gas_parts)
    if solid_parts is not None and solid_parts.size:
        parts.append(solid_parts)
    if energy is not None and has_temperature_var(cell):
        parts.append(np.array([float(energy)], dtype=np.float64))
    if not parts:
        return np.zeros(0, dtype=np.float64)
    out = np.concatenate(parts)
    assert out.size == n_var(cell), f"stencil projection size {out.size} != n_var {n_var(cell)}"
    return out


def project_stencil_balance_residual(cell: Cell) -> np.ndarray:
    """Hamel inventory stencil: zu/rez/in − holdup, without R or inter-phase exchange.

    Ref: Eq. (2-4)–(2-6) linear transport/inventory terms only.
    """
    kind = cell_kind(cell)
    if kind in {"cyclone", "return_leg"}:
        gas_parts = None
        solid_parts = None
        if n_solid_var(cell) > 0:
            solid_res = calc_solid_balance_residual(
                m_solid_zu=cell.m_solid_zu,
                m_solid_rez=cell.m_solid_rez,
                m_solid_in=cell.m_solid_in,
                m_solid_auf_in=cell.m_solid_auf_in,
                m_solid_ab_in=cell.m_solid_ab_in,
                R_solid=np.zeros_like(cell.m_solid),
                size_migration=np.zeros_like(cell.m_solid),
                m_solid=cell.m_solid,
                K_solid_auf=cell.K_solid_auf,
                K_solid_ab=cell.K_solid_ab,
                solid_state_model=str(cell.solid_state_model),
                k0_anchor=str(getattr(cell, "holdup_k0_anchor", "engineering")),
            )
            solid_parts = solid_res.reshape(-1)[solid_flat_indices_for_nr(cell)]
        energy = float(cell._calc_side_temperature_closure_residual())
        return _project_balance_parts(cell, gas_parts=gas_parts, solid_parts=solid_parts, energy=energy)

    if kind == "freeboard" and str(cell.solid_state_model) == "freeboard_closure":
        gas_res = calc_gas_balance_residual(
            N_zu_d=cell.N_zu_d,
            N_rez_d=cell.N_rez_d,
            N_d_in=cell.N_d_in,
            R_gas_d=np.zeros_like(cell.N_d),
            N_d=cell.N_d,
            N_ex=np.zeros_like(cell.N_d),
            N_zu_b=cell.N_zu_b,
            N_rez_b=cell.N_rez_b,
            N_b_in=cell.N_b_in,
            R_gas_b=np.zeros_like(cell.N_b),
            N_b=cell.N_b,
        )
        gas_idx = main_nr_gas_idx_for_cell(cell)
        gas_combined = gas_res[: gas_res.size // 2] + gas_res[gas_res.size // 2 :]
        gas_parts = gas_combined[gas_idx]
        return _project_balance_parts(cell, gas_parts=gas_parts, solid_parts=None, energy=cell.calc_energy_balance())

    gas_res = calc_gas_balance_residual(
        N_zu_d=cell.N_zu_d,
        N_rez_d=cell.N_rez_d,
        N_d_in=cell.N_d_in,
        R_gas_d=np.zeros_like(cell.N_d),
        N_d=cell.N_d,
        N_ex=np.zeros_like(cell.N_d),
        N_zu_b=cell.N_zu_b,
        N_rez_b=cell.N_rez_b,
        N_b_in=cell.N_b_in,
        R_gas_b=np.zeros_like(cell.N_b),
        N_b=cell.N_b,
    )
    gas_idx = main_nr_gas_idx_for_cell(cell)
    gas_parts = np.concatenate((gas_res[: gas_res.size // 2][gas_idx], gas_res[gas_res.size // 2 :][gas_idx]))

    solid_parts = None
    if n_solid_var(cell) > 0:
        solid_res = calc_solid_balance_residual(
            m_solid_zu=cell.m_solid_zu,
            m_solid_rez=cell.m_solid_rez,
            m_solid_in=cell.m_solid_in,
            m_solid_auf_in=cell.m_solid_auf_in,
            m_solid_ab_in=cell.m_solid_ab_in,
            R_solid=np.zeros_like(cell.m_solid),
            size_migration=np.zeros_like(cell.m_solid),
            m_solid=cell.m_solid,
            K_solid_auf=cell.K_solid_auf,
            K_solid_ab=cell.K_solid_ab,
            solid_state_model=str(cell.solid_state_model),
            k0_anchor=str(getattr(cell, "holdup_k0_anchor", "engineering")),
        )
        solid_parts = solid_res.reshape(-1)[solid_flat_indices_for_nr(cell)]

    energy = float(cell.calc_energy_balance())
    return _project_balance_parts(cell, gas_parts=gas_parts, solid_parts=solid_parts, energy=energy)


def _holdup_solid_self_derivative(cell: Cell, flat_idx: int) -> float | None:
    """Unscaled ∂(holdup stencil residual)/∂m for one NR solid DOF.

    Mirrors ``calc_solid_balance_residual`` anchor branches (R=0, no migration)
    for char/ash transport buckets only.
    """
    flat_indices = solid_flat_indices_for_nr(cell)
    pos = int(np.where(flat_indices == flat_idx)[0][0])
    nk, n_comp = cell.m_solid.shape
    nk_idx = pos // int(n_comp)
    comp_idx = pos % int(n_comp)
    transport_mask = holdup_transport_char_ash_mask((nk, n_comp))
    if float(transport_mask[nk_idx, comp_idx]) <= 0.0:
        return None
    m_solid_auf = cell.m_solid_auf_in * transport_mask
    m_solid_ab = cell.m_solid_ab_in * transport_mask
    k_auf = cell.K_solid_auf * transport_mask
    k_ab = cell.K_solid_ab * transport_mask
    k_sum = k_auf + k_ab
    net_source = (
        cell.m_solid_zu
        + cell.m_solid_rez
        + cell.m_solid_in
        + m_solid_auf
        + m_solid_ab
        + np.zeros_like(cell.m_solid)
    )
    k_val = float(k_sum[nk_idx, comp_idx])
    net_val = float(net_source[nk_idx, comp_idx])
    abs_net = abs(net_val)
    if str(getattr(cell, "holdup_k0_anchor", "engineering")) == HOLDUP_K0_ANCHOR_NONE:
        return -k_val
    if abs(k_val) <= 1e-14 and abs_net <= 1e-14:
        return -1.0
    if abs(k_val) <= 1e-14 and abs_net > 1e-14:
        return -1.0
    if k_val > 1e-14:
        return -k_val
    return None


def inventory_self_derivative_unscaled(
    cell: Cell,
    *,
    row_local: int,
    col_local: int,
) -> float | None:
    """Unscaled ∂F_row/∂x_col for Hamel inventory / side-closure diagonal."""
    if row_local != col_local:
        return None
    kind = cell_kind(cell)
    ng = n_gas_var(cell)
    gas_mode = gas_mode_for_nr(cell)

    if kind in {"cyclone", "return_leg"}:
        if has_temperature_var(cell):
            try:
                t_off = temperature_offset(cell)
            except ValueError:
                return None
            if row_local == t_off:
                return -calc_side_temperature_closure_cp_dot(cell)
        return None

    if kind == "bed" and has_temperature_var(cell):
        try:
            t_off = temperature_offset(cell)
        except ValueError:
            return None
        if row_local == col_local == t_off:
            return -calc_bed_energy_outflow_cp_dot(cell)

    if gas_mode == "two_phase" and kind == "bed":
        gas_half = ng // 2
        if col_local < gas_half and row_local == col_local:
            return -1.0
        if gas_half <= col_local < ng and row_local == col_local:
            return -1.0
        sol_off = ng
        nv_sol = n_solid_var(cell)
        if col_local >= sol_off and row_local == col_local and nv_sol > 0:
            local = col_local - sol_off
            if local >= nv_sol:
                return None
            flat = solid_flat_indices_for_nr(cell)
            if str(cell.solid_state_model) == "legacy_stream":
                return -1.0
            if str(cell.solid_state_model) == "holdup_transport":
                return _holdup_solid_self_derivative(cell, int(flat[local]))
        return None

    if gas_mode == "single" and col_local < ng and row_local == col_local:
        return -1.0

    nv_sol = n_solid_var(cell)
    if nv_sol > 0:
        sol_off = ng
        if col_local >= sol_off and row_local == col_local:
            local = col_local - sol_off
            if local >= nv_sol:
                return None
            flat = solid_flat_indices_for_nr(cell)
            if str(cell.solid_state_model) == "legacy_stream":
                return -1.0
            if str(cell.solid_state_model) == "holdup_transport" and kind == "bed":
                return _holdup_solid_self_derivative(cell, int(flat[local]))
    return None


def side_temperature_reference_K(cell: Cell) -> float:
    if cell_kind(cell) == "cyclone" and float(np.sum(np.maximum(cell.N_b_in + cell.N_d_in, 0.0))) > 1e-12:
        return float(cell.T_in_gas)
    solid_ab = float(np.sum(np.maximum(cell.m_solid_ab_in, 0.0)))
    solid_auf = float(np.sum(np.maximum(cell.m_solid_auf_in, 0.0)))
    if solid_ab > 1e-12:
        return float(cell.T_solid_ab_in)
    if solid_auf > 1e-12:
        return float(cell.T_solid_auf_in)
    return float(cell.T_in_solid)


def calc_bed_energy_outflow_cp_dot(cell: Cell) -> float:
    """Effective Cp·ṁ [W/K] for Eq. (2-7) outlet enthalpy H_aus(T).

    Uses current ``N_b+N_d`` and ``_solid_outflow_rates()`` with frozen flows w.r.t. T.
    Ref: ``Cell.calc_energy_balance`` → ``calc_energy_balance_residual`` H_out term.
    """
    if cell_kind(cell) != "bed":
        return 0.0

    def _h_out(T: float) -> float:
        gas_out = np.maximum(cell.N_b + cell.N_d, 0.0)
        solid_out = np.maximum(cell._solid_outflow_rates(), 0.0)
        return float(cell._calc_gas_enthalpy_flow(gas_out, T) + cell._calc_solid_enthalpy_flow(solid_out, T))

    T0 = float(cell.T)
    dT = 1.0e-5
    return max((_h_out(T0 + dT) - _h_out(T0)) / dT, 1.0e-12)


def calc_side_temperature_closure_cp_dot(cell: Cell) -> float:
    """Effective Cp·ṁ [W/K] for side-element temperature closure."""
    T_ref = side_temperature_reference_K(cell)
    gas_in = np.maximum(cell.N_b_in + cell.N_d_in, 0.0)
    solid_in = np.maximum(cell.m_solid_in + cell.m_solid_auf_in + cell.m_solid_ab_in, 0.0)
    dT = 1.0
    cp_dot = 0.0
    if float(np.sum(gas_in)) > 1e-12:
        cp_dot += abs(
            calc_gas_enthalpy_flow(gas_in, T_ref + dT, cell._h_cache)
            - calc_gas_enthalpy_flow(gas_in, T_ref, cell._h_cache)
        ) / dT
    if float(np.sum(solid_in)) > 1e-12:
        cp_dot += abs(
            calc_solid_enthalpy_flow(
                solid_in,
                T_ref + dT,
                ash_dry_wt=cell.solid.ash_dry_wt,
                VM_daf=cell.solid.VM_daf,
                h_f_dry=cell.solid.h_f_dry,
                cp_char_fn=cp_char,
                cp_ash_fn=cp_ash,
                cp_sand_fn=cp_sand,
                solid_enthalpy_mode=cell.solid_enthalpy_mode,
            )
            - calc_solid_enthalpy_flow(
                solid_in,
                T_ref,
                ash_dry_wt=cell.solid.ash_dry_wt,
                VM_daf=cell.solid.VM_daf,
                h_f_dry=cell.solid.h_f_dry,
                cp_char_fn=cp_char,
                cp_ash_fn=cp_ash,
                cp_sand_fn=cp_sand,
                solid_enthalpy_mode=cell.solid_enthalpy_mode,
            )
        ) / dT
    return max(float(cp_dot), 1.0e3)


def analytic_packed_derivative(
    cell: Cell,
    *,
    row_local: int,
    col_local: int,
) -> float | None:
    return inventory_self_derivative_unscaled(cell, row_local=row_local, col_local=col_local)


def _holdup_solid_band_coupling_derivative(
    row_cell: Cell,
    col_cell: Cell,
    *,
    row_local: int,
    col_local: int,
    block_delta: int,
) -> float | None:
    """Unscaled ∂F_row/∂m_col for holdup_transport axial K·m inflow stencil.

    Ref: ``_set_bed_solid_auf_in_from_source`` / ``_set_bed_solid_ab_in_from_bed_source``;
    Hamel (1999) Eq. (2-6) with frozen ``K_auf``/``K_ab`` (not packed NR DOFs).
    """
    if str(row_cell.solid_state_model) != "holdup_transport":
        return None
    if str(col_cell.solid_state_model) != "holdup_transport":
        return None
    row_ng = n_gas_var(row_cell)
    col_ng = n_gas_var(col_cell)
    if row_local < row_ng or col_local < col_ng or row_local != col_local:
        return None
    row_nv = n_solid_var(row_cell)
    col_nv = n_solid_var(col_cell)
    row_sol = row_local - row_ng
    col_sol = col_local - col_ng
    if row_nv <= 0 or col_nv <= 0 or row_sol != col_sol:
        return None
    if row_sol < 0 or row_sol >= row_nv or col_sol < 0 or col_sol >= col_nv:
        return None
    flat_idx = int(solid_flat_indices_for_nr(row_cell)[row_sol])
    nk_idx = flat_idx // int(N_SOLID_COMP)
    comp_idx = flat_idx % int(N_SOLID_COMP)
    if comp_idx not in (S_CHAR, S_ASH):
        return None
    m_col = float(np.maximum(col_cell.m_solid[nk_idx, comp_idx], 0.0))
    if m_col <= 0.0:
        return 0.0
    if block_delta == 1:
        return float(col_cell.K_solid_auf[nk_idx, comp_idx])
    if block_delta == -1:
        return float(col_cell.K_solid_ab[nk_idx, comp_idx])
    return None


def inventory_inflow_coupling_derivative(
    row_cell: Cell,
    col_cell: Cell,
    *,
    row_local: int,
    col_local: int,
    row_block: int,
    col_block: int,
) -> float | None:
    """Unscaled band coupling for axial gas (+1) and solid inlet stencil.

    Ref: ``routing._apply_bed_upstream_inlet`` — ``N_d_in = prev.N_d``;
    legacy ``m_solid_in = upper_frac * above + lower_frac * prev``;
    holdup ``m_solid_auf_in = K_auf_col * m_col``, ``m_solid_ab_in = K_ab_col * m_col``.
    """
    block_delta = int(row_block) - int(col_block)
    if block_delta not in (-1, 1):
        return None
    row_kind = cell_kind(row_cell)
    col_kind = cell_kind(col_cell)
    if row_kind not in {"bed", "freeboard"} or col_kind not in {"bed", "freeboard"}:
        return None

    row_ng = n_gas_var(row_cell)
    col_ng = n_gas_var(col_cell)
    if (
        block_delta == 1
        and row_local == col_local
        and row_local < row_ng
        and col_local < col_ng
    ):
        row_mode = gas_mode_for_nr(row_cell)
        col_mode = gas_mode_for_nr(col_cell)
        if row_mode == col_mode:
            if row_mode == "two_phase" and row_kind == "bed":
                return 1.0
            if row_mode == "single":
                return 1.0

    if (
        row_kind == "bed"
        and col_kind == "bed"
        and str(row_cell.solid_state_model) == "legacy_stream"
        and row_local == col_local
        and row_local >= row_ng
        and col_local >= col_ng
    ):
        row_nv = n_solid_var(row_cell)
        col_nv = n_solid_var(col_cell)
        row_sol = row_local - row_ng
        col_sol = col_local - col_ng
        if row_nv <= 0 or col_nv <= 0 or row_sol != col_sol:
            return None
        if row_sol < 0 or row_sol >= row_nv or col_sol < 0 or col_sol >= col_nv:
            return None
        lower = float(getattr(row_cell, "_nr_solid_lower_inlet_frac", 0.0))
        lower = float(np.clip(lower, 0.0, 1.0))
        if block_delta == 1:
            return lower
        return 1.0 - lower

    holdup = _holdup_solid_band_coupling_derivative(
        row_cell,
        col_cell,
        row_local=row_local,
        col_local=col_local,
        block_delta=block_delta,
    )
    if holdup is not None:
        return holdup
    return None


def _apply_single_analytic_overlay(
    out: np.ndarray,
    *,
    row_local: int,
    raw: float | None,
    dF_local: np.ndarray,
    eq_scale_local: np.ndarray,
    dF_stencil_local: np.ndarray | None,
) -> None:
    if raw is None:
        return
    scale = float(eq_scale_local[row_local])
    if scale <= 0.0:
        return
    if dF_stencil_local is not None:
        stencil_val = float(dF_stencil_local[row_local])
        analytic_scaled = float(raw) / scale
        ref = max(abs(stencil_val), abs(analytic_scaled), 1.0e-12)
        if abs(analytic_scaled - stencil_val) > 0.05 * ref:
            return
        out[row_local] = analytic_scaled + (float(dF_local[row_local]) - stencil_val)
    else:
        out[row_local] = float(raw) / scale


def apply_analytic_overlay_to_band_column(
    dF_local: np.ndarray,
    row_cell: Cell,
    col_cell: Cell,
    *,
    col_local: int,
    eq_scale_local: np.ndarray,
    dF_stencil_local: np.ndarray | None,
    row_block: int,
    col_block: int,
) -> np.ndarray:
    """Hybrid overlay for immediate downstream band blocks (+1 gas inflow stencil)."""
    out = np.asarray(dF_local, dtype=np.float64).copy()
    for row_local in range(n_var(row_cell)):
        raw = inventory_inflow_coupling_derivative(
            row_cell,
            col_cell,
            row_local=row_local,
            col_local=col_local,
            row_block=row_block,
            col_block=col_block,
        )
        _apply_single_analytic_overlay(
            out,
            row_local=row_local,
            raw=raw,
            dF_local=dF_local,
            eq_scale_local=eq_scale_local,
            dF_stencil_local=dF_stencil_local,
        )
    return out


def apply_analytic_overlay_to_column(
    dF_local: np.ndarray,
    cell: Cell,
    *,
    col_local: int,
    eq_scale_local: np.ndarray,
    dF_stencil_local: np.ndarray | None = None,
) -> np.ndarray:
    """Hybrid overlay: analytic inventory / side-T diagonal + (full FD − stencil FD)."""
    out = np.asarray(dF_local, dtype=np.float64).copy()
    nv = n_var(cell)
    if col_local < 0 or col_local >= nv:
        return out
    raw = inventory_self_derivative_unscaled(cell, row_local=col_local, col_local=col_local)
    _apply_single_analytic_overlay(
        out,
        row_local=col_local,
        raw=raw,
        dF_local=dF_local,
        eq_scale_local=eq_scale_local,
        dF_stencil_local=dF_stencil_local,
    )
    return out
