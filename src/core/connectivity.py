"""Thesis-mode connectivity and routing for Hamel-style reactor blocks."""

from __future__ import annotations

from typing import Any
from dataclasses import asdict, dataclass

import numpy as np

from src.core.cell import Cell, S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, cp_ash, cp_char, cp_sand

_F_W = 0.25
_ZETA_W = 0.40


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

    This preserves each primary species' total molar flow and only changes the
    dense/bubble partition.  It prevents outer Vorabrechnung refresh from moving
    the bottom inlet split while leaving the accepted state on an obsolete split.
    """
    if not cells:
        return
    cell = cells[0]
    dense_frac = resolve_bottom_gas_inlet_dense_fraction(cell, cfg)
    idx = GAS_SPECIES_INDEX
    for sp in ("O2", "H2O", "N2"):
        j = idx[sp]
        total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
        cell.N_d[j] = max(total * dense_frac, 1e-12)
        cell.N_b[j] = max(total * (1.0 - dense_frac), 1e-12)


def preproject_bottom_primary_gas_state_to_exchange_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    improvement_factor: float = 0.75,
) -> bool:
    """Preproject bed0 primary gas split against the local K_bd exchange closure.

    Hamel gives the Vorabrechnung module responsibility for hydrodynamics and
    the outer Abgleich with the cell model, but does not prescribe an explicit
    start-value projection.  This bridge is therefore deliberately narrow:
    only O2/H2O/N2 in the bottom bed cell are repartitioned, each species'
    total molar flow is conserved, and reaction products are left untouched.
    """
    if not cells:
        return False
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    idxs = [GAS_SPECIES_INDEX[sp] for sp in ("O2", "H2O", "N2")]
    totals = np.array([float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in idxs], dtype=np.float64)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 1e-10):
        return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()

    def _invalidate_gas_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _set_dense_values(values: np.ndarray) -> None:
        clipped = np.clip(np.asarray(values, dtype=np.float64), 1e-12, totals - 1e-12)
        for value, total, j in zip(clipped, totals, idxs):
            cell.N_d[j] = float(value)
            cell.N_b[j] = float(total - value)
        _invalidate_gas_state_cache()

    def _primary_dense_residual(values: np.ndarray) -> np.ndarray:
        _set_dense_values(values)
        res = np.asarray(cell.residuals(), dtype=np.float64)
        return res[idxs].copy()

    x0 = np.array([float(cell.N_d[j]) for j in idxs], dtype=np.float64)
    initial_res = _primary_dense_residual(x0)
    initial_norm = float(np.linalg.norm(initial_res, ord=2))
    if not np.isfinite(initial_norm) or initial_norm <= 1e-12:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        return False

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            _primary_dense_residual,
            x0=np.clip(x0, 1e-12, totals - 1e-12),
            bounds=(np.full_like(totals, 1e-12), totals - 1e-12),
            max_nfev=12,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        candidate = np.asarray(sol.x, dtype=np.float64)
        final_res = _primary_dense_residual(candidate)
    except Exception:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        cell.residuals()
        return False

    final_norm = float(np.linalg.norm(final_res, ord=2))
    if np.isfinite(final_norm) and final_norm <= improvement_factor * initial_norm:
        _set_dense_values(candidate)
        cell.residuals()
        return True

    cell.N_d[:] = original_d
    cell.N_b[:] = original_b
    _invalidate_gas_state_cache()
    cell.residuals()
    return False


def preproject_bottom_major_gas_state_to_local_balance(
    cells: list[Cell],
    cfg: Any,
    *,
    improvement_factor: float = 0.75,
) -> bool:
    """Preproject bed0 major gas x0 against the local two-phase gas balances.

    This is a start-value bridge for the bottom oxidation/devolatilization zone:
    primary species totals (O2/H2O/N2) remain fixed, while major kinetic products
    (CO/H2/CH4/CO2) may receive dense and bubble outlet support.  Tar species are
    intentionally excluded because the current tar source is a separate lumped
    pyrolysis carrier and including it worsened the O2/H2O split in static probes.
    """
    if not cells:
        return False
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    primary_names = ("O2", "H2O", "N2")
    product_names = ("CO", "H2", "CH4", "CO2")
    primary_idxs = [GAS_SPECIES_INDEX[sp] for sp in primary_names]
    product_idxs = [GAS_SPECIES_INDEX[sp] for sp in product_names]
    objective_idxs = primary_idxs + product_idxs
    objective_rows = objective_idxs + [j + cell.N_d.shape[0] for j in objective_idxs]
    primary_totals = np.array(
        [float(max(cell.N_d[j] + cell.N_b[j], 0.0)) for j in primary_idxs],
        dtype=np.float64,
    )
    if np.any(~np.isfinite(primary_totals)) or np.any(primary_totals <= 1e-10):
        return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()

    def _invalidate_gas_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _set_values(values: np.ndarray) -> None:
        x = np.asarray(values, dtype=np.float64)
        k = 0
        for value, total, j in zip(x[k : k + len(primary_idxs)], primary_totals, primary_idxs):
            dense = float(np.clip(value, 1e-12, total - 1e-12))
            cell.N_d[j] = dense
            cell.N_b[j] = float(total - dense)
        k += len(primary_idxs)
        for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
            cell.N_d[j] = max(float(value), 1e-12)
        k += len(product_idxs)
        for value, j in zip(x[k : k + len(product_idxs)], product_idxs):
            cell.N_b[j] = max(float(value), 1e-12)
        _invalidate_gas_state_cache()

    def _residual_objective(values: np.ndarray) -> np.ndarray:
        _set_values(values)
        res = np.asarray(cell.residuals(), dtype=np.float64)
        return res[objective_rows].copy()

    initial_res = np.asarray(cell.residuals(), dtype=np.float64)
    initial_gas_max = float(
        max(
            np.max(np.abs(initial_res[: cell.N_d.shape[0]])),
            np.max(np.abs(initial_res[cell.N_d.shape[0] : 2 * cell.N_d.shape[0]])),
        )
    )
    if not np.isfinite(initial_gas_max) or initial_gas_max <= 1e-12:
        return False

    x0 = np.concatenate(
        (
            np.array([float(cell.N_d[j]) for j in primary_idxs], dtype=np.float64),
            np.array([max(float(cell.N_d[j]), 1e-9) for j in product_idxs], dtype=np.float64),
            np.array([max(float(cell.N_b[j]), 1e-9) for j in product_idxs], dtype=np.float64),
        )
    )
    lower = np.full_like(x0, 1e-12)
    upper = np.concatenate(
        (
            primary_totals - 1e-12,
            np.full(len(product_idxs), 80.0, dtype=np.float64),
            np.full(len(product_idxs), 80.0, dtype=np.float64),
        )
    )

    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            _residual_objective,
            x0=np.clip(x0, lower, upper),
            bounds=(lower, upper),
            max_nfev=20,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        candidate = np.asarray(sol.x, dtype=np.float64)
        _set_values(candidate)
        final_res = np.asarray(cell.residuals(), dtype=np.float64)
    except Exception:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        _invalidate_gas_state_cache()
        cell.residuals()
        return False

    final_gas_max = float(
        max(
            np.max(np.abs(final_res[: cell.N_d.shape[0]])),
            np.max(np.abs(final_res[cell.N_d.shape[0] : 2 * cell.N_d.shape[0]])),
        )
    )
    if np.isfinite(final_gas_max) and final_gas_max <= improvement_factor * initial_gas_max:
        return True

    cell.N_d[:] = original_d
    cell.N_b[:] = original_b
    _invalidate_gas_state_cache()
    cell.residuals()
    return False


def preproject_bottom_total_gas_and_temperature_to_energy_closure(
    cells: list[Cell],
    cfg: Any,
    *,
    gas_relax: float = 1.0,
    temperature_bounds_K: tuple[float, float] = (650.0, 1600.0),
    improvement_factor: float = 0.75,
) -> bool:
    """Joint bed0 x0 projection for total gas closure and energy accessibility.

    This is an initialization-only bridge for the bottom devolatilization /
    oxidation stiffness.  It does not modify reaction sources and does not put
    drying/pyrolysis gas into the bubble source.  Instead:

    - major gas dense/bubble outlet values are solved with bounded least squares
      against total gas closure, phase residuals, and energy closure;
    - TAR/NH3 source support is added only to dense outlet variables because the
      direct pyrolysis source remains a suspension-phase source;
    - temperature is bounded to keep the projection in a plausible x0 range.
    """
    if not cells:
        return False
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return False

    original_d = cell.N_d.copy()
    original_b = cell.N_b.copy()
    original_T = float(cell.T)

    def _invalidate_state_cache() -> None:
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False

    def _gas_max(res: np.ndarray) -> float:
        n = cell.N_d.shape[0]
        return float(max(np.max(np.abs(res[:n])), np.max(np.abs(res[n : 2 * n]))))

    def _combined_gas_max(res: np.ndarray) -> float:
        n = cell.N_d.shape[0]
        return float(np.max(np.abs(res[:n] + res[n : 2 * n])))

    initial_res = np.asarray(cell.residuals(), dtype=np.float64)
    initial_energy_abs = abs(float(initial_res[-1]))
    initial_gas_max = _gas_max(initial_res)
    initial_combined_gas_max = _combined_gas_max(initial_res)
    if not np.isfinite(initial_energy_abs) or initial_energy_abs <= 1e-9:
        return False

    n_gas = cell.N_d.shape[0]
    projection_omega = float(np.clip(gas_relax, 0.0, 1.0))
    major_names = ("CO", "CO2", "H2", "H2O", "CH4", "O2", "N2")
    major_idxs = [GAS_SPECIES_INDEX[sp] for sp in major_names]
    dense_only_source_names = ("NH3", "TAR1", "TAR2")

    fixed_dense_support = original_d.copy()
    fixed_bubble_support = original_b.copy()
    for sp in dense_only_source_names:
        j = GAS_SPECIES_INDEX[sp]
        dense_support = max(float(initial_res[j]), 0.0)
        fixed_dense_support[j] = max(float(fixed_dense_support[j] + projection_omega * dense_support), 0.0)
        # Keep bubble untouched: direct drying/pyrolysis source remains dense-only.
        fixed_bubble_support[j] = max(float(fixed_bubble_support[j]), 0.0)

    def _h_gas(flow: np.ndarray, temperature: float) -> float:
        return calc_gas_enthalpy_flow(np.asarray(flow, dtype=np.float64), float(temperature), cell._h_cache)

    def _h_solid(flow: np.ndarray, temperature: float) -> float:
        return calc_solid_enthalpy_flow(
            np.asarray(flow, dtype=np.float64),
            float(temperature),
            ash_dry_wt=cell.solid.ash_dry_wt,
            VM_daf=cell.solid.VM_daf,
            h_f_dry=cell.solid.h_f_dry,
            cp_char_fn=cp_char,
            cp_ash_fn=cp_ash,
            cp_sand_fn=cp_sand,
        )

    attr_lo = float(getattr(cell, "_nr_temperature_min_K", temperature_bounds_K[0]))
    attr_hi = float(getattr(cell, "_nr_temperature_max_K", temperature_bounds_K[1]))
    if not np.isfinite(attr_lo):
        attr_lo = float(temperature_bounds_K[0])
    if not np.isfinite(attr_hi):
        attr_hi = float(temperature_bounds_K[1])
    t_lo = max(float(temperature_bounds_K[0]), attr_lo, 300.0)
    t_hi = min(float(temperature_bounds_K[1]), attr_hi, 2500.0)
    t_hi = max(t_hi, t_lo + 1.0)

    def _set_candidate(x: np.ndarray) -> None:
        arr = np.asarray(x, dtype=np.float64)
        cell.N_d[:] = fixed_dense_support
        cell.N_b[:] = fixed_bubble_support
        for value, j in zip(arr[: len(major_idxs)], major_idxs):
            cell.N_d[j] = max(float(value), 1e-12)
        for value, j in zip(arr[len(major_idxs) : 2 * len(major_idxs)], major_idxs):
            cell.N_b[j] = max(float(value), 1e-12)
        cell.T = float(np.clip(arr[-1], t_lo, t_hi))
        _invalidate_state_cache()

    ref_gas = max(
        float(
            np.sum(
                np.maximum(
                    cell.N_zu_d
                    + cell.N_zu_b
                    + cell.N_d_in
                    + cell.N_b_in
                    + cell.N_rez_d
                    + cell.N_rez_b,
                    0.0,
                )
            )
        ),
        1.0,
    )
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)

    def _objective(x: np.ndarray) -> np.ndarray:
        _set_candidate(x)
        res = np.asarray(cell.residuals(), dtype=np.float64)
        combined = (res[major_idxs] + res[[n_gas + j for j in major_idxs]]) / ref_gas
        phase = np.concatenate((res[major_idxs], res[[n_gas + j for j in major_idxs]])) / ref_gas
        energy = np.array([float(res[-1]) / ref_energy], dtype=np.float64)
        # Total closure is the priority; phase rows remain in the objective so
        # the projection cannot "solve" total gas by destroying two-phase balance.
        return np.concatenate((combined, 0.4 * phase, 3.0 * energy))

    x0 = np.concatenate(
        (
            np.array([max(float(original_d[j]), 1e-12) for j in major_idxs], dtype=np.float64),
            np.array([max(float(original_b[j]), 1e-12) for j in major_idxs], dtype=np.float64),
            np.array([float(np.clip(original_T, t_lo, t_hi))], dtype=np.float64),
        )
    )
    lower = np.full_like(x0, 1e-12)
    upper = np.full_like(x0, 120.0)
    lower[-1] = t_lo
    upper[-1] = t_hi
    try:
        from scipy.optimize import least_squares

        sol = least_squares(
            _objective,
            x0=np.clip(x0, lower, upper),
            bounds=(lower, upper),
            max_nfev=40,
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
        )
        raw_candidate = np.asarray(sol.x, dtype=np.float64)
        # Use the bounded LSQ result as a direction, not as a hidden local solve.
        # A full projection closes bed0 aggressively but can simply push an
        # enthalpy shock into bed1/bed2.  Relaxing the state move keeps this as a
        # start-value preconditioner.
        relaxed_candidate = np.asarray(x0 + projection_omega * (raw_candidate - x0), dtype=np.float64)
        _set_candidate(relaxed_candidate)
        final_res = np.asarray(cell.residuals(), dtype=np.float64)
    except Exception:
        cell.N_d[:] = original_d
        cell.N_b[:] = original_b
        cell.T = original_T
        _invalidate_state_cache()
        cell.residuals()
        return False

    final_energy_abs = abs(float(final_res[-1]))
    final_combined_gas_max = _combined_gas_max(final_res)
    final_gas_max = _gas_max(final_res)

    energy_improved = np.isfinite(final_energy_abs) and final_energy_abs <= improvement_factor * initial_energy_abs
    gas_improved = (
        np.isfinite(final_combined_gas_max)
        and final_combined_gas_max <= improvement_factor * max(initial_combined_gas_max, 1e-9)
        and np.isfinite(final_gas_max)
        and final_gas_max <= max(initial_gas_max, 1e-9)
    )
    if energy_improved and gas_improved:
        cell._bottom_joint_preprojection_diag = {
            "initial_energy_abs_W": float(initial_energy_abs),
            "final_energy_abs_W": float(final_energy_abs),
            "initial_gas_max_mol_s": float(initial_gas_max),
            "final_gas_max_mol_s": float(final_gas_max),
            "initial_combined_gas_max_mol_s": float(initial_combined_gas_max),
            "final_combined_gas_max_mol_s": float(final_combined_gas_max),
            "candidate_temperature_K": float(cell.T),
        }
        return True

    cell.N_d[:] = original_d
    cell.N_b[:] = original_b
    cell.T = original_T
    _invalidate_state_cache()
    cell.residuals()
    return False


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
            # Propagate the remaining resident reactive stock, not the original
            # fresh/source support. Reusing m_solid_zu + m_solid_in would count
            # the same VM/moisture budget again in every lower-bed cell.
            reactive_source = np.maximum(cell.m_solid, 0.0)
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
    top_downflow_components: np.ndarray | None = None,
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

    m_ab_from_above_total = float(max(top_downflow_total, 0.0))
    m_ab_from_above_components = np.zeros(2, dtype=np.float64)
    if isinstance(top_downflow_components, np.ndarray) and top_downflow_components.shape[0] >= 2:
        m_ab_from_above_components[:] = np.maximum(top_downflow_components[:2], 0.0)
    elif m_ab_from_above_total > 0.0:
        m_ab_from_above_components[:] = 0.5 * m_ab_from_above_total

    for i in range(len(bed_cells) - 1, -1, -1):
        cell = bed_cells[i]
        up_in = float(m_auf_totals[i - 1]) if i > 0 else 0.0
        ext_in_char = float(np.sum(np.maximum(cell.m_solid_zu[:, S_CHAR] + cell.m_solid_rez[:, S_CHAR] + cell.m_solid_in[:, S_CHAR], 0.0)))
        ext_in_ash = float(np.sum(np.maximum(cell.m_solid_zu[:, S_ASH] + cell.m_solid_rez[:, S_ASH] + cell.m_solid_in[:, S_ASH], 0.0)))
        reaction_char = float(np.sum(np.asarray(cell.R_solid[:, S_CHAR], dtype=np.float64)))
        reaction_ash = float(np.sum(np.asarray(cell.R_solid[:, S_ASH], dtype=np.float64)))
        # Bed axial upflow remains hydrodynamics-driven (shared ``K_auf``), while the
        # top-down closure for ``K_ab`` keeps char/ash source terms separated.
        m_ab_total = max(
            up_in - float(m_auf_totals[i]) + m_ab_from_above_total + ext_in_char + ext_in_ash + reaction_char + reaction_ash,
            0.0,
        )
        m_ab_char = max(up_in - float(m_auf_totals[i]) + m_ab_from_above_components[0] + ext_in_char + reaction_char, 0.0)
        m_ab_ash = max(up_in - float(m_auf_totals[i]) + m_ab_from_above_components[1] + ext_in_ash + reaction_ash, 0.0)
        k_ab_char = m_ab_char / max(float(holdups[i]), 1e-12) if holdups[i] > 1e-12 else 0.0
        k_ab_ash = m_ab_ash / max(float(holdups[i]), 1e-12) if holdups[i] > 1e-12 else 0.0

        cell.K_solid_auf.fill(0.0)
        cell.K_solid_ab.fill(0.0)
        # Thesis holdup transport: only char/ash are transported across cells.
        cell.K_solid_auf[:, S_CHAR] = float(k_auf_vals[i])
        cell.K_solid_auf[:, S_ASH] = float(k_auf_vals[i])
        cell.K_solid_ab[:, S_CHAR] = float(k_ab_char)
        cell.K_solid_ab[:, S_ASH] = float(k_ab_ash)
        m_ab_from_above_total = m_ab_total
        m_ab_from_above_components[0] = m_ab_char
        m_ab_from_above_components[1] = m_ab_ash


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
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)

    for i, cell in enumerate(bed_cells):
        if i > 0:
            below = bed_cells[i - 1]
            cell.m_solid_auf_in[:, :] = below.K_solid_auf * np.maximum(below.m_solid, 0.0)
            cell.T_solid_auf_in = float(below.T)
        if i + 1 < len(bed_cells):
            above = bed_cells[i + 1]
            cell.m_solid_ab_in[:, :] = above.K_solid_ab * np.maximum(above.m_solid, 0.0)
            cell.T_solid_ab_in = float(above.T)
        elif top_above_cell is not None:
            cell.m_solid_ab_in[:, :] = np.maximum(top_above_cell._solid_downflow_rates(), 0.0)
            cell.T_solid_ab_in = float(top_above_cell.T)
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


def _set_degenerate_single_phase_gas_inlet(prev: Cell, cell: Cell) -> None:
    """Route upstream total gas into the suspension channel of degenerate segments."""
    prev_kind = str(getattr(prev, "cell_type", "bed")).strip().lower()
    prev_total = prev.N_d if prev_kind in {"freeboard", "cyclone", "return_leg"} else prev.N_d + prev.N_b
    cell.N_d_in[:] = np.maximum(prev_total, 0.0)
    cell.N_b_in.fill(0.0)


def _set_explicit_freeboard_inlet_from_prev(prev: Cell, cell: Cell) -> None:
    """Refresh one explicit freeboard cell inlet from its upstream source cell."""
    _set_degenerate_single_phase_gas_inlet(prev, cell)
    cell.T_in_gas = float(prev.T)
    cell.T_in_solid = float(prev.T)
    cell.m_solid_in.fill(0.0)
    # Freeboard-closure solids remain outside the NR solid residual, but their
    # entrained char/ash stream still carries sensible/formation enthalpy into
    # the freeboard energy residual.
    cell.m_solid_in[:] = propagated_solid_stream_from_cell(prev, include_reactive=False)
    active_mask = getattr(cell, "_freeboard_active_char_ash_mask", None)
    if (
        str(getattr(cell, "solid_state_model", "legacy_stream")) != "freeboard_closure"
        and isinstance(active_mask, np.ndarray)
        and active_mask.shape[0] == cell.solid.n_size_classes
    ):
        mask = np.maximum(active_mask, 0.0).reshape(-1)
        cell.m_solid_in[:, S_CHAR] *= mask
        cell.m_solid_in[:, S_ASH] *= mask
    if str(getattr(cell, "solid_state_model", "legacy_stream")) == "holdup_transport":
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
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(bottom_below_cell.T)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(below.T)
        if i + 1 < len(freeboard_cells):
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
            cell.T_solid_ab_in = float(above.T)
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
        cell.T_solid_auf_in = float(cell.T)
        cell.T_solid_ab_in = float(cell.T)
        if str(getattr(cell, "solid_state_model", "legacy_stream")) == "freeboard_closure":
            continue
        if i == 0:
            cell.m_solid_auf_in[:, :] = np.maximum(bottom_below_cell._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(bottom_below_cell.T)
        else:
            below = freeboard_cells[i - 1]
            cell.m_solid_auf_in[:, :] = np.maximum(below._solid_upflow_rates(), 0.0)
            cell.T_solid_auf_in = float(below.T)
        if i + 1 < n:
            above = freeboard_cells[i + 1]
            cell.m_solid_ab_in[:, :] = np.maximum(above._solid_downflow_rates(), 0.0)
            cell.T_solid_ab_in = float(above.T)
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


def seed_holdup_from_inflows(cell: Cell) -> bool:
    """Seed one holdup state from current transport inflows and ``K`` coefficients.

    For thesis ``holdup_transport`` states, upper bed cells depend on lower-cell
    ``K_auf * m_solid`` transport. During initialization this creates a one-way
    dependency chain: bed0 is known from fresh feed, bed1 depends on bed0, bed2
    depends on bed1, and so on. This helper raises under-filled holdup states
    toward their current steady ``inflow / K`` support without lowering existing
    inventory, so it is safe as an initialization sweep before NR updates.
    """
    if str(cell.solid_state_model) not in {"holdup_transport", "freeboard_closure"}:
        return False
    inflow = np.maximum(
        cell.m_solid_zu
        + cell.m_solid_rez
        + cell.m_solid_in
        + cell.m_solid_auf_in
        + cell.m_solid_ab_in,
        0.0,
    )
    # Eq.2-6 steady holdup must also carry positive local solid formation terms
    # such as char generated from devolatilization. Negative reaction terms remain
    # sinks and are intentionally not used to seed resident inventory.
    inflow = inflow + np.maximum(cell.R_solid, 0.0)
    K_total = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    seeded = np.zeros_like(cell.m_solid)
    mask = K_total > 1e-12
    if not np.any(mask):
        return False
    seeded[mask] = inflow[mask] / K_total[mask]
    underfilled_mask = mask & (np.maximum(cell.m_solid, 0.0) < seeded) & (seeded > 0.0)
    if not np.any(underfilled_mask):
        return False
    cell.m_solid[underfilled_mask] = seeded[underfilled_mask]
    return True


def seed_bed_holdup_chain_from_transport(
    bed_cells: list[Cell],
    *,
    top_above_cell: Cell | None = None,
) -> bool:
    """Propagate thesis bed holdup seeding bottom→top from current transport states.

    Hamel Eq. 2.6 transport uses neighboring solid holdup to define axial inflow
    terms. After switching x0 to active-solid-only semantics, upper bed cells can
    remain empty unless we explicitly walk this dependency chain once during
    initialization. This sweep seeds only empty cells and then refreshes inflow
    terms from the newly created resident holdups.
    """
    if not bed_cells:
        return False
    changed_any = False
    for _ in range(max(3 * len(bed_cells), 1)):
        update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
        changed = False
        for cell in bed_cells:
            changed |= seed_holdup_from_inflows(cell)
        changed_any |= changed
        if not changed:
            break
    update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
    return changed_any


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

    dense_frac = resolve_bottom_gas_inlet_dense_fraction(cell, cfg)
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
        # Reactive component seeding: if above has no VM/moisture but seeding is needed,
        # pull from lower cell. This is a smooth update (no hard switch on solid amount).
        if reactive_inlet_allowed:
            above_reactive_empty = float(
                np.sum(np.maximum(above_solid[:, S_VM] + above_solid[:, S_MOISTURE], 0.0))
            ) <= 1e-12
            if above_reactive_empty:
                above_solid[:, S_VM] = lower_solid[:, S_VM]
                above_solid[:, S_MOISTURE] = lower_solid[:, S_MOISTURE]
        # Temperature: use above cell T only if above supplies any char/ash, else use prev T.
        # This avoids a hard binary switch on total solid that caused Jacobian discontinuity.
        above_char_ash = float(np.sum(np.maximum(above_solid[:, S_CHAR] + above_solid[:, S_ASH], 0.0)))
        above_T = above.T if above_char_ash > 1e-12 else prev.T
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
    _set_degenerate_single_phase_gas_inlet(top, cyclone_cell)
    cyclone_cell.T_in_gas = float(top.T)
    cyclone_cell.m_solid_in.fill(0.0)
    cyclone_cell.m_solid_auf_in[:, :] = propagated_solid_stream_from_cell(top, include_reactive=False)
    cyclone_cell.m_solid_ab_in.fill(0.0)
    cyclone_cell.T_in_solid = float(top.T)
    cyclone_cell.T_solid_auf_in = float(top.T)
    cyclone_cell.T_solid_ab_in = float(cyclone_cell.T)
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
    return_leg_cell.T_solid_auf_in = float(return_leg_cell.T)
    return_leg_cell.T_solid_ab_in = float(cyclone_cell.T)
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
    top_downflow_components: np.ndarray | None = None
    if freeboard_cells:
        top_downflow_matrix = np.maximum(freeboard_cells[0]._solid_downflow_rates(), 0.0)
        top_downflow_total = float(np.sum(top_downflow_matrix))
        top_downflow_components = np.array(
            [
                float(np.sum(top_downflow_matrix[:, S_CHAR])),
                float(np.sum(top_downflow_matrix[:, S_ASH])),
            ],
            dtype=np.float64,
        )
    if not apply_frozen_bed_solid_transport_coefficients(cells):
        update_bed_solid_transport_coefficients(
            cells,
            cfg,
            top_downflow_total=top_downflow_total,
            top_downflow_components=top_downflow_components,
        )
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
    cyclone_cell.N_d[:] = np.maximum(cyclone_cell.N_d_in, 0.0)
    cyclone_cell.N_b.fill(0.0)
    cyclone_cell.T = float(cyclone_cell.T_in_gas)
    side_transport_frozen = apply_frozen_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell)
    if not side_transport_frozen:
        update_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell, cfg)
    seed_holdup_from_inflows(cyclone_cell)
    return_leg_cell.m_solid_ab_in[:, :] = np.maximum(cyclone_cell._solid_downflow_rates(), 0.0)
    return_leg_cell.T_in_solid = float(cyclone_cell.T)
    return_leg_cell.T_solid_ab_in = float(cyclone_cell.T)
    return_leg_cell.T = float(cyclone_cell.T)
    if not side_transport_frozen:
        update_side_block_solid_transport_coefficients(cyclone_cell, return_leg_cell, cfg)
    seed_holdup_from_inflows(return_leg_cell)
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
    if int(changed_cell_idx) == 0:
        for i in range(len(cells)):
            propagate_upstream(cells, cfg, i)
        if not apply_frozen_bed_solid_transport_coefficients(cells):
            update_bed_solid_transport_coefficients(cells, cfg)
        update_bed_solid_transport_inflows(cells)
        for cell in cells:
            project_zero_source_solid_components(cell)
        return
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
