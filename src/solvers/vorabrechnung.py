"""Vorabrechnung（预算）：在全局 NR 前生成与轴向温度一致的化学初值 x₀。

对应 Hamel (1999) Module 2 / §2.2 思想：先估计轴向温度与热解进度，再构造元素大致平衡的
气相摩尔流初值，降低 F(x₀) 与内层 NR 难度。

注：此为工程近似，非论文 Fortran 逐行复现。
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from src.core.cell import Cell, S_CHAR, S_VM, S_MOISTURE, S_ASH
from src.core.connectivity import cell_total_solid_holdup
from src.core.species import GAS_SPECIES_INDEX, get_atom_count, gibbs_molar
from src.thermal.devolatilization import devolatilization_rate_for_cell
from src.thermal.drying import solve_drying_CN
from src.thermodynamics.gibbs_hamel_reduced import ReducedHamelGibbsSolver
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer


class _MajorSpeciesDBAdapter:
    """Thermo adapter for major-species Gibbs initialization in Vorabrechnung."""

    def get_mu0(self, species: str, T: float) -> float:
        return gibbs_molar(species, T)

    def get_atom_count(self, species: str, element: str) -> int:
        return get_atom_count(species, element)


_MAJOR_GIBBS_MINIMIZER: GibbsMinimizer | None = None
_MAJOR_GIBBS_REDUCED_SOLVER: ReducedHamelGibbsSolver | None = None
_MAJOR_GIBBS_CANDIDATES: list[str] = ["CO2", "CO", "CH4", "H2", "H2O", "O2", "N2"]
_MAJOR_GIBBS_WARMSTART_MAX_RESIDUAL = 1e-8
_MAJOR_GIBBS_WARMSTART_MAX_CLOSURE_REL = 1e-6


def _get_major_gibbs_minimizer() -> GibbsMinimizer:
    global _MAJOR_GIBBS_MINIMIZER
    if _MAJOR_GIBBS_MINIMIZER is None:
        _MAJOR_GIBBS_MINIMIZER = GibbsMinimizer(_MajorSpeciesDBAdapter())
    return _MAJOR_GIBBS_MINIMIZER


def _get_major_gibbs_reduced_solver() -> ReducedHamelGibbsSolver:
    global _MAJOR_GIBBS_REDUCED_SOLVER
    if _MAJOR_GIBBS_REDUCED_SOLVER is None:
        _MAJOR_GIBBS_REDUCED_SOLVER = ReducedHamelGibbsSolver(_MajorSpeciesDBAdapter())
    return _MAJOR_GIBBS_REDUCED_SOLVER


def _legacy_major_seed(
    *,
    nC_vm: float,
    nH_vm: float,
    nO_vm: float,
    H2O_feed: float,
    N2_feed: float,
    o2_remaining: float,
    frac_height: float,
) -> dict[str, float]:
    """Legacy heuristic major-species seed (kept as robust fallback)."""
    co2_frac = max(0.0, 1.0 - frac_height)
    co_frac = min(1.0, frac_height + 0.2)

    n_CO2 = nC_vm * co2_frac * 0.5
    n_CO = nC_vm * co_frac * 0.5
    n_CH4 = nC_vm * 0.05
    n_H2 = max(nH_vm - 2.0 * (nC_vm * co2_frac * 0.1), 0.01 * max(nH_vm, 1e-12))
    n_H2O = H2O_feed + (o2_remaining * 0.0) - nO_vm * 0.5
    n_H2O = max(n_H2O, 0.01 * H2O_feed)
    n_N2 = N2_feed
    n_O2 = max(o2_remaining, 0.0)
    return {
        "CO2": float(max(n_CO2, 0.0)),
        "CO": float(max(n_CO, 0.0)),
        "CH4": float(max(n_CH4, 0.0)),
        "H2": float(max(n_H2, 0.0)),
        "H2O": float(max(n_H2O, 0.0)),
        "N2": float(max(n_N2, 0.0)),
        "O2": float(max(n_O2, 0.0)),
    }


def _major_elements_from_feeds(
    *,
    nC_vm: float,
    nH_vm: float,
    nO_vm: float,
    H2O_feed: float,
    N2_feed: float,
    o2_remaining: float,
) -> dict[str, float]:
    """Build major-element atom pool (C/H/O/N) for A1-style Gibbs x0."""
    return {
        "C": float(max(nC_vm, 0.0)),
        "H": float(max(2.0 * H2O_feed + 2.0 * nH_vm, 0.0)),
        "O": float(max(2.0 * o2_remaining + H2O_feed + 2.0 * nO_vm, 0.0)),
        "N": float(max(2.0 * N2_feed, 0.0)),
    }


def _normalize_major_gibbs_solver_mode(mode: str | None) -> str:
    """Normalize major-gibbs solver mode to one of augmented/hamel_reduced/shadow_compare."""
    m = str(mode or "augmented").strip().lower()
    aliases = {
        "reduced": "hamel_reduced",
        "hamel_reduced": "hamel_reduced",
        "augmented": "augmented",
        "aug": "augmented",
        "shadow_compare": "shadow_compare",
        "shadow": "shadow_compare",
    }
    out = aliases.get(m)
    if out is None:
        raise ValueError(
            f"Unsupported major_gibbs_solver_mode={mode!r}; "
            "expected one of 'augmented', 'hamel_reduced', 'shadow_compare'"
        )
    return out


def _solve_major_gibbs_seed(
    *,
    T: float,
    P: float,
    elements: dict[str, float],
    lambda0: np.ndarray | None,
    ln_N0: float | None,
    solver_mode: str = "augmented",
) -> tuple[dict[str, float], dict[str, Any]]:
    """Major-species Gibbs seed with mode dispatch (augmented / reduced / shadow compare)."""
    mode = _normalize_major_gibbs_solver_mode(solver_mode)
    if mode == "augmented":
        minimizer = _get_major_gibbs_minimizer()
        res, diag = minimizer.solve(
            T=float(T),
            P=float(P),
            elements=elements,
            candidates=_MAJOR_GIBBS_CANDIDATES,
            lambda0=lambda0,
            ln_N0=ln_N0,
            return_diag=True,
        )
        diag = dict(diag)
        diag["solver_mode"] = "augmented"
        out = {sp: float(max(res.get(sp, 0.0), 0.0)) for sp in _MAJOR_GIBBS_CANDIDATES}
        return out, diag

    if mode == "hamel_reduced":
        reduced = _get_major_gibbs_reduced_solver()
        res, diag = reduced.solve(
            T=float(T),
            P=float(P),
            elements=elements,
            candidates=_MAJOR_GIBBS_CANDIDATES,
            lambda0=lambda0,
            tol=1e-8,
            max_iter=120,
            return_diag=True,
        )
        diag = dict(diag)
        diag["solver_mode"] = "hamel_reduced"
        out = {sp: float(max(res.get(sp, 0.0), 0.0)) for sp in _MAJOR_GIBBS_CANDIDATES}
        return out, diag

    # shadow_compare: run both, prefer reduced as thesis-target primary output,
    # but keep augmented diagnostics for side-by-side audit.
    reduced = _get_major_gibbs_reduced_solver()
    red_res, red_diag = reduced.solve(
        T=float(T),
        P=float(P),
        elements=elements,
        candidates=_MAJOR_GIBBS_CANDIDATES,
        lambda0=lambda0,
        tol=1e-8,
        max_iter=120,
        return_diag=True,
    )
    minimizer = _get_major_gibbs_minimizer()
    aug_res, aug_diag = minimizer.solve(
        T=float(T),
        P=float(P),
        elements=elements,
        candidates=_MAJOR_GIBBS_CANDIDATES,
        lambda0=lambda0,
        ln_N0=ln_N0,
        return_diag=True,
    )

    red_out = {sp: float(max(red_res.get(sp, 0.0), 0.0)) for sp in _MAJOR_GIBBS_CANDIDATES}
    aug_out = {sp: float(max(aug_res.get(sp, 0.0), 0.0)) for sp in _MAJOR_GIBBS_CANDIDATES}
    red_credible = _major_gibbs_seed_is_credible(guess=red_out, diag=red_diag, elements=elements)
    selected = "hamel_reduced" if red_credible else "augmented"
    primary_out = red_out if selected == "hamel_reduced" else aug_out
    primary_diag = dict(red_diag if selected == "hamel_reduced" else aug_diag)
    primary_diag["solver_mode"] = "shadow_compare"
    primary_diag["selected_solver"] = selected
    primary_diag["shadow_compare"] = {
        "hamel_reduced": {
            "converged": bool(red_diag.get("converged", False)),
            "final_residual": float(red_diag.get("final_residual", np.inf)),
            "sum_guess": float(np.sum(list(red_out.values()))),
        },
        "augmented": {
            "converged": bool(aug_diag.get("converged", False)),
            "final_residual": float(aug_diag.get("final_residual", np.inf)),
            "sum_guess": float(np.sum(list(aug_out.values()))),
        },
    }
    return primary_out, primary_diag


def _major_gibbs_seed_is_credible(
    *,
    guess: dict[str, float],
    diag: dict[str, Any],
    elements: dict[str, float],
) -> bool:
    """Validate that a major-species Gibbs x0 looks like a usable molar-flow seed."""
    vals = np.array([float(guess.get(sp, 0.0)) for sp in _MAJOR_GIBBS_CANDIDATES], dtype=np.float64)
    if not np.all(np.isfinite(vals)) or np.any(vals < 0.0):
        return False
    if not bool(diag.get("converged", False)):
        return False
    total_guess = float(np.sum(vals))
    total_elements = max(float(np.sum([max(v, 0.0) for v in elements.values()])), 1e-12)
    return total_guess <= 10.0 * total_elements


def _major_gibbs_diag_allows_warmstart(diag: dict[str, Any]) -> bool:
    """Allow axial warm-start only from well-converged previous cell solutions."""
    if not bool(diag.get("converged", False)):
        return False
    residual = diag.get("final_residual")
    if residual is None:
        return False
    try:
        residual_f = float(residual)
    except (TypeError, ValueError):
        return False
    if not (np.isfinite(residual_f) and residual_f <= _MAJOR_GIBBS_WARMSTART_MAX_RESIDUAL):
        return False

    solver_mode = str(diag.get("solver_mode", "")).strip().lower()
    if solver_mode == "shadow_compare":
        selected = str(diag.get("selected_solver", "")).strip().lower()
        if selected != "hamel_reduced":
            return False

    if str(diag.get("converged_by", "")).strip().lower() == "closure_relaxed":
        return False

    closure_rel = diag.get("closure_rel")
    if closure_rel is not None:
        try:
            closure_rel_f = float(closure_rel)
        except (TypeError, ValueError):
            return False
        if not (np.isfinite(closure_rel_f) and closure_rel_f <= _MAJOR_GIBBS_WARMSTART_MAX_CLOSURE_REL):
            return False
    return True


def estimate_axial_T_profile(
    n_cells: int,
    T_inlet: float,
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    C_dry: float,
    H_dry: float,
    moisture_wt: float,
    P: float,
) -> np.ndarray:
    """估计轴向温度剖面（底高顶低），用于 Vorabrechnung 初值。

    简化绝热温升 + 线性轴向衰减；顶温做下限使剖面落在常见床层/出口区间。
    """
    moisture_frac = moisture_wt / 100.0
    dry_fuel = fuel_feed_kg_s * (1.0 - moisture_frac)

    M_C = 12.011e-3
    n_C = dry_fuel * (C_dry / 100.0) / M_C

    dH_comb = 393_500.0  # J/mol O2（量级近似）
    total_mol_s = O2_feed + H2O_feed + N2_feed + max(n_C * 0.5, 1e-6)
    Cp_mix = 32.0  # J/(mol·K)

    dT_adiabatic = (O2_feed * dH_comb) / max(total_mol_s * Cp_mix, 1.0)
    # LU/HTW 稳定分支通常落在约 1100–1250 K；过热的 Vorabrechnung 初值
    # 会把 GS/NR 更容易推向下部热支。
    T_bottom = float(np.clip(T_inlet + dT_adiabatic, 1050.0, 1325.0))
    T_top = float(np.clip(T_bottom * 0.86, 920.0, 1180.0))

    if n_cells <= 1:
        return np.array([T_bottom], dtype=np.float64)

    T_profile = np.array(
        [
            T_bottom + (T_top - T_bottom) * (i / max(n_cells - 1, 1))
            for i in range(n_cells)
        ],
        dtype=np.float64,
    )
    return T_profile


def vorabrechnung_tau_for_cell(cell: Cell) -> float:
    """Vorabrechnung 使用的固相停留时间尺度 [s]。

    使用全床高度除以 u_mf 作为固相停留时间估计，确保干燥和热解在
    Vorabrechnung 预计算时接近完全（Hamel 单次预算语义）。
    非均匀床层网格必须使用显式 ``_vorab_bed_height``；``n_cells * dh``
    仅作为旧 uniform-mesh fallback。
    """
    bed_height = getattr(cell, "_vorab_bed_height", None)
    if bed_height is None:
        n = int(getattr(cell, "_n_vorab_cells", 1))
        bed_height = float(cell.geo.dh) * max(n, 1)
    return float(max(float(bed_height), float(cell.geo.dh)) / max(float(cell.u_mf), 1e-3))


def refresh_cell_vorabrechnung(cell: Cell, *, force: bool = False) -> None:
    """按 Hamel 外层 Abgleich 语义刷新单格 Vorabrechnung。

    顺序约束：
    1. 先重算 hydrodynamics（u_mf, u_b, eps_b, K_bd...）
    2. 再按最新流体力学停留时间计算 drying/pyrolysis 缓存
    3. inner NR 期间只读这些缓存，不在残差扰动中重复重算
    """
    cell.calc_hydrodynamics()
    if force:
        cell.invalidate_vorabrechnung_cache()
    cell.compute_vorabrechnung(vorabrechnung_tau_for_cell(cell))


def refresh_cell_hydrodynamics_for_frozen_inner(cell: Cell) -> None:
    """Refresh only hydrodynamics/transport and snapshot for inner-NR freeze.

    This path intentionally keeps already prepared drying/pyrolysis source caches intact.
    """
    cell.calc_hydrodynamics()
    cell._snapshot_vorabrechnung_hydrodynamics()


def initialize_fixed_vorabrechnung_sources_for_cells(
    cells: List[Cell],
    *,
    T_reference: float,
) -> None:
    """Initialize fixed drying/pyrolysis Vorabrechnung sources at a single reference temperature.

    Hamel-style strict mode: drying/DAEM are budgeted once (Vorabrechnung) and then
    treated as fixed sources for subsequent outer/inner iterations.
    """
    T_ref = float(max(T_reference, 1.0))
    for cell in cells:
        T_saved = float(cell.T)
        cell.T = T_ref
        cell.invalidate_vorabrechnung_cache()
        refresh_cell_vorabrechnung(cell, force=False)
        cell.T = T_saved
        refresh_cell_hydrodynamics_for_frozen_inner(cell)


def refresh_vorabrechnung_for_cells(
    cells: List[Cell],
    *,
    force: bool = False,
    refresh_sources: bool = True,
) -> None:
    """批量刷新 Vorabrechnung（outer Abgleich 调用入口）。

    Parameters
    ----------
    refresh_sources:
        True: refresh hydrodynamics + drying/pyrolysis sources.
        False: refresh hydrodynamics/snapshots only, keep drying/pyrolysis fixed.
    """
    if not refresh_sources:
        for cell in cells:
            refresh_cell_hydrodynamics_for_frozen_inner(cell)
        return
    for cell in cells:
        refresh_cell_vorabrechnung(cell, force=force)


def generate_initial_x0(
    cells: List[Cell],
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    moisture_wt: float,
    ash_dry_wt: float,
    VM_daf: float,
    T_profile: np.ndarray,
    fuel_type: str = "coal",
    use_hamel_major_gibbs_x0: bool = True,
    strict_hamel_major_gibbs_x0: bool = False,
    major_gibbs_solver_mode: str = "augmented",
) -> None:
    """按轴向高度与估计热解进度，写入各 cell 的 N_d/N_b/T/m_solid 初值（就地修改）。

    当 ``use_hamel_major_gibbs_x0=True`` 时，主组分初值采用 A1 风格 Gibbs 最小化
    （带 lambda/lnN warm-start；可选 augmented/reduced/shadow_compare）。
    - 非 strict：失败时自动回退到 legacy 启发式分配（兼容旧流程）。
    - strict：失败即抛错，避免 thesis 主线静默 fallback。
    """
    n = len(cells)
    idx = GAS_SPECIES_INDEX

    total_inlet = O2_feed + H2O_feed + N2_feed

    moisture_frac = moisture_wt / 100.0
    ash_frac = ash_dry_wt / 100.0
    vm_daf_frac = VM_daf / 100.0
    dry_feed = fuel_feed_kg_s * (1.0 - moisture_frac)
    daf_feed = dry_feed * (1.0 - ash_frac)

    M_C, M_H, M_O = 12.011e-3, 1.00794e-3, 15.999e-3
    c_dry = C_dry / 100.0
    h_dry = H_dry / 100.0
    o_dry = O_dry / 100.0
    to_daf = 1.0 / max(1.0 - ash_frac, 1e-9)
    c_daf = c_dry * to_daf
    h_daf = h_dry * to_daf
    o_daf = o_dry * to_daf

    w_vm_dry = vm_daf_frac * (1.0 - ash_frac)
    w_char_dry = 1.0 - ash_frac - w_vm_dry

    daem_fuel = "brown_coal"
    if fuel_type in ("wood", "biomass"):
        daem_fuel = "wood"
    elif fuel_type in ("coal", "brown_coal"):
        daem_fuel = "brown_coal"

    major_lambda_ws: np.ndarray | None = None
    major_lnN_ws: float | None = None

    for i, cell in enumerate(cells):
        T_i = float(T_profile[i])
        frac_height = (i + 0.5) / max(n, 1)
        solid_T_init = max(float(getattr(cell, "T_in_solid", T_profile[0])), 293.15)

        o2_consumed_frac = min(1.0 - np.exp(-5.0 * frac_height), 1.0)
        o2_remaining = O2_feed * max(1.0 - o2_consumed_frac, 0.0)

        tau_est = cell.geo.dh / max(0.05, cell.u_mf if cell.u_mf > 0 else 0.05)
        X_vm, _ = devolatilization_rate_for_cell(
            T_bed=T_i,
            tau_cell=tau_est,
            T_init=solid_T_init,
            VM_daf=vm_daf_frac,
            fuel_type=daem_fuel,
        )

        X_vm_cum = min(X_vm * (frac_height + 0.1), 1.0)
        m_vm_released = daf_feed * vm_daf_frac * X_vm_cum

        nC_vm = m_vm_released * c_daf / M_C
        nH_vm = m_vm_released * h_daf / M_H / 2.0
        nO_vm = m_vm_released * o_daf / M_O / 2.0

        major_guess: dict[str, float]
        if use_hamel_major_gibbs_x0:
            try:
                elems = _major_elements_from_feeds(
                    nC_vm=float(nC_vm),
                    nH_vm=float(nH_vm),
                    nO_vm=float(nO_vm),
                    H2O_feed=float(H2O_feed),
                    N2_feed=float(N2_feed),
                    o2_remaining=float(o2_remaining),
                )
                major_guess, diag = _solve_major_gibbs_seed(
                    T=T_i,
                    P=float(cell.P),
                    elements=elems,
                    lambda0=major_lambda_ws,
                    ln_N0=major_lnN_ws,
                    solver_mode=major_gibbs_solver_mode,
                )
                if not _major_gibbs_seed_is_credible(guess=major_guess, diag=diag, elements=elems):
                    if strict_hamel_major_gibbs_x0:
                        raise RuntimeError(
                            "major Gibbs x0 produced non-credible seed in strict mode "
                            f"(converged={bool(diag.get('converged', False))}, "
                            f"residual={diag.get('final_residual')})"
                        )
                    else:
                        raise RuntimeError("major Gibbs x0 produced non-credible molar-flow seed")
                if _major_gibbs_diag_allows_warmstart(diag):
                    lam_diag = diag.get("lambda")
                    if isinstance(lam_diag, np.ndarray) and np.all(np.isfinite(lam_diag)):
                        major_lambda_ws = np.array(lam_diag, dtype=np.float64)
                    else:
                        major_lambda_ws = None
                    ln_diag = diag.get("ln_N")
                    if isinstance(ln_diag, float) and np.isfinite(ln_diag):
                        major_lnN_ws = float(ln_diag)
                    else:
                        major_lnN_ws = None
                else:
                    major_lambda_ws = None
                    major_lnN_ws = None
            except Exception:
                if strict_hamel_major_gibbs_x0:
                    raise
                # Failed/low-quality cell solution must not poison downstream axial warm-start.
                major_lambda_ws = None
                major_lnN_ws = None
                major_guess = _legacy_major_seed(
                    nC_vm=float(nC_vm),
                    nH_vm=float(nH_vm),
                    nO_vm=float(nO_vm),
                    H2O_feed=float(H2O_feed),
                    N2_feed=float(N2_feed),
                    o2_remaining=float(o2_remaining),
                    frac_height=float(frac_height),
                )
        else:
            major_guess = _legacy_major_seed(
                nC_vm=float(nC_vm),
                nH_vm=float(nH_vm),
                nO_vm=float(nO_vm),
                H2O_feed=float(H2O_feed),
                N2_feed=float(N2_feed),
                o2_remaining=float(o2_remaining),
                frac_height=float(frac_height),
            )

        n_CO2 = float(major_guess["CO2"])
        n_CO = float(major_guess["CO"])
        n_CH4 = float(major_guess["CH4"])
        n_H2 = float(major_guess["H2"])
        n_H2O = float(major_guess["H2O"])
        n_N2 = float(major_guess["N2"])
        n_O2 = float(major_guess["O2"])

        total = n_CO2 + n_CO + n_CH4 + n_H2 + n_H2O + n_N2 + n_O2
        total = max(total, total_inlet * 0.5)

        cell.N_d.fill(0.0)
        cell.N_b.fill(0.0)
        # Initialize gas split using initialized hydrodynamics instead of fixed 70/30.
        # Bubble-gas share is proportional to visible gas holdup in bubble vs dense phase.
        eps_b = float(getattr(cell, "eps_b", np.nan))
        eps_d_void = float(getattr(cell, "eps_d_voidage", np.nan))
        if np.isfinite(eps_b) and np.isfinite(eps_d_void):
            gas_holdup_total = max(eps_b + eps_d_void, 1e-12)
            bubble_share = float(np.clip(eps_b / gas_holdup_total, 0.05, 0.95))
        else:
            bubble_share = 0.30
        dense_share = float(1.0 - bubble_share)

        for sp, val in [
            ("CO2", n_CO2),
            ("CO", n_CO),
            ("CH4", n_CH4),
            ("H2", n_H2),
            ("H2O", n_H2O),
            ("N2", n_N2),
            ("O2", n_O2),
        ]:
            cell.N_d[idx[sp]] = max(val * dense_share, 1e-12)
            cell.N_b[idx[sp]] = max(val * bubble_share, 1e-12)

        if i == 0:
            bottom_dense_frac = float(getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", np.nan))
            if np.isfinite(bottom_dense_frac):
                bottom_dense_frac = float(np.clip(bottom_dense_frac, 0.0, 1.0))
                for sp in ("O2", "H2O", "N2"):
                    j = idx[sp]
                    total_sp = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                    cell.N_d[j] = max(total_sp * bottom_dense_frac, 1e-12)
                    cell.N_b[j] = max(total_sp * (1.0 - bottom_dense_frac), 1e-12)

        cell.T = T_i
        # 初始化 4 组分固相猜测 [char, vm, moisture, ash]
        nk = cell.solid.n_size_classes
        m_total_est = fuel_feed_kg_s * (0.3 + 0.7 * (1.0 - frac_height)) * 0.5
        if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
            cell.calc_hydrodynamics()
            m_total_est = max(cell_total_solid_holdup(cell), 1e-12)
            dry_diag = solve_drying_CN(
                cell.solid.d_p,
                T_i,
                solid_T_init,
                moisture_wt,
                max(float(tau_est), 0.05),
                Nr=12,
                Nt=80,
                pressure_pa=float(cell.P),
                return_history=True,
            )
            x_dry_cum = min(float(dry_diag["X_dry"][-1]) * (frac_height + 0.1), 1.0)
            active_seed = np.maximum(cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in, 0.0)
            char_seed = float(np.sum(active_seed[:, S_CHAR]))
            vm_seed = float(np.sum(active_seed[:, S_VM]))
            moist_seed = float(np.sum(active_seed[:, S_MOISTURE]))
            ash_seed = float(np.sum(active_seed[:, S_ASH]))

            def _class_fraction(comp_idx: int) -> np.ndarray:
                col = np.maximum(active_seed[:, comp_idx], 0.0)
                total = float(np.sum(col))
                if total <= 1e-12:
                    return np.full(nk, 1.0 / max(nk, 1), dtype=np.float64)
                return col / total

            char_frac = _class_fraction(S_CHAR)
            vm_frac = _class_fraction(S_VM)
            moist_frac = _class_fraction(S_MOISTURE)
            ash_frac_classes = _class_fraction(S_ASH)

            # Hamel distinguishes active fuel solids from the inert bed inventory.
            # Without an explicit inert-material state in ``m_solid``, seed only the
            # active fuel inventory from mapped solid inflows; keep total bed holdup
            # for hydrodynamics / transport coefficients only.
            moist_mass = max(moist_seed * (1.0 - x_dry_cum), 0.0)
            vm_mass = max(vm_seed * (1.0 - X_vm_cum), 0.0)
            ash_mass = max(ash_seed, 0.0)
            char_mass = max(char_seed, 0.0)

            cell.m_solid[:, S_CHAR] = char_mass * char_frac
            cell.m_solid[:, S_VM] = vm_mass * vm_frac
            cell.m_solid[:, S_MOISTURE] = moist_mass * moist_frac
            cell.m_solid[:, S_ASH] = ash_mass * ash_frac_classes
        else:
            cell.m_solid[:, S_CHAR] = max(m_total_est * (1.0 - moisture_frac) * w_char_dry / nk, 1e-12)
            cell.m_solid[:, S_VM] = max(m_total_est * (1.0 - moisture_frac) * w_vm_dry / nk, 1e-12)
            cell.m_solid[:, S_MOISTURE] = max(m_total_est * moisture_frac / nk, 1e-12)
            cell.m_solid[:, S_ASH] = max(m_total_est * (1.0 - moisture_frac) * ash_frac / nk, 1e-12)
