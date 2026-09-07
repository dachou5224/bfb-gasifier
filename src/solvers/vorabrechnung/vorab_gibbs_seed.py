"""Vorabrechnung — vorab gibbs seed segment (OPT-004)."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Dict, List

import numpy as np

from src.core.cell import Cell, S_CHAR, S_VM, S_MOISTURE, S_ASH
from src.core.connectivity import cell_total_solid_holdup
from src.core.species import GAS_SPECIES_INDEX, get_atom_count, gibbs_molar
from src.thermal.devolatilization import (
    devolatilization_rate_for_cell,
    vm_devolatilization_zone_cumulative_fraction,
    vm_devolatilization_zone_increment_fraction,
)
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
    """构造 major-species A1 元素池 b_i（Hamel Anhang A1 (A.2)/(A.3)）。

    (A.2) 元素守恒：求解器内强制 Σ_j a_ij N_j = b_i（见 ``GibbsMinimizer`` / ``ReducedHamelGibbsSolver``）。

    (A.3) 已知源汇 N_Q 修正：本函数在 seed 层把进料与 VM 释放的原子流并入 b_i，
    再交给 A1 Newton 求平衡摩尔数 N_j。对应关系：

    - C ← VM 碳 ``nC_vm``
    - H ← 2·H2O_feed + 2·nH_vm``（进料 H₂O 与 VM 氢）
    - O ← 2·o2_remaining + H2O_feed + 2·nO_vm``（剩余 O₂、H₂O 氧、VM 氧）
    - N ← 2·N2_feed``（进料 N₂）

    ``o2_remaining`` / ``H2O_feed`` / ``N2_feed`` 即 Vorab 侧已知气相源项 N_Q；
    VM 元素来自干燥/热解链，不属于 Gibbs 内未知量。

    Source: Hamel (1999) Anhang A1 (A.2), (A.3); ``vorab_x0`` major-species x0 路径
    """
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
    candidates: list[str] | None = None,
) -> tuple[dict[str, float], dict[str, Any]]:
    """Major-species Gibbs seed with mode dispatch (augmented / reduced / shadow compare)."""
    species = list(candidates or _MAJOR_GIBBS_CANDIDATES)
    mode = _normalize_major_gibbs_solver_mode(solver_mode)
    if mode == "augmented":
        minimizer = _get_major_gibbs_minimizer()
        res, diag = minimizer.solve(
            T=float(T),
            P=float(P),
            elements=elements,
            candidates=species,
            lambda0=lambda0,
            ln_N0=ln_N0,
            return_diag=True,
        )
        diag = dict(diag)
        diag["solver_mode"] = "augmented"
        out = {sp: float(max(res.get(sp, 0.0), 0.0)) for sp in species}
        return out, diag

    if mode == "hamel_reduced":
        reduced = _get_major_gibbs_reduced_solver()
        res, diag = reduced.solve(
            T=float(T),
            P=float(P),
            elements=elements,
            candidates=species,
            lambda0=lambda0,
            tol=1e-8,
            max_iter=120,
            return_diag=True,
        )
        diag = dict(diag)
        diag["solver_mode"] = "hamel_reduced"
        out = {sp: float(max(res.get(sp, 0.0), 0.0)) for sp in species}
        return out, diag

    # shadow_compare: run both, prefer reduced as thesis-target primary output,
    # but keep augmented diagnostics for side-by-side audit.
    reduced = _get_major_gibbs_reduced_solver()
    red_res, red_diag = reduced.solve(
        T=float(T),
        P=float(P),
        elements=elements,
        candidates=species,
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
        candidates=species,
        lambda0=lambda0,
        ln_N0=ln_N0,
        return_diag=True,
    )

    red_out = {sp: float(max(red_res.get(sp, 0.0), 0.0)) for sp in species}
    aug_out = {sp: float(max(aug_res.get(sp, 0.0), 0.0)) for sp in species}
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


# 不含 C₂H₄：§4.3 例举物种，但 Ch.5 无动力学；微量且无可承接路径，不纳入模型
_PYROLYSIS_GIBBS_CANDIDATES: list[str] = ["CO2", "CO", "CH4", "H2", "H2O"]
_PYROLYSIS_GIBBS_ELEM_NAMES: list[str] = ["C", "H", "O"]
_PYROLYSIS_GIBBS_TOL: float = 1e-10
_PYROLYSIS_GIBBS_MAX_ITER: int = 250
# 兼容旧 re-export；Hamel §4.3 路径不再用工程 closure 覆盖
_PYROLYSIS_GIBBS_CLOSURE_REL_TOL: float = 1e-4
_PYROLYSIS_GIBBS_RESIDUAL_REL_TOL: float = 1e-8


def _pyrolysis_elements_closed_system(
    *,
    nC_mol: float,
    nH_mol: float,
    nO_mol: float,
) -> dict[str, float]:
    """Hamel §4.3 / Anhang A1：热解气封闭元素池（仅 VM 释放 C/H/O）。

    Source: Hamel (1999) §4.3 pyrolysis gas composition; Anhang A1
    """
    return {
        "C": float(max(nC_mol, 0.0)),
        "H": float(max(2.0 * nH_mol, 0.0)),
        "O": float(max(2.0 * nO_mol, 0.0)),
    }


def _hamel_teer_gas_pool_strict_an_eq_b(
    *,
    c_gas: float,
    h_gas: float,
    o_gas: float,
) -> tuple[dict[str, float], dict[str, float]]:
    """Koks/Teer 已扣完后的气相元素池：全量进 ``b``，并检验 ``A n = b`` 可行。

    Hamel §4.3：Koks 与 Teer 产量先定，剩余元素全部归气相，再 ``Σ a_ij n_j = b_i``。
    不再把「吃不掉的 C」二次划入 Koks；若当前候选在 ``n≥0`` 下无法吃完全部
    ``(c_gas,h_gas,o_gas)``，直接报错（应回查上游 Koks/Teer 或候选集）。

    Returns
    -------
    elems : dict
        气相元素池（原子流 [mol/s]），等于输入全量
    mole_seed : dict
        某一可行非负基本解（供 min-G 初值）

    Source: Hamel (1999) §4.3 Eq.(4.26)–(4.30)
    """
    from scipy.optimize import linprog

    reduced = _get_major_gibbs_reduced_solver()
    elem_names = _PYROLYSIS_GIBBS_ELEM_NAMES
    candidates = _PYROLYSIS_GIBBS_CANDIDATES
    a = reduced._build_atom_matrix(elem_names, candidates)
    c0 = max(float(c_gas), 0.0)
    h0 = max(float(h_gas), 0.0)
    o0 = max(float(o_gas), 0.0)
    b_full = np.array([c0, h0, o0], dtype=np.float64)

    if float(np.sum(b_full)) <= 1.0e-20:
        z = {e: 0.0 for e in elem_names}
        return z, {sp: 0.0 for sp in candidates}

    eq = linprog(
        np.zeros(len(candidates), dtype=np.float64),
        A_eq=a,
        b_eq=b_full,
        bounds=(0.0, None),
        method="highs",
    )
    if not bool(eq.success):
        raise RuntimeError(
            "Hamel §4.3: Teer/Koks 后气相元素池对候选物种 An=b 不可行 "
            f"(c={c0:.6g}, h={h0:.6g}, o={o0:.6g}); "
            "应回查上游 Koks/Teer 产量或候选集，禁止二次把 C 划回 Koks"
        )
    n_seed = np.asarray(eq.x, dtype=np.float64)
    elems = {elem_names[i]: float(b_full[i]) for i in range(len(elem_names))}
    mole_seed = {candidates[j]: float(max(n_seed[j], 0.0)) for j in range(len(candidates))}
    return elems, mole_seed


def _hamel_gas_pool_with_surplus_c_as_coke(
    *,
    c_gas: float,
    h_gas: float,
    o_gas: float,
) -> tuple[dict[str, float], dict[str, float], float]:
    """已废弃：原「过剩 C→Koks」延伸。现等价于严格 ``An=b``（surplus 恒为 0）。"""
    elems, mole_seed = _hamel_teer_gas_pool_strict_an_eq_b(
        c_gas=c_gas, h_gas=h_gas, o_gas=o_gas
    )
    return elems, mole_seed, 0.0


def _hamel_project_pyrolysis_element_pool_for_a1(
    *,
    c_gas: float,
    h_gas: float,
    o_gas: float,
) -> tuple[dict[str, float], dict[str, float]]:
    """兼容旧名：Teer 后气相全量 ``b`` + 可行 mole_seed（严格 ``An=b``）。"""
    return _hamel_teer_gas_pool_strict_an_eq_b(c_gas=c_gas, h_gas=h_gas, o_gas=o_gas)


def _legacy_fixed_ratio_pyrolysis_increment_not_hamel(
    *,
    dC: float,
    dH: float,
    dO: float,
) -> dict[str, float]:
    """非 Hamel：固定比例热解增量（仅测试/审计对照，生产路径禁用）。"""
    return {
        "CO": 0.40 * dC,
        "CO2": 0.30 * dC,
        "CH4": 0.05 * dC,
        "H2": max(0.15 * dH, 1.0e-12),
        "H2O": max(0.10 * dO, 0.0),
        "O2": 0.0,
        "N2": 0.0,
    }


def _pyrolysis_guess_element_totals(
    guess: dict[str, float],
    *,
    candidates: list[str] | None = None,
) -> dict[str, float]:
    """Gibbs/回退 guess → C/H/O 原子流 [mol/s]。"""
    from src.core.elemental_ledger import gas_element_molar_rates
    from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX

    species = list(candidates or _PYROLYSIS_GIBBS_CANDIDATES)
    rates = np.zeros(len(GAS_SPECIES), dtype=np.float64)
    for sp in species:
        rates[GAS_SPECIES_INDEX[sp]] = float(max(guess.get(sp, 0.0), 0.0))
    return gas_element_molar_rates(rates, fuel_type="coal")


def _pyrolysis_cap_gibbs_guess_to_element_pool(
    guess: dict[str, float],
    elements: dict[str, float],
    *,
    candidates: list[str] | None = None,
) -> dict[str, float]:
    """等比缩 Gibbs 主组分，使 C/H/O 不超出 Teer 分出后的气相元素池。"""
    species = list(candidates or _PYROLYSIS_GIBBS_CANDIDATES)
    base = {sp: float(max(guess.get(sp, 0.0), 0.0)) for sp in species}
    out = _pyrolysis_guess_element_totals(base, candidates=species)
    scale = 1.0
    for el in ("C", "H", "O"):
        budget = float(elements.get(el, 0.0))
        used = float(out.get(el, 0.0))
        if budget <= 1.0e-20:
            if used > 1.0e-12:
                return {sp: 0.0 for sp in species}
            continue
        if used > budget * (1.0 + 1.0e-9):
            scale = min(scale, budget / used)
    if scale >= 1.0 - 1.0e-12:
        return base
    return {sp: base[sp] * scale for sp in species}


def _pyrolysis_guess_respects_element_pool(
    guess: dict[str, float],
    elements: dict[str, float],
    *,
    candidates: list[str] | None = None,
    atol: float = 1.0e-6,
) -> bool:
    """Guess 主组分不得超出封闭 C/H/O 池。"""
    species = list(candidates or _PYROLYSIS_GIBBS_CANDIDATES)
    out = _pyrolysis_guess_element_totals(guess, candidates=species)
    for el in ("C", "H", "O"):
        if float(out.get(el, 0.0)) > float(elements.get(el, 0.0)) + atol:
            return False
    total = float(sum(float(guess.get(sp, 0.0)) for sp in species))
    return total > 1.0e-20


def _solve_pyrolysis_gibbs_an_eq_b(
    *,
    T_K: float,
    P_Pa: float,
    elements: dict[str, float],
    mole_seed: dict[str, float] | None = None,
) -> tuple[dict[str, float], dict[str, Any]]:
    """§4.3 封闭热解气：``min G`` s.t. ``A n = b``, ``n ≥ 0``。

    Hamel §4.3 / Anhang A1 的物理目标是元素守恒下的 Gibbs 平衡。
    现有 lambda-only 降维 Newton 偶发伪收敛且 ``A n ≠ b``；
    此处用约束优化直接强制 ``A n = b``（与论文 (4.30) 一致）。

    Source: Hamel (1999) §4.3 Eq.(4.26)–(4.31); Anhang A1 物理目标
    """
    from scipy.optimize import minimize

    from src.core.constants import P0, Rg

    species = list(_PYROLYSIS_GIBBS_CANDIDATES)
    elem_names = list(_PYROLYSIS_GIBBS_ELEM_NAMES)
    reduced = _get_major_gibbs_reduced_solver()
    a = reduced._build_atom_matrix(elem_names, species)
    b = np.array([float(max(elements.get(e, 0.0), 0.0)) for e in elem_names], dtype=np.float64)
    if float(np.sum(b)) <= 1.0e-20:
        z = {sp: 0.0 for sp in species}
        return z, {
            "converged": True,
            "solver_mode": "an_eq_b_min_G",
            "max_abs_element_closure": 0.0,
            "pyrolysis_gibbs_elements": {e: 0.0 for e in elem_names},
        }

    mu0 = np.array(
        [float(reduced.species_db.get_mu0(sp, float(T_K))) for sp in species],
        dtype=np.float64,
    )
    log_p = float(np.log(float(P_Pa) / P0))
    rt = float(Rg * float(T_K))

    def _gibbs(n_vec: np.ndarray) -> float:
        n_pos = np.maximum(n_vec, 0.0)
        n_tot = float(np.sum(n_pos))
        if n_tot <= 1.0e-30:
            return 0.0
        y = n_pos / n_tot
        return float(n_pos @ (mu0 + rt * (np.log(np.maximum(y, 1.0e-300)) + log_p)))

    if mole_seed is None:
        x0 = np.full(len(species), 1.0e-12, dtype=np.float64)
    else:
        x0 = np.array(
            [float(max(mole_seed.get(sp, 0.0), 0.0)) for sp in species],
            dtype=np.float64,
        )
        if float(np.sum(x0)) <= 1.0e-20:
            x0 = np.full(len(species), 1.0e-12, dtype=np.float64)
    x0 = np.maximum(x0, 1.0e-14)

    cons = [
        {
            "type": "eq",
            "fun": lambda n, i=i: float(a[i] @ np.maximum(n, 0.0) - b[i]),
        }
        for i in range(len(elem_names))
    ]
    opt = minimize(
        _gibbs,
        x0,
        method="SLSQP",
        bounds=[(0.0, None)] * len(species),
        constraints=cons,
        options={"ftol": float(_PYROLYSIS_GIBBS_TOL), "maxiter": int(_PYROLYSIS_GIBBS_MAX_ITER), "disp": False},
    )
    n_sol = np.maximum(np.asarray(opt.x, dtype=np.float64), 0.0)
    closure = a @ n_sol - b
    max_abs = float(np.max(np.abs(closure), initial=0.0))

    # SLSQP 偶发 success=False 但解已满足 An=b；回退到可行 mole_seed
    if max_abs > 1.0e-6:
        n_seed = np.maximum(x0.copy(), 0.0)
        clos_seed = a @ n_seed - b
        if float(np.max(np.abs(clos_seed), initial=0.0)) <= 1.0e-8:
            n_sol = n_seed
            closure = clos_seed
            max_abs = float(np.max(np.abs(closure), initial=0.0))
        else:
            # 再试：均匀扰动初值
            for scale in (0.5, 2.0, 0.1):
                x_try = np.maximum(x0 * scale, 1.0e-14)
                opt2 = minimize(
                    _gibbs,
                    x_try,
                    method="SLSQP",
                    bounds=[(0.0, None)] * len(species),
                    constraints=cons,
                    options={
                        "ftol": float(_PYROLYSIS_GIBBS_TOL),
                        "maxiter": int(_PYROLYSIS_GIBBS_MAX_ITER),
                        "disp": False,
                    },
                )
                n2 = np.maximum(np.asarray(opt2.x, dtype=np.float64), 0.0)
                c2 = a @ n2 - b
                m2 = float(np.max(np.abs(c2), initial=0.0))
                if m2 < max_abs:
                    n_sol, closure, max_abs, opt = n2, c2, m2, opt2
                if max_abs <= 1.0e-6:
                    break

    out = {species[j]: float(n_sol[j]) for j in range(len(species))}
    diag: dict[str, Any] = {
        "converged": max_abs <= 1.0e-6,
        "converged_by": "an_eq_b_min_G",
        "solver_mode": "an_eq_b_min_G",
        "message": str(getattr(opt, "message", "")),
        "opt_success": bool(getattr(opt, "success", False)),
        "max_abs_element_closure": max_abs,
        "element_closure": {elem_names[i]: float(closure[i]) for i in range(len(elem_names))},
        "pyrolysis_gibbs_elements": {elem_names[i]: float(b[i]) for i in range(len(elem_names))},
        "G": float(_gibbs(n_sol)),
    }
    return out, diag


def _solve_pyrolysis_gibbs_reduced(
    *,
    T_K: float,
    P_Pa: float,
    elements: dict[str, float],
    mole_seed: dict[str, float] | None = None,
) -> tuple[dict[str, float], dict[str, Any]]:
    """§4.3 热解池：严格 ``A n = b`` 下 ``min G``（见 ``_solve_pyrolysis_gibbs_an_eq_b``）。

    Source: Hamel (1999) §4.3, Anhang A1 物理目标
    """
    return _solve_pyrolysis_gibbs_an_eq_b(
        T_K=T_K,
        P_Pa=P_Pa,
        elements=elements,
        mole_seed=mole_seed,
    )


def solve_pyrolysis_gas_gibbs_composition(
    *,
    T_K: float,
    P_Pa: float,
    nC_mol: float,
    nH_mol: float,
    nO_mol: float,
    solver_mode: str = "hamel_reduced",
    pyrolysis_gibbs_temperature_floor_K: float | None = None,
) -> dict[str, float]:
    """Hamel §4.3：封闭热解气包 Gibbs 组成 [mol/s]（不含 O₂/N₂ 进料）。

    Source: Hamel (1999) §4.3 Eq.(4.26)–(4.31), Anhang A1

    ``nC_mol``：碳原子摩尔流 [mol/s]。
    ``nH_mol`` / ``nO_mol``：氢/氧 **半原子** 摩尔流 [mol/s]（即 H、O 原子流 ÷ 2），
    与 Vorab ``nH_vm_mol_s`` / ``nO_vm_mol_s`` 及 ``vm_element_molar_rates`` 的 H/O ÷2 口径一致。

    候选物种：CO, CO₂, CH₄, H₂, H₂O（NH₃/H₂S 在 Teer 分出前处理）。
    §4.3 文中 z.B. 曾列 C₂H₄，但 Ch.5 无动力学承接，模型不纳入。
    输入 ``nC/nH/nO`` 须已是 **Koks+Teer 扣完后** 的气相元素；本函数对全量 ``b``
    做严格 ``A n = b``（元素用完），不再二次把 C 划回 Koks。
    ``pyrolysis_gibbs_temperature_floor_K`` 已废弃（Hamel 无温度地板），保留仅兼容旧调用。
    """
    _ = (solver_mode, pyrolysis_gibbs_temperature_floor_K)
    dC = float(max(nC_mol, 0.0))
    dH = float(max(nH_mol, 0.0))
    dO = float(max(nO_mol, 0.0))
    if dC <= 1.0e-20 and dH <= 1.0e-20 and dO <= 1.0e-20:
        return {sp: 0.0 for sp in (*_PYROLYSIS_GIBBS_CANDIDATES, "O2", "N2")}

    elems_gas = _pyrolysis_elements_closed_system(nC_mol=dC, nH_mol=dH, nO_mol=dO)
    elems_a1, mole_seed = _hamel_teer_gas_pool_strict_an_eq_b(
        c_gas=dC,
        h_gas=2.0 * dH,
        o_gas=2.0 * dO,
    )
    if sum(elems_a1.values()) <= 1.0e-20:
        return {sp: 0.0 for sp in (*_PYROLYSIS_GIBBS_CANDIDATES, "O2", "N2")}

    # Hamel §4.3：在给定局部 T,P 下求平衡（无温度地板/升档续算）
    guess, diag = _solve_pyrolysis_gibbs_reduced(
        T_K=float(T_K),
        P_Pa=float(P_Pa),
        elements=elems_a1,
        mole_seed=mole_seed,
    )
    if not bool(diag.get("converged", False)):
        raise RuntimeError("Hamel §4.3 pyrolysis Gibbs (An=b) did not converge")
    if float(diag.get("max_abs_element_closure", np.inf)) > 1.0e-4:
        raise RuntimeError("Hamel §4.3 pyrolysis Gibbs element closure too loose")

    # 全量元素须用完：禁止再 cap 缩产物
    used = _pyrolysis_guess_element_totals(guess)
    for el in ("C", "H", "O"):
        budget = float(elems_gas.get(el, 0.0))
        got = float(used.get(el, 0.0))
        if abs(got - budget) > max(1.0e-6, 1.0e-6 * max(budget, 1.0)):
            raise RuntimeError(
                f"Hamel §4.3 pyrolysis Gibbs did not consume all {el}: "
                f"used={got:.6g} budget={budget:.6g}"
            )
    if not _pyrolysis_guess_respects_element_pool(guess, elems_gas):
        raise RuntimeError("Hamel §4.3 pyrolysis Gibbs exceeds Teer gas element pool")

    out = {sp: 0.0 for sp in ("O2", "N2")}
    for sp in _PYROLYSIS_GIBBS_CANDIDATES:
        out[sp] = float(max(guess.get(sp, 0.0), 0.0))
    return out

