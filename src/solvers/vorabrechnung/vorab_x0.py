"""Vorabrechnung — vorab x0 segment (OPT-004)."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Dict, List

import numpy as np

from src.core.cell import Cell, S_CHAR, S_VM, S_MOISTURE, S_ASH
from src.core.connectivity import cell_total_solid_holdup, propagate_upstream
from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX, get_atom_count, gibbs_molar
from src.kinetics.char_reactions import R4_M_C, phi_c
from src.thermal.devolatilization import (
    devolatilization_rate_for_cell,
    vm_devolatilization_zone_cumulative_fraction,
    vm_devolatilization_zone_increment_fraction,
)
from src.thermal.drying import solve_drying_CN
from src.thermodynamics.gibbs_hamel_reduced import ReducedHamelGibbsSolver
from src.thermodynamics.gibbs_minimizer import GibbsMinimizer
from src.solvers.vorabrechnung.vorab_budgets import (
    VorabrechnungCellBudget,
    _bed_cell_indices,
    _resolve_product_bed_cell_indices,
    build_vorabrechnung_cell_budgets,
    normalize_x0_holdup_mode,
    vorabrechnung_tau_for_cell,
)
from src.solvers.vorabrechnung.vorab_cell_mapping import (
    _hydrodynamic_bubble_share_for_cell,
    _incremental_pyrolysis_major_mol_s,
)
from src.solvers.vorabrechnung.vorab_gibbs_seed import (
    _legacy_major_seed,
    _major_elements_from_feeds,
    _major_gibbs_diag_allows_warmstart,
    _major_gibbs_seed_is_credible,
    _solve_major_gibbs_seed,
    solve_pyrolysis_gas_gibbs_composition,
)

PYROLYSIS_DENSE_SEED_SPECIES: tuple[str, ...] = ("NH3", "TAR1", "TAR2")
MAJOR_PRODUCT_SEED_SPECIES: tuple[str, ...] = ("CO", "CO2", "H2", "CH4")
_VORAB_MAJOR_PRODUCT_SOURCE_FLOOR_MOL_S = 1.0e-3


def _major_product_init_gas_source_mol_s(cell: Cell, species_index: int) -> float:
    """Init 产物 holdup 闭合用的固定气相源项 [mol/s]。

    优先 ``_vm_gas_source_cache``（Vorab 热解/干燥 single-shot）；无 Vorab 增量时
    回退 ``max(R_gas,0)``（床层下游 char 气化等）。避免 init 态 R5/R8 净消耗使
    CH4 等 ``max(R_gas,0)=0`` 而 Check2 rel 爆炸。

    Ref: Hamel (1999) Vorabrechnung; Eq.2-4/2-5 init holdup bridge
    """
    j = int(species_index)
    vorab_source = 0.0
    if bool(getattr(cell, "_vm_cache_valid", False)):
        gas_src = getattr(cell, "_vm_gas_source_cache", None)
        if gas_src is not None:
            vorab_source = float(max(float(np.asarray(gas_src, dtype=np.float64)[j]), 0.0))
    if vorab_source > _VORAB_MAJOR_PRODUCT_SOURCE_FLOOR_MOL_S:
        return vorab_source
    return float(max(float(cell.R_gas_d[j]) + float(cell.R_gas_b[j]), 0.0))


def _resolve_product_holdup_bed_cell_limit(cfg: Any, bed_cell_limit: int) -> int:
    """VM 分区启用时，产物 seed/Check2 默认仅覆盖分区 bed 格。

    ``bed_cell_limit<=0`` 在无 VM 分区时表示全床层；有 VM 分区时表示分区宽度。
    """
    limit = int(bed_cell_limit)
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    if n_zone > 0:
        if limit <= 0:
            return n_zone
        return int(max(limit, n_zone))
    return limit


def seed_major_product_holdup_from_vorab_gas_balance(
    cells: List[Cell],
    cfg: Any,
    *,
    bed_cell_limit: int = 2,
    max_passes: int = 3,
    blend_factor: float = 0.35,
    apply_bc_fn: Any | None = None,
) -> dict[str, float]:
    """将 bed 区产物 holdup 温和对齐固定 Vorab 源项 + 轴向 inflow（Eq.2.4/2.5）。

    在 ``_initialize_thesis_single_shot_vorab_sources`` 之后调用：仅调整
    CO/CO2/H2/CH4 总量 toward ``N_in + R_gas``，主气化剂 holdup 不动；
    两相按本地水力份额分配。采用阻尼混合，避免一次性吃掉完整残差导致 holdup 爆炸。

    Ref: Hamel (1999) Eq.2-4/2-5; Vorabrechnung single-shot source bridge
    """
    from src.core.connectivity import propagate_upstream

    idx = GAS_SPECIES_INDEX
    effective_limit = _resolve_product_holdup_bed_cell_limit(cfg, bed_cell_limit)
    bed_indices = _resolve_product_bed_cell_indices(cells, effective_limit)
    if not bed_indices:
        return {"cells_seeded": 0.0, "closure_passes": 0.0, "closure_total_mol_s": 0.0}

    n_vm_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    all_bed = _bed_cell_indices(cells)
    vm_zone_bed_indices = set(all_bed[: max(n_vm_zone, 0)])

    alpha = float(np.clip(blend_factor, 0.05, 1.0))
    n_passes = int(max(max_passes, 1))
    cells_touched = 0
    closure_total = 0.0

    for _ in range(n_passes):
        touched_this_pass = 0
        for ki, i in enumerate(bed_indices):
            cell = cells[i]
            if str(getattr(cell, "cell_type", "bed")) != "bed":
                continue
            if apply_bc_fn is not None:
                apply_bc_fn()
            cell.residuals()
            cell.calc_hydrodynamics()
            bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
            dense_share = float(1.0 - bubble_share)
            updated = False
            for sp in MAJOR_PRODUCT_SEED_SPECIES:
                j = int(idx[sp])
                inflow_total = float(
                    max(
                        cell.N_zu_d[j]
                        + cell.N_zu_b[j]
                        + cell.N_d_in[j]
                        + cell.N_b_in[j]
                        + cell.N_rez_d[j]
                        + cell.N_rez_b[j],
                        0.0,
                    )
                )
                in_vm_zone = bool(n_vm_zone > 0 and i in vm_zone_bed_indices)
                if in_vm_zone:
                    if not bool(getattr(cell, "_vm_cache_valid", False)):
                        continue
                    gas_src = getattr(cell, "_vm_gas_source_cache", None)
                    vorab_only = float(
                        max(float(np.asarray(gas_src, dtype=np.float64)[j]), 0.0)
                    ) if gas_src is not None else 0.0
                    if vorab_only <= _VORAB_MAJOR_PRODUCT_SOURCE_FLOOR_MOL_S:
                        continue
                    source_total = vorab_only
                else:
                    source_total = _major_product_init_gas_source_mol_s(cell, j)
                target_total = max(inflow_total + source_total, 1e-12)
                current_total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                rel_gap = abs(target_total - current_total) / max(target_total, 1e-12)
                if rel_gap <= 0.02:
                    continue
                blended_total = current_total + alpha * (target_total - current_total)
                delta = abs(blended_total - current_total)
                if delta <= 1.0e-18:
                    continue
                cell.N_d[j] = max(blended_total * dense_share, 1e-12)
                cell.N_b[j] = max(blended_total * bubble_share, 1e-12)
                closure_total += delta
                updated = True
            if updated:
                cell._thermo_cache_valid = False
                cell._hydro_cache_valid = False
                touched_this_pass += 1
            if ki == 0 and len(bed_indices) > 1:
                propagate_upstream(cells, cfg, bed_indices[1])
        for j in range(bed_indices[0] + 1, len(cells)):
            propagate_upstream(cells, cfg, j)
        cells_touched = max(cells_touched, touched_this_pass)

    return {
        "cells_seeded": float(cells_touched),
        "closure_passes": float(n_passes),
        "closure_total_mol_s": float(closure_total),
        "bed_cells_targeted": float(len(bed_indices)),
    }


def seed_pyrolysis_gas_dense_holdup_for_nr_x0(
    cells: List[Cell],
    *,
    holdup_scale: float = 1.0,
    residual_closure_passes: int = 2,
) -> dict[str, float]:
    """Seed NH3/TAR dense-phase holdup from fixed Vorabrechnung pyrolysis sources.

    Direct drying/pyrolysis gas enters the suspension phase only.  Bubble outlets
    for these species stay at zero so phase-split rows are not polluted by a
    hidden bubble source.

    The seed uses signed dense-balance closure moves (not ``tau * source``) so
    bottom cells with ``N_d=0`` but positive pyrolysis source do not overshoot.

    Ref: connectivity ``preproject_bed_gas_temperature_to_energy_closure`` dense
    support; Hamel Vorabrechnung start-value bridge.
    """
    idx = GAS_SPECIES_INDEX
    scale = float(max(holdup_scale, 0.0))
    n_passes = int(max(residual_closure_passes, 1))
    cells_touched = 0
    closure_total = 0.0

    for _ in range(n_passes):
        for cell in cells:
            if str(getattr(cell, "cell_type", "bed")) != "bed":
                continue
            if not bool(getattr(cell, "_vm_cache_valid", False)):
                continue
            gas_src = np.asarray(getattr(cell, "_vm_gas_source_cache", None), dtype=np.float64)
            pyro_source = sum(max(float(gas_src[int(idx[sp])]), 0.0) for sp in PYROLYSIS_DENSE_SEED_SPECIES)
            if pyro_source <= 1e-20:
                continue
            res = np.asarray(cell.residuals(), dtype=np.float64)
            tau_val = max(float(cell.geo.dh) / max(float(cell.u_mf), 1e-3), 1e-4)
            updated = False
            for sp in PYROLYSIS_DENSE_SEED_SPECIES:
                j = int(idx[sp])
                src_rate = max(float(gas_src[j]), 0.0)
                current = float(cell.N_d[j])
                dense_gap = scale * float(res[j])
                if src_rate > 1.0e-20:
                    tau_target = scale * src_rate * tau_val
                    if current <= 1.0e-12:
                        dense_gap = tau_target - current
                    elif dense_gap < 0.0 and current <= max(tau_target, 1.0e-12) * 1.05:
                        dense_gap = 0.0
                    elif dense_gap > 0.0:
                        dense_gap = min(dense_gap, max(tau_target - current, 0.0))
                elif abs(dense_gap) <= 1e-20:
                    continue
                if abs(dense_gap) <= 1e-20:
                    continue
                cell.N_d[j] = max(current + dense_gap, 1e-12)
                cell.N_b[j] = 0.0
                closure_total += abs(dense_gap)
                updated = True
            if updated:
                cell._thermo_cache_valid = False
                cell._hydro_cache_valid = False
                cells_touched += 1

    return {
        "cells_seeded": float(cells_touched),
        "tau_seed_total_mol_s": 0.0,
        "residual_closure_total_mol_s": float(closure_total),
    }


def _cell_gas_inflow_mol_s(cell: Cell, species: str) -> float:
    """单格某物种 gas inflow 总量 [mol/s]（zu + in + rez）。"""
    j = int(GAS_SPECIES_INDEX[species])
    return float(
        max(
            cell.N_zu_d[j]
            + cell.N_zu_b[j]
            + cell.N_d_in[j]
            + cell.N_b_in[j]
            + cell.N_rez_d[j]
            + cell.N_rez_b[j],
            0.0,
        )
    )


def _zero_vorab_budget() -> VorabrechnungCellBudget:
    return VorabrechnungCellBudget(
        cell_index=0,
        axial_xi=0.0,
        tau_cumulative_s=0.0,
        t_budget_K=300.0,
        x_dry_cumulative=0.0,
        x_vm_cumulative=0.0,
        o2_remaining_mol_s=0.0,
        nC_vm_mol_s=0.0,
        nH_vm_mol_s=0.0,
        nO_vm_mol_s=0.0,
    )


def _apply_gas_mol_flows_to_cell_holdup(cell: Cell, flows: dict[str, float]) -> None:
    """按本地水力份额将总摩尔流写入 N_d/N_b。"""
    idx = GAS_SPECIES_INDEX
    bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
    dense_share = float(1.0 - bubble_share)

    cell.N_d.fill(0.0)
    cell.N_b.fill(0.0)
    for sp, val in flows.items():
        if sp not in idx:
            continue
        total = float(max(val, 0.0))
        if total <= 1.0e-20:
            continue
        j = int(idx[sp])
        cell.N_d[j] = max(total * dense_share, 1e-12)
        cell.N_b[j] = max(total * bubble_share, 1e-12)


def _apply_fast_oxidation_stoich_to_phase_holdup(
    N: np.ndarray,
    *,
    fuel_type: str = "coal",
    n_passes: int = 3,
) -> dict[str, float]:
    """对单相 holdup 做快氧化化学计量投影（Startwert，非 F(x) 限幅）。

    每遍：R12（H₂）→ R5（CO）→ R6（CH₄）→ R10（焦油）。
    R10 产 H₂/CO，故多遍迭代直到 O₂ 或燃料侧耗尽。

    Ref: Hamel (1999) §2.1 Startwertwahl；Kap.5 R5/R6/R10/R12 计量
    """
    from src.kinetics.tar_reactions import get_tar_component_stoichiometry

    idx = GAS_SPECIES_INDEX
    n = np.asarray(N, dtype=np.float64)
    # 按 TAR1/TAR2 组分分别用 R10 计量，避免 lumped 配比与库存不一致留下残余
    sto10_by_comp: dict[str, dict[str, float]] = {}
    for comp in ("TAR1", "TAR2"):
        try:
            sto10_by_comp[comp] = get_tar_component_stoichiometry("R10", str(fuel_type), comp)
        except ValueError:
            continue

    xi12 = xi5 = xi6 = xi10 = 0.0
    for _ in range(max(int(n_passes), 1)):
        d12 = float(min(max(n[idx["H2"]], 0.0), 2.0 * max(n[idx["O2"]], 0.0)))
        n[idx["H2"]] = max(float(n[idx["H2"]]) - d12, 0.0)
        n[idx["O2"]] = max(float(n[idx["O2"]]) - 0.5 * d12, 0.0)
        n[idx["H2O"]] = float(n[idx["H2O"]]) + d12
        xi12 += d12

        d5 = float(min(max(n[idx["CO"]], 0.0), 2.0 * max(n[idx["O2"]], 0.0)))
        n[idx["CO"]] = max(float(n[idx["CO"]]) - d5, 0.0)
        n[idx["O2"]] = max(float(n[idx["O2"]]) - 0.5 * d5, 0.0)
        n[idx["CO2"]] = float(n[idx["CO2"]]) + d5
        xi5 += d5

        d6 = float(min(max(n[idx["CH4"]], 0.0), max(n[idx["O2"]], 0.0) / 1.5))
        n[idx["CH4"]] = max(float(n[idx["CH4"]]) - d6, 0.0)
        n[idx["O2"]] = max(float(n[idx["O2"]]) - 1.5 * d6, 0.0)
        n[idx["CO"]] = float(n[idx["CO"]]) + d6
        n[idx["H2O"]] = float(n[idx["H2O"]]) + 2.0 * d6
        xi6 += d6

        for comp, sto in sto10_by_comp.items():
            nu_o2 = abs(float(sto.get("O2", 0.0)))
            if nu_o2 <= 1e-15:
                continue
            d10 = float(min(max(float(n[idx[comp]]), 0.0), max(float(n[idx["O2"]]), 0.0) / nu_o2))
            if d10 <= 0.0:
                continue
            n[idx[comp]] = max(float(n[idx[comp]]) - d10, 0.0)
            n[idx["O2"]] = max(float(n[idx["O2"]]) - nu_o2 * d10, 0.0)
            n[idx["CO"]] = float(n[idx["CO"]]) + float(sto.get("CO", 0.0)) * d10
            n[idx["H2"]] = float(n[idx["H2"]]) + float(sto.get("H2", 0.0)) * d10
            xi10 += d10

    # 末遍 R10 产 H₂/CO：再做一轮轻质氧化闭合
    d12 = float(min(max(n[idx["H2"]], 0.0), 2.0 * max(n[idx["O2"]], 0.0)))
    n[idx["H2"]] = max(float(n[idx["H2"]]) - d12, 0.0)
    n[idx["O2"]] = max(float(n[idx["O2"]]) - 0.5 * d12, 0.0)
    n[idx["H2O"]] = float(n[idx["H2O"]]) + d12
    xi12 += d12
    d5 = float(min(max(n[idx["CO"]], 0.0), 2.0 * max(n[idx["O2"]], 0.0)))
    n[idx["CO"]] = max(float(n[idx["CO"]]) - d5, 0.0)
    n[idx["O2"]] = max(float(n[idx["O2"]]) - 0.5 * d5, 0.0)
    n[idx["CO2"]] = float(n[idx["CO2"]]) + d5
    xi5 += d5
    d6 = float(min(max(n[idx["CH4"]], 0.0), max(n[idx["O2"]], 0.0) / 1.5))
    n[idx["CH4"]] = max(float(n[idx["CH4"]]) - d6, 0.0)
    n[idx["O2"]] = max(float(n[idx["O2"]]) - 1.5 * d6, 0.0)
    n[idx["CO"]] = float(n[idx["CO"]]) + d6
    n[idx["H2O"]] = float(n[idx["H2O"]]) + 2.0 * d6
    xi6 += d6
    # R6 可能再产 CO；吃掉残余
    d5b = float(min(max(n[idx["CO"]], 0.0), 2.0 * max(n[idx["O2"]], 0.0)))
    n[idx["CO"]] = max(float(n[idx["CO"]]) - d5b, 0.0)
    n[idx["O2"]] = max(float(n[idx["O2"]]) - 0.5 * d5b, 0.0)
    n[idx["CO2"]] = float(n[idx["CO2"]]) + d5b
    xi5 += d5b

    N[:] = n
    return {"xi12": xi12, "xi5": xi5, "xi6": xi6, "xi10": xi10}


def _r1_alpha_from_cell(cell: Cell) -> float:
    """Hamel Eq.5.9–5.10：``α = 1/φ_c``，按炭质量加权粒径。"""
    T = float(max(getattr(cell, "T", 1150.0), 300.0))
    d_p_classes = np.asarray(getattr(cell.solid, "d_p_classes", [0.001]), dtype=np.float64)
    m_char = np.maximum(np.asarray(cell.m_solid[:, S_CHAR], dtype=np.float64), 0.0)
    if d_p_classes.size <= 0:
        d_c = 1.0e-3  # [m]
    elif m_char.size == d_p_classes.size and float(np.sum(m_char)) > 1e-15:
        d_c = float(np.average(d_p_classes, weights=m_char))
    else:
        d_c = float(max(d_p_classes[0], 1e-9))
    phi = float(phi_c(d_c, T))
    return 1.0 / max(phi, 1e-30)


def _apply_char_r1_stoich_to_holdup(
    N: np.ndarray,
    m_char: np.ndarray,
    *,
    alpha: float,
) -> dict[str, float]:
    """气相残氧 + 炭库存的 R1 Startwert 计量（Hamel Eq.5.8 链）。

    未知量混用：``N`` 为摩尔流 [mol/s]，``m_char`` 为库存 [kg]。投影把 1 s 的
    氧流接到炭库存上，使底格不再处于「残氧 ∩ 空合成气」的 R12 刚性格点。

    计量（``ξ`` = 炭摩尔 [mol]，``α = 1/φ_c``）::

        O₂  -= α ξ
        CO  += 2 (1-α) ξ
        CO₂ += (2α-1) ξ
        C   -= ξ

    Ref: Hamel (1999) Eq. (5.8)–(5.10)；``src/core/cell_kinetics.py`` R1 源项
    """
    idx = GAS_SPECIES_INDEX
    n = np.asarray(N, dtype=np.float64)
    m = np.maximum(np.asarray(m_char, dtype=np.float64), 0.0)
    alpha = float(np.clip(alpha, 0.5, 1.0))
    o2 = float(max(n[idx["O2"]], 0.0))
    c_mol = float(np.sum(m)) / max(float(R4_M_C), 1e-30)  # [mol] 1 s 当量
    if o2 <= 1e-15 or c_mol <= 1e-15:
        return {"xi_c": 0.0, "o2_used": 0.0, "alpha": alpha}
    xi = float(min(c_mol, o2 / alpha))
    if xi <= 1e-15:
        return {"xi_c": 0.0, "o2_used": 0.0, "alpha": alpha}
    n[idx["O2"]] = max(o2 - alpha * xi, 0.0)
    n[idx["CO"]] = float(n[idx["CO"]]) + 2.0 * (1.0 - alpha) * xi
    n[idx["CO2"]] = float(n[idx["CO2"]]) + (2.0 * alpha - 1.0) * xi
    take_kg = xi * float(R4_M_C)
    total_kg = float(np.sum(m))
    if total_kg > 1e-15:
        m *= max(total_kg - take_kg, 0.0) / total_kg
    N[:] = n
    m_char[:] = m
    return {"xi_c": xi, "o2_used": alpha * xi, "alpha": alpha}


# Species whose totals change under fast-oxidation stoichiometry (R12/R5/R6/R10).
_FAST_OX_PARTICIPATING_SPECIES: tuple[str, ...] = (
    "H2",
    "CO",
    "CH4",
    "TAR1",
    "TAR2",
    "O2",
    "H2O",
    "CO2",
)


def _assign_participating_totals_preserving_phase_fractions(
    cell: Cell,
    n_tot_after: np.ndarray,
    n_d_before: np.ndarray,
    n_b_before: np.ndarray,
    *,
    dense_share_fallback: float,
) -> None:
    """Write only participating species; preserve each species' prior dense fraction.

    Non-participating holdups (N₂, H₂S, …) stay untouched so P is locally continuous:
    tiny fuel–O₂ overlap ⇒ O(overlap) displacement, not a full-vector redistribute.
    """
    idx = GAS_SPECIES_INDEX
    dense_share_fallback = float(np.clip(dense_share_fallback, 0.0, 1.0))
    for sp in _FAST_OX_PARTICIPATING_SPECIES:
        j = idx[sp]
        total_new = float(max(n_tot_after[j], 0.0))
        total_old = float(max(n_d_before[j], 0.0) + max(n_b_before[j], 0.0))
        if total_old > 1e-15:
            frac_d = float(max(n_d_before[j], 0.0) / total_old)
        else:
            frac_d = dense_share_fallback
        frac_d = float(np.clip(frac_d, 0.0, 1.0))
        cell.N_d[j] = total_new * frac_d
        cell.N_b[j] = total_new * (1.0 - frac_d)


def project_bed_holdup_fast_oxidation_closure(
    cells: List[Cell],
    cfg: Any | None = None,
    *,
    bed_cell_limit: int | None = None,
    phase_update: str = "full_hydro",
) -> dict[str, float]:
    """床层 holdup 快氧化投影：消除可燃气体/焦油与 O₂ 共存；残氧再按 R1 吃炭。

    ``transport_from_vorab`` 把 ``o2_remaining`` 与热解产物直接叠加，未做氧化闭合；
    同构裸动力学下 R12/R10 会把 F(x) 撑到天文量级。气相燃料耗尽后若仍有 O₂，
    再按 Hamel R1 计量把残氧接到炭库存上（Kap.7：氧在底区很快耗尽），避免
    底格落在「残氧 ∩ 空合成气、源项仍产 H₂」的刚性格点。只改 holdup，不改 R(x)。

    Parameters
    ----------
    phase_update:
        ``full_hydro`` — x₀ Startwert：合计氧化后按水力份额重分整条气体向量。
        ``local_participating`` — 线搜索旧路径：合计氧化后保留各物种原相分数。
        ``per_phase`` — 各相独立做气相快氧化；气泡 O₂ 与悬浮相 H₂ 可以共存。
          不跑 R1 置零。Newton 试探若把 O₂ 搬进有 H₂ 的相，只在该相按计量烧掉。

    Ref: Hamel (1999) §2.1 Startwertwahl
    """
    if cfg is not None and not bool(
        getattr(cfg, "vorab_transport_x0_fast_oxidation_closure_thesis", True)
    ):
        return {
            "enabled": 0.0,
            "cells_touched": 0.0,
            "xi12_total": 0.0,
            "xi_r1_total": 0.0,
        }

    mode = str(phase_update).strip().lower()
    if mode not in {"full_hydro", "local_participating", "per_phase"}:
        raise ValueError(
            f"Unsupported phase_update={phase_update!r}; "
            "expected 'full_hydro', 'local_participating', or 'per_phase'"
        )

    fuel_type = str(getattr(cfg, "fuel_type", "coal") if cfg is not None else "coal")
    bed_indices = _bed_cell_indices(cells)
    if bed_cell_limit is not None:
        bed_indices = bed_indices[: max(int(bed_cell_limit), 0)]

    xi12_total = 0.0
    xi5_total = 0.0
    xi6_total = 0.0
    xi10_total = 0.0
    xi_r1_total = 0.0
    touched = 0
    idx = GAS_SPECIES_INDEX

    def _accumulate_gas_fastox(n_vec: np.ndarray) -> dict[str, float]:
        d = _apply_fast_oxidation_stoich_to_phase_holdup(n_vec, fuel_type=fuel_type)
        return d

    apply_r1 = bool(
        getattr(cfg, "vorab_transport_x0_char_oxidation_o2_closure_thesis", True)
    ) if cfg is not None else True
    for i in bed_indices:
        cell = cells[i]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        n_d_before = np.asarray(cell.N_d, dtype=np.float64).copy()
        n_b_before = np.asarray(cell.N_b, dtype=np.float64).copy()
        n_tot = n_d_before + n_b_before
        fuel_before = float(
            n_tot[idx["H2"]]
            + n_tot[idx["CO"]]
            + n_tot[idx["CH4"]]
            + n_tot[idx["TAR1"]]
            + n_tot[idx["TAR2"]]
        )
        o_before = float(n_tot[idx["O2"]])
        char_kg = float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0)))
        if o_before <= 1e-15:
            continue
        if fuel_before <= 1e-15 and not (apply_r1 and char_kg > 1e-15 and mode != "per_phase"):
            continue
        if mode == "per_phase":
            n_d = n_d_before.copy()
            n_b = n_b_before.copy()
            diag = {"xi12": 0.0, "xi5": 0.0, "xi6": 0.0, "xi10": 0.0}
            for n_vec in (n_d, n_b):
                fuel_p = float(
                    n_vec[idx["H2"]]
                    + n_vec[idx["CO"]]
                    + n_vec[idx["CH4"]]
                    + n_vec[idx["TAR1"]]
                    + n_vec[idx["TAR2"]]
                )
                if float(n_vec[idx["O2"]]) > 1e-15 and fuel_p > 1e-15:
                    d = _accumulate_gas_fastox(n_vec)
                    for key in diag:
                        diag[key] = float(diag[key]) + float(d[key])
            cell.N_d[:] = np.maximum(n_d, 0.0)
            cell.N_b[:] = np.maximum(n_b, 0.0)
            xi12_total += float(diag["xi12"])
            xi5_total += float(diag["xi5"])
            xi6_total += float(diag["xi6"])
            xi10_total += float(diag["xi10"])
            cell._thermo_cache_valid = False
            touched += 1
            continue
        # 合计相上做氧化，避免分相独立投影制造虚假 phase-split
        diag = {"xi12": 0.0, "xi5": 0.0, "xi6": 0.0, "xi10": 0.0}
        if fuel_before > 1e-15:
            diag = _accumulate_gas_fastox(n_tot)
        o_left = float(max(n_tot[idx["O2"]], 0.0))
        if o_left > 1e-15 and char_kg > 1e-15 and apply_r1:
            r1d = _apply_char_r1_stoich_to_holdup(
                n_tot,
                cell.m_solid[:, S_CHAR],
                alpha=_r1_alpha_from_cell(cell),
            )
            xi_r1_total += float(r1d["xi_c"])
            if float(max(n_tot[idx["O2"]], 0.0)) > 1e-15:
                d2 = _accumulate_gas_fastox(n_tot)
                for key in ("xi12", "xi5", "xi6", "xi10"):
                    diag[key] = float(diag[key]) + float(d2[key])
        bubble_share = float(_hydrodynamic_bubble_share_for_cell(cell))
        dense_share = float(1.0 - bubble_share)
        if mode == "local_participating":
            _assign_participating_totals_preserving_phase_fractions(
                cell,
                n_tot,
                n_d_before,
                n_b_before,
                dense_share_fallback=dense_share,
            )
        else:
            cell.N_d[:] = np.maximum(n_tot * dense_share, 0.0)
            cell.N_b[:] = np.maximum(n_tot * bubble_share, 0.0)
            if i == 0 and cfg is not None:
                from src.core.connectivity.gas_inlet import resolve_bottom_gas_inlet_dense_fraction

                bottom_dense_frac = float(
                    getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", np.nan)
                )
                if not np.isfinite(bottom_dense_frac):
                    bottom_dense_frac = float(resolve_bottom_gas_inlet_dense_fraction(cell, cfg))
                bottom_dense_frac = float(np.clip(bottom_dense_frac, 0.0, 1.0))
                for sp in ("O2", "H2O", "N2"):
                    j = idx[sp]
                    total_sp = float(max(n_tot[j], 0.0))
                    cell.N_d[j] = total_sp * bottom_dense_frac
                    cell.N_b[j] = total_sp * (1.0 - bottom_dense_frac)
        xi12_total += float(diag["xi12"])
        xi5_total += float(diag["xi5"])
        xi6_total += float(diag["xi6"])
        xi10_total += float(diag["xi10"])
        cell._thermo_cache_valid = False
        touched += 1

    return {
        "enabled": 1.0,
        "cells_touched": float(touched),
        "xi12_total": float(xi12_total),
        "xi5_total": float(xi5_total),
        "xi6_total": float(xi6_total),
        "xi10_total": float(xi10_total),
        "xi_r1_total": float(xi_r1_total),
        "phase_update_local": float(mode == "local_participating"),
        "phase_update_per_phase": float(mode == "per_phase"),
    }


_BED0_EQ24_SPECIES: tuple[str, ...] = ("O2", "H2O", "N2", "CO", "CO2", "H2", "CH4")
_BED0_SYNGAS_SPECIES: tuple[str, ...] = ("CO", "H2", "CH4")
_BED0_INERT_EXCHANGE_SPECIES: tuple[str, ...] = ("N2", "H2O", "CO2")


def _split_bed0_inert_exchange_species(cell: Cell, *, n_passes: int) -> None:
    """按线性化 Eq.2.2 重分 N₂/H₂O/CO₂；不碰 O₂ 与合成气。

    ``n_tot = zu + R`` 贴 Eq.2.4 总和；``n_d = (zu_d + R_d + α n_tot) / (1+α+β)``
    给 N_ex 驱动力。R / 总量取当前态（种子与按相 fastox 之后），禁止把
    投影后再算的 R 喂回 O₂/合成气。

    Ref: Hamel (1999) Eq.2.2, 2.4, 2.5; handoff §5.56.
    """
    idx = GAS_SPECIES_INDEX
    n_pass = int(max(n_passes, 1))
    for _ in range(n_pass):
        cell.residuals()
        n_d_tot = float(max(np.sum(np.maximum(cell.N_d, 0.0)), 1e-12))
        n_b_tot = float(max(np.sum(np.maximum(cell.N_b, 0.0)), 1e-12))
        k_ex = (
            float(cell.K_bd)
            * float(cell.V_b)
            * float(cell.P)
            / (Rg * max(float(cell.T), 1.0))
        )
        alpha = float(np.clip(k_ex / n_b_tot, 0.0, 50.0))
        beta = float(np.clip(k_ex / n_d_tot, 0.0, 50.0))
        denom = 1.0 + alpha + beta
        for sp in _BED0_INERT_EXCHANGE_SPECIES:
            j = int(idx[sp])
            zu_d = float(cell.N_zu_d[j] + cell.N_d_in[j] + cell.N_rez_d[j])
            zu_b = float(cell.N_zu_b[j] + cell.N_b_in[j] + cell.N_rez_b[j])
            zu = zu_d + zu_b
            r_d = float(cell.R_gas_d[j])
            r_b = float(cell.R_gas_b[j])
            r_sum = float(np.clip(r_d + r_b, -1.0e3, 1.0e3))
            cap = 8.0 * max(abs(zu), abs(r_sum), 1.0)
            n_tot = float(np.clip(zu + r_sum, 1e-12, cap))
            n_d = float((zu_d + float(np.clip(r_d, -cap, cap)) + alpha * n_tot) / max(denom, 1e-12))
            n_d = float(np.clip(n_d, 1e-12, max(n_tot - 1e-12, 1e-12)))
            cell.N_d[j] = n_d
            cell.N_b[j] = float(max(n_tot - n_d, 1e-12))
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
    cell.calc_hydrodynamics()
    cell.residuals()


def apply_bed0_eq24_two_phase_startwert(
    cells: List[Cell],
    cfg: Any | None = None,
    *,
    n_passes: int = 1,
) -> dict[str, float]:
    """底格 Startwert：按 Eq.2.4/2.5 **总和** 闭合 N，再按相分工分相。

    不改 R(x) / K_bd / Eq.2.7。气相 fastox 之后调用：

    - 总量 ``N = N_zu + N_in + N_rez + R``（N_ex 在两相残差里抵消，总和不含交换）
    - O₂ 留在气泡相、热解合成气留在悬浮相，避免同相 C_H2·C_O2 点燃裸 R12
    - N₂/H₂O/CO₂ 用线性化 Eq.2.2 分配，使 y_b ≠ y_d 时 N_ex 有驱动力
    - 只用 fastox 态上冻结的 R 做投影；禁止把投影后再算出的 R 喂回 N
    - 种子后须再跑 ``apply_bed0_inert_exchange_split_startwert``（handoff §5.56）

    Ref: Hamel (1999) Eq.2.2, 2.4, 2.5; Kap.5 异相仅悬浮相
    """
    if cfg is not None and not bool(
        getattr(cfg, "vorab_bed0_eq24_two_phase_startwert_thesis", False)
    ):
        return {"enabled": 0.0, "passes": 0.0, "o2_res_sum": 0.0, "n_ex_o2": 0.0}
    if not cells:
        return {"enabled": 0.0, "passes": 0.0, "o2_res_sum": 0.0, "n_ex_o2": 0.0}
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return {"enabled": 0.0, "passes": 0.0, "o2_res_sum": 0.0, "n_ex_o2": 0.0}

    idx = GAS_SPECIES_INDEX
    n_pass = int(max(n_passes, 1))
    cell.calc_hydrodynamics()
    cell.residuals()
    frozen_r_d = np.asarray(cell.R_gas_d, dtype=np.float64).copy()
    frozen_r_b = np.asarray(cell.R_gas_b, dtype=np.float64).copy()
    n_d_tot = float(max(np.sum(np.maximum(cell.N_d, 0.0)), 1e-12))
    n_b_tot = float(max(np.sum(np.maximum(cell.N_b, 0.0)), 1e-12))
    k_ex = (
        float(cell.K_bd)
        * float(cell.V_b)
        * float(cell.P)
        / (Rg * max(float(cell.T), 1.0))
    )
    alpha = float(np.clip(k_ex / n_b_tot, 0.0, 50.0))
    beta = float(np.clip(k_ex / n_d_tot, 0.0, 50.0))
    denom = 1.0 + alpha + beta
    for _ in range(n_pass):
        for sp in _BED0_EQ24_SPECIES:
            j = int(idx[sp])
            zu_d = float(cell.N_zu_d[j] + cell.N_d_in[j] + cell.N_rez_d[j])
            zu_b = float(cell.N_zu_b[j] + cell.N_b_in[j] + cell.N_rez_b[j])
            zu = zu_d + zu_b
            r_d = float(frozen_r_d[j])
            r_b = float(frozen_r_b[j])
            r_sum = float(np.clip(r_d + r_b, -1.0e3, 1.0e3))
            cap = 8.0 * max(abs(zu), abs(r_sum), 1.0)
            n_tot = float(np.clip(zu + r_sum, 1e-12, cap))
            if sp == "O2":
                n_tot = float(max(n_tot, 0.25 * max(zu, 0.0), 1e-12))
                n_d = 0.0
                n_b = float(max(n_tot, 0.0))
            elif sp in _BED0_SYNGAS_SPECIES:
                n_d = float(max(n_tot, 0.0))
                n_b = 0.0
            else:
                r_d_c = float(np.clip(r_d, -cap, cap))
                n_d = float((zu_d + r_d_c + alpha * n_tot) / max(denom, 1e-12))
                n_d = float(np.clip(n_d, 1e-12, max(n_tot - 1e-12, 1e-12)))
                n_b = float(max(n_tot - n_d, 1e-12))
            cell.N_d[j] = n_d
            cell.N_b[j] = n_b
        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
    cell.calc_hydrodynamics()
    cell.residuals()
    j_o2 = int(idx["O2"])
    last_nex_o2 = float(cell.N_ex[j_o2])
    last_o2_res = float(
        cell.N_zu_d[j_o2]
        + cell.N_zu_b[j_o2]
        + cell.N_d_in[j_o2]
        + cell.N_b_in[j_o2]
        + cell.N_rez_d[j_o2]
        + cell.N_rez_b[j_o2]
        + cell.R_gas_d[j_o2]
        + cell.R_gas_b[j_o2]
        - cell.N_d[j_o2]
        - cell.N_b[j_o2]
    )

    y_d = cell._mole_fractions("d")
    y_b = cell._mole_fractions("b")
    return {
        "enabled": 1.0,
        "passes": float(n_pass),
        "o2_res_sum": float(last_o2_res),
        "n_ex_o2": float(last_nex_o2),
        "y_o2_d": float(y_d[int(idx["O2"])]),
        "y_o2_b": float(y_b[int(idx["O2"])]),
        "n_h2": float(cell.N_d[int(idx["H2"])] + cell.N_b[int(idx["H2"])]),
        "n_o2": float(cell.N_d[int(idx["O2"])] + cell.N_b[int(idx["O2"])]),
    }


def apply_bed0_r1_oxidizer_seed(
    cells: List[Cell],
    cfg: Any | None = None,
) -> dict[str, float]:
    """底格 Startwert：从气泡挪 O₂ 进悬浮相，按相 fastox 后留下点亮 R1 的残氧。

    §5.49 分相让 R1=0。Newton 第一步的 ``dN_d,O2≈0``，放开上段 ``dO2+`` 又会
    把线搜索打到 rms~10²。这里只改 x₀：把 ``stoich(H₂/CO/CH₄)+seed`` 从
    ``N_b,O2`` 转到 ``N_d,O2``（总量守恒），再按相 fastox。剩余 seed 点亮 R1。

    ``vorab_bed0_r1_oxidizer_seed_mol_s_thesis<=0`` 时关闭。

    Source: Hamel (1999) Eq.2.4 Startwertwahl; Kap.5 异相仅悬浮相; handoff §5.52.
    """
    seed = float(getattr(cfg, "vorab_bed0_r1_oxidizer_seed_mol_s_thesis", 0.0) or 0.0)
    if cfg is None or seed <= 0.0 or not cells:
        return {"enabled": 0.0, "moved_o2": 0.0, "leftover_o2_d": 0.0, "seed_mol_s": 0.0}
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return {"enabled": 0.0, "moved_o2": 0.0, "leftover_o2_d": 0.0, "seed_mol_s": seed}

    idx = GAS_SPECIES_INDEX
    j_o2 = int(idx["O2"])
    stoich = float(
        0.5 * max(float(cell.N_d[idx["H2"]]), 0.0)
        + 0.5 * max(float(cell.N_d[idx["CO"]]), 0.0)
        + 1.5 * max(float(cell.N_d[idx["CH4"]]), 0.0)
    )
    available = float(max(cell.N_b[j_o2], 0.0))
    moved = float(min(max(stoich + seed, 0.0), available))
    if moved <= 0.0:
        return {
            "enabled": 1.0,
            "moved_o2": 0.0,
            "leftover_o2_d": float(max(cell.N_d[j_o2], 0.0)),
            "seed_mol_s": seed,
            "stoich_o2_d": stoich,
        }
    cell.N_b[j_o2] = float(max(float(cell.N_b[j_o2]) - moved, 0.0))
    cell.N_d[j_o2] = float(max(float(cell.N_d[j_o2]), 0.0) + moved)
    cell._thermo_cache_valid = False
    project_bed_holdup_fast_oxidation_closure(
        cells,
        cfg,
        bed_cell_limit=1,
        phase_update="per_phase",
    )
    leftover = float(max(cell.N_d[j_o2], 0.0))
    return {
        "enabled": 1.0,
        "moved_o2": moved,
        "leftover_o2_d": leftover,
        "seed_mol_s": seed,
        "stoich_o2_d": stoich,
    }


def apply_bed0_inert_exchange_split_startwert(
    cells: List[Cell],
    cfg: Any | None = None,
    *,
    n_passes: int = 3,
) -> dict[str, float]:
    """底格 Startwert：种子/fastox 之后只重分 N₂/H₂O/CO₂。

    §5.49 分相用冻结 R 与分相前总量。R1 种子 + 按相 fastox 改了 y 与
    H₂O/CO₂ 总和后，N_ex 仍按旧分相，N₂ 两相残差对称到 Ê≈±0.21。
    这里重跑线性化 Eq.2.2，不改 O₂/合成气、不改 k / K_bd。

    门控与 ``vorab_bed0_eq24_two_phase_startwert_thesis`` 相同，无新开关。

    禁止对 O₂ 套本公式（handoff §5.57）：按 ``zu+R`` 抬 ``N_b,O2`` 会把 N_ex
    打到 10¹；带 ``N_d`` 下限的线性化能收 O₂ 帽子，但 bed1 焓变差。Newton
    已能把 O₂ Ê 收到 ~0.02，不必再改 x₀。

    Ref: Hamel (1999) Eq.2.2, 2.4, 2.5; handoff §5.56–5.57.
    """
    if cfg is not None and not bool(
        getattr(cfg, "vorab_bed0_eq24_two_phase_startwert_thesis", False)
    ):
        return {"enabled": 0.0, "passes": 0.0, "n2_res_d": 0.0, "n_ex_n2": 0.0}
    if not cells:
        return {"enabled": 0.0, "passes": 0.0, "n2_res_d": 0.0, "n_ex_n2": 0.0}
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return {"enabled": 0.0, "passes": 0.0, "n2_res_d": 0.0, "n_ex_n2": 0.0}

    n_pass = int(max(n_passes, 1))
    leftover_o2 = float(max(cell.N_d[int(GAS_SPECIES_INDEX["O2"])], 0.0))
    _split_bed0_inert_exchange_species(cell, n_passes=n_pass)
    idx = GAS_SPECIES_INDEX
    j_n2 = int(idx["N2"])
    zu_d = float(cell.N_zu_d[j_n2] + cell.N_d_in[j_n2] + cell.N_rez_d[j_n2])
    n2_res_d = float(zu_d + cell.R_gas_d[j_n2] - cell.N_d[j_n2] + cell.N_ex[j_n2])
    return {
        "enabled": 1.0,
        "passes": float(n_pass),
        "n2_res_d": n2_res_d,
        "n_ex_n2": float(cell.N_ex[j_n2]),
        "leftover_o2_d": leftover_o2,
        "leftover_o2_d_after": float(max(cell.N_d[int(idx["O2"])], 0.0)),
    }


def _build_transport_holdup_gas_flows(
    cell: Cell,
    *,
    cell_index: int,
    o2_remaining: float,
    H2O_feed: float,
    N2_feed: float,
    pyro_increment: dict[str, float],
) -> dict[str, float]:
    """Hamel §2.1 transport x₀：主气 bed0 来自进料，上游继承 inflow；产物 inflow+热解增量。"""
    product_species = ("CO", "CO2", "H2", "CH4")
    flows: dict[str, float] = {}

    if cell_index == 0:
        flows["O2"] = float(max(o2_remaining, 0.0))
        flows["H2O"] = float(max(H2O_feed, 0.0))
        flows["N2"] = float(max(N2_feed, 0.0))
    else:
        for sp in ("O2", "H2O", "N2"):
            flows[sp] = _cell_gas_inflow_mol_s(cell, sp)

    for sp in product_species:
        inflow = _cell_gas_inflow_mol_s(cell, sp)
        inc = float(max(pyro_increment.get(sp, 0.0), 0.0))
        flows[sp] = max(inflow + inc, inflow, 1.0e-12 if inc > 0.0 else 0.0)

    for sp in ("H2O",):
        if cell_index > 0:
            inflow = _cell_gas_inflow_mol_s(cell, sp)
            inc = float(max(pyro_increment.get(sp, 0.0), 0.0))
            flows[sp] = max(inflow + inc, flows.get(sp, inflow))

    return flows


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
    cell_budgets: list[VorabrechnungCellBudget] | None = None,
    major_gibbs_temperature_floor_K: float = 0.0,
    x0_holdup_mode: str = "",
    cfg: Any | None = None,
    use_pyrolysis_gibbs_composition: bool = True,
    pyrolysis_gibbs_solver_mode: str = "hamel_reduced",
    pyrolysis_gibbs_temperature_floor_K: float = 750.0,
) -> None:
    """按轴向高度与估计热解进度，写入各 cell 的 N_d/N_b/T/m_solid 初值（就地修改）。

    当 ``cell_budgets`` 非空时采用 A 档 Vorabrechnung 预算（``build_vorabrechnung_cell_budgets``）；
    否则保留 legacy 启发式（``exp(-5·xi)``、``tau=dh/u_mf``）。

    x₀ holdup 模式（``x0_holdup_mode`` / ``normalize_x0_holdup_mode``）：
    - ``transport_from_vorab``（Hamel §2.1 默认）：主气沿 inflow/进料，产物 inflow+§4.3 热解增量；
      **禁止**逐格 major-Gibbs holdup（避免 ~1105 K CO/CO₂ 分岔）。
    - ``major_gibbs_per_cell``：legacy 扩展，每格 A1 Gibbs（``use_hamel_major_gibbs_x0=True`` 等价）。
    - ``legacy_heuristic``：固定比例 ``_legacy_major_seed``。

    Source: Hamel (1999) §2.1 Startwertwahl, §4.3 pyrolysis gas composition
    """
    n = len(cells)
    idx = GAS_SPECIES_INDEX

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
    n_daf_per_c = (float(getattr(cfg, "nitrogen_fraction", 0.0)) * 12.011e-3) / max(
        float(C_dry) * 14.0067e-3, 1e-12
    )
    s_daf_per_c = (float(getattr(cfg, "S_dry", 0.0)) * 12.011e-3) / max(
        float(C_dry) * 32.065e-3, 1e-12
    )
    sulfur_volatile_frac = float(getattr(cells[0].solid, "sulfur_volatile_frac", 0.5)) if cells else 0.5
    from src.core.species import fuel_tar_target_hc_ratio
    from src.thermal.hegermann_yield import xc_koks_hegermann, xc_teer_hegermann

    pyrolysis_tar_frac_fixed = (
        float(getattr(cells[0].solid, "pyrolysis_tar_carbon_frac", 0.2)) if cells else 0.2
    )
    pyrolysis_tar_yield_mode = (
        str(getattr(cells[0].solid, "pyrolysis_tar_yield_mode", "hegermann")) if cells else "hegermann"
    )
    target_tar_hc = None
    if C_dry > 1e-12 and H_dry >= 0.0:
        target_tar_hc = fuel_tar_target_hc_ratio(C_dry, H_dry)

    w_vm_dry = vm_daf_frac * (1.0 - ash_frac)
    w_char_dry = 1.0 - ash_frac - w_vm_dry

    daem_fuel = "brown_coal"
    if fuel_type in ("wood", "biomass"):
        daem_fuel = "wood"
    elif fuel_type in ("coal", "brown_coal"):
        daem_fuel = "brown_coal"

    major_lambda_ws: np.ndarray | None = None
    major_lnN_ws: float | None = None
    x0_mode = normalize_x0_holdup_mode(x0_holdup_mode, use_hamel_major_gibbs_x0=use_hamel_major_gibbs_x0)
    use_transport_x0 = x0_mode == "transport_from_vorab"
    use_per_cell_gibbs = x0_mode == "major_gibbs_per_cell"
    budget_prev = _zero_vorab_budget()

    for i, cell in enumerate(cells):
        budget = cell_budgets[i] if cell_budgets is not None and i < len(cell_budgets) else None
        T_i = float(T_profile[i] if i < len(T_profile) else (budget.t_budget_K if budget is not None else 900.0))
        frac_height = float(budget.axial_xi if budget is not None else (i + 0.5) / max(n, 1))
        solid_T_init = max(float(getattr(cell, "T_in_solid", T_profile[0])), 293.15)

        if budget is not None:
            o2_remaining = float(budget.o2_remaining_mol_s)
            tau_est = float(budget.tau_cumulative_s)
            X_vm_cum = float(budget.x_vm_cumulative)
            nC_vm = float(budget.nC_vm_mol_s)
            nH_vm = float(budget.nH_vm_mol_s)
            nO_vm = float(budget.nO_vm_mol_s)
            x_dry_cum_prefetch = float(budget.x_dry_cumulative)
        else:
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
            x_dry_cum_prefetch = None

        curr_budget = budget if budget is not None else VorabrechnungCellBudget(
            cell_index=i,
            axial_xi=float(frac_height),
            tau_cumulative_s=float(tau_est),
            t_budget_K=float(T_i),
            x_dry_cumulative=float(x_dry_cum_prefetch or 0.0),
            x_vm_cumulative=float(X_vm_cum),
            o2_remaining_mol_s=float(o2_remaining),
            nC_vm_mol_s=float(nC_vm),
            nH_vm_mol_s=float(nH_vm),
            nO_vm_mol_s=float(nO_vm),
        )

        if use_transport_x0:
            if i > 0 and cfg is not None:
                from src.core.connectivity import propagate_upstream

                propagate_upstream(cells, cfg, i)
            cell.calc_hydrodynamics()
            if not bool(use_pyrolysis_gibbs_composition):
                raise ValueError("Hamel strict: generate_initial_x0 requires vorab_pyrolysis_gibbs_composition_thesis=True")
            if (
                pyrolysis_tar_yield_mode == "hegermann"
                and str(fuel_type).strip().lower() == "coal"
            ):
                pyrolysis_tar_frac = float(xc_teer_hegermann(float(cell.P), float(T_i)))
                pyrolysis_koks_frac = float(xc_koks_hegermann(float(cell.P), float(T_i)))
            else:
                pyrolysis_tar_frac = pyrolysis_tar_frac_fixed
                pyrolysis_koks_frac = 0.0
            pyro_inc = _incremental_pyrolysis_major_mol_s(
                curr_budget,
                budget_prev,
                T_K=float(T_i),
                P_Pa=float(cell.P),
                pyrolysis_gibbs_solver_mode=str(pyrolysis_gibbs_solver_mode),
                pyrolysis_gibbs_temperature_floor_K=float(pyrolysis_gibbs_temperature_floor_K),
                fuel_type=str(fuel_type),
                pyrolysis_tar_carbon_frac=pyrolysis_tar_frac,
                pyrolysis_koks_carbon_frac=pyrolysis_koks_frac,
                n_daf_per_c=n_daf_per_c,
                s_daf_per_c=s_daf_per_c,
                sulfur_volatile_frac=sulfur_volatile_frac,
                target_tar_hc_ratio=target_tar_hc,
                pyrolysis_ch4_to_h2_co_split_frac=float(
                    getattr(cfg, "pyrolysis_ch4_to_h2_co_split_frac", 0.0)
                )
                if cfg is not None
                else 0.0,
                pyrolysis_co_reduce_h2_inject_frac=float(
                    getattr(cfg, "pyrolysis_co_reduce_h2_inject_frac", 0.0)
                )
                if cfg is not None
                else 0.0,
            )
            flows = _build_transport_holdup_gas_flows(
                cell,
                cell_index=i,
                o2_remaining=float(o2_remaining),
                H2O_feed=float(H2O_feed),
                N2_feed=float(N2_feed),
                pyro_increment=pyro_inc,
            )
            _apply_gas_mol_flows_to_cell_holdup(cell, flows)
            if i == 0:
                bottom_dense_frac = float(getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", np.nan))
                if np.isfinite(bottom_dense_frac):
                    bottom_dense_frac = float(np.clip(bottom_dense_frac, 0.0, 1.0))
                    for sp in ("O2", "H2O", "N2"):
                        j = idx[sp]
                        total_sp = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                        cell.N_d[j] = max(total_sp * bottom_dense_frac, 1e-12)
                        cell.N_b[j] = max(total_sp * (1.0 - bottom_dense_frac), 1e-12)
            budget_prev = curr_budget
        elif use_per_cell_gibbs:
            major_guess: dict[str, float]
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
                    T=float(max(T_i, major_gibbs_temperature_floor_K)),
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
            _apply_gas_mol_flows_to_cell_holdup(cell, major_guess)
            if i == 0:
                bottom_dense_frac = float(getattr(cell, "_vorab_bottom_gas_inlet_dense_frac", np.nan))
                if np.isfinite(bottom_dense_frac):
                    bottom_dense_frac = float(np.clip(bottom_dense_frac, 0.0, 1.0))
                    for sp in ("O2", "H2O", "N2"):
                        j = idx[sp]
                        total_sp = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                        cell.N_d[j] = max(total_sp * bottom_dense_frac, 1e-12)
                        cell.N_b[j] = max(total_sp * (1.0 - bottom_dense_frac), 1e-12)
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
            _apply_gas_mol_flows_to_cell_holdup(cell, major_guess)
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
            if x_dry_cum_prefetch is not None:
                x_dry_cum = x_dry_cum_prefetch
            else:
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
                    feed_frac = np.asarray(cell.solid.mass_fractions, dtype=np.float64)
                    feed_sum = float(np.sum(feed_frac))
                    if feed_sum > 1e-12:
                        return feed_frac / feed_sum
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
