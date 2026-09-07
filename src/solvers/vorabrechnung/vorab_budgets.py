"""Vorabrechnung — vorab budgets segment (OPT-004)."""
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
def normalize_x0_holdup_mode(mode: str | None, *, use_hamel_major_gibbs_x0: bool) -> str:
    """解析 x₀ holdup 模式：transport_from_vorab | major_gibbs_per_cell | legacy_heuristic。"""
    m = str(mode or "").strip().lower()
    aliases = {
        "transport": "transport_from_vorab",
        "transport_from_vorab": "transport_from_vorab",
        "hamel_transport": "transport_from_vorab",
        "major_gibbs": "major_gibbs_per_cell",
        "major_gibbs_per_cell": "major_gibbs_per_cell",
        "gibbs_per_cell": "major_gibbs_per_cell",
        "legacy": "legacy_heuristic",
        "legacy_heuristic": "legacy_heuristic",
    }
    resolved = aliases.get(m)
    if resolved is not None:
        return resolved
    if use_hamel_major_gibbs_x0:
        return "major_gibbs_per_cell"
    return "legacy_heuristic"


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
    T_bottom_floor_K: float = 1050.0,
    T_bottom_cap_K: float = 1325.0,
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
    T_bottom = float(
        np.clip(
            T_inlet + dT_adiabatic,
            float(max(T_bottom_floor_K, 300.0)),
            float(max(T_bottom_cap_K, T_bottom_floor_K + 1.0)),
        )
    )
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


def march_vorabrechnung_bed_temperature_from_budgets(
    *,
    T_inlet: float,
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    moisture_wt: float,
    daf_feed_kg_s: float,
    vm_daf_frac: float,
    budgets: list["VorabrechnungCellBudget"],
    bed_indices: list[int],
    T_cap_K: float = 1200.0,
    heat_loss_frac: float = 0.05,
    dH_comb_J_per_mol_O2: float = 393_500.0,
    drying_enthalpy_J_per_kg_moisture: float = 2.26e6,
    pyrolysis_enthalpy_J_per_kg_daf_vm: float = 1.0e6,
    Cp_gas_J_per_mol_K: float = 32.0,
) -> np.ndarray:
    """沿床层由 Vorab 预算增量做焓步进，生成 bed 温度剖面 [K]。

    自床底 ``T_inlet`` 起，按相邻格间 O₂ 消耗（放热）、干燥与热解（吸热）逐步更新；
    仅作 Startwertwahl 初值，不读取 Fig 7.4 验证曲线。

    Source: Hamel (1999) Vorabrechnung 沿程预算 → 焓平衡初温；thesis 工程实现
    """
    if not bed_indices or not budgets:
        return np.array([], dtype=np.float64)

    moisture_frac = float(moisture_wt) / 100.0
    T_floor = float(max(T_inlet, 300.0))
    T_cap = float(max(T_cap_K, T_floor + 1.0))
    alpha_exo = float(np.clip(1.0 - heat_loss_frac, 0.0, 1.0))
    base_gas_mol_s = float(max(O2_feed + H2O_feed + N2_feed, 1e-6))

    marched = np.zeros(len(bed_indices), dtype=np.float64)
    prev_o2_consumed = 0.0
    prev_x_dry = 0.0
    prev_x_vm = 0.0
    T_prev = float(T_inlet)

    for j, bi in enumerate(bed_indices):
        budget = budgets[bi]
        o2_consumed = float(max(O2_feed - budget.o2_remaining_mol_s, 0.0))
        o2_inc = float(max(o2_consumed - prev_o2_consumed, 0.0))
        dry_inc = float(max(budget.x_dry_cumulative - prev_x_dry, 0.0))
        vm_inc = float(max(budget.x_vm_cumulative - prev_x_vm, 0.0))

        Q_exo = alpha_exo * o2_inc * float(dH_comb_J_per_mol_O2)
        Q_dry = dry_inc * moisture_frac * float(fuel_feed_kg_s) * float(drying_enthalpy_J_per_kg_moisture)
        Q_pyro = vm_inc * float(daf_feed_kg_s) * float(vm_daf_frac) * float(pyrolysis_enthalpy_J_per_kg_daf_vm)
        Q_net = Q_exo - Q_dry - Q_pyro

        pyro_gas_mol_s = float(max(budget.nC_vm_mol_s + budget.nH_vm_mol_s * 0.5, 0.0))
        n_gas = float(max(base_gas_mol_s + pyro_gas_mol_s, 1e-6))
        dT = Q_net / max(n_gas * float(Cp_gas_J_per_mol_K), 1.0)

        if j == 0:
            T_cell = float(np.clip(T_inlet + dT, T_floor, T_cap))
        else:
            T_cell = float(np.clip(T_prev + dT, T_floor, T_cap))

        marched[j] = T_cell
        T_prev = T_cell
        prev_o2_consumed = o2_consumed
        prev_x_dry = float(budget.x_dry_cumulative)
        prev_x_vm = float(budget.x_vm_cumulative)

    return marched


def _vorab_march_enthalpy_delta_T(
    *,
    o2_inc_mol_s: float,
    dry_inc: float,
    vm_inc: float,
    fuel_feed_kg_s: float,
    moisture_wt: float,
    daf_feed_kg_s: float,
    vm_daf_frac: float,
    pyro_gas_mol_s: float,
    base_gas_mol_s: float,
    heat_loss_frac: float,
    dH_comb_J_per_mol_O2: float,
    drying_enthalpy_J_per_kg_moisture: float,
    pyrolysis_enthalpy_J_per_kg_daf_vm: float,
    Cp_gas_J_per_mol_K: float,
) -> float:
    """单格 Vorab 预算增量 → 气相焓步进温升 [K]。"""
    moisture_frac = float(moisture_wt) / 100.0
    alpha_exo = float(np.clip(1.0 - heat_loss_frac, 0.0, 1.0))
    Q_exo = alpha_exo * float(o2_inc_mol_s) * float(dH_comb_J_per_mol_O2)
    Q_dry = float(dry_inc) * moisture_frac * float(fuel_feed_kg_s) * float(drying_enthalpy_J_per_kg_moisture)
    Q_pyro = (
        float(vm_inc) * float(daf_feed_kg_s) * float(vm_daf_frac) * float(pyrolysis_enthalpy_J_per_kg_daf_vm)
    )
    Q_net = Q_exo - Q_dry - Q_pyro
    n_gas = float(max(base_gas_mol_s + pyro_gas_mol_s, 1e-6))
    return float(Q_net / max(n_gas * float(Cp_gas_J_per_mol_K), 1.0))


def _compute_vorabrechnung_cell_budget(
    *,
    cell: Cell,
    cell_index: int,
    axial_xi: float,
    tau_i: float,
    T_i: float,
    solid_T_init: float,
    moisture_wt: float,
    vm_daf_frac: float,
    daem_fuel: str,
    daf_feed: float,
    c_daf: float,
    h_daf: float,
    o_daf: float,
    O2_feed: float,
    stoich: float,
    M_C: float,
    M_H: float,
    M_O: float,
) -> VorabrechnungCellBudget:
    X_vm, _ = devolatilization_rate_for_cell(
        T_bed=float(T_i),
        tau_cell=float(tau_i),
        T_init=float(solid_T_init),
        VM_daf=float(vm_daf_frac),
        fuel_type=daem_fuel,
    )
    x_vm_cum = float(np.clip(X_vm, 0.0, 1.0))

    dry_diag = solve_drying_CN(
        cell.solid.d_p,
        float(T_i),
        float(solid_T_init),
        moisture_wt,
        float(tau_i),
        Nr=12,
        Nt=80,
        pressure_pa=float(cell.P),
        return_history=True,
    )
    x_dry_cum = float(np.clip(float(dry_diag["X_dry"][-1]), 0.0, 1.0))

    m_vm_released = float(daf_feed) * float(vm_daf_frac) * x_vm_cum
    nC_vm = m_vm_released * c_daf / M_C
    nH_vm = m_vm_released * h_daf / M_H / 2.0
    nO_vm = m_vm_released * o_daf / M_O / 2.0

    o2_consumed = min(float(O2_feed), stoich * float(nC_vm))
    o2_remaining = max(float(O2_feed) - o2_consumed, 0.0)

    return VorabrechnungCellBudget(
        cell_index=int(cell_index),
        axial_xi=float(axial_xi),
        tau_cumulative_s=float(tau_i),
        t_budget_K=float(T_i),
        x_dry_cumulative=x_dry_cum,
        x_vm_cumulative=x_vm_cum,
        o2_remaining_mol_s=o2_remaining,
        nC_vm_mol_s=float(nC_vm),
        nH_vm_mol_s=float(nH_vm),
        nO_vm_mol_s=float(nO_vm),
    )


@dataclass(frozen=True)
class VorabrechnungCellBudget:
    """单格 Vorabrechnung 预算（A 档 Startwertwahl 语义）。

    沿床层累积停留时间 ``tau_cumulative_s = xi * tau_bed`` 驱动干燥/DAEM；
    ``o2_remaining_mol_s`` 由累积 VM 碳释放化学计量消耗，取代 ``exp(-5·xi)`` 启发式。

    Ref: Hamel (1999) Vorabrechnung / Startwertwahl 思想；thesis-aligned 工程实现。
    """

    cell_index: int
    axial_xi: float  # [-] 0=床底中心 … 1=床顶/自由段
    tau_cumulative_s: float  # [s]
    t_budget_K: float  # [K] 沿床焓步进物理预算；march 模式下 NR 初温见 derive_nr_temperature_seed_from_vorab_budgets
    x_dry_cumulative: float  # [-]
    x_vm_cumulative: float  # [-]
    o2_remaining_mol_s: float  # [mol/s]
    nC_vm_mol_s: float  # [mol/s] 累积至该高度的 VM 碳
    nH_vm_mol_s: float  # [mol/s]
    nO_vm_mol_s: float  # [mol/s]


def _bed_cell_indices(cells: List[Cell]) -> list[int]:
    return [i for i, cell in enumerate(cells) if str(getattr(cell, "cell_type", "bed")) == "bed"]


def _resolve_product_bed_cell_indices(cells: List[Cell], bed_cell_limit: int) -> list[int]:
    """bed_cell_limit<=0 表示全床层；>0 表示自 bed0 起前 N 格。"""
    bed_indices = _bed_cell_indices(cells)
    limit = int(bed_cell_limit)
    if limit <= 0:
        return bed_indices
    return bed_indices[:limit]


def _rebuild_vorab_budget_vm_state(
    budget: VorabrechnungCellBudget,
    *,
    x_vm_cumulative: float,
    daf_feed: float,
    vm_daf_frac: float,
    c_daf: float,
    h_daf: float,
    o_daf: float,
    O2_feed: float,
    stoich: float,
    M_C: float,
    M_H: float,
    M_O: float,
) -> VorabrechnungCellBudget:
    """按累积 VM 释放率刷新预算中的元素/O₂ 状态。"""
    x_vm = float(np.clip(x_vm_cumulative, 0.0, 1.0))
    m_vm_released = float(daf_feed) * float(vm_daf_frac) * x_vm
    nC_vm = m_vm_released * c_daf / M_C
    nH_vm = m_vm_released * h_daf / M_H / 2.0
    nO_vm = m_vm_released * o_daf / M_O / 2.0
    o2_remaining = max(float(O2_feed) - min(float(O2_feed), stoich * float(nC_vm)), 0.0)
    return replace(
        budget,
        x_vm_cumulative=x_vm,
        o2_remaining_mol_s=float(o2_remaining),
        nC_vm_mol_s=float(nC_vm),
        nH_vm_mol_s=float(nH_vm),
        nO_vm_mol_s=float(nO_vm),
    )


def _apply_vm_devolatilization_zone_spread_to_budgets(
    budgets: list[VorabrechnungCellBudget],
    bed_indices: list[int],
    *,
    n_zone_cells: int,
    profile: str,
    daf_feed: float,
    vm_daf_frac: float,
    c_daf: float,
    h_daf: float,
    o_daf: float,
    O2_feed: float,
    stoich: float,
    M_C: float,
    M_H: float,
    M_O: float,
) -> list[VorabrechnungCellBudget]:
    """将 DAEM 总 VM 释放沿床底脱挥发区分区展开。

    Ref: Hamel (1999) Kap.2.1 — Freisetzung … in Abhängigkeit der Reaktorhöhe，
    再映射到 jeder Zelle；Kap.4 DAEM 决定总量。
    """
    n_zone = int(max(n_zone_cells, 0))
    if n_zone <= 0 or not bed_indices:
        return budgets
    x_vm_total = max(float(budgets[i].x_vm_cumulative) for i in bed_indices)
    bed_rank = {int(bi): int(j) for j, bi in enumerate(bed_indices)}
    bed_set = set(bed_indices)
    out: list[VorabrechnungCellBudget] = []
    for i, budget in enumerate(budgets):
        if i not in bed_set:
            out.append(budget)
            continue
        rank = bed_rank[int(i)]
        if rank < n_zone:
            x_new = vm_devolatilization_zone_cumulative_fraction(
                rank,
                n_zone,
                x_vm_total,
                profile=profile,
            )
        else:
            x_new = x_vm_total
        out.append(
            _rebuild_vorab_budget_vm_state(
                budget,
                x_vm_cumulative=x_new,
                daf_feed=float(daf_feed),
                vm_daf_frac=float(vm_daf_frac),
                c_daf=float(c_daf),
                h_daf=float(h_daf),
                o_daf=float(o_daf),
                O2_feed=float(O2_feed),
                stoich=float(stoich),
                M_C=float(M_C),
                M_H=float(M_H),
                M_O=float(M_O),
            )
        )
    return out


def assign_vm_release_increment_caps_from_budgets(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
) -> dict[str, float]:
    """按 Vorab VM 分区预算为 bed cell 设置单格热解释放上限。"""
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    profile = str(getattr(cfg, "vorab_vm_devolatilization_zone_profile_thesis", "linear"))
    bed_indices = _bed_cell_indices(cells)
    for cell in cells:
        cell._vm_release_fraction_cap = None
    if n_zone <= 0 or not bed_indices or not cell_budgets:
        return {"vm_release_caps_assigned": 0.0}
    x_vm_total = max(float(cell_budgets[i].x_vm_cumulative) for i in bed_indices if i < len(cell_budgets))
    assigned = 0.0
    for rank, bi in enumerate(bed_indices):
        if bi >= len(cells) or rank >= n_zone:
            continue
        cap = vm_devolatilization_zone_increment_fraction(
            rank,
            n_zone,
            x_vm_total,
            profile=profile,
        )
        cells[bi]._vm_release_fraction_cap = float(cap)
        if bool(getattr(cells[bi], "_vorab_drying_pyro_sources_frozen", False)):
            assigned += 1.0
            continue
        cells[bi].invalidate_vorabrechnung_cache()
        assigned += 1.0
    return {"vm_release_caps_assigned": float(assigned)}


def _fuel_vm_mass_kg_from_bed0(cells: List[Cell]) -> float:
    """床底进料 VM 总质量 [kg/s]（zu 流）。"""
    bed_indices = _bed_cell_indices(cells)
    if not bed_indices:
        return 0.0
    bed0 = cells[bed_indices[0]]
    return float(max(np.sum(np.maximum(bed0.m_solid_zu[:, S_VM], 0.0)), 0.0))


def _vm_zone_axial_vm_seed_kg(
    *,
    cell_rank: int,
    total_vm_kg: float,
    prev_x_vm_cumulative: float,
) -> float:
    """VM 脱挥发分区内：该格入口 VM 存量 = 总进料 × (1 − 上游累积释放率)。"""
    if total_vm_kg <= 1.0e-20:
        return 0.0
    x_prev = float(np.clip(prev_x_vm_cumulative, 0.0, 1.0))
    return float(max(total_vm_kg * (1.0 - x_prev), 0.0))


def _expected_vm_holdup_kg_vm_zone(
    *,
    cell_rank: int,
    total_vm_kg: float,
    prev_x_vm_cumulative: float,
    x_vm_cumulative: float,
) -> float:
    """VM 分区语义下该格 reactive VM holdup [kg/s]。"""
    vm_seed = _vm_zone_axial_vm_seed_kg(
        cell_rank=cell_rank,
        total_vm_kg=total_vm_kg,
        prev_x_vm_cumulative=prev_x_vm_cumulative,
    )
    x_vm = float(np.clip(x_vm_cumulative, 0.0, 1.0))
    return float(max(vm_seed * (1.0 - x_vm), 0.0))


def refresh_vm_zone_fixed_pyrolysis_sources(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
) -> dict[str, float]:
    """仅刷新 VM 分区 bed 格的固定热解源项缓存（single-shot 对齐，非全床 refresh）。

    thesis single-shot 已冻结源项时跳过（outer/init 不得重算 DAEM）。

    Ref: Hamel Vorabrechnung 脱挥发分区；thesis VM zone ↔ solid reactive 对齐
    """
    from src.solvers.vorabrechnung.vorab_policy import vorab_outer_may_refresh_drying_pyro_sources

    if not vorab_outer_may_refresh_drying_pyro_sources(cfg):
        return {"vm_zone_sources_refreshed": 0.0}
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    bed_indices = _bed_cell_indices(cells)
    if n_zone <= 0 or not bed_indices or not cell_budgets:
        return {"vm_zone_sources_refreshed": 0.0}

    from src.solvers.vorabrechnung.vorab_refresh import (
        refresh_cell_hydrodynamics_for_frozen_inner,
        refresh_cell_vorabrechnung,
    )

    assign_vm_release_increment_caps_from_budgets(cells, cell_budgets, cfg)
    refreshed = 0.0
    for rank, bi in enumerate(bed_indices):
        if rank >= n_zone or bi >= len(cells):
            break
        cell = cells[bi]
        cell.invalidate_vorabrechnung_cache()
        refresh_cell_vorabrechnung(cell, force=True)
        refresh_cell_hydrodynamics_for_frozen_inner(cell)
        refreshed += 1.0
    return {"vm_zone_sources_refreshed": float(refreshed)}


def align_vm_devolatilization_zone_vorab_state(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Any | None = None,
    reseed_holdup_chain_fn: Any | None = None,
    refresh_zone_pyrolysis_sources: bool = False,
) -> dict[str, float]:
    """VM 脱挥发分区：solid reactive holdup + 单格 cap +（可选）分区源项刷新。

    将 Vorab 累积 ``x_vm_cumulative`` 与 bed0 总 VM 进料链式映射到 bed0..n−1，
    供 init Check2 与 outer Abgleich 共用。

    Ref: Hamel (1999) Kapitel 4 脱挥发分区；Vorabrechnung ↔ Zellenmodell Check2
    """
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    if n_zone <= 0 or not cell_budgets:
        return {
            "vm_zone_aligned": 0.0,
            "bed_solid_reactive_cells_profiled": 0.0,
            "vm_release_caps_assigned": 0.0,
            "vm_zone_sources_refreshed": 0.0,
        }

    from src.solvers.vorabrechnung.vorab_cell_mapping import (
        apply_vorabrechnung_bed_solid_reactive_profile_from_budget,
    )

    solid_diag = apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
        cells,
        cell_budgets,
        cfg,
        apply_bc_fn=None,
        reseed_holdup_chain_fn=None,
    )
    cap_diag = assign_vm_release_increment_caps_from_budgets(cells, cell_budgets, cfg)
    if reseed_holdup_chain_fn is not None:
        reseed_holdup_chain_fn()
    if apply_bc_fn is not None:
        apply_bc_fn()

    source_diag = {"vm_zone_sources_refreshed": 0.0}
    if refresh_zone_pyrolysis_sources:
        source_diag = refresh_vm_zone_fixed_pyrolysis_sources(cells, cell_budgets, cfg)
        if apply_bc_fn is not None:
            apply_bc_fn()

    return {
        "vm_zone_aligned": 1.0,
        "bed_solid_reactive_cells_profiled": float(solid_diag.get("bed_solid_reactive_cells_profiled", 0.0)),
        "vm_release_caps_assigned": float(cap_diag.get("vm_release_caps_assigned", 0.0)),
        "vm_zone_sources_refreshed": float(source_diag.get("vm_zone_sources_refreshed", 0.0)),
    }


def derive_nr_temperature_seed_from_vorab_budgets(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    *,
    O2_feed: float,
    cfg: Any | None = None,
) -> np.ndarray:
    """由 Vorab 预算链导出 NR 可反应窗初温剖面（与 ``t_budget_K`` 焓步进解耦）。

    沿程进度由 ``axial_xi`` 与 O₂ 消耗份额混合；床底/床顶锚于配置的反应温度窗，
    不读取 validation JSON，也不把冷态 ``t_budget_K`` 直接当作 NR 初温。

    Source: Hamel (1999) Startwertwahl；thesis 工程：预算物理 / NR seed 双轨
    """
    n = len(cells)
    if n <= 0 or not cell_budgets:
        return np.array([], dtype=np.float64)

    T_bottom = float(getattr(cfg, "vorab_nr_temperature_seed_bottom_K_thesis", 750.0) if cfg else 750.0)
    T_top = float(getattr(cfg, "vorab_nr_temperature_seed_top_K_thesis", 1100.0) if cfg else 1100.0)
    progress_power = float(
        getattr(cfg, "vorab_nr_temperature_seed_progress_power_thesis", 0.65) if cfg else 0.65
    )
    o2_weight = float(
        getattr(cfg, "vorab_nr_temperature_seed_o2_progress_weight_thesis", 0.35) if cfg else 0.35
    )
    T_bottom = float(max(T_bottom, 300.0))
    T_top = float(max(T_top, T_bottom + 1.0))
    progress_power = float(np.clip(progress_power, 0.05, 3.0))
    o2_weight = float(np.clip(o2_weight, 0.0, 1.0))
    o2_feed = float(max(O2_feed, 1e-12))

    bed_indices = _bed_cell_indices(cells)
    T_seed = np.zeros(n, dtype=np.float64)
    last_bed_T = T_bottom
    bed0_cold_anchor = bool(
        getattr(cfg, "vorab_nr_temperature_seed_bed0_cold_anchor_thesis", True) if cfg else True
    )
    bed0_soft = float(
        getattr(cfg, "vorab_nr_temperature_seed_bed0_soft_progress_thesis", 0.0) if cfg else 0.0
    )
    bed0_soft = float(np.clip(bed0_soft, 0.0, 1.0))

    for j, bi in enumerate(bed_indices):
        if bi >= len(cell_budgets):
            continue
        budget = cell_budgets[bi]
        p_xi = float(np.clip(budget.axial_xi, 0.0, 1.0)) ** progress_power
        o2_consumed_frac = float(
            np.clip((o2_feed - float(budget.o2_remaining_mol_s)) / o2_feed, 0.0, 1.0)
        ) ** progress_power
        progress = float(np.clip((1.0 - o2_weight) * p_xi + o2_weight * o2_consumed_frac, 0.0, 1.0))
        if bed0_cold_anchor and j == 0:
            # 空气/分级燃烧：bed0 钉冷入口窗；可选 soft 进度把部分 O₂/ξ 热抬入 bed0，
            # 避免硬钉 700 K 与上层 O₂ 剖面脱节。上层床始终跟 O₂/ξ（不再用床序号跃升）。
            T_cell = float(T_bottom + (T_top - T_bottom) * progress * bed0_soft)
        else:
            T_cell = float(T_bottom + (T_top - T_bottom) * progress)
        T_seed[bi] = T_cell
        last_bed_T = T_cell

    for i in range(n):
        if i not in bed_indices:
            T_seed[i] = last_bed_T

    return T_seed


def _map_bed_temperature_profile_to_cells(
    cells: List[Cell],
    T_bed_profile: np.ndarray,
) -> np.ndarray:
    """床层 T 剖面 → 全 reactor cell 向量（非床层取床顶值）。"""
    bed_indices = _bed_cell_indices(cells)
    n = len(cells)
    T_full = np.zeros(n, dtype=np.float64)
    if not bed_indices:
        fill = float(T_bed_profile[-1]) if len(T_bed_profile) else 900.0
        T_full.fill(fill)
        return T_full
    for i in range(n):
        if i in bed_indices:
            j = bed_indices.index(i)
            T_full[i] = float(T_bed_profile[j]) if j < len(T_bed_profile) else float(T_bed_profile[-1])
        else:
            T_full[i] = float(T_bed_profile[min(len(bed_indices) - 1, len(T_bed_profile) - 1)])
    return T_full


def build_nr_temperature_seed_for_march_mode(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    *,
    T_inlet: float,
    O2_feed: float,
    H2O_feed: float,
    N2_feed: float,
    fuel_feed_kg_s: float,
    C_dry: float,
    H_dry: float,
    moisture_wt: float,
    P: float,
    cfg: Any | None = None,
    T_bottom_cap_K: float = 1325.0,
) -> np.ndarray:
    """march 模式下 NR 初温 seed（与冷态 ``t_budget_K`` 预算链解耦）。

    ``adiabatic_anchor``：``estimate_axial_T_profile`` 绝热锚定（legacy 同款，非 Fig 7.4）
    ``budget_progress``：``derive_nr_temperature_seed_from_vorab_budgets``
    ``blend``：两者按 ``vorab_nr_temperature_seed_adiabatic_blend_thesis`` 混合

    Source: Hamel (1999) Startwertwahl；thesis 双轨
    """
    mode = str(getattr(cfg, "vorab_nr_temperature_seed_mode_thesis", "adiabatic_anchor") if cfg else "adiabatic_anchor")
    bed_indices = _bed_cell_indices(cells)
    n_bed = max(len(bed_indices), 1)

    t_floor_K = float(
        getattr(cfg, "vorab_bottom_temperature_floor_K_thesis", 1050.0) if cfg is not None else 1050.0
    )
    t_cap_K = float(
        getattr(cfg, "vorab_bottom_temperature_cap_K_thesis", T_bottom_cap_K)
        if cfg is not None
        else T_bottom_cap_K
    )
    t_floor_K = float(max(t_floor_K, 300.0))
    t_cap_K = float(max(t_cap_K, t_floor_K + 1.0))
    T_adiabatic_bed = estimate_axial_T_profile(
        n_cells=n_bed,
        T_inlet=float(T_inlet),
        O2_feed=float(O2_feed),
        H2O_feed=float(H2O_feed),
        N2_feed=float(N2_feed),
        fuel_feed_kg_s=float(fuel_feed_kg_s),
        C_dry=float(C_dry),
        H_dry=float(H_dry),
        moisture_wt=float(moisture_wt),
        P=float(P),
        T_bottom_floor_K=t_floor_K,
        T_bottom_cap_K=t_cap_K,
    )
    T_adiabatic = _map_bed_temperature_profile_to_cells(cells, T_adiabatic_bed)

    if mode == "budget_progress":
        return derive_nr_temperature_seed_from_vorab_budgets(
            cells,
            cell_budgets,
            O2_feed=float(O2_feed),
            cfg=cfg,
        )

    if mode == "blend":
        T_budget = derive_nr_temperature_seed_from_vorab_budgets(
            cells,
            cell_budgets,
            O2_feed=float(O2_feed),
            cfg=cfg,
        )
        blend = float(getattr(cfg, "vorab_nr_temperature_seed_adiabatic_blend_thesis", 0.85) if cfg else 0.85)
        blend = float(np.clip(blend, 0.0, 1.0))
        return T_adiabatic * blend + T_budget * (1.0 - blend)

    return T_adiabatic


def _axial_xi_for_bed_cell(
    cells: List[Cell],
    bed_indices: list[int],
    cell_index: int,
    *,
    H_bed: float,
) -> float:
    """床层单元中心无量纲高度 xi ∈ [0, 1]。"""
    if cell_index not in bed_indices or H_bed <= 1e-12:
        return 1.0
    j = bed_indices.index(cell_index)
    z_center = sum(float(cells[bed_indices[k]].geo.dh) for k in range(j)) + 0.5 * float(
        cells[cell_index].geo.dh
    )
    return float(np.clip(z_center / H_bed, 0.0, 1.0))


def build_vorabrechnung_cell_budgets(
    cells: List[Cell],
    *,
    T_inlet: float,
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
    fuel_type: str = "coal",
    P: float = 101325.0,
    T_bottom_floor_K: float = 1050.0,
    T_bottom_cap_K: float = 1325.0,
    o2_per_mol_c_released: float = 0.5,
    cfg: Any | None = None,
) -> tuple[np.ndarray, list[VorabrechnungCellBudget]]:
    """构建与连接矩阵解耦的 Vorabrechnung 单格预算链（A 档 Startwertwahl）。

    1. 床层 T 剖面：``vorab_bed_temperature_march_thesis`` 时先用 ``T_inlet`` 算干燥/热解，
       再沿床焓步进；否则 ``estimate_axial_T_profile`` 绝热锚定
    2. 全床 ``tau_bed = H_bed / u_mf``（``vorabrechnung_tau_for_cell``）
    3. 沿程 ``tau_i = xi_i * tau_bed`` 驱动 CN 干燥与 DAEM 热解
    4. O2 剩余 = feed − min(feed, stoich × 累积 nC_VM)

    Source: Hamel (1999) Vorabrechnung; specs Startwertwahl 语义
    """
    n = len(cells)
    bed_indices = _bed_cell_indices(cells)
    n_bed = max(len(bed_indices), 1)
    H_bed = float(sum(float(cells[i].geo.dh) for i in bed_indices)) if bed_indices else float(
        sum(float(c.geo.dh) for c in cells)
    )

    use_march = cfg is not None and bool(getattr(cfg, "vorab_bed_temperature_march_thesis", False))
    t_floor_K = float(
        getattr(cfg, "vorab_bottom_temperature_floor_K_thesis", T_bottom_floor_K)
        if cfg is not None
        else T_bottom_floor_K
    )
    t_cap_K = float(
        getattr(cfg, "vorab_bottom_temperature_cap_K_thesis", T_bottom_cap_K)
        if cfg is not None
        else T_bottom_cap_K
    )
    t_floor_K = float(max(t_floor_K, 300.0))
    t_cap_K = float(max(t_cap_K, t_floor_K + 1.0))
    if use_march:
        T_bed_profile = np.full(n_bed, float(T_inlet), dtype=np.float64)
    else:
        T_bed_profile = estimate_axial_T_profile(
            n_cells=n_bed,
            T_inlet=float(T_inlet),
            O2_feed=float(O2_feed),
            H2O_feed=float(H2O_feed),
            N2_feed=float(N2_feed),
            fuel_feed_kg_s=float(fuel_feed_kg_s),
            C_dry=float(C_dry),
            H_dry=float(H_dry),
            moisture_wt=float(moisture_wt),
            P=float(P),
            T_bottom_floor_K=t_floor_K,
            T_bottom_cap_K=t_cap_K,
        )
    T_profile = np.zeros(n, dtype=np.float64)
    for i in range(n):
        if i in bed_indices:
            T_profile[i] = float(T_bed_profile[bed_indices.index(i)])
        else:
            T_profile[i] = float(T_bed_profile[-1])

    tau_bed = float(vorabrechnung_tau_for_cell(cells[bed_indices[0]] if bed_indices else cells[0]))

    moisture_frac = float(moisture_wt) / 100.0
    ash_frac = float(ash_dry_wt) / 100.0
    vm_daf_frac = float(VM_daf) / 100.0
    dry_feed = float(fuel_feed_kg_s) * (1.0 - moisture_frac)
    daf_feed = dry_feed * (1.0 - ash_frac)

    M_C, M_H, M_O = 12.011e-3, 1.00794e-3, 15.999e-3
    c_dry = float(C_dry) / 100.0
    h_dry = float(H_dry) / 100.0
    o_dry = float(O_dry) / 100.0
    to_daf = 1.0 / max(1.0 - ash_frac, 1e-9)
    c_daf = c_dry * to_daf
    h_daf = h_dry * to_daf
    o_daf = o_dry * to_daf

    daem_fuel = "brown_coal"
    if fuel_type in ("wood", "biomass"):
        daem_fuel = "wood"
    elif fuel_type in ("coal", "brown_coal"):
        daem_fuel = "brown_coal"

    stoich = float(max(o2_per_mol_c_released, 0.0))
    budgets: list[VorabrechnungCellBudget] = []
    base_gas_mol_s = float(max(O2_feed + H2O_feed + N2_feed, 1e-6))

    march_params = {}
    if use_march and cfg is not None:
        march_params = {
            "heat_loss_frac": float(getattr(cfg, "vorab_bed_temperature_march_heat_loss_frac_thesis", 0.05)),
            "dH_comb_J_per_mol_O2": float(
                getattr(cfg, "vorab_bed_temperature_march_dH_comb_J_per_mol_O2_thesis", 393_500.0)
            ),
            "drying_enthalpy_J_per_kg_moisture": float(
                getattr(cfg, "vorab_bed_temperature_march_drying_enthalpy_J_per_kg_thesis", 2.26e6)
            ),
            "pyrolysis_enthalpy_J_per_kg_daf_vm": float(
                getattr(
                    cfg,
                    "vorab_bed_temperature_march_pyrolysis_enthalpy_J_per_kg_daf_vm_thesis",
                    1.0e6,
                )
            ),
            "Cp_gas_J_per_mol_K": float(
                getattr(cfg, "vorab_bed_temperature_march_Cp_gas_J_per_mol_K_thesis", 32.0)
            ),
            "T_floor_K": float(max(T_inlet, 300.0)),
            "T_cap_K": float(
                max(
                    getattr(cfg, "vorab_bottom_temperature_cap_K_thesis", T_bottom_cap_K),
                    max(T_inlet, 300.0) + 1.0,
                )
            ),
        }

    if use_march and bed_indices:
        prev_o2_consumed = 0.0
        prev_x_dry = 0.0
        prev_x_vm = 0.0
        T_entry = float(T_inlet)
        T_floor = float(march_params["T_floor_K"])
        T_cap = float(march_params["T_cap_K"])

        for j, bi in enumerate(bed_indices):
            cell = cells[bi]
            xi = _axial_xi_for_bed_cell(cells, bed_indices, bi, H_bed=H_bed)
            tau_i = float(max(xi * tau_bed, float(cell.geo.dh) / max(float(cell.u_mf), 1e-3)))
            solid_T_init = float(max(T_inlet, 293.15))
            budget = _compute_vorabrechnung_cell_budget(
                cell=cell,
                cell_index=bi,
                axial_xi=xi,
                tau_i=tau_i,
                T_i=T_entry,
                solid_T_init=solid_T_init,
                moisture_wt=float(moisture_wt),
                vm_daf_frac=float(vm_daf_frac),
                daem_fuel=daem_fuel,
                daf_feed=float(daf_feed),
                c_daf=c_daf,
                h_daf=h_daf,
                o_daf=o_daf,
                O2_feed=float(O2_feed),
                stoich=stoich,
                M_C=M_C,
                M_H=M_H,
                M_O=M_O,
            )
            o2_consumed = float(max(O2_feed - budget.o2_remaining_mol_s, 0.0))
            o2_inc = float(max(o2_consumed - prev_o2_consumed, 0.0))
            dry_inc = float(max(budget.x_dry_cumulative - prev_x_dry, 0.0))
            vm_inc = float(max(budget.x_vm_cumulative - prev_x_vm, 0.0))
            pyro_gas_mol_s = float(max(budget.nC_vm_mol_s + budget.nH_vm_mol_s * 0.5, 0.0))
            dT = _vorab_march_enthalpy_delta_T(
                o2_inc_mol_s=o2_inc,
                dry_inc=dry_inc,
                vm_inc=vm_inc,
                fuel_feed_kg_s=float(fuel_feed_kg_s),
                moisture_wt=float(moisture_wt),
                daf_feed_kg_s=float(daf_feed),
                vm_daf_frac=float(vm_daf_frac),
                pyro_gas_mol_s=pyro_gas_mol_s,
                base_gas_mol_s=base_gas_mol_s,
                **{k: march_params[k] for k in march_params if k not in ("T_floor_K", "T_cap_K")},
            )
            T_out = float(np.clip(T_entry + dT, T_floor, T_cap))
            budget = replace(budget, t_budget_K=T_out)
            T_profile[bi] = T_out
            budgets.append(budget)
            T_entry = T_out
            prev_o2_consumed = o2_consumed
            prev_x_dry = float(budget.x_dry_cumulative)
            prev_x_vm = float(budget.x_vm_cumulative)

        budget_by_index = {budget.cell_index: budget for budget in budgets}
        last_bed_T = float(T_profile[bed_indices[-1]])
        for i, cell in enumerate(cells):
            if i in bed_indices:
                continue
            xi = 1.0
            tau_i = float(max(tau_bed, float(cell.geo.dh) / max(float(cell.u_mf), 1e-3)))
            solid_T_init = max(float(getattr(cell, "T_in_solid", last_bed_T)), 293.15)
            budget = _compute_vorabrechnung_cell_budget(
                cell=cell,
                cell_index=i,
                axial_xi=xi,
                tau_i=tau_i,
                T_i=last_bed_T,
                solid_T_init=solid_T_init,
                moisture_wt=float(moisture_wt),
                vm_daf_frac=float(vm_daf_frac),
                daem_fuel=daem_fuel,
                daf_feed=float(daf_feed),
                c_daf=c_daf,
                h_daf=h_daf,
                o_daf=o_daf,
                O2_feed=float(O2_feed),
                stoich=stoich,
                M_C=M_C,
                M_H=M_H,
                M_O=M_O,
            )
            budget_by_index[i] = replace(budget, t_budget_K=last_bed_T)
            T_profile[i] = last_bed_T
        budgets = [budget_by_index[i] for i in range(n)]
    else:
        for i, cell in enumerate(cells):
            xi = _axial_xi_for_bed_cell(cells, bed_indices, i, H_bed=H_bed) if i in bed_indices else 1.0
            tau_i = float(max(xi * tau_bed, float(cell.geo.dh) / max(float(cell.u_mf), 1e-3)))
            T_i = float(T_profile[i])
            solid_T_init = max(float(getattr(cell, "T_in_solid", T_profile[0])), 293.15)
            budgets.append(
                _compute_vorabrechnung_cell_budget(
                    cell=cell,
                    cell_index=i,
                    axial_xi=xi,
                    tau_i=tau_i,
                    T_i=T_i,
                    solid_T_init=solid_T_init,
                    moisture_wt=float(moisture_wt),
                    vm_daf_frac=float(vm_daf_frac),
                    daem_fuel=daem_fuel,
                    daf_feed=float(daf_feed),
                    c_daf=c_daf,
                    h_daf=h_daf,
                    o_daf=o_daf,
                    O2_feed=float(O2_feed),
                    stoich=stoich,
                    M_C=M_C,
                    M_H=M_H,
                    M_O=M_O,
                )
            )

    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)) if cfg is not None else 0
    if n_zone > 0 and bed_indices:
        profile = str(getattr(cfg, "vorab_vm_devolatilization_zone_profile_thesis", "linear"))
        budgets = _apply_vm_devolatilization_zone_spread_to_budgets(
            budgets,
            bed_indices,
            n_zone_cells=n_zone,
            profile=profile,
            daf_feed=float(daf_feed),
            vm_daf_frac=float(vm_daf_frac),
            c_daf=float(c_daf),
            h_daf=float(h_daf),
            o_daf=float(o_daf),
            O2_feed=float(O2_feed),
            stoich=float(stoich),
            M_C=float(M_C),
            M_H=float(M_H),
            M_O=float(M_O),
        )

    return T_profile, budgets



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

