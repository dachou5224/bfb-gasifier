"""Vorabrechnung — vorab cell mapping segment (OPT-004)."""
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
from src.core.connectivity import align_bottom_primary_gas_state_to_inlet_split
from src.solvers.vorabrechnung.vorab_budgets import (
    VorabrechnungCellBudget,
    _bed_cell_indices,
    _expected_vm_holdup_kg_vm_zone,
    _fuel_vm_mass_kg_from_bed0,
    _vm_zone_axial_vm_seed_kg,
)
from src.solvers.vorabrechnung.vorab_gibbs_seed import (
    _legacy_fixed_ratio_pyrolysis_increment_not_hamel,
    solve_pyrolysis_gas_gibbs_composition,
)

def _hydrodynamic_bubble_share_for_cell(cell: Cell) -> float:
    """可见气相 holdup 中的气泡份额（Vorab 水力学 → 两相分配）。"""
    eps_b = float(getattr(cell, "eps_b", np.nan))
    eps_d_void = float(getattr(cell, "eps_d_voidage", np.nan))
    if np.isfinite(eps_b) and np.isfinite(eps_d_void):
        gas_holdup_total = max(eps_b + eps_d_void, 1e-12)
        return float(np.clip(eps_b / gas_holdup_total, 0.05, 0.95))
    return 0.30


_POST_FASTOX_SYNGAS_PASSTHROUGH_SPECIES: tuple[str, ...] = (
    "CO",
    "CO2",
    "H2",
    "CH4",
    "TAR1",
    "TAR2",
)


def _apply_post_fastox_post_staged_passthrough(
    cells: List[Cell],
    cfg: Any,
    species: tuple[str, ...],
) -> dict[str, float]:
    """无 O₂ 补气床格：把指定物种 holdup 贴到 inflow+zu+rez，再按水力两相分。

    跳过 bed0 与 ``zu_O2>0`` 的格，避免冲掉 fastox / R1 种子。
    不改 k / K_bd / Eq.2.7。

    Ref: Hamel (1999) Eq.2.1–2.2 Startwertwahl; handoff §5.62–5.63.
    """
    from src.core.connectivity import propagate_upstream

    bed_indices = _bed_cell_indices(cells)
    if not bed_indices:
        return {"enabled": 1.0, "cells_updated": 0.0}
    idx = GAS_SPECIES_INDEX
    updated = 0.0
    for i in bed_indices[1:]:
        propagate_upstream(cells, cfg, i)
        cell = cells[i]
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        zu_o2 = float(max(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]], 0.0))
        if zu_o2 > 1.0e-15:
            continue
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)
        wrote = False
        for sp in species:
            j_sp = int(idx[sp])
            total = float(
                max(
                    cell.N_d_in[j_sp]
                    + cell.N_b_in[j_sp]
                    + cell.N_zu_d[j_sp]
                    + cell.N_zu_b[j_sp]
                    + cell.N_rez_d[j_sp]
                    + cell.N_rez_b[j_sp],
                    0.0,
                )
            )
            if total <= 1.0e-15:
                continue
            cell.N_d[j_sp] = max(total * dense_share, 1e-12)
            cell.N_b[j_sp] = max(total * bubble_share, 1e-12)
            wrote = True
        if not wrote:
            continue
        cell._thermo_cache_valid = False
        updated += 1.0
        for j in range(i + 1, len(cells)):
            propagate_upstream(cells, cfg, j)
    return {"enabled": 1.0, "cells_updated": float(updated)}


def apply_post_fastox_post_staged_h2o_passthrough(
    cells: List[Cell],
    cfg: Any,
) -> dict[str, float]:
    """fastox 之后再贴一次无 O₂ 补气床格的 H₂O 透传。

    流映射在 fastox 之前：bed2 氧化后 ``N_H2O`` 下降，bed3+ 仍停在轴向
    ``max(inflow, current)`` 的偏高值（flag 开着时 init gap 仍 ≈+2.3 mol/s）。
    仅 ``vorab_a_tier_post_staged_h2o_passthrough_thesis`` 打开时生效。
    不改 k / K_bd / Eq.2.7。

    Ref: Hamel (1999) Eq.2.1–2.2 Startwertwahl; handoff §5.28, §5.62.
    """
    if not bool(getattr(cfg, "vorab_a_tier_post_staged_h2o_passthrough_thesis", False)):
        return {"enabled": 0.0, "cells_updated": 0.0}
    return _apply_post_fastox_post_staged_passthrough(cells, cfg, ("H2O",))


def apply_post_fastox_post_staged_syngas_passthrough(
    cells: List[Cell],
    cfg: Any,
) -> dict[str, float]:
    """fastox 之后再贴一次无 O₂ 补气床格的 CO/CO₂/H₂/CH₄/TAR。

    ``seed_major_product_holdup`` 在 fastox 之前，且 zone=3 时 limit 只到 bed2。
    fastox 改 bed2 后，bed3+ 仍停在映射时的 Gibbs 拷贝（CO₂ 偏高、CH₄/TAR 归零）。
    仅 ``vorab_a_tier_post_staged_syngas_passthrough_thesis`` 打开时生效。
    不改 H₂O / O₂ / N₂。不改 k / K_bd / Eq.2.7。勿与 H₂O 透传叠开验收。

    Ref: Hamel (1999) Eq.2.1–2.2 Startwertwahl; handoff §5.63.
    """
    if not bool(getattr(cfg, "vorab_a_tier_post_staged_syngas_passthrough_thesis", False)):
        return {"enabled": 0.0, "cells_updated": 0.0}
    return _apply_post_fastox_post_staged_passthrough(
        cells, cfg, _POST_FASTOX_SYNGAS_PASSTHROUGH_SPECIES
    )


def _incremental_pyrolysis_major_mol_s(
    budget_curr: VorabrechnungCellBudget,
    budget_prev: VorabrechnungCellBudget,
    *,
    T_K: float | None = None,
    P_Pa: float | None = None,
    pyrolysis_gibbs_solver_mode: str = "hamel_reduced",
    pyrolysis_gibbs_temperature_floor_K: float = 750.0,
    fuel_type: str = "coal",
    pyrolysis_tar_carbon_frac: float = 0.2,
    pyrolysis_koks_carbon_frac: float = 0.0,
    pyrolysis_tar_yield_mode: str = "hegermann",
    n_daf_per_c: float = 0.0,
    s_daf_per_c: float = 0.0,
    sulfur_volatile_frac: float = 0.5,
    target_tar_hc_ratio: float | None = None,
    pyrolysis_ch4_to_h2_co_split_frac: float = 0.0,
    pyrolysis_co_reduce_h2_inject_frac: float = 0.0,
) -> dict[str, float]:
    """相邻床格间 Vorab 热解增量 → 主组分摩尔流 [mol/s]（§4.3 Gibbs 累积差分）。

    coal + ``hegermann`` 时用 Eq.4.21–4.22 的 ``x_C,Koks`` / ``x_C,Teer``，
    与 ``calc_drying_pyrolysis_sources`` 同口径，避免固定 tar=0.2、koks=0
    使 Teer 后气相 ``An=b`` 不可行。
    """
    dC = float(max(budget_curr.nC_vm_mol_s - budget_prev.nC_vm_mol_s, 0.0))
    dH = float(max(budget_curr.nH_vm_mol_s - budget_prev.nH_vm_mol_s, 0.0))
    dO = float(max(budget_curr.nO_vm_mol_s - budget_prev.nO_vm_mol_s, 0.0))
    if dC <= 1.0e-20 and dH <= 1.0e-20 and dO <= 1.0e-20:
        return {sp: 0.0 for sp in ("CO", "CO2", "H2", "CH4", "H2O", "O2", "N2")}
    if T_K is None or P_Pa is None:
        raise ValueError("Hamel strict: pyrolysis increment requires T_K and P_Pa for §4.3 Gibbs")
    from src.core.cell_pyrolysis import incremental_pyrolysis_products_from_budgets

    tar_frac = float(pyrolysis_tar_carbon_frac)
    koks_frac = float(pyrolysis_koks_carbon_frac)
    mode = str(pyrolysis_tar_yield_mode).strip().lower()
    if mode == "hegermann" and str(fuel_type).strip().lower() == "coal":
        from src.thermal.hegermann_yield import xc_koks_hegermann, xc_teer_hegermann

        tar_frac = float(xc_teer_hegermann(float(P_Pa), float(T_K)))
        koks_frac = float(xc_koks_hegermann(float(P_Pa), float(T_K)))

    delta = incremental_pyrolysis_products_from_budgets(
        budget_curr,
        budget_prev,
        T_K=float(T_K),
        P_Pa=float(P_Pa),
        fuel_type=str(fuel_type),
        pyrolysis_tar_carbon_frac=tar_frac,
        pyrolysis_koks_carbon_frac=koks_frac,
        n_daf_per_c=float(n_daf_per_c),
        s_daf_per_c=float(s_daf_per_c),
        sulfur_volatile_frac=float(sulfur_volatile_frac),
        target_tar_hc_ratio=target_tar_hc_ratio,
        pyrolysis_gibbs_solver_mode=str(pyrolysis_gibbs_solver_mode),
        pyrolysis_gibbs_temperature_floor_K=float(pyrolysis_gibbs_temperature_floor_K),
        pyrolysis_ch4_to_h2_co_split_frac=float(pyrolysis_ch4_to_h2_co_split_frac),
        pyrolysis_co_reduce_h2_inject_frac=float(pyrolysis_co_reduce_h2_inject_frac),
    )
    return {
        sp: float(max(delta.get(sp, 0.0), 0.0))
        for sp in ("CO", "CO2", "H2", "CH4", "H2O", "O2", "N2")
    }


def _seed_bed_solid_reactive_holdup_from_vorab_budget(
    cell: Cell,
    budget: VorabrechnungCellBudget,
    *,
    cfg: Any | None = None,
    cell_rank: int | None = None,
    bed_indices: list[int] | None = None,
    cell_budgets: list[VorabrechnungCellBudget] | None = None,
    total_vm_kg: float | None = None,
) -> bool:
    """按 Vorab 预算刷新单格 reactive solid（VM/水分）holdup，不动 char/ash。

    char/ash 由 ``seed_bed_holdup_chain_from_transport`` 维护；此处只闭合
    干燥/热解进度 ``x_dry_cumulative`` / ``x_vm_cumulative`` 与 active inflow 预算。

    Ref: Hamel (1999) Vorabrechnung 固相 zu/in 与干燥/热解分配
    """
    if str(getattr(cell, "solid_state_model", "legacy_stream")) not in {
        "holdup_transport",
        "freeboard_closure",
    }:
        return False

    nk = int(cell.solid.n_size_classes)
    active_seed = np.maximum(cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in, 0.0)
    vm_seed = float(np.sum(active_seed[:, S_VM]))
    moist_seed = float(np.sum(active_seed[:, S_MOISTURE]))

    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)) if cfg is not None else 0
    use_vm_zone = (
        n_zone > 0
        and cell_rank is not None
        and bed_indices is not None
        and cell_budgets is not None
        and total_vm_kg is not None
        and int(cell_rank) < n_zone
    )
    if not use_vm_zone and vm_seed + moist_seed <= 1.0e-20:
        return False

    def _class_fraction(comp_idx: int) -> np.ndarray:
        col = np.maximum(active_seed[:, comp_idx], 0.0)
        total = float(np.sum(col))
        if total <= 1.0e-12:
            feed_frac = np.asarray(cell.solid.mass_fractions, dtype=np.float64)
            feed_sum = float(np.sum(feed_frac))
            if feed_sum > 1.0e-12:
                return feed_frac / feed_sum
            return np.full(nk, 1.0 / max(nk, 1), dtype=np.float64)
        return col / total

    vm_frac = _class_fraction(S_VM)
    moist_frac = _class_fraction(S_MOISTURE)
    x_dry_cum = float(np.clip(budget.x_dry_cumulative, 0.0, 1.0))
    x_vm_cum = float(np.clip(budget.x_vm_cumulative, 0.0, 1.0))

    if use_vm_zone:
        rank = int(cell_rank)
        prev_x_vm = (
            0.0
            if rank <= 0
            else float(np.clip(cell_budgets[bed_indices[rank - 1]].x_vm_cumulative, 0.0, 1.0))
        )
        vm_seed = _vm_zone_axial_vm_seed_kg(
            cell_rank=rank,
            total_vm_kg=float(total_vm_kg),
            prev_x_vm_cumulative=prev_x_vm,
        )
        vm_mass = _expected_vm_holdup_kg_vm_zone(
            cell_rank=rank,
            total_vm_kg=float(total_vm_kg),
            prev_x_vm_cumulative=prev_x_vm,
            x_vm_cumulative=x_vm_cum,
        )
        if vm_seed <= 1.0e-20 and moist_seed <= 1.0e-20:
            return False
        vm_frac = np.full(nk, 1.0 / max(nk, 1), dtype=np.float64)
    else:
        vm_mass = max(vm_seed * (1.0 - x_vm_cum), 0.0)

    moist_mass = max(moist_seed * (1.0 - x_dry_cum), 0.0)

    cell.m_solid[:, S_VM] = vm_mass * vm_frac
    cell.m_solid[:, S_MOISTURE] = moist_mass * moist_frac
    return True


def apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Any | None = None,
    reseed_holdup_chain_fn: Any | None = None,
) -> dict[str, float]:
    """床层 reactive solid（VM/水分）holdup 对齐 Vorab 干燥/热解预算剖面。"""
    if not bool(getattr(cfg, "vorab_a_tier_bed_solid_reactive_profile_thesis", True)):
        return {"bed_solid_reactive_cells_profiled": 0.0}

    bed_indices = _bed_cell_indices(cells)
    profiled = 0.0
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    total_vm_kg = _fuel_vm_mass_kg_from_bed0(cells) if n_zone > 0 else 0.0
    for rank, i in enumerate(bed_indices):
        if i >= len(cell_budgets):
            continue
        if _seed_bed_solid_reactive_holdup_from_vorab_budget(
            cells[i],
            cell_budgets[i],
            cfg=cfg,
            cell_rank=rank,
            bed_indices=bed_indices,
            cell_budgets=cell_budgets,
            total_vm_kg=total_vm_kg,
        ):
            profiled += 1.0

    if profiled > 0.0 and reseed_holdup_chain_fn is not None:
        reseed_holdup_chain_fn()
    if apply_bc_fn is not None:
        apply_bc_fn()
    return {"bed_solid_reactive_cells_profiled": float(profiled)}


def apply_vorabrechnung_bed_stream_mapping(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    upstream_closure: bool | None = None,
) -> dict[str, float]:
    """将 Vorab zu/ab 流映射到床层气相 holdup（Hamel zu/ab P0 原型）。

    默认（``upstream_closure=False``）：bed0 主气化剂总量锚定 ``N_zu_d+N_zu_b``，
    两相分配按本地水力学 ``eps_b/(eps_b+eps_d_void)``（非对称复制 ``N_d=N_zu_d``）。
    实验路径（``upstream_closure=True``）：上游床格主气沿 ``N_d_in+N_b_in`` 链式映射
    并做 major phase-split 预投影；产物 holdup 保留 ``generate_initial_x0`` 结果。

    Source: Hamel (1999) Kap.2 p.16 Vorabrechnung → Zellen zu/ab 映射
    """
    from src.core.connectivity import (
        preproject_bed_major_gas_phase_split_to_exchange_closure,
        propagate_upstream,
        set_bottom_cell_feeds,
    )
    from src.core.connectivity.gas_inlet import resolve_bottom_o2_dense_fraction

    bed_indices = _bed_cell_indices(cells)
    if not bed_indices or not cell_budgets:
        return {"bed_stream_cells_mapped": 0.0, "upstream_closure_passes": 0.0}

    if upstream_closure is None:
        upstream_closure = bool(getattr(cfg, "vorab_a_tier_bed_stream_upstream_closure_thesis", False))

    set_bottom_cell_feeds(cells, cfg)
    idx = GAS_SPECIES_INDEX
    primary_species = ("O2", "H2O", "N2")

    bed0 = cells[bed_indices[0]]

    def _assign_primary_holdup_with_hydro(cell: Cell, *, include_inflow: bool) -> None:
        """主气 holdup = 本格 ``N_zu``（+ 可选上游 ``N_in+N_rez``）按水力两相分配。

        分级 O2 补气格保持 ``include_inflow=False``（仅 ``N_zu``），与 axial O2 blend 分工；
        H2O/N2 累积闭合由 axial / inert passthrough 负责（handoff §5.27）。
        """
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)
        o2_dense = float(resolve_bottom_o2_dense_fraction(cell, cfg))
        for sp in primary_species:
            j = int(idx[sp])
            total = float(max(cell.N_zu_d[j] + cell.N_zu_b[j], 0.0))
            if include_inflow:
                total += float(
                    max(
                        cell.N_d_in[j]
                        + cell.N_b_in[j]
                        + cell.N_rez_d[j]
                        + cell.N_rez_b[j],
                        0.0,
                    )
                )
            if total <= 1.0e-15:
                total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
            if total <= 1.0e-15:
                continue
            d_frac = o2_dense if sp == "O2" else dense_share
            cell.N_d[j] = max(total * d_frac, 1e-12)
            cell.N_b[j] = max(total * (1.0 - d_frac), 1e-12)
        cell._thermo_cache_valid = False

    post_staged_h2o = bool(
        getattr(cfg, "vorab_a_tier_post_staged_h2o_passthrough_thesis", False)
    )

    def _assign_post_staged_inert_passthrough(cell: Cell) -> None:
        """无本地 O2 补气的床格：N2 透传；H2O 仅在 thesis 旗标下透传。

        N2：§5.27 默认开。H2O：§5.28 opt-in（``vorab_a_tier_post_staged_h2o_passthrough_thesis``），
        须配 LSQ2+延迟焓透传；否则易落入 CO≈0.21 坏吸引子（§5.29）。
        """
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)
        species = ("N2", "H2O") if post_staged_h2o else ("N2",)
        for sp in species:
            j = int(idx[sp])
            total = float(
                max(
                    cell.N_d_in[j]
                    + cell.N_b_in[j]
                    + cell.N_zu_d[j]
                    + cell.N_zu_b[j]
                    + cell.N_rez_d[j]
                    + cell.N_rez_b[j],
                    0.0,
                )
            )
            if total <= 1.0e-15:
                continue
            cell.N_d[j] = max(total * dense_share, 1e-12)
            cell.N_b[j] = max(total * bubble_share, 1e-12)
        cell._thermo_cache_valid = False

    # 自下而上：先传播入口。有 O2 补气的格仍只锚本地 ``N_zu``（氧化剖面由 axial O2
    # blend 负责；勿改成 inflow+zu，否则中段 O2 holdup 暴涨、塌出口 CO）。
    # 无 O2 补气的格：N2 透传（§5.27）；H2O 见 ``post_staged_h2o`` 旗标。
    staged_anchors = 0.0
    _assign_primary_holdup_with_hydro(bed0, include_inflow=False)
    if float(bed0.N_zu_d[idx["O2"]] + bed0.N_zu_b[idx["O2"]]) > 1.0e-15:
        staged_anchors += 1.0
    for bi in bed_indices[1:]:
        propagate_upstream(cells, cfg, bi)
        cell = cells[bi]
        zu_o2 = float(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]])
        if zu_o2 > 1.0e-15:
            _assign_primary_holdup_with_hydro(cell, include_inflow=False)
            staged_anchors += 1.0
        else:
            _assign_post_staged_inert_passthrough(cell)

    if staged_anchors <= 0.0:
        _assign_primary_holdup_with_hydro(bed0, include_inflow=False)

    for i in range(1, len(cells)):
        propagate_upstream(cells, cfg, i)

    axial_diag = apply_vorabrechnung_bed_axial_primary_profile_from_budget(
        cells,
        cell_budgets,
        cfg,
    )

    if not upstream_closure:
        return {
            "bed_stream_cells_mapped": float(max(staged_anchors, 1.0)),
            "staged_primary_gas_anchors": float(staged_anchors),
            "upstream_closure_passes": 0.0,
            **axial_diag,
        }

    closure_passes = 0
    product_species = ("CO", "CO2", "H2", "CH4")
    for k, i in enumerate(bed_indices[1:], start=1):
        cell = cells[i]
        saved_products = {
            sp: (
                float(cell.N_d[int(idx[sp])]),
                float(cell.N_b[int(idx[sp])]),
            )
            for sp in product_species
        }
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)

        for sp in primary_species:
            j = int(idx[sp])
            inflow_total = float(max(cell.N_d_in[j] + cell.N_b_in[j], 0.0))
            total = max(inflow_total, 1e-12)
            cell.N_d[j] = max(total * dense_share, 1e-12)
            cell.N_b[j] = max(total * bubble_share, 1e-12)

        for sp in product_species:
            j = int(idx[sp])
            cell.N_d[j] = max(saved_products[sp][0], 1e-12)
            cell.N_b[j] = max(saved_products[sp][1], 1e-12)

        cell._thermo_cache_valid = False

        if upstream_closure:
            if preproject_bed_major_gas_phase_split_to_exchange_closure(
                cells,
                cfg,
                cell_index=i,
                preserve_primary_inlet_totals=True,
            ):
                closure_passes += 1

        for j in range(i + 1, len(cells)):
            propagate_upstream(cells, cfg, j)

    return {
        "bed_stream_cells_mapped": float(len(bed_indices)),
        "upstream_closure_passes": float(closure_passes),
        **axial_diag,
    }


def apply_vorabrechnung_bed_axial_primary_profile_from_budget(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    bed_cell_limit: int = 0,
    blend_factor: float | None = None,
) -> dict[str, float]:
    """沿床层将主气化剂 holdup 对齐 Vorab 预算剖面（A 档 zu/ab 扩展）。

    - bed0：O₂/H₂O/N₂ 保持 ``N_zu_d/N_zu_b`` 锚定
    - bed i>0：O₂ 总量温和混合 toward ``o2_remaining_mol_s``（默认 35% 步长）；
      H₂O = ``max(inflow, current)``；若 ``post_staged_h2o_passthrough`` 且无 O₂ 补气则 inflow+zu+rez；
      N₂ = inflow + ``N_zu`` + ``N_rez``（禁止 ``max(inflow, current)`` 保留过填）

    Source: Hamel (1999) Vorabrechnung 沿程 O₂ 消耗预算 → Zellen holdup
    """
    from src.core.connectivity import propagate_upstream

    if not bool(getattr(cfg, "vorab_a_tier_bed_axial_primary_profile_thesis", True)):
        return {"bed_axial_primary_cells_profiled": 0.0}

    bed_indices = _bed_cell_indices(cells)
    if len(bed_indices) <= 1 or not cell_budgets:
        return {"bed_axial_primary_cells_profiled": 0.0}

    limit_cfg = int(getattr(cfg, "vorab_a_tier_bed_axial_primary_profile_cells_thesis", 0))
    if bed_cell_limit > 0:
        limit_cfg = bed_cell_limit if limit_cfg <= 0 else min(limit_cfg, bed_cell_limit)
    n_target = len(bed_indices) if limit_cfg <= 0 else min(limit_cfg, len(bed_indices))
    target_indices = bed_indices[1:n_target] if n_target > 1 else []
    if not target_indices:
        return {"bed_axial_primary_cells_profiled": 0.0}

    if blend_factor is None:
        blend_factor = float(getattr(cfg, "vorab_a_tier_bed_axial_primary_blend_thesis", 0.35))
    alpha = float(np.clip(blend_factor, 0.05, 1.0))

    idx = GAS_SPECIES_INDEX
    profiled = 0.0

    for i in range(1, len(cells)):
        propagate_upstream(cells, cfg, i)

    for i in target_indices:
        if i >= len(cell_budgets):
            continue
        cell = cells[i]
        budget = cell_budgets[i]
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)

        j_o2 = int(idx["O2"])
        o2_target = float(max(budget.o2_remaining_mol_s, 1e-12))
        o2_current = float(max(cell.N_d[j_o2] + cell.N_b[j_o2], 0.0))
        o2_blended = max(o2_current + alpha * (o2_target - o2_current), 1e-12)
        cell.N_d[j_o2] = max(o2_blended * dense_share, 1e-12)
        cell.N_b[j_o2] = max(o2_blended * bubble_share, 1e-12)

        # H2O：默认 max(inflow, current)。仅当 post-staged H2O thesis 且无 O2 补气时
        # 改为 inflow+zu+rez（§5.28–5.29）。
        # N2：一律 inflow+zu+rez（§5.27）。
        j_h2o = int(idx["H2O"])
        inflow_h2o = float(max(cell.N_d_in[j_h2o] + cell.N_b_in[j_h2o], 0.0))
        local_zu_h2o = float(max(cell.N_zu_d[j_h2o] + cell.N_zu_b[j_h2o], 0.0))
        local_rez_h2o = float(max(cell.N_rez_d[j_h2o] + cell.N_rez_b[j_h2o], 0.0))
        zu_o2 = float(max(cell.N_zu_d[idx["O2"]] + cell.N_zu_b[idx["O2"]], 0.0))
        post_h2o = bool(
            getattr(cfg, "vorab_a_tier_post_staged_h2o_passthrough_thesis", False)
        )
        if post_h2o and zu_o2 <= 1.0e-15:
            total_h2o = max(inflow_h2o + local_zu_h2o + local_rez_h2o, 1e-12)
        else:
            current_h2o = float(max(cell.N_d[j_h2o] + cell.N_b[j_h2o], 0.0))
            total_h2o = max(inflow_h2o, current_h2o, 1e-12)
        cell.N_d[j_h2o] = max(total_h2o * dense_share, 1e-12)
        cell.N_b[j_h2o] = max(total_h2o * bubble_share, 1e-12)

        j_n2 = int(idx["N2"])
        inflow_n2 = float(max(cell.N_d_in[j_n2] + cell.N_b_in[j_n2], 0.0))
        local_zu_n2 = float(max(cell.N_zu_d[j_n2] + cell.N_zu_b[j_n2], 0.0))
        local_rez_n2 = float(max(cell.N_rez_d[j_n2] + cell.N_rez_b[j_n2], 0.0))
        total_n2 = max(inflow_n2 + local_zu_n2 + local_rez_n2, 1e-12)
        cell.N_d[j_n2] = max(total_n2 * dense_share, 1e-12)
        cell.N_b[j_n2] = max(total_n2 * bubble_share, 1e-12)

        cell._thermo_cache_valid = False
        profiled += 1.0
        for j in range(i + 1, len(cells)):
            propagate_upstream(cells, cfg, j)

    return {
        "bed_axial_primary_cells_profiled": float(profiled),
        "bed_axial_primary_blend_factor": float(alpha),
    }


def apply_vorabrechnung_bed01_phase_split_after_vorab_closure(
    cells: List[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Any | None = None,
) -> dict[str, float]:
    """产物/源项闭合后对 bed0–bed1 做 major phase-split 预投影（A 档 bed0–bed1 联合）。"""
    from src.core.connectivity import (
        preproject_bed_major_gas_phase_split_to_exchange_closure,
        propagate_upstream,
    )

    if not bool(getattr(cfg, "vorab_a_tier_bed01_phase_split_after_vorab_closure_thesis", True)):
        return {"bed01_phase_split_accepted": 0.0, "bed01_phase_split_cells": 0.0}

    bed_indices = _bed_cell_indices(cells)
    max_bed_index = int(max(getattr(cfg, "vorab_a_tier_bed01_phase_split_max_bed_index_thesis", 1), 0))
    target_indices = bed_indices[: max(max_bed_index + 1, 1)]
    if not target_indices:
        return {"bed01_phase_split_accepted": 0.0, "bed01_phase_split_cells": 0.0}

    accepted = 0.0
    bed0_index = bed_indices[0] if bed_indices else 0
    for i in target_indices:
        preserve_primary = bool(i == bed0_index)
        ok = preproject_bed_major_gas_phase_split_to_exchange_closure(
            cells,
            cfg,
            cell_index=i,
            preserve_primary_inlet_totals=preserve_primary,
            apply_bc_fn=apply_bc_fn,
        )
        if ok:
            accepted += 1.0
        for j in range(i + 1, len(cells)):
            propagate_upstream(cells, cfg, j)
    return {
        "bed01_phase_split_accepted": float(accepted),
        "bed01_phase_split_cells": float(len(target_indices)),
    }


def apply_vorabrechnung_bed_cell_mapping(
    cells: List[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
) -> dict[str, float]:
    """将 Vorab 预算映射到床层 cell x0（Hamel Vorab→Zellenmodell P0 原型）。

    - 每床格 ``T`` 取自 ``VorabrechnungCellBudget.t_budget_K``
    - 反应产物（CO/CO2/H2/CH4）两相分配按本地 ``eps_b/(eps_b+eps_d_void)``
    - bed0 主气化剂（O2/H2O/N2）在映射后强制对齐 Vorab 入口 dense/bubble 份额
      （``align_bottom_primary_gas_state_to_inlet_split``，Kap.3.7–3.10）

    Source: Hamel (1999) Kap.2 p.16 Vorabrechnung abbilden auf Zellenmodell
    """
    from src.core.connectivity import align_bottom_primary_gas_state_to_inlet_split

    if not cells or not cell_budgets:
        return {"bed_cells_mapped": 0.0}

    idx = GAS_SPECIES_INDEX
    product_species = ("CO", "CO2", "H2", "CH4")
    mapped = 0

    for i, cell in enumerate(cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        if i >= len(cell_budgets):
            continue
        budget = cell_budgets[i]
        cell.T = float(budget.t_budget_K)
        cell.calc_hydrodynamics()
        bubble_share = _hydrodynamic_bubble_share_for_cell(cell)
        dense_share = float(1.0 - bubble_share)

        if i == 0:
            for sp in product_species:
                j = int(idx[sp])
                total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                if total <= 1e-20:
                    continue
                cell.N_d[j] = max(total * dense_share, 1e-12)
                cell.N_b[j] = max(total * bubble_share, 1e-12)
        else:
            major_and_products = ("O2", "H2O", "N2", *product_species)
            for sp in major_and_products:
                j = int(idx[sp])
                total = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
                if total <= 1e-20:
                    continue
                cell.N_d[j] = max(total * dense_share, 1e-12)
                cell.N_b[j] = max(total * bubble_share, 1e-12)

        cell._thermo_cache_valid = False
        cell._hydro_cache_valid = False
        mapped += 1

    if mapped > 0:
        align_bottom_primary_gas_state_to_inlet_split(cells, cfg)

    return {"bed_cells_mapped": float(mapped)}


