"""Init 阶段 Vorabrechnung ↔ Zellenmodell 对齐（Hamel Check2 / Node-B 桥接）。

对应流程图 Node B：PRECALC 结果映射到 cell x0 后，在进入首次 NR 前评估
预算链与 cell 状态是否一致；未对齐时执行窄域 corrective bridge（非完整
回退 Node A）。

Ref: Hamel (1999) Bild 2.2 Check2; docs/hamel_engineering_bridge_notes.md §7–8
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

from src.core.cell import Cell, S_MOISTURE, S_VM
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.vorabrechnung import (
    MAJOR_PRODUCT_SEED_SPECIES,
    VorabrechnungCellBudget,
    _VORAB_MAJOR_PRODUCT_SOURCE_FLOOR_MOL_S,
    _bed_cell_indices,
    _expected_vm_holdup_kg_vm_zone,
    _fuel_vm_mass_kg_from_bed0,
    _major_product_init_gas_source_mol_s,
    _resolve_product_bed_cell_indices,
    _resolve_product_holdup_bed_cell_limit,
    align_vm_devolatilization_zone_vorab_state,
    apply_vorabrechnung_bed_axial_primary_profile_from_budget,
    apply_vorabrechnung_bed_solid_reactive_profile_from_budget,
    apply_vorabrechnung_bed_stream_mapping,
    seed_major_product_holdup_from_vorab_gas_balance,
)


def _solid_reactive_budget_rel_max(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
) -> float:
    """床层 VM/水分 holdup 与 Vorab 干燥/热解预算的最大相对偏差。"""
    rel_max = 0.0
    bed_indices = _bed_cell_indices(cells)
    n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    total_vm_kg = _fuel_vm_mass_kg_from_bed0(cells) if n_zone > 0 else 0.0
    for rank, i in enumerate(bed_indices):
        if i >= len(cell_budgets):
            continue
        cell = cells[i]
        if str(getattr(cell, "solid_state_model", "legacy_stream")) not in {
            "holdup_transport",
            "freeboard_closure",
        }:
            continue
        budget = cell_budgets[i]
        active_seed = np.maximum(cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in, 0.0)
        vm_seed = float(np.sum(active_seed[:, S_VM]))
        moist_seed = float(np.sum(active_seed[:, S_MOISTURE]))
        if vm_seed + moist_seed <= 1.0e-20 and not (n_zone > 0 and rank < n_zone):
            continue
        x_dry = float(np.clip(budget.x_dry_cumulative, 0.0, 1.0))
        x_vm = float(np.clip(budget.x_vm_cumulative, 0.0, 1.0))
        if n_zone > 0 and rank < n_zone and total_vm_kg > 1.0e-20:
            prev_x_vm = (
                0.0
                if rank <= 0
                else float(np.clip(cell_budgets[bed_indices[rank - 1]].x_vm_cumulative, 0.0, 1.0))
            )
            exp_vm = _expected_vm_holdup_kg_vm_zone(
                cell_rank=rank,
                total_vm_kg=total_vm_kg,
                prev_x_vm_cumulative=prev_x_vm,
                x_vm_cumulative=x_vm,
            )
            ref = max(exp_vm, total_vm_kg, 1.0e-12)
        else:
            exp_vm = max(vm_seed * (1.0 - x_vm), 0.0)
            ref = max(exp_vm + max(moist_seed * (1.0 - x_dry), 0.0), vm_seed + moist_seed, 1.0e-12)
        exp_moist = max(moist_seed * (1.0 - x_dry), 0.0)
        act_vm = float(np.sum(cell.m_solid[:, S_VM]))
        act_moist = float(np.sum(cell.m_solid[:, S_MOISTURE]))
        if not (n_zone > 0 and rank < n_zone):
            ref = max(exp_vm + exp_moist, vm_seed + moist_seed, 1.0e-12)
        rel_max = max(
            rel_max,
            abs(act_vm - exp_vm) / ref,
            abs(act_moist - exp_moist) / ref,
        )
    return float(rel_max)


def assess_vorabrechnung_cell_alignment(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None] | None = None,
) -> dict[str, float | bool]:
    """评估 A 档 Vorab 预算与当前 cell 状态的对齐度（init Check2-lite）。"""
    if apply_bc_fn is not None:
        apply_bc_fn()

    idx = GAS_SPECIES_INDEX
    bed_indices = _bed_cell_indices(cells)
    ref_primary = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_o2 = max(float(cfg.O2_feed), 1.0)

    primary_rel_max = 0.0
    if bed_indices:
        c0 = cells[bed_indices[0]]
        for sp in ("O2", "H2O", "N2"):
            j = int(idx[sp])
            hold = float(max(c0.N_d[j] + c0.N_b[j], 0.0))
            zu = float(max(c0.N_zu_d[j] + c0.N_zu_b[j], 0.0))
            if zu > 1.0e-12:
                primary_rel_max = max(primary_rel_max, abs(hold - zu) / zu)

    temperature_max_delta_K = 0.0
    o2_budget_rel_max = 0.0
    for i, budget in enumerate(cell_budgets):
        if i >= len(cells):
            break
        if i in bed_indices:
            temperature_max_delta_K = max(
                temperature_max_delta_K,
                abs(float(cells[i].T) - float(budget.t_budget_K)),
            )
            j = int(idx["O2"])
            o2_hold = float(max(cells[i].N_d[j] + cells[i].N_b[j], 0.0))
            o2_budget = float(max(budget.o2_remaining_mol_s, 0.0))
            o2_budget_rel_max = max(o2_budget_rel_max, abs(o2_hold - o2_budget) / ref_o2)

    product_rel_max = 0.0
    n_vm_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    vm_zone_bed_indices = set(_bed_cell_indices(cells)[: max(n_vm_zone, 0)])
    product_bed_limit = _resolve_product_holdup_bed_cell_limit(
        cfg,
        int(getattr(cfg, "vorab_a_tier_product_holdup_source_closure_bed_cells_thesis", 0)),
    )
    product_cells = _resolve_product_bed_cell_indices(cells, product_bed_limit)
    for i in product_cells:
        cell = cells[i]
        cell.residuals()
        in_vm_zone = bool(n_vm_zone > 0 and i in vm_zone_bed_indices)
        for sp in MAJOR_PRODUCT_SEED_SPECIES:
            j = int(idx[sp])
            hold = float(max(cell.N_d[j] + cell.N_b[j], 0.0))
            inflow = float(
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
            if in_vm_zone:
                # VM 分区格必须用固定 Vorab 源；cache 未就绪时跳过，避免 R_gas≈0 使 rel 爆炸。
                if not bool(getattr(cell, "_vm_cache_valid", False)):
                    continue
                gas_src = getattr(cell, "_vm_gas_source_cache", None)
                vorab_only = float(
                    max(float(np.asarray(gas_src, dtype=np.float64)[j]), 0.0)
                ) if gas_src is not None else 0.0
                if vorab_only <= _VORAB_MAJOR_PRODUCT_SOURCE_FLOOR_MOL_S:
                    continue
                source = vorab_only
            else:
                source = _major_product_init_gas_source_mol_s(cell, j)
            target = max(inflow + source, 1.0e-12)
            product_rel_max = max(product_rel_max, abs(hold - target) / target)

    tol_primary = float(max(getattr(cfg, "vorab_a_tier_init_abgleich_primary_rel_tol_thesis", 0.05), 0.0))
    tol_T = float(max(getattr(cfg, "vorab_a_tier_init_abgleich_temperature_K_thesis", 30.0), 0.0))
    tol_product = float(max(getattr(cfg, "vorab_a_tier_init_abgleich_product_rel_tol_thesis", 0.40), 0.0))
    tol_o2_profile = float(max(getattr(cfg, "vorab_a_tier_init_abgleich_o2_profile_rel_tol_thesis", 0.85), 0.0))
    tol_solid = float(max(getattr(cfg, "vorab_a_tier_init_abgleich_solid_rel_tol_thesis", 0.25), 0.0))

    solid_reactive_rel_max = _solid_reactive_budget_rel_max(cells, cell_budgets, cfg)

    aligned = (
        primary_rel_max <= tol_primary
        and temperature_max_delta_K <= tol_T
        and product_rel_max <= tol_product
        and o2_budget_rel_max <= tol_o2_profile
        and solid_reactive_rel_max <= tol_solid
    )

    return {
        "aligned": bool(aligned),
        "bed0_primary_holdup_rel_max": float(primary_rel_max),
        "temperature_max_delta_K": float(temperature_max_delta_K),
        "product_source_holdup_rel_max": float(product_rel_max),
        "o2_budget_profile_rel_max": float(o2_budget_rel_max),
        "solid_reactive_budget_rel_max": float(solid_reactive_rel_max),
        "tol_primary_rel": tol_primary,
        "tol_temperature_K": tol_T,
        "tol_product_rel": tol_product,
        "tol_o2_profile_rel": tol_o2_profile,
        "tol_solid_reactive_rel": tol_solid,
    }


def _product_holdup_seed_bed_cell_limit(cfg: Any) -> int:
    return _resolve_product_holdup_bed_cell_limit(
        cfg,
        int(getattr(cfg, "vorab_a_tier_product_holdup_source_closure_bed_cells_thesis", 2)),
    )


def _init_allows_bottom_gas_preprojection(cfg: Any) -> bool:
    """Init Check2 桥接是否允许 bed0 gas LSQ（默认关，outer pre-inner 仍可用 thesis 桥接）。"""
    return bool(
        getattr(cfg, "nr_init_bottom_primary_exchange_preprojection_thesis", False)
        or getattr(cfg, "nr_init_bottom_major_local_balance_preprojection_thesis", False)
        or getattr(cfg, "nr_init_bottom_major_phase_split_preprojection_thesis", False)
    )


def _run_vorab_cell_abgleich_bridge_pass(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    reseed_holdup_chain_fn: Callable[[], None] | None = None,
    allow_pyrolysis_source_refresh: bool = True,
    allow_bottom_gas_preprojection: bool = True,
) -> None:
    """单遍 Vorab ↔ cell 窄域修正（bed0 流映射、产物/固相预算、bed0 major 预投影）。"""
    from src.core.connectivity import (
        preproject_bottom_major_gas_state_to_local_balance,
        preproject_bottom_primary_gas_state_to_exchange_closure,
    )
    from src.solvers.vorabrechnung.vorab_policy import vorab_outer_may_refresh_drying_pyro_sources

    if not allow_pyrolysis_source_refresh:
        refresh_pyro_sources = False
    else:
        refresh_pyro_sources = vorab_outer_may_refresh_drying_pyro_sources(cfg)

    n_vm_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
    product_closure = bool(getattr(cfg, "vorab_a_tier_product_holdup_source_closure_thesis", True))
    # VM 分区：预投影会冲掉 bed0 CH4 等 holdup，产物 seed 放到 pass 末段。
    product_seed_before_preproj = product_closure and n_vm_zone <= 0

    apply_vorabrechnung_bed_stream_mapping(
        cells,
        cell_budgets,
        cfg,
        upstream_closure=False,
    )
    if bool(getattr(cfg, "vorab_a_tier_bed_axial_primary_profile_thesis", False)):
        apply_vorabrechnung_bed_axial_primary_profile_from_budget(
            cells,
            cell_budgets,
            cfg,
        )
    if bool(getattr(cfg, "vorab_a_tier_bed_solid_reactive_profile_thesis", True)):
        n_zone = int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0))
        if n_zone > 0:
            align_vm_devolatilization_zone_vorab_state(
                cells,
                cell_budgets,
                cfg,
                apply_bc_fn=apply_bc_fn,
                reseed_holdup_chain_fn=reseed_holdup_chain_fn,
                refresh_zone_pyrolysis_sources=refresh_pyro_sources,
            )
        else:
            apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
                cells,
                cell_budgets,
                cfg,
                apply_bc_fn=apply_bc_fn,
                reseed_holdup_chain_fn=reseed_holdup_chain_fn,
            )
    if product_seed_before_preproj:
        seed_major_product_holdup_from_vorab_gas_balance(
            cells,
            cfg,
            bed_cell_limit=_product_holdup_seed_bed_cell_limit(cfg),
            max_passes=2,
            apply_bc_fn=apply_bc_fn,
        )
    apply_bc_fn()
    if allow_bottom_gas_preprojection:
        preproject_bottom_primary_gas_state_to_exchange_closure(cells, cfg)
        preproject_bottom_major_gas_state_to_local_balance(cells, cfg)
    from src.solvers.vorabrechnung import apply_vorabrechnung_bed01_phase_split_after_vorab_closure

    if bool(getattr(cfg, "vorab_a_tier_bed01_phase_split_after_vorab_closure_thesis", False)):
        apply_vorabrechnung_bed01_phase_split_after_vorab_closure(
            cells,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )
    if product_closure and n_vm_zone > 0:
        apply_bc_fn()
        seed_major_product_holdup_from_vorab_gas_balance(
            cells,
            cfg,
            bed_cell_limit=_product_holdup_seed_bed_cell_limit(cfg),
            max_passes=4,
            blend_factor=0.55,
            apply_bc_fn=apply_bc_fn,
        )
        apply_bc_fn()
    if reseed_holdup_chain_fn is not None:
        reseed_holdup_chain_fn()
        apply_bc_fn()


def run_vorab_cell_abgleich_bridge(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    reseed_holdup_chain_fn: Callable[[], None] | None = None,
    max_passes: int | None = None,
    temperature_blend: bool = True,
    allow_pyrolysis_source_refresh: bool | None = None,
    allow_bottom_gas_preprojection: bool = True,
) -> dict[str, Any]:
    """Vorabrechnung ↔ Zellenmodell Check2 窄域桥接（init / outer 共用）。

    Hamel Bild 2.2：Check2 未对齐时的 corrective pass，非完整 Node-A PRECALC 回退。
    """
    from src.solvers.vorabrechnung.vorab_policy import vorab_outer_may_refresh_drying_pyro_sources

    if allow_pyrolysis_source_refresh is None:
        allow_pyrolysis_source_refresh = vorab_outer_may_refresh_drying_pyro_sources(cfg)
    if max_passes is None:
        max_passes = int(max(getattr(cfg, "vorab_a_tier_init_abgleich_max_passes_thesis", 2), 0))
        if int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)) > 0:
            max_passes = max(max_passes, 3)

    history: list[dict[str, float | bool]] = []
    metrics = assess_vorabrechnung_cell_alignment(
        cells,
        cell_budgets,
        cfg,
        apply_bc_fn=apply_bc_fn,
    )
    history.append(dict(metrics))
    bridge_passes = 0

    if bool(metrics["aligned"]) or max_passes <= 0:
        return {
            "aligned": bool(metrics["aligned"]),
            "bridge_passes": float(bridge_passes),
            "history": history,
            "final": metrics,
        }

    blend_T = (
        float(np.clip(getattr(cfg, "vorab_a_tier_init_abgleich_temperature_blend_thesis", 0.25), 0.0, 1.0))
        if temperature_blend
        else 0.0
    )
    if bool(getattr(cfg, "vorab_bed_temperature_march_thesis", False)):
        blend_T = float(
            np.clip(getattr(cfg, "vorab_a_tier_init_abgleich_temperature_blend_march_thesis", 0.0), 0.0, 1.0)
        )
    bed_indices = _bed_cell_indices(cells)

    for _ in range(max_passes):
        if bool(metrics["aligned"]):
            break
        bridge_passes += 1

        _run_vorab_cell_abgleich_bridge_pass(
            cells,
            cell_budgets,
            cfg,
            apply_bc_fn=apply_bc_fn,
            reseed_holdup_chain_fn=reseed_holdup_chain_fn,
            allow_pyrolysis_source_refresh=allow_pyrolysis_source_refresh,
            allow_bottom_gas_preprojection=allow_bottom_gas_preprojection,
        )
        if blend_T > 0.0 and bed_indices:
            for i in bed_indices:
                if i >= len(cell_budgets):
                    continue
                t_ref = float(cell_budgets[i].t_budget_K)
                cells[i].T = float((1.0 - blend_T) * float(cells[i].T) + blend_T * t_ref)
                cells[i]._thermo_cache_valid = False
                cells[i]._hydro_cache_valid = False
        apply_bc_fn()
        metrics = assess_vorabrechnung_cell_alignment(
            cells,
            cell_budgets,
            cfg,
            apply_bc_fn=apply_bc_fn,
        )
        history.append(dict(metrics))

    return {
        "aligned": bool(metrics["aligned"]),
        "bridge_passes": float(bridge_passes),
        "history": history,
        "final": metrics,
    }


def run_init_vorab_cell_abgleich_bridge(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    max_passes: int | None = None,
    reseed_holdup_chain_fn: Callable[[], None] | None = None,
) -> dict[str, Any]:
    """Init 末段 Check2-lite：未对齐则重跑 Vorab ↔ cell 窄域桥接。"""
    if max_passes is None:
        max_passes = int(max(getattr(cfg, "vorab_a_tier_init_abgleich_max_passes_thesis", 2), 0))
        if not _init_allows_bottom_gas_preprojection(cfg):
            max_passes = 0
    return run_vorab_cell_abgleich_bridge(
        cells,
        cell_budgets,
        cfg,
        apply_bc_fn=apply_bc_fn,
        reseed_holdup_chain_fn=reseed_holdup_chain_fn,
        max_passes=max_passes,
        temperature_blend=True,
        allow_bottom_gas_preprojection=_init_allows_bottom_gas_preprojection(cfg),
    )


def run_outer_vorab_cell_abgleich_bridge(
    cells: list[Cell],
    cell_budgets: list[VorabrechnungCellBudget],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None],
    reseed_holdup_chain_fn: Callable[[], None] | None = None,
) -> dict[str, Any]:
    """Outer Check2 失败后的单遍 Vorab ↔ cell 修正（Hamel Node-B 轻量回退）。"""
    if not bool(getattr(cfg, "vorab_a_tier_outer_abgleich_bridge_thesis", True)):
        return {"enabled": False, "aligned": False, "bridge_passes": 0.0, "history": [], "final": {}}
    max_passes = int(max(getattr(cfg, "vorab_a_tier_outer_abgleich_max_passes_thesis", 1), 0))
    result = run_vorab_cell_abgleich_bridge(
        cells,
        cell_budgets,
        cfg,
        apply_bc_fn=apply_bc_fn,
        reseed_holdup_chain_fn=reseed_holdup_chain_fn,
        max_passes=max_passes,
        temperature_blend=False,
    )
    result["enabled"] = True
    return result


def measure_outer_abgleich_gas_metrics(
    cells: list[Cell],
    cfg: Any,
    *,
    apply_bc_fn: Callable[[], None] | None = None,
) -> dict[str, float]:
    """Outer Check2 气相扩展：bed major phase-split RMS 与 bed0 主气 zu 漂移。"""
    if apply_bc_fn is not None:
        apply_bc_fn()

    from src.solvers.global_nr_solver import (
        _residual_gas_phase_metrics,
        build_equation_scales,
        global_residual,
        pack_reactor,
    )

    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)

    def _bc() -> None:
        if apply_bc_fn is not None:
            apply_bc_fn()

    f_raw = np.asarray(global_residual(pack_reactor(cells), cells, _bc), dtype=np.float64)
    scale = np.asarray(
        build_equation_scales(cells, ref_gas, ref_solid, ref_energy),
        dtype=np.float64,
    )
    f_hat = f_raw / scale
    gas_phase = _residual_gas_phase_metrics(f_hat, cells)

    bed0_primary_rel = 0.0
    bed_indices = _bed_cell_indices(cells)
    if bed_indices:
        c0 = cells[bed_indices[0]]
        idx = GAS_SPECIES_INDEX
        for sp in ("O2", "H2O", "N2"):
            j = int(idx[sp])
            hold = float(max(c0.N_d[j] + c0.N_b[j], 0.0))
            zu = float(max(c0.N_zu_d[j] + c0.N_zu_b[j], 0.0))
            if zu > 1.0e-12:
                bed0_primary_rel = max(bed0_primary_rel, abs(hold - zu) / zu)

    return {
        "gas_phase_split_rms": float(gas_phase["gas_phase_split_rms"]),
        "gas_phase_split_max_abs": float(gas_phase["gas_phase_split_max_abs"]),
        "bed0_primary_rel_drift": float(bed0_primary_rel),
    }


def outer_convergence_from_config(cfg: Any) -> "OuterConvergence":
    """Build Check2 criterion from reactor config (thesis A-tier gates)."""
    from src.solvers.convergence import OuterConvergence

    return OuterConvergence(
        rel_tol=float(getattr(cfg, "nr_outer_check2_rel_tol", 1e-4)),
        abs_tol_T=float(getattr(cfg, "nr_outer_check2_abs_tol_T_K", 1.0)),
        enable_gas_phase_split_gate=bool(
            getattr(cfg, "nr_outer_check2_gas_phase_split_gate_thesis", True)
        ),
        max_gas_phase_split_rms=float(
            getattr(cfg, "nr_outer_check2_max_gas_phase_split_rms_thesis", 0.05)
        ),
        enable_bed0_primary_gate=bool(
            getattr(cfg, "nr_outer_check2_bed0_primary_gate_thesis", False)
        ),
        max_bed0_primary_rel_drift=float(
            getattr(cfg, "nr_outer_check2_max_bed0_primary_rel_drift_thesis", 0.10)
        ),
    )
