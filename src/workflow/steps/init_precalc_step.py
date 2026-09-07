"""流程图：赋初值 + Vorabrechnung 快照（轴向 T、x0、水动力冻结）。

对应 mmd：InitVal → PRECALC（在进入 NR/Abgleich 前完成初始化段）。
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from time import perf_counter
from typing import TYPE_CHECKING, Any

import numpy as np

from src.core.constants import Rg
from src.core.species import GAS_SPECIES_INDEX

if TYPE_CHECKING:
    from src.core.reactor import Reactor


def _reconcile_upper_bed_solid_holdup_if_enabled(
    reactor: Reactor,
    *,
    allow_reduce: bool | None = None,
) -> None:
    """bed1+ char/ash holdup 对齐当前 ``K_auf`` 输运支撑。"""
    if reactor._reconcile_upper_bed_solid_holdup_for_nr(allow_reduce=allow_reduce):
        reactor._apply_all_bc_for_nr()


def _snap_bed_holdup_to_local_support(
    reactor: Reactor,
    *,
    bed_indices: tuple[int, ...] = (1,),
    max_passes: int = 3,
) -> None:
    """对指定床格把 char/ash 贴到当前 Eq.2.6 support（每次后刷新 BC）。

    全床同时 ``allow_reduce`` 会对抢 auf/ab 入流；能量主导的 bed1 须单独收敛。
    Ref: Hamel (1999) Eq.2.6–2.7; handoff §5.55.
    """
    if not bool(getattr(reactor.config, "nr_init_align_upper_bed_solid_holdup_thesis", True)):
        return
    from src.core.connectivity.solid_transport import align_bed_solid_holdup_to_transport_inflow_support

    cells = reactor.cells
    for _ in range(max(int(max_passes), 1)):
        changed = False
        for i in bed_indices:
            if i < 0 or i >= len(cells):
                continue
            cell = cells[i]
            if str(getattr(cell, "cell_type", "bed")) != "bed":
                continue
            changed |= bool(
                align_bed_solid_holdup_to_transport_inflow_support(cell, allow_reduce=True)
            )
        if not changed:
            break
        reactor._apply_all_bc_for_nr()
    # 末次 BC 会微移 support；再贴一次且不再刷新，让 NR 看到的 hold≈support。
    for i in bed_indices:
        if 0 <= i < len(cells) and str(getattr(cells[i], "cell_type", "bed")) == "bed":
            align_bed_solid_holdup_to_transport_inflow_support(cells[i], allow_reduce=True)


def _o2_syngas_overlap_mol2_s2(cell: Any) -> float:
    """O₂ 与 H₂/CO/CH₄ holdup 乘积 [ (mol/s)² ]：R12 刚性格点探针。"""
    idx = GAS_SPECIES_INDEX
    n = np.asarray(cell.N_d, dtype=np.float64) + np.asarray(cell.N_b, dtype=np.float64)
    o2 = float(max(n[idx["O2"]], 0.0))
    fuel = float(max(n[idx["H2"]], 0.0) + max(n[idx["CO"]], 0.0) + max(n[idx["CH4"]], 0.0))
    return o2 * fuel


def _run_bottom_cell_startwert_solve(reactor: Reactor) -> None:
    """快氧化后对 bed0 做单格守恒求解；变差或冷根则回滚。

    Ref: Hamel (1999) §2.1 Startwertwahl（生成初值）；冷根见 handoff §5.30。
    """
    cfg = reactor.config
    if not bool(getattr(cfg, "vorab_init_solve_bottom_cell_thesis", False)):
        reactor._vorab_init_bottom_cell_solve_diag = {"enabled": 0.0, "accepted": 0.0}
        return
    cells = getattr(reactor, "cells", None) or []
    if not cells:
        reactor._vorab_init_bottom_cell_solve_diag = {"enabled": 1.0, "accepted": 0.0, "empty": 1.0}
        return
    cell = cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        reactor._vorab_init_bottom_cell_solve_diag = {"enabled": 1.0, "accepted": 0.0, "not_bed": 1.0}
        return

    from src.solvers.cell_solver import _pack_state, _unpack_state, evaluate_cell_state, solve_cell

    T_floor = float(getattr(cfg, "vorab_init_solve_bottom_cell_T_floor_K_thesis", 800.0))
    snap = _pack_state(cell).copy()
    before = evaluate_cell_state(cell)
    overlap_before = _o2_syngas_overlap_mol2_s2(cell)
    info = solve_cell(
        cell,
        max_iter=12,
        stiff_stabilization=True,
        skip_homotopy=True,
        skip_multistart=True,
        temperature_floor_K=T_floor,
        verbose=False,
    )
    after = evaluate_cell_state(cell)
    overlap_after = _o2_syngas_overlap_mol2_s2(cell)
    T_after = float(cell.T)
    rms_before = float(before["rms_scaled"])
    rms_after = float(after["rms_scaled"])
    accept = (
        np.isfinite(rms_after)
        and rms_after <= 0.98 * max(rms_before, 1e-16)
        and T_after >= T_floor - 1e-6
        and overlap_after <= max(overlap_before, 1e-12) * 1.05
    )
    if not accept:
        _unpack_state(snap, cell)
        cell._thermo_cache_valid = False
    reactor._vorab_init_bottom_cell_solve_diag = {
        "enabled": 1.0,
        "accepted": float(accept),
        "rms_before": rms_before,
        "rms_after": rms_after if accept else rms_before,
        "rms_attempt": rms_after,
        "T_after": float(cell.T),
        "T_attempt": T_after,
        "overlap_before": overlap_before,
        "overlap_after": overlap_after if accept else overlap_before,
        "overlap_attempt": overlap_after,
        "solver_converged": float(bool(info.get("converged"))),
        "solver_rms": float(info.get("rms_scaled", np.nan)),
    }


@dataclass(frozen=True)
class GlobalNRInitPrecalcResult:
    """INIT + 预计算阶段计时与中间量（供结果汇总与审计）。"""

    resolved_init_strategy: str
    nr_init_s_total: float  # [s]
    nr_vorabrechnung_s: float  # [s]
    T_est: np.ndarray  # [K]


@dataclass
class _InitPrecalcState:
    """五阶段 init 编排共享状态（顺序敏感，勿跳步）。"""

    resolved: str
    T_est: np.ndarray  # [K]
    cell_budgets: Any | None
    use_a_tier_budget: bool
    use_cell_mapping: bool
    use_bed_stream_mapping: bool
    use_march: bool
    t_floor_K: float  # [K]
    t_cap_K: float  # [K]
    init_started: float  # [s]
    vorab_started: float  # [s]


def _prepare_macro_hydrodynamics_seed(
    reactor: Reactor,
    T_est: np.ndarray,
) -> list[float | None]:
    """Seed Vorabrechnung hydrodynamics from macro gas throughput, not cell inventories.

    Hamel's Vorabrechnung requires hydrodynamics to be evaluated before the
    chemistry-consistent x0 exists. To keep that ordering while avoiding
    dependence on stale ``N_d/N_b`` state, we temporarily:
    - override ``u0_target`` from the reactor-level primary gas feed, and
    - clear the bed-cell gas inventories that would otherwise back-drive
      ``calc_cell_hydrodynamics``.
    """
    cfg = reactor.config
    original_u0_targets = [cell.u0_target for cell in reactor.cells]
    if cfg.u0_target is not None:
        for cell in reactor.cells:
            cell.N_d.fill(0.0)
            cell.N_b.fill(0.0)
        return original_u0_targets

    area = math.pi / 4.0 * float(cfg.D_bed) ** 2
    primary_mol_s = float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed)
    for i, cell in enumerate(reactor.cells):
        cell.N_d.fill(0.0)
        cell.N_b.fill(0.0)
        cell.u0_target = (
            (primary_mol_s * Rg * float(T_est[i]) / float(cfg.P)) / max(area, 1e-12)
            if primary_mol_s > 0.0
            else 0.0
        )
    return original_u0_targets


def _restore_macro_hydrodynamics_seed(
    reactor: Reactor,
    original_u0_targets: list[float | None],
) -> None:
    """Restore user/configured hydrodynamics targets after Vorabrechnung x0 generation."""
    for cell, original_u0 in zip(reactor.cells, original_u0_targets):
        cell.u0_target = original_u0


def _vorab_budget_gas_mol_s_at_cell(
    cfg,
    budget: Any | None,
) -> float:
    """Vorab 预算在该格高度上的总气量 [mol/s]（主气 + 热解增量代理）。"""
    base_gas_mol_s = float(max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 0.0))
    if budget is None:
        return base_gas_mol_s
    pyro_gas_mol_s = float(
        max(float(budget.nC_vm_mol_s) + 0.5 * float(budget.nH_vm_mol_s), 0.0)
    )
    return float(max(base_gas_mol_s + pyro_gas_mol_s, 0.0))


def _apply_init_bed_vorab_gas_u0_targets(
    reactor: Reactor,
    *,
    T_est: np.ndarray,
    cell_budgets: Any | None,
) -> list[float | None]:
    """Init 阶段 bed 水力学：用 Vorab 气量 throughput 固定 ``u0_target`` 并刷新 ``eps_b/K``。

    Transport x0 的 ``N_d/N_b`` 在 bed0 常低于 Vorab 预算的总气量（尤其 O₂/Steam）；
    直接 ``calc_hydrodynamics`` 会得到极低 ``eps_b`` 与 ``K_auf``，进而触发
    ``inflow/K`` holdup 爆炸。Hamel Vorabrechnung 在 chemistry-consistent x0 之前
    先用宏观气量算水力学（``_prepare_macro_hydrodynamics_seed``）；此处将同一语义
    延伸到 init BC/K 刷新，并在 NR 开始前恢复 ``u0_target=None``。

    Ref: Hamel (1999) Kapitel 2.1 Vorabrechnung; Eq.3.24 visible bubble fraction
    """
    cfg = reactor.config
    saved = [cell.u0_target for cell in reactor.cells]
    if not bool(getattr(cfg, "vorab_init_bed_vorab_gas_u0_target_thesis", True)):
        return saved
    if cfg.u0_target is not None:
        return saved

    area = math.pi / 4.0 * float(cfg.D_bed) ** 2
    for i, cell in enumerate(reactor.cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        budget = None
        if cell_budgets is not None and i < len(cell_budgets):
            budget = cell_budgets[i]
        gas_mol_s = _vorab_budget_gas_mol_s_at_cell(cfg, budget)
        t_i = float(T_est[i] if i < len(T_est) else cell.T)
        cell.u0_target = (
            (gas_mol_s * Rg * t_i / float(cfg.P)) / max(area, 1e-12)
            if gas_mol_s > 0.0
            else 0.0
        )
        cell._thermo_cache_valid = False
        cell.calc_hydrodynamics()
    return saved


def _vorab_bottom_temperature_bounds_K(cfg) -> tuple[float, float]:
    """Vorab 绝热 clip 床底温度上下界 [K]。"""
    floor_K = float(getattr(cfg, "vorab_bottom_temperature_floor_K_thesis", 1050.0))
    cap_K = float(getattr(cfg, "vorab_bottom_temperature_cap_K_thesis", 1325.0))
    floor_K = float(max(floor_K, 300.0))
    cap_K = float(max(cap_K, floor_K + 1.0))
    return floor_K, cap_K


def _install_temperature_fence_from_vorabrechnung(reactor: Reactor, T_est: np.ndarray) -> None:
    """按 Vorabrechnung T 估计给每格挂 NR 温度围栏。

    底格可用 ``nr_temperature_fence_bed0_lower_margin_K_thesis`` 单独加宽下沿
    （handoff §5.58）；上段可用 ``nr_temperature_fence_upper_bed_lower_margin_K_thesis``
    （handoff §5.64，CORE 仍 None）。未覆盖的格仍用全局 ``lower_margin``。
    """
    cfg = reactor.config
    enabled = bool(getattr(cfg, "nr_temperature_fence_enabled_thesis", True))
    use_march = bool(getattr(cfg, "vorab_bed_temperature_march_thesis", False))
    lower_margin = float(max(getattr(cfg, "nr_temperature_fence_lower_margin_K_thesis", 150.0), 0.0))
    if use_march:
        lower_margin = float(
            min(
                lower_margin,
                float(getattr(cfg, "vorab_nr_temperature_fence_lower_margin_K_thesis", 80.0)),
            )
        )
    upper_margin = float(max(getattr(cfg, "nr_temperature_fence_upper_margin_K_thesis", 350.0), 0.0))
    march_upper_anchor = float(
        getattr(cfg, "vorab_nr_temperature_fence_upper_anchor_K_thesis", 1400.0)
    )
    for i, cell in enumerate(reactor.cells):
        if i >= len(T_est) or not enabled:
            for attr in ("_vorab_temperature_estimate_K", "_nr_temperature_min_K", "_nr_temperature_max_K"):
                if hasattr(cell, attr):
                    delattr(cell, attr)
            continue
        t_ref = float(T_est[i])
        cell._vorab_temperature_estimate_K = t_ref
        abs_min_K = float(max(getattr(cfg, "nr_temperature_fence_absolute_min_K_thesis", 300.0), 300.0))
        cell_lower_margin = lower_margin
        if i == 0:
            bed0_margin = getattr(cfg, "nr_temperature_fence_bed0_lower_margin_K_thesis", None)
            if bed0_margin is not None:
                cell_lower_margin = float(max(float(bed0_margin), 0.0))
        else:
            upper_bed_margin = getattr(
                cfg, "nr_temperature_fence_upper_bed_lower_margin_K_thesis", None
            )
            if upper_bed_margin is not None:
                cell_lower_margin = float(max(float(upper_bed_margin), 0.0))
        rel_min_K = float(t_ref - cell_lower_margin)
        lower_mode = str(getattr(cfg, "nr_temperature_fence_lower_mode_thesis", "relative"))
        if lower_mode == "absolute":
            cell._nr_temperature_min_K = abs_min_K
        else:
            cell._nr_temperature_min_K = float(max(abs_min_K, rel_min_K))
        if use_march:
            cell._nr_temperature_max_K = float(
                min(2500.0, max(t_ref + upper_margin, march_upper_anchor))
            )
        else:
            cell._nr_temperature_max_K = float(min(2500.0, t_ref + upper_margin))


def _apply_bed0_nr_startwert_T(reactor: Reactor) -> None:
    """底格 NR Startwert T：不改围栏中心 T_est，只改 cell.T。

    须在 generate_initial_x0 / 末次 snap 之后调用（那些步骤会写回 Vorab T）。
    Ref: Hamel (1999) Kap.2 Startwertwahl; handoff §5.59.
    """
    t0 = getattr(reactor.config, "vorab_bed0_nr_startwert_T_K_thesis", None)
    if t0 is None or not reactor.cells:
        return
    cell = reactor.cells[0]
    if str(getattr(cell, "cell_type", "bed")) != "bed":
        return
    lo = float(getattr(cell, "_nr_temperature_min_K", 300.0))
    hi = float(getattr(cell, "_nr_temperature_max_K", 2500.0))
    if not np.isfinite(lo):
        lo = 300.0
    if not np.isfinite(hi):
        hi = 2500.0
    lo = max(lo, 300.0)
    hi = min(max(hi, lo + 1.0), 2500.0)
    cell.T = float(np.clip(float(t0), lo, hi))
    cell._thermo_cache_valid = False
    if hasattr(cell, "_h_cache") and cell._h_cache is not None:
        cell._h_cache.clear()


def _resolve_init_strategy_and_seed_T(
    reactor: Reactor,
    *,
    init_strategy: str | None,
) -> _InitPrecalcState:
    """阶段 1：解析 init 策略、轴向 T 估计、温度 fence、进料/上游传播。"""
    from src.core.reactor import _resolve_nr_init_strategy
    from src.solvers.vorabrechnung import estimate_axial_T_profile

    cfg = reactor.config
    resolved = _resolve_nr_init_strategy(init_strategy)
    init_started = perf_counter()
    vorab_started = perf_counter()
    use_a_tier_budget = bool(getattr(cfg, "vorab_a_tier_startwert_budget_thesis", True))
    use_cell_mapping = use_a_tier_budget and bool(getattr(cfg, "vorab_a_tier_cell_mapping_thesis", True))
    use_bed_stream_mapping = use_a_tier_budget and bool(
        getattr(cfg, "vorab_a_tier_bed_stream_mapping_thesis", True)
    )
    use_march = bool(getattr(cfg, "vorab_bed_temperature_march_thesis", False))
    t_floor_K, t_cap_K = _vorab_bottom_temperature_bounds_K(cfg)
    if use_march:
        T_est = np.full(len(reactor.cells), float(cfg.T_inlet), dtype=np.float64)
    else:
        T_est = estimate_axial_T_profile(
            n_cells=cfg.n_cells,
            T_inlet=cfg.T_inlet,
            O2_feed=cfg.O2_feed,
            H2O_feed=cfg.H2O_feed,
            N2_feed=cfg.N2_feed,
            fuel_feed_kg_s=cfg.fuel_feed,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            moisture_wt=cfg.moisture_wt,
            P=cfg.P,
            T_bottom_floor_K=t_floor_K,
            T_bottom_cap_K=t_cap_K,
        )
    for i, cell in enumerate(reactor.cells):
        cell.T = float(T_est[i])
    _install_temperature_fence_from_vorabrechnung(reactor, T_est)
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    return _InitPrecalcState(
        resolved=resolved,
        T_est=T_est,
        cell_budgets=None,
        use_a_tier_budget=use_a_tier_budget,
        use_cell_mapping=use_cell_mapping,
        use_bed_stream_mapping=use_bed_stream_mapping,
        use_march=use_march,
        t_floor_K=t_floor_K,
        t_cap_K=t_cap_K,
        init_started=init_started,
        vorab_started=vorab_started,
    )


def _run_vorabrechnung_snapshot(reactor: Reactor, state: _InitPrecalcState) -> _InitPrecalcState:
    """阶段 2：Vorabrechnung x0 / A-tier budget / stream mapping（macro hydro 临时 seed）。"""
    from src.core.connectivity import preproject_bed_major_gas_phase_split_to_exchange_closure
    from src.solvers.vorabrechnung import (
        apply_vorabrechnung_bed_stream_mapping,
        build_vorabrechnung_cell_budgets,
        generate_initial_x0,
    )

    cfg = reactor.config
    T_est = state.T_est
    cell_budgets = state.cell_budgets
    original_u0_targets = _prepare_macro_hydrodynamics_seed(reactor, T_est)
    try:
        for cell in reactor.cells:
            cell.calc_hydrodynamics()
        reactor._snapshot_bottom_gas_inlet_split_from_vorabrechnung()
        reactor._set_bottom_cell_feeds()
        for i in range(len(reactor.cells)):
            reactor._propagate_upstream(i)

        if state.use_a_tier_budget:
            from src.solvers.vorabrechnung import build_nr_temperature_seed_for_march_mode

            t_floor_K, t_cap_K = _vorab_bottom_temperature_bounds_K(cfg)
            T_marched, cell_budgets = build_vorabrechnung_cell_budgets(
                reactor.cells,
                T_inlet=cfg.T_inlet,
                O2_feed=cfg.O2_feed,
                H2O_feed=cfg.H2O_feed,
                N2_feed=cfg.N2_feed,
                fuel_feed_kg_s=cfg.fuel_feed,
                C_dry=cfg.C_dry,
                H_dry=cfg.H_dry,
                O_dry=cfg.O_dry,
                moisture_wt=cfg.moisture_wt,
                ash_dry_wt=cfg.ash_dry_wt,
                VM_daf=cfg.VM_daf,
                fuel_type=cfg.fuel_type,
                P=cfg.P,
                T_bottom_floor_K=t_floor_K,
                T_bottom_cap_K=t_cap_K,
                cfg=cfg,
            )
            seed_mode = str(getattr(cfg, "vorab_nr_temperature_seed_mode_thesis", "adiabatic_anchor"))
            if state.use_march:
                T_est = build_nr_temperature_seed_for_march_mode(
                    reactor.cells,
                    cell_budgets,
                    T_inlet=float(cfg.T_inlet),
                    O2_feed=float(cfg.O2_feed),
                    H2O_feed=float(cfg.H2O_feed),
                    N2_feed=float(cfg.N2_feed),
                    fuel_feed_kg_s=float(cfg.fuel_feed),
                    C_dry=float(cfg.C_dry),
                    H_dry=float(cfg.H_dry),
                    moisture_wt=float(cfg.moisture_wt),
                    P=float(cfg.P),
                    cfg=cfg,
                    T_bottom_cap_K=t_cap_K,
                )
                reactor._vorab_marched_t_budget_K = np.array(T_marched, dtype=np.float64)
            elif seed_mode in ("budget_progress", "blend"):
                T_est = build_nr_temperature_seed_for_march_mode(
                    reactor.cells,
                    cell_budgets,
                    T_inlet=float(cfg.T_inlet),
                    O2_feed=float(cfg.O2_feed),
                    H2O_feed=float(cfg.H2O_feed),
                    N2_feed=float(cfg.N2_feed),
                    fuel_feed_kg_s=float(cfg.fuel_feed),
                    C_dry=float(cfg.C_dry),
                    H_dry=float(cfg.H_dry),
                    moisture_wt=float(cfg.moisture_wt),
                    P=float(cfg.P),
                    cfg=cfg,
                    T_bottom_cap_K=t_cap_K,
                )
            else:
                T_est = T_marched
            for i, cell in enumerate(reactor.cells):
                cell.T = float(T_est[i])
            _install_temperature_fence_from_vorabrechnung(reactor, T_est)

        generate_initial_x0(
            cells=reactor.cells,
            O2_feed=cfg.O2_feed,
            H2O_feed=cfg.H2O_feed,
            N2_feed=cfg.N2_feed,
            fuel_feed_kg_s=cfg.fuel_feed,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            O_dry=cfg.O_dry,
            moisture_wt=cfg.moisture_wt,
            ash_dry_wt=cfg.ash_dry_wt,
            VM_daf=cfg.VM_daf,
            T_profile=T_est,
            fuel_type=cfg.fuel_type,
            use_hamel_major_gibbs_x0=bool(cfg.vorab_major_gibbs_x0),
            strict_hamel_major_gibbs_x0=bool(cfg.thesis_mode),
            major_gibbs_solver_mode=str(cfg.major_gibbs_solver_mode),
            cell_budgets=cell_budgets,
            major_gibbs_temperature_floor_K=(
                float(getattr(cfg, "vorab_major_gibbs_min_temperature_K_thesis", 750.0))
                if str(getattr(cfg, "vorab_nr_temperature_seed_mode_thesis", "adiabatic_anchor"))
                != "adiabatic_anchor"
                else 0.0
            ),
            x0_holdup_mode=str(getattr(cfg, "vorab_x0_holdup_mode_thesis", "transport_from_vorab")),
            cfg=cfg,
            use_pyrolysis_gibbs_composition=bool(
                getattr(cfg, "vorab_pyrolysis_gibbs_composition_thesis", True)
            ),
            pyrolysis_gibbs_solver_mode=str(
                getattr(cfg, "vorab_pyrolysis_gibbs_solver_mode_thesis", "hamel_reduced")
            ),
            pyrolysis_gibbs_temperature_floor_K=float(
                getattr(cfg, "vorab_major_gibbs_min_temperature_K_thesis", 750.0)
            ),
        )
        if state.use_bed_stream_mapping and cell_budgets is not None:
            reactor._vorab_a_tier_bed_stream_mapping_diag = apply_vorabrechnung_bed_stream_mapping(
                reactor.cells,
                cell_budgets,
                cfg,
                upstream_closure=bool(
                    getattr(cfg, "vorab_a_tier_bed_stream_upstream_closure_thesis", False)
                ),
            )
            if bool(getattr(cfg, "vorab_a_tier_bed1_phase_split_after_stream_mapping_thesis", False)):
                bed_indices = [
                    i
                    for i, c in enumerate(reactor.cells)
                    if str(getattr(c, "cell_type", "bed")) == "bed"
                ]
                bed1_split_ok = False
                if len(bed_indices) >= 2:
                    bed1_split_ok = preproject_bed_major_gas_phase_split_to_exchange_closure(
                        reactor.cells,
                        cfg,
                        cell_index=bed_indices[1],
                        preserve_primary_inlet_totals=True,
                    )
                    for i in range(bed_indices[1] + 1, len(reactor.cells)):
                        reactor._propagate_upstream(i)
                reactor._vorab_a_tier_bed1_phase_split_diag = {
                    "bed1_index": float(bed_indices[1]) if len(bed_indices) >= 2 else -1.0,
                    "accepted": bool(bed1_split_ok),
                }
            else:
                reactor._vorab_a_tier_bed1_phase_split_diag = None
        else:
            reactor._vorab_a_tier_bed_stream_mapping_diag = None
            reactor._vorab_a_tier_bed1_phase_split_diag = None
        if state.use_cell_mapping and cell_budgets is not None:
            from src.solvers.vorabrechnung import apply_vorabrechnung_bed_cell_mapping

            reactor._vorab_a_tier_cell_mapping_after_x0_diag = apply_vorabrechnung_bed_cell_mapping(
                reactor.cells,
                cell_budgets,
                cfg,
            )
    finally:
        _restore_macro_hydrodynamics_seed(reactor, original_u0_targets)

    state.T_est = T_est
    state.cell_budgets = cell_budgets
    return state


def _run_post_vorabrechnung_holdup_and_reactive(reactor: Reactor, state: _InitPrecalcState) -> None:
    """阶段 3：BC 刷新、holdup 链 seed、A-tier 固相 reactive profile。"""
    from src.solvers.vorabrechnung import (
        align_vm_devolatilization_zone_vorab_state,
        apply_vorabrechnung_bed_solid_reactive_profile_from_budget,
    )

    cfg = reactor.config
    cell_budgets = state.cell_budgets
    reactor._apply_all_bc_for_nr()
    reactor._seed_initialized_holdup_chain_for_nr()
    if state.use_a_tier_budget and cell_budgets is not None and bool(
        getattr(cfg, "vorab_a_tier_bed_solid_reactive_profile_thesis", True)
    ):
        if int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)) > 0:
            reactor._vorab_a_tier_bed_solid_reactive_profile_diag = align_vm_devolatilization_zone_vorab_state(
                reactor.cells,
                cell_budgets,
                cfg,
                apply_bc_fn=reactor._apply_all_bc_for_nr,
                reseed_holdup_chain_fn=reactor._seed_initialized_holdup_chain_for_nr,
                refresh_zone_pyrolysis_sources=False,
            )
        else:
            reactor._vorab_a_tier_bed_solid_reactive_profile_diag = (
                apply_vorabrechnung_bed_solid_reactive_profile_from_budget(
                    reactor.cells,
                    cell_budgets,
                    cfg,
                    apply_bc_fn=reactor._apply_all_bc_for_nr,
                    reseed_holdup_chain_fn=reactor._seed_initialized_holdup_chain_for_nr,
                )
            )
    else:
        reactor._vorab_a_tier_bed_solid_reactive_profile_diag = None
    if state.use_a_tier_budget:
        reactor._vorab_a_tier_cell_budgets = cell_budgets


def _run_bottom_zone_preprojections(reactor: Reactor, state: _InitPrecalcState) -> None:
    """阶段 4：床底/zone 预投影（primary、major、product seed、phase split、solid energy）。"""
    from src.core.connectivity import (
        preproject_bottom_major_gas_state_to_local_balance,
        preproject_bottom_primary_gas_state_to_exchange_closure,
        run_bottom_zone_solid_energy_preprojection,
    )
    from src.solvers.vorabrechnung import (
        apply_vorabrechnung_bed01_phase_split_after_vorab_closure,
        seed_major_product_holdup_from_vorab_gas_balance,
        _resolve_product_holdup_bed_cell_limit,
    )

    cfg = reactor.config
    cell_budgets = state.cell_budgets
    reactor._apply_all_bc_for_nr()
    run_primary_preproj = not (
        state.use_cell_mapping
        and not bool(getattr(cfg, "nr_bottom_primary_exchange_preprojection_when_cell_mapped_thesis", False))
    )
    if run_primary_preproj and bool(
        getattr(cfg, "nr_init_bottom_primary_exchange_preprojection_thesis", False)
    ):
        preproject_bottom_primary_gas_state_to_exchange_closure(reactor.cells, cfg)
    if bool(getattr(cfg, "nr_init_bottom_major_local_balance_preprojection_thesis", False)):
        preproject_bottom_major_gas_state_to_local_balance(reactor.cells, cfg)
    reactor._refresh_explicit_freeboard_transport_from_closure()
    reactor._initialize_explicit_side_block_states()
    # 产物 holdup seed 依赖固定 Vorab 源；VM 分区时须先 single-shot。
    reactor._initialize_thesis_single_shot_vorab_sources()
    if state.use_a_tier_budget and bool(getattr(cfg, "vorab_a_tier_product_holdup_source_closure_thesis", True)):
        product_seed_diag = seed_major_product_holdup_from_vorab_gas_balance(
            reactor.cells,
            cfg,
            bed_cell_limit=_resolve_product_holdup_bed_cell_limit(
                cfg,
                int(getattr(cfg, "vorab_a_tier_product_holdup_source_closure_bed_cells_thesis", 2)),
            ),
            max_passes=3,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
        )
        reactor._vorab_a_tier_product_holdup_seed_diag = product_seed_diag
        reactor._apply_all_bc_for_nr()
        reactor._vorab_a_tier_bed01_phase_split_diag = apply_vorabrechnung_bed01_phase_split_after_vorab_closure(
            reactor.cells,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
        )
        reactor._apply_all_bc_for_nr()
    else:
        reactor._vorab_a_tier_product_holdup_seed_diag = None
        reactor._vorab_a_tier_bed01_phase_split_diag = None
    reactor._apply_all_bc_for_nr()
    # bed0 phase-split 改在 init 末段（align/abgleich 之后）；见 _run_init_bed0_phase_split_closure
    if bool(getattr(cfg, "nr_init_bottom_zone_solid_energy_preprojection_thesis", False)):
        reactor._bottom_zone_solid_energy_preprojection_diag = run_bottom_zone_solid_energy_preprojection(
            reactor.cells,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        )
    else:
        reactor._bottom_zone_solid_energy_preprojection_diag = None
    reactor._apply_all_bc_for_nr()


def _run_bed0_eq24_two_phase_startwert(reactor: Reactor) -> None:
    """气相 fastox 之后：底格 Eq.2.4/2.5 总量闭合 + 两相分相 Startwert。"""
    from src.solvers.vorabrechnung.vorab_x0 import (
        apply_bed0_eq24_two_phase_startwert,
        apply_bed0_inert_exchange_split_startwert,
        apply_bed0_r1_oxidizer_seed,
    )

    reactor._vorab_bed0_eq24_two_phase_diag = apply_bed0_eq24_two_phase_startwert(
        reactor.cells,
        reactor.config,
    )
    reactor._vorab_bed0_r1_oxidizer_seed_diag = apply_bed0_r1_oxidizer_seed(
        reactor.cells,
        reactor.config,
    )
    # 种子/fastox 改了 y 与 H2O/CO2 总和；N2/H2O/CO2 须按当前态重贴 Eq.2.2。
    reactor._vorab_bed0_inert_split_diag = apply_bed0_inert_exchange_split_startwert(
        reactor.cells,
        reactor.config,
    )


def _run_init_bed0_phase_split_closure(reactor: Reactor) -> None:
    """Init 末段：bed0–bed1 major phase-split（rm=1），须在 align/abgleich 之后。"""
    from src.core.connectivity import preproject_bed_major_gas_phase_split_to_exchange_closure

    cfg = reactor.config
    if not bool(getattr(cfg, "nr_init_bottom_major_phase_split_preprojection_thesis", False)):
        reactor._init_bed0_phase_split_closure_diag = None
        return

    max_bed_index = 0
    if bool(getattr(cfg, "nr_init_bed1_phase_split_preprojection_thesis", False)):
        max_bed_index = 1

    bed_indices = [
        i
        for i, cell in enumerate(reactor.cells)
        if str(getattr(cell, "cell_type", "bed")) == "bed" and i <= max_bed_index
    ]
    accepted: dict[str, bool] = {}
    for bed_index in bed_indices:
        accepted[f"bed{bed_index}_accepted"] = bool(
            preproject_bed_major_gas_phase_split_to_exchange_closure(
                reactor.cells,
                cfg,
                cell_index=bed_index,
                improvement_factor=0.85 if bed_index == 0 else 0.95,
                preserve_primary_inlet_totals=(bed_index == 0),
                apply_bc_fn=reactor._apply_all_bc_for_nr,
            )
        )
        reactor._apply_all_bc_for_nr()
    reactor._init_bed0_phase_split_closure_diag = accepted


def _run_init_abgleich_and_joint_preprojection(reactor: Reactor, state: _InitPrecalcState) -> None:
    """阶段 5：cell mapping 对齐、init abgleich、pyrolysis seed、hydro freeze、joint 预投影。"""
    from src.core.connectivity import (
        align_bottom_primary_gas_state_to_inlet_split,
        run_bottom_joint_start_value_preprojection,
    )
    from src.solvers.vorabrechnung import seed_pyrolysis_gas_dense_holdup_for_nr_x0

    cfg = reactor.config
    cell_budgets = state.cell_budgets
    if state.use_cell_mapping and cell_budgets is not None:
        align_bottom_primary_gas_state_to_inlet_split(reactor.cells, cfg)
        reactor._vorab_a_tier_cell_mapping_diag = {
            "bed_cells_mapped": float(sum(1 for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed")),
            "final_primary_inlet_split_align": True,
        }
        reactor._apply_all_bc_for_nr()
    else:
        reactor._vorab_a_tier_cell_mapping_diag = None
    # Hamel Kap.2.1 / thesis single-shot：先冻结干燥/热解源，再做 Check2 abgleich。
    # 否则 VM 分区格 `_vm_cache_valid=False` 时产物 rel 会对 R_gas≈0 爆炸。
    reactor._initialize_thesis_single_shot_vorab_sources()
    if state.use_a_tier_budget and cell_budgets is not None and bool(
        getattr(cfg, "vorab_a_tier_init_abgleich_bridge_thesis", True)
    ):
        from src.solvers.vorab_cell_abgleich import run_init_vorab_cell_abgleich_bridge

        reactor._vorab_a_tier_init_abgleich_diag = run_init_vorab_cell_abgleich_bridge(
            reactor.cells,
            cell_budgets,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            reseed_holdup_chain_fn=reactor._seed_initialized_holdup_chain_for_nr,
        )
    else:
        reactor._vorab_a_tier_init_abgleich_diag = None
    if bool(getattr(cfg, "nr_seed_pyrolysis_gas_dense_holdup_x0_thesis", False)):
        seed_diag = seed_pyrolysis_gas_dense_holdup_for_nr_x0(
            reactor.cells,
            holdup_scale=float(getattr(cfg, "nr_pyrolysis_dense_holdup_seed_scale_thesis", 1.0)),
            residual_closure_passes=2,
        )
        reactor._pyrolysis_dense_holdup_seed_diag = seed_diag
    else:
        reactor._pyrolysis_dense_holdup_seed_diag = None
    reactor._vorab_a_tier_cell_budgets = cell_budgets if state.use_a_tier_budget else None
    reactor._initialize_nr_hydrodynamics_freeze_cache()
    run_joint_preproj = bool(getattr(cfg, "nr_init_bottom_joint_preprojection_thesis", False))
    if state.use_cell_mapping and not bool(
        getattr(cfg, "nr_bottom_joint_preprojection_when_cell_mapped_thesis", False)
    ):
        run_joint_preproj = False
    if run_joint_preproj:
        bottom_joint = run_bottom_joint_start_value_preprojection(
            reactor.cells,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
        )
    else:
        bottom_joint = {
            "joint": False,
            "joint_gas_relax": 0.0,
            "joint_global_gate_pass": False,
            "joint_rms_scaled_before": 0.0,
            "joint_rms_scaled_after": 0.0,
            "joint_merit_scaled_before": 0.0,
            "joint_merit_scaled_after": 0.0,
            "joint_split_rms_before": 0.0,
            "joint_split_rms_after": 0.0,
            "joint_split_aware_gate": False,
            "bc_aware_bed0_bed1_joint": False,
            "bed1_guard_tripped": False,
            "bed1_energy_abs_before_W": 0.0,
            "bed1_energy_abs_after_W": 0.0,
            "skipped_for_vorab_cell_mapping": bool(state.use_cell_mapping),
        }
    reactor._bottom_joint_preprojection_run_diag = bottom_joint
    _run_init_bed0_phase_split_closure(reactor)
    _reconcile_upper_bed_solid_holdup_if_enabled(reactor)
    reactor._apply_all_bc_for_nr()
    # 放在 abgleich / phase-split / holdup 之后：否则预算重播会冲掉快氧化 Startwert
    if bool(getattr(cfg, "vorab_transport_x0_fast_oxidation_closure_thesis", True)):
        from src.solvers.vorabrechnung.vorab_x0 import project_bed_holdup_fast_oxidation_closure

        reactor._vorab_transport_x0_fast_oxidation_diag = project_bed_holdup_fast_oxidation_closure(
            reactor.cells, cfg
        )
        reactor._apply_all_bc_for_nr()
    else:
        reactor._vorab_transport_x0_fast_oxidation_diag = None
    _run_bed0_eq24_two_phase_startwert(reactor)
    reactor._apply_all_bc_for_nr()
    _run_bottom_cell_startwert_solve(reactor)
    if bool(getattr(cfg, "vorab_transport_x0_fast_oxidation_closure_thesis", True)):
        reactor._apply_all_bc_for_nr()
        reactor._refresh_explicit_freeboard_transport_from_closure()
        reactor._initialize_explicit_side_block_states()
        reactor._apply_all_bc_for_nr()
    # fastox / Eq.2.4 / R1 种子会扰动 K 与固相入流。再对齐一次，否则 bed1/bed9
    # 欠库存 → m_in>m_out 焓洞（Eq.2.7，handoff §5.54）。允许降库存，把过高格拉回 support。
    _reconcile_upper_bed_solid_holdup_if_enabled(reactor, allow_reduce=True)
    reactor._apply_all_bc_for_nr()
    # 全床联立对齐后 bed1 char 仍欠 ~10%；对其单独 snap（handoff §5.55）。
    _snap_bed_holdup_to_local_support(reactor, bed_indices=(1,), max_passes=3)
    # snap/BC 后再贴一次惰性分相，让 NR 看到的 N2/H2O/CO2 与当前 N_ex 一致（§5.56）。
    from src.solvers.vorabrechnung.vorab_x0 import apply_bed0_inert_exchange_split_startwert

    reactor._vorab_bed0_inert_split_final_diag = apply_bed0_inert_exchange_split_startwert(
        reactor.cells,
        reactor.config,
    )
    from src.solvers.vorabrechnung.vorab_cell_mapping import (
        apply_post_fastox_post_staged_h2o_passthrough,
        apply_post_fastox_post_staged_syngas_passthrough,
    )

    reactor._apply_all_bc_for_nr()
    reactor._vorab_post_fastox_h2o_passthrough_diag = apply_post_fastox_post_staged_h2o_passthrough(
        reactor.cells,
        reactor.config,
    )
    reactor._vorab_post_fastox_syngas_passthrough_diag = (
        apply_post_fastox_post_staged_syngas_passthrough(
            reactor.cells,
            reactor.config,
        )
    )
    reactor._apply_all_bc_for_nr()
    _apply_bed0_nr_startwert_T(reactor)


def run_init_and_precalc_for_global_nr(
    reactor: Reactor,
    *,
    init_strategy: str | None,
    gs_warmup_steps: int | None,
) -> GlobalNRInitPrecalcResult:
    """执行轴向温度估计、进料/上游传播、初值 x0、水动力冻结与 thesis 单次源项初始化。

    编排顺序（OPT-013）：
    1. ``_resolve_init_strategy_and_seed_T``
    2. ``_run_vorabrechnung_snapshot``
    3. ``_run_post_vorabrechnung_holdup_and_reactive``
    4. ``_run_bottom_zone_preprojections``
    5. ``_run_init_abgleich_and_joint_preprojection``
    """
    del gs_warmup_steps  # 保留接口；Global NR 路径不使用 GS warmup

    state = _resolve_init_strategy_and_seed_T(reactor, init_strategy=init_strategy)
    state = _run_vorabrechnung_snapshot(reactor, state)
    saved_u0_targets = _apply_init_bed_vorab_gas_u0_targets(
        reactor,
        T_est=state.T_est,
        cell_budgets=state.cell_budgets,
    )
    try:
        _run_post_vorabrechnung_holdup_and_reactive(reactor, state)
    finally:
        _restore_macro_hydrodynamics_seed(reactor, saved_u0_targets)
    _run_bottom_zone_preprojections(reactor, state)
    _run_init_abgleich_and_joint_preprojection(reactor, state)

    nr_vorabrechnung_s = perf_counter() - state.vorab_started
    nr_init_s_total = perf_counter() - state.init_started

    return GlobalNRInitPrecalcResult(
        resolved_init_strategy=state.resolved,
        nr_init_s_total=float(nr_init_s_total),
        nr_vorabrechnung_s=float(nr_vorabrechnung_s),
        T_est=np.asarray(state.T_est, dtype=np.float64),
    )
