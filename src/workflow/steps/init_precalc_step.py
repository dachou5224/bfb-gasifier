"""流程图：赋初值 + Vorabrechnung 快照（轴向 T、x0、水动力冻结）。

对应 mmd：InitVal → PRECALC（在进入 NR/Abgleich 前完成初始化段）。
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from time import perf_counter
from typing import TYPE_CHECKING

import numpy as np

from src.core.constants import Rg

if TYPE_CHECKING:
    from src.core.reactor import Reactor


@dataclass(frozen=True)
class GlobalNRInitPrecalcResult:
    """INIT + 预计算阶段计时与中间量（供结果汇总与审计）。"""

    resolved_init_strategy: str
    nr_init_s_total: float  # [s]
    nr_vorabrechnung_s: float  # [s]
    T_est: np.ndarray  # [K]


def _prepare_macro_hydrodynamics_seed(
    reactor: "Reactor",
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
    reactor: "Reactor",
    original_u0_targets: list[float | None],
) -> None:
    """Restore user/configured hydrodynamics targets after Vorabrechnung x0 generation."""
    for cell, original_u0 in zip(reactor.cells, original_u0_targets):
        cell.u0_target = original_u0


def run_init_and_precalc_for_global_nr(
    reactor: "Reactor",
    *,
    init_strategy: str | None,
    gs_warmup_steps: int | None,
) -> GlobalNRInitPrecalcResult:
    """执行轴向温度估计、进料/上游传播、初值 x0、水动力冻结与 thesis 单次源项初始化。"""
    from src.core.reactor import _resolve_nr_init_strategy
    from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0

    cfg = reactor.config
    resolved = _resolve_nr_init_strategy(init_strategy, gs_warmup_steps, allow_legacy_gs=False)

    init_started = perf_counter()
    vorab_started = perf_counter()
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
    )
    for i, cell in enumerate(reactor.cells):
        cell.T = float(T_est[i])
    reactor._set_bottom_cell_feeds()
    for i in range(len(reactor.cells)):
        reactor._propagate_upstream(i)
    original_u0_targets = _prepare_macro_hydrodynamics_seed(reactor, T_est)
    try:
        for cell in reactor.cells:
            cell.calc_hydrodynamics()

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
        )
    finally:
        _restore_macro_hydrodynamics_seed(reactor, original_u0_targets)
    # The explicit freeboard closure consumes bed-top entrained solid rates via
    # ``top._solid_upflow_rates()``. Under holdup_transport these rates are
    # zero until the thesis NR boundary path has populated ``K_solid_auf/K_ab``.
    # Refresh BC first so particle-trace solids become real freeboard holdup
    # instead of being silently dropped during the first closure sync.
    reactor._apply_all_bc_for_nr()
    reactor._refresh_explicit_freeboard_transport_from_closure()
    reactor._initialize_explicit_side_block_states()
    reactor._initialize_nr_hydrodynamics_freeze_cache()
    reactor._initialize_thesis_single_shot_vorab_sources()
    nr_vorabrechnung_s = perf_counter() - vorab_started
    nr_init_s_total = perf_counter() - init_started

    return GlobalNRInitPrecalcResult(
        resolved_init_strategy=resolved,
        nr_init_s_total=float(nr_init_s_total),
        nr_vorabrechnung_s=float(nr_vorabrechnung_s),
        T_est=np.asarray(T_est, dtype=np.float64),
    )
