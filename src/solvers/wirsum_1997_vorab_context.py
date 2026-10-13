"""阶段 D：内层冻结的 Vorabrechnung 上下文。

把水力、``K_solid_*``、干燥/热解 cache 收成同一份 ``q``。内层只读；
外层按 ``vorab_outer_source_policy`` 决定是否完整重算 ``q = P(x)``。
本模块不改生产默认：``auto`` 仍是工程 single-shot。

自由板轨迹目前不属于外层：显式图跳过 closure 同步，非显式图调用后立即返回。

Source: Hamel (1999) §2.1 Bild 2.2; Wirsum (1997) §2.2;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 D
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

from src.core.cell import Cell
from src.core.cell_state_cache import enable_inner_nr_vorabrechnung_freeze
from src.solvers.global_nr_solver import global_residual, pack_reactor
from src.solvers.vorabrechnung.vorab_policy import (
    resolve_vorab_outer_source_policy,
    vorab_tau_model_label,
)
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    MAPPING_EXPR,
    example_lu_like_feed,
    two_outer_maps,
)
from src.solvers.vorabrechnung.vorab_refresh import (
    backup_vorabrechnung_hydro_freeze,
    restore_vorabrechnung_hydro_freeze,
)

VORAB_CONTEXT_VERSION: int = 1

HydroBackup = list[tuple[dict[str, float], np.ndarray, np.ndarray, bool, bool]]
ResidualApplyBc = Callable[[], None]


@dataclass
class VorabFrozenContext:
    """一份冻结的 ``q``。version 变了就不能拿旧快照复评。"""

    version: int
    hydro_backup: HydroBackup
    drying_pyro: list[dict[str, Any]]
    source_policy: str
    tau_model: str
    freeboard_trajectory_owner: str
    outer_updates_freeboard_trajectory: bool
    n_cells: int


def freeboard_trajectory_contract(reactor: Any | None = None) -> dict[str, Any]:
    """当前生产外层对自由板轨迹的所有权。"""
    explicit = False
    if reactor is not None:
        explicit = bool(reactor._use_explicit_freeboard_solver_graph())
    return {
        "explicit_solver_graph": explicit,
        "outer_updates_freeboard_trajectory": False,
        "owner": "init_precalc",
        "note": (
            "显式图跳过 closure 同步，以免覆盖 NR 气/T；"
            "非显式图调用 refresh 后因不是显式图立即返回。"
            "固体已从 NR 未知量删除。不要简单取消全部轨迹刷新。"
        ),
    }


def _snapshot_drying_pyro(cell: Cell) -> dict[str, Any]:
    return {
        "sources_frozen": bool(getattr(cell, "_vorab_drying_pyro_sources_frozen", False)),
        "vm_cache_valid": bool(getattr(cell, "_vm_cache_valid", False)),
        "vm_gas_source": np.array(cell._vm_gas_source_cache, dtype=np.float64, copy=True),
        "vm_solid_sink": np.array(cell._vm_solid_sink_cache, dtype=np.float64, copy=True),
        "T": float(cell._vm_cache_T),
        "tau": float(cell._vm_cache_tau),
        "T_init": float(cell._vm_cache_T_init),
        "m_vm_in": float(cell._vm_cache_m_vm_in),
        "m_moist_in": float(cell._vm_cache_m_moist_in),
    }


def _restore_drying_pyro(cell: Cell, snap: dict[str, Any]) -> None:
    cell._vorab_drying_pyro_sources_frozen = bool(snap["sources_frozen"])
    cell._vm_cache_valid = bool(snap["vm_cache_valid"])
    cell._vm_gas_source_cache[:] = snap["vm_gas_source"]
    cell._vm_solid_sink_cache[:] = snap["vm_solid_sink"]
    cell._vm_cache_T = float(snap["T"])
    cell._vm_cache_tau = float(snap["tau"])
    cell._vm_cache_T_init = float(snap["T_init"])
    cell._vm_cache_m_vm_in = float(snap["m_vm_in"])
    cell._vm_cache_m_moist_in = float(snap["m_moist_in"])


def snapshot_vorab_frozen_context(
    cells: list[Cell],
    *,
    reactor: Any | None = None,
    cfg: Any | None = None,
) -> VorabFrozenContext:
    """把当前 live 水力/K 与干燥热解 cache 收成内层只读上下文。"""
    if not cells:
        raise ValueError("snapshot_vorab_frozen_context 需要至少一个 cell")
    for cell in cells:
        cell._snapshot_vorabrechnung_hydrodynamics()
    fb = freeboard_trajectory_contract(reactor)
    policy_cfg = cfg if cfg is not None else getattr(reactor, "config", None)
    return VorabFrozenContext(
        version=VORAB_CONTEXT_VERSION,
        hydro_backup=backup_vorabrechnung_hydro_freeze(cells),
        drying_pyro=[_snapshot_drying_pyro(cell) for cell in cells],
        source_policy=(
            resolve_vorab_outer_source_policy(policy_cfg)
            if policy_cfg is not None
            else "unspecified"
        ),
        tau_model=vorab_tau_model_label(),
        freeboard_trajectory_owner=str(fb["owner"]),
        outer_updates_freeboard_trajectory=bool(fb["outer_updates_freeboard_trajectory"]),
        n_cells=int(len(cells)),
    )


def install_vorab_frozen_context(cells: list[Cell], context: VorabFrozenContext) -> None:
    """复原冻结 ``q`` 并打开内层 freeze。线搜索 / Jacobian / 回退共用这一份。"""
    if int(context.version) != VORAB_CONTEXT_VERSION:
        raise ValueError(
            f"VorabFrozenContext version {context.version} != {VORAB_CONTEXT_VERSION}"
        )
    if len(cells) != int(context.n_cells):
        raise ValueError(
            f"context n_cells={context.n_cells} 与 cells={len(cells)} 不一致"
        )
    restore_vorabrechnung_hydro_freeze(cells, context.hydro_backup)
    for cell, snap in zip(cells, context.drying_pyro):
        _restore_drying_pyro(cell, snap)
        enable_inner_nr_vorabrechnung_freeze(cell)


def eval_F_with_frozen_context(
    x: npt.ArrayLike,
    cells: list[Cell],
    apply_bc_fn: ResidualApplyBc,
    context: VorabFrozenContext,
) -> npt.NDArray[np.float64]:
    """``F(x; q_frozen)``。每次评残差前复原同一份上下文。"""
    install_vorab_frozen_context(cells, context)
    return global_residual(np.asarray(x, dtype=np.float64), cells, apply_bc_fn)


def _seed_reduced_holdup_reactor() -> Any:
    from src.solvers.wirsum_1997_holdup_startwert import apply_km_holdup_startwert_axial_sweep
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    apply_km_holdup_startwert_axial_sweep(list(reactor.cells), reactor._apply_all_bc_for_nr)
    return reactor


def build_stage_d_mapping_payload() -> dict[str, Any]:
    """两次外层映射的原子/焓账。不改生产 single-shot，不重写 N。"""
    feed = example_lu_like_feed()
    tau = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    common = dict(
        tau_s=tau,
        occupancy_source="named_test_age_bins",
        fuel_type="coal",
    )
    out = two_outer_maps(
        feed,
        state_a={**common, "P_Pa": 2.5e6, "T_K": 1100.0, "X_vm": np.array([0.0, 0.3, 0.7, 1.0])},
        state_b={**common, "P_Pa": 2.5e6, "T_K": 1200.0, "X_vm": np.array([0.0, 0.4, 0.8, 1.0])},
    )
    return {
        "source": "Wirsum alignment stage D pyrolysis mapping ledger",
        "solved": False,
        "mapping_expr": MAPPING_EXPR,
        "is_hamel_closed_form": False,
        "rewrites_unknown_N": out["rewrites_unknown_N"],
        "default_policy_unchanged": out["default_policy_unchanged"],
        "atom_C_gap_a": out["map_a"]["atom_C_gap"],
        "atom_C_gap_b": out["map_b"]["atom_C_gap"],
        "H_map_gap_a_W": out["map_a"]["H_map_gap_W"],
        "H_map_gap_b_W": out["map_b"]["H_map_gap_W"],
        "gas_hf_counted_once": out["map_a"]["gas_hf_counted_once"],
        "vm_pool_basis": out["map_a"]["vm_pool_basis"],
        "not_startwert": True,
    }


def build_stage_energy_outer_payload() -> dict[str, Any]:
    """能量行 + 约化外层 P。不改生产 single-shot。"""
    from src.solvers.wirsum_1997_coupled_f_energy import (
        build_stage_energy_outer_payload as _payload,
    )

    return _payload()


def build_stage_cd_coupled_payload() -> dict[str, Any]:
    """C 两级分 + D 热解源接入联立 F 的诊断包。不改生产默认。"""
    from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
        build_stage_cd_coupled_payload as _payload,
    )

    return _payload()


def build_stage_d_purity_payload() -> dict[str, Any]:
    """holdup Startwert 后：未冻结可历史依赖；冻结后 ``F(x),F(x)`` 与 ``F(x),F(y),F(x)`` 一致。"""
    live = _seed_reduced_holdup_reactor()
    x_live = pack_reactor(list(live.cells))
    f1 = global_residual(x_live.copy(), list(live.cells), live._apply_all_bc_for_nr)
    f2 = global_residual(x_live.copy(), list(live.cells), live._apply_all_bc_for_nr)
    unfrozen_delta = float(np.max(np.abs(f2 - f1)))

    frozen = _seed_reduced_holdup_reactor()
    cells = list(frozen.cells)
    x0 = pack_reactor(cells)
    ctx = snapshot_vorab_frozen_context(cells, reactor=frozen, cfg=frozen.config)
    fa = eval_F_with_frozen_context(x0, cells, frozen._apply_all_bc_for_nr, ctx)
    fb = eval_F_with_frozen_context(x0, cells, frozen._apply_all_bc_for_nr, ctx)
    y = x0 + 1.0e-6
    _ = eval_F_with_frozen_context(y, cells, frozen._apply_all_bc_for_nr, ctx)
    fc = eval_F_with_frozen_context(x0, cells, frozen._apply_all_bc_for_nr, ctx)
    return {
        "source": "Wirsum alignment stage D vorab frozen context",
        "solved": False,
        "path": "reduced_3cell_x0_plus_km_holdup_axial_sweep",
        "not_a_convergence_claim": True,
        "units": {"residual": "mixed SI packed F", "holdup": "kg"},
        "n_unknowns": int(x0.size),
        "source_policy": ctx.source_policy,
        "tau_model": ctx.tau_model,
        "freeboard_trajectory_owner": ctx.freeboard_trajectory_owner,
        "outer_updates_freeboard_trajectory": ctx.outer_updates_freeboard_trajectory,
        "unfrozen_max_abs_F_xx": unfrozen_delta,
        "frozen_max_abs_F_xx": float(np.max(np.abs(fb - fa))),
        "frozen_max_abs_F_xyx": float(np.max(np.abs(fc - fa))),
        "note": (
            "未冻结路径：apply_bc 用上一评留下的 R_solid/水力重算 K。"
            "冻结后内层只读 q。默认生产外层仍是 engineering_single_shot。"
        ),
    }
