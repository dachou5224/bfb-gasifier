"""固相库存 Startwert：对 K>0 的组分反演 ``M = m_in / K``。

reduced x0 在 ``m_zu`` 还是 0 时写 holdup，库存保持 0。
Eq. 2.7 的 ``H_aus,solid = K·M`` 因此为 0，而 ``H_zu`` 仍带入燃料焓。

本核只处理 ``K_auf+K_ab > 0`` 的组分（LU 上是 char / ash）。
``K_VM = K_moisture = 0`` 时 Eq. (2-6) 不约束库存（P1），不发明出流，
也不把 VM 补进能量出口。

不是 residual 锚（``holdup_k0_anchor``），也不是 ``source − M``。

Source: Hamel (1999) §2.1 Startwertwahl; Eq. 2.6 出流闭合 ``m = K M``
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.solvers.wirsum_1997_energy_decomp import decompose_reactor_eq27

_SOLID_NAMES = ("char", "vm", "moisture", "ash")
_SOLID_IDX = (S_CHAR, S_VM, S_MOISTURE, S_ASH)
_K_ATOL = 1.0e-12  # [1/s]


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def inlet_solid_streams(cell: Any) -> npt.NDArray[np.float64]:
    """当前格全部固相入口流 [kg/s]。"""
    return np.maximum(
        np.asarray(cell.m_solid_zu, dtype=np.float64)
        + np.asarray(cell.m_solid_rez, dtype=np.float64)
        + np.asarray(cell.m_solid_in, dtype=np.float64)
        + np.asarray(cell.m_solid_auf_in, dtype=np.float64)
        + np.asarray(cell.m_solid_ab_in, dtype=np.float64),
        0.0,
    )


def k_total_1_per_s(cell: Any) -> npt.NDArray[np.float64]:
    """``K_auf + K_ab`` [1/s]。"""
    return np.maximum(
        np.asarray(cell.K_solid_auf, dtype=np.float64)
        + np.asarray(cell.K_solid_ab, dtype=np.float64),
        0.0,
    )


def km_support_holdup(
    m_in: npt.ArrayLike,
    K: npt.ArrayLike,
    *,
    k_atol: float = _K_ATOL,
) -> npt.NDArray[np.float64]:
    """``M = m_in / K``，仅 ``K > k_atol`` [kg]。K=0 处保持 0。"""
    inflow = np.asarray(m_in, dtype=np.float64)
    k = np.asarray(K, dtype=np.float64)
    if inflow.shape != k.shape:
        raise ValueError(f"m_in shape {inflow.shape} != K shape {k.shape}")
    out = np.zeros_like(inflow)
    mask = k > float(k_atol)
    out[mask] = inflow[mask] / k[mask]
    return out


def apply_km_holdup_startwert(cell: Any) -> dict[str, Any]:
    """就地把 K>0 组分写成 ``M = m_in / K``。不改 F(x)，不写 VM 若 K=0。"""
    before = np.asarray(cell.m_solid, dtype=np.float64).copy()
    m_in = inlet_solid_streams(cell)
    k = k_total_1_per_s(cell)
    target = km_support_holdup(m_in, k)
    seeded = np.asarray(cell.m_solid, dtype=np.float64).copy()
    writable = k > _K_ATOL
    seeded[writable] = target[writable]
    cell.m_solid[:] = seeded
    if hasattr(cell, "_thermo_cache_valid"):
        cell._thermo_cache_valid = False
    per_comp = {}
    for name, idx in zip(_SOLID_NAMES, _SOLID_IDX):
        per_comp[name] = {
            "m_in_kg_s": _finite(float(np.sum(m_in[:, idx]))),
            "K_max_1_per_s": _finite(float(np.max(k[:, idx])) if k.size else 0.0),
            "M_before_kg": _finite(float(np.sum(before[:, idx]))),
            "M_after_kg": _finite(float(np.sum(seeded[:, idx]))),
            "writable": bool(np.any(writable[:, idx])),
        }
    return {
        "solid_state_model": str(getattr(cell, "solid_state_model", "")),
        "xi_M_kg": _finite(float(np.sum(np.abs(seeded - before)))),
        "per_comp": per_comp,
    }


def apply_km_holdup_startwert_to_cells(cells: list[Any]) -> dict[str, Any]:
    """逐格做 ``M = m_in / K``。只动 ``holdup_transport``。

    一次扫完全床：上层若此时还没有 ``m_auf_in``，库存仍为 0。
    轴向接力用 ``apply_km_holdup_startwert_axial_sweep``。
    """
    per_cell: list[dict[str, Any]] = []
    touched = 0
    for i, cell in enumerate(cells):
        if str(getattr(cell, "solid_state_model", "")) != "holdup_transport":
            per_cell.append({"cell": i, "skipped": True})
            continue
        diag = apply_km_holdup_startwert(cell)
        diag["cell"] = i
        diag["skipped"] = False
        if float(diag["xi_M_kg"] or 0.0) > 0.0:
            touched += 1
        per_cell.append(diag)
    return {
        "name": "km_holdup_startwert",
        "units": {"holdup": "kg", "K": "1/s", "inlet": "kg/s"},
        "cells_touched": touched,
        "per_cell": per_cell,
    }


def apply_km_holdup_startwert_axial_sweep(
    cells: list[Any],
    apply_bc_fn: Any,
) -> dict[str, Any]:
    """自下而上：种第 i 格后立刻刷新 BC，让第 i+1 格看见 ``K M`` 上行流。

    仍只反演 K>0 的 char/ash。不给 VM 加 K，不改 residual。
    """
    if apply_bc_fn is None:
        raise ValueError("axial sweep requires apply_bc_fn so upflow reaches the next cell")
    per_cell: list[dict[str, Any]] = []
    touched = 0
    for i, cell in enumerate(cells):
        if str(getattr(cell, "solid_state_model", "")) != "holdup_transport":
            per_cell.append({"cell": i, "skipped": True, "pass": i})
            apply_bc_fn()
            continue
        diag = apply_km_holdup_startwert(cell)
        apply_bc_fn()
        diag["cell"] = i
        diag["skipped"] = False
        if float(diag["xi_M_kg"] or 0.0) > 0.0:
            touched += 1
        per_cell.append(diag)
    return {
        "name": "km_holdup_startwert_axial_sweep",
        "units": {"holdup": "kg", "K": "1/s", "inlet": "kg/s"},
        "cells_touched": touched,
        "per_cell": per_cell,
    }


def build_holdup_startwert_ab_payload() -> dict[str, Any]:
    """reduced 3-cell：关 / 开 ``M = m_in / K`` 后的 Eq. 2.7。不跑 Newton。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    before = decompose_reactor_eq27(list(reactor.cells))
    applied = apply_km_holdup_startwert_to_cells(list(reactor.cells))
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    after = decompose_reactor_eq27(list(reactor.cells))
    return {
        "source": "Wirsum alignment km holdup Startwert A/B",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "km_holdup_startwert",
        "not": [
            "holdup_k0_anchor",
            "vm_in_eq27_outlet",
            "temperature_fence",
            "source_minus_M",
        ],
        "units": {"enthalpy": "W", "holdup": "kg", "K": "1/s"},
        "applied": applied,
        "before": before,
        "after": after,
        "note": (
            "generate_initial_x0 时 m_zu=0，holdup 被写成 0。"
            "BC 补上给料后按 M=m_in/K 只填 char/ash。"
            "VM 的 K=0，库存与 H_aus 仍不闭合；不把 VM 补进 Eq.2.7。"
        ),
    }


def _prepare_lu_cells():
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    return reactor


def build_holdup_sweep_ab_payload() -> dict[str, Any]:
    """一次扫描 vs 自下而上接力。不跑 Newton，不给 VM 加 K。"""
    from src.solvers.wirsum_1997_energy_leftover import attribute_eq27_leftover

    one = _prepare_lu_cells()
    apply_km_holdup_startwert_to_cells(list(one.cells))
    one._apply_all_bc_for_nr()
    for cell in one.cells:
        cell.residuals()
    one_eq = decompose_reactor_eq27(list(one.cells))

    sweep = _prepare_lu_cells()
    applied = apply_km_holdup_startwert_axial_sweep(
        list(sweep.cells), sweep._apply_all_bc_for_nr
    )
    for cell in sweep.cells:
        cell.residuals()
    sweep_eq = decompose_reactor_eq27(list(sweep.cells))
    return {
        "source": "Wirsum alignment km holdup axial sweep A/B",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "km_holdup_startwert_axial_sweep",
        "not": [
            "K_vm",
            "vm_in_eq27_outlet",
            "holdup_k0_anchor",
            "temperature_fence",
        ],
        "units": {"enthalpy": "W", "holdup": "kg"},
        "applied": applied,
        "one_shot": {
            "eq27": {
                "bed0_residual_W": one_eq["bed0_residual_W"],
                "per_cell_residual_W": [row["residual_W"] for row in one_eq["per_cell"]],
                "holdup_zero": [row["holdup_zero"] for row in one_eq["per_cell"]],
            },
            "leftover": [attribute_eq27_leftover(row) for row in one_eq["per_cell"]],
        },
        "axial_sweep": {
            "eq27": {
                "bed0_residual_W": sweep_eq["bed0_residual_W"],
                "per_cell_residual_W": [row["residual_W"] for row in sweep_eq["per_cell"]],
                "holdup_zero": [row["holdup_zero"] for row in sweep_eq["per_cell"]],
                "m_holdup_kg": [row["m_holdup_kg"] for row in sweep_eq["per_cell"]],
            },
            "leftover": [attribute_eq27_leftover(row) for row in sweep_eq["per_cell"]],
        },
        "note": (
            "一次扫描时 bed1 还没有 m_auf_in，other_solid_in 是空库存。"
            "接力后 char/ash 上行进入 bed1/2 的 M=m_in/K。"
            "bed0 的 VM 桶不变；不给 VM 加 K。"
        ),
    }
