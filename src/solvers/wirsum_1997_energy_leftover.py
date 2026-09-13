"""holdup Startwert 之后，把 Eq. 2.7 残差拆成可加总的桶。

``M = m_in / K`` 填上 char/ash 后 bed0 仍约 −5 MW。本核只归因，
不给 VM 加 K，不把 VM 写入 H_aus，不加 T clip。

桶之和必须等于 ``H_in − H_out − Q``。

Source: Hamel (1999) Eq. 2.7–2.8
"""

from __future__ import annotations

from typing import Any

import numpy as np

from src.solvers.wirsum_1997_energy_decomp import decompose_reactor_eq27
from src.solvers.wirsum_1997_holdup_startwert import apply_km_holdup_startwert_to_cells


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def attribute_eq27_leftover(row: dict[str, Any]) -> dict[str, Any]:
    """把一格 Eq. 2.7 残差拆成互斥桶 [W]。"""
    hin = row["H_in_parts_W"]
    h_zu = row["H_zu_by_comp_W"]
    h_sout = row["H_solid_out_by_comp_W"]
    gas_in = (
        float(hin["H_gas_axial_in"] or 0.0)
        + float(hin["H_gas_zu"] or 0.0)
        + float(hin["H_gas_rez"] or 0.0)
    )
    buckets = {
        "char_H_zu_minus_H_sout": float(h_zu["char"] or 0.0) - float(h_sout["char"] or 0.0),
        "vm_H_zu_minus_H_sout": float(h_zu["vm"] or 0.0) - float(h_sout["vm"] or 0.0),
        "moisture_H_zu_minus_H_sout": float(h_zu["moisture"] or 0.0)
        - float(h_sout["moisture"] or 0.0),
        "ash_H_zu_minus_H_sout": float(h_zu["ash"] or 0.0) - float(h_sout["ash"] or 0.0),
        "other_solid_in": (
            float(hin["H_solid_rez"] or 0.0)
            + float(hin["H_solid_axial_in"] or 0.0)
            + float(hin["H_solid_auf_in"] or 0.0)
            + float(hin["H_solid_ab_in"] or 0.0)
        ),
        "gas_H_in_minus_H_out": gas_in - float(row["H_gas_out_W"] or 0.0),
        "minus_Q": -float(row["Q_loss_W"] or 0.0),
    }
    residual = float(row["residual_W"] or 0.0)
    summed = float(sum(buckets.values()))
    abs_sum = sum(abs(v) for v in buckets.values())
    share = {
        k: (abs(v) / abs_sum if abs_sum > 0.0 else 0.0) for k, v in buckets.items()
    }
    dominant = max(buckets, key=lambda k: abs(buckets[k]))
    return {
        "cell": row.get("cell"),
        "residual_W": _finite(residual),
        "buckets_W": {k: _finite(v) for k, v in buckets.items()},
        "abs_share": {k: _finite(v) for k, v in share.items()},
        "bucket_sum_W": _finite(summed),
        "sum_matches_residual": bool(abs(summed - residual) <= 1.0),
        "dominant_bucket": dominant,
        "holdup_zero": bool(row.get("holdup_zero")),
        "vm_in_not_in_energy_out": bool(row.get("vm_in_not_in_energy_out")),
    }


def build_energy_leftover_payload() -> dict[str, Any]:
    """reduced 3-cell：先 ``M=m_in/K``，再拆剩余能量残差。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    apply_km_holdup_startwert_to_cells(list(reactor.cells))
    reactor._apply_all_bc_for_nr()
    for cell in reactor.cells:
        cell.residuals()
    eq27 = decompose_reactor_eq27(list(reactor.cells))
    per_cell = [attribute_eq27_leftover(row) for row in eq27["per_cell"]]
    bed0 = per_cell[0]
    return {
        "source": "Wirsum alignment Eq.2.7 leftover after km holdup Startwert",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "name": "eq27_energy_leftover",
        "not": [
            "K_vm",
            "vm_in_eq27_outlet",
            "temperature_fence",
            "energy_clip",
        ],
        "units": "W",
        "eq27": {
            "bed0_residual_W": eq27["bed0_residual_W"],
            "bed0_holdup_zero": eq27["bed0_holdup_zero"],
        },
        "per_cell": per_cell,
        "bed0_dominant_bucket": bed0["dominant_bucket"],
        "note": (
            "char/ash 用 M=m_in/K 闭合后，剩余残差应能加总。"
            "VM 桶是 K_VM=0 的结构项，不是调参对象。"
            "bed1/2 若仍 holdup_zero，是单次扫描未接力，不是给 VM 加 K。"
            "不用实验 CO。"
        ),
    }
