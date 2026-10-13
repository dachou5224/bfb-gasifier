#!/usr/bin/env python3
"""bed9 char holdup 闭合诊断 + 修复 A/B（handoff §5.21）。

1) gate+hot 末态逐床 solid / support / flux
2) 同态 posthoc：allow_reduce align / 含 R 的 m* / bed9 solid LSQ
3) 全量 A/B：postinner reconcile、solid_max=9、bed-top solid LSQ 开关

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed9_holdup_closure_ab.py
"""

from __future__ import annotations

from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)
from audit_bed9_fb_solid_enthalpy import STACK, _bridge_audit
from audit_local2_ps_abs_rms_stack import _topk_residuals
from audit_rate_gate_rms_stack_ab import ISOMORPH, PS, _bed2_R

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.cell import S_ASH, S_CHAR
from src.core.connectivity.preprojection import preproject_bed_solid_holdup_to_transport_closure
from src.core.connectivity.solid_transport import (
    _transport_inflow_support_holdup,
    align_bed_solid_holdup_to_auf_upflow_closure,
    align_bed_solid_holdup_to_transport_inflow_support,
    max_upper_bed_solid_auf_flux_ratio_error,
    max_upper_bed_solid_stream_imbalance,
)


def _bed_solid_row(cell: Any, idx: int) -> dict[str, float]:
    cell.residuals()
    m = np.maximum(np.asarray(cell.m_solid, dtype=np.float64), 0.0)
    support = _transport_inflow_support_holdup(cell)
    m_in = (
        np.maximum(cell.m_solid_auf_in, 0.0)
        + np.maximum(cell.m_solid_ab_in, 0.0)
        + np.maximum(cell.m_solid_in, 0.0)
        + np.maximum(cell.m_solid_zu, 0.0)
        + np.maximum(cell.m_solid_rez, 0.0)
    )
    m_out = cell._solid_outflow_rates()
    R = np.asarray(cell.R_solid, dtype=np.float64)
    K = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    # 含反应汇的稳态目标（负 R 也计入）
    net = m_in + R
    m_star = np.zeros_like(m)
    mask = K > 1e-12
    m_star[mask] = np.maximum(net[mask], 0.0) / K[mask]
    res = cell.calc_solid_balance()
    return {
        "bed": float(idx),
        "m_char": float(np.sum(m[:, S_CHAR])),
        "m_ash": float(np.sum(m[:, S_ASH])),
        "support_char": float(np.sum(support[:, S_CHAR])),
        "m_star_char": float(np.sum(m_star[:, S_CHAR])),
        "m_in_char": float(np.sum(m_in[:, S_CHAR])),
        "m_out_char": float(np.sum(m_out[:, S_CHAR])),
        "R_char": float(np.sum(R[:, S_CHAR])),
        "res_char": float(np.sum(res[:, S_CHAR])),
        "res_abs_max": float(np.max(np.abs(res))),
        "over_support_char": float(np.sum(m[:, S_CHAR] - support[:, S_CHAR])),
        "over_mstar_char": float(np.sum(m[:, S_CHAR] - m_star[:, S_CHAR])),
        "K_sum_char_mean": float(np.mean(K[:, S_CHAR])),
    }


def _chain_table(reactor: Any) -> list[dict[str, float]]:
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    reactor._apply_all_bc_for_nr()
    return [_bed_solid_row(c, i) for i, c in enumerate(beds)]


def _metrics(reactor: Any, result: dict[str, Any] | None = None) -> dict[str, Any]:
    topk, groups, rms_hat, _ = _topk_residuals(reactor, k=6)
    exit_gas = dict((result or {}).get("bed_exit_gas") or {})
    if not exit_gas and result is None:
        # best-effort from top bed
        top = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"][-1]
        n = top.N_d + top.N_b
        tot = float(np.sum(n))
        from src.core.species import GAS_SPECIES_INDEX

        j = int(GAS_SPECIES_INDEX["CO"])
        exit_gas = {"CO": float(n[j] / max(tot, 1e-30))}
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]
    bridge = _bridge_audit(reactor)
    return {
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float((result or {}).get("rms_scaled_final", rms_hat)),
        "groups": {k: float(v) for k, v in groups.items()},
        "topk": topk,
        "bed2_RH2O": float(_bed2_R(reactor)["R_H2O"]),
        "Fhat_T_bed9": float(bridge["Fhat_T_bed9"]),
        "res_E_bed9": float(bridge["bed9_enthalpy"]["res_Eout"]),
        "dm_bed9": float(bridge["bed9_enthalpy"]["dm_K_vs_in"]),
        "flux_mismatch": float(max_upper_bed_solid_auf_flux_ratio_error(beds, min_bed_index=1)),
        "stream_imbalance": float(max_upper_bed_solid_stream_imbalance(beds, min_bed_index=1)),
        "bed9": _bed_solid_row(reactor.cells[9], 9),
    }


def _snapshot_m(reactor: Any) -> list[np.ndarray]:
    return [np.asarray(c.m_solid, dtype=np.float64).copy() for c in reactor.cells]


def _restore_m(reactor: Any, snap: list[np.ndarray]) -> None:
    for c, m in zip(reactor.cells, snap):
        c.m_solid[:] = m
        c._thermo_cache_valid = False
        c._hydro_cache_valid = False


def _posthoc_probes(reactor: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"before": _metrics(reactor)}
    snap = _snapshot_m(reactor)
    beds = [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]

    # A: allow_reduce align all upper
    _restore_m(reactor, snap)
    reactor._apply_all_bc_for_nr()
    changed = False
    for c in beds[1:]:
        changed |= align_bed_solid_holdup_to_transport_inflow_support(c, allow_reduce=True)
    reactor._apply_all_bc_for_nr()
    out["align_reduce_upper"] = {**_metrics(reactor), "changed": bool(changed)}

    # B: auf upflow closure on bed9
    _restore_m(reactor, snap)
    reactor._apply_all_bc_for_nr()
    changed = align_bed_solid_holdup_to_auf_upflow_closure(reactor.cells[9], allow_reduce=True)
    reactor._apply_all_bc_for_nr()
    out["auf_closure_bed9"] = {**_metrics(reactor), "changed": bool(changed)}

    # C: reaction-aware m* on bed9 char/ash
    _restore_m(reactor, snap)
    reactor._apply_all_bc_for_nr()
    cell = reactor.cells[9]
    cell.residuals()
    m_in = (
        np.maximum(cell.m_solid_auf_in, 0.0)
        + np.maximum(cell.m_solid_ab_in, 0.0)
        + np.maximum(cell.m_solid_in, 0.0)
        + np.maximum(cell.m_solid_zu, 0.0)
        + np.maximum(cell.m_solid_rez, 0.0)
    )
    R = np.asarray(cell.R_solid, dtype=np.float64)
    K = np.maximum(cell.K_solid_auf + cell.K_solid_ab, 0.0)
    for comp in (S_CHAR, S_ASH):
        mask = K[:, comp] > 1e-12
        net = m_in[:, comp] + R[:, comp]
        cell.m_solid[mask, comp] = np.maximum(net[mask], 0.0) / K[mask, comp]
    cell._thermo_cache_valid = False
    reactor._apply_all_bc_for_nr()
    out["mstar_R_bed9"] = _metrics(reactor)

    # D: solid LSQ bed9
    _restore_m(reactor, snap)
    reactor._apply_all_bc_for_nr()
    ok = preproject_bed_solid_holdup_to_transport_closure(
        reactor.cells,
        reactor.config,
        cell_index=9,
        apply_bc_fn=reactor._apply_all_bc_for_nr,
        apply_local_bc_fn=reactor._apply_local_bc_for_nr,
        skip_rms_scaled=0.0,
        improvement_factor=0.5,
    )
    out["lsq_bed9"] = {
        **_metrics(reactor),
        "accepted": bool(ok),
        "diag": dict(getattr(reactor.cells[9], "_bottom_zone_solid_preprojection_diag", {}) or {}),
    }

    _restore_m(reactor, snap)
    reactor._apply_all_bc_for_nr()
    return out


def _run(label: str, overrides: dict[str, object], max_iter: int = 8) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    met = _metrics(reactor, result)
    chain = _chain_table(reactor)
    hist = list(result.get("nr_outer_history") or [])
    holdup_flags = [
        {
            "outer": i,
            "preinner_holdup_reconciled": h.get("preinner_holdup_reconciled"),
            "preinner_run_holdup_reconcile": h.get("preinner_run_holdup_reconcile"),
            "postinner_holdup_reconciled": h.get("postinner_holdup_reconciled"),
            "solid_preprojections_accepted": ((h.get("zone_diag") or {}) if isinstance(h.get("zone_diag"), dict) else {}).get(
                "solid_preprojections_accepted"
            ),
        }
        for i, h in enumerate(hist)
        if isinstance(h, dict)
        and (
            "preinner_holdup_reconciled" in h
            or "postinner_holdup_reconciled" in h
            or "preinner_run_holdup_reconcile" in h
        )
    ]
    print(
        f"[{label}] CO={met['exit_CO']:.4f} rms={met['rms']:.4f} "
        f"Fhat9={met['Fhat_T_bed9']:+.3f} res_char9={met['bed9']['res_char']:+.3e} "
        f"over_sup={met['bed9']['over_support_char']:+.3e} "
        f"flux={met['flux_mismatch']:.3f} stream={met['stream_imbalance']:.3f} "
        f"wall={wall_s:.0f}s",
        flush=True,
    )
    return {
        "label": label,
        "wall_s": float(wall_s),
        "metrics": met,
        "chain": chain,
        "holdup_flags": holdup_flags,
        "overrides": dict(overrides),
        "_reactor": reactor,  # stripped before write
        "_result": result,
    }


def main() -> None:
    base = _run("gate_hot_base", {})
    reactor = base.pop("_reactor")
    base.pop("_result", None)
    posthoc = _posthoc_probes(reactor)
    for name, ph in posthoc.items():
        if name == "before":
            continue
        b9 = ph.get("bed9") or {}
        print(
            f"  [posthoc {name}] CO={ph.get('exit_CO', float('nan')):.4f} "
            f"Fhat9={ph.get('Fhat_T_bed9', float('nan')):+.3f} "
            f"res_char={b9.get('res_char', float('nan')):+.3e} "
            f"accepted={ph.get('accepted', ph.get('changed', ''))}",
            flush=True,
        )

    arms = [base]
    for label, ov in (
        (
            "postinner_holdup_reconcile",
            {"nr_outer_postinner_holdup_reconcile_thesis": True},
        ),
        (
            "solid_energy_maxbed9",
            {"nr_outer_preinner_solid_energy_max_bed_index_thesis": 9},
        ),
        (
            "bed_top_solid_lsq",
            {"nr_outer_preinner_bed_top_solid_holdup_thesis": True},
        ),
        (
            "bed_top_solid_lsq_after_y_co",
            {
                "nr_outer_preinner_bed_top_solid_holdup_thesis": True,
                "nr_outer_preinner_bed_top_solid_holdup_after_y_co_thesis": True,
            },
        ),
    ):
        arm = _run(label, ov)
        arm.pop("_reactor", None)
        arm.pop("_result", None)
        arms.append(arm)

    out = data_path("phase2_bed9_holdup_closure_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed9 char holdup closure + repair A/B; §5.21",
            "arms": arms,
            "posthoc": posthoc,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
