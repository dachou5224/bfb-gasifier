#!/usr/bin/env python3
"""bed9 固相 dm 分解 + 延迟床顶焓透传 A/B（handoff §5.20）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed9_dm_deferred_passthrough.py
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

from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import holdup_transport_char_ash_mask

COMP = ("char", "VM", "moisture", "ash")
COMP_IDX = (S_CHAR, S_VM, S_MOISTURE, S_ASH)


def _sum_comp(arr: np.ndarray) -> dict[str, float]:
    out = {name: float(np.sum(arr[:, idx])) for name, idx in zip(COMP, COMP_IDX)}
    out["total"] = float(np.sum(arr))
    return out


def _bed9_mass_decomp(reactor: Any) -> dict[str, Any]:
    cell = reactor.cells[9]
    reactor._apply_all_bc_for_nr()
    cell.residuals()
    mig = np.asarray(cell._calc_size_migration(), dtype=np.float64)
    R = np.asarray(cell.R_solid, dtype=np.float64)
    m = np.maximum(np.asarray(cell.m_solid, dtype=np.float64), 0.0)
    K_auf = np.maximum(np.asarray(cell.K_solid_auf, dtype=np.float64), 0.0)
    K_ab = np.maximum(np.asarray(cell.K_solid_ab, dtype=np.float64), 0.0)
    m_out = (K_auf + K_ab) * m
    m_auf = np.maximum(np.asarray(cell.m_solid_auf_in, dtype=np.float64), 0.0)
    m_ab = np.maximum(np.asarray(cell.m_solid_ab_in, dtype=np.float64), 0.0)
    m_in = np.maximum(np.asarray(cell.m_solid_in, dtype=np.float64), 0.0)
    m_zu = np.maximum(np.asarray(cell.m_solid_zu, dtype=np.float64), 0.0)
    m_rez = np.maximum(np.asarray(cell.m_solid_rez, dtype=np.float64), 0.0)
    m_in_all = m_auf + m_ab + m_in + m_zu + m_rez
    mask = holdup_transport_char_ash_mask(m.shape)
    m_out_masked = (K_auf * mask + K_ab * mask) * m
    m_in_masked = (m_auf + m_ab) * mask + m_in + m_zu + m_rez
    res = np.asarray(cell.calc_solid_balance(), dtype=np.float64)
    # 稳态：m_out ≈ m_in + R + mig（char/ash 输运支）
    budget = m_in_all + R + mig - m_out
    return {
        "m_in_all": _sum_comp(m_in_all),
        "m_auf_in": _sum_comp(m_auf),
        "m_ab_in": _sum_comp(m_ab),
        "m_solid_in": _sum_comp(m_in),
        "m_zu": _sum_comp(m_zu),
        "m_rez": _sum_comp(m_rez),
        "m_out_Km": _sum_comp(m_out),
        "m_out_masked": _sum_comp(m_out_masked),
        "m_in_masked_like_balance": _sum_comp(m_in_masked),
        "R_solid": _sum_comp(R),
        "size_migration": _sum_comp(mig),
        "dm_out_minus_in": _sum_comp(m_out - m_in_all),
        "R_plus_mig": _sum_comp(R + mig),
        "budget_in_R_mig_minus_out": _sum_comp(budget),
        "solid_res_sum": _sum_comp(res),
        "solid_res_abs_max": float(np.max(np.abs(res))),
        "holdup": _sum_comp(m),
        "K_auf_mean_char": float(np.mean(K_auf[:, S_CHAR])),
        "K_ab_mean_char": float(np.mean(K_ab[:, S_CHAR])),
        "K_auf_mean_vm": float(np.mean(K_auf[:, S_VM])),
        "K_ab_mean_vm": float(np.mean(K_ab[:, S_VM])),
    }


def _arm_history(result: dict[str, Any]) -> list[dict[str, Any]]:
    hist = list(result.get("nr_outer_history") or [])
    out: list[dict[str, Any]] = []
    for i, h in enumerate(hist):
        if not isinstance(h, dict):
            continue
        keys = (
            "bed_top_energy_passthrough_armed",
            "bed_top_energy_passthrough_newly_armed",
            "bed_top_energy_passthrough_y_co_top",
            "trace_o2_rate_gate_armed",
            "trace_o2_rate_gate_newly_armed",
        )
        if any(k in h for k in keys):
            row = {"outer": i}
            for k in keys:
                if k in h:
                    row[k] = h[k]
            out.append(row)
    # also dig preinner diags if nested
    for i, h in enumerate(hist):
        if not isinstance(h, dict):
            continue
        pre = h.get("preinner") or h.get("preinner_diag") or {}
        if isinstance(pre, dict) and "bed_top_energy_passthrough_armed" in pre:
            out.append({"outer": i, **{k: pre[k] for k in pre if "passthrough" in k or "y_co" in k}})
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
    bridge = _bridge_audit(reactor)
    decomp = _bed9_mass_decomp(reactor)
    topk, groups, rms_hat, _ = _topk_residuals(reactor, k=8)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", rms_hat)),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "bed2_RH2O": float(_bed2_R(reactor)["R_H2O"]),
        "groups": {k: float(v) for k, v in groups.items()},
        "topk": topk,
        "bridge": bridge,
        "bed9_mass_decomp": decomp,
        "passthrough_armed": bool(getattr(reactor.cells[9], "_solid_energy_passthrough", False)),
        "arm_hist": _arm_history(result),
        "overrides": dict(overrides),
    }
    e = bridge["bed9_enthalpy"]
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} "
        f"Fhat9={bridge['Fhat_T_bed9']:+.3f} dm={decomp['dm_out_minus_in']['total']:+.3e} "
        f"R+mig={decomp['R_plus_mig']['total']:+.3e} resE={e['res_Eout']:.3e} "
        f"pt={out['passthrough_armed']} wall={wall_s:.0f}s",
        flush=True,
    )
    print(
        f"  dm by comp: char={decomp['dm_out_minus_in']['char']:+.3e} "
        f"ash={decomp['dm_out_minus_in']['ash']:+.3e} "
        f"VM={decomp['dm_out_minus_in']['VM']:+.3e} "
        f"mois={decomp['dm_out_minus_in']['moisture']:+.3e}",
        flush=True,
    )
    print(
        f"  R: char={decomp['R_solid']['char']:+.3e} ash={decomp['R_solid']['ash']:+.3e} "
        f"mig_tot={decomp['size_migration']['total']:+.3e} "
        f"solid_res_max={decomp['solid_res_abs_max']:.3e}",
        flush=True,
    )
    return out


def main() -> None:
    arms = [
        _run("gate_hot_base", {}),
        _run(
            "passthrough_after_y_co",
            {
                "nr_bed_top_solid_energy_passthrough_thesis": True,
                "nr_bed_top_solid_energy_passthrough_after_y_co_thesis": True,
                "nr_bed_top_solid_energy_passthrough_y_co_threshold_thesis": 0.14,
            },
        ),
    ]
    out = data_path("phase2_bed9_dm_deferred_passthrough_o8.json")
    write_json(
        out,
        {
            "note": "bed9 solid dm decomp + deferred bed-top energy passthrough; §5.20",
            "arms": arms,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
