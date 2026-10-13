#!/usr/bin/env python3
"""bed9↔fb0 固相质量/焓连续审计 + 床顶焓透传 A/B（handoff §5.19）。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed9_fb_solid_enthalpy.py
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
from audit_local2_ps_abs_rms_stack import _topk_residuals
from audit_rate_gate_rms_stack_ab import GATE_AFTER, HOT, ISOMORPH, PS, _bed2_R

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.cell import S_ASH, S_CHAR
from src.core.cell_balances import (
    calc_energy_balance_H_in,
    calc_gas_enthalpy_flow,
    calc_solid_enthalpy_flow,
)
from src.core.connectivity.preprojection import _preproj_cell_energy_residual_abs
from src.core.species import cp_ash, cp_char, cp_sand
from src.solvers.nr_indexing import has_temperature_var, n_var, temperature_offset

STACK = {**GATE_AFTER, **HOT}


def _solid_kw(cell: Any) -> dict[str, Any]:
    return dict(
        ash_dry_wt=cell.solid.ash_dry_wt,
        VM_daf=cell.solid.VM_daf,
        h_f_dry=cell.solid.h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
        solid_enthalpy_mode=cell.solid_enthalpy_mode,
    )


def _flux_snap(cell: Any) -> dict[str, float]:
    cell.residuals()
    up = cell._solid_upflow_rates()
    down = cell._solid_downflow_rates()
    e_out = cell._solid_energy_outlet_rates()
    m_in = (
        np.maximum(cell.m_solid_in, 0.0)
        + np.maximum(cell.m_solid_auf_in, 0.0)
        + np.maximum(cell.m_solid_ab_in, 0.0)
        + np.maximum(cell.m_solid_zu, 0.0)
        + np.maximum(cell.m_solid_rez, 0.0)
    )
    return {
        "T": float(cell.T),
        "model": str(cell.solid_state_model),
        "m_up": float(np.sum(up)),
        "m_down": float(np.sum(down)),
        "m_e_out": float(np.sum(e_out)),
        "m_in_all": float(np.sum(m_in)),
        "m_auf_in": float(np.sum(cell.m_solid_auf_in)),
        "m_ab_in": float(np.sum(cell.m_solid_ab_in)),
        "m_hold": float(np.sum(np.maximum(cell.m_solid, 0.0))),
        "char_up": float(np.sum(up[:, S_CHAR])),
        "ash_up": float(np.sum(up[:, S_ASH])),
        "E_abs": float(_preproj_cell_energy_residual_abs(cell)),
        "passthrough": bool(getattr(cell, "_solid_energy_passthrough", False)),
    }


def _enthalpy_gap(cell: Any) -> dict[str, float]:
    cell.residuals()
    H_in = float(
        calc_energy_balance_H_in(
            N_b_in=cell.N_b_in,
            N_d_in=cell.N_d_in,
            T_in_gas=cell.T_in_gas,
            N_zu_b=cell.N_zu_b,
            N_zu_d=cell.N_zu_d,
            T_zu_gas=cell.T_zu_gas,
            N_rez_b=cell.N_rez_b,
            N_rez_d=cell.N_rez_d,
            T_rez_gas=cell.T_rez_gas,
            m_solid_rez=cell.m_solid_rez,
            T_rez_solid=cell.T_rez_solid,
            m_solid_in=cell.m_solid_in,
            T_in_solid=cell.T_in_solid,
            m_solid_zu=cell.m_solid_zu,
            T_zu_solid=cell.T_zu_solid,
            h_cache=cell._h_cache,
            m_solid_auf_in=cell.m_solid_auf_in,
            T_solid_auf_in=cell.T_solid_auf_in,
            m_solid_ab_in=cell.m_solid_ab_in,
            T_solid_ab_in=cell.T_solid_ab_in,
            ash_dry_wt=cell.solid.ash_dry_wt,
            VM_daf=cell.solid.VM_daf,
            h_f_dry=cell.solid.h_f_dry,
            solid_enthalpy_mode=cell.solid_enthalpy_mode,
        )
    )
    H_gas = float(calc_gas_enthalpy_flow(cell.N_b + cell.N_d, cell.T, cell._h_cache))
    H_sol_e = float(calc_solid_enthalpy_flow(cell._solid_energy_outlet_rates(), cell.T, **_solid_kw(cell)))
    H_sol_k = float(calc_solid_enthalpy_flow(cell._solid_outflow_rates(), cell.T, **_solid_kw(cell)))
    m_pass = (
        np.maximum(cell.m_solid_in, 0.0)
        + np.maximum(cell.m_solid_auf_in, 0.0)
        + np.maximum(cell.m_solid_ab_in, 0.0)
        + np.maximum(cell.m_solid_zu, 0.0)
        + np.maximum(cell.m_solid_rez, 0.0)
    )
    H_sol_pass = float(calc_solid_enthalpy_flow(m_pass, cell.T, **_solid_kw(cell)))
    return {
        "H_in": H_in,
        "H_gas_out": H_gas,
        "H_sol_Kout": H_sol_k,
        "H_sol_Eout": H_sol_e,
        "H_sol_pass": H_sol_pass,
        "res_Kout": H_in - H_gas - H_sol_k,
        "res_Eout": H_in - H_gas - H_sol_e,
        "res_pass": H_in - H_gas - H_sol_pass,
        "dm_K_vs_in": float(np.sum(cell._solid_outflow_rates()) - np.sum(m_pass)),
        "dm_E_vs_in": float(np.sum(cell._solid_energy_outlet_rates()) - np.sum(m_pass)),
    }


def _Fhat_T(reactor: Any, bed: int) -> float:
    from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor

    cells = reactor.cells
    cfg = reactor.config
    apply = reactor._apply_all_bc_for_nr
    apply()
    ref_gas = max(
        float(getattr(cfg, "O2_feed", 0.0) + getattr(cfg, "H2O_feed", 0.0) + getattr(cfg, "N2_feed", 0.0)),
        1.0,
    )
    ref_solid = max(abs(float(getattr(cfg, "fuel_feed", 1.0))), 1.0e-9)
    ref_energy = max(abs(float(getattr(cfg, "fuel_feed", 1.0))) * 20.0e6, 1.0e6)
    F = np.asarray(global_residual(pack_reactor(cells), cells, apply), dtype=np.float64)
    scale = np.asarray(build_equation_scales(cells, ref_gas, ref_solid, ref_energy), dtype=np.float64)
    F_hat = F / np.maximum(scale, 1e-30)
    off = 0
    for i, cell in enumerate(cells):
        nv = int(n_var(cell))
        if i == bed and has_temperature_var(cell):
            return float(F_hat[off + int(temperature_offset(cell))])
        off += nv
    return float("nan")


def _bridge_audit(reactor: Any) -> dict[str, Any]:
    bed9 = reactor.cells[9]
    bed8 = reactor.cells[8]
    fb0 = reactor.freeboard_cells[0]
    reactor._apply_all_bc_for_nr()
    return {
        "bed8": _flux_snap(bed8),
        "bed9": _flux_snap(bed9),
        "fb0": _flux_snap(fb0),
        "bed9_enthalpy": _enthalpy_gap(bed9),
        "continuity": {
            "bed8_up_minus_bed9_auf": float(
                np.sum(bed8._solid_upflow_rates()) - np.sum(bed9.m_solid_auf_in)
            ),
            "bed9_up_minus_fb0_auf_or_in": float(
                np.sum(bed9._solid_upflow_rates())
                - np.sum(np.maximum(fb0.m_solid_auf_in, 0.0) + np.maximum(fb0.m_solid_in, 0.0))
            ),
            "fb0_down_minus_bed9_ab": float(
                np.sum(fb0._solid_downflow_rates()) - np.sum(bed9.m_solid_ab_in)
            ),
        },
        "Fhat_T_bed9": _Fhat_T(reactor, 9),
        "Fhat_T_bed0": _Fhat_T(reactor, 0),
    }


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
    topk, groups, rms_hat, max_abs = _topk_residuals(reactor, k=8)
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
        "overrides": dict(overrides),
    }
    e = bridge["bed9_enthalpy"]
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} E={out['energy']:.4f} "
        f"Fhat9={bridge['Fhat_T_bed9']:+.3f} "
        f"res_E={e['res_Eout']:.3e} res_pass={e['res_pass']:.3e} "
        f"dm={e['dm_K_vs_in']:+.3e} wall={wall_s:.0f}s",
        flush=True,
    )
    c = bridge["continuity"]
    print(
        f"  continuity: b8up-b9auf={c['bed8_up_minus_bed9_auf']:+.3e} "
        f"b9up-fb_in={c['bed9_up_minus_fb0_auf_or_in']:+.3e} "
        f"fbdown-b9ab={c['fb0_down_minus_bed9_ab']:+.3e}",
        flush=True,
    )
    return out


def main() -> None:
    arms = [
        _run("gate_hot_base", {}),
        _run(
            "bed_top_energy_passthrough",
            {"nr_bed_top_solid_energy_passthrough_thesis": True},
        ),
    ]
    # counterfactual on base end-state: force passthrough attribute and remeasure
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    timed_phase2_solve(reactor, **phase2_solve_kwargs(max_global_iter=8))
    before = _bridge_audit(reactor)
    reactor.cells[9]._solid_energy_passthrough = True
    after = _bridge_audit(reactor)
    print(
        f"[posthoc passthrough] res_E {before['bed9_enthalpy']['res_Eout']:.3e}→"
        f"{after['bed9_enthalpy']['res_Eout']:.3e} "
        f"Fhat9 {before['Fhat_T_bed9']:+.3f}→{after['Fhat_T_bed9']:+.3f}",
        flush=True,
    )
    out = data_path("phase2_bed9_fb_solid_enthalpy_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed9↔fb solid mass/enthalpy + bed-top energy passthrough A/B; §5.19",
            "arms": arms,
            "posthoc_passthrough": {
                "before": before,
                "after": after,
            },
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
