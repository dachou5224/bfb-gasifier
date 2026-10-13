#!/usr/bin/env python3
"""bed9/bed0 能量残差：投影/postinner A/B（handoff §5.18）。

基线栈：ps + 延迟 rate-gate + hot_restart。
试验：postinner fb T 投影；末态 oneshot bed T→能量闭合；扩大 solid/energy preinner 上界。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_bed9_energy_closure_ab.py
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
from audit_gate_hot_topk_residuals import _run as _run_stack  # noqa: F401 — not used
from audit_local2_ps_abs_rms_stack import _topk_residuals
from audit_rate_gate_rms_stack_ab import GATE_AFTER, HOT, ISOMORPH, PS, _bed2_R

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.connectivity.preprojection import (
    preproject_bed_temperature_to_local_energy_closure,
)
from src.core.species import GAS_SPECIES_INDEX
from src.solvers.nr_indexing import has_temperature_var, n_gas_var, temperature_offset


STACK = {**GATE_AFTER, **HOT}


def _energy_snapshot(reactor: Any, beds: tuple[int, ...] = (0, 9)) -> dict[str, Any]:
    topk, groups, rms_hat, max_abs = _topk_residuals(reactor, k=20)
    out: dict[str, Any] = {
        "rms_hat": float(rms_hat),
        "max_abs": float(max_abs),
        "groups": {k: float(v) for k, v in groups.items()},
        "T": {},
        "F_hat_T": {},
    }
    # locate T residuals
    apply = reactor._apply_all_bc_for_nr
    from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor

    cells = reactor.cells
    cfg = reactor.config
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
    offset = 0
    for i, cell in enumerate(cells):
        from src.solvers.nr_indexing import n_var

        nv = int(n_var(cell))
        if i in beds and has_temperature_var(cell):
            tloc = int(temperature_offset(cell))
            out["T"][f"bed{i}"] = float(cell.T)
            out["F_hat_T"][f"bed{i}"] = float(F_hat[offset + tloc])
        offset += nv
    # top energy-like
    out["topk_energyish"] = [
        r
        for r in topk
        if r.get("kind") in ("energy", "solid_or_other")
        or (r.get("kind") == "gas" and r.get("species") in ("N2", "H2O"))
    ][:8]
    out["topk"] = topk[:8]
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
    snap = _energy_snapshot(reactor)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", snap["rms_hat"])),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "converged": bool(result.get("converged", False)),
        "bed2_RH2O": float(_bed2_R(reactor)["R_H2O"]),
        "snap": snap,
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} E={out['energy']:.4f} "
        f"bed9_T_Fhat={snap['F_hat_T'].get('bed9', float('nan')):+.3f} "
        f"bed0_T_Fhat={snap['F_hat_T'].get('bed0', float('nan')):+.3f} "
        f"T9={snap['T'].get('bed9', float('nan')):.1f} wall={wall_s:.0f}s",
        flush=True,
    )
    return out


def _oneshot_bed_T_project(reactor: Any, beds: tuple[int, ...]) -> dict[str, Any]:
    before = _energy_snapshot(reactor)
    n_ok = 0
    for i in beds:
        ok = preproject_bed_temperature_to_local_energy_closure(
            reactor.cells,
            reactor.config,
            cell_index=int(i),
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            apply_local_bc_fn=getattr(reactor, "_apply_local_bc_for_nr", None),
            temperature_bounds_K=(650.0, 1600.0),
        )
        if ok:
            n_ok += 1
    reactor._apply_all_bc_for_nr()
    after = _energy_snapshot(reactor)
    return {"n_ok": n_ok, "before": before, "after": after}


def main() -> None:
    arms = [
        _run("gate_hot_base", {}),
        _run(
            "postinner_fb_T_proj",
            {"nr_outer_postinner_freeboard_energy_T_projection_thesis": True},
        ),
        _run(
            "solid_energy_maxbed9",
            {"nr_outer_preinner_solid_energy_max_bed_index_thesis": 9},
        ),
    ]

    # oneshot on base stack end-state
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(STACK)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall = timed_phase2_solve(reactor, **phase2_solve_kwargs(max_global_iter=8))
    oneshot = {
        "bed9": _oneshot_bed_T_project(reactor, (9,)),
        # restore? re-solve too expensive; do bed9+bed0 on fresh? use current after bed9
    }
    # re-get fresh for bed0+bed9 joint
    reactor2 = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor2)
    timed_phase2_solve(reactor2, **phase2_solve_kwargs(max_global_iter=8))
    oneshot["bed0_and_bed9"] = _oneshot_bed_T_project(reactor2, (0, 9))
    co = float((result.get("bed_exit_gas") or {}).get("CO", float("nan")))
    print(
        f"[oneshot bed9] Fhat9 {oneshot['bed9']['before']['F_hat_T'].get('bed9'):+.3f}→"
        f"{oneshot['bed9']['after']['F_hat_T'].get('bed9'):+.3f} "
        f"rms_hat {oneshot['bed9']['before']['rms_hat']:.4f}→{oneshot['bed9']['after']['rms_hat']:.4f}",
        flush=True,
    )
    print(
        f"[oneshot bed0+9] Fhat9 {oneshot['bed0_and_bed9']['before']['F_hat_T'].get('bed9'):+.3f}→"
        f"{oneshot['bed0_and_bed9']['after']['F_hat_T'].get('bed9'):+.3f} "
        f"Fhat0 {oneshot['bed0_and_bed9']['before']['F_hat_T'].get('bed0'):+.3f}→"
        f"{oneshot['bed0_and_bed9']['after']['F_hat_T'].get('bed0'):+.3f} "
        f"rms_hat {oneshot['bed0_and_bed9']['before']['rms_hat']:.4f}→"
        f"{oneshot['bed0_and_bed9']['after']['rms_hat']:.4f}",
        flush=True,
    )

    out = data_path("phase2_bed9_energy_closure_ab_o8.json")
    write_json(
        out,
        {
            "note": "bed9/bed0 energy closure A/B on gate+hot stack; §5.18",
            "base_CO_from_oneshot_run": co,
            "arms": arms,
            "oneshot": {
                "bed9": {
                    "n_ok": oneshot["bed9"]["n_ok"],
                    "before_F": oneshot["bed9"]["before"]["F_hat_T"],
                    "after_F": oneshot["bed9"]["after"]["F_hat_T"],
                    "before_rms_hat": oneshot["bed9"]["before"]["rms_hat"],
                    "after_rms_hat": oneshot["bed9"]["after"]["rms_hat"],
                    "before_groups": oneshot["bed9"]["before"]["groups"],
                    "after_groups": oneshot["bed9"]["after"]["groups"],
                },
                "bed0_and_bed9": {
                    "n_ok": oneshot["bed0_and_bed9"]["n_ok"],
                    "before_F": oneshot["bed0_and_bed9"]["before"]["F_hat_T"],
                    "after_F": oneshot["bed0_and_bed9"]["after"]["F_hat_T"],
                    "before_rms_hat": oneshot["bed0_and_bed9"]["before"]["rms_hat"],
                    "after_rms_hat": oneshot["bed0_and_bed9"]["after"]["rms_hat"],
                    "before_groups": oneshot["bed0_and_bed9"]["before"]["groups"],
                    "after_groups": oneshot["bed0_and_bed9"]["after"]["groups"],
                },
            },
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
