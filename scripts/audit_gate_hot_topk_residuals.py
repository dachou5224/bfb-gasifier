#!/usr/bin/env python3
"""gate+hot 末态 top 残差分解（handoff §5.17）：谁挡 fully converge。

用法：
  PYTHONUNBUFFERED=1 python3 scripts/audit_gate_hot_topk_residuals.py
"""

from __future__ import annotations

from typing import Any

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


def _run(label: str, overrides: dict[str, object], max_iter: int) -> dict[str, Any]:
    cfg = dict(ISOMORPH)
    cfg.update(PS)
    cfg.update(overrides)
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg)
    run_phase2_init(reactor)
    result, wall_s = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=int(max_iter))
    )
    assert isinstance(result, dict)
    topk, groups, rms_hat, max_abs = _topk_residuals(reactor, k=15)
    exit_gas = dict(result.get("bed_exit_gas") or {})
    bed2 = _bed2_R(reactor)
    # 按 cell / kind 汇总 |F_hat|
    by_cell: dict[str, float] = {}
    by_kind: dict[str, float] = {}
    for row in topk:
        ckey = f"cell{row.get('cell', '?')}"
        by_cell[ckey] = by_cell.get(ckey, 0.0) + float(row["abs"])
        k = str(row.get("kind", "?"))
        if k == "gas":
            k = f"gas:{row.get('species', '?')}"
        by_kind[k] = by_kind.get(k, 0.0) + float(row["abs"])
    out = {
        "label": label,
        "wall_s": float(wall_s),
        "exit_CO": float(exit_gas.get("CO", float("nan"))),
        "rms": float(result.get("rms_scaled_final", rms_hat)),
        "split": float(result.get("rms_scaled_gas_phase_split_final", float("nan"))),
        "energy": float(result.get("rms_scaled_energy_final", float("nan"))),
        "converged": bool(result.get("converged", False)),
        "n_acc": sum(
            1 for x in (result.get("nr_accepted_lambda_history") or []) if x is not None
        ),
        "rms_hat": float(rms_hat),
        "max_abs": float(max_abs),
        "groups": {k: float(v) for k, v in groups.items()},
        "topk": topk,
        "topk_by_cell": dict(sorted(by_cell.items(), key=lambda kv: -kv[1])),
        "topk_by_kind": dict(sorted(by_kind.items(), key=lambda kv: -kv[1])),
        "bed2": bed2,
        "overrides": dict(overrides),
    }
    print(
        f"[{label}] CO={out['exit_CO']:.4f} rms={out['rms']:.4f} split={out['split']:.4f} "
        f"max|F|={max_abs:.3f} bed2_RH2O={bed2['R_H2O']:+.2f}",
        flush=True,
    )
    for row in topk[:8]:
        loc = (
            f"c{row.get('cell')}.{row.get('phase', '')}."
            f"{row.get('species', row.get('kind'))}"
        )
        print(f"  {loc:28s} F_hat={row['F_hat']:+.4f}", flush=True)
    return out


def main() -> None:
    arms = [
        _run("ps_o8", {}, 8),
        _run("gate_after_yco_o8", dict(GATE_AFTER), 8),
        _run("hot_restart_o8", dict(HOT), 8),
        _run("gate_plus_hot_o8", {**GATE_AFTER, **HOT}, 8),
    ]
    out = data_path("phase2_gate_hot_topk_residuals_o8.json")
    write_json(
        out,
        {
            "note": "top residuals after gate/hot stack; who blocks fully converge; §5.17",
            "arms": arms,
        },
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
