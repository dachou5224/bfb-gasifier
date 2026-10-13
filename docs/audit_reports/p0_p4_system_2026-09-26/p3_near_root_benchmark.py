"""P3 近根扰动。每项 120 s。超时记失败，不以靠近原答案为通过。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q

_ROOT = Path(__file__).resolve().parents[3]
_X = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
_OUT = _ROOT / "docs/audit_reports/p0_p4_system_2026-09-26/p3_near_root_benchmark.json"
_BUDGET_S = 120.0


def _guesses() -> list[tuple[str, np.ndarray]]:
    base = np.asarray(_X, dtype=np.float64)
    hot = base.copy()
    hot[-2:] += 50.0
    cold = base.copy()
    cold[-2:] -= 50.0
    half = base.copy()
    pos = half[:-2] > 0.0
    half[:-2][pos] *= 0.5
    double = base.copy()
    double[:-2][pos] *= 2.0
    return [("T_plus_50K", hot), ("T_minus_50K", cold), ("positive_half", half), ("positive_double", double)]


def _one(name: str, x0: np.ndarray) -> dict:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    t0 = time.monotonic()
    status = "timeout"
    blocks = None
    worst = None
    x = np.asarray(x0, dtype=np.float64)
    x_out = x
    outers = 0
    try:
        while time.monotonic() - t0 < _BUDGET_S:
            remain = _BUDGET_S - (time.monotonic() - t0)
            q = outer_P_named(x, nq)
            x = _e8_log_scaled_correct(q, x, deadline_s=remain)
            outers += 1
            q = outer_P_named(x, nq)
            f = coupled_residual_named(x, q, rate_scale=1.0)
            blocks = named_endpoint_blocks(f)
            i = int(np.argmax(np.abs(f)))
            worst = {"index": i, "value": float(f[i])}
            x_out = x
            if bool(blocks["passed"]):
                status = "passed"
                break
            status = "failed"
    except TimeoutError:
        status = "timeout" if blocks is None else "failed"
    elapsed = time.monotonic() - t0
    if elapsed > _BUDGET_S and status != "passed":
        status = "timeout"
    row = {
        "name": name,
        "status": status,
        "elapsed_s": elapsed,
        "budget_s": _BUDGET_S,
        "outer_refreshes": outers,
        "blocks": None if blocks is None else {k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in blocks.items()},
        "worst_row": worst,
        "T_bed": float(x_out[-2]),
        "T_fb": float(x_out[-1]),
    }
    print(name, status, round(elapsed, 2), flush=True)
    return row


def main() -> None:
    rows = [_one(name, x0) for name, x0 in _guesses()]
    _OUT.write_text(json.dumps({"cases": rows}, indent=2), encoding="utf-8")
    print("wrote", _OUT)


if __name__ == "__main__":
    main()
