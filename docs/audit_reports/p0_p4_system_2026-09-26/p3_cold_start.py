"""冷启动。初值来自入口和冻结 q，不读 E8/P3 端点。预算 600 s。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import edges_from_x
from src.solvers.wirsum_1997_named_f import (
    coupled_residual_named,
    default_x0_named,
    outer_P_named,
)
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/p0_p4_system_2026-09-26/p3_cold_start.json"
_BUDGET_S = 600.0


def _initial() -> np.ndarray:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = default_x0_named(nq)
    if float(x[42]) <= 0.0:
        x[42] = 1.2e-5
    if float(x[43]) <= 0.0:
        x[43] = 1.8e-5
    q = outer_P_named(x, nq)
    ed = edges_from_x(x[:39], q.energy.species)
    tau = np.maximum(np.asarray(q.hydro.tau_s, dtype=np.float64), 1.0e-6)
    x[46:48] = np.maximum(np.asarray(ed["m_eject"], dtype=np.float64) * tau, 1.0e-12)
    return x


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = _initial()
    t0 = time.monotonic()
    status = "timeout"
    blocks = None
    worst = None
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
        "name": "cold_start",
        "reads_endpoint": False,
        "status": status,
        "elapsed_s": elapsed,
        "budget_s": _BUDGET_S,
        "outer_refreshes": outers,
        "blocks": None if blocks is None else {
            k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in blocks.items()
        },
        "worst_row": worst,
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
    }
    _OUT.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(json.dumps({"status": status, "elapsed_s": elapsed, "outers": outers, "T_bed": row["T_bed"], "T_fb": row["T_fb"]}), flush=True)


if __name__ == "__main__":
    main()
