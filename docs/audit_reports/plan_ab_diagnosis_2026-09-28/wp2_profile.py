"""减半初值上冻结 q 修正 10 s，数真实残差调用。不改门限。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp2_profile.json"


def main() -> None:
    x = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"].copy()
    pos = x[:-2] > 0.0
    x[:-2][pos] *= 0.5
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    q = outer_P_named(x, nq)
    calls = {"F": 0}
    real = coupled_residual_named

    def counted(z, *args, **kwargs):
        calls["F"] += 1
        return real(z, *args, **kwargs)

    import src.solvers.wirsum_1997_named_f as named
    import src.solvers.wirsum_1997_named_f_e8 as e8

    named.coupled_residual_named = counted
    e8.coupled_residual_named = counted
    t0 = time.perf_counter()
    try:
        _e8_log_scaled_correct(q, x, deadline_s=10.0, max_nfev=10**9, ftol=1.0e-12)
        status = "returned"
    except TimeoutError:
        status = "timeout"
    finally:
        named.coupled_residual_named = real
        e8.coupled_residual_named = real
    elapsed = time.perf_counter() - t0
    row = {
        "status": status,
        "elapsed_s": elapsed,
        "F_calls": calls["F"],
        "F_per_s": calls["F"] / max(elapsed, 1.0e-9),
        "note": "Same half-state frozen-q corrector as the 10 s diagnosis. Char balance is assembled once per F.",
    }
    _OUT.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
