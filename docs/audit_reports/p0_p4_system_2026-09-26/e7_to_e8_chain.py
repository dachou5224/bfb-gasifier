"""从快照里的 E7 x 走现有 E8 延拓。600 s 到点即停。不手改未知量。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f_e8 import solve_named_e8_from_e7
from src.solvers.wirsum_1997_named_f_endpoint_io import load_named_e8_snapshot

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/p0_p4_system_2026-09-26/e7_to_e8_chain.json"


def main() -> None:
    snap = load_named_e8_snapshot()
    e7 = snap["e7"]
    t0 = time.monotonic()
    out = solve_named_e8_from_e7(e7=e7, deadline_s=600.0, max_iter=8)
    elapsed = time.monotonic() - t0
    hist = [
        {
            "hydro_scale": row["hydro_scale"],
            "fate_scale": row["fate_scale"],
            "stage_passed": bool(row["stage_passed"]),
        }
        for row in out.get("continuation", [])
    ]
    cert = out.get("endpoint_cert", {})
    summary = {
        "elapsed_s": elapsed,
        "budget_s": 600.0,
        "timed_out": bool(out.get("timed_out")),
        "endpoint_certified": bool(out.get("endpoint_certified")),
        "T_K": out.get("T_K"),
        "gas_mol_s_l2": cert.get("gas_mol_s_l2"),
        "energy_W": cert.get("energy_W"),
        "n_scales": len(hist),
        "last_scale": hist[-1] if hist else None,
        "continuation": hist,
    }
    _OUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("elapsed_s", "timed_out", "endpoint_certified", "last_scale", "n_scales")}))


if __name__ == "__main__":
    main()
