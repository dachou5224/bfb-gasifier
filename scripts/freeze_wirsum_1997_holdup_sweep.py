"""写出一次扫描 vs 轴向接力 holdup Startwert。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_holdup_startwert import build_holdup_sweep_ab_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_holdup_sweep.json"


def _slim_leftover(row: dict) -> dict:
    return {
        "cell": row["cell"],
        "residual_W": row["residual_W"],
        "dominant_bucket": row["dominant_bucket"],
        "holdup_zero": row["holdup_zero"],
        "buckets_W": row["buckets_W"],
    }


def main() -> None:
    payload = build_holdup_sweep_ab_payload()
    slim = {
        "source": payload["source"],
        "case": payload["case"],
        "solved": payload["solved"],
        "name": payload["name"],
        "not": payload["not"],
        "units": payload["units"],
        "note": payload["note"],
        "applied": {
            "cells_touched": payload["applied"]["cells_touched"],
            "per_cell": [
                {
                    "cell": row.get("cell"),
                    "xi_M_kg": row.get("xi_M_kg"),
                    "char_M_kg": (row.get("per_comp") or {}).get("char", {}).get("M_after_kg"),
                    "vm_M_kg": (row.get("per_comp") or {}).get("vm", {}).get("M_after_kg"),
                }
                for row in payload["applied"]["per_cell"]
            ],
        },
        "one_shot": {
            "per_cell_residual_W": payload["one_shot"]["eq27"]["per_cell_residual_W"],
            "holdup_zero": payload["one_shot"]["eq27"]["holdup_zero"],
            "leftover": [_slim_leftover(row) for row in payload["one_shot"]["leftover"]],
        },
        "axial_sweep": {
            "per_cell_residual_W": payload["axial_sweep"]["eq27"]["per_cell_residual_W"],
            "holdup_zero": payload["axial_sweep"]["eq27"]["holdup_zero"],
            "m_holdup_kg": payload["axial_sweep"]["eq27"]["m_holdup_kg"],
            "leftover": [_slim_leftover(row) for row in payload["axial_sweep"]["leftover"]],
        },
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("one_shot_E", slim["one_shot"]["per_cell_residual_W"])
    print("sweep_E", slim["axial_sweep"]["per_cell_residual_W"])
    print("sweep_holdup_zero", slim["axial_sweep"]["holdup_zero"])


if __name__ == "__main__":
    main()
