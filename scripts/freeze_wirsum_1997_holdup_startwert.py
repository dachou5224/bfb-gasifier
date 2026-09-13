"""写出 ``M = m_in / K`` Startwert A/B。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_holdup_startwert import build_holdup_startwert_ab_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_holdup_startwert.json"


def _slim_cell(row: dict) -> dict:
    return {
        "cell": row["cell"],
        "T_K": row["T_K"],
        "residual_W": row["residual_W"],
        "H_in_W": row["H_in_W"],
        "H_solid_zu_W": row["H_in_parts_W"]["H_solid_zu"],
        "H_gas_out_W": row["H_gas_out_W"],
        "H_solid_out_W": row["H_solid_out_W"],
        "Q_loss_W": row["Q_loss_W"],
        "m_holdup_kg": row["m_holdup_kg"],
        "m_energy_out_kg_s": row["m_energy_out_kg_s"],
        "H_zu_by_comp_W": row["H_zu_by_comp_W"],
        "H_solid_out_by_comp_W": row["H_solid_out_by_comp_W"],
        "holdup_zero": row["holdup_zero"],
        "H_solid_out_zero": row["H_solid_out_zero"],
        "vm_in_not_in_energy_out": row["vm_in_not_in_energy_out"],
    }


def main() -> None:
    payload = build_holdup_startwert_ab_payload()
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
                    "skipped": row.get("skipped"),
                    "xi_M_kg": row.get("xi_M_kg"),
                    "char_M_kg": (row.get("per_comp") or {}).get("char", {}).get("M_after_kg"),
                    "ash_M_kg": (row.get("per_comp") or {}).get("ash", {}).get("M_after_kg"),
                    "vm_M_kg": (row.get("per_comp") or {}).get("vm", {}).get("M_after_kg"),
                }
                for row in payload["applied"]["per_cell"]
            ],
        },
        "before": {
            "bed0_residual_W": payload["before"]["bed0_residual_W"],
            "per_cell": [_slim_cell(row) for row in payload["before"]["per_cell"]],
        },
        "after": {
            "bed0_residual_W": payload["after"]["bed0_residual_W"],
            "per_cell": [_slim_cell(row) for row in payload["after"]["per_cell"]],
        },
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("before_E", slim["before"]["bed0_residual_W"])
    print("after_E", slim["after"]["bed0_residual_W"])
    print("after_H_sout", slim["after"]["per_cell"][0]["H_solid_out_W"])
    print("after_holdup", slim["after"]["per_cell"][0]["m_holdup_kg"])


if __name__ == "__main__":
    main()
