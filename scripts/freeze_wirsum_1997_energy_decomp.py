"""写出 Eq. 2.7 能量分项（reduced 3-cell，R12 Startwert 前/后）。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_energy_decomp import build_eq27_energy_decomp_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_energy_decomp.json"


def _slim_cell(row: dict) -> dict:
    return {
        "cell": row["cell"],
        "T_K": row["T_K"],
        "heat_loss_mode": row["heat_loss_mode"],
        "H_in_W": row["H_in_W"],
        "H_in_parts_W": row["H_in_parts_W"],
        "H_gas_out_W": row["H_gas_out_W"],
        "H_solid_out_W": row["H_solid_out_W"],
        "Q_loss_W": row["Q_loss_W"],
        "residual_W": row["residual_W"],
        "decomp_matches_live": row["decomp_matches_live"],
        "dominant_H_in_term": row["dominant_H_in_term"],
        "H_zu_by_comp_W": row["H_zu_by_comp_W"],
        "m_in_kg_s": row["m_in_kg_s"],
        "m_energy_out_kg_s": row["m_energy_out_kg_s"],
        "m_holdup_kg": row["m_holdup_kg"],
        "K_auf_max_1_per_s": row["K_auf_max_1_per_s"],
        "holdup_zero": row["holdup_zero"],
        "H_solid_out_zero": row["H_solid_out_zero"],
        "vm_in_not_in_energy_out": row["vm_in_not_in_energy_out"],
        "solid_state_model": row["solid_state_model"],
    }


def main() -> None:
    payload = build_eq27_energy_decomp_payload()
    slim = {
        "source": payload["source"],
        "case": payload["case"],
        "solved": payload["solved"],
        "name": payload["name"],
        "not": payload["not"],
        "units": payload["units"],
        "note": payload["note"],
        "before": {
            "bed0_residual_W": payload["before"]["bed0_residual_W"],
            "bed0_dominant_H_in_term": payload["before"]["bed0_dominant_H_in_term"],
            "bed0_holdup_zero": payload["before"]["bed0_holdup_zero"],
            "per_cell": [_slim_cell(row) for row in payload["before"]["per_cell"]],
        },
        "after": {
            "bed0_residual_W": payload["after"]["bed0_residual_W"],
            "bed0_dominant_H_in_term": payload["after"]["bed0_dominant_H_in_term"],
            "bed0_holdup_zero": payload["after"]["bed0_holdup_zero"],
            "per_cell": [_slim_cell(row) for row in payload["after"]["per_cell"]],
        },
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    b0 = slim["after"]["per_cell"][0]
    print("bed0_E", b0["residual_W"])
    print("H_in", b0["H_in_W"], "H_sout", b0["H_solid_out_W"], "Q", b0["Q_loss_W"])
    print("holdup_zero", b0["holdup_zero"], "H_zu", b0["H_in_parts_W"]["H_solid_zu"])


if __name__ == "__main__":
    main()
