"""写出同相 R12 Startwert A/B（reduced 3-cell x0）。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_startwert_r12 import build_startwert_r12_ab_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_startwert_r12.json"


def main() -> None:
    payload = build_startwert_r12_ab_payload()
    slim = {
        "source": payload["source"],
        "case": payload["case"],
        "solved": payload["solved"],
        "name": payload["name"],
        "not": payload["not"],
        "units": payload["units"],
        "note": payload["note"],
        "applied": {
            "xi12_total_mol_s": payload["applied"]["xi12_total_mol_s"],
            "cells_touched": payload["applied"]["cells_touched"],
            "per_cell": [
                {
                    "cell": row["cell"],
                    "xi12_mol_s": row["xi12_mol_s"],
                    "dense_xi12": row["dense"]["xi12_mol_s"],
                    "bubble_xi12": row["bubble"]["xi12_mol_s"],
                }
                for row in payload["applied"]["per_cell"]
            ],
        },
        "before": {
            "max_xi12_possible_mol_s": payload["before"]["max_xi12_possible_mol_s"],
            "max_xi5_possible_mol_s": payload["before"]["max_xi5_possible_mol_s"],
            "max_xi6_possible_mol_s": payload["before"]["max_xi6_possible_mol_s"],
            "max_r12_extent_mol_s": payload["before"]["max_r12_extent_mol_s"],
            "residual_l2": payload["before"]["residual_l2"],
            "residual_linf": payload["before"]["residual_linf"],
            "overlap": payload["before"]["overlap"],
        },
        "after": {
            "max_xi12_possible_mol_s": payload["after"]["max_xi12_possible_mol_s"],
            "max_xi5_possible_mol_s": payload["after"]["max_xi5_possible_mol_s"],
            "max_xi6_possible_mol_s": payload["after"]["max_xi6_possible_mol_s"],
            "max_r12_extent_mol_s": payload["after"]["max_r12_extent_mol_s"],
            "residual_l2": payload["after"]["residual_l2"],
            "residual_linf": payload["after"]["residual_linf"],
            "overlap": payload["after"]["overlap"],
        },
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("before_xi12_holdup", slim["before"]["max_xi12_possible_mol_s"])
    print("after_xi12_holdup", slim["after"]["max_xi12_possible_mol_s"])
    print("before_r12_extent", slim["before"]["max_r12_extent_mol_s"])
    print("after_r12_extent", slim["after"]["max_r12_extent_mol_s"])
    print("before_l2", slim["before"]["residual_l2"])
    print("after_l2", slim["after"]["residual_l2"])
    print("left_xi5", slim["after"]["max_xi5_possible_mol_s"])


if __name__ == "__main__":
    main()
