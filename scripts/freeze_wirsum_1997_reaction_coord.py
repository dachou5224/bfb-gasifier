"""写出阶段 F3 离开炭 R1 协调。不跑全炉。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_reaction_coord import build_stage_f3_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_f3_reaction_coord.json"


def main() -> None:
    payload = build_stage_f3_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("leave", payload["m_leave_char_kg_s"])
    print("newton", payload["newton"])
    print("e1_root", payload["root"]["e1_mol_s"], "e1_inlet", payload["inlet_extent_e1_mol_s"])
    print("C_gap", payload["root"]["carbon_gas_minus_e1"])


if __name__ == "__main__":
    main()
