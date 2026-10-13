"""写出阶段 F 后半协调核。不跑全炉、不接生产 simulate_freeboard。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_bed_freeboard_coord import build_stage_f2_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_f2_bed_freeboard_coord.json"


def main() -> None:
    payload = build_stage_f2_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("T_Q0", payload["inner_Q0"]["T"], "T_Qw", payload["inner_Qwall"]["T"])
    print("d50_Q0", payload["outer_Q0"]["d50_m"], "d50_Qw", payload["inner_Qwall"]["d50_m"])
    print("eta_leave", payload["outer_Q0"]["eta"], "eta_at_old_cut", payload["inner_Qwall"]["eta_at_Q0_cut"])
    print("const_minus_lapple", payload["constant_eta_minus_lapple_Q0"])


if __name__ == "__main__":
    main()
