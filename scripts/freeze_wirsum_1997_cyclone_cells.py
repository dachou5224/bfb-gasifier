"""写出 Hamel Bild 7.7 旋风 cell 核。不跑全炉。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_cyclone_cells import build_stage_f_cyclone_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_f_cyclone_cells.json"


def main() -> None:
    payload = build_stage_f_cyclone_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("leave", payload["m_leave_char_kg_s"])
    print("closed", payload["recycle_closed"])
    print("open_dF", payload["recycle_open"]["dF_bed_dMret"], "F_bed", payload["recycle_open"]["F_bed"])


if __name__ == "__main__":
    main()
