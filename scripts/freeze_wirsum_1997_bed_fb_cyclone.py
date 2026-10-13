"""写出阶段 F5 床–自由空间–旋风交接。不跑全炉。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_bed_fb_cyclone import build_stage_f5_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_f5_bed_fb_cyclone.json"


def main() -> None:
    payload = build_stage_f5_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("fine_off", payload["fine_recycle_closed"])
    print("fine_on_rez", payload["fine_recycle_open"]["m_rez_kg_s"], "F_bed", payload["fine_recycle_open"]["F_bed"])
    print("coarse_leave", payload["coarse_recycle_closed"]["m_leave_char_kg_s"])


if __name__ == "__main__":
    main()
