"""写出阶段 F 床–自由空间界面恒等式。不跑 Newton，不接生产 simulate_freeboard。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_bed_freeboard import build_stage_f_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_f_bed_freeboard.json"


def main() -> None:
    payload = build_stage_f_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    closed = payload["recycle_closed_coarse"]
    print("wirsum_eject", closed["wirsum_1997_exact"]["m_eject_kg_s"])
    print("hamel_eject", closed["hamel_extended"]["m_eject_kg_s"])
    print("eject_delta", payload["eject_hamel_minus_wirsum_kg_s"])
    print("fate_coarse", closed["hamel_extended"]["fate"])
    print("fate_fine", payload["recycle_open_fine_hamel"]["fate"])
    print("recycle_fine", payload["recycle_open_fine_hamel"]["m_recycle_kg_s"])


if __name__ == "__main__":
    main()
