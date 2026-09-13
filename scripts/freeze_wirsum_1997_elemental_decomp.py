"""写出 P5 源项元素分解（reduced 3-cell x0）。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_elemental_decomp import build_p5_elemental_decomp_payload

OUT = ROOT / "data" / "wirsum_alignment_p5_elemental_decomp.json"


def main() -> None:
    payload = build_p5_elemental_decomp_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("max_r12", payload["max_r12_extent_mol_s"])
    print("max_HO_after_r12", payload["max_HO_after_r12"])


if __name__ == "__main__":
    main()
