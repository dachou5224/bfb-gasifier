"""写出 P5 LU 重诊清单（reduced 3-cell x0）。

输出 ``data/wirsum_alignment_p5_lu_reentry.json``。不跑 Newton / full LU。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_lu_reentry import build_p5_lu_reentry_payload

OUT = ROOT / "data" / "wirsum_alignment_p5_lu_reentry.json"


def main() -> None:
    payload = build_p5_lu_reentry_payload(include_jacobian=False)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("open_items:", payload["open_items"])


if __name__ == "__main__":
    main()
