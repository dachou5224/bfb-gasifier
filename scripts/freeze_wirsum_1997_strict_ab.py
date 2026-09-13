"""写出 P2.1 Wirsum strict Newton 四组 A/B 表。

输出 ``data/wirsum_1997_strict_ab.json``。不跑 full LU。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_strict import build_reduced_holdup_ab_payload

OUT = ROOT / "data" / "wirsum_1997_strict_ab.json"


def main() -> None:
    payload = build_reduced_holdup_ab_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
