"""写出 P4 Vorabrechnung outer map A/B。

输出 ``data/outer_iteration_history.json``。不跑 full LU。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_outer_diagnostics import build_outer_map_ab_payload

OUT = ROOT / "data" / "outer_iteration_history.json"


def main() -> None:
    payload = build_outer_map_ab_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
