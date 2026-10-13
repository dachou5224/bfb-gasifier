"""写出阶段 A 两条基线的上下文契约。只初始化，不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_baseline_contract import build_stage_a_baseline_payload

OUT = ROOT / "data" / "wirsum_alignment_stage_a_baselines.json"


def main() -> None:
    payload = build_stage_a_baseline_payload()
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    for name, snap in payload["paths"].items():
        print(
            name,
            "thesis",
            snap["knobs"]["thesis_mode"],
            "n_unknowns",
            snap["n_unknowns"],
            "vm_cache",
            snap["vm_cache_valid"],
        )


if __name__ == "__main__":
    main()
