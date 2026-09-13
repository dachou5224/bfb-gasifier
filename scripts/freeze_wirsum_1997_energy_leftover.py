"""写出 holdup Startwert 之后的 Eq. 2.7 剩余桶。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_energy_leftover import build_energy_leftover_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_energy_leftover.json"


def main() -> None:
    payload = build_energy_leftover_payload()
    slim = {
        "source": payload["source"],
        "case": payload["case"],
        "solved": payload["solved"],
        "name": payload["name"],
        "not": payload["not"],
        "units": payload["units"],
        "note": payload["note"],
        "eq27": payload["eq27"],
        "bed0_dominant_bucket": payload["bed0_dominant_bucket"],
        "per_cell": [
            {
                "cell": row["cell"],
                "residual_W": row["residual_W"],
                "dominant_bucket": row["dominant_bucket"],
                "sum_matches_residual": row["sum_matches_residual"],
                "holdup_zero": row["holdup_zero"],
                "vm_in_not_in_energy_out": row["vm_in_not_in_energy_out"],
                "buckets_W": row["buckets_W"],
                "abs_share": row["abs_share"],
            }
            for row in payload["per_cell"]
        ],
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    b0 = slim["per_cell"][0]
    print("bed0_E", b0["residual_W"], "dom", b0["dominant_bucket"])
    print("buckets", {k: round(v / 1e6, 3) for k, v in b0["buckets_W"].items()})


if __name__ == "__main__":
    main()
