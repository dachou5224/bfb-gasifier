"""写出阶段 C 约化输运恒等式。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.core.solid_transport_reduced import (
    engineering_closure_jump_kg,
    mixed_denominator_gap,
    species_own_identity,
)

OUT = ROOT / "data" / "wirsum_alignment_stage_c_solid_transport.json"


def main() -> None:
    M = np.array([[4.08, 1.13], [4.79, 1.32], [4.63, 1.28]], dtype=np.float64)
    M_closure = np.array([5.204, 191.261, 184.872], dtype=np.float64)
    K_up = np.array([0.079, 0.079, 0.032], dtype=np.float64)
    feed = np.array([[0.322, 0.089], [0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    payload = {
        "source": "Wirsum alignment stage C solid transport",
        "solved": False,
        "engineering_15pct_jump": engineering_closure_jump_kg(187.958),
        "species_own": species_own_identity(M, K_up, feed),
        "engineering_mixed_denominator": mixed_denominator_gap(M, M_closure, K_up, feed),
        "note": "默认生产路径仍是 engineering_bottom_threshold；species_own 为 opt-in / strict。",
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("jump_ratio", payload["engineering_15pct_jump"]["jump_ratio"])
    print("mixed_gap", payload["engineering_mixed_denominator"]["max_abs_gap_kg_s"])
    print("species_km_gap", payload["species_own"]["km_vs_down_gap_kg_s"])


if __name__ == "__main__":
    main()
