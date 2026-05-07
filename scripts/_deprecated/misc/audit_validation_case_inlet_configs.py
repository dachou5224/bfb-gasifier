#!/usr/bin/env python3
"""审计 validation_cases.json 中各 case 的 inlet / recycle / topology 配置完备度。"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tests.validation_case_utils import extract_case_inlet_topology_hints, load_validation_json_root


def _fmt_num(val: object, nd: int = 3) -> str:
    if isinstance(val, (int, float)):
        return f"{float(val):.{nd}f}"
    return "-"


def _readiness(h: dict[str, object]) -> str:
    secondary_ready = h["secondary_injection_xi"] is not None and h["secondary_air_Nm3_h"] is not None
    geometry_ready = h["reactor_height_m"] is not None and h["bed_diameter_m"] is not None and h["freeboard_diameter_m"] is not None
    recycle_ready = bool(h["recirculation"]) and bool(h["recycle_topology"] or h["cyclone_present"])
    if secondary_ready and geometry_ready:
        return "builder-ready"
    if geometry_ready and (h["secondary_injection_xi"] is not None or h["recirculation"]):
        return "hint-ready"
    if recycle_ready or geometry_ready:
        return "partial"
    return "sparse"


def main() -> int:
    data = load_validation_json_root()
    print("=" * 176)
    print("Validation case inlet/topology readiness audit")
    print("=" * 176)
    print(
        f"{'case':<28} {'plant':<24} {'P[MPa]':>7} {'sec_xi':>7} {'sec_air':>10} {'sec_O2':>10} "
        f"{'steam':>8} {'recycle':>8} {'cyclone':>8} {'bedD':>7} {'fbD':>7} {'H':>7} {'bedH':>7} {'status':>14}"
    )
    print("-" * 176)

    for key, node in data.items():
        if not key.startswith("CASE_"):
            continue
        hints = extract_case_inlet_topology_hints(node)
        op = node.get("inputs", {}).get("operating_conditions", {})
        print(
            f"{key:<28} {str(hints['plant'])[:24]:<24} "
            f"{_fmt_num(op.get('pressure_MPa'), 2):>7} "
            f"{_fmt_num(hints['secondary_injection_xi'], 3):>7} "
            f"{_fmt_num(hints['secondary_air_Nm3_h'], 1):>10} "
            f"{_fmt_num(hints['secondary_air_O2_mol_s'], 3):>10} "
            f"{_fmt_num(hints['steam_feed_kg_h'], 1):>8} "
            f"{str(bool(hints['recirculation'])):>8} "
            f"{str(bool(hints['cyclone_present'])):>8} "
            f"{_fmt_num(hints['bed_diameter_m'], 3):>7} "
            f"{_fmt_num(hints['freeboard_diameter_m'], 3):>7} "
            f"{_fmt_num(hints['reactor_height_m'], 2):>7} "
            f"{_fmt_num(hints['bed_height_m'], 2):>7} "
            f"{_readiness(hints):>14}"
        )

    print("-" * 176)
    print("note:")
    print("  - builder-ready: geometry + secondary injection location/flow 都已结构化，可直接驱动 freeboard inlet config。")
    print("  - hint-ready: 已有 geometry 和 secondary/recycle hint，但 secondary flow 或其他关键量仍缺。")
    print("  - partial: 只有部分 reactor/recycle/topology 提示。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
