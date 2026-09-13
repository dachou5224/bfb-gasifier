"""写出 R12 Startwert 前后 packed 残差族分解。不跑 Newton。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.solvers.wirsum_1997_residual_decomp import build_residual_decomp_payload

OUT = ROOT / "data" / "wirsum_alignment_p6_residual_decomp.json"


def _slim_cell(row: dict) -> dict:
    return {
        "cell": row["cell"],
        "gas_residual_l2": row["gas_residual_l2"],
        "gas_linf": row["gas_linf"],
        "solid_residual_l2": row["solid_residual_l2"],
        "energy_residual_W": row["energy_residual_W"],
        "r5_extent_mol_s": row["r5"]["extent_mol_s"],
        "r5_looks_like_pair": row["r5_looks_like_pair"],
        "r12_extent_mol_s": row["r12"]["extent_mol_s"],
    }


def _slim_side(side: dict) -> dict:
    return {
        "residual_l2": side["residual_l2"],
        "residual_linf": side["residual_linf"],
        "linf_label": side["linf_label"],
        "linf_family": side["linf_family"],
        "linf_over_l2": side["linf_over_l2"],
        "family_l2": side["family_l2"],
        "dominant_family": side["dominant_family"],
        "max_r5_extent_mol_s": side["max_r5_extent_mol_s"],
        "max_r12_extent_mol_s": side["max_r12_extent_mol_s"],
        "per_cell": [_slim_cell(row) for row in side["per_cell"]],
    }


def main() -> None:
    payload = build_residual_decomp_payload()
    slim = {
        "source": payload["source"],
        "case": payload["case"],
        "solved": payload["solved"],
        "name": payload["name"],
        "not": payload["not"],
        "units": payload["units"],
        "note": payload["note"],
        "before": _slim_side(payload["before"]),
        "after": _slim_side(payload["after"]),
    }
    OUT.write_text(json.dumps(slim, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("before_dom", slim["before"]["dominant_family"], slim["before"]["linf_label"])
    print("after_dom", slim["after"]["dominant_family"], slim["after"]["linf_label"])
    print("after_energy_l2", slim["after"]["family_l2"]["energy"])
    print("after_gas_l2", slim["after"]["family_l2"]["gas"])
    print("after_r5", slim["after"]["max_r5_extent_mol_s"])


if __name__ == "__main__":
    main()
