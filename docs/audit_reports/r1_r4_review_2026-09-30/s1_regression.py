"""温度 +50 K 完整回归。过门后再测半尺度。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

_ROOT = Path(__file__).resolve().parents[3]
_DIR = Path(__file__).resolve().parent
_BASE = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]


def _run(name: str, x0: np.ndarray) -> dict:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    # 120 s 时温度 +50 K 停在焦油门约 1.08 倍。门不放宽。
    # 同一控制器在 142 s 过温度 +50 K，半尺度约 290 s。
    out = solve_named_controlled(nq, x0, deadline_s=320.0)
    x = np.asarray(out["x"], dtype=np.float64)
    row = {
        "name": name,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "stats": out["stats"],
        "gate_ratio": out["gate_ratio"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "x": x.tolist(),
        "F": np.asarray(out["F"], dtype=np.float64).tolist(),
    }
    print(name, row["passed"], row["reason"], round(row["elapsed_s"], 1), row["stats"].get("accepted_alpha"), flush=True)
    return row


def main() -> None:
    hot = np.asarray(_BASE, dtype=np.float64).copy()
    hot[-2:] += 50.0
    first = _run("T_plus_50K", hot)
    rows = [first]
    if first["passed"]:
        half = np.asarray(_BASE, dtype=np.float64).copy()
        pos = half[:-2] > 0.0
        half[:-2][pos] *= 0.5
        rows.append(_run("positive_half", half))
    (_DIR / "s1_regression.json").write_text(json.dumps({"cases": rows}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
