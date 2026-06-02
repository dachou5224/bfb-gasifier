#!/usr/bin/env python3
"""LU 温度剖面审计：将当前模型结果与 result-profile.png 提取点做近似对比。

说明
----
- 该数据源来自 `docs/charts/result-profile.png` 的手工数字化，仅适合作为
  reactor-level 温度 profile 的软验证。
- 由于原图是位图且左图浓度曲线被图例/注释遮挡，本脚本当前只审计温度 profile。
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_reactor_config


DATA_PATH = REPO_ROOT / "data" / "validation_profile_result_profile.json"


def _load_profile_data() -> dict[str, Any]:
    with open(DATA_PATH, encoding="utf-8") as f:
        return json.load(f)


def _interp_profile(x_query: np.ndarray, x_data: np.ndarray, y_data: np.ndarray) -> np.ndarray:
    return np.interp(x_query, x_data, y_data, left=y_data[0], right=y_data[-1])


def run_audit() -> dict[str, Any]:
    cfg = build_phase1_htw_lu_reactor_config()
    result = Reactor(cfg).solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    xi_model = np.array([float(c.geo.h_center / cfg.H_bed) for c in Reactor(cfg).cells], dtype=float)
    T_model_C = np.array([float(t - 273.15) for t in result["T_profile"]], dtype=float)

    data = _load_profile_data()
    exp_points = data["temperature_profile_experiment"]
    sim_points = data["temperature_profile_simulation"]

    xi_exp = np.array([float(p["xi"]) for p in exp_points], dtype=float)
    T_exp = np.array([float(p["temperature_C"]) for p in exp_points], dtype=float)
    xi_png_sim = np.array([float(p["xi"]) for p in sim_points], dtype=float)
    T_png_sim = np.array([float(p["temperature_C"]) for p in sim_points], dtype=float)

    T_model_at_exp = _interp_profile(xi_exp, xi_model, T_model_C)
    T_model_at_png_sim = _interp_profile(xi_png_sim, xi_model, T_model_C)

    err_exp = T_model_at_exp - T_exp
    err_png_sim = T_model_at_png_sim - T_png_sim

    rows_exp: list[dict[str, float]] = []
    for xi, t_ref, t_mod, err in zip(xi_exp, T_exp, T_model_at_exp, err_exp):
        rows_exp.append(
            {
                "xi": float(xi),
                "T_ref_C": float(t_ref),
                "T_model_C": float(t_mod),
                "error_C": float(err),
            }
        )

    rows_png_sim: list[dict[str, float]] = []
    for xi, t_ref, t_mod, err in zip(xi_png_sim, T_png_sim, T_model_at_png_sim, err_png_sim):
        rows_png_sim.append(
            {
                "xi": float(xi),
                "T_png_sim_C": float(t_ref),
                "T_model_C": float(t_mod),
                "error_C": float(err),
            }
        )

    mae_exp = float(np.mean(np.abs(err_exp)))
    rmse_exp = float(math.sqrt(float(np.mean(err_exp**2))))
    mae_png_sim = float(np.mean(np.abs(err_png_sim)))
    rmse_png_sim = float(math.sqrt(float(np.mean(err_png_sim**2))))

    return {
        "source_image": data["source_image"],
        "confidence": data["confidence"],
        "solver": PHASE1_HTW_LU_SOLVE_KWARGS.get("solver"),
        "solve_kwargs": dict(PHASE1_HTW_LU_SOLVE_KWARGS),
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "model_profile_C": [
            {"xi": float(xi), "T_model_C": float(t)}
            for xi, t in zip(xi_model, T_model_C)
        ],
        "vs_experiment": {
            "n_points": int(len(rows_exp)),
            "mae_C": mae_exp,
            "rmse_C": rmse_exp,
            "rows": rows_exp,
        },
        "vs_png_simulation": {
            "n_points": int(len(rows_png_sim)),
            "mae_C": mae_png_sim,
            "rmse_C": rmse_png_sim,
            "rows": rows_png_sim,
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="审计当前 LU 温度剖面 vs PNG 提取点")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="当 vs_experiment 的 MAE > 120°C 或 RMSE > 140°C 时返回非零",
    )
    args = ap.parse_args()

    out = run_audit()

    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        print("LU temperature profile audit (digitized from result-profile.png)")
        print(
            f"  converged={out['converged']} n_iter={out['n_iter']} "
            f"solver={out['solver']} confidence={out['confidence']['temperature_profile']}"
        )
        print(
            f"  vs experiment: MAE={out['vs_experiment']['mae_C']:.1f}°C "
            f"RMSE={out['vs_experiment']['rmse_C']:.1f}°C"
        )
        print(
            f"  vs png simulation: MAE={out['vs_png_simulation']['mae_C']:.1f}°C "
            f"RMSE={out['vs_png_simulation']['rmse_C']:.1f}°C"
        )
        print("  xi    T_ref(°C)  T_model(°C)  err(°C)")
        for row in out["vs_experiment"]["rows"]:
            print(
                f"  {row['xi']:>4.2f}   {row['T_ref_C']:>8.1f}    "
                f"{row['T_model_C']:>10.1f}  {row['error_C']:>8.1f}"
            )

    if args.strict:
        if out["vs_experiment"]["mae_C"] > 120.0 or out["vs_experiment"]["rmse_C"] > 140.0:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
