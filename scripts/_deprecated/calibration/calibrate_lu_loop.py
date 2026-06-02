"""LU 工况简化标定循环（Phase-1）。

扫描少量高杠杆参数，按 KPI 相对误差打分并输出最优组合。
仅用于离线标定，不在主路径自动启用。
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import src.kinetics.char_reactions as cr
import src.kinetics.gas_reactions as gr
from src.core.reactor import Reactor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_reactor_config,
    load_validation_case_node,
)


@dataclass
class RunRow:
    dense_frac: float
    heat_loss: float
    r5_scale: float
    r4_scale: float
    T_exit: float
    Xc: float
    CO: float
    CO2: float
    H2: float
    CH4: float
    score: float


def rel_err(sim: float, ref: float) -> float:
    if ref <= 0.0:
        return 0.0
    return abs(sim - ref) / ref


def main() -> int:
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref["exit_gas_dry_mol_frac"]
    t_ref = float(ref["exit_temperature_K"])
    xc_ref = float(ref["carbon_conversion_pct"]) / 100.0

    # 记录原始常数
    base_r5 = gr.R5_suspension_k_T05
    base_r4 = cr.R4_kf_k0

    dense_grid = [0.20, 0.35, 0.50, 0.65]
    hl_grid = [0.10, 0.12]
    r5_grid = [1.0, 0.3, 0.1]
    r4_grid = [1.0, 3.0]

    rows: list[RunRow] = []

    print("=" * 88)
    print("LU calibration loop: dense_frac / heat_loss / R5_scale / R4_scale")
    print("=" * 88)

    n_total = len(dense_grid) * len(hl_grid) * len(r5_grid) * len(r4_grid)
    n_done = 0

    try:
        for dense_frac in dense_grid:
            for heat_loss in hl_grid:
                for r5_scale in r5_grid:
                    for r4_scale in r4_grid:
                        n_done += 1
                        gr.R5_suspension_k_T05 = base_r5 * r5_scale
                        cr.R4_kf_k0 = base_r4 * r4_scale

                        cfg = build_phase1_htw_lu_reactor_config()
                        cfg.gas_inlet_dense_frac = dense_frac
                        cfg.heat_loss_frac = heat_loss

                        try:
                            res = Reactor(cfg).solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
                        except Exception as e:
                            print(
                                f"[{n_done:03d}/{n_total}] dense={dense_frac:.2f} hl={heat_loss:.2f} "
                                f"r5={r5_scale:.2f} r4={r4_scale:.1f} | FAILED: {type(e).__name__}: {e}"
                            )
                            continue
                        y = res["exit_gas_dry"]
                        t = float(res["T_profile"][-1])
                        xc = float(res["carbon_conv"])
                        co = float(y.get("CO", 0.0))
                        co2 = float(y.get("CO2", 0.0))
                        h2 = float(y.get("H2", 0.0))
                        ch4 = float(y.get("CH4", 0.0))

                        # 加权评分：优先温度+碳转化，同时约束四主组分
                        e_t = rel_err(t, t_ref)
                        e_x = rel_err(xc, xc_ref)
                        e_co = rel_err(co, float(ref_dry["CO"]))
                        e_co2 = rel_err(co2, float(ref_dry["CO2"]))
                        e_h2 = rel_err(h2, float(ref_dry["H2"]))
                        e_ch4 = rel_err(ch4, float(ref_dry["CH4"]))
                        score = (
                            1.2 * e_t
                            + 1.2 * e_x
                            + 1.4 * e_co
                            + 1.0 * e_co2
                            + 1.0 * e_h2
                            + 0.8 * e_ch4
                        )

                        rows.append(
                            RunRow(
                                dense_frac=dense_frac,
                                heat_loss=heat_loss,
                                r5_scale=r5_scale,
                                r4_scale=r4_scale,
                                T_exit=t,
                                Xc=xc,
                                CO=co,
                                CO2=co2,
                                H2=h2,
                                CH4=ch4,
                                score=score,
                            )
                        )

                        print(
                            f"[{n_done:03d}/{n_total}] "
                            f"dense={dense_frac:.2f} hl={heat_loss:.2f} r5={r5_scale:.2f} r4={r4_scale:.1f} | "
                            f"T={t:.0f} Xc={xc:.3f} CO={co*100:.2f}% CO2={co2*100:.2f}% H2={h2*100:.2f}% CH4={ch4*100:.2f}% | "
                            f"score={score:.3f}"
                        )
    finally:
        gr.R5_suspension_k_T05 = base_r5
        cr.R4_kf_k0 = base_r4

    rows.sort(key=lambda x: x.score)

    print("\n" + "=" * 88)
    print("Top 10 parameter sets")
    print("=" * 88)
    for i, r in enumerate(rows[:10], start=1):
        print(
            f"#{i:02d} score={r.score:.3f} | dense={r.dense_frac:.2f} hl={r.heat_loss:.2f} r5={r.r5_scale:.2f} r4={r.r4_scale:.1f} | "
            f"T={r.T_exit:.0f} Xc={r.Xc:.3f} CO={r.CO*100:.2f}% CO2={r.CO2*100:.2f}% H2={r.H2*100:.2f}% CH4={r.CH4*100:.2f}%"
        )

    best = rows[0]
    print("\nBEST:")
    print(best)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
