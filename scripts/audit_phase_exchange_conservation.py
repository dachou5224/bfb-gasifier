#!/usr/bin/env python3
"""专项审计：bubble / dense 相间交换与两相守恒闭合。

覆盖两层：
1. 纯代数不变量：`N_ex` 在 dense / bubble 残差里严格反号，系统净贡献为零
2. 真实 LU solved cell：在实际 cell 状态下复核上述性质，并输出主导交换组分
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.cell_balances import calc_gas_balance_residual
from src.core.cell_hydrodynamics import calc_phase_exchange
from src.core.species import GAS_SPECIES, N_GAS
from src.core.reactor import Reactor
from tests.validation_case_utils import PHASE1_HTW_LU_SOLVE_KWARGS, build_phase1_htw_lu_reactor_config


def _audit_algebraic_invariants() -> dict[str, Any]:
    c_b = np.array([2.0, 1.0, 0.5] + [0.0] * (N_GAS - 3), dtype=np.float64)
    c_d = np.array([1.0, 1.5, 0.5] + [0.0] * (N_GAS - 3), dtype=np.float64)
    n_ex = calc_phase_exchange(K_bd=2.5, V_b=0.2, C_b=c_b, C_d=c_d)
    zero = np.zeros(N_GAS, dtype=np.float64)

    res = calc_gas_balance_residual(
        N_zu_d=zero,
        N_rez_d=zero,
        N_d_in=zero,
        R_gas_d=zero,
        N_d=zero,
        N_ex=n_ex,
        N_zu_b=zero,
        N_rez_b=zero,
        N_b_in=zero,
        R_gas_b=zero,
        N_b=zero,
    )

    res_d = res[:N_GAS]
    res_b = res[N_GAS:]
    dense_match = float(np.max(np.abs(res_d - n_ex)))
    bubble_match = float(np.max(np.abs(res_b + n_ex)))
    net_zero = float(np.max(np.abs(res_d + res_b)))

    n_ex_zero = calc_phase_exchange(K_bd=3.0, V_b=0.1, C_b=c_b, C_d=c_b)
    equal_phase_zero = float(np.max(np.abs(n_ex_zero)))

    return {
        "dense_match_max_abs": dense_match,
        "bubble_match_max_abs": bubble_match,
        "net_zero_max_abs": net_zero,
        "equal_phase_zero_max_abs": equal_phase_zero,
        "pass": (
            dense_match <= 1e-12
            and bubble_match <= 1e-12
            and net_zero <= 1e-12
            and equal_phase_zero <= 1e-12
        ),
    }


def _top_exchange_species(n_ex: np.ndarray, k: int = 5) -> list[dict[str, float | str]]:
    order = np.argsort(-np.abs(n_ex))
    rows: list[dict[str, float | str]] = []
    for j in order[:k]:
        rows.append({"species": GAS_SPECIES[int(j)], "N_ex_mol_s": float(n_ex[int(j)])})
    return rows


def _audit_lu_cells() -> dict[str, Any]:
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    selected = sorted(set([0, 1, 2, len(reactor.cells) - 1]))
    rows: list[dict[str, Any]] = []
    worst_dense = 0.0
    worst_bubble = 0.0
    worst_net = 0.0

    for i in selected:
        cell = reactor.cells[i]
        cell.calc_hydrodynamics()
        cell.calc_exchange()
        cell.calc_reactions()

        res_with = calc_gas_balance_residual(
            N_zu_d=cell.N_zu_d,
            N_rez_d=cell.N_rez_d,
            N_d_in=cell.N_d_in,
            R_gas_d=cell.R_gas_d,
            N_d=cell.N_d,
            N_ex=cell.N_ex,
            N_zu_b=cell.N_zu_b,
            N_rez_b=cell.N_rez_b,
            N_b_in=cell.N_b_in,
            R_gas_b=cell.R_gas_b,
            N_b=cell.N_b,
        )
        res_noex = calc_gas_balance_residual(
            N_zu_d=cell.N_zu_d,
            N_rez_d=cell.N_rez_d,
            N_d_in=cell.N_d_in,
            R_gas_d=cell.R_gas_d,
            N_d=cell.N_d,
            N_ex=np.zeros_like(cell.N_ex),
            N_zu_b=cell.N_zu_b,
            N_rez_b=cell.N_rez_b,
            N_b_in=cell.N_b_in,
            R_gas_b=cell.R_gas_b,
            N_b=cell.N_b,
        )

        exch_d = res_with[:N_GAS] - res_noex[:N_GAS]
        exch_b = res_with[N_GAS:] - res_noex[N_GAS:]
        dense_match = float(np.max(np.abs(exch_d - cell.N_ex)))
        bubble_match = float(np.max(np.abs(exch_b + cell.N_ex)))
        net_zero = float(np.max(np.abs(exch_d + exch_b)))
        combined_residual_noex = res_with[:N_GAS] + res_with[N_GAS:]
        reaction_scale = float(np.max(np.abs(cell.R_gas_d + cell.R_gas_b))) if np.any(cell.R_gas_d + cell.R_gas_b) else 0.0
        exchange_scale = float(np.max(np.abs(cell.N_ex))) if np.any(cell.N_ex) else 0.0

        worst_dense = max(worst_dense, dense_match)
        worst_bubble = max(worst_bubble, bubble_match)
        worst_net = max(worst_net, net_zero)

        rows.append(
            {
                "cell": i,
                "xi": float(cell.geo.h_center / cfg.H_bed),
                "T_K": float(cell.T),
                "dense_match_max_abs": dense_match,
                "bubble_match_max_abs": bubble_match,
                "net_zero_max_abs": net_zero,
                "combined_residual_max_abs": float(np.max(np.abs(combined_residual_noex))),
                "exchange_scale_max_abs": exchange_scale,
                "reaction_scale_max_abs": reaction_scale,
                "top_exchange_species": _top_exchange_species(cell.N_ex),
            }
        )

    return {
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0)),
        "worst_dense_match_max_abs": worst_dense,
        "worst_bubble_match_max_abs": worst_bubble,
        "worst_net_zero_max_abs": worst_net,
        "rows": rows,
        "pass": worst_dense <= 1e-10 and worst_bubble <= 1e-10 and worst_net <= 1e-10,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="审计 bubble/dense 相间交换与两相守恒闭合")
    ap.add_argument("--json", action="store_true", help="输出 JSON")
    args = ap.parse_args()

    algebra = _audit_algebraic_invariants()
    lu = _audit_lu_cells()
    out = {"algebraic_invariants": algebra, "lu_cells": lu}

    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        print("=" * 96)
        print("Phase exchange & two-phase conservation audit")
        print("=" * 96)
        print(
            f"algebraic: pass={algebra['pass']} "
            f"dense={algebra['dense_match_max_abs']:.3e} "
            f"bubble={algebra['bubble_match_max_abs']:.3e} "
            f"net={algebra['net_zero_max_abs']:.3e} "
            f"equal_phase_zero={algebra['equal_phase_zero_max_abs']:.3e}"
        )
        print(
            f"LU cells: pass={lu['pass']} converged={lu['converged']} n_iter={lu['n_iter']} "
            f"worst_dense={lu['worst_dense_match_max_abs']:.3e} "
            f"worst_bubble={lu['worst_bubble_match_max_abs']:.3e} "
            f"worst_net={lu['worst_net_zero_max_abs']:.3e}"
        )
        print("-" * 96)
        print(
            f"{'cell':>4} {'xi':>6} {'T(K)':>8} {'|ex|_max':>10} {'|R|_max':>10} "
            f"{'dense_err':>10} {'bubble_err':>11} {'net_err':>10} {'comb_res':>10}"
        )
        for row in lu["rows"]:
            print(
                f"{row['cell']:>4d} {row['xi']:>6.3f} {row['T_K']:>8.1f} "
                f"{row['exchange_scale_max_abs']:>10.3e} {row['reaction_scale_max_abs']:>10.3e} "
                f"{row['dense_match_max_abs']:>10.3e} {row['bubble_match_max_abs']:>11.3e} "
                f"{row['net_zero_max_abs']:>10.3e} {row['combined_residual_max_abs']:>10.3e}"
            )
            tops = ", ".join(
                f"{r['species']}={r['N_ex_mol_s']:.3e}" for r in row["top_exchange_species"]
            )
            print(f"      top N_ex: {tops}")

    return 0 if (algebra["pass"] and lu["pass"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
