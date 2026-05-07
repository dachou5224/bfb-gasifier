"""LU 工况 GS 迭代历史审计。

输出每次全局迭代的 dT、max_res、T_exit、carbon_conv，
用于定位“看似收敛但精度失真”或“状态发散”模式。
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config


def main() -> int:
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    result = reactor.solve(max_global_iter=20, tol_global=1.0, solver="gauss_seidel", verbose=False)

    print("=" * 76)
    print("LU GS history audit")
    print("=" * 76)
    print(f"converged={result.get('converged')} n_iter={result.get('n_iter')} rms_scaled_gs={result.get('rms_scaled_gs'):.3e}")
    print(f"T_exit={result['T_profile'][-1]:.2f} carbon_conv={result['carbon_conv']:.4f}")
    print("-" * 76)
    print(
        f"{'iter':>4s} {'dT_max':>12s} {'max_res':>12s} {'max_rms':>12s} {'T_exit':>10s} "
        f"{'prof_pen':>10s} {'runaway':>10s} {'o2_reb':>10s} {'vm_reb':>10s} {'moist_reb':>10s} {'o2_slip':>10s}"
    )
    print("-" * 76)

    for row in result.get("history", []):
        print(
            f"{row['iter']:4d} {row['dT_max']:12.4e} {row['max_res']:12.4e} "
            f"{row['max_rms']:12.4e} {row['T_exit']:10.2f} "
            f"{row.get('profile_penalty', float('nan')):10.4e} "
            f"{row.get('runaway_penalty', float('nan')):10.4e} "
            f"{row.get('o2_rebound', float('nan')):10.4e} "
            f"{row.get('vm_rebound', float('nan')):10.4e} "
            f"{row.get('moist_rebound', float('nan')):10.4e} "
            f"{row.get('o2_slip', float('nan')):10.4e}"
        )

    print("=" * 76)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
