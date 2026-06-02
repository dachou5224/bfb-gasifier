"""HTW LU 基线流体力学快照审计。

目的：
1. 复用 `tests/validation_case_utils.py` 的 Phase 1 基线配置；
2. 直接查看底部 cell 在给定气体进料下的 `u0/u_mf/d_b/u_b/epsilon_b/K_bd`；
3. 对比“由当前状态反推 u0”与“显式固定 u0_target”两种模式。

用法：
    python scripts/audit_htw_hydrodynamics_snapshot.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config, load_case_LU


def seed_bottom_cell_with_inlet_gas(reactor: Reactor) -> None:
    cfg = reactor.config
    cell = reactor.cells[0]
    idx = GAS_SPECIES_INDEX

    # 复制 reactor._set_bottom_cell_feeds() 的主气体分配逻辑，
    # 但这里直接写到状态向量 N_b / N_d，便于单独审 hydrodynamics。
    cell.N_b.fill(0.0)
    cell.N_d.fill(0.0)

    cell.N_d[idx["O2"]] = cfg.O2_feed * 0.2
    cell.N_d[idx["H2O"]] = cfg.H2O_feed * 0.2
    cell.N_d[idx["N2"]] = cfg.N2_feed * 0.2

    cell.N_b[idx["O2"]] = cfg.O2_feed * 0.8
    cell.N_b[idx["H2O"]] = cfg.H2O_feed * 0.8
    cell.N_b[idx["N2"]] = cfg.N2_feed * 0.8

    # 取 HTW 代表性床层温度而不是冷态入口温度，便于看操作区间内水力学。
    cell.T = 1100.0
    cell.P = cfg.P


def run_snapshot(u0_target: float | None) -> None:
    case = load_case_LU()
    cfg = build_phase1_htw_lu_reactor_config(case)
    cfg.u0_target = u0_target
    reactor = Reactor(cfg)
    seed_bottom_cell_with_inlet_gas(reactor)

    cell = reactor.cells[0]
    cell.calc_hydrodynamics()

    print("-" * 72)
    mode = "inferred-from-state" if u0_target is None else f"fixed-u0={u0_target:.3f} m/s"
    print(f"Mode: {mode}")
    print(f"Case: {case['_source_case_key']}  P={cfg.P/1e6:.2f} MPa  T={cell.T:.1f} K")
    print(f"Feeds: O2={cfg.O2_feed:.3f}  H2O={cfg.H2O_feed:.3f}  N2={cfg.N2_feed:.3f} mol/s")
    print()
    print(f"u0       = {cell.u0:.5f} m/s")
    print(f"u_mf     = {cell.u_mf:.5f} m/s")
    print(f"d_b      = {cell.d_b:.5f} m")
    print(f"u_b      = {cell.u_b:.5f} m/s")
    print(f"epsilon_b= {cell.eps_b:.5f}")
    print(f"epsilon_d= {cell.eps_d:.5f}")
    print(f"K_bd     = {cell.K_bd:.5f} 1/s")

    ratio = cell.u0 / max(cell.u_mf, 1e-12)
    print(f"u0/u_mf  = {ratio:.3f}")

    print("Quick checks:")
    print(f"  u_mf in sanity band [0.02, 0.08]? {'YES' if 0.02 <= cell.u_mf <= 0.08 else 'NO'}")
    print(f"  d_b  in sanity band [0.05, 0.30]? {'YES' if 0.05 <= cell.d_b <= 0.30 else 'NO'}")
    print(f"  K_bd in sanity band [1, 15]?      {'YES' if 1.0 <= cell.K_bd <= 15.0 else 'NO'}")


if __name__ == "__main__":
    print("HTW LU hydrodynamics snapshot audit")
    run_snapshot(u0_target=None)
    run_snapshot(u0_target=1.0)
