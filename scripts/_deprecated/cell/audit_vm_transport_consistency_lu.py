#!/usr/bin/env python3
"""LU 工况下 VM 传播与热解源项一致性审计。

目标：
- 明确 shared ``global_nr`` 解里每个 cell 的 VM 库存、VM 入口和热解源项入口；
- 识别“state 里有 VM，但热解源项入口为 0”的 orphan VM cells；
- 证明当前 shared NR 下 `VM` 几乎只从底部 fresh feed 进入热解计算。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import S_MOISTURE, S_VM
from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def _sum_pos(arr: np.ndarray) -> float:
    return float(np.sum(np.maximum(arr, 0.0)))


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS),
        check_x0=True,
    )

    orphan_cells: list[int] = []
    active_py_cells: list[int] = []

    print("=" * 170)
    print("LU VM transport consistency audit (shared global NR)")
    print("=" * 170)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e} "
        f"Texit={result['T_profile'][-1]:.1f}K"
    )
    print("-" * 170)
    print(
        f"{'cell':>4} {'xi':>5} {'T[K]':>8} {'state_VM':>10} {'zu_VM':>10} {'in_VM':>10} "
        f"{'pyro_VM_in':>11} {'vm_rel':>10} {'x_vm':>8} {'state_H2O':>11} {'zu_H2O':>10} {'in_H2O':>10} {'orphan_VM':>10}"
    )
    print("-" * 170)

    for i, cell in enumerate(reactor.cells):
        cell.calc_hydrodynamics()
        tau = float(cell.geo.dh / max(cell.u_mf, 1e-3))
        zu_vm = _sum_pos(cell.m_solid_zu[:, S_VM])
        in_vm = _sum_pos(cell.m_solid_in[:, S_VM])
        state_vm = _sum_pos(cell.m_solid[:, S_VM])
        pyro_vm_in = _sum_pos(cell.m_solid_zu[:, S_VM] + cell.m_solid_in[:, S_VM])

        zu_h2o = _sum_pos(cell.m_solid_zu[:, S_MOISTURE])
        in_h2o = _sum_pos(cell.m_solid_in[:, S_MOISTURE])
        state_h2o = _sum_pos(cell.m_solid[:, S_MOISTURE])

        bundle = calc_drying_pyrolysis_sources(
            tau=tau,
            T=float(cell.T),
            P=float(cell.P),
            d_p=cell.solid.d_p,
            moisture_wt=cell.solid.moisture_wt,
            ash_dry_wt=cell.solid.ash_dry_wt,
            C_dry=cell.solid.C_dry,
            H_dry=cell.solid.H_dry,
            O_dry=cell.solid.O_dry,
            nitrogen_fraction=cell.solid.nitrogen_fraction,
            sulfur_fraction=cell.solid.sulfur_fraction,
            sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
            pyrolysis_tar_carbon_frac=cell.solid.pyrolysis_tar_carbon_frac,
            fuel_type=cell.fuel_type,
            m_vm_in=pyro_vm_in,
            m_moist_in=_sum_pos(cell.m_solid_zu[:, S_MOISTURE] + cell.m_solid_in[:, S_MOISTURE]),
            solid_shape=cell.R_solid.shape,
            char_index=0,
            vm_index=S_VM,
            moisture_index=S_MOISTURE,
        )
        vm_rel = float(max(-np.sum(bundle.solid_sink[:, S_VM]), 0.0))

        orphan_vm = bool(state_vm > 1e-8 and pyro_vm_in <= 1e-12)
        if orphan_vm:
            orphan_cells.append(i)
        if vm_rel > 1e-12:
            active_py_cells.append(i)

        print(
            f"{i:>4d} {float(cell.geo.h_center / cfg.H_bed):>5.2f} {float(cell.T):>8.1f} "
            f"{state_vm:>10.4e} {zu_vm:>10.4e} {in_vm:>10.4e} {pyro_vm_in:>11.4e} "
            f"{vm_rel:>10.4e} {bundle.x_vm:>8.3f} {state_h2o:>11.4e} {zu_h2o:>10.4e} {in_h2o:>10.4e} {str(orphan_vm):>10}"
        )

    print("-" * 170)
    print(f"active pyrolysis cells: {active_py_cells}")
    print(f"orphan VM cells       : {orphan_cells}")
    print("Notes:")
    print("  - `pyro_VM_in = m_solid_zu[VM] + m_solid_in[VM]`，是当前热解源项真正看到的 VM 入口。")
    print("  - 若某 cell `state_VM > 0` 但 `pyro_VM_in = 0`，说明状态里的 VM 与热解源项链路已脱钩。")
    print("  - 当前 `_propagated_solid_stream()` 只传播 `CHAR + ASH`，因此上部 cell 通常不会得到新的 VM 入口。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
