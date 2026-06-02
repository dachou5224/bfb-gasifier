#!/usr/bin/env python3
"""LU 工况热解产物分配审计。

目标：
- 量化 shared global NR 解上各 cell 的 `x_vm` / `x_dry` / VM 气源；
- 审计 `allocate_pyrolysis_products_elemental()` 对 `CO/CH4/TAR1/TAR2` 的分配模式；
- 明确 `CH4` 是否在源头就偏弱，以及 `TAR2` 是否被完全闲置。
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    print("=" * 150)
    print("LU pyrolysis allocator audit (shared global NR)")
    print("=" * 150)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"init={result.get('nr_init_strategy')} jacobian={result.get('nr_jacobian_strategy')} "
        f"Texit={result['T_profile'][-1]:.1f}K n_iter={result.get('n_iter')} "
        f"rms={float(result.get('rms_scaled_final', float('nan'))):.3e}"
    )
    print("-" * 150)
    print(
        f"{'cell':>4} {'xi/H':>6} {'T[K]':>8} {'tau[s]':>8} {'m_VM_in':>10} {'m_Moist':>10} "
        f"{'x_vm':>7} {'x_dry':>7} {'CO':>8} {'CH4':>8} {'TAR1':>8} {'TAR2':>8} {'H2O':>8} {'NH3':>8}"
    )
    print("-" * 150)

    total_tar1 = 0.0
    total_tar2 = 0.0
    total_ch4 = 0.0
    total_co = 0.0
    active_cells = 0

    for i, cell in enumerate(reactor.cells):
        cell.calc_hydrodynamics()
        tau = cell.geo.dh / max(cell.u_mf, 1e-3)
        m_vm_in = float(cell.m_solid_zu[:, 1].sum() + cell.m_solid_in[:, 1].sum())
        m_moist_in = float(cell.m_solid_zu[:, 2].sum() + cell.m_solid_in[:, 2].sum())
        bundle = calc_drying_pyrolysis_sources(
            tau=tau,
            T=cell.T,
            P=cell.P,
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
            m_vm_in=m_vm_in,
            m_moist_in=m_moist_in,
            solid_shape=cell.R_solid.shape,
            char_index=0,
            vm_index=1,
            moisture_index=2,
        )

        src = bundle.gas_source
        co = float(src[idx["CO"]])
        ch4 = float(src[idx["CH4"]])
        tar1 = float(src[idx["TAR1"]])
        tar2 = float(src[idx["TAR2"]])
        h2o = float(src[idx["H2O"]])
        nh3 = float(src[idx["NH3"]])

        if max(co, ch4, tar1, tar2, h2o, nh3) > 1e-12:
            active_cells += 1
        total_co += co
        total_ch4 += ch4
        total_tar1 += tar1
        total_tar2 += tar2

        print(
            f"{i:>4d} {float(cell.geo.h_center/cfg.H_bed):>6.2f} {cell.T:>8.1f} {tau:>8.3f} "
            f"{m_vm_in:>10.4e} {m_moist_in:>10.4e} {bundle.x_vm:>7.3f} {bundle.x_dry:>7.3f} "
            f"{co:>8.4f} {ch4:>8.4f} {tar1:>8.4f} {tar2:>8.4f} {h2o:>8.4f} {nh3:>8.4f}"
        )

    print("-" * 150)
    print(
        f"totals: active_cells={active_cells} CO={total_co:.4f} CH4={total_ch4:.4f} "
        f"TAR1={total_tar1:.4f} TAR2={total_tar2:.4f}"
    )
    print("Notes:")
    print("  - 若 active_cells 仅为 1，说明当前 shared NR 下 VM 释放几乎只发生在 cell0。")
    print("  - 若 TAR2 总为 0，则当前 allocator 实际只在使用单 tar surrogate。")
    print("  - 若 CH4 总量在 pyrolysis 源头就很小，后续再调 R7 也无法把出口 CH4 拉回文献量级。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
