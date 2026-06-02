"""LU 工况单轮 GS 扫描逐 cell 诊断。

目标：定位哪些 cell 在一次扫描中没有被真正求稳，及其与温度 / O2 / VM 的关系。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.validation_case_utils import build_phase1_htw_lu_reactor_config
from src.core.reactor import Reactor
from src.core.cell import S_CHAR, S_VM, S_MOISTURE
from src.core.species import GAS_SPECIES_INDEX as idx
from src.solvers.cell_solver import solve_cell


def main() -> int:
    cfg = build_phase1_htw_lu_reactor_config()
    r = Reactor(cfg)
    r._set_bottom_cell_feeds()

    print("=" * 120)
    print("LU GS single-sweep cell trace")
    print("=" * 120)
    print("cell  h(m)   T_pre   res_max   rms_sc   conv  opt_ok   O2_out   CO_out  CO2_out  H2_out   CH4_out  mVM      mMoist   mChar")
    print("-" * 120)

    max_res = 0.0
    for i in range(len(r.cells)):
        r._propagate_upstream(i)
        cell = r.cells[i]
        if np.sum(np.maximum(cell.N_d + cell.N_b, 0.0)) < 1e-9:
            cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + cell.N_rez_d, 0.0)
            cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + cell.N_rez_b, 0.0)
            cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez, 0.0)

        cell.calc_hydrodynamics()
        cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))

        for _pc in range(2):
            cell.calc_exchange()
            cell.calc_reactions()
            mig = cell._calc_size_migration()
            cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez + cell.R_solid + mig, 0.0)
            cell.calc_hydrodynamics()

        t_pre = cell.T
        res_cell = solve_cell(cell, stiff_stabilization=True, verbose=False)
        max_res = max(max_res, float(res_cell.get("residual", 0.0)))

        cell.calc_exchange()
        cell.calc_reactions()
        mig = cell._calc_size_migration()
        cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez + cell.R_solid + mig, 0.0)

        N = np.maximum(cell.N_d + cell.N_b, 0.0)
        h = (i + 0.5) * cfg.H_bed / cfg.n_cells
        print(
            f"{i:>3d}  {h:>4.1f}  {t_pre:>6.1f}  {res_cell.get('residual', 0.0):>8.2e}  {res_cell.get('rms_scaled', 0.0):>7.2e}  "
            f"{str(bool(res_cell.get('physically_converged', False))):>5s}  {str(bool(res_cell.get('optimizer_success', False))):>6s}  {N[idx['O2']]:>7.3f}  {N[idx['CO']]:>7.3f}  "
            f"{N[idx['CO2']]:>8.3f}  {N[idx['H2']]:>7.3f}  {N[idx['CH4']]:>8.3f}  "
            f"{np.sum(cell.m_solid[:, S_VM]):>7.5f}  {np.sum(cell.m_solid[:, S_MOISTURE]):>8.5f}  {np.sum(cell.m_solid[:, S_CHAR]):>7.5f}"
        )

    print("-" * 120)
    print(f"sweep max residual = {max_res:.3e}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
