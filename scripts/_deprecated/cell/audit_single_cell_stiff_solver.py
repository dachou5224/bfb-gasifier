"""单 cell 刚性稳定化对比审计。

对同一初值，比较 solve_cell 在：
1) 默认模式（stiff_stabilization=False）
2) 刚性稳定化模式（stiff_stabilization=True）
下的 residual / rms_scaled / 关键组分。
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.validation_case_utils import build_phase1_htw_lu_reactor_config
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX as idx
from src.solvers.cell_solver import solve_cell


def _prepare_cell(cell_index: int):
    cfg = build_phase1_htw_lu_reactor_config()
    r = Reactor(cfg)
    r._set_bottom_cell_feeds()

    for i in range(cell_index + 1):
        r._propagate_upstream(i)
        c = r.cells[i]
        if np.sum(np.maximum(c.N_d + c.N_b, 0.0)) < 1e-9:
            c.N_d[:] = np.maximum(c.N_d_in + c.N_zu_d + c.N_rez_d, 0.0)
            c.N_b[:] = np.maximum(c.N_b_in + c.N_zu_b + c.N_rez_b, 0.0)
            c.m_solid[:] = np.maximum(c.m_solid_in + c.m_solid_zu + c.m_solid_rez, 0.0)

        c.calc_hydrodynamics()
        c.compute_vorabrechnung(c.geo.dh / max(c.u_mf, 1e-3))
        for _ in range(2):
            c.calc_exchange()
            c.calc_reactions()
            mig = c._calc_size_migration()
            c.m_solid[:] = np.maximum(c.m_solid_in + c.m_solid_zu + c.m_solid_rez + c.R_solid + mig, 0.0)
            c.calc_hydrodynamics()

    return r.cells[cell_index]


def _run_once(cell, stiff: bool):
    c = copy.deepcopy(cell)
    out = solve_cell(c, stiff_stabilization=stiff, verbose=False)
    N = np.maximum(c.N_d + c.N_b, 0.0)
    return out, {
        "O2": float(N[idx["O2"]]),
        "CO": float(N[idx["CO"]]),
        "CO2": float(N[idx["CO2"]]),
        "H2": float(N[idx["H2"]]),
        "CH4": float(N[idx["CH4"]]),
        "T": float(c.T),
    }


def main() -> int:
    print("=" * 80)
    print("Single-cell stiff solver audit")
    print("=" * 80)

    for cell_index in (0, 1, 3):
        base_cell = _prepare_cell(cell_index)
        off, y_off = _run_once(base_cell, stiff=False)
        on, y_on = _run_once(base_cell, stiff=True)

        print(f"\n[cell {cell_index}]")
        print(
            f"  default:  opt_ok={off.get('optimizer_success')} rms={off.get('rms_scaled'):.3e} "
            f"res={off.get('residual'):.3e}"
        )
        print(
            f"  stiff-on: opt_ok={on.get('optimizer_success')} rms={on.get('rms_scaled'):.3e} "
            f"res={on.get('residual'):.3e} attempted={on.get('stiff_attempted')} accepted={on.get('stiff_accepted')}"
        )
        print(
            f"  gas(default) O2={y_off['O2']:.3f} CO={y_off['CO']:.3f} CO2={y_off['CO2']:.3f} "
            f"H2={y_off['H2']:.3f} CH4={y_off['CH4']:.3f} T={y_off['T']:.1f}"
        )
        print(
            f"  gas(stiff)   O2={y_on['O2']:.3f} CO={y_on['CO']:.3f} CO2={y_on['CO2']:.3f} "
            f"H2={y_on['H2']:.3f} CH4={y_on['CH4']:.3f} T={y_on['T']:.1f}"
        )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
