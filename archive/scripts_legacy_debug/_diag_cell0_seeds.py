"""Scan T-seeds for cell0 in both scenarios to map the convergence landscape."""
import sys
sys.path.insert(0, '.')

from scripts.audit_single_cell_oxygen_convergence import _prepare_cell, kinetics_scaling, Scenario
from src.solvers.cell_solver import _pack_state, _unpack_state, _build_residual_scales, _residual_wrapper
import numpy as np
from scipy.optimize import least_squares

SCENARIOS = [
    ('stable_window', dict(dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)),
    ('peaky_branch',  dict(dense=0.35, heat_loss=0.12, r5s=1.00, r6s=1.00, r4s=1.00)),
]

T_SEEDS = [800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400, 1500]

for sc_name, kw in SCENARIOS:
    sc = Scenario(sc_name, **kw)
    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        # check cell index 0 only (the persistent failure)
        for ci in [0]:
            cell = _prepare_cell(sc, ci)
            x0_raw = _pack_state(cell)
            scales = _build_residual_scales(cell)
            n = len(x0_raw)
            lb = np.full(n, 0.0); lb[-1] = 300.0
            ub = np.full(n, 1e4); ub[-1] = 2500.0

            from src.core.species import GAS_SPECIES_INDEX as IDX
            O2_d = float(cell.N_zu_d[IDX['O2']])
            O2_b = float(cell.N_zu_b[IDX['O2']])
            print(f"\n=== {sc_name} cell{ci} (h={cell.geo.h_center:.3f}m) T_init={x0_raw[-1]:.0f}K  O2_zu_d={O2_d:.3e}  O2_zu_b={O2_b:.3e} ===")

            best_rms = 1e9
            best_T   = None
            for t0 in T_SEEDS:
                x_t = x0_raw.copy(); x_t[-1] = float(t0)
                r = least_squares(
                    _residual_wrapper, x_t, args=(cell,), bounds=(lb, ub),
                    x_scale=np.ones(n), jac='2-point', loss='linear',
                    ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=300
                )
                xo = np.clip(r.x, lb, ub)
                _unpack_state(xo, cell)
                res = cell.residuals()
                rms = float(np.sqrt(np.mean((res / scales) ** 2)))
                flag = " *** PHYS" if rms < 0.10 else ""
                print(f"  T0={t0:5d}K => T_final={cell.T:7.1f}K  rms={rms:.4f}{flag}")
                if rms < best_rms:
                    best_rms = rms
                    best_T = cell.T
            print(f"  --> best rms={best_rms:.4f}  T={best_T:.1f}K")
