"""Deep analysis: what's the best achievable rms for cell0 stable_window?
Uses x_best (after homotopy) as base, scans temperature seeds with generous budget.
"""
import sys
sys.path.insert(0, '.')

from scripts.audit_single_cell_oxygen_convergence import _prepare_cell, kinetics_scaling, Scenario
from src.solvers.cell_solver import (
    _pack_state, _unpack_state, _build_residual_scales,
    _residual_wrapper, _residual_wrapper_scaled, _homotopy_warmup, solve_cell
)
from src.core.cell import N_SOLID_COMP
from src.core.species import N_GAS, GAS_SPECIES_INDEX as IDX
import numpy as np
from scipy.optimize import least_squares

SCENARIOS = [
    ('stable_window', dict(dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)),
    ('peaky_branch',  dict(dense=0.35, heat_loss=0.12, r5s=1.00, r6s=1.00, r4s=1.00)),
]

T_SEEDS = [850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1300, 1400]

def best_solve(cell, x_base, t_seed, lb, ub, scales, nfev=600, verbose=False):
    """Try to minimize residuals starting from x_base with T=t_seed.
    Uses x_scale='jac' for automatic scaling of ill-conditioned problems.
    """
    n = len(x_base)
    diag = np.ones(n)
    # Better x_scale: use actual variable magnitudes, floor at 1e-3
    diag[:2*N_GAS] = np.maximum(np.abs(x_base[:2*N_GAS]), 1e-3)
    diag[-1] = 1000.0
    
    x_try = np.array(x_base, dtype=float)
    x_try[-1] = float(np.clip(t_seed, lb[-1], ub[-1]))
    
    # First pass: soft_l1 to get out of bad local areas
    r1 = least_squares(
        _residual_wrapper_scaled, x_try, args=(cell, scales), bounds=(lb, ub),
        x_scale='jac', jac='2-point', loss='soft_l1', f_scale=0.5,
        ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=nfev//2
    )
    x1 = np.clip(r1.x, lb, ub)
    
    # Second pass: linear to polish
    r2 = least_squares(
        _residual_wrapper, x1, args=(cell,), bounds=(lb, ub),
        x_scale='jac', jac='2-point', loss='linear',
        ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=nfev//2
    )
    x2 = np.clip(r2.x, lb, ub)
    _unpack_state(x2, cell)
    res = cell.residuals()
    rms = float(np.sqrt(np.mean((res/scales)**2)))
    return x2, rms, cell.T

for sc_name, kw in SCENARIOS:
    sc = Scenario(sc_name, **kw)
    with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
        cell = _prepare_cell(sc, 0)
        x0_raw = _pack_state(cell)
        scales = _build_residual_scales(cell)
        n = len(x0_raw)
        lb = np.zeros(n); lb[-1] = 300.0
        ub = np.full(n, 1e4); ub[-1] = 3000.0

        # First: run solve_cell to get x_best (homotopy-warmed)
        info = solve_cell(cell, stiff_stabilization=True, verbose=False)
        x_best = _pack_state(cell)
        rms_solve = info['rms_scaled']
        print(f"\n{'='*60}")
        print(f"{sc_name} cell0: solve_cell rms={rms_solve:.4f}  T={cell.T:.1f}K")
        print(f"  N_d non-zero species: {np.sum(cell.N_d > 1e-6)}/{N_GAS}")
        print(f"  N_b non-zero species: {np.sum(cell.N_b > 1e-6)}/{N_GAS}")
        
        # Now try T-seeds from x_best using better solver strategy
        print(f"\n  T-seed scan from x_best (x_scale='jac', 600 nfev each):")
        best_overall_rms = rms_solve
        best_overall_T = cell.T
        for t0 in T_SEEDS:
            x2, rms, T_final = best_solve(cell, x_best, t0, lb, ub, scales, nfev=600)
            flag = " *** PHYS" if rms < 0.10 else ""
            print(f"    T0={t0:5d}K => T_final={T_final:7.1f}K  rms={rms:.4f}{flag}")
            if rms < best_overall_rms:
                best_overall_rms = rms
                best_overall_T = T_final
        
        print(f"\n  --> Overall best: rms={best_overall_rms:.4f}  T={best_overall_T:.1f}K")
        
        # Also try: run a fresh homotopy from T=900K
        print(f"\n  Homotopy-from-900K test:")
        cell2 = _prepare_cell(sc, 0)
        x_low = _pack_state(cell2)
        x_low[-1] = 900.0
        _unpack_state(x_low, cell2)
        _homotopy_warmup(cell2, n_steps=7)
        x_hom = _pack_state(cell2)
        scales2 = _build_residual_scales(cell2)
        r_hom = least_squares(
            _residual_wrapper, x_hom, args=(cell2,), bounds=(lb, ub),
            x_scale='jac', jac='2-point', loss='linear',
            ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=1000
        )
        x_hom2 = np.clip(r_hom.x, lb, ub)
        _unpack_state(x_hom2, cell2)
        res_hom = cell2.residuals()
        rms_hom = float(np.sqrt(np.mean((res_hom/scales2)**2)))
        print(f"    Homotopy from 900K => T_final={cell2.T:.1f}K  rms={rms_hom:.4f}")
        
        # Residual breakdown for best solution
        _unpack_state(x_best, cell)  # restore best
        res_best = cell.residuals()
        species_names = ['H2','CO','CO2','N2','O2','H2O','CH4','H2S','NH3','Tar','N2O']
        print(f"\n  Scaled residual breakdown (solve_cell best, rms={rms_solve:.4f}):")
        for i, nm in enumerate(species_names[:N_GAS]):
            sr_d = res_best[i] / scales[i]
            sr_b = res_best[N_GAS+i] / scales[N_GAS+i]
            if abs(sr_d) > 0.05 or abs(sr_b) > 0.05:
                print(f"    {nm:6s}: dense={sr_d:+.3f}  bubble={sr_b:+.3f}")
        sr_T = res_best[-1] / scales[-1]
        print(f"    Energy: {sr_T:+.4f}")
