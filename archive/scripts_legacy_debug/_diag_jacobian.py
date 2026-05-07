"""Inspect residual Jacobian structure for cell0 stable_window."""
import sys
sys.path.insert(0, '.')

from scripts.audit_single_cell_oxygen_convergence import _prepare_cell, kinetics_scaling, Scenario
from src.solvers.cell_solver import _pack_state, _unpack_state, _build_residual_scales, _residual_wrapper
from src.core.cell import N_SOLID_COMP
from src.core.species import N_GAS, GAS_SPECIES_INDEX as IDX
import numpy as np

sc = Scenario('stable_window', dense=0.40, heat_loss=0.10, r5s=0.50, r6s=1.00, r4s=1.50)
with kinetics_scaling(sc.r5s, sc.r6s, sc.r4s):
    cell = _prepare_cell(sc, 0)
    x0_raw = _pack_state(cell)
    scales = _build_residual_scales(cell)
    n = len(x0_raw)
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    print(f"n_vars={n}  N_GAS={N_GAS}  nk={nk}  N_SOLID_COMP={N_SOLID_COMP}  nv_sol={nv_sol}")
    print(f"T_init={x0_raw[-1]:.1f}")
    print(f"N_d species: {x0_raw[:N_GAS]}")
    print(f"N_b species: {x0_raw[N_GAS:2*N_GAS]}")
    print(f"m_solid sum: {x0_raw[2*N_GAS:2*N_GAS+nv_sol].sum():.4e}")

    # Residuals at T=800K with original species
    x_t = x0_raw.copy(); x_t[-1] = 800.0
    _unpack_state(x_t, cell)
    res_800 = cell.residuals()
    rms_800 = np.sqrt(np.mean((res_800/scales)**2))
    print(f"\nAt T=800K:")
    print(f"  rms={rms_800:.4f}  max={np.max(np.abs(res_800)):.4e}")
    print(f"  any NaN={np.any(np.isnan(res_800))}  any Inf={np.any(np.isinf(res_800))}")
    
    # Check dres/dT finite difference
    x_plus = x_t.copy(); x_plus[-1] = 808.0
    _unpack_state(x_plus, cell)
    res_808 = cell.residuals()
    dT = (res_808 - res_800) / 8.0
    print(f"  dres/dT: any NaN={np.any(np.isnan(dT))}  any Inf={np.any(np.isinf(dT))}  max={np.max(np.abs(dT)):.4e}")
    
    # Compute full Jacobian numerically
    print("\nComputing full Jacobian numerically...")
    _unpack_state(x_t, cell)
    f0 = cell.residuals()
    J = np.zeros((len(f0), n))
    for j in range(n):
        h = max(abs(x_t[j]) * 1e-5, 1e-8)
        x_p = x_t.copy(); x_p[j] += h
        # clip to bounds
        x_p[:-1] = np.maximum(x_p[:-1], 0.0)
        x_p[-1] = np.clip(x_p[-1], 300.0, 3000.0)
        _unpack_state(x_p, cell)
        f_p = cell.residuals()
        J[:, j] = (f_p - f0) / h
    
    # Check Jacobian quality
    print(f"  J any NaN: {np.any(np.isnan(J))}  any Inf: {np.any(np.isinf(J))}")
    print(f"  J max abs: {np.max(np.abs(J)):.4e}")
    
    # Singular values
    try:
        sv = np.linalg.svd(J, compute_uv=False)
        print(f"  Singular values (top 5): {sv[:5]}")
        print(f"  Singular values (bottom 5): {sv[-5:]}")
        print(f"  Condition number: {sv[0]/sv[-1]:.4e}")
        n_zero_sv = np.sum(sv < 1e-10)
        print(f"  Near-zero singular values (< 1e-10): {n_zero_sv}")
    except Exception as e:
        print(f"  SVD failed: {e}")
    
    # Check which columns are all-zero (stuck at bounds?)
    col_norms = np.linalg.norm(J, axis=0)
    zero_cols = np.where(col_norms < 1e-14)[0]
    print(f"\n  Zero-norm Jacobian columns: {len(zero_cols)} → indices {zero_cols[:20]}")
    print(f"  Variables at lower bound (0): {np.sum(x_t[:-1] <= 1e-15)}")
    print(f"  Species names: {list(IDX.keys())}")
