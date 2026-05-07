"""Diagnostic: trace exactly what happens in cell 1 during first GS iteration."""
import sys, numpy as np
sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES_INDEX as idx, GAS_SPECIES
from src.solvers.cell_solver import solve_cell, _pack_state

cfg = ReactorConfig(
    n_cells=10, H_bed=6.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=293.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.337, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.0, heat_loss_frac=0.08,
)
reactor = Reactor(cfg)
reactor._set_bottom_cell_feeds()
cell0 = reactor.cells[0]

# Pre-fill cell 0 and solve it exactly as GS does
cell0.N_d[:] = np.maximum(cell0.N_d_in + cell0.N_zu_d + cell0.N_rez_d, 0.0)
cell0.N_b[:] = np.maximum(cell0.N_b_in + cell0.N_zu_b + cell0.N_rez_b, 0.0)
cell0.m_solid[:] = np.maximum(cell0.m_solid_in + cell0.m_solid_zu + cell0.m_solid_rez, 0.0)
cell0.calc_hydrodynamics()
cell0.compute_vorabrechnung(cell0.geo.dh / max(cell0.u_mf, 1e-3))
cell0.calc_reactions()
cell0.N_d[:] = np.maximum(cell0.N_d_in + cell0.N_zu_d + 0.2 * cell0.R_gas_d, 0.0)
cell0.N_b[:] = np.maximum(cell0.N_b_in + cell0.N_zu_b + 0.2 * cell0.R_gas_b, 0.0)
cell0.m_solid[:] = np.maximum(cell0.m_solid_in + cell0.m_solid_zu + 0.2 * cell0.R_solid, 0.0)
solve_cell(cell0)

print("=== Cell 0 output ===")
for j, sp in enumerate(GAS_SPECIES):
    if cell0.N_d[j] > 0.001 or cell0.N_b[j] > 0.001:
        print(f"  {sp:5s}: N_d={cell0.N_d[j]:.4f} N_b={cell0.N_b[j]:.4f}")
print(f"  T={cell0.T:.1f}K, sum N_d={np.sum(cell0.N_d):.4f}, sum N_b={np.sum(cell0.N_b):.4f}")

# Set up cell 1
reactor._propagate_upstream(1)
cell1 = reactor.cells[1]
print()
print("=== Cell 1 inputs ===")
for j, sp in enumerate(GAS_SPECIES):
    if cell1.N_d_in[j] > 0.001 or cell1.N_b_in[j] > 0.001:
        print(f"  {sp:5s}: N_d_in={cell1.N_d_in[j]:.4f} N_b_in={cell1.N_b_in[j]:.4f}")
print(f"  sum N_d_in={np.sum(cell1.N_d_in):.4f}, sum N_b_in={np.sum(cell1.N_b_in):.4f}")
print(f"  m_solid_in sum: {np.sum(cell1.m_solid_in):.6f}")
print(f"  N_zu_d sum: {np.sum(cell1.N_zu_d):.6f}, N_zu_b sum: {np.sum(cell1.N_zu_b):.6f}")

# Pre-fill cell 1
cell1.N_d[:] = np.maximum(cell1.N_d_in + cell1.N_zu_d + cell1.N_rez_d, 0.0)
cell1.N_b[:] = np.maximum(cell1.N_b_in + cell1.N_zu_b + cell1.N_rez_b, 0.0)
cell1.m_solid[:] = np.maximum(cell1.m_solid_in + cell1.m_solid_zu + cell1.m_solid_rez, 0.0)
cell1.calc_hydrodynamics()
print(f"  V_d={cell1.V_d:.4f} V_b={cell1.V_b:.4f} K_bd={cell1.K_bd:.6e}")
cell1.compute_vorabrechnung(cell1.geo.dh / max(cell1.u_mf, 1e-3))
cell1.calc_reactions()
print(f"  R_gas_d: max={np.max(np.abs(cell1.R_gas_d)):.4e}")
print(f"  vm_gas_source CH4={cell1._vm_gas_source_cache[idx['CH4']]:.4e}")
for j, sp in enumerate(GAS_SPECIES):
    if abs(cell1.R_gas_d[j]) > 0.001:
        print(f"  R_gas_d[{sp}]={cell1.R_gas_d[j]:.6f}")

cell1.N_d[:] = np.maximum(cell1.N_d_in + cell1.N_zu_d + 0.2 * cell1.R_gas_d, 0.0)
cell1.N_b[:] = np.maximum(cell1.N_b_in + cell1.N_zu_b + 0.2 * cell1.R_gas_b, 0.0)
cell1.m_solid[:] = np.maximum(cell1.m_solid_in + cell1.m_solid_zu + 0.2 * cell1.R_solid, 0.0)

print()
print("=== Cell 1 initial guess for solve_cell ===")
for j, sp in enumerate(GAS_SPECIES):
    if cell1.N_d[j] > 0.001 or cell1.N_b[j] > 0.001:
        print(f"  {sp:5s}: N_d={cell1.N_d[j]:.4f} N_b={cell1.N_b[j]:.4f}")

# Check residuals at initial guess
print()
print("=== Residuals at initial guess ===")
res_pre = cell1.residuals()
print(f"  max |res| = {np.max(np.abs(res_pre)):.4e}")
for j, sp in enumerate(GAS_SPECIES):
    if abs(res_pre[j]) > 0.1 or abs(res_pre[j + 11]) > 0.1:
        print(f"  {sp:5s}: res_d={res_pre[j]:.4e} res_b={res_pre[j+11]:.4e}")
print(f"  energy res = {res_pre[-1]:.4e}")

# Now solve
print()
print("=== Solving cell 1 ===")
# Reset to initial guess
cell1.N_d[:] = np.maximum(cell1.N_d_in + cell1.N_zu_d + 0.2 * cell1.R_gas_d, 0.0)
cell1.N_b[:] = np.maximum(cell1.N_b_in + cell1.N_zu_b + 0.2 * cell1.R_gas_b, 0.0)
cell1.m_solid[:] = np.maximum(cell1.m_solid_in + cell1.m_solid_zu + 0.2 * cell1.R_solid, 0.0)
r = solve_cell(cell1)
print(f"  solve_cell result: {r}")
print()
print("=== Cell 1 output ===")
for j, sp in enumerate(GAS_SPECIES):
    if cell1.N_d[j] > 0.001 or cell1.N_b[j] > 0.001:
        print(f"  {sp:5s}: N_d={cell1.N_d[j]:.4f} N_b={cell1.N_b[j]:.4f}")
print(f"  T={cell1.T:.1f}K, sum N_d={np.sum(cell1.N_d):.4f}, sum N_b={np.sum(cell1.N_b):.4f}")
print()
print("=== Residuals at solution ===")
res_post = cell1.residuals()
print(f"  max |res| = {np.max(np.abs(res_post)):.4e}")
for j, sp in enumerate(GAS_SPECIES):
    if abs(res_post[j]) > 0.1 or abs(res_post[j + 11]) > 0.1:
        print(f"  {sp:5s}: res_d={res_post[j]:.4e} res_b={res_post[j+11]:.4e}")
print(f"  energy res = {res_post[-1]:.4e}")
