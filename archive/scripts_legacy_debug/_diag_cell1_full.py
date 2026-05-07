"""Trace cell 1 in full GS run."""
import sys, numpy as np
sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES_INDEX as idx, GAS_SPECIES
from src.solvers.cell_solver import solve_cell

cfg = ReactorConfig(
    n_cells=10, H_bed=6.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=293.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.337, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.08,
)
reactor = Reactor(cfg)

# Manually replicate the GS first iteration
reactor._set_bottom_cell_feeds()
reactor.cells[0].T_rez_gas = reactor.cells[-1].T

for i in range(len(reactor.cells)):
    reactor._propagate_upstream(i)
    cell = reactor.cells[i]
    old_T = cell.T

    cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + cell.N_rez_d, 0.0)
    cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + cell.N_rez_b, 0.0)
    cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + cell.m_solid_rez, 0.0)

    cell.calc_hydrodynamics()
    cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
    cell.calc_reactions()
    
    if i == 1:
        print(f"=== Cell 1 before second pre-fill ===")
        print(f"  R_gas_d max={np.max(np.abs(cell.R_gas_d)):.4e}")
        print(f"  R_gas_d[CH4]={cell.R_gas_d[idx['CH4']]:.4e}")
        print(f"  R_gas_d[H2]={cell.R_gas_d[idx['H2']]:.4e}")
        print(f"  R_gas_d[CO]={cell.R_gas_d[idx['CO']]:.4e}")
        print(f"  m_solid sum={np.sum(cell.m_solid):.6e}")
        print(f"  N_d[N2]={cell.N_d[idx['N2']]:.4f}, N_b[N2]={cell.N_b[idx['N2']]:.4f}")
        print(f"  N_d sum={np.sum(cell.N_d):.4f}")

    cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + 0.2 * cell.R_gas_d, 0.0)
    cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + 0.2 * cell.R_gas_b, 0.0)
    cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + 0.2 * cell.R_solid, 0.0)

    if i == 1:
        print(f"\n=== Cell 1 initial guess for solve_cell ===")
        for j, sp in enumerate(GAS_SPECIES):
            if cell.N_d[j] > 0.001 or cell.N_b[j] > 0.001:
                print(f"  {sp:5s}: N_d={cell.N_d[j]:.4f} N_b={cell.N_b[j]:.4f}")

    r = solve_cell(cell)
    cell.T = old_T + 0.4 * (cell.T - old_T)

    if i == 1:
        print(f"\n=== Cell 1 after solve_cell ===")
        print(f"  solve result: converged={r['converged']}, residual={r['residual']:.4e}")
        for j, sp in enumerate(GAS_SPECIES):
            if cell.N_d[j] > 0.001 or cell.N_b[j] > 0.001:
                print(f"  {sp:5s}: N_d={cell.N_d[j]:.4f} N_b={cell.N_b[j]:.4f}")
        print(f"  T={cell.T:.1f}K, N2_d={cell.N_d[idx['N2']]:.4f}, N2_b={cell.N_b[idx['N2']]:.4f}")

print("\n=== Final state after all cells (1 GS pass) ===")
for i, cell in enumerate(reactor.cells):
    n2_d = cell.N_d[idx['N2']]
    n2_b = cell.N_b[idx['N2']]
    ch4_d = cell.N_d[idx['CH4']]
    total = np.sum(cell.N_d) + np.sum(cell.N_b)
    print(f"  cell{i}: T={cell.T:.0f}K N2_d={n2_d:.3f} N2_b={n2_b:.3f} CH4_d={ch4_d:.3f} total={total:.3f}")
