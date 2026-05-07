"""Debug recirculation state."""
import sys
import numpy as np

sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig

cfg = ReactorConfig(
    n_cells=10, H_bed=5.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=300.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.25, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.1,
)
reactor = Reactor(cfg)

# Run 1 iteration manually to see what happens
reactor._set_bottom_cell_feeds()
if cfg.recirculation_frac > 0:
    reactor.cells[0].T_rez_gas = reactor.cells[-1].T

# Step through cell 0 and 1
for i in range(2):
    reactor._propagate_upstream(i)
    cell = reactor.cells[i]
    print(f"\nCell {i} BEFORE solve:")
    print(f"  T={cell.T:.1f}K")
    print(f"  N_rez_d[O2]={cell.N_rez_d[1]:.6f}")
    print(f"  m_solid_zu[char]={np.sum(cell.m_solid_zu[:, 0]):.6f}")
    
    from src.solvers.cell_solver import solve_cell
    cell.calc_hydrodynamics()
    cell.compute_vorabrechnung(cell.geo.dh/max(cell.u_mf, 1e-3))
    cell.calc_reactions()
    cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + 0.2 * cell.R_gas_d, 0.0)
    cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + 0.2 * cell.R_gas_b, 0.0)
    cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + 0.2 * cell.R_solid, 0.0)
    res_cell = solve_cell(cell)
    
    print(f"\nCell {i} AFTER solve:")
    print(f"  T={cell.T:.1f}K")
    print(f"  N_rez_d[O2]={cell.N_rez_d[1]:.6f}")
    print(f"  m_solid_zu[char]={np.sum(cell.m_solid_zu[:, 0]):.6f}")
    print(f"  residual={res_cell.get('residual', 0.0):.3e}")
