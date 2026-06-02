"""诊断脚本：整炉轴向剖面详细分析。"""
import sys
import numpy as np
sys.path.insert(0, ".")

from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES
from src.solvers.cell_solver import _build_residual_scales

cfg = ReactorConfig(
    n_cells=10, H_bed=5.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=300.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.25, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.1,
)
reactor = Reactor(cfg)
res = reactor.solve(max_global_iter=10, tol_global=1e-3)

print("=== 轴向温度剖面 ===")
idx = reactor.cells[0].geo
from src.core.species import GAS_SPECIES_INDEX as IDX

for i, cell in enumerate(reactor.cells):
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    cell.calc_reactions()
    raw_res = cell.residuals()
    scales = _build_residual_scales(cell)
    rms = float(np.sqrt(np.mean((raw_res / scales) ** 2)))
    N_total = float(np.sum(cell.N_b + cell.N_d))
    o2_frac = float((cell.N_b[IDX["O2"]] + cell.N_d[IDX["O2"]]) / max(N_total, 1e-9))
    co_frac = float((cell.N_b[IDX["CO"]] + cell.N_d[IDX["CO"]]) / max(N_total, 1e-9))
    h2_frac = float((cell.N_b[IDX["H2"]] + cell.N_d[IDX["H2"]]) / max(N_total, 1e-9))
    m_char = float(np.sum(cell.m_solid[:, 0]))  # S_CHAR=0
    eps_b = cell.eps_b
    print(
        f"  cell{i:2d}: T={cell.T:7.1f}K  rms={rms:.3e}  eps_b={eps_b:.3f}"
        f"  O2={o2_frac*100:.2f}%  CO={co_frac*100:.2f}%  H2={h2_frac*100:.2f}%"
        f"  m_char={m_char:.4f}kg/s"
    )

print()
print("=== 出口气体 (干基) ===")
dg = res["exit_gas_dry"]
for sp in ["CO", "CO2", "H2", "H2O", "CH4", "O2", "N2"]:
    if sp in dg:
        print(f"  {sp}: {dg[sp]*100:.2f}%")

print()
print(f"  T_exit={res['T_profile'][-1]:.1f}K   (目标~1173K)")
print(f"  Xc={res['carbon_conv']:.3f}   (炭转化率)")
print(f"  converged={res['converged']}   n_iter={res['n_iter']}")
