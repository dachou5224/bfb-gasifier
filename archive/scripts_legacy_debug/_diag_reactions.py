"""诊断各反应速率——追踪H2/CO2缺失的根因。"""
import sys, numpy as np
sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX as IDX

cfg = ReactorConfig(
    n_cells=10, H_bed=5.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=300.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.25, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.1,
)
reactor = Reactor(cfg)
res = reactor.solve(max_global_iter=5, tol_global=1e-3)

print("=== 各 cell 关键气体组分（干基 mol/s）===")
print(f"  {'cell':>4}  {'T(K)':>8}  {'CO':>10}  {'CO2':>10}  {'H2':>10}  {'H2O':>10}  {'CH4':>10}  {'N2':>10}")
print("  " + "-"*78)
for i, cell in enumerate(reactor.cells):
    Nd = cell.N_d
    Nb = cell.N_b
    Ntot = Nd + Nb
    co  = Ntot[IDX["CO"]]
    co2 = Ntot[IDX["CO2"]]
    h2  = Ntot[IDX["H2"]]
    h2o = Ntot[IDX["H2O"]]
    ch4 = Ntot[IDX["CH4"]] if "CH4" in IDX else 0.0
    n2  = Ntot[IDX["N2"]]
    print(f"  {i:>4}  {cell.T:>8.1f}  {co:>10.4f}  {co2:>10.6f}  {h2:>10.4f}  {h2o:>10.4f}  {ch4:>10.6f}  {n2:>10.4f}")

print()
print("=== 底层 cell 0 详细反应速率 ===")
cell0 = reactor.cells[0]
cell0.calc_hydrodynamics()
cell0.compute_vorabrechnung(cell0.geo.dh / max(cell0.u_mf, 1e-3))
cell0.calc_reactions()

# 检查 R_gas_d (dense-phase 气体生成速率，mol/s)
print(f"  T_cell0 = {cell0.T:.1f} K")
print(f"  R_gas_d (mol/s 各组分变化):")
for j, sp in enumerate(GAS_SPECIES):
    v = float(np.sum(cell0.R_gas_d[j])) if cell0.R_gas_d.ndim > 1 else float(cell0.R_gas_d[j])
    if abs(v) > 1e-9:
        print(f"    {sp:8s}: {v:+12.6f}")

print()
print("=== 中间 cell 5 详细反应速率 ===")
cell5 = reactor.cells[5]
cell5.calc_hydrodynamics()
cell5.compute_vorabrechnung(cell5.geo.dh / max(cell5.u_mf, 1e-3))
cell5.calc_reactions()
print(f"  T_cell5 = {cell5.T:.1f} K")
print(f"  R_gas_d (mol/s 各组分变化):")
for j, sp in enumerate(GAS_SPECIES):
    v = float(np.sum(cell5.R_gas_d[j])) if cell5.R_gas_d.ndim > 1 else float(cell5.R_gas_d[j])
    if abs(v) > 1e-9:
        print(f"    {sp:8s}: {v:+12.6f}")

print()
print("=== 进口条件核查 ===")
print(f"  O2_feed  = {cfg.O2_feed:.4f} mol/s")
print(f"  H2O_feed = {cfg.H2O_feed:.4f} mol/s  (steam/O2 = {cfg.steam_to_o2_molar:.2f})")
print(f"  N2_feed  = {cfg.N2_feed:.4f} mol/s")
c0 = reactor.cells[0]
print(f"  N_zu_d[H2O] = {c0.N_zu_d[IDX['H2O']]:.4f}  N_zu_b[H2O] = {c0.N_zu_b[IDX['H2O']]:.4f}")
print(f"  N_zu_d[O2]  = {c0.N_zu_d[IDX['O2']]:.4f}  N_zu_b[O2]  = {c0.N_zu_b[IDX['O2']]:.4f}")
print(f"  N_zu_d[N2]  = {c0.N_zu_d[IDX['N2']]:.4f}  N_zu_b[N2]  = {c0.N_zu_b[IDX['N2']]:.4f}")
