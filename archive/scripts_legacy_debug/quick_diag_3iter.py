"""快速诊断：3次迭代后各cell状态。"""
import sys, numpy as np
sys.path.insert(0, '/Users/liuzhen/AI-projects/bfb-gasifier')

from tests.validation_case_utils import build_phase1_htw_lu_reactor_config
from src.core.reactor import Reactor
from src.core.cell import S_CHAR, S_VM
from src.core.species import GAS_SPECIES_INDEX as idx, gas_diffusivity_correlation
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R4_effective

cfg = build_phase1_htw_lu_reactor_config()
r = Reactor(cfg)
res = r.solve(max_global_iter=3, tol_global=1.0, solver='gauss_seidel', verbose=False)

print(f"n_iter={res['n_iter']} residual={res['residual']:.3e}")
print(f"T_exit={res['T_profile'][-1]:.1f} carbon_conv={res['carbon_conv']:.4f}")

bot, top = r.cells[0], r.cells[-1]
m_char_zu = bot.m_solid_zu[:,S_CHAR].sum()
m_char_out = top.m_solid[:,S_CHAR].sum()
print(f"bot.m_solid_zu[CHAR]={m_char_zu:.5f}  top.m_solid[CHAR]={m_char_out:.5f}")

print('--- 各cell T m_char area r1 r2 r4 o2_sup h2o_sup ---')
for i, c in enumerate(r.cells):
    c.calc_hydrodynamics(); c.calc_exchange()
    D_g = gas_diffusivity_correlation(c.T, c.P)
    Cd = c._concentrations('d')
    X = c._compute_char_conversion()
    dc = d_core_from_spm_char_conversion(X, c.solid.d_p)
    areas = c._calc_char_surface_area_per_class()
    r1_v, alpha = rate_R1(c.T, Cd[idx['O2']], c.solid.d_p, D_g, dc, 'coal')
    r2_v = rate_R2(c.T, Cd[idx['H2O']], c.solid.d_p, D_g, dc)
    r4_v = rate_R4_effective(c.T, max(Cd[idx['CO2']]*8.314*c.T,0), max(Cd[idx['CO']]*8.314*c.T,0), c.solid.d_p, D_g, dc)
    o2_sup = c.N_zu_d[idx['O2']] + c.N_d_in[idx['O2']] + c.N_zu_b[idx['O2']] + c.N_b_in[idx['O2']]
    h2o_sup = c.N_zu_d[idx['H2O']] + c.N_d_in[idx['H2O']] + c.N_zu_b[idx['H2O']] + c.N_b_in[idx['H2O']]
    print(f"cell{i}: T={c.T:.0f} area={areas.sum():.2e} r1={r1_v:.2e} r2={r2_v:.2e} r4={r4_v:.2e} o2={o2_sup:.3e} h2o={h2o_sup:.3e} m_char={c.m_solid[:,S_CHAR].sum():.5f} tau_d={c.V_d/max(c.u0,c.u_mf,1e-3):.2f}s")

print('\n--- 出口干基气相组成 ---')
dy = res['exit_gas_dry']
for sp in ['CO','CO2','H2','CH4','H2O','N2']:
    print(f"  {sp}: {dy.get(sp, 0)*100:.3f}%")
