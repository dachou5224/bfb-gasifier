import json, copy, math
import numpy as np
from tests.global_nr_fixtures import build_initialized_lu_reactor
from src.solvers.global_nr_solver import global_residual,pack_reactor
from src.solvers.wirsum_1997_holdup_startwert import apply_km_holdup_startwert_axial_sweep
from src.core.species import GAS_SPECIES_INDEX as I,enthalpy_molar,cp_char,cp_ash
from src.core.cell import S_CHAR,S_VM,S_MOISTURE,S_ASH
from src.core.energy_balance_defaults import H_FORMATION_CHAR_J_PER_KG
from src.core.connectivity.solid_transport import cell_total_solid_holdup,cell_resident_holdup_for_transport_closure,update_bed_solid_transport_coefficients
from src.kinetics.gas_reactions import rate_R12
from src.kinetics.arrhenius import k_jensen_h2
from src.solvers.vorabrechnung.vorab_policy import vorab_outer_may_refresh_drying_pyro_sources
r=build_initialized_lu_reactor(3); c=r.cells[0]; r._apply_all_bc_for_nr(); [v.residuals() for v in r.cells]
out={'case':'current worktree reduced3; no Newton','config':{k:getattr(r.config,k) for k in ['thesis_mode','thesis_vorab_sources_single_shot','nr_solver_mode','holdup_k0_anchor']}}
out['outer_may_refresh_sources']=vorab_outer_may_refresh_drying_pyro_sources(r.config)
thermo=c._get_local_thermo_bundle(); a=float(thermo['C_d'][I['H2']]); b=float(thermo['C_d'][I['O2']]); T=c.T
old=rate_R12(T,a,b,c.P,thermo['y_d']); literal=k_jensen_h2(1.631e9,3420.,T)*a**1.5*b
out['R12']={'T':T,'CH2':a,'CO2_oxygen':b,'current_rate':old,'pdf_literal_rate':literal,'literal_over_current':literal/old,'Vd':c.V_d}
# Standard reaction enthalpies implied by actual species properties.
hc=H_FORMATION_CHAR_J_PER_KG*.012011
h=lambda s:enthalpy_molar(s,298.15)
out['thermochemistry']={'char_hf_kJ_mol':hc/1000,'C_to_CO2_current_kJ_mol':(h('CO2')-h('O2')-hc)/1000,'C_plus_H2O_to_CO_H2_current_kJ_mol':(h('CO')+h('H2')-h('H2O')-hc)/1000,'C_plus_CO2_to_2CO_current_kJ_mol':(2*h('CO')-h('CO2')-hc)/1000,'water_vapor_hf_kJ_mol':h('H2O')/1000,'liquid_water_hf_J_kg':-285830/.01801528}
m=np.zeros_like(c.m_solid); m[0,S_CHAR]=1
out['thermochemistry']['char_1kg_hf_J']=c._calc_solid_enthalpy_flow(m,298.15)
m[:]=0;m[0,S_MOISTURE]=.01801528
out['thermochemistry']['liquid_water_1mol_hf_J']=c._calc_solid_enthalpy_flow(m,298.15)
w_ash=c.solid.ash_dry_wt/100; w_vm=c.solid.VM_daf/100*(1-w_ash);w_char=1-w_ash-w_vm
m[:]=0;m[0]=[w_char,w_vm,0,w_ash]
out['thermochemistry'].update({'dryfuel_hf_from_HHV_J_kg':c.solid.h_f_dry,'dryfuel_components_hf_J_kg':c._calc_solid_enthalpy_flow(m,298.15),'dryfuel_mismatch_J_kg':c._calc_solid_enthalpy_flow(m,298.15)-c.solid.h_f_dry,'char_fraction_dry':w_char})
apply_km_holdup_startwert_axial_sweep(r.cells,r._apply_all_bc_for_nr)
x=pack_reactor(r.cells); fs=[];ks=[]
for _ in range(3):
 fs.append(global_residual(x,r.cells,r._apply_all_bc_for_nr).copy());ks.append(np.stack([v.K_solid_ab.copy() for v in r.cells]))
out['seeded_repeat']={'F_linf_deltas':[float(np.max(abs(fs[i]-fs[i-1]))) for i in [1,2]],'K_linf_deltas':[float(np.max(abs(ks[i]-ks[i-1]))) for i in [1,2]]}
# interleaved evaluate perturbed T then return x
xx=x.copy();xx[-1]+=10
_ = global_residual(xx,r.cells,r._apply_all_bc_for_nr)
fback=global_residual(x,r.cells,r._apply_all_bc_for_nr)
out['seeded_repeat']['F_x_y_x_linf_delta']=float(np.max(abs(fback-fs[-1])))
out['inventories']=[{'M_geo':cell_total_solid_holdup(v),'M_active':float(v.m_solid[:,[S_CHAR,S_ASH]].sum()),'closure_M':cell_resident_holdup_for_transport_closure(v,is_bottom_bed_cell=i==0)} for i,v in enumerate(r.cells)]
# freeze physics and transport at point, then repeat order check
for v in r.cells:
 v.compute_vorabrechnung(v.geo.dh/max(v.u_mf,.001));v._snapshot_vorabrechnung_hydrodynamics();v._freeze_vorabrechnung_inner_nr=True
f0=global_residual(x,r.cells,r._apply_all_bc_for_nr);_=global_residual(xx,r.cells,r._apply_all_bc_for_nr);f1=global_residual(x,r.cells,r._apply_all_bc_for_nr)
out['frozen_F_x_y_x_linf_delta']=float(np.max(abs(f1-f0)))
# Bottom inventory branch one-sided limits
v=r.cells[0]; geo=cell_total_solid_holdup(v)
branches=[]
for f in [.149999,.15,.150001]:
 v.m_solid[:]=0;v.m_solid[0,S_CHAR]=f*geo
 branches.append({'fraction':f,'closure_mass':cell_resident_holdup_for_transport_closure(v,is_bottom_bed_cell=True)})
out['bottom_inventory_branch']=branches
out['R12']['O2_concentration_for_1mol_s_current']=(1/(k_jensen_h2(1.631e9,3420.,T)*a*out['R12']['Vd']))**2
out['R12']['O2_concentration_for_1mol_s_pdf_literal']=1/(k_jensen_h2(1.631e9,3420.,T)*a**1.5*out['R12']['Vd'])
print(json.dumps(out,indent=2));open('/tmp/bfb_secondary_evidence.json','w').write(json.dumps(out,indent=2))
