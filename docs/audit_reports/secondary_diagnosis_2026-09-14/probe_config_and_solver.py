import json,numpy as np
from tests.global_nr_fixtures import build_initialized_lu_reactor
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton,apply_wirsum_1997_strict_config
from src.solvers.wirsum_1997_startwert_r12 import apply_same_phase_r12_startwert_to_cells
from src.solvers.wirsum_1997_holdup_startwert import apply_km_holdup_startwert_axial_sweep
from src.solvers.vorabrechnung.vorab_policy import vorab_outer_may_refresh_drying_pyro_sources
from src.solvers.global_nr_solver import global_residual,pack_reactor
from src.core.cell import S_MOISTURE
import src.core.cell_balances as cb
r=build_initialized_lu_reactor(3);r._apply_all_bc_for_nr();[c.residuals() for c in r.cells]
apply_same_phase_r12_startwert_to_cells(r.cells);apply_km_holdup_startwert_axial_sweep(r.cells,r._apply_all_bc_for_nr)
[c.residuals() for c in r.cells]
before=[c.calc_energy_balance() for c in r.cells]; orig=cb.calc_solid_enthalpy_flow;charhf=cb.H_FORMATION_CHAR_J_PER_KG
cb.H_FORMATION_CHAR_J_PER_KG=0.
def diag_solid(m,T,**kw):
 return orig(m,T,**kw)+float(np.sum(np.asarray(m)[:,S_MOISTURE]))*(-285830/.01801528)
cb.calc_solid_enthalpy_flow=diag_solid
after=[c.calc_energy_balance() for c in r.cells]
cb.calc_solid_enthalpy_flow=orig;cb.H_FORMATION_CHAR_J_PER_KG=charhf
cfg=build_phase2_htw_lu_freeboard_reactor_config(); phase2={k:getattr(cfg,k) for k in ['thesis_mode','thesis_vorab_sources_single_shot','nr_solver_mode','holdup_k0_anchor','nr_refresh_hydrodynamics_on_accepted_step_thesis']};phase2['outer_refresh_pyro_allowed']=vorab_outer_may_refresh_drying_pyro_sources(cfg)
apply_wirsum_1997_strict_config(cfg);phase2['after_strict_outer_refresh_pyro_allowed']=vorab_outer_may_refresh_drying_pyro_sources(cfg)
J=np.diag([1e16,1.]);res=wirsum_1997_strict_newton(lambda x:J@x-np.ones(2),lambda x:J,np.zeros(2))
false=wirsum_1997_strict_newton(lambda x:1e22*x+1e35*x*x+1e8,lambda x:np.diag(1e22+2e35*x),np.zeros(1))
out={'energy_same_state_diagnostic_only':{'before_W':before,'char0_water_liquid_hf_W':after,'delta_W':(np.array(after)-before).tolist()},'phase2_config':phase2,'strict_rank_diagonal':{'failure':res.failure_reason,'rank':res.rank_last,'true_det':float(np.linalg.det(J))},'strict_small_step':{'converged':false.converged,'F':false.F.tolist(),'x':false.x.tolist(),'iters':false.iterations}}
print(json.dumps(out,indent=2));open('/tmp/bfb_secondary_followup.json','w').write(json.dumps(out,indent=2))
