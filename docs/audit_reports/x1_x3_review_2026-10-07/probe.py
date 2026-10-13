"""Fixed-point objective consistency and path smoothness; no path solves."""
import json
from pathlib import Path
import numpy as np
import scipy.optimize as opt
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,original_residual_from_tar_rows
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct,residual_block_scale
from src.solvers.wirsum_1997_named_f_control import progress_phi,progress_ratios,step_improves_progress,certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
r=Path(__file__).parent;p=Path('docs/audit_reports/w1_w3_review_2026-10-07');a=json.loads((p/'x2_s032.json').read_text());b=json.loads((p/'x2_s048.json').read_text())
x=np.array(a['x']);smooth=[]
for s in [.319,.3195,.32,.3205,.321]:
 q=outer_P_named(x,q_on_wesseling1_path(s));f=coupled_residual_named(x,q)
 smooth.append({'s':s,'E':f[-2:].tolist(),'phi':progress_phi(f),'fuel_kg_s':q.energy.feed.m_dry_kg_s,'inlet_O2':float(q.energy.species.N_in_b[5])})
q=outer_P_named(x,q_on_wesseling1_path(.32));f=coupled_residual_named(x,q);basecert=certify_controlled_state(x,f,q.hydro);ledger=external_imbalance(x,q)
x=np.array(b['x']);q=outer_P_named(x,q_on_wesseling1_path(.48));f=coupled_residual_named(x,q);sc=residual_block_scale(50);gate=np.full(50,1e-6);gate[-2:]=1
orig=opt.lsq_linear;capture={}
def inverse_matrix(A):return np.column_stack([original_residual_from_tar_rows(A[:,j]*sc)/gate for j in range(A.shape[1])])
rows=[];cache={}
for mode in ['current','gate_aligned','gate_energy_tangent']:
 def fn(A,b,**kw):
  if mode in ['gate_aligned','gate_energy_tangent']:
   A=inverse_matrix(A);b=original_residual_from_tar_rows(b*sc)/gate
  if mode=='gate_energy_tangent':
   logs=[i for i in np.flatnonzero(x>0) if i not in range(40,46)]
   ti=[logs.index(48),logs.index(49)];ri=[j for j in range(A.shape[1]) if j not in ti]
   V=np.zeros((A.shape[1],len(ri)));V[ri,:]=np.eye(len(ri));V[ti,:]=-np.linalg.solve(A[-2:,ti],A[-2:,ri])
   lo,hi=kw['bounds'];kw2=dict(kw);kw2['bounds']=(lo[ri],hi[ri])
   fit=orig(A@V,b,**kw2);v=V@fit.x
   factor=min([1.]+[hi[j]/v[j] if v[j]>0 else lo[j]/v[j] for j in ti if v[j]!=0])
   fit.x=v*factor
  else:fit=orig(A,b,**kw)
  capture['cost']=float(np.linalg.norm(A@fit.x-b));return fit
 opt.lsq_linear=fn
 try:z=_e8_log_scaled_correct(q,x,log_span=.001953125,ratio_span=.001953125,stats=cache)
 finally:opt.lsq_linear=orig
 ff=coupled_residual_named(z,outer_P_named(z,q_on_wesseling1_path(.48)));pred=np.array(cache['predicted_raw'])
 rows.append({'mode':mode,'phi0':progress_phi(f),'predicted_phi':progress_phi(pred),'actual_phi':progress_phi(ff),'ratios':progress_ratios(ff),'accepted_by_progress_rule':step_improves_progress(f,ff),'linear_status':cache['linear_status'],'physical_step_max':float(np.max(np.abs(z-x)))})
out={'s032_recert':basecert,'s032_cho':{k:ledger['gap_mol_s'][k] for k in ['C','H','O']},'smoothness_at_fixed_s032_x':smooth,'s048_same_matrix_comparison':rows}
(r/'evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
