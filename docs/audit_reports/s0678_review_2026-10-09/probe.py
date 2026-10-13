"""Cold recertification and actual corrector Jacobian columns. No solve."""
import json,time
from pathlib import Path
import numpy as np
import src.solvers.wirsum_1997_named_f_e8 as e8
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,original_residual_from_tar_rows
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));x=np.array(r['x']);nq=q_on_wesseling1_path(r['s']);t=time.monotonic();q=outer_P_named(x,nq);f=coupled_residual_named(x,q);cap={};orig=e8._gate_progress_step
# Capture the matrix exactly where the production corrector would solve it.
def grab(matrix,entrance_row,scale,lo,hi,quality=None,**kwargs):
 cap.update(matrix=matrix.copy(),scale=scale.copy());return np.zeros(matrix.shape[1]),original_residual_from_tar_rows(entrance_row*scale),0
e8._gate_progress_step=grab
try:e8._e8_log_scaled_correct(q,x,stats={},temperature_radius_scale=4,deadline_s=15)
finally:e8._gate_progress_step=orig
z0=e8._seed_downstream_o2(q,x.copy());ids=[int(i) for i in e8._e8_positive_log_ids(z0) if i not in range(40,46)];rows=[]
for i in (48,49,16,19):
 j=ids.index(i);used=original_residual_from_tar_rows(cap['matrix'][:,j]*cap['scale']);vals=[]
 for h in (1e-5,1e-6):
  xp=z0.copy();xm=z0.copy();xp[i]*=np.exp(h);xm[i]*=np.exp(-h);qp=outer_P_named(xp,nq);qm=outer_P_named(xm,nq)
  deriv=(coupled_residual_named(xp,qp)-coupled_residual_named(xm,qm))/(2*h)
  tauprime=(np.array(qp.hydro.tau_s)-np.array(qm.hydro.tau_s))/(2*h)
  vals.append({'h':h,'full_coarse_fb':deriv[47],'used_coarse_fb':used[47],'missing_coarse_fb':deriv[47]-used[47],'tau_log_derivative':tauprime.tolist(),'tau_term_expected':float(x[47]/q.hydro.tau_s[1]**2*tauprime[1]),'full_energy':deriv[-2:].tolist(),'used_energy':used[-2:].tolist(),'max_raw_difference':float(np.max(np.abs(deriv-used)))})
 rows.append({'state_index':i,'samples':vals})
out={'source_sha256':named_kernel_source_sha256(),'s':r['s'],'cert':certify_controlled_state(x,f,q.hydro),'ledger':external_imbalance(x,q),'tau':np.array(q.hydro.tau_s).tolist(),'columns':rows,'wall_s':time.monotonic()-t}
def plain(o):
 if isinstance(o,dict):return {str(k):plain(v) for k,v in o.items()}
 if isinstance(o,(list,tuple)):return [plain(v) for v in o]
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 return o
Path(__file__).with_name('evidence.json').write_text(json.dumps(plain(out),indent=2));print(json.dumps(plain({k:v for k,v in out.items() if k!='ledger'}),indent=2))
