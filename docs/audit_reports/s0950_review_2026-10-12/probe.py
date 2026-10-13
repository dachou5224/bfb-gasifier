import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));t=json.load(open('docs/audit_reports/equation_loop_2026-10-08/advance_0.952.json'))
x=np.array(t['trial_x']);nq=q_on_wesseling1_path(t['s']);q=outer_P_named(x,nq);f=coupled_residual_named(x,q)
xr=np.array(r['x']);qr=outer_P_named(xr,q_on_wesseling1_path(r['s']));fr=coupled_residual_named(xr,qr)
scans={};roots={};cols={}
for i in [5,21]:
 def row(v):
  z=x.copy();z[i]=10.**v;return float(coupled_residual_named(z,q)[i])
 scans[i]=[[v,row(v)] for v in range(-70,-3,2)]
 for (a,fa),(b,fb) in zip(scans[i][:-1],scans[i][1:]):
  if fa*fb<0:
   lr=brentq(row,a,b,xtol=1e-12);roots[i]={'old':float(x[i]),'new':10**lr,'log_growth':float(np.log(10**lr/x[i]))};break
 h=1e-4;xp=x.copy();xm=x.copy();xp[i]*=np.exp(h);xm[i]*=np.exp(-h)
 c=(coupled_residual_named(xp,outer_P_named(xp,nq))-coupled_residual_named(xm,outer_P_named(xm,nq)))/(2*h)
 cols[i]={'self':float(c[i]),'norm':float(np.linalg.norm(c))}
z=x.copy()
for i,v in roots.items():z[i]=v['new']
fz=coupled_residual_named(z,outer_P_named(z,nq))
def plain(o):
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 raise TypeError(type(o))
out={'fingerprint':named_kernel_source_sha256(),'saved_fingerprint':r['source_sha256'],'root_s':r['s'],'root_cert':certify_controlled_state(xr,fr,qr.hydro),'root_ledger':external_imbalance(xr,qr),'trial_saved_F_difference':float(max(abs(f-np.array(t['F'])))),'oxygen_scans':scans,'frozen_scalar_roots':roots,'full_log_columns':cols,'reseed_only_F':fz.tolist(),'reseed_only_cert':certify_controlled_state(z,fz,outer_P_named(z,nq).hydro)}
Path(__file__).with_name('evidence.json').write_text(json.dumps(out,indent=2,default=plain))
print(json.dumps({k:v for k,v in out.items() if k not in ['oxygen_scans','reseed_only_F','root_ledger']},indent=2,default=plain))
