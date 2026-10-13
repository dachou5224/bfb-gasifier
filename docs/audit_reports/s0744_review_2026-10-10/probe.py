"""Read-only equilibrium accuracy probe; no continuation or production edits."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from src.solvers import vorabrechnung
import src.solvers.vorabrechnung.vorab_gibbs_seed as g
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.core.constants import Rg,P0
saved=json.load(open('docs/audit_reports/equation_loop_2026-10-08/advance_0.746.json'))
x=np.array(saved['trial_x']);nq=q_on_wesseling1_path(.746)
orig=g._solve_pyrolysis_gibbs_an_eq_b
records=[]
use_refined=False
def wrapper(**kw):
 out,diag=orig(**kw); species=list(g._PYROLYSIS_GIBBS_CANDIDATES); db=g._get_major_gibbs_reduced_solver();A=db._build_atom_matrix(list(g._PYROLYSIS_GIBBS_ELEM_NAMES),species);Z=null_space(A)
 n=np.array([out[s] for s in species]);n0=n.copy();T=kw['T_K'];c=np.array([db.species_db.get_mu0(s,T) for s in species])/(Rg*T)+np.log(kw['P_Pa']/P0)
 def grad(n):return Z.T@(c+np.log(n/n.sum()))
 before=grad(n)
 for it in range(20):
  v=grad(n)
  if max(abs(v))<1e-12:break
  H=Z.T@(np.diag(1/n)-np.ones((len(n),len(n)))/n.sum())@Z
  dn=Z@np.linalg.solve(H,-v);alpha=min(1.,min((-0.9*n[dn<0]/dn[dn<0]),default=1.))
  while np.linalg.norm(grad(n+alpha*dn))>np.linalg.norm(v) and alpha>1e-8:alpha*=.5
  n+=alpha*dn
 records.append({'T':T,'species':species,'n':n0.tolist(),'refined_n':n.tolist(),'correction_L2':float(np.linalg.norm(n-n0)),'stationarity_before':before.tolist(),'stationarity_after':grad(n).tolist(),'diag':diag})
 return ({sp:float(n[i]) for i,sp in enumerate(species)} if use_refined else out),diag
g._solve_pyrolysis_gibbs_an_eq_b=wrapper
for dt in [0,1e-6,1e-5,5e-5,1e-4,.01,.1]:
 z=x.copy();z[-2]+=dt;outer_P_named(z,nq)
g._solve_pyrolysis_gibbs_an_eq_b=orig
r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));xr=np.array(r['x']);qr=outer_P_named(xr,q_on_wesseling1_path(r['s']));fr=coupled_residual_named(xr,qr)
scan_records=list(records)
use_refined=True
g._solve_pyrolysis_gibbs_an_eq_b=wrapper
qr2=outer_P_named(xr,q_on_wesseling1_path(r['s']));fr2=coupled_residual_named(xr,qr2)
qt2=outer_P_named(x,nq);ft2=coupled_residual_named(x,qt2)
g._solve_pyrolysis_gibbs_an_eq_b=orig
result={'refined_root_cert':certify_controlled_state(xr,fr2,qr2.hydro),'refined_trial_cert':certify_controlled_state(x,ft2,qt2.hydro),'fingerprint':named_kernel_source_sha256(),'saved_fingerprint':r['source_sha256'],'root_s':r['s'],'root_cert':certify_controlled_state(xr,fr,qr.hydro),'samples':scan_records}
Path(__file__).with_name('evidence.json').write_text(json.dumps(result,indent=2,default=lambda o:o.item() if isinstance(o,np.generic) else o.tolist()))
for a in records:print(a['T'],a['correction_L2'],a['stationarity_before'],a['diag']['opt_success'],a['diag']['message'])
print('root',result['root_cert'])
