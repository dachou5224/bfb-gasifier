"""Read-only diagnosis; no coupled solve."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_p2 import gas_velocity_and_density
from src.physics.wirsum_1997_freeboard import hamel_eq386_asymptotic_z
from src.core.constants import g
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
root=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));fail=json.load(open('docs/audit_reports/equation_loop_2026-10-08/advance_0.698.json'))
x=np.array(root['x']);target=fail['s'];nq=q_on_wesseling1_path(target)
def state(s,z):
 q=outer_P_named(z,q_on_wesseling1_path(s));f=coupled_residual_named(z,q);b=q.reference
 ug,rho,_,_=gas_velocity_and_density(z[16:23],np.zeros(7),T=z[-1],P=b.P,area=b.area)
 d=b.d_p[1];km=18*b.mu_g/(b.rho_p*d*d);kz=.75*hamel_eq386_asymptotic_z(b.phi_s)*rho/(b.rho_p*d);ge=g*(1-rho/b.rho_p);vt=2*ge/(km+np.sqrt(km*km+4*kz*ge));tau=q.hydro.tau_s[1]
 return {'s':s,'ug':ug,'rho':rho,'vt':vt,'margin':vt-ug,'tau':tau,'escape_rate':1/tau,'w_leave':q.hydro.w_leave.tolist(),'T':z[-2:].tolist(),'coarse_M':z[47],'coarse_F':f[47],'cert':certify_controlled_state(z,f,q.hydro)}
rows={'root':state(root['s'],x),'target_at_old_x':state(target,x),'prediction':state(target,np.array(fail['prediction_x'])),'trial':state(target,np.array(fail['trial_x']))}
derivs=[]
for i in (49,16,19):
 h=1e-6;xp=x.copy();xm=x.copy();xp[i]*=np.exp(h);xm[i]*=np.exp(-h);p=state(root['s'],xp);m=state(root['s'],xm)
 derivs.append({'i':i,**{k:(p[k]-m[k])/(2*h) for k in ('ug','vt','margin','tau','escape_rate')}})
hist=fail['history'];summary={'n':len(hist),'accepted':sum(bool(r.get('accepted')) for r in hist),'first':hist[0],'last':hist[-1]}
q=outer_P_named(x,q_on_wesseling1_path(root['s']));out={'source_sha256':named_kernel_source_sha256(),'states':rows,'log_derivatives':derivs,'history':summary,'ledger':external_imbalance(x,q)}
def plain(o):
 if isinstance(o,dict):return {str(k):plain(v) for k,v in o.items()}
 if isinstance(o,(tuple,list)):return [plain(v) for v in o]
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 return o
Path(__file__).with_name('evidence.json').write_text(json.dumps(plain(out),indent=2));print(json.dumps(plain({k:v for k,v in out.items() if k!='ledger'}),indent=2))
