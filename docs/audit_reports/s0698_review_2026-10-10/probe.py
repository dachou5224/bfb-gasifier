"""Cold checks and a frozen cyclone block calculation, not a coupled solve."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,char_face_flows_kg_s
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
root=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/advance_0.700.json'));nq=q_on_wesseling1_path(.7);x=np.array(r['trial_x']);q=outer_P_named(x,nq);f=coupled_residual_named(x,q)
ka=float(q.energy.species.k.K_auf);kb=float(q.energy.species.k.K_ab);leave=char_face_flows_kg_s(x,q)[2];A=np.array([[ka+kb,-kb],[-ka,ka+kb]]);mass=np.linalg.solve(A,[sum(leave),0]);z=x.copy();z[37:39]=mass;fz=coupled_residual_named(z,outer_P_named(z,nq))
columns=[]
for i in (37,38):
 h=1e-7;xp=x.copy();xp[i]+=h;fp=coupled_residual_named(xp,outer_P_named(xp,nq));physical=(fp-f)/h;columns.append({'i':i,'physical_rows37_38':physical[37:39].tolist(),'log_exact_rows37_38':(physical[37:39]*x[i]).tolist()})
xr=np.array(root['x']);qr=outer_P_named(xr,q_on_wesseling1_path(root['s']));fr=coupled_residual_named(xr,qr)
out={'source_sha256':named_kernel_source_sha256(),'saved_root_source':root['source_sha256'],'root_cert':certify_controlled_state(xr,fr,qr.hydro),'root_ledger':external_imbalance(xr,qr),'root_w_leave':qr.hydro.w_leave.tolist(),'root_M_cyclone':xr[37:39].tolist(),'trial_saved_F_max_difference':float(max(abs(f-np.array(r['F'])))),'trial_w_leave':q.hydro.w_leave.tolist(),'trial_tau':q.hydro.tau_s.tolist(),'trial_solid_rows':{i:float(f[i]) for i in (14,15,37,38,39,46,47)},'cyclone':{'ka':ka,'kb':kb,'incoming_kg_s':float(sum(leave)),'old_mass':x[37:39].tolist(),'frozen_equilibrium_mass':mass.tolist(),'log_growth_needed':np.log(mass/x[37:39]).tolist(),'columns':columns},'candidate_F_solid':{i:float(fz[i]) for i in (14,15,37,38,39,46,47)},'candidate_energy':fz[-2:].tolist(),'candidate_is_root':False}
def plain(o):
 if isinstance(o,dict):return {str(k):plain(v) for k,v in o.items()}
 if isinstance(o,(list,tuple)):return [plain(v) for v in o]
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 return o
Path(__file__).with_name('evidence.json').write_text(json.dumps(plain(out),indent=2));print(json.dumps(plain({k:v for k,v in out.items() if k!='root_ledger'}),indent=2))
