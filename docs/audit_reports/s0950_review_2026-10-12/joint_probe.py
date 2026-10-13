"""One bounded joint solve from a diagnostically reseeded state; no official checkpoint writes."""
import json,time
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_control import solve_named_controlled,certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
p=Path(__file__).parent;e=json.loads((p/'evidence.json').read_text());t=json.load(open('docs/audit_reports/equation_loop_2026-10-08/advance_0.952.json'));x=np.array(t['trial_x'])
for i,v in e['frozen_scalar_roots'].items():x[int(i)]=v['new']
nq=q_on_wesseling1_path(t['s']);tm=time.monotonic();r=solve_named_controlled(nq,x,deadline_s=20,temperature_radius_scale=4);elapsed=time.monotonic()-tm
print('keys',list(r));
def plain(o):
 if isinstance(o,np.ndarray):return o.tolist()
 if isinstance(o,np.generic):return o.item()
 raise TypeError(type(o))
# Return schema is inspected before further certification.
out={k:v for k,v in r.items() if k not in ['q']};out['audit_elapsed_s']=elapsed
(p/'joint_evidence.json').write_text(json.dumps(out,indent=2,default=plain))
print('elapsed',elapsed)
for k,v in r.items():
 if k not in ['x','q','F','history']:print(k,str(v)[:1000])
xc=np.array(r['x']);qc=outer_P_named(xc,nq);fc=coupled_residual_named(xc,qc)
out['cold_cert']=certify_controlled_state(xc,fc,qc.hydro);out['cold_ledger']=external_imbalance(xc,qc)
(p/'joint_evidence.json').write_text(json.dumps(out,indent=2,default=plain))
