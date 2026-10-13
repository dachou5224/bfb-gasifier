"""Recertify paired checkpoints and probe freeboard tar redistribution only."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import coupled_residual_named,outer_P_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state,progress_phi
from src.solvers.external_state_ledger import external_imbalance
r=Path(__file__).parent;old=Path('docs/audit_reports/v1_v3_review_2026-10-05')
w=json.loads((old/'w2_once.json').read_text())['s_0_16'];c=json.loads((old/'w2_from_008.json').read_text());nq=q_on_wesseling1_path(.16)
def assess(x):
 q=outer_P_named(x,nq);f=coupled_residual_named(x,q);cert=certify_controlled_state(x,f,q.hydro);l=external_imbalance(x,q)
 cho={k:l['gap_mol_s'][k] for k in ['C','H','O']}
 return dict(cert=cert,cho=cho,enthalpy_gap_W=l['enthalpy_W']['gap'],phi=progress_phi(f),x=x.tolist(),F=f.tolist())
x=np.array(c['x']);base=assess(x);f=np.array(base['F']);out={'warm':assess(np.array(w['x'])),'cold':base,'redistributions':[]}
for delta in [1e-7,1e-6,-f[44],(f[45]-f[44])/2]:
 z=x.copy();z[44]+=delta;z[45]-=delta;row=assess(z);row['delta']=float(delta);row['major_F_max_change']=float(np.max(np.abs(np.array(row['F'])[:40]-f[:40])));out['redistributions'].append(row)
(r/'evidence.json').write_text(json.dumps(out,indent=2))
for key in ['warm','cold']:print(key,out[key]['cert']['passed'],out[key]['cert']['gate_ratio'],out[key]['cho'])
for row in out['redistributions']:print('delta',row['delta'],'pass',row['cert']['passed'],'gates',row['cert']['gate_ratio'],'cho',row['cho'],'gas_change',row['major_F_max_change'])
