import json,numpy as np
from pathlib import Path
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,daem_X_at_T,coupled_residual_named
from src.solvers.wirsum_1997_pyrolysis_ledger import map_one_state
from src.solvers.external_state_ledger import external_imbalance
from src.core.elemental_ledger import H2O_MOLAR_MASS_KG_PER_MOL
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
p=json.load(open('docs/audit_reports/r1_r4_review_2026-09-30/s3_path_step.json'))
rows=[]
def visit(o):
 if isinstance(o,dict):
  if 's' in o and 'x' in o and 'F' in o: rows.append(o)
  for v in o.values():visit(v)
 elif isinstance(o,list):
  for v in o:visit(v)
visit(p)
print('saved',[(r['s'],r.get('passed')) for r in rows])
out=[]
for s in [.02,.22,.23]:
 matches=[r for r in rows if abs(float(r['s'])-s)<1e-8]
 if not matches:continue
 r=matches[-1];x=np.array(r['x']);nq=q_on_wesseling1_path(s);q=outer_P_named(x,nq)
 m=map_one_state(q.energy.feed,P_Pa=q.reference.P,T_K=x[-2],X_vm=np.array([0.,daem_X_at_T(x[-2],q.tau_daem_s)]),tau_s=np.array([q.tau_daem_s]),occupancy_source=q.energy.occupancy_source)
 ledger=external_imbalance(x,q);prod=m['cells'][0]['n_py_mol_s'];f=coupled_residual_named(x,q)
 out.append({'s':s,'feed_water_mol_s':q.energy.feed.m_moist_kg_s/H2O_MOLAR_MASS_KG_PER_MOL,'equation_water_mol_s':float(q.energy.species.n_moist[3]),'blocks':named_endpoint_blocks(f),'tar_rows':f[40:46].tolist(),'tar_x':x[40:46].tolist(),'pyro_atom_gaps':{k:v for k,v in m.items() if k.startswith('atom_')},'omitted_minor_H':3*prod.get('NH3',0)+2*prod.get('H2S',0),'external_gap':ledger['gap_mol_s'],'inlet_total':ledger['inlet_mol_s']['total']})
end=[]
for s in [0,.999999,1]:
 q=q_on_wesseling1_path(s);end.append({'s':s,'N_in_d':q.energy.species.N_in_d.tolist(),'V_d':q.reference.V_d,'V_b':q.reference.V_b,'trajectory_continue_s':q.reference.trajectory_continue_s})
result={'states':out,'path_endpoints':end}
Path('docs/audit_reports/s1_s3_review_2026-09-30/evidence.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
