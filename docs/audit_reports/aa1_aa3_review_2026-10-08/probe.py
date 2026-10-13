"""Read-only cold recertification; no continuation."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wirsum_1997_named_f import coupled_residual_named,outer_P_named
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
r=json.loads(Path('docs/audit_reports/z1_z3_review_2026-10-08/aa1_aa2.json').read_text()); rows=[]
for p in r['aa2']['points']:
 x=np.array(p['x']);q=outer_P_named(x,q_on_wesseling1_path(p['s']));f=coupled_residual_named(x,q);c=certify_controlled_state(x,f,q.hydro);led=external_imbalance(x,q)
 rows.append({'s':p['s'],'cert':c,'cho':{k:led['gap_mol_s'][k] for k in ('C','H','O')},'enthalpy_gap_W':led['enthalpy_W']['gap'],'saved_row_F_max_diff':float(np.max(np.abs(f-np.array(p['F'])))),'temperatures_K':x[-2:].tolist()})
last=rows[-1];x=np.array(r['certified_x']);f=coupled_residual_named(x,outer_P_named(x,q_on_wesseling1_path(r['certified_s'])))
out={'current_source_sha256':named_kernel_source_sha256(),'rows':rows,'top_level_certified_F_max_diff':float(np.max(np.abs(f-np.array(r['certified_F'])))),'saved_top_level_energy':r['certified_F'][-2:],'recomputed_energy':f[-2:].tolist(),'top_level_F_matches_AA1_energy':bool(np.allclose(r['certified_F'][-2:],r['aa1']['energy_rows_W'],rtol=0,atol=1e-12))}
Path(__file__).with_name('evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
