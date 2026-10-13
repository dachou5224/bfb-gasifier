"""Cold endpoint and default initial-state audit; no coupled solve."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,default_x0_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.named_kernel_sources import named_kernel_source_sha256
r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));q=q_on_wesseling1_path(1.);x=np.array(r['x']);p=outer_P_named(x,q);f=coupled_residual_named(x,p);x0=default_x0_named(q);p0=outer_P_named(x0,q);f0=coupled_residual_named(x0,p0)
out={'source_sha256':named_kernel_source_sha256(),'saved_source_sha256':r['source_sha256'],'saved_s':r['s'],'cert':certify_controlled_state(x,f,p.hydro),'ledger':external_imbalance(x,p),'default_x':x0.tolist(),'default_cert':certify_controlled_state(x0,f0,p0.hydro),'root_x':x.tolist(),'default_largest_rows':sorted(enumerate(f0.tolist()),key=lambda a:abs(a[1]),reverse=True)[:10]}
Path(__file__).with_name('evidence.json').write_text(json.dumps(out,indent=2,default=lambda o:o.tolist() if isinstance(o,np.ndarray) else o.item()))
print('fingerprints',out['source_sha256'],out['saved_source_sha256']);print('cert',out['cert']);print('gaps',out['ledger']['gap_mol_s']);print('initial cert',out['default_cert']);print('largest',out['default_largest_rows'])
print('root vs default',[(i,float(x[i]),float(x0[i])) for i in [14,15,37,38,39,46,47,48,49]])
