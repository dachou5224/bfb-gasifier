"""Deterministic protocol probe of continuation; no kinetic solve."""
from unittest.mock import patch
import json
from pathlib import Path
import src.physics.wirsum_1997_freeboard as fb
calls=[]
def initial(*args,**kwargs):
    calls.append(kwargs['cell_height'])
    return ('undecided',5.0,0.9,0.2)
def continued(*args,**kwargs):
    return ('reaches_next_cell',1.0,True)
with patch.object(fb,'integrate_particle_fate_with_time',initial),patch.object(fb,'continue_undecided_event',continued):
    result=fb.integrate_particle_through_cells(1.0,[1.0,1.0],u_g=.3,rho_g=1.,rho_p=1000.,d_p=.001,mu_g=1e-5,z=1.,t_max=5.,continue_budget_s=15.)
out={'synthetic_events':['initial integration undecided in cell 0','refined continuation reaches internal upper face at +1 s'],'actual_result':result,'cell_integrations':len(calls),'finding':'No integration in cell 1 follows the accepted internal-face event. This is a protocol reproduction, not a Wesseling trajectory result.'}
Path('docs/audit_reports/review_abcd_2026-09-28/evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
