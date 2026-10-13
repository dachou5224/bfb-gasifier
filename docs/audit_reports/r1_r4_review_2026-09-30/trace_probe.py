import json
from pathlib import Path
from unittest.mock import patch
import src.physics.wirsum_1997_freeboard as fb
from src.solvers.wesseling_two_zone import build_wesseling1_two_zone_q,wesseling1_physical_x0
orig=fb.continue_undecided_event
integ=fb.integrate_particle_fate_with_time
rows=[]
class Done(Exception):pass
def wrap(*args,**kwargs):
    if kwargs['d_p']<.00029:return orig(*args,**kwargs)
    ends=[]
    def logged(*a,**kw):
        r=integ(*a,**kw);ends.append({'dt':kw['dt'],'fate':r[0],'elapsed':r[1],'height':r[2],'velocity':r[3]});return r
    with patch.object(fb,'integrate_particle_fate_with_time',logged):r=orig(*args,**kwargs)
    rows.append({'initial_u_h':args[:2],'cell_height':kwargs['cell_height'],'extra_s':kwargs['extra_s'],'results':ends,'accepted':r[-1]})
    if not r[-1]:raise Done()
    return r
with patch.object(fb,'continue_undecided_event',wrap):
    try:wesseling1_physical_x0(build_wesseling1_two_zone_q())
    except Done:pass
out=Path('docs/audit_reports/r1_r4_review_2026-09-30/trace_evidence.json');out.write_text(json.dumps(rows,indent=2));print(out.read_text())
