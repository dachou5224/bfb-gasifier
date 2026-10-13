"""Recreate two saved-target predictors only; do not run correctors."""
import json,time
from pathlib import Path
import numpy as np
from src.solvers.wirsum_1997_named_f import coupled_residual_named,outer_P_named
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wirsum_1997_named_f_control import progress_ratios
r=json.loads(Path('docs/audit_reports/z1_z3_review_2026-10-08/aa1_aa2.json').read_text());points=r['aa2']['points'];out=[]
for ix in (2,4):
 prev=points[ix-1];target=points[ix];s=prev['s'];x=np.array(prev['x']);h=1e-4
 def full(t):return coupled_residual_named(x,outer_P_named(x,q_on_wesseling1_path(t)))
 offset=(full(s+h)-full(s-h))/(2*h)*target['ds'];stats={};t=time.monotonic()
 z=_e8_log_scaled_correct(outer_P_named(x,q_on_wesseling1_path(s)),x,log_span=.5,ratio_span=1,stats=stats,residual_offset=offset,deadline_s=15)
 f=coupled_residual_named(z,outer_P_named(z,q_on_wesseling1_path(target['s'])))
 out.append({'s':target['s'],'wall_s':time.monotonic()-t,'solved':stats.get('subproblem_solved'),'predictor_temperatures':z[-2:].tolist(),'root_minus_predictor_K':(np.array(target['x'])[-2:]-z[-2:]).tolist(),'gate_ratios':progress_ratios(f),'predictor_x':z.tolist()})
Path(__file__).with_name('predict_evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps([{k:v for k,v in p.items() if k!='predictor_x'} for p in out],indent=2))
