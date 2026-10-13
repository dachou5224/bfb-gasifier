"""Read-only fixed-state and bounded timeout audit."""
import json,runpy
from pathlib import Path
from dataclasses import replace
import numpy as np
from src.solvers.wesseling_two_zone import build_wesseling1_two_zone_q,wesseling1_physical_x0
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_control import reactivate_zero_unknowns
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q,_e8_log_scaled_correct
q0=build_wesseling1_two_zone_q();x=wesseling1_physical_x0(q0);q=outer_P_named(x,q0)
v=q0.reference.area*q0.reference.H_bed
cold=runpy.run_path('docs/audit_reports/p0_p4_system_2026-09-26/p3_cold_start.py')['_initial']()
cq=outer_P_named(cold,replace(freeze_e8_q(),stage=11));z,opened=reactivate_zero_unknowns(cq,cold);f=coupled_residual_named(z,cq)
stats={}
try:_e8_log_scaled_correct(cq,z,deadline_s=0.02,stats=stats)
except TimeoutError:status='timeout'
else:status='returned'
d={'wesseling':{'bed_geometric_volume_m3':v,'V_d':q0.energy.species.geo_ex.V_d,'V_b':q0.energy.species.geo_ex.V_b,'sum_phase_to_bed_volume':(q0.energy.species.geo_ex.V_d+q0.energy.species.geo_ex.V_b)/v,'w_undecided':q.hydro.w_undecided.tolist(),'tau_s':q.hydro.tau_s.tolist(),'effective_tau_floor_s':np.maximum(q.hydro.tau_s,1e-6).tolist()},'cold':{'opened':opened,'remaining_zero_rows':[{'index':int(i),'F':float(f[i])} for i in np.flatnonzero(z==0)]},'timeout':{'status':status,'stats_after_call':stats}}
p=Path('docs/audit_reports/wp1_wp4_review_2026-09-28/evidence.json');p.write_text(json.dumps(d,indent=2));print(p.read_text())
