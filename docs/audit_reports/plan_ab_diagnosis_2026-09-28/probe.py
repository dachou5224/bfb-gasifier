"""Bounded audit only: one 10 s frozen-q correction; no model edits."""
import cProfile, io, json, pstats, runpy, time
from pathlib import Path
from dataclasses import replace
import numpy as np
import src.solvers.wirsum_1997_named_f_e8 as e8
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
out=Path('docs/audit_reports/plan_ab_diagnosis_2026-09-28')
nq=replace(e8.freeze_e8_q(Q_W=0),stage=11)
cold=runpy.run_path('docs/audit_reports/p0_p4_system_2026-09-26/p3_cold_start.py')['_initial']()
cq=outer_P_named(cold,nq); seeded=e8._seed_downstream_o2(cq,cold); cf=coupled_residual_named(seeded,cq)
zeros=[{'index':int(i),'F':float(cf[i])} for i in np.flatnonzero(seeded==0)]
x=np.load('docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz')['x'].copy(); x[:-2]*=.5
q=outer_P_named(x,nq)
orig=e8.coupled_residual_named
stats={'F_calls':0,'F_seconds':0.,'best_merit':float('inf')}
def counted(z,*a,**kw):
    t=time.perf_counter(); f=orig(z,*a,**kw); stats['F_seconds']+=time.perf_counter()-t;stats['F_calls']+=1
    merit=e8.residual_block_merit(f)
    if merit<stats['best_merit']:
        stats['best_merit']=merit;stats['best_evaluated_blocks']=named_endpoint_blocks(f)
    return f
e8.coupled_residual_named=counted
prof=cProfile.Profile(); t=time.perf_counter(); prof.enable()
try:
    e8._e8_log_scaled_correct(q,x,deadline_s=10); stats['status']='returned'
except TimeoutError: stats['status']='timeout'
finally: prof.disable();e8.coupled_residual_named=orig
stats['elapsed_seconds']=time.perf_counter()-t
s=io.StringIO();pstats.Stats(prof,stream=s).sort_stats('cumulative').print_stats(20)
(out/'profile.txt').write_text(s.getvalue())
(out/'evidence.json').write_text(json.dumps({'cold_zero_coordinates_after_seed':zeros,'half_state_frozen_q_profile':stats,'note':'Profiling overhead included. Best evaluated state can be a finite-difference probe, not accepted iterate. No P calls inside corrector.'},indent=2))
print((out/'evidence.json').read_text()); print(s.getvalue())
