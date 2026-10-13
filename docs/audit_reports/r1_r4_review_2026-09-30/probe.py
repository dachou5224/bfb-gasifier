"""One saved-state corrector + directional probes, no endpoint writes."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q,_e8_log_scaled_correct,residual_block_merit as merit
from src.physics.wirsum_1997_freeboard import accept_refined_event
nq=replace(freeze_e8_q(),stage=11)
x=np.array(json.load(open('docs/audit_reports/review_abcd_2026-09-28/r2_tplus50.json'))['x'])
q=outer_P_named(x,nq); stats={}
p=_e8_log_scaled_correct(q,x,deadline_s=10,max_nfev=8,max_j=4,ftol=1e-8,temperature_factor=1.5,log_span=4.,stats=stats)
rows=[]
for a in [0,1e-4,.001,.01,.1,.25,.5,1]:
    y=x+a*(p-x);qy=outer_P_named(y,nq)
    rows.append({'alpha':a,'frozen':merit(coupled_residual_named(y,q)),'composite':merit(coupled_residual_named(y,qy))})
# Hybrid q distinguishes pyrolysis refresh from hydrodynamic refresh at same candidate.
qy=outer_P_named(p,nq)
pyro_only=replace(q,energy=replace(qy.energy,species=replace(qy.energy.species,K_eject=q.energy.species.K_eject,w_leave=q.energy.species.w_leave,w_return=q.energy.species.w_return,geo_fb=q.energy.species.geo_fb,geo_ex=q.energy.species.geo_ex,k=q.energy.species.k,D_g=q.energy.species.D_g,rho_g_fb_ref=q.energy.species.rho_g_fb_ref)),hydro=q.hydro)
out={'stats':stats,'direction':rows,'candidate_pyro_refresh_only_merit':merit(coupled_residual_named(p,pyro_only)),'both_undecided_same_time_accepted':accept_refined_event('undecided',15.,'undecided',15.),'note':'Same saved accepted point; single new corrector proposal, not the historical rejected proposal. Hybrid q is diagnostic, not a model endpoint.'}
Path('docs/audit_reports/r1_r4_review_2026-09-30/evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
