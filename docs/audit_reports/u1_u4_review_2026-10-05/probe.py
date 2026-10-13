"""Single fixed-state derivative audit; no solver/source modifications."""
import json,time,hashlib
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,refresh_pyrolysis_named,apply_tar_reaction_transport_rows,tar_consumption_fractions
from src.solvers.wirsum_1997_named_f_e8 import residual_block_scale,residual_block_merit,_e8_positive_log_ids,_e8_log_scaled_correct
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
root=Path(__file__).parent
d=json.loads(Path('docs/audit_reports/t1_t4_review_2026-09-30/u3_path.json').read_text())
x=np.array(d['s_0_16_controlled']['x']); nq=q_on_wesseling1_path(.16); q=outer_P_named(x,nq); scale=residual_block_scale(len(x)); ids=_e8_positive_log_ids(x); a1,a2=tar_consumption_fractions()
def inv(y):
    z=y.copy()
    for i in [40,42,44]:
        z[i]=a1*y[i]+y[i+1];z[i+1]=a2*y[i]-y[i+1]
    return z
count={'F':0,'P':0}; und=[]
def eval(z,full):
    qq=outer_P_named(z,nq) if full else refresh_pyrolysis_named(z,q)
    count['P']+=int(full);count['F']+=1
    if full:und.append(float(np.max(qq.hydro.w_undecided)))
    f=coupled_residual_named(z,qq)
    return apply_tar_reaction_transport_rows(z,qq,f)/scale
f0=coupled_residual_named(x,q); g0=eval(x,True); mats=[[],[]];t0=time.monotonic()
for k,i in enumerate(ids):
    if time.monotonic()-t0>90: raise RuntimeError('fixed probe budget exceeded')
    h=1e-7 if i>=48 else 1e-4
    xp=x.copy();xm=x.copy();xp[i]*=np.exp(h);xm[i]*=np.exp(-h)
    for j,full in enumerate([False,True]):mats[j].append((eval(xp,full)-eval(xm,full))/(2*h))
    if k%10==0:print('columns',k,'elapsed',time.monotonic()-t0,flush=True)
Ji,Jh=[np.column_stack(m) for m in mats]
# Equivalent raw residual derivatives from stable transformed rows.
Ri,Rh=[inv(J) for J in [Ji,Jh]]
st={}; prop=_e8_log_scaled_correct(q,x,max_j=4,max_nfev=8,ftol=1e-4,log_span=4,temperature_factor=1.5,stats=st)
p=(prop[ids]-x[ids])/x[ids]
y0=f0/scale
slope={'inner':float(y0@(Ri@p)),'full':float(y0@(Rh@p)),'feedback':float(y0@((Rh-Ri)@p))}
colterms=y0[:,None]*(Rh-Ri)*p[None,:]
# A single equilibrated linear solve per matrix, bounded in log displacement by rescaling.
trials=[]
for label,J in [('inner',Ji),('full',Jh)]:
    cn=np.maximum(np.linalg.norm(J,axis=0),1e-12);A=J/cn
    v,res,rank,sv=np.linalg.lstsq(A,-g0,rcond=1e-10);v=v/cn
    v*=min(1.,.1/max(np.max(np.abs(v)),1e-300))
    rj=inv(J);pred=float(y0@(rj@v))
    rr=[]
    for alpha in [.001,.01,.1,1.]:
        z=x.copy();z[ids]*=np.exp(alpha*v);gg=eval(z,True);f=inv(gg)*scale
        rr.append({'alpha':alpha,'merit':residual_block_merit(f),'gas_l2':float(np.linalg.norm(f[np.r_[0:14,16:37]])),'tar_l2':float(np.linalg.norm(f[40:46]))})
    trials.append({'matrix':label,'rank':int(rank),'singular_values':sv.tolist(),'max_log_step':float(np.max(np.abs(v))),'predicted_raw_merit_squared_half_slope':pred,'trials':rr})
# Check old corrector-direction derivatives directly at two small displacements.
direct=[]
for eps in [1e-3,1e-4]:
    z=x+eps*(prop-x)
    g=eval(z,True);direct.append({'epsilon':eps,'raw_merit':float(np.linalg.norm(inv(g))),'slope_half_merit_squared':float((np.linalg.norm(inv(g))**2-np.linalg.norm(y0)**2)/(2*eps))})
np.savez(root/'jacobians.npz',x=x,ids=ids,Ji=Ji,Jh=Jh,p=p,g0=g0)
out={'state':'U3 new s=0.16 failed point','baseline_merit':residual_block_merit(f0),'active_ids':ids.tolist(),'zero_ids':np.flatnonzero(x==0).tolist(),'jacobian_elapsed_plus_trials_s':time.monotonic()-t0,'calls':count,'max_undecided':max(und),'inner_corrector_stats':st,'corrector_slopes':slope,'feedback_by_column':sorted([{'index':int(i),'slope':float(colterms[:,k].sum())} for k,i in enumerate(ids)],key=lambda r:abs(r['slope']),reverse=True),'feedback_by_row':sorted([{'index':i,'slope':float(t)} for i,t in enumerate(colterms.sum(axis=1))],key=lambda r:abs(r['slope']),reverse=True),'direct_check':direct,'linear_trials':trials}
(root/'evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
