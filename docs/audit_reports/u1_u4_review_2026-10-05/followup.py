"""Reuse saved Jacobians, bounded linear subproblem and P-group attribution."""
import json,time
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.optimize import lsq_linear
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,tar_consumption_fractions
from src.solvers.wirsum_1997_named_f_e8 import residual_block_scale,residual_block_merit
r=Path(__file__).parent; d=np.load(r/'jacobians.npz');x=d['x'];ids=d['ids'];p=d['p'];nq=q_on_wesseling1_path(.16);q=outer_P_named(x,nq);sc=residual_block_scale(50);a1,a2=tar_consumption_fractions();f=coupled_residual_named(x,q);y=f/sc
out={'repeated_merit':[residual_block_merit(coupled_residual_named(x,outer_P_named(x,nq))) for _ in range(3)],'bounded':[],'feedback_groups':[]}
def inv(J):
    z=J.copy()
    for i in [40,42,44]:z[i]=a1*J[i]+J[i+1];z[i+1]=a2*J[i]-J[i+1]
    return z
for key in ['Ji','Jh']:
    J=d[key];lo=np.full(len(ids),-.1);hi=-lo;lo[ids>=48]=-.001;hi[ids>=48]=.001
    fit=lsq_linear(J,-d['g0'],bounds=(lo,hi),method='bvls',tol=1e-10,max_iter=100)
    v=fit.x;rows=[]
    for alpha in [.001,.01,.1,1]:
        z=x.copy();z[ids]*=np.exp(alpha*v);ff=coupled_residual_named(z,outer_P_named(z,nq))
        rows.append({'alpha':alpha,'merit':residual_block_merit(ff),'energy_W':ff[-2:].tolist(),'solid_l2':float(np.linalg.norm(ff[[14,15,37,38,39,46,47]])),'tar_l2':float(np.linalg.norm(ff[40:46])),'gas_l2':float(np.linalg.norm(ff[np.r_[0:14,16:37]]))})
    out['bounded'].append({'matrix':key,'linear_status':int(fit.status),'linear_residual_norm':float(np.linalg.norm(J@v+d['g0'])),'predicted_raw_slope':float(y@(inv(J)@v)),'steps':rows})
for h in [1e-3,1e-4]:
    dx=np.zeros_like(x);dx[ids]=x[ids]*p
    qp=outer_P_named(x+h*dx,nq);qm=outer_P_named(x-h*dx,nq)
    def hybrid(qq,group):
        sp=q.energy.species;ss=qq.energy.species
        if group=='Kbd':sp=replace(sp,geo_ex=replace(sp.geo_ex,K_bd=ss.geo_ex.K_bd))
        if group=='Dg':sp=replace(sp,D_g=ss.D_g)
        if group=='eject':sp=replace(sp,K_eject=ss.K_eject)
        if group=='trajectory_geometry':sp=replace(sp,geo_fb=ss.geo_fb,w_leave=ss.w_leave,w_return=ss.w_return)
        return replace(q,energy=replace(q.energy,species=sp),hydro=qq.hydro if group=='trajectory_geometry' else q.hydro)
    for group in ['Kbd','Dg','eject','trajectory_geometry']:
        df=(coupled_residual_named(x,hybrid(qp,group))-coupled_residual_named(x,hybrid(qm,group)))/(2*h*sc)
        out['feedback_groups'].append({'h':h,'group':group,'slope':float(y@df),'derivative_norm':float(np.linalg.norm(df))})
(r/'followup.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
