"""Saved-point audit only; no continuation or production modifications."""
import json
from pathlib import Path
import numpy as np
import scipy.optimize as opt
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,apply_tar_reaction_transport_rows,tar_consumption_fractions
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance
r=Path(__file__).parent;d=json.loads(Path('docs/audit_reports/x1_x3_review_2026-10-07/y2_continue.json').read_text());p=d['points'][-1];s=p['s'];x=np.array(p['x'])
def raw(s,x):return coupled_residual_named(x,outer_P_named(x,q_on_wesseling1_path(s)))
q=outer_P_named(x,q_on_wesseling1_path(s));f=raw(s,x);cert=certify_controlled_state(x,f,q.hydro);led=external_imbalance(x,q)
h=1e-4;der=(raw(s+h,x)-raw(s-h,x))/(2*h);a1,a2=tar_consumption_fractions()
def T(v):
 out=v.copy()
 for i in [40,42,44]:out[i]=v[i]+v[i+1];out[i+1]=a2*v[i]-a1*v[i+1]
 return out
zero=apply_tar_reaction_transport_rows(x,q,np.zeros(50));off=der*.002
bad=apply_tar_reaction_transport_rows(x,q,off);good=T(off)
orig=opt.minimize;records=[]
def capture(fun,v,**kw):
 fit=orig(fun,v,**kw)
 records.append({'success':bool(fit.success),'status':int(fit.status),'message':str(fit.message),'nit':int(fit.nit),'initial_cost':float(fun(v)),'final_cost':float(fun(fit.x)),'constraint_values':[float(c['fun'](fit.x)) for c in kw['constraints']]})
 return fit
opt.minimize=capture
scans=[];stats={}
try:
 for ds in [.02,.002]:
  z=_e8_log_scaled_correct(q,x,log_span=.5,ratio_span=1.,stats=stats,residual_offset=der*ds)
  scans.append({'ds':ds,'linear_status':stats['linear_status'],'max_state_change':float(np.max(np.abs(z-x))),'predicted_energy':np.array(stats['predicted_raw'])[-2:].tolist(),'zero_step_energy':(f+der*ds)[-2:].tolist(),'optimizer':records[-1]})
finally:opt.minimize=orig
# Re-solve the identical ds=.002 local subproblem with scaled variables/objective
# and exact derivatives of its existing affine/quadratic constraints.
norm_records=[]
def normalized(fun,start,**kw):
 if fun.__name__ != 'cost':return orig(fun,start,**kw)
 D=np.maximum(np.max(np.abs(np.array(kw['bounds'])),axis=1),1e-12)
 norm=max(float(fun(start)),1.)
 jac0=kw['jac'];cons=[]
 for c in kw['constraints']:
  defaults=c['fun'].__defaults__
  if len(defaults)==3:
   row,level,limit=defaults;cs=max(float(limit),1.)
   sign=-1. if c['fun'](start+np.ones_like(start)*1e-7)<c['fun'](start) else 1.
   # Determine affine gradient without cancellation from closure's expression.
   names=c['fun'].__code__.co_consts
   # Two paired constraints have opposite gradients; preserve their order.
   sign=-1. if len(cons)%2==0 else 1.
   g=lambda y,row=row,sign=sign,cs=cs:sign*row*D/cs
  else:
   bb,bj=defaults;cs=1e-12
   g=lambda y,bb=bb,bj=bj,cs=cs:-2.*(bj.T@(bb+bj@(D*y)))*D/cs
  cf=c['fun'];cons.append({'type':'ineq','fun':lambda y,cf=cf,cs=cs:cf(D*y)/cs,'jac':g})
 fit=orig(lambda y:fun(D*y)/norm,start/D,jac=lambda y:jac0(D*y)*D/norm,bounds=[(l/d,h/d) for (l,h),d in zip(kw['bounds'],D)],constraints=cons,method='SLSQP',options={'maxiter':200,'ftol':1e-12})
 fit.x=fit.x*D;norm_records.append({'success':bool(fit.success),'status':int(fit.status),'message':str(fit.message),'nit':int(fit.nit)})
 return fit
opt.minimize=normalized
try:
 zz=_e8_log_scaled_correct(q,x,log_span=.5,ratio_span=1.,stats=stats,residual_offset=der*.002)
 normscan={'solver':norm_records[-1],'predicted_energy':np.array(stats['predicted_raw'])[-2:].tolist(),'linear_status':stats['linear_status'],'max_state_change':float(np.max(np.abs(zz-x)))}
finally:opt.minimize=orig
out={'normalized_same_subproblem':normscan,'s':s,'cert':cert,'cho':{k:led['gap_mol_s'][k] for k in ['C','H','O']},'enthalpy_gap_W':led['enthalpy_W']['gap'],'zero_transform_ratio_rows':zero[[41,43,45]].tolist(),'offset_bad_ratio_rows':bad[[41,43,45]].tolist(),'offset_correct_ratio_rows':good[[41,43,45]].tolist(),'offset_transform_error_norm':float(np.linalg.norm(bad-good)),'scans':scans,'material_feasible_tolerance_example':{'block_norm':1e-5,'constraint':1e-12-(1e-5)**2,'accepted_by_internal_feasibility':1e-12-(1e-5)**2>=-1e-8}}
(r/'evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
