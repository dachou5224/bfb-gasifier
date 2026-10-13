from pathlib import Path
_ART = Path(__file__).resolve().parent
import json, numpy as np
from scipy.optimize import least_squares
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q,blend_hydro_named_q,_balance_freeboard_o2,_e8_log_scaled_correct,residual_e7_plus_scaled_e8
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
p='docs/audit_reports/p2_tar_start_2026-09-23/bfb_p2_hydro_march.npz'
d=np.load(p);x=d['x'];s=.962
q7,q8=freeze_e7_q(),freeze_e8_q()
q=blend_hydro_named_q(outer_P_named(x,q7),outer_P_named(x,q8),hydro_scale=s,fate_scale=0.)
z=_e8_log_scaled_correct(q,x)
scale=np.full(45,.01);scale[[14,15,37,38,39]]=.001;scale[-1]=1.
def f(z):return coupled_residual_named(z,q)/scale
print('baseline',named_endpoint_blocks(f(z)*scale),flush=True)
# Existing positive oxygen from scalar balance; do not seed from residual or bed oxygen.
assert z[21]>0
ids=np.flatnonzero(z>0);base=z[ids].copy()
def unpack(v):
 zz=z.copy();zz[ids]=base*np.exp(v);return zz
def fun(v):return f(unpack(v))
def jac(v):
 out=np.zeros((45,len(ids)))
 for j,i in enumerate(ids):
  h=1e-7 if i==44 else 1e-4
  vp=v.copy();vm=v.copy();vp[j]+=h;vm[j]-=h
  out[:,j]=(fun(vp)-fun(vm))/(2*h)
 return out
J=jac(np.zeros(len(ids)));r=fun(np.zeros(len(ids)));g=J.T@r;j=int(np.where(ids==21)[0][0]);col=J[:,j]
print('fb_o2',z[21],'logcol_norm',np.linalg.norm(col),'gradient',g[j], 'column_16_17_21',col[[16,17,21]],flush=True)
fit=least_squares(fun,np.zeros(len(ids)),jac=jac,max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds=(-40.,40.))
y=unpack(fit.x);ff=f(y)*scale;h=residual_e7_plus_scaled_e8(y,q7,q8,hydro_scale=s,fate_scale=0.)
print('joint',fit.nfev,fit.message,'optimality',fit.optimality,'frozen',named_endpoint_blocks(ff),'composite',named_endpoint_blocks(h),flush=True)
np.savez(str(_ART / 'bfb_p2_0924_joint.npz'),x=y,F=ff,H=h,hydro_scale=s,fate_scale=0.)
json.dump({'baseline':named_endpoint_blocks(f(z)*scale),'fb_o2':float(z[21]),'column_norm':float(np.linalg.norm(col)),'omitted_gradient':float(g[j]),'joint':named_endpoint_blocks(ff),'composite':named_endpoint_blocks(h),'nfev':fit.nfev,'message':fit.message,'optimality':float(fit.optimality)},open(str(_ART / 'bfb_p2_0924_probe.json'),'w'),indent=2)
