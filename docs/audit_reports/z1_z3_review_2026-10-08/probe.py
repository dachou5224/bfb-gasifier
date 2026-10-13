"""Fixed-point predictor audit; no continuation or correction solve."""
import json, time
from pathlib import Path
import numpy as np
import scipy.optimize as opt
from src.solvers.wirsum_1997_named_f import outer_P_named, coupled_residual_named
from src.solvers.wirsum_1997_named_f_control import progress_phi, progress_ratios, certify_controlled_state
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.external_state_ledger import external_imbalance
out=Path(__file__).parent
saved=json.loads(Path('docs/audit_reports/x1_x3_review_2026-10-07/y2_continue.json').read_text())['points'][-1]
s=float(saved['s']); x=np.array(saved['x']); ds=.002
full=lambda t,z:coupled_residual_named(z,outer_P_named(z,q_on_wesseling1_path(t)))
f=full(s,x); h=1e-4; offset=(full(s+h,x)-full(s-h,x))/(2*h)*ds
capture={}; orig=opt.minimize

def wrapper(fun,*args,**kw):
    fit=orig(fun,*args,**kw)
    if fun.__name__=='cost_y' and kw.get('method')=='SLSQP':
        y=fit.x; g=kw['jac'](y); cols=[]; margins=[]; labels=[]
        for i,c in enumerate(kw['constraints']):
            val=float(c['fun'](y))
            if val<1e-7:
                cols.append(c['jac'](y)); margins.append(val); labels.append('constraint_'+str(i))
        for i,(lo,hi) in enumerate(kw['bounds']):
            e=np.zeros(y.size);e[i]=1
            if y[i]-lo<1e-7: cols.append(e);margins.append(y[i]-lo);labels.append('lower_'+str(i))
            if hi-y[i]<1e-7: cols.append(-e);margins.append(hi-y[i]);labels.append('upper_'+str(i))
        A=np.array(cols).T
        lam,_=opt.nnls(A,g)
        r=g-A@lam
        capture.update(active=labels,multipliers=lam.tolist(),stationarity_norm=float(np.linalg.norm(r)),gradient_norm=float(np.linalg.norm(g)),relative_stationarity=float(np.linalg.norm(r)/max(np.linalg.norm(g),1e-30)),complementarity_max=float(np.max(np.abs(lam*np.array(margins)))),normalized_objective=float(fun(y)))
    return fit
opt.minimize=wrapper
stats={};t=time.monotonic()
try: z=_e8_log_scaled_correct(outer_P_named(x,q_on_wesseling1_path(s)),x,log_span=.5,ratio_span=1,stats=stats,residual_offset=offset)
finally: opt.minimize=orig
ft=full(s+ds,z); before=full(s+ds,x);pred=np.array(stats['predicted_raw'])
ledger=external_imbalance(x,outer_P_named(x,q_on_wesseling1_path(s)))
r={'source_sha256':named_kernel_source_sha256(),'wall_s':time.monotonic()-t,'s040_cert':certify_controlled_state(x,f,outer_P_named(x,q_on_wesseling1_path(s)).hydro),'s040_ledger':ledger,'kkt':capture,'target_entrance_phi':progress_phi(before),'predicted_phi':progress_phi(pred),'actual_phi':progress_phi(ft),'actual_ratios':progress_ratios(ft),'predicted_ratios':progress_ratios(pred),'target_rho':(progress_phi(before)-progress_phi(ft))/(progress_phi(before)-progress_phi(pred)),'largest_actual_rows':sorted(enumerate(ft.tolist()),key=lambda p:abs(p[1]),reverse=True)[:8],'x_pred':z.tolist(),'F_target':ft.tolist(),'new_root':False}
def plain(v):
    if isinstance(v,dict):return {str(k):plain(w) for k,w in v.items()}
    if isinstance(v,(list,tuple)):return [plain(w) for w in v]
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    return v
(out/'evidence.json').write_text(json.dumps(plain(r),indent=2))
print(json.dumps(plain({k:v for k,v in r.items() if k not in ('x_pred','F_target','s040_ledger')}),indent=2))
