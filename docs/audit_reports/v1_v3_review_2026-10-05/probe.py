"""Fixed-state trust model check, no path rerun or source edits."""
import json
from pathlib import Path
import numpy as np
import scipy.optimize as opt
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,tar_consumption_fractions
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct,residual_block_scale
from src.solvers.wirsum_1997_named_f_control import material_gate_ratios,step_reduces_material_violation
root=Path(__file__).parent
x=np.array(json.loads(Path('docs/audit_reports/t1_t4_review_2026-09-30/u3_path.json').read_text())['s_0_16_controlled']['x'])
nq=q_on_wesseling1_path(.16);q=outer_P_named(x,nq); f=coupled_residual_named(x,q);ids=np.flatnonzero(x>0);cap={};orig=opt.lsq_linear
# Capture an existing production call without writing or replacing production files.
def capture(A,b,**kw):
 cap.update(A=A.copy(),b=b.copy(),kw=kw);return orig(A,b,**kw)
opt.lsq_linear=capture
try:_e8_log_scaled_correct(q,x,log_span=.5,temperature_factor=1.5)
finally:opt.lsq_linear=orig
A=cap['A'];b=cap['b'];sc=residual_block_scale(50);a1,a2=tar_consumption_fractions()
def inv(g):
 r=g.copy()
 for i in [40,42,44]:r[i]=a1*g[i]+g[i+1];r[i+1]=a2*g[i]-g[i+1]
 return r
rows=[]
for radius in [.5,.05,.005,.0005]:
 lo=np.full(len(ids),-radius);hi=-lo
 # Limit temperature consistently with radius in this diagnostic only.
 lo[ids>=48]=-min(radius,.001);hi[ids>=48]=min(radius,.001)
 fit=orig(A,b,bounds=(lo,hi),method='bvls',tol=1e-10,max_iter=100);v=fit.x
 for alpha in [1,.1,.01]:
  z=x.copy();z[ids]*=np.exp(alpha*v);qq=outer_P_named(z,nq);ff=coupled_residual_named(z,qq)
  pred=inv(-b+alpha*A@v)*sc;gr=material_gate_ratios(ff);pg=material_gate_ratios(pred)
  rows.append(dict(radius=radius,alpha=alpha,ratios=gr,predicted_ratios=pg,aggregate=float(np.linalg.norm(gr)),maximum=max(gr),componentwise_accept=step_reduces_material_violation(f,ff),energy_W=ff[-2:].tolist(),closed=bool(np.max(qq.hydro.w_undecided)<=1e-12),linear_status=int(fit.status)))
out=dict(state='U3 saved state; V2 final vector was not saved',baseline=material_gate_ratios(f),baseline_aggregate=float(np.linalg.norm(material_gate_ratios(f))),rows=rows)
(root/'evidence.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
