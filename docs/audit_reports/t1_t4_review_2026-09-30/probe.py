"""Fixed-state audit only: no solves of path, no production edits."""
import json
from pathlib import Path
import numpy as np
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named,coupled_residual_named,refresh_pyrolysis_named
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct,residual_block_merit
p=Path('docs/audit_reports/s1_s3_review_2026-09-30/t3_corrected_path.json')
d=json.loads(p.read_text())
r=[r for r in d['toward_022'] if abs(r['s']-.16)<1e-12 and not r['passed_block_gates']][-1]
x=np.array(r['x']); nq=q_on_wesseling1_path(.16); q=outer_P_named(x,nq)
f=lambda z:coupled_residual_named(z,q)
f0=f(x); ids=[42,43]; results=[]
for h in [1e-3,1e-4,1e-5]:
    cols=[]
    for i in ids:
        a=x.copy(); b=x.copy(); a[i]*=np.exp(h); b[i]*=np.exp(-h)
        cols.append((f(a)-f(b))/(2*h))
    J=np.column_stack(cols); v=np.linalg.lstsq(J[ids,:],-f0[ids],rcond=1e-8)[0]; z=x.copy(); z[ids]*=np.exp(v); ff=f(z)
    results.append(dict(h=h,J_tar=J[ids,:].tolist(),singular_values=np.linalg.svd(J[ids,:],compute_uv=False).tolist(),v=v.tolist(),new_tar=z[ids].tolist(),tar_before=f0[ids].tolist(),tar_after=ff[ids].tolist(),merit=residual_block_merit(ff),gas_delta_norm=float(np.linalg.norm((ff-f0)[:14]))))
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry
from src.solvers.wirsum_1997_named_f import R10_FUEL_TYPE
sto=get_lumped_tar_stoichiometry('R10',R10_FUEL_TYPE); sto11=get_lumped_tar_stoichiometry('R11',R10_FUEL_TYPE)
a1=-sto['TAR1']; a2=-sto['TAR2']
contrast=dict(a1=a1,a2=a2,r11_a1=-sto11['TAR1'],r11_a2=-sto11['TAR2'],feed=dict(q.energy.n_tar),x_contrast=[float(a2*x[i]-a1*x[i+1]) for i in [40,42,44]],F_contrast=[float(a2*f0[i]-a1*f0[i+1]) for i in [40,42,44]])
st={}; z=_e8_log_scaled_correct(q,x,max_nfev=8,max_j=4,ftol=1e-4,temperature_factor=1.5,log_span=4,stats=st)
def inner(z): return coupled_residual_named(z,refresh_pyrolysis_named(z,q))
def outer(z): return coupled_residual_named(z,outer_P_named(z,nq))
steps=[]
for a in [0,1e-3,.01,.1,1]:
    y=x+a*(z-x)
    steps.append(dict(alpha=a,stale=residual_block_merit(f(y)),inner=residual_block_merit(inner(y)),composite=residual_block_merit(outer(y))))
minor=q.energy.species.n_NH3_bypass_mol_s+q.energy.species.n_H2S_bypass_mol_s
out=dict(contrast=contrast,tar_probe=results,corrector_stats=st,step_comparison=steps,T_change=(z[-2:]-x[-2:]).tolist(),minor_mol_s=minor,minor_fraction_fb=float(minor/(minor+sum(x[16:23]))))
Path(__file__).with_name('evidence.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
