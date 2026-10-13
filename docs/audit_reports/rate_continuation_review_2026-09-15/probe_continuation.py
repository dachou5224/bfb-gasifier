"""Read-only model probe: endpoint residual and projected-Jacobian consistency."""
import json
import numpy as np
import src.solvers.wirsum_1997_named_f_r10 as r10
from src.solvers.wirsum_1997_named_f_e6 import residual_named_e6
from src.solvers.wirsum_1997_named_f_r12_r5 import solve_named_r12_then_r5

records=[]
original=r10.wirsum_1997_strict_newton

def wrapper(res,jac,x,**kw):
    stats={'negative_evaluations':0,'worst_negative':0.0}
    def watched(z):
        low=float(np.min(z[:-1]))
        if low<0:
            stats['negative_evaluations']+=1
            stats['worst_negative']=min(stats['worst_negative'],low)
        return res(z)
    sol=original(watched,jac,x,**kw)
    stats.update(stop=sol.failure_reason,raw_final_min=float(np.min(sol.x[:-1])),raw_final_residual=float(np.linalg.norm(sol.F)),accepted_steps=len(sol.accepted_alphas))
    records.append(stats)
    return sol
r10.wirsum_1997_strict_newton=wrapper
out={}
a=solve_named_r12_then_r5(max_iter=20)
out['r5']={k:a[k] for k in ['physical_convergence','residual_l2','T_K','continuation']}
b=r10.solve_named_e6_then_r10(max_iter=20)
x=np.array(b['x']);q0=r10.freeze_e6_q();q1=r10.freeze_e6_r10_q()
f=lambda z:r10.residual_e6_plus_scaled_r10(z,q0,q1,r10_scale=1.0)
out['r10']={k:b[k] for k in ['physical_convergence','residual_l2','T_K','continuation','x_nonneg_projection']}
out['r10']['fresh_direct_full_residual_l2']=float(np.linalg.norm(residual_named_e6(x,q1)))
out['r10']['endpoint_identity_max_gap']=float(np.max(np.abs(f(x)-residual_named_e6(x,q1))))
out['r10']['returned_min']=float(np.min(x[:-1]));out['r10']['inner_stats']=records
# Entire negative half-line is flat after projection; implementation still uses boundary Jacobian.
i=int(np.argmin(x[:-1]));z=x.copy();z[i]=-1e-5
zp=z.copy();zp[i]+=1e-6
fd=(f(zp)-f(z))/1e-6
J=r10.jacobian_e6_plus_scaled_r10(z,q0,q1,r10_scale=1.0)
out['projection_counterexample']={'index':i,'x_i':float(z[i]),'domain_accepts':r10._in_domain(z),'fd_column_norm':float(np.linalg.norm(fd)),'implemented_column_norm':float(np.linalg.norm(J[:,i])),'F_equals_F_clipped':bool(np.array_equal(f(z),f(r10._clip_nonneg(z))))}
print(json.dumps(out,indent=2))
