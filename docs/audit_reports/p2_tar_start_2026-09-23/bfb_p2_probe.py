from pathlib import Path
_ART = Path(__file__).resolve().parent
import json, numpy as np
from src.solvers.wirsum_1997_named_f import *
from src.solvers.wirsum_1997_named_f_e8 import as_e8_state, freeze_e8_q, _fd_jacobian_trace, _e8_active_step
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_endpoint_io import recertify_named_e7
D=json.load(open('data/named_e8_endpoint_2026-09-17.json'))
x7=np.array(D['e7']['x']); x=as_e8_state(x7)
q7=outer_P_named(x7,freeze_e7_q())
# At s=0, stage 8 retains its source model with stage-7 transport.
q=replace(q7,stage=8,reference=freeze_e8_q().reference)
f=lambda z:coupled_residual_named(z,q)
F=f(x)
print('E7',recertify_named_e7(x7))
print('F0',np.linalg.norm(F),F[-1], 'x_tar', x[40:44], 'F_b',F[7:14], 'F_tar',F[40:44])
for name,fun in [('R10',r10_major_and_tar_source),('R11',r11_major_and_tar_source)]:
 a=fun(x7,q7,rate_scale=1); b=fun(x,q,rate_scale=1)
 print('source_delta',name,'dense',b[0]-a[0],'bubble',b[1]-a[1])
J=_fd_jacobian_trace(f,x)
p,_,_=_e8_active_step(J,F,x)
print('p clipped',np.linalg.norm(J@p+F),'dot',F@(J@p),'T',p[-1])
print('tar secant columns')
for h in [1e-20,1e-16,1e-12,1e-8,1e-5]:
 y=x.copy();y[42]+=h; col=(f(y)-F)/h
 print(h, col[42],col[7],col[10],col[-1])
print('row diff full',[(i,float(v)) for i,v in enumerate(F) if abs(v)>1e-6])
np.savez(str(_ART / 'bfb_p2_probe.npz'),x=x,F=F,J=J,p=p)
