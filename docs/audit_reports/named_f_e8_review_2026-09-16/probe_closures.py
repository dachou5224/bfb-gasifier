"""Local closure probes; no full solve or E8 continuation chain."""
import json,numpy as np
from src.solvers.wirsum_1997_named_f import default_x0_named,outer_P_named,daem_X_at_T,r10_c_tar,_conc_major
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q,residual_named_e8,jacobian_named_e8
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e6 import freeze_e6_q
from src.solvers.wirsum_1997_named_f_r10 import freeze_e6_r10_q
from src.solvers.wirsum_1997_named_f_r11 import freeze_e6_r10_r11_q
from src.core.constants import Rg
q=freeze_e8_q();x=default_x0_named(q);out=[]
for T in [800.,1420.]:
 z=x.copy();z[-1]=T
 p=outer_P_named(z,q);pp=outer_P_named(z,p)
 J=jacobian_named_e8(z,freeze_e6_q(),freeze_e6_r10_q(),freeze_e6_r10_r11_q(),freeze_e7_q(),q)
 f=residual_named_e8(z,q);cols=[]
 for i in [14,15,37,38,39]:
  zp=z.copy();h=1e-6;zp[i]+=h
  fd=(residual_named_e8(zp,q)-f)/h
  # Solid rows avoid subtraction against enormous gas reaction sources.
  rows=[14,15,37,38,39]
  cols.append({'column':i,'solid_rows':rows,'fd':fd[rows].tolist(),'implemented':J[rows,i].tolist(),'max_gap':float(np.max(abs(fd[rows]-J[rows,i])))})
 Ctar=r10_c_tar(z[:7],z[7:14],z[40]+z[41],T,q.energy.species.geo_ex.P)
 Cmajor=_conc_major(z[:7],T,q.energy.species.geo_ex.P)
 out.append({'T':T,'w_leave':p.energy.species.w_leave.tolist(),'K':p.energy.species.K_eject.tolist(),'rho_once':p.energy.species.geo_fb.rho_g,'rho_twice':pp.energy.species.geo_fb.rho_g,'X_DAEM':daem_X_at_T(T,q.tau_daem_s),'unreleased_C_mol_s':q.energy.feed.nC_vm*(1-daem_X_at_T(T,q.tau_daem_s)),'major_plus_tar_over_P_RT':float((sum(Cmajor)+Ctar)/(q.energy.species.geo_ex.P/(Rg*T))),'jacobian_columns':cols})
print(json.dumps(out,indent=2))
