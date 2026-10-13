"""Compare analytic R12 boundary derivative with a one-sided difference."""
import json
import numpy as np
from src.solvers.wirsum_1997_named_f_r12 import _r12_r_and_dC
from src.solvers.wirsum_1997_stiffness_full_gas import IDX
C=np.zeros(7); C[IDX['H2']]=1.0
r,J=_r12_r_and_dC(1100.0,C)
rows=[]
for h in [1e-8,1e-10,1e-12]:
    Cp=C.copy(); Cp[IDX['O2']]=h
    rp,_=_r12_r_and_dC(1100.0,Cp)
    rows.append({'h_mol_m3':h,'one_sided_derivative_per_s':float((rp-r)/h)})
print(json.dumps({'T_K':1100.0,'C_H2_mol_m3':1.0,'C_O2_mol_m3':0.0,'analytic_derivative_per_s':float(J[IDX['O2']]),'differences':rows},indent=2))
