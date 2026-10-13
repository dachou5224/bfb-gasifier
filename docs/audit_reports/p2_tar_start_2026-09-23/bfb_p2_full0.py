from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_probe.py')).read().split("print('E7'")[0])
from scipy.optimize import least_squares
z0=np.load(str(_ART / 'bfb_p2_small.npz'))['x'];nq7=freeze_e7_q();nq8=freeze_e8_q()
def res(z):
 qq=replace(outer_P_named(z,nq7),stage=8,reference=nq8.reference)
 return coupled_residual_named(z,qq)
# Interior positive coordinates, never eliminate oxygen. Existing exact structural zeros stay fixed.
print('zero slots',np.where(z0[:-1]==0)[0],flush=True)
for i in [21,28,35]:
 if z0[i]<=0:z0[i]=1e-24
ids=np.where(z0>0)[0];base=z0[ids].copy()
scale=np.full(45,0.01);scale[[14,15,37,38,39]]=0.001;scale[-1]=1.
def unpack(v):
 z=z0.copy();z[ids]=base*np.exp(v);return z
def fun(v):return res(unpack(v))/scale
r=least_squares(fun,np.zeros(len(ids)),jac='3-point',diff_step=1e-5,max_nfev=80,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds=(-40,40))
z=unpack(r.x);ff=res(z)
gasids=list(range(14))+list(range(16,37));si=[14,15,37,38,39]
print('RESULT',r.nfev,r.message,'gas',np.linalg.norm(ff[gasids]),'solid',np.linalg.norm(ff[si]),'tar',np.linalg.norm(ff[40:44]),'energy',ff[-1],'T',z[-1],'optimality',r.optimality,flush=True)
print('F',ff,flush=True)
np.savez(str(_ART / 'bfb_p2_full0.npz'),x=z,F=ff)
