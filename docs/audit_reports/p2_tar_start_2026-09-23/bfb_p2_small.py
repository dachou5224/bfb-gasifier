from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_probe.py')).read().split("print('E7'")[0])
from scipy.optimize import least_squares
ids=np.array(list(range(14))+[40,41,42,43])
# Bed plus four tar variables; keep zero inflow species fixed only if corresponding row is zero.
ids=np.array([i for i in ids if x[i]>0 or abs(F[i])>1e-12])
z0=x.copy();z0[42:44]=[1.2e-5,1.8e-5]
print('ids',ids,'x',z0[ids],flush=True)
base=z0[ids].copy()
def unpack(v):
 z=z0.copy();z[ids]=base*np.exp(v);return z
def fun(v):return f(unpack(v))[ids]
r=least_squares(fun,np.zeros(len(ids)),jac='3-point',diff_step=1e-4,max_nfev=80,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds=(-40,40))
z=unpack(r.x); ff=f(z)
print('RESULT',r.nfev,r.message,'bedtar',np.linalg.norm(ff[ids]),'energy',ff[-1],'allgas',np.linalg.norm(ff[:37]),'min',min(z),'optimality',r.optimality)
print('bed',ff[:14], 'tar',ff[40:44], 'x_tar',z[40:44])
np.savez(str(_ART / 'bfb_p2_small.npz'),x=z,F=ff)
