from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_probe.py')).read().split("print('E7'")[0])
from scipy.optimize import least_squares
from src.solvers.wirsum_1997_named_f_e8 import blend_hydro_named_q,residual_e7_plus_scaled_e8
import time
nq7=freeze_e7_q();nq8=freeze_e8_q(); z=np.load(str(_ART / 'bfb_p2_verified0.npz'))['x']
scale=np.full(45,0.01);scale[[14,15,37,38,39]]=0.001;scale[-1]=1.
def metrics(ff):return dict(gas=float(np.linalg.norm(ff[list(range(14))+list(range(16,37))])),solid=float(np.linalg.norm(ff[[14,15,37,38,39]])),tar=float(np.linalg.norm(ff[40:44])),energy=float(ff[-1]))
def passed(m):return m['gas']<1e-6 and m['solid']<1e-6 and m['tar']<1e-6 and abs(m['energy'])<1
history=[]
for s in [0.1,0.5,1.]:
 ok=False
 for k in range(8):
  q=blend_hydro_named_q(outer_P_named(z,nq7),outer_P_named(z,nq8),hydro_scale=s,fate_scale=0.)
  z0=z.copy();ids=np.where(z0>0)[0];base=z0[ids].copy()
  def unpack(v):
   zz=z0.copy();zz[ids]=base*np.exp(v);return zz
  def fun(v):return coupled_residual_named(unpack(v),q)/scale
  def jac(v):
   out=np.zeros((45,len(ids)))
   for j,i in enumerate(ids):
    h=1e-7 if i==44 else 1e-4
    vp=v.copy();vm=v.copy();vp[j]+=h;vm[j]-=h
    out[:,j]=(fun(vp)-fun(vm))/(2*h)
   return out
  r=least_squares(fun,np.zeros(len(ids)),jac=jac,max_nfev=80,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds=(-40,40))
  trial=unpack(r.x);ff=residual_e7_plus_scaled_e8(trial,nq7,nq8,hydro_scale=s,fate_scale=0.);m=metrics(ff)
  entry=dict(s=s,picard=k,nfev=r.nfev,T=float(trial[-1]),**m);history.append(entry);print(entry,flush=True)
  z=trial
  if passed(m):ok=True;break
 if not ok:break
 np.savez(str(_ART / 'bfb_p2_hydro_accepted.npz'),x=z,F=ff,hydro_scale=s,fate_scale=0.)
json.dump(history,open(str(_ART / 'bfb_p2_hydro_history.json'),'w'),indent=2)
print('END',s,ok,flush=True)
