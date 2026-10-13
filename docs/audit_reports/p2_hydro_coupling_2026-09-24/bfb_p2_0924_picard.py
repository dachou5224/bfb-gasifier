from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_0924_probe.py')).read().split("z=_e8_log_scaled_correct")[0])
z=np.load(str(_ART / 'bfb_p2_0924_joint.npz'))['x']
scale=np.full(45,.01);scale[[14,15,37,38,39]]=.001;scale[-1]=1.
history=[]
for k in range(6):
 q=blend_hydro_named_q(outer_P_named(z,q7),outer_P_named(z,q8),hydro_scale=s,fate_scale=0.)
 ids=np.flatnonzero(z>0);base=z[ids].copy();z0=z.copy()
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
 fit=least_squares(fun,np.zeros(len(ids)),jac=jac,max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,bounds=(-40.,40.))
 z=unpack(fit.x);h=residual_e7_plus_scaled_e8(z,q7,q8,hydro_scale=s,fate_scale=0.)
 m=named_endpoint_blocks(h);history.append(dict(iteration=k,nfev=fit.nfev,**m));print(history[-1],flush=True)
 if m['passed']:break
np.savez(str(_ART / 'bfb_p2_0924_picard.npz'),x=z,H=h,hydro_scale=s,fate_scale=0.)
json.dump(history,open(str(_ART / 'bfb_p2_0924_picard.json'),'w'),indent=2)
