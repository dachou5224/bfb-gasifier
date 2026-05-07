import sys, time
sys.path.insert(0, '.')
from src.core.reactor import Reactor, ReactorConfig

cfg = ReactorConfig(ER=0.25)
r = Reactor(cfg)
t0 = time.time()
res = r.solve(max_global_iter=1, tol_global=1e-3)
print('ok', res['n_iter'], f"T={res['T_profile'][-1]:.1f}", f"rms={res['rms_scaled_gs']:.3e}", f"dt={time.time()-t0:.1f}s")
