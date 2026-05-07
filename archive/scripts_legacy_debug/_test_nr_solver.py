"""Test the global NR solver with rms tracking."""
import sys
import numpy as np
import time

sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig

cfg = ReactorConfig(
    n_cells=10, H_bed=5.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=300.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.25, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.1,
)
t0 = time.time()
reactor = Reactor(cfg)
res = reactor.solve(max_global_iter=20, tol_global=1e-3, solver="global_nr")
dt = time.time() - t0

print(f"T_exit={res['T_profile'][-1]:.1f}K conv={res['converged']} Xc={res['carbon_conv']:.3f} time={dt:.1f}s")
print(f"  converged_outer={res.get('converged_outer')}")
print(f"  converged_inner_nr={res.get('converged_inner_nr')}")
print(f"  rms_scaled_final={res.get('rms_scaled_final'):.4e}")
print(f"  n_iter={res['n_iter']}")
if 'norm_history' in res:
    print(f"  history:")
    for i, rms in enumerate(res['norm_history'], 1):
        print(f"    iter={i} rms={rms:.4e}")
