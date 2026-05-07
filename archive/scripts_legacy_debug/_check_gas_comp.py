"""检查出口气体组分并与 Lu et al. 文献对比。"""
import sys, numpy as np
sys.path.insert(0, ".")
from src.core.reactor import Reactor, ReactorConfig

cfg = ReactorConfig(
    n_cells=10, H_bed=6.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=293.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.337, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.08,
)
reactor = Reactor(cfg)
res = reactor.solve(max_global_iter=5, tol_global=1e-3)

print("=== 出口气体（湿基摩尔分数）===")
for sp, v in sorted(res["exit_gas"].items(), key=lambda x: -x[1]):
    if v > 1e-6:
        print(f"  {sp:6s}: {v*100:7.3f}%")

print()
print("=== 出口气体（干基摩尔分数）===")
for sp, v in sorted(res["exit_gas_dry"].items(), key=lambda x: -x[1]):
    if v > 1e-6:
        print(f"  {sp:6s}: {v*100:7.3f}%")

# Hamel (1999) HTW Wesseling Sim 1 — Fig. 7.4 干基读数
ref = {"CO": 15.7, "CO2": 13.3, "H2": 14.5, "CH4": 3.4, "N2": 52.0, "O2": 0.0}

print()
print("=== 与 Hamel HTW Wesseling Sim1 文献对比（干基）===")
print(f"  {'组分':6s}  {'模型':>8s}  {'文献':>8s}  {'绝对偏差':>10s}  {'状态':>6s}")
print("  " + "-"*52)
dry = res["exit_gas_dry"]
total_err = 0.0
for sp in ["CO", "CO2", "H2", "CH4", "N2", "O2"]:
    model_val = dry.get(sp, 0.0) * 100
    ref_val = ref.get(sp, 0.0)
    diff = model_val - ref_val
    total_err += diff**2
    flag = "✓" if abs(diff) < 3.0 else ("△" if abs(diff) < 8.0 else "✗")
    print(f"  {sp:6s}  {model_val:8.2f}%  {ref_val:8.2f}%  {diff:+10.2f}%  {flag:>6s}")

rms_comp = (total_err / len(ref))**0.5
print()
print(f"  组分 RMS 偏差 = {rms_comp:.2f} pp")
print()
print(f"  T_exit    = {res['T_profile'][-1]:.1f} K   (目标: ~1100 K,  偏差 {res['T_profile'][-1]-1100:+.1f} K)")
print(f"  Xc        = {res['carbon_conv']:.3f}       (目标: ~0.95,    偏差 {res['carbon_conv']-0.95:+.3f})")
print(f"  converged = {res['converged']}   n_iter = {res['n_iter']}")
