"""完整输入一致性审查：将模型配置与 Hamel(1999) HTW Wesseling Sim1 文献参数逐项对比。"""
import sys, json, numpy as np
sys.path.insert(0, ".")

from src.core.reactor import ReactorConfig
from src.core.feed_inlet import compute_gas_feeds_mol_s, stoichiometric_o2_mol_s

# ── 文献参数（CASE_HTW_WESSELING_1, Hamel 1999 Table 7.1） ──
REF = {
    # 燃料
    "fuel_type":        "Rhenish Brown Coal",
    "fuel_feed_kg_h":   3377.4,
    "fuel_feed_kg_s":   3377.4 / 3600,
    "moisture_wt_pct":  16.9,
    "ash_dry_wt_pct":   11.41,
    "VM_daf_pct":       53.42,
    "C_dry":            61.5,
    "H_dry":            4.1,
    "O_dry":            21.8,
    "N_dry":            0.68,
    "S_dry":            0.51,
    # 运行条件
    "P_MPa":            2.5,
    "ER":               0.337,
    "T_inlet_K":        293.0,
    "H_bed_m":          6.0,
    "H_total_m":        14.5,
    "D_bed_m":          0.6,
    # 进料流量（文献计算值）
    "O2_actual_kg_h":   1660.0,
    "O2_Nm3_h":         1162.7,
    "air_actual_Nm3_h": 5536.9,
    # 粒径
    "d_p_mm_range":     [1.5, 3.0],
    "n_size_classes":   10,
    # 热损失
    "heat_loss_pct":    "5-10%",
    # 蒸汽/氧摩尔比（文献值需反推）
    # air_actual_Nm3_h / O2_Nm3_h = (O2+N2)/O2 = 4.76 → N2/O2=3.76 √
    # steam: Hamel Table 7.1 gives H2O/O2 ~ 0.74 (from steam_kg_h/O2_kg_h * MW_O2/MW_H2O)
}

# ── 当前模型配置 ──
cfg = ReactorConfig(
    n_cells=10, H_bed=6.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=293.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, S_dry=0.51, HHV_dry=22.0,
    ER=0.337, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.1, heat_loss_frac=0.075,
    nitrogen_fraction=0.68,
)

SEP = "=" * 68

def check(label, model_val, ref_val, unit="", tol_pct=2.0, fmt=".4g"):
    if ref_val is None or ref_val == "__MISSING__":
        tag = "？"
        pct = float("nan")
    else:
        pct = 100.0 * (float(model_val) - float(ref_val)) / (abs(float(ref_val)) + 1e-30)
        tag = "✓" if abs(pct) <= tol_pct else ("△" if abs(pct) <= 10.0 else "✗")
    m_str = format(float(model_val), fmt) if isinstance(model_val, (int, float)) else str(model_val)
    r_str = format(float(ref_val), fmt) if isinstance(ref_val, (int, float)) else str(ref_val)
    pct_str = f"{pct:+.1f}%" if not np.isnan(pct) else "  —"
    print(f"  {tag}  {label:<35s} 模型={m_str:>12s}{unit}  文献={r_str:>12s}{unit}  偏差={pct_str}")


print(SEP)
print("  [1] 燃料进料量")
print(SEP)
ref_kg_s = REF["fuel_feed_kg_h"] / 3600
check("fuel_feed [kg/s]", cfg.fuel_feed, ref_kg_s, " kg/s", tol_pct=0.5, fmt=".5f")

print()
print(SEP)
print("  [2] 燃料工业分析 / 元素分析")
print(SEP)
check("moisture_wt [%]",  cfg.moisture_wt,   REF["moisture_wt_pct"],  " %", tol_pct=0.1)
check("ash_dry_wt  [%]",  cfg.ash_dry_wt,    REF["ash_dry_wt_pct"],   " %", tol_pct=0.1)
check("VM_daf      [%]",  cfg.VM_daf,        REF["VM_daf_pct"],        " %", tol_pct=0.1)
check("C_dry       [%]",  cfg.C_dry,         REF["C_dry"],             " %", tol_pct=0.1)
check("H_dry       [%]",  cfg.H_dry,         REF["H_dry"],             " %", tol_pct=0.1)
check("O_dry       [%]",  cfg.O_dry,         REF["O_dry"],             " %", tol_pct=0.1)
check("N_dry       [%]",  cfg.nitrogen_fraction, REF["N_dry"],         " %", tol_pct=0.5)
check("S_dry       [%]",  cfg.S_dry,         REF["S_dry"],             " %", tol_pct=0.5)

# 元素闭合（干基 wt% 之和）
elem_sum = cfg.C_dry + cfg.H_dry + cfg.O_dry + cfg.nitrogen_fraction + cfg.S_dry + cfg.ash_dry_wt
print(f"\n  元素闭合（干基，含灰）: {elem_sum:.2f}%  (应≈100%,  差 {elem_sum-100:.2f}%)")

print()
print(SEP)
print("  [3] 反应器几何")
print(SEP)
check("H_bed [m]",   cfg.H_bed,          REF["H_bed_m"],    " m", tol_pct=1.0)
check("D_bed [m]",   cfg.D_bed,          REF["D_bed_m"],    " m", tol_pct=0.1)
check("P [MPa]",     cfg.P/1e6,          REF["P_MPa"],      " MPa", tol_pct=0.1)
check("T_inlet [K]", cfg.T_inlet,        REF["T_inlet_K"],  " K", tol_pct=0.5)
check("ER [-]",      cfg.ER,             REF["ER"],         "",    tol_pct=0.5)

print()
print(SEP)
print("  [4] 气化剂流量（由 ER 反推）")
print(SEP)
# O2 实际流量反推
O2_ref_mol_s   = REF["O2_actual_kg_h"] / 3600.0 / 0.032
N2_ref_mol_s   = O2_ref_mol_s * (79.0 / 21.0)
# 文献蒸汽：反推 steam_to_o2
# Hamel Table 7.1: 水蒸汽 ~5.86 kg/s... 需从文献重新计算
# 用文献 O2 计算化学计量量反推实际 steam
O2_stoich = stoichiometric_o2_mol_s(cfg.fuel_feed, cfg.moisture_wt, cfg.C_dry,
                                     cfg.H_dry, cfg.O_dry, cfg.S_dry)
O2_actual = REF["ER"] * O2_stoich
H2O_from_ratio = 0.8 * O2_actual  # 文献未给出准确蒸汽量，用典型0.8估算

check("O2_stoich [mol/s]",  O2_stoich,              None,              " mol/s", fmt=".4f")
print(f"  ℹ  O2_stoich(计算) = {O2_stoich:.4f} mol/s  (文献: {O2_ref_mol_s:.4f} mol/s,  偏差 {100*(O2_stoich-O2_ref_mol_s)/O2_ref_mol_s:+.1f}%)")
check("O2_feed  [mol/s]",   cfg.O2_feed,            O2_actual,         " mol/s", tol_pct=1.0, fmt=".4f")
check("N2_feed  [mol/s]",   cfg.N2_feed,            N2_ref_mol_s,      " mol/s", tol_pct=2.0, fmt=".4f")
print(f"  ℹ  N2/O2 比 = {cfg.N2_feed/cfg.O2_feed:.3f}  (空气应为 79/21={79/21:.3f})")
print(f"  ℹ  H2O/O2 比 = {cfg.H2O_feed/cfg.O2_feed:.3f}  (模型 steam_to_o2={cfg.steam_to_o2_molar})")

print()
print(SEP)
print("  [5] 粒径 / 固体属性")
print(SEP)
d_p_mean_mm = cfg.d_p * 1000.0
check("d_p [mm] (mean)",  d_p_mean_mm, 2.25, " mm", tol_pct=20.0)  # 文献1.5-3.0mm，取均值2.25
print(f"  ℹ  文献范围 1.5–3.0 mm，当前 d_p = {d_p_mean_mm:.2f} mm")
print(f"  ℹ  n_age_classes = {cfg.n_age_classes}  (文献: 10 discrete size classes)")
check("rho_s [kg/m³]",  cfg.rho_s,   1400.0, " kg/m³", tol_pct=5.0)
check("eps_mf [-]",     cfg.eps_mf,  0.45,   "",       tol_pct=5.0)

print()
print(SEP)
print("  [6] 热损失 / 操作条件")
print(SEP)
check("heat_loss_frac [-]",   cfg.heat_loss_frac,    0.075, "", tol_pct=20.0)
check("recirculation_frac [-]", cfg.recirculation_frac, None, "")
print(f"  ℹ  文献: via Cyclone（旋风分离返料），recirculation=True")

print()
print(SEP)
print("  [7] 汇总：高风险不一致项")
print(SEP)
issues = []
if abs(cfg.fuel_feed - ref_kg_s) / ref_kg_s > 0.005:
    issues.append(f"  ✗  fuel_feed: {cfg.fuel_feed:.5f} vs {ref_kg_s:.5f} kg/s ({100*(cfg.fuel_feed-ref_kg_s)/ref_kg_s:+.2f}%)")
if cfg.n_age_classes != 10:
    issues.append(f"  ✗  n_age_classes={cfg.n_age_classes} (文献要求 10 粒径离散类)")
if abs(d_p_mean_mm - 2.25) / 2.25 > 0.3:
    issues.append(f"  ✗  d_p={d_p_mean_mm:.2f}mm 偏离文献范围 1.5-3.0mm")
if cfg.S_dry == 0.0:
    issues.append(f"  △  S_dry=0.0% (文献 0.51%，影响 H2S 生成与 O2 消耗)")
if cfg.nitrogen_fraction == 0.0:
    issues.append(f"  △  N_dry=0.0% (文献 0.68%，影响 NH3 生成)")
if abs(cfg.steam_to_o2_molar - 0.8) > 0.1:
    issues.append(f"  △  steam_to_o2={cfg.steam_to_o2_molar:.2f} (文献未明确，典型0.5-1.0)")

if issues:
    for iss in issues:
        print(iss)
else:
    print("  ✓  无高风险不一致项")
