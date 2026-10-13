# Hamel Parity 审计 — 收敛不确定性消歧 — 2026-07-14

**证据链**：`data/hamel_1999_extracted.txt`（§2.1–2.3, Tab.2.2, Bild 2.2）
→ `docs/validation_cases_hamel_provenance.md`
→ `docs/audit_reports/hamel_parity_energy_2026-07-06.md` / `hamel_parity_global_nr_2026-07-06.md`
→ 当前 Phase2 LU 工程行为（seed900 + K-refresh）

**范围**：针对「seed900 / Tin / Check1 skip Vorab / final-E 主导格 / split 物种」五类不确定点，
用 Hamel 原文裁决**哪些是方程缺口、哪些是 Startwert/工程桥接、哪些只能数值判定**。

---

## 差异表

| 不确定点 | 分册/论文位置 | 论文式（摘要） | 代码/工程位置 | 判定 | 证据 |
|----------|---------------|----------------|---------------|------|------|
| Tin=293 K | Tab.2.2；`validation_cases_hamel_provenance.md` | 进气 T/P 为 Betriebsdaten；验证表环境 293 K | `case_lu_loader` / `T_inlet` | **一致（边界）** | 禁止抬 Tin 当修复 |
| bed0 seed 850/900 | §2.1 Startwertwahl | Startwerte 由用户生成或沿用旧解；**无**床温闭式 | `vorab_nr_temperature_seed_bottom_K_thesis=900` | **有意偏离（Startwert B 档）** | 非 Eq.2.7 parity |
| Check1 未过不刷 Vorab | Bild 2.2 | Zellenmodell NR 至「Lösung gefunden？」；Abgleich 失败才回 Vorabrechnung | `nr_strict_check1_before_refresh_thesis` | **部分一致** | 流程图支持；中途是否刷 K 论文未写 |
| accept-step K-refresh | §2.1 Vorabrechnung | Vorab 负责水力/流；Zellenmodell 解 inventory | `refresh_cells_hydrodynamics_for_frozen_inner` | **有意偏离（B 档）** | 修 stale K；非 Bild 2.2 明文 |
| energy-regress 强制 refresh | Bild 2.2 | 无「能量回退则立即 Vorab」条款 | `nr_check1_energy_regress_refresh_ratio_thesis=0` | **有意偏离（默认关）** | A/B 伤 split；保持 OFF |
| continuation ΔT 收紧 | Eq.2.9 Wirsum 阻尼 | 远 Startwert 用阻尼；**无**绝对 K 预算 | `nr_inner_abs_dT_budget_*` | **有意偏离（B 档）** | 工程强化；continuation 默认 0 |
| final E 主导格 | Eq.2.7 | 整格焓平衡；**不规定**轴向谁主导残差 | live residual 剖析 | **论文无空间归因** | 须数值环 |
| split 物种尖峰 | Eq.2.4–2.5 | bubble/dense + `ex_bd`；无「N2 应主导」 | gas_phase_split metrics | **论文无物种归因** | 须数值环 |
| Eq.2.7 热损形式 | Eq.2.7 `Q_W`/`Q_WÜ` | 壁面 + 换热器；autotherm 仅环境损失 | `hamel_eq27_explicit` | **一致（见 2026-07-06 审计）** | 能量地板更像 Startwert/流 |

**判定枚举**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离（B 档）` | `论文无归因`

---

## 统计

| 判定 | 条数 |
|------|------|
| 一致 | 2（Tin；Eq.2.7 热损口径） |
| 部分一致 | 1（Bild 2.2 Check1/Check2） |
| 有意偏离（B 档） | 4（seed；K-refresh；energy-regress OFF；ΔT budget） |
| 论文无归因 | 2（final-E 格；split 物种） |
| 不一致 | 0 |

---

## 约束（后续刀不得违反）

1. **不改** `T_inlet`（保持 LU/Hamel 293 K）。
2. **不默认打开** energy-regress refresh / continuation-ΔT 收紧。
3. seed900 在文档中标为 **Startwertwahl 工程选择**，验收看收敛/出口，不宣称温度同构。
4. 下一刀优先：**seed900 下数值重测 final E/S 主导格**；仅当怀疑 `ex_bd`/`K_bd` 时再开 exchange-cell parity。

---

## 建议下一步

1. 跑 seed900 默认的 residual dominance 探针（final + O1）。
2. 若 final E 仍 bed0 → 继续 Startwert/固相焓流出刀；若换位到中上层 → 追该格 H_in/H_out。
3. 可选：把本表链入 `docs/hamel_engineering_bridge_notes.md` Startwert 节。


---

## 数值续测（seed900 默认，2026-07-14）

数据：`data/phase2_seed900_residual_dominance.json`

| 阶段 | E_rms | S_rms | 能量主导格 | 气相尖峰 |
|------|-------|-------|------------|----------|
| pre-O1 | 0.184 | 0.043 | **bed0** | bed1 CO/H2；bed0 O₂ |
| post-O1 | 0.256 | 0.113 | bed3 / bed0 | 混杂；上层出现 N₂ |
| final（8 outer） | **0.093** | **0.047** | **bed1 → bed0**（非 bed9） | bed0 CO/H₂O/O₂ + bed9 H₂O/H₂/CO |

**结论**：文献无法指定主导格；数值表明 **final 能量地板仍在底部 bed0/bed1**（T≈1250/1355 K，Tin 仍 293 K），不是先前 seed850 轨迹里看到的 bed9 优先。

**下一刀**：底部床加热后的焓闭合（冷进气 + 固相 H_out），而非先追自由板/上层。

---

## bed0/bed1 焓分解与温度围栏 A/B（续，同日）

数据：
- `data/phase2_bed01_energy_components.json`（焓项）
- `data/phase2_bed0_fence_hit_ab_v2.json` / `data/phase2_bed0_only_fence_ab.json`
- `data/phase2_fence_lift_energy_dominance.json`

### 证伪表

| 假设 | 结果 | 证据 |
|------|------|------|
| H_Qloss（缺壁面热损） | **否证** | Phase2 已用 `hamel_eq27_explicit`；Q≪\|E\| |
| H_opposite_couple（床间固体焓对消） | **确认** | bed0 E≈+3.8e6 W，bed1 E≈−4.6e6 W；抬 T0 +10 K → dE0≈−dE1 |
| H_fence（seed+350 卡住床温） | **确认卡住；否证为修复** | 见下表 |
| 冷 Tin=293 不可改 | **保持** | 文献 Betriebsdaten；抬 Tin 伤 split |

### 围栏几何

- `tmax = T_vorab + nr_temperature_fence_upper_margin`（march OFF）
- seed900 + margin350 → **bed0 tmax=1250**；final bed0/bed1 **双双顶格**

### A/B（8 outer，单路径 `reactor.solve`）

| 变体 | bed0 T | bed0 顶格? | bed0 E_W | E_rms | S_rms |
|------|--------|------------|----------|-------|-------|
| baseline margin350 | 1250 | 是 | +3.81e6 | **0.093** | **0.047** |
| bed0-only +150（tmax1400） | 1400 | 是 | +5.8e5 | 0.090 | 0.077 |
| bed0-only +250（tmax1500） | 1432 | **否** | ~3e4 | **0.104** | **0.088** |
| 全局 margin500 | 1400 | 是 | +5.8e5 | 0.090 | 0.077 |

抬开围栏后 bed 局部峰值从 ~4.5e6 W 降到 ~7e5 W，主导格从 bed1/bed0 换到 bed2/bed1/bed7…；但 **全局 E_rms 不降、S 明显变差**。当前 1250 围栏对 (E,S) 是 Pareto 保护，不是该默认抬开。

### 约束增补

5. **不默认抬** `nr_temperature_fence_upper_margin` / bed0-only 围栏放宽（已 A/B 否证为修复）。

### 下一刀（已执行 → split；见下节）

---

## Phase-split 物种刀（seed900 默认，同日续）

数据：
- `data/phase2_bed0_split_species_ab.json`
- `data/phase2_bed0_split_term_decomp.json`
- `data/phase2_bed0_split_r5_df_ab.json`

### 基线地图（S=0.047，maxS=0.195）

| 格 | split_rms | 最差物种 | 备注 |
|----|-----------|----------|------|
| **bed9** | **0.105** | H₂O（=maxS） | 与 bed0 **同量级** |
| bed0 | 0.101 | CO | 次之：H₂O、N₂、CO₂、O₂ |
| bed8 | 0.019 | H₂O | |
| bed1 | 0.008 | H₂O | |

此前 top_gas 里的 bed0 CO/H₂O/O₂ 仍在，但 **全局 maxS 实际在 bed9 H₂O**。

### 项分解（Hamel Eq.2.4/2.5 残差）

**bed0 CO**（raw split≈−14）：
- `d_zu=0`（df=0.5 进料对称）
- `d_R≈+19`：dense 产 CO（R_d≈+8.9），bubble 耗 CO（R_b≈−10.3，主要为 R5）
- `2·N_ex≈−29`：交换大力补偿，仍留残差
- → **反应相不对称 + 交换未完全抹平**（非缺 K_bd；K_bd≈13.6 在工作）

**bed9 H₂O/H₂/CO**：
- 反应弱；主导是 **轴向 `d_in` 相不对称 + `2·N_ex`**
- → 上层入口相分裂口问题，不是 bed0 R5

### A/B 证伪

| 旋钮 | S | E | 判定 |
|------|---|---|------|
| df=0.30 | **0.030** | 0.117 | bed0 split≈0；S↓但 E↑；split 上移 bed3/1 |
| df=0.35 / 0.40 / 0.70 | ≥基线 | ↑ | **否证**为稳妥修复 |
| `gas_inlet_o2_dense_frac=0.70` | 0.51 | 0.65 | **否证**（灾难） |
| `r5_bubble_scale_bed0=0.5/0` | ↑ | ↑ | **否证**（T 塌、S 恶化） |

### 约束增补

6. **不默认改** `gas_inlet_dense_frac` / `gas_inlet_o2_dense_frac` / `r5_bubble_scale_bed0`（已 A/B；df0.30 仅作已知 Pareto 点，未升默认）。

### 下一刀

1. **bed9 轴向入口相分裂**（`N_d_in`/`N_b_in` 来自 bed8 出口如何分相）——直接打 maxS。  
2. bed0 仅在有 Hamel R5/`ex_bd` 同构疑点时再开 parity；禁止再盲调 R5 scale。  
3. df=0.30 若要升默认，须更长 ladder 证明 E 可回收。

---

## bed9 轴向入口 / 交换残差刀（同日续）

数据：
- `data/phase2_bed9_axial_split_audit.json`
- `data/phase2_bed9_lower_fence_ab.json`

### 证伪表

| 假设 | 判定 | 证据 |
|------|------|------|
| H_axial_BC（`N_in≠prev.N`） | **否证** | bed7–9 `err_d/b=0`；`routing._apply_bed_upstream_inlet` 同构 Hamel |
| H_bed8_out_split（入口相分比异常） | **否证为修复** | `f_d≈0.70` 沿床连续；强制 50/50 入口使 H₂O split −15→−22 **更差** |
| H_ex_bed9（交换主导 maxS） | **确认** | bed9 H₂O：`split_noex≈−0.83`，`2·N_ex≈−14.7`（≈95%） |
| H_bed9_T_floor（下围栏） | **确认卡住；否证为修复** | 见下表 |

### 机制摘要

- 轴向 BC 正确：`N_d_in=prev.N_d`，`N_b_in=prev.N_b`。
- bed9 H₂O/H₂/CO 的 scaled split 几乎全是 **`2·N_ex = 2 K_bd V_b (C_b−C_d)` 留在残差里**；对流+反应已近平衡。
- bed9 **T=951.4 = tmin = vorab−150**（下围栏顶格）；与 bed0 上围栏对称。

### 下围栏 A/B（8 outer）

| 变体 | T9 | at_lo | 2Nex H₂O | E | S | maxS |
|------|----|-------|----------|---|---|------|
| lo150（默认） | 951 | 是 | −14.7 | **0.093** | **0.047** | **0.195** |
| lo80（抬高地板） | 1021 | 是 | **−0.69** | 0.155 | 0.107 | 0.407 |
| lo250/300（放冷） | 875 | 否 | −167 | 0.464 | 0.471 | 2.18 |

抬高地板能抹平 bed9 `N_ex`，但全局 S/E 变差；放冷则气泡 H₂O→0、全面崩溃。当前 lo150 仍是 Pareto 护栏。

### 约束增补

7. **不默认改** `nr_temperature_fence_lower_margin_K_thesis`（已 A/B）。
8. **不改** 轴向气相 `N_in←prev.N` 接线（与 Hamel 一致）。

### 下一刀（已执行 → C_b/C_d / 自由板能量；见下节）

---

## bed9 `|dy|` / 自由板能量假根（2026-07-16）

数据：
- `data/phase2_bed9_CdCb_Vb_audit.json`
- `data/phase2_bed9_freeboard_T300_audit.json`
- `data/phase2_fb0_energy_false_root.json`
- `data/phase2_fb_energy_passthrough_postfix.json`

### 证伪 / 确认

| 假设 | 判定 | 证据 |
|------|------|------|
| H_Vb_hydro（V_b 暴涨） | **否证** | bed7→9：V_b 略降；K_bd·V_b 仅 1.86→2.20 |
| H_CdCb_thermo（浓度公式错） | **否证** | `N_ex = K_bd V_b c_fac (y_b−y_d)` 数值闭合 |
| H_dy_growth（\|dy\| 主导） | **确认** | bed9 \|dy_H2O\|≈0.011 ≫ bed7 的 0.0005 |
| H_fb_T300_false_root | **确认并修复** | 见下 |

### 根因链

1. `freeboard_closure` 固体不在 NR；能量却用 `K·m_solid` 作出流 → **`m_in≠m_out` 固体焓洞 ~2 MW**。  
2. NR 把自由板 T 压到 ~300 K，用冷气出口焓「抵消」该洞（fb0 在 Tin=951 时 F_T≈0）。  
3. `T_solid_ab_in=300 K` 回流冷却 bed9 → 顶住下围栏 → \|dy\| / `N_ex` 爆炸 → maxS。

### 最小修复

`Cell._solid_energy_outlet_rates()`：`freeboard_closure` 下 Eq.2.7 `H_aus` 固体载体改为**入口固体透传**（在 cell T 下出流）。  
运输/返料仍用 `_solid_downflow_rates()`（不改）。  
测试：`tests/test_energy_balance_parity.py::test_freeboard_closure_energy_uses_solid_inlet_passthrough`

### 修复后（8 outer，单路径 solve）

| 指标 | 修复前 | 修复后 |
|------|--------|--------|
| E_rms | 0.093 | **0.015** |
| S_rms | 0.047 | **0.028** |
| maxS | 0.195 | **0.098** |
| fb@300K | 3/8 | **0/8** |
| bed9 下围栏 | 顶格 | **离开**（T≈978） |
| Tab→bed9 | 300 K | **945 K** |
| dy_H2O / N_ex | 0.011 / −7.3 | **0.0002 / −0.13** |

### 约束增补

9. freeboard_closure **能量**固体出口必须与入口质量闭合；禁止再用 `K·m` 当 Eq.2.7 载体。  
10. 运输返料质量仍走 trajectory/`_solid_downflow_rates`，与能量透传分离。

### 下一刀（已执行 → ladder20 + 剩余 maxS 地图；见下节）

---

## 修复后 ladder20 + 剩余 maxS 地图（2026-07-16）

数据：`data/phase2_postfix_ladder20_split_map.json`

### 稳定性

| outer | E | S | maxS | converged | fb@300 | bed9 下围栏 |
|------:|---|---|------|-----------|--------|-------------|
| 8 | 0.0151 | 0.0280 | 0.0983 | **True** | 0 | 否 |
| 20 | 0.0151 | 0.0280 | 0.0983 | **True** | 0 | 否 |
| Δ(8→20) | **0** | **0** | **0** | — | — | — |

自由板能量透传修复后，**8 outer 已定态**；加长不回弹。

### 剩余 split 地图（主导已离开 bed9）

| 排名 | 格 | split_rms | 最差物种 | 备注 |
|-----:|----|-----------|----------|------|
| 1 | **bed0** | 0.0606 | **CO₂** (−0.098=maxS) | 仍顶 **上围栏 T=1250** |
| 2 | **bed1** | 0.0585 | **H₂** (−0.096) | 次之 H₂O/CO |
| 3 | bed2 | 0.0202 | H₂O | |
| 4+ | bed3… | ≤0.011 | — | bed9 不再进前排 |

顶物种：bed0 CO₂/H₂O/H₂/O₂ + bed1 H₂/H₂O/CO。

### 结论

- 自由板假根链已闭合；当前地板回到 **底部 bed0/bed1 相分裂**（与 bed0 上围栏耦合）。  
- 全局 `converged=True`（solver 口径）；**勿**等同 Fig.7 / 出口 LU 验收。

### Outer preinner bed1 T2（2026-07-16，单路径 `reactor.solve`）

| 配置 | E | S | maxS(map) | merit | bed0_T |
|------|--:|--:|----------:|------:|-------:|
| bed1 OFF（默认，修复前） | **0.015** | 0.028 | 0.061 | 0.028 | 1250 |
| bed1 ON + hotter-than-bed0 门控 | 0.026 | **0.023** | **0.037** | **0.026** | 1138 |
| bed1 ON 每轮强制 T2 | 0.031 | 0.028 | 0.037 | 0.031 | 1062 |

- 修复：`include_bed1=True` 在 `max_bed_index=0` 时也能把 bed1 列入候选（此前为死旗标）。  
- **不**把 Phase2 默认改为 ON：merit 略好但 E 近似翻倍。

### bed0 CO₂ `N_ex` 缺口（同日闭合）

**根因链**（非 Hamel 方程错误）：

1. 事后 bed0 T2 可将 `N_ex(CO₂)` −0.05→**+3.78**、`dy`→−0.015，床0 split→0。  
2. 晚期 outer（O7）全局 S≈0.028 **略低于** `tau_split=0.03` → **跳过 T2**，而 `max|split|≈0.098` 仍在 bed0 CO₂。  
3. 中期（O6）T2 常被 global merit gate 回滚（抬 E 超过降 S）。  
4. NR 落到 combined≈0 / equal-y 吸引子，留下相分裂地板。

**修复**：`resolve_preinner_tier_plan` 增加 `gas_phase_split_max_abs ≥ 0.05` 亦可触发 T2  
（`nr_outer_preinner_split_tier_min_max_abs_thesis`）。

| 8 outer | E | S | bed0 split | Nex CO₂ | merit | hot |
|---------|--:|--:|-----------:|--------:|------:|-----|
| 修复前 | 0.015 | 0.028 | 0.061 | −0.05 | 0.028 | bed0 |
| +max_abs T2 | 0.025 | 0.021 | **≈0** | **+3.78** | 0.025 | bed1 |
| +bed1 late follow-up | 0.025 | **0.012** | ≈0 | +3.78 | 0.025 | bed2 |

探针：`data/phase2_split_maxabs_trigger_o8.json`、`data/phase2_bed1_hotspot_followup_gated_o8.json`。

### bed1 H₂ late-floor follow-up（同日）

**根因**：选床在 bed0 T2 **之前**；bed0 仍热时 bed1 进不了候选；bed0 清掉后 bed1（H₂）成唯一热点。

**否证**：中期无条件 bed1 T2 → E↑、T 离围栏（merit 0.031）。  
**修复**：`nr_outer_preinner_bed1_hotspot_followup_thesis` —— 仅当  
`bed0_split < 0.01` **且** `global S < 0.04` **且** `bed1 ≥ 0.05` 时做 bed1 T2。

O7：`b0=true, b1=true`；O4–6：`b1=false`。事后探针 E 几乎不变、S 0.021→0.012。

### bed2 follow-up 否证（2026-07-16）

事后 bed2 T2：S 0.012→0.010，E≈不变，**merit 仍 E 主导（0.025）** → `worth_for_merit=false`。  
不落地 bed2 late follow-up。探针：`data/phase2_bed2_t2_post_probe.json`。

### 当前收敛地板 vs LU 出口（同日）

| 量 | 值 | 备注 |
|----|-----|------|
| merit / E / S | 0.025 / 0.025 / 0.012 | E 主导；bed0 顶 1250 K，E₀≈1.71 MW |
| bed_exit 干基 CO/CO₂ | **0.126 / 0.150** | vs LU 0.13 / 0.11（床顶尚可） |
| exit 干基 CO/CO₂ | **0.048 / 0.209** | 来自 **cyclone / NR 显式自由板** |
| 解析 FB（终态床顶重算）CO | **≈0.109** | 几乎不烧 CO；与 NR 自由板分叉 |

**自由板诊断修复**：`result_exit_summary` 不再使用中途缓存的 `_last_explicit_freeboard_closure` 填剖面（曾出现虚假 profile[0] O₂≈0.06）；改为始终用终态床顶重跑 `simulate_freeboard`。

**NR 显式自由板轴向（终态）**：fb0 CO≈0.109 → fb7 CO≈0.045，同时 H₂↑（WGS 型），O₂≈0。二次风 O₂ 仅 0.078 mol/s（xi=0.6），不是主因。  
解析闭包同入口几乎平坦 → **下一刀是 NR 自由板 vs Hamel 解析闭包的气相化学同构**。

探针：`data/phase2_energy_outlet_detail_o8.json`、`data/phase2_fb_inlet_vs_profile0.json`、`data/phase2_freeboard_co_ab_o8.json`。

### 下一刀

1. **NR 自由板 WGS/均相 vs `simulate_freeboard`**：同床顶入口为何 NR 把 CO→CO₂+H₂、解析几乎不动。  
2. bed0 能量地板（围栏下 E≈1.7 MW）——merit 要降必须动 E。  
3. 禁止抬 Tin / 默认抬围栏；**勿**默认全量 bed1/bed2 T2。
