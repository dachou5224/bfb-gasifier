# 交接文档：反应源项 Hamel 同构下的 NR 收敛推进

**日期：** 2026-08-06  
**状态：** 未 fully converge；平台期已定位到 merit 景观，非阻尼系数  
**语言：** 接手 agent 请用简体中文回复用户  
**勿改：** 用户明确要求勿改计划 `.plan.md`（若存在）

---

## 0. 给接手 agent 的 60 秒摘要

主线是床层 **关掉 extent 限幅**（Hamel 同构）后让 HTW LU Phase2 全局 NR 收敛。  
稳定默认栈（gate+hot + LSQ n=2，H2O 透传关，`PYTHONHASHSEED=0`）在 **max-merit 时代**约 CO≈0.159 / rms≈0.046。  
**Wirsum Eq.2.23 已落地默认**（关 split/energy max-merit）。  
**Hamel+Wirsum 核心栈**（关 LS fastox / syngas / PS / GATE / HOT / LSQ）：o8 **未收敛**，**n_acc=0**，CO≈0.012；证据 `data/phase2_hamel_wirsum_core_o4_o8.json`（§5.43）。  
此前 CO≈0.116/0.159 是叠在平台桥上的结果，不能当同构验收。  
方法：逐 cell T/组分 → 追因。底格反应/传质/传热审计见 §5.48–5.49：根因是 Startwert 不在 F=0 流形上。  
**§5.49–5.64**：两相 Startwert + 按相 fastox + 分相切向 clip + 底格 R1 残氧种子 + 每步 |ΔT| 帽 40 K + 种子后固相对齐 + bed1 局部 snap + 种子后 N₂/H₂O/CO₂ 重分相 + 底格围栏下沿 T_ref−350≈800 K + **底格 Startwert T=820 K**。CORE o8 rms **≈0.012**，底格 T≈**840 K**。λ 历史长度 = inner budget（o8=40），其中 **非 None 接受步 ≈25**。§5.61：bed0 焓洞是 VM zone=3 + `K_VM=0`。§5.62：H₂O 透传砸栏。§5.63：产物透传填 E3 但 o8 rms 变差。§5.64：840 K 是底格 T 吸引子；o8 能量平台是 **budget 用尽**（线搜索仍在 λ=1 接受）；o16 energy_rms 0.045→0.035；bed7–9 贴上段 1000 K 栏，加宽到 800 只换来小数点后第四位，CORE 不上段围栏。勿用 850，勿开 H₂O/产物透传当 CORE 默认，勿把 VM 补进 Eq.2.7。勿开合计 LS fastox / PTC / L4，勿开「允许 O₂ 进燃料相」全床 clip。勿把全局 Tmin 放到 600 K 以下。勿开 `auf_upflow_closure`。勿对 O₂ 套惰性 Eq.2.2。

---

## 1. 用户目标与硬约束

| 项 | 内容 |
|----|------|
| 目标 | 反应源项 Hamel 同构（去 extent 限幅）下推进收敛 |
| 床层限幅 | **不恢复** `extent_limiters_enabled_thesis` |
| 可用杠杆 | Wirsum 阻尼、初值、自由板 dt 界；床层/自由板政策可分离 |
| **诊断方法** | **先**全炉逐 cell 查 T/组分分布异常，**再**对异常格做源汇/守恒追因；**禁止**以抬单组分（如 CO）或压 rms 为唯一成功标准 |
| 语言 | 简体中文 |
| 项目规则 | 先读 `docs/CLAUDE.md`；SI；Arrhenius 仅工厂函数；勿滥加 logging/抽象 |

相关技能（按需）：

- `.cursor/skills/bfb-diagnosing-bugs/SKILL.md` — 先建反馈环再假设  
- `.cursor/skills/bfb-hamel-parity-audit/SKILL.md` — 方程级 parity  
- `.cursor/skills/bfb-numerical-guard/SKILL.md` — 求解器数值清单  

---

## 2. 已落地改动（生产默认方向）

### 2.1 同构与 Startwert

| 开关 / 位置 | 含义 |
|-------------|------|
| `ReactorConfig.extent_limiters_enabled_thesis=False` | 床层限幅关 |
| `freeboard_extent_dt_bounds_thesis=True` | 自由板 dt 界（与床层解耦） |
| `vorab_transport_x0_fast_oxidation_closure_thesis=True` | 快氧化 x₀：R12→R5→R6→R10；合计 holdup 氧化再按水力重分 |
| `vorab_transport_x0_char_oxidation_o2_closure_thesis=False`（CORE） | 不再把 holdup 氧置零（§5.48） |
| `vorab_bed0_eq24_two_phase_startwert_thesis=True`（CORE） | 底格 Eq.2.4 总量 + O₂/合成气分相（§5.49） |
| `nr_line_search_per_phase_fastox_thesis=True`（CORE） | 线搜索按相 fastox，非合计（§5.50） |
| `nr_clip_same_phase_oxidizer_fuel_step_thesis=True`（CORE） | 禁止同相增加 O₂∩合成气（§5.51） |
| `nr_clip_allow_oxidizer_into_fuel_thesis=False` | 放开全床 `dN_O2>0` 会把试探打到 6.7/132（§5.52） |
| `vorab_bed0_r1_oxidizer_seed_mol_s_thesis=2.0`（CORE） | 底格气泡→悬浮相残氧种子，点亮 R1（§5.52） |
| `nr_inner_t_step_default_cap_K_thesis=40.0`（CORE；生产 600） | 每步 |ΔT| 帽，避免底格一刀砸进围栏（§5.53） |
| `nr_temperature_fence_bed0_lower_margin_K_thesis=350`（CORE；生产 None） | 仅底格 Tmin=T_ref−350≈800 K（§5.58） |
| `nr_temperature_fence_upper_bed_lower_margin_K_thesis=None`（CORE） | 上段 Tmin；opt-in 350 让 bed7–9 离开 1000 K 栏，o8 增益很小（§5.64） |
| `vorab_bed0_nr_startwert_T_K_thesis=820`（CORE；生产 None） | 底格 NR 初温 820 K（§5.60；880 见 §5.59） |
| `vorab_vm_devolatilization_zone_cells_thesis=3`（CORE 钉死） | Hamel Kap.2.1 轴向 VM 释放；zone=1 填焓洞但 o8 变差（§5.61） |
| `vorab_a_tier_post_staged_h2o_passthrough_thesis=False`（CORE） | bed3 H₂O 透传；opt-in 在 fastox 后再贴（§5.62）。打开会砸进 800 K 栏 |
| `vorab_a_tier_post_staged_syngas_passthrough_thesis=False`（CORE） | bed3+ CO/CO₂/CH₄/TAR 透传；opt-in 在 fastox 后再贴（§5.63）。填 E3 但 o8 rms 变差 |
| 种子后 `reconcile_upper_bed_solid_holdup`（`allow_reduce=True`） | fastox/R1 后重贴 Eq.2.6，堵 bed1/bed9 焓洞（§5.54） |
| `_snap_bed_holdup_to_local_support(bed1)` | 全床联立后再单独贴 bed1 char/ash（§5.55） |
| `apply_bed0_inert_exchange_split_startwert` | 种子/snap 后只重分 N₂/H₂O/CO₂（§5.56；同一 eq24 开关） |
| `vorab_init_solve_bottom_cell_thesis=False` | 底格 `solve_cell`；默认关（残氧+R12 会卡住） |
| 实现 | `src/solvers/vorabrechnung/vorab_x0.py`；在 abgleich **之后**调用（`init_precalc_step.py`） |
| 测试 | `tests/test_vorabrechnung_hamel_x0_transport.py` |

### 2.2 NR 线搜索与 Pre-Inner

| 开关 | 含义 |
|------|------|
| `nr_line_search_gas_phase_split_merit_thesis=False` | **Wirsum Eq.2.23**：默认不用 split 进 max-merit |
| `nr_line_search_energy_merit_thesis=False` | **Wirsum Eq.2.23**：默认不用 energy 进 max-merit；接受标量 = 缩放残差 RMS |
| `nr_line_search_fast_oxidation_projection_thesis=False` | LS 试探快氧化投影（工程桥；**默认关**。Vorab x₀ 快氧化仍开） |
| `nr_line_search_syngas_collapse_guard_thesis=False` | 拒绝合成气伪下降（工程桥；**默认关**） |
| 接线 | `global_nr_iteration.py` / `nr_inner_step.py` / `reactor.py`；ladder `p0/p1/p2` 已对齐关 max-merit |
| Phase2 | `nr_outer_preinner_bed1_phase_split_thesis=True`，`max_bed_index=2`（见 `tests/validation/config_phase2_freeboard.py`） |
| 全局 split 守卫 | 单格 LSQ 若抬高全局 split → 回滚（`preprojection.py`） |
| 注意 | 接受步若水力刷新 → 丢弃缓存 J；split/energy max-merit 仅 opt-in |

### 2.3 关键代码锚点

- 线搜索 / 阻尼：`src/solvers/global_nr_iteration.py`（约 398–620 行）  
- merit / 接受谓词：`src/solvers/global_nr_solver.py` — `_line_search_merit`、`_line_search_trial_acceptable`  
- Wirsum 旋钮：`src/solvers/wirsum_damping.py`  
- 配置默认：`src/core/reactor.py`（`nr_*_thesis`）  
- LU harness：`scripts/_lu_harness.py`  

---

## 3. 当前数值状态（HTW LU Phase2，同构限幅关）

| 阶段 | o8 量级 | converged |
|------|---------|-----------|
| 仅去限幅 / 裸 R12 | 天文 / stall | False |
| +快氧化 x₀ | rms ~0.72→0.14 | False |
| +R10+重分 | split~0.07, energy~0.076 | False |
| +LS fastox | split~0.05–0.07, energy~0.027 | False |
| +bed1 preinner（现默认） | **split~0.037**, energy best~**0.013** | **False**（平台） |

**注意：** JSON 里 `rms_scaled_energy_final` 常是 **best-so-far**；晚期活状态 energy 可回弹到 ~**0.040** 并主导 merit。验收时分开报告。

证据文件（仓库 `data/`）：

- `phase2_extent_raw_damping_vs_lsq_ab_o8.json`  
- `phase2_extent_raw_damping_step_model_ab_o8.json`  
- `phase2_extent_raw_ls_accept_probe_o8.json`  
- `phase2_extent_raw_ls_fail_root_o8.json`  
- `phase2_extent_raw_ls_merit_decomp_o8.json`  
- `phase2_extent_raw_energy_desc_accept_ab_o8.json`  
- 更早：`phase2_extent_raw_bed1_preinner_o8.json`、`*_o20.json`、`*_fast_ox_*` 等  

复现风格：用 `_lu_harness.build_phase2_reactor` + `run_phase2_init` + `phase2_solve_kwargs(max_global_iter=8)`；SEED≈1250 K bed0 fence（见历史 throwaway harness）。

---

## 4. 收敛诊断结论（已钉死）

### 4.1 阻尼 A/B：无效

关 LSQ 后，`prefer_full_step=False`、`λ₀=0.25/0.125`、更多 halvings、关 split-merit、`lm`/`ptc` → 末态**锁同一点**（split≈0.0376, energy≈0.070）。  
原因：LS **零接受步**，不是 λ 太小。

LSQ on + soft damp：split 略降（~0.0356），energy 变差（~0.028）。

### 4.2 线搜索失败机制

- baseline：失败率约 20/32；晚期 best_trial merit≈0.045 **>** 当前≈0.037  
- 准确插桩（`_line_search_trial_acceptable`）：拒绝原因 **100% `merit_not_decreased`**  
- 典型交易：

```text
当前:  merit≈0.040 (energy 主导), split≈0.037, energy≈0.040
最佳试探: merit≈0.045 (split 主导), split≈0.045, energy≈0.022, rms≈0.017
→ energy↓ 但 split↑ → max-merit 上升 → 拒
```

### 4.3 「energy 降即可接受」已证伪

门控政策（energy 主导且 energy_trial < energy_F 则接受；含 split 软帽）→ split/energy **相对 baseline 变差**。  
`nr_refresh_hydrodynamics_on_accepted_step_thesis=True` **不能**自动修复抬高的 split。  
与 `_line_search_trial_acceptable` 注释中的 post-F3 病理一致。

### 4.4 LSQ 定位

- **非 Hamel**：论文无 holdup LSQ 处方  
- **工程上仍有效**：压 energy、提供可走轨迹  
- 勿再扩大 maxbed；保持全局 split 守卫  

---

## 5. 建议的下一刀（接手优先序）

1. **压 rms / 收敛** — §5.11–5.12：CO 桥 OK；rms 主导=bed2 H2O/H2；**下一刀：床2 耦合或避开 bed6 local2**  
2. **平台期 cell×物种 / N2 / 相源 / absorb / clip** — §5.2–5.6 已钉死  
3. **数值 A/B（refresh-aware / split 约束）** — 已证伪为默认，见 §5.1  
4. **文献：** 优先 Wirsum 1998；次选 Wozny 1983、1999 Final Report  
5. **验收：** 分开 best-so-far vs 末态活指标；勿恢复限幅当物理修复  

**不要做：** 盲目扩大 LSQ / maxbed9 冻总量；宣称「阻尼同构已对齐」；无证据改 merit 权重当最终解；默认开 local2x2。

### 5.1 Refresh-aware / split 约束 A/B（已证伪为默认）

开关（`reactor.py`，**默认 False**；工程桥）：

- `nr_line_search_refresh_aware_accept_thesis`
- `nr_line_search_split_constrained_step_thesis`

实现：`_refresh_aware_joint_acceptable` / `_zero_gas_holdup_newton_step`；LS 接线在 `global_nr_iteration.py`。

证据：`data/phase2_refresh_aware_split_constrain_ab_o8.json`（当前 T-rescale 栈，o8）

| mode | CO | split | energy | 备注 |
|------|---:|------:|-------:|------|
| baseline | 0.041 | 0.050 | 0.042 | 4 接受 / 21 LS fail |
| refresh_aware | 0.041 | 0.050 | 0.042 | **61 probe / 0 accept** |
| split_constrain（冻气相 dx） | 0.011 | 0.032 | 0.058 | CO 锁 init；有害 |
| both | 0.011 | 0.032 | 0.058 | 同 split_constrain |

插桩：61 次 refresh 探针中 **merit_post 无一 ≤ merit_F**（最好仍 +4.6e-4）；counterfactual「只看 merit」亦为 0。  
结论与 §4.3 一致：水力刷新**不能**修复 energy–split 交易；冻全气相 Newton 步会毁掉 CO 轨迹。  
勿开生产默认。下一刀回到 **cell×物种 split 主导项**（或文献）。

### 5.2 平台期 cell×物种 split 主导项（2026-08-07）

脚本：`scripts/audit_split_plateau_cell_species.py`  
证据：`data/phase2_split_plateau_cell_species_o8.json`（同构栈 o8，CO≈0.041）

**物种方差（split 项 sq 份额）**

| 物种 | var 份额 | 备注 |
|------|--------:|------|
| **N2** | **0.45** | 第一主导（惰性相分配） |
| **CO** | **0.37** | 与 N2 合计 ≈0.82 |
| H2 | 0.14 | 第三；床层 H₂ holdup 仍近 0 |
| H2O/O2/CO2/CH4 | <0.03 each | 非平台主因 |

**轴向方差**

| 区域 | var 份额 | 最差物种 |
|------|--------:|----------|
| bed6–9（上段） | **0.64** | N2 / bed9 为 CO |
| bed9 单格 | 0.23 | CO≈N2（|Δ|≈0.129） |
| bed0–2（下段） | 0.13 | bed1 仅 0.026（H2） |

**与旧叙事的差异**

- 旧平台 dissect / bed1 preinner 假设热点在 **bed0/1 × CO/H2**。  
- **当前同构平台** split 质量在 **上段床 × N2/CO**；继续扩 bed1 LSQ **打不中** 主方差。  
- 末态活 merit 由 **rms/gas** 主导（energy≈0.041 < split≈0.050）；handoff §4.2 的 energy↔split 交易在此点很弱（Newton+fastox λ 梯仅 ΔE~1e-6，却抬 split）。

**建议下一刀**（已落地 §5.3–5.4）

1. ~~反应物种相分配~~ → 见 §5.4：密相源 + 库存未吸收。  
2. 勿从 split merit 剔除 N2（§5.3）。  
3. 文献仍缺 Wirsum/Wozny 时，保持工程桥标注。

### 5.3 N2 大 split：不是误装反应，是惰性放大相组成差（2026-08-07）

证据：

- `data/phase2_upper_bed_n2_co_split_term_decomp.json`
- `data/phase2_n2_split_mole_fraction_decomp.json`
- `data/phase2_n2_correct_closure_and_gate_ab_o8.json`

| 检查 | 结果 |
|------|------|
| `R_gas[N2]` | **全床 = 0**（未进 R1–R12 / 焦油计量） |
| 浓度定义 | `C_i = y_i P/(R T)`（非 N/V） |
| 上段 N2 | `y_d−y_b ≈ +0.009`，`\|2 N_ex\| ≈ \|res_split\|`，`res_sum≈0` |
| bed9 Δy 主导 | N2 (+0.009) > H2O (−0.004) > CO (−0.003) |

机制：反应气（CO/H2O 等）两相分配不均 → 密相 N2 摩尔分数略高 → 大流量惰性 × `K_bd V_b` 把小 Δy 放大成大绝对 split。N2 在守恒式中正确；根因在**反应物种相分配**，不是 N2 误入反应方程。

A/B（求解期从 split 门禁去掉 N2）：CO 仍 ≈0.041，full-gate split 仍 ≈0.050——只改报表口径，不推进轨迹。

下一刀：上段 **CO/H2O 相间交换与源项相定位**（非删 N2）。

### 5.4 上段 CO/H2：密相源 + 库存未吸收（2026-08-07）

脚本：`scripts/audit_upper_bed_co_h2_phase_source.py`  
证据：`data/phase2_upper_bed_co_h2_phase_source_o8.json`

| 量 | bed6–9 |
|----|--------|
| split \|rs\|² 份额 | N2 0.48 / **CO 0.40** / H2 0.12 / H2O≈0.003 |
| CO `R_d/(R_d+R_b)` | **≈1.00**（几乎全在密相；气泡源≈0） |
| H2 同口径 | **≈1.00**（fastox 后 H₂ holdup≈0，`N_ex≈0`） |

**CO split 恒等式**（`res_d−res_b = in_diff + R_diff − out_diff + 2 N_ex`）均值：

| 项 | 均值 [mol/s] | 角色 |
|----|-------------:|------|
| `R_diff` | **+4.3** | 密相净产 CO |
| `2 N_ex` | **+4.1** | 气泡略富 CO → 交换进密相，**与 R_diff 同号叠加** |
| in/out_diff | ~0.09 / 0.07 | 可忽略 |

同时 `res_sum ≈ R_net`（总量未闭合）：密相产 CO **未进** `N_d` 库存/出流。  
组成悖论：`y_d,CO − y_b,CO ≈ −0.004`（气泡略富），与「密相产 CO」相反——是库存滞后的症状，不是把异相反应误接到气泡。

H2：split ≈ `R_diff`（纯密相源；`N_ex≈0`）。  
H2O：密相大汇，但 split 部件对消，**不是**上段方差主因。

体积份额重分源的反事实（诊断，非 Hamel）：CO split 可降约一半，H2 近消——证明代数杠杆在 **R 相定位 / 密相库存吸收**；异相 R1–R4 / R7–R8 **必须留在密相**，禁止当生产补丁。

**建议下一刀**（已做 §5.5）

1. ~~密相 holdup 吸收探针~~ → 见 §5.5。  
2. 勿扩 bed1；勿把炭反应搬进气泡。

### 5.5 上段 dense-CO 吸收 A/B（2026-08-07）

脚本：`scripts/audit_upper_dense_syngas_absorb_ab.py`  
证据：`data/phase2_upper_dense_syngas_absorb_ab_o8.json`  
实现（默认关）：`nr_outer_preinner_upper_dense_syngas_absorb_thesis` + `oneshot_bed_dense_species_sum_absorb`  
说明：非 bed0 phase-split LSQ **冻住** `N_d+N_b`，代数上无法关 `res_sum`。

| arm | exit CO | rms | split | energy | 备注 |
|-----|--------:|----:|------:|-------:|------|
| baseline | 0.0327 | 0.0502 | 0.0502 | 0.0417 | |
| absorb 全量 CO+H2 | 0.0327 | 0.0502 | 0.0502 | 0.0417 | 局部门禁 **全程拒收**（oneshot 抬 split→2.9） |
| absorb_soft_CO α=0.1 | 0.0327 | 0.0502 | 0.0502 | 0.0417 | 局部门禁仍拒收 |
| **soft_CO α=0.1 + skip local gate** | **0.0381** | **0.0479** | **0.0410** | 0.0417 | 2/10 outer 经全局门禁留下；上段 \|res_split\| 8.5→3.3 |
| maxbed9（冻总量 LSQ） | 0.0041 | 0.0629 | 0.0333 | 0.0629 | split↓ 但 **CO 塌缩**；有害 |

结论：

1. **全量吸收 / 带局部 split 门禁**：证伪（不改轨迹）。  
2. **maxbed9 扩 LSQ**：证伪为生产路径（杀 CO）。  
3. **软吸收 CO×0.1 且跳过局部门禁、靠全局 merit**：弱正号（CO/split 略好），仍远低于旧平台 CO≈0.15；`res_sum` 仍大。工程桥，**勿默认开**。  
4. 下一刀：查为何 NR 不自行抬 `N_d,CO`（Jacobian / fastox / clip）→ §5.6。

### 5.6 为何 NR 不抬上段 `N_d,CO`（2026-08-07）

脚本 / 证据：

- `scripts/audit_upper_dense_co_newton_block.py` → `data/phase2_upper_dense_co_newton_block_o8.json`
- `scripts/audit_solid_ftb_relative_floor_ab.py` → `data/phase2_solid_ftb_relative_floor_ab_o8.json`
- 实现（默认关）：`nr_clip_solid_ftb_relative_floor_thesis`

**平台态一步解剖（同构 o8）**

| 量 | bed6–9 dense CO |
|----|-----------------|
| 理想吸收 `res_sum` | 3.2–6.4 mol/s |
| Newton `dx_raw` | 0.15–0.18（仅 **~4%** 理想） |
| legacy `dx_clip` | **~7e-5**（α≈**4.7e-4**） |
| fastox 对 CO | Δ≈0（非阻塞） |

**α 屠杀根因**：bed1 某固相类库存 ~**1.38e-6 kg**，Newton 给负步 → 固相 FTB `α = 0.8·ξ/|dx|` 进**全局** α，把全部气相步（含上段 dense-CO）缩到近零。  
不是 fastox；不是 dense-CO 列奇异（`J_dd≈-0.12`，有信号）。

**相对地板 opt-in**（`solid_ftb_relative_floor`：地板 `max(1e-6, 1e-4·max m)`，对齐气相）

| | α | dense-CO `dx_clip` | 全炉 o8 |
|--|--:|--:|--|
| legacy | 4.7e-4 | ~7e-5 | CO 0.0327 / split 0.050 |
| relative floor | ~0.14 | ~0.02–0.026 | CO 0.0327 / split 0.050（**轨迹几乎不变**） |

平台线搜索（relative）：λ≤0.03 可降 merit，但 split 几乎不动；λ≥0.125 因 split 拒收。

结论：

1. **Clip 全局 α 被微量固相掐死** = 真 bug / 工程缺陷；相对地板可解除平台 oneshot，**勿默认开**（全炉 A/B 无出口 CO 收益）。  
2. **更深阻塞**：即便放开 clip，Newton 只请求 ~4% 理想吸收，且 LS 在较大 λ 上仍被 split 挡住——与 §5.4「`R_diff`+`2 N_ex` 同号」一致。  
3. 下一刀优先：**联合密相库存 + 气泡相 CO（局部 2×2）**，见 §5.7；单改 clip 不够。

### 5.7 上段 CO 局部 2×2（J 低估 vs 联合闭合）（2026-08-07）

脚本 / 证据：

- 平台解剖：`data/phase2_upper_co_joint_nex_oneshot_o8.json`、`data/phase2_upper_co_local2x2_soft_o8.json`
- 全炉 A/B：`scripts/audit_upper_co_local2x2_ab.py` → `data/phase2_upper_co_local2x2_ab_o8.json`
- 实现（**默认关**）：`nr_outer_preinner_upper_co_local2x2_thesis`（+ `min_bed_index` / `fraction`）  
  → `run_outer_preinner_upper_co_local2x2`（`preprojection.py`），挂入 `run_outer_preinner_start_value_bridge`；`outer_loop` 记 `upper_co_local2x2_*`

**平台态：全系统 Newton vs 局部 CO 2×2（bed6–9）**

| | dense-CO 步 / 理想吸收 |
|--|--:|
| 全系统 Newton `dx_full` | **~4%** |
| 局部 CO 2×2 `dx2`（FD 列与 J 一致） | **~50–56%** |
| 单变量 `dx1`（只动 dense） | ~0.55 mol/s 量级，不足 |

朴素「joint Nex」迭代（按 `Nex_tgt` 强拧相份额）：**发散**（`Nb→0`，`Nex→−55`，split→0.54）——失败路径，勿再试。

平台 oneshot soft apply（不动全炉 NR）：

| frac | merit | split | \|res_split\| mean | 门禁 |
|-----:|------:|------:|------------------:|:----:|
| 0（基线） | 0.0545 | 0.0502 | ~8.5 | — |
| **1.0** | **0.0538** | **0.0359** | **~1.0** | ✓ |
| 0.5 | 0.0787 | 0.0412 | ~4.6 | ✗（非线性；半步更差） |

**全炉 o8 A/B（max_global_iter=8）**

| arm | CO | rms | split | energy | local2 接受 | n_acc |
|-----|---:|----:|------:|-------:|------------:|------:|
| baseline | 0.0327 | 0.0502 | 0.0502 | 0.0417 | 0/0 | 8 |
| **local2x2** | **0.0705** | 0.0632 | 0.0603 | 0.0513 | **2/10**（仅 outer 1–2） | 6 |
| local2x2+rel_ftb | 0.0705 | 0.0633 | 0.0603 | 0.0513 | 2/10 | 6 |
| local2x2+soft_CO α=0.1 skip-local | **0.0726** | 0.0626 | **0.0584** | 0.0491 | 2/10 | 7 |

解读：

1. **J/耦合压矮是真的**：局部 2×2 步长约为全系统 dense-CO 步的 **10×+**；FD 与组装 J 列一致 → 不是列奇异，是全局耦合把 CO 步压矮。  
2. **工程桥有 CO 正号**：出口 CO 0.033→**0.070**（对照平台仍 ~0.154）；+soft absorb 再微抬到 **0.073**、split 略好于纯 local2x2，仍差于 baseline split。  
3. **代价**：末态 split/rms/energy 略差；门禁（`split` 与 `merit` 均 ≤1.001×）只放行前两轮 Pre-Inner，之后提议被回滚。  
4. 平台 oneshot 降 split、全炉轨迹升 split → **后续 Newton/LS 把局部门禁收益冲掉或推向另一盆地**；勿仅凭 oneshot 默认开启。  
5. **勿默认开**；`frac` 默认 1.0（半步已证伪于平台）；rel_ftb 叠加无增量。

**建议下一刀** → 已做 §5.8。

### 5.8 outer≥3 拒因 = split LSQ 爆炸；post-gate local2；CO 缺口是库存（2026-08-07）

脚本 / 证据：

- `scripts/audit_upper_co_local2x2_reject.py` → `data/phase2_upper_co_local2x2_reject_o8.json`
- 逐层探针：`data/phase2_upper_co_local2x2_tier_probe_o6.json`
- 全局门禁：`data/phase2_upper_co_local2x2_global_gate_probe_o8.json`
- 复测 A/B：`scripts/audit_upper_co_local2x2_ab.py` → `data/phase2_upper_co_local2x2_ab_o8.json`

**拒因分类（local2x2 bundle 旧路径）**

| reject_class | 次数 | 含义 |
|---|---:|---|
| accepted | 2（其中 outer2 为虚报） | 真正 sticky 主要是 outer1 |
| **global_rollback** | **8** | local2 本地门禁过，全局回滚 |
| local_gate | 0 | 不是本地 split/merit 拒 |

逐层 RMS（bridge#3+）：

```text
start   rms≈0.063
split   rms≈13.5   ← phase-split LSQ 引爆
local2  rms≈13.5   （略降 split，本地仍 accept）
→ global gate 回滚整捆
```

关全局门禁对照：CO→**0.132**、上段 `|res_sum|→0.35`、`|res_split|→0.06`，但 **rms→16.6**（无效稳态）——证明「多吃几步 local2」能抬 CO，但不能裸关全局门禁。

**CO 0.03→0.15：库存闭合，不是缺产率**

| arm | exit CO | R_net CO | N_top CO | upper \|res_sum\| |
|-----|--------:|---------:|---------:|------------------:|
| baseline | 0.033 | 35.6 | 3.07 | 4.51 |
| local2 bundle | 0.070 | 34.9 | 6.92 | 3.58 |
| no_global_gate | 0.132 | 20.7 | 13.8 | 0.35 |
| 平台对照 | ~0.154 | — | — | — |

`R_net` 并未短缺；出口 CO 低是因为上段 `res_sum` 未吸入密相库存。闭合库存即抬 `N_top`。

**工程修复（默认：local2 开启时 post-gate=True）**

| 改动 | 位置 |
|------|------|
| `nr_outer_preinner_upper_co_local2x2_post_gate_thesis`（默认 True） | `reactor.py` |
| local2 移到全局门禁/fallback **之后**单独跑 | `run_outer_preinner_start_value_bridge` |
| local2 本地门禁加 `rms_ok` | `run_outer_preinner_upper_co_local2x2` |
| `split_only_fallback` 不再虚报 local2 accepted | 同上 bridge |

**全炉 o8 复测**

| arm | CO | rms | split | energy | local2 真接受 |
|-----|---:|----:|------:|-------:|-------------:|
| baseline | 0.033 | 0.050 | 0.050 | 0.042 | 0 |
| local2 **bundle**（post_gate=False） | 0.070 | 0.063 | 0.060 | 0.051 | 1/10 |
| **local2 post_gate** | **0.113** | 0.087 | **0.052** | 0.087 | 2/8 |
| post_gate+soft_CO | **0.114** | 0.086 | **0.051** | 0.086 | 2/8 |

结论：

1. outer≥3 拒因是 **bed1/底部 phase-split LSQ 在 local2 改写上段后引爆 rms**，不是 local2 本地门禁。  
2. post-gate 解耦后 CO 0.033→**0.113**（对照 ~0.154），split 接近 baseline；rms/energy 仍偏高——**勿默认开**整条 local2 链。  
3. 剩余 CO 缺口仍是库存 sticky 步不够 + 可能的 split 爆炸根因；**不是**再加反应速率。  
4. 单测：`tests/test_preinner_merit_driven.py` 15 passed。

**建议下一刀** → 已做 §5.9。

### 5.9 split LSQ 爆炸 = bed0 bubble/H2O；全局残差守卫（2026-08-07）

脚本 / 证据：

- `scripts/audit_preinner_split_lsq_detonation.py` → `data/phase2_preinner_split_lsq_detonation_o6.json`
- 守卫后复测：`data/phase2_split_residual_guard_postcheck_o6.json`
- A/B：`data/phase2_split_residual_guard_ab_o8.json`
- 单测：`tests/test_preinner_split_residual_guard.py`

**爆炸定位（12/12 次一致）**

| 项 | 值 |
|----|-----|
| 触发床 | **bed0** phase-split LSQ（`ok=True` 局部接受） |
| 全局 argmax | **cell0 / bubble / H2O**，`F_hat≈148` |
| top5 | bubble+dense **H2O ≈ +148**；**H2 ≈ −148**；bubble **O2 ≈ −74** |
| 旧守卫漏洞 | 只看 `gas_phase_split_rms` 不升；局部 split 可降而 **NR rms/max_abs 爆炸** |

**修复（默认开）**

`nr_outer_preinner_split_global_residual_guard_thesis=True`：每格 LSQ 后若

- `rms > running_rms × 1.25`，或  
- `max_abs > max(running_max×1.25, running_max+0.5)`  

则回滚该格 `N_d/N_b`（记 `global_residual_rollbacks`）。bed1 followup 同步守卫。

守卫后 split 出口：`max_abs` 保持 ~0.9（**0 次** post-split 爆炸）；bridge#3+ 上 `bed0_accepted=False`，`res_rb=1`。

**全炉 o8 A/B**

| arm | CO | rms | split | energy |
|-----|---:|----:|------:|-------:|
| baseline（守卫开，默认） | 0.033 | **0.047** | **0.047** | 0.043 |
| baseline_noguard | 0.033 | 0.050 | 0.050 | 0.042 |
| local2 post_gate + 守卫 | **0.113** | **0.063** | **0.051** | **0.054** |
| local2 post_gate 无守卫 | 0.113 | 0.087 | 0.052 | 0.087 |

解读：

1. 爆炸源是 **bed0 H2O/H2 相分配 LSQ**，不是上段 CO local2。  
2. 残差守卫默认开：无 local2 时 split/rms 略好；有 local2 时 **rms/energy 从 0.087→0.063/0.054**，CO 不变。  
3. 守卫是症状拦截；根因仍是 bed0 LSQ 目标（局部 split 向量）与全系统 NR 残差不一致——更深刀可改 bed0 目标/界限，非扩 maxbed9。

**建议下一刀** → 已做 §5.10。

### 5.10 bed0 凭空造 CO → phase-only majors（2026-08-07）

脚本 / 证据：

- `data/phase2_bed0_lsq_h2o_twist_o8.json`（LSQ 后、守卫回滚前状态）
- A/B：`data/phase2_bed0_phase_only_ab_o8.json`
- 复测：`data/phase2_bed0_phase_only_postcheck_o6.json`（blown=0）
- 单测：`tests/test_preinner_bed0_phase_only.py`

**状态轨迹（爆炸瞬间）**

| 物种 | Δ总量 | 相变化 |
|------|------:|--------|
| O2/H2O/N2 | **0**（冻总量） | H2O 仅 ~0.13 mol/s 相间挪动 |
| H2 | ~0 | 仍≈0 |
| **CO** | **+0.85**（0→0.85） | 密相凭空出现 |
| CO2 | −0.30 | 总量被改 |

结论：argmax 虽在 bubble/H2O，**触发物是产物自由 DOF**（`preserve_primary_inlet_totals` 下产物 `N_d/N_b` 独立 ≤80），不是 H2O 大挪动。造出的 CO 经 BC/反应耦合把 H2O/H2 残差打到 ~150。

**修复（默认开）**

`nr_outer_preinner_bed0_phase_only_majors_thesis=True`：bed0 Pre-Inner 与上层一样，**只优化各 major 相份额、冻物种总量**。  
残差守卫仍保留作安全带。

**全炉 A/B（o8；守卫开）**

| arm | CO | rms | split | energy | local2 接受 |
|-----|---:|----:|------:|-------:|------------:|
| baseline | 0.033 | 0.048 | 0.048 | 0.042 | 0 |
| local2 + **legacy** bed0 造产物 | 0.113 | 0.063 | 0.051 | 0.054 | 2/10 |
| **local2 + phase-only bed0** | **0.122** | **0.058** | **0.040** | **0.050** | **3/10** |
| local2+soft absorb | 0.095 | 0.057 | 0.053 | 0.057 | 2/10 |
| local2 o12 | 0.122 | 0.058 | 0.040 | 0.050 | 3/12（无增量） |

postcheck：`blown=0`；首轮 bed0 LSQ 可 `accepted=True` 且 `max_abs` 下降（0.22→0.19）。

解读：

1. bed0 造 CO 是 §5.9 爆炸的根因；phase-only 是根治，守卫是兜底。  
2. 与 local2 post-gate 叠加后：CO 0.033→**0.122**（对照 ~0.154），**split 反优于 baseline**（0.040）。  
3. soft absorb / 加长到 o12 **无进一步 CO 收益**——剩余缺口需别的杠杆（更多有效库存步、动力学/下游、或 Newton 耦合）。  
4. local2 与 bed0 phase-only 仍属工程桥组合；**local2 本身勿默认开**；bed0 phase-only 默认开（修 bug）。

**建议下一刀** → 已做 §5.11。

### 5.11 CO 缺口是上段库存；phase_share 吸收越过对照（2026-08-07）

脚本 / 证据：

- `scripts/audit_co_gap_to_platform.py` → `data/phase2_co_gap_to_platform_o8.json`
- 末态 oneshot：`data/phase2_upper_co_endstate_absorb_o8.json`、`data/phase2_upper_co_endstate_joint_absorb_o8.json`
- 全炉 A/B：`data/phase2_phase_share_absorb_gate_ab_o8.json`
- 单测：`tests/test_phase_share_absorb.py`

**逐 bed 预算（local2 vs baseline）**

| | baseline | local2 |
|--|--------:|-------:|
| exit CO | 0.033 | 0.122 |
| N_top CO | 3.06 | 12.66 |
| 对照 0.154 所需 N_top | 14.5 | 16.0 |
| ΔN vs 对照 | −11.4 | **−3.3** |
| R_net 下段/上段 | 18.5 / 17.3 | 15.9 / 13.6 |
| 上段正 res_sum | 18.1 | **5.5** |
| y_if 吸尽上段 res_sum | 0.225 | **0.175** |

解读：缺口**不是**缺产率（`R_net` 甚至略降）；是上段 `res_sum` 未吸入 holdup。local2 主要抬 bed6–9 的 `N`、压 `res_sum`。再吸完剩余 ~5.5 mol/s 即可超过 0.154。

**末态 absorb 模式**

| mode | frac=1 y | split | 门禁 |
|------|--------:|------:|:----:|
| dense_only | 0.145 | 0.040→**0.112** | ✗ |
| **phase_share** | **0.145** | 0.040→**0.039** | ✓ |
| equal_phase | 0.145 | 0.040→0.042 | ✗（frac=1） |

密相单吸抬 split；**按当前相份额分给两相**不抬 split。

**工程接线（均默认关，除既有 bed0 phase-only / 残差守卫）**

| 开关 | 含义 |
|------|------|
| `nr_outer_preinner_upper_dense_absorb_mode_thesis="phase_share"` | 相份额吸收 |
| `nr_outer_preinner_upper_dense_absorb_post_gate_thesis=True` | 门禁/local2 之后跑 |
| phase_share 局部门禁 | `split_ok ∧ rms_ok`（不因 energy merit 微升拒收） |

**全炉 o8**

| arm | CO | rms | split | energy |
|-----|---:|----:|------:|-------:|
| baseline | 0.033 | 0.048 | 0.048 | 0.042 |
| local2 only | 0.122 | 0.058 | 0.040 | 0.050 |
| **local2 + ps_abs f=1** | **0.163** | 0.060 | **0.039** | 0.060 |
| local2 + ps_abs f=0.6 | **0.159** | 0.060 | 0.039 | 0.060 |
| ps_abs + skip_local | 0.174 | **0.147** | 0.087 | 0.087 |

结论：

1. 对照 CO≈0.154 **已被 phase_share 库存闭合越过**（0.163）；勿再加动力学调参当主路径。  
2. skip_local 会毁掉 rms/split——禁止。  
3. 组合仍为工程桥：**local2 + phase_share absorb 勿默认开**；bed0 phase-only / split 残差守卫保持默认开。  
4. 新代价：rms/energy ~0.060（高于 baseline 0.048/0.042）——下一刀压残差/收敛，不是抬 CO。

**建议下一刀** → 已做 §5.12。

### 5.12 rms 税 = bed2 H2O/H2 反应翻转（2026-08-10）

脚本 / 证据：

- `scripts/audit_local2_ps_abs_rms_stack.py` → `data/phase2_local2_ps_abs_rms_stack_o20.json`
- bed2 探针：`data/phase2_bed2_h2o_hotspot_vs_local2_o8.json`
- min_bed A/B：`data/phase2_local2_minbed_ps_abs_ab_o12.json`

**收敛验收（o20）**

| arm | CO | rms | split | energy | converged | n_acc |
|-----|---:|----:|------:|-------:|:---------:|------:|
| baseline | 0.033 | 0.048 | 0.048 | 0.042 | no | 4 |
| local2 only | 0.122 | 0.058 | 0.040 | 0.050 | no | 5 |
| ps_stack f=1 | **0.163** | 0.059 | **0.039** | 0.055 | no | 5 |
| ps_stack f=0.3 | **0.156** | 0.059 | 0.039 | 0.055 | no | 5 |

- 出口 CO 在 f∈[0.3,1] 稳定于 **0.15–0.16**（对照 ~0.154）。  
- **未 fully converge**；best-so-far merit 卡在 outer≈2（~0.058）。  
- 软化 absorb **几乎不改变 rms** → rms 税不是 absorb 分数。

**残差主导（ps / local2 末态 top）**

| 排名 | 位置 | \|F_hat\| |
|-----:|------|--------:|
| 1–2 | **bed2 dense/bubble H2O** | ~0.43 |
| 3–4 | **bed2 dense/bubble H2** | ~0.40 |
| 5–6 | bed2 O2 | ~0.22 |

baseline 同点 H2O `|F|~0.01`；local2 后 **反应源翻转**：

| bed2 | baseline | local2 (min_bed=6) |
|------|--------:|-------------------:|
| R_H2O | −2.8 | **+33** |
| R_H2 | +3.9 | **−32** |
| R_O2 | −0.3 | **−18** |
| F_H2O | ~0 | **~+0.43**（密相+气泡同号 → 总量/源项问题，非 split） |

**min_bed 试验**

| arm | CO | rms | bed2 H2O 热点 |
|-----|---:|----:|:-------------:|
| local2 min6 + ps f1 | 0.163 | 0.059 | 有 |
| local2 **min7** + ps f1 | 0.138 | **0.085** | 无（但 rms 另坏） |
| local2 min7 + abs 自 bed6 | 0.067 | 0.052 | 无（CO 塌） |
| local2 min8 + ps f1 | 0.123 | 0.062 | 无 |
| absorb only | 0.093 | 0.054 | 无 |

解读：

1. rms/energy 上升的主因是 **bed2 H2O↔H2/O2 反应在 local2 含 bed6 时被轴向耦合拧爆**，不是上段 CO 方程本身。  
2. 避开 bed6（min7+）可消 bed2 热点，但目前 **保不住 CO≈0.15 且压 rms**——需另找床2闭合或更温和的 bed6 步。  
3. 工程桥现状：`local2(min6)+phase_share(f≈0.3–1)` 可交付对照 CO，附带 **~0.01 rms 税** 与未收敛；**勿默认开**。  
4. o12/o20 加长迭代 **无收敛收益**（n_acc 仍~5）。

**建议下一刀**

1. 解剖 local2(bed6) → bed2 的耦合路径（组成回流 / T / 速率分支）。  
2. 或设计「CO 进入 0.15 带宽后冻结 upper 桥、专修 bed2」的两阶段 outer。  
3. 文献侧仍优先 Wirsum。

### 5.13 bed6→bed2 翻转路径 + in-solve 冻结（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed6_local2_to_bed2_coupling.py` → `data/phase2_bed6_local2_to_bed2_coupling_o4.json`
- 轨迹：`data/phase2_bed2_flip_traj_twostage_o12.json`
- 冻结 A/B：`scripts/audit_freeze_upper_bridge_ab.py` → `data/phase2_freeze_upper_bridge_ab_o8.json`

**耦合假设**

| 假设 | 结果 | 证据 |
|------|------|------|
| H1：oneshot 改 bed6 CO + BC refresh → bed2 `R_H2O` 立即翻转 | **证伪** | oneshot 后 bed2 `R_H2O` 不变 |
| H2：翻转发生在 NR 轨迹（outer2→3 Newton） | **确认** | br#1–2 仍 ≈−2.9；br#3 进入 preinner 前已是 **+67** |

结论：local2(bed6) **不直接**拧爆 bed2；它把全炉状态推到一条 Newton 会放大 bed2 H2O/H2 源项的分支。

**两阶段陷阱**

- 两次完整 `reactor.solve()` 不能做「先抬 CO、再冻桥」：每次 solve 重跑 Vorab/`run_init_and_precalc`，CO 假塌缩。  
- 必须用同一次 solve 内冻结，或 `_solve_global_nr(..., skip_precalc=True)` 续算。

**已接线（默认关）**

```text
nr_outer_freeze_upper_co_bridge_after_y_co_thesis = False
nr_outer_freeze_upper_co_bridge_y_co_threshold_thesis = 0.14
nr_outer_freeze_upper_co_bridge_after_local2_accepts_thesis = 0   # >0 才启用
```

运行时旗标：`_upper_co_bridge_frozen_runtime`、`_upper_co_local2_accept_count_runtime`（每次 `solve` 清除）。  
冻结后关掉 **local2 + phase_share absorb**；split/zone preinner 仍跑。  
`nr_outer_history` 现含 `upper_co_bridge_frozen` / `y_co_top`。

**冻结 A/B（同构 + ps_f0.3，o8）**

| arm | CO | rms | split | bed2 R_H2O | freeze@ | local2 |
|-----|---:|----:|------:|-----------:|--------:|--------|
| ps_f0.3（对照） | **0.153** | 0.061 | 0.039 | **+65.5** | — | 3/10 |
| freeze y_CO≥0.14 | 0.141 | 0.059 | 0.040 | +65.5 | 6 | 3/5 |
| freeze y_CO≥0.10 | 0.130 | 0.059 | 0.040 | +65.4 | 4 | 3/3 |
| freeze after local2×1 | 0.068 | **0.049** | 0.037 | +49.8 | 1 | 1/1 |
| freeze after local2×2 | 0.089 | 0.058 | 0.044 | +65.2 | 2 | 2/2 |

解读：

1. **y_CO 阈值冻结太晚**：翻转已在 outer2→3 Newton 完成；之后停桥只略压 rms（~0.002），**修不了 bed2**，还把 CO 从 0.15 拉回 0.13–0.14。  
2. **更早冻（local2×1）**：CO 塌到 ~0.07，bed2 仍翻（+50），只是幅度略小——**单次 local2 接受已足以铺好翻转轨迹**。  
3. in-solve 冻结是可用的工程旋钮，但 **不是 bed2 闭合解**；**勿默认开**。  
4. 旧结论「skip_precalc 续算塌 CO」已修正：在 **不重跑 Vorab** 的前提下续算可保 CO≈0.15 且自愈 bed2（见下方）。

**窗口精定位（已做）**

脚本：`scripts/audit_bed2_flip_newton_window.py` → `data/phase2_bed2_flip_newton_window_o4.json`  
补充：`data/phase2_bed2_flip_inner2_rates.json`、`data/phase2_bed2_flip_inner2_cellvars.json`

| 事件 | bed2 R_H2O |
|------|-----------:|
| refresh#1 / preinner#1 / inner#1 / postinner#1 | ≈ −2.86 |
| preinner#2（local2 再接受，无 absorb） | ≈ −2.86 |
| **inner#2**（2 Newton iter；Check1 continuation，**无** refresh#2） | **−2.86 → +73~+90** |
| postinner#2 | 保持大正 |

要点：

1. 翻转不在 preinner/refresh/postinner，而在 **第 2 次 inner NR**。  
2. bed2 的 T / `N_H2O` / 水力几乎不动（ΔN_H2O~0.02），但 `R_gas` 自身爆炸。  
3. 这是 **速率面跳跃**（extent 限幅关闭下的动力学敏感）；local2 只是把 NR 推进该盆地。

**Δx 掩码：触发 DOF = bed2 微量 O₂**

脚本：`scripts/audit_bed2_flip_inner2_dx_mask.py` → `data/phase2_bed2_flip_inner2_dx_mask.json`  
ε 扫描：`data/phase2_bed2_flip_o2_epsilon_sweep.json`

| 掩码 | bed2 R_H2O | 解读 |
|------|-----------:|------|
| `none` / 非 bed2 | ≈ −2.86 | 不翻 |
| `cell2_bed_only` / `gas_only` / `sp_O2_allcells` | ≈ +90 | 翻 |
| **`bed2_O2`（‖Δx‖≈1e-16）** | ≈ +90 | **充分条件** |
| `gas_except_O2` / `all_except_cell2` | ≈ −2.86 | O₂ 必要 |

ε 扫描（仅改 bed2 `N_O2`，其余冻结在 inner#2 入口）：

- `N_O2=0` → 好分支；`1e-18` 开始拧号；**`1e-16` 全翻**；更大 ε 近似线性（`R_H2O ~ 8e8 · N_O2`）。  
- Newton 把 O₂ 从精确 0 抬到 ~1e-16 灰尘 → 裸 R1/R12 速率面爆炸。

**微量 O₂ 粘滞（已接线，默认关）**

```text
nr_snap_trace_o2_holdup_thesis = False
nr_snap_trace_o2_atol_mol_s_thesis = 1e-12
nr_snap_trace_o2_min_bed_index_thesis = 0   # 2=仅 bed2+
```

- residual 前 snap + `_clip_dx` 禁止贫氧格引入微量 O₂（有入口 O₂ 时不挡）。  
- A/B：`data/phase2_snap_trace_o2_ab_o8.json` — **从 solve 开头就开** → bed2 修好（R≈−2.5）但 **CO 塌到 ~0.05**。  
- 结论：抬 CO 的 NR 轨迹与 bed2 O₂ 病态盆地耦合；**不能**默认开 snap。

**自愈路径：skip_precalc 二次 outer（关键）**

`data/phase2_snap_o2_twostage_skip_precalc.json`、`data/phase2_cont_bridge_on_off_heal.json`

| 阶段 | CO | rms | bed2 R_H2O |
|------|---:|----:|-----------:|
| A：ps_f0.3 o8 | 0.153 | 0.061 | +65 |
| B：`_solve_global_nr(skip_precalc=True)` 再 o8（桥开或关均可） | **0.153–0.156** | **0.048–0.050** | **−2.6** |

- in-solve freeze（即便 o20）**不能**自愈 bed2（仍 +65）——卡在同一 outer/Check1 轨迹。  
- **重新进入 outer abgleich**（skip_precalc，不重跑 Vorab）即可在保住 CO 的同时把 bed2 拉回好分支。  
- 桥开/关对自愈几乎无差别；自愈主因是 **outer 回路复位**，不是 snap/冻桥。

**建议下一刀**（已由 §5.14 接手）

1. ~~在单次 `solve` 内做「CO 带宽后 outer 回路热重启」~~ → 见 §5.14。  
2. snap_o2 / freeze 保留为诊断旋钮，**勿默认开**。  
3. 文献侧仍优先 Wirsum。

### 5.14 单次 solve 内 outer 热重启（2026-08-10）

脚本 / 证据：

- `scripts/audit_outer_hot_restart_ab.py` → `data/phase2_outer_hot_restart_ab_o8.json`
- 补证：`data/phase2_outer_hot_restart_nested_fbproj_o8.json`

**已接线（默认关）**

```text
nr_outer_hot_restart_after_y_co_thesis = False
nr_outer_hot_restart_y_co_threshold_thesis = 0.14
nr_outer_hot_restart_bonus_inner_budget_thesis = 0   # >0 → cont_max_iter≈bonus/5
```

实现：首轮 outer 结束后若床顶 `y_CO≥阈值`，先做 **finalize 自由板 T 投影**，再嵌套一次  
`_solve_global_nr(skip_precalc=True)`（重建 inner 闭包；`skip_precalc` 路径不再套娃热重启）。  
结果字段：`nr_outer_hot_restart`；续算段 history 标 `outer_hot_restarted`。

**A/B（同构 + ps_f0.3）**

| arm | CO | rms | split | bed2 R_H2O | converged | 解读 |
|-----|---:|----:|------:|-----------:|:---------:|------|
| A：ps o8 | 0.153 | 0.061 | 0.039 | **+65.5** | False | 基线翻 |
| B：o8 + 热重启（嵌套 skip_precalc + 续算前 fb T 投影） | **0.156** | **0.050** | 0.042 | **−2.61** | False | **对齐 oracle** |
| C：ps o16 无热重启 | 0.156 | 0.059 | 0.039 | +65.6 | False | 加预算不够 |
| D：A 后显式 `_solve_global_nr(skip_precalc=True)` o8 | 0.156 | 0.050 | 0.042 | **−2.61** | False | oracle |

**证伪路径（勿再走）**

| 路径 | 结果 |
|------|------|
| CO 带宽后 **中途**复位 refresh/merit（同 outer 循环内） | bed2 仍 +65；CO 保住 |
| 同 `_solve_inner` 闭包二次 `_run_outer` | bed2 仍 +65 |
| 嵌套 skip_precalc **但无**续算前 freeboard T 投影 | bed2 仍 +65 |

**结论**

1. 热重启有效，且与显式两段 skip_precalc **数值对齐**（CO≈0.156、rms≈0.050、bed2 R_H2O≈−2.6）。  
2. **必要件**：续算前 `project_explicit_freeboard_temperatures_to_energy_closure`（与 `solve()` finalize 同款）。  
3. **勿默认开**——工程桥，非 Hamel 正文；需 COS 带宽触发，且加倍 wall≈180–250s/o8。  
4. ~~近零 O₂ 速率保护~~ → 见 §5.15。

### 5.15 微量 O₂ 速率门控（2026-08-10）

脚本 / 证据：

- `scripts/audit_trace_o2_rate_gate_ab.py` → `data/phase2_trace_o2_rate_gate_ab_o8.json`
- posthoc：`data/phase2_trace_o2_rate_gate_posthoc.json`
- 单测：`tests/test_trace_o2_rate_gate.py`

**动机**

snap holdup（§5.13）从 solve 开头启用 → bed2 修好但 **CO 塌到 ~0.05**。  
需要「不动 holdup、只挡速率面」且 **先抬 CO、再门控**。

**实现（默认关）**

```text
nr_trace_o2_rate_gate_thesis = False
nr_trace_o2_rate_gate_atol_mol_s_thesis = 1e-12
nr_trace_o2_rate_gate_min_bed_index_thesis = 2
nr_trace_o2_rate_gate_after_y_co_thesis = False
nr_trace_o2_rate_gate_y_co_threshold_thesis = 0.14
```

- `gate_trace_o2_for_rate_eval`：局地 `N_O2`/`C_O2` 低于阈值时，速率评价用的 C/y 中 O₂ 置 0。  
- **不看 `N_in`**：bed2 常有上游微量 O₂ 入口，按入口跳过会让门控失效（已踩坑）。  
- `after_y_co`：preinner 前检查床顶 `y_CO`，达阈值才武装（`arm_trace_o2_rate_gate_if_y_co`）。

**A/B（同构 + ps_f0.3，o8）**

| arm | CO | rms | split | bed2 R_H2O | wall |
|-----|---:|----:|------:|-----------:|-----:|
| ps_f0.3 | 0.153 | 0.061 | 0.039 | **+65.5** | ~104s |
| rate_gate **立即** bed2+ | 0.052 | 0.050 | 0.019 | −2.48 | ~101s |
| rate_gate **after y_CO≥0.14** bed2+ | **0.153** | 0.061 | 0.040 | **−2.69** | ~101s |
| hot_restart（§5.14） | **0.156** | **0.050** | 0.042 | −2.61 | ~213s |

**posthoc**：ps 末态仅武装门控、不改状态 → `R_H2O +65.5→−2.71`，CO 仍 0.153。

**结论**

1. **延迟 rate-gate** 可在单次 solve、**无热重启税**下同时保 CO≈0.15 与 bed2 好分支。  
2. 立即门控与 snap 一样塌 CO——抬 CO 的 NR 轨迹与 bed2 O₂ 病态盆地仍耦合。  
3. rms 仍停在 ~0.061（未吃到 hot_restart 的 0.050）；若要压 rms 仍可用热重启，或在武装后再加短续算。  
4. **勿默认开**（工程桥；改写速率评价，非 Hamel 正文）。

**建议下一刀**（已由 §5.16 接手）

1. ~~延迟门控武装后短续算 / 与热重启组合~~ → 见 §5.16。  
2. 或压 split/energy 平台冲 fully converge。  
3. 文献侧仍优先 Wirsum。

### 5.16 延迟 rate-gate 压 rms：加长 outer / 组合热重启（2026-08-10）

脚本 / 证据：`scripts/audit_rate_gate_rms_stack_ab.py` → `data/phase2_rate_gate_rms_stack_ab.json`

**A/B（同构 + ps_f0.3）**

| arm | CO | rms | split | energy | bed2 R_H2O | n_acc | wall |
|-----|---:|----:|------:|-------:|-----------:|------:|-----:|
| ps o8 | 0.153 | 0.061 | 0.039 | 0.061 | **+65.5** | 5 | ~98s |
| gate after y_CO o8 | 0.153 | 0.061 | 0.040 | 0.061 | −2.69 | 5 | ~109s |
| gate after y_CO **o16** | **0.156** | 0.057 | 0.040 | 0.057 | −2.62 | 5 | ~169s |
| hot_restart o8 | 0.156 | 0.050 | 0.042 | 0.049 | −2.61 | 7 | ~239s |
| **gate after y_CO + hot o8** | **0.156** | **0.049** | **0.040** | **0.049** | −2.61 | 5 | **~185s** |

**解读**

1. 仅加长 outer（o16）在延迟门控上略压 rms（0.061→0.057），**吃不完** hot_restart 的 ~0.050。  
2. **门控 + 热重启** 略优于单热重启：rms 0.049、split 0.040（单热 0.042），且 wall **更短**（~185s vs ~239s）——门控先把 bed2 拉回好分支，续算少做无用功。  
3. 均未 fully converge；n_acc 仍~5–7。  
4. 工程推荐栈（均默认关）：`local2+phase_share` → `rate_gate after y_CO` → 可选 `hot_restart`。  
5. **勿默认开**。

**建议下一刀**（已由 §5.17 接手）

1. ~~残差分解：gate+hot 末态 top 残差~~ → 见 §5.17。  
2. 或压 split/energy 平台。  
3. 文献侧仍优先 Wirsum。

### 5.17 gate+hot 末态 top 残差：能量/床顶 T 挡收敛（2026-08-10）

脚本 / 证据：`scripts/audit_gate_hot_topk_residuals.py` → `data/phase2_gate_hot_topk_residuals_o8.json`

**对照**

| arm | rms | gas_rms | energy_rms | max\|F̂\| | #1 残差 |
|-----|----:|--------:|-----------:|--------:|---------|
| ps o8 | 0.061 | **0.066** | 0.079 | 0.428 | **bed2 H2O/H2**（翻转税） |
| gate after y_CO o8 | 0.061 | 0.023 | 0.079 | 0.210 | bed9 **T** |
| hot_restart o8 | 0.050 | 0.025 | 0.069 | 0.179 | bed9 **T** |
| gate+hot o8 | **0.049** | **0.022** | **0.069** | 0.179 | bed9 **T** |

**gate+hot top6（解码后）**

| 排名 | 位置 | F̂ |
|-----:|------|--:|
| 1 | **bed9 T（能量）** | +0.179 |
| 2 | bed3 dense N2 | −0.126 |
| 3 | bed0 T | −0.081 |
| 4 | bed3 bubble H2O | −0.067 |
| 5–6 | bed5 N2 dense/bubble | ±0.064 |

**解读**

1. gate/hot 已消掉 bed2 H2O/H2/O2 热点；**气相 rms 从 0.066→0.022**。  
2. 未收敛主因切换为 **能量支**：`energy_rms≈0.069` 主导合成 rms；床顶 **bed9 T** 为最大单项。  
3. 次级：中段 **N2 相分配**（bed3/5）与 **bed0 T**——偏 split/底区能量，不再是 bed2 动力学。  
4. 固相残差可忽略（`solid_rms≈0.004`）。

**建议下一刀**（已由 §5.18 接手）

1. ~~专攻 bed9/bed0 能量闭合~~ → 见 §5.18。  
2. 或收紧 bed3–5 N2/H2O split（慎扩 LSQ）。  
3. 文献侧仍优先 Wirsum。

### 5.18 bed9 能量：非局部 T 根，而是床顶焓结构不平衡（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed9_energy_closure_ab.py` → `data/phase2_bed9_energy_closure_ab_o8.json`
- 项分解：`data/phase2_bed9_energy_terms.json`
- T 扫描：`data/phase2_bed9_T_sweep_bc.json`、`data/phase2_bed9_energy_residual_probe.json`

**基线**：gate+hot（§5.16）→ CO≈0.156、rms≈0.049、bed9 `|F̂_T|≈0.179`。

**A/B（在 gate+hot 上）**

| arm | CO | rms | bed9 F̂_T | 结论 |
|-----|---:|----:|----------:|------|
| gate+hot | 0.156 | 0.049 | +0.179 | 对照 |
| postinner fb T 投影 | 0.158 | **0.052** | +0.179 | **变差**；勿开 |
| solid/energy preinner max_bed=9 | 0.156 | 0.049 | +0.179 | 无效 |
| oneshot `preproject_bed_temperature…` bed9(/+0) | — | — | +0.179 | **拒收/无效** |

**物理分解（bed9，Eq.2.7）**

| 项 | 值 |
|----|---:|
| T / T_in_gas | 1150.82 / 1150.77 K |
| H_in | −6.21×10⁷ W |
| H_gas_out | −6.60×10⁶ W |
| H_sol_out | −5.89×10⁷ W |
| Q_loss | 0 |
| **res = H_in−H_out** | **+3.47×10⁶ W** |

- 邻床 bed7/8 同量级 |res| 仅 ~0.3 MW；**床顶 bed9 ~3.5 MW** 突出。  
- 局地 T 扫描（有/无 BC）：1140–1160 K 内 **无** 使 |E|→0 的根（最佳仍 >2.7 MW）。  
- 故不是「T 偏了 1 K」的局部投影问题，而是 **床顶固相/气相出流焓与 auf/ab/自由板耦合的结构不平衡**。

**结论**

1. 自由板 postinner T 投影、扩大 solid-energy preinner、单格 T→能量 LSQ **均不能**压 bed9 F̂_T。  
2. 未收敛能量支的根在 **bed-top ↔ freeboard 固相焓载体**（`_solid_energy_outlet_rates` / 扬析回流），不是再拧 NR 的 T。  
3. **勿**默认打开 `nr_outer_postinner_freeboard_energy_T_projection_thesis`（已有注释：会拉偏床层）。

**建议下一刀**（已由 §5.19 接手）

1. ~~审计 bed9↔fb0 固相质量/焓连续~~ → 见 §5.19。  
2. ~~床顶焓出口透传 A/B~~ → 见 §5.19。  
3. 并行可收 bed3 N2 split；文献侧仍优先 Wirsum。

### 5.19 bed9↔fb 质量连续完好；床顶焓透传压 F̂_T 但塌 CO（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed9_fb_solid_enthalpy.py` → `data/phase2_bed9_fb_solid_enthalpy_ab_o8.json`
- 开关：`nr_bed_top_solid_energy_passthrough_thesis`（默认 **False**）
- 实现：`Cell._solid_energy_outlet_rates` 认 `_solid_energy_passthrough`；`Reactor.solve` 挂床顶床
- 单测：`tests/test_bed_top_solid_energy_passthrough.py`

**界面质量连续（gate+hot 末态）**

| 差 | 值 |
|----|---:|
| bed8 `m_up` − bed9 `m_auf_in` | **0** |
| bed9 `m_up` − fb0 固体入口 | **0** |
| fb0 `m_down` − bed9 `m_ab_in` | **0** |

界面传质本身闭合；问题不在「床↔自由板断链」。

**床顶格内焓载体（同态）**

| 量 | 值 |
|----|---:|
| bed9 `m_e_out` − `m_in_all`（K·m 出口 vs 入口） | **+0.268 kg/s** |
| Eq.2.7 `res`（现用 E-out） | **+3.45 MW** |
| 同态改入口透传后的 `res_pass` | **−0.078 MW** |
| posthoc 强制透传：F̂_T | +0.179 → **−0.009** |

**A/B（gate+hot 栈）**

| arm | CO | rms | bed9 F̂_T | bed2 R(H₂O) | 结论 |
|-----|---:|----:|----------:|------------:|------|
| gate+hot | 0.156 | 0.049 | +0.179 | −2.61 | 对照 |
| +床顶焓透传 | **0.100** | 0.046 | **−0.015** | **+2.91** | F̂_T 修好，**CO 塌 + bed2 再翻** |

**结论**

1. §5.18 诊断成立：bed9 能量洞来自 **格内** `m_out(K·m) ≠ m_in` 的焓载体，不是界面漏质量。  
2. 类 `freeboard_closure` 的床顶透传在数学上能抹平该洞，但作为 **全程 NR 工程桥** 会把出口 CO 从 ~0.156 打到 ~0.10，并扰动 bed2 O₂ 分支——**勿默认开**。  
3. 下一刀应针对 **为何床顶 holdup 平衡下 K·m 出口与入口差 0.27 kg/s**（反应汇/扬析/尺寸迁移），或做 **门控/晚启透传**（仅 CO 已高后），而不是无条件透传。

**建议下一刀**（已由 §5.20 接手）

1. ~~分解 bed9 固相 `dm`~~ → 见 §5.20。  
2. ~~延迟床顶焓透传~~ → 见 §5.20。  
3. 并行：bed3 N2 / bed0 T；文献侧仍优先 Wirsum。

### 5.20 bed9 `dm`＝未收敛固相残差；延迟焓透传保 CO 且压 F̂_T（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed9_dm_deferred_passthrough.py` → `data/phase2_bed9_dm_deferred_passthrough_o8.json`
- 开关：`nr_bed_top_solid_energy_passthrough_after_y_co_thesis`（默认 **False**；须同时开 `nr_bed_top_solid_energy_passthrough_thesis`）
- 实现：`install_bed_top_solid_energy_passthrough_on_cells` / `arm_bed_top_solid_energy_passthrough_if_y_co`（`global_nr_solver.py`；preinner 武装）
- 单测：`tests/test_bed_top_solid_energy_passthrough.py`

**bed9 质量分解（gate+hot 末态）**

| 项 | char | ash | total |
|----|-----:|----:|------:|
| `m_in_all` | 3.800 | 1.090 | 4.890 |
| `m_out` (K·m) | 4.039 | 1.119 | 5.158 |
| **dm = out−in** | **+0.239** | **+0.029** | **+0.268** |
| `R_solid` | −0.049 | 0 | −0.049 |
| size_migration | 0 | 0 | 0 |
| **solid residual Σ** | **−0.288** | **−0.029** | **−0.317** |

恒等式：`dm = −(R+mig) − res_solid`（此处 `0.268 ≈ 0.049 + 0.317`）。  
故焓洞的质量侧主因是 **bed9 固相 NR 未收敛**（`|res|_max≈0.33 kg/s`，几乎全是 char），不是反应汇或尺寸迁移，也不是界面断链（§5.19 已证界面差=0）。

**A/B（gate+hot 栈）**

| arm | CO | rms | energy_rms | bed9 F̂_T | 结论 |
|-----|---:|----:|-----------:|----------:|------|
| gate+hot | 0.156 | 0.049 | 0.069 | +0.179 | 对照；#1=bed9 T |
| 全程透传（§5.19） | **0.100** | 0.046 | — | ≈0 | 塌 CO，勿用 |
| **透传 after y_CO≥0.14** | **0.156** | **0.0395** | **0.039** | **−0.009** | 保 CO；bed9 T 出榜 |

延迟武装后末态气/固与基线几乎同轨（T9、solid_res 同值），能量行因载体改写而跌落——属 **工程桥抹平能量表象**，**不**消除 char holdup 残差。top1 转为 bed3 N2（F̂≈−0.126），其次 bed0 T。

**结论**

1. 修 bed9 能量的正途是压 **固相 holdup 残差**（或弄清为何 bed-top `K·m` 持续 > 入口+R），不是无条件焓透传。  
2. `after_y_co` 透传可作 **可选工程桥**：保平台 CO、压 energy_rms / F̂_T；**勿默认开**。  
3. 全程透传仍禁止（塌 CO）。

**建议下一刀**（已由 §5.21 接手）

1. ~~专攻 bed9 char holdup~~ → 见 §5.21。  
2. 并行收 bed3 N2 split。  
3. 文献侧仍优先 Wirsum。

### 5.21 bed9 holdup：reconcile 跑了但 NR 回灌；延迟床顶 solid LSQ 有效（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed9_holdup_closure_ab.py` → `data/phase2_bed9_holdup_closure_ab_o8.json`
- 开关：`nr_outer_preinner_bed_top_solid_holdup_thesis` + `…_after_y_co_thesis`（均默认 **False**）
- 实现：`run_outer_preinner_start_value_bridge` post-gate 调 `preproject_bed_solid_holdup_to_transport_closure(top)`
- 单测：`tests/test_bed_top_solid_holdup_preinner.py`

**末态链（gate+hot）**

| bed | res_char | over_support_char | over_m*_char | 备注 |
|----:|---------:|------------------:|-------------:|------|
| 7 | −0.038 | +0.007 | +0.056 | 近闭合 |
| 8 | +0.002 | −0.054 | −0.002 | 闭合 |
| **9** | **−0.288** | **+0.430** | **+0.518** | 独突出 |

- `flux_mismatch≈0.42`（>0.15），`stream_imbalance≈0.05`（<0.15）→ 仅 auf 比失衡，stream 因 `K_ab` 尚可。  
- `_transport_inflow_support_holdup` **只加正 R**，不含炭消耗 → `m*`（含负 R）比 support 更低；即使贴 support 仍留 ~`R` 量级残差。  
- preinner `run_holdup_reconcile=True` 约 **10/20** outer 且 `holdup_reconciled` 记 True → **不是漏作用域**，而是 **inner NR 回灌** bed9 holdup。

**同态 posthoc（不重跑 NR）**

| 探针 | res_char | F̂_T | 结论 |
|------|---------:|-----:|------|
| align_reduce 上段 | −0.020 | −0.029 | 有效 |
| auf_upflow 只 bed9 | **−3.19** | **+2.41** | 有害 |
| m*=(in+R)/K bed9 | **−0.004** | −0.042 | 最佳同态 |
| solid LSQ bed9 | −0.037 | −0.015 | 有效 |

**全量 A/B（gate+hot 栈）**

| arm | CO | rms | F̂_T | res_char9 | 结论 |
|-----|---:|----:|----:|----------:|------|
| base | 0.156 | 0.049 | +0.179 | −0.288 | 对照 |
| postinner holdup reconcile | 0.158 | **0.055** | +0.202 | −0.317 | **变差** |
| solid_energy max_bed=9 | 0.156 | 0.049 | +0.179 | −0.288 | 无效（与 §5.18 一致） |
| bed-top solid LSQ（立即） | **0.132** | 0.063 | −0.019 | −0.043 | 修 holdup，**塌 CO** |
| **bed-top solid LSQ after y_CO** | **0.155** | **0.049** | **−0.017** | **−0.038** | 保 CO + 压 F̂_T/res |

**结论**

1. bed9 焓洞的质量根是 **过填 holdup**（相对 support / m*），不是界面断链。  
2. 现有 upper-bed reconcile **会跑**，但压不住床顶；扩 `solid_energy_max_bed=9`、postinner reconcile **无效或有害**。  
3. **勿**对床顶单独做 auf_upflow closure。  
4. 可选工程桥：`bed_top_solid_holdup` + `after_y_co`（与延迟焓透传同类；**勿默认开**）。  
5. 更深正途：support 计入负 R，或阻止 inner 回灌床顶 holdup（缩放/冻结），而非再扩 bottom-zone LSQ。

**建议下一刀**（已由 §5.22 接手）

1. ~~负 R support~~ → 见 §5.22。  
2. ~~LSQ + 延迟焓透传组合~~ → 见 §5.22。  
3. 并行 bed3 N2；文献侧仍优先 Wirsum。

### 5.22 负 R support 有害；延迟焓透传压 rms 最优，LSQ 修 holdup 但不叠加（2026-08-10）

脚本 / 证据：

- `scripts/audit_bed9_negR_support_combo_ab.py` → `data/phase2_bed9_negR_support_combo_ab_o8.json`
- 开关：`nr_holdup_support_include_negative_R_thesis`（默认 **False**）
- 实现：`_transport_inflow_support_holdup(..., include_negative_R=)`；`reconcile_upper_bed_solid_holdup_for_init` 读 cfg
- 单测：`tests/test_holdup_support_negative_R.py`

**A/B（gate+hot 栈，o8）**

| arm | CO | rms | energy_rms | F̂_T9 | res_char9 | 结论 |
|-----|---:|----:|-----------:|-----:|----------:|------|
| base | 0.156 | 0.049 | 0.069 | +0.179 | −0.288 | 对照；#1=bed9 T |
| **negR support** | **0.138** | 0.053 | 0.070 | +0.169 | −0.277 | **塌 CO + bed2 翻**；勿开 |
| bed-top LSQ after y_CO | 0.155 | 0.049 | 0.069 | **−0.017** | **−0.038** | 修 holdup；#1→**bed8 T** |
| **energy PT after y_CO** | 0.156 | **0.0395** | **0.039** | **−0.009** | −0.288 | **rms 最优**；holdup 未修 |
| LSQ + energy PT | 0.155 | 0.051 | 0.068 | −0.010 | −0.045 | **不叠加**；逊于单开 PT |
| negR + LSQ / full combo | 0.138 | 0.053 | 0.070 | +0.169 | −0.277 | 被 negR 毒化（CO 未达门控） |

**解读**

1. reconcile 贴「含负 R 的 m*」会在 **CO 抬升前**过度削 holdup → bed2 O₂/H₂O 分支翻、出口 CO 卡在 ~0.14 下，延迟门控无法武装。  
2. 延迟焓透传仍是压 `energy_rms` / 合成 rms 的最强工程桥（同 §5.20），但 **不碰** char 残差。  
3. 延迟床顶 solid LSQ 把洞从 bed9 挪到 **bed8 T**（F̂≈+0.157），合成 rms 几乎不变。  
4. LSQ+PT 同时开互相干扰，rms 回升到 ~0.051——**不要默认叠**。

**结论**

1. **`nr_holdup_support_include_negative_R_thesis` 勿开**（已实现仅供对照）。  
2. 当前可选最佳工程栈（均勿默认）：  
   - 要 rms：`gate+hot` + `energy_passthrough after_y_co` → rms≈0.0395  
   - 要 holdup 物理：`gate+hot` + `bed_top_solid_holdup after_y_co` → res_char≈−0.04，但能量热点上移 bed8  
3. 正途下一刀应跟 **bed8 能量/holdup**（LSQ 修 bed9 后的新 #1），或 bed3 N2；勿再拧负 R reconcile。

**建议下一刀**（已由 §5.23 接手）

1. ~~审计 bed8~~ → 见 §5.23。  
2. ~~并行 bed3 N2~~ → 见 §5.23（split_max3 弱效）。  
3. 文献侧仍优先 Wirsum。

### 5.23 bed8 是 LSQ 推上去的过填；双床顶 LSQ（n=2）同时压 holdup 与 rms（2026-08-11）

脚本 / 证据：

- `scripts/audit_bed8_after_top_lsq_ab.py` → `data/phase2_bed8_after_top_lsq_ab_o8.json`
- 开关：`nr_outer_preinner_bed_top_solid_holdup_n_beds_thesis`（默认 **1**；试 2=bed8+bed9）
- 单测：`tests/test_bed_top_solid_holdup_preinner.py`（含 n_beds=2）

**迁移证据（gate+hot → LSQ n=1）**

| 量 | base bed8 | LSQ n=1 bed8 | LSQ n=1 bed9 |
|----|----------:|-------------:|-------------:|
| F̂_T | −0.021 | **+0.157** | −0.017 |
| res_char | +0.002 | **−0.225** | −0.038 |
| dm (out−in) | −0.045 | **+0.209** | −0.008 |
| res_E | −0.32 MW | **+3.03 MW** | ~0 |

床顶 LSQ 修好 bed9 后，**同类过填/焓洞上移到 bed8**（非新物理机制）。  
同态：bed8 solid LSQ 有效（F8→+0.022）；单格 T→能量投影 **拒收**。

**全量 A/B**

| arm | CO | rms | energy_rms | F8 / F9 | res_char 8/9 | 结论 |
|-----|---:|----:|-----------:|--------:|-------------:|------|
| base | 0.156 | 0.049 | 0.069 | −0.02 / +0.18 | ~0 / −0.29 | 对照 |
| LSQ n=1 | 0.155 | 0.049 | 0.069 | **+0.16** / −0.02 | −0.23 / −0.04 | 热点→bed8 |
| **LSQ n=2** | **0.156** | **0.040** | **0.057** | **+0.01 / −0.02** | **−0.04 / −0.04** | **双闭合；#1=bed3 N2** |
| LSQ n=1 + energy PT | 0.155 | 0.051 | 0.068 | +0.15 / −0.01 | −0.22 / −0.05 | 不叠加 |
| LSQ n=1 + split_max3 | 0.156 | 0.052 | 0.069 | +0.16 / −0.02 | −0.23 / −0.03 | N2 \|F̂\| 0.126→0.102，整体变差 |

**结论**

1. bed8 热点是 **上游连锁过填**，不是独立 T 根。  
2. `n_beds=2`（床顶两格 solid LSQ，after y_CO）是目前 **holdup 物理 + rms** 兼顾最好的工程桥（rms≈0.040，接近单开延迟焓透传的 0.0395，且真正压 dm）。**勿默认开**。  
3. 扩 phase_split 到 bed3 **不能**作为主路径（rms↑）。  
4. 未收敛新 #1：**bed3 N2**（F̂≈−0.126）；其次 bed3/bed7 能量量级项。

**建议下一刀**（已由 §5.24 接手）

1. ~~专攻 bed3 N2~~ → 见 §5.24。  
2. ~~LSQ n=2 vs 焓透传对比~~ → 见 §5.24。  
3. 文献侧仍优先 Wirsum。

### 5.24 bed3 N2：惰性相分配 + 总量洞；targeted split 只挪热点（2026-08-11）

脚本 / 证据：

- `scripts/audit_bed3_n2_on_lsq_n2_ab.py` → `data/phase2_bed3_n2_on_lsq_n2_ab_o8.json`
- 开关：`nr_outer_preinner_targeted_phase_split_beds_thesis` + `…_after_y_co_thesis`（默认空/关）
- 单测：`tests/test_targeted_phase_split_preinner.py`

**LSQ n=2 末态 bed3 N2 分解**

| 量 | 值 | 解读 |
|----|---:|------|
| R_d / R_b | **0 / 0** | 确为惰性（同 §5.3） |
| y_d − y_b | **+0.011** | 小组成差 |
| N_ex | −3.80 mol/s | 相间交换放大 |
| res_split | **+10.1** | 相残差主导表象 |
| **res_sum** | **−14.9** | **总量也未闭合**（非纯 split） |

同态 bed3 phase-split LSQ：`dy→0`，`res_split→−0.64`，接受——但 **总量洞仍在**。

**A/B（均在 gate+hot 上）**

| arm | CO | rms | energy_rms | #1 | dy(N2) | 结论 |
|-----|---:|----:|-----------:|----|-------:|------|
| LSQ n=2 | 0.156 | 0.040 | 0.057 | bed3 N2 −0.126 | +0.011 | holdup 好；N2 #1 |
| LSQ n=2 + max_bed=3 | 0.157 | **0.049** | 0.057 | bed3 T | ≈0 | **rms 变差** |
| LSQ n=2 + target bed3 | 0.156 | 0.040 | 0.057 | **bed3 T +0.107** | ≈0 | N2↓，热点→能量 |
| LSQ n=2 + target 3+5 | 0.156 | 0.040 | 0.057 | bed3 T | ≈0 | 与单 bed3 同 |
| **energy PT after y_CO** | 0.156 | **0.0395** | **0.039** | bed3 N2 | +0.011 | **rms 最优** |
| LSQ n=2 + energy PT | 0.154 | 0.040 | 0.041 | bed3 N2 | +0.011 | 略伤 CO；N2 未消 |

**推荐工程栈（均勿默认开）**

| 目标 | 栈 | rms | 代价 |
|------|----|----:|------|
| **稳定默认** | gate+hot + LSQ n=2（H2O 透传**关**） | ≈0.046 | CO≈0.159；§5.29 |
| 冲 rms（opt-in） | LSQ n=2 + H2O 透传 + 延迟焓透传；`PYTHONHASHSEED=0` + OMP=1 | ≈0.039 | CO≈0.151；有双吸引子风险 |
| 压 N2 表象 | ~~targeted split / 链式 sum~~ | — | **已由 §5.27 init N2 透传取代** |

**结论**

1. 扩 `max_bed=3` 仍有害（复现 §5.23）。  
2. targeted bed3 split **能**抹平 Δy，但把 #1 换成 bed3 能量，且 **res_sum≈−15** 说明根因含 **N2 总量/交换闭合**，不是再拧一次相 LSQ。  
3. LSQ n=2 与延迟焓透传 **择一**：叠开略掉 CO 且不消 N2。  
4. 下一正途：bed3 **总量**（`res_sum`）+ 能量，或接受双策略分治。

**建议下一刀**（已由 §5.25 接手）

1. ~~分解 bed3 N2 `res_sum` / 总量吸收~~ → 见 §5.25。  
2. 文献侧仍优先 Wirsum。

### 5.25 bed3 N2 是 holdup 过填；双向放出闭合本格但热点上移 bed4（2026-08-11）

脚本 / 证据：

- `scripts/audit_bed3_n2_sum_closure_ab.py` → `data/phase2_bed3_n2_sum_closure_ab_o8.json`
- 开关：`nr_outer_preinner_targeted_species_sum_*`（默认空；`bidirectional` 默认 True）
- `oneshot_bed_dense_species_sum_absorb(..., bidirectional=True)` 支持负 `res_sum` 放出
- 单测：`tests/test_species_sum_bidirectional.py`

**预算恒等式（LSQ n=2 末态，R=0）**

| 量 | 值 |
|----|---:|
| N_in | 46.31 mol/s |
| N_hold | 61.26 mol/s |
| **gap = hold−in** | **+14.95** |
| res_sum | **−14.95** |
| N_ex | −3.80（在 sum 中对消，不贡献总量洞） |

→ 根因是 **N2 holdup 过填**，不是交换项；旧 absorb（只吸正残差）对本格 **无效**。

**同态 posthoc**

| 探针 | res_sum | dy | 结论 |
|------|--------:|---:|------|
| absorb 单向 | −14.95 | +0.011 | 无效 |
| **双向放出** | **0** | +0.012 | 总量闭合；Δy 仍在 |
| 放出 + split | **0** | ≈0 | 总量+相都闭合 |
| bed3 T→能量 | −14.95 | +0.011 | 拒收 |

**全量 A/B（gate+hot + LSQ n=2）**

| arm | CO | rms | #1 | bed3 res_sum | 结论 |
|-----|---:|----:|----|-------------:|------|
| LSQ n=2 | 0.156 | 0.0402 | bed3 N2 −0.126 | −14.95 | 对照 |
| N2 sum bed3 | 0.156 | 0.0407 | **bed4 N2 −0.164** | **0** | 洞上移 bed4 |
| **sum + split bed3** | 0.156 | **0.0394** | bed4 N2 −0.141 | **0** | rms≈焓透传；热点→bed4 |
| sum f=0.5 | 0.156 | 0.0407 | bed4 N2 −0.164 | ≈0 | 同 cascade |

**结论**

1. bed3 N2 总量洞 = **过填 holdup**（与床顶 char 过填同构：局部闭合→上游/下游邻格冒头）。  
2. 双向 `phase_share` 放出是正确局部手术；**勿**期望单格结束全炉 N2。  
3. `sum+split` 把合成 rms 压到 **0.0394**（与延迟焓透传 0.0395 持平），但仍以 bed4 N2/T 为 #1。  
4. bed3 T 投影无效；下一刀若追 N2，应对 **bed3–5 链式 sum**（慎，易伤 CO），或接受分治栈。

**建议下一刀**（已由 §5.26 接手）

1. ~~链式 N2 sum / cascade 通盘~~ → 见 §5.26。  
2. 文献侧仍优先 Wirsum。

### 5.26 N2 通盘：不是逐格独立 bug，是 bed3 过填抬高整段平台（2026-08-11）

脚本 / 证据：

- `scripts/audit_n2_axial_cascade.py` → `data/phase2_n2_axial_cascade_o8.json`
- 基线栈：gate+hot + 床顶 solid LSQ `n_beds=2`

**全床 N2 剖面（LSQ n=2 末态）**

| bed | N_in | N_hold | gap(=hold−in) | dy | 备注 |
|----:|-----:|-------:|--------------:|---:|------|
| 0 | 42.3 | 42.3 | ≈0 | −0.013 | 闭合 |
| 1 | 51.7 | 45.6 | **−6.11** | −0.007 | **欠填** |
| 2 | 47.7 | 46.3 | −1.36 | −0.003 | 略欠 |
| **3** | **46.3** | **61.3** | **+14.95** | +0.011 | **唯一过填尖峰** |
| 4–9 | 61.3 | 61.3 | ≈0 | +0.005~0.013 | 局部闭合，但绝对值被 bed3 抬高 |

- `gap_sum ≈ +7.48 = (+14.95) + (−6.11) + (−1.36)`：系统净多 ~7.5 mol/s N2。  
- bed4–9 的 `gap≈0` **具有误导性**：它们只是与错误平台 `N=61.3` 自洽（`N_out≡N_hold` 的 CSTR 链），不是「没问题」。  
- 轴向：`bed2.N_out = bed3.N_in = 46.3`，但 `bed3.N_hold = 61.3` → **断裂发生在 bed3 格内**，不是 bed2→3 界面传质断链。

**同态自下而上逐格放出（证实 cascade）**

| 步骤 | peak gap 位置 | 说明 |
|------|--------------|------|
| base | bed3 +14.95 | |
| bleed bed1（补欠填） | bed3 | 峰值不动 |
| bleed bed2 | bed3（gap 降为 +7.48） | 欠填被吃掉，净过剩显性化 |
| bleed bed3 | **bed4 +7.48** | 出流下降 → 下一格入口降、holdup 未动 |
| bleed bed4…8 | peak 逐格上移 | 典型推扫 |
| bleed bed9 | gap_sum→0 | 链末端才收口 |

机制（CSTR）：`N_out,k = N_hold,k = N_in,k+1`。压 bed k 的过填 ⇒ 降 `N_in,k+1` ⇒ 若 `N_hold,k+1` 仍停在旧平台，则 gap 原样上移。

**全量链式 sum A/B（after y_CO）**

| arm | CO | rms | N2 峰值 | 结论 |
|-----|---:|----:|---------|------|
| LSQ n=2 | **0.156** | 0.040 | bed3 +14.95 | 对照 |
| sum beds 3–5 | 0.156 | 0.041 | **bed6 +14.95** | 尖峰上移，不根治 |
| sum beds 3–9 | **0.203** | 0.040 | 上段 gap=0；bed1 欠填仍在 | **塌平台 CO** |
| sum beds 1–9 | **0.178** | 0.0395 | 全 gap≈0 | CO 仍偏高 |

**与 §5.3 的关系**

- §5.3（早期平台）：上段 N2 以 **Δy×交换放大** 为主、`res_sum≈0`。  
- 当前 gate+hot+LSQ n=2：中段另出现 **bed3 总量过填 → 抬高 bed4–9 平台**；Δy 仍在（+0.005~0.013），但 #1 残差已被总量洞主导。  
- 两者叠加：惰性相分配问题仍在；总量 cascade 是后加的中段结构缺陷。

**结论**

1. 用户观察正确：**问题会沿 cell 上移**——因为是一条抬高的 N2 流量平台，不是 10 个独立 bug。  
2. 手术点应在 **bed3 为何相对入口多囤 ~15 mol/s**（及 bed1 为何欠填 ~6），而不是对 bed4–9 逐个放出。  
3. 全链 sum 能抹平 gap，但 **破坏出口 CO（0.156→0.18~0.20）**——勿默认开。  
4. 推荐仍分治：rms→延迟焓透传；holdup→LSQ n=2；N2 深追先查 **bed2/3 过渡 + bed1 欠填**，勿盲扩链式放出。

**建议下一刀**（已由 §5.27 接手）

1. ~~专攻 bed3 N2 过填来源~~ → 见 §5.27（init / axial `max(inflow,current)`）。  
2. 文献侧仍优先 Wirsum。

### 5.27 bed3 N2 过填根因在 Vorab axial；N2 按 inflow 闭合（2026-08-11）

证据 / 产物：

- `data/phase2_n2_inert_inflow_fix_o8.json`
- 改动：`src/solvers/vorabrechnung/vorab_cell_mapping.py`
  - axial：`N2 = inflow + N_zu + N_rez`（禁止 `max(inflow, current)`）
  - stream mapping：无 O2 补气床格做 **N2 透传**；H2O / O2 路径不动
- 单测：`test_stream_mapping_clears_inflated_n2_holdup_on_non_injection_beds`

**根因（init，不是 NR）**

1. 仅 `run_phase2_init`（未跑 NR）时 bed3 N2 gap 已达 **+23**（NR 末态约 +15 更轻）。  
2. 分级进气只改写有 `N_zu(O2)` 的 bed0–2；bed3+ 保留抬高的 transport x0。  
3. `apply_vorabrechnung_bed_axial_primary_profile_from_budget` 对 N2 用 `max(inflow, current)` → **保留过填并向上传播**，形成 §5.26 平台。  
4. bed1 欠填：分级锚点只写本地 `N_zu`，旧 axial 取 `max(inflow, zu)≈inflow`，丢掉本格 `N_zu` 累积。

**试错（勿再开）**

| 尝试 | 结果 |
|------|------|
| 分级 O2 也改 `inflow+zu` | 中段 O2 holdup 暴涨，CO→~0.10 |
| H2O 也改严格 `inflow+zu` | 抹热解水分，CO→~0.11、rms↑ |
| 链式 N2 sum（§5.26） | 尖峰上移或塌 CO |

**修复后（gate+hot，o8）**

| arm | CO | rms | init/end bed3 N2 gap | 上段 N_hold |
|-----|---:|----:|---------------------:|------------:|
| LSQ n=2 after y_CO | **0.159** | 0.046 | ≈0 / ≈0 | **53.8**（旧 61.3） |
| 延迟焓透传 after y_CO | 0.206 | 0.066 | ≈0 / ≈0 | 53.8 |

- N2 轴向 gap 全床 ≈0；#1 残差不再是 bed3 N2 总量洞。  
- CO 仍贴平台（对照 ~0.154）；rms 略高于旧「带 N2 过填」的 0.040——那是错误平台上的表象最优。  
- 延迟焓透传在 N2 闭合后 **不再** 优于 LSQ n=2（CO/rms 双差）；holdup 物理栈优先 LSQ n=2。

**结论**

1. bed3 N2 洞是 **Vorab axial / stream 惰性闭合 bug**，不是联立 NR「把残差上推」。  
2. 只修 N2；O2 仍 zu-only + axial blend；H2O 仍 `max(inflow, current)`。  
3. 勿默认开链式 N2 sum；勿把 H2O/O2 一并改成严格透传。

**建议下一刀**（已由 §5.28 接手）

1. ~~在 N2 已闭合的 init 上重扫 energy / split 热点~~ → 见 §5.28。  
2. 文献侧仍优先 Wirsum。

### 5.28 后 N2 残差重扫；post-staged H2O 透传 + LSQ2（±焓透传）压 rms（2026-08-11）

证据：

- `data/phase2_post_n2_rms_rescan_ab.json`（N2-only 基线上的 o8 重扫）
- `data/phase2_post_n2_h2o_passthrough_o8.json` / `data/phase2_post_n2_h2o_stack_confirm_o8.json`
- 改动：`vorab_cell_mapping` — 无 O2 补气格 **H2O+N2** 透传；axial H2O 在无 O2 格用 inflow+zu+rez

**N2-only 基线重扫（gate+hot）**

| arm | CO | rms | 结论 |
|-----|---:|----:|------|
| LSQ n=2 | 0.161 | 0.046 | 对照；E≈0.065 主导；#1=bed3/0 能量 |
| LSQ + 延迟焓透传 | 0.161 | 0.046 | 几乎无效 |
| LSQ + phase-split max=3 | 0.205 | 0.054 | **有害**（复现） |
| LSQ o12 | 0.161 | 0.046 | 与 o8 相同 → **已平台** |

- posthoc `T→energy` 对 bed0/bed3 **拒收**；能量洞不是单格 LSQ 能拧的。  
- init 已见 bed3 **H2O gap≈+8**（与已修的 N2 过填同类：`max(inflow,current)` 保留 x0）。

**post-staged H2O 透传 A/B（确认跑）**

| arm | CO | rms | 结论 |
|-----|---:|----:|------|
| LSQ n=2 | **0.151** | **0.0392** | 保平台；能量仍 #1（bed0） |
| 仅延迟焓透传 | 0.187 | 0.044 | 无 LSQ → 塌 CO |
| **LSQ n=2 + 延迟焓透传** | **0.151** | **0.0387** | **当前 rms 最优** |

注意：曾出现 LSQ-only 坏吸引子（CO≈0.22 / 全床 T≈1150）；同栈复跑落到 0.151。工程上 **LSQ2+焓透传** 更稳。

**结论**

1. N2 闭合后 rms 瓶颈是 **能量（bed0/bed3）+ bed3 H2O 过填**，不是再加 outer。  
2. bed3+ H2O 与 N2 同源（axial `max(inflow,current)`）；**只对无 O2 补气格**清 H2O，分级段仍 `max(inflow,current)` 保热解水。  
3. 推荐栈（均勿默认开）：`gate+hot` + **床顶 solid LSQ n=2 after y_CO** + **延迟焓透传 after y_CO** → CO≈0.151 / rms≈0.039。  
4. 勿扩 phase-split max=3；勿单开焓透传。

**建议下一刀**（已由 §5.29 接手）

1. ~~bed0 能量 / 吸引子稳定性~~ → 见 §5.29。  
2. 文献侧仍优先 Wirsum。

### 5.29 bed0 能量钉死；双吸引子与线程/hash；H2O 改 opt-in（2026-08-11）

证据：

- `data/phase2_bed0_energy_post_h2o_diag.json`
- `data/phase2_bed0_energy_seed_epre_ab_o8.json` / `phase2_bed0_energy_seed_epre_omp1_ab_o8.json`
- `data/phase2_attractor_stability_o8.json` / `phase2_n2only_lsq2_stability_o8.json`
- 改动：
  - `vorab_a_tier_post_staged_h2o_passthrough_thesis`（默认 **False**）
  - `_lu_harness` 强制 BLAS=1 线程 + `PYTHONHASHSEED` 警告

**bed0 能量**

| 观察 | 值 |
|------|-----|
| seed | R12 开 → bottom=1150、cold_anchor=False；fence≈[1000,1500] |
| 末态 T0 | ≈1150（贴 seed，几乎不动） |
| F̂_T(T) 扫描（坏吸引子末态） | 零点约 **1330 K**；1150 处 F̂_T&lt;0（要升温） |
| posthoc T→energy LSQ | bed0/bed3 **拒收** |

seed/Epre A/B（即使钉线程）均未优于对照，且常有害：

| arm | CO | rms | 结论 |
|-----|---:|----:|------|
| base LSQ2+Ept | 0.15 或 0.21 | — | 双吸引子 |
| Epre bed0 | 0.118 | 0.044 | **塌 CO** |
| seed1250 | 0.174 | 0.063 | 变差 |
| seed1300 | 0.226 | 0.079 | T→1450，更差 |

**双吸引子（同栈）**

| 吸引子 | CO | rms | 特征 |
|--------|---:|----:|------|
| 好 | **0.151** | **0.039** | rate-gate@iter3；T 仍≈1150 |
| 坏 | **0.214** | **0.058** | 全床 T≈1150；N2/能量上榜 |

触发因素（已证实相关）：

1. **OpenBLAS/OMP 多线程**（未钉 1 时易进坏支）  
2. **`PYTHONHASHSEED` 未固定**（须启动前 `PYTHONHASHSEED=0`；进程内设置无效）  
3. post-staged **H2O 透传**放大分叉面（N2-only 更稳）

**推荐栈（更新）**

| 目标 | 栈 | CO / rms | 备注 |
|------|----|---------|------|
| **稳定默认** | gate+hot + LSQ n=2；**H2O 透传关** | ≈0.159 / 0.046 | N2 已闭合 |
| 冲 rms（opt-in） | 上 + H2O 透传 + 延迟焓透传；`OMP=1` + `PYTHONHASHSEED=0` | ≈0.151 / 0.039 | 仍可能偶发坏支，复跑 |

**结论**

1. bed0 能量洞是 **NR 不走温度自由度**（seed 钉死），不是 fence 无根；抬 seed / 开 Epre **无效或有害**。  
2. §5.28 的 rms 最优依赖 H2O 透传，但引入双吸引子 → H2O 改 **thesis opt-in**。  
3. 审计必须：`PYTHONHASHSEED=0 OMP_NUM_THREADS=1 python3 scripts/...`（harness 已强制 BLAS=1）。

**建议下一刀**（已由 §5.30 接手）

1. ~~查 NR 为何不更新 bed0 T~~ → 见 §5.30。  
2. 文献侧仍优先 Wirsum。

### 5.30 bed0 T：Newton 大步降温被线搜索回退（2026-08-12）

证据：

- `data/phase2_bed0_T_nr_trace_o8.json`
- `data/phase2_bed0_T_newton_ab_o8.json`
- `data/phase2_bed0_T_ls_reject_probe.json`
- 栈：gate+hot + LSQ n=2，H2O 透传关，`PYTHONHASHSEED=0`

**Outer 轨迹（N2-only LSQ2）**

| 阶段 | dT_outer | 说明 |
|------|---------:|------|
| pre-hot o1–o2 | 0.3 → 46 | 有温变 |
| pre-hot o3–o10 | **0** | 全炉 T 冻结；energy_rms 仍缓降（preinner/气相） |
| hot-restart o1 | 161 | 他床主导；bed0 净变≈0 |
| 末态 | T0≈1150 | 贴 seed；F̂_T≈−0.13 |

**Newton / clip（`_clip_dx` 探针）**

- 多次提出 **raw dT0 ≈ −347 K**（降温），且 raw≈clip → **不是** split-T-cap=150 或 fence 上界在掐步。  
- 与 dT_budget=350 同量级（几乎吃满预算）。  
- 全求解净 ΔT0≈−0.3 K。

**线搜索探针**

| step | clip dT0 | 试探 T | 结果 |
|-----:|---------:|--------|------|
| 1 | −346 | 1150→**1000**（下界）→回 1150 | merit≈E≈0.105 主导 |
| 4–8 | −222~−350 | 同样摸到 ~1000 再回退 | 接受后 T 仍≈1150 |

结论：NR **有**温度自由度；方向是**大幅降温**，LS 试满步后 **merit 不降而回退**，表象为「钉在 seed」。

**A/B（勿当默认）**

| arm | CO | rms | ΔT0 | 结论 |
|-----|---:|----:|----:|------|
| base | 0.161 | 0.046 | −0.4 | 对照 |
| split-T-cap→600 | 0.161 | 0.046 | −0.4 | 无效（大步本就未按 150 限） |
| dT_budget=0 | 0.154 | 0.048 | −0.3 | 略贴平台 CO；T 仍钉 |
| **关 energy merit** | **0.116** | 0.030 | −0.1 | **塌 CO**；勿关 |
| Tcap600+无 budget | 0.154 | 0.048 | −0.3 | 同无 budget |

**与单变量 T 扫描的矛盾**

- 固定组成下 F̂_T(T) 在 ~1330 K 过零（§5.29 坏吸引子扫描）暗示「要升温」。  
- 耦合 Newton 却稳定给出 **dT0&lt;0**。  
- 根因属 **气–能耦合下的能量假根 / 线性化方向**，不是简单 seed/fence（呼应早期 `phase2_bed0_energy_false_root_*`）。

**结论**

1. bed0 T 钉死 = **LS 拒绝（回退）Newton 降温大步**，不是 DOF 被关掉。  
2. 抬 seed / 松 T-cap / 关 energy merit **都不能**根治；关 energy merit 会塌出口 CO。  
3. 下一正途：能量方程与气相耦合的可解性（Wirsum 阻尼 / 假根），或 composition-aware 的 T 预投影；勿再拧 seed。

**建议下一刀**

1. 分解「满步 T→1000」时 merit 中 gas/split/energy 谁上升（LS 回退的直接原因）。  
2. 文献侧仍优先 Wirsum。

### 5.31 LS 回退分解：gas/split 升、energy 降；freeze-gas 可接受但塌 CO（2026-08-12）

证据：

- `data/phase2_bed0_T_ls_merit_decomp_o8.json` — 大 `|dT0|` 线搜索试探的 merit 分量  
- `data/phase2_bed0_T_ls_projection_ab.json` — 首个 `|dT0|>100` 步的投影 A/B  
- `data/phase2_bed0_T_split_constrained_ab_o8.json` — 打开既有 `zero_gas` 门控的 o8  
- 栈：gate+hot + LSQ n=2，H2O 关，`PYTHONHASHSEED=0`

**满步回退时谁在升（o8，n=40 大降温试探）**

| 指标 | reject 上均值 Δ | 主导计数 |
|------|----------------:|----------|
| energy | **−0.01**（略降） | 0 |
| split | +1.6（早期局部 +30） | 2 |
| gas / rms | **+7e4 / +5e4** | **37** |

- 基点常已是 energy 主导 merit（E≈merit≈0.105）；试探时 **E 继续略降**，但 **gas/split 爆炸** → `merit_trial` 升 → 回退。  
- §5.30「merit≈E」描述的是**基点地板**；回退的直接原因是 **trial 上 gas/split 升**，不是 energy 升。

**投影 A/B（首个 dT0≈−347 截获）**

| arm | LS 可接受？ | 现象 |
|-----|:-----------:|------|
| full | 否 | split/gas ~1e9 |
| gas×0.1 / 压 bed0 T | 否 | 仍炸 |
| zero_T（只动气） | 否 | 气方向本身有毒 |
| **zero_gas**（冻气相 holdup） | **是** | Δmerit≈ΔE≈−0.005；T→1000 |

**工程门控 A/B（`nr_line_search_split_constrained_step_thesis`）**

| arm | CO | rms | T0 | 说明 |
|-----|---:|----:|---:|------|
| 关（默认） | **0.159** | 0.046 | 1149 | 平台；大降温步被 LS 挡住 |
| 开 | **0.094** | 0.068 | **1000** | 触发 33 次 freeze-gas；接受冷步 → **塌 CO** |

结论：

1. LS 回退在保护出口 CO：挡住的是「E 略降 + 毒气体 Newton 步」的耦合假根方向。  
2. 现成 `zero_gas`（energy 主导时冻气相）能使冷步 **被接受**，但会把 bed0 送进 **T 下界坏支**，勿默认开。  
3. 毒在 **气相 Newton 分量**（`zero_T` 仍炸）；单纯拧 T-cap/seed/关 energy merit 仍非正途。  
4. 下一正途仍是 **改线性化方向**（Wirsum 阻尼细节 / composition-aware 升温预投影），使能量下降与升温支一致，而不是让 freeze-gas 放行降温。

**建议下一刀**

1. 文献：Wirsum (1998) 阻尼 Newton 对 T/组成步的处理。  
2. 代码探针（opt-in）：仅当 `|dT0|` 大且 trial 升温方向时才允许能量步；或对 bed0 做 composition-frozen 单变量 T 校正（与耦合 Newton 对照），验证能否离开 1150 而不塌 CO。

### 5.32 composition-frozen T：预投影无效；F_E(T) 路径依赖（2026-08-12）

证据：

- `data/phase2_bed0_T_comp_frozen_preproj_ab_o8.json`
- `data/phase2_bed0_T_heat_override_ab_o8.json`
- `data/phase2_bed0_T_frozen_vs_coupled_dT_o4.json`
- `data/phase2_bed0_T_local_slope_o8.json`
- `data/phase2_bed0_T_clean_scan_o8.json`
- 栈同 §5.31（gate+hot + LSQ n=2，H2O 关，`PYTHONHASHSEED=0`）

**A. 粗扫 / 预投影 / 强制升温（不可当根）**

| 臂 | CO | rms | T0 末 | 备注 |
|----|---:|----:|------:|------|
| base | 0.159 | 0.046 | 1149 | 平台 |
| post_init / mid_o2 `preproject_bed_temperature…` | ≈0.16 | ≈0.047 | 1150 | **拒收**（界内无法过 improvement gate） |
| force T0=1500 后 o8 | **0.146** | 0.047 | **1150** | NR 拉回 seed；CO 变差 |

粗扫表象：`|F̂_T0|` 在 1500 K 最小、界内无变号 →「要升温」。但见 C。

**B. 升温改写 Newton dT0（opt-in 探针）**

- 条件：耦合 `dT0<-80` 且冻组成 secant `dT_f>+5` → 改为升温 capped 80 K。  
- 结果：见 18 次大降温，**0 次改写**；CO/rms/T0 与 base 相同。  
- 原因（o4 截获）：`F̂≈−0.145` 但局部 `∂F̂/∂T≈−1.7e−4` → 冻组成 secant 也给出 **dT_f≈−850 K**（与耦合同向降温）。

**C. 关键：F_E(T) 非路径无关**

每点 `unpack(x0)→设 T→清 cache→residual` 仍不稳定：

| T | 同点两次 F̂ |
|---:|-------------|
| 1000 | −0.108 / −0.162 |
| 1150 | −0.160 / −0.133 |
| 1500 | −0.008 / **+0.031** |

同 T 三次 reps 有时稳、扫完再评 F0 可从 −0.133 跳到 +0.11。  
结论：**仅 pack 的 holdup+T 不足以冻结能量残差路径**（`apply_bc` / `global_residual` 改写未打包隐藏态）。此前「单变量零点 ~1330/1500」与升温预投影叙事**不可作为闭环证据**。

**结论**

1. 现成 bed0 能量 T 预投影 / 强制抬 T **不能**离开 1150 且保 CO。  
2. 「冻组成要升温、耦合要降温」在**可信局部线性化**上并不成立：两者在 NR 大步点都指向降温；粗扫升温像是路径污染假象。  
3. bed0 T 钉死的正途仍须：**(a)** 可复现的能量残差探针（深快照），或 **(b)** Wirsum 全文阻尼细节；勿再基于脏 T-scan 加工程桥。

**建议下一刀**

1. 深快照（或 `cell` 深拷贝）后重做干净 F_E(T)；若仍无界内根，则能量债不在 T DOF。  
2. 文献：Wirsum (1998)。  
3. 并行可回到 rms 主项（bed2 等），bed0 T 在无干净探针前降优先级。

### 5.33 深快照 F_E(T)：无界内根；「要升温」已证伪（2026-08-12）

证据：

- `data/phase2_bed0_T_hidden_state_diff_o4.json` — pack 恢复后的隐藏态 diff  
- `data/phase2_bed0_T_deep_snap_scan_o8.json` — 深快照路径无关扫描  
- 栈同 §5.31

**污染源（§5.32 复现）**

热扰动 T=1500 后再 `unpack(x0)`：F̂ 从 −0.133 → **+0.111**。  
残留差异主要是 **bed0 `K_solid_auf`**、**bed1 `m_solid_auf_in`**（及 `_work_res`）。  
深快照恢复全部 ndarray/标量后，F̂ **精确回到** −0.133（含再 apply）。

**干净扫描（o8 末态，每点深恢复）**

| 项 | 结果 |
|----|------|
| 同 T 三次 + 热污染后再评 | **repro_ok** |
| [1000,1500] 变号根 | **无** |
| min \|F̂\| | **1000 K**（−0.108）；1500 K 最差（−0.197） |
| 局部 ∂F̂/∂T（±25 K） | ≈ **−1.7×10⁻⁴** |
| 冻组成 secant dT | ≈ **−770 K**（降温） |
| 耦合 Newton dT0 | ≈ −350 K（同向） |

F̂(T) 近似单调：升温 → F̂ 更负、\|F̂\| 更大。无 BC refresh 的冻入口扫描同样无界内根。

**结论**

1. 旧叙事「单变量要升温 / 零点 ~1330–1500」是 **扬析系数未恢复** 造成的假象，作废。  
2. 在可信冻组成下，能量残差与耦合 Newton **一致要降温**；LS 挡住冷步是因为 **气相爆炸**（§5.31），不是 T 方向反了。  
3. 界内 **不存在** 仅靠 bed0 T 消掉的能量根；能量债不在「把 T 拧到某种子」上。  
4. freeze-gas 放行降温能降 E、但塌 CO——气–能必须联合，不能拆成纯 T 预投影。

**建议下一刀**

1. 气–能联合：小步冷 T + 同步气相修正（非全 `zero_gas`），或文献 Wirsum 阻尼。  
2. 否则降优先级 bed0 T，转 rms/split 主项（bed2 等）。  
3. 凡 bed0 T 扫描必须深快照（至少 `K_solid_auf` / 下游 `m_solid_auf_in`）。

### 5.34 气–能联合软桥：小步冷仍塌 CO；缩气不够（2026-08-12）

证据：`data/phase2_bed0_T_gas_energy_joint_ab_o8.json`（栈同 §5.31）

**截获首步（dT0≈−347）投影**

| arm | LS ok? | 备注 |
|-----|:------:|------|
| full / gas×0.01 / gas×0.001 | 否 | split/gas 仍炸（×0.001 仍 Δmerit~7e4） |
| zero_bed0_gas + dT20 | 否 | 他床气相分量足够毒 |
| **zero_gas**（含 dT cap 10–80） | **是** | ΔE≈Δmerit≈−0.0056；split 微升可接受 |
| T_only_dT20 | 是 | 同 zero_gas_dT20 |

**o8 策略 A/B**（`|dT0|>80` 时改写 clip）

| mode | CO | rms | T0 | energy | n_pol |
|------|---:|----:|---:|-------:|------:|
| base | **0.159** | **0.046** | 1149 | 0.065 | 0 |
| zero_gas_dT20 | **0.114** | 0.068 | **1000** | 0.096 | 40 |
| zero_gas_dT40 | **0.114** | 0.059 | **1000** | 0.039 | 22 |
| gas×0.01_dT20 | 0.159 | 0.051 | 1146 | 0.072 | 7 |
| zero_bed0_gas_dT20 | 0.160 | 0.046 | 1150 | 0.065 | 3 |

**结论**

1. LS 只接受 **全冻气相**；任何残留气相 Newton 分量（含 ×0.001、仅 bed0）仍足以炸 merit。  
2. 全冻气 + 小步冷 T **不能**停在中间温度：多 outer 仍滑到 **T 下界**，CO≈0.11（与 §5.31 全步 freeze-gas 同病）。  
3. 缩气 / 仅冻 bed0 气：保 CO，但 **几乎不松 bed0 T**，能量不改善。  
4. bed0 T 钉死在当前栈上 **无安全的 clip 级工程桥**；暂停该类旋钮。

**建议下一刀**

1. **转主线**：压 rms/split（bed2 H2O/H2 等历史主项），接受 bed0 T≈1150 为平台特征直至有 Wirsum/新线性化。  
2. 文献：Wirsum (1998) 阻尼细节（若可获取）。  
3. 勿再开 soft `zero_gas`* 进默认。

### 5.35 平台残差重定位：bed2 已愈；压 rms 仍靠 H2O 透传（2026-08-12）

证据：

- `data/phase2_rms_landscape_current_stack_o8.json`
- `data/phase2_rms_press_ab_o8.json`（含 absorb / bed3 / H2O 透传）
- 栈：gate+hot + LSQ n=2，H2O 透传默认关，`PYTHONHASHSEED=0`

**景观（相对 §5.12 bed2 热点时代）**

| 项 | 现状 |
|----|------|
| bed2 `R_H2O` | **≈−3.06**（好分支；gate+hot 已愈） |
| 全局 | rms≈**energy_rms≈0.046**；split≈0.031；gas_rms≈0.062 |
| topk | **bed0 能量**、**bed3 能量**、bed3 dense **H2O**、bed9/8 dense **H2(/H2O)** |
| bed3 H2O | N≈21.5，`res_sum≈−11`（过填） |

**压 rms A/B**

| arm | CO | rms | bed3 N_H2O | 结论 |
|-----|---:|----:|-----------:|------|
| base | **0.159** | 0.046 | 21.5 | 稳定默认 |
| **H2O 透传** | 0.151 | **0.039** | 16.2 | 唯一明显压 rms；仍双吸引子风险（§5.29） |
| absorb +H2 / +H2O | 0.143 | 0.049 | ~21 | 伤 CO，rms 变差 |
| post_init bed3 H2O 双向 absorb | 0.160 | 0.046 | 21.5 | init 速率未就绪可致 **N 爆炸**（~1e13），NR 拉回原状 |
| mid_o2 bed3 H2O 双向 absorb | 0.155 | 0.050 | 20.6 | 暂降 N→9 后 **NR 再过填**；rms 更差 |
| targeted phase-split bed3 | 0.159 | 0.046 | 21.5 | 无效 |

**结论**

1. 历史 bed2 H2O/H2 翻转 **不再是** 当前平台 rms 主因。  
2. rms 地板由 **能量（bed0/3）** 与 **bed3 H2O 过填** 共撑；上段扩 absorb 有害。  
3. 已知唯一压 rms 旋钮仍是 **`vorab_a_tier_post_staged_h2o_passthrough_thesis`**（opt-in）；oneshot/targeted-split **不能**替代透传。  
4. bed3 过填会在 NR 中 **重建**——根在轴向/Vorab 持水政策，不在末态 LSQ。

**建议下一刀**

1. 解剖 bed3 H2O：为何 mid absorb 后 NR 把 N 从 ~9 拉回 ~21（与 N2 `max(inflow,current)` 同类？）。  
2. 或接受 rms≈0.046 平台，专攻能量主项（有深快照约束，§5.33）。  
3. H2O 透传保持 opt-in + `PYTHONHASHSEED=0`；勿默认、勿扩 upper absorb。

### 5.36 bed3 H2O：「NR 再膨胀」= 重跑 Vorab 假象（2026-08-12）

证据：

- `data/phase2_bed3_h2o_reinflate_trace_o8.json`（逐 outer；含误用二次 `solve`）
- `data/phase2_bed3_h2o_reinflate_skip_ab.json`（`skip_precalc` 对照 + init gap）
- 代码：`vorab_cell_mapping.py` H2O 默认 `max(inflow, current)`；透传旗标同 §5.28

**Init 已过填（根因）**

| 旗标 | bed3 N_H2O | inflow | gap N−in |
|------|----------:|-------:|---------:|
| H2O 透传 **关**（默认） | 20.6 | 12.5 | **+8.0** |
| H2O 透传 **开** | 14.8 | 12.5 | +2.3 |

单次 o8（无 absorb）：init 20.6 → 末态 21.5（gap 仍 ≈+7.4）。过填不是末态 NR「造」出来的，是 **x0 政策**。

**「再膨胀」复现拆解**

| 步骤 | N_H2O |
|------|------:|
| o2 后 | 20.6 |
| bidirectional absorb | **9.0** |
| 再调用 `reactor.solve()`（默认重跑 init/Vorab） | **立刻 20.6** |
| absorb 后 `_solve_global_nr(skip_precalc=True)` ×o6 | **仍 9.0** |

结论：§5.35 所见「NR 把 N 拉回 21」主要是 **续算未设 `skip_precalc`，Vorab axial 再次 `max(inflow,current)` 写回**；真 NR 在 skip 路径下 **不回填**。

**absorb 不能当透传替代**

skip 保住 N≈9 后：CO≈0.155、**rms≈0.061**（差于 base 0.046）。  
N&lt;inflow 留下正 `res_sum`，不是透传目标态（透传 init≈14.8）。

**结论**

1. bed3 H2O 过填 = **Vorab `max(inflow,current)`**（与已修的 N2 同类；H2O 仍 opt-in）。  
2. 诊断续算必须 `skip_precalc=True`，否则会误判「NR 再膨胀」。  
3. oneshot absorb ≠ 透传；勿当压 rms 默认桥。  
4. 默认栈可接受 bed3 gap≈+8 与 rms≈0.046；要压 gap/rms 仍只有透传（双吸引子风险，§5.29）。

**建议下一刀**

1. 接受当前平台，转其他主项；或文献 Wirsum。  
2. 若再动 H2O init：只考虑与 N2 同构的非 O2 格 `inflow+zu+rez`，并强制 `PYTHONHASHSEED=0` 双吸引子回归——**勿默认开**。  
3. 审计脚本续算一律 `skip_precalc`。

### 5.37 能量地板：bed0/bed3 Hin−Hout 对撞；焓桥无效（2026-08-12）

证据：

- `data/phase2_energy_term_decomp_platform_o8.json`
- `data/phase2_energy_press_ab_o8.json`
- 栈同 §5.35（gate+hot + LSQ n=2，H2O 关）

**项分解（o8 末态）**

| bed | E [W] | F̂_T | Hin−Hout [W] | Q [W] | dm(E−K) |
|----:|------:|-----:|-------------:|------:|--------:|
| 0 | **−2.50e6** | −0.133 | −2.42e6 | +8.5e4 | **0** |
| 3 | **+2.49e6** | +0.133 | +2.57e6 | +8.5e4 | **0** |
| 2 | −1.12e6 | −0.060 | −1.04e6 | +8.5e4 | 0 |
| 9 | −5.5e5 | −0.030 | −4.7e5 | +8.5e4 | 0 |

- 全局 `energy_rms≈0.046` ≈ `|F̂_T0|` 量级；**Q ≪ |E|**。  
- 固相能量出口与输运出口 **已重合**（`dm_E−K=0`）——与早期 bed9 E/K 分裂不同。  
- bed0 与 bed3 的 E **近等幅反号**（同 T≈1150）：轴向焓级联不平衡，不是单床 T 种子问题（呼应 §5.33 无界内 T 根）。

**A/B**

| arm | CO | rms | F̂_T0 | 结论 |
|-----|---:|----:|------:|------|
| base | 0.159 | 0.0461 | −0.133 | 对照 |
| 延迟床顶焓透传 | 0.159 | 0.0457 | −0.135 | **几乎无增益**（E/K 已齐） |
| solid-energy preinner maxbed9 | 0.159 | 0.0461 | −0.133 | 无效 |

**结论**

1. 当前能量地板 = **bed0/bed3 焓流闭合差**，不是热损、不是床顶 E/K、不是可拧的单变量 T。  
2. §5.20–5.22 的延迟焓透传在本栈 **失去杠杆**（前提 dm≠0 已消失）。  
3. 与 §5.33–5.34 一致：勿再期望 clip/焓桥消掉 ±2.5 MW。

**建议下一刀**

1. **冻结**能量/焓工程桥试验；平台验收以 CO≈0.15–0.16 / rms≈0.046 为基线。  
2. 文献：Wirsum 阻尼是否改变气–能联合步。  
3. 可选旁路：H2O 透传 opt-in 压 rms（已知代价，§5.29）。

### 5.38 能量分项：进料主导 bed0；对消对；dm×hs 假杠杆（2026-08-12）

证据：`data/phase2_energy_hin_parts_axial_o8.json`（栈同 §5.35，CO=0.159 / rms=0.046）。

**H_in 分项**

| bed | 主导入流焓 | 量级 |
|----:|------------|------|
| 0 | **仅** `gas_zu` + `solid_zu`（进料） | −2.2 / −8.2 MW |
| 3 | `solid_auf` + `solid_ab`（级联） | −14.3 / −10.4 MW |

**轴向固相焓连续（H_out,s vs 下格 H_auf）**

- bed0→1、1→2：**gap = 0**（质量与焓均对齐）。  
- bed2 以上 gap 变大：因上段 **auf+ab 双向**，总 `H_out,s ≠ H_auf,next`，属指标假象，非断裂。

**dm×hs 反事实（把固相质量亏缺焓塞进气出流）**

| | E [MW] |
|--|-------:|
| bed0 基线 | −2.50 |
| +dm·h_s,out | **+3.25**（翻号） |
| +dm·h_s,in | **+1.88**（仍翻号） |

气化产物已在 `H_out,gas`（生成焓）中部分体现；朴素 dm×hs **双重计算**，不是可修的记账洞。

**全局**

| 量 | 值 |
|----|---:|
| bed0 + bed3 | **−0.016 MW**（对消） |
| ΣE_bed | **−1.88 MW** |
| ΣQ_bed | +0.85 MW |
| ΣE_fb | 0 |

净地板在 bed2（−1.12）、bed5（−0.71）、bed9（−0.55）等，不是单点 T/焓桥能拧掉的。

**结论**

1. §5.37「对撞」= **局部对消对**；全局仍余 ≈−1.9 MW。  
2. 非固相 E/K、非 bed0→1 连续断裂、非 dm→气焓漏记。  
3. 与 §5.33–5.37 一致：**冻结能量工程桥**。

**建议下一刀**

1. 平台基线验收（CO≈0.15–0.16 / rms≈0.046）。  
2. 文献 Wirsum（§6.4：**全文已入库**；笔记 `docs/wirsum_1997_nr_solver_notes.md`）。  
3. 勿再开焓透传 / soft-freeze / dm×hs 类桥。

### 5.39 平台验收：全炉逐 cell 剖面与异常（2026-08-12）

证据：

- `data/phase2_platform_cell_profile_o8.json`
- Canvas：工作区 `canvases/platform-cell-profile.canvas.tsx`
- 栈同 §5.35；CO=0.159 / rms=0.046 / converged=False

**出口 vs LU Fig.7.4**

| 物种 | 平台湿基 | LU 湿基 | 判定 |
|------|--------:|--------:|------|
| CO | 0.159 | 0.13 | 可作 NR 基线 |
| H₂ | **≈0** | 0.12 | **临界异常** |
| CH₄ | **≈0** | 0.028 | **临界异常** |
| CO₂ | 0.040 | 0.11 | 高严重度偏低 |
| H₂O | 0.229 | 0.175 | 偏高 |

**轴向要点**

- bed0：y_O₂=0.093，CO=0，dm=+0.50（氧化+气化，预期内）；E=−2.5 MW  
- bed1：O₂ 耗尽  
- bed3–5：CO 谷 ~0.020；bed3 E=+2.5 MW（与 bed0 对消）  
- bed6–9：CO 回升→0.159；H₂ 始终≈0；上段 topk= dense H₂/H₂O  
- T≈1150 K 全床平直（无文献吸热凹）

**结论**

1. **数值平台可冻结**（CO/rms）。  
2. **物理合成气不可验收**——H₂/CH₄ 塌缩是下一主线。  
3. 能量地板已知，不阻塞 H₂ 调查。

**建议下一刀**

1. H₂ 预算分解（产：R2/挥发分/R8；耗：R12/R6/交换）；先 bed0–2 与 bed8–9。  
2. 对照限幅开/关或旧验证快照（曾有湿 H₂≈0.087）做差分。  
3. 勿回到能量工程桥。

### 5.40 方法纠偏 + 逐格 H₂：产得出来、存不进去（2026-08-12）

**方法（用户确认，写入 §1）**

1. 全炉逐 cell 查 **T + 多组分** 分布，标异常格。  
2. 对异常格做源汇/守恒追因。  
3. **禁止**以抬单一出口组分（如 CO）或压 rms 作为唯一成功标准——§5.39 已证明：CO/rms「平台」可与 H₂/CH₄ 物理崩溃共存。

证据：`data/phase2_platform_h2_cell_budget_o8.json`（同平台栈 o8）。

**逐格 H₂（摘要）**

| bed | y_H₂ | N_H₂ [mol/s] | R_H₂,net | 主产项 |
|----:|-----:|-------------:|---------:|--------|
| 0 | 0 | **0** | **+2.58** | R2 1.28 + VM 1.30 |
| 1 | ~1e-5 | 7.7e-4 | +4.28 | R2 2.87 + VM 1.30 |
| 3–5 | ~3e-5 | ~2.5e-3 | +3.8–4.1 | R2 ~3.7–4.0 |
| 9 | ~9e-6 | 8.0e-4 | **+6.79** | R2 5.82 + R8 0.97 |

全炉化学净产 H₂ ≈ **+44 mol/s**（R2≈38，VM≈3.5，R8≈2.5；R12/R3≈0）。  
出口 N_H₂ ≈ **8×10⁻⁴ mol/s**。

**根因定位（本刀）**

1. **不是**「R2/挥发分没开」——底部与全床都在产 H₂。  
2. **不是**「R12 氧化吃光」——extent R12≈0，且 `R_gas` 与化学计量净产一致为正。  
3. **是**：稳态应有 `N ≈ in + R`（约数 mol/s），实际 N 钉在 ~10⁻³；H₂ 方程残差巨大（与 topk dense H₂ 一致），**NR 吸引子不积累 H₂ 库存**。  
4. 首异常格 **bed0**：R=+2.58 且 N=0（产源已在、库存为零）。

**建议下一刀（仍按逐格）**

1. bed0→bed1：H₂ 是否在 NR 未知量 / 是否被 clip、fastox、syngas-guard、线搜索 merit 打回。  
2. 单格：强制抬 `N_H2` 后残差是否下降（证伪「方程不要 H₂」）。  
3. 同步扫 CH₄（同类症状）。

### 5.41 bed0 O₂–H₂ 刚性；上段可强制闭合（2026-08-12）

证据：`data/phase2_platform_h2_nh2_lock_probe_o8.json`（平台 o8 末态探针）。

**DOF**

- H₂ **在** NR 未知量内（与 CO/CO₂/H₂O/CH₄/O₂/N₂ + 热解气同列）。排除「没进 Jacobian」。

**bed0（首异常格）**

| 探针 | 结果 |
|------|------|
| 基线 | N=0，R≈+2.58，F≈+2.58（dense F̂≈0.032） |
| 强制 N←in+R，无 fastox | R → **−10¹³** 量级（R12 与残留 O₂ 刚性爆炸） |
| 强制后再跑 LS fastox | H₂ **2.58→0**；O₂ 仍剩 ~4.8 mol/s |

结论：bed0 在 O₂ 未耗尽时 **不能** 稳态积 H₂；fastox（R12 优先）与裸动力学一致地把 H₂ 打回零。N=0 时 R2/VM 仍往 `R_gas` 记产氢 → 残差结构性存在，直到 O₂ 区消失或计量闭合改写。

**bed1–9（O₂≈0）**

| 探针 | 结果 |
|------|------|
| 强制 N←in+R，再 fastox | H₂ **保持** 4–7 mol/s（fastox 不烧） |
| 出口 y_H₂ | **≈0.067**（对照 LU 湿基 0.12；基线 ≈0） |

结论：上段 **不是** fastox/动力学禁区；方程允许有限 H₂ 库存。平台吸引子 N_H₂~10⁻³ = **NR 路径未走到该盆地**（步进/初值/merit），不是「产不出」。

**分层**

1. **bed0**：O₂–H₂ 分区问题（氧化带）。  
2. **bed1+**：求解器未积累（可强制闭合证伪）。

**建议下一刀**

1. 上段：轴向 H₂ seed / 产物 holdup 对齐 `in+R`（仅 O₂≈0 格），或短 NR 续算看是否锁住 y_H₂≈0.07。  
2. bed0：勿幻想在破氧前堆 H₂；查 O₂ 消耗与产氢区是否错位（过宽氧化带）。  
3. A/B 慎关全床 fastox（bed0 可能更炸）；若关，仅作对照且限 iter。

### 5.42 Wirsum ‖F‖₂ 接受落地（关 max-merit 默认）（2026-08-20）

代码：`reactor.py` / `_lu_harness` LADDER p0–p2 / `_line_search_merit` docstring；getattr 回退 `False`。  
证据：`data/phase2_wirsum_rms_merit_ab_o4_o8.json`（`PYTHONHASHSEED=0`，栈 ISOMORPH+PS+GATE_AFTER+HOT+LSQ_n2，H2O 透传关）。

| arm | CO | H₂ | rms | split | energy | n_acc |
|-----|-----|-----|-----|-------|--------|-------|
| **wirsum_rms o4**（新默认） | **0.116** | ~4e-6 | **0.025** | 0.033 | 0.106 | 1 |
| legacy max-merit o4 | 0.161 | ~8e-6 | 0.065 | 0.031 | 0.065 | 3 |
| **wirsum_rms o8**（新默认） | **0.116** | ~4e-6 | **0.020** | 0.029 | 0.076 | 2 |
| 旧平台 o8（max-merit） | 0.159 | ~9e-6 | 0.046 | — | — | — |

**结论**

1. 同构接受准则已默认：仅缩放残差 RMS/‖F̂‖₂ 下降。  
2. 相对旧平台：**rms 明显下降，出口 CO 回落到 ~0.116**（与历史「关 energy merit 塌 CO」一致，§5.x）；H₂ 仍近零。  
3. fastox / syngas guard **默认已关**（工程桥）。需要旧 CO 平台时显式 opt-in L4 + max-merit。  
4. 数值「平台冻」指标须按新默认重标定；物理验收仍未过。

### 5.43 Hamel+Wirsum 核心栈（去掉平台桥）（2026-08-20）

命名栈：`scripts/_lu_harness.py` `HAMEL_WIRSUM_CORE`；`reactor.py` 同步关 LS fastox / syngas guard。  
证据：`data/phase2_hamel_wirsum_core_o4_o8.json`（`PYTHONHASHSEED=0`）。

**保留（同构/Startwert）**

- 床层 `extent_limiters_enabled_thesis=False`
- Vorab 快氧化 x₀（`vorab_transport_x0_fast_oxidation_closure_thesis=True`）
- Wirsum ‖F̂‖₂/RMS 接受；outer-fixed hydro
- Phase2 builder：`exact_hamel`、bed1 preinner（config_phase2，Startwert 基础设施）

**关闭（L4 / LS 工程桥）**

- LS fastox、syngas collapse guard
- local2 / phase_share absorb、rate-gate、hot restart、床顶 solid LSQ、H2O 透传

| arm | CO | H₂ | rms | n_acc | converged |
|-----|-----|-----|-----|-------|-----------|
| **core o4** | **0.012** | 0 | 0.020 | **0** | False |
| **core o8** | **0.012** | 0 | 0.020 | **0** | False |
| core + LS fastox o4 | **0.012** | 0 | 0.020 | **0** | False |
| Wirsum 叠在 L4 平台 o8 | 0.116 | ~4e-6 | 0.020 | 2 | False |
| 旧 max-merit 平台 o8 | 0.159 | ~9e-6 | 0.046 | — | False |

**结论**

1. 纯 Hamel 方程 + Wirsum 接受 **没有 fully converge**；o4≡o8，Inner **一步都没接受**。  
2. 只加 LS fastox（仍无 L4）**同样 n_acc=0**（`data/phase2_hamel_wirsum_core_ls_fastox_o4_o8.json`）。先前能走步的是 PS/GATE/HOT/LSQ，不是 Wirsum。  
3. 残差未天文爆炸（Vorab x₀ 仍在），但出口合成气塌到 init 附近。  
4. 禁止把 L4 整包写回「同构默认」。下一刀在核心栈上查 **为何所有 Newton 试探被拒**（逐 cell / clip_history），而不是加回 absorb。

### 5.44 核心栈 n_acc=0：J 好、步太刚、Wirsum 回溯不够（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_ls_reject_o1.json` 及 T-only / deep-LS / uniform-α 对照。

**现象（o1，CORE）**

- `line_search_failures=1`，12 次回溯全拒，`n_acc=0`
- 未剪裁 Newton：`||J dx + F̂|| / ||F̂|| ≈ 6×10⁻⁴`（线性步合格，Wirsum Eq.2.21）
- 现有 holdup FTB 剪裁后线性缺陷升到 **0.78**（方向被毁），但**不是**拒步主因：去掉 FTB 仍全拒
- 最小试探 λ=1/2048 时 trial RMS 仍 ~10³–10⁵；24 次回溯到 λ≈1.2×10⁻⁷ 时 best trial RMS≈**0.66**，仍高于起点 **0.020**
- 末态 = 起步态：bed0 y_O₂≈0.09、CO/H₂=0；床顶 y_CO≈0.012、H₂=0；T 全钉在 seed≈1150 K。主导残差是 **bed 能量**（|F̂|≈0.15）和 bed9 CO/H₂

**结论**

1. 核心栈不是「收敛了但不对」，而是 **Inner 根本迈不出 Wirsum 步**。  
2. Jacobian 线性化是好的；非线性（刚性能量/氧化带）让 λ=1…10⁻⁷ 都升 ‖F‖。  
3. 出口 CO≈0.012 是 **Startwert/Pre-Inner 后的态**，不是 NR 走出来的。L4 平台的 CO≈0.12 来自改写初值后的另一次轨迹。  
4. 下一刀应在核心栈上改 **Startwert（能量/T–组成）或正则化 Newton（LM/更小信任域且保持方向）**，不要加回 absorb/local2。

### 5.45 PTC 能接住 Wirsum 步，但步长无物理位移（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_lm_ptc_o1.json`、`data/phase2_hamel_wirsum_core_ptc_alpha_o1.json`。

| arm | n_acc | best trial RMS | λ | CO | 备注 |
|-----|-------|----------------|---|-----|------|
| Newton CORE o1 | 0 | ~10³–10⁵ | — | 0.012 | §5.44 |
| LM μ₀=1e-4 o1 | 0 | 536 | 1/2048 | 0.012 | 线性缺陷 0.29，FTB 后 935 |
| PTC α=1、12 次回溯 | 0 | **0.044** | 1/2048 | 0.012 | 已靠近起点 0.020 |
| **PTC α=1、24 次回溯 o1** | **1** | **0.020409** | **3.8e-6** | 0.012 | 首次 CORE 接受 |
| 同上 o4 | 2 | — | 4e-6, 1e-6 | 0.012 | rms 0.02041→0.02017；energy 仍 0.0685 |

**结论**

1. 把方向缩到足够小（PTC）并把 Wirsum 二分加到 ~2⁻¹⁸，‖F‖₂ **可以**严格下降。  
2. 接受步不改变出口 CO/H₂，能量 RMS 不动。这是数值上的 ε 步，不是离开 Startwert。  
3. PTC 不是 Wirsum 原文（原文是 \(F'Δx=-F\) 再回溯）；**不要**把 PTC 写成同构默认。  
4. 真正的下一刀仍是 **Startwert**：1150 K + bed0 残氧 + 合成气≈0 的能量不平衡，NR 用 ε 步走不出去。

### 5.46 初值确认：拧 T 种子几乎不动能量/O₂（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_init_seed_scan.json`（CORE，只 init，无 NR）。

| 种子 | T_bed0 | y_O₂ bed0 | y_CO 顶 | energy RMS |
|------|--------|-----------|---------|------------|
| R12 政策 1150（现默认） | 1150 | 0.091 | 0.012 | 0.108 |
| budget 900 + 冷锚 | 900 | 0.089 | 0.001 | 0.108 |
| adiabatic_anchor | 1180 | 0.095 | 0.000 | 0.111 |

**结论：** 是初值，但不是「把 1150 改回 900」就能解。三种 T 下 bed0 都是 **残氧 ~9%、CO/H₂=0**，能量 RMS 几乎相同。下一步应查 Vorab 快氧化 x₀ 与能量闭合，而不是再扫底温。见 §5.47。

### 5.47 bed0 本身不可行：fastox 把合成气抽空，氧与源项还在（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_bed0_init_anatomy.json`（CORE，只 init，无 NR；`PYTHONHASHSEED=0`）。

**fastox 开（现默认）bed0**

| 量 | 值 | 含义 |
|----|----|------|
| T / T_zu_gas | 1150 K / **293 K** | 冷一次风打进热格子 |
| y_O₂ / N_O₂ | 0.091 / 5.97 mol/s | 进料 O₂=11.23，动力学只吃 3.25 |
| N_CO / N_H₂ / N_CH₄ | **0** | 快氧化抽空 holdup |
| R_CO / R_H₂ / R_O₂ | +2.24 / +2.51 / −3.25 mol/s | 炭/气化仍在产合成气、耗氧 |
| E | **−2.47 MW** | Eq.2.7 远未闭合 |
| fastox ξ | ξ12=37.5，ξ5=27.8（全床合计） | 只吃气相燃料，**不**吃固体炭 |

**fastox 关：** holdup 里有 CO/H₂/CH₄，但 R12 裸速率把 R_O₂、R_H₂ 撑到 ~10¹³ mol/s。这就是为何 x₀ 必须做快氧化壳。

**轴向（fastox 开）**

- **只有 bed0 残氧**（y_O₂=0.091）；bed1+ y_O₂=0。  
- **全床 y_H₂=0**，同时每格 R_H₂≈+3 mol/s：holdup 在快氧化流形上，源项不在。  
- 能量洞不独属于 bed0：bed1 −2.52 MW，bed9 −4.82 MW。bed0 的独特点是 **残氧 ∩ 正在产 H₂/CO**。

**结论**

1. 用户判断对：CORE 下 **bed0 的 Startwert 自己就不在准稳态流形上**，不是再拧全局阻尼能走出去的。  
2. 机制：Hamel §2.1 快氧化只消除「H₂/CO/CH₄/焦油 与 O₂ 共存」；氧有剩 → 合成气被置零，炭/蒸汽气化立刻把 H₂/CO 源打回残差。Newton 步一旦把 H₂ 填进 bed0，R12 就把试探残差打到 10³–10⁵（§5.44）。  
3. 关 fastox 更糟（裸 R12）。拧 T 无效（§5.46）。  
4. 下一刀应改 **bed0 的 x₀ 闭合**（氧耗尽 vs 合成气库存 vs 能量），不是加回 L4、不是 PTC 默认。见 §5.48。

### 5.48 底格 Startwert：R1 投影关掉残氧，但 O₂ 流残差未闭合（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_bed0_char_r1_init.json`、`phase2_hamel_wirsum_core_bed0_char_r1_o1.json` / `_o4.json`。

气相 fastox 之后补 Hamel R1 计量（`α=1/φ_c`）：残氧接到炭库存，产物进 CO/CO₂。开关 `vorab_transport_x0_char_oxidation_o2_closure_thesis`（默认开，随 fastox）。

| 量 | 仅气相 fastox（§5.47） | +R1 投影 |
|----|------------------------|-----------|
| bed0 y_O₂ / N_O₂ | 0.091 / 5.97 mol/s | **0 / 0** |
| E | −2.47 MW | **−0.34 MW** |
| R_O₂ | −3.25 | **0**（格子里没有 C_O₂，R1 熄火） |
| N_zu O₂ | 11.23 | 11.23 |
| 单格 rms_scaled | — | 0.069 |

**O₂ 方程：** 投影后 `N_out(O₂)=0` 且 `R=0`，残差 ≈ 整股进料 11.23 mol/s，比投影前 ~2 mol/s **更差**。能量变好是因为 CO₂ 生成焓进了 H_out，不是动力学吃掉了进气氧。

**单格 `solve_cell`：** 在 R1 之后（氧已零）rms 0.0686→0.0695，拒收。在气相 fastox 后、残氧仍在时开求解，裸 R12 把 least_squares **卡住**（>3 min，已杀）。故 `vorab_init_solve_bottom_cell_thesis` **默认关**；若打开则 `skip_homotopy` + `max_iter=12`，冷根/重叠变差回滚。

**全局 NR：** CORE o1/o4 的 `accepted_lambda_history` 全是 `None`（线搜索失败，不是接受步）。出口 CO 仍 0.012。R1 投影**没有**打开 Newton 路径。

**结论**

1. Hamel Kap.7「氧在第二格耗尽」是**解的轴向结果**，不能靠把 holdup 氧直接置零来冒充；置零会熄灭 R1。  
2. 底格要求 `R_O₂ ≈ N_out − N_zu`，必须在 **C_O₂>0** 时让动力学真正吃进气。  
3. 下一刀不是再置零氧，而是让 N、R(C)、N_ex 同流形；见 §5.49。

### 5.49 底格两相 Startwert：O₂ 气泡 / 合成气悬浮（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_bed0_eq24_two_phase_init.json`、`..._o1.json`。

审计结论（未改方程）：水力与 K_bd 正常；N_ex≈0 是因为 fastox `full_hydro` 把 y_b=y_d；R1 置零熄火；VM 只进 R 不进 N。三者耦合，改一处会牵动其余。

**落地（只改 x₀）**

- CORE：`vorab_transport_x0_char_oxidation_o2_closure_thesis=False`
- CORE：`vorab_bed0_eq24_two_phase_startwert_thesis=True`（生产默认仍关）
- `apply_bed0_eq24_two_phase_startwert`：冻结 fastox 态 R，一次投影 `N≈N_zu+R`；O₂ 全在气泡、CO/H₂/CH₄ 全在悬浮；N₂/H₂O/CO₂ 线性化 Eq.2.2。密集相 O₂ 必须**精确 0**（1e-12 仍被 R12 的 C_O₂^0.5 放大到 10⁷ mol/s）。

| 量 | R1 置零（§5.48） | 两相 Startwert |
|----|------------------|----------------|
| bed0 N_O₂ / y_O₂,b | 0 / 0 | **8.04 / 0.209**（气泡） |
| N_H₂ 悬浮 | 0 | **2.49** |
| E | −0.34 MW | **−0.038 MW** |
| 同相 O₂·H₂ | 0 | 0（分相） |
| N_ex(O₂) | 0 | **+47 mol/s** |
| 单格 rms_scaled | 0.069 | **0.32**（相残差） |
| CORE o1 n_acc / y_CO,top | 0 / 0.012 | **0 / 0.012** |
| LS best trial rms | ~10³–10⁵ | **3610**（λ=1/2048） |

**新根因：** R1 只在悬浮相。O₂ 全在气泡 ⇒ C_d,O₂=0 ⇒ R1=0 ⇒ Eq.2.2 仍把 O₂ 以 N_ex≈47 打进悬浮相残差。Newton 要消这条残差，只能把 O₂ 搬进悬浮相，于是与已种子的 H₂ 同相，裸 R12 爆炸。  
所以「不置零氧 + 分相 + VM 进 N」在**赋值点**上自洽，但 **Newton 第一步的方向就是 R12 悬崖**。

下一刀（仍不改 Hamel 方程）：让线搜索试探不跨相点燃 R12。见 §5.50。

### 5.50 线搜索按相 fastox：压住 R12 悬崖，但仍非下降方向（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_per_phase_fastox_o1.json`。

合计 `local_participating` fastox 会把气泡 O₂ 与悬浮 H₂ 当成一锅烧，毁掉 §5.49 分相。新模式 `phase_update="per_phase"`：各相独立做 R12/R5/R6/R10 计量，**不跑 R1 置零**。

开关 `nr_line_search_per_phase_fastox_thesis`（生产默认关，CORE 开）。**不是** `nr_line_search_fast_oxidation_projection_thesis`（合计 LS fastox 仍关）。

| 试探 | best trial rms（λ=1/2048） |
|------|---------------------------|
| 无 LS 投影（§5.49） | 3610 |
| 按相 fastox 仅 bed0 | 3610（上段格 Newton 带进微量 O₂，裸 R12） |
| 按相 fastox 全 bed | **6.74** |

CORE o1：12 次按相投影，`n_acc=0`，y_CO,top 仍 0.012。投影让 F 有限，但投影点不在 JΔx 的下降方向上（6.74 ≫ 当前 0.047）。

下一刀：在无同相 overlap 的流形上取下降步。见 §5.51。

### 5.51 分相切向 clip：禁止同相增加 O₂∩合成气（2026-08-20）

证据：`data/phase2_hamel_wirsum_core_manifold_dx_probe.json`、`data/phase2_hamel_wirsum_core_same_phase_clip_o1.json`、`data/phase2_hamel_wirsum_core_same_phase_clip_o4.json`。

**诊断（未改方程）**

第一步 clipped Newton 的最大分量是温度（|dT|~100 K），气相泄漏很小但致命：

- bed1 `dN_d,O2=+0.022 mol/s` 进入已有 H₂/CO 的悬浮相
- bed0 `dN_b,CO=+8.5×10⁻⁵ mol/s` 进入已有 8 mol/s O₂ 的气泡

λ=1/2048 时：无 fastox rms=3610（5 格同相 overlap）；按相 fastox 后 rms=6.74，最差是 **bed0 H₂ 残差**（fastox 烧掉气泡里那点 CO 后把 N_ex 打歪）。单格加 O₂ 再 fastox 解释不了 6.74，必须看**整步 Δx**。

在 raw Δx 上把「同相增加 overlap」的分量置零，再 `_clip_dx`，再按相 fastox：

| 试探 | rms_scaled | n_acc | y_CO,top |
|------|------------|-------|----------|
| 当前点（init） | 0.0476 | 0 | 0.012 |
| clip+fastox λ=1/2048（§5.50） | 6.74 | 0 | 0.012 |
| 切向 clip + fastox λ=1（probe） | **0.0365** | — | — |
| CORE o1 | **0.0373** | **1**（λ=1） | 0.0065 |
| CORE o4 | **0.0345** | **4**（全 λ=1） | 0.0086 |

o4 底格仍 `N_d,O2=0`，`N_b,O2` 8.04→6.08，`N_d,H2` 2.49→0.26。下降步在走，但 R1 仍未点亮。y_CO,top 先掉后略回，远低于 LU 0.13——不要为此开平台桥。

**落地（只改步长流形）**

- `nr_clip_same_phase_oxidizer_fuel_step_thesis`（生产默认关，CORE 开）
- `_zero_same_phase_oxidizer_fuel_newton_step`：某相已有合成气则禁止 `dN_O2>0`；某相已有 O₂ 则禁止 `dN_{H2,CO,CH4}>0`。负步保留。
- 接线在 `_clip_dx` **之前**，避免 overlap 分量去抢气相 trust region。
- 仍开按相 fastox；仍关合计 LS fastox。

这不能单独点亮 R1（悬浮相仍无 O₂），只是让线搜索看到真实下降方向。O₂ 进悬浮相吃炭仍要靠后续可接受步慢慢走，或 N_ex/气泡侧闭合。

下一刀：见 §5.52。

### 5.52 底格 R1 残氧种子：从气泡挪 O₂ 进悬浮相（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_r1_o2_leak_probe.json`、`data/phase2_hamel_wirsum_core_r1_seed_o4.json`。

**诊断**

- 裸 O₂∩H₂（无 fastox）即使 1e−8 mol/s 也把 R12 打到 10⁸。
- 只在 holdup 上加悬浮相 O₂ 再按相 fastox：加到超过计量后 R1 点亮，rms 不炸（init 加 8 mol/s：0.048→0.036，Rd,O₂≈−4）。
- 第一步 Newton 的 `dN_d,O2≈0`（底格并不想进氧）；上段 `dN_d,O2=+0.085` 是爆炸源。把 clip 改成「允许 O₂ 进燃料相」后 CORE o4 **退回 n_acc=0、试探 rms=6.74**。该开关保持默认关。

**落地（只改 x₀）**

- `vorab_bed0_r1_oxidizer_seed_mol_s_thesis`（生产默认 0，CORE=2.0 mol/s）
- `apply_bed0_r1_oxidizer_seed`：在 §5.49 分相之后，把 `stoich(H₂/CO/CH₄)+seed` 从 `N_b,O2` 转到 `N_d,O2`，再按相 fastox。不跑 R1 置零。
- 仍开 §5.51 对称切向 clip；`nr_clip_allow_oxidizer_into_fuel_thesis=False`。

| 量 | §5.51 o4 | §5.52 init | §5.52 o4 |
|----|----------|------------|----------|
| rms_scaled | 0.0345 | **0.0314** | **0.0162** |
| n_acc | 4（λ=1） | — | **3**（0.5/1/0.5，第 4 步拒） |
| split | 0.124 | — | **0.028** |
| bed0 N_d,O2 | 0 | **1.84** | 3.04 |
| bed0 Rd,O2 | 0 | **−1.64** | **−1.76** |
| N_ex,O2 | 34.6 | 10.9 | 1.43 |
| y_CO,top | 0.0086 | — | 0.0067 |

R1 已点亮，O₂ 交换从 47 收到 1.4。y_CO 仍远低于 LU 0.13；底格 T 从 1150 落到 1000 K。下一刀：见 §5.53。

### 5.53 能量/温度：不要升温 Startwert，改每步 |ΔT| 帽（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_energy_t_probe.json`、`data/phase2_hamel_wirsum_core_energy_t_tcap_ab.json`、`data/phase2_hamel_wirsum_core_energy_t_tcap40_o4.json`、`data/phase2_hamel_wirsum_core_energy_t_tcap40_o8.json`。

**诊断**

- 围栏：`T_ref−150` → bed0 Tmin=1000.23 K。§5.52 默认帽 600 K，第一步 `|dT|` 顶到围栏宽度（350 K），T 一刀落到下沿，第 4 步拒。
- Eq.2.7 在冻结组成、只动 T 时**不是**「E>0 就升温」。干净 scan：T=1100 时 |E| 与 rms 最好（E=−0.33 MW，rms=0.0297）；升到 1300–1500 气相残差变差（0.034→0.078）。既有 `preproject_bed_temperature_to_local_energy_closure` 因 `improvement_factor=0.75` 拒收。
- 强制 +80 K Startwert：o4 仍 n_acc=3、T=1000、rms=0.016。升温不是刀。
- `nr_inner_abs_dT_budget_K_thesis=80`：每个 inner 用尽预算后后续步 `|dT|≈0`，o4 只接受 2 步后拒。

**落地（只改步长流形）**

- `nr_inner_t_step_default_cap_K_thesis`：生产默认 600 K，**CORE=40 K**。
- 不接能量 T 投影，不加宽围栏，不改 k / K_bd / Q。

| 量 | §5.52 o4（帽 600） | §5.53 o4（帽 40） | §5.53 o8（帽 40） |
|----|-------------------|------------------|------------------|
| rms_scaled | 0.0162 | **0.0179** | **0.0138** |
| n_acc | 3（第 4 拒） | **4（全 λ=1）** | 29 |
| bed0 T | 1000.23（围栏） | **1070** | 1000.23（围栏） |
| bed0 E | −0.31 MW | **+0.19 MW** | −0.57 MW |
| bed0 N_d,O2 | 3.04 | 3.14 | 3.35 |

帽 40 把围栏猛砸推迟：o4 仍在 1070 K、能量根未翻到 ~900 K，四步都能走。o8 组成演化后根仍下移，T 会再贴围栏，但 rms 继续降到 0.014。下一刀：见 §5.54。

### 5.54 组成–焓：fastox/R1 种子后重贴固相 Eq.2.6（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_energy_enthalpy_decomp_init.json`、`data/phase2_hamel_wirsum_core_energy_solid_gap_init.json`、`data/phase2_hamel_wirsum_core_energy_holdup_realign_o4.json`。

**诊断**

- energy_rms 主导格是 **bed1（E≈−6.5 MW，Ê=−0.34）** 和 **bed9（−4.8 MW）**，底格只排第四。bed1 的 E(T) 在 1000–1200 K 几乎平坦，升温/降温 Startwert 无效。
- Eq.2.7 固相出口用 `K·m`。bed1/bed9 在 fastox + Eq.2.4 + R1 种子之后 **欠库存**（hold/support≈0.84），`m_in−m_out≈0.41 kg/s` → 结构焓洞。对齐函数原先只在这些步骤**之前**跑一次。
- 对 bed1/bed9 做 T Newton（帽 50 K）几乎不降 |E|；o4 轨迹与 §5.53 相同。
- 种子后再 `reconcile_upper_bed_solid_holdup(..., allow_reduce=True)`：init energy_rms 0.101→**0.074**，rms 0.031→**0.029**；bed9 E −4.8→−1.3。`include_negative_R` 会把能量弄得更差，保持关。

**落地（只改 x₀ 顺序）**

- `_run_bed0_eq24_two_phase_startwert` 与 fastox **之后**再调 `_reconcile_upper_bed_solid_holdup_if_enabled(allow_reduce=True)`。
- 早段对齐仍只升不降。无新 CORE 开关；沿用 `nr_init_align_upper_bed_solid_holdup_thesis`。
- 不改 Eq.2.7，不加宽围栏，不改 k / K_bd。

| 量 | §5.53 init / o4 | §5.54 init / o4 |
|----|-----------------|-----------------|
| rms_scaled | 0.0314 / 0.0179 | **0.0286** / 0.0190 |
| energy_rms | 0.101 / 0.056 | **0.074** / 0.061 |
| n_acc o4 | 4（全 λ=1） | 4（全 λ=1） |
| bed0 T o4 | 1070 | 1070 |
| bed1 E init | −6.45 MW | **−5.42 MW** |
| bed9 E init | −4.82 MW | **−1.30 MW** |

下一刀：见 §5.55。

### 5.55 bed1 局部 snap：全床联立对齐贴不紧 char（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_energy_holdup_bed1_snap_o4.json`。

**诊断**

- §5.54 后 bed1 仍 E≈−5.4 MW。E(T) 平坦。`auf_upflow_closure` 把 E 打到 −12 / rms 0.073，**禁止当 init 默认**。
- 剩余洞拆开：VM 入流 `K=0` 约 **−1.82 MW**；char/ash 在末次 BC 后又欠约 10%（约 1.5 MW）。气相最大项是底格 **N₂ 相分裂** Ê=+0.215/−0.215（合计≈0），不是守恒错误。
- 全床再跑一轮 `allow_reduce` 会对抢 auf/ab，bed1 char 仍贴不拢。只对 **bed1** 做「align → BC」2–3 次后 char/ash = support。

**落地（只改 x₀）**

- `_snap_bed_holdup_to_local_support(bed_indices=(1,), max_passes=3)`，接在种子后全床 reconcile 之后。
- 不改 Eq.2.7，不加宽围栏，不改 k / K_bd。

| 量 | §5.54 init / o4 | §5.55 init / o4 |
|----|-----------------|-----------------|
| rms_scaled | 0.0286 / 0.0190 | **0.0280** / **0.0170** |
| energy_rms | 0.074 / 0.061 | **0.066** / **0.050** |
| n_acc o4 | 4（全 λ=1） | 4（全 λ=1） |
| bed0 T o4 | 1070 | 1070 |
| bed1 E | −5.42 / −3.20 | **−4.16** / **−1.78** |

下一刀：见 §5.56。

### 5.56 底格惰性相重分：种子后 N_ex 与 N₂/H₂O/CO₂ 脱节（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_energy_inert_split_o4.json`。

**诊断**

- §5.55 后 gas 最大项是底格 N₂ 相分裂 Ê=+0.215/−0.215（合计≈0）。CO₂/H₂O 同类对打。总和 Eq.2.4 已闭合，不是守恒错误。
- §5.49 分相用**冻结 R** 与分相前 `n_d_tot/n_b_tot`。R1 种子 + 按相 fastox 改了 y 与 H₂O/CO₂ 总和后，N_ex 仍按旧分相。
- 探针：只对 N₂/H₂O/CO₂ 重跑线性化 Eq.2.2（`n_tot=zu+R`，`n_d=(zu_d+R_d+α n_tot)/(1+α+β)`），不碰 O₂/合成气。init rms 0.028→0.018；N₂ Ê 0.215→0.03。

**落地（只改 x₀）**

- `apply_bed0_inert_exchange_split_startwert`：门控 = `vorab_bed0_eq24_two_phase_startwert_thesis`，**无新 CORE 开关**。
- 接在 R1 种子之后，以及 bed1 snap 之后（末次，不再 `apply_bc`）。
- 不改 k / K_bd / Eq.2.7，不加宽围栏，不把 split merit 写回接受准则。

| 量 | §5.55 init / o4 | §5.56 init / o4 |
|----|-----------------|-----------------|
| rms_scaled | 0.0280 / 0.0170 | **0.0178** / **0.0174** |
| gas_rms | 0.0344 / 0.0191 | **0.0189** / **0.0189** |
| energy_rms | 0.066 / 0.050 | **0.059** / **0.055** |
| n_acc o4 | 4（全 λ=1） | 4（全 λ=1） |
| bed0 T o4 | 1070 | 1070 |
| bed1 E | −4.16 / −1.78 | **−2.12** / **−0.80** |
| N₂ Ê_dense init | +0.215 | **−0.031** |
| O₂ Ê_dense o4 | （init 0.16） | **0.018** |
| N_d,O2 | 1.84 / 3.10 | 1.84 / 3.11 |

o4 rms 与 §5.55 同量级，但 Startwert 已在流形上：Newton 改去收 O₂ 而不是对打 N₂。底格 T 仍一刀贴围栏（1000 K 下沿）。

下一刀：见 §5.57。

### 5.57 O₂ Startwert 刀全部回到同一吸引子（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_energy_o2_knife_o8.json`。

**诊断（只改 x₀，不改 k / K_bd / Eq.2.7）**

底格 O₂ 总和未闭合（N=5.05 vs zu+R=9.24，gap=+4.19），Ê_dense≈0.125。这不是守恒写错，是种子/fastox 后冻结 R 与活 R 脱节。试了五刀：

| 刀 | init 结果 | o4 |
|----|-----------|-----|
| 总和闭合、保留 N_d | N_b→7.4，N_ex→25，rms **0.029** | — |
| 下限线性化 Eq.2.2 | O₂ Ê 0.125→**0.005**，但 E1 −2.1→**−4.1**，ene 0.067 | — |
| 底格 T₀=1070 | init 几乎不变 | **与 §5.56 o4 逐位相同** |
| 合成气进气泡 + 全部 O₂ 进悬浮相 | init rms **0.039** | 回到同一吸引子（T=1070，rms 0.017） |
| 上段过量 H₂O 按 zu+R 缩 | 能量在格间挪，ene_rms 不动 | — |

冻结组成下 E₀(T) 几乎平坦（1000 K：−2.21；1200 K：−2.83）。主导项是燃料 `H_szu≈−8.2 MW` vs `H_sout≈−4.8 MW`（VM `K=0` 不进出口），不是 O₂ 分相。

**o8（§5.56 栈，无新开关）**

- rms **0.0157**，gas 0.0166，ene 0.052
- 底格 T=**1000.23 K**（围栏下沿 = T_ref−150）
- E0 仍 **−2.67 MW**；E1 −2.12→**−0.44 MW**
- 前 8 步全 λ=1，随后大量极小 λ / 拒步

Newton 自己会把 O₂ Ê 收到 ~0.018，再把 T 砸进围栏。加宽围栏只会得到更冷的假根，E0 不降。

**落地**

- **不改 CORE**。`apply_bed0_inert_exchange_split_startwert` 注明禁止对 O₂ 套同一公式。
- 勿把 VM 补进 Eq.2.7 出口当默认（那是工程桥，不是 Hamel 同构）。

下一刀：见 §5.58（用户要求下探底格围栏）。

### 5.58 底格围栏下沿从 1000 K 放到 800 K（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed0_tmin_ab.json`。

**动机。** §5.57 o8 贴在 T_ref−150=1000.23 K，Newton 还在往下走。冻结组成扫描（unpack 不夹 T）：800 K 附近 rms 略降，650 K 气相不连续，不能当默认。

**实验（仅改底格 Tmin，上段仍 150 K；tcap=40；`PYTHONHASHSEED=0`）**

| 底格 margin / Tmin | 迭代 | T0 | rms | gas | ene | E0 MW |
|--------------------|------|----|-----|-----|-----|-------|
| 150 / 1000（§5.57） | o8 | **1000.2**（贴栏） | 0.0157 | 0.0166 | 0.0525 | −2.67 |
| 250 / 900 | o8 | **900.2**（贴栏） | **0.0130** | 0.0135 | 0.0445 | −2.20 |
| 350 / 800 | o8 | **877.2**（未贴栏） | **0.0127** | 0.0134 | 0.0427 | −2.10 |
| 350 / 800 | o12 | 878.8 | 0.0140 | 0.0137 | 0.0513 | −0.08 |
| 450 / 700 | o12 | 878.8 | 与上一行相同 | | | |

o8：900 仍贴栏，800 才让 Newton 停在 ~878 K。再放到 700 无增量。o12 把 E0 收到近 0，但 energy_rms 回升、总 rms 变差。

**落地**

- 新旋钮 `nr_temperature_fence_bed0_lower_margin_K_thesis`（生产默认 `None`=全局 150）。
- **CORE=350**（Tmin≈800 K）。上段围栏不动。
- 实现：`init_precalc_step._install_temperature_fence_from_vorabrechnung`；测试 `test_bed0_fence_lower_margin_only_widens_bottom_cell`。
- 勿把全局/绝对下沿放到 600 K 以下；勿改 k / K_bd / Eq.2.7。

下一刀：见 §5.59。

### 5.59 底格 Startwert T=880 K：o8 不再先花在降温上（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed0_tstart_ab.json`。

**动机。** §5.58 o8 从 Vorab T₀=1150 降到 877 K，8 步大多在降温。O₂ Ê 仍≈0.10。围栏中心保持 T_est（Tmin≈800），只把 NR 见到的 `cell.T` 放到 880 K。

| 底格 Startwert T | o8 T0 | rms | gas | ene | O₂ Ê | NdO2 |
|------------------|-------|-----|-----|-----|------|------|
| 1150（Vorab） | 877.2 | 0.0127 | 0.0134 | 0.0427 | 0.103 | 2.67 |
| **880** | **842.8** | **0.0110** | **0.0100** | 0.0424 | **0.032** | 3.37 |

从 880 出发，Newton 仍继续降温到 843 K（未贴 800 K 栏）。能量仍由 bed0 E≈−2.1 MW 主导（Ê≈−0.11），ene_rms 几乎不动。

**落地**

- `vorab_bed0_nr_startwert_T_K_thesis`（生产 `None`）。在惰性分相之后调用 `_apply_bed0_nr_startwert_T`，贴已装围栏。
- **CORE=880 K**。不改 k / K_bd / Eq.2.7，不加宽全局 Tmin。

下一刀：见 §5.60。

### 5.60 Startwert T=820 K：850 会砸栏，820 升回 ~840（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed0_tstart_850_820.json`。

**动机。** §5.59 从 880 降到 843。再靠近吸引子时，850 与 820 行为相反。

| Startwert T | o8 T0 | rms | gas | ene | O₂ Ê | 备注 |
|-------------|-------|-----|-----|-----|------|------|
| 880（§5.59） | 842.8 | 0.0110 | 0.0100 | 0.0424 | 0.032 | 未贴栏 |
| 850 | **800.7** | 0.0100 | 0.0095 | 0.0371 | 0.033 | **贴 800 K 栏**（tcap 40，两步砸到底） |
| **820** | **839.9** | **0.0100** | 0.0100 | **0.0353** | 0.029 | 从 820 **升温**到 ~840 |

850 的 800 K 是 tcap 假冷根，不是能量吸引子。880 与 820 都停在 ~840 K。

**落地**

- **CORE Startwert T=820 K**（生产仍 `None`）。
- 勿把 CORE 设成 850。勿再靠加宽围栏去「允许 800 K」。

下一刀：见 §5.61。

### 5.61 bed0 焓洞是 VM zone=3 + K_VM=0；zone=1 填洞但 o8 变差（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed0_vm_zone_ab_o8.json`。

**结构。** Eq.2.7 固相出口是 `K·m`。`K` 只给 char/ash（B-tier：VM 不走 auf/ab）。燃料 VM 的 `H_szu≈−3.09 MW` 进底格，但 `H_sout` 不含 VM。热解产物只能走 `H_gout`。Phase2 / Hamel Kap.2.1 把脱挥发分到 **zone=3**：bed0 cap≈1/3，`R_vm≈−0.123 kg/s` vs `m_zu,VM=0.369 kg/s`。剩余 2/3 的释放在 bed1/bed2，而 VM 质量并不随 `K` 上行。single-shot 源项冻结后，NR 把 T 从 820 升到 840 **也不会**提高 bed0 的 `x_vm`（cap 钉死）。

char 几乎已经闭合：`H_szu,char≈−5.11` vs `H_sout≈−5.04`。把 bed0 char/ash 再 snap 到 Eq.2.6 support，init `|E0|` 只动 ~0.04 MW，不值得落地。

**A/B（只改 Vorab 分区，不改 Eq.2.7 / k / K_bd）。** 同 `HAMEL_WIRSUM_CORE`，`PYTHONHASHSEED=0`，o8。

| arm | init E0 | init rms / ene | o8 T0 | o8 rms | o8 gas / ene | o8 H₂ Ê | 备注 |
|-----|---------|----------------|-------|--------|--------------|---------|------|
| zone=3（Phase2 默认） | −1.65 MW | 0.0185 / 0.061 | **840** | **0.0116** | 0.0106 / 0.045 | 0.016 | `x_vm,eff=0.33`；`n_acc=0` |
| **zone=1** | **−0.80 MW** | 0.0174 / 0.051 | 1014 | **0.0187** | 0.0208 / 0.056 | **0.053** | `x_vm,eff≈1`；H₂/CO Ê ~3× |
| zone=0 | −47 MW | 3.80 / 18.1 | — | — | — | — | `H_sout≈−451 MW`，holdup 炸 |

zone=1 把全部热解气源堆进 bed0。两相 Startwert + 按相 fastox 仍把 `N_d,H₂/CO` 打成 0（避免与残氧同相），于是 `R_H₂` 大、holdup 为 0 → 气相残差主导，o8 比 zone=3 **差**。zone=0 关掉分区映射，上层床 holdup 失支。

**落地**

- CORE **钉死** `vorab_vm_devolatilization_zone_cells_thesis=3`（与 Phase2 一致，防止以后有人用 zone=1「修能量」）。
- **不**开 `_solid_energy_passthrough` 当床层默认，**不**给 VM 加 `K`，**不**把 VM 补进 Eq.2.7 出口。
- 勿 snap bed0（相对 bed1 局部 snap 无增量）。

下一刀：见 §5.63。

### 5.62 bed3 +2.3 MW 是 H₂O 过填；透传填洞但 o8 砸栏（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed3_h2o_passthrough_o8.json`。

**结构。** CORE init 能量主导不只是 bed0 VM 洞（§5.61）。bed3 `E≈+2.28 MW` 时 char/ash 已贴 Eq.2.6（`in−out≈0`）。气相：`H_gout≈−4.27 MW` vs 上游 bed2 `H_gout≈−1.96 MW`（同 T≈1150 K）。组分跳变是 **H₂O 20.6 vs inflow 12.5**（gap +8.0 mol/s，与 §5.36 同一过填），外加 CO₂/CH₄/TAR 的 Gibbs x0 拷贝。bed3–bed9 的 `H_gout` 被复制成同一份。

对 bed2/3/7/9 再做局部 char snap：bed7/bed9 焓能降，但会把洞推到邻居（bed8 变差）；bed3 几乎不动。不是 holdup 刀。

**H₂O 透传 A/B**（不改 Eq.2.7 / k / K_bd）。`PYTHONHASHSEED=0`，o8。

| arm | init E3 | init H₂O gap | init rms / ene | o8 T0 | o8 rms | o8 E3 | 备注 |
|-----|---------|--------------|----------------|-------|--------|-------|------|
| CORE（旗标关） | +2.28 | +8.04 | 0.0185 / 0.061 | **840** | **0.0116** | +1.56 | `n_acc_all=40` |
| 旗标开、仅映射时 | +1.09 | +2.29 | 0.0174 / 0.056 | — | — | — | fastox 后 bed2 的 H₂O 已降，bed3 未再贴 |
| **旗标开 + fastox 后再贴** | **+0.29** | **0** | 0.0191 / 0.059 | **800.2** | **0.0125** | +0.43 | **贴 800 K 栏** |

fastox 后再贴：`apply_post_fastox_post_staged_h2o_passthrough`（门控仍是既有 `vorab_a_tier_post_staged_h2o_passthrough_thesis`）。这是 §5.28 透传在 fastox 之后的正确性，不是新物理方程。

o8 把 T0 砸进围栏，总 rms 变差。与 §5.29 双吸引子是另一条路，但 CORE（无 LSQ）同样不能靠它压能量。

**n_acc 误读。** §5.61 表里的 `n_acc=0` 读的是**末次 inner** 的 `nr_accepted_lambda_history`。外层累计 `n_acc_all≈40`（8 个 outer）。T 从 820 升到 840 就是这些接受步，不是围栏外的暗改。

**落地**

- 保留 fastox 后 H₂O 再贴，**CORE 旗标仍 False**。
- 勿把 H₂O 透传写成 Wirsum/Hamel 同构默认。勿对 bed3+ 做全床 char snap。

下一刀：见 §5.63。

### 5.63 bed3 产物 x0 跳变与 H₂O 同族；透传填 E3 但 o8 rms 变差（2026-08-24）

证据：`data/phase2_hamel_wirsum_core_bed3_syngas_passthrough_o8.json`。

**结构。** `seed_major_product_holdup_from_vorab_gas_balance` 在 fastox **之前**，且 zone=3 时 `_resolve_product_holdup_bed_cell_limit` 只覆盖 bed0–2。fastox 改 bed2 的 CO/CO₂/CH₄/TAR 后，bed3+ 仍停在映射时的 Gibbs 拷贝（CORE init：CO₂ 3.67 vs inflow 1.87，CH₄/TAR=0 vs 0.86/0.11）。与 §5.62 的 H₂O 过填是**同一时序洞**，物种不同。

**产物透传 A/B**（不含 H₂O，不改 Eq.2.7 / k / K_bd）。`PYTHONHASHSEED=0`，o8。`solve()` 会重跑 init，须走配置旗标，不能事后改 x0。

| arm | init E3 | init CO₂ / CH₄ / TAR2 | init rms / ene | o8 T0 | o8 rms | o8 E3 | 备注 |
|-----|---------|----------------------|----------------|-------|--------|-------|------|
| CORE（旗标关） | +2.28 | 3.67 / 0 / 0 | 0.0185 / 0.061 | **840** | **0.0116** | +1.56 | 非 None 接受 25；λ 条目 40 |
| **产物透传（H₂O 仍关）** | **+1.38** | **1.87 / 0.86 / 0.11** | 0.0186 / 0.060 | **834** | **0.0141** | **+0.05** | 贴 inflow；未贴栏；gas_rms 变差；外层 5 步早停 |

Init 上能量 rms 几乎不动：剩下的 E3 仍是 H₂O 20.6 vs 12.5。o8 把 bed3 焓几乎填平，但总 rms 与气相 rms 都变差，T0 略冷、未贴 800 K 栏。与 §5.62 H₂O 透传「填洞但砸栏」不同，这条是「填洞但走歪」。

**n_acc 计数。** 外层拼起来的 `nr_accepted_lambda_history` 长度 ≈40（含 None）。非 None 接受步 CORE=25、产物臂=30。§5.62 表里的 40 是列表长度，不是接受次数。

**落地**

- `apply_post_fastox_post_staged_syngas_passthrough`；旗标 `vorab_a_tier_post_staged_syngas_passthrough_thesis` 默认 **False**。
- 接在惰性分相 / H₂O 再贴之后、Startwert T 之前。跳过 bed0 与 `zu_O2>0` 的格。
- **CORE 旗标仍 False。** 勿与 H₂O 透传叠开验收。勿写成 Wirsum/Hamel 同构默认。

下一刀：见 §5.64。

### 5.64 T0=840 是底格吸引子；o8 能量平台是 inner budget 用尽（2026-08-24）

证据：

- `data/phase2_hamel_wirsum_core_ls_plateau_o4.json`（短探针，**不能**当平台证据）
- `data/phase2_hamel_wirsum_core_ls_plateau_o8_o16.json`
- `data/phase2_hamel_wirsum_core_upper_tmin350_o8.json`

**短探针陷阱。** `max_global_iter<=4` 走 `short_probe_mode`：inner budget=`max_iter`（o4=4），不是 `5×max_iter`。o4 全 λ=1、每步 max|ΔT|=40 K，T0 只到 826 K。平台诊断必须用 o8。

**o8 vs o16**（CORE，`PYTHONHASHSEED=0`；energy_rms 为事后 `global_residual` 分组）。`solve()` 的 `rms_scaled_energy_final` 数值更低，对比时勿混用。

| 臂 | budget used/total | n_acc / n_fail | n_outer | T0 | rms | energy_rms | E0 / E3 |
|----|-------------------|----------------|---------|----|-----|------------|---------|
| o8 | **40/40** | 25 / 15 | 8 | **840.0** | 0.0116 | **0.0451** | −2.18 / +1.56 |
| o16 | **80/80** | 55 / 25 | 13 | **840.2** | **0.0105** | **0.0354** | −1.99 / +0.96 |

T0 几乎不动：840 K 是底格温度吸引子，不是「还没走完」。o8 末段仍在 λ=1 接受；停是因为 thesis inner budget=`max_iter×5=40` 用尽，不是线搜索 merit 死掉。多给预算，组成 Newton 继续压 energy（E3 1.56→0.96），T0 仍钉在 840。o16 在 outer 13 早停（`nr_outer_merit_early_stop_min_iters=8`），budget 仍用满 80。

每步 `max_abs_dT_step=40` 来自上段格，不是底格。o8 末态 bed7–9 贴 **T_est−150≈1000.8 K**（全局下沿）。主导能量格：bed0（VM 洞，§5.61）、bed3（H₂O 过填，§5.62）、以及贴栏的 bed7/bed9。

**上段围栏 A/B**（o8；`solve()` 的 rms/energy 字段）。T0 仍 ~841 K，未砸 800 K 栏。

| arm | T7–9 | rms | ene | E9 |
|-----|------|-----|-----|----|
| CORE 上段 margin=150 | **1000.8 / 1000.9 / 1000.9**（贴栏） | 0.01126 | 0.0362 | +1.09 |
| 上段 margin=350（Tmin≈800） | **979 / 973 / 982** | 0.01084 | 0.0344 | +0.86 |

离开 1000 K 栏，但增益在小数点后第四位，且小于 o16 多跑组成步的收益。§5.58 明确「上段围栏不动」仍然成立。

**落地**

- 新旋钮 `nr_temperature_fence_upper_bed_lower_margin_K_thesis`（生产 / CORE 均为 **None**）。
- 实现：`init_precalc_step._install_temperature_fence_from_vorabrechnung`；测试 `test_upper_bed_fence_lower_margin_does_not_change_bed0`。
- **不**把 CORE o8 改成 o16；验收栈仍是 o8。多预算是组成收敛，不是新物理方程。
- 勿把全局/上段 Tmin 放到 600 K 以下；勿改 k / K_bd / Eq.2.7。

下一刀：T0 钉在 840 之后的**组成**残差（bed0 VM 焓洞、bed3 H₂O 过填仍在）。勿再叠 H₂O/产物透传当 CORE，勿靠加宽上段围栏压 rms。

---

## 6. Hamel / Siegen 文献求证状态

### 6.1 谱系文档（必读）

**`docs/siegen_cell_model_lineage.md`** — 完整书目、传承图、已读摘要、出版物清单索引。  
`docs/README.md` 已挂 L0 索引。

### 6.2 Hamel 正文已澄清

- §2.3：Wozny NR + 块三对角 + Nebenelemente；Eq. 2.9；阻尼见 **Wirsum (1998)**  
- **无** energy/split 对抗、`merit=max`、单支接受规则  
- 本地 Diss. PDF：`docs/_2001_VDI-Dis._…Brennstoffe.pdf`（gitignore）

### 6.3 已读全文（算法价值有限）

| 文献 | 路径 / 结论 |
|------|-------------|
| Rajan & Wen (1980) | 用户 `Documents/AIChE Journal - July 1980 - Rajan - ….pdf`；细胞/两相骨架；**非**全局 NR |
| Wozny et al. CIT (1987) | 用户 `Documents/Chemie Ingenieur Technik - February 1987 - Wozny - ….pdf`；2 页；细则回指 1983 Habilitation |
| Hamel & Krumm PT (2001) | 用户 Work 路径下 Powder Technology PDF；求解指 **1999 Final Report [10]** |
| Weil et al. ISCC (2003) | 用户 siegon 文献夹；平衡 vs 细胞；细节回指 Hamel Diss. |
| 锡根出版物清单 | 网页打印 PDF；检索索引，无早期 Diss. |
| **Wirsum (1997/98) VDI-6/383** | 本地 `docs/_1997_VDI-Dis._Wirsum_Klaerschlamm_blasenbildende_Wirbelschichtfeuerungen.pdf`；§2.3.3 已摘 → `docs/wirsum_1997_nr_solver_notes.md` |

### 6.4 仍缺全文 / 部分同构仍阻塞

| 文献 | 为何关键 | 现状 |
|------|----------|------|
| **Wozny (1983)** Habilitation | 块三对角 NR **消元/组装**细节 | 仍缺 |
| **Hartleben (1983)** | 细胞 + Wozny → 鼓泡床的桥 | 仍缺 |
| **Hamel/Fett/Krumm (1999)** Final Report IV A4-9021.1 | PT 文称含 mathematical solution procedure | 仍缺 |
| Wirsum 步长因子 \(2^{-i}\) 印刷字形 | OCR 将因子印残 | **接受准则已硬核实**；因子为高置信解读 |

**Wirsum 已闭合（细节表见 §6.6）：** 接受准则 = **‖F‖₂ 严格下降**（Eq. 2.23）+ 回溯阻尼；线性 Newton 步 \(F'\Delta x=-F\)（Eq. 2.21）；**不是** `merit=max(rms,split,energy)`。后者及 fastox/syngas/LSQ 门控继续标 **工程桥**。

### 6.5 Rajan vs 本模型骨架差异（摘要）

共用：轴向分室 + 两相气。  
分叉：固相仅悬浮相（Hamel）vs backflow cell；能量整 cell 一条 vs 分室；水力 Hilligardt 链 vs 燃烧管束；全局 NR+外层 Abgleich vs elutriation 迭代；气化动力学 vs 燃烧/石灰石/NOₓ。详见谱系文档读后摘要。

### 6.6 handoff 宣称 × Wirsum 全文核查表（2026-08-20）

对照：`docs/_1997_VDI-Dis._Wirsum_…pdf` §2.3.3（OCR：`data/wirsum_1997_extract/page_020–021.txt`）；笔记 `docs/wirsum_1997_nr_solver_notes.md`。

#### A. 现已可用 Wirsum 全文 **确认** 的点

| # | handoff / 既有说法 | Wirsum 原文依据 | 结论 |
|---|-------------------|-----------------|------|
| A1 | 远初值用 **gedämpfte Newton**（Hamel 回指 Wirsum） | §2.3.3：近解收敛快、远解常发散 → 采用 gedämpfte 变体 | **确认** |
| A2 | 求解属 **Wozny NR** 族 | 「基于 Wozny [149] 的 Newton-Raphson」 | **确认** |
| A3 | 雅可比呈 **块三对角**（邻格耦合） | 「blocktridiagonalen Funktionalmatritzen」；\(F_i\) 通常只依赖本格+邻格 | **确认** |
| A4 | 非邻格流可破坏纯三对角（Nebenelemente / CFB 等） | 明文：非邻格物质/能量流入为例外（CFB、非邻格颗粒对流热） | **确认**（与 Hamel Nebenelemente 叙事同族） |
| A5 | 接受准则是 **残差范数下降**，不是分量 max-merit | Eq. **(2.23)**：\(\|F(x^{\nu+1})\|_2 < \|F(x^{\nu})\|_2\) | **确认** |
| A6 | handoff §6.2 / §4：**无** energy↔split 对抗处方、**无** `merit=max`、**无** 单支 energy-only 接受 | §2.3.3 全文未见 split/energy 分量门控 | **确认「Hamel/Wirsum 正文无此物」** |
| A7 | Newton 线性步 \(F'\Delta x=-F\) | Eq. **(2.21)** OCR 清晰：`F'(x0) Ax = -F(x0)` | **确认**（原笔记「待校对」可降级） |
| A8 | 阻尼 = 在满 Newton 步上做 **回溯试探** 直至 (2.23) | 步骤 4：求 \(j\) 使 \(\|F(x+\,\cdot\,\Delta x)\|_2\) 下降；步骤 5 更新 | **确认「有回溯」**；因子字形 OCR 残，**解读为 \(2^{-i}\)**（高置信，非 E1 字形） |
| A9 | 实用终止：迭代上限 / \(\|\Delta x\|\) / \(\|F\|\) | §2.3.3 三条 Abbruchbedingungen | **确认** |
| A10 | 「只拧 `λ` / halvings、仍用 max-merit」≠ 测 Wirsum 同构 | 原文标量是 \(\|F\|_2\)，不是 max(split,energy) | **确认**：§4.1 阻尼旋钮 A/B **无效**与「未测对 Wirsum 准则」不矛盾 |
| A11 | fastox / syngas-guard / LSQ preinner / absorb 非 Wirsum | §2.3.3 未见 | **确认：工程桥** |
| A12 | 主题是污泥 **燃烧** cell 模型，非气化动力学 | 书名/前言 | **确认**；求解语法仍可被 Hamel 气化复用 |

#### B. handoff 中曾指望 Wirsum、但全文 **给不出** 的点

| # | 曾指望 / 提问 | Wirsum 核查结果 |
|---|---------------|-----------------|
| B1 | 阻尼如何处理 **T 与组成联合步**（§5.30–5.34） | **无**：只有整体 \(\|F\|_2\)；无按变量类（T/气/固）的分策 |
| B2 | 是否规定 **composition-frozen T** 或气–能软桥 | **无** |
| B3 | 是否规定 **缩放 / 无量纲 residual** | **无**（本仓库 `F_hat` 缩放仍属工程实现细节） |
| B4 | 是否规定 Vorab / Startwert **配方**（快氧化、H₂ seed 等） | §2.3.3 **无**；Startwert 只作为「需阻尼」的动机 |
| B5 | bed0 O₂–H₂ 刚性、上段 H₂ 库存（§5.40–5.41） | **超出**本书求解段落；属物理/初值，非 Wirsum LS 细则 |
| B6 | Wozny **块消元算法**伪代码 | 只引用 Wozny；**消元细节仍缺** 1983 Habilitation |
| B7 | Hamel Eq. (2.9) 印刷形式 | 在 **Hamel** Diss.，不在 Wirsum；Wirsum 对应 (2.15)–(2.23) |

#### C. 对 handoff 主线的直接含义

1. **§4.2 energy↓/split↑ → max-merit 拒步**：是 **本仓库工程 merit** 的机制；Wirsum 若只看 \(\|F\|_2\)，同一试探可能接受或拒绝——**不能**再用「等 Wirsum 全文」回避「要不要改 merit」的工程决策，只能标清：改/不改都是工程选择。  
2. **§4.1 阻尼旋钮扫不动平台**：与 Wirsum **不冲突**——未改接受标量时，改 λ 只是同景观下的步长。  
3. **「文献对齐」最小实验（已落地默认）**：`nr_line_search_gas_phase_split_merit_thesis=False` 与 `nr_line_search_energy_merit_thesis=False`；接受准则 = **残差 RMS/‖F̂‖₂ 下降**（保留缩放）。ladder `p0/p1/p2` 同步。旧 max-merit 栈仅作 opt-in A/B。  
4. **勿宣称**：max-merit、fastox LS、syngas guard、refresh-aware accept、L4（local2/gate/hot/LSQ）=「Wirsum 同构」（上述工程桥默认关）。

**落地记录（2026-08-20）：** `reactor.py` 默认关；`scripts/_lu_harness.py` LADDER p0–p2；`getattr` 回退 `False`；`_line_search_merit` docstring → Wirsum Eq.2.23。

---

## 7. 配置与复现备忘

常用 thesis 旋钮（默认值以 `reactor.py` 为准）：

```text
extent_limiters_enabled_thesis = False
vorab_transport_x0_fast_oxidation_closure_thesis = True
vorab_transport_x0_char_oxidation_o2_closure_thesis = True  # R1 残氧；§5.48
vorab_init_solve_bottom_cell_thesis = False               # 底格 fsolve，默认关
nr_line_search_fast_oxidation_projection_thesis = False  # 工程桥，默认关
nr_line_search_syngas_collapse_guard_thesis = False      # 工程桥，默认关
# 平台桥（local2/absorb/gate/hot/LSQ）默认关；验收用 HAMEL_WIRSUM_CORE
nr_refresh_hydrodynamics_on_accepted_step_thesis = False
nr_line_search_gas_phase_split_merit_thesis = False      # Wirsum Eq.2.23
nr_line_search_energy_merit_thesis = False               # Wirsum Eq.2.23
nr_prefer_full_step_thesis = True
nr_lambda_init_thesis = 1.0
nr_inner_t_step_default_cap_K_thesis = 600.0  # CORE=40；§5.53
nr_temperature_fence_bed0_lower_margin_K_thesis = None  # CORE=350；§5.58
nr_temperature_fence_upper_bed_lower_margin_K_thesis = None  # CORE 关；§5.64 上段 350 增益很小
vorab_bed0_nr_startwert_T_K_thesis = None  # CORE=820；§5.60（880 见 §5.59）
vorab_vm_devolatilization_zone_cells_thesis = 3  # CORE 钉死；§5.61 勿改 1
vorab_a_tier_post_staged_h2o_passthrough_thesis = False  # CORE 关；§5.62 打开砸栏
vorab_a_tier_post_staged_syngas_passthrough_thesis = False  # CORE 关；§5.63 填 E3 但 o8 rms 变差
nr_outer_preinner_phase_split_*  (Phase2: bed1, max_bed≈2)
```

诊断插桩手法（历史 throwaway，未合入生产）：

- wrap `_line_search_merit` / `_line_search_trial_acceptable`  
- 读 `nr_accepted_lambda_history`、`nr_clip_history`、`nr_outer_history`  

---

## 8. 交接检查清单

- [ ] 读 `docs/CLAUDE.md` + 本文 + `docs/siegen_cell_model_lineage.md`  
- [ ] 复现一次 Phase2 o8（Wirsum ‖F‖₂ 接受默认），记录 CO/rms/y_H₂；未 fully converge 仍预期  
- [ ] 勿改 `.plan.md`；勿恢复床层 extent 限幅  
- [ ] 新数值结论写入 `data/phase2_*.json` 并在此文档追加路径  
- [ ] 若临时恢复 max-merit：显式 opt-in，勿写回默认  

---

## 9. 一句话给接手者

**验收用 HAMEL_WIRSUM_CORE。§5.60：T=820→o8≈840 K。§5.61：VM zone=3。§5.62：H₂O 透传砸栏。§5.63：产物透传填 E3 但 o8 rms 变差。§5.64：840 是 T 吸引子；o8 能量平台是 inner budget 用尽；上段 1000 K 栏勿进 CORE。勿把 VM 补进 Eq.2.7，勿用 850，勿改 k、K_bd。**
