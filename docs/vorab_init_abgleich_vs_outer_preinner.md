# Init Abgleich 与 Outer Pre-Inner 职责划分

更新时间：2026-06-04

本文档说明 Hamel 流程图中 **Node B（PRECALC → Zellenmodell 映射）** 与 **外层 Abgleich 每轮 inner NR 前 Startwert 桥接** 在 Python 实现中的分工，避免两套逻辑职责重叠或互相破坏。

流程图对照：`docs/gasifier_model_flowchart_trilingual.mmd`（Check1 / Check2 / Node A / Node B）。

---

## 1. 一句话定位

| 机制 | 阶段 | 驱动信息 | 目标 |
|------|------|----------|------|
| **Init Vorab–Cell Abgleich** | 首次 NR 前（init 末段） | **冻结** Vorab 预算链 `cell_budgets` | PRECALC 语义 ↔ x₀ 一致（Check2-lite） |
| **Outer Pre-Inner Bridge** | 每个 outer 轮：refresh **之后**、inner NR **之前** | **当前** cell 状态 + 水力 refresh 后的残差 | 降低 inner NR 起步 split/solid/energy 刚性 |
| **Outer Check2 Bridge** | 每个 outer 轮：inner Check1 **通过** 且 **气相** Check2 **失败** 后 | 同 Init，复用冻结 `cell_budgets` | 轻量 Node-B 修正，**非**完整 Node-A PRECALC |

三者不是同义重复：Init Abgleich 面向 **Vorab 预算一致性**；Pre-Inner 面向 **NR 方程残差可解性**；Outer Check2 Bridge 是 outer 环上 **按需** 的 Vorab 窄域回退。

---

## 2. 时间线（默认 thesis A 档）

### 2.1 Init（`run_init_and_precalc_for_global_nr`）

```
T 剖面 / build_vorabrechnung_cell_budgets
  → generate_initial_x0（Hamel transport x₀ 或 legacy per-cell Gibbs）
  → apply_vorabrechnung_bed_stream_mapping（bed0 zu 锚定）
  → 固相 reactive 剖面 / holdup chain
  → thesis single-shot 源项冻结
  → 产物 holdup seed（阻尼 blend）
  → [opt-in] bed01 phase-split
  → bottom primary / major / joint / zone 预投影链
  → 【Init Abgleich Bridge】← 若 Check2-lite 未通过，最多 N 遍
  → 冻结 _vorab_a_tier_cell_budgets
  → 首次 global NR
```

**Init Abgleich 入口**：`src/solvers/vorab_cell_abgleich.py::run_init_vorab_cell_abgleich_bridge`  
**调用点**：`src/workflow/steps/init_precalc_step.py`（init 预投影链之后）

### 2.2 Outer 每轮（`run_global_nr_outer_abgleich`）

```
Vorabrechnung refresh（水力 / 源项 snapshot 刷新）
  → 【Outer Pre-Inner Bridge】
  → inner global NR（Check1）
  → Check2：ΔT + 可选气相门控
  → 【Outer Check2 Bridge】← 仅气相门控失败时，单遍
  → 下一轮 refresh 或结束
```

**Pre-Inner 入口**：`src/core/connectivity.py::run_outer_preinner_start_value_bridge`  
**Outer Check2 Bridge 入口**：`src/solvers/vorab_cell_abgleich.py::run_outer_vorab_cell_abgleich_bridge`  
**编排**：`src/solvers/outer_loop.py` + `src/core/reactor.py::_solve_global_nr`

---

## 3. 职责矩阵

| 能力 | Init Abgleich | Outer Pre-Inner | Outer Check2 Bridge |
|------|:-------------:|:-----------------:|:-------------------:|
| 读 Vorab `cell_budgets` | ✅ 冻结预算 | ❌ | ✅ 冻结预算 |
| bed0 主气 holdup ↔ `N_zu` | ✅ | ❌ | ✅ |
| 产物 holdup ↔ inflow + `R_gas` | ✅ 阻尼 seed | ❌ | ✅ |
| 固相 VM/水分 ↔ 干燥/热解预算 | ✅ | ❌（zone 走 NR 残差） | ✅ |
| 轴向 O₂ 剖面（opt-in） | ✅ 仅 flag 开 | ❌ | ✅ 仅 flag 开 |
| bed major **phase-split** 预投影 | ✅ bottom + opt-in bed01 | ✅ bed0（默认） | ✅ 同 Init 窄域 |
| bottom **solid/energy** 预投影 | ✅ init 链已有 zone | ✅ **主职责** | ❌ |
| 温度 blend → `t_budget_K` | ✅ init 专用 | ❌ | ❌ |
| 全局 merit gate 回滚 | ❌ | ✅ | ❌ |
| 完整 PRECALC / Node A 回退 | ❌ | ❌ | ❌ |

**设计原则**

1. **Vorab 预算项**（zu/ab、产物源项闭合、固相 reactive 剖面）→ Init Abgleich / Outer Check2 Bridge。  
2. **NR 残差几何**（split RMS、bottom zone solid/energy）→ Outer Pre-Inner。  
3. Pre-Inner **不得**重跑 `build_vorabrechnung_cell_budgets` 或改写 single-shot 源项。  
4. Init Abgleich **不得**假设水力每 outer 已 refresh；它只在 init 末态对齐 x₀。

---

## 4. Init Abgleich 细节

### 4.1 Check2-lite 评估

函数：`assess_vorabrechnung_cell_alignment`

| 指标 | 含义 | 默认阈值（config） |
|------|------|-------------------|
| `bed0_primary_holdup_rel_max` | bed0 O₂/H₂O/N₂ holdup vs `N_zu` | `vorab_a_tier_init_abgleich_primary_rel_tol_thesis` = 0.05 |
| `temperature_max_delta_K` | 床层 T vs `t_budget_K` | 30 K |
| `product_source_holdup_rel_max` | 产物 vs inflow + `R_gas` | 0.40 |
| `o2_budget_profile_rel_max` | 床层 O₂ vs `o2_remaining_mol_s` | 0.85 |
| `solid_reactive_budget_rel_max` | VM/水分 vs `x_vm` / `x_dry` 预算 | 0.25 |

全部通过 → `aligned=True`，bridge 零遍。

### 4.2 单遍修正顺序（`_run_vorab_cell_abgleich_bridge_pass`）

1. `apply_vorabrechnung_bed_stream_mapping`（`upstream_closure=False`）  
2. [opt-in] `apply_vorabrechnung_bed_axial_primary_profile_from_budget`  
3. `apply_vorabrechnung_bed_solid_reactive_profile_from_budget` + holdup chain reseed  
4. `seed_major_product_holdup_from_vorab_gas_balance`  
5. `preproject_bottom_primary_gas_state_to_exchange_closure`  
6. `preproject_bottom_major_gas_state_to_local_balance`  
7. [opt-in] `apply_vorabrechnung_bed01_phase_split_after_vorab_closure`  
8. [init only] 温度向 `t_budget_K` blend（`vorab_a_tier_init_abgleich_temperature_blend_thesis`）

诊断字段：`reactor._vorab_a_tier_init_abgleich_diag`

---

## 5. Outer Pre-Inner 细节

函数：`run_outer_preinner_start_value_bridge`

### 5.1 默认序列

1. **Gas phase-split**（`run_outer_preinner_gas_phase_split_preprojection`）  
   - Reactor 默认仅 **bed0**（`max_bed_index=0`）；HTW Phase2 staged-primary 口径为 **max=2**（含 bed2）  
   - bed1 需显式 `nr_outer_preinner_bed1_phase_split_thesis=True`  
   - bed0：`preserve_primary_inlet_totals=True`；bed1：False  

2. **Bottom zone solid/energy**（`run_bottom_zone_solid_energy_preprojection`）  
   - 默认 bed0–bed2（`nr_outer_preinner_solid_energy_max_bed_index_thesis=2`）  

3. **bed0 resplit**（split + solid/energy 均开时）  

4. **Global merit gate**：若 `merit_scaled` 变差 → 回滚 bottom 快照  

### 5.2 与 Init 链的关系

- Init 阶段已跑过的 **joint / zone / bottom split** 给出 x₀；Pre-Inner 在 **refresh 改变 eps_b、K_bd、源项 snapshot** 后，针对 **新线性化点** 再投影。  
- Pre-Inner 使用 `_measure_nr_init_scaled_metrics`（与 init 诊断同口径），但 **不读取** `cell_budgets`。

诊断字段：`reactor._outer_preinner_last_diag`；`outer_history` 中 `preinner_*` 键。

### 5.3 Check1 续跑

`strict_check1_before_refresh_thesis=True` 时，若上一轮 Check1 失败，可跳过 refresh 但仍允许 Pre-Inner（`nr_outer_preinner_on_check1_continuation_thesis`）。若 continuation 内线搜索无进展（RMS 平台、`n_iter≤1`），下一轮恢复 Vorab refresh（`src/solvers/outer_loop.py`，单测 `test_strict_check1_reallows_refresh_when_continuation_stalls`）。

---

## 6. Outer Check2 Bridge 细节

触发条件（`outer_loop.py`）：

- inner Check1 **已通过**；  
- Check2 **未**通过；  
- **且** 气相门控失败（split RMS 或 bed0 主气漂移），**非** 仅 ΔT 失配。

原因：仅 ΔT 失败时做 Vorab holdup 重映射会破坏已收敛内层状态（probe 已验证）。

配置：

- `vorab_a_tier_outer_abgleich_bridge_thesis` = True  
- `vorab_a_tier_outer_abgleich_max_passes_thesis` = 1  
- 无温度 blend（与 Init 区分）

---

## 7. 配置速查

### Init Abgleich

| Config | 默认 | 说明 |
|--------|------|------|
| `vorab_a_tier_init_abgleich_bridge_thesis` | True | 总开关 |
| `vorab_a_tier_init_abgleich_max_passes_thesis` | 2 | 最大修正遍数 |
| `vorab_a_tier_init_abgleich_*_tol_thesis` | 见 §4.1 | Check2-lite 阈值 |
| `vorab_a_tier_init_abgleich_temperature_blend_thesis` | 0.25 | init 专用 T blend |

### Outer Pre-Inner

| Config | 默认 | 说明 |
|--------|------|------|
| `nr_outer_preinner_phase_split_preprojection_thesis` | True | split 预投影 |
| `nr_outer_preinner_phase_split_max_bed_index_thesis` | 0 | 0=bed0 |
| `nr_outer_preinner_bed1_phase_split_thesis` | False | bed1 split |
| `nr_outer_preinner_solid_energy_preprojection_thesis` | True | solid/energy zone |
| `nr_outer_preinner_solid_energy_max_bed_index_thesis` | 2 | 至 bed2 |
| `nr_outer_preinner_global_merit_gate_thesis` | True | merit 回滚 |

### Outer Check2

| Config | 默认 | 说明 |
|--------|------|------|
| `nr_outer_check2_gas_phase_split_gate_thesis` | True | split RMS 门控 |
| `nr_outer_check2_max_gas_phase_split_rms_thesis` | 0.05 | 阈值 |
| `vorab_a_tier_outer_abgleich_bridge_thesis` | True | Check2 窄域 bridge |

---

## 8. 反模式（已排除或禁止）

| 错误做法 | 后果 |
|----------|------|
| Pre-Inner 里重算 Vorab 预算 / 改 single-shot 源项 | 破坏 thesis 语义；与 refresh 职责冲突 |
| Init Abgleich 每 outer 无条件运行 | 覆盖 NR 已接受状态；outer 收敛变差 |
| Check2 失败（仅 ΔT）时跑 Vorab bridge | holdup 被重置；probe 回归恶化 |
| Init 与 Pre-Inner 同时默认 bed1 split（不同 preserve 语义） | bed1 N₂/O₂ split 残差放大 |
| 用 Pre-Inner merit gate 替代 Init Check2-lite | Vorab 预算漂移不可见 |

---

## 9. 代码索引

| 文档概念 | 实现 |
|----------|------|
| Init Check2-lite + bridge | `src/solvers/vorab_cell_abgleich.py` |
| Init 编排 | `src/workflow/steps/init_precalc_step.py` |
| Outer Pre-Inner | `src/core/connectivity.py`（`run_outer_preinner_*`） |
| Outer 循环 + Check2 | `src/solvers/outer_loop.py` |
| Reactor 接线 | `src/core/reactor.py::_solve_global_nr` |
| 配置定义 | `src/core/reactor.py` `ReactorConfig` |
| x₀ transport + §4.3 热解 Gibbs | `src/solvers/vorabrechnung.py` |
| 流程图映射 | `src/workflow/simulation_runner.py` |

---

## 10. 维护规则

1. 新增 Vorab→cell 映射能力 → 优先放入 **Init Abgleich / Outer Check2 Bridge** 共享路径（`run_vorab_cell_abgleich_bridge`），不要塞进 Pre-Inner。  
2. 新增 NR 残差预投影 → 优先放入 **Pre-Inner**，需 merit gate。  
3. 修改任一 bridge 顺序 → 同步更新本文档 §3–§5 与相关 probe 脚本说明。  
4. 与 Hamel 原文差异处标注「论文证据 / 工程补缝」（见 `docs/hamel_engineering_bridge_notes.md`）。

---

## 11. VTT B 档验收与 Check2 分工（2026-06）

VTT 单回路（`vtt_minimal_loop.json`）**B 档 Jacobian 验收**已通过：`unexpected_blocks==0`、内层 NR 收敛。

在默认 A 档配置下 `reactor.solve()` 仍可能 `converged_fully=False`，实测门控分解如下：

| 门控 | VTT 实测 | 阻塞？ |
|------|----------|--------|
| Pre-inner `global_gate_pass` | True | 否 |
| Check2 气相 phase-split | False（RMS ≈0.63 > 0.05） | **是** |
| Check2 **ΔT**（`abs_tol_T=0.5 K`） | False（首轮 ~300 K） | **是** |
| Outer inner Check1 | False（`inner_rms≈0.63 > tol_global=0.15`） | **是** |
| Outer Check2 Bridge | 未触发（需 `inner_converged` 且气相失败） | — |

因此：**Pre-Inner 不是 VTT 外层不收敛根因**；短探针（`max_global_iter≤4`）下 Check1+Check2 双重失败是预期。full outer + `nr_inner_tol_rms_thesis=0.02` + Check1 stall recovery 可达成 `converged_fully=True`（见 `test_vtt_minimal_loop_outer_a_tier_convergence`）。详见 `docs/hamel_b_tier_thesis_extraction.md` §10。

自动化门控断言：`tests/test_global_nr_solver.py::test_vtt_minimal_loop_outer_check2_gate_diagnosis`（短探针）、`test_vtt_minimal_loop_outer_a_tier_convergence`（A 档 full outer）。

---

## 12. Hamel 对齐 x₀ 与 §4.3 热解释放（2026-06-04）

### 12.1 原文边界

| Hamel 段落 | 语义 | Python 实现 |
|------------|------|-------------|
| §2.1 Startwertwahl | 配置阶段赋**全部** bilanz 变量；可沿用上一轮 | `generate_initial_x0` + init 链 |
| Vorabrechnung | 水力 + 干燥/热解 → zu/ab 与源项 | `build_vorabrechnung_cell_budgets` |
| §4.3 / Anhang A1 | **热解气包**封闭 Gibbs（非全反应器进料） | `solve_pyrolysis_gas_gibbs_composition` |
| Zellenmodell | 主反应走动力学 | global NR 残差 |

**禁止**：把 A1 用于每格 holdup 全进料 Gibbs（legacy `major_gibbs_per_cell`）——会在 ~1105 K 产生 CO/CO₂ 分岔，导致 bed4+ `N_CO≈0` 而 upstream inflow 仍 ~15 mol/s CO。

### 12.2 x₀ holdup 模式（`vorab_x0_holdup_mode_thesis`）

| 模式 | 行为 | 默认 |
|------|------|------|
| `transport_from_vorab` | 主气 bed0=进料/预算；上游=轴向 inflow；产物=`inflow + §4.3 热解增量` | **thesis 默认** |
| `major_gibbs_per_cell` | 每格 A1 Gibbs（legacy 扩展） | opt-in（`vorab_major_gibbs_x0=True`） |
| `legacy_heuristic` | `_legacy_major_seed` 固定比例 | 无预算 fallback |

### 12.3 热解组成与产物 seed（Cell + Vorab 统一 §4.3）

**2026-06-04 更新**：Cell NR 与 Vorabrechnung 共用同一热解产物链，消除 legacy 双路径分裂。

| 模块 | 入口 | 行为 |
|------|------|------|
| Cell 干燥/热解 | `calc_drying_pyrolysis_sources` → `allocate_pyrolysis_products_hamel` | Teer/NH₃/H₂S 分出 → §4.3 Gibbs |
| Vorab 增量 | `incremental_pyrolysis_products_from_budgets` | **累积 VM 预算差分**（避免 bed0 单格大块 VM Gibbs 塌缩为纯 CO） |
| Gibbs 核心 | `solve_pyrolysis_gas_gibbs_composition` | NNLS 可行投影 + `infer_lambda_from_mole_guess` λ₀ → A1（`tol=1e-10`, `max_iter=250`） |

要点：

- Teer 分出后的 `(c_gas,h_gas,o_gas)` 对五主组分常 **不可行**（`A n=b` 无非负解）；用 `_hamel_project_pyrolysis_element_pool_for_a1` 在 `An≤b, n≥0` 下 NNLS/SLSQP 投影得到 A1 可行子预算（**非**固定优先级启发式）。
- `ReducedHamelGibbsSolver.infer_lambda_from_mole_guess`：由 NNLS 投影摩尔 guess 的对数比反推 λ 初值（Hamel A1-8）。
- 热解专用容差：`closure_rel_tol=1e-6`、`omega=0.85`；收敛后仍对**原始** Teer 后池做 cap 校验。
- A1 不收敛或超出原始池 → `RuntimeError`（无启发式回退）。
- `_incremental_pyrolysis_major_mol_s`：相邻床格 VM 增量 → 同上 Gibbs 链（须传 `T_K`, `P_Pa`）。
- Gibbs **不**直接写全量 holdup，只决定热解化学计量增量。
- `vorab_a_tier_product_holdup_source_closure_bed_cells_thesis=0`：产物 holdup seed 覆盖**全部** bed。

单测：`tests/test_cell_pyrolysis.py`（Gibbs H₂、元素不超限）、`tests/test_vorabrechnung_hamel_x0_transport.py`。

### 12.4 配置速查

| Config | 默认 | 说明 |
|--------|------|------|
| `vorab_x0_holdup_mode_thesis` | `transport_from_vorab` | x₀ holdup 策略 |
| `vorab_major_gibbs_x0` | False | per-cell Gibbs（legacy） |
| `vorab_pyrolysis_gibbs_composition_thesis` | True | §4.3 热解 Gibbs |
| `vorab_a_tier_product_holdup_source_closure_bed_cells_thesis` | 0 | 0=全 bed 产物 seed |

### 12.5 测试

`tests/test_vorabrechnung_hamel_x0_transport.py`
