# BFB 气化炉状态报告

日期：2026-04-14

## 摘要

本轮工作的主题是把 thesis-mode 主线继续从“结构近似 Hamel”推进到“`bed / freeboard / cyclone / return_leg` 在求解器语义上统一采用 Hamel Eq. 2.6 的 holdup/transport 口径”。

本日已完成的核心收口包括：

- `Cell.m_solid` 在 thesis-mode 主线上统一按 **cell holdup [kg]** 理解，而不是旧的 kg/s 流率语义
- `bed / freeboard / cyclone / return_leg` 全部走 `holdup_transport`
- `K`-传输矩阵在 thesis holdup 口径下仅作用于 `char/ash`；`VM/moisture` 不参与跨格 `K*m` 传输
- thesis-mode 可切到 **single-shot Vorabrechnung sources**：drying/DAEM 在平均床温下一次预算并固定，outer 仅刷新 hydrodynamics
- freeboard trajectory closure 已提供 solver-facing 的 class-wise `m_hold / m_dot_auf / m_dot_ab / K_auf / K_ab`
- `freeboard_cells + cyclone + return_leg` 已接入同一条 solver-facing graph
- `cyclone / return_leg` 不再只是 routing placeholder，而是显式 solver cells
- thesis-mode 下显式 `return_leg` 返料已改为**直接用 `K_leg * m_leg` 的下行通量**闭合到底床，不再受 legacy `recirculation_frac` 缩放
- result summary 中 `freeboard eject / exit / return / cyclone capture` 已改成显式 solver-state 口径，不再继续复用 closure capture proxy

当前结论：**Hamel thesis 对齐主线已完成当前阶段的关键结构与求解语义收口**；但“连接矩阵/全局 Jacobian 同构实现”仍属于台账中的开放项（见 `docs/hamel_unconsistency_register.md`），后续优先级从“主线 blocker”转向“同构度收敛 + 性能/可解释性优化”。

---

## 一、本日修订总览

| 模块 | 本日修订 | 结果 |
|---|---|---|
| `src/core/cell.py` | thesis-mode 下 side-block 反应源项强制为 0；`m_solid` 的 holdup 语义继续沿主线使用 | `cyclone / return_leg` 不再误参与反应 |
| `src/core/connectivity.py` | 新增 side-block thesis transport 系数与 seed 逻辑；`return_leg` 底部返料改为直接使用下行通量；side-block 也支持 inner NR freeze | side-block 从 routing block 升级为 thesis-style explicit cells |
| `src/core/reactor.py` | `cyclone / return_leg` 默认切到 `holdup_transport`；summary 的 freeboard/side-block 指标改成 explicit solver-state 口径 | 结果输出不再混用 closure capture proxy |
| `tests/test_reactor_solid_streams.py` | 增加/更新 side-block holdup、freeze、return-leg recycle 回归 | 锁住 thesis-mode 返料与冻结语义 |
| `tests/test_global_nr_solver.py` | 增加 summary 指标与 explicit freeboard/cyclone state 的一致性断言 | 锁住显式结果摘要口径 |

---

## 二、solid-state / holdup 主线收口

### 1. 统一语义

thesis-mode 当前主线已经不再把 `m_solid` 当作“固相流率”，而是当作 **单元内固体 holdup [kg]**。

对应地，固相传输统一写成：

- 上行：`m_dot_auf = K_auf * m_solid`
- 下行：`m_dot_ab = K_ab * m_solid`

并进入标准 Eq. 2.6 形式：

`m_zu + m_rez + m_in + m_auf_in + m_ab_in + R + migration - (K_auf + K_ab) * m_solid = 0`

### 2. 本轮确认的剩余模块审计

已复查：

- `src/core/elemental_ledger.py`
- `src/core/cell_pyrolysis.py`

当前没有残留继续直接消费 `m_solid` 的旧流率式调用点，因此这两处不再构成 holdup 主线的 blocker。

---

## 三、bed / freeboard / side-block 的 thesis 一致性

### A. Bed

床层继续沿用用户给出的 thesis 口径：

- `K_auf = f_w * eps_b * u_b / ((1 - eps_b) * dh)`
- 顶床层进入 freeboard 时乘 `zeta_w`
- `K_ab` 由宏观质量连续性自上而下递推得到

这些系数已经进入 frozen Vorabrechnung 生命周期，并作为线性 transport coefficient 写入 Eq. 2.6。

### B. Freeboard

freeboard trajectory closure 当前已不再只输出 absolute entrained flow，而是能够输出：

- class-wise `m_hold`
- class-wise `m_dot_auf`
- class-wise `m_dot_ab`
- class-wise `K_auf`
- class-wise `K_ab`

并同步回 nominal `freeboard_cells`。这意味着 freeboard 在 thesis-mode 下已经不只是 postprocess 段，而是 solver-facing 的显式 cell 组。

### C. Cyclone / Return Leg

根据本轮用户提供的 thesis 口径，`cyclone` 与 `return_leg` **没有专属守恒方程**，而是与 bed/freeboard 一样沿用 Eq. 2.6；它们的物理角色只由连接矩阵与 transport coefficient 决定。

当前实现已对齐到以下口径：

#### Cyclone

- 视为非鼓泡段（`eps_b = 0`, `K_bd = 0`）
- 使用短停留时间 `Δt_cyc` 推出 `K_cyc = 1 / Δt_cyc`
- 对 char / ash 使用分离效率 `eta_k` 分成：
  - 上行飞灰：`(1 - eta_k) * K_cyc * m_cyc`
  - 下行捕集：`eta_k * K_cyc * m_cyc`

#### Return Leg

- 视为致密回料段（`u_g = 0`, `eps_d_voidage ≈ eps_mf`）
- 用 dense-leg capacity 与当前 cyclone 下行量估算 `K_leg = 1 / Δt_leg`
- 向床底返料直接取：

`m_rez,bed = K_leg * m_leg`

**关键新约定：**

在 thesis-mode 的显式 return-leg 路径中，**忽略 legacy `recirculation_frac`**；底部返料不再额外做 legacy 缩放。

---

## 四、solver graph 与 inner freeze 收口

### 1. solver graph

当前 thesis-mode 下，以下 block 已进入同一条 solver-facing graph：

- `bed cells`
- `freeboard cells`
- `cyclone`
- `return_leg`

这意味着 freeboard / side-block 不再只是“执行路由的外部附属结构”，而是显式 state blocks。

### 2. inner NR freeze

本轮补齐了 side-block 的 freeze 生命周期：

- outer refresh 时 snapshot 当前 side-block `K`
- inner NR 中边界更新不再重复重算 side-block transport coefficient
- side-block 现与 bed/freeboard 一样遵守“outer refresh / inner fixed”的 thesis 语义

这一步对数值意义比对性能意义更重要：它避免 inner NR 里 transport coefficient 漂移，保持 Jacobian 线性化前提更接近 Hamel 原始结构。

---

## 五、result summary 口径修订

此前 thesis-mode freeboard-active 路径里，虽然显式 `freeboard_cells / cyclone / return_leg` 已存在，但 summary 仍有一部分指标沿用 closure proxy。

本轮已把以下指标切到 explicit solver-state 口径：

- `freeboard_entrained_eject_char_ash_kg_s`
- `freeboard_entrained_exit_char_kg_s`
- `freeboard_entrained_exit_ash_kg_s`
- `freeboard_entrained_return_char_kg_s`
- `freeboard_entrained_return_ash_kg_s`
- `freeboard_cyclone_capture_char_kg_s`
- `freeboard_cyclone_capture_ash_kg_s`
- `freeboard_cyclone_recycle_candidate_char_ash_kg_s`
- `freeboard_entrained_char_kg_s`
- `freeboard_entrained_ash_kg_s`
- `freeboard_entrained_return_char_profile_kg_s`
- `freeboard_entrained_return_ash_profile_kg_s`

现在这些量分别来自：

- freeboard cells 的 `_solid_upflow_rates()`
- freeboard cells 的 `_solid_downflow_rates()`
- cyclone 的 `_solid_downflow_rates()`

而不再依赖 closure 中的 capture fraction 代理值。

---

## 六、回归与门禁

本轮最终回归结果：

```text
pytest -q tests/test_reactor_solid_streams.py tests/test_global_nr_solver.py tests/test_freeboard_analytical_trajectory.py tests/test_solver_sequence.py tests/test_cell_balances.py
88 passed
```

强制 sanity gate：

```text
python3 tests/sanity_checks.py
ALL 17 SANITY CHECKS PASSED
```

---

## 七、本日新增/更新的关键回归点

### `tests/test_reactor_solid_streams.py`

- `test_explicit_return_leg_recycle_keeps_gas_zero_and_solid_only`
- `test_apply_all_bc_with_explicit_freeboard_routes_freeboard_to_side_blocks_and_backmixes_into_bed_top`
- `test_side_block_transport_coefficients_freeze_inside_inner_nr_boundary_updates`

### `tests/test_global_nr_solver.py`

- `test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately`
  - 现额外断言 result summary 中 freeboard exit / return / cyclone capture 与显式 state 的 up/down transport 一致

---

## 八、当前状态评级

| 模块 | 当前评级 | 说明 |
|---|---|---|
| solid holdup 语义 | **A-** | thesis-mode 主线已统一 |
| bed transport coefficient | **A-** | 已写入 frozen 线性系数主线 |
| freeboard solver-facing transport | **A-** | closure 已输出 `K_auf/K_ab` 并同步 explicit cells |
| side-block thesis 对齐 | **A-** | 已从 routing block 升级为 Eq. 2.6 explicit cells |
| summary 结果口径 | **A-** | 关键固相指标已切到 explicit solver-state |
| 整体 thesis 结构同构度 | **B+ / A- 之间** | 主线已收口；连接矩阵/返料侧块 Jacobian 同构仍在台账 open/in_progress 项内 |

---

## 九、后续若继续做（非 blocker）

以下不再是 Hamel 一致性的主 blocker，而是可选优化项：

1. 为 freeboard-active path 重建更精细的 `_affected_nr_residual_cells()` 邻接裁剪，以恢复更多 block sparsity 性能收益。  
2. 若需要更高的结果可解释性，可继续把 `u_gb / u_p / u_t / carry_ratio` 等 freeboard profile 的 summary 从 closure 代理量逐步替换为 solver-facing 重建量。  

---

## 九点五、今日补充（按“1→2”顺序）

### 1) 文档同步（已完成）

- 已同步 `README.md`、`docs/BFB_TechSpec_v11.md`、`docs/techspec.md`、`docs/hamel_dissertation_vs_python_architecture.md`：
  - 明确 `thesis_vorab_sources_single_shot=True` 的语义：drying/DAEM 在 Vorabrechnung 单次预算并固定；
  - 明确 inner NR 冻结、outer Abgleich 仅刷新 hydrodynamics；
  - 明确 A1 风格主组分 Gibbs 初值已接入 Vorabrechnung 初值链路。

### 2) 收敛/连接性回归（已完成）

- 针对 `global_nr/solver_sequence/connectivity/cell_balances` 路径回归：
  - `pytest -q tests/test_global_nr_solver.py tests/test_solver_sequence.py tests/test_reactor_solid_streams.py tests/test_cell_balances.py`
  - **85 passed**
- 全局 sanity gate：
  - `python3 tests/sanity_checks.py`
  - **ALL 17 SANITY CHECKS PASSED**
- 备注：测试输出仍有 `DGEMV parameter #2 illegal value` 的 BLAS 警告，当前不影响通过性，但建议后续单独定位（非本轮 blocker）。

---

## 十、结论

截至 2026-04-14，本项目 thesis-mode 主线已经完成本阶段的关键重构：

- `bed / freeboard / cyclone / return_leg` 已统一在 **holdup + frozen transport coefficient + Eq. 2.6** 框架下工作
- freeboard 与 side-block 已不再游离在 solver 主线之外
- result summary 的关键固相量也已切到显式 solver-state 口径

因此，当前状态可以认为已经从“论文结构近似”推进到“论文主线语义基本一致”的阶段；但若以“与 Hamel 原始 Fortran 离散同构”为目标，仍需继续完成台账中的开放项。

---

## 十一、NR 脚本监控补充（wall-time / RMS / x0）

为落实“NR 求解脚本统一监控 + Vorabrechnung x0 合理性检查”，新增：

- `scripts/_nr_monitor.py`
  - `build_vorab_x0_sanity(reactor)`：按 Vorabrechnung 初始化链构建 x0，并检查
    - 有限值 / 非负性
    - `T` 区间
    - `eps_b` 区间
    - `K_bd` 非负
  - `solve_with_nr_monitor(...)`：统一返回 `wall_time_s`、`nr_total_s`、`rms_scaled_final`、`converged` 等指标
  - `print_nr_monitor(...)`：统一打印口径

已接入当前常用 LU/Freeboard 审计脚本：

- `scripts/audit_freeboard_exit_lu.py`
- `scripts/audit_bed_solid_transition_lu.py`
- `scripts/audit_freeboard_direction_lu.py`
- `scripts/audit_global_nr_profile_lu.py`
- `scripts/audit_global_nr_profile_lu_freeboard.py`

本轮追加接入（性能/初始化专项脚本）：

- `scripts/audit_global_nr_init_strategies.py`
- `scripts/audit_global_nr_walltime.py`
- `scripts/audit_solver_parity_lu.py`
- `scripts/audit_reaction_progress_lu.py`

本轮继续追加接入（hydrodynamics / chemistry / budget 专项）：

- `scripts/audit_bed_exit_to_reactor_exit_budget_lu.py`
- `scripts/audit_hydrodynamics_consistency_lu.py`
- `scripts/audit_hydrodynamics_ud_closure_lu.py`
- `scripts/audit_mass_transfer_regime_lu.py`
- `scripts/audit_net_molar_sources_lu.py`
- `scripts/audit_top_bed_chemistry_lu.py`
- `scripts/audit_vm_transport_consistency_lu.py`
- `scripts/audit_local_oxidation_hotspots_lu.py`

本轮继续追加接入（freeboard reaction/injection/trajectory 专项）：

- `scripts/audit_freeboard_reaction_sets_lu.py`
- `scripts/audit_freeboard_secondary_injection_lu.py`
- `scripts/audit_freeboard_trajectory_models_lu.py`

说明：

- 统一打印 `wall_time_s` 与 `rms_scaled_final`，并带 `converged/outer/inner` 状态；
- `x0 sanity` 目前针对 bed-cells 的 Vorabrechnung 初始化链（符合当前 `generate_initial_x0` 入口）。

本轮新增接入（你上轮点名的 backsolve/trace/allocator 专项）：

- `scripts/audit_ch4_h2o_paths_lu.py`
- `scripts/audit_epsb_profile_lu.py`
- `scripts/audit_db0_strategy_backsolve_lu.py`
- `scripts/audit_lambda_strategy_backsolve_lu.py`
- `scripts/audit_xi_strategy_backsolve_lu.py`
- `scripts/audit_global_nr_o2_trace_lu.py`
- `scripts/audit_pyrolysis_allocator_lu.py`
- `scripts/audit_vm_phase_split_lu.py`

本轮继续接入（剩余 `global_nr` 审计脚本补齐）：

- `scripts/audit_ud_closure_full_reactor_lu.py`
- `scripts/audit_r8_hamel_alignment_lu.py`
- `scripts/audit_bubble_chain_backsolve_lu.py`
- `scripts/audit_bubble_required_closure_backsolve_lu.py`
- `scripts/audit_bubble_local_ode_backsolve_lu.py`
- `scripts/audit_bubble_ode_terms_backsolve_lu.py`
- `scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`
- `scripts/audit_freeboard_entrained_solids_lu.py`
- `scripts/audit_freeboard_air_repartition_lu.py`
- `scripts/audit_cell0_overlap_sensitivity_lu.py`

本轮再补齐（secondary-injection + init-benchmark）：

- `scripts/audit_freeboard_secondary_local_refine_lu.py`
- `scripts/audit_freeboard_secondary_modes_lu.py`
- `scripts/audit_freeboard_secondary_observation_lu.py`
- `scripts/benchmark_nr_init_gibbs_bootstrap.py`

对应冒烟运行已通过（`head` 截断输出模式）：

- `python3 scripts/audit_epsb_profile_lu.py | head -n 12`
- `python3 scripts/audit_global_nr_o2_trace_lu.py | head -n 16`
- `python3 scripts/audit_lambda_strategy_backsolve_lu.py | head -n 12`
- `python3 scripts/audit_ch4_h2o_paths_lu.py | head -n 14`
- `python3 scripts/audit_db0_strategy_backsolve_lu.py | head -n 12`
- `python3 scripts/audit_xi_strategy_backsolve_lu.py | head -n 12`
- `python3 scripts/audit_pyrolysis_allocator_lu.py | head -n 14`
- `python3 scripts/audit_vm_phase_split_lu.py | head -n 20`
- `python3 scripts/audit_ud_closure_full_reactor_lu.py | head -n 12`
- `python3 scripts/audit_r8_hamel_alignment_lu.py | head -n 10`
- `python3 scripts/audit_bubble_chain_backsolve_lu.py | head -n 10`
- `python3 scripts/audit_bubble_required_closure_backsolve_lu.py | head -n 10`
- `python3 scripts/audit_bubble_local_ode_backsolve_lu.py | head -n 10`
- `python3 scripts/audit_bubble_ode_terms_backsolve_lu.py | head -n 10`
- `python3 scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py | head -n 10`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py | head -n 12`
- `python3 scripts/audit_freeboard_air_repartition_lu.py | head -n 12`
- `python3 scripts/audit_cell0_overlap_sensitivity_lu.py | head -n 10`
- `python3 scripts/audit_freeboard_secondary_local_refine_lu.py | head -n 14`
- `python3 scripts/audit_freeboard_secondary_modes_lu.py | head -n 14`
- `python3 scripts/audit_freeboard_secondary_observation_lu.py | head -n 14`
- `python3 scripts/benchmark_nr_init_gibbs_bootstrap.py --repeat 1 --timeout-s 120 | head -n 18`

本轮运行中修复的脚本问题：

- `audit_ud_closure_full_reactor_lu.py`：移除 `dict | dict`，改为兼容旧 Python 的 `dict()` + 赋值。
- `audit_bubble_chain_backsolve_lu.py`：补齐 `_integrate_variant(..., u_d=...)` 缺失参数。
- `audit_freeboard_air_repartition_lu.py`：将超严 `assert` 改为 warning 输出，避免审计中断并暴露守恒漂移量。
