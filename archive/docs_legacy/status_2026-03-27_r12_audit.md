# Status Snapshot — 2026-03-27

> **编号更正（2026-04-11）**：本页讨论 `R12`、`R6`、`R7` 时采用的是当时的 **legacy implementation labels**。按 Hamel (1999) thesis 编号应解释为：`impl R12 -> thesis R6 (H2 oxidation)`，`impl R6 -> thesis R7 (CH4 oxidation)`，`impl R7 -> thesis R9 (methane reforming)`。详见 `docs/reaction_numbering.md`。

## Completed in this cycle

- Rewired char reactions in cell reaction kernel:
  - Activated `R2`, `R3`, `R4` in `calc_reactions()`.
  - Updated char sink to include `R1+R2+R3+R4`.
- Added runtime toggle for `R12`:
  - `Cell.enable_r12` (default `True`).
  - `ReactorConfig.enable_r12` propagated to all cells.
- Added LU sensitivity script:
  - `scripts/audit_r12_toggle_lu.py` for ON/OFF comparison with Phase-1 solver settings.
- Stabilized conservation audit logic:
  - `scripts/audit_cell_conservation.py` now checks exchange-only closure and reaction elemental closure.

## Latest validation results

### 1) Reaction activation audit

- Imported and active: `R1-R12` pathways are all reachable in `calc_reactions()`.
- Spec note remains: `R12` is outside strict `R1-R11` baseline.

### 2) Single-cell conservation audit

- Case A (no reaction): exchange residual = `0.000e+00`.
- Case B (reaction-only): elemental source closure max residual = `1.388e-17`.
- Script conclusion: `BASIC STRUCTURE OK`.

### 3) LU end-to-end status (still off target)

- Baseline run remains non-converged at 20 iterations.
- Exit temperature remains too high (~1589 K vs target ~1100 K).
- Carbon conversion remains too low (~0% vs target ~95%).
- Dry gas composition still deviates strongly from validation targets.

### 4) R12 ON/OFF KPI comparison (same Phase-1 settings)

- OFF: `T_exit=1417.7 K`, `carbon_conv=0.0000`, `CO_dry=0.0011`, `CO2_dry=0.0199`, `H2_dry=0.0173`, `CH4_dry=0.7548`.
- ON : `T_exit=1626.9 K`, `carbon_conv=0.0000`, `CO_dry=0.0022`, `CO2_dry=0.0293`, `H2_dry=0.0111`, `CH4_dry=0.6565`.
- Both: non-converged at 20 iterations.

Relative error vs validation target (OFF / ON):

- `T_exit`: `28.88% / 47.90%`
- `carbon_conv`: `100.00% / 100.00%`
- `CO_dry`: `99.18% / 98.33%`
- `CO2_dry`: `81.87% / 73.38%`
- `H2_dry`: `85.56% / 90.75%`
- `CH4_dry`: `2595.62% / 2244.59%`

## Current interpretation

- Structural chemistry wiring gap (`R2/R3/R4`) has been fixed.
- Core conservation bookkeeping in isolated cell tests is internally consistent.
- Main remaining issue appears at full-reactor/solver coupling level (global behavior), not only local source-term closure.

## Solver-path notes (current)

- `global_nr_solver` uses FD Jacobian + damped line search and appears structurally consistent, but convergence remains weak under current initialization/scales for LU.
- `reactor._solve_gauss_seidel()` still includes strong operator splitting / pre-reaction pre-update; this likely dominates drift when chemistry is stiff.
- Recirculation is still treated as outer-loop coupling, not fully embedded as a simultaneously solved implicit block with all upstream/downstream constraints.
- Current symptom pattern (`carbon_conv≈0`, high `CH4_dry`) suggests global-path mismatch dominates over local stoichiometric closure.

## Implementation update (reactor-level debug started)

### A) 已实施项

- 对齐了 global NR 入口与 Phase-1 共享基线：
  - `tests/test_table2_LU_global_nr.py` 改为复用 `build_phase1_htw_lu_reactor_config()`。
  - 新增可控参数：`enable_r12`, `u0_target`。
- 新增同参对照脚本：`scripts/audit_solver_parity_lu.py`（GS vs global_nr）。
- 调整 NR 有限差分步长下限（避免 `x≈0` 时步长过小）。
- 修复返料耦合一致性：底部返料同时回注 `N_rez_d` 与 `N_rez_b`，并在底部进料初始化时清零两相返料缓存。

### B) 最新对照结果（同参，`enable_r12=False`）

- GS:
  - `T_exit=1499.8 K`, `carbon_conv=0.0000`
  - dry: `CO=0.1374`, `CO2=0.4519`, `H2=0.1223`, `CH4=0.0164`
- global_nr:
  - `T_exit=1044.2 K`, `carbon_conv=0.0000`, `rms_scaled_final=1.23e5`
  - dry: `CO=0.0086`, `CO2=0.0354`, `H2=0.1364`, `CH4=0.3153`

结论：两条路径都未达“可接受收敛+KPI”状态，但“失真模式”显著不同，说明路径/耦合问题仍是主矛盾。

### C) 最新 LU 端到端（`tests/test_table2_LU.py`）

- `converged=False`, `n_iter=20`
- `T_exit=1789.2 K`（相对偏差 59.75%）
- `carbon_conv=99.4%`（接近目标）
- 干基主组分仍全面超差：
  - `CO` 误差 94.76%
  - `CO2` 误差 56.50%
  - `H2` 误差 39.35%
  - `CH4` 误差 50.54%

结论：碳转化已可被驱动，但组分分配与温度场仍明显错误。

## Runtime optimization update (solver too slow issue)

已实施的提速改动：

- `global_nr_solver`：移除默认 Jacobian 构造过程中的大量 debug 输出（仅 `verbose=True` 时打印）。
- `reactor._solve_global_nr`：减少前置 GS warm-up（5→最多2），并将 NR 内外层迭代预算与 `max_iter` 对齐，避免隐式超预算。
- `cell.compute_vorabrechnung`：增加基于 `(T, tau)` 的缓存门控（温度变化小于阈值时复用），避免重复执行干燥/热解昂贵计算。
- `cell._calc_gas_enthalpy_flow`：增加气体焓向量小缓存（按温度键），减少反复 NASA 多项式求值。
- `cell_solver`：收紧 `least_squares` 兜底 `max_nfev` 预算，降低单格求解极端长尾。

测速结果（同一 LU 脚本）：

- `tests/test_table2_LU.py` 运行时间由此前数百秒级降至约 `151 s`（`real 151.30`）。

备注：当前仍属于“可跑但不够快”；global NR 仍受 FD-Jacobian O(n²) 成本主导，后续需继续做数值路径级优化。

## Convergence & Accuracy Plan (execution)

### Phase 1 — 固定基线与可比性（立即执行）

目标：所有结论建立在同一工况与同一参数上。

1. 固定基线：`enable_r12=False`（严格 R1-R11），`n_cells=10`，`recirculation_frac=0.1`，`heat_loss_frac=0.1`。
2. 使用同参对照脚本 `scripts/audit_solver_parity_lu.py` 作为主入口。

通过判据：

- GS/NR 输出字段齐全（`converged`, `n_iter`, `residual`, `rms_scaled_final`, KPI）。
- 同参重复运行偏差稳定（不出现数量级跳变）。

### Phase 2 — GS 路径收敛形态修正

目标：减少分裂策略引起的假收敛与组分漂移。

1. 审计并收敛化 `reactor._solve_gauss_seidel()` 中“预叠加反应 + 温度松弛”策略。
2. 引入轻量 early-stop（温度与残差联合判据），避免无效迭代拖时。

通过判据：

- `max_global_iter=20` 下运行时间进一步下降。
- `CO2_dry`、`CH4_dry` 不再出现异常振荡型漂移。

### Phase 3 — NR 路径质量提升（在可接受耗时内）

目标：提高 `global_nr` 的有效步比例，而不是盲目堆迭代。

1. 继续优化 FD-Jacobian 数值质量（步长/缩放/线搜索触发）。
2. 控制外层与内层预算总量，保证 `max_iter` 可预期。

通过判据：

- `rms_scaled_final` 持续下降。
- 相比 GS，NR 至少在 `T_exit` 或 `CO2_dry` 上具备可解释优势。

### Phase 4 — 精度拉齐（KPI 逐项消差）

目标：从“可跑”进入“接近验证值”。

1. 先温度与碳转化，再组分：`T_exit` / `carbon_conv` → `CO,CO2,H2,CH4`。
2. 固化回归表（每次提交记录误差）。

阶段目标（LU）：

- `T_exit` 相对误差 ≤ 20%（先达阶段线，再冲 15%）。
- `carbon_conv` 相对误差 ≤ 20%。
- 干基主组分中至少 2 项先进入容差线（CO/CO2/H2 ≤15%，CH4 ≤20%）。

### Next run command (single-source KPI)

- `/Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python scripts/audit_solver_parity_lu.py`

## Zone tracing update (pyrolysis / combustion / gasification)

新增诊断脚本：

- `scripts/audit_zone_profile_lu.py`

诊断内容：逐 cell 输出 `T`, `O2_out`, `O2_cons`, `m_VM`, `m_Moisture`, `m_Char`，并定位：

1. O2 耗尽位置
2. VM 释放完成位置
3. Moisture 释放完成位置

最新观察（3 外迭代快照）：

- `O2` 在上部接近耗尽（约 `cell 9` 才到 1% 阈值）；
- `VM` 与 `Moisture` 在床层内未完成释放（未达到 1% 阈值）；
- `cell 6~9` 出现 `m_VM`/`m_Moisture` 近似平台（`out/in ≈ 1`），提示热解/干燥源项在上部驱动不足。

与此对应的实现修复：

- 修复了 `calc_reactions()` 缓存命中路径下 `R_solid` 丢失 `VM/MOISTURE` 扣减的问题（缓存气源恢复时同步恢复固相扣减项），避免“气相有释放、固相不扣减”的口径不闭合。

当前结论：

- 分区追踪证据支持“上部热解/干燥驱动不足 + O2 过晚耗尽”并存；
- 下一步应在不破坏守恒的前提下，优先提升区带可分离性（底部更强燃烧消氧、床中段气化主导、上部热解收尾），再进入参数细调。

## Phase-2 execution notes (convergence/accuracy)

本轮已执行 GS 路径修正与诊断增强，结论如下：

1. **速度层面**：
  - 通过前期优化（缓存、预算、日志）已把 LU 典型运行压缩到分钟级（约 50–90s 区间）。

2. **收敛形态层面**：
  - 新增 `scripts/audit_gs_history_lu.py`，可输出每轮 `dT_max / max_res / T_exit / carbon_conv`。
  - 观察到明显“多稳态/跳变”模式：前几轮残差下降，随后可能突跳到高温/高残差分支。

3. **精度层面**：
  - 仍未进入可接受区间，尤其 `carbon_conv`、`H2`、`CH4` 误差仍大。
  - 当前主问题已从“局部守恒是否正确”转为“全局迭代轨道稳定性与分支选择”。

4. **当前最佳调试方向**：
  - 继续使用 `audit_gs_history_lu.py` 做轨道诊断；
  - 下一步优先做“失败 cell 精细重启策略 + 分支抑制（避免后期跳变）”，而非继续盲目调反应参数。
## Phase-2 incremental implementation (latest)

本轮新增/调整：

- 固体返料语义收紧（按 TechSpec 守恒口径统一）：
  - `m_solid_rez` 仅回流 `CHAR + ASH`；
  - 不再回流 `VM / MOISTURE`；
  - 理由：旋风返料在 Phase-1 视为已反应后固体，不应再次作为 fresh-feed 参与干燥/热解源项。

- `reactor._solve_gauss_seidel`：
  - 保留“失败 cell 小步显式回退”；
  - 增加“高残差 cell 入口重启重试一次”；
  - 增加“按最小残差回退到最佳迭代态（best-iterate restore）”；
  - 增加大跳变拒绝机制（抑制明显异常分支）。
- `cell_solver.solve_cell`：
  - 切换为有界 `least_squares` 单路径（移除 root + fallback 双路径），提高轨道可重复性。
- 新增历史诊断脚本：`scripts/audit_gs_history_lu.py`。

最新观测：

- GS 轨道比之前更可读（可看到前几轮残差下降，再进入分支竞争）。
- 仍未达到验证精度，尤其组分分配仍偏离明显。
- 当前结论：问题核心仍在“全局分支稳定性 + 反应器级耦合”，不是单元守恒或单个反应式是否启用。

## Single-cell solver stability update (2026-03-30)

目标：先把单 cell 求解稳定性与“优化器收敛 ≠ 物理收敛”的问题拆开。

已实施：

- `src/solvers/cell_solver.py`
  - 增加方程尺度归一包装：`_residual_wrapper_scaled(...)`。
  - 保留原始 `least_squares` 作为主路径（默认行为不变）。
  - 新增可选刚性稳定化路径：`stiff_stabilization=True` 时，触发“scaled → raw”两段式重求解，并按 `rms_scaled/max_res` 选择更优解。
  - 返回新增诊断字段：`stiff_attempted`、`stiff_accepted`。

- 新增脚本：`scripts/audit_single_cell_stiff_solver.py`
  - 对同一初值比较 `stiff_stabilization=False/True`。
  - 输出 `residual`、`rms_scaled`、关键气相组分与温度。

最新观测：

- `cell 0`（最难点）默认 `rms_scaled≈4.768`，开启刚性稳定化后 `rms_scaled≈4.767`，改进极小。
- `cell 1/3` 为中等难度，默认已相对稳定（`rms_scaled≈3.3e-1`），不触发稳定化。
- 说明当前“单 cell 难解”主要仍集中在底部第一格，且局部最小值问题尚未被两段式策略根本解除。

阶段结论：

- 已建立“可开关的刚性稳定化机制 + 可复现实验脚本”，后续可继续迭代而不污染默认主路径。
- 下一步应优先引入更强的 continuation（反应源项同伦 / 分块求解顺序），专攻 `cell 0` 的高残差盆地切换。

## Bottom-zone refinement + source homotopy update (2026-03-30, later)

本轮按“先改善底部燃烧区离散，再做源项延拓”的策略实现并验证：

### 1) 非均匀轴向网格（底部加密）

- `src/core/reactor.py`
  - `ReactorConfig` 新增可选字段 `dh_list`。
  - 当提供 `dh_list` 时，`n_cells=len(dh_list)`，`H_bed=sum(dh_list)`；否则保持旧的均匀网格 `H_bed/n_cells`。
  - `Reactor._build_cells()` 改为按每个 cell 的 `dh_i` 和累计 `h_center_i` 构建几何。

- `tests/validation_case_utils.py`
  - 新增 `build_phase1_htw_lu_refined_config(case, n_fine=3, dh_fine=0.15)`。
  - 默认生成：底部 `3 × 0.15 m` + 上部 `10 × 0.555 m`，总床高保持 `6.0 m`。

### 2) 反应源项同伦（source homotopy）

- `src/core/cell.py`
  - `calc_reactions(rate_multiplier: float = 1.0)`：新增反应速率比例因子。
  - 在 O2/H2O 限速后统一缩放反应项（气相 + 固相），默认 `1.0` 与旧行为一致。

- `src/solvers/cell_solver.py`
  - 新增 `_homotopy_warmup(cell, n_steps=5)`：将反应源项从小到大逐步提升（0→1），每步做小预算 `least_squares`，并把结果作为下一步初值。
  - 在 `solve_cell(..., stiff_stabilization=True)` 时先执行该 warm-up，再进入主求解。

### 3) 对比脚本与结果

- 新增脚本：`scripts/audit_refined_mesh_lu.py`
  - 对比 baseline（10×0.6m）与 refined（3×0.15m + 10×0.555m）。

关键结果（本轮运行）：

- Baseline（标准审计脚本）：
  - `T_exit=1208.8 K`（通过温度容差）
  - `carbon_conv=26.6%`（未通过）
  - 干基 `CO/CO2/H2/CH4` 全部未通过

- Refined mesh：
  - `T_exit=1209.0 K`（通过）
  - `carbon_conv=44.5%`（较 baseline 明显提升，但仍未达标）
  - 干基：`CO2`、`H2` 已进入容差；`CO`、`CH4` 仍偏差大

- 单格 homotopy 效果（cell0）：
  - 同一状态下 `rms_scaled` 由约 `0.713` 降到约 `0.294`（改善明显）

### 4) 阶段结论

- 底部加密 + 源项同伦组合**有效改善了可解性与部分 KPI**（尤其 `H2`、`CO2`）。
- 当前主要短板已收敛到：
  1) `CO` 偏低，
  2) `CH4` 近零（重整/氧化过强或区带分布仍不合理），
  3) 碳转化仍低于目标（44.5% vs 95%）。

下一步建议：

- 在 refined 网格上做“分区选择性同伦”：仅对 cell0–cell2 强化 continuation，减少上部 cell 过度氧化副作用。
- 在不改守恒结构前提下，针对 `CO/CH4` 做有限范围动力学敏感性扫描（R5/R6/R7/R4）。

## Temperature-profile tuning investigation (2026-03-30)

已新增脚本：

- `scripts/audit_temperature_profile_tuning.py`

覆盖场景：

1. baseline uniform（10×0.6m）
2. refined default（3×0.15m + 10×0.555m）
3. refined + heat_loss=0.12
4. refined + gas_inlet_dense_frac=0.35
5. refined + dense=0.35 + heat_loss=0.12
6. #5 + R5×0.3 + R4×3
7. #6 + R6×0.5

关键温度剖面结论：

- `refined_default`：温度曲线平滑，`T_exit≈1200K`，但 `Xc` 偏低（~0.025）。
- `dense=0.35`：温度仍可控（`T_exit≈1200K`），`Xc` 显著上升（~0.784），但 O2 仅在上部接近耗尽（cell 12）。
- `dense=0.35 + heat_loss=0.12`：出现强烈多峰/跳变（最高约 `1825K`），判定为不稳定分支，不建议继续沿此路径调参。
- 在不稳定组合上叠加 `R5↓, R4↑` 后：
  - 可把温度拉回可控区（`T_exit≈1224K` 或 `1200K`），
  - 但 `Xc` 在不同子组合下大幅波动（约 0.088 ~ 1.000），说明参数耦合敏感，存在“温度可稳但转化/组分漂移”的风险。

当前建议的温度调优窗口：

- 以 `refined + dense≈0.35` 作为主线（温度剖面稳定且有较高转化潜力）；
- 暂不叠加更大热损（`heat_loss=0.12`）作为主配置，避免诱发高温分支；
- 后续在该窗口内小步扫描 `R5/R6/R4`，优先约束 `CO` 与 `CH4`，并同步盯 `T_profile` 单调性/平滑性（避免出现局部尖峰 > 1400K）。

## Focused temperature-window scan (2026-03-30, continued)

新增脚本：

- `scripts/audit_temperature_profile_window.py`

扫描窗口：

- `dense_frac ∈ {0.30, 0.35, 0.40}`
- `R5_scale ∈ {0.5, 0.7, 1.0}`
- `R6_scale ∈ {0.5, 0.7, 1.0}`
- `R4_scale ∈ {1.0, 1.5, 2.0}`
- 固定 `refined mesh + heat_loss=0.10`

结论摘要：

- 以“温度剖面平滑 + 无过热尖峰 + KPI 方向性”为综合准则，较优候选集中在：
  - `dense≈0.40, R5≈0.5, R6≈1.0, R4≈1.5`
  - `dense≈0.35, R5≈0.5, R6≈1.0, R4≈1.0`
  - `dense≈0.40, R5≈0.5, R6≈0.5, R4≈1.5`

- 这些候选的共同特征：
  - `T_exit` 维持在约 `1200–1240 K`
  - `T_peak` 大多控制在 `1250–1280 K` 左右（明显优于此前 1400–1800 K 跳峰分支）
  - `carbon_conv` 可升高，但 `CO/CH4` 仍偏低，说明当前主矛盾已从“温度剖面是否稳定”转为“气相分配路径错误”

### R7（蒸汽重整）敏感性补充

在候选窗口上进一步测试 `R7_scale`：

- `R7_scale=1.0`：可出现较高 `CH4`（~0.078），但 `Xc` 明显下降（~0.265），且并不改善主组分整体匹配。
- `R7_scale=0.5~0.3`：温度较稳，但 `CH4` 仍很低（~0.001）。
- `R7_scale=0.1`：`CO/H2/Xc` 被强烈推高，但 `CO2` 过高、`T_exit` 掉到 ~1088 K，偏离另一侧分支。

阶段判断：

- 单独通过 `R7` 调整无法同时修复 `CH4` 与整体组分误差；
- 当前更像是 **区带反应竞争（R5/R6/R7/R4）+ O2 耗尽位置过晚** 的组合问题，而非单一重整速率常数问题；
- 下一步应优先继续盯：
  1. O2 前移耗尽（底部 2–5 cell），
  2. 保持 `T_peak < 1300–1350 K`，
  3. 再看 `CH4` 是否在中上部被保留下来。

## Oxygen reaction trace + peak-cause first analysis (2026-03-31)

按“先查原因，不先硬限温峰”的原则，新增：

- `scripts/audit_oxygen_reaction_trace.py`

该脚本逐 cell 输出：

- `O2_supply`（两相入口+回流）
- `O2_out`
- `limit_factor_o2`
- `bubble O2 demand`（R5b + 1.5*R6b）
- `dense/char/tar O2 demand`（R5d/R6d/R12/R1/R9/R10）
- 温峰邻域（peak cell ±1）反应分解

### 本轮观察（稳定候选与参考分支）

1. **温峰附近 bubble 相氧相关反应占比偏高**：
  - 多个温峰邻域 cell 出现 bubble O2 需求占有效 O2 需求较高（甚至极高）的情况。

2. **当前 O2 limiter 覆盖不完整（结构性现象）**：
  - 现实现中 limiter 主要作用于 dense/char/tar 聚合组；
  - bubble `R5b/R6b` 未被同一 O2 limiter 直接约束，可能在局部放大放热峰。

3. **温峰并非单一机理**：
  - 某些温峰点由 `R1`（炭燃烧）贡献明显；
  - 某些温峰点由 bubble 氧化贡献更突出；
  - 即“char 燃烧 + bubble 氧化 + 供氧分配”共同塑形。

4. **阶段判断（与用户策略一致）**：
  - 目前不应先做“强硬限温峰”约束；
  - 应先把温峰成因分解清楚，再决定限速/分配策略。

### 下一步（原因导向）

- 在同一工况下继续做 peak-neighborhood 追踪：
  1. 先定位温峰由 `R1` 主导还是 bubble `R5/R6` 主导；
  2. 再看 O2 供给分配（dense/bubble）是否与主导机理匹配；
  3. 最后才评估是否需要把 bubble 氧化纳入统一 O2 限制框架。

### 补充：VM / 挥发分释放链路（本轮重点）

根据用户建议，本轮已把“VM 释放 → 快速均相氧化”的链路显式加入氧追踪：

- 新增诊断量：`VM_O2eq`
  - 用热解缓存气源中的 `CO / CH4 / H2 / H2S / TAR` 换算成“若被快速均相氧化所需的 O2 当量”。

- 新增诊断量：`homogeneous_O2`
  - 把 `R5 / R6 / R9 / R10 / R12`（含 bubble + dense + tar oxidation）作为一个统一的“快速均相氧化桶”。

最新观察：

1. 在温峰邻域，`homogeneous_O2` 往往占有效 O2 需求的主导部分；
2. `VM_O2eq` 在温峰邻域并不为零，说明刚释放出的挥发分确实为快速均相氧化提供了直接燃料；
3. 因此，用户判断是正确的：**温峰很大概率与挥发分释放后被快速均相氧化直接耦合有关**，而不是只看 char combustion 就够。

当前工作结论：

- 后续关于 O2 burn-out 的分析，应把 `R5/R6/R10/R12` 视为第一优先级，`R1` 作为并行解释项；
- 若后面需要改模型，应优先考虑“氧分配与快速均相氧化耦合”的改法，而不是先加一个生硬的温峰上限。

### 单 cell 优先诊断（2026-03-31）

根据用户要求，本轮新增：

- `scripts/audit_single_cell_oxygen_convergence.py`

执行策略：

- 先逐个传播上游入口到目标 cell；
- 再对目标 cell 单独调用 `solve_cell()`；
- 只有当单 cell 的 `rms_scaled / physically_converged` 过关后，才认为它适合继续用于整炉传播分析。

最新结果（稳定候选分支，下部 cells 0–5）：

- 所有 tested single-cell 都是 `optimizer_success=True`，但 `physically_converged=False`；
- 典型 `rms_scaled` 约在 `0.15 ~ 0.48`，仍高于当前物理收敛线；
- 说明：**下部单 cell 自身尚未完全收敛，当前很多整炉行为不能直接归咎于传播路径本身。**

更重要的是：

- 单 cell 氧耗主导项明显被 `R10 tar oxidation` 占据；
- 在下部单 cell 中，`R10` 对 O2 的需求量级远高于其它 O2 相关项；
- 同时单 cell 中确实存在非零 tar 存量与非零 `VM_O2eq`，支持“VM 释放后快速进入 tar/均相氧化链路”的判断。

阶段判断：

- 当前最该先修的是 **单 cell 级别的快速均相氧化口径/量级问题**；
- 在单 cell 未做到物理收敛前，不宜把主要注意力放在整炉传播策略上。
