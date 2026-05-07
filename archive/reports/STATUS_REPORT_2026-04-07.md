# BFB 气化炉状态补充报告

日期：2026-04-07

> **编号更正（2026-04-11）**：本报告早期段落大量使用的是 **legacy implementation labels**，不是 Hamel (1999) thesis 原始编号。后续阅读时请按以下映射解释：`impl R3 -> thesis R4`，`impl R4 -> thesis R3`，`impl R6 -> thesis R7`，`impl R7 -> thesis R9`，`impl R12 -> thesis R6`，`impl R9 -> H2S placeholder（非 thesis 主编号）`。权威总表见 `docs/reaction_numbering.md`。

## 本轮结论

- 已重新读取现有状态文档并核对当前实现，确认仓库主线仍以 `stable_window_candidate` 为稳定基线。
- 当前稳定基线再次验证通过：
  - `python3 -m pytest tests/test_cell_kinetics.py tests/test_phase_gate_20.py tests/test_phase_gate_30.py tests/test_table2_LU.py -q`
  - 结果：`10 passed, 1 xfailed`
- `python3 scripts/audit_oxygen_reaction_trace.py` 复核后，主问题仍集中在：
  - `bubble` 相 O2 需求在多个 cell 中异常主导；
  - `cell 2-4` 的 `R5b` 与 `R10b` 无限速率量级远高于局部 O2 预算；
  - 现有分相 O2 limiter 只是把问题显性化，没有解释其来源。

## 本轮新增证据

新增脚本：

- `python3 scripts/audit_bubble_timescales.py`
- `python3 scripts/audit_bubble_rate_components.py`
- `python3 scripts/audit_bubble_species_budget.py`

该脚本在稳定基线下输出 `bubble` 相的：

- `tau_b = dh / u_b`
- `tau_ex = 1 / K_bd`
- `tau_r5 = C_O2 / r5b_vol`
- `F_hydro = V_b / tau_b`
- `F_state = N_b_total * R * T / P`
- `ratio = F_state / F_hydro`

稳定基线下的关键观察：

- `cell 2-4` 中，`tau_r5 ≈ 1.7e-4 ~ 3.2e-4 s`
- 同区段 `tau_b ≈ 0.18 ~ 0.20 s`
- 同区段 `tau_ex ≈ 0.97 ~ 1.23 s`
- 因此 `tau_r5 << tau_b << tau_ex`
- `F_state / F_hydro` 在多数 cell 约为 `0.44 ~ 0.49`

解释：

- 以当前口径，`bubble R5` 在数值上已经处于“相对流动/交换几乎瞬时”的极限区；
- 同时，状态变量 `N_b` 反推的 bubble 体积流持续小于 hydrodynamics 给出的 bubble 体积吞吐。

### `R5b` / `R10b` 主导权切换的直接原因

`audit_bubble_rate_components.py` 在稳定基线下对 `cell 0-4` 输出了：

- `C_CO`, `C_O2`, `C_tar`
- `Qp/Keq`, `drv5`
- `k5`, `k10`
- `phi5 = C_CO * sqrt(C_O2)`
- `phi10 = sqrt(C_tar) * C_O2`
- `r5v`, `r10v`

关键结论：

- `cell1` 中：
  - `C_CO ≈ 0`
  - `Qp/Keq ≈ 2.83e3`
  - 因此 `drv5 = 0`
  - 所以 `R5b = 0` 不是因为 `R10b` 异常增强，而是 `R5b` 被 **CO 缺失 + Gibbs clamp** 一起关掉。
- `cell1` 中 `C_tar` 仍非零，因此 `R10b` 还能保持有限值。
- 到 `cell2-4`：
  - `drv5` 恢复到 `~1`
  - `phi5` 与 `phi10` 已是同量级
  - 但 `k5/k10` 已约为 `7-10x`
  - 所以 `R5b` 会自然压过 `R10b`

一句话总结：

- `cell1` 是 “`R5b` 被关掉，所以 `R10b` 显得主导”；
- `cell2-4` 是 “`R5b` 恢复后，由于 `k5` 明显更大，因此迅速接管主导权”。

### 更关键的新发现：整炉接受态 vs 单 cell 收敛态存在明显分支差距

`audit_bubble_species_budget.py` 在稳定基线下比较了：

- 整炉 GS 接受态（full-reactor accepted state）
- 在相同上游条件下重新做的单 cell 收敛态（single-cell resolved state）

结果：

- `cell1`
  - full-reactor accepted:
    - `T ≈ 1309.5 K`
    - `Nb(CO) ≈ 0.0000`
    - `Nb(O2) ≈ 1.2577`
    - 局部 `rms ≈ 1.213e-01`
  - single-cell resolved:
    - `T ≈ 1265.4 K`
    - `Nb(CO) ≈ 5.6612`
    - `Nb(O2) ≈ 0.0001`
    - `rms ≈ 1.268e-03`

- `cell2`
  - full-reactor accepted:
    - `T ≈ 1294.8 K`
    - `Nb(CO) ≈ 0.3777`
    - `Nb(O2) ≈ 1.0804`
    - 局部 `rms ≈ 9.950e-02`
  - single-cell resolved:
    - `T ≈ 1229.2 K`
    - `Nb(CO) ≈ 9.1889`
    - `Nb(O2) ≈ 0.0001`
    - `rms ≈ 2.963e-03`

这说明：

- `cell1-2` 的异常不只是“本地动力学表达式太猛”；
- 更像是 **整炉 GS 接受态停在了偏氧化、低 CO 的局部分支**；
- 而在相同上游条件下，单 cell 本身存在一个 **低 O2、高 CO、低残差** 的可收敛分支。

### 单 cell 物种预算含义

在单 cell 收敛态下：

- `cell1` 的 bubble `CO` 预算为：
  - inflow `≈ 0.6445`
  - reaction `≈ -0.2945`
  - exchange_to_bubble `≈ +5.3207`
  - outflow `≈ 5.6612`

- `cell1` 的 bubble `O2` 预算为：
  - inflow `≈ 0.5148`
  - reaction `≈ -0.5426`
  - exchange_to_bubble `≈ +0.0305`
  - outflow `≈ 0.0001`

解释：

- 在局部收敛支上，`cell1` 的 bubble `CO` 并不是被 `R5b` 大量吃掉；
- 相反，它主要是通过 `exchange` 从 dense 相补入，然后作为 bubble 出口继续向上传播；
- `O2` 才是在该 cell 内被基本耗尽。

因此，当前最值得怀疑的是：

- GS 扫描/接受逻辑在 `cell1-2` 附近保留了偏氧化支；
- 而不是简单的 `R5b` 或 `R10b` 常数本身就必须立刻修改。

## 已证伪路径

我做过一个实验性分支：

- 用 `N_b` 反推的实际 bubble 体积流去缩放 `bubble` 均相反应有效体积；
- 并同步更新了相关审计脚本口径。

结果：

- 该改动会让稳定候选的温峰更高，`oxygen trace` 表现更差；
- 因此本轮已将这条实验分支**完整回退**，没有保留在当前实现中。

结论：

- “简单的 bubble 体积吞吐一致性缩放”不是正确修正方向；
- 现阶段不应直接把 throughput mismatch 当作根因补丁写入模型。

## GS lower-pair 分支修正试验

本轮还验证了一条 solver-side 路径：在 GS 整轮扫描结束后，对 `cell1-2` 做一次额外的 full local re-solve，尝试把接受态拉回单 cell 已知存在的低 `O2` / 高 `CO` 分支。

### 正证据

直接在稳定基线最终态上手动调用 lower-pair full solve 时，确实能把局部状态拉向更合理的支路：

- `cell1`
  - 从 `T ≈ 1309.5 K, Nb(CO) ≈ 0.0000, Nb(O2) ≈ 1.2577`
  - 拉到 `T ≈ 1302.1 K, Nb(CO) ≈ 3.4965, Nb(O2) ≈ 0.0373`
  - 局部 `rms` 从 `≈ 1.213e-01` 降到 `≈ 8.750e-02`

- `cell2`
  - 从 `T ≈ 1294.8 K, Nb(CO) ≈ 0.3777, Nb(O2) ≈ 1.0804`
  - 拉到 `T ≈ 1244.5 K, Nb(CO) ≈ 7.2234, Nb(O2) ≈ 0.0003`
  - 局部 `rms` 到 `≈ 1.0e-04`

这再次说明：

- solver 并不是“根本到不了”那条局部分支；
- 问题在于它一旦嵌回整炉 GS 外循环，分支选择会重新偏到别处。

### 负证据

把这条 lower-pair full solve 真正嵌回 `src/core/reactor.py` 的 GS 主流程后，出现了两种都不能接受的结果：

1. **不加额外 profile 约束时**
   - lower-pair 修正确实会被接受；
   - 但稳定候选会变成 `Texit=1139.7 K, Tpeak=1326.5 K@cell1`；
   - 虽然 `cell1-2` 的 `bubble O2` 被显著压低，`CO` 也回升，但温峰比基线 `1309.5 K` 更高，`oxygen trace` 口径更差。

2. **给 `cell1` 加 profile-aware temperature cap 后**
   - lower-pair 修正又退化成“不被接受”；
   - 整炉行为基本回到原稳定基线；
   - 说明这条路径目前缺少一个既能保留 cooler branch、又不破坏 GS 全局接受逻辑的稳定嵌入方式。

### 当前处置

因此，这条 solver-side lower-pair 修正路径目前也视为**已证伪/未成熟**：

- 实验代码已从当前实现中回退；
- 当前仓库重新回到原稳定基线口径；
- 本轮保留的是证据，而不是这条修正本身。

回退后重新验证：

- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU.py tests/test_cell_kinetics.py -q`
  - `8 passed, 1 xfailed`
- `python3 scripts/audit_oxygen_reaction_trace.py`
  - 重新回到 `Texit=1139.7 K`
  - `Tpeak=1309.5 K@cell1`

## 文档算法口径 vs 当前实现

本轮还重新核对了 `techspec / docs / code` 的算法层口径，结论比较明确：

- `docs/BFB_TechSpec_v11.md` 与 `docs/CLAUDE.md` 记载的**论文/目标架构**是：
  - **外层**：全局 Newton-Raphson
  - **耦合结构**：分块三对角 Jacobian
  - **单 cell / 局部子问题**：作为全局方程组中的局部残差块，而不是整炉主算法本身

- 当前**默认验证口径**并不是上面这套，而是：
  - `tests/validation_case_utils.py` 中 `PHASE1_HTW_LU_SOLVE_KWARGS["solver"] = "gauss_seidel"`
  - `src/core/reactor.py` 默认也走 `solver="gauss_seidel"`

- 当前**单 cell 求解器**也已经不是文档旧口径里写的 `fsolve`：
  - 实际实现见 `src/solvers/cell_solver.py`
  - 主体是**有界 `least_squares`**
  - 并带有 `homotopy warm-up`、scaled residual stage、temperature-seeded multistart 等稳健化机制

所以现状应理解为：

- **文档目标**：全局 NR + block-tridiagonal Jacobian
- **当前主实现**：GS outer loop + bounded `least_squares` single-cell solver
- **当前技术债**：文档、项目目标与默认验证路径之间仍有结构性偏差

## 当前 global NR 状态

仓库里已经有一条 `global_nr` 路径（`src/solvers/global_nr_solver.py`），而且从实现说明上看，确实是在朝论文架构靠拢：

- 全炉状态打包成全局向量
- 稀疏有限差分 Jacobian
- 阻尼 Newton 步 + line search
- 循环返料通过 `apply_bc_fn` 施加到全局残差前

但在当前 LU Phase-1 稳定候选口径上直接运行：

- `solver="global_nr"`
- `max_global_iter=20`
- 同样的 `R5=0.50, R6=1.00, R4=1.50`

得到的是：

- `converged=False`
- `rms_scaled_final ≈ 1.20e-02`
- `Texit ≈ 1172.5 K`
- `Tpeak ≈ 1352.9 K@cell0`

这说明：

- `global_nr` 方向**在架构上更接近论文**；
- 但它**还不是当前可直接替换默认 GS 的成熟路径**；
- 现阶段如果要做算法层推进，优先级应该是把 `global_nr` 变成可验证、可对比、可回归的主路径，而不是继续给 GS 加更多局部 branch patch。

## bottom temperature cap 试验

你提到 `Tpeak > 1300 K` 未必一定不可接受，这个担心是合理的。

我做了一个**无代码 monkeypatch 实验**：只把当前 GS 中 `cell0` full candidate ladder 的 `temperature_cap_K` 去掉，不改其它逻辑。

结果对比：

- baseline
  - `final_rms ≈ 0.1500`
  - `Texit ≈ 1139.7 K`
  - `Tpeak ≈ 1309.5 K@cell1`
  - iter1 `cell0` 走 `full`，`full_cap = 1400.0`

- no bottom cap
  - `final_rms ≈ 0.1391`
  - `Texit ≈ 1139.7 K`
  - `Tpeak ≈ 1310.2 K@cell1`
  - iter1 `cell0` 仍走 `full`，只是 `full_cap = None`

含义：

- 当前 bottom `temp cap` 更像是一个**branch steering heuristic**；
- 它不是一个已经被证实的“硬物理上限”；
- 去掉它后，温峰只小幅上升了 `~0.7 K`，但整体 `rms` 反而略有改善。

因此，后续如果继续沿 GS 路径试验，`temp cap` 不应再被默认视为不可碰的物理约束，更合理的定位是：

- 可调的 continuation / branch-selection 工具；
- 是否保留，应该由整炉回归与审计结果决定，而不是仅凭 `1300+ K` 先验。

补充验证：

- 直接把 bottom `temp cap` 从实现里拿掉后，`tests/test_table2_LU.py::test_phase1_lu_window_prefers_convergence_aware_candidate`
  会失败：
  - baseline `dense=0.30` 不再优于 `dense=0.35`
- 因此本轮**没有**把“去掉 temp cap”留在主实现里；
- 当前代码仍保留原 `cell0` continuation cap，只把“可去掉”视为后续 global NR 主线上的开放问题，而不是 GS 路径上的已落地改动。

## 本轮已落地代码变更

本轮真正保留在代码里的改动，集中在 **global NR 入口收敛**，而不是直接切换所有默认 solver：

- 在 `tests/validation_case_utils.py` 中新增了
  - `PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS = {"max_global_iter": 20, "tol_global": 1.0, "solver": "global_nr"}`

- `tests/test_table2_LU_global_nr.py`
  - 改为复用这组 shared kwargs；
  - 不再在文件内部散写另一套 NR solve 参数

- `scripts/audit_solver_parity_lu.py`
  - 改为 GS 走 `PHASE1_HTW_LU_SOLVE_KWARGS`
  - NR 走 `PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS`
  - 使 LU parity 审计开始以 shared config 为准

这意味着当前项目已经有两套明确口径：

- **稳定门控基线**：GS shared kwargs
- **论文一致性开发主线**：global NR shared kwargs

后续若要真正“切到 global_nr”，可以优先围绕这组 shared NR kwargs 做收敛与回归，而不必再从散落脚本里回收参数。

## global NR warmup 敏感性

本轮还做了一个直接针对 `_solve_global_nr()` 的敏感性实验：固定同一 LU refined config，仅改变进入 NR 前的 GS warmup 步数。

结果：

- `warmup=0`
  - `Texit ≈ 1189.5 K`
  - `Tpeak ≈ 1375.0 K@cell0`
  - `rms_scaled_final ≈ 2.66e-02`

- `warmup=1`
  - `Texit ≈ 1139.7 K`
  - `Tpeak ≈ 1309.5 K@cell1`
  - `rms_scaled_final ≈ 1.89e-02`
  - 温度剖面几乎贴回当前 GS stable branch

- `warmup=2`
  - `Texit ≈ 1172.5 K`
  - `Tpeak ≈ 1352.9 K@cell0`
  - `rms_scaled_final ≈ 1.20e-02`

这说明：

- 当前 `global_nr` 对 warm-start basin **极敏感**；
- 进入 NR 之前的 GS warmup 不是一个无害细节，而是在实质性决定它最后走向哪条温度分支；
- 若要把项目主线迁到 `global_nr`，优先级很高的一件事是把
  - GS warmup 步数
  - Vorabrechnung 初值
  - outer/inner convergence contract
  这三者固定成可复现实验口径。

## 下一步建议

下一轮优先做“来源分解”而不是“直接限速”：

1. 拆开 `bubble` 相 `R5b` 与 `R10b` 的来源，按 `cell` 输出 `CO/O2/TAR` 浓度、驱动力与速率常数。
2. 下一轮重点转向 `bubble CO` 的来源与恢复机制：
   - 为什么 `cell1` 的 `bubble CO` 被抽空；
   - 为什么到 `cell2` 又迅速抬升。
3. 优先检查 GS 接受态为何没有回到单 cell 可收敛分支：
   - 是否需要在整轮 GS 扫描后，对 `rms > 0.1` 的下部 cell 做一次本地 re-solve / branch correction；
   - 是否是 `_evaluate_current_gs_state()` 只做评估、不做纠偏，导致 `cell1` 停留在偏氧化支。
4. 在没有更多证据前，不继续尝试“按体积流/停留时间硬缩放 bubble 反应”的路径，也不优先改 `R10/R5` 常数。

## 本轮续补：global NR 显式 warmup policy

基于上面的结论，本轮没有继续沿 GS acceptance / stable-window 逻辑加补丁，而是把
`global_nr` 的初始化 contract 正式显式化。

### 已落地改动

- `src/core/reactor.py`
  - `Reactor.solve(...)` 新增 `nr_gs_warmup_steps` 显式参数
  - `_solve_global_nr()` 不再只依赖内部硬编码 warmup，而是：
    - 默认仍保留旧逻辑 `min(2, max_iter // 10)` 作为 fallback
    - 允许调用方明确指定 warmup 步数
  - `global_nr` 返回结果新增元数据：
    - `nr_gs_warmup_steps`
    - `nr_outer_max`
    - `nr_outer_iters`
    - `nr_inner_tol_rms`

- `tests/validation_case_utils.py`
  - `PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS` 更新为：
    - `{"max_global_iter": 20, "tol_global": 1.0, "solver": "global_nr", "nr_gs_warmup_steps": 1}`
  - 这意味着当前 LU 的 shared NR 开发口径已经显式锁定到 `warmup=1`

- `tests/test_table2_LU_global_nr.py`
  - 新增 slow 回归，确认 shared NR policy 会把 `nr_gs_warmup_steps=1` 写入结果

- `tests/test_lu_gs_vs_global_nr.py`
  - NR integration path 也改为显式传入 `nr_gs_warmup_steps=1`

- `scripts/audit_solver_parity_lu.py`
  - 输出中新增 `nr_gs_warmup_steps`
  - 便于后续直接审计“同一 solver 参数下的 branch / KPI 差异”

### 本轮验证

- `python3 -m py_compile src/core/reactor.py tests/validation_case_utils.py tests/test_table2_LU_global_nr.py tests/test_lu_gs_vs_global_nr.py scripts/audit_solver_parity_lu.py`
  - 通过

- `python3 -m pytest tests/test_table2_LU.py::test_phase1_lu_window_prefers_convergence_aware_candidate -q`
  - `1 passed`

- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_policy_reports_explicit_warmup -q`
  - `1 passed`

- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU.py tests/test_cell_kinetics.py -q`
  - `8 passed, 1 xfailed`

- `python3 scripts/audit_oxygen_reaction_trace.py`
  - GS 基线未变：
    - `Texit=1139.7 K`
    - `Tpeak=1309.5 K@cell1`

### 当前 shared NR 数值口径

`python3 scripts/audit_solver_parity_lu.py` 输出：

- GS baseline
  - `Texit=1139.7 K`
  - `carbon_conv=0.9070`
  - `CO_dry=0.1990`

- shared global NR (`nr_gs_warmup_steps=1`)
  - `Texit=1162.8 K`
  - `carbon_conv=0.9512`
  - `CO_dry=0.2566`
  - `rms_scaled_final=2.57e-02`
  - `residual=4.932e+05`
  - `converged=False`

解释：

- 现在的 `global_nr` 结果终于是**可复现的 solver policy 输出**，而不是“内部默认 warmup 恰好落在哪条 branch”；
- 但它仍然**不是 GS 的直接替代品**，因为 KPI 与全局残差还没有一起收敛到理想状态；
- 更重要的是，既然主线准备迁到 `global_nr`，那后续的好坏判断就不应再首先依赖 GS 的
  `stable_window_candidate / convergence-aware raw score`，
  而应转向：
  - NR 全局残差
  - 轴向温度/组分 profile 自洽性
  - 与文献 KPI 的偏差

### 更新后的下一步

如果继续沿论文算法主线推进，优先级建议调整为：

1. 固定 `global_nr` 评估口径，不再拿 GS stable-window ranking 作为首要 reject 条件。
2. 继续比较 `warmup=0/1/2` 三条分支，但重点比较的是：
   - 全局残差
   - `Texit / carbon_conv / CO-CO2-H2`
   - 温度峰值是否被上部吸热气化与传播自然拉回
3. 若 hotter `cell0` branch 在全局 profile 上自洽，则不应仅因 `Tpeak > 1300 K` 自动判坏。

## 本轮续补：global NR wall-time 瓶颈与 block-tridiag Jacobian

本轮把注意力从“NR 走哪条分支”继续推进到“NR 为什么这么慢”，并专门核对了
`global_nr` 与 block-tridiagonal Jacobian 的实现瓶颈。

### 直接证据：当前主瓶颈不是线性求解，而是 Jacobian 装配

对 shared LU NR 口径（`warmup=1`）做 `cProfile` 后，热点非常明确：

- 整体约 `95 s`
- `src/solvers/global_nr_solver.py::build_jacobian_fd`
  - 累计约 `70.7 s`
- `src/core/reactor.py::_solve_gauss_seidel`（即 NR 前的 GS warmup）
  - 约 `22.2 s`
- `splu`
  - 只有 `~0.012 s`

这说明：

- 现阶段**不是** sparse linear solve 在拖慢 NR；
- 真正昂贵的是：
  - Jacobian 列有限差分时反复重算整炉 residual
  - 以及进入 NR 前的 GS warmup / outer 空转

### 结构判断：物理上是 block-tridiag，旧实现却付了近 dense-FD 的代价

从当前 reactor 耦合关系看：

- `cell i` residual 依赖：
  - `x_i`（本格）
  - `x_{i-1}`（气相上游入口）
  - `x_{i+1}`（上部固相回落入口）
  - `x_top -> cell0`（recycle 角块）

因此 Jacobian 结构应是：

- 主块 `B_i`
- 下邻块 `A_i`
- 上邻块 `C_i`
- 再加一个顶格回到底格的 recycle corner block

也就是说，它本质上是 **block-tridiagonal + recycle corner**。

但旧版 `build_jacobian_fd()` 的做法是：

- 每扰动一列变量
- 都 `unpack_reactor + apply_bc_fn`
- 然后对**所有 cell**调用 `residuals()`

对于 LU shared config：

- `n_cells=10`
- `n_var_per_cell=27`
- `n_total=270`

于是每次 Jacobian 装配都要做：

- `270` 次列扰动
- 每次扰动重算 `10` 个 cell residual
- 即 `2700` 次 cell-level residual 调用

### 已落地优化 1：按 block-tridiag 受影响块装配 FD Jacobian

本轮在 `src/solvers/global_nr_solver.py` 中保留了：

- `dense_fd`

并新增默认策略：

- `block_tridiag_fd`

核心做法：

- 每个全局变量先映射回所属 `cell j`
- 只重算受该扰动影响的 residual 块：
  - `cell j-1`
  - `cell j`
  - `cell j+1`
  - 若 `j` 为顶格，再额外加 `cell0` recycle 块

因此 Jacobian 的非零模式保持不变，但装配时不再为远距离零块白做 residual 评估。

同时，`global_nr` 现在返回额外元数据：

- `nr_jacobian_strategy`
- `nr_timing`
- `nr_counts`

其中会记录：

- Jacobian 装配耗时
- line search residual 耗时
- residual 调用次数
- Jacobian nnz

### 已落地优化 2：裁掉 outer loop 末段空转

本轮还给 `src/core/reactor.py::_solve_global_nr()` 加了一个保守的 stall break：

- 若某次 outer solve 中
  - inner `n_iter <= 1`
  - 且 `rms_scaled_final` 相比上一 outer 没有实质改善
- 则直接停 outer loop

原因是审计显示后半段 outer 已进入：

- 第 0 个 inner NR 步就 line-search fail
- `rms` 基本不再变化

这部分不再提供新分支，只是在烧时间。

### 数值一致性验证

新增测试：

- `tests/test_global_nr_solver.py`
  - 比较 `dense_fd` 与 `block_tridiag_fd` 在初始化 LU state 上的 Jacobian
  - 最大绝对差 `< 1e-8`
  - 且 block-tridiag residual 调用数更少

说明这轮改动是：

- **同一个 FD Jacobian 的结构化装配优化**
- 不是换了不同数值模型

### 当前 wall-time 结果

`python3 scripts/audit_global_nr_walltime.py` 当前输出：

- `dense_fd`
  - `wall_s = 52.29 s`
  - `outer_iters = 3`
  - `n_iter = 9`
  - `gs_warmup = 13.77 s`
  - `inner_solve_total = 38.32 s`
  - `last-inner jacobian_build = 3.72 s`
  - `jacobian_cell_residual_calls = 2700`

- `block_tridiag_fd`
  - `wall_s = 22.90 s`
  - `outer_iters = 3`
  - `n_iter = 9`
  - `gs_warmup = 12.30 s`
  - `inner_solve_total = 10.41 s`
  - `last-inner jacobian_build = 1.09 s`
  - `jacobian_cell_residual_calls = 783`

速度收益：

- 端到端 NR wall-time：`2.28x` 加速
- Jacobian build：`3.40x` 加速
- residual cell 调用数：`2700 -> 783`

更关键的是，两条策略的结果完全一致：

- `Texit = 1162.8 K`
- `Tpeak = 1372.5 K@cell0`
- `carbon_conv = 0.9512`
- `rms_scaled_final = 2.570e-02`

### 这轮之后的瓶颈顺序

现在 `block_tridiag_fd` 落地后，shared LU NR 口径下的主要耗时顺序已经变成：

1. `GS warmup`
2. 多次 outer NR solve 总和
3. 单次 inner NR 的 Jacobian build
4. `splu`（几乎可忽略）

一句话总结：

- **block-tridiag Jacobian 装配已经证明是高收益项，并已落地**
- 下一轮若继续加速，优先级不再是 linear solver，而是：
  - 减少 warmup 成本
  - 减少 outer / inner 空转次数

### 当前建议的下一步

如果继续按“论文算法主线 + 最少迭代次数”推进，下一轮优先做：

1. 审计并重构 GS warmup：
   - 现在 `warmup=1` 物理上有用，但 wall-time 仍显著；
   - 应继续拆出“真正需要的 warmup 内容”和“只是借用了昂贵的 GS cell solve”。
2. 为 outer loop 增加更明确的 basin / stagnation 诊断：
   - 现在已有 `nr_outer_history`
   - 可以进一步据此判断是否该直接减少 outer 次数或换更轻的 refresh。
3. 暂时**不要**优先投入 full block-tridiagonal 线性求解器：
   - 因为当前 `splu` 不是瓶颈；
   - 在 Jacobian 装配和迭代次数还没压下来前，先换线性代数 backend 的收益会很小。

## 本轮续补：global NR 改为显式 init strategy，并切到 paper-style Vorabrechnung

本轮把 `global_nr` 的初始化从“隐式 GS warmup 副作用”拆成了显式策略。

代码变更：

- `src/core/reactor.py`
  - `Reactor.solve(...)` 新增 `nr_init_strategy`
  - `_solve_global_nr()` 现在先解析初始化策略：
    - `vorabrechnung`
    - `gs_warmup`
  - `nr_init_strategy`、`nr_init_s_total`、`nr_vorabrechnung_s` 会写回结果元数据
- `tests/validation_case_utils.py`
  - shared LU NR 口径从
    - `nr_gs_warmup_steps=1`
  - 改成
    - `nr_init_strategy="vorabrechnung"`
    - `nr_jacobian_strategy="block_tridiag_fd"`
- 新增 `scripts/audit_global_nr_init_strategies.py`
  - 统一对比 `vorabrechnung` / `gs_warmup(1)` / `gs_warmup(2)`

### 为什么这么改

前面已经确认：

- 论文里的 `Vorabrechnung` 是 cheap bootstrap / chemically consistent `x0`
- 不是昂贵的逐 cell `least_squares` mini-solve

所以把 shared NR 开发口径继续绑在 `GS warmup=1` 上，算法上不再合理；它更像 branch steering，而不是论文式初始化。

### init-strategy 审计结果

`python3 scripts/audit_global_nr_init_strategies.py` 当前输出：

- `vorabrechnung`
  - `wall_s = 13.41 s`
  - `init_total = 0.01 s`
  - `outer_iters = 3`
  - `n_iter = 11`
  - `rms_scaled_final = 3.397e-02`
  - `Texit = 1179.5 K`
  - `Tpeak = 1365.0 K@cell0`
  - `carbon_conv = 0.9950`

- `gs_warmup(steps=1)`
  - `wall_s = 24.33 s`
  - `init_total = 13.27 s`
  - `gs_warmup = 13.25 s`
  - `outer_iters = 3`
  - `n_iter = 9`
  - `rms_scaled_final = 2.570e-02`
  - `Texit = 1162.8 K`
  - `Tpeak = 1372.5 K@cell0`
  - `carbon_conv = 0.9512`

- `gs_warmup(steps=2)`
  - `wall_s = 22.28 s`
  - `init_total = 18.35 s`
  - `gs_warmup = 18.34 s`
  - `outer_iters = 2`
  - `n_iter = 3`
  - `rms_scaled_final = 3.123e-02`
  - `Texit = 1144.8 K`
  - `Tpeak = 1311.2 K@cell1`
  - `carbon_conv = 0.9499`

结论：

- `Vorabrechnung` 初始化本身几乎免费，`~0.01 s`
- 原先的 wall-time 大头确实是 `GS warmup`
- `warmup=1/2` 仍然明显在做 branch steering，不只是“提供 Jacobian 初值”
- 因而 shared NR 口径改成 `nr_init_strategy="vorabrechnung"` 是合理的论文一致性选择

### shared NR 新 wall-time

切到 `Vorabrechnung` shared policy 后，再跑 `python3 scripts/audit_global_nr_walltime.py`：

- `dense_fd`
  - `wall_s = 37.76 s`
  - `init_total = 0.01 s`
  - `inner_solve_total = 37.63 s`
  - `jacobian = 3.39 s`

- `block_tridiag_fd`
  - `wall_s = 17.24 s`
  - `init_total = 0.01 s`
  - `inner_solve_total = 17.10 s`
  - `jacobian = 1.09 s`

速度收益：

- 端到端 wall-time：`2.19x`
- Jacobian build：`3.11x`
- residual cell 调用：`2700 -> 783`

而两条 Jacobian 策略结果保持一致：

- `Texit = 1179.5 K`
- `Tpeak = 1365.0 K@cell0`
- `carbon_conv = 0.9950`
- `rms_scaled_final = 3.397e-02`

### 这轮之后的判断

现在可以更明确地下结论：

1. `GS warmup` 不是 `global_nr` 必需的算法部件。
2. 若按论文算法主线推进，shared NR 应该以 `Vorabrechnung` 为初始化口径。
3. 现阶段 NR 的主要问题已不是“初始化太重”，而是：
   - basin 选择
   - 全局 residual / KPI 漂移
   - inner / outer 总步数仍偏多
4. `block_tridiag_fd` 仍然值得保留，因为即便去掉重型 warmup，它也继续提供稳定的 `2x+` wall-time 收益。

### 当前建议的下一步

下一轮最该做的是：

1. 直接审计 `vorabrechnung -> global_nr` 这条 shared 路径为何落到 `Texit≈1179.5 K, Xc≈0.995` 的 hotter / over-converted branch。
2. 优先减少 `n_iter=11` 这类 inner NR 步数，而不是再回头恢复 `GS warmup`。
3. 若要加 bootstrap，也应是 cheap chemical-consistency refresh，而不是单 cell `least_squares`。

## 本轮续补：shared global NR 的轴向 profile / 出口组分审计

按当前决定，`global_nr + vorabrechnung + block_tridiag_fd` 现在视为固定开发主线；
因此这轮不再把 `cell0` hotter / higher-conversion branch 当成首要问题，而是转向：

- 轴向温度 profile
- 轴向主气相 profile
- 出口干基 syngas 组成

### 新增脚本

- `scripts/audit_global_nr_profile_lu.py`

功能：

- 固定 `PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS`
- 对照 `validation_cases.json -> CASE_HTW_WESSELING_1.outputs.axial_profiles`
- 审计：
  - `T_K`
  - `CO/CO2/H2/CH4/O2/H2O` 湿基轴向 profile
  - `CO/CO2/H2/CH4` 出口干基组成

### 当前结果

`python3 scripts/audit_global_nr_profile_lu.py` 输出：

- 求解状态
  - `init=vorabrechnung`
  - `jacobian=block_tridiag_fd`
  - `outer_iters=3`
  - `n_iter=11`
  - `rms_scaled_final=3.397e-02`
  - `nr_total_s=12.70 s`

- 出口温度
  - `model = 1179.5 K`
  - `ref = 1100.0 K`
  - `abs_err = 79.5 K`
  - `rel_err = 7.23%`

- 出口干基主气相
  - `CO`
    - `model = 0.2345`
    - `ref = 0.1300`
    - `rel_err = 80.40%`
  - `CO2`
    - `model = 0.0438`
    - `ref = 0.1100`
    - `rel_err = 60.14%`
  - `H2`
    - `model = 0.1356`
    - `ref = 0.1200`
    - `rel_err = 12.97%`
  - `CH4`
    - `model = 0.0010`
    - `ref = 0.0280`
    - `rel_err = 96.59%`

- 轴向 profile 误差摘要
  - `T_K`
    - `MAE = 203.5 K`
    - `RMSE = 271.0 K`
    - 最差点：`xi=0.00`, `+615 K`
  - `CO_mol_wet`
    - `MAE = 0.0440`
    - 最差点：`xi=1.00`, `+0.0880`
  - `CO2_mol_wet`
    - `MAE = 0.0457`
    - 最差点：`xi=0.05`, `+0.0727`
  - `H2_mol_wet`
    - `MAE = 0.0182`
    - 整体相对最好
  - `CH4_mol_wet`
    - `MAE = 0.0199`
    - 最差点：`xi=1.00`, `-0.0271`
  - `O2_mol_wet`
    - `MAE = 0.0087`
    - 量级尚可，但 `xi=0.10~0.30` 仍偏高
  - `H2O_mol_wet`
    - `MAE = 0.1262`
    - 最差点：`xi=0.05`, `-0.1877`

### 读数后的判断

这轮最重要的新结论不是“温峰高”，而是：

1. `H2O` profile 全程明显偏低。
2. `CH4` 在全床层几乎被吃光。
3. 出口 `CO` 明显偏高，同时 `CO2` 明显偏低。
4. `H2` 反而已经相对接近参考。

这组组合指向的主矛盾更像是：

- 含蒸汽/含甲烷的气相反应链太强或分配不对；
- `CO ↔ CO2` 的氧化/转化平衡仍偏向 `CO` 保留；
- 后续主线应优先查：
  - `H2O` 消耗路径
  - `CH4` 消耗路径
  - 上部 cell 到出口段的 `CO/CO2` 再分配

而不是再优先盯 `cell0` hotter branch 本身。

### 当前建议的下一步

下一轮优先做：

1. 直接审计 top-half / exit 附近的 `H2O` 与 `CH4` 反应去向。
2. 将 `CO/CO2/H2/H2O/CH4` 的 profile 偏差与具体反应源项对应起来。
3. 以 profile / exit gas 为主目标，而不是再把 `Tpeak` 当主判据。

## 本轮续补：关键反应动力学 / Keq 沿程审计

按最新方向，本轮补了一个专门追 `R4/R7/R8` 的脚本：

- `scripts/audit_reaction_progress_lu.py`

固定口径：

- `global_nr`
- `nr_init_strategy="vorabrechnung"`
- `nr_jacobian_strategy="block_tridiag_fd"`

输出内容：

- `R4` Boudouard
  - `rate_R4_effective`
  - `Q4proxy = P_CO^2 / P_CO2`
  - Langmuir-Hinshelwood 分母项 `1 + kb*PCO2 + kc*PCO`
- `R7` steam reforming
  - `K_eq`
  - `Q_p`
  - `drv = 1 - Q_p/K_eq`
  - 净速率与 forward/reverse split
- `R8` WGSR
  - `K_eq`
  - `Q_p`
  - `drv = 1 - Q_p/K_eq`
  - 净速率与 forward/reverse split

说明：

- 当前代码库**没有** `R4` 的显式热力学 `K_eq` 实现；
- 因此这轮对 `R4` 不伪造平衡判据，只追动力学量级与组成代理。

### 当前读数

`python3 scripts/audit_reaction_progress_lu.py` 表明：

1. `R8 WGSR`
   - 下部 `cell 0-3`：
     - `Q8 > K8`
     - `drv8 < 0`
     - 净速率为负
     - 即局部组成推动 **reverse WGSR**
   - 中上部 `cell 4-9`：
     - `Q8 < K8`
     - `drv8 > 0`
     - 净速率转正，并逐格增强
     - 例如 `cell9`：
       - `K8 = 7.543e-01`
       - `Q8 = 2.719e-01`
       - `drv8 = 6.40e-01`
       - `r8 = 2.469e+02`
   - 这与前面看到的 top-half 净源项一致：
     - 上部在持续吃 `H2O`
     - 生 `H2/CO2`

2. `R7 steam reforming`
   - 几乎全床层 `Q7 >> K7`
   - `drv7` 被裁到 `-1e6`
   - 净 `r7` 大多为负
   - 也就是当前局部组成下，`R7` 处在 **reverse / methanation-like** 一侧
   - 这很关键：
     - 当前 `CH4` 偏低**不是**因为上部还在强 `steam reforming`
     - 更像是前段没有把 `CH4` 建起来，后段逆向恢复也不够

3. thesis `R3 Boudouard`（implementation `R4eff`）
   - `R4eff` 自底向顶持续减小：
     - `cell0 ≈ 1.0e-03`
     - `cell9 ≈ 4.2e-05`
   - `Q4proxy = P_CO^2 / P_CO2` 自底向顶快速升高：
     - `7.0e5 -> 4.1e6`
   - LH 分母项大约在 `4.6 ~ 6.3`
   - 读数表明当前 `R4` 更像是：
     - 在高 `CO` / 低 `CO2` 环境下逐步失去驱动力空间
     - 但由于缺少显式 `K_eq` 实现，还不能把它严格定性为“近平衡”

### 这轮之后的判断

这次沿程审计把优先级进一步收窄了：

1. `H2O deficit` 的确是关键问题，因为 `R8` 在上部持续正向推进。
2. `CH4` 偏低的首要嫌疑不再是上部 `R7` 过强，而更像是前段生成不足 / 早段分配错误。
3. `CO_dry` 偏高也不应简单归因于上部继续生 `CO`：
   - 上部 `R8` 其实在往 `CO2` 方向推；
   - 更可能是前段组成已经推偏，再叠加 `H2O` 赤字放大 dry-basis `CO`。

### 当前建议的下一步

下一轮优先做：

1. 直接审计前半床层的 `CH4` 来源与去向：
   - `VM -> CH4`
   - `R6`
   - `R7`
   - `R3`
2. 直接审计 `H2O` 赤字是由哪些源项主导：
   - `R2`
   - `R7`
   - `R8`
   - tar steam reforming (`R11`)
3. 若后续需要对 `R4` 做更严格热力学审计，再单独补 `Boudouard` 的 `K_eq/Qp` 实现，而不是在主模型里先拍脑袋加修正。

## 本轮续补：CH4 / H2O 路径拆分与实现对齐检查

### 实现对齐检查

按 `bfb-scientific-audit` 口径先做了两项快速核查：

1. `grep -R "np.exp" -n src/kinetics src/thermodynamics src/core`
   - 没发现新的 Arrhenius 违规写法；
   - 反应速率主干仍走 `k_hobbs / k_standard / k_jensen_r7`；
   - 剩余 `np.exp` 主要在：
     - `R8` 的经验修正项
     - `K_eq`/Gibbs 相关稳定实现
2. `python3 tests/sanity_checks.py`
   - `11/11 PASS`

同时发现一个需要单独盯住的实现/文档不一致：

- `src/kinetics/gas_reactions.py` 中 `R8` 现在使用 `P_bar = P / 1e5`
- 但 `docs/CLAUDE.md` 和早期计划文本里写的是 `atm`
- `docs/source_units_audit.md` 又记成了 `bar`

这说明：

- `R8` 压力单位口径目前**文档未统一**
- 这轮先不直接改主模型
- 但它已经被明确记录成下一轮必须核对的 source-of-truth 问题

### 新增脚本

- `scripts/audit_ch4_h2o_paths_lu.py`

作用：

- 固定 shared `global_nr` 口径
- 逐 cell 拆分 `CH4` / `H2O` 的净贡献项
- 与 `build_reaction_sources()` 相同 limiter 口径重算

覆盖项：

- `pyrolysis`
- `R2`
- `R3`
- `R6b / R6d`
- `R7`
- `R8`
- `R11b / R11d`
- `R12b / R12d`

### 当前结论

`python3 scripts/audit_ch4_h2o_paths_lu.py` 的结果非常有帮助：

1. `CH4`
   - `cell0` 主要来自 `pyrolysis`
     - `CH4_py ≈ +0.628`
     - `CH4_R6 ≈ -0.083`
     - `CH4_R7 ≈ -0.020`
     - `CH4_net ≈ +0.529`
   - `cell1-9`
     - `pyrolysis` 基本为 `0`
     - `R3` 基本也为 `0`
     - `R7` 由于当前是逆向，反而在**生成** `CH4`
     - `R11` 也在生成少量 `CH4`
   - 这说明：
     - 当前 `CH4` 偏低的主因**不是**后段 `R7` 过强把它吃掉；
     - 更像是**前段只在 cell0 建了一次 CH4 库存，之后没有持续来源**

2. `H2O`
   - `cell0`
     - `H2O_py ≈ +0.636`
     - `H2O_R2 ≈ -12.784`
     - `H2O_net ≈ -9.661`
     - 底部最主要赤字来自 `R2`
   - `cell6-9`
     - `H2O_R8` 逐渐成为最大负项
     - `cell9` 约：
       - `H2O_R2 ≈ -1.727`
       - `H2O_R7 ≈ +1.003`（逆向）
       - `H2O_R8 ≈ -6.594`
       - `H2O_R11 ≈ -0.053`
       - `H2O_net ≈ -7.365`
   - 这说明：
     - 上部 `H2O deficit` 的主因是 **WGSR (`R8`)**
     - 底部 `H2O deficit` 的主因是 **char-steam gasification (`R2`)**

### 现在更明确的优先级

这一轮之后，主线已经可以进一步收窄成：

1. `CH4` 偏低
   - 优先审 `pyrolysis -> CH4` 分配是否过弱
   - 而不是先怀疑上部 `R7` 过强
2. `H2O` 偏低
   - 底部优先审 `R2`
   - 上部优先审 `R8`
3. `R8` 的压力单位口径
   - 必须与 Hamel 原文重新对齐
   - 这是当前最可疑的“实现未完全对齐论文算法”项之一

## 本轮续补：热解源头审计与 R8 文档/实现冲突确认

这轮继续按“实现必须尽量贴 Hamel 原文”推进，但仍然**不改主模型**，先补证据。

### 新增脚本

- `scripts/audit_pyrolysis_allocator_lu.py`

功能：

- 固定 shared `global_nr` 口径
- 逐 cell 输出：
  - `m_VM_in`
  - `m_Moist_in`
  - `x_vm`
  - `x_dry`
  - `CO/CH4/TAR1/TAR2/H2O/NH3` 热解源项

### 当前读数

`python3 scripts/audit_pyrolysis_allocator_lu.py` 显示：

- 仅 `cell0` 有非零 VM 源项：
  - `m_VM_in ≈ 3.6895e-01`
  - `x_vm ≈ 0.208`
  - `CO ≈ 1.1817`
  - `CH4 ≈ 0.6283`
  - `TAR1 ≈ 0.1480`
  - `TAR2 = 0.0000`
  - `H2O ≈ 0.6361`
- `cell1-9`
  - `m_VM_in = 0`
  - 所有 VM 气源基本为 `0`

汇总：

- `active_cells = 1`
- `total CH4 from pyrolysis = 0.6283`
- `total TAR1 = 0.1480`
- `total TAR2 = 0.0000`

### 这轮之后更明确的实现判断

1. 当前 shared NR 下，VM 释放几乎只发生在 `cell0`。
2. 当前 `allocate_pyrolysis_products_elemental()` 实际只在使用 `TAR1`：
   - `TAR2` 在热解阶段始终为零；
   - 这与项目其余部分的“双 tar 代理”口径并不完全一致。
3. 因此当前 `CH4` 偏低的一个非常现实的嫌疑是：
   - 源头只在 `cell0` 给了一次 `CH4`
   - 后续没有持续 VM 补给
   - 也没有第二 tar surrogate 参与分配

### R8 单位口径：冲突已确认

本轮再次核对本地文档与实现：

- `docs/CLAUDE.md`
  - 明写：`R8` 压力单位应为 `atm`
- `src/kinetics/gas_reactions.py`
  - 当前实现：`P_bar = P / 1e5`
- `docs/source_units_audit.md`
  - 记录成 `bar`

所以这里不是“我怀疑”，而是已经确认存在：

- **文档与实现不一致**
- 且这个不一致正落在当前最敏感的 `R8/WGSR` 链上

### 当前建议的下一步

下一轮优先级现在可以非常明确地定成：

1. 核对 Hamel 原文 `R8 Eq.5.45` 的压力单位与完整速率式；
2. 审 `allocate_pyrolysis_products_elemental()` 是否应该：
   - 允许 `TAR2` 参与热解源分配；
   - 或至少不要把 VM 释放几乎全部压缩在 `cell0`；
3. 若要真正修出口 `CH4`，优先查热解源头，而不是先调上部 `R7`。

### 补充快照：top-half / exit 段净气相源项

为避免误判“出口 `CO` 偏高是不是上部还在继续生 `CO`”，我又直接看了一次当前 shared NR 解上的
`R_gas_d + R_gas_b` 净源项快照。

上部几格（`cell 7-9`）结果：

- `cell 7`
  - `R_CO = -2.473`
  - `R_CO2 = +5.154`
  - `R_H2 = +5.986`
  - `R_H2O = -7.355`
- `cell 8`
  - `R_CO = -5.374`
  - `R_CO2 = +6.520`
  - `R_H2 = +5.468`
  - `R_H2O = -7.375`
- `cell 9`
  - `R_CO = -5.877`
  - `R_CO2 = +6.709`
  - `R_H2 = +5.383`
  - `R_H2O = -7.365`

这个读数很关键：

- 上半床层到出口段**并没有**继续净生成 `CO`
- 相反，它在：
  - 持续消耗 `CO`
  - 生成 `CO2`
  - 生成 `H2`
  - 大量消耗 `H2O`

因此当前 shared NR 的出口偏差更准确地说是：

1. `H2O` profile 赤字很重，且这个赤字会一直延续到出口；
2. 出口 `CO_dry` 偏高不应简单解读为“顶部还在强生 CO”，而更可能是：
   - 前段已建立的组成偏差
   - 再叠加 `H2O` 偏低带来的 dry-basis 放大效应

这进一步支持下一轮优先级：

- 先审 `H2O` deficit 的来源；
- 再审 `CH4` 为什么始终起不来；
- 然后回头看前半床层如何把 `CO/CO2` 初始分配推偏。

---

## 2026-04-08 补充：R8 source-of-truth 与 VM transport consistency

这轮继续按“先核实现口径，再决定是否动主模型”推进，仍然**不改 solver / kinetics**，只补证据。

### 新增脚本

- `scripts/audit_r8_hamel_alignment_lu.py`
- `scripts/audit_vm_transport_consistency_lu.py`

并再次跑了：

- `python3 tests/sanity_checks.py`

结果仍是 `11/11 PASS`。

### R8：当前主敏感项不是 `atm/bar`，而是 `a_R8`

`python3 scripts/audit_r8_hamel_alignment_lu.py` 在当前 shared `global_nr` 解上比较了 4 种 R8 口径：

1. 当前代码：`P_bar + a_R8=1.0`
2. 仅改催化因子：`P_bar + a_R8=0.02`
3. 仅改压力基准：`P_atm + a_R8=1.0`
4. 两者都改：`P_atm + a_R8=0.02`

积分到 dense-phase cell source 的汇总结果：

- `current (bar, a=1.0) ≈ 3.9173e+01`
- `bar, a=0.02 ≈ 7.8346e-01`
- `atm, a=1.0 ≈ 3.9131e+01`
- `atm, a=0.02 ≈ 7.8262e-01`

相对当前实现的比值：

- `a=0.02 / current = 2.0000e-02`
- `atm / current = 9.9893e-01`
- `atm+a=0.02 / current = 1.9979e-02`

这轮得到的判断非常明确：

1. 在 LU / `P=2.5 MPa`（约 `25 bar`）工况下，`atm` vs `bar` 的差异只有 `~0.1%` 量级；
2. 因此 `atm/bar` 更像是**文档/实现一致性问题**，不是当前 WGSR 偏强的主因；
3. 当前真正可能把 `R8` 放大到过强的，是 `src/kinetics/gas_reactions.py` 里：
   - `R8_a_R8 = 1.0`
   - 但本地 Hamel 提取记录多处写的是煤灰/焦场景 `a_R8 = 0.02`

所以后续若继续对齐 Hamel，`R8_a_R8` 的 source-of-truth 优先级已经高于 `atm/bar`。

### shared NR 下已确认存在 `orphan VM`

`python3 scripts/audit_vm_transport_consistency_lu.py` 现在把每个 cell 的：

- `state_VM`
- `m_solid_zu[VM]`
- `m_solid_in[VM]`
- `pyro_VM_in = m_solid_zu[VM] + m_solid_in[VM]`
- `vm_release`

放到了同一张表里。

核心结果：

- `cell0`
  - `state_VM ≈ 2.4836e-01`
  - `pyro_VM_in ≈ 3.6895e-01`
  - `vm_release ≈ 7.6828e-02`
- `cell1-9`
  - `state_VM` 仍为非零，约 `2.66e-02 ~ 7.11e-02`
  - 但 `pyro_VM_in = 0`
  - `vm_release = 0`

脚本汇总：

- `active pyrolysis cells = [0]`
- `orphan VM cells = [1,2,3,4,5,6,7,8,9]`

这说明当前 shared NR 下已经不是简单的“VM 只在 cell0 释放”而已，而是更具体的结构问题：

1. `_propagated_solid_stream()` 明确只传播 `CHAR + ASH`；
2. `compute_vorabrechnung()` / `calc_drying_pyrolysis_sources()` 实际只看：
   - `m_solid_zu[VM] + m_solid_in[VM]`
3. 但 `global_nr` 最终状态向量里，上部 cell 仍保留了非零 `state_VM`
4. 因而形成了：
   - **状态里有 VM**
   - **热解源项入口却是 0**
   - 即 `orphan VM`

这意味着当前 `CH4 / tar` 偏差不只是“参数可能不对”，还可能包含：

- **状态变量与热解源项链路脱钩**

### 当前结论后的优先级重排

下一轮若继续按 Hamel 对齐，我认为优先级应明确改成：

1. 先核 `R8_a_R8` 的 Hamel source-of-truth，而不是继续纠缠 `atm/bar`；
2. 先解决 `orphan VM` / `VM state` 与 `pyrolysis source` 脱钩问题；
3. 之后再判断是否需要改 `TAR2` 热解分配或重新设计轴向 VM release 口径。

---

## 2026-04-08 续报：消除 orphan VM，并按 Hamel 口径下调 `R8_a_R8`

这轮做了两类实现改动，并都重新验证：

1. **NR 边界投影**
   - 在 `src/core/reactor.py` 的 `_apply_all_bc_for_nr()` 中加入投影：
     - 若某 cell 从 `zu + in + rez` 接收到的 `VM` 为 0，则强制 `m_solid[:, VM] = 0`
     - 若某 cell 从 `zu + in + rez` 接收到的 `MOISTURE` 为 0，则强制 `m_solid[:, MOISTURE] = 0`
   - 在 `src/solvers/global_nr_solver.py` 中让 `x` 在每次 `global_residual()` 后重新与投影后的 `cells` 对齐，避免 line-search / Jacobian 继续拿陈旧状态向量。

2. **`R8_a_R8` 口径对齐**
   - 将 `src/kinetics/gas_reactions.py` 中 `R8_a_R8` 从 `1.0` 调整为 `0.02`
   - 依据：本地 Hamel 提取记录与 `docs/chat-record.md` 多处都指向煤灰/焦场景 `a_R8 = 0.02`

### orphan VM 已消除

`python3 scripts/audit_vm_transport_consistency_lu.py` 现在显示：

- `active pyrolysis cells = [0]`
- `orphan VM cells = []`

也就是说，这轮之后 shared `global_nr` 下已经不再出现：

- 上部 cell `state_VM > 0`
- 但 `pyro_VM_in = 0`

的脱钩状态。

### NR wall-time / residual

`python3 scripts/audit_global_nr_walltime.py`：

- `dense_fd`
  - `wall_s ≈ 70.73s`
  - `n_iter = 15`
  - `rms_scaled_final ≈ 2.576e-02`
- `block_tridiag_fd`
  - `wall_s ≈ 17.82s`
  - `n_iter = 15`
  - `rms_scaled_final ≈ 2.576e-02`

对比：

- `wall-time speedup ≈ 3.97x`
- `jacobian-build speedup ≈ 3.03x`
- `jacobian cell residual calls: 2700 -> 783`

说明：

- `orphan VM` 投影没有破坏 shared NR；
- `block_tridiag_fd` 仍然是明确有效的加速手段；
- 但 `R8_a_R8=0.02` 之后 inner NR 迭代数从 `11 -> 15`，因此总 wall-time 比上一版更长。

### LU parity：出口组分变化方向

`python3 scripts/audit_solver_parity_lu.py` 当前结果：

- `gauss_seidel`
  - `Texit = 1139.7 K`
  - `Xc = 0.8969`
  - `CO = 0.1759`
  - `CO2 = 0.1176`
  - `H2 = 0.1776`
  - `CH4 = 0.0300`
- `global_nr`
  - `Texit = 1192.0 K`
  - `Xc = 1.0000`
  - `CO = 0.1777`
  - `CO2 = 0.0841`
  - `H2 = 0.1509`
  - `CH4 = 0.0021`
  - `rms_scaled_final = 2.576e-02`

相对于这轮之前：

- `global_nr`
  - `CO` 明显下降（更接近参考）
  - `CO2` 明显上升（更接近参考）
  - `CH4` 有回升，但仍远低于目标
  - `H2` 仍偏高
  - `Texit` 进一步升高到 `1192 K`

同时要注意：

- `R8_a_R8` 是机制级参数，所以 **GS 基线也一起被带动**
- 新 GS 出口现在 `CO/CO2/CH4` 更接近参考，但 `H2` 偏高更明显

### 当前 profile 结论

`python3 scripts/audit_global_nr_profile_lu.py`：

- exit dry gas:
  - `CO = 0.1777` vs `0.1300` → `36.68%`
  - `CO2 = 0.0841` vs `0.1100` → `23.51%`
  - `H2 = 0.1509` vs `0.1200` → `25.76%`
  - `CH4 = 0.0021` vs `0.0280` → `92.36%`
- profile MAE:
  - `CO_wet` 改善到 `0.0312`
  - `CO2_wet` 约 `0.0435`
  - `H2_wet` 约 `0.0287`
  - `CH4_wet` 仍约 `0.0198`
  - `H2O_wet` 仍最差，约 `0.1397`

这说明：

1. `R8_a_R8=0.02` 把 `CO/CO2` 往更合理方向拉回了；
2. 但 `CH4` 仍然几乎起不来，所以热解/甲烷源头问题依旧存在；
3. `H2O deficit` 仍是全局最大的 profile 偏差；
4. 当前 shared NR 仍更像“偏热、偏完全转化”的 branch，只是 syngas composition 已比前一版更接近文献。

### 本轮验证

- `python3 -m py_compile src/core/reactor.py src/solvers/global_nr_solver.py src/kinetics/gas_reactions.py tests/test_global_nr_solver.py`
- `python3 tests/sanity_checks.py` → `11/11 PASS`
- `python3 -m pytest tests/test_global_nr_solver.py tests/test_phase_gate_30.py tests/test_table2_LU.py tests/test_cell_kinetics.py -q`
  - `11 passed, 1 xfailed`

### 下一步建议

这轮之后，优先级可以进一步收敛成：

1. 继续查 `CH4` 源头：
   - `TAR2` 仍未参与热解分配
   - `VM` 仍只在 `cell0` 形成热解源
2. 继续查 `H2O deficit`：
   - 现在 `R8` 已按 Hamel 口径减弱，若 `H2O` 仍显著偏低，则更该回头查 `R2 / pyrolysis / steam feed` 口径
3. 不再回头纠缠 `orphan VM`，因为这条 solver consistency 问题已经被消掉

---

## 2026-04-08 续报：热解双 tar 分配已接通，但 `CH4` 仍在 `cell0` 被 `R6` 当场吃掉

这轮继续沿“先核 source，再看 sink”推进。

### 实现改动

修改了 `src/core/cell_pyrolysis.py`：

1. `allocate_pyrolysis_products_elemental()` 不再只把 tar 分到 `TAR1`
2. 现在按 `calc_tar_surrogate_fractions(fuel_type)` 的双代理 H/C 配比，把 tar 分到：
   - `coal`: `TAR1=C6H6`, `TAR2=C10H8`
   - `biomass`: `TAR1=C10H8`, `TAR2=C16H34`
3. 同时修正了 `carbon_to_char` 计算，改为把 `TAR1 + TAR2` 的碳都计入气相侧

对应测试/门控也一起更新：

- `tests/test_cell_pyrolysis.py`
- `tests/test_thermal_models.py`
- `tests/sanity_checks.py`

验证：

- `python3 -m pytest tests/test_cell_pyrolysis.py tests/test_thermal_models.py -q` → `8 passed`
- `python3 tests/sanity_checks.py` → `11/11 PASS`

### 热解源头：现在确实是双 tar 了

`python3 scripts/audit_pyrolysis_allocator_lu.py` 现在显示：

- `cell0`
  - `CO ≈ 0.9994`
  - `CH4 ≈ 0.5619`
  - `TAR1 ≈ 0.0237`
  - `TAR2 ≈ 0.0609`
  - `H2O ≈ 0.7919`
- totals:
  - `CH4 = 0.5619`
  - `TAR1 = 0.0237`
  - `TAR2 = 0.0609`

也就是说：

1. 之前“热解只生成 `TAR1`”这条实现错误已经修掉；
2. 现在煤工况下 tar 的分配比例与 `coal` 目标 H/C 口径一致；
3. 但 `active_cells` 仍只有 `1`，所以 VM release 的轴向分布问题依旧存在。

### parity / profile：`CO / CO2 / Xc` 继续改善，但 `CH4` 仍没有回来

`python3 scripts/audit_solver_parity_lu.py`：

- `global_nr`
  - `Texit = 1184.5 K`
  - `carbon_conv = 0.9507`
  - `CO = 0.1684`
  - `CO2 = 0.0893`
  - `H2 = 0.1539`
  - `CH4 = 0.0000`
  - `rms_scaled_final = 2.552e-02`

这轮相对上一版：

- `Xc` 从 `1.0000` 回到了 `0.9507`，几乎贴到文献 `95%`
- `CO/CO2` 继续朝参考改善
- `Texit` 也从 `1192.0K` 降到 `1184.5K`
- 但 `CH4` 出口仍几乎为 `0`

`python3 scripts/audit_global_nr_profile_lu.py` 也一致：

- `CO dry` 误差降到 `29.51%`
- `CO2 dry` 误差降到 `18.85%`
- `H2 dry` 仍偏高 `28.24%`
- `CH4 dry` 仍是 `100%` 级偏差
- `H2O wet` 仍是最大的 profile MAE（`~0.1490`）

### `CH4` 的直接 sink 已定位：不是 `R7`，而是 `cell0` 的 `R6`

`python3 scripts/audit_ch4_h2o_paths_lu.py` 现在给出的关键读数：

- `cell0`
  - `CH4_py ≈ +0.562`
  - `CH4_R6 ≈ -1.530`
  - `CH4_R7 ≈ +0.069`
  - `CH4_R11 ≈ 0`
  - `CH4_net ≈ -0.898`

而 `cell1-9`：

- `CH4_py = 0`
- `CH4_R6 ≈ 0`
- `CH4_R7 > 0`（逆向回补少量 `CH4`）
- 但回补量远不足以恢复出口甲烷

这意味着：

1. 当前 shared NR 下，`CH4` 不是在上部被继续重整掉；
2. 相反，它是在 `cell0` 一释放出来，就被 `R6 CH4 oxidation` 强烈吃掉；
3. 而这又和当前结构事实直接相关：
   - `VM` 仍只在 `cell0` 形成热解源
   - `cell0` 同时又是 `O2` 最富、温度最高的区段

### `R7 / R8` 当前语义

`python3 scripts/audit_reaction_progress_lu.py` 显示：

- `R7` 几乎全床层仍是逆向净速率（`drv7 < 0`, `r7 < 0`）
  - 所以它在当前解里更像是**回补 CH4**，不是继续把 CH4 往下吃
- `R8` 下部仍偏逆向，上部趋向正向
  - `H2O deficit` 仍明显存在，但 `CH4≈0` 已不能再怪到 `R8`

### 当前判断

到这一步，`CH4` 的问题主链已经很具体了：

1. 热解源头现在已不再是“单 tar bug”
2. 但 VM/CH4 源项仍集中在 `cell0`
3. `cell0` 的 `R6` 会把这部分 `CH4` 当场烧掉
4. 所以上部虽然 `R7` 有逆向回补，也恢复不了出口甲烷

换句话说：

- **当前 `CH4` 偏低的更根本问题不是 `R7`**
- 而是 **“VM 释放时空位置” 与 “底部高 O2 / 高 T combustion zone” 强重叠**

### 下一步建议

如果继续，我建议优先做：

1. 一个轻量灵敏度探针：
   - 只审 `cell0` 的 `R6 / O2 / VM` overlap
   - 看 `CH4` 对 `R6_scale` 或底部 `O2` 相分配的响应
2. 或者更贴物理地继续查：
   - 为什么当前 shared NR / Vorabrechnung 口径下，VM release 仍几乎只压在 `cell0`
   - 是否该把论文里“固定 VM yields”重新翻译成不同的 reactor-level source 装配，而不是当前这种 bottom-only release

---

## 2026-04-08 续报：`cell0` overlap 灵敏度指向底部 O2 相分配，而不是简单下调 `R6`

这轮先不再直接改主模型，而是做了几组 shared `global_nr` 的无代码灵敏度试验，
专门看：

- `r6_scale`
- `gas_inlet_dense_frac`
- 两者联动

同时新增脚本：

- `scripts/audit_cell0_overlap_sensitivity_lu.py`

用于固化这类开发审计。

### 关键试验结果

#### A. 仅改 `R6`

`dense=0.40` 下：

- `r6_scale=1.00`
  - `Texit ≈ 1191.5 K`
  - `Xc ≈ 0.9222`
  - `CO ≈ 0.1844`
  - `CO2 ≈ 0.0869`
  - `H2 ≈ 0.1948`
  - `CH4 ≈ 0.0000`
- `r6_scale=0.50`
  - `Texit ≈ 1254.5 K`
  - `Xc ≈ 1.0000`
  - `CO ≈ 0.1490`
  - `CO2 ≈ 0.1200`
  - `H2 ≈ 0.0804`
  - `CH4 ≈ 0.0001`
- `r6_scale=0.25`
  - `Texit ≈ 1259.5 K`
  - `Xc ≈ 1.0000`
  - `CO ≈ 0.1497`
  - `CO2 ≈ 0.1197`
  - `H2 ≈ 0.0743`
  - `CH4 ≈ 0.0001`

判断：

1. 下调 `R6` 确实会把 `CO/CO2` 拉向参考；
2. 但几乎**救不回 `CH4`**；
3. 同时会把 `Texit` 明显抬高，并把 `H2` 压得过低；
4. 所以“直接砍 `R6`”更像粗暴调参，不像主控旋钮。

#### B. 仅改底部 O2 相分配

`r6_scale=1.00` 下：

- `dense=0.30`（当前 shared NR 口径）
  - `Texit ≈ 1184.5 K`
  - `Xc ≈ 0.9507`
  - `CO ≈ 0.1684`
  - `CO2 ≈ 0.0893`
  - `H2 ≈ 0.1539`
  - `CH4 ≈ 0.0000`
  - `cell0 wet O2 ≈ 0.0601`
- `dense=0.20`
  - `Texit ≈ 1194.5 K`
  - `Xc ≈ 1.0000`
  - `CO ≈ 0.1214`
  - `CO2 ≈ 0.1218`
  - `H2 ≈ 0.1504`
  - `CH4 ≈ 0.0008`
  - `cell0 wet O2 ≈ 0.0700`
- `dense=0.40`
  - `Texit ≈ 1191.5 K`
  - `Xc ≈ 0.9222`
  - `CO ≈ 0.1844`
  - `CO2 ≈ 0.0869`
  - `H2 ≈ 0.1948`
  - `CH4 ≈ 0.0000`
  - `cell0 wet O2 ≈ 0.0430`

判断：

1. `dense=0.20` 比当前 `dense=0.30` 更能把 `CO/CO2` 同时拉到接近 LU 参考；
2. `CH4` 虽仍很低，但开始出现微弱回升；
3. 这说明当前 shared NR 下，**底部 O2 的相分配** 比“硬砍 `R6`”更像主控旋钮。

#### C. 联动：`r6=0.5` + `dense=0.30`

- `Texit ≈ 1261.1 K`
- `Xc ≈ 1.0000`
- `CO ≈ 0.1150`
- `CO2 ≈ 0.1404`
- `H2 ≈ 0.0673`
- `CH4 ≈ 0.0000`

判断：

- 这组虽然把 `CO` 拉得更低，但整体已经明显过头：
  - `CO2` 过高
  - `H2` 过低
  - `Texit` 过高
- 所以它不适合直接作为 shared NR 新口径。

### 这轮后的总判断

到目前为止，`CH4≈0` 的链路可以进一步收缩成：

1. 热解源头已不再是单 tar bug；
2. 但 VM/CH4 仍几乎只在 `cell0` 释放；
3. `cell0` 同时仍处于高 `O2` / 高 `T` overlap 区；
4. 因而 `R6` 会把 `CH4` 当场吃掉；
5. 而改善这件事最敏感的旋钮，当前看更像：
   - **底部 O2 相分配**
   - 而不是单独修改 `R6 kinetics`

### 下一步建议

若继续推进 global NR 主线，我建议优先顺序改成：

1. 先把 shared `global_nr` 开发口径单独审到 `gas_inlet_dense_frac ≈ 0.20`
   - 仅作为 NR 开发线，不回写 GS 稳定基线
2. 然后再看 `CH4` 是否仍被 `cell0` 的 `R6` 当场清零
3. 若仍然如此，再去做更细的 `cell0` local overlap / limiter 审计

---

## 2026-04-08 续报：shared `global_nr` 已与 GS 基线解耦，NR 开发口径切到 `dense=0.20`

根据上一轮 `cell0 overlap` 灵敏度，当前更像主控旋钮的是底部 `O2` 相分配，
而不是直接修改 `R6 kinetics`。这轮我把它正式落成了 **NR-only 开发口径**，
同时保持 GS 基线不变。

### 实现改动

在 `tests/validation_case_utils.py` 中新增：

- `build_phase1_htw_lu_global_nr_reactor_config()`

口径：

- 基于现有 `build_phase1_htw_lu_reactor_config()`
- 但仅对 shared `global_nr` 开发线覆盖：
  - `gas_inlet_dense_frac = 0.20`

这样现在项目里有两条明确分离的 LU 配置：

1. **GS 基线**
   - 仍使用 `build_phase1_htw_lu_reactor_config()`
   - 保持 `dense=0.30`
2. **shared global NR 开发口径**
   - 使用 `build_phase1_htw_lu_global_nr_reactor_config()`
   - 固定 `dense=0.20`

我也把 NR 相关脚本/测试入口切到了新 builder：

- `tests/test_table2_LU_global_nr.py`
- `scripts/audit_global_nr_profile_lu.py`
- `scripts/audit_global_nr_walltime.py`
- `scripts/audit_reaction_progress_lu.py`
- `scripts/audit_ch4_h2o_paths_lu.py`
- `scripts/audit_pyrolysis_allocator_lu.py`
- `scripts/audit_r8_hamel_alignment_lu.py`
- `scripts/audit_vm_transport_consistency_lu.py`
- `scripts/audit_cell0_overlap_sensitivity_lu.py`

同时修正了 `scripts/audit_solver_parity_lu.py`，现在它明确是：

- **GS baseline vs shared global NR dev config**

不再误把 GS 也切到 `dense=0.20`。

### 测试 / 门控

- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_policy_reports_explicit_init_strategy tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_config_uses_separate_dense_fraction tests/test_table2_LU.py::test_phase1_lu_config_uses_current_tuned_window -q`
  - `3 passed`

### 当前 shared NR（`dense=0.20`）profile 结果

`python3 scripts/audit_global_nr_profile_lu.py`：

- `Texit = 1194.5 K`
- `rms_scaled_final = 2.423e-02`
- exit dry gas:
  - `CO = 0.1214` vs `0.1300` → `6.60%`
  - `CO2 = 0.1218` vs `0.1100` → `10.73%`
  - `H2 = 0.1504` vs `0.1200` → `25.37%`
  - `CH4 = 0.0008` vs `0.0280` → `97.22%`

相比上一版 shared NR（`dense=0.30`）：

- `CO` 明显更接近参考
- `CO2` 也更接近参考
- `H2O wet` profile MAE 略有改善（`~0.1490 -> ~0.1348`）
- `H2` 仍偏高
- `CH4` 仍然接近零
- `Texit` 反而略高

### GS baseline vs NR dev config 对照

`python3 scripts/audit_solver_parity_lu.py`：

- `gauss_seidel`（baseline, `dense=0.30`）
  - `Texit = 1139.7 K`
  - `Xc = 0.8988`
  - `CO = 0.1946`
  - `CO2 = 0.1023`
  - `H2 = 0.1908`
  - `CH4 = 0.0327`
- `global_nr`（shared dev config, `dense=0.20`）
  - `Texit = 1194.5 K`
  - `Xc = 1.0000`
  - `CO = 0.1214`
  - `CO2 = 0.1218`
  - `H2 = 0.1504`
  - `CH4 = 0.0008`
  - `rms_scaled_final = 2.423e-02`

这里的意义是：

1. shared NR 的 `CO/CO2` 已明显优于 GS baseline
2. 但 `CH4` 仍远差于 GS
3. `H2` 仍偏高，`Texit` 仍偏高
4. 所以当前 NR 开发主线已经从“整体都不对”收缩到了：
   - **syngas 里 `CO/CO2` 已基本对上**
   - 剩下主要是 `CH4/H2/Texit`

### 当前判断

这轮之后，我认为可以把 shared NR 主线的优先级进一步收窄成：

1. `dense=0.20` 值得保留为 NR 开发口径
   - 因为它改善的是 `CO/CO2/H2O` 这类更关键的整体 profile
   - 不是单个反应的硬调参
2. `CH4≈0` 依然没有解决
   - 所以主问题继续指向 `cell0` 的 `VM/O2/R6 overlap`
3. 但这一步至少证明：
   - 继续沿 global NR 主线推进是有意义的
   - 旧 GS 判据和旧 GS 最优窗口不应再直接拿来限制 NR

## VM 气源分相灵敏度补充审计（2026-04-08）

这轮继续追 `cell0 VM/O2 overlap`，但先不改主模型，而是用 dev-only monkeypatch 做了一个更窄的科学审计：

- 新增脚本：`python3 scripts/audit_vm_phase_split_lu.py`
- 固定 shared `global_nr` 口径：
  - `nr_init_strategy="vorabrechnung"`
  - `nr_jacobian_strategy="block_tridiag_fd"`
  - NR dev config `gas_inlet_dense_frac=0.20`
- 只改变一件事：
  - 把热解气源 `gas_src_vm` 中的 `0% / 50% / 100%` 分别写入 `bubble`
  - dense 相拿剩余部分

### 先确认的算法/物理口径

本地 Hamel 提取与已有说明并没有给出“热解气必须先进 bubble”这样的证据。相反，现有口径更支持当前默认实现是合理起点：

- 论文把 DAEM / 热解看成颗粒级过程，再把 `CO/H2/CH4/Teer` 作为释放出的挥发分产物汇总。
- `bubble` 在 Hamel 口径下是 **solid-free** 的均相反应区；
- `suspension / emulsion` 才是颗粒、床料和灰催化都在的区域。

因此，“颗粒内热解 → 产物先进入 dense / suspension，再通过 `N_ex` 与 bubble 交换”在物理上是合理推断；这里的 phase-split 审计只是验证这件事在当前 NR 主线下是否值得继续偏离。

### 灵敏度结果

`python3 scripts/audit_vm_phase_split_lu.py` 的关键结果：

- `vm->bubble = 0.00`
  - `Texit = 1194.5 K`
  - `Tpeak = 1359.4 K @ cell1`
  - `CO = 0.1214`
  - `CO2 = 0.1218`
  - `H2 = 0.1504`
  - `CH4 = 0.0008`
  - `rms_scaled_final = 2.423e-02`
  - `n_iter = 14`

- `vm->bubble = 0.50`
  - `Texit = 1369.5 K`
  - `Tpeak = 1526.1 K @ cell1`
  - `CO = 0.0512`
  - `CO2 = 0.1744`
  - `H2 = 0.0441`
  - `CH4 = 0.0000`
  - `rms_scaled_final = 9.934e-03`
  - `n_iter = 51`

- `vm->bubble = 1.00`
  - `Texit = 1379.5 K`
  - `Tpeak = 1538.7 K @ cell1`
  - `CO = 0.0497`
  - `CO2 = 0.1782`
  - `H2 = 0.0427`
  - `CH4 = 0.0000`
  - `rms_scaled_final = 9.964e-03`
  - `n_iter = 53`

### 结论

这条路现在可以视为**不值得继续**，而且原因有两层：

1. **物理结果更差**
   - 把 VM 气源往 `bubble` 挪，不但没有救回 `CH4`
   - 反而把 shared NR 推到更热、更氧化的 branch：
     - `Texit` 直接从 `1194.5 K` 抬到 `~1370-1380 K`
     - `CO` 大幅掉到 `~0.05`
     - `CO2` 被推高到 `~0.175`
     - `H2` 也被压得过低

2. **算法代价显著更高**
   - `n_iter` 从 `14` 暴涨到 `51-53`
   - 对当前 `global_nr` 主线来说，这不符合“尽量减少外层/总迭代次数来加速 solver”的目标

一句话总结：

- **“把 VM 气源分到 bubble 来缓解 `cell0` overlap” 这条想法在当前 Hamel-aligned NR 主线下被证伪了。**
- 当前默认的“热解气先入 dense / suspension”不仅更符合颗粒所在相的物理图景，也明显更有利于收敛与 wall-time。

### 对下一步的影响

因此，下一步不该继续在 `VM phase split` 上烧时间，而应回到更可能影响出口 syngas / profile 且不显著增加迭代成本的路径：

1. 继续审 `cell0` 的 `O2` 相分配 / `N_ex` / `bubble through-flow`
2. 继续审 `CH4` 源头与 `R6` 消耗链，但不要再通过把热解气硬塞进 `bubble` 来做
3. 保持 `vorabrechnung + global_nr + block_tridiag_fd` 为 shared NR 主线，不回退到重型 warmup

## shared global NR 前几格 O2 追踪补充（2026-04-08）

为了避免只凭 `CH4≈0` 做猜测，这轮还把当前 shared NR 口径下 `cell0-2` 的 O2 预算直接拆出来了：

- 新增脚本：`python3 scripts/audit_global_nr_o2_trace_lu.py`
- 复用 `audit_oxygen_reaction_trace.py` 的单 cell O2 分解逻辑
- 但固定到当前 shared `global_nr` 开发口径（`dense=0.20`）

关键读数：

- `cell0`
  - `O2_supply_bubble = 11.502`
  - `O2_supply_dense = 4.220`
  - `bubble_o2_unlimited = 7.918`
  - `dense_o2_unlimited = 189.579`
  - `limit_factor_o2_bubble = 1.000`
  - `limit_factor_o2_dense = 0.022`
  - `O2_shortfall_bubble = 0.000`
  - `O2_shortfall_dense = 185.380`
  - `O2_transfer_potential_bd = 3.526`
  - `vm_oxidizable_o2eq = 1.661`
  - `ch4_net_proxy = -64.731`

- `cell1`
  - `bubble_o2_unlimited = 1.131`
  - `dense_o2_unlimited = 15.260`
  - `limit_factor_o2_bubble = 1.000`
  - `limit_factor_o2_dense = 0.179`
  - `O2_shortfall_dense = 12.533`

- `cell2`
  - `bubble_o2_unlimited = 0.449`
  - `dense_o2_unlimited = 10.861`
  - `limit_factor_o2_bubble = 1.000`
  - `limit_factor_o2_dense = 0.137`
  - `O2_shortfall_dense = 9.368`

### 这组数的含义

现在 shared NR 主线下，前几格的主矛盾已经很清楚：

1. **不是 bubble O2 不够**
   - `cell0-2` 的 `limit_factor_o2_bubble` 都基本是 `1.0`
   - `bubble shortfall` 也几乎为零

2. **而是 dense 相 O2 预算严重短缺**
   - `cell0` 最极端，`dense` 无限制需氧量比预算高出近两个数量级
   - `cell1-2` 仍然明显是 `dense shortfall >> bubble shortfall`

3. **bubble->dense 交换潜力有，但远不足以填平 dense 端缺口**
   - `cell0` 的 `O2_transfer_potential_bd ≈ 3.526`
   - 但 dense 端 shortfall 是 `≈185.380`

4. **因此当前 shared NR 的优先问题不是 bubble kinetics 本身**
   - 而是 `cell0-2` 的 `dense O2 budget / N_ex / phase partition`
   - 以及它们如何把底部 `VM` 与 `O2` 在 dense 端过度重叠

一句话总结：

- 在 shared `global_nr` 主线下，下一步最值得继续查的是 **dense 相 O2 shortfall 的来源**，而不是再去优先怀疑 bubble O2 limiter。

## Hamel-aligned hydrodynamics / freeboard 接入（2026-04-08）

这轮开始把“床层出口 != 反应器出口”落实到代码结构里，并先处理了最接近 Hamel 原文、且能在当前求解器中稳定落地的两块：

1. **hydrodynamics 口径前移**
2. **freeboard-aware reactor exit**

### 1. 床内 hydrodynamics 口径修正

本轮修改了：

- `src/physics/phase_fractions.py`
- `src/core/cell_hydrodynamics.py`
- `src/core/cell.py`
- `src/physics/freeboard.py`

#### 1.1 `eps_b` 不再只用旧的 excess-gas 简化式

当前 `cell_hydrodynamics` 已改成：

- 用 Hamel Chapter 3 的 **visible bubble fraction** 形式估算 `eps_b`
- 用 `u_d = max(u_mf/eps_mf, min(u0, u_b))` 作为 suspension 气体特征速度
- `u_br` 也改为直接吃这个 `u_d`，不再硬绑 `u_mf/eps_mf`

同时补充了：

- `n_rz`
- `u_d`
- `eps_d_voidage`

其中 `eps_d_voidage` 用 Richardson-Zaki 关系：

$$\epsilon_{d,\mathrm{void}} = \epsilon_{mf} \left(\frac{u_d}{u_{mf}}\right)^{1/n_{RZ}}$$

注意：

- 当前 `Cell.eps_d` 仍然是**相体积分率**（用于 `V_d = eps_d * V_cell`）
- 新增的 `eps_d_voidage` 才是更接近 Hamel 文本里的 **emulsion porosity**
- 这样做是为了先和论文口径接轨，同时不把现有两相守恒结构直接打碎

#### 1.2 `beta_A` 改回本地 Hamel 提取口径

`src/physics/freeboard.py` 的 `calc_beta_a(...)` 已从原先的工程近似改成：

$$\beta_A = \frac{1}{3 \cdot d_{b,ws}}$$

这和本地 Hamel 提取记录（`docs/chat-record.md`，Chapter 3 / Eq. 3.75）一致。

### 2. freeboard-aware reactor exit

本轮新增：

- `src/core/freeboard_segment.py`
- `scripts/audit_freeboard_exit_lu.py`
- `tests/validation_case_utils.py::build_phase2_htw_lu_freeboard_reactor_config`

实现方式是：

- **不**把 freeboard 硬塞进当前 bed 内两相 `Cell` 主求解链
- 而是在 bed solve 完成后，追加一个 **gas-only / thermal freeboard axial segment**
- 用 `u_gb` / `beta_A` 给出自由板段停留时间
- 在 freeboard 段中继续计算 gas-only 反应
  - `R5, R6, R7, R8, R10, R11b, R12`
- 并通过 gas enthalpy 守恒更新温度

这样当前结果对象里已经明确区分：

- `bed_T_profile`
- `bed_exit_gas` / `bed_exit_gas_dry`
- `freeboard_T_profile`
- `reactor_exit_T`
- `exit_gas` / `exit_gas_dry`（现在是 reactor exit，不再默认为 bed exit）

### 3. freeboard-aware LU 审计结果

`python3 scripts/audit_freeboard_exit_lu.py`：

- `freeboard_active = True`
- `H_freeboard = 8.50 m`
- `n_freeboard_cells = 8`
- `beta_A = 0.7817`

bed vs reactor exit：

- `bed_exit_T = 1182.0 K`
- `reactor_exit_T = 1181.7 K`
- `ref_exit_T = 1100.0 K`

干基主气相：

- `CO`
  - bed exit: `0.1445`
  - reactor exit: `0.1640`
  - ref: `0.1300`
- `CO2`
  - bed exit: `0.1008`
  - reactor exit: `0.0857`
  - ref: `0.1100`
- `H2`
  - bed exit: `0.1550`
  - reactor exit: `0.1361`
  - ref: `0.1200`
- `CH4`
  - bed exit: `0.0000`
  - reactor exit: `0.0013`
  - ref: `0.0280`

这至少证明了两件事：

1. **当前代码结构里，bed exit 和 reactor exit 已经被真正区分开**
2. **freeboard 段确实会继续改写出口组成**

但也要明确：

- 这还是“最小 freeboard 集成”
- 不是 Hamel 全部连接矩阵 / cyclone / duct / particle-ejection 的完整实现

### 4. 对现有 shared NR 的影响

在 `H_freeboard=0` 的 shared NR bed-only 口径下，新的 hydrodynamics 修正已经改变了 baseline 数字：

`python3 scripts/audit_global_nr_profile_lu.py`：

- `Texit = 1182.0 K`
- `CO = 0.1445`
- `CO2 = 0.1008`
- `H2 = 0.1550`
- `CH4 = 0.0000`
- `n_iter = 12`
- `rms_scaled_final = 4.083e-02`

相比前一版 shared NR：

- `CO/CO2` 仍在可讨论区间，但方向已经变化
- `CH4` 仍然没有被救回
- `H2O` / `O2` profile 仍需继续审

### 5. 验证与回归

本轮已验证：

- `python3 -m py_compile src/physics/phase_fractions.py src/core/cell_hydrodynamics.py src/physics/freeboard.py src/core/freeboard_segment.py src/core/reactor.py tests/validation_case_utils.py scripts/audit_freeboard_exit_lu.py`
- `python3 -m pytest tests/test_cell_balances.py tests/test_global_nr_solver.py -q`
  - `11 passed`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU.py tests/test_table2_LU_global_nr.py tests/test_cell_kinetics.py -q`
  - `10 passed, 1 xfailed`
- `python3 tests/sanity_checks.py`
  - `11/11 PASS`

### 6. 当前最值得继续的方向

这轮之后，下一步最该继续的不是再纠缠“是否要把 freeboard 接入”，因为这件事已经开始做了；而是：

1. 继续核对 **freeboard 反应集合** 是否需要再收窄/重排，以更贴近 Hamel 的 gas-only / thermal 区段
2. 用新的 `bed_exit` vs `reactor_exit` 双口径，重新定义 LU 验证基准
3. 回到 `cell0-2` 的 **dense O2 shortfall / N_ex**，因为 freeboard 接入并没有自动解决底部 `CH4≈0`

## freeboard-aware validation 基准补充（2026-04-08）

这轮进一步把“验证口径也必须跟出口定义一致”补齐了：

- 新增 `scripts/audit_global_nr_profile_lu_freeboard.py`
- 新增 `tests/validation_case_utils.py::build_phase2_htw_lu_freeboard_reactor_config`
- 新增 builder 测试与 freeboard-aware result 测试

### 关键确认

`validation_cases.json` 中 LU case 的轴向 `xi` 明显覆盖整炉：

- `xi = [0.0, ..., 0.45, 0.6, 0.75, 0.8, 0.9, 1.0]`

也就是说：

- `0.45` 附近是 bed top
- `0.6+` 已经在 freeboard
- 以前用 bed-only 模型直接对这组 profile / exit，本来就在结构上不一致

### freeboard-aware profile 审计结果

`python3 scripts/audit_global_nr_profile_lu_freeboard.py`：

- `freeboard_active = True`
- `n_bed = 10`
- `n_freeboard = 8`
- `n_iter = 12`
- `rms_scaled_final = 4.083e-02`
- `nr_total_s = 17.86 s`

温度：

- `bed_exit_T = 1182.0 K`
- `reactor_exit_T = 1181.7 K`
- `ref_exit_T = 1100.0 K`
- 相对偏差 `7.43%`

干基出口：

- `CO = 0.1640` vs `0.1300` → `26.17%`
- `CO2 = 0.0857` vs `0.1100` → `22.14%`
- `H2 = 0.1361` vs `0.1200` → `13.39%`
- `CH4 = 0.0013` vs `0.0280` → `95.33%`

轴向 profile MAE：

- `T_K = 172.2 K`
- `CO_wet = 0.0258`
- `CO2_wet = 0.0201`
- `H2_wet = 0.0197`
- `CH4_wet = 0.0200`
- `O2_wet = 0.0285`
- `H2O_wet = 0.1247`

### 这组结果说明什么

现在可以更明确地把问题拆开：

1. **出口定义已经更一致了**
   - 现在比较的是更接近文献含义的 `reactor exit`
   - 不再把 `bed exit` 直接当最终出口

2. **但当前最小 freeboard 集成还不是最终正确口径**
   - `H2` 相对确实更接近文献
   - 但 `CO/CO2` 被 freeboard 段继续推偏
   - 说明现在的 freeboard 反应集合 / 停留时间 / 热损仍需继续校正

3. **下一步不该回退 freeboard**
   - 而是要继续把 freeboard 本身做对
   - 同时保留 `bed exit` 与 `reactor exit` 双口径，避免再混淆

## freeboard 默认反应集收紧（2026-04-08）

这轮把 freeboard 的默认反应集从：

- `("R5", "R6", "R7", "R8", "R10", "R11", "R12")`

收紧为：

- `("R5", "R6", "R7", "R10", "R11", "R12")`

也就是默认 **不再在 gas-only freeboard 里启用 `R8 WGSR`**。原因很直接：

- 当前 `R8` 实现对应的是煤灰/焦催化 `WGSR`
- 但最小 freeboard 段当前是 gas-only 近似
- 把 `R8` 直接带进去，会在 LU 工况上系统性把 `CO/CO2` 往错误方向推偏

为避免默认值漂移，这条口径已经同时固化到：

- `src/core/reactor.py::ReactorConfig.freeboard_enabled_reactions`
- `src/core/freeboard_segment.py::_freeboard_reaction_step(...)`
- `tests/validation_case_utils.py::build_phase2_htw_lu_freeboard_reactor_config(...)`

### freeboard reaction-set 审计

`python3 scripts/audit_freeboard_reaction_sets_lu.py`：

- `default_no_r8`
  - `Texit = 1198.6 K`
  - `CO = 0.1427`
  - `CO2 = 0.1023`
  - `H2 = 0.1513`
  - `CH4 = 0.0015`
  - `sum|R8| = 0.0000`
- `legacy_with_r8`
  - `Texit = 1181.7 K`
  - `CO = 0.1640`
  - `CO2 = 0.0857`
  - `H2 = 0.1361`
  - `CH4 = 0.0013`
  - `sum|R8| = 1.7160`

对 LU 文献出口（`CO=0.1300, CO2=0.1100, H2=0.1200, CH4=0.0280`）而言：

- 去掉 `R8` 后，`CO` 相对误差从 `26.17%` 降到 `9.79%`
- `CO2` 相对误差从 `22.14%` 降到 `6.96%`
- `CH4` 仍几乎为零，只是从 `0.0013` 到 `0.0015`
- `H2` 仍偏高，而且 `Texit` 被抬到了 `1198.6 K`

### freeboard-aware reactor-exit 审计

`python3 scripts/audit_freeboard_exit_lu.py`：

- `bed_exit_T = 1182.0 K`
- `reactor_exit_T = 1198.6 K`
- `ref_exit_T = 1100.0 K`

干基出口：

- `CO: 0.1445 -> 0.1427`（bed → reactor）
- `CO2: 0.1008 -> 0.1023`
- `H2: 0.1550 -> 0.1513`
- `CH4: 0.0000 -> 0.0015`

说明现在的 freeboard 默认口径已经不再把 `CO/CO2` 往更坏的方向推，但它仍然：

- 把 `Texit` 抬高
- 没有解决 `CH4≈0`
- 对 `H2` 的拉低仍不够

### freeboard-aware profile 审计

`python3 scripts/audit_global_nr_profile_lu_freeboard.py`：

- `converged = False`
- `n_iter = 12`
- `rms_scaled_final = 4.083e-02`
- `nr_total_s = 17.80 s`

出口干基相对误差：

- `CO = 9.79%`
- `CO2 = 6.96%`
- `H2 = 26.06%`
- `CH4 = 94.70%`

轴向 profile MAE：

- `T_K = 177.7 K`
- `CO_wet = 0.0210`
- `CO2_wet = 0.0153`
- `H2_wet = 0.0243`
- `CH4_wet = 0.0200`
- `O2_wet = 0.0285`
- `H2O_wet = 0.1294`

### 当前判断

这轮之后，freeboard 主线已经更清楚了：

1. `R8` 不该作为当前 gas-only freeboard 的默认反应
2. 当前 reactor-exit 的剩余主问题不再是 `CO/CO2`，而是 `Texit / H2 / CH4`
3. 下一步应优先审：
   - `freeboard_heat_loss_frac`
   - `beta_A / tau` 对停留时间的影响
   - `CH4` 在 bed-top 之前就已接近零这一结构性问题

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/validation_case_utils.py tests/test_table2_LU_global_nr.py scripts/audit_freeboard_reaction_sets_lu.py scripts/audit_freeboard_exit_lu.py scripts/audit_global_nr_profile_lu_freeboard.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately -q`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_kinetics.py -q`
- `python3 tests/sanity_checks.py`

## freeboard entrained-solids 依据补充（2026-04-08）

这轮把“论文里的 freeboard 是否真的含固体”进一步压成了可量化证据。

### 原始算法口径

本地 Hamel 提取记录已经足够明确地表明：

- freeboard 不是默认纯气相段
- 床层表面存在由 bubble bursting 触发的颗粒抛射 / 夹带
- 这些颗粒随后进入 freeboard / cyclone / recirculation 拓扑

关键证据：

- `docs/chat-record.md` 问答 48 明确给出 bed→freeboard 的初始固体通量：
  - `F0r`（roof mechanism）
  - `F0w`（wake mechanism）
  - `F0 = F0r + F0w`
- 同一提取中还明确写了：
  - particles are "thrown" upward at the bed/freeboard interface
  - cyclone output is coupled back to the primary bed cell with separation efficiency `~90–95%`

这说明 Hamel 的原始算法口径应理解为：

1. freeboard 有 entrained solids basis
2. 因而不能先验把 freeboard 简化成“只剩 homogeneous gas reactions”
3. 但 freeboard 的 hydrodynamics 也不是床内那套 dense bubbling bed 两相结构，而是独立的 trajectory / ghost-bubble / cyclone 拓扑

### 新增审计脚本

- `scripts/audit_freeboard_entrained_solids_lu.py`

脚本固定用当前 shared `global_nr + freeboard` LU 口径，在 bed-top 处用 Hamel 提取的 `F0r/F0w` 公式估算初始 entrained-solids 量级。

### LU 审计结果

`python3 scripts/audit_freeboard_entrained_solids_lu.py`：

- `eps_b = 0.0100`
- `eps_d_void = 0.9900`
- `u_b = 2.4929 m/s`
- `d_b = 0.4264 m`
- `F0r = 0.0009 kg/m²/s`
- `F0w = 0.0249 kg/m²/s`
- `F0 = 0.0258 kg/m²/s`
- `m_dot_eject0 = 0.0073 kg/s`

对照当前顶格固体流：

- `top_out_char = 0.0288 kg/s`
- `top_out_ash = 0.0061 kg/s`
- `top_out_char_ash = 0.0349 kg/s`
- `m_dot_eject0 / top_out_char_ash = 0.209`

### 这组数字的意义

即使按当前保守 LU 口径，床顶界面的初始抛射固体也约等于顶格 `char+ash` 流的 `~21%`。

这已经足够说明：

- Hamel 口径下 freeboard 不应默认视为纯气相
- 当前 `src/core/freeboard_segment.py` 的 gas-only freeboard 只能视为最小开发近似
- 下一步若继续对齐论文，应优先补的是：
  - entrained `char/ash/fines` 状态
  - freeboard 内对应的 heterogeneous / catalytic reaction enablement
  - cyclone separation + recycle closure

### 当前代码缺口

现有实现中：

- `src/core/freeboard_segment.py` 仍只携带气相状态
- `src/core/reactor.py::_recycled_solid_stream()` / `_propagated_solid_stream()` 仍仅处理 bed 内 `CHAR + ASH`
- `top_solid_inlet_frac` 只是顶格边界 surrogate，不等价于 freeboard entrainment model

所以这轮之后，freeboard 主线的优先级已经更明确：

1. 不再争论“是否该纯气相”
2. 直接转向 entrained-solid-aware freeboard 结构设计

### 验证

- `python3 -m py_compile scripts/audit_freeboard_entrained_solids_lu.py`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py`

## entrained-solid-aware freeboard 骨架（2026-04-08）

这轮没有直接上完整 cyclone/freeboard 联立，而是先把最小可执行的 entrained-solids 骨架接进了当前 freeboard 段。

### 实现

核心改动：

- `src/core/freeboard_segment.py`
- `src/core/reactor.py`
- `tests/test_global_nr_solver.py`

当前最小骨架包含三件事：

1. **freeboard 状态里显式保存 entrained `char/ash` profile**
   - `m_dot_char`
   - `m_dot_ash`
   - `rho_solid`
   - `rho_cat`

2. **bed-top 初始抛射量不再只存在于审计脚本里**
   - `simulate_freeboard(...)` 内部直接用 Hamel 提取的 `F0r/F0w` 计算 `m_dot_eject0`
   - 结果对象新增：
     - `freeboard_entrained_eject_flux_kg_m2_s`
     - `freeboard_entrained_eject_char_ash_kg_s`
     - `freeboard_entrained_char_kg_s`
     - `freeboard_entrained_ash_kg_s`
     - `freeboard_entrained_solid_density_kg_m3`
     - `freeboard_entrained_catalyst_density_kg_m3`

3. **`R11` 的悬浮相催化路径已接入 freeboard**
   - `rate_R11_suspension(...)` 现在可由 freeboard 内的 `rho_cat` surrogate 驱动
   - `freeboard_reaction_diag` 额外拆出了：
     - `R11b`
     - `R11d`

### 当前 LU 结果

这一步的作用目前更像“把论文缺口变成显式状态”，而不是马上重写出口：

- `python3 scripts/audit_freeboard_entrained_solids_lu.py`
  - `m_dot_eject0 = 0.0073 kg/s`
  - `eject0 / top_out_char_ash = 0.209`
  - `rho_cat[max] = 1.03e-03 kg/m3`
  - `rho_cat[avg] = 2.93e-04 kg/m3`
  - `sum|R11d| = 0.0000e+00`

这说明：

1. freeboard 内确实已经有非零 entrained-solids basis
2. 但按当前 LU 口径，这个最小 `rho_cat` surrogate 还不足以真正点亮 `R11d`
3. 所以下一步若要让 freeboard 里的“固体参与”真正改变出口，需要继续补：
   - 更真实的 entrained hold-up / trajectory
   - cyclone separation / recycle closure
   - 可能还包括 freeboard 内 char 异相路径，而不只 tar catalytic path

### 对现有 reactor-exit 的影响

这轮接入后，`freeboard-aware` 的 reactor-exit 数字没有被带偏：

- `python3 scripts/audit_freeboard_exit_lu.py`
  - `bed_exit_T = 1182.0 K`
  - `reactor_exit_T = 1198.6 K`
  - `CO = 0.1427`
  - `CO2 = 0.1023`
  - `H2 = 0.1513`
  - `CH4 = 0.0015`

- `python3 scripts/audit_global_nr_profile_lu_freeboard.py`
  - `n_iter = 12`
  - `rms_scaled_final = 4.083e-02`
  - `nr_total_s = 16.05 s`
  - reactor-exit 相对误差仍为：
    - `CO = 9.79%`
    - `CO2 = 6.96%`
    - `H2 = 26.06%`
    - `CH4 = 94.70%`

也就是说，当前这一步主要是**把 freeboard 里的固体基底显式化**，还不是出口结果的主调参步骤。

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_kinetics.py -q`
- `python3 scripts/audit_freeboard_exit_lu.py`
- `python3 scripts/audit_global_nr_profile_lu_freeboard.py`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py`
- `python3 tests/sanity_checks.py`

## freeboard cyclone closure 指标接入（2026-04-08）

这轮继续把 freeboard 从“有 entrained basis”推进到“有 closure 指标”，但仍然**没有**把它强耦合回主求解器。

### 实现

新增 / 更新：

- `src/core/reactor.py`
- `src/core/freeboard_segment.py`
- `tests/validation_case_utils.py`
- `tests/test_table2_LU_global_nr.py`
- `tests/test_global_nr_solver.py`
- `scripts/audit_freeboard_exit_lu.py`
- `scripts/audit_freeboard_entrained_solids_lu.py`

新增配置：

- `freeboard_cyclone_capture_char_frac = 0.90`
- `freeboard_cyclone_capture_ash_frac = 0.95`

新增结果字段：

- `freeboard_entrained_exit_char_kg_s`
- `freeboard_entrained_exit_ash_kg_s`
- `freeboard_cyclone_capture_char_kg_s`
- `freeboard_cyclone_capture_ash_kg_s`
- `freeboard_cyclone_recycle_candidate_char_ash_kg_s`

注意：这些量当前仍是 **post-solve closure accounting**，还没有反写进 `bottom recycle` 或 `global_nr` 主残差。

### LU 审计结果

`python3 scripts/audit_freeboard_exit_lu.py`：

- `entrained eject = 0.0073 kg/s`
- `freeboard exit char = 0.0000 kg/s`
- `freeboard exit ash = 0.0000 kg/s`
- `cyclone capture char = 0.0000 kg/s`
- `cyclone capture ash = 0.0000 kg/s`
- `recycle candidate = 0.0000 kg/s`

`python3 scripts/audit_freeboard_entrained_solids_lu.py`：

- `m_dot_eject0 = 0.0073 kg/s`
- `eject0 / top_out_char_ash = 0.209`
- `entrained_exit_char = 0.0000 kg/s`
- `entrained_exit_ash = 0.0000 kg/s`
- `cyclone_capture_char = 0.0000 kg/s`
- `cyclone_capture_ash = 0.0000 kg/s`
- `recycle_candidate = 0.0000 kg/s`

### 这组数字说明什么

这轮最有价值的结论不是“cyclone 现在已经有回料”，而是：

1. **初始抛射基底是非零的**
   - bed/freeboard 界面存在 `O(1e-3 ~ 1e-2 kg/s)` 的 entrained basis

2. **但在当前最小 freeboard 轨迹近似里，固体沿程衰减过快**
   - `beta_A` 主导的指数衰减把 entrained solids 在 freeboard 末端压到了接近 0
   - 所以 `cyclone_capture_efficiency` 现在还不是主控旋钮

3. **下一步真正该查的是 trajectory / hold-up decay，而不是 capture efficiency**
   - 当前 closure accounting 已经把问题定位清楚：瓶颈在 `u_gb / beta_A / tau / solid hold-up`

### 当前判断

这轮之后，freeboard 主线的优先级应改写为：

1. 先审 entrained solids 在 freeboard 中的衰减是否过快
2. 再决定是否需要把 cyclone closure 真正并入 `bottom recycle`
3. 在此之前，不值得急着把 `cyclone_capture_efficiency` 继续精调

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/validation_case_utils.py tests/test_global_nr_solver.py tests/test_table2_LU_global_nr.py scripts/audit_freeboard_exit_lu.py scripts/audit_freeboard_entrained_solids_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 scripts/audit_freeboard_exit_lu.py`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py`

## freeboard hydrodynamics 灵敏度（2026-04-08）

这轮按你建议把 freeboard 的优先级彻底验证了一遍：先查 **A. hydrodynamics / solids movement**，再考虑 **B. char reactions**。

### 实现

新增参数化但**不改变默认口径**：

- `freeboard_u_gb_scale`
- `freeboard_beta_a_scale`

位置：

- `src/core/reactor.py`
- `src/core/freeboard_segment.py`
- `tests/validation_case_utils.py`

新增审计脚本：

- `scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`

### 审计结果

`python3 scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`

扫描了：

- `u_gb_scale = 1.00 / 1.53 / 2.00`
- `beta_A_scale = 0.25 / 0.50 / 1.00 / 1.50`

结果最关键的模式是：

1. **`u_gb_scale` 几乎不影响 reactor-exit**
   - 各组下 `Texit` 都保持 `1198.6 K`
   - `CO / CO2 / H2 / CH4` 基本不变
   - 说明当前 freeboard 停留时间变化并不是主控出口偏差的第一来源

2. **`beta_A_scale` 几乎单独控制 entrained-solids survival**
   - `beta_A_scale = 0.25` → `survive ≈ 0.1899`
   - `beta_A_scale = 0.50` → `survive ≈ 0.0361`
   - `beta_A_scale = 1.00` → `survive ≈ 0.0013`
   - `beta_A_scale = 1.50` → `survive ≈ 0.0000`

3. **而 reactor-exit 气相几乎完全不随这组扫参变化**
   - 说明当前 freeboard 里“固体参与”还没真正进入主反应链
   - 现阶段主要是 hydrodynamics 把 solids 太快衰减掉了

### 结论

这轮已经把优先级顺序锁死了：

1. **先做 A：freeboard hydrodynamics / trajectory / hold-up**
2. **再做 B：char reactions in freeboard**

因为在当前实现里：

- `B` 其实还没有真正展开
- freeboard 内没有显式 `R1–R4` char 路径
- 现有与固体相关的只有一个很弱的 `R11d` surrogate
- solids 还没等到进入反应主链，就先被 `beta_A` 控制的衰减压没了

所以若现在直接上 char reactions，得到的结论大概率还是假象：不是反应把 solids 吃掉，而是 hydrodynamics 先把 solids 送没了。

### 下一步建议

下一步应该优先做：

1. 把 freeboard solids 的 trajectory / hold-up 从当前单指数 surrogate 提升到更接近 Hamel 的轨迹口径
2. 等 solids 在 freeboard 末端不再被“数值上清零”后，再接 `R1–R4` 或其子集

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/validation_case_utils.py tests/test_table2_LU_global_nr.py scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`

## freeboard trajectory model 对比（2026-04-08）

这轮把 freeboard solids trajectory 从单一指数 surrogate 扩成了**可切换的双模型**：

- `surrogate_exp`
- `force_balance`

### 实现

改动位置：

- `src/core/freeboard_segment.py`
- `src/core/reactor.py`
- `tests/validation_case_utils.py`
- `tests/test_global_nr_solver.py`
- `tests/test_table2_LU_global_nr.py`

新增内容：

1. `ReactorConfig.freeboard_trajectory_model`
   - 默认仍为 `surrogate_exp`
   - 先不直接切默认口径，避免把新的 trajectory 假设硬写进现有 freeboard baseline

2. `force_balance` 轨迹分支
   - 用 `Haider-Levenspiel C_D`
   - 按 `gravity + buoyancy + drag` 更新颗粒速度
   - 保留 Gaussian 起始速度分布
   - `beta_A` 只继续作用于 `u_gb(z)`，不再直接当作 solids 质量衰减因子

3. 新增审计脚本
   - `scripts/audit_freeboard_trajectory_models_lu.py`

### 审计结果

`python3 scripts/audit_freeboard_trajectory_models_lu.py`

关键读数：

- `surrogate_exp, betaA=1.00` → `survive = 0.0016`
- `surrogate_exp, betaA=0.50` → `survive = 0.0404`
- `surrogate_exp, betaA=0.25` → `survive = 0.1863`

对比：

- `force_balance, betaA=1.00` → `survive = 1.0000`
- `force_balance, betaA=0.50` → `survive = 1.0000`
- `force_balance, betaA=0.25` → `survive = 1.0000`

而且这 6 组 case 的 reactor-exit 几乎完全相同：

- `Texit = 1198.6 K`
- `CO = 0.1427`
- `CO2 = 0.1023`
- `H2 = 0.1513`
- `CH4 = 0.0015`

### 结论

这轮给出了一个很强的结构性信号：

1. 旧 `surrogate_exp` 的确把 `beta_A` 用成了**过强的 solids mass-decay surrogate**
   - 这解释了为什么前面会出现 `survival ≈ 0`

2. 新 `force_balance` 又走到了另一端
   - 在当前最小模型下几乎变成 `all survive`
   - 说明仅有 `gravity/buoyancy/drag` 还不够，需要补更真实的 **deposition / return / fall-back criterion**

3. 所以当前最合理的下一步不是切默认模型，而是继续完善 `force_balance`
   - 让它既摆脱旧的过强指数衰减
   - 又避免把所有 ejecta 都无损送到 cyclone

### 下一步建议

优先继续做：

1. 在 `force_balance` 下补 **stall / settling / return-to-bed** 判据
2. 让 freeboard 末端保留“非零但非全保留”的 solids basis
3. 然后再决定是否把 `R1-R4` 接入 freeboard entrained-solid path

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/validation_case_utils.py tests/test_global_nr_solver.py tests/test_table2_LU_global_nr.py scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py scripts/audit_freeboard_trajectory_models_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_global_nr_solver.py::test_freeboard_force_balance_trajectory_mode_smoke tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_kinetics.py -q`
- `python3 tests/sanity_checks.py`
- `python3 scripts/audit_freeboard_trajectory_models_lu.py`

## freeboard hydrodynamics class-aware bundle（2026-04-09）

这轮继续往 Hamel freeboard hydrodynamics 靠，重点不是再调化学，而是把 freeboard 里的 solids movement 结构补完整。

### 实现

改动位置：

- `src/core/freeboard_segment.py`
- `src/core/reactor.py`
- `tests/test_global_nr_solver.py`
- `scripts/audit_freeboard_exit_lu.py`
- `scripts/audit_freeboard_entrained_solids_lu.py`
- `scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`
- `scripts/audit_freeboard_trajectory_models_lu.py`

本轮具体补了三件事：

1. freeboard trajectory 从单一 `d_p_bar` 升级成 **size-class × velocity-bin bundle**
   - 入口不再只吃 `d_p_bar / m_char_total / m_ash_total`
   - 改为直接吃 bed-top 的：
     - `d_p_classes`
     - `m_char[:,k]`
     - `m_ash[:,k]`
   - 这更接近 Hamel 的颗粒群轨迹口径，而不是单颗粒平均化

2. 显式输出 `return-to-bed`
   - `force_balance` 分支下，若 bundle 在某 freeboard 段内无法穿越该段，则整股 bundle 记为 `return-to-bed`
   - 结果对象新增：
     - `freeboard_entrained_return_char_kg_s`
     - `freeboard_entrained_return_ash_kg_s`
     - 分段 profile 版本

3. 审计脚本同步看三条 hydrodynamics 去向
   - `eject`
   - `exit/to-cyclone`
   - `return-to-bed`

### 审计结果

`python3 scripts/audit_freeboard_entrained_solids_lu.py`

新的关键信号是：

- `n_size_classes = 1`
- `d_p_bar = 5.0000e-04 m`
- `m_dot_eject0 = 0.0073 kg/s`
- `entrained_exit_char = 0.0000 kg/s`
- `return_char = 0.0000 kg/s`

`python3 scripts/audit_freeboard_trajectory_models_lu.py`

结果变成：

- `surrogate_exp`
  - `betaA=1.00` → `survive = 0.0016`, `return = 0.0000`
  - `betaA=0.50` → `survive = 0.0404`, `return = 0.0000`
  - `betaA=0.25` → `survive = 0.1863`, `return = 0.0000`

- `force_balance`
  - `betaA=1.00` → `survive = 1.0000`, `return = 0.0000`
  - `betaA=0.50` → `survive = 1.0000`, `return = 0.0000`
  - `betaA=0.25` → `survive = 1.0000`, `return = 0.0000`

而 reactor-exit 仍几乎完全不变：

- `Texit = 1198.6 K`
- `CO = 0.1427`
- `CO2 = 0.1023`
- `H2 = 0.1513`
- `CH4 = 0.0015`

### 结论

这轮把 freeboard hydrodynamics 的根因进一步压缩清楚了：

1. 当前 `force_balance` 之所以还是 `all survive`
   - **不是** 因为 `return-to-bed` 逻辑没实现
   - 而是因为当前 LU 口径下 **bed-top 只有 1 个 size class**
   - 代表粒径还是单一 `0.5 mm`

2. 因此 freeboard 现在虽然已经支持 class-aware hydrodynamics
   - 但上游 bed solver 还没有给它真正的颗粒粒径分布
   - 所以它还不可能出现 Hamel 风格的“大颗粒回落 / 小颗粒上带”的中间态

3. 这也意味着下一步优先级已经改变：
   - 现在最该补的是 **bed/top solids size discretization**
   - 不是继续在当前单一粒径 basis 上调 freeboard chemistry

### 下一步建议

接下来优先做：

1. 把 bed/top solids 的多粒径离散真正接回主求解
   - 让 `top.solid.d_p_classes` 不再全部相同
   - 让 `top.m_solid[:, CHAR/ASH]` 真正带尺寸分布

2. 然后再复查 freeboard hydrodynamics
   - 看 `return-to-bed` 是否自然出现
   - 看 `to-cyclone` 是否落在合理量级

3. 等 solids hydrodynamics basis 成立后，再评估 freeboard char chemistry 的方向

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py scripts/audit_freeboard_exit_lu.py scripts/audit_freeboard_entrained_solids_lu.py scripts/audit_freeboard_trajectory_models_lu.py scripts/audit_freeboard_hydrodynamics_sensitivity_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_global_nr_solver.py::test_freeboard_force_balance_trajectory_mode_smoke tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_kinetics.py -q`
- `python3 tests/sanity_checks.py`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py`
- `python3 scripts/audit_freeboard_exit_lu.py`

## 2026-04-09: secondary inlet 观测面诊断

这轮没有改 freeboard 主物理，只新增了 `secondary inlet` 的诊断观测面。

### 代码改动

- 在 [src/core/freeboard_segment.py](src/core/freeboard_segment.py) 给 secondary-injection slice 挂出三组状态：
  - `pre_mix`
  - `post_mix_pre_rxn`
  - `post_rxn`
- 在 [src/core/reactor.py](src/core/reactor.py) 把这些观测字段抬到 `result`
- 在 [tests/test_global_nr_solver.py](tests/test_global_nr_solver.py) 新增 `test_freeboard_secondary_observation_exposes_post_mix_pre_rxn_state`
- 新增审计脚本 [scripts/audit_freeboard_secondary_observation_lu.py](scripts/audit_freeboard_secondary_observation_lu.py)

### 关键结果

`python3 scripts/audit_freeboard_secondary_observation_lu.py`

LU、`sec_air_frac=0.05`、`refine=8`、`distributed_uniform` 下：

- `max O2 post-mix = 9.1e-4`
- `max O2 post-rxn = 0.0`

注入窗口内每个 local slice 都呈现同一模式：

- `O2pre = 0`
- `O2mix ≈ 9e-4`
- `O2post = 0`

同时：

- `COmix < COpre`，`COpost < COmix`
- `H2mix < H2pre`，`H2post ≈ H2mix or slightly recovers`
- `Tpost > Tmix`

这说明当前 secondary inlet 段内的**局部快氧化是存在的**，但它只在 `post-mix/pre-rxn` 观测面可见；若只输出 slice endpoint，就会把局部 `O2 spike` 完全洗掉。

### 结论

到这一步可以把问题进一步收缩为：

1. `secondary air` 不是简单“额外加空气”问题，total-air repartition 路径是对的
2. `O2 spike` 至少有一大部分是**采样面 / operator 定义**问题，而不只是网格粗细问题
3. 后续若要和 Hamel 的 profile 更一致，freeboard 注入段应明确区分：
   - 混合后未反应状态
   - 反应后状态

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py scripts/audit_freeboard_secondary_observation_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_secondary_local_refine_expands_profile_near_injection tests/test_global_nr_solver.py::test_freeboard_secondary_distributed_mode_reports_metadata tests/test_global_nr_solver.py::test_freeboard_secondary_observation_exposes_post_mix_pre_rxn_state -q`
- `python3 scripts/audit_freeboard_secondary_observation_lu.py`

## 2026-04-09: HTW Wesseling temperature anchors re-confirmed

用户重新回到 Hamel 原文，补提取了 HTW Wesseling 两个工况的关键温度锚点。  
这轮不改求解器，只修正 `validation` 口径的证据等级。

### 更新

- [data/validation_cases.json](data/validation_cases.json)
  - `CASE_HTW_WESSELING_1.outputs.axial_profiles._observations.temperature_curve_confidence`
  - `CASE_HTW_WESSELING_2.outputs.axial_profiles._observations.temperature_curve_confidence`

### 修正后的口径

此前我把 `CASE_HTW_WESSELING_1` 的 `ξ≈0.6` 温峰降成了较低置信。  
现在应修正为：

- **主温度锚点是可信的**
  - Air/Steam: `ξ=0.00/0.10/0.20/0.45/0.60/0.80/1.00`
    对应约 `750/1020/1100/1060/1240/1130/1100 K`
  - O2/Steam: `ξ=0.00/0.10/0.20/0.45/0.60/0.80/1.00`
    对应约 `800/1150/1200/1150/1380/1300/1220 K`
- **仍较低置信的是锚点之间的细部曲率**
  - 尤其是 bed 内是否存在更细的局部起伏、谷值与次峰
  - 这部分不应由当前稀疏 digitized 点硬编码断言

也就是说，当前正确的验证语义应是：

1. `ξ≈0.6` 的 secondary-air / secondary-O2 thermal peak **存在**
2. 但 bed 内更细的 `T(x)` 起伏形状仍需更高分辨率重读图
3. 后续 local-audit 里，`T@0.6` 可以保留为主锚点之一，但不应用它单独替代整条曲线的形状约束

### 验证

- `python3 - <<'PY' ... json.loads(Path('data/validation_cases.json').read_text()) ... PY`

## 2026-04-09: HTW Wesseling staged-flow semantics corrected

用户补充了 Hamel Table 7.1 / §7.1.2 中关于 Wesseling staged flow 的更精确信息。  
这轮修正的是 **validation data semantics**，不是求解器本身。

### 修正内容

- [data/validation_cases.json](data/validation_cases.json)
  - `CASE_HTW_WESSELING_1.inputs.gasification_agent`
    - `primary_air_Nm3_h` 不再伪装成显式表格值
    - 新增 `total_air_Nm3_h = 5536.9`
    - 新增 `secondary_air_Nm3_h = 30.1`
  - `CASE_HTW_WESSELING_2.inputs.gasification_agent`
    - 新增 `total_O2_Nm3_h = 1921.6`
    - 新增 `secondary_O2_Nm3_h = 37.5`
    - `secondary_injection_agent` 明确为 `O2`
- [tests/validation_case_utils.py](tests/validation_case_utils.py)
  - `extract_case_inlet_topology_hints()` 现在支持：
    - `total_air_Nm3_h`
    - `total_O2_Nm3_h`
    - `secondary_O2_Nm3_h`
    - generic `secondary_agent_*`
  - 若 `primary` 未直给，但 `total - secondary` 可导，则输出 derived primary flow
- [tests/test_table2_LU_global_nr.py](tests/test_table2_LU_global_nr.py)
  - 更新 LU case 期望
  - 新增 `CASE_HTW_WESSELING_2` 的 `secondary O2` 解析测试

### 现在的正确语义

对 Wesseling：

1. `ER` / `lambda_O2` 约束 **总 oxidant**
2. `secondary stream` 有明确数值：
   - Case 1: `secondary air = 30.1 Nm3/h`
   - Case 2: `secondary O2 = 37.5 Nm3/h`
3. `primary distributor flow` **不是表格显式给值**，只能：
   - 保留为 qualitative source note
   - 或在需要时由 `total - secondary` 做 derived estimate

这比旧口径更准确，因为旧数据把 `5536.9 Nm3/h` 误写成了 `primary_air_Nm3_h`，而它实际对应的是 **total air budget**。

### 验证

- `python3 -m py_compile tests/validation_case_utils.py tests/test_table2_LU_global_nr.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json tests/test_table2_LU_global_nr.py::test_extract_case_inlet_topology_hints_converts_numeric_secondary_air tests/test_table2_LU_global_nr.py::test_extract_case_inlet_topology_hints_supports_wesseling_secondary_o2 -q`

## 2026-04-09: LU freeboard rerun with explicit case secondary air = 30.1 Nm3/h

这轮把 `Phase 2` freeboard builder 正式切到 **case-defined secondary stream**：

- [tests/validation_case_utils.py](tests/validation_case_utils.py)
  - 新增 `apply_case_secondary_stream()`
  - `build_phase2_htw_lu_freeboard_reactor_config()` 现在默认按 case source-of-truth 把 `secondary stream` 从 primary 中扣出，再挂到 freeboard `xi≈0.6`
- [scripts/audit_freeboard_secondary_observation_lu.py](scripts/audit_freeboard_secondary_observation_lu.py)
  - 不再使用假设 `sec_air_frac`
  - 直接使用 `CASE_HTW_WESSELING_1` 的 `secondary air = 30.1 Nm3/h`

### 关键结果

`python3 scripts/audit_global_nr_profile_lu_freeboard.py`

- `bed_exit_T = 1182.0 K`
- `reactor_exit_T = 1199.1 K`（ref `1100 K`，`+9.0%`）
- 出口干基：
  - `CO = 0.1423`（ref `0.1300`，`+9.5%`）
  - `CO2 = 0.1029`（ref `0.1100`，`-6.5%`）
  - `H2 = 0.1550`（ref `0.1200`，`+29.2%`）
  - `CH4 = 0.0000`（ref `0.0280`）

`python3 scripts/audit_freeboard_secondary_observation_lu.py`

- `secondary_agent = Air`
- `secondary_Nm3_h = 30.1`
- `max O2 post-mix = 1.0e-4`
- `max O2 post-rxn = 0.0`

注入窗口内只出现很弱的局部效应：

- `T` 大约只从 `1198.3 K` 缓慢抬到 `1205.6 K`
- `CO` 只从 `0.1335` 缓慢降到 `0.1318`
- `H2` 几乎不动

### 解释

这说明一件很重要的事：

1. 当 secondary stream 改成 **Table 7.1 / §7.1.2 给出的显式值 `30.1 Nm3/h`** 后，
2. 当前 1D freeboard 口径下的局部 `secondary-air` 扰动变得非常弱，
3. 已经不再接近之前用“假设分流比例”扫出来的强局部峰值。

换句话说，先前较强的 `xi≈0.6` 局部响应，确实有一部分来自**人为假设了更大的 secondary split**。  
现在回到文献显式 secondary flow 后，模型更像是在说：

- 这股 `30.1 Nm3/h` secondary air 对 1D fully-mixed slice 来说太弱
- 若 Hamel 图上仍有明显 `O2/T` 局部峰，则还需要额外机制解释：
  - 局部 jet / incomplete cross-sectional mixing
  - multi-level staging topology
  - 或 HTW 图中 secondary stream 的有效局部浓度并不等同于简单的 1D 全截面平均

### 验证

- `python3 -m py_compile tests/validation_case_utils.py tests/test_table2_LU_global_nr.py scripts/audit_global_nr_profile_lu_freeboard.py scripts/audit_freeboard_secondary_observation_lu.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json tests/test_table2_LU_global_nr.py::test_apply_secondary_air_repartition_preserves_total_air_budget tests/test_table2_LU_global_nr.py::test_apply_case_secondary_stream_uses_explicit_case_budget tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 scripts/audit_global_nr_profile_lu_freeboard.py`
- `python3 scripts/audit_freeboard_secondary_observation_lu.py`

## 2026-04-09: bed exit already dominates reactor exit under explicit case secondary air

新增桥接型审计脚本 [scripts/audit_bed_exit_to_reactor_exit_budget_lu.py](scripts/audit_bed_exit_to_reactor_exit_budget_lu.py)，
把顶床层 3 个 cell 的主反应与 freeboard 净再分配放到同一张表里看。

### 关键结果

`python3 scripts/audit_bed_exit_to_reactor_exit_budget_lu.py`

在 `secondary air = 30.1 Nm3/h` 的正式 case 输入下：

- `bed_exit_T = 1182.0 K`
- `reactor_exit_T = 1199.1 K`

bed exit -> reactor exit 的湿基净变化非常小：

- `CO: 0.13553 -> 0.13272`（`-0.00281`）
- `CO2: 0.09368 -> 0.09593`（`+0.00226`）
- `H2: 0.14498 -> 0.14459`（`-0.00039`）
- `H2O: 0.06747 -> 0.06735`（`-0.00012`）
- `CH4: 0.00000 -> 0.00000`

freeboard 反应积分也很弱：

- `R5 sum ≈ 0.126`
- `R12 sum ≈ 0.004`
- `R6 sum ≈ 9.2e-5`
- `R7 sum ≈ -9.2e-5`
- `R8/R10/R11 ≈ 0`

反过来看顶床层 3 个 cell（`xi/H = 0.75, 0.85, 0.95`）：

- `R7` 都是显著逆向（`-11.0, -13.8, -17.7`）
- `R8` 都是弱逆向（`-0.043, -0.076, -0.116`）
- `R4_area` 仍为正，但已很小（`0.043 -> 0.021`）
- `CH4` 在进入 freeboard 前已经是 `0`
- `H2O_net` 在顶床层仍是显著负值（约 `-6.8 ~ -7.2`）

### 结论

这轮证据把优先级进一步锁死了：

1. 在 case-defined `30.1 Nm3/h secondary air` 下，freeboard 只是**小修饰项**
2. `reactor exit` 组成主要已经在 **bed top** 被决定
3. 后续主线应从 secondary-air / freeboard tuning 转回 **bed chemistry**
4. 尤其要继续查：
   - 为什么顶床层 `R7` 持续强逆向
   - 为什么 `CH4` 在 bed 内已经接近 `0`
   - 为什么 `H2O` 到 bed top 仍维持显著赤字

### 验证

- `python3 -m py_compile scripts/audit_bed_exit_to_reactor_exit_budget_lu.py`
- `python3 scripts/audit_bed_exit_to_reactor_exit_budget_lu.py`

## 2026-04-09: top-bed chemistry is dominated by R2, not by freeboard continuation

新增专项脚本 [scripts/audit_top_bed_chemistry_lu.py](scripts/audit_top_bed_chemistry_lu.py)，只看 `xi/H >= 0.65` 的顶床层。

### 关键结果

`python3 scripts/audit_top_bed_chemistry_lu.py`

顶床层 `cell 6-9`：

- `R2_area = 199.1, 173.8, 161.1, 148.6`
- `R4_area = 0.058, 0.043, 0.031, 0.021`
- `R7_raw = -8.73, -10.93, -13.65, -17.58`
- `R8_raw = -0.002, -0.040, -0.077, -0.118`
- `CH4_net = +0.350, +0.497, +0.657, +0.928`
- `H2O_net = -7.533, -7.297, -6.981, -6.774`

积分后更直观：

- `top-bed sums:`
  - `R2_area = 682.597`
  - `R4_area = 0.153`
  - `R7_raw = -50.883`
  - `R8_raw = -0.237`
  - `CH4_net = +2.431`
  - `H2O_net = -28.586`

### 结论

这把 bed 主线问题进一步锁死了：

1. 顶床层 `H2O deficit` 的主控项不是 `R8`，而是 **`R2` char-steam gasification**
2. 顶床层 `R7` 持续为负，说明当前组成稳定推动**逆向甲烷化**，不是 steam reforming
3. 虽然 `CH4_net` 在顶床层是正的，但 `yCH4≈0`，说明这点逆向回补量远远不够把 `CH4` 拉回到文献层级
4. `R4` 到顶床层时已经弱到可以忽略，不再是主控旋钮

所以后续优先级应继续集中在：

- bed 内 `VM -> CH4` 源头为什么太弱
- `R2` 为什么把蒸汽持续拉成显著赤字
- 上部床层组成为什么会稳定把 `R7` 推到逆向

### 验证

- `python3 -m py_compile scripts/audit_top_bed_chemistry_lu.py`
- `python3 scripts/audit_top_bed_chemistry_lu.py`

## 2026-04-09: bed hydrodynamics consistency audit points to solved-state `u_d` loop as the first root cause

新增专项脚本 [scripts/audit_hydrodynamics_consistency_lu.py](scripts/audit_hydrodynamics_consistency_lu.py)，
把当前主路径与 Hamel 本地提取逐项并排：

- `u_d`
- `eps_d_raw / eps_d_voidage`
- `alpha_b`
- `d_b`（Mori-Wen vs ODE path）
- `lambda_b`（current helper vs Hamel extract）

### 关键结果

`python3 scripts/audit_hydrodynamics_consistency_lu.py`

1. **初始化 top-bed hydrodynamics 是正常量级**

- top-bed initial snapshot:
  - `u_d ≈ 0.097`
  - `eps_b = 0.01`
  - `eps_d ≈ 0.555`

2. **solved-state hydrodynamics 把整床都推成了极端稀相**

- all bed cells:
  - `u0 ≈ 1.40`
  - `u_mf ≈ 0.043 - 0.046`
  - `u_d/u_mf ≈ 31 - 34`
  - `eps_d_raw ≈ 1.07`
  - `eps_d_voidage = 0.99`（被 clip）
  - `dense solid fraction = 0.01`

也就是说，当前 top-bed “过稀”并不是局部渐变问题，而是 solved-state hydrodynamics **整床都塌成了 99% voidage**。

3. **当前 bed 主路径与 Hamel 仍有多处结构差异**

- `d_b`：主路径仍用 **Mori-Wen**，不是 Hilligardt ODE
- `lambda_b`：
  - current helper: `d_b/(0.5*u_b) * (P/P0)^(-0.2)` → `~0.09 - 0.18 s`
  - Hamel extract: `280*u_mf/g * (P/P0)^(-0.7)` → `~0.13 - 0.14 s`
  - 二者量级接近，但公式并不一致
- `alpha_b`：
  - current specs/techspec 主口径：`u_b/u_mf ≈ 44 - 56`
  - chat-record 另一提取：`(u_b/u_d)*eps_b ≈ 0.014 - 0.018`
  - 两者冲突，需回到 Hamel 原文定 source-of-truth
- `K_bd`：
  - 当前主路径始终用 Eq. 3.50 mixed form
  - fast/slow bubble branch logic 并未真正进入主路径

### 当前优先级判断

这轮审计把 hydrodynamics 修正顺序也锁出来了：

1. **第一优先级：`u_d` / `eps_d` solved-state 闭环**
   - 因为它已经单独足以把整床打成 `eps_d=0.99`
2. **第二优先级：`alpha_b` / fast-slow regime source-of-truth**
   - 当前 specs 与 chat-record 冲突
3. **第三优先级：`d_b(h)` ODE / `lambda_b` / `K_bd` branch**
   - 这些仍重要，但在 `u_d` 修正之前不宜先调 chemistry

### 验证

- `python3 -m py_compile scripts/audit_hydrodynamics_consistency_lu.py`
- `python3 scripts/audit_hydrodynamics_consistency_lu.py`

## 2026-04-09: validation JSON 显式化 `secondary/recycle/topology` 提示

### 背景

当前 freeboard / secondary-air / recycle 相关开发已经明显依赖 `validation_cases.json`
里的 case-specific 拓扑信息，但此前这些信息并没有被统一显式化：

- `secondary air` 在 HTW Wesseling Sim 1/2 中只体现在输出 profile 注释里
- `recirculation via cyclone` 只有一个布尔值和简短 note
- `bed/freeboard diameter` 在 HTW case 中只有单一 `diameter_m`

这会导致后续脚本和 builder 继续手写假设，而不是从数据源读取。

### 本轮改动

#### 1. 在 `data/validation_cases.json` 中补显式提示字段

对 `CASE_HTW_WESSELING_1` 与 `CASE_HTW_WESSELING_2` 增补：

- `inputs.gasification_agent.primary_air_Nm3_h`
- `inputs.gasification_agent.secondary_air_Nm3_h`
- `inputs.gasification_agent.secondary_injection_xi_reactor`
- `inputs.gasification_agent.secondary_injection_agent`
- `inputs.gasification_agent.secondary_injection_note`
- `inputs.reactor.bed_diameter_m`
- `inputs.reactor.freeboard_diameter_m`
- `inputs.operating_conditions.recycle_topology`
- `inputs.operating_conditions.cyclone_present`

注意：

- 对 HTW Sim 1/2，`secondary_air_Nm3_h` 仍是 `__MISSING__`
- 这轮**没有**虚构二次风定量流量
- 只是把文献已经明确给出的 topology / location hint 变成结构化字段

#### 2. 在 `tests/validation_case_utils.py` 正式读取这些字段

`load_case_LU()` 现在会额外暴露：

- `recirculation_note`
- `recycle_topology`
- `cyclone_present`
- `bed_diameter_m`
- `freeboard_diameter_m`
- `primary_air_Nm3_h`
- `secondary_air_Nm3_h`
- `secondary_injection_xi`
- `secondary_injection_agent`
- `secondary_injection_note`

同时：

- `build_phase1_htw_lu_*_reactor_config()` 改为使用 `bed_diameter_m`
- 若 `recirculation=false`，则 `recirculation_frac=0` 且 `recycle_gas=False`
- `build_phase2_htw_lu_freeboard_reactor_config()` 会从 case 中继承
  `freeboard_secondary_injection_xi`

也就是说，freeboard builder 现在至少在**注入位置**上，不再需要脚本里额外硬编码 `0.60` 才能与 case 对齐。

### 当前结论

这轮没有直接改变求解结果，但把 case source-of-truth 明显收紧了：

1. `secondary injection @ xi≈0.6`
2. `cyclone-based recirculation`
3. `bed/freeboard geometry`

现在都已经成为可程序化读取的 validation input hint，而不再只是散落在注释与 profile observations 里的文本信息。

### 新增/更新验证

- `tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json`
- `tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment`

这些测试确认：

- `load_case_LU()` 能正确暴露 `secondary/recycle/topology` 提示
- `Phase 2` freeboard builder 会继承 `secondary_injection_xi = 0.60`

### 验证

- `python3 -m py_compile tests/validation_case_utils.py tests/test_table2_LU_global_nr.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment tests/test_global_nr_solver.py::test_freeboard_secondary_injection_xi_uses_global_reactor_coordinate -q`

## 2026-04-09: 通用 validation case inlet/topology 抽取层

### 目标

把 `validation_cases.json` 里跨 case 复用的这些提示统一抽出来：

- `secondary air`
- `primary air`
- `steam feed`
- `recirculation / cyclone`
- `bed/freeboard geometry`

避免后续 freeboard / staging / recycle 开发继续只围绕 `CASE_HTW_WESSELING_1`
手写特例。

### 实施

在 `tests/validation_case_utils.py` 新增：

- `nm3_h_to_mol_s_stp()`
- `air_nm3_h_to_species_mol_s()`
- `extract_case_inlet_topology_hints()`

这个抽取层现在会统一给出：

- `secondary_air_Nm3_h`
- `secondary_air_O2_mol_s`
- `secondary_air_N2_mol_s`
- `primary_air_Nm3_h`
- `primary_air_O2_mol_s`
- `primary_air_N2_mol_s`
- `steam_feed_kg_h`
- `secondary_injection_xi`
- `recirculation`
- `cyclone_present`
- `recycle_topology`
- `reactor_height_m`
- `bed_height_m`
- `bed_diameter_m`
- `freeboard_diameter_m`

并支持 `reactor._ref = "_WSV400_SHARED"` 这类 shared reactor block 的解析。

### 新增审计

新增：

- `scripts/audit_validation_case_inlet_configs.py`

它会给所有 `CASE_*` 输出一张 readiness 表，按：

- `builder-ready`
- `hint-ready`
- `partial`
- `sparse`

分类当前 case 是否已足够驱动 freeboard inlet / recycle builder。

### 当前审计结论

`python3 scripts/audit_validation_case_inlet_configs.py`

结果显示：

- `CASE_HTW_WESSELING_1/2`：`hint-ready`
  - 已有 `secondary_injection_xi=0.60`
  - 已有 `cyclone recirculation`
  - 缺的是二次风定量流量

- `CASE_VTT_PRESSURIZED_PEAT_13`：`hint-ready`
  - 已有 `secondary_air_Nm3_h=41.6`
  - 已有 `steam_feed_kg_h=7.6`
  - 已有 `recirculation=true`
  - 缺的是显式 `secondary_injection_xi`

- `CASE_VTT_PRESSURIZED_PEAT_6` / `CASE_VTT_PRESSURIZED_SAWDUST_14`：`partial`
  - 有 `secondary_air_Nm3_h`
  - 但无 recirculation / explicit xi

- `CASE_WSV400_A-F`：现在都至少是 `partial`
  - `WSV400_B-F` 已通过 `_WSV400_SHARED` 继承到 `bed/freeboard geometry`
  - 不再被错误判成 `sparse`

### 结论

这轮的意义不在于立刻改 solver，而在于把“哪些 case 已经具备 freeboard/staging/recycle 配置基础”明确量化了。

按当前 readiness，下一步如果要从 LU 外扩，我会优先看：

1. `VTT_PRESSURIZED_PEAT_13`
2. `VTT_PRESSURIZED_PEAT_6`
3. `HTW_WESSELING_2`

因为它们比 WSV 系列更接近当前 freeboard + secondary-air 开发主线。

### 验证

- `python3 -m py_compile tests/validation_case_utils.py tests/test_table2_LU_global_nr.py scripts/audit_validation_case_inlet_configs.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json tests/test_table2_LU_global_nr.py::test_extract_case_inlet_topology_hints_converts_numeric_secondary_air tests/test_table2_LU_global_nr.py::test_extract_case_inlet_topology_hints_resolves_wsv_shared_reactor_block tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 scripts/audit_validation_case_inlet_configs.py`

## 2026-04-09: LU case 改成 total-air repartition 审计

### 背景

在重新核对本地文档后，当前更合理的 HTW Sim 1 口径是：

- `secondary` 应理解为 `secondary air`
- `ER / λ_O2` 给的是**总空气预算**
- 本地 source-of-truth 没有给出 `primary/secondary air split`

这意味着先前“把 freeboard secondary 当成额外加空气”的试法并不理想，
因为它会在固定 `ER` 之外再引入额外氧化剂。

### 实施

在 `tests/validation_case_utils.py` 新增：

- `apply_secondary_air_repartition()`

它会：

- 保持总 `O2/N2` 不变
- 只把总空气在 `primary` 和 `secondary` 之间重分配
- 不改变总蒸汽 `H2O_feed`

并新增专项脚本：

- `scripts/audit_freeboard_air_repartition_lu.py`

扫描：

- `secondary air fraction = 0.00, 0.05, 0.10, 0.15, 0.20, 0.25`

同时保持：

- `O2_total = O2_primary + O2_secondary = const`
- `N2_total = N2_primary + N2_secondary = const`

### 审计结果

`python3 scripts/audit_freeboard_air_repartition_lu.py`

得到：

- `sec_air_frac=0.00`
  - `T@0.60 = 1198.6 K`
  - `CO@0.60 = 0.133`
  - `CO2@0.60 = 0.095`
  - `H2@0.60 = 0.141`
  - `Texit = 1198.6 K`

- `sec_air_frac=0.05`
  - `T@0.60 = 1205.0 K`
  - `CO@0.60 = 0.145`
  - `H2@0.60 = 0.145`
  - `Texit = 1276.2 K`

- `sec_air_frac=0.10`
  - `T@0.60 = 1209.9 K`
  - `CO@0.60 = 0.150`
  - `H2@0.60 = 0.150`
  - `Texit = 1360.3 K`

- `sec_air_frac=0.20`
  - `T@0.60 = 1219.8 K`
  - `CO@0.60 = 0.166`
  - `H2@0.60 = 0.160`
  - `Texit = 1527.8 K`

核心现象：

1. 即使**总空气守恒**，当前 coarse freeboard (`n_freeboard_cells=8`) 下仍然**留不住**文献的 `O2 spike`
2. 更糟的是，随着更多空气改到 freeboard，`CO/H2` 在 `xi≈0.6` 不降反升
3. `Texit` 仍然被显著抬高

### 结论

这说明问题已经进一步收缩：

- 主矛盾不再是“是不是额外加了总空气”
- 而是当前 freeboard 中：
  - `secondary air` 的同段瞬时完全混合
  - 同段内 lumped reaction consumption
  - 以及过粗的轴向离散

共同把文献里的局部 `O2 spike + CO/H2 dip` 洗平了

换句话说：

**把总空气改成 repartition 是必要修正，但它本身还不足以复现 Fig. 7.4 的局部 secondary-air 指纹。**

### 直接下一步

下一步最值得做的是：

1. 不再继续大范围扫 `secondary_air_frac`
2. 直接把 freeboard 注入段拆成更细的 sub-step / finer axial discretization
3. 重点检查：
   - `secondary air` 注入后是否被同段内瞬时完全吃掉
   - `O2` 是否因为反应与混合都 lump 在同一段而看不到局部峰

### 验证

- `python3 -m py_compile tests/validation_case_utils.py tests/test_table2_LU_global_nr.py scripts/audit_freeboard_air_repartition_lu.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_apply_secondary_air_repartition_preserves_total_air_budget tests/test_table2_LU_global_nr.py::test_load_case_lu_exposes_secondary_and_recycle_hints_from_json -q`
- `python3 scripts/audit_freeboard_air_repartition_lu.py`

## 2026-04-09: primary / secondary inlet 局部热点单元内诊断

### 目标

用户指出一个关键方向：如果要捕捉快速氧化反应，就不能只盯整体出口，
而必须直接盯：

1. `cell0` 的 `primary inlet`
2. freeboard 的 `secondary inlet`

所以这轮没有去改 solver，而是先补了一个局部热点诊断脚本，看当前实现是否在
单元内就把 `O2` 峰洗掉了。

### 新增脚本

- `scripts/audit_local_oxidation_hotspots_lu.py`

功能：

- 对 `cell0` 复用现有 `_cell_o2_trace()`，输出 `O2 supply / demand / limiter / shortfall`
- 对 freeboard `secondary inlet` 段，重放当前实现里的 `12` 个 reaction substeps
- 逐步打印：
  - `T`
  - `O2`
  - `CO`
  - `CO2`
  - `H2`
  - `CH4`
  - `R5 / R6 / R12`

审计场景：

- `baseline`
- `sec_air_05`（总空气重分配 5% 到 secondary）
- `sec_air_10`（总空气重分配 10% 到 secondary）

### 关键结果

#### 1. `cell0` 仍然是极强的 dense-side fast oxidation hotspot

`baseline` 下：

- `O2sup_b = 11.502`
- `O2sup_d = 3.423`
- `O2dem_b = 2.083`
- `O2dem_d = 10144.031`
- `lim_b = 1.000`
- `lim_d = 0.000`

主要快氧化项：

- `R6d = 9359.727`
- `R10 = 420.741`
- `R1 = 312.091`
- `R5d = 51.686`

这再次说明：

- `cell0` 的主矛盾仍是 **dense-side O2 shortfall + fast oxidation overlap**
- 不是 bubble O2 不够

#### 2. freeboard `secondary inlet` 段的 `O2` 在 `step1-step2` 内就被吃光

以 `sec_air_05` 为例：

- `step1`: `T=1290.6K`, `O2=0.0000`, `R5=0.7016`, `R6=0.0025`, `R12=0.0199`
- `step2`: `T=1276.2K`, `O2=0.0000`

以 `sec_air_10` 为例：

- `step1`: `T=1379.7K`, `O2=0.0001`, `R5=1.4032`, `R6=0.0066`, `R12=0.0346`
- `step2`: `T=1360.3K`, `O2=0.0000`

这条证据非常关键：

**当前 coarse freeboard 口径下，secondary inlet 段并不是“没有局部快氧化”，而是“局部快氧化太快，以至于在同一 cell 的前 1–2 个 substeps 内就把 `O2` 清空了”。**

也就是说，文献里的 `O2 spike` 很可能不是“方向错”，而是被当前数值处理：

- 同段瞬时完全混合
- 同段 lumped reaction
- 过粗轴向离散

共同洗平了。

### 结论

这轮基本把下一步锁死了：

1. 不该再继续扫出口 KPI 或大范围扫 `secondary_air_frac`
2. 应优先改 freeboard `secondary inlet` 段的数值处理
3. 候选方向是：
   - 注入段拆成更细的 sub-step
   - 注入段局部网格加密
   - 或把“混合”和“反应”拆成顺序子过程

### 验证

- `python3 -m py_compile scripts/audit_local_oxidation_hotspots_lu.py`
- `python3 scripts/audit_local_oxidation_hotspots_lu.py`

## 2026-04-09: freeboard secondary inlet 局部加密与注入模式审计

### 背景

上一轮已经确认：

- freeboard `secondary inlet` 段里，`O2` 在 `step1-step2` 内就被吃到接近 0

所以这轮我先做两件事，但都**不改变默认主线**：

1. 给 `secondary inlet` 段加一个可选的 `local refine`
2. 给 `secondary air` 加一个可选的 `injection mode`

目标不是立刻改 baseline，而是看：

- 仅靠局部加密能不能把文献里的 `O2 spike` 放出来
- 若不能，问题是不是在当前 `lumped injection operator`

### 新增能力

在 `src/core/reactor.py` / `src/core/freeboard_segment.py` 新增：

- `freeboard_secondary_local_refine`
- `freeboard_secondary_injection_mode`

其中：

- 默认 `freeboard_secondary_local_refine = 1`
- 默认 `freeboard_secondary_injection_mode = "lumped"`

可选模式：

- `"lumped"`
- `"distributed_uniform"`

含义：

- `lumped`：secondary air 在 injection 段一次性并入
- `distributed_uniform`：secondary air 均匀分摊到该段的局部细分 slice

### 回归

新增测试：

- `tests/test_global_nr_solver.py::test_freeboard_secondary_local_refine_expands_profile_near_injection`
- `tests/test_global_nr_solver.py::test_freeboard_secondary_distributed_mode_reports_metadata`

两条都通过，说明：

- 默认口径未漂移
- 打开局部加密后，freeboard profile 长度会按预期扩展
- 新的 `distributed_uniform` 元数据已正确挂出

### 审计 1：局部加密是否足够

`python3 scripts/audit_freeboard_secondary_local_refine_lu.py`

条件：

- `secondary air fraction = 0.05`
- 比较 `refine = 1 / 4 / 8`

结果：

- `refine=1`：`max O2 near xi≈0.6 = 0.0000`
- `refine=4`：仍为 `0.0000`
- `refine=8`：仍为 `0.0000`

虽然 `refine=8` 后 `xi≈0.64-0.70` 出现了更细的温度与组分渐变，但 `O2` 仍在每个局部 slice 末端被吃到 0。

结论：

**单靠局部加密，不足以恢复文献里的 `O2 spike`。**

### 审计 2：注入模式是否足够

`python3 scripts/audit_freeboard_secondary_modes_lu.py`

条件：

- `refine = 8`
- `secondary air fraction = 0.05`
- 比较 `lumped` vs `distributed_uniform`

结果：

- `lumped`: `maxO2_window = 0.00000`
- `distributed_uniform`: `maxO2_window = 0.00000`

但 `distributed_uniform` 确实带来了更平滑、逐步升温的局部路径：

- `xi≈0.638`：`T≈1212 K`
- `xi≈0.666`：`T≈1241 K`
- `xi≈0.702`：`T≈1279 K`

也就是说：

- `distributed_uniform` 改善了局部 `T/CO/H2` 的演化形态
- 但**仍不足以保留非零 `O2` 到 slice 末端**

### 当前最强结论

到这一步，问题已经收缩成：

1. 不是“网格太粗” alone
2. 也不只是“注入一次性 lumped” alone
3. 更像是当前数值定义里：
   - `post-mix -> immediate reaction within same slice`
   - 且观测点落在 `reaction-updated slice state`

导致任何短寿命 `O2` 峰都在输出之前被吃掉

### 直接下一步

下一步最值得做的不是继续调 `refine` 或 `mode`，而是查：

1. `slice` 内的观测点定义
2. 是否需要显式暴露：
   - `post-mix / pre-reaction` 状态
   - 或 `mixing` 与 `reaction` 的 operator split

如果文献中的 `O2` 采样更接近局部混合后而非完全反应后的截面平均，那么当前输出口径天然会低估 `O2 spike`。

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py scripts/audit_freeboard_secondary_local_refine_lu.py scripts/audit_freeboard_secondary_modes_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_secondary_local_refine_expands_profile_near_injection tests/test_global_nr_solver.py::test_freeboard_secondary_distributed_mode_reports_metadata -q`
- `python3 scripts/audit_freeboard_secondary_local_refine_lu.py`
- `python3 scripts/audit_freeboard_secondary_modes_lu.py`

## 2026-04-09: freeboard `secondary air @ xi≈0.6` 审计与坐标修正

### 背景

本地文档和 `validation_cases.json` 对 HTW Wesseling Sim 1 的指纹是明确的：

- `xi≈0.6` 有 `secondary air`
- 应出现局部 `O2` 回升
- `T` 在该点附近抬到约 `1240 K`
- `CO/H2` 在该点下探，随后 toward exit 略回升

但当前实现里，`freeboard_secondary_injection_xi` 最初是按 **freeboard 局部坐标** 解释的，和文献/验证口径里的**整炉归一化坐标**不一致。

### 实施

我先修了坐标语义：

- `src/core/freeboard_segment.py`
- `src/core/reactor.py`
- `tests/test_global_nr_solver.py`

现在 `secondary_injection_xi` 优先按 **global reactor xi** 解释；`Phase 2 LU` 口径下，`xi=0.60` 会落到正确的 freeboard 段，而不再是 freeboard 内部的 `0.60`。

同时新增审计脚本：

- `scripts/audit_freeboard_secondary_injection_lu.py`

该脚本固定 shared freeboard dev config，扫描：

- `secondary_O2 = 0.00, 0.05, 0.10, 0.15, 0.20 × primary O2`
- 注入位置 `xi = 0.60`
- 注入组分按 `air = O2 + 3.76 N2`

并直接对照文献在 `xi=0.45 / 0.60 / 1.00` 的局部指纹。

### 8 段 freeboard 扫描结论

`python3 scripts/audit_freeboard_secondary_injection_lu.py`

在当前 shared `n_freeboard_cells=8` 口径下：

- `secondary_O2 = 0` 时：
  - `O2@0.60 = 0.000`
  - `T@0.60 = 1198.6 K`
  - `dCO(0.60-0.45) ≈ 0`
  - `dH2(0.60-0.45) ≈ 0`

- 随着二次风从 `0.05 -> 0.20 × primary O2` 增加：
  - `T@0.60` 只小幅抬到 `1210.5 K`
  - `CO@0.60` 只从 `0.133 -> 0.130`
  - `H2@0.60` 只从 `0.141 -> 0.140`
  - `O2@0.60` 始终是 `0.000`
  - 但 `Texit` 会被一路抬高到 `1488.3 K`

这说明在 **8 段 coarse freeboard** 下，二次风量增加主要表现为“把后段整体推热”，而不是复现文献那种尖锐的 `xi≈0.6` 局部 spike/dip。

### 加密 freeboard 的只读补充审计

我又做了一轮只读加密检查，不改主实现，只把 freeboard 提高到 `n_freeboard_cells=16`，并测试：

- `secondary_O2 = 0.10 × primary O2`
- `secondary_injection_xi = 0.60`

结果出现了更接近文献方向的局部响应：

- injection segment 落在 `xi ≈ 0.615`
- `T` 在该点从约 `1198.6 K` 直接跳到 `1348.9 K`
- `CO` 从约 `0.1329` 跳降到 `0.1007`
- `H2` 从约 `0.1408` 跳降到 `0.1367`

但 `O2` 仍然在该段内被瞬时吃光，profile 上保持 `0.0`。

### 当前判断

这轮证据把问题收得更清楚了：

1. `secondary air` 这件事本身是对的，且位置现在已按正确 `global xi` 对齐
2. 当前模型已经能在更细 freeboard 离散下复现 **局部 partial-combustion heat release**，即 `T` 跳升和 `CO/H2` 下探
3. 仍然缺的是文献里的 **残余 O2 spike**
4. 所以主问题更像是：
   - 二次风在单段内混合过快
   - 同段内反应过于 lumped
   - `n_freeboard_cells=8` 对这种局部 injection 过粗

而不是“要不要 secondary air”这个方向本身有误。

### 下一步

下一步最值得做的不是继续盲扫二次风量，而是：

1. 提高 freeboard 局部离散，至少用于 `secondary-air` 审计口径
2. 审 `secondary air` 注入后的同段瞬时完全混合假设
3. 视需要把注入段拆成：
   - injection / partial-mixing sub-step
   - downstream reaction sub-step

这样才有机会把 `O2 spike` 从“被当段内瞬时吃掉”改成“在 profile 上可见”。

### 验证

- `python3 -m py_compile src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py scripts/audit_freeboard_secondary_injection_lu.py`
- `python3 -m pytest tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_global_nr_solver.py::test_freeboard_force_balance_trajectory_mode_smoke tests/test_global_nr_solver.py::test_freeboard_secondary_injection_xi_uses_global_reactor_coordinate tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment -q`
- `python3 scripts/audit_freeboard_secondary_injection_lu.py`

## freeboard chemistry direction audit（2026-04-09）

在 shared freeboard dev 口径切到 `force_balance` 之后，这轮继续查的不是参数，而是：

- freeboard 里到底哪些反应真的在发生
- 当前出口偏差是否还值得继续追 `R10/R11`
- `R8` 是否该回到默认集合

### 审计

新增脚本：

- `scripts/audit_freeboard_direction_lu.py`

同时复跑：

- `scripts/audit_freeboard_reaction_sets_lu.py`

### reaction-set 结论

`python3 scripts/audit_freeboard_reaction_sets_lu.py`

在新的 `force_balance` 基线上：

- `default_no_r8`
  - `Texit = 1198.6 K`
  - `CO = 0.1427`
  - `CO2 = 0.1023`
  - `H2 = 0.1513`
  - `CH4 = 0.0015`
  - `sum|R7| = 0.1399`
  - `sum|R8| = 0.0000`
  - `sum|R11| = 0.0000`

- `legacy_with_r8`
  - `Texit = 1181.7 K`
  - `CO = 0.1640`
  - `CO2 = 0.0857`
  - `H2 = 0.1361`
  - `CH4 = 0.0013`
  - `sum|R8| = 1.7160`

这说明：

1. 一旦在 freeboard 中打开 `R8`，它会显著重分配 `CO/CO2/H2`
2. 但当前 source-of-truth 仍指向 `R8` 是 ash/char-catalyzed WGSR
3. 因此默认 freeboard 反应集**不应**把 `R8` 放回去，它应继续只作为显式审计项

### direction audit 结论

`python3 scripts/audit_freeboard_direction_lu.py`

当前 shared freeboard dev 口径（`force_balance + default_no_r8`）下：

- `CO`: `0.134809 -> 0.132892`（微降）
- `CO2`: `0.094009 -> 0.095290`（微升）
- `H2`: `0.144535 -> 0.140844`（微降）
- `H2O`: `0.067317 -> 0.068933`（微升）
- `CH4`: `0.000000 -> 0.001382`（微升）
- `O2`: `0.000493 -> 0.000000`（被吃尽）
- `TAR1/TAR2`: 全程 `≈ 0`

对应的 freeboard 反应积分：

- `R5 = +0.049315`
- `R7 = -0.139937`（逆向）
- `R12 = +0.001556`
- `R8 = 0`
- `R10 = 0`
- `R11 = 0`

### 物理解释

这轮已经把 freeboard chemistry 方向压得很清楚：

1. **当前 freeboard 不是 tar-limited chemistry**
   - `TAR1/TAR2 ≈ 0`
   - `R10/R11 ≈ 0`
   - 所以继续优先追 `R11` 没有价值

2. **当前 freeboard 主导的是轻微后氧化 + 逆向 R7 修复**
   - 少量残余 `O2` 被 `R5/R12` 吃掉
   - `R7` 以逆向形式微量回补 `CH4`
   - 这也是为什么 `CH4` 从 bed top 的 `0` 抬到 reactor exit 的 `0.0014`

3. **CH4 主问题仍不在 freeboard**
   - 因为 freeboard 只是在做小幅回补
   - 当前 `CH4` 低的主因仍在 bed 末端之前

4. **R8 继续保持关闭是合理的**
   - 它一旦打开，确实强烈改写出口
   - 但从论文/提取口径看，它仍不是 gas-only freeboard 的自然默认项

### 当前结论

到这一步，freeboard chemistry 的优先级已经可以重排：

1. `R10/R11` 不是当前主矛盾
2. `R8` 继续保持审计项，不回默认
3. 后续若继续做 freeboard，本体上更值得查的是：
   - 是否要引入 **secondary oxygen / staging** 口径
   - 以及在工业 HTW case 下，freeboard 是否需要更接近论文的分段进气/旋风/连接管拓扑

### 验证

- `python3 -m py_compile scripts/audit_freeboard_direction_lu.py`
- `python3 scripts/audit_freeboard_reaction_sets_lu.py`
- `python3 scripts/audit_freeboard_direction_lu.py`
- `python3 scripts/audit_freeboard_trajectory_models_lu.py`

## freeboard dev 口径切到 force_balance（2026-04-09）

这轮继续沿“先把 freeboard hydrodynamics 建对，再看化学方向”推进，结论已经足够强，可以把 shared freeboard 开发口径从 `surrogate_exp` 切到 `force_balance`。

### 新证据

我给 freeboard 结果对象补了终端沉降/携带诊断：

- `u0`
- `u_gb`
- `u_p_mean`
- `u_t_mean`
- `carry_ratio = u_gb / u_t`

相关实现：

- `src/physics/freeboard.py`
- `src/core/freeboard_segment.py`
- `src/core/reactor.py`
- `scripts/audit_freeboard_entrained_solids_lu.py`
- `scripts/audit_freeboard_exit_lu.py`

`python3 scripts/audit_freeboard_entrained_solids_lu.py` 给出的关键读数是：

- `u0[avg] = 1.4255 m/s`
- `u_gb[avg] = 1.6557 m/s`
- `u_p[avg] = 1.7337 m/s`
- `u_t[avg] = 0.8019 m/s`
- `carry_ratio[min/max] = 1.784 / 3.066`

也就是说，在当前 LU 的**单粒径** freeboard 口径下，`u_gb > u_t` 是持续成立的。  
这意味着：

1. `force_balance` 下几乎全 carry-over 是**物理上一致**的结果
2. 旧 `surrogate_exp` 仍把 solids 压成接近 0，则已经不只是“简化”，而是与当前 carry criterion **冲突**

### 口径切换

因此我把 shared Phase 2 freeboard 开发 builder 改成：

- `freeboard_trajectory_model = "force_balance"`

位置：

- `tests/validation_case_utils.py`
- `tests/test_table2_LU_global_nr.py`

注意：

- `ReactorConfig` 的库级默认值仍保留 `surrogate_exp`
- 但**shared freeboard dev config** 已切到 `force_balance`
- 后续 freeboard chemistry 审计应以这条口径为准

### 切换后的 LU 结果

`python3 scripts/audit_freeboard_entrained_solids_lu.py`

- `entrained_exit_char = 0.0060 kg/s`
- `entrained_exit_ash = 0.0013 kg/s`
- `cyclone_capture_char = 0.0054 kg/s`
- `cyclone_capture_ash = 0.0012 kg/s`
- `recycle_candidate = 0.0066 kg/s`

`python3 scripts/audit_freeboard_exit_lu.py`

- `u0_end = 1.4274 m/s`
- `u_gb_end = 1.4305 m/s`
- `u_p_end = 0.6287 m/s`
- `u_t_end = 0.8019 m/s`
- `carry_ratio_end = 1.784`

出口气仍基本没变：

- `Texit = 1198.6 K`
- `CO = 0.1427`
- `CO2 = 0.1023`
- `H2 = 0.1513`
- `CH4 = 0.0015`

说明这次切换主要是在**把 hydrodynamics basis 校正到更一致**，而不是立刻重写 reactor-exit chemistry。

### 结论

到这一步可以把 freeboard hydrodynamics 的当前结论锁定为：

1. 在 **当前单粒径 LU 口径** 下，非零 carry-over / cyclone candidate 是合理的
2. `surrogate_exp` 已不再适合作为 freeboard 开发主线
3. 后续 freeboard chemistry direction audit 应基于 `force_balance`

### 下一步建议

下一步不该再围绕 `surrogate_exp` 做比较，而应直接进入：

1. freeboard reaction activation / deactivation 审计
2. 判断在“有 entrained solids carry-over”口径下，哪些反应该继续开启
3. 重点看 tar / CO / H2 / CH4 在 freeboard 段内的净方向

### 验证

- `python3 -m py_compile src/physics/freeboard.py src/physics/__init__.py src/core/freeboard_segment.py src/core/reactor.py tests/test_global_nr_solver.py scripts/audit_freeboard_entrained_solids_lu.py scripts/audit_freeboard_exit_lu.py`
- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_phase2_lu_freeboard_config_restores_reactor_height_segment tests/test_global_nr_solver.py::test_freeboard_aware_result_exposes_bed_and_reactor_exit_separately tests/test_global_nr_solver.py::test_freeboard_force_balance_trajectory_mode_smoke -q`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_kinetics.py -q`
- `python3 tests/sanity_checks.py`
- `python3 scripts/audit_freeboard_entrained_solids_lu.py`
- `python3 scripts/audit_freeboard_exit_lu.py`

---

## 2026-04-09: `u_d` closure sensitivity confirms bed hydrodynamics is the first-order root cause

### 背景

在 LU shared `global_nr` 口径下，之前已经确认：

- solved-state 全床几乎都落到
  - `u_d ~= u0`
  - `eps_b = 0.01` floor
  - `eps_d_voidage = 0.99` clip
- 这会把 top-bed 直接推成 `~1%` solid fraction

本轮目标不是立刻改主路径，而是把 `u_d` 从硬编码实现改成**显式可审计 closure**，量化它对 hydrodynamics 和 top-bed chemistry 的影响。

### 实现

新增最小开发入口：

- `ReactorConfig.hydrodynamics_u_d_closure = "current"`（默认不变）
- `Cell.u_d_closure`
- `calc_cell_hydrodynamics(..., u_d_closure=...)`

当前支持三条口径：

1. `current`
   - `u_d = max(u_mf/eps_mf, min(u0, u_b))`
2. `umf_over_epsmf`
   - `u_d = u_mf / eps_mf`
3. `backsolve_visible_epsb`
   - 先用当前 `eps_b = (u0-u_mf)/u_b`
   - 再由 visible-bubble 公式
     - `eps_b = (u0-u_d)/(u_b + (n_b-1)u_d)`
   - 反解 `u_d`

注意：

- 第 3 条只是把当前仓库里**已经存在的两条相分率关系**拼成一个显式 closure
- 它不是已经确认的 Hamel 原文 source-of-truth

### 新审计

新增：

- `scripts/audit_hydrodynamics_ud_closure_lu.py`

审计结果：

| closure | Texit [K] | n_iter | rms | eps_b(top) | eps_d(top) | solid(top) | u_d/u_mf(top) | R2(top) | R7(top) | R8(top) | CH4(top) | H2O(top) | COdry | CO2dry | H2dry |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `current` | 1182.0 | 12 | `4.083e-02` | 0.0100 | 0.9900 | 0.0100 | 31.694 | 682.597 | -50.883 | -0.237 | 2.431 | -28.586 | 0.1445 | 0.1008 | 0.1550 |
| `umf_over_epsmf` | 1372.0 | 52 | `9.951e-03` | 0.5317 | 0.5480 | 0.4520 | 2.222 | 0.000 | -0.007 | 1.854 | 0.007 | -1.847 | 0.0573 | 0.1747 | 0.0412 |
| `backsolve_visible_epsb` | 1204.5 | 18 | `2.803e-02` | 0.5596 | 0.4500 | 0.5500 | 0.512 | 42.773 | -19.332 | -4.288 | 10.333 | 0.753 | 0.0830 | 0.1442 | 0.1633 |

### 读数解释

这轮结论已经足够强，可以锁两点：

1. **当前 `u_d ~= u0` 的 closure 本身就是第一顺位问题**
   - 只要把 `u_d` 换成别的合理 closure，top-bed `eps_d_voidage` 就能从 `0.99` 直接回到 `0.45~0.55`
   - 说明“top-bed 过稀”首先不是 chemistry 推出来的，而是 hydrodynamics closure 推出来的

2. **bed chemistry 对这条 closure 极其敏感**
   - `R2/R7/R8/CH4/H2O` 都会跟着 `u_d` 大幅漂移
   - 因此当前 `current` closure 下对 top-bed chemistry 的物理解读不能直接当真

其中更有价值的一条是：

- `backsolve_visible_epsb`
  - 不会像 `umf_over_epsmf` 一样把炉况推到极热的另一条支路
  - 但它已经足以把 top-bed 固相分率恢复到 `~55%`
  - 同时把 `R2` 从 `682.6` 压到 `42.8`
  - 把 `R8` 从弱逆向改成显著逆向
  - 把 `CH4_net` 从 `+2.43` 拉到 `+10.33`

这说明：

- 当前 top-bed `H2O deficit` 很大一部分可能不是“真实 chemistry 必然结果”
- 而是错误 hydrodynamics 把 char-steam contact basis 放大了

### 代码与测试

新增/更新：

- `src/core/cell_hydrodynamics.py`
- `src/core/cell.py`
- `src/core/reactor.py`
- `tests/test_cell_balances.py`
- `tests/test_global_nr_solver.py`
- `scripts/audit_hydrodynamics_ud_closure_lu.py`

验证：

- `python3 -m py_compile src/core/cell_hydrodynamics.py src/core/cell.py src/core/reactor.py scripts/audit_hydrodynamics_ud_closure_lu.py tests/test_cell_balances.py tests/test_global_nr_solver.py`
- `python3 -m pytest tests/test_cell_balances.py::test_ud_closure_alternatives_prevent_visible_bubble_collapse tests/test_global_nr_solver.py::test_reactor_propagates_ud_closure_to_cells -q`
- `python3 scripts/audit_hydrodynamics_ud_closure_lu.py`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_table2_LU_global_nr.py tests/test_global_nr_solver.py tests/test_cell_balances.py tests/test_cell_kinetics.py -q`

结果：

- targeted tests: `2 passed`
- audit 脚本已跑通
- broader regression: `31 passed`

### 当前建议

下一步不该先回头调 chemistry，也不该继续围绕 freeboard 小修小补。  
最该做的是：

1. 把 `u_d` 的 source-of-truth 再核到 Hamel 原文
2. 在确认前，先用 `backsolve_visible_epsb` 作为 **audit branch**
3. 复查 top-bed `eps_b/eps_d/R2/R7/R8/CH4/H2O`
4. 再决定 `alpha_b / K_bd / d_b(h)` 哪一条是下一个主偏差

### 补充：`alpha_b` / `K_bd` regime 审计

新增：

- `scripts/audit_mass_transfer_regime_lu.py`

目标：

- 在同一 LU 工况下，对比两套 `alpha_b` 定义：
  - `alpha_b = u_b / u_mf`
  - `alpha_b = (u_b / u_d) * eps_b`
- 同时并排看三种 `K_bd` 量级：
  - 当前实现 `Eq.3.50` mixed form
  - fast-bubble harmonic branch（`Kbc/Kcd`）
  - slow-bubble Preto branch（`chi=1`）

#### 结果 1：`current` closure 下两套 `alpha_b` 定义全床完全冲突

`current` 口径：

- `alpha_umf`：全床 `44 ~ 56`，全部判成 `fast`
- `alpha_ud = (u_b/u_d)*eps_b`：全床 `0.014 ~ 0.018`，全部判成 `slow`

也就是说：

- **10/10 cells** 出现 regime disagreement

这是一个很强的信号，说明当前 `u_d ~= u0` 不只是让 `eps_b/eps_d` 失真，还会把 `alpha_b` 判别本身拆成两套彼此矛盾的结果。

#### 结果 2：`backsolve_visible_epsb` 下 regime 冲突消失

`backsolve_visible_epsb` 口径：

- `alpha_umf`：全床 `45 ~ 57`，`fast`
- `alpha_ud`：全床 `60 ~ 72`，同样 `fast`

即：

- **0/10 cells** disagreement

这说明至少在“判别是否 fast bubble”这件事上，`backsolve_visible_epsb` 比当前口径更自洽。

#### 结果 3：`K_bd` 量级被 `u_d` closure 改写了 1–2 个数量级

top-bed（cells 6–9）：

`current`

- `Kbd_cur ≈ 8.37 ~ 9.05 1/s`
- `Kbd_fast ≈ 0.066 ~ 0.074 1/s`
- `Kbd_slow ≈ 0.482 ~ 0.497 1/s`

`backsolve_visible_epsb`

- `Kbd_cur ≈ 0.211 ~ 0.224 1/s`
- `Kbd_fast ≈ 0.066 ~ 0.075 1/s`
- `Kbd_slow ≈ 0.473 ~ 0.488 1/s`

含义很清楚：

1. 当前 `current` 实现的 `K_bd` 远高于 branch-based 估计
2. 仅仅把 `u_d` 改成更自洽的 closure，`Kbd_cur` 就会从 `~9` 直接降到 `~0.21`
3. 因此当前主偏差不只是 regime classification，**through-flow basis (`u_d -> u_br`) 本身**也在放大相间交换

### 这一步的结论

到这里可以把 bed hydrodynamics 的问题再往前锁一步：

1. 第一根因仍是 `u_d` closure
2. 但它的后果已经扩展到：
   - `eps_b / eps_d_voidage`
   - `alpha_b` regime 判别
   - `u_br`
   - `K_bd`
3. 所以当前 `R2/R7/R8/CH4/H2O` 的偏差，不应先当成单纯 kinetics 问题

### 验证

- `python3 -m py_compile scripts/audit_mass_transfer_regime_lu.py`
- `python3 scripts/audit_mass_transfer_regime_lu.py`

### 补充：`backsolve_visible_epsb` 对整炉 KPI 的净效果

新增：

- `scripts/audit_ud_closure_full_reactor_lu.py`

目的：

- 不只看局部 hydrodynamics / top-bed chemistry
- 直接比较 `current` 与 `backsolve_visible_epsb` 在 **bed exit / reactor exit / axial profiles** 上的净效果

#### 结果

| closure | Texit [K] | Tbed [K] | n_iter | rms | CO | CO2 | H2 | CH4 | MAE_T | MAE_CO | MAE_H2O |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `current` | 1199.1 | 1182.0 | 12 | `3.987e-02` | 0.1423 | 0.1029 | 0.1550 | 0.0000 | 178.25 | 0.0209 | 0.1297 |
| `backsolve_visible_epsb` | 1216.1 | 1197.2 | 16 | `3.002e-02` | 0.1183 | 0.1265 | 0.1599 | 0.0014 | 191.60 | 0.0311 | 0.1463 |

出口干基相对误差：

`current`

- `CO`: `9.46%`
- `CO2`: `6.49%`
- `H2`: `29.19%`
- `CH4`: `100%`

`backsolve_visible_epsb`

- `CO`: `9.00%`
- `CO2`: `14.97%`
- `H2`: `33.28%`
- `CH4`: `95.09%`

#### 解读

这一步的结论很重要：

1. `backsolve_visible_epsb` **确实修正了 hydrodynamics 自洽性**
   - `eps_b / eps_d / alpha_b / K_bd` 都比 `current` 更像一个内部一致的两相模型

2. 但它**不能单独升级成正式主路径**
   - 虽然 `rms` 略有改善
   - 但整炉 `Texit`、`MAE_T`、`CO2/H2/H2O` 并没有一起变好

3. 因此当前偏差不是“只要修 `u_d` 就结束”
   - 更可能是 `u_d`
   - `K_bd`
   - `d_b(h)`
   - `lambda_b`
   这几项在当前实现里存在**协同偏差**

也就是说，当前最合理的判断是：

- `u_d` 仍然是**第一根因入口**
- 但下一步不能只拿 `backsolve_visible_epsb` 直接替代默认公式
- 而应把它作为 **audit branch**，继续联动审 `K_bd / d_b / lambda_b`

### 验证

- `python3 -m py_compile scripts/audit_ud_closure_full_reactor_lu.py`
- `python3 scripts/audit_ud_closure_full_reactor_lu.py`

### 补充：`backsolve_visible_epsb` 下的 bubble-chain 审计

新增：

- `scripts/audit_bubble_chain_backsolve_lu.py`

目的：

- 固定 `backsolve_visible_epsb` 这条更自洽的相分率口径
- 只拆 `d_b(h)` / `xi_b` / `lambda_b` 三条 bubble-growth 链
- 不再引入出口 KPI 作为主判据

对比了四种 `d_b(h)` 口径：

1. `MW`
   - 当前主路径 `Mori-Wen`
2. `ODE_fix_cur`
   - `xi_b = 0.35`
   - `lambda_b = current = d_b/(0.5u_b)*(P/P0)^(-0.2)`
3. `ODE_dyn_cur`
   - `xi_b = 0 for fast bubble`
   - `lambda_b = current`
4. `ODE_dyn_hamel`
   - `xi_b = 0 for fast bubble`
   - `lambda_b = 280*u_mf/g*(P/P0)^(-0.7)`

#### 结果

在 `backsolve_visible_epsb` 分支上：

- `MW` 与 `ODE_fix_cur` 仍在同一量级
  - top-bed `d_b ≈ 0.40 ~ 0.43 m`
- 但一旦把 `xi_b` 切成 dynamic fast-bubble 口径（`xi_b=0`）
  - 不论用 `current lambda` 还是 `hamel lambda`
  - `d_b` 都几乎塌到极小值
  - `ODE_dyn_cur ≈ 1e-4 m`
  - `ODE_dyn_hamel ≈ 0.039 ~ 0.054 m`

top-bed 代表值（cells 6–9）：

| cell | `d_b cell` | `MW` | `ODE_fix_cur` | `ODE_dyn_cur` | `ODE_dyn_hamel` |
|---|---:|---:|---:|---:|---:|
| 6 | 0.4027 | 0.3934 | 0.0001 | 0.0001 | 0.0393 |
| 7 | 0.4156 | 0.4184 | 0.0001 | 0.0001 | 0.0389 |
| 8 | 0.4234 | 0.4184 | 0.0001 | 0.0001 | 0.0389 |
| 9 | 0.4299 | 0.4322 | 0.0001 | 0.0001 | 0.0388 |

注：

- `ODE_fix_cur` 的第一格与 `MW` 对齐，但后面很快塌到 `1e-4`
- 这说明 bubble-growth ODE 里的参数组合目前并不稳定

#### 解释

这轮结论很关键：

1. `u_d` 不是唯一剩余问题
   - 修完 `u_d` 之后，bubble-growth ODE 的 `xi_b / lambda_b` 立刻浮出来

2. 当前 `Mori-Wen` 主路径之所以还维持在 `0.4 m` 量级
   - 更像是一个经验稳定替代
   - 不是已经证明与 Hamel 的 ODE fully aligned

3. 当前真正需要核的 source-of-truth 已经进一步收窄到：
   - `xi_b` 是否真的应在 fast bubble 下取 `0`
   - `lambda_b` 到底应采用哪一条原文口径
   - ODE 中 `eps_b` 应该用哪一条定义参与 growth term

也就是说，下一步 bed hydrodynamics 的主线已经不是再回头追 `u_d`，而是：

- **核 `xi_b / lambda_b / eps_b(ODE)` 三件事**

### 验证

- `python3 -m py_compile scripts/audit_bubble_chain_backsolve_lu.py`
- `python3 scripts/audit_bubble_chain_backsolve_lu.py`

### 补充：bubble ODE 各项拆解

新增：

- `scripts/audit_bubble_ode_terms_backsolve_lu.py`

目的：

- 在 `backsolve_visible_epsb` 分支上
- 不再只看积分后的 `d_b(h)`
- 而是直接拆：
  - growth term
  - decay term
  - `eps_b` 定义
  - `xi_b`
  - `lambda_b`

#### 关键结果

top-bed 与全床趋势一致：

- `eps_vis` 与 `eps_exc` 几乎相同
  - 例如 cell 6–9：
    - `eps_vis ≈ 0.553 ~ 0.566`
    - `eps_exc ≈ 0.553 ~ 0.566`

- 因此：
  - `g_dyn_vis` 与 `g_dyn_exc` 也几乎完全相同
  - 说明 ODE 崩塌**主因不在 `eps_b` 定义**

同时，growth / decay 的量级关系是：

- `g_fix ≈ 0.090 ~ 0.102`
- `g_dyn ≈ 0.058 ~ 0.063`
- `dec_cur ≈ 0.316`（几乎全床常数）
- `dec_ham ≈ 0.208 ~ 0.416`

于是：

- `dd_fix = g_fix - dec_cur` 全床都为负
- `dd_dyn_ham = g_dyn - dec_ham` 全床也都为负，而且更负

#### 这一步最重要的新发现

当前 `lambda_b(current)` 会把 decay term 变成几乎**与 `d_b`、`u_b` 无关的常数**：

因为代码中

- `lambda_b = d_b / (0.5*u_b) * (P/P0)^(-0.2)`

代回

- `decay = d_b / (3*lambda_b*u_b)`

可化成

- `decay = (0.5/3) * (P/P0)^(0.2)`

也就是说，当前 ODE 里的 split/decay 项在给定压力下几乎就是一个固定值。  
在 LU 的 `2.5 MPa` 工况下，这个固定值约 `0.316`，明显高于当前 growth term（`0.058 ~ 0.102`），于是 `dd_b/dh` 自然全程为负。

这使得问题进一步收窄：

1. **`eps_b` 定义不是当前 ODE 崩塌主因**
2. **`xi_b` 从 0.35 切到 dynamic fast-bubble 0 会进一步压缩 growth term**
3. **当前 `lambda_b(current)` 的代数形式会把 decay 锁成高压下的强常数项**

### 当前判断

因此下一步最该核的 source-of-truth 已经非常明确：

1. `lambda_b` 是否真的应采用当前这条会导致常数 decay 的形式
2. `xi_b` 在 fast bubble 下是否真的取 `0`
3. 若两者都取 Hamel 提取口径，growth/decay 是否还能给出非负或接近平衡的 `dd_b/dh`

### 验证

- `python3 -m py_compile scripts/audit_bubble_ode_terms_backsolve_lu.py`
- `python3 scripts/audit_bubble_ode_terms_backsolve_lu.py`

### 补充：`lambda_b` source-of-truth 已上升为下一优先级

本轮做了两件事：

1. 把 `lambda_b` 做成显式策略，但默认行为不变
   - `src/physics/bubble_dynamics.py`
   - 支持：
     - `current`
     - `hamel_280`
2. 新增隔离审计
   - `scripts/audit_lambda_strategy_backsolve_lu.py`

#### 新结果

在 `backsolve_visible_epsb` 分支、固定 `xi_b=0.35` 的前提下：

| cell | `d_b cell` | `ODE_cur` | `ODE_ham` | `lam_cur` | `lam_ham` |
|---|---:|---:|---:|---:|---:|
| 6 | 0.4027 | 0.0001 | 0.0723 | 0.1710 | 0.1323 |
| 7 | 0.4156 | 0.0001 | 0.0719 | 0.1747 | 0.1336 |
| 8 | 0.4234 | 0.0001 | 0.0719 | 0.1776 | 0.1351 |
| 9 | 0.4299 | 0.0001 | 0.0719 | 0.1797 | 0.1367 |

top-bed 均值：

- `current` ODE: `0.0001 m`
- `hamel_280` ODE: `0.0719 m`

#### 解释

这一步的意义很明确：

1. 在固定 `xi_b=0.35` 时，**单独**把 `lambda_b` 从 `current` 换成 `hamel_280`
   - 已经能把 bubble ODE 从几乎塌到 0 的状态拉回到非零量级

2. 它虽然还没回到当前 `MW / cell d_b ≈ 0.4 m` 的量级
   - 但已经足够说明：
   - **`lambda_b` 本身就是下一优先级 source-of-truth**

3. 因而当前 bubble-chain 的剩余问题可进一步拆成：
   - `lambda_b` 把 ODE 从 `1e-4` 拉回 `~0.07`
   - `xi_b` / 经验 `Mori-Wen` / 其余 growth-term 口径 再决定为什么还没到 `~0.4`

也就是说，下一步最值得做的不是继续扫 `eps_b`，而是继续围绕：

- `lambda_b`
- `xi_b`
- `MW vs ODE growth`

### 代码与测试

新增/更新：

- `src/physics/bubble_dynamics.py`
- `tests/test_bubble_dynamics.py`
- `scripts/audit_lambda_strategy_backsolve_lu.py`

验证：

- `python3 -m py_compile src/physics/bubble_dynamics.py tests/test_bubble_dynamics.py scripts/audit_lambda_strategy_backsolve_lu.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q`
- `python3 scripts/audit_lambda_strategy_backsolve_lu.py`
- `python3 tests/sanity_checks.py`
- `python3 -m pytest tests/test_phase_gate_30.py tests/test_cell_kinetics.py tests/test_bubble_dynamics.py -q`

结果：

- `tests/test_bubble_dynamics.py`: `3 passed`
- `tests/sanity_checks.py`: `ALL 11 SANITY CHECKS PASSED`
- targeted regression: `7 passed`

### 补充：`xi_b` 影响已被隔离，优先级落后于 `lambda_b`

本轮把 `xi_b` 也做成了显式策略，但默认行为不变：

- `fixed_035`
- `hamel_regime`

对应代码：

- `src/physics/bubble_dynamics.py`
- `src/physics/__init__.py`
- `tests/test_bubble_dynamics.py`

新增隔离审计：

- `scripts/audit_xi_strategy_backsolve_lu.py`

#### 结果

固定 `lambda_b = hamel_280` 后，比较：

- `xi_b = 0.35`
- `xi_b = 0 for fast bubbles`（`hamel_regime`，在当前 LU 分支上全床都是 fast）

top-bed（cells 6–9）：

| cell | `d_b cell` | `ODE_fix` | `ODE_dyn` |
|---|---:|---:|---:|
| 6 | 0.4027 | 0.0723 | 0.0393 |
| 7 | 0.4156 | 0.0719 | 0.0389 |
| 8 | 0.4234 | 0.0719 | 0.0389 |
| 9 | 0.4299 | 0.0719 | 0.0388 |

均值：

- `fixed_035`: `0.0719 m`
- `hamel_regime`: `0.0389 m`

#### 解读

这一步已经足够把剩余优先级排清楚：

1. `lambda_b` 的影响更大
   - 它能把 ODE top-bed `d_b` 从 `~1e-4` 拉到 `~0.072`

2. `xi_b` 也有实质影响
   - 在 `lambda_b=hamel_280` 固定后，再把 `xi_b` 切成 `hamel_regime`
   - 会把 `d_b` 从 `~0.072` 进一步压到 `~0.039`

3. 因而当前 bubble-chain 的优先级可排序为：
   - **`lambda_b`**
   - **`xi_b`**
   - 然后才是其余 growth-term 细节

也就是说，当前最合理的判断是：

- `lambda_b` 是 bubble ODE 的下一主偏差
- `xi_b` 是紧随其后的次级偏差
- 两者都还没回到能解释当前 `MW / cell d_b ≈ 0.4 m` 的程度

### 验证

- `python3 -m py_compile src/physics/bubble_dynamics.py src/physics/__init__.py tests/test_bubble_dynamics.py scripts/audit_xi_strategy_backsolve_lu.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q`
- `python3 scripts/audit_xi_strategy_backsolve_lu.py`

结果：

- `tests/test_bubble_dynamics.py`: `4 passed`

### 补充：`d_b0` 不是当前 bubble-gap 的主因

新增隔离审计：

- `scripts/audit_db0_strategy_backsolve_lu.py`

目的：

- 在 `backsolve_visible_epsb` 分支上
- 固定更接近论文的 bubble ODE 审计口径
  - `lambda_b = hamel_280`
  - `xi_b = fixed_035`
- 只替换 `d_b0`

对比了三种 `d_b0`：

1. `current(Darton)`
2. `Hamel(q0/orifice)`
   - `d_b0 = 1.3 * (Vdot0^2/g)^0.2`
   - `Vdot0` 以 **总表观体积流 / 孔数** 推断
3. `Hamel(qex/orifice)`
   - 同式
   - 但 `Vdot0` 以 **excess gas = (u0-u_mf)` / 孔数** 推断
   - 仅作 audit，不是已锁定 source-of-truth

#### 结果

`d_b0` 本身有差别：

- `current(Darton) = 0.11313 m`
- `Hamel(q0/orifice) = 0.09096 m`
- `Hamel(qex/orifice) = 0.08984 m`

但积分到 top-bed 后，三条几乎重合：

| cell | `d_b cell` | `ODE_cur` | `ODE_q0` | `ODE_qex` |
|---|---:|---:|---:|---:|
| 6 | 0.4027 | 0.0723 | 0.0720 | 0.0720 |
| 7 | 0.4156 | 0.0719 | 0.0719 | 0.0719 |
| 8 | 0.4234 | 0.0719 | 0.0719 | 0.0719 |
| 9 | 0.4299 | 0.0719 | 0.0718 | 0.0718 |

top-bed 均值：

- `current = 0.0719 m`
- `hamel_q0 = 0.0719 m`
- `hamel_qex = 0.0719 m`

#### 解读

这一步已经足够把 `d_b0` 从主嫌疑里基本排除：

1. `d_b0` 的不同取法只影响床底最初几格
2. 到 top-bed 后，解几乎完全收敛到同一条 `~0.072 m` 轨道
3. 因此当前 `ODE (~0.07 m)` 与 `MW / cell (~0.4 m)` 的巨大差距，**主因不在 `d_b0`**

于是 bubble-chain 的剩余优先级可以再收紧为：

1. `lambda_b`
2. `xi_b`
3. growth/decay 主体本身（而不是 `d_b0`）

### 验证

- `python3 -m py_compile scripts/audit_db0_strategy_backsolve_lu.py`
- `python3 scripts/audit_db0_strategy_backsolve_lu.py`

### 补充：常参 ODE vs 局部 ODE 的口径差异很小

新增审计：

- `scripts/audit_bubble_local_ode_backsolve_lu.py`

目的：

- 验证当前 `ODE (~0.07 m)` 与 `MW/cell (~0.4 m)` 的 gap
- 是否只是因为前面一直拿 **床底固定 `u0/u_mf` 的常参 ODE** 去对比 **各 cell 局部条件下的 MW / cell d_b**

做法：

- 固定 bubble ODE 审计口径：
  - `hydrodynamics_u_d_closure = backsolve_visible_epsb`
  - `lambda_b = hamel_280`
  - `xi_b = fixed_035`
- 对比：
  1. `constant-bottom ODE`
  2. `piecewise-local ODE`（每个 cell 用本地 `u0/u_mf` 走一段）

#### 结果

| cell | `d_b cell` | `ODE_const` | `ODE_local` |
|---|---:|---:|---:|
| 6 | 0.4027 | 0.0721 | 0.0727 |
| 7 | 0.4156 | 0.0719 | 0.0734 |
| 8 | 0.4234 | 0.0719 | 0.0739 |
| 9 | 0.4299 | 0.0719 | 0.0747 |

top-bed 均值：

- `constant-bottom ODE = 0.0719 m`
- `piecewise-local ODE = 0.0737 m`
- `cell / MW = 0.4179 m`

#### 解读

这一步已经足够说明：

1. `piecewise-local ODE` 只比 `constant-bottom ODE` 略大一点
2. 因而当前 `ODE (~0.07 m)` 与 `MW/cell (~0.4 m)` 的 gap
   - **基本不在 `u0/u_mf` 沿程变化**
3. 换句话说，之前的常参 ODE 对比虽然更粗，但**不是主误差来源**

于是剩余主问题可以继续收紧为：

- `lambda_b`
- `xi_b`
- growth/decay 主体本身

而不是：

- `d_b0`
- 或 `u0/u_mf` 的沿程更新

### 验证

- `python3 -m py_compile scripts/audit_bubble_local_ode_backsolve_lu.py`
- `python3 scripts/audit_bubble_local_ode_backsolve_lu.py`

### 补充：反推所需 closure 后，`lambda_b` 优先级进一步坐实

新增审计：

- `scripts/audit_bubble_required_closure_backsolve_lu.py`

目的：

- 在 `backsolve_visible_epsb` 分支上
- 直接反推：
  - 若保持 `xi_b = 0.35`，为了让当前 `cell d_b` 满足 `dd_b/dh ≈ 0`，需要多大的 `lambda_b`
  - 若保持当前 `lambda_b`，为了满足 `dd_b/dh ≈ 0`，需要多大的 `xi_b`

#### 结果

代表性的 top-bed（cells 6–9）：

| cell | `lam_cur` | `lam_ham` | `lam_req` | `req/cur` | `req/ham` | `xi_req(cur)` | `xi_req(ham)` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | 0.1710 | 0.1323 | 0.5924 | 3.46 | 4.48 | 0.794 | 0.835 |
| 7 | 0.1747 | 0.1336 | 0.6075 | 3.48 | 4.55 | 0.796 | 0.839 |
| 8 | 0.1776 | 0.1351 | 0.6210 | 3.50 | 4.60 | 0.800 | 0.843 |
| 9 | 0.1797 | 0.1367 | 0.6300 | 3.51 | 4.61 | 0.802 | 0.845 |

全床中位数：

- `lam_req / current ≈ 3.42`
- `lam_req / hamel_280 ≈ 4.25`
- `xi_req(current) ≈ 0.786`
- `xi_req(hamel_280) ≈ 0.821`

#### 解读

这一步把优先级进一步锁死了：

1. 若坚持 `xi_b = 0.35`
   - 那么要支撑当前 `cell d_b≈0.4 m`
   - `lambda_b` 需要普遍放大到当前的 **3.4x**
   - 相对 `hamel_280` 甚至要 **4.2x**

2. 反过来，若坚持当前 `lambda_b`
   - 为了让 ODE 平衡到当前 `cell d_b`
   - 需要的 `xi_b` 会接近 `0.8`
   - 这已经远高于当前项目默认的 `0.35`

3. 因而当前最合理的优先级判断是：
   - **主偏差更像 `lambda_b` 太小**
   - 而不是仅仅 `xi_b` 太低

也就是说，bubble-chain 剩余问题现在可以再收敛成：

1. `lambda_b`
2. `xi_b`
3. 然后才是更细的 growth-term 结构

### 验证

- `python3 -m py_compile scripts/audit_bubble_required_closure_backsolve_lu.py`
- `python3 scripts/audit_bubble_required_closure_backsolve_lu.py`

### 补充：hydrodynamics pressure correlations 逐项核对

新增审计：

- `scripts/audit_hydrodynamics_pressure_correlations.py`

目的：

- 把 hydrodynamics 相关 pressure-sensitive 项统一过一遍
- 区分“公式正确且与本地 source-of-truth 一致”与“仍 unresolved”
- 量化 `101300` / `101325` 的基准混用量级

#### 结果

在 `P_ref = 101325 Pa`、`P_LU = 2.5 MPa` 下：

- `rho_g`
  - code ratio = `24.6731`
  - 与理想气体 `rho_g ~ P` 完全一致

- `D_g`
  - code ratio = `0.04053`
  - 与 `D_g ~ 1/P` 一致
  - 但公式内部用 `101300 Pa`
  - 而项目常数 `P0 = 101325 Pa`
  - 两者偏差仅 `0.0247%`

- `u_br`
  - code ratio = `0.618253`
  - 与 `u_br ~ (P/P0)^(-0.15)` 完全一致

- `lambda_b`
  - `current` ratio = `0.52669`（`(P/P0)^(-0.2)`）
  - `hamel_280` ratio = `0.106034`（`(P/P0)^(-0.7)`）
  - 说明实现层面两条口径都能正确算
  - 但 **source-of-truth 仍未统一**

- `K_bd`
  - 对流项 ratio = `0.618253`
  - 扩散项 ratio = `0.201321`
  - 总体 ratio = `0.499448`
  - 符合其分别继承 `u_br` 与 `sqrt(D_g)` 的 pressure path

- `beta_A / u_gb`
  - 当前都 **没有显式 pressure exponent**
  - 仅通过上游 `d_b`、`rho_g`、terminal-velocity 链间接受压强影响

#### 解读

这轮 pressure audit 把问题收缩成一条主线：

1. **正确且基本一致**
   - `rho_g`
   - `u_br`
   - `K_bd` 的 pressure path

2. **公式正确，但常数基准未完全统一**
   - `D_g` 仍在公式里写 `101300`
   - 而项目常数 `P0 = 101325`
   - 量级影响极小，但从 scientific-audit 口径看仍值得后续统一

3. **真正 unresolved 的 pressure correction**
   - 仍是 `lambda_b`

4. **当前没有本地证据支持**
   - `beta_A` 或 `u_gb` 还缺一条额外显式 pressure correlation

也就是说，hydrodynamics 相关 pressure correlations 里：

- **最大的不一致源不是 `u_br` 或 `K_bd`**
- **而仍然是 `lambda_b` 的 pressure correction source-of-truth**

### 验证

- `python3 -m py_compile scripts/audit_hydrodynamics_pressure_correlations.py`
- `python3 scripts/audit_hydrodynamics_pressure_correlations.py`
- `python3 tests/sanity_checks.py`

### 补充：统一 `D_g` 的 pressure baseline

这轮做了一处可以安全落地的小修正：

- [src/core/species.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/core/species.py)
  - `gas_diffusivity_correlation()` 不再硬编码 `101300/P`
  - 改为统一使用 `core.constants.P0 / P`
- [tests/test_phase1.py](/Users/liuzhen/AI-projects/bfb-gasifier/tests/test_phase1.py)
  - 对应测试改成 `P0` / `2*P0`
- [docs/source_units_audit.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/source_units_audit.md)
  - 说明同步更新
- [scripts/audit_hydrodynamics_pressure_correlations.py](/Users/liuzhen/AI-projects/bfb-gasifier/scripts/audit_hydrodynamics_pressure_correlations.py)
  - 现在会把 `D_g` 判成“已统一”

#### 解读

这不会改变 `D_g ~ 1/P` 的公式形状，也不会实质改变数值结果；它只是把原来 `101300` 与项目常数 `P0=101325` 的双基准清掉。  
这样现在 hydrodynamics 相关 pressure correlations 里：

- `rho_g`：一致
- `D_g`：一致且基准已统一
- `u_br`：一致
- `K_bd`：一致
- **剩余真正 unresolved 的 pressure correction 只剩 `lambda_b`**

#### 验证

- `python3 -m py_compile src/core/species.py tests/test_phase1.py scripts/audit_hydrodynamics_pressure_correlations.py`
- `python3 -m pytest tests/test_phase1.py -q` → `37 passed`
- `python3 scripts/audit_hydrodynamics_pressure_correlations.py`
- `python3 tests/sanity_checks.py`

### 补充：`lambda_b` pressure-window 审计

新增审计：

- `scripts/audit_lambda_b_pressure_window.py`

目的：

- 在 `0.1–2.5 MPa` 的代表性压力窗口内
- 直接比较 `lambda_b(current)` 与 `lambda_b(hamel_280)`
- 同时比较它们在 ODE 衰减项 `d_b/(3*lambda_b*u_b)` 上造成的量级差异

#### 结果

固定：

- `d_b = 0.12 m`
- `u_b = 1.2 m/s`
- `u_mf = 0.05 m/s`

得到：

| P [MPa] | `lam_cur` [s] | `lam_ham` [s] | `ham/cur` | `dec_cur` [1/m] | `dec_ham` [1/m] | `dec_ham/dec_cur` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.101 | 0.2000 | 1.4271 | 7.14 | 0.1667 | 0.0234 | 0.140 |
| 0.500 | 0.1453 | 0.4669 | 3.21 | 0.2294 | 0.0714 | 0.311 |
| 1.000 | 0.1265 | 0.2874 | 2.27 | 0.2635 | 0.1160 | 0.440 |
| 1.500 | 0.1167 | 0.2164 | 1.85 | 0.2857 | 0.1541 | 0.539 |
| 2.000 | 0.1102 | 0.1769 | 1.61 | 0.3026 | 0.1884 | 0.623 |
| 2.500 | 0.1053 | 0.1513 | 1.44 | 0.3164 | 0.2203 | 0.696 |

#### 解读

这一步给了一个很重要的新约束：

1. **不能简单以指数判断强弱**
   - 虽然 `hamel_280` 的 pressure exponent 是更强的 `-0.7`
   - 但因为它在常压下的基值 `280*u_mf/g` 很大
   - 所以在 `2.5 MPa` 下它给出的 `lambda_b` **仍然比 current 更大**
   - 只是差距已从常压的 `7.1x` 缩到 `1.44x`

2. 也就是说：
   - 若 `ODE_dyn_hamel` 仍然显著压小 `d_b`
   - **主因不能只归咎于 pressure exponent `-0.7`**
   - 还必须联动看 `u_mf` 基值、`xi_b`、以及 growth/decay 主体

3. 因而当前对 `lambda_b` 的更精确判断应改写为：
   - unresolved 点不只是“指数是 `-0.2` 还是 `-0.7`”
   - 而是**整条 `lambda_b` 构造式的 source-of-truth**

### 验证

- `python3 -m py_compile scripts/audit_lambda_b_pressure_window.py`
- `python3 scripts/audit_lambda_b_pressure_window.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q` → `4 passed`

### 补充：`lambda_b` source-of-truth 已锁定并切进默认实现

基于对 Hamel (1999) 原文的进一步核对：

- Eq. 3.35, p.30：
  - `lambda_b = 280 * u_mf / g`
- Eq. 3.42, p.32（after Heinbockel 1995）：
  - `lambda_b = 280 * u_mf / g * (P/P0)^(-0.7)`
  - `P0 = 101300 Pa`

并且论文正文明确说明：

- `lambda_b` 是**外部经验关联**
- **不是**由 `d_b/u_b` 推导出的内部变量

#### 已落实到实现的改动

- [src/core/constants.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/core/constants.py)
  - 新增 `P0_HAMEL = 101300 Pa`
- [src/physics/bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/bubble_dynamics.py)
  - `bubble_lifetime()` 默认策略切为 `hamel_280`
  - 旧 `current` 保留为 legacy/audit 对照
- [src/physics/mass_transfer.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/mass_transfer.py)
  - `calc_u_br()` 默认参考压切到 `P0_HAMEL`
- [src/core/species.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/core/species.py)
  - `D_g` 的 pressure baseline 也统一到 `P0_HAMEL`
- [specs/02_hydrodynamics.md](/Users/liuzhen/AI-projects/bfb-gasifier/specs/02_hydrodynamics.md)
  - 明确补入 Eq. 3.35 / Eq. 3.42
- [docs/techspec.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/techspec.md)
  - 同步补入论文页码与方程号

#### 解读

这一步之后，hydrodynamics 相关 pressure correlations 的状态变成：

- `rho_g`：一致
- `D_g`：一致
- `u_br`：一致
- `K_bd`：一致
- `lambda_b`：**已按 Hamel 原文锁定**

也就是说，之前 “`lambda_b` 是 hydrodynamics 中唯一 unresolved pressure correction” 这条判断现在已经失效。  
后面如果 `bubble ODE` 仍与当前 `MW/cell d_b` 差距很大，问题就该继续回到：

1. `xi_b`
2. growth/decay 主体
3. `MW` 与 `Hilligardt ODE` 的模型结构差异

而不该再把 `lambda_b` 当 source-of-truth 未决项。

#### 验证

- `python3 -m py_compile src/core/constants.py src/core/species.py src/physics/mass_transfer.py src/physics/bubble_dynamics.py tests/test_bubble_dynamics.py tests/test_phase1.py scripts/audit_hydrodynamics_pressure_correlations.py scripts/audit_lambda_b_pressure_window.py scripts/audit_hydrodynamics_consistency_lu.py scripts/audit_bubble_chain_backsolve_lu.py`
- `python3 -m pytest tests/test_phase1.py tests/test_bubble_dynamics.py -q` → `42 passed`
- `python3 tests/sanity_checks.py` → `ALL 11 SANITY CHECKS PASSED`
- `python3 scripts/audit_hydrodynamics_pressure_correlations.py`
- `python3 scripts/audit_lambda_b_pressure_window.py`

### 补充：`xi_b` 原文定义已锁定并切进 ODE 口径

进一步核对 Hamel (1999) Chapter 3.1.2, p.26–27, Eq.3.34 后，`xi_b` 的 source-of-truth 已明确：

- `alpha_b = u_b / u_d`
- `0 < alpha_b < 1`：
  - `xi_b = 1 - alpha_b^3`
- `alpha_b > 1`：
  - `xi_b = 0`

也就是说：

- **不是** `u_b/u_mf`
- 也**不是**连续软过渡到别的值
- 而是按 fast/slow bubble 的**硬切换**

#### 已落实到实现的改动

- [src/physics/bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/bubble_dynamics.py)
  - `classify_bubble_regime()` 改为按 `u_b/u_d`
  - `bubble_diameter_ode()` / `integrate_bubble_diameter()` 新增 `u_d`
  - ODE 默认 `xi_strategy` 切为 `hamel_regime`
- [specs/02_hydrodynamics.md](/Users/liuzhen/AI-projects/bfb-gasifier/specs/02_hydrodynamics.md)
  - 明确补入 Eq.3.34 的分段定义
- [docs/techspec.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/techspec.md)
  - 修正 `alpha_b` 定义为 `u_b/u_d`
- [docs/BFB_TechSpec_v11.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/BFB_TechSpec_v11.md)
  - 同步修正
- [tests/test_bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/tests/test_bubble_dynamics.py)
  - 新增 `classify_bubble_regime_uses_ud_ratio`

#### 新结果

复跑 [audit_xi_strategy_backsolve_lu.py](/Users/liuzhen/AI-projects/bfb-gasifier/scripts/audit_xi_strategy_backsolve_lu.py) 后：

- `u0_ref = 1.4344`
- `u_mf_ref = 0.0435`
- `u_d_ref = 0.0196`

因此当前 `backsolve_visible_epsb` 分支里：

- 全床 `alpha_b = u_b/u_d ≈ 99 ~ 114`
- 所以 `xi_b` 在 `hamel_regime` 下全床都被硬切到 `0`

对应 ODE top-bed：

- `fixed_035`：`d_b ≈ 0.0719 m`
- `hamel_regime`：`d_b ≈ 0.0388 m`

#### 解读

这一步把一个重要结论钉死了：

1. 之前把 `xi_b` 建在 `u_b/u_mf` 上，确实是错口径
2. 按论文原文改成 `u_b/u_d` 后
   - 当前 LU / `backsolve_visible_epsb` 分支会把全床都判成 **fast bubble**
   - 因而 `xi_b = 0` 成为主导结果
3. 这会把 ODE growth term 再往下压一层

也就是说，若后面 `Hilligardt ODE` 仍然显著小于 `MW / cell d_b`
那就更不该回头怀疑 `xi_b` 定义本身了；问题会继续集中到：

- `u_d` 的物理量级
- ODE growth/decay 结构
- `MW` 与 `Hilligardt ODE` 的模型结构差异

#### 验证

- `python3 -m py_compile src/physics/bubble_dynamics.py tests/test_bubble_dynamics.py scripts/audit_xi_strategy_backsolve_lu.py scripts/audit_bubble_chain_backsolve_lu.py scripts/audit_hydrodynamics_consistency_lu.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q` → `6 passed`
- `python3 scripts/audit_xi_strategy_backsolve_lu.py`
- `python3 tests/sanity_checks.py`

### 更正：`chi` 不属于 Eq.3.33，而属于 Preto 传质模型

进一步核对后，上一轮关于 `chi` 的判断已被纠正：

- `xi_b`：属于 bubble-growth ODE（Eq.3.33）
- `chi`：属于 Section 3.1.4 的 Preto 传质模型

也就是说：

- Eq.3.33 的 bubble-growth ODE **不含 `chi`**
- `chi` 不应接入 `dd_b/dh` 的 growth term

#### 已落实到实现

- [src/physics/bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/bubble_dynamics.py)
  - 已把误接入 ODE 的 `chi_b` 完整移除
- [tests/test_bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/tests/test_bubble_dynamics.py)
  - 已移除对应错误测试
- [specs/02_hydrodynamics.md](/Users/liuzhen/AI-projects/bfb-gasifier/specs/02_hydrodynamics.md)
  - ODE 结构已改回不含 `chi`
- [docs/techspec.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/techspec.md)
  - 已改回不含 `chi`
- [docs/BFB_TechSpec_v11.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/BFB_TechSpec_v11.md)
  - 已改回不含 `chi`

#### 正确口径

- `chi` 仍然重要，但它应该出现在 **Preto 的 `K_bd` 模型**
- 不应与 `xi_b` 混入同一条 bubble-growth ODE

所以当前最需要继续从原文确认的 `chi` 问题，应当放在：

- **mass-transfer / Preto branch**

而不是：

- **bubble-growth ODE**

#### 验证

- `python3 -m py_compile src/physics/bubble_dynamics.py tests/test_bubble_dynamics.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q` → `7 passed`
- `python3 tests/sanity_checks.py` → `ALL 11 SANITY CHECKS PASSED`

### 补充：`psi_b` 与 Eq.3.33 主路径属性已锁定

进一步核对 Hamel 原文后，以下几项现在可以视为已确认：

#### `psi_b`

Hamel (1999) p.28：

- Eq.3.14
  - technical gas distributor：`psi_b = 0.76`
  - porous plate：`psi_b = 0.67`
- Eq.3.15（Wein 1992）
  - `psi_b = 0.17 * u_mf^(-0.33)`

项目当前默认仍取 technical gas distributor 的 `0.76`，这与 HTW 工业分布板口径一致。

#### Eq.3.33 / Eq.3.41 的地位

论文在 Chapter 3.1.3 明确说明：

- 通过修改后的 Hilligardt 方程，可扩展到加压床层
- 也就是说，Eq.3.33 / pressure-modified Eq.3.41 是论文的**主 bubble-growth 模型**
- **不是**被 `Darton` / `Mori-Wen` 经验式取代

#### 已落实到实现

- [src/physics/bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/bubble_dynamics.py)
  - 新增 `bubble_interaction_factor()`
  - 支持 `technical_distributor / porous_plate / wein_1992`
- [src/physics/__init__.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/physics/__init__.py)
  - 导出 `bubble_interaction_factor`
- [tests/test_bubble_dynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/tests/test_bubble_dynamics.py)
  - 新增 `test_bubble_interaction_factor_supports_hamel_sources`
- [specs/02_hydrodynamics.md](/Users/liuzhen/AI-projects/bfb-gasifier/specs/02_hydrodynamics.md)
  - 补入 `psi_b` 三条来源
  - 明确写出 Eq.3.33/3.41 是论文主路径
- [docs/techspec.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/techspec.md)
  - 同步补入
- [docs/BFB_TechSpec_v11.md](/Users/liuzhen/AI-projects/bfb-gasifier/docs/BFB_TechSpec_v11.md)
  - 同步补入

#### 解读

这一步把一个更大的结构性结论钉死了：

- 论文主路径是 **Hilligardt ODE**
- 当前项目 bed 主路径仍是 [cell_hydrodynamics.py](/Users/liuzhen/AI-projects/bfb-gasifier/src/core/cell_hydrodynamics.py) 里的 `Mori-Wen`

所以现在最大的 hydrodynamics 结构差距已经不是单个参数了，而是：

1. `u_d`
2. `xi_b`
3. `lambda_b`
4. `chi`

### 补充：Hilligardt ODE 已接入 bed hydrodynamics audit branch

本轮没有切默认主路径；`cell_hydrodynamics` 仍默认走 `Mori-Wen`。但我已把论文口径的 bubble-growth ODE 作为 **可切换 audit branch** 正式接入：

- `src/core/cell_hydrodynamics.py`
  - 新增 `bubble_diameter_model`
  - 支持：
    - `mori_wen`（默认，保持当前行为）
    - `hilligardt_ode`（审计分支）
  - `hilligardt_ode` 会吃：
    - `psi_b_strategy`
    - `lambda_strategy`
    - `xi_strategy`
    - `u_d_closure`
  - 采用短固定点迭代，把 `u_d <-> u_b <-> d_b(h)` 串回同一条 hydrodynamics closure
- `src/core/reactor.py`
  - `ReactorConfig` 新增：
    - `hydrodynamics_bubble_diameter_model`
    - `hydrodynamics_psi_b_strategy`
    - `hydrodynamics_lambda_strategy`
    - `hydrodynamics_xi_strategy`
- `src/core/cell.py`
  - 上述配置已能从 `ReactorConfig -> Cell -> calc_cell_hydrodynamics()` 贯通

对应验证：

- `tests/test_cell_balances.py`
  - 新增 `test_hilligardt_ode_bubble_branch_is_available_for_hydrodynamics_audit`
- `tests/test_global_nr_solver.py`
  - `test_reactor_propagates_ud_closure_to_cells` 现同时检查 ODE hydrodynamics 配置传播
- `python3 -m pytest tests/test_cell_balances.py::test_hilligardt_ode_bubble_branch_is_available_for_hydrodynamics_audit tests/test_global_nr_solver.py::test_reactor_propagates_ud_closure_to_cells tests/test_bubble_dynamics.py -q`
  - 结果：`9 passed`
- `python3 tests/sanity_checks.py`
  - 结果：`ALL 11 SANITY CHECKS PASSED`

新增审计：

- `scripts/audit_bed_bubble_model_lu.py`

它现在按 **Vorabrechnung initialized LU state** 比较 `Mori-Wen` 与 `Hilligardt ODE` 两条 bed hydrodynamics 链，避免求解器 branch-selection 干扰。当前结果：

- `Mori-Wen`
  - `d_b` 自底到顶约 `0.15 -> 0.40 m`
  - `eps_b` 约 `0.68 -> 0.51`
  - `K_bd` 约 `0.65 -> 0.23 1/s`
- `Hilligardt ODE`（`backsolve_visible_epsb + lambda_b=hamel_280 + xi_b=hamel_regime`）
  - `d_b` 自底到顶约 `0.081 -> 0.037 m`
  - `eps_b` **全床被顶到 `0.70` clip**
  - `K_bd` 约 `1.34 -> 3.32 1/s`

当前最重要的判断：

1. 论文主路径 ODE 已能进入 bed 主链，不再只是离线脚本比较。
2. 一旦把 ODE 接进主链，当前剩余偏差已不只是 `u_d`。
3. **`visible bubble fraction / eps_b closure` 已上升为下一主偏差**：
   - 在 initialized LU state 上，ODE 分支会把全床 `eps_b` 顶到 `0.70` clip。
   - 这说明当前 `eps_b = calc_visible_bubble_fraction(u0, u_b, u_d)` 与 ODE/`alpha_b=u_b/u_d` 联动后，系统会被推到另一个极端。
4. 因此下一步最该继续审的，不是 `lambda_b` 或 `chi`，而是：
   - `eps_b` 的 source-of-truth
   - 它在 Hilligardt ODE 里是否应直接由 `u0/u_b/u_d` 这样闭合
   - 以及 `visible bubble fraction` 与床内两相体积分率的关系是否被简化错位

### 补充：`eps_b` 的当前主疑点已收窄到 `u_d` 反算闭环，而不是公式二选一

新增脚本：

- `scripts/audit_eps_b_closure_lu.py`

它在 **Vorabrechnung initialized LU bed state** 上，直接并排计算：

- `eps_excess = (u0 - u_mf) / u_b`
- `eps_visible = (u0 - u_d) / (u_b + (n_b-1)u_d)`

当前结果非常明确：

- 在 `backsolve_visible_epsb` 分支上，无论是 `Mori-Wen` 还是 `Hilligardt ODE`
- 全床 `eps_excess == eps_visible`（数值误差仅 `1e-16` 量级）

这说明：

1. 当前主疑点已经不是“代码到底该用 excess 还是 visible 公式”
2. 而是 **`u_d` closure 本身正在反算并强制满足这条等式**
3. 因而一旦 `u_b` 被 ODE 压小，`eps_b` 就会被整体抬高到接近 `0.9`，再撞上 cell-level clip `0.7`

所以后续最该核的原文问题应改写为：

- `u_d` 是否应该由 visible bubble fraction 反算得到
- 还是 `u_d` / `eps_b` / `eps_d_voidage` 应分别由不同 closure 独立给出

### 补充：`Eq.3.11 / Eq.3.12` 的 `u_d` closure 已接入 audit branch

本轮把原文给出的两条 `u_d` 关联式接入了 `src/core/cell_hydrodynamics.py`：

- `hilligardt_eq311`
  - `u_d = u_mf + (u0 - u_mf)/3`
- `wein_1992_eq312`
  - `u_d = 1.45 * u_mf`

默认主路径仍未切换；当前只是把它们作为 **audit branch** 暴露出来。

对应验证：

- `tests/test_cell_balances.py::test_ud_closure_alternatives_prevent_visible_bubble_collapse`
  - 已新增对 `hilligardt_eq311` / `wein_1992_eq312` 的公式检查
- `python3 -m pytest tests/test_cell_balances.py::test_ud_closure_alternatives_prevent_visible_bubble_collapse -q`
  - 结果：`1 passed`

新增脚本：

- `scripts/audit_ud_closure_initialized_lu.py`

它在 **Vorabrechnung initialized LU state** 上比较五条 `u_d` closure：

- `current`
- `umf_over_epsmf`
- `hilligardt_eq311`
- `wein_1992_eq312`
- `backsolve_visible_epsb`

当前 top-bed (`xi/H >= 0.65`) 结果：

- `current`
  - `eps_b_top = 0.010`
  - `eps_d_voidage_top = 0.990`
  - `solid_top = 0.010`
  - `u_d/u_mf = 27.06`
- `umf_over_epsmf`
  - `eps_b_top = 0.460`
  - `eps_d_voidage_top = 0.551`
  - `solid_top = 0.449`
  - `u_d/u_mf = 2.22`
- `hilligardt_eq311`
  - `eps_b_top = 0.261`
  - `eps_d_voidage_top = 0.801`
  - `solid_top = 0.199`
  - `u_d/u_mf = 9.69`
- `wein_1992_eq312`
  - `eps_b_top = 0.486`
  - `eps_d_voidage_top = 0.495`
  - `solid_top = 0.506`
  - `u_d/u_mf = 1.45`
- `backsolve_visible_epsb`
  - `eps_b_top = 0.519`
  - `eps_d_voidage_top = 0.450`
  - `solid_top = 0.550`
  - `u_d/u_mf = 0.531`

当前最重要的新判断：

1. `current` 口径现在可以进一步降级，它明显把 top-bed 直接打成 `0.99` voidage。
2. `hilligardt_eq311` 虽然来自原文，但在当前 initialized LU state 上仍然让 top-bed 过稀：
   - `u_d/u_mf ≈ 9.7`
   - `eps_d_voidage_top ≈ 0.80`
3. `wein_1992_eq312` 在 initialized hydrodynamics 上给出的量级最像一个仍属于 bubbling bed 的 top-bed：
   - `eps_b_top ≈ 0.49`
   - `eps_d_voidage_top ≈ 0.49`
4. 这说明下一步不该再把 `u_d` 当成单一未决项，而要继续查：
   - Hamel 在 HTW 加压床层主算例里，`u_d` 最终到底采用 Eq.3.11 还是 Eq.3.12
   - Eq.3.11 是否只适用于某一类分布板/操作窗口
   - `backsolve_visible_epsb` 虽然给出更“饱满”的 dense phase，但它不符合原文逻辑顺序，不能直接当 source-of-truth

### 补充：Heinbockel 高压 bubble chain 已接入 audit branch

根据 Hamel (1999) Chapter 3.1.3 的原文，Heinbockel 高压修正并不是“只给 `lambda_b` 加一个指数”，而是**整条 bubble chain 的全局替换**：

- Eq.3.41：pressurized bubble-growth ODE
- Eq.3.42：bubble lifetime
- Eq.3.43：bubble rise velocity
- Eq.3.44：bubble through-flow

本轮代码实现：

- `src/physics/bubble_dynamics.py`
  - `bubble_rise_velocity(..., strategy=...)`
    - `hilligardt_eq313`（默认）
    - `heinbockel_eq343`
  - `bubble_diameter_ode(..., ode_strategy=...)`
    - `hilligardt_eq333`（默认）
    - `heinbockel_eq341`
- `src/core/reactor.py`
  - `ReactorConfig` 新增：
    - `hydrodynamics_bubble_velocity_strategy`
    - `hydrodynamics_bubble_ode_strategy`
- `src/core/cell.py` / `src/core/cell_hydrodynamics.py`
  - 上述配置已贯通到 bed hydrodynamics 主链

默认行为未改；当前仍只作为 audit branch 使用。

对应验证：

- `tests/test_bubble_dynamics.py`
  - 新增：
    - `test_bubble_rise_velocity_supports_heinbockel_pressure_form`
    - `test_integrate_bubble_diameter_supports_heinbockel_ode_branch`
- `tests/test_global_nr_solver.py`
  - `test_reactor_propagates_ud_closure_to_cells` 现同时检查：
    - `hydrodynamics_bubble_velocity_strategy`
    - `hydrodynamics_bubble_ode_strategy`
- `python3 -m pytest tests/test_bubble_dynamics.py tests/test_global_nr_solver.py::test_reactor_propagates_ud_closure_to_cells -q`
  - 结果：`10 passed`

新增审计：

- `scripts/audit_heinbockel_initialized_lu.py`

它固定：

- `u_d = wein_1992_eq312`
- `bubble_diameter_model = hilligardt_ode`

只比较：

- `hilligardt_atm_chain`
  - `velocity_strategy = hilligardt_eq313`
  - `ode_strategy = hilligardt_eq333`
- `heinbockel_pressurized_chain`
  - `velocity_strategy = heinbockel_eq343`
  - `ode_strategy = heinbockel_eq341`

在 **initialized LU state** 上，结果非常关键：

- `hilligardt_atm_chain`
  - `d_b ≈ 0.081 -> 0.039 m`
  - `u_b ≈ 1.32 -> 1.56 m/s`
  - `eps_b = 0.70`（全床撞 clip）
  - `K_bd ≈ 2.67 -> 6.57 1/s`
- `heinbockel_pressurized_chain`
  - `d_b ≈ 0.113 -> 0.207 m`
  - `u_b ≈ 18.8 -> 19.3 m/s`
  - `eps_b ≈ 0.061`
  - `K_bd ≈ 3.03 -> 1.40 1/s`

当前最重要的新判断：

1. 只要把 `u_d` 固定到 `Eq.3.12 (Wein)`，再启用 Heinbockel 高压链，`eps_b` 就不再爆到 `0.7~0.9`。
2. 这说明我们前面在 ODE 分支上看到的 `eps_b` 极端值，确实不是“少了 pressure exponent”这么简单，而是：
   - `u_d` closure
   - atmospheric ODE
   - visible bubble fraction 闭环
   三者组合出来的错误 basin。
3. Heinbockel 分支当前又出现了另一个需要继续核对的量级问题：
   - `u_b ≈ 19 m/s` 在 HTW 2.5 MPa 下非常激进
   - 这虽然把 `eps_b` 拉回了低值，但可能又把 rise velocity 推得过高

因此下一步最该继续核的，不再是“有没有 pressure correction”，而是：

- Eq.3.43 的括号和系数在实现上是否还需进一步核实
- `u_b` 在高压下如此之大的量级是否与 Fig.3.3–3.5 的实验窗口一致
- `Eq.3.41 + Eq.3.43` 接入后，`K_bd` 的量级是否仍落在 Hamel Fig.3.6 可接受区间
5. 以及 **主路径仍未从 MW 切到 Hilligardt ODE**

也就是说，后面若要进一步对齐 Hamel，最终不可避免要面对：

- 是否把 bed 主 `d_b(h)` 从 `Mori-Wen` 切到 ODE

#### 验证

- `python3 -m py_compile src/physics/bubble_dynamics.py src/physics/__init__.py tests/test_bubble_dynamics.py`
- `python3 -m pytest tests/test_bubble_dynamics.py -q` → `8 passed`
- `python3 tests/sanity_checks.py` → `ALL 11 SANITY CHECKS PASSED`

### 2026-04-10: Heinbockel OCR correction pass on Eq.3.41 / Eq.3.42 / Eq.3.43

根据用户提供的德文论文第25页与第31页高清 OCR，本轮把 Heinbockel 高压链里的三处 source-of-truth 重新校正：

- `Eq.3.42` 继续锁定为：
  - `lambda_b = 280 * u_mf / g * (P/P0)^(-0.7)`
  - 也就是此前代码里的 `u_mf/g` 口径是对的，不切到任何 `epsilon_b^(1/3)/d_b^(1/3)` 变形。
- `Eq.3.41` 的分母不应再误用 `n_b = 2.7`，而应使用：
  - `n_B = 5`
  - 含义是同时聚并的气泡数量（Briens et al. 1988），不同于 `Eq.3.44` 里 `u_br = n_b * u_d * (P0/P)^0.15` 的 `n_b = 2.7`。
- `Eq.3.43` 当前新增了一个按 Hamel 第31页 OCR 直接实现的 audit branch：
  - `u_b = psi_b * (u0-u_mf) * ((P/P0)^0.2 + 1) + u_b,i`
  - 同时保留旧的 `2.14*(P/P0)^0.7` 版本为 `legacy` 对照，不再作为唯一 Heinbockel 实现。

新增审计：

- `scripts/audit_heinbockel_initialized_lu.py`

当前在 `u_d = wein_1992_eq312` 的 initialized LU state 上，对比结果如下：

- `hilligardt_atm_chain`
  - `u_b ≈ 1.32 -> 1.56 m/s`
  - `eps_b = 0.70`（全床仍撞 clip）
- `heinbockel_pressurized_chain_hamel_ocr`
  - `u_b ≈ 1.78 -> 2.75 m/s`
  - `d_b ≈ 0.18 -> 0.76 m`
  - `eps_b ≈ 0.41 -> 0.63`
  - `K_bd ≈ 1.09 -> 0.25 1/s`
- `heinbockel_pressurized_chain_legacy_2p14_p07`
  - `u_b ≈ 18.9 -> 19.3 m/s`
  - `eps_b ≈ 0.061`

当前新判断：

1. 旧 `2.14*(P/P0)^0.7` 分支给出的 `u_b ~ 19 m/s` 已可视为过激 legacy branch。
2. 按用户二次高清 OCR 校正为 `[(P/P0)^0.2 + 1]` 后，`Eq.3.43` 已能在 `P=P0` 时退化回常压 `2.0` 系数口径，不再存在之前的退化性张力。
3. 因此当前 `Eq.3.43` 的主未决项已从“括号符号是否正确”收缩为更细的符号层问题，例如原文中的 `v_{b,s}` 是否应显式保留在 `u_{b,i}` 前。

#### 验证

- `python3 -m py_compile src/physics/bubble_dynamics.py scripts/audit_heinbockel_initialized_lu.py tests/test_bubble_dynamics.py`
- `python3 scripts/audit_heinbockel_initialized_lu.py`

### 2026-04-10: Eq.3.15 source-of-truth locked

根据用户对德文原著第23页的像素级校对，`Eq.3.15` 现已正式锁定为：

- `psi_b = 0.17 * u_mf^(-0.33)`
- 适用区间：`0.014 m/s < u_mf < 0.18 m/s`

因此此前为了 OCR 分歧临时保留的 `wein_1992_eq315_055` audit branch 已从代码中移除，不再保留 `-0.55` 作为候选实现。

当前 `psi_b` 相关实现状态：

- `technical_distributor = 0.76`
- `porous_plate = 0.67`
- `wein_1992 = 0.17 * u_mf^(-0.33)`（已锁定）

新增审计：

- `scripts/audit_psi_b_strategy_initialized_lu.py`

它现在仅比较：

- `technical_distributor`
- `wein_1992`

在当前 Heinbockel initialized LU chain 下，top-bed 对比为：

- `technical_distributor`
  - `u_b_top = 4.1569 m/s`
  - `eps_b_top = 0.2726`
  - `K_bd_top = 0.4292 1/s`
- `wein_1992`
  - `u_b_top = 3.1519 m/s`
  - `eps_b_top = 0.3565`
  - `K_bd_top = 0.4191 1/s`

这说明：

1. `Eq.3.15` 的指数歧义已经结束，后续不再需要围绕 `-0.33 / -0.55` 打转。
2. `technical_distributor=0.76` 与 `wein_1992` 在 initialized Heinbockel chain 上会明显改变 `u_b / eps_b`，但对 `K_bd` 影响较小。

### 2026-04-10: shared global-NR builder switched to Hamel/Wein hydrodynamics chain

本轮将 shared `global_nr` 开发口径正式切换到论文一致性的 hydrodynamics 主链：

- `u_d = Eq.3.12 (Wein 1992)`
- `psi_b = Eq.3.15 = 0.17 * u_mf^(-0.33)`
- `d_b(h)` = `Hilligardt ODE`
- 加压修正 = `Heinbockel Eq.3.41 / 3.42 / 3.43 / 3.44`
- `xi_b` = `Eq.3.34`

具体已接入：

- `tests/validation_case_utils.py`
  - `build_phase1_htw_lu_global_nr_reactor_config()`
  - `build_phase2_htw_lu_freeboard_reactor_config()`

对应回归：

- `tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_config_uses_hamel_wein_hydrodynamics_chain`

#### chemistry-facing impact after hydrodynamics switch

1. top-bed chemistry audit  
`scripts/audit_top_bed_chemistry_lu.py`

结果：

- `Texit = 1185.1 K`
- top-bed (`cell 6-9`) sums:
  - `R2_area = 482.233`
  - `R4_area = 0.111`
  - `R7_raw = -24.930`
  - `R8_raw = -1.045`
  - `CH4_net = +1.025`
  - `H2O_net = -16.284`

结论：

- 顶床层 `H2O deficit` 仍由 `R2` 主控；
- `R7` 仍持续逆向，但绝对值较旧 hydrodynamics 口径已明显减弱；
- `R8` 仍只是弱逆向；
- `CH4` 仍在 bed 内过低，但其“完全被顶床层 chemistry 决定”的程度下降了。

2. freeboard-aware full profile audit  
`scripts/audit_global_nr_profile_lu_freeboard.py`

结果：

- `n_iter = 13`
- `rms = 3.024e-02`
- `nr_total_s = 88.58 s`
- `bed_exit_T = 1189.5 K`
- `reactor_exit_T = 1212.1 K`
- dry exit gas:
  - `CO = 0.1381`
  - `CO2 = 0.1101`
  - `H2 = 0.1474`
  - `CH4 = 0.0012`

相对文献出口偏差：

- `CO`: `6.25%`
- `CO2`: `0.06%`
- `H2`: `22.85%`
- `CH4`: `95.86%`
- `Texit`: `10.19%`

结论：

1. hydrodynamics 主链切换后，`CO/CO2` 已经非常接近参考值；
2. 主要剩余偏差收缩到：
   - `H2` 偏高
   - `CH4` 近零
   - `H2O` 轴向 profile 明显偏低
3. 这说明后续 chemistry 审计应全部建立在新的 Hamel/Wein hydrodynamics 上，旧 hydrodynamics 下关于 `R2/R7/R8/CH4/H2O` 的定量结论不再宜直接沿用。

#### 验证

- `python3 -m pytest tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_config_uses_separate_dense_fraction tests/test_table2_LU_global_nr.py::test_phase1_lu_global_nr_shared_config_uses_hamel_wein_hydrodynamics_chain tests/test_global_nr_solver.py::test_reactor_defaults_psi_b_to_wein_eq315 -q` → `3 passed`
- `python3 scripts/audit_top_bed_chemistry_lu.py`
- `python3 scripts/audit_global_nr_profile_lu_freeboard.py`

### 2026-04-10: bed solid-transition audit after hydrodynamics switch

新增审计：

- `scripts/audit_bed_solid_transition_lu.py`

目的：

- 直接检查当前 shared NR 新主线下，床内固相占比是否从 bed 向 freeboard 呈合理过渡
- 明确区分：
  - `dense_share = 1 - eps_b`
  - `dense_solid = 1 - eps_d_voidage`
  - `bulk_solid = (1 - eps_b) * (1 - eps_d_voidage)`

结果：

- bottom `bulk_solid = 0.2669`
- mid `bulk_solid = 0.3006`
- top `bulk_solid = 0.3127`

也就是 **bulk solid holdup 不是向床顶下降，而是反而上升**。

当前判断：

1. 这说明 hydrodynamics 虽然已经比旧主链更接近 Hamel source-of-truth，但 **bed-to-freeboard 的固相过渡仍未完全物理化**。
2. 当前更像是：
   - `eps_b` 自底向顶下降
   - `eps_d_voidage` 几乎保持常数 `~0.494`
   - 从而导致 `(1-eps_b)*(1-eps_d_voidage)` 反而上升
3. 因此 hydrodynamics 现在可以说是：
   - **major source-of-truth constants/closures 已基本锁定**
   - 但 **bed axial solids-transition behavior 仍未收口**

这也意味着：

- chemistry 审计现在应以新 hydrodynamics 主线为基础继续做；
- 但凡涉及 “床顶边界条件是否已物理可信” 的结论，都仍需保留这条 bed-solid-transition 风险说明。

#### 验证

- `python3 scripts/audit_bed_exit_to_reactor_exit_budget_lu.py`
- `python3 scripts/audit_freeboard_direction_lu.py`
- `python3 scripts/audit_bed_solid_transition_lu.py`

### 2026-04-10: 主线默认 `psi_b` 切换到 Eq.3.15

根据用户最终确认，`Eq.3.15` 已作为 `psi_b` 的主 source-of-truth：

- `psi_b = 0.17 * u_mf^(-0.33)`

因此本轮已将项目主线默认值从：

- `technical_distributor = 0.76`

切换为：

- `wein_1992`

涉及默认值变更的位置：

- `src/core/cell.py`
- `src/core/cell_hydrodynamics.py`
- `src/core/reactor.py`

保留项：

- `technical_distributor = 0.76`
- `porous_plate = 0.67`

它们仍作为可选配置保留，但不再是默认主路径。

### 2026-04-10 14:15 CST: interim handoff before pause

本轮新增与修正：

- 新增 `src/physics/phase_fractions.py::calc_bulk_solid_holdup()`
- 新增 `scripts/audit_epsb_profile_lu.py`
- 更新 `scripts/audit_bed_solid_transition_lu.py`，显式输出 `eps_avg`
- 补充 `tests/test_phase1.py` 相分率层级回归
- 修正 `tests/test_cell_balances.py` 中已过时的 `u_d_closure` 断言

#### 当前已确认的结论

1. `bulk solid` 公式本身没有符号错误：

   - `bulk_solid = (1 - eps_b) * (1 - eps_d_voidage)`
   - 等价于 `1 - [eps_b + (1 - eps_b) * eps_d_voidage]`

2. 当前真正的风险不是 `bulk_solid` 公式，而是 **`eps_b(h)` 的轴向趋势错误**：

   - 在当前 shared NR 主线下，`eps_b` 从 `0.4729 -> 0.3815`
   - `eps_avg` 从 `0.7331 -> 0.6873`
   - `bulk_solid` 从 `0.2669 -> 0.3127`

   这代表模型当前预测“床层越往上越密”，与 Hamel 对鼓泡床 bed 内部“空隙率向上增加”的物理图景相反。

3. `Eq.3.24` 的 profile 审计已经把直接原因拆清：

   - 分子 `u0 - u_d` 沿床层几乎不变（约 `1.34 ~ 1.44`）
   - 分母 `u_b + (n_b - 1) * u_d` 沿床层持续增大（`2.84 -> 3.57`）
   - 因此 `eps_b` 被结构性压低

4. 因而，当前 bed 内 hydrodynamics 的主要未决点已经收缩为：

   - `Eq.3.12 (u_d)`、`Eq.3.41/3.43 (d_b/u_b)`、`Eq.3.24 (eps_b)` 这组组合在 HTW 工况下的相互作用
   - 特别是：为什么 `d_b(h)` / `u_b(h)` 的上升，没有带来 Hamel 预期的 `eps_avg(h)` 上升

5. chemistry 方面当前只能给出保守判断：

   - freeboard 仍是小修饰项
   - 主要组成偏差仍来自 bed
   - 但在 bed 内 solids-transition 尚未物理化之前，不宜对 bed-top chemistry 做过强定量结论

#### 当前最可靠的验证快照

- `python3 -m pytest tests/test_phase1.py tests/test_cell_balances.py -q` → `47 passed`
- `python3 scripts/audit_bed_solid_transition_lu.py`
- `python3 scripts/audit_epsb_profile_lu.py`

#### 休息后恢复工作的优先顺序

1. **继续 bed hydrodynamics，不先回 chemistry**

   目标：
   - 让 bed 内 `eps_avg(h)` 恢复为随高度上升
   - 让 `bulk_solid(h)` 恢复为随高度下降
   - 到 bed surface 再与 freeboard 发生模型切换

2. **直接围绕 `Eq.3.24` 做 closure 审计**

   重点检查：
   - `u0` 是否应在 bed 内按局部状态更新，而不是近似常值
   - `u_d` 在 HTW 主算例里与 `Eq.3.24` 的组合是否仍有解释偏差
   - `u_b` 的高压修正是否仍偏强，导致分母沿程增长过快

3. **只有在 solids-transition 修正后，才重开 chemistry audit**

   届时再重新量化：
   - `R2` 对 `H2O deficit` 的主导强度
   - `R7` 逆向与 `CH4` 近零之间的关系
   - `H2` 偏高是否仍然存在
