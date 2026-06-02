# Status Report（Session Handoff）

**更新时间**：2026-06-01T00:00:00+08:00  
**仓库**：`dachou5224/bfb-gasifier`  
**目标主线**：在不放弃 baseline Newton + line-search 的前提下，提升 phase2（bed+freeboard）收敛质量与效率。

## 1) 本次会话关键结论

1. **bed-only 入口局部加密有效**：残差显著下降，说明入口细化对局部刚性/热点有帮助。  
2. **替代步进（LM/PTC/equil_newton）在当前问题上整体不优于 baseline**：多数情况下残差更差或 walltime 更长。  
3. **non-monotone line-search 仅在“窄窗+轻微放松”下有局部收益**，但经常带来额外 walltime。  
4. **最有效路线是 baseline 内优化**：line-search trial cap + lambda_seed 自适应 + only-on-failure non-monotone，在 bed-only 上实现了近似不增耗时的残差改进。  
5. **phase2 纳入 freeboard 后收敛平台仍明显**，仅靠 line-search 微调未突破。  
6. **根因指向 Jacobian 结构组装剪枝漏耦合**：structured/pruned 相对 dense 存在显著缺失耦合；关闭 pruning 后缺失比例明显下降。  
7. **已落地修复方向**：`band_plus_side_elements_structured` 统一禁用 local-BC/affected-row pruning，改为 full-BC 组装（包含 explicit freeboard 场景）。

## 2) 已完成的重要代码改动

### A. 网格与配置能力
- `src/core/reactor.py`
  - 新增 `bed_dh_profile`（床层非均匀网格）。
  - 扩展 thesis 求解配置：`nr_step_model_thesis`（`newton|lm|ptc|equil_newton`）及相关参数。
  - 增加 `nr_line_search_max_trials_thesis` 等 line-search 调优参数。

- `tests/validation_case_utils.py`
  - 新增 `build_phase2_htw_lu_local_refined_damped_config(...)`，用于入口局部加密+阻尼联动实验。

### B. NR 求解器策略
- `src/solvers/global_nr_solver.py`
  - 新增/接入 `LM`、`PTC`、`equil_newton` 步进模型。
  - 支持 non-monotone 接受准则（可配置）。
  - 增加 line-search 试探预算上限（trial cap）与统计（`line_search_cap_hits`）。
  - 优化 `lambda_seed` 自适应。
  - 将 non-monotone 限制为失败后触发。
  - **最新修复**：`band_plus_side_elements_structured` 下总是禁用 local callbacks（full-BC），避免 recycle/非局部边界耦合被 pruning 漏掉。

- `src/workflow/steps/nr_inner_step.py`
  - 将上述新配置透传到 `solve_global_nr`；扩展模型白名单和参数校验。

### C. 回归测试
- `tests/test_reactor_solid_streams.py`
  - 新增 `bed_dh_profile` 正常/异常测试。
- `tests/test_global_nr_solver.py`
  - 覆盖 step model 参数透传与未知模型拒绝。
  - **最新同步**：将 explicit freeboard structured 测试改为验证“忽略 local callbacks 且矩阵一致”。

## 3) 关键实验观察（可供后续 agent 直接复用）

- **bed-only（baseline 内微调）**：可在近似不增 walltime 下显著降低 RMS 与 component max。  
- **phase2（含 freeboard）**：tuned 组合主要表现为提速，收敛质量未明显优于 baseline，残差平台仍在。  
- **dense vs structured 审计**：pruned structured 的缺失耦合在 freeboard 场景下明显，关闭 pruning 可显著降低 missing fraction；缺失集中在 bed 长程耦合与 recycle 链路，不是单一 freeboard 边界项问题。

## 4) 当前状态（接手重点）

1. Jacobian 结构修复已实现并有针对性单测。  
2. 需要基于该修复**重跑 phase2 A/B**，确认是否突破 `~3.6e-2` 残差平台，并评估 `line_search_failures / line_search_evaluations / walltime`。  
3. 若质量改善成立，再做小步调参（优先 `line_search_max_trials` 与 failure-only non-monotone 窄窗），目标是“质量不退化 + 耗时下降”。

## 5) 已确认的门控结果（本轮）

- `pytest -q tests/test_global_nr_solver.py -k 'band_plus_side_elements_structured or block_tridiag_uses_local_callback'` 通过。  
- `python3 tests/sanity_checks.py` 通过（`ALL 17 SANITY CHECKS PASSED`）。

## 6) 注意事项

- 当前工作区有大量历史改动（含未跟踪文件），接手时请**避免误回滚**无关文件。  
- 若继续做 Jacobian 结构相关实验，优先聚焦：
  - `src/solvers/global_nr_solver.py`
  - `src/solvers/structured_jacobian.py`
  - `tests/test_global_nr_solver.py`
  - `tests/validation_case_utils.py`

## 7) 2026-05-16 续查结论（freeboard 收敛 root cause）

### A. 超时纪律

后续所有完整 phase2 求解实验必须使用硬超时包装，例如 `subprocess.run(..., timeout=420)`。
若 420s 内无结果，直接记录 `TIMEOUT`，不要继续等待。长实验前先确认无后台 Python：
`ps -o pid,ppid,etime,command`。

### B. 已排除/降权的假设

1. **当前真实 inner NR 下 structured Jacobian 仍漏耦合**：降权。
   - 在真实 freeze 上下文中复测 dense vs `band_plus_side_elements_structured`，`diff nnz=0`。
   - 非 freeze 下确有长程 bed coupling 差异，但那不是 inner NR 的实际线性化状态。

2. **温度 Newton 步长上限过大**：排除为主因。
   - 将 T step limit 临时压到 300/100/50 K 均使 phase2 变差。
   - baseline 600 K 反而能到当前较低平台。

3. **freeboard reaction stiffness 主导**：排除为主因。
   - 关闭 freeboard reactions 后 RMS 约 `1.79e-2`，比 baseline `1.45e-2` 差。

4. **secondary air 局部注入突变主导**：排除为主因。
   - 将 secondary air 搬到底部或改 `distributed_uniform + local_refine=4` 均更差。

5. **recycle/cyclone capture 主导**：排除为主因。
   - 关闭 recycle/capture 后结果几乎不变。

6. **outer refresh 未同步 freeboard closure transport**：排除为直接修复。
   - `preserve_gas_state=True` 强制刷新 freeboard closure/transport 后 RMS 变差到约 `2.21e-2`。

7. **freeboard/side-block 单气相入口合并即可解决**：排除为单独充分修复，但保留为结构修正。
   - Hamel 口径下 freeboard 使用 ghost-bubble/颗粒轨迹解析闭合，side elements 主要承载固体分离/返料；二者都不应在入口继续强制保留床层 `N_d/N_b` 两相分配。
   - 已将进入 explicit freeboard / cyclone 的上游气体改为 `N_d_in = prev.N_d + prev.N_b`、`N_b_in = 0`。该修正保持总气相守恒，同时避免不存在 bubble phase 的退化段承担独立 bubble inlet residual。
   - 同一 3-bed/2-freeboard 初始态 A/B：`gas_combined_rms` 不变（`6.37e-2`），但 `gas_phase_split_rms` 从 `1.26e-1` 降到 `7.48e-2`，总 RMS 从 `6.12e-2` 降到 `5.02e-2`。
   - 标准 phase2 长跑（420s 硬超时包装，实际 397s 返回）：`rms_scaled_final ≈ 1.4729e-2`，`energy RMS ≈ 2.91e-2`，`gas RMS ≈ 1.89e-2`，说明 energy 改善但 gas/combined 残差仍未突破；入口合并是必要但不充分。

### C. 当前最可信 root cause

phase2 平台不是总守恒失败，而是**床顶/自由板接口附近的 dense-bubble 相分配残差与 freeboard 能量残差互相拉扯**。

关键证据：
- baseline phase2：`rms_scaled_final ≈ 1.4458e-2`，`energy RMS ≈ 3.36e-2`，`gas RMS ≈ 1.79e-2`。
- 最大逐项残差来自 `bed[9]` 顶部床层 gas split：
  - `bed[9] gas_comb_max ≈ 1.8e-3`（总气相守恒很好）
  - `bed[9] gas_split_max ≈ 9.18e-1`（dense/bubble 分配很差）
  - CO/H2O/N2 等 dense 与 bubble residual 成对反号，combined 几乎抵消。
- freeboard 多个 cell 的 gas/energy residual 仍显著，例如 `freeboard[16] gas_comb_max ≈ 9.25e-2`。
- 手动保持 `N_d+N_b` 不变、投影修正 `N_d/N_b` 会让非线性残差爆炸，说明不能做后处理投影；需要在 Newton 方程尺度/merit function/变量参数化层面处理。

### D. 下一步建议

1. 已增加一个低成本诊断输出：最终结果透传
   `rms_scaled_gas_combined_final` / `rms_scaled_gas_phase_split_final` 以及对应 max_abs，
   避免只看普通 gas RMS。
2. 优先考虑 Hamel-compatible 的变量/尺度改造：
   - bed cells 保留两相方程，但 line-search merit 同时报 `combined/split`，避免总守恒好而 split 大时误判；
   - freeboard cells 物理上应更接近 single-gas phase，下一步应只对 freeboard 采用 single-gas unknown/combined residual，或给 freeboard bubble residual 做结构性消元；单纯入口合并已验证不够。
   - cyclone/return-leg 等 side elements 也要继续审查：它们当前无反应/退化水动力，但仍可能带着 full gas/T residual rows 进入 NR。
3. 暂不建议继续投入 LM/PTC/non-monotone：
   - non-monotone `relax=1.005` 只能小幅改善到约 `1.39e-2`，耗时约 20min 且不稳定。
   - `relax=1.02` 更差。

### E. 2026-05-16 接口结构修正续查

1. 已在 residual projection 层加入退化单相段处理：
   - `cell_type in {"freeboard", "cyclone", "return_leg"}` 时，dense gas residual row 改为 `gas_d + gas_b` 的 total gas balance。
   - bubble gas residual row 改为 `N_b` anchor，驱动退化段 bubble holdup 回到 0。
   - 诊断中的 `gas_phase_split` 只统计真实 bed two-phase cells，避免 freeboard/cyclone 的 anchor row 污染 split 指标。
2. 3-bed/2-freeboard 小网格有限迭代验证：
   - `rms_scaled_final ≈ 3.004e-3`，`gas_combined_rms ≈ 2.736e-3`，`energy_rms ≈ 5.018e-3`。
   - 指标语义更干净，但最终平台基本不变；说明 row 消元是结构正确性修复，不是完整收敛突破。
3. 小网格 cell hotspot：
   - bed[0] 仍是最大 phase split 来源：`gas_split_max ≈ 2.86e-2`。
   - bed[1]/bed[2] split 明显较小：`≈5.0e-3 / 4.7e-3`。
   - freeboard[0]/[1] 已无 split residual，`gas_comb_max ≈3.9e-3 / 1.25e-3`。
   - cyclone/return_leg gas residual 基本为 0；side elements 不是当前小网格平台主因。
4. 下一步建议：
   - 优先审查 bed[0] 与外部 feed/recycle/phase inlet 的两相分配结构，而不是继续盯 freeboard。
   - 若继续验证标准 phase2，必须使用更短交互策略：先 residual-only / 小网格，再标准长跑；长跑必须硬超时且不阻塞交互。

### F. 2026-05-17 性能/freeze 续查

1. 后台状态：
   - 多次 `ps -o pid,ppid,etime,command` 均未发现残留 Python/pytest/NR 求解进程。
   - freeze 主要来自 agent 交互层等待长命令返回，而不是后台进程泄漏。
2. 标准 phase2 `max_global_iter=5` 探针：
   - 35s 硬超时仍未完成第一次内层 Jacobian。
   - 超时栈位于 `assemble_jacobian_fd_structured -> cell.residuals() -> reactions/solid_balance`。
3. 已落地性能结构修正：
   - explicit freeboard 的 side-pair 结构从 head×tail 保守网收窄为实际影响拓扑。
   - 标准 phase2：`side_pairs 215 -> 40`。
   - 理论 cell residual calls / Jacobian：`7475 -> 2710`。
   - 退化段 bubble anchor columns 用解析 Jacobian，进一步估算到 `2283` 次 cell residual calls / Jacobian。
4. 仍未达标：
   - 35s 探针仍超时，说明单个 cell residual 成本较高；只收窄图还不够。
   - freeboard cell 数不是主要因素：nfb=1..8 估算调用数约 `2099 -> 2283`，变化不大。
   - bed cell 数是主要线性成本源：nb=3/5/7/10 时估算调用数约 `759/1167/1575/2187`。
5. 下一步性能方向：
   - 优先降低 Jacobian rebuild 次数（例如 phase2 默认 `nr_jacobian_lag_steps` 做受控 A/B），但第一轮 Jacobian 仍需优化。
   - 更深层需要优化 bed residual pipeline：缓存温度相关 Arrhenius/thermo、避免每个 FD scalar perturbation 全量重算与该变量无关的反应块。
   - 在 agent 操作层面，标准全模型求解应以 `timeout<=100s` 为硬上限；探针命令 `yield_time_ms<=10s`，不返回则轮询一次并停止等待。

## 8) 2026-05-29 续记：Hamel 结构桥接工程记录

已新增 `docs/hamel_engineering_bridge_notes.md`，用于记录这轮解决过的 band+side-element、bed/freeboard 缝合规则，特别是 Hamel 原文没有逐项给出实现细节的部分。

核心记录点：

1. **cell-type-specific NR layout**：
   - bed 保留 dense/bubble 两相 gas + solid holdup + T。
   - freeboard 使用 total single gas + T；solid 由 trajectory/closure 管理，不进入 NR unknown。
   - cyclone/return_leg 作为 solid side elements，不再伪装成 bubbling-bed cell。

2. **bed -> freeboard 显式投影**：
   - `N_g,fb,in = N_d,bed,out + N_b,bed,out`。
   - bubble split residual 只存在于 bed 内，不跨入 freeboard。
   - freeboard 对 dense/bubble 纯重分配不敏感，只响应 total gas 改变。

3. **side-element Jacobian 缝合**：
   - `bed -> freeboard -> cyclone -> return_leg -> bed` 作为显式拓扑进入 residual/Jacobian。
   - `band_plus_side_elements_structured` 禁用 local pruning，避免 recycle/freeboard/side coupling 被剪掉。

4. **closure refresh 不突变原则**：
   - freeboard closure refresh 可作为诊断，但默认 solver path 不强刷。
   - result finalizer 只报告 gap，不回写突变 NR 已接受状态。

5. **当前收敛状态**：
   - Phase2 LU 默认 `nr_jacobian_lag_steps=4`。
   - 10-outer probe：`rms_scaled_final≈0.00182`，`max_abs_scaled_final≈0.01383`。
   - RMS 已基本收敛，剩余严格失败点是少数 max residual rows。

## 9) 2026-05-30 续记：bottom inlet split 归入 Vorabrechnung

Hamel 的求解流程没有直接给出 `N_zu,d / N_zu,b` 固定比例；它只把
`Startwertwahl` 放在 configuration，同时明确 `Vorabrechnung` 负责气固流、
孔隙率、相间交换与停留时间。基于这个证据，Phase1/Phase2 global-NR 主线
已从 legacy fixed split 转为 `gas_inlet_split_strategy="precalc_hydrodynamic_flux"`。

当前实现：

1. `fixed` 策略仍保留，`gas_inlet_dense_frac` 作为 legacy/fallback。
2. Hamel-aligned 策略使用预计算水力学通量：
   `dense_frac = u_d * (1 - eps_b) / u0`。
3. init/precalc 和 outer refresh 在水力学刷新后会重新投影 bottom gas feed，
   避免把旧 fixed split 当作最终边界条件。
4. 新增单测覆盖 fixed、hydrodynamic、precalc 缺失时 fallback，以及 feed rows 投影。

### 9.1 短 probe 收敛检查

同一当前代码、Phase2 LU、`max_global_iter=6` 对照：

| split strategy | resolved dense frac | rms_scaled_final | max_abs_scaled_final | component max | wall |
|---|---:|---:|---:|---:|---:|
| `fixed` | 0.2000 | 0.002453 | 0.018946 | 0.006102 | 43 s |
| `precalc_hydrodynamic_flux` | 0.2227 | 0.004819 | 0.044168 | 0.007291 | 33 s |

结论：目前没有收敛改善；预计算通量 split 作为 Hamel 结构修正是合理的，
但直接替换 legacy `0.20` 会让短 probe 的 bed0/bed1 gas residual 变差。

新策略 hotspot：

- `bed`：RMS≈0.00701，max≈0.0515。
- `freeboard`：RMS≈0.000606，max≈0.00119。
- `cyclone`：RMS≈0.000993，max≈0.00132。
- `return_leg`：近零。

最差 rows 集中在 `bed[0]` / `bed[1]` 的 dense gas：`TAR2`、`CO`、`CO2`、
`N2/H2` split、`CH4/NH3/H2O`。这表明下一步不是 freeboard/side 缝合，
而是 bottom bed 两相入口、Kbd/exchange 与快速反应源之间的局部刚性。

### 9.2 bottom bed 刚性续查

已修正一个实现语义问题：`precalc_hydrodynamic_flux` 现在在
Vorabrechnung 水力学刷新后 snapshot 到 `bed[0]`，inner NR 中
`set_bottom_cell_feeds()` 读取冻结值。此前它会从 live `u0/u_d/eps_b`
重新解析，容易把 bottom inlet split 变成滞后的移动边界。

修正后短 probe：

- frozen/resolved dense frac：0.222725885。
- live `u_d(1-eps_b)/u0` 与 frozen 值一致，不再漂移。
- 但 `max_global_iter=6` 的 residual 仍为
  `rms_scaled_final≈0.004819`、`max_abs_scaled_final≈0.04417`。

项级诊断（fixed 0.20 vs precalc 0.2227）说明：

- 小幅入口 split 变化会触发 bed0 dense reaction terms 的数 mol/s 级重排；
  CO、CO2、CH4、H2O 的 dense source/sink 是主放大器。
- N2 total residual 近零，但 phase split residual 变差，说明两相入口、
  `N_ex/Kbd` 和 outlet phase composition 没有形成平滑闭合。
- 继续方向应检查 bottom cell 两相气体 balance 中，bubble bulk inlet、
  throughflow/exchange、以及快速 homogeneous/pyrolysis gas source 的职责是否重复或错位。

### 9.3 damping 与 reaction continuation 诊断

当前 global NR 已有 Hamel-style damped Newton：

- thesis path 默认 full Newton step first，再二分降低 `lambda`。
- 失败步最多试探 `nr_line_search_max_trials_thesis` 次；Phase2 当前 cap 为 12。
- `dx` 先经过 gas/solid/T 物理 clipping；thesis path 不启用 GD fallback。
- 可配置 `lm` / `ptc` / `equil_newton`，但短 probe 显示它们不是当前 bed0 刚性的直接解。

Phase2 LU、`max_global_iter=6` 的 line-search 证据：

### 9.4 bed0 primary gas exchange 预投影（2026-06-01）

新增一个只作用于初始化阶段的窄口径预投影：
`preproject_bottom_primary_gas_state_to_exchange_closure()`。

设计边界：

- 只作用于 `bed[0]`，且只调整 `O2/H2O/N2`。
- 严格保持每个 primary species 的总量 `N_d + N_b` 不变。
- 不移动 `CO/CO2/H2/CH4` 等反应产物，避免把 dense 快速反应产物硬塞入 bubble。
- 使用当前 Vorabrechnung 水力学下的 cell residual，解 3 变量 bounded least-squares；
  只有 primary dense residual 范数至少下降到原来的 `0.75` 以下才接受，否则回滚。
- 只在 `run_init_and_precalc_for_global_nr()` 的 x0/BC 初始化后执行一次；
  outer refresh 仍只做 bottom inlet split/state split 对齐，不做 exchange-aware state rewrite。

静态审计结果（不进入 global NR）：

- 修复前典型 bed0：N2 dense/bubble residual 为 `+10.58/-10.58 mol/s`，总和约 0，
  属于纯相分配刚性。
- 默认初始化后：N2 dense residual `≈ -7e-15`，bubble residual `≈ 3e-10`。
- bed0 gas max residual 从约 `13.6 mol/s` 降到约 `9.52 mol/s`。
- 剩余主导项转为 `CO` dense residual `≈ 9.52 mol/s` 与 `CO2` exchange/反应项，
  因此下一步焦点应从 N2 phase split 转向 bottom fast reaction / CO outlet 初值闭合。

单测门禁：

- `tests/test_bottom_gas_inlet_split.py::test_bottom_primary_exchange_preprojection_reduces_bed0_primary_phase_stiffness`
  覆盖 primary 总量守恒、N2 相刚性下降、产品组分不被预投影改写。

### 9.5 bed0 major product gas x0 局部投影（2026-06-01）

继续追查发现：primary exchange 预投影后，bed0 最大 gas residual 转到
`CO/H2/CH4` dense rows。项级拆解显示：

- `CO` dense residual 约 `9.52 mol/s`，主要来自 `R_gas_d[CO]`。
- 其中约 `5.67 mol/s` 直接来自 Vorabrechnung VM/pyrolysis cache。
- `CH4` dense source 几乎全部来自 VM/pyrolysis。
- `generate_initial_x0()` 的 Hamel major Gibbs seed 在 bed0 强氧化入口下给
  `CO/H2/CH4 ≈ 1e-12`，因此与后续固定 kinetic/pyrolysis source cache 脱节。
- 关闭 major Gibbs x0、回退 legacy seed 会让局部 two-phase residual 更差，
  因为 legacy 组成与当前 `K_bd/N_ex`、O2/H2O 限制器不一致。

新增第二个只用于初始化的局部投影：
`preproject_bottom_major_gas_state_to_local_balance()`。

设计边界：

- 仍只作用于 `bed[0]`。
- `O2/H2O/N2` 只允许相分配变化，总量保持不变。
- `CO/H2/CH4/CO2` 允许 dense/bubble 出口初值获得局部 source support。
- `TAR1/TAR2` 暂不纳入；静态试验显示纳入 tar 会放大 O2/H2O split residual。
- 接受准则：全气相 max residual 必须至少降到原来的 `0.75` 以下，否则回滚。

静态审计结果（不进入 global NR）：

- primary 预投影后：bed0 gas max residual `≈ 9.52 mol/s`，主导为 `CO` dense。
- major gas 局部投影后：bed0 gas max residual `≈ 3.19 mol/s`。
- `CO` total residual 从约 `9.52 mol/s` 降至约 `2.54 mol/s`。
- 剩余主导项为 `CH4` dense source `≈ 3.19 mol/s`，说明下一步应检查
  VM/pyrolysis CH4 source 与 bottom oxidation/steam reforming 初值之间的闭合。

单测门禁：

- `tests/test_bottom_gas_inlet_split.py::test_bottom_major_gas_preprojection_reduces_bed0_product_residual_without_moving_primary_totals`
  覆盖 major gas 投影降低 bed0 gas max、primary 总量守恒、CO/H2 获得非零出口支撑。

### 9.6 CH4 oxidation/reforming supply limiter（2026-06-01）

继续追 CH4 后发现：默认 major gas 投影后，`CH4` 仍接近零，但
`R_gas_d[CH4]≈3.19 mol/s`，几乎全部来自 VM/pyrolysis source。手动扫 CH4
出口初值时，旧实现中只要给 `CH4 total≈0.001 mol/s`，R6/R7 就会在 O2 充足时
允许数十 mol/s 级 CH4 消耗，导致 residual 从 `+3.19` 瞬间跳到大负值。

修复：

- 在 `build_reaction_sources()` 中为 CH4 消耗增加 supply limiter。
- bubble 相限制 `r6b_lim` 不超过 bubble 可用 CH4。
- dense 相限制 `r6d_lim + ext7_lim` 不超过 dense 可用 CH4。
- dense 可用量包括 `N_zu_d/N_d_in/N_rez_d`、bubble→dense exchange、VM/pyrolysis
  CH4 source，以及 R3 hydrogasification 的 CH4 生成。
- 该 limiter 仿照已有 CO bubble supply limiter，属于守恒型数值护栏，不改动力学常数。

静态效果：

- `CH4 total=0.001 mol/s` 时，CH4 residual 不再从正源项跳到 `~-8.75 mol/s`；
  dense/bubble CH4 residual 保持在小量级。
- 默认初始化后的 bed0 gas max 仍约 `3.19 mol/s`，因为默认 CH4 出口初值仍在
  `~1e-12`；下一步不是再加 limiter，而是研究 VM/pyrolysis CH4 source 是否应在
  x0 中获得一个极小但非零的 residence/outlet support，且不能放大 O2/H2O split。

单测门禁：

- `tests/test_bottom_gas_inlet_split.py::test_bed0_ch4_oxidation_is_limited_by_available_ch4_supply`
  覆盖微量 CH4 初值下，R6/R7 消耗不超过本相可用 CH4 supply。

### 9.7 Hamel 证据下的 pyrolysis 两相分配原则（2026-06-01）

Hamel 原文/本地证据没有给出“bed0 初始 pyrolysis 产物如何在 dense/bubble
两相分配”的显式规则。可确认的是：

- `Vorabrechnung` 负责按高度计算水分与热解产物释放，并映射到 cell 模型。
- 气相守恒方程分为 `Suspensionsphase` `(2-4)` 与 `Blasenphase` `(2-5)`。
- 气泡相不含固体；两相之间通过同一个 `ex_bd` 相间交换项耦合。
- 4.3 节只说明热解气组成由 Gibbs/free enthalpy minimization 得到，并未说明
  热解气直接进入 bubble。

因此当前项目原则固定为：

- 干燥/热解气源项先进 `dense/suspension`：`gas_src_vm -> R_gas_d`。
- `R_gas_b` 不接收直接 pyrolysis source。
- bubble 中的 pyrolysis-derived species 只能通过 `N_ex/K_bd` exchange 与 bubble
  outlet balance 获得。
- 初始化投影中允许 `CO/H2/CH4/CO2` 的 bubble outlet 变量移动，只表示为满足
  exchange/outlet balance 的初值支撑，不代表把 pyrolysis source 直接分配到 bubble。

新增单测：

- `tests/test_cell_kinetics.py::test_pyrolysis_gas_source_enters_dense_phase_not_bubble_phase`
  明确锁住 `gas_src_vm` 只进入 dense source、bubble source 为零。

短 probe 注意：

- 2026-06-01 尝试 `max_global_iter=2` 的 Phase2 全局求解，外层 `timeout=30s`
  触发 `TIMEOUT_30S`。当前不要把全局求解当作快速门禁；若要评估整体收敛改善，
  应单独做带硬超时的 profile，并拆分 init/Jacobian/line-search wall time。

| path | rms_scaled_final | max_abs_scaled_final | line search failures | backtracks | observation |
|---|---:|---:|---:|---:|---|
| `precalc_hydrodynamic_flux` | 0.004819 | 0.044168 | 9 | 105 | 后期 `lambda≈1e-3`，30%-50% 变量被 clipping |
| legacy fixed `0.20` | 0.002453 | 0.018946 | 8 | 124 | residual 更低，后期 clipping 仅 1-2 个温度变量 |
| `lm` step model | 0.009565 | 0.104889 | 11 | 142 | clipping 降低但 residual 更差 |
| `ptc` step model | 0.031809 | 0.397715 | 6 | 44 | 很快卡在高 residual |

结论：已有 damping，而且正在大量介入；但它只是缩小全局步长，不能修复
`bed0` 中 bottom inlet split、`K_bd/N_ex` 相间交换、dense 相快速反应源项
共同形成的局部刚性。

新增了默认关闭的 reaction-source continuation 诊断钩子：

- `Cell.nr_reaction_rate_multiplier` 可作为 global NR residual 的默认反应倍率；
  显式传入 `cell.residuals(rate_multiplier=...)` 的旧路径保持优先。
- `ReactorConfig.nr_reaction_continuation_stages_thesis`、
  `nr_reaction_continuation_height_m_thesis`、
  `nr_reaction_continuation_inner_iter_cap_thesis` 可让下部 bed cell 先低反应倍率 warm-start，
  最终阶段强制回到 `1.0` 完整物理方程。
- 该功能默认 `(1.0,)`、height `0.0`，不会改变 Phase2 默认验证路径。

短 probe：对下部 1 m 使用 `(0.25, 0.5, 1.0)`、warmup cap 2，
最终完整方程反而变差到 `rms_scaled_final≈0.01234`、`max_abs_scaled_final≈0.0965`。
因此不建议把 reaction continuation 作为默认收敛修正；它目前只作为诊断工具保留。

更有希望的方向是下部床层局部加密：`build_phase2_htw_lu_local_refined_damped_config`
短 probe 没有 line-search failure，accepted `lambda` 维持在 `0.25 -> 0.0625`，
说明缩小 bottom control volume 比单纯增加 damping 更接近问题本质。

### 9.4 非均匀床层 mesh 的 Vorabrechnung tau 修正

检查局部加密后发现一个实现假设：`vorabrechnung_tau_for_cell()` 使用
`cell.geo.dh * n_cells / u_mf` 估计固相停留时间，这只在 uniform bed mesh 下等于
`H_bed / u_mf`。对于 bottom-refined mesh：

- 细底层 cell 的 drying/pyrolysis tau 被低估。
- 粗上部 cell 的 tau 被高估。
- 这会把 Vorabrechnung source 分布扭曲，并放大 solid holdup residual。

已修正：

- `Reactor._build_cells()` 给 bed cell 写入 `_vorab_bed_height = H_bed`。
- `vorabrechnung_tau_for_cell()` 优先使用 `_vorab_bed_height / u_mf`，
  `n_cells * dh` 只作为旧 uniform fallback。
- 新增单测确认非均匀 bed mesh 下所有 bed cell 的 Vorabrechnung tau 都等于
  `H_bed / u_mf`。

修正后 refined short probe：

| config | max_global_iter | rms_scaled_final | max_abs_scaled_final | line search failures | wall |
|---|---:|---:|---:|---:|---:|
| refined bottom mesh | 4 | 0.012887 | 0.189787 | 0 | 19.5 s |
| refined bottom mesh | 6 | 0.012374 | 0.136207 | 19 | 86.8 s |

结论：tau bug 必须修，但局部加密本身仍不能直接作为默认收敛修正。
它改善了早期步长接受性，却把 solid holdup/transport residual 和 freeboard residual
放大。下一步应回到 bottom bed 两相 gas balance 与 solid holdup transport 的尺度/闭合，
尤其是 `K_auf/K_ab` 与 char/ash holdup residual 是否和非均匀 `dh` 一致。

### 9.8 底部 bed energy balance 初步审计

新增快速诊断脚本：

- `scripts/audit_bottom_bed_energy_balance_lu.py`

该脚本只做静态审计，不运行全局 NR 长求解。当前 LU Phase2 初始化后：

| cell | T [K] | energy residual [W] | current-stream T root | gas-closed T root | main observation |
|---|---:|---:|---:|---:|---|
| bed0 | 1325.0 | -6.106e6 | none in 250-2500 K | 531 K | 当前出口气相未带足 gas-balance closure 的焓支撑 |
| bed1 | 1304.4 | -4.359e6 | 747 K | 812 K | 受 bed0 出口与 VM/moisture 轴向传递影响 |

bed0 焓流分解：

- fresh gas: `-2.793e6 W`
- fresh solid/fuel: `-5.607e6 W`
- hot solid downflow from bed1: `+2.126e6 W`
- outlet gas: `-3.057e6 W`
- outlet solid: `+2.927e6 W`

解释：

- 这不是单纯的 temperature 初值问题；在当前出口物流下，bed0 能量方程没有物理温度根。
- 若先用当前 gas residual 关闭 total gas outlet，bed0 出现约 `531 K` 的能量根，说明能量残差与
  gas-balance 初值耦合很强。
- 当前主导问题更像是 bottom bed 的 `gas x0 / pyrolysis products / oxidation support / T x0`
  没有做联合预投影，而不是 freeboard/side-element 拼接错误。

下一步建议：

- 不再先改 side elements。
- 审查 `generate_initial_x0()` 中 1050-1325 K 的 bottom temperature seed 是否过热。
- 尝试一个 bounded bottom-bed joint preprojection：先使 bed0 total gas 接近 balance closure，
  再将 `T` 投影到能量方程可达根附近，但必须保持 pyrolysis gas 仍只进入 dense phase。

### 9.9 bottom-bed bounded joint preprojection 实施结果

已新增初始化阶段的 bed0 joint preprojection：

- 函数：`preproject_bottom_total_gas_and_temperature_to_energy_closure()`
- 接入点：`run_init_and_precalc_for_global_nr()`
- 变量：bed0 的 7 个主气体 dense outlet、7 个主气体 bubble outlet、`T`
- 目标：total gas closure、两相 gas residual、energy residual
- 约束：`T` bounded in `650-1600 K`
- dense-only 支撑：`NH3/TAR1/TAR2` 只允许加到 dense outlet；bubble outlet 不作为
  pyrolysis direct source。

局部效果：

| metric | before | after |
|---|---:|---:|
| bed0 T | 1325 K | ~1078 K |
| bed0 energy residual | ~-6.1 MW | ~-0.05 MW |
| bed0 current-stream T root | none | ~1071 K |
| bed0 gas-closed T root | ~531 K | ~1095 K |

新增测试：

- `tests/test_bottom_gas_inlet_split.py::test_bottom_joint_gas_temperature_preprojection_reduces_energy_and_total_gas_residual`

重要副作用：

- bed0 的问题被消除后，bed1 接收新的 bed0 outlet enthalpy，成为新的 energy hotspot：
  `bed1 energy residual ≈ -10.6 MW`。
- 全局静态初始 residual 因此没有改善：
  `scaled RMS ≈ 0.034`，`scaled max ≈ 0.565`，主导项转移到 bed1 energy。

结论：

- bed0 bounded joint preprojection 的局部数学目标成立。
- 但作为全局收敛修正还不够；它证明 bottom-bed 能量问题不能只按单 cell 修，
  需要把 bed0 -> bed1 的下游能量响应纳入同一个 bottom-zone bridge，或者顺序投影后
  同时检查相邻 cell residual，避免把 MW 级缺口上移。

### 9.10 Vorabrechnung temperature fence

新增“possible solution”温度围栏，防止 NR 或初始化投影掉入物理上不可信的低温数学根：

- 配置：
  - `nr_temperature_fence_enabled_thesis = True`
  - `nr_temperature_fence_lower_margin_K_thesis = 150`
  - `nr_temperature_fence_upper_margin_K_thesis = 350`
- 初始化阶段用 `estimate_axial_T_profile()` 的 `T_est[i]` 给每个 bed cell 写入：
  - `_vorab_temperature_estimate_K`
  - `_nr_temperature_min_K = max(300, T_est - lower_margin)`
  - `_nr_temperature_max_K = min(2500, T_est + upper_margin)`
- 应用路径：
  - `nr_indexing.unpack_cell()` 对所有 NR trial/unpack 后温度夹紧。
  - `global_nr_solver` line-search trial clipping 使用 per-cell fence，而不是固定
    `300-2500 K`。
  - bottom joint preprojection 的 temperature bounds 与该 fence 取交集。

LU Phase2 当前效果：

- bed0 `T_est = 1325 K`，因此下限为 `1175 K`。
- 之前 joint projection 会把 bed0 推到约 `1078 K`；现在被限制到约 `1175 K`。
- bed0 energy residual 保持在可接受的低 MW 以下量级：约 `-0.39 MW`，不再追逐
  `~1070 K` 的低温数学根。
- bed1 仍然是全局 hotspot，约 `-10.2 MW`，说明围栏解决“不可行低温支路”问题，
  但不解决 bottom-zone 轴向能量联动。

新增/更新测试：

- `tests/test_global_nr_solver.py::test_unpack_reactor_respects_vorabrechnung_temperature_fence`
- `tests/test_bottom_gas_inlet_split.py::test_bottom_joint_gas_temperature_preprojection_reduces_energy_and_total_gas_residual`
