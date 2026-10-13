# Hamel Parity 审计 — a1-gibbs — 2026-07-06

**证据链**：`docs/hamel_submodels/07_appendix_a1_gibbs_solver.md`（PDF pp. 91–94 / `data/hamel_pages_png/page_91.png`–`page_94.png`）→ `src/thermodynamics/gibbs_minimizer.py` → `src/thermodynamics/gibbs_hamel_reduced.py` → `src/solvers/vorabrechnung/vorab_gibbs_seed.py`

**调用入口**：Vorab x0 / major-species seed（`major_gibbs_solver_mode`: `augmented` | `hamel_reduced`）；Cell 微量组分（`minor_species.py` → `GibbsMinimizer`）。

**判定口径**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离`

---

## 1. 方程编号对照总表

| Hamel 编号 | 论文含义（摘要） | 实现模块 | 代码锚点 | 判定 | 备注 |
|------------|------------------|----------|----------|------|------|
| **(A.1)** | 目标函数 \(G=\sum_j \mu_j N_j\) | — | 无显式 \(G\) 求值 | **缺失（有意）** | 拉格朗日重写后直接解 KKT；不构造 \(G\) 标量，与 Hamel 算法链等价 |
| **(A.2)** | 元素守恒 \(\sum_j a_{ij}N_j = b_i\) | 调用方 | `vorab_gibbs_seed._major_elements_from_feeds`；solver `elements` 入参 | **部分一致** | 守恒在残差中强制；(A.2) docstring 已标注；\(N_Q\) 在 seed 层并入 \(b_i\) |
| **(A.3)** | \(N_Q\) 修正后的元素量 \(b_i\) | 调用方 | `vorab_gibbs_seed._major_elements_from_feeds`；pyrolysis `_hamel_project_pyrolysis_element_pool_for_a1` | **部分一致** | (A.3) 进料 O₂/H₂O/N₂ + VM 元素 → `b_i` 已在 docstring 逐项说明 |
| **(A.5)–(A.11)** | 理想气体化学势 + 浓度约束 + \(\lambda_i\) 中间式 | `GibbsMinimizer` | `_calc_moles`：`c=\mu^0/(RT)+\ln(P/P_0)`，`n_j=N\exp(a^T\lambda-c)` | **部分一致** | 合并为增广 \(NE{+}1\) 系统，未逐步标注 (A.5)–(A.11) |
| **(A.14)–(A.17)** | 继续拉格朗日展开 / 消元中间量 | `GibbsMinimizer` | 同上增广形式 | **部分一致** | 结构等价，编号未逐式标注 |
| **(A.18)** | 元素方程组（消元前） | `GibbsMinimizer` | `_residuals`：`f[:N_E] = a @ n - b` | **一致** | 增广系统第一块 |
| **(A.20)** | 总摩尔归一 \(\sum n_j/N = 1\) | `GibbsMinimizer` | `_residuals`：`f[N_E] = sum(n)/N - 1` | **一致** | 以 `ln N` 为第 \(NE{+}1\) 未知量 |
| **(A.21)** | 消去 \(z\)：\(z=A_{ref}/b_{ref}\) | `ReducedHamelGibbsSolver` | `_recover_moles_from_ratio`：`z=a_ref/b_ref`，`n=q/z` | **一致** | 参考元素消元 |
| **(A.23)** | 降阶指数项 \(q_j=\exp(a^T\lambda-c)\) | `ReducedHamelGibbsSolver` | `_calc_q` | **一致** | |
| **(A.25)** | 元素加权和 \(A_i=\sum_j a_{ij}q_j\) | `ReducedHamelGibbsSolver` | `_calc_aq` | **一致** | |
| **(A.26)** | 降阶方程 \(F(\lambda)=0\) | `ReducedHamelGibbsSolver` | `_reduced_residual`：`b_ref·A_i - b_i·A_ref=0` | **一致** | 参考行消元后的 \(NE{-}1\) 或 \(NE\) 维系统 |
| **(A.27)** | Newton 线性化 \(J\Delta\lambda=-F\) | 两者 | `GibbsMinimizer._jacobian` + `solve(J,-f)`；`ReducedHamelGibbsSolver._reduced_jacobian` + `solve` | **一致** | 增广 vs 降阶两种 Jacobian 组装 |
| **(A.28)** | \(NE\) 个线性方程求 \(\Delta\lambda\) | 两者 | 同上 `np.linalg.solve` / `lstsq` | **一致** | 降阶 solver 固定 \(\lambda_{ref}=0\) |
| **(A.29)** | (A.28) 整理后的中间式 | — | 未单独实现 | **缺失（可接受）** | 直接实现 (A.26)+(A.28)，跳过显式 (A.29) |
| **(A.30)** | 阻尼更新 \(\lambda+\eta\Delta\lambda\) | 两者 | `GibbsMinimizer.solve` 回溯 `eta`；`ReducedHamelGibbsSolver.solve` 同理 | **一致** | 默认 `feasible_backtracking`，非论文唯一策略 |
| **(A.31)** | 阻尼细节 / 步长约束 | 两者 | `omega0`、`eta*=0.5` 循环 | **部分一致** | 工程回溯，未逐字复现 (A.31) |
| **(A.33)** | 总摩尔数 \(N\) 回代 | `GibbsMinimizer` | `N=exp(ln_N)`；`ReducedHamelGibbsSolver` | `_recover_moles_from_ratio` 后 `sum(n)` | **一致** | 两种参数化 |
| **(A.34)** | Holub & Vonka 收敛判据 | 两者 | `tol_mode="holub_vonka_a34"`（增广 `n_j` / 降阶 `q_j`）；默认 `l2_residual` + 降阶 `closure_rel_tol` | **部分一致** | (A.34) 分式已实现为可选开关；默认仍为工程 L2/松弛 |
| **(A.35)–(A.40)** | 起始 \(y_j\) / 泰勒 + 守恒联立 Newton | `GibbsMinimizer` | `hamel_a40_initial_guess` + `initial_guess_mode="hamel_a40"` | **部分一致** | (A.35) y^0 + (A.38) N^0 + 线性化 (A.40) 单步；默认仍为 `taylor_linearized` |
| **(A.35)–(A.40)** | 起始值（降阶路径） | `ReducedHamelGibbsSolver` | `hamel_a40` → λ 种子 + `_coarse_lambda_seed` 网格 | **部分一致** | (A.40) λ 种子可选；网格搜索仍为工程扩展 |

---

## 2. 两套 Solver 分工

| 路径 | 配置 | 方程形态 | 典型用途 |
|------|------|----------|----------|
| **增广 Newton** | `major_gibbs_solver_mode=augmented` | (A.18)+(A.20)，\(NE{+}1\) 未知量 \((\lambda, \ln N)\) | `GibbsMinimizer`；Cell 微量组分 `minor_species.py` |
| **降阶 Newton** | `major_gibbs_solver_mode=hamel_reduced`（thesis 默认） | (A.21)–(A.26) 消元后仅 \(\lambda\) | `ReducedHamelGibbsSolver`；Vorab major / pyrolysis Gibbs seed |

---

## 3. 旧编号 → Hamel 编号（代码迁移对照）

历史 docstring 使用 `A1-x` 内部编号，与论文 `(A.x)` **非一一对应**。迁移对照：

| 代码旧标号 | 建议 Hamel 锚点 | 函数 |
|------------|-----------------|------|
| `A1-8` | **(A.15)–(A.18)** / **(A.23)** | `GibbsMinimizer._calc_moles`；`ReducedHamelGibbsSolver._calc_q` |
| `A1-9` | **(A.18)** | `GibbsMinimizer._residuals` 元素块 |
| `A1-10` | **(A.20)** | `GibbsMinimizer._residuals` 归一块 |
| `A1-12`, `A1-13` | **(A.27)–(A.28)** | `GibbsMinimizer._jacobian` |
| `A1-14` | **(A.35)–(A.40)** | `GibbsMinimizer._initial_guess` |

---

## 4. 工程扩展项（非 Hamel 编号）

| 项 | 位置 | 判定 | 说明 |
|----|------|------|------|
| `feasible_backtracking` / `legacy_half` | `GibbsMinimizer.solve` | **有意偏离** | 可行域回溯；默认非 Holub 原文唯一阻尼 |
| `_coarse_lambda_seed` 网格 | `ReducedHamelGibbsSolver` | **有意偏离** | `hamel_a40` 种子后仍保留网格鲁棒性 |
| `_least_squares_refine` (TRF) | `ReducedHamelGibbsSolver` | **有意偏离** | Newton 失败后的 trust-region 兜底 |
| `_reduced_jacobian_fd` | `ReducedHamelGibbsSolver` | **有意偏离** | 解析 Jacobian 奇异时 FD 回退 |
| `closure_rel_tol` / `residual_rel_tol` | `ReducedHamelGibbsSolver.solve` | **有意偏离** | 仅 `tol_mode=l2_residual` 时启用 `converged_by=closure_relaxed` |
| 残差行缩放 `_residual_row_scale` | `ReducedHamelGibbsSolver` | **有意偏离** | 大数量级元素池数值稳定 |
| `infer_lambda_from_mole_guess` | `ReducedHamelGibbsSolver` | **工程扩展** | Vorab NNLS 初值反推；物理等价于 (A.35) 邻域 |

---

## 5. 统计与审计结论

| 判定 | 数量（含重复路径） |
|------|-------------------|
| 一致 | 10 |
| 部分一致 | 8 |
| 缺失（可接受 / 有意） | 2（(A.1) 目标不显式；(A.29) 中间式跳过） |
| 有意偏离（工程扩展） | 7 |

**总体结论（a1-gibbs scope）**：

- **算法主干同构**：元素守恒 + 理想气体化学势 + Newton–阻尼迭代在两条 solver 路径上均已实现；thesis 主路径 `hamel_reduced` 对应 Hamel 降阶链 **(A.21)–(A.28)**。
- **未完全同构**：工程扩展（网格种子、TRF 兜底、残差行缩放）仍保留；**生产默认已为 Hamel 初值/判据**（2026-07-08）。
- **P2 已闭合（2026-07-08）**：生产默认 **`initial_guess_mode=hamel_a40`** + **`tol_mode=holub_vonka_a34`**（`hamel_parity_defaults.py` + 两 solver）；`holub_equilibrium_residual` / `closure_relaxed` 为数值回退。
- **争议状态**：由「部分一致、中间式未标注」降为 **「主干同构 + 初值/判据默认 Hamel（有证据）」**。

---

## 6. 建议下一步（P1 剩余）

1. ~~将 docstring 全面改用 `(A.x)`~~（已完成）
2. ~~(A.34) `tol_mode="holub_vonka_a34"`~~（**生产默认**，2026-07-08）
3. ~~`vorab_gibbs_seed.py` 标注 (A.2)/(A.3)~~（已完成）
4. 回归：`pytest tests/test_gibbs_minimizer.py tests/test_gibbs_hamel_reduced.py` + `sanity_checks.py`
5. ~~可选：初值 (A.35)–(A.40) 分步同构~~（**生产默认 `hamel_a40`**，2026-07-08）
6. 可选：非 thesis 水力链 parity（`u_d` 默认路径 / Mori-Wen vs Heinbockel）

**自动化**：`python3 scripts/audit_hamel_consistency.py` → PASS（不含 A1 数值锚点）
