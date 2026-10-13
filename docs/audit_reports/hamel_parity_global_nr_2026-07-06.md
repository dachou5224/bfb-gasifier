# Hamel Parity 审计 — global-nr / Jacobian — 2026-07-06（2026-07-07 P2 台账闭合）

## 范围

- Hamel (1999) §2.3 NR 矩阵结构（块三对角 + Nebenelemente 侧块）
- B 档 P0-2：`side_block_pairs == f(connectivity_graph)`、NR 邻域、结构验收
- B 档 P2：hybrid 半解析 Jacobian（inventory stencil + FD 补反应/交换）
- 实现：`src/solvers/structured_jacobian.py`、`src/solvers/jacobian_analytic.py`、`src/solvers/jacobian_transport_freeze.py`、`src/solvers/global_nr_solver.py`、`src/core/connectivity_graph.py`

**交叉引用**：`docs/hamel_b_tier_thesis_extraction.md` §4、§P0-2、§10

---

## 台账闭合口径（2026-07-07）

| 层级 | 含义 | 判定 |
|---|---|---|
| **P0 结构** | 块三对角 + 图驱动 Nebenelemente + Woodbury 线性步 + `unexpected_blocks==0` | **已闭合 ✅** |
| **P2 数值 hybrid** | 论文无 Fortran 逐项导数；hybrid = 解析 inventory/侧块 T + `(J_fd−J_stencil)` 补 R/ex/hydro | **B 档验收闭合 ✅**（见下表与测试） |
| **生产默认** | 默认 `hybrid_*`（`nr_hybrid_jacobian_thesis=True`）；`False` 回退 structured FD | **2026-07-07 已切换** |
| **Accepted gap** | 侧块固相 DOF 全 FD；间接 hydro/ex 无 Hamel 闭式；外层 Check2 与 Jacobian 分离 | **记录接受，不阻塞 P0/P2 台账** |

Hamel 原文未给出全局 NR Jacobian 的完整解析式；**「完全逐项解析同构」不是本项目 B 档目标**。B 档 Jacobian parity 闭合定义为：**结构同构 + hybrid 与 structured FD 在声明稀疏图上一致（容差内）**。

---

## P0-2 Nebenelemente 判定总表

| 论文/台账项 | Hamel 含义 | 实现 | 判定 |
|---|---|---|---|
| 主链块三对角 | 相邻 cell 方程–变量耦合 | `build_solver_graph` + `band_block_pairs`（\|i−j\|≤1） | **一致（结构）** |
| Nebenelemente 侧块 | 非邻接 rez/zu 入流产生 \|i−j\|>1 块 | `jacobian_side_block_pairs_from_connectivity_graph` → `build_jacobian_structure` | **一致（结构）** |
| 侧块进入线性步 | Woodbury：`J = J_band + U V^T` | `StructuredJacobian.to_uv` + `solve_structured_linear_step` | **一致（离散布局）** |
| Side-block `T` 进入 NR | Nebenelemente 能量闭包行 | `nr_indexing.has_temperature_var` 含 cyclone/return_leg | **一致** |
| band+side 全 BC 组装 | 回收耦合对 local-BC 剪枝敏感 | `band_plus_side_elements_structured` 禁用 local callback | **工程闭合** |
| Jacobian 数值（hybrid） | 论文无逐项闭式 | `jacobian_analytic.py` + `hybrid_*` + `jacobian_transport_freeze.py` | **B 档 hybrid 闭合**（见 §P2） |

---

## P2 Hybrid 解析覆盖（2026-07-07）

公式：`J = J_analytic_stencil/scale + (J_fd − J_stencil_fd)`；5% stencil guard 不一致时回退该行 FD。

| 残差组 | Hamel 引用 | 解析项 | 状态 |
|---|---|---|---|
| 气相 dense/bubble 对角 | Eq. (2-4)/(2-5) | `−1` inventory | ✅ hybrid |
| 气相轴向 band 入流 | Eq. (2-4)/(2-5) | 下游对上游 `+1` | ✅ hybrid |
| 固相 legacy band | Eq. (2-6) | `lower_frac` / `upper_frac` | ✅ hybrid |
| 固相 holdup 对角 | Eq. (2-6) | `−(K_auf+K_ab)` 或 anchor `−1` | ✅ hybrid（需 K 冻结，见下） |
| 固相 holdup band | Eq. (2-6) | `+K_auf_col` / `+K_ab_col` | ✅ hybrid（需 K 冻结） |
| 床层能量 `T` | Eq. (2-7) | `−∂H_out/∂T`（`calc_bed_energy_outflow_cp_dot`） | ✅ hybrid |
| 侧块能量 `T` | Eq. (2-7) Nebenelement | `−cp_dot` | ✅ hybrid |
| 气/固 R、ex_bd、迁移 | Eq. (2-4)–(2-6) | 无闭式 | **FD 补差**（accepted） |
| 床层 `T` 的 hydro/间接耦合 | Eq. (2-7) | 无单独闭式 | **FD 补差**（accepted） |
| 侧块固相 `m` | Eq. (2-6) side | 无闭式 | **全 FD**（accepted） |

### `K_solid` 列扰动冻结

`K_auf`/`K_ab` 不在 packed NR 向量中；列 FD 时 `apply_bc` 若重算 `K(m)` 会破坏 holdup 解析 band/对角。

- 模块：`src/solvers/jacobian_transport_freeze.py`
- 接入：`assemble_jacobian_fd_structured`、`build_jacobian_fd_dense`、`build_jacobian_fd_block_tridiag`
- BC：`nr_boundary._bed_transport_coefficients_are_frozen` / `_side_transport_coefficients_are_frozen`
- 每列扰动前：`unpack_reactor(x) + apply_bc_fn()` 重置 packed 状态

---

## Hybrid 策略与默认

| 策略 | 说明 |
|---|---|
| `band_plus_side_elements_structured` | 纯 FD fallback（`nr_hybrid_jacobian_thesis=False`） |
| `hybrid_band_plus_side_elements_structured` | **生产默认**（VTT/侧块图） |
| `block_tridiag_structured` / `hybrid_block_tridiag_structured` | 纯 bed 链；默认 hybrid |

Reactor 默认：`_default_nr_jacobian_strategy()` → `hybrid_*`（2026-07-07）；设 `nr_hybrid_jacobian_thesis=False` 回退 structured FD。

---

## B 档验收测试

| 门控 | 测试 | 容差 / 条件 | 状态 |
|---|---|---|---|
| `unexpected_blocks == 0` | `test_vtt_minimal_loop_jacobian_structure_has_no_unexpected_blocks` | — | ✅ |
| 内层 NR 收敛 | `test_vtt_minimal_loop_inner_nr_convergence_acceptance`（`@slow`） | merit `rms_scaled < 0.02` | ✅ |
| 图侧块 ⊆ 结构侧块 | `test_jacobian_side_pairs_from_graph_match_expected_vtt_wraparound` | — | ✅ |
| structured vs dense FD | `test_structured_bed_recycle_side_block_matches_dense_fd_on_initialized_lu_state` | `< 1e-8` | ✅ |
| VTT hybrid vs FD | `test_vtt_hybrid_jacobian_matches_fd_and_keeps_structure` | max \|Δ\| `< 1e-3`；独立 reactor | ✅ |
| VTT 内层 NR（hybrid） | `test_vtt_minimal_loop_inner_nr_convergence_with_hybrid_jacobian`（`@slow`） | merit `rms_scaled < 0.02` | ✅ |
| LU 组装 benchmark | `test_jacobian_assembly_benchmark_lu_structured_hybrid_dense_parity`（`@slow`） | hybrid≈structured≈dense | ✅ |
| LU holdup hybrid vs FD | `test_holdup_lu_hybrid_jacobian_matches_fd_and_keeps_structure` | max \|Δ\| `< 1e-8` | ✅ |
| K 冻结 | `test_jacobian_transport_k_stays_frozen_during_structured_column_perturbation` | — | ✅ |
| 解析单元 | `tests/test_jacobian_analytic.py`（16 项） | — | ✅ |

配置：`data/connectivity/vtt_minimal_loop.json`；VTT 策略默认 `hybrid_band_plus_side_elements_structured`（`nr_hybrid_jacobian_thesis=True`）。

---

## Accepted gap（不纳入 B 档 Jacobian 阻塞项）

1. **侧块固相 Jacobian**：cyclone/return_leg 的 `m` 行仍为 structured FD（Hamel 无侧块 holdup 解析表）。
2. **生产默认 FD 回退**：设 `nr_hybrid_jacobian_thesis=False` 可恢复纯 structured FD。
3. **外层 Check2**：`test_vtt_minimal_loop_outer_check2_gate_diagnosis` 与 Jacobian 结构验收分离（§10 VTT 诊断）。
4. **dense / structured / hybrid walltime benchmark**：`scripts/benchmark_jacobian_assembly.py`；结果 `data/jacobian_assembly_benchmark_latest.json`。
5. **VTT hybrid 双 reactor 容差 `1e-3`**：初始化微差；LU 单链可达 `1e-8`。

---

## 代码锚点

| 功能 | 位置 |
|---|---|
| 图 → 侧块对 | `connectivity_graph.jacobian_side_block_pairs_from_connectivity_graph` |
| 图 → NR 邻域 | `connectivity_graph.affected_nr_residual_cells` |
| 结构声明 | `structured_jacobian.build_jacobian_structure` |
| FD 组装 + unexpected 检测 | `structured_jacobian.assemble_jacobian_fd_structured` |
| Hybrid 解析 | `jacobian_analytic.py` |
| K 冻结 | `jacobian_transport_freeze.py` |
| Hybrid 组装 | `assemble_jacobian_hybrid_structured` |
| 策略路由 | `global_nr_solver.build_jacobian_fd` |
| 导数来源表 | `jacobian_analytic.JACOBIAN_DERIVATIVE_SOURCES` |

---

## 自动化

```bash
python3 scripts/audit_hamel_consistency.py
pytest tests/test_jacobian_analytic.py -q
pytest tests/test_global_nr_solver.py -k "jacobian or vtt_hybrid or holdup_lu_hybrid" -q
pytest tests/test_global_nr_solver.py -k "vtt_minimal_loop_inner_nr" -q  # @slow
```
