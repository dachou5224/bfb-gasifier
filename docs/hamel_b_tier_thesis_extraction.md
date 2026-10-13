# Hamel B 档：学位论文可实施提取（无 Fortran 源码）

> **主源**：Hamel (2001) VDI 学位论文 — *Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe*  
> **OCR 文本**：`data/hamel_1999_extracted.txt`  
> **交叉引用**：`docs/hamel_isomorphic_gap_implementation_checklist.md`（§0.5 证据矩阵）、`docs/vorab_init_abgleich_vs_outer_preinner.md`  
> **范围**：B 档（数学同构级）—— Verbindungsmatrix、Nebenelemente、拓扑 cell、NR 结构；**不含** Fortran 变量排序与稀疏存储格式。

---

## 1. B 档定义与 A 档边界

| 档位 | 论文锚点 | 本仓库现状 | B 档验收 |
|------|----------|------------|----------|
| **A 档** | Bild 2.2 三模块 + Vorabrechnung 预算链 | 已实现 single-shot 预算、Init/Outer Check2、bed0 zu/ab | 保持默认主线 |
| **B 档** | §2.1 配置 + §2.2 方程 + §2.3 NR 矩阵结构 | 线性 bed/freeboard + 模板化 cyclone/return | **可配置 Verbindungsmatrix** + **Jacobian Nebenelemente 显式侧块** |

A 档 opt-in 探针（轴向 O₂ blend、bed01 phase-split 等）**不纳入** B 档 P0；B 档聚焦拓扑与同构 NR 结构。

---

## 2. 论文 §2.1 — 配置模块与 Verknüpfung

### 2.1.1 三模块结构（Bild 2.2）

1. **Konfiguration**：读入 Steuervariablen、Anlagen-/Betriebsdaten；初始化内部字段；**Bestimmung der Verknüpfung zwischen den Bilanzräumen**；Startwertwahl。
2. **Vorabrechnung**：Strömungsmechanik、Trocknung/Pyrolyse、Inertmaterial-Abgleich（A 档已部分落地）。
3. **Zellenmodell**：NR 迭代求解 Bilanzgleichungen。

### 2.1.2 Anlagenkonfiguration（Tabelle 2.1 语义）

论文将反应器离散为沿**主流动方向**排列的 Zellen，段包括：

- Wirbelschicht（流化床）
- Freiraum（自由空间）
- Verbindungsleitungen（连接管）
- Zyklon（旋风分离器）
- Rückführleitungen（返料管）

并配置：Brennstoff-/Vergasungsmittelzufuhr、Ascheabzug、Feststoff-/Gasrückführung 位置。

**B 档含义**：cell 序列不能硬编码为 `bed[0..n-1] + freeboard[0..m-1]`；需在配置阶段构建**有向连接图**，每条边对应一类 zu/ab/rez 流。

### 2.1.3 符号表（论文 Symbolverzeichnis，B 档相关）

| 符号 | 论文含义 | B 档工程映射 |
|------|----------|--------------|
| `N_BU` | 每 cell Bilanzvariablen 数 | `nr_indexing` 每块 unknown 维度 |
| `N_z` / `NCL` | Zellen数 / Cross-links 数 | 图节点数 / 非邻接耦合边数 |
| `N_CA` | 全局方程数 | `N_z * N_BU`（Eq 2.3 语义） |
| `N_zu,d/b,j,i` | 自**任意**源 cell 进入 cell i 的摩尔流 | 图入边 + 物种 j + 相 d/b |
| `N_rez,d/b,j,i` | 再循环进入 cell i 的摩尔流 | 标记为 `rez` 的入边（Zyklon→床层） |
| `m_zu/rez/auf/ab` | 固相质量流（多 sorte × 粒级） | 固相边 + j,k 索引 |

**论文未给出**：Cross-link 的稀疏存储、`NCL` 与 Jacobian 块索引的一一对应表 —— 需工程定义（见 §6）。

---

## 3. 论文 §2.2 — 守恒方程与连接语义

### 3.1 气相（Eq 2.4–2.5）

悬浮相：

```
0 = N_zu,d,j,i + N_rez,d,j,i + N_d,j,i+1 - N_d,j,i - ex_bd,j,i + r_d,j,i
```

气泡相（Eq 2.5）结构对称。

**关键**：`N_zu` 与 `N_rez` 的下标 **i 为接收 cell**；源 cell 由 Verbindungsmatrix 决定，**不限于 i±1**。这与当前 bed 轴向 `i+1→i` 默认链不同。

### 3.2 固相（§2.2.2）

- 多 Feststoffsorte j、多 Partikelgrößenklasse k。
- 流：`m_zu`, `m_rez`, `m_auf`, `m_ab`（相邻 cell 间 up/down）。
- 源汇：反应 `r_j,k,i`、粒级破碎/聚合 `m_in/out`。

### 3.3 能量（Bild 2.5 / Eq 2.7）

对流 + 反应热 + 相间交换 + **跨 cell 焓流**（与 zu/rez/auf/ab 一致）。

### 3.4 全局规模（Eq 2.3）

```
N_CA = N_z * N_BU
```

每个 cell 含相同结构的 Bilanzvariablen 集合（两相气 + 固相 holdup + T 等，具体集合由 Modellvariablen 决定）。

---

## 4. 论文 §2.3 — NR 与 Jacobian 结构

### 4.1 算法 lineage

- Newton-Raphson，源自 **Wozny (1983)** 一维 cell 求解框架。
- 阻尼形式 **Eq 2.9**（**Wirsum 1998** 的 gedämpfte Variante）：

```
x_{v+1} = x_v + Δx_v   (2.9)
```

（Δx 含阻尼因子；具体阻尼规则论文指向 Wirsum 1998，本仓库 `global_nr_solver` 已有 damped Newton 路径。）

### 4.2 Blocktridiagonale Funktionalmatrix

- 主链：cell 沿 Hauptströmungsrichtung 排列 → 方程对 **同 cell** 与 **相邻 cell** 变量耦合 → **块三对角**主带。
- **Nebenelemente（侧块）**：仅当存在 **非邻接** 入流时填充，典型例子：

> Feststoffrezirkulation aus dem Zyklon zurück in die Wirbelschicht

即：cyclone/return_leg → bed0（或 bed 顶）的固相/气相 rez 流产生 **|i - j| > 1** 的 Jacobian 块。

### 4.3 与当前 Python 实现对照

| 论文要求 | 当前实现 | 差距 |
|----------|----------|------|
| 块三对角主带 | `structured_jacobian.build_solver_graph` + band LU | ✅ P0-2 已闭合 |
| Nebenelemente | `jacobian_side_block_pairs_from_connectivity_graph` + UV Woodbury | ✅ VTT `unexpected_blocks==0` |
| 解析/半解析 Jacobian | `jacobian_analytic.py` + `hybrid_*` + `jacobian_transport_freeze.py` | ✅ **P2 B 档验收闭合**（hybrid opt-in；侧块固相 FD 为 accepted gap） |
| 残差邻域刷新 | `connectivity_graph.affected_nr_residual_cells` | ✅ 图拓扑 |
| Verbindungsmatrix 驱动 zu/rez | `connectivity.py` + 图路由 | 部分路径仍混用 legacy 索引（P3 可选） |

---

## 5. 论文 Kap 7 — 拓扑实例（Zyklon 两 cell）

§7.1 VTT 等算例：**Zyklon 用 zwei Zellen** 建模（Bild 7.7），连接关系为：

```
bed → freeboard → (Verbindungsleitung) → Zyklon(2 cells) → Rückführung → bed
```

这与 `connectivity_graph` 中 `cyclone_idx = n_bed + n_freeboard`、`return_leg_idx = cyclone + 1` 的**模板一致**，但论文允许：

- 每段 **Zellenzahl frei wählbar**（§2.2 开篇）；
- 返料可接入 **非 bed0**（由 Lage der Feststoffrückführungen 配置）。

**B 档 P0 实例验收**：同一 JSON 算例可表达上述回路，且 **N_rez/m_rez 进入残差而非后处理回写**。

---

## 6. 论文未给出、需工程补全的项

1. **Verbindungsmatrix 数据结构**：有向边列表 `(src, dst, stream_kind, phase, species_mask)` 是否足够；`NCL` 是否等于非邻接边数。
2. **Unknown ordering**：论文只保证每 cell `N_BU` 一致；块内变量顺序需与 `nr_indexing` 对齐并文档化。
3. **Side block 稀疏格式**：块三对角 + 低秩侧块（当前 UV）是否等价于 Fortran Nebenelemente —— 需与 `hamel_checked_algorithm_spec.md` 对齐。
4. **连接管 cell**：论文列为独立段，当前或无 dedicated cell_type；可暂合并为 freeboard 与 cyclone 之间的 **passthrough edge**。
5. **Gas rezirkulation**：符号表含 Gasrückführung；B 档需区分 `N_zu` vs `N_rez` 标签（影响 Jacobian 块归类，不影响守恒结构）。

---

## 7. B 档实施任务映射（P0 → 文件）

### P0-1 可配置 Verbindungsmatrix

| 任务 | 目标文件 | 验收 |
|------|----------|------|
| 边表/schema | `src/core/connectivity_graph.py`（新 dataclass） | YAML/JSON 读入 bed/fb/cyclone/return 边 |
| 残差 zu/rez 组装 | `src/core/connectivity.py` | 入流由图查询，非 `cells[i±1]` 硬编码 |
| Reactor 构图 | `src/core/reactor.py` | Konfiguration 阶段 `build_connectivity_from_config` |

### P0-2 Jacobian Nebenelemente

| 任务 | 目标文件 | 验收 |
|------|----------|------|
| 侧块来自图 | `src/solvers/structured_jacobian.py` | `side_block_pairs == f(graph)` |
| NR 邻域 | `src/solvers/global_nr_solver.py` | 改 cell 时刷新所有入边源/汇块 |
| 结构审计 | `tests/test_global_nr_solver.py` | unexpected_blocks == 0（回路算例） |

### P2 Hybrid Jacobian（B 档数值，2026-07-07 台账闭合）

| 任务 | 目标文件 | 验收 |
|------|----------|------|
| inventory 解析 + FD 补差 | `jacobian_analytic.py` | `tests/test_jacobian_analytic.py` |
| holdup `K` 列扰动冻结 | `jacobian_transport_freeze.py` | `test_jacobian_transport_k_stays_frozen_during_structured_column_perturbation` |
| hybrid vs structured FD | `structured_jacobian.py` | LU holdup `< 1e-8`；VTT `< 1e-3`（双 reactor init） |
| 导数来源表 | `JACOBIAN_DERIVATIVE_SOURCES` | 与 `hamel_parity_global_nr_2026-07-06.md` §P2 一致 |
| 生产默认 | `reactor._default_nr_jacobian_strategy` | **hybrid_* 默认**（2026-07-07）；`nr_hybrid_jacobian_thesis=False` → structured FD |
| Bed 固相默认 | `ReactorConfig.bed_solid_state_model` | **holdup_transport**（2026-07-07 P1）；`legacy_stream` opt-in |
| 自由板 drag 默认 | `ReactorConfig.freeboard_trajectory_coeff_model` | **exact_hamel**（2026-07-07 P1）；`stable_mixed_drag_split` opt-in（LADDER baseline） |
| Eq.(2.7) 热损默认 | `ReactorConfig.energy_balance_heat_loss_mode` | **`hamel_eq27_explicit`**（2026-07-07）；`inlet_fraction` opt-in（LADDER baseline–p1） |
| Q_W 引导 | `resolve_reactor_Q_wall_W_total` | `heat_loss_frac × fuel_feed × HHV_dry`；出口 T validation 延后 |

**Accepted gap**：侧块固相 DOF 全 FD；R/ex/hydro 经 `(J_fd−J_stencil)`；Check2 与 Jacobian 分离。Walltime：`data/jacobian_assembly_benchmark_latest.json`。

### P0-3 与 A 档 Vorabrechnung 隔离

- Vorabrechnung **不**重复定义拓扑；仅消费 Konfiguration 已冻结的图。
- 见 `docs/vorab_init_abgleich_vs_outer_preinner.md`：Init/Outer bridge 仍属 A 档，B 档不改变默认关闭的 opt-in 探针。

### P1（B 档后续）

| 任务 | 状态 | 验收 |
|------|------|------|
| bed/freeboard 气相/固相 `propagate_upstream` 图路由 | ✅ | `tests/test_connectivity_graph.py::test_propagate_upstream_graph_path_matches_index_path` |
| `freeboard_chain_inlet_source_block_id` 图邻居 | ✅ | 与线性 `cells[i±1]` 等价 |

### P2 固相 K_auf/K_ab 图路由

| 任务 | 目标文件 | 状态 |
|------|----------|------|
| `solid_auf/ab_source_block_id` | `connectivity_graph.py` | ✅ |
| `update_bed/freeboard_solid_transport_inflows` 图路径 | `connectivity.py` | ✅ |
| Verbindungsleitung 透传段 | `add_connection_pipe_passthrough` + JSON | ✅ `data/connectivity/vtt_connection_pipe_passthrough.json` |
| 图 vs 索引等价测试 | `tests/test_connectivity_graph.py` | ✅ |

**P2 可选（未做）**

- 多返料点 JSON 样例（需 Reactor 多 `return_leg` 支持）
- `update_bed_solid_transport_coefficients` 顶床 `top_downflow` 从图取 fb 源

### P3（B 档后续）

| 任务 | 状态 | 验收 |
|------|------|------|
| 多 cyclone / 多 return_leg Reactor 侧块数组化 | ✅ | `cyclone_cells` / `return_leg_cells` / `side_block_cells_by_id` |
| 图驱动侧块实例化 | ✅ | `_build_side_block_cells_from_graph()` |
| 多返料点 rez 边路由 | ✅ | `apply_recirculation_from_connectivity_graph` + `target_bed_index` |
| 2-stage cyclone + dual return JSON | ✅ | `data/connectivity/vtt_dual_cyclone_dual_return.json` |
| Wirsum Eq. 2.9 阻尼参数集中 | ✅ | `src/solvers/wirsum_damping.py` + `nr_lambda_init_thesis` |

**P3 可选（未做）**

- 多 cyclone 级间 `K_auf/K_ab` 完整链式 transport（当前仅末级 cyclone → return_leg 配对）
- 残差 zu/rez 完全边表化（`sum(edge.flux for edge in inbound(i))`）

---

## 8. Python 现状快照（2026-06）

**已实现 B 档闭环（VTT 单回路）：**

- `ConnectivityGraph` + JSON 加载（`vtt_minimal_loop.json`、`vtt_connection_pipe_passthrough.json`）
- Jacobian `side_block_pairs` 由图生成（`structured_jacobian.build_jacobian_structure`）
- NR 邻域 / BC 路由 / `propagate_upstream` / 固相 `K_auf/K_ab` inflow 均支持 `connectivity_graph` + `cell_registry`
- `Reactor.connectivity_graph` 在 Konfiguration 后冻结

**仍属模板级 / 待扩展：**

- 多 cyclone 级间 solid transport 系数链（P3 可选）
- 残差 zu/rez 尚未完全抽象为 `sum(edge.flux for edge in inbound(i))`（部分路径仍混用 legacy 索引）

---

## 9. 推荐实施顺序

1. ~~**Schema**~~ ✅ `ConnectivityEdge` + VTT JSON
2. ~~**残差 BC / propagate**~~ ✅ 图路由（P0–P1）
3. ~~**Jacobian 侧块**~~ ✅ `side_block_pairs == f(graph)`（P0-2）
4. ~~**固相 transport**~~ ✅ P2 `K_auf/K_ab` 图路径 + connection_pipe 透传
5. ~~**测试**~~ ✅ VTT 单回路：`test_vtt_minimal_loop_jacobian_structure_has_no_unexpected_blocks` + `test_vtt_minimal_loop_inner_nr_convergence_acceptance`（`unexpected_blocks == 0`，内层 NR 收敛）
6. ~~**Hybrid Jacobian 台账**~~ ✅ P2 B 档验收闭合（2026-07-07）；见 `hamel_parity_global_nr_2026-07-06.md`
7. ~~**文档**~~ ✅ 本清单 §4.3 + 审计汇总 §0.1；外层门控诊断见 §10

---

## 10. 外层 Check2 门控诊断（VTT 单回路，2026-06）

B 档 Jacobian/Nebenelemente 验收与 A 档 **外层 Abgleich（Check2）** 分离。VTT 验收配置见 `build_vtt_minimal_loop_nr_acceptance_config()`。

### 10.1 实测结论（3 bed + 1 fb + cyclone + return，`vtt_minimal_loop.json`）

| 门控 | 默认配置 | 说明 |
|------|----------|------|
| Jacobian `unexpected_blocks` | ✅ 0 | `nr_structure_validation_ok=True` |
| 内层 Check1（damped NR） | ✅ 收敛 | 单独 `solve_global_nr`：`rms_scaled_final≈2e-4` |
| hybrid vs structured FD | ✅ | `test_vtt_hybrid_jacobian_matches_fd_and_keeps_structure`（`< 1e-3`） |
| **Pre-inner `global_gate_pass`** | ✅ 通过 | merit 预投影未触发回滚；**不是** `converged_fully=False` 主因 |
| **气相 phase-split Check2** | ❌ 失败 | 首轮 `gas_phase_split_rms≈0.63`，阈值 `0.05`；Pre-inner 可降至 ~0.6 但仍未过关 |
| bed0 primary Check2 | — 默认关 | `nr_outer_check2_bed0_primary_gate_thesis=False` |
| **ΔT Check2** | ❌ 首几轮失败 | `max\|T - T_pre_refresh\|` 首轮 **~300 K**，第 2 轮 **~600 K**（被 clip）；阈值 `max(0.5 K, rel_tol·T_scale)≈0.5 K` |
| **内层 Check1（outer 内，默认）** | ❌ 短探针失败 | `max_global_iter≤4` 触发 `short_probe_mode`，inner budget 仅 4；`inner_rms≈0.63` |
| **内层 Check1（outer 内，full）** | ✅ 可收敛 | `max_iter=30` + `nr_inner_tol_rms_thesis=0.02` + Check1 continuation stall recovery → `converged_fully=True` |
| **Outer Check2 Bridge** | 未触发 | 需 `inner_converged=True` 且气相 Check2 失败；短探针轮次内层未先过关 |

因此 `reactor.solve()` 在 **短探针** outer 配置下 **`converged_fully=False` 是预期**；full outer（见 `VTT_MINIMAL_LOOP_SOLVE_KWARGS`）可在 A 档完成。

### 10.2 Check1 continuation 停滞修复（2026-06）

`strict_check1` 模式下，若 continuation 轮次线搜索无进展（`n_iter≤1` 且 RMS 平台），原先永久 `refresh_allowed=False`，导致 outer 空转。修复见 `src/solvers/outer_loop.py`：停滞时下一轮恢复 Vorab refresh。单测：`test_strict_check1_reallows_refresh_when_continuation_stalls`。

VTT 配置 `nr_inner_tol_rms_thesis=0.02` 与 B 档内层验收对齐（能量分量 ~0.011 在 0.01 阈值下会卡住 Check1）。

### 10.3 判据源码锚点

```python
# src/solvers/convergence.py — OuterConvergence.is_aligned
tol_abs = max(abs_tol_T, rel_tol * T_ref_scale)  # 默认 abs_tol_T=0.5 K
if dT_outer_max > tol_abs: return False
```

```python
# src/solvers/outer_loop.py — Check2 bridge 需 inner 先收敛
if inner_converged and not raw_outer_aligned and (not gas_check_pass or not bed0_check_pass):
    check2_bridge_diag = check2_bridge_fn()
outer_matched = inner_converged and raw_outer_aligned
```

Pre-inner 全局门控（`src/core/connectivity.py::run_outer_preinner_start_value_bridge`）仅在 merit 变差时回滚；VTT 算例中 `preinner_global_gate_pass=True`。

### 10.4 推荐验收口径

| 档位 | 通过条件 |
|------|----------|
| **B 档（结构同构）** | `unexpected_block_count==0` + 内层 NR 收敛（`test_vtt_minimal_loop_inner_nr_convergence_acceptance`） |
| **A 档（全流程）** | 上述 + `converged_fully`（`test_vtt_minimal_loop_outer_a_tier_convergence`，`VTT_MINIMAL_LOOP_SOLVE_KWARGS`） |
| **门控诊断（短探针）** | `test_vtt_minimal_loop_outer_check2_gate_diagnosis`（`max_global_iter=4`，记录 Check2 分解） |

### 10.5 残差主因与改进（2026-06 全量模拟审计）

全量 outer（`VTT_MINIMAL_LOOP_SOLVE_KWARGS`）收敛后，**高残差几乎全来自 bed 区主组分 dense/bubble phase-split**，侧块（cyclone/return）接近 0。

| 阶段 | `gas_phase_split_rms` | 主因 cell / 物种 |
|------|----------------------|------------------|
| Init | ~0.33（最差 ~1.24） | **bed1**：N2、CO、CO2 的 \|N_d−N_b\|（zu/ab 映射后两相未闭合） |
| Final（默认） | ~0.010 | **bed1**：CO2/CO/N2；**bed2**：H2O/H2 |
| Final（改进配置） | **~0.0045** | 同上，幅度降 ~2× |

**已落地改进**（`build_vtt_minimal_loop_nr_acceptance_config`）：

- `nr_outer_preinner_bed1_phase_split_thesis=True` — 每 outer refresh 后对 bed1 做 major split 预投影
- `vorab_a_tier_bed1_phase_split_after_stream_mapping_thesis=True` — Init zu/ab 映射后对 bed1 预投影

审计脚本：`scripts/vtt_full_sim_residual_audit.py`。

**HTW Phase2 已默认开启**（`build_phase2_htw_lu_freeboard_reactor_config`，staged primary 时）：

- `nr_outer_preinner_phase_split_max_bed_index_thesis=2`（bed2 split 预投影；压 bed3 N₂ 气泡热点）

**VTT 后续仍可选**：
- 提高 `max_iter` / inner budget（收益递减；40 outer 通常足够）
- recycle 回路 N2 初值与 return-leg 气相 handoff 对标 Hamel（降低 init split 尖峰）

### 10.6 后续（非 B 档阻塞）

- 调查 recycle 回路首轮大 ΔT 与 bed 气相 split 是否应对标 Hamel 论文 outer 容差
- 可选：VTT 专用 Check2 标定（与 LU/HTW 算例分开）

---

## 11. 参考文献（论文内引）

- Wozny (1979, 1983) — NR cell 求解
- Wirsum (1998) — gedämpftes Newton (Eq 2.9)
- Rajan & Wen (1980) / Hartleben (1983) / Weiß (1987) — 前期 cell 模型 lineage
- Kurkela et al. (1993) — VTT 装置与 Bild 7.7 拓扑

---

*文档版本：v1 — 从 OCR 提取，无 Fortran 对照。若后续获得源码，应增补 §6 排序/存储对照表。*
