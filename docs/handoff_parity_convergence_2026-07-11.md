# Parity Audit + Convergence 诊断 — 交接文档

> **读者**：接手 Hamel parity 审计、Phase2 LU 收敛、outer/inner NR 诊断的 agent。  
> **日期**：2026-07-11  
> **前置阅读**：[`docs/CLAUDE.md`](CLAUDE.md) · [`docs/vorab_init_abgleich_vs_outer_preinner.md`](vorab_init_abgleich_vs_outer_preinner.md) · [`docs/optimization_handoff_summary.md`](optimization_handoff_summary.md)

---

## 1. 一句话现状

| 维度 | 状态 |
|------|------|
| **Hamel 方程 parity（P0/P1/P2 台账）** | 结构/Jacobian/A1 Gibbs/水力主链 **已闭合**；出口 T / Fig.7 对标 **仍开放** |
| **Convergence（Phase2 p2 LU）** | **O1 refresh split shock 已修复**（Hamel 语义）；末态 merit 从锁死 0.223 降至 ~0.15–0.18，**仍未收敛** |
| **诊断基础设施** | Init / refresh / outer / inner 全链路 trace 脚本就绪；数据在 `data/phase2_*` |
| **调参** | `regression_rel` 等 outer 早停阈值 **暂未动**（用户明确要求） |

---

## 2. Parity Audit 进展

### 2.1 台账与自动化

| 资源 | 用途 |
|------|------|
| [`docs/audit_reports/hamel_parity_audit_summary_2026-07-04.md`](audit_reports/hamel_parity_audit_summary_2026-07-04.md) | 七 scope 汇总 + 2026-07-06/07/08 闭合更新 |
| [`docs/audit_reports/hamel_parity_global_nr_2026-07-06.md`](audit_reports/hamel_parity_global_nr_2026-07-06.md) | Jacobian P0/P2 hybrid |
| [`docs/hamel_b_tier_thesis_extraction.md`](hamel_b_tier_thesis_extraction.md) | B 档拓扑、VTT Check2 诊断 §10 |
| [`docs/hamel_engineering_bridge_notes.md`](hamel_engineering_bridge_notes.md) | 工程桥接 vs 论文证据边界（§8 bottom split） |
| `python3 scripts/audit_hamel_consistency.py` | 方程锚点扫描 |
| `python3 tests/sanity_checks.py` | 17/17 门禁 |

### 2.2 已闭合（不阻塞收敛工作）

- **P0 Jacobian 结构**：块三对角 + 图驱动 Nebenelemente + Woodbury；VTT `unexpected_blocks==0`
- **P2 hybrid Jacobian**：生产默认 `nr_hybrid_jacobian_thesis=True`
- **Bed Eq.(2-6)**：默认 `holdup_transport`；`legacy_stream` 为 audit opt-in
- **自由板**：默认 `exact_hamel`；ladder baseline/p0/p1 仍用 `stable_mixed_drag_split`
- **A1 Gibbs**：默认 `hamel_a40` + `holub_vonka_a34`
- **Vorab single-shot**：干燥/DAEM 源项 init 后冻结；outer 仅 refresh 水力

### 2.3 Accepted gap（记录接受，勿当 bug 反复修）

- 侧块固相 Jacobian 全 FD；R/ex/hydro 经 `(J_fd − J_stencil)` 补差
- **外层 Check2 与 Jacobian 验收分离**（VTT 上 Check2 可失败而 inner NR 可收敛）
- holdup char/ash mask、holdup guardrails（B-tier）
- 出口温度 / Fig.7 parity **延后**（待能量/反应链闭合）

### 2.4 新发现的不一致点（2026-07-10，**已修代码**）

**问题**：`reactor._snapshot_bottom_gas_inlet_split_from_vorabrechnung()` 在 **每次 outer refresh** 调用 `align_bottom_primary_gas_state_to_inlet_split`，强行改 cell0 的 `N_d/N_b`。

**Hamel 原文（Bild 2.2 / §2.1）**：

- Vorabrechnung → zu/ab **边界流** + 水力系数
- Zellenmodell → NR 未知量 `N_d/N_b`，靠 ex_bd 闭合
- Check2 失败 → **回到 Vorabrechnung**，不是 mid-outer 改 inventory

**修复**（`src/core/reactor.py`）：refresh 路径 **只 snapshot 水力 dense_frac**，不再 align。`align` 仅保留 Init / `apply_vorabrechnung_bed_cell_mapping`。

**证据**：

| 指标 | 修复前 | 修复后 |
|------|--------|--------|
| 1st refresh split | 0.044 → **0.223** | 0.044 → **0.044** |
| 全量 p40 O1 merit | 0.223（1 iter 锁死） | **0.103**（inner 8 iter） |
| 全量 p40 末态 merit | 0.223 | **~0.15–0.18** |

数据：`data/phase2_refresh_split_shock.json` vs `data/phase2_refresh_split_shock_post_fix.json`  
全量：`data/phase2_full_sim_merit_trace_p40_post_fix.json`

---

## 3. Convergence 诊断进展

### 3.1 诊断脚本（本阶段新增/强化）

| 脚本 | 作用 | 默认输出 |
|------|------|----------|
| [`scripts/audit_full_sim_merit_trace.py`](../scripts/audit_full_sim_merit_trace.py) | Init（独立 reactor）+ 标准 `solve()`；outer/inner history | `data/phase2_full_sim_merit_trace.json` |
| [`scripts/audit_refresh_split_shock.py`](../scripts/audit_refresh_split_shock.py) | init → 1st refresh → pre-inner；逐步 cell split 分解 | `data/phase2_refresh_split_shock.json` |
| [`scripts/audit_outer_residual_decay.py`](../scripts/audit_outer_residual_decay.py) | outer/inner 残差衰减、T-cap / merit ES / bed2 等 CLI | `data/phase2_outer_residual_decay_audit_*.json` |
| [`scripts/audit_init_residual_breakdown.py`](../scripts/audit_init_residual_breakdown.py) | init 分项 RMS、top 方程族 | — |
| [`scripts/audit_bed_phase_split_contribution.py`](../scripts/audit_bed_phase_split_contribution.py) | 逐 bed split 方差贡献 | `data/phase2_bed_phase_split_contribution.json` |

**Harness SSOT**：[`scripts/_lu_harness.py`](../scripts/_lu_harness.py) — `LADDER_LEVELS["p2"]`、`phase2_solve_kwargs()`、`build_phase2_reactor()`。

**重要**：全量测试必须用 **`reactor.solve()` 单路径**。`run_phase2_init()` + `solve()` 双 init 会污染轨迹；早期 trace 曾因此误报「好结果 ~0.044」。

### 3.2 Init 典型画像（p2 LU，post-fix）

| 指标 | 值 | 主导 |
|------|-----|------|
| merit | ~0.104 | **energy** |
| gas_phase_split_rms | ~0.044 | — |
| top 方程 | cell10 `T:energy_balance` | F̂ ≈ −0.42 |

### 3.3 Outer / Inner 机制（post-fix p40）

```
O1: refresh 后 split 仍 ~0.044；inner 8 iter → merit ~0.103（split 主导）
O2–O8: 缓升 ~0.103 → ~0.119（energy/split 交替）
O9: +0.03 split 跳变（refresh_skipped + preinner global_gate 回滚）
→ plateau break，9 outer，inner budget 39/200
```

**仍存在的瓶颈**（修复 align 之后）：

1. **Pre-inner global gate 常失败**：bed0 T2 局部接受但 holdup/联合门控回滚（`holdup_flux_mismatch` ~0.72）
2. **Inner NR 难压 split**：O1 后 split merit ~0.10+，T-cap=150K 已生效但不够
3. **Outer plateau 早停**：`no_inner_progress && rms_stalled` 在 ~9–11 outer 退出（非 merit early_stop flag）
4. **bed2** 占 split 方差 ~54%（历史审计）；selective bed2 T2 易 destabilize，默认仍 bed0-only

### 3.4 已实现的收敛机制（代码在树中，部分未充分验证）

| 机制 | 位置 | 默认 | 备注 |
|------|------|------|------|
| Inner split-dominant **T-cap** | `global_nr_solver.resolve_inner_t_max_step_K` | cap=150K, ratio=0.85 | A/B：cap=0 复现 split 追逐；cap=150 抑制 inner RMS 抬高 |
| **Composite merit** outer 早停 | `convergence.OuterMeritEarlyStop` + `outer_loop` | enabled, min_outer=8, reg_rel=**0.25** | 用户说**暂不调** regression_rel |
| T2 pre-inner selective bed2 | `preprojection._select_outer_preinner_phase_split_bed_indices` | max_bed_index=**0** | bed2 全量 T2 有害 |
| T2-only fallback on global gate fail | `preprojection` | strict gate 下极少触发 | |
| `inner_group_history` 钩子 | `global_nr_iteration` + `nr_inner_record_group_history_thesis` | audit 脚本打开 | |

测试：[`tests/test_outer_loop_convergence.py`](../tests/test_outer_loop_convergence.py)（18 过）、[`tests/test_preinner_merit_driven.py`](../tests/test_preinner_merit_driven.py)。

### 3.5 历史 benchmark 参考（align 修复前，仅供对比）

| 运行 | outer | 末态 RMS/split | 备注 |
|------|-------|----------------|------|
| p20 + T-cap v2 | 17 | **0.036** | 最佳工作点（**可能含 double-init 或旧轨迹**；需 post-fix 重跑） |
| p40 无保护 | 40 | 0.107 | outer~25 refresh rebound |
| p40 align 锁死 | 40 | 0.223 | **根因：refresh align**（已修） |

**待办**：用 post-fix 代码 **重跑** `benchmark_phase2_convergence_ladder.py --levels p2 --outer-max 20` 建立新 baseline。

---

## 4. 推荐命令（接手第一天）

```bash
cd bfb-gasifier

# 门禁
python3 tests/sanity_checks.py
python3 scripts/audit_hamel_consistency.py

# Refresh shock 快检（~2 min init + refresh）
python3 scripts/audit_refresh_split_shock.py

# 全量 trace（~7–8 min，p40）
python3 scripts/audit_full_sim_merit_trace.py --max-global-iter 40 \
  --output data/phase2_full_sim_merit_trace_latest.json

# Outer 衰减 + 可选 knob
python3 scripts/audit_outer_residual_decay.py --max-global-iter 20

# 相关单测
pytest tests/test_bottom_gas_inlet_split.py tests/test_outer_loop_convergence.py -q
pytest tests/test_vorabrechnung_a_tier_budget.py -q   # 见 §5 已知失败
```

**Ladder SSOT**：`python3 scripts/benchmark_phase2_convergence_ladder.py --levels p2 --outer-max 20`

---

## 5. 开口项目（按优先级）

### P0 — 收敛主链

- [ ] **Post-fix baseline**：重跑 convergence ladder p20/p40，更新 `data/phase2_convergence_ladder_*.json`
- [ ] **O9+ split 缓升**：refresh_skipped continuation 下 pre-inner gate / holdup 联合门控 — 读 `joint_gate.py` + O9 `outer_rows` in trace JSON
- [ ] **Init energy 0.104**：cell10 `T:energy_balance` 是否与热损模式 / x₀ 一致（`audit_init_residual_breakdown.py`）
- [ ] **确认测试**：`test_march_vorabrechnung_widens_nr_temperature_fence_upper_bound` **1 failed / 39 passed**（2026-07-10）；判断是否回归或 flaky

### P1 — Hamel 语义对齐

- [ ] **p2 改用 `gas_inlet_split_strategy=precalc_hydrodynamic_flux`**（Phase1/LU 已是；p2 仍 `fixed` 0.5）— 见 `config_phase2_freeboard.py`、`hamel_engineering_bridge_notes.md` §8
- [ ] **Check2 失败路径**：Vorab zu/ab remapping（`vorab_cell_abgleich` outer bridge）vs 仅 pre-inner T2
- [ ] **SCRIPTS_INDEX**：登记 `audit_full_sim_merit_trace.py`、`audit_refresh_split_shock.py`、`audit_outer_residual_decay.py`

### P2 — 优化/实验（用户暂禁或待讨论）

- [ ] Outer **composite merit 早停** `regression_rel` 标定（用户：**暂时不用**）
- [ ] **Selective bed2** pre-inner T2（split 方差主因 bed2，但全量/early 有害）
- [ ] 延长 outer 到 40+ **无** split shock 后的真实曲线（post-fix 仅跑到 9 outer）

### P3 — Parity 出口

- [ ] Hamel Fig.7 / Table 7 出口 CO/CO₂/T 对标（`data/hamel_fig7_parity_digitized.json`）
- [ ] `calibrate_phase2_co_co2_window.py` 在 **收敛态** 下重跑（非 plateau 态）

---

## 6. 反模式（勿重复踩坑）

| 反模式 | 后果 |
|--------|------|
| `run_phase2_init()` 后再 `solve()` 当标准测试 | 轨迹异常（8 outer / 0.36 或误报 0.044） |
| outer refresh 里 `align` inventory | split 0.044→0.223（**已修**） |
| 单因素压 RMS 判断 outer 停/行 | 应用 **composite merit** = max(rms, split, energy) |
| Init Abgleich 每 outer 无条件跑 | 破坏 NR 已接受状态 |
| bed2 全量 T2 pre-inner | destabilize；仅 selective 可选 |
| 无 pre-inner 保护地推到 outer 40 | 历史：0.107 split 平台（align 修复前为 0.223 锁死） |

---

## 7. 关键文件索引

| 主题 | 路径 |
|------|------|
| Outer 循环 | `src/solvers/outer_loop.py` |
| Inner 迭代 + T-cap | `src/solvers/global_nr_iteration.py`, `global_nr_solver.resolve_inner_t_max_step_K` |
| Pre-inner | `src/core/connectivity/preprojection.py`, `joint_gate.py` |
| Refresh | `src/core/reactor.py::_refresh_vorabrechnung_sources_for_nr` |
| Bottom split | `src/core/connectivity/gas_inlet.py` |
| Init 链 | `src/workflow/steps/init_precalc_step.py` |
| Check2 bridge | `src/solvers/vorab_cell_abgleich.py` |
| p2 配置 | `tests/validation/config_phase2_freeboard.py`, `_lu_harness.LADDER_LEVELS["p2"]` |
| 诊断技能 | `.cursor/skills/bfb-diagnosing-bugs/`, `bfb-hamel-parity-audit/` |

---

## 8. 数据文件 SSOT（2026-07-10/11）

| 文件 | 含义 |
|------|------|
| `data/phase2_refresh_split_shock.json` | align 修复 **前** refresh shock |
| `data/phase2_refresh_split_shock_post_fix.json` | align 修复 **后** |
| `data/phase2_full_sim_merit_trace_p40_v2.json` | 标准 solve，align 修复前（0.223 锁死） |
| `data/phase2_full_sim_merit_trace_p40_post_fix.json` | 标准 solve，align 修复后 |
| `data/phase2_bed_phase_split_contribution.json` | bed0/1/2 split 方差 |
| `data/phase2_outer_residual_decay_audit_p20_*.json` | T-cap / merit ES 实验 |

---

## 9. 建议接手顺序

1. 读本文 + `vorab_init_abgleich_vs_outer_preinner.md` §5–6  
2. 跑 `audit_refresh_split_shock.py` + `audit_full_sim_merit_trace.py` 确认环境  
3. 重跑 `benchmark_phase2_convergence_ladder.py` 建立 post-fix baseline  
4. 针对 O9 split 跳变做 pre-inner / gate 审计（JSON 里 `outer_rows` + `preinner_bridge_diag`）  
5. 再议 p2 `precalc_hydrodynamic_flux` 与 Check2 remapping  

---

*文档由 2026-07-10 parity/convergence 会话整理；align 修复尚未单独 commit，接手前请 `git diff src/core/reactor.py src/core/connectivity/gas_inlet.py`。*
