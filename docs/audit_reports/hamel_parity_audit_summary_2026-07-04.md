# Hamel 同构 Parity 审计汇总报告

**审计日期**：2026-07-04（初版）；**2026-07-06 闭合更新**见 §0；**2026-07-07 Jacobian P2 台账闭合**见 §0.1；**2026-07-08 文档/ A1 默认**见 §0.4  
**审计范围**：全部 7 个可覆盖 scope（`hydro-core`、`freeboard`、`exchange-cell`、`global-nr`、`drying-pyro`、`kinetics`、`a1-gibbs`）  
**证据链**：`docs/hamel_submodels/00–07` → `specs/hamel_checked_algorithm_spec.md` → `src/` → `docs/hamel_isomorphic_gap_implementation_checklist.md` §0.5  
**自动化锚点**：`python3 scripts/audit_hamel_consistency.py` → **PASS**  
**最小测试门禁**：`python3 tests/sanity_checks.py` → **ALL 17 SANITY CHECKS PASSED**

---

## 0. 2026-07-06 闭合更新（相对 7/4 初版）

| 项 | 原判定 | 现判定 | 证据 |
|---|---|---|---|
| Eq. (3.50) `K_bd` 对流项 | 部分一致（OCR 误识 `u_br/V_b`） | **一致** | PDF page_27；`(3/2)u_{b,r}/d_b`；`mass_transfer.py` |
| Eq. (3.46)–(3.49) 替代 `K_bd` | 缺失 | **已实现** | `calc_kbd_kunii_levenspiel` 等 |
| Eq. (3.43) `v_b` | 部分一致（未暴露） | **一致（默认 v_b=1）** | `bubble_rise_velocity(..., v_b=)` + `ReactorConfig.hydrodynamics_v_b` |
| Eq. (4.3) Stefan 前沿 | 部分一致（湿基分母） | **一致** | `drying.py` 干基 `w_{0,tr}` |
| Eq. (4.5)/(4.7)/(4.8) | 缺失 | **已实现** | `drying.py` |
| Eq. (3.70)–(3.71)/(3.96)/(3.97) | 缺失 | **已实现** | `freeboard.py`, `trajectory_common.py` |
| Eq. (5.45) R8 压力指数 | 部分一致（裁剪） | **一致** | 恢复 `p^((0.5-p/250))` |
| Eq. (5.49) R8 SI | 缺失 | **已实现** | `rate_R8_si()` |
| Eq. (2-6) `legacy_stream` | 部分一致（缺 auf/ab） | **默认 holdup_transport（2026-07-07 P1）** | `bed_solid_state_model=holdup_transport`；`legacy_stream` 为 opt-in 审计分支 |
| `thesis_mode` 自由板 drag | 部分一致 | **全局默认 `exact_hamel`（2026-07-07 P1）** | `hamel_parity_defaults.py` + `ReactorConfig` |
| Vorab single-shot 隔离 | 部分一致 | **P1 已落地** | `vorab_policy.py` + 源项冻结 + init 顺序调整 |
| A1 Gibbs 对照表 | 部分一致 | **P1 已闭合（有证据）** | `docs/audit_reports/hamel_parity_a1_gibbs_2026-07-06.md` |
| Eq. (2.7)/(2.8) 能量/焓 | 部分一致 | **部分一致（P1 已落地）** | `hamel_parity_energy_2026-07-06.md`；Eq.(2.8) 默认一致；Eq.(2.7) 可选显式 Q |
| P0-2 Jacobian Nebenelemente | 部分一致（FD） | **结构一致（P0 已闭合）** | `hamel_parity_global_nr_2026-07-06.md`；VTT `unexpected_blocks==0` |
| 非 thesis 水力链默认 | 部分一致（Cell/calc 默认 mori_wen） | **一致（Hamel 主链为默认）** | `hydrodynamics_defaults.py`；见 `hamel_parity_hydro_non_thesis_2026-07-06.md` |

### 0.1 2026-07-07 Jacobian P2 台账闭合

| 项 | 判定 | 证据 |
|---|---|---|
| P0 Jacobian 结构 | **已闭合 ✅** | 块三对角 + 图驱动 Nebenelemente + Woodbury；VTT `unexpected_blocks==0` |
| P2 hybrid 数值 | **B 档验收闭合 ✅** | 床层气/固 inventory + band；holdup `K·m`（`jacobian_transport_freeze`）；侧块 `T` `−cp_dot`；LU hybrid vs FD `< 1e-8` |
| 生产默认 Jacobian | **hybrid（2026-07-07）** | `nr_hybrid_jacobian_thesis=True` 默认；`False` → structured FD |
| Accepted gap | **记录接受** | 侧块固相全 FD；R/ex/hydro 经 FD 补差；Check2 与 Jacobian 分离 |

### 0.2 2026-07-07 P1：固相 + 自由板默认

| 项 | 判定 | 证据 |
|---|---|---|
| Bed Eq.(2-6) 默认模型 | **holdup_transport（生产默认）** | `ReactorConfig.bed_solid_state_model`；`src/core/hamel_parity_defaults.py` |
| `legacy_stream` | **有意偏离（audit opt-in）** | `bed_solid_state_model=legacy_stream` |
| 自由板轨迹阻力 | **exact_hamel（生产默认）** | `ReactorConfig.freeboard_trajectory_coeff_model`；`segment.simulate_freeboard` |
| `stable_mixed_drag_split` | **有意偏离（audit opt-in）** | `_lu_harness` LADDER baseline/p0/p1 仍显式启用 |
| 自动化 | **PASS** | `audit_hamel_consistency.py` + `tests/test_hamel_parity_defaults.py` |

仍开放（**非阻塞**）：holdup 数值 guardrails（B-tier accepted gap，已文档化）；**出口 T Hamel validation**（待反应/能量校正）。

### 0.2.1 2026-07-07 P1：粒径迁移 + holdup guardrails

| 项 | 判定 | 证据 |
|---|---|---|
| 粒径迁移库存基准 | **per_class_holdup（生产默认）** | `calc_size_migration(inventory_mode=...)`；`ReactorConfig.size_migration_inventory_mode` |
| `total_bed_proxy` | **有意偏离（audit opt-in）** | `size_migration_inventory_mode=total_bed_proxy` |
| holdup char/ash mask | **B-tier accepted gap** | `holdup_transport_char_ash_mask()`；VM/moisture 不参与 K 输运 |
| holdup NR guardrails | **B-tier accepted gap** | `_apply_holdup_transport_numerical_guardrails()`；`disconnected` / `source_without_transport` |

### 0.3 2026-07-07 Eq.(2-7) 热损模式探针

| 项 | 结论 | 证据 |
|---|---|---|
| Phase2 fast 探针 | **H_in 常为负 → `f·max(H_in,0)` 实际热损≈0** | `scripts/audit_energy_heat_loss_mode_compare.py` |
| 显式校准 | **`abs_H_in_audit_bridge`**：`f·|H_in|` 映射至 `Q_wall_W_total` | `resolve_inlet_fraction_calibration_W()` |
| 生产默认 | **`hamel_eq27_explicit`（2026-07-07）** | `DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE`；`inlet_fraction` opt-in |
| Q_W 初值 | **`heat_loss_frac × fuel_feed × HHV_dry`**（未显式设 `Q_wall_W_total` 时） | `resolve_reactor_Q_wall_W_total()` |
| 出口 T validation | **延后** | 待反应/能量守恒校正后再对标 Hamel Fig.7 / Table 7 |

### 0.4 2026-07-08 A1 Gibbs 生产默认

| 项 | 判定 | 证据 |
|---|---|---|
| 初值 | **`hamel_a40`（生产默认）** | `HAMEL_GIBBS_INITIAL_GUESS_MODE`；`GibbsMinimizer` / `ReducedHamelGibbsSolver` |
| 收敛判据 | **`holub_vonka_a34`（生产默认）** | `HAMEL_GIBBS_TOL_MODE`；KKT 平衡时 `holub_equilibrium_residual` 回退 |
| Legacy | **opt-in** | `taylor_linearized` + `l2_residual` |
| 自动化 | **PASS** | `audit_hamel_consistency.py` + `tests/test_gibbs_*.py` |

---

## 1. 各 Scope 判定统计

| Scope | 分册 | 一致 | 部分一致 | 不一致 | 缺失 | 有意偏离 | 主要代码模块 |
|-------|------|------|----------|--------|------|----------|--------------|
| `hydro-core` | `03_hydrodynamics_core_chain.md` | 4 | 6 | 0 | 4 | 0 | `src/physics/minimum_fluidization.py`, `bubble_dynamics.py`, `phase_fractions.py`, `mass_transfer.py` |
| `freeboard` | `04_freeboard_particle_trajectory_analytic.md` | 15 | 4 | 0 | 4 | 0 | `src/physics/freeboard.py`, `src/core/freeboard/` |
| `exchange-cell` | `02_cell_balances_and_exchange.md` | 4 | 2 | 0 | 0 | 0 | `src/core/cell_balances.py`, `cell.py`, `connectivity_graph.py` |
| `global-nr` | `01_global_nr_and_connectivity.md` | 4 | 2 | 0 | 0 | 0 | `src/solvers/global_nr_solver.py`, `structured_jacobian.py`, `connectivity_graph.py` |
| `drying-pyro` | `05_drying_and_daem_pyrolysis.md` | 5 | 6 | 0 | 3 | 0 | `src/thermal/drying.py`, `devolatilization.py`, `core/cell_pyrolysis.py` |
| `kinetics` | `06_kinetics_r1_r11_and_equilibrium_driving.md` | 17 | 6 | 0 | 1 | 0 | `src/kinetics/gas_reactions.py`, `char_reactions.py`, `tar_reactions.py`, `arrhenius.py`, `src/thermodynamics/equilibrium.py` |
| `a1-gibbs` | `07_appendix_a1_gibbs_solver.md` | 4 | 6 | 0 | 0 | 0 | `src/thermodynamics/gibbs_minimizer.py`, `gibbs_hamel_reduced.py` |
| **合计** | — | **53** | **36** | **0** | **12** | **0** | — |

**判定口径**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离（见 gap §x）`

---

## 2. 关键发现与证据

### 2.1 一致项亮点（可直接采信）

| 方程/表 | Scope | 代码位置 | 说明 |
|---------|-------|----------|------|
| Eq. 3.44 `u_{b,r}` | hydro-core | `src/physics/mass_transfer.py:24-34` | 压力修正指数 `-0.15`、`n_b=2.7`、Pa 输入均与论文一致；audit 数值锚点 PASS |
| Eq. 3.42 `λ_b` | hydro-core | `src/physics/bubble_dynamics.py:275-278` | 系数 `280`、指数 `-0.7`、`u_mf/g` 结构一致 |
| Eq. 3.50 `K_bd` 扩散项 | hydro-core | `src/physics/mass_transfer.py:37-61` | `sqrt(144·D_g·ε_mf·u_b/(π·d_b^3))` 与论文一致 |
| Eq. 3.74–3.76 `u_gb` / `β_A` | freeboard | `src/physics/freeboard.py:13-70` | ghost-bubble 衰减链完整实现 |
| Eq. 3.88–3.105 颗粒轨迹 | freeboard | `src/core/freeboard/trajectory_common.py` | Riccati 型 ODE、Δ>0/Δ<0 双分支、高度/速度反函数均实现 |
| Eq. 2-4 / 2-5 气相守恒 | exchange-cell | `src/core/cell_balances.py:27-64` | 双相残差与 `ex_bd` 符号正确 |
| Verbindungsmatrix | exchange-cell / global-nr | `src/core/connectivity_graph.py` | `zu/rez/auf/ab` 边类；VTT 单回路 `unexpected_blocks==0` |
| Eq. 2.9 NR 更新 | global-nr | `src/solvers/global_nr_solver.py` | 阻尼 Newton 更新式与 Wirsum (1998) 对应 |
| Nebenelemente 侧块 | global-nr | `src/solvers/structured_jacobian.py` | 侧块显式进入 Jacobian 与 Woodbury 校正求解 |
| Eq. 4.9 `Nu_P` | drying-pyro | `src/thermal/drying.py:71-89` | `2 + 1.2·Re^0.5·Pr^(1/3)` 完全按式实现 |
| Eq. 4.10–4.12 DAEM | drying-pyro | `src/thermal/devolatilization.py` | Gauss-Hermite 积分、径向体积分、Table 4.1 参数一致 |
| Eq. 5.36 R5 气泡 | kinetics | `src/kinetics/gas_reactions.py:82-106` | `k_standard` + 浓度幂次；audit 数值锚点 PASS |
| Eq. 5.38–5.39 R12 H₂ | kinetics | `src/kinetics/gas_reactions.py:261-279` | `k_jensen_h2`；audit 数值锚点 PASS |
| Eq. 5.40–5.41 R6 CH₄ | kinetics | `src/kinetics/gas_reactions.py:143-153` | `k_jensen_r7`；audit 数值锚点 PASS |
| Eq. 5.51 R7 重整 | kinetics | `src/kinetics/gas_reactions.py:160-179` | SI 换算 `k0=6.113e-7`、E=137327 J/mol 一致 |
| Eq. 5.59 R10 焦油氧化 | kinetics | `src/kinetics/tar_reactions.py:138-154` | `k_hobbs` + T·P^0.3·C_tar^0.5·C_O2；Pa 口径 |
| Eq. 5.60–5.61 R10 参数 | kinetics | `src/kinetics/tar_reactions.py:129-132` | Olefine/Aromaten 参数与论文一致 |
| Eq. A.27–A.28 Newton 步 | a1-gibbs | `src/thermodynamics/gibbs_minimizer.py`, `gibbs_hamel_reduced.py` | 两个 solver 均实现阻尼 Newton |

### 2.2 部分一致项（当前仍开放）

> **注**：§0 / §0.2 已闭合项（K_bd、R8、干燥 4.x、自由板 3.70–3.97、Jacobian、水力默认、固相/自由板 P1 默认）不再重复列出。

| 方程/表 | Scope | 代码位置 | 偏离说明 |
|---------|-------|----------|----------|
| Eq. (2-6) holdup guardrails | exchange-cell | `cell_balances.py` | 默认 `holdup_transport` 已启用；char/ash mask + NR guardrails 为 **B-tier accepted gap**（已文档化） |
| 粒径迁移 | exchange-cell | `calc_size_migration` | 默认 **`per_class_holdup`**；`total_bed_proxy` opt-in |
| Eq. (2-7) 热损 | exchange-cell | `calc_energy_balance_residual` | 默认 **`hamel_eq27_explicit`**；`inlet_fraction` opt-in |
| 侧块能量 | exchange-cell | side-block 温度贯通 | 非完整 Eq.(2-7) |
| Eq. 4.1–4.4 干燥 | drying-pyro | `drying.py` | Crank-Nicolson 数值替代；**`w_{0,tr}` 干基已按 page_36 复核** |
| Eq. 5.47 R8 平衡常数 | kinetics | `gas_reactions.py` | 修正 Keq，非 OCR 字面量 |
| A1 标注/判据/初值 | a1-gibbs | `gibbs_*.py` | **默认 `hamel_a40` + `holub_vonka_a34`（2026-07-08）**；Taylor/L2 opt-in |
| Vorab outer 刷新 | global-nr | `vorab_policy.py` | single-shot 源项已隔离；outer 水力刷新策略仍属 A 档调优 |
| 外层 Check2 | global-nr | `outer_loop.py` | 气相 split / ΔT 与 B 档 Jacobian 验收分离 |

### 2.3 已闭合项（7/4 初版曾标为缺失/部分一致）

| 方程/表 | 现状态 | 证据 |
|---------|--------|------|
| Eq. 3.46–3.49 替代 `K_bd` | **已实现（opt-in model）** | `mass_transfer.py` |
| Eq. 3.70–3.71 / 3.96 / 3.97 | **已实现** | `freeboard.py`, `trajectory_common.py` |
| Eq. 4.5 / 4.7 / 4.8 | **已实现** | `drying.py` |
| Eq. 5.45 / 5.49 R8 | **一致** | `gas_reactions.py` |
| Eq. 3.50 对流项 | **一致（PDF page_27）** | `(3/2)u_{b,r}/d_b` |
| 默认轨迹系数 | **exact_hamel（P1）** | `hamel_parity_defaults.py` |
| Bed 固相默认 | **holdup_transport（P1）** | `ReactorConfig.bed_solid_state_model` |

### 2.4 仍缺失 / 低优先级

| 方程/表 | Scope | 说明 |
|---------|-------|------|
| Eq. 4.8 独立 Bi 标量 | drying-pyro | 仅网格 Biot，无 standalone `Bi` 辅助函数 |
| K_bd 第三项（若论文有） | hydro-core | 主链 (3.50) 两项式已 PDF 锁定；第三项未在 page_27 等价式中出现 |

---

## 3. 与 Gap Checklist 的对比

当前 `docs/hamel_isomorphic_gap_implementation_checklist.md` §0.5–§0.6 已与本次审计结论**同步（2026-07-08）**。历史 diff 如下：

| 条目 | 7/4 gap §0.5（旧） | 现判定（§0.5/§0.6） |
|------|---------------------|---------------------|
| R10 焦油氧化 | 不一致 | **一致（P1 已闭合）** |
| 默认水力学配置 | 部分一致（Cell legacy） | **一致（Hamel 主链默认）** |
| Jacobian 侧块 | 部分一致（FD） | **P0 结构闭合 + P2 hybrid 生产默认** |
| Bed 固相 / 迁移 | 部分一致 | **holdup_transport + per_class_holdup 默认** |
| 能量 Eq.(2-7) | 部分一致 | **hamel_eq27_explicit 默认** |
| 自由板 drag | 部分一致 | **exact_hamel 默认** |
| A1 Gibbs | 部分一致（opt-in Hamel 分支） | **hamel_a40 + holub_vonka_a34 默认** |
| 干燥 w 干湿基 | 待 PDF 复核 | **page_36 已复核** |

---

## 4. 建议下一步动作

### P0（文档）

1. ~~**更新 audit summary §2.2/2.3**~~ ✅ 2026-07-07 已同步 §0.2 / §2.2–2.4
2. ~~**更新 gap checklist §0.5**~~ ✅ 2026-07-08 已同步 §0.5–§0.6 + 生产默认表
3. **统一方程编号引用**：`02_cell_balances` docstring 与 Hamel `(2-4)…(2-8)` 对齐

### P1（中优先级）

4. ~~**固相默认 holdup_transport**~~ ✅ `hamel_parity_defaults.py` + `ReactorConfig.bed_solid_state_model`
5. ~~**自由板默认 exact_hamel**~~ ✅ `ReactorConfig` + `simulate_freeboard` 默认
6. ~~**Eq.(2-7) 热损模式**~~ ✅ 默认 `hamel_eq27_explicit`（2026-07-07）；出口 T validation **待反应/能量校正**
7. ~~**干燥 `w` 干湿基**~~ ✅ `initial_moisture_dry_basis_from_wt_percent()` + page_36 锚点测试

### P2（低优先级 / 研究项）

8. ~~**A1 同构数值回归**~~ ✅ 默认 `hamel_a40` + `holub_vonka_a34`（2026-07-08）
9. ~~**解析 Jacobian / 生产默认 hybrid**~~ ✅ 2026-07-07
10. **Vorabrechnung strict single-shot 默认**：A 档 outer 策略仍待标定
11. **Jacobian walltime benchmark**：`data/jacobian_assembly_benchmark_latest.json` 已落盘

---

## 5. 完成标准自检

- [x] 全部 7 个 scope 的方程清单均已逐条审计并填入表格（含 `缺失` 项）。
- [x] 所有 `部分一致` 均有文件:行或分册段落证据。
- [x] 未用代码反向修改分册；文档变更建议单独列出并标注需 PDF 依据。
- [x] 自动化锚点 `audit_hamel_consistency.py` 与 `tests/sanity_checks.py` 均通过。
- [ ] 用户确认 scope 内无遗漏方程（或明确记录「本次跳过」项及原因）—— 待用户审阅。

---

## 6. 附录：各 Scope 详细报告路径

| Scope | 报告文件 |
|-------|----------|
| `hydro-core` | 由 agent `parity-hydro-core` 产出（见会话记录） |
| `freeboard` | 由 agent `parity-freeboard` 产出（见会话记录） |
| `exchange-cell` | `docs/audit_reports/hamel_parity_exchange_cell_2026-07-04.md` |
| `global-nr` | `docs/audit_reports/hamel_parity_global_nr_2026-07-06.md` |
| `drying-pyro` | 由 agent `parity-drying-pyro` 产出（见会话记录） |
| `kinetics` | 由 agent `parity-kinetics` 产出（见会话记录） |
| `a1-gibbs` | 由 agent `parity-a1-gibbs` 产出（见会话记录） |

---

*本报告由 `bfb-hamel-parity-audit` skill 执行生成，遵循证据链：论文 → 分册 → checked_spec → 代码 → gap checklist。*
