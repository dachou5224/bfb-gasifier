# BFB 症状 → 反馈环映射

Agent 在 Phase 1 按用户描述的**主症状**选环；多条症状并存时，先建**最窄、最快**的那条。

## 符号

- **seam**：可独立调用的公共边界（单速率、单格、单模块 pytest）
- **tight**：秒级、确定性、断言具体物理量（非仅「未崩溃」）
- **slow**：分钟级；仅作确认，不作内层假设-验证循环

---

## A. 单元 / 方程级（优先，最快）

| 主症状 | 首选反馈环 | seam |
|--------|------------|------|
| 某 Hamel 方程/Table 与代码不符 | `python3 scripts/audit_hamel_consistency.py` | `arrhenius`、`*_reference` 模块 |
| u_mf / u_br / K_bd 水力链 | `pytest tests/test_bubble_dynamics.py tests/test_mass_transfer.py -q` + `audit_hydrodynamics_consistency_lu.py` | `calc_u_br`, `calc_kbd` |
| R5/R6/R8/R10 速率或指数 | `pytest tests/test_cell_kinetics.py tests/test_tar_reactions.py -q` | `rate_R5_*`, `rate_R12`, `get_K_eq` |
| R10 压力项 Pa/bar 混用 | `python3 scripts/audit_r10_pressure_correction.py` | R10 速率函数 |
| 单格 fsolve 不收敛 / infeasible | `python3 scripts/audit_single_cell_oxygen_convergence.py` | `solve_cell` |
| 单格守恒残差 | `python3 scripts/audit_cell_conservation.py` | `Cell.residuals` |

## B. 床层 / 自由板耦合

| 主症状 | 首选反馈环 | 备注 |
|--------|------------|------|
| 炭质量不守恒 / bed-top gap | `python3 scripts/audit_char_mass_conservation_lu.py --mode phase2 --max-global-iter 1` | 区分 stale closure vs 算法 |
| CO/CO₂ 轴向预算异常 | `python3 scripts/audit_bed_co_co2_budget_lu.py` | 输出 JSON 到 `data/` |
| O₂ 床层分配 / bed0-bed1 split | `python3 scripts/audit_bed01_o2_inlet_split_lu.py` | |
| 自由板方向 / 夹带固相 | `audit_freeboard_direction_lu.py`, `audit_freeboard_entrained_solids_lu.py` | |
| 自由板反应集影响 | `audit_freeboard_reaction_sets_lu.py` | 关反应对照 |
| 均相/气泡 R 分工 | `audit_homogeneous_bubble_reduction_lu.py` | |

## C. 全局 NR / 外层 Abgleich

| 主症状 | 首选反馈环 | 备注 |
|--------|------------|------|
| init 残差高 / x₀ 不可行 | `audit_init_residual_breakdown.py`, `audit_init_staged_w1_w2.py` | 读 `vorab_init_abgleich_vs_outer_preinner.md` |
| outer 不收敛 / rms_scaled 平台 | `benchmark_phase2_convergence_ladder.py --level p2 --max-global-iter 20` | **timeout 420s** |
| Check1 vs Check2 语义混淆 | 查 `docs/vorab_init_abgleich_vs_outer_preinner.md` §2–3 | 再跑 `audit_global_nr_profile_lu.py` |
| Jacobian / structured 差异 | `pytest tests/test_global_nr_solver.py -k 'band_plus or block_tridiag' -q` | 配合 `bfb-numerical-guard` |
| 假收敛（dT 小但残差大） | `audit_global_nr_profile_lu.py --strict` | 看 `converged_outer` vs inner |

## D. 验证 / 文献对标（慢，非 root-cause 首选）

| 主症状 | 首选反馈环 | 前提 |
|--------|------------|------|
| 出口 CO/CO₂ vs Table 2 LU | `audit_phase1_htw_lu.py` | **parity 已闭合** |
| Fig.7 parity | `pytest tests/test_fig7_parity_digitization.py -q` | 数字化 JSON 已更新 |
| 多 case 趋势 | `benchmark_phase2_multi_case_matrix.py` | 勿用 LU 默认外推 |

## E. 回归 / 对比

| 场景 | 反馈环 |
|------|--------|
| 某次提交后变差 | `git bisect run` + 上表中最窄 pytest/audit |
| 两 config 开关对比 | `benchmark_phase2_convergence_ladder.py` baseline vs p2 |
| 两次 JSON 实验 diff | diff `data/phase2_*_latest.json` 中 `rms_scaled_final`、物种出口 |

---

## 模块审计入口（不确定症状时）

```bash
cd bfb-gasifier && python3 scripts/run_module_audits.py --module physics
cd bfb-gasifier && python3 scripts/run_module_audits.py --list
```

各 module 含 `triage_hint` 指向 `src/` 子树。

---

## F. 脚本目录维护（非运行时诊断）

| 主症状 | 入口 |
|--------|------|
| scripts 膨胀、重复 audit、未登记脚本 | `.cursor/skills/bfb-scripts-hygiene/` + `scripts/SCRIPTS_INDEX.md` |
