# Deprecated Scripts

本目录存放**已退出主链**的调试、对照与一次性审计脚本。默认 CI / `run_test_sequence.py` / `run_module_audits.py` **不会**调用它们。

## 使用规则

- 需要复盘历史推理时可查阅，但不要作为当前实现依据。
- 多数脚本的 docstring 仍引用旧路径 `scripts/*.py`；运行前请改用 `_deprecated/...` 完整路径。
- `audit_sequential_stack.py` 等组合脚本可能包含已失效的子进程路径。

## 近期迁入（瘦身计划 2026-06）

| 原路径 | 新路径 | 原因 |
|--------|--------|------|
| `audit_gs_history_lu.py` | `_deprecated/nr/` | NR-only 后 `solver=gauss_seidel` 不可用 |
| `audit_solver_parity_lu.py` | `_deprecated/nr/` | GS vs NR 对照已无意义 |
| `benchmark_structured_jacobian.py` | `_deprecated/nr/` | 一次性 Jacobian benchmark |
| `audit_bottom_bed_energy_balance_lu.py` | `_deprecated/thermal/` | 孤立 Phase2 能量诊断 |
| `audit_phase2_freeboard_coeff_compare.py` | `_deprecated/freeboard/` | 零 gate 引用 sensitivity |
| `audit_psi_b_strategy_initialized_lu.py` | `_deprecated/hydrodynamics/` | 与 heinbockel 审计重复 |

## 2026-07-04 Tier B 归档

| 原路径 | 新路径 | 原因 |
|--------|--------|------|
| `audit_char_reduction_sensitivity_lu.py` | `_deprecated/cell/` | 一次性 char 还原灵敏度实验 |
| `audit_bed_bottom_co_co2_profile.py` | `_deprecated/cell/` | 零引用；可合并进 `audit_bed_co_co2_budget_lu.py` |
| `audit_init_holdup_seed_variants_w2.py` | `_deprecated/nr/` | 一次性 init holdup seed 变体实验 |
| `benchmark_phase2_temperature_chain_lu.py` | `_deprecated/nr/` | 一次性 temperature chain sweep |
| `benchmark_phase2_vm_zone_iter_sweep.py` | `_deprecated/nr/` | 一次性 VM zone 迭代 sweep |
| `benchmark_phase2_vm_zone_r5_mitigation_sweep.py` | `_deprecated/nr/` | 一次性 R5 mitigation sweep |
| `audit_temperature_fence_secondary_lu.py` | `_deprecated/thermal/` | 一次性 secondary injection 温度围栏实验 |
| `calibrate_phase2_fixed50_r5_composition_grid.py` | `_deprecated/calibration/` | 一次性 fixed50 R5 组成网格标定 |
| `sweep_phase2_gas_inlet_dense_frac_lu.py` | `_deprecated/calibration/` | 一次性气相入口 dense fraction sweep |
