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
