# Scripts Index (Cleaned)

本索引用于区分当前推荐脚本与历史脚本，避免重复审计链并降低维护成本。

## Policy

- `scripts/*.py` 根目录只保留当前主链推荐入口与高频审计脚本。
- 历史调试/重复脚本统一迁移到 `scripts/_deprecated/`。
- 新增脚本前应先检查是否可复用现有脚本，避免功能重复。

## Active Scripts (Recommended)

- `_audit_inputs.py`
- `_nr_monitor.py`
- `analyze_k10_units.py`
- `audit_bed_exit_to_reactor_exit_budget_lu.py`
- `audit_cell_conservation.py` (legacy-kept)
- `audit_freeboard_direction_lu.py`
- `audit_freeboard_entrained_solids_lu.py`
- `audit_freeboard_reaction_sets_lu.py`
- `audit_global_nr_profile_lu.py` (legacy-kept)
- `audit_gs_history_lu.py`
- `audit_hydrodynamics_consistency_lu.py`
- `audit_phase1_htw_lu.py`
- `audit_phase_exchange_conservation.py` (legacy-kept)
- `audit_r10_pressure_correction.py`
- `audit_r10_tar_oxidation.py`
- `audit_single_cell_oxygen_convergence.py`
- `audit_solver_parity_lu.py` (legacy-kept)
- `audit_zone_profile_lu.py`
- `benchmark_nr_init_gibbs_bootstrap.py`
- `run_phase_gate_00.py`
- `run_phase_gate_10.py`
- `run_phase_gate_20.py`
- `run_phase_gate_30.py`
- `run_phase_gate_40.py`
- `run_module_audits.py`
- `run_test_sequence.py`
- `validation_performance_report.py`
- `verify_r10_correction.py`

## Deprecated Scripts (Quarantine)

已迁移脚本请从以下分组查找，不再作为默认执行入口：

- `scripts/_deprecated/backsolve/`
- `scripts/_deprecated/freeboard/`
- `scripts/_deprecated/hydrodynamics/`
- `scripts/_deprecated/nr/`
- `scripts/_deprecated/cell/`
- `scripts/_deprecated/thermal/`
- `scripts/_deprecated/calibration/`
- `scripts/_deprecated/misc/`

## Migration Notes

- 本次迁移按 evidence-based 分级执行：文档绑定 + 主链必要 + 回归价值。
- 迁移目标是减少并行脚本链，统一审计口径到主链 `src/core` + `src/solvers` + `src/workflow`。
