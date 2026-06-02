# Scripts Index (Cleaned)

本索引用于区分当前推荐脚本与历史脚本，避免重复审计链并降低维护成本。

## Policy

- `scripts/*.py` 根目录只保留当前主链推荐入口与高频审计脚本。
- 历史调试/重复脚本统一迁移到 `scripts/_deprecated/`（见 `_deprecated/README.md`）。
- 新增脚本前应先检查是否可复用现有脚本，避免功能重复。

## Active Scripts (Recommended)

- `_audit_inputs.py`
- `_nr_monitor.py`
- `analyze_k10_units.py`
- `audit_bed_exit_to_reactor_exit_budget_lu.py`
- `audit_cell_conservation.py` (legacy-kept)
- `audit_char_mass_conservation_lu.py`（被 phase1 audit / tests import）
- `audit_flowchart_mapping.py`（`run_module_audits` workflow 模块）
- `audit_freeboard_direction_lu.py`
- `audit_freeboard_entrained_solids_lu.py`
- `audit_freeboard_reaction_sets_lu.py`
- `audit_global_nr_profile_lu.py` (legacy-kept)
- `audit_hamel_consistency.py`（Hamel 锚点一致性）
- `audit_hydrodynamics_consistency_lu.py`
- `audit_phase1_htw_lu.py`
- `audit_phase_exchange_conservation.py` (legacy-kept)
- `audit_r10_pressure_correction.py`
- `audit_r10_tar_oxidation.py`
- `audit_single_cell_oxygen_convergence.py`
- `audit_zone_profile_lu.py`
- `benchmark_nr_init_gibbs_bootstrap.py`
- `run_phase_gate_00.py` … `run_phase_gate_40.py`
- `run_module_audits.py`
- `run_test_sequence.py`
- `validate_slimming_plan.py`
- `validation_performance_report.py`
- `verify_r10_correction.py`

## Deprecated Scripts (Quarantine)

已迁移脚本请从以下分组查找，不再作为默认执行入口：

- `scripts/_deprecated/backsolve/`
- `scripts/_deprecated/freeboard/`（含 `audit_phase2_freeboard_coeff_compare.py`）
- `scripts/_deprecated/hydrodynamics/`（含 `audit_psi_b_strategy_initialized_lu.py`）
- `scripts/_deprecated/nr/`（含 GS 历史审计、`benchmark_structured_jacobian.py`）
- `scripts/_deprecated/cell/`
- `scripts/_deprecated/thermal/`（含 `audit_bottom_bed_energy_balance_lu.py`）
- `scripts/_deprecated/calibration/`
- `scripts/_deprecated/misc/`

## Migration Notes

- 主链 gate 编排：`run_test_sequence.py`、`run_module_audits.py`
- 瘦身假设回归：`tests/test_codebase_slimming_plan.py`、`scripts/validate_slimming_plan.py`
- 历史 GS 求解路径已自 `src/core/reactor.py` 移除；相关脚本见 `_deprecated/nr/`
