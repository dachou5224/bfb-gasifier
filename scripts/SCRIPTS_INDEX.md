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
- `audit_bed_co_co2_budget_lu.py`（CO/CO₂ 床层来源–去向审计；SYMPTOM-MAP 主链）
- `audit_phase2_coco2_structural_rca.py`（出口 CO/CO₂ 结构偏差多因子 RCA：床内+自由板+裁定）
- `audit_bed_exit_to_reactor_exit_budget_lu.py`
- `audit_bed01_o2_inlet_split_lu.py`（bed0-bed1 O₂ 入口分配审计）
- `audit_split_plateau_cell_species.py`（平台期 cell×物种 phase-split 主导项；handoff P0）
- `audit_upper_bed_co_h2_phase_source.py`（上段床 CO/H2 相源定位与 split 项分解；handoff §5.4）
- `audit_upper_dense_syngas_absorb_ab.py`（上段密相 CO/H2 总量吸收 A/B；handoff §5.5）
- `audit_upper_dense_co_newton_block.py`（上段 dense-CO Newton/clip/fastox 阻塞；handoff §5.6）
- `audit_solid_ftb_relative_floor_ab.py`（固相 FTB 相对地板 A/B；handoff §5.6）
- `audit_upper_co_local2x2_ab.py`（上段 CO 局部 2×2 Pre-Inner A/B；handoff §5.7–5.8）
- `audit_upper_co_local2x2_reject.py`（local2x2 拒因分类 + CO 预算；handoff §5.8）
- `audit_preinner_split_lsq_detonation.py`（Pre-Inner split LSQ 爆炸 bed×物种；handoff §5.9）
- `audit_co_gap_to_platform.py`（逐 bed CO 预算 vs 平台 y≈0.154；handoff §5.11）
- `audit_local2_ps_abs_rms_stack.py`（local2+phase_share rms/energy 栈与 top 残差；handoff §5.12）
- `audit_bed6_local2_to_bed2_coupling.py`（local2→bed2 耦合 H1/H2；handoff §5.13）
- `audit_freeze_upper_bridge_ab.py`（同一次 solve 内冻结 upper CO 桥 A/B；handoff §5.13）
- `audit_bed2_flip_newton_window.py`（bed2 翻转 outer 窗口：refresh/preinner/inner/postinner；§5.13）
- `audit_bed2_flip_inner2_dx_mask.py`（inner#2 Δx 分组掩码 → bed2 O₂ 灰尘触发；§5.13）
  - 相关：`data/phase2_bed2_flip_o2_epsilon_sweep.json`、`data/phase2_snap_trace_o2_ab_o8.json`、`data/phase2_cont_bridge_on_off_heal.json`、`data/phase2_freeze_heal_bed2_o20.json`、`data/phase2_snap_o2_twostage_skip_precalc.json`
- `audit_outer_hot_restart_ab.py`（CO 带宽后嵌套 skip_precalc 热重启 A/B；§5.14）
  - 相关：`data/phase2_outer_hot_restart_ab_o8.json`、`data/phase2_outer_hot_restart_nested_fbproj_o8.json`
- `audit_trace_o2_rate_gate_ab.py`（微量 O₂ 速率门控：立即 vs CO 带宽后武装；§5.15）
  - 相关：`data/phase2_trace_o2_rate_gate_ab_o8.json`、`data/phase2_trace_o2_rate_gate_posthoc.json`
- `audit_rate_gate_rms_stack_ab.py`（延迟门控 + o16 / 热重启压 rms；§5.16）
  - 相关：`data/phase2_rate_gate_rms_stack_ab.json`
- `audit_gate_hot_topk_residuals.py`（gate/hot 末态 top 残差；§5.17）
  - 相关：`data/phase2_gate_hot_topk_residuals_o8.json`
- `audit_bed9_energy_closure_ab.py`（bed9/bed0 能量投影 A/B；§5.18）
  - 相关：`data/phase2_bed9_energy_closure_ab_o8.json`、`data/phase2_bed9_energy_terms.json`、`data/phase2_bed9_T_sweep_bc.json`
- `audit_bed9_fb_solid_enthalpy.py`（bed9↔fb 质量/焓连续 + 床顶焓透传 A/B；§5.19）
  - 相关：`data/phase2_bed9_fb_solid_enthalpy_ab_o8.json`
- `audit_bed9_dm_deferred_passthrough.py`（bed9 固相 dm 分解 + 延迟焓透传 A/B；§5.20）
  - 相关：`data/phase2_bed9_dm_deferred_passthrough_o8.json`
- `audit_bed9_holdup_closure_ab.py`（bed9 char holdup / support / 床顶 solid LSQ A/B；§5.21）
  - 相关：`data/phase2_bed9_holdup_closure_ab_o8.json`
- `audit_bed9_negR_support_combo_ab.py`（负 R support vs LSQ/焓透传组合；§5.22）
  - 相关：`data/phase2_bed9_negR_support_combo_ab_o8.json`
- `audit_bed8_after_top_lsq_ab.py`（LSQ 后 bed8 热点 + n_beds/split A/B；§5.23）
  - 相关：`data/phase2_bed8_after_top_lsq_ab_o8.json`
- `audit_bed3_n2_on_lsq_n2_ab.py`（LSQ n=2 上 bed3 N2 + targeted split A/B；§5.24）
  - 相关：`data/phase2_bed3_n2_on_lsq_n2_ab_o8.json`
- `audit_bed3_n2_sum_closure_ab.py`（bed3 N2 总量过填 + 双向放出 A/B；§5.25）
  - 相关：`data/phase2_bed3_n2_sum_closure_ab_o8.json`
- `audit_n2_axial_cascade.py`（全床 N2 轴向 cascade 通盘；§5.26）
  - 相关：`data/phase2_n2_axial_cascade_o8.json`
- N2 inert inflow 闭合验证（§5.27；inline / 复跑 `_lu_harness`）
  - 相关：`data/phase2_n2_inert_inflow_fix_o8.json`
- 后 N2 残差重扫 / post-staged H2O 栈确认（§5.28）
  - 相关：`data/phase2_post_n2_rms_rescan_ab.json`、`data/phase2_post_n2_h2o_stack_confirm_o8.json`
- bed0 能量 / 双吸引子（§5.29；审计请 `PYTHONHASHSEED=0`）
  - 相关：`data/phase2_attractor_stability_o8.json`、`data/phase2_n2only_lsq2_stability_o8.json`、`data/phase2_bed0_energy_seed_epre_omp1_ab_o8.json`
- bed0 T Newton/LS 回退（§5.30）
  - 相关：`data/phase2_bed0_T_nr_trace_o8.json`、`data/phase2_bed0_T_newton_ab_o8.json`、`data/phase2_bed0_T_ls_reject_probe.json`
- bed0 T LS merit 分解 / freeze-gas 投影（§5.31；inline + `PYTHONHASHSEED=0`）
  - 相关：`data/phase2_bed0_T_ls_merit_decomp_o8.json`、`data/phase2_bed0_T_ls_projection_ab.json`、`data/phase2_bed0_T_split_constrained_ab_o8.json`
- bed0 composition-frozen T / 路径依赖（§5.32；inline）
  - 相关：`data/phase2_bed0_T_comp_frozen_preproj_ab_o8.json`、`data/phase2_bed0_T_heat_override_ab_o8.json`、`data/phase2_bed0_T_frozen_vs_coupled_dT_o4.json`、`data/phase2_bed0_T_local_slope_o8.json`、`data/phase2_bed0_T_clean_scan_o8.json`
- bed0 深快照 F_E(T)（§5.33；inline）
  - 相关：`data/phase2_bed0_T_hidden_state_diff_o4.json`、`data/phase2_bed0_T_deep_snap_scan_o8.json`
- bed0 气–能联合软桥 A/B（§5.34；inline）
  - 相关：`data/phase2_bed0_T_gas_energy_joint_ab_o8.json`
- 平台 rms 景观 / 压 rms A/B（§5.35；inline）
  - 相关：`data/phase2_rms_landscape_current_stack_o8.json`、`data/phase2_rms_press_ab_o8.json`
- bed3 H2O 再膨胀 / skip_precalc 对照（§5.36；inline）
  - 相关：`data/phase2_bed3_h2o_reinflate_trace_o8.json`、`data/phase2_bed3_h2o_reinflate_skip_ab.json`
- 能量项分解 / 焓桥 A/B（§5.37；inline）
  - 相关：`data/phase2_energy_term_decomp_platform_o8.json`、`data/phase2_energy_press_ab_o8.json`
- 能量 H_in 分项 / ΣE / dm×hs 反事实（§5.38；inline）
  - 相关：`data/phase2_energy_hin_parts_axial_o8.json`
- 平台全炉逐 cell 剖面 + 异常（§5.39；inline）
  - 相关：`data/phase2_platform_cell_profile_o8.json`
- 逐格 H₂ 源汇预算（§5.40；inline）
  - 相关：`data/phase2_platform_h2_cell_budget_o8.json`
- H₂ DOF / 强制闭合 / fastox 探针（§5.41；inline）
  - 相关：`data/phase2_platform_h2_nh2_lock_probe_o8.json`
- `audit_bed_phase_split_contribution.py`（逐 bed split 方差贡献 / T2 探针）
  - 相关证据：`data/phase2_bed0_lsq_h2o_twist_o8.json`、`data/phase2_bed0_phase_only_ab_o8.json`（§5.10）；`data/phase2_phase_share_absorb_gate_ab_o8.json`（§5.11）；`data/phase2_bed2_h2o_hotspot_vs_local2_o8.json`、`data/phase2_local2_minbed_ps_abs_ab_o12.json`（§5.12）；`data/phase2_bed6_local2_to_bed2_coupling_o4.json`、`data/phase2_bed2_flip_traj_twostage_o12.json`、`data/phase2_freeze_upper_bridge_ab_o8.json`（§5.13）
- `audit_cell_conservation.py` (legacy-kept)
- `audit_char_mass_conservation_lu.py`（被 phase1 audit / tests import）
- `audit_flowchart_mapping.py`（`run_module_audits` workflow 模块）
- `audit_freeboard_direction_lu.py`
- `audit_freeboard_entrained_solids_lu.py`
- `audit_freeboard_reaction_sets_lu.py`
- `audit_global_nr_profile_lu.py` (legacy-kept)
- `audit_homogeneous_bubble_reduction_lu.py`（均相/气泡 R 分工审计）
- `audit_hamel_consistency.py`（Hamel 锚点一致性）
- `audit_energy_heat_loss_mode_compare.py`（Phase2：`inlet_fraction` vs `hamel_eq27_explicit` 热损对比）
- `audit_hydrodynamics_consistency_lu.py`
- `audit_init_residual_breakdown.py`（init 残差 / x₀ 不可行诊断）
- `audit_init_staged_w1_w2.py`（staged init w1/w2 对比审计）
- `audit_phase1_htw_lu.py`
- `audit_phase_exchange_conservation.py` (legacy-kept)
- `audit_r10_pressure_correction.py`
- `audit_r10_tar_oxidation.py`
- `audit_single_cell_oxygen_convergence.py`
- `audit_zone_profile_lu.py`
- `benchmark_nr_init_gibbs_bootstrap.py`
- `freeze_wirsum_alignment_baseline.py`（Wirsum 对齐 P0：冻结 reduced/LU 残差与 Jacobian 指标）
- `benchmark_jacobian_assembly.py`（structured / hybrid / dense 组装 walltime + parity）
- `benchmark_opt012_walltime_compare.py`（opt012 walltime 对比；optimization docs）
- `benchmark_phase2_convergence_ladder.py`（P2 收敛 ladder SSOT；`_lu_harness.LADDER_LEVELS`）
- `benchmark_phase2_multi_case_matrix.py`（多 case 矩阵趋势）
- `benchmark_preinner_perf_experiments.py`（preinner 性能实验；可选前后对比）
- `build_validation_cases_master_table.py`（生成 15 Case 主表与 JSON）
- `calibrate_phase2_co_co2_window.py`（diagnostic-only `r*_scale` 扫描）
- `compare_low_nr_init_temperature_lu.py`（低 NR init 温度对比）
- `compare_vorab_march_lu.py`（Vorab march 对比；`--scenario staged|weights|both`）
- `enrich_validation_cases_from_hamel.py`（从 Hamel 文献 enrich validation cases）
- `merge_fig7_parity_outputs.py`（合并 Fig7 parity 输出）
- `render_hamel_fig7_parity_pages.py`（渲染 parity 底图）
- `run_phase_gate.py`（统一 gate runner；替代已删除的 `run_phase_gate_00–40.py` 薄包装）
- `sensitivity_ud_closure_o2_co_co2.py`（UD closure O₂/CO/CO₂ 灵敏度）
- `run_module_audits.py`
- `run_test_sequence.py`
- `validate_slimming_plan.py`
- `validation_performance_report.py`
- `vtt_full_sim_residual_audit.py`（VTT 全模拟残差审计）
- `verify_r10_correction.py`

## Infrastructure（非入口，勿删）

- `_lu_harness.py` — Phase2 LU 公共 harness（`LADDER_LEVELS` SSOT）
- `_audit_common.py` — audit 脚本共享纯函数（Tier M 合并产物）

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
- 瘦身流程 skill：`.cursor/skills/bfb-scripts-hygiene/`（含 `SCRIPT-TIERS.md` 分级台账）
- 历史 GS 求解路径已自 `src/core/reactor.py` 移除；相关脚本见 `_deprecated/nr/`
- **2026-07-04 Tier C 清理**：删除薄包装 `run_phase_gate_00–40.py`，统一使用 `python3 scripts/run_phase_gate.py <00|10|20|30|40|all>`；删除 `compare_staged_march_weights_lu.py` / `compare_vorab_march_staged_lu.py`，统一使用 `python3 scripts/compare_vorab_march_lu.py --scenario staged|weights|both`
- **2026-07-04 Tier B 归档**：将 9 个零引用实验脚本迁入 `scripts/_deprecated/{cell|nr|thermal|calibration}/`，详见 `scripts/_deprecated/README.md`
