# Scripts 分级台账

**状态：2026-07-04 Tier B/C/S 已落地**（Copilot + 收尾验收通过）

## 当前统计

| 指标 | 瘦身前 (2026-07-03) | 当前 |
|------|---------------------|------|
| 根目录 `scripts/*.py` | 61 | **46** |
| `SCRIPTS_INDEX` Active | 26 | **44** |
| 未登记脚本 | 34 | **0** |
| `_deprecated/` | 67 | **76** |

验收：

```bash
pytest tests/test_codebase_slimming_plan.py -m slimming_audit -q
python3 scripts/validate_slimming_plan.py --skip-pytest-full
```

---

## ✅ Tier C — 已完成（删除薄包装）

| 原脚本 | 替代 |
|--------|------|
| `run_phase_gate_00–40.py` | `python3 scripts/run_phase_gate.py 00\|10\|…\|all` |
| `compare_staged_march_weights_lu.py` | `compare_vorab_march_lu.py --scenario weights` |
| `compare_vorab_march_staged_lu.py` | `compare_vorab_march_lu.py --scenario staged` |

编码于 `tests/test_codebase_slimming_plan.py` → `REMOVED_SCRIPTS`。

---

## ✅ Tier B — 已完成（归档 `_deprecated/`）

| 原路径 | 新路径 |
|--------|--------|
| `audit_char_reduction_sensitivity_lu.py` | `_deprecated/cell/` |
| `audit_bed_bottom_co_co2_profile.py` | `_deprecated/cell/` |
| `audit_init_holdup_seed_variants_w2.py` | `_deprecated/nr/` |
| `benchmark_phase2_temperature_chain_lu.py` | `_deprecated/nr/` |
| `benchmark_phase2_vm_zone_iter_sweep.py` | `_deprecated/nr/` |
| `benchmark_phase2_vm_zone_r5_mitigation_sweep.py` | `_deprecated/nr/` |
| `audit_temperature_fence_secondary_lu.py` | `_deprecated/thermal/` |
| `calibrate_phase2_fixed50_r5_composition_grid.py` | `_deprecated/calibration/` |
| `sweep_phase2_gas_inlet_dense_frac_lu.py` | `_deprecated/calibration/` |

详见 `scripts/_deprecated/README.md` §2026-07-04。

---

## ✅ Tier S/A — 已完成（登记 `SCRIPTS_INDEX.md`）

Gate 与文档工具脚本已全部登记，含：

- `benchmark_phase2_convergence_ladder.py`、`benchmark_phase2_multi_case_matrix.py`
- `audit_bed_co_co2_budget_lu.py`、`audit_init_*`、`run_phase_gate.py`
- `build_validation_cases_master_table.py`、`merge_fig7_parity_outputs.py` 等

Infrastructure（不列入 Active 入口列表）：

- `_lu_harness.py`、`_audit_common.py` — 见 `SCRIPTS_INDEX.md` §Infrastructure

---

## 🔄 Tier M — 进行中 / 可选后续

| 簇 | 状态 | 下一步 |
|----|------|--------|
| bed CO/CO₂ audit | 部分 | `_audit_common.py` 已抽取 `prepare_cell_for_analysis`；可继续合并 bottom profile 逻辑 |
| Init 诊断三件套 | 未合并 | `audit_init_*` → 单 CLI `--mode staged\|holdup\|residual` |
| Vorab compare | ✅ | 已合并为 `compare_vorab_march_lu.py` |
| Phase2 benchmark sweep | 部分 | sweep 类已 deprecated；主链 ladder + multi_case 保留 |

**原则**：Tier M 每次只合并一簇；每步跑 `validate_slimming_plan.py`。

---

## 禁止回退

- ❌ 把 `_deprecated/` 脚本移回根目录而不更新 slimming 测试
- ❌ 恢复 `run_phase_gate_XX.py` 薄包装
- ❌ 新增根目录脚本而不更新 `SCRIPTS_INDEX.md`（`test_all_root_scripts_indexed_or_infrastructure` 会 FAIL）
