# OPT-002 → OPT-001 实施计划与 QA 门禁

> 基于 [`codebase_optimization_index.md`](codebase_optimization_index.md) 推荐顺序；OPT-007/008/012 已完成。

## 总览

| 阶段 | 任务 | 交付物 | QA Gate |
|------|------|--------|---------|
| **A** | OPT-002 Phase 1 | `preprojection_engine.py` + solid holdup 迁移 | **QA-A** ✅ |
| **B** | OPT-002 Phase 2 | phase-split + T-local 迁移 | **QA-B** ✅ |
| **C** | OPT-002 Phase 3 | joint relax 模板 + gas/T LSQ 迁移 | **QA-C** ✅ |
| **D** | OPT-001 | `connectivity/` 分包 + re-export | **QA-D** ✅ |
| **E** | OPT-013 | `init_precalc_step` 阶段函数化 | **QA-E** ✅ |
| **F** | OPT-005 | `global_nr_iteration.py` 主循环提取 | **QA-F** ✅ |
| **G** | OPT-006 | `result_exit_summary.py` 分段汇总 | **QA-G** ✅ |
| **H** | OPT-004 | `vorabrechnung/` 包拆分 | **QA-H** ✅ |

后续（本计划外、按需并行）：OPT-009。

---

## 阶段 A — OPT-002 Phase 1（当前）

**目标：** 提取 LSQ + improvement gate + restore 骨架；迁移最简函数 `preproject_bed_solid_holdup_to_transport_closure`。

**新建模块：** `src/core/preprojection_engine.py`

| 组件 | 职责 |
|------|------|
| `preproj_cell_residual_vector` / `preproj_cell_energy_residual_abs` | 从 connectivity 迁入（OPT-012 partial residual） |
| `PreprojectionBcContext` | trial `refresh` / accept-reject `finalize` |
| `LsqPreprojectionSpec` | bounds、max_nfev、tol |
| `run_improvement_gated_lsq(...)` | try/LSQ/accept/restore/sync 统一路径 |

**connectivity 兼容：** 模块级 `_preproj_*` 别名保留，benchmark legacy shim 无需改动。

### QA-A（必须通过才能进 Phase 2）

```bash
pytest tests/test_preprojection_engine.py tests/test_bottom_gas_inlet_split.py -q
python3 tests/sanity_checks.py
```

**通过准则：**

- [ ] engine 单元测试：mock LSQ accept/reject/exception 三条路径
- [ ] `test_bottom_gas_inlet_split` 全绿（含 partial residual slice 测试）
- [ ] sanity_checks 17/17
- [ ] solid holdup 行为不变：diag 字段、`improvement_factor` 门控、skip_rms 逻辑

---

## 阶段 B — OPT-002 Phase 2

**迁移：**

1. `preproject_bed_major_gas_phase_split_to_exchange_closure`
2. `preproject_bed_temperature_to_local_energy_closure`

**风险：** primary inlet totals 分支、T bounds 与 BC refresh 顺序。

### QA-B

```bash
pytest tests/test_bottom_gas_inlet_split.py::test_init_precalc_records_guarded_joint_preprojection_relax \
  tests/test_bottom_gas_inlet_split.py::test_partial_cell_residual_groups_match_full_vector_slices -q
python3 scripts/benchmark_opt012_walltime_compare.py --repeats 2
```

**通过准则：**

- [ ] joint relax 诊断字段数值与 Phase A baseline JSON diff 在 rtol=0（逻辑路径未改）
- [ ] preinner_bridge median 不劣于 OPT-012 基线（±10% 噪声带）

---

## 阶段 C — OPT-002 Phase 3

**迁移：**

1. `preproject_bed_gas_temperature_to_energy_closure`（大函数，neighbor couple 保留局部逻辑）
2. `run_relax_schedule(...)` — joint / bc-aware alpha 循环
3. `run_bottom_joint_start_value_preprojection` 改用 relax 引擎

### QA-C

```bash
pytest tests/test_bottom_gas_inlet_split.py tests/test_vorabrechnung_a_tier_budget.py -q
python3 scripts/run_phase_gate.py 20 --json
```

**通过准则：**

- [ ] gate 20 structural parity 通过
- [ ] vorab A-tier budget 测试无回归
- [ ] connectivity.py 预投影区行数减少 ≥300 行（相对 Phase A 前）

---

## 阶段 D — OPT-001 connectivity 分包

**目录：**

```
src/core/connectivity/
  __init__.py          # 对外 API re-export（import 路径不变）
  preprojection.py
  solid_transport.py
  gas_inlet.py
  nr_boundary.py
  joint_gate.py
```

`preprojection_engine.py` 暂留 `src/core/`，Phase D 移入 `connectivity/preprojection_engine.py`。

### QA-D

```bash
pytest tests/test_connectivity_graph.py tests/test_bottom_gas_inlet_split.py tests/test_global_nr_solver.py -q
npx gitnexus@latest impact _bottom_joint_gate_merit   # 人工确认 fan-in 无断链
```

---

## 阶段 E — OPT-013 init 阶段函数化

依赖 OPT-002 引擎稳定后拆分 `run_init_and_precalc_for_global_nr`（381 行 → 5 阶段函数）。

### QA-E

```bash
pytest tests/test_vorabrechnung_*.py tests/test_bottom_gas_inlet_split.py -q
python3 scripts/run_phase_gate.py all --strict
```

---

## 阶段 F — OPT-005 NR 主循环提取

**交付：**

- `src/solvers/global_nr_iteration.py` — `execute_global_nr_solve`
- `global_nr_solver.solve_global_nr` — 薄 wrapper；迭代体经 `_gnr` 模块别名调用以兼容测试 monkeypatch

### QA-F

```bash
pytest tests/test_global_nr_solver.py -q -k "line_search or hydrodynamics or retries"
python3 tests/sanity_checks.py
pytest tests/test_flowchart_mapping_contract.py tests/test_outer_loop_convergence.py -q
```

---

## 阶段 H — OPT-004 vorabrechnung 包拆分

**交付：** `src/solvers/vorabrechnung/` 五子模块 + `__init__.py` 66 项 re-export

### QA-H

```bash
pytest tests/test_vorabrechnung_*.py tests/test_init_precalc_step_stages.py -q
python3 tests/sanity_checks.py
```

---

| 指标 | 命令 |
|------|------|
| 预投影 walltime | `python3 scripts/benchmark_opt012_walltime_compare.py --repeats 3` |
| 收敛 ladder | `python3 scripts/benchmark_phase2_convergence_ladder.py --levels p2 --outer-max 10` |
| 科学审计 | bfb-scientific-audit skill on 改动文件 |

---

## 回滚策略

- 每阶段独立 commit（用户请求时）
- engine 未接管的函数保持原实现；单函数回滚 = revert 该函数调用点
- QA 失败：停下一阶段，fix 当前阶段后再过 Gate
