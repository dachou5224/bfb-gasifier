# Equation Map（代码 ↔ 论文）

> 与 [`docs/gasifier_model_flowchart_trilingual.mmd`](../../../../docs/gasifier_model_flowchart_trilingual.mmd) 及
> `src.workflow.simulation_runner.describe_hamel_flowchart_mapping()` 一致：**主求解实现**在 `src/solvers/`，
> `src/gasifier/solvers/` 仅为同名符号重导出（见 `src/gasifier/solvers/__init__.py`）。

## 流程与全局 NR（Check1 / Check2）

- **Vorabrechnung / 初值**：`src/solvers/vorabrechnung.py`（`estimate_axial_T_profile`、`generate_initial_x0`、`refresh_vorabrechnung_for_cells`）
- **全域 Newton–Raphson**：`src/solvers/global_nr_solver.py`（`solve_global_nr`）
- **外层 Abgleich（预计算 ↔ 单元模型）**：`src/solvers/outer_loop.py`（`run_global_nr_outer_abgleich`）
- **结果汇总**：`src/solvers/result_builder.py`（`finalize_global_nr_result`）
- **逐格 fsolve（legacy）**：`src/solvers/cell_solver.py`（`solve_cell`）

## Hydrodynamics

- `src/gasifier/models/hydrodynamics.py`（薄封装）→ 实际计算：`src/core/cell_hydrodynamics.py`、`src/physics/*`
- `src/gasifier/models/correlations/drag_models.py`
  - Source: specs/02_hydrodynamics.md（待细化到方程号）

## Heat Transfer

- `src/gasifier/models/heat_transfer.py`（薄封装）→ 关联式：`src/gasifier/models/correlations/heat_transfer_coeff.py`；主链热损失分配见 `src/core/reactor.py`

## Reaction Model

- `src/gasifier/models/reaction_model.py`
  - Source: specs/04_kinetics.md（R1–R11）
- `src/gasifier/models/correlations/reaction_kinetics.py`
  - Source: 速率倍率护栏 / 与 `src/kinetics` 对齐
