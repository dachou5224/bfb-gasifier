---
name: bfb-numerical-guard
description: BFB 数值稳定性监控。针对 fsolve（cell_solver）与 global_nr（global_nr_solver），检查残差缩放、np.exp 保护、Jacobian 稀疏结构等。在修改求解器、残差或边界条件后启用。
---

# bfb-numerical-guard（数值稳定性监控）

## 何时启用

- 修改 **`src/solvers/cell_solver.py`**（逐格 `fsolve`）、**`src/solvers/global_nr_solver.py`**（全局阻尼 NR）、或 **`Cell.residuals`** / 边界条件闭包。
- 调整 **残差缩放**、**未知量打包顺序**、**收敛判据** 后。
- 出现 **Newton–Raphson 不收敛**、**残差范数爆炸**、或 `sanity_checks` 通过但端到端异常时（可与 `gstack-investigate` 配合）。

## 1. 残差缩放（Scaling）

### global_nr

- 阅读 `global_nr_solver.py` 中 **`build_equation_scales`** 与 **RMS ‖F̂‖** 的定义；缩放应使用 **按方程类型的物理参考量**（气相 mol/s、固相 kg/s、能量 W 等），避免仅用当前 ‖F‖ 动态缩放导致 NR「只看某一类方程」。
- 修改 `ref_gas_mol_s`、`ref_solid_kg_s`、`ref_energy_W` 的传入（见 `reactor._solve_global_nr`）时，确认与工况量级一致。

### fsolve / cell

- 若对残差做手工缩放或未知量归一化，需在注释中说明目的，并检查是否破坏 **雅可比条件数**（避免一边倒）。

## 2. `np.exp` 与刚性项

- 对 **Arrhenius、平衡常数、DAEM** 等路径中的 `np.exp(±x)`：
  - 检查 `x` 是否可能超出浮点安全范围；必要时使用 **`np.clip(x, -limit, limit)`** 或项目已有惯例（见 `.cursor/rules/scientific-modeling-core.mdc` 示例）。
- 禁止在无界温度/压力试验点上直接堆叠多个 `exp` 而无防护。

## 3. Jacobian 与稀疏结构（global_nr）

- **未知量顺序** 必须与 `n_var` / `pack_cell` / `unpack_cell` / `cell_offsets` 一致。
- 有限差分 Jacobian 使用 **`scipy.sparse`** 组装时，确认：
  - 块结构是否覆盖 **格间耦合**（边界条件闭包后顶格↔底格是否进入残差）。
  - 非零模式是否与 `n_var` 维度匹配；修改变量数后是否需更新图案（pattern）。

## 4. 收敛判据与「假收敛」

- **`reactor._solve_global_nr`**：外层收敛必须依赖 **内层 NR 的 `converged`** 与温度变化阈值，不得仅凭 `dT` 过小而忽略巨大残差（见 `docs/validation_gap_analysis.md` 与近期 reactor 逻辑）。
- 报告结果时区分 **`converged_outer` / `converged_inner_nr`**（若返回 dict 中含这些键）。

## 5. 快速自检命令

```bash
cd bfb-gasifier && python3 tests/sanity_checks.py
```

针对求解器的专项测试（若存在）：

```bash
pytest tests/test_lu_gs_vs_global_nr.py -q
```

## 相关文件

- `src/solvers/cell_solver.py`
- `src/solvers/global_nr_solver.py`
- `src/core/reactor.py`（`_solve_global_nr`、`_solve_gauss_seidel`）
- `src/core/cell.py`（`residuals`）
- `.cursor/rules/scientific-modeling-core.mdc`（§数值稳定性）
