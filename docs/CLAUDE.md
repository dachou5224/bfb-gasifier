# BFB 气化炉 — Agent 开发约定

算法规格与 Hamel 方程转录见 [`docs/README.md`](README.md)（`hamel_submodels/00–07` 为方程 SoT）。本文档为**实现约定**与 **Phase 顺序**。

## 全局编程规则

### 单位制

- 全程 **SI**：Pa, K, m, mol, kg, s
- 唯一例外：R8 WGSR 方程内部将 Pa 换算为 atm，**外部接口仍传 Pa**

### Arrhenius 形式

- **只能**调用 `src/kinetics/arrhenius.py` 工厂函数
- R1–R4：`k_hobbs`；标准形式：`k_standard`；R7：`k_jensen_r7`
- ❌ 禁止在其他文件手写 `k = A * exp(-E/RT)`

### 关键常数（不得在各函数内硬编码）

```python
Rg = 8.314      # J/(mol·K)
g  = 9.81       # m/s²
P0 = 101325.0   # Pa
n_b = 2.7       # Hilligardt 气体交换因子
```

优先使用 `src/core/constants.py`。

### 禁止事项

- ❌ R8 中直接用 Pa 代入压力项
- ❌ 假设 Boudouard k_b/k_c 指数为负
- ❌ 使用常数 `n_RZ = 4.65`（须按 Re_s 分段）
- ❌ 忽略 R5 悬浮相 `C_H2O^0.5` 项
- ❌ 添加 logging、设计模式、抽象基类

### 方程注释

每个核心函数 docstring 须注明方程编号，例如：

```python
"""
Source: docs/hamel_submodels/03_hydrodynamics_core_chain.md; Hamel (1999) Eq.3.50
"""
```

**勿再引用**已清理的 `specs/00–04`、`specs/species.md`。

## 开发顺序（严格按此顺序，不跳步）

| Phase | 内容 |
|-------|------|
| 1 | `species.py`、`arrhenius.py`、`validation_cases.json` |
| 2 | 流体力学（u_mf、气泡、相含率、自由板） |
| 3 | 动力学（炭/气/tar 反应） |
| 4 | `drying.py`、`devolatilization.py` |
| 5 | `cell.py`、`cell_solver.py`、`reactor.py`、`sanity_checks.py` |
| 6 | `test_table2_LU.py`、`app.py` |

## 一致性优先

1. 先核对 **代码 vs `docs/hamel_submodels/00–07`**
2. 不一致时回到 **Hamel 论文 PDF**（页码 / 方程号 / 表格）
3. 禁止用当前实现反向补文档
4. 一致性未闭合前，不得进入调参

## 验收

## 相关文档

| 文档 | 用途 |
|------|------|
| [`docs/README.md`](README.md) | Hamel 算法文档索引 |
| [`specs/hamel_checked_algorithm_spec.md`](../specs/hamel_checked_algorithm_spec.md) | 已核验摘要（本地-only） |
| [`docs/module_dependency_overview.md`](module_dependency_overview.md) | 代码模块依赖 |
