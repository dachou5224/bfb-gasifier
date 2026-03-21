# 实现范围摘要（与论文离散不等价处）

> 完整论证见 **`docs/validation_gap_analysis.md` §0`**。

## 仍采用论文中的 **控制方程形式**

- 两相 cell、$\dot{N}_{ex,bd}$、气相/固相/能量守恒（见 `01_conservation_equations.md`）。
- 水力学、动力学、干燥/热解子模型方程（`02`–`04`）。

## 与 Hamel **程序结构** 不一致之处（摘要）

1. **求解器**：论文为 **全炉联立 NR + 分块三对角 Jacobian**；代码为 **`Reactor.solve` 扫描 + 每格 `fsolve`**（Gauss–Seidel 式）。
2. **Vorabrechnung**：论文有预算初值；代码 **无**，`_build_cells()` 初值简单。
3. **干燥/DAEM**：论文可在预算中固定产率；代码在 **`residuals()`** 中 **每次** 调用 CN+DAEM（计算代价与 Jacobian 扰动行为不同）。
4. **循环返料**：论文在 **Jacobian 侧边块**闭合；代码在 **扫描外**更新 `N_rez` 等。
5. **拓扑**：论文 **Verbindungsmatrix**（旋风、连接管 cell）；代码 **仅床层线性 cell**，无 ξ>1 段（与 `H_freeboard=0` 等配置一致）。

## 阅读顺序建议

1. `specs/01`…`04` — 方程是否正确实现进 `Cell`。
2. `validation_gap_analysis.md` §0 — 为何与 Fig. 7.4 / 论文曲线**数值上**仍可能差一截。
