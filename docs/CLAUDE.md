# BFB Gasifier 1D Model — Agent Instructions (CLAUDE.md)
**基于**：Hamel & Krumm (2001) + Hamel (1999) 技术说明书 v10.0

---

## 项目概述

一维稳态鼓泡流化床（BFB）气化炉模型，Python 实现。
沿轴向离散为串联 cell，每个 cell 分气泡相和乳化相，耦合流体力学、干燥、热解、化学动力学。

**验证目标**（Table 2 工况 LU，HTW 加压炉，褐煤，P=2.5 MPa）：
- 出口温度误差 < ±10%（实验约 900°C）
- 碳转化率误差 < ±10%（实验约 85%）
- CO 摩尔分数误差 < ±15%（实验约 0.25）

---

## 论文图示（Bild 2.2 / 2.3）

- **程序结构（Bild 2.2）** 与 **单格气泡/悬浮摩尔衡算（Bild 2.3）** 的 **Mermaid** 图源已归档在 **`README.md`**（「Hamel 原文图示」节）、**`docs/BFB_TechSpec_v11.md` §0**、**`docs/hamel_dissertation_vs_python_architecture.md` §0**。
- 与 `specs/01_conservation_equations.md` 中 Eq. 3.50（$K_{bd}$）及两相源项符号一致。

---

## 项目结构

```
gasifier_1d/
├── CLAUDE.md                  ← 本文件（Agent 每次必读）
├── specs/
│   ├── 01_conservation_equations.md
│   ├── 02_hydrodynamics.md
│   ├── 03_drying_devolatilization.md
│   └── 04_kinetics.md
├── src/
│   ├── physics/
│   │   ├── __init__.py
│   │   ├── minimum_fluidization.py
│   │   ├── bubble_dynamics.py
│   │   ├── phase_fractions.py
│   │   ├── mass_transfer.py
│   │   └── freeboard.py
│   ├── kinetics/
│   │   ├── __init__.py
│   │   ├── arrhenius.py       ← 三种 Arrhenius 形式的工厂函数（必须用这里的）
│   │   ├── char_reactions.py  ← R1–R4
│   │   ├── gas_reactions.py   ← R5–R9
│   │   └── tar_reactions.py   ← R10–R11
│   ├── thermal/
│   │   ├── drying.py
│   │   └── devolatilization.py
│   ├── core/
│   │   ├── cell.py
│   │   ├── reactor.py
│   │   └── species.py
│   └── solvers/
│       └── cell_solver.py
├── tests/
│   ├── sanity_checks.py       ← 每次修改后必须运行（见下文「修订后工作流」）
│   └── test_table2_LU.py      ← 最终验证
├── data/
│   ├── validation_cases.json  ← **验证工况主数据**（嵌套结构；Table 2 LU ↔ `CASE_HTW_WESSELING_1`）
│   └── test_cases.json        ← 兼容：扁平 `CASE_LU` 快照，与上键等价；**新增/改工况请以 validation_cases.json 为准**
└── app.py                     ← Streamlit 界面（最后实现）
```

---

## 修订后工作流（Agent / 开发者）

每次对 **`src/`**、**`tests/`** 或影响数量级的 **`specs/`**、**`docs/`** 做实质性修改后，在收尾前**必须**运行：

```bash
cd bfb-gasifier && python3 tests/sanity_checks.py
```

- **期望输出**：`ALL 11 SANITY CHECKS PASSED`（含 Chapter 4 干燥/热解相关项）。
- **自动化**：可在**后台**执行，**无需**每次征得用户同意；若沙箱内 NumPy 异常，改用本机 `python3` 或请求完整权限后重试。
- 失败时不得宣称「任务完成」，应修复或说明阻塞原因。

**pytest（可选）**：

- `ReactorConfig.n_age_classes` 默认 **1**，映射到 `SolidProps.n_size_classes`（多粒径/龄期类未启用）。
- `tests/test_table2_LU.py`：`CASE_HTW_WESSELING_1` 与 `data/validation_cases.json` **默认数值对照**（slow，约数十秒～1 分钟）；标定前可设 **`BFB_RELAX_VALIDATION=1`** 跳过与 JSON 的严格断言。
- 单位/量纲约定见 **`docs/source_units_audit.md`**（与 `missing_parameters_summary.md` 交叉引用）。

---

## ⚠️ 全局编程规则（必须遵守）

### 1. 单位制
- **全程 SI**：Pa, K, m, mol, kg, s
- 唯一例外：R8 WGSR 速率方程（Eq.5.45）中 P 单位为 **atm**
  ```python
  P_atm = P_Pa / 101325.0   # 在 R8 函数内部换算，外部仍传 Pa
  ```

### 2. Arrhenius 形式——**禁止混用**，只能调用 `arrhenius.py` 中的工厂函数

```python
# src/kinetics/arrhenius.py 中定义，其他文件 import 使用

def k_hobbs(k0, E, T):
    """形式 A：R1–R4（炭气化）。k = k0 * T * exp(-E/(Rg*T))"""
    return k0 * T * np.exp(-E / (8.314 * T))

def k_standard(k0, E, T):
    """形式 C：R6 等标准形式。k = k0 * exp(-E/(Rg*T))"""
    return k0 * np.exp(-E / (8.314 * T))

def k_jensen_r7(A, E_T, T):
    """形式 B：R7 专用。k = (A/T) * exp(-E_T/T)，E_T=E/Rg 单位 K"""
    return (A / T) * np.exp(-E_T / T)
```

### 3. 关键物理常数（不要硬编码在各函数内）

```python
Rg   = 8.314      # J/(mol·K)
g    = 9.81       # m/s²
P0   = 101325.0   # Pa（标准大气压，加压修正基准）
n_b  = 2.7        # 气体交换因子（Hilligardt 常数）
```

### 4. 方程编号注释
每个函数必须在 docstring 中注明方程编号：
```python
def bubble_rise_velocity(u0, u_mf, d_b, nu_d, psi_b=0.76):
    """
    气泡上升速度。
    Source: specs/02_hydrodynamics.md §2.1, Hamel (1999) Eq.4.4/4.5
    """
```

### 5. 禁止事项
- ❌ 禁止使用标准 `k = A*exp(-E/RT)` 形式实现 R1–R4
- ❌ 禁止在 R8 中直接用 Pa 代入压力项
- ❌ 禁止假设 k_b 和 k_c（Boudouard）的指数为负（两者为正，代表吸附热）
- ❌ 禁止使用 `n_RZ = 4.65`（常数），必须按 Re_s 分段计算
- ❌ 禁止忽略 R5 悬浮相中的 `C_H2O^0.5` 项
- ❌ 禁止添加 logging、设计模式、抽象基类（保持简单）

---

## 开发顺序（严格按此顺序，不跳步）

```
Phase 1（基础，无依赖）
  [1.1] src/core/species.py        — 气体热力学属性、生成焓
  [1.2] src/kinetics/arrhenius.py  — 三种 Arrhenius 工厂函数
  [1.3] data/validation_cases.json — 验证工况数据集（含 HTW/VTT 等；Table 2 LU 见 `CASE_HTW_WESSELING_1`）

Phase 2（流体力学）
  [2.1] src/physics/minimum_fluidization.py  — Ar, Re_mf, u_mf
  [2.2] src/physics/bubble_dynamics.py       — u_b, ODE, K_bd
  [2.3] src/physics/phase_fractions.py       — n_RZ, epsilon_d, epsilon_b
  [2.4] src/physics/freeboard.py             — beta_A, u_gb, C_D_haider

Phase 3（动力学）
  [3.1] src/kinetics/char_reactions.py       — R1(k_l串联), R2, R3, R4(Weeda)
  [3.2] src/kinetics/gas_reactions.py        — R5(两相), R6, R7, R8, R9
  [3.3] src/kinetics/tar_reactions.py        — R10, R11

Phase 4（热解）
  [4.1] src/thermal/drying.py                — Agarwal（Crank-Nicolson FD）
  [4.2] src/thermal/devolatilization.py      — DAEM（Gauss-Hermite）

Phase 5（求解器）
  [5.1] src/core/cell.py                     — 单 cell 守恒方程
  [5.2] src/solvers/cell_solver.py             — scipy.optimize.fsolve
  [5.3] src/core/reactor.py                  — 多 cell 扫描迭代（外循环）
  [5.4] tests/sanity_checks.py               — 数量级验证

> **与 Hamel 论文离散差异**：论文为**全局 Newton–Raphson**、**分块三对角 Jacobian**、**Vorabrechnung** 与**连接矩阵**（旋风/连接管 cell）；本实现为 **Gauss–Seidel 式底→顶扫描**，**无**预算模块与拓扑矩阵。详见 **`docs/validation_gap_analysis.md` §0**。

Phase 6（验证与界面）
  [6.1] tests/test_table2_LU.py              — 工况 LU 端到端验证
  [6.2] app.py                               — Streamlit
```

---

## 数量级验证（sanity_checks.py 中实现）

**每次代码/规格修订后**运行 `python3 tests/sanity_checks.py`（见上文「修订后工作流」）。以下为各检验含义：

| 检验 | 输入条件 | 期望范围 |
|------|---------|---------|
| u_mf | 褐煤 d_p=1mm, T=900K, P=2.5MPa | 0.02–0.08 m/s |
| d_b(H_bed) | u0=0.3 m/s, P=2.5 MPa, H=5 m | 0.05–0.3 m |
| K_bd | d_b=0.1m, P=2.5MPa, T=1000K | 1–15 s⁻¹ |
| R_boudouard | T=1073K, P_CO2=5e4 Pa | ~1e-5 mol/(m²·s) |
| R8 WGSR | T=1073K, P=2.5MPa, y_CO=0.25 | 符号正确（朝平衡方向） |
| DAEM | 升温至 900°C，积分 100s | 挥发分释放率 > 80% |

---

## 工况 LU 输入参数（Table 2，Hamel 2001）

**权威来源**：`data/validation_cases.json` → **`CASE_HTW_WESSELING_1`**（Table 7.1 Sim Nr.1 = 文献 Table 2 LU，HTW Wesseling Air/Steam）。

- 代码与测试（`tests/test_table2_LU.py`）从该 JSON 读取并**展平**为扁平字段；不再以 `test_cases.json` 为唯一依据。
- `data/test_cases.json` 中的 `CASE_LU` 仅为与上键等价的扁平快照，便于人工对照；若数值冲突，**以 `validation_cases.json` 为准**。

```python
# 展平后与下列量级一致（详见 JSON inputs.*）
CASE_LU_FLAT ≈ {
    'reactor': 'HTW_pressurised',
    'fuel': 'brown_coal_RB',
    'P': 2.5e6,            # Pa（inputs.operating_conditions.pressure_MPa）
    'T_inlet': 293.0,      # K（inputs.operating_conditions.T_inlet_K）
    'fuel_feed': 3377.4,   # kg/h（inputs.fuel.feed_rate_kg_h）
    'primary_agent': 'air_steam',  # Luft/Dampf
    'ER': 0.337,
    'recirculation': True,
    'H_bed': 14.5,         # m（inputs.reactor.height_m）
    'D_freeboard': 0.6,    # m
    'moisture_wt': 16.9,   # wt%
    'ash_dry_wt':  11.41,  # wt%
    'VM_daf':      53.42,  # wt%
    'C_dry':       61.5,   # wt%
    'H_dry':       4.1,    # wt%
    'O_dry':       21.8,   # wt%
    'N_dry':       0.68,   # wt%
    'S_dry':       0.51,   # wt%
}
```

---

## 当前已知缺口（不阻塞核心实现）

1. R10/R11：已实现 Hamel Eq.5.59 / Table 5.4（气泡/悬浮 R11 分 Serio/Corella）；若与实验对标仍可调 `SolidProps.catalyst_solid_fraction` 等
2. Agarwal 解析解：用 Crank-Nicolson FD 替代，标注 `NOTE`
3. 泥炭/锯屑 x_FI 参考点（T_ref, P_ref）：暂用 T_ref=1073K, P_ref=101325Pa
4. 自由板区速度分布标准差 σ：暂忽略分布，仅用均值 u_p0=1.53*u_b
