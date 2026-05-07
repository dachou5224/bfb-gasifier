# BFB 气化炉 1D 模型 — 全景理解笔记

## 本文件用途
记录对整个项目建模方法、实现架构、当前状态和已知缺口的全面理解。

---

## 1. 项目概述

基于 **Hamel & Krumm (2001)** 和 **Hamel (1999) 博士论文**的一维稳态鼓泡流化床（BFB）气化炉模型，Python 实现。

**验证目标**：HTW Wesseling 加压炉（Table 2 LU 工况），褐煤，P=2.5 MPa
- 出口温度误差 < ±10%（实验 ~1120K）
- 碳转化率误差 < ±10%（实验 ~95%）  
- CO 摩尔分数误差 < ±15%（实验 ~0.25）

**当前验证结果**（稳定配置 dense=0.40）：
- T_exit: 1093K vs 1120K (0.64% 误差) ✅
- CO₂ (dry): 11.75% vs 11% (6.8%) ✅
- 碳转化率: 100% vs 95% (5.26%) ✅

---

## 2. 两相模型核心架构

### 2.1 Cell 离散化
反应器沿轴向离散为串联 cell（默认 n_cells=10），每个 cell 包含：
- **气泡相 (B)**：纯气体，仅发生均相反应 R5-R11，高速 u_b
- **乳化/悬浮相 (D)**：气体+固体，异相 R1-R4 + 均相 R5-R11，最小流化速度 u_mf

### 2.2 相间传质（核心耦合方程 Eq. 3.50）
$$\dot{N}_{ex,bd,j,i} = K_{bd,i} \cdot V_{b,i} \cdot (C_{j,b,i} - C_{j,d,i})$$

K_bd 由对流项 + 扩散项组成：
$$K_{bd} = \frac{3 u_{br}}{2 d_b} + \sqrt{\frac{144 D_g \epsilon_{mf} u_b}{\pi d_b^3}}$$

压力效应：$u_{br} = n_b \cdot u_d \cdot (P/P_0)^{-0.15}$

### 2.3 守恒方程

**气相摩尔守恒（悬浮相）**：
$$0 = \dot{N}_{zu,d} + \dot{N}_{rez,d} + \dot{N}_{d,i+1} + \dot{N}_{r,d} - \dot{N}_{d,i} - \dot{N}_{ex,bd}$$

**气相摩尔守恒（气泡相）**：
$$0 = \dot{N}_{zu,b} + \dot{N}_{rez,b} + \dot{N}_{b,i+1} + \dot{N}_{r,b} - \dot{N}_{b,i} + \dot{N}_{ex,bd}$$

**固相质量守恒**（含粒径分级）：
$$0 = \dot{m}_{zu} + \dot{m}_{rez} - \dot{m}_{aus} + \dot{m}_{auf,i+1} - \dot{m}_{auf,i} + \dot{m}_{ab,i-1} - \dot{m}_{ab,i} + \dot{m}_r$$

**全局能量守恒**：
$$0 = \dot{H}_{zu} + \dot{H}_{rez} + \dot{H}_{auf,i+1} + \dot{H}_{ab,i-1} - \dot{H}_{auf,i} - \dot{H}_{ab,i} - \dot{H}_{aus} - \dot{Q}_W$$
（焓值必须包含标准生成焓，反应热自动进入）

---

## 3. 反应体系 (R1-R11)

### 3.1 异相炭反应 (R1-R4)
| 编号 | 反应 | 模型 | Arrhenius 形式 |
|------|------|------|---------------|
| R1 | C + O₂ → CO₂ | SPM (k_l 串联) | k_hobbs |
| R2 | C + H₂O → CO + H₂ | SCM | k_hobbs |
| R3 | C + H₂ → CH₄ | SCM | k_hobbs |
| R4 | C + CO₂ → 2CO (Boudouard) | Weeda L-H | k_hobbs |

### 3.2 均相气相反应 (R5-R9)
| 编号 | 反应 | 特殊处理 |
|------|------|---------|
| R5 | CO + ½O₂ → CO₂ | 气泡/悬浮不同表达式；悬浮含 C_H₂O^0.5 |
| R6 | CH₄ + 2O₂ → CO₂ + 2H₂O | de Souza-Santos 形式 |
| R7 | CH₄ + H₂O → 3H₂ + CO | 可逆，(1-Qp/Keq) |
| R8 | CO + H₂O ⇌ CO₂ + H₂ (WGSR) | P 用 atm，可逆 |
| R9 | H₂S + ½O₂ → SO₂ | 可被 Gibbs 替代 |

### 3.3 焦油反应 (R10-R11)
- **R10**：焦油氧化（气泡/悬浮相同动力学）
- **R11**：焦油裂解（气泡: Serio 热裂解；悬浮: Corella 催化裂解）

### 3.4 Gibbs-动力学耦合 (Eq. 5.45)
$$R_{net} = R_{kinetic} \cdot (1 - Q_p/K_{eq})$$
- R5: max(0, 1-Qp/Keq)（不允许逆反应）
- R7, R8: (1-Qp/Keq)（允许逆反应）
- 微量物种 (H₂S, NH₃等): Gibbs 自由能最小化

### 3.5 Arrhenius 工厂函数（必须使用）
```python
k_hobbs(k0, E, T)      # R1-R4: k = k0·T·exp(-E/(Rg·T))
k_standard(k0, E, T)   # R6 等: k = k0·exp(-E/(Rg·T))
k_jensen_r7(A, E_T, T) # R7: k = (A/T)·exp(-E_T/T)
```

---

## 4. 流体力学

### 4.1 最小流化速度 u_mf
Ergun 方程（Ar-Re_mf 关系）求解

### 4.2 气泡动力学
- **上升速度**：u_b = ψ_b(u₀ - u_mf) + u_{b,i}，ψ_b=0.76
- **气泡直径**：当前用 Mori-Wen 代数模型（非 Hilligardt ODE）
- **慢泡/快泡判别**：α_b = u_b/u_mf

### 4.3 相体积分数
n_RZ 必须按 Re_s 分段计算（禁止使用常数 4.65）

---

## 5. 干燥与热解

### 5.1 干燥模型
球形颗粒：外干壳 + 内湿润区，蒸发前沿半径 r_e
- Eq. 4.4: 修正蒸发焓 h_v' = h_v + (c_w + c_s/w₀_tr)(T_e - T₀)
- Eq. 4.9: Nu_p = 2 + 1.2 Re^{1/2} Pr^{1/3}
- 实现: Crank-Nicolson FD（替代 Agarwal 解析解）

### 5.2 热解 DAEM
分布活化能模型 + Gauss-Hermite 积分
- Eq. 4.10: 挥发分释放方程
- Eq. 4.12: 径向体积积分
- 元素分配: C/H/O 守恒 → CO/H₂/CH₄/TAR + CO₂/H₂O fallback

---

## 6. 求解策略

### 6.1 论文方法 vs 代码实现（关键差异！）

| 方面 | 论文 (Hamel) | 代码实现 |
|------|-------------|---------|
| **求解器** | 全局 Newton-Raphson + 块三对角 Jacobian | Gauss-Seidel 扫描 + 逐 cell fsolve |
| **预算 (Vorabrechnung)** | 平均温度一次预算 → 化学一致 x₀ | 缺失；硬编码初始猜测 |
| **干燥/DAEM** | 预算中一次计算，固定产率 | 每次外迭代重算（用当前 cell T） |
| **循环** | Jacobian 侧块 (Sherman-Morrison) | 后扫描显式赋值，需额外外迭代 |
| **拓扑** | 连接矩阵（旋风、管道 cell） | 硬编码线性序列，无旋风/管道 |

### 6.2 两条求解路径
- `solver="gauss_seidel"`（默认，生产可用）：~4s 收敛
- `solver="global_nr"`（调试中）：接近论文，但初值敏感

### 6.3 数值稳定性措施
- 残差归一化（能量 vs 组分量级差 8 个数量级）
- 反应预步进（Gauss-Seidel 预热）
- 供给限制（O₂ 消耗 ≤ 95% 供给）
- np.exp 溢出保护

---

## 7. 气体物种

11 组分固定物种列表：
CO, CO₂, H₂, H₂O, CH₄, O₂, N₂, H₂S, NH₃, TAR1, TAR2

焦油代理分子：
- 煤：TAR1=C₆H₆, TAR2=C₁₀H₈
- 生物质：TAR1=C₁₀H₈, TAR2=C₁₆H₃₄

热力学：NASA 7 系数多项式

---

## 8. 代码结构

```
src/
├── core/
│   ├── cell.py          # 单 cell 守恒残差方程
│   ├── reactor.py       # 多 cell 扫描迭代（外循环）
│   ├── species.py       # 气体热力学、NASA 系数
│   ├── constants.py     # 物理常数 Rg, g, P0
│   ├── composition.py   # 组分处理
│   └── feed_inlet.py    # 进料计算
├── physics/
│   ├── minimum_fluidization.py  # Ar, Re_mf, u_mf
│   ├── bubble_dynamics.py       # u_b, d_b (Mori-Wen)
│   ├── phase_fractions.py       # n_RZ, ε_d, ε_b
│   ├── mass_transfer.py         # K_bd (Eq. 3.50)
│   └── freeboard.py             # 自由板区（未集成）
├── kinetics/
│   ├── arrhenius.py       # 三种工厂函数
│   ├── char_reactions.py  # R1-R4
│   ├── gas_reactions.py   # R5-R9
│   └── tar_reactions.py   # R10-R11
├── thermal/
│   ├── drying.py            # Crank-Nicolson 干燥
│   └── devolatilization.py  # DAEM
├── thermodynamics/
│   ├── equilibrium.py      # K_eq, Q_p, 驱动力
│   ├── gibbs_minimizer.py  # Gibbs 自由能最小化
│   └── minor_species.py    # S/N 微量物种
└── solvers/
    ├── cell_solver.py       # 单 cell fsolve
    ├── global_nr_solver.py  # 全局 NR（调试中）
    └── vorabrechnung.py     # 预算模块
```

---

## 9. 已知缺口与阻塞项

### 🔴 高优先级
1. **u₀ 强耦合**：u₀ 从化学计算反算 → 流体力学-化学不可分离 → 振荡
2. **求解器架构**：GS vs 全局 NR 决策未最终确定
3. **O₂ 预算不完整**：气泡相 R5b/R6b 未进统一限制器

### 🟡 中优先级
4. **R3/R4 集成状态不明**：可能未激活
5. **WGSR 双公式**：两种动力学表达式并存
6. **气泡直径模型**：用 Mori-Wen 而非规范要求的 Hilligardt ODE
7. **R10 k₀ 待文献验证**：临时降低 100 倍

### 🟢 低优先级
8. **自由板区未集成**：模块存在但反应器未使用
9. **干燥/DAEM 校准**：物性参数非燃料特异
10. **连接矩阵缺失**：无旋风/管道 cell

---

## 10. 源代码实现详情

### 代码规模
- **总计 ~6,400 行 Python**
- Phase 1-6 全部完成，HTW LU 验证通过

### 模块清单

| 模块 | 行数 | 状态 | 关键内容 |
|------|------|------|---------|
| core/ | ~2700 | ✅ | cell.py(21KB 单cell残差), reactor.py(43KB 多cell扫描), species.py(NASA多项式), cell_balances/kinetics/hydrodynamics/pyrolysis 辅助模块 |
| physics/ | ~650 | ✅ | minimum_fluidization(Ergun), bubble_dynamics(Darton/Mori-Wen), mass_transfer(K_bd), phase_fractions(n_RZ分段) |
| kinetics/ | ~900 | ✅ | arrhenius(三工厂), char_reactions(R1-R4 SPM), gas_reactions(R5-R9,R12), tar_reactions(R10-R11两组分代理) |
| thermal/ | ~850 | ✅ | devolatilization(DAEM 10点GH积分), drying(Crank-Nicolson径向) |
| solvers/ | ~850 | ✅ | cell_solver(least_squares+homotopy), global_nr_solver(FD Jacobian+稀疏LU), vorabrechnung |
| thermodynamics/ | ~700 | ✅ | equilibrium(K_eq,Qp,驱动力), gibbs_minimizer(Lagrange乘子法), minor_species(SO₂/COS/HCN/NO) |
| tests/ | ~2500 | ✅ | 16个测试文件 + sanity_checks(11项门控) |
| app.py | ~400 | ✅ | Streamlit 交互式仪表盘 |

### Cell 状态变量 (x 向量)
```
x = [N_d(11), N_b(11), m_solid(nk×4), T]  ≈ 27 维
```
- N_d, N_b: 悬浮/气泡相 11 组分摩尔流量 [mol/s]
- m_solid[nk,4]: 粒径分级固体 [炭,挥发分,水分,灰] [kg/s]
- T: cell 温度 [K]

### 内部依赖图
```
reactor.py
  ├→ cell.py → cell_balances/hydrodynamics/kinetics/pyrolysis
  ├→ cell_solver.py / global_nr_solver.py
  ├→ physics/ (u_mf, d_b, u_b, K_bd, ε_b)
  ├→ kinetics/ (R1-R12 速率)
  ├→ thermal/ (干燥, DAEM)
  └→ thermodynamics/ (K_eq, Gibbs)
```

### 反应实现完整清单 (R1-R12)

| R# | 反应 | 相 | 速率函数 | Arrhenius 形式 |
|----|------|---|---------|---------------|
| R1 | C + (1/Φ_c)O₂ → CO/CO₂ | 悬浮 | rate_R1 | k_hobbs |
| R2 | C + H₂O → CO + H₂ | 悬浮 | rate_R2 | k_hobbs |
| R3 | C + 2H₂ → CH₄ | 悬浮 | rate_R3 | k_hobbs |
| R4 | C + CO₂ → 2CO | 悬浮 | rate_R4_effective | k_hobbs (L-H) |
| R5 | 2CO + O₂ → 2CO₂ | B/D 不同 | rate_R5_bubble/suspension | k_standard |
| R6 | CH₄ + O₂ → 产物 | B/D | rate_R6 | k_standard |
| R7 | CH₄ + H₂O → CO + 3H₂ | B/D | rate_R7 | k_jensen_r7 |
| R8 | CO + H₂O ⇌ CO₂ + H₂ | B/D | rate_R8 | k_standard |
| R9 | H₂S + ½O₂ → SO₂ | B/D | rate_R9 | 占位/Gibbs |
| R10 | TAR + O₂ → CO + H₂ | B/D | rate_R10 | k_hobbs(P^0.3) |
| R11 | TAR + H₂O → CO+H₂+CH₄ | B:Serio/D:Corella | rate_R11 | k_standard |
| R12 | 2H₂ + O₂ → 2H₂O | B/D | rate_R12 | k_standard |

---

## 11. 开发规则总结

1. **单位制**：全程 SI，唯一例外 R8 内部 atm 转换
2. **Arrhenius**：必须用 arrhenius.py 工厂函数
3. **常数**：从 constants.py 导入，禁止硬编码
4. **注释**：每个函数 docstring 标注 Hamel 方程编号
5. **验证**：每次修改后运行 `python3 tests/sanity_checks.py` → "ALL 11 SANITY CHECKS PASSED"
6. **语言**：项目交流统一简体中文
