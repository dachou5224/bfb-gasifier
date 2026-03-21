# BFB 模型缺失参数与实现缺口清单

> 按子模型分类 | 含文献检索建议 | **源码审查更新于 2026-03-21**

**相关文档**：
- `docs/source_units_audit.md` — **`src/` 单位与量纲约定**（流率 vs 持料量、焓参考、压力 Pa、干燥/热解源项等）。
- `docs/validation_gap_analysis.md` — **出口温度/气相组成与 `validation_cases.json` 偏差过大的原因**（几何、u0 闭环、能量闭合、收敛、标定等），与本文「待确认参数」交叉对照。

---

## 检索约定：与 Hamel (1999) 原著对齐（供文献 / 聊天机器人检索）

以下记号、方程编号与 **Hamel, S. (1999)** 博士论文（Dissertation）正文一致，便于在 PDF 中全文检索或向 AI 提供**锚点字符串**。

| 检索用关键词（建议原样复制） | 含义 |
|------------------------------|------|
| **Gleichung 3.44** / **Gl. 3.44**（**p. 32**） | **\(u_{b,r}\)**（bubble through-flow）；**Heinbockel (1995)** 对 **Hilligardt (1986)** 的加压修正；**\(n_b\)** 见 **Gl. 3.23**（p. 27）；代码 **`u_br`** |
| **Gleichung 3.50** / **Gl. 3.50**（**p. 35**）；亦见 **(4-12)** | 气泡–悬浮相**总传质系数** **\(K_{bd}\)** [1/s]；**对流 + 扩散** 混合（**Sit & Grace (1981)** 渗透理论）；流体力学与化学摩尔衡算的主接口 |
| **(4-2)** | Archimedes / Ergun–Re 关系（最小流化） |
| **(4-4)** | Hilligardt 气泡上升速度 \(u_b\)（含 \(\psi_b\)） |
| **(4-6)** | Hilligardt 气泡直径沿高微分方程 \(\mathrm{d}d_b/\mathrm{d}h\) |
| **Gleichung 5.20** / **Gl. 5.20** | 气体二元扩散系数 \(D_g\) 的关联式 |
| **§5.2.6**；**Gleichung 5.59** / **Gl. 5.59** | 焦油氧化 R10 速率 |
| **Tabelle 5.4** | 焦油 R10/R11 前置因子与活化能（芳香 / 烯烃；Serio；Corella） |
| **Gleichung 5.45** / **Gl. 5.45**（约 p. 94） | 动力学速率与平衡驱动力耦合（\(1-Q_p/K_{eq}\) 类修正） |
| **Gleichung 6.11**、**Gleichung 6.13**；**Tabelle 6.3** | 炭表面反应、有效动力学参数（Braunkohle / Holz / Hausmüll 等） |
| **Gleichung 4.2**（**p. 52–53**） | 颗粒外表面 **\(r=R_0\)** 热流边界：导热通量 = 对流换热（\(\lambda_s\,\mathrm{d}T/\mathrm{d}r=\alpha\,(T_{ws}-T_s)\)） |
| **Gleichung 4.4**（**p. 52–53**） | 修正蒸发焓 **\(h_v'\)**（潜热 + 湿芯自 **\(T_0\)** 升至 **\(T_e\)** 的显热，含 **\(w_{0,tr}\)**） |
| **Gleichung 4.6**（**p. 52–53**） | 干壳导热解析解温度边界：**\(T|_{r=R_0}=T_s\)**，**\(T|_{r=r_e}=T_e\)**（蒸发前沿） |
| **Gleichung 4.9**、**4.10–4.12** | Nu、DAEM 与径向积分（**Kapitel 4**） |
| **Tabelle 2** | 验证工况（如 LU 等）；与 Hamel & Krumm (2001) **Table 2** 对应 |
| **\(R_g\)** | 通用气体常数；文中与代码 `Rg` 对应 **[J/(mol·K)]** |

**小数与指数（德文排版）**：原著中指数常写作 **\(P^{0{,}3}\)**、**\(C^{0{,}5}\)**（逗号为小数点）；检索时亦可试 **0.3** / **0,3** 两种写法。  
**编号勿混**：**Kapitel 4 干燥** 用 **Gleichung 4.2 / 4.4 / 4.6**（pp. 52–53）；表中 **(4-4)、(4-6)** 指 **流体力学 Hilligardt** 括号式编号，与干燥 **Gl. 4.4、Gl. 4.6** 不是同一式。  
**符号并列**：下文在关键处同时给出 **Gleichung 编号** 与代码变量名，便于对照。

---

## 0. 源码审查摘要（`src/`）

以下为对 `bfb-gasifier/src` 全量 Python 模块的**实现状态**归纳，用于与本文“待确认参数”对照。

| 模块 | 路径 | 状态 | 说明 |
|------|------|------|------|
| 常数 | `core/constants.py` | ✅ | `Rg, g, P0, n_b, T_REF` 等 |
| 物种与物性 | `core/species.py` | ⚠️ 部分 | NASA `Cp/H/S` 已覆盖主气相；`TAR1/TAR2` 需 `configure_tar_components_by_fuel()` 配 MW，**完整 NASA 需 `register_tar_component_properties()`**；理想气密度、幂律粘度、**Gl. 5.20**（\(D_g\)）扩散已实现 |
| Cell | `core/cell.py` | ⚠️ 部分 | `calc_hydrodynamics/exchange/reactions/gas_balance/solid_balance/energy_balance/residuals` 已实现；**干燥/热解已作为源项接入（先干燥后热解）**；**d_core** 由 SPM+炭转化率；**粒径类迁移 `m_left/m_right` 未实现**；`mu_g,20=1.8e-5` 固定（TODO 混合气）；`enthalpy_molar(TAR*)` 未注册时焓为 0 |
| Reactor | `core/reactor.py` | ⚠️ 算法 | 多 cell **Gauss–Seidel 扫描** + `solve_cell`；**非**论文全局 NR / 分块三对角 Jacobian；**无 Vorabrechnung、无 Verbindungsmatrix**；循环在 `fsolve` **外**显式更新。与 Hamel 程序差异见 **`docs/validation_gap_analysis.md` §0**。其余：`ER` + `primary_agent` + 元素分析 时由 **`feed_inlet.compute_gas_feeds_mol_s`** 填 **O2/H2O/N2_feed**；`configure_tar_components_by_fuel` 已调用 |
| 单 cell 求解 | `solvers/cell_solver.py` | ✅ | `fsolve` + 数值 Jacobian |
| Arrhenius | `kinetics/arrhenius.py` | ✅ | `k_hobbs/k_standard/k_jensen_r7` + clip |
| 炭反应 | `kinetics/char_reactions.py` | ✅ | R1–R4、`phi_c`、`D_d_A`、`k_l` 串联；**R1 为 Hamel 有效动力学 + C_O2^0.5**；**d_core** 由 `d_core_from_spm_char_conversion` + 转化率 |
| 气相反应 | `kinetics/gas_reactions.py` | ⚠️ 部分 | R5–R9 已实现；R5/R7/R8 用 `thermodynamics.equilibrium` 驱动力；**R9 仍为简化占位参数** |
| 焦油反应 | `kinetics/tar_reactions.py` | ✅ | R10 **Gl. 5.59** / **Tabelle 5.4**（两相同一热力学，芳香/烯烃加权）；R11 **气泡 Serio**、**悬浮 Corella×ρ_cat**；`cell` 中 **C_tar = TAR1+TAR2**；`rho_cat ≈ ρ_s(1−ε_mf)·f_cat` |
| 热力学平衡 | `thermodynamics/equilibrium.py` | ✅ | `K_eq`、反应商、`calc_gibbs_driving_force`；R8 平衡用 **Benson 拟合**，与 Hamel 指数式可并存待统一 |
| 微量物种 Gibbs | `thermodynamics/minor_species.py`、`gibbs_minimizer.py` | ⚠️ | `Cell.use_gibbs_minor` 关闭时走简化 R9 |
| 最小流化 | `physics/minimum_fluidization.py` | ✅ | Ergun / `compute_u_mf` |
| 气泡 | `physics/bubble_dynamics.py` | ✅ | Hilligardt `u_b`、Mori-Wen / Darton / ODE 可选；**`cell.calc_hydrodynamics` 使用 `mori_wen_bubble_diameter`** |
| 相分率 | `physics/phase_fractions.py` | ✅ | `calc_epsilon_b/d`、分段 `n_RZ` |
| 传质 | `physics/mass_transfer.py` | ✅ | `calc_u_br`、`calc_kbd` |
| 自由板 | `physics/freeboard.py` | ✅ 函数已实现 | **`reactor` / 多 cell 顶部未集成**自由板轨迹 |
| 干燥 | `thermal/drying.py` | ⚠️ | Crank-Nicolson 框架；已实现 **Gl. 4.4** 修正蒸发焓与 **Gl. 4.9** Nu，并**接入 `Cell` 源项**；物性参数仍需标定 |
| 热解 DAEM | `thermal/devolatilization.py` | ⚠️ | DAEM + GH 积分 + **Gl. 4.12** 径向体积分；已接入 `Cell`；**元素守恒分配已实现（CO/H2/CH4/TAR 为主，必要时 CO2/H2O 兜底）** |
| 数量级检验 | `tests/sanity_checks.py` | ✅ | 已扩展为 11 项（含 **Gl. 4.4**、**Gl. 4.9**、**Gl. 4.12**、元素守恒分配、干燥-热解源项接入） |

**与旧版文档的差异**：干燥/热解现已耦合进 `Cell` 反应源项（先干燥后热解，含径向 DAEM 体积分），但参数标定与严格式闭合（分层粒径迁移、完整 tar 热力学）仍有缺口。

---

## 概述

- **已通过**：`python tests/sanity_checks.py` 中 **11/11** 门控（原 6 项 + Chapter 4 干燥/热解耦合检查）。
- **与 `validation_cases.json` 对比（默认）**：`tests/test_table2_LU.py` 中 slow 用例**按 JSON**（`CASE_HTW_WESSELING_1.outputs`）断言温度、干基主气相、碳转化率容差；命令行可运行 `python3 tests/test_table2_LU.py` 打印逐项对比。未标定前可设 **`BFB_RELAX_VALIDATION=1`** 跳过严格断言；当前该用例带 **`pytest.mark.xfail`**，标定对齐后应**删除 xfail** 作为硬门控。
- **仍需用户或文献锁定**：动力学中的**燃料专用有效参数**、干燥/DAEM 与 **Hamel 原文一致**的闭合、以及 **reactor 边界条件**（ER、蒸汽/氧比等）。

---

## A. 异相炭反应动力学（`char_reactions.py`）

> 影响：碳转化率、出口 CO/CO₂ 比、温度剖面  
> 优先级：🔴 高  
> **原著锚点**：**Kapitel 5.1**、**Kapitel 6**；**Tabelle 6.3**；**Gleichung 6.11**、**Gleichung 6.13**（表面通量与有效因子）

### 源码备注

- **R1**：已实现 **Hamel 有效动力学**（`k_hobbs` + **Tabelle 6.3**：`coal`→Braunkohle \(k_0=39\)、\(E=50\,\mathrm{kJ/mol}\)；`biomass`→Holz \(k_0=48\)、\(E=53\,\mathrm{kJ/mol}\)），表面通量与 **Gl. 6.11** 串联形式一致：**\(r = k_l\,C_{O_2}^{0{,}5}\)**（代码幂次：`0.5`）。
- **d_core**：**SPM** — \(d_{\mathrm{Kern}} = d_p\,(1-X)^{1/3}\)（与收缩核表述一致），其中 **\(X\)** 由 `Cell._compute_char_conversion()`（稳态炭质量流：\(1 - \dot m_{\mathrm{C,aus}}/\dot m_{\mathrm{C,ein}}\)）得到；R2/R3 共用同 `d_core`。

### R1 炭燃烧 — **Tabelle 6.3**（与代码择一或扩展）

| 燃料（Hamel 命名） | \(k_0\) | 单位（与论文一致） | \(E\) [kJ/mol] | 反应级数 \(n\) | 模型 |
|--------------------|---------|---------------------|----------------|----------------|------|
| Lignite A (Hobbs et al. 引用) | 1,22 | m/(s·K) | 85,6 | 1,0 | Hobbs et al. (1992) |
| **Braunkohle** | 39,0 | \(\mathrm{(kmol/(s^2\cdot m))^{0{,}5}}\) | 50,0 | 0,5 | Shrinking Core |
| **Holz / Wood** | 48,0 | \(\mathrm{(kmol/(s^2\cdot m))^{0{,}5}}\) | 53,0 | 0,5 | Shrinking Particle |
| **Hausmüll** | 40,0 | \(\mathrm{(kmol/(s^2\cdot m))^{0{,}5}}\) | 50,0 | 0,5 | Shrinking Particle |

（表中数字按德文习惯用**逗号作小数点**书写，便于与 PDF 对照；代码中仍为 float 点号。）

### R2–R4、灰层、\(\Phi_c\)、R8 催化（文献数值）

（原 §A 中表格与公式仍适用，作为**标定目标**；实现已覆盖 **Gl. 6.11** / **Gl. 6.13**、\(\Phi_c\) 分段。）

### 仍待确认 / 改进

| 项目 | 当前代码状态 | 说明 |
|------|-------------|------|
| R1 Waste 燃料 | 未单独分支 | Table 6.3 Hausmüll 可与 `fuel` 扩展第三分支 |
| `d_core` 初值 | 无炭流时 X=0 | 与进料/初值猜测一致时需审视收敛 |

---

## B. 均相气相反应（`gas_reactions.py` + `thermodynamics/equilibrium.py`）

> 优先级：🔴 高

### 源码备注

- **R5/R7/R8**：已实现；可逆部分依赖 **`gibbs_molar` + `get_K_eq` + 驱动力**（**Gl. 5.45** 类形式：\(R_{\mathrm{net}} = R_{\mathrm{kin}}\,(1 - Q_p/K_{eq})\)）。
- **R8 WGSR（水煤气变换）**  
  - **原著**：检索 **Gl. 5.45**（约 p. 94）及 **Chen et al. (1987)** 在 Hamel 中的转述；平衡常数经验式常写作 **\(\exp(-36{,}893 + 4019/T)\)** 一类（指数中 **4019**、**36,893** 为常见检索数字；**以你手中论文 PDF 为准**）。  
  - **代码**：`equilibrium.py` 中 WGSR 另用 **Benson** 型 \(\exp(4577{,}8/T - 4{,}33)\) 拟合；**与上式系数不同**，需在验证算例中**统一为一种**并记录依据。
- **R6 CH₄ 氧化**：已实现（de Souza-Santos 形式）。
- **R9 H₂S 氧化**：**仍为简化占位**（`R9_k0`, `R9_E`），待文献或实验标定。

---

## C. 焦油反应（`tar_reactions.py` + `cell.py`）

> 优先级：🔴 高

### 源码备注

- **化学计量**：双代理 `TAR1/TAR2`，`get_lumped_tar_stoichiometry` 与 `specs/species.md` 一致。
- **R10**（**§5.2.6**；**Gl. 5.59**；**Tabelle 5.4**）— 原著常用形式（检索用）：  
  $$R_{C_mH_n} = -k_{10}\,\exp\!\left(-\frac{E}{R_g\,T}\right)\,T\,P^{0{,}3}\,C_{C_mH_n}^{0{,}5}\,C_{O_2}$$  
  代码：`rate_R10(T, C_tar, C_O2, P, fuel_type)`；**气泡/悬浮相同** \(k_{10}\)、\(E/R_g\)，按代理 **芳香（Aromaten）** vs **烯烃/烷烃** 摩尔分率加权；浓度 **\(C_{\mathrm{tar}} = C_{\mathrm{TAR1}} + C_{\mathrm{TAR2}}\)**。  
  **Tabelle 5.4** 中典型常数（示例，便于检索）：芳香 \(k_{10}=20{,}700\)，\(E/R_g = 9\,650\,\mathrm{K}\)；烯烃/烷烃 \(k_{10}=59{,}8\)，\(E/R_g = 12\,200\,\mathrm{K}\)（**Siminski (1972)** 经 Hamel）。
- **R11**  
  - **气泡相**（**Serio et al. (1987)**，**Tabelle 5.4**）：\(R_{\mathrm{Teer}} = -k_{0,R11}\,\exp(-E/(R_g T))\,C_{\mathrm{Teer}}\)，\(k_0 = 5{,}42\times 10^{4}\,\mathrm{s^{-1}}\)，\(E = 100{,}5\,\mathrm{kJ/mol}\) → `rate_R11_bubble`。  
  - **悬浮相**（**Corella et al. (1991)**，**Tabelle 5.4**）：前置因子 **\(0{,}7\,\mathrm{m^3/(s\cdot kg_{Kat})}\)**，\(E = 63{,}1\,\mathrm{kJ/mol}\)，须乘局部催化剂质量密度 \(\rho_{\mathrm{Kat}}\) → `rate_R11_suspension`（`Cell._catalyst_bulk_density()`）。
- **缺口 / 标定**：`SolidProps.catalyst_solid_fraction`（灰/砂催化份额）默认 1.0；若需区分床料与灰分可再拆参数。

---

## D. 干燥（`drying.py`）

> 优先级：🟡 中  
> **原著锚点**：**Kapitel 4**（**pp. 52–53**）；**Gl. 4.2**、**Gl. 4.4**、**Gl. 4.6**、**Gl. 4.9**

### D.1 **Gleichung 4.2**：外表面热流边界（\(r=R_0\)）

颗粒外表面（**\(r=R_0\)**）处，**传入颗粒的导热热流**与**床层对流换热**相等：
$$\lambda_s \left. \frac{\mathrm{d}T}{\mathrm{d}r} \right|_{r=R_0} = \alpha \cdot (T_{ws} - T_s) = \dot{q}(t)$$

| 符号 | 含义 |
|------|------|
| **\(\lambda_s\)** | 干壳导热系数 |
| **\(\alpha\)** | 床层与颗粒间对流换热系数（代码中由 **Gl. 4.9** Nu 得 `convective_htc_from_nusselt`） |
| **\(T_{ws}\)** | 流化床温度（代码中如 `T_bed`） |
| **\(T_s\)** | 颗粒表面温度 |
| **\(\dot q(t)\)** | 表面热流密度 |

### D.2 **Gleichung 4.4**：修正蒸发焓 \(h_v'\)

有效蒸发所需能量 **\(h_v'\)**：潜热 **\(h_v\)** 加上将**湿芯**自初温 **\(T_0\)** 升至蒸发温度 **\(T_e\)** 的显热（干基初含水 **\(w_{0,tr}\)**）：
$$h_v' = h_v + \left( c_w + \frac{c_s}{w_{0,tr}} \right) \cdot (T_e - T_0)$$

| 符号 | 含义 |
|------|------|
| **\(c_w\)** | 水的比热容（代码 `C_WATER`） |
| **\(c_s\)** | 固体比热容（代码 `CP_W` 等） |
| **\(w_{0,tr}\)** | 干基初始含水率（代码参数 `w0_tr`） |

**实现**：`thermal/drying.py` → `corrected_evaporation_enthalpy()` 与上式一致。

### D.3 **Gleichung 4.6**：干壳区温度边界

干壳（**蒸发前沿 \(r_e\)** 与外表面 **\(R_0\)** 之间）导热方程的解析/数值解用下列温度限：
$$T|_{r=R_0} = T_s \quad \text{und} \quad T|_{r=r_e} = T_e$$

| 符号 | 含义 |
|------|------|
| **\(T|_{r=R_0}\)** | 外表面温度 = **\(T_s\)** |
| **\(T|_{r=r_e}\)** | 内部蒸发前沿温度 = **\(T_e\)**（沸点/蒸发温度） |

**实现**：Crank–Nicolson 径向求解（`solve_drying_CN`）中以外部床温与蒸发前沿约束体现该结构。

### D.4 参数与耦合

| 参数 / 项 | 说明 |
|-----------|------|
| `LAMBDA_W`, `RHO_W`, `CP_W`, `H_EVAP`, `T_EVAP`, `h_conv` | 仍为占位或典型值；**\(\lambda_s\)** 与代码湿壳导热需区分时需单独标定 |
| **与 Cell 耦合** | 已通过 `Cell.calc_reactions` 接入干燥源项（先干燥后热解）；**Gl. 4.4** 已实现；**4.2 / 4.6** 与 CN 边界一致；**Gl. 4.9** 用于 \(\alpha\) |

---

## E. 热解 DAEM（`devolatilization.py`）

> 优先级：🟡 中  
> **原著锚点**：**Gleichung 4.10**（挥发分释放）、**4.11**（活化能分布）、**4.12**（径向体积分）；**Kapitel 8**（DAEM 参数与文献）

| 参数 / 项 | 说明 |
|-----------|------|
| `A_DAEM`, `E0_DAEM`, `SIGMA_DAEM` | 代码有默认值；**待与 Hamel §8 / Kapitel 8 或燃料对齐** |
| **tar/gas/char 分配** | 已实现 C/H/O 元素守恒分配（主产物 CO/H2/CH4/TAR，必要时 CO2/H2O 兜底） |
| **与 Cell 耦合** | 已注入 `calc_reactions` 的固/气相源项链（含 **Gl. 4.12** 径向 DAEM 体积分） |

---

## F. 流体力学（`bubble_dynamics.py` / `phase_fractions.py`）

> 优先级：🟡 中

### F.1 已实现与调用链

| 环节 | 实现位置 | 原著锚点（检索） |
|------|----------|------------------|
| 最小流化速度 | `Cell.calc_hydrodynamics` → `compute_u_mf` | **(4-2)** 型 Ergun–\(Re_{mf}\) 关系 |
| 气泡直径 **\(d_b\)** | **`mori_wen_bubble_diameter(h, u0, u_mf, D_bed)`** | Mori–Wen / Darton 代数式；**非** **(4-6)** Hilligardt \(\mathrm{d}d_b/\mathrm{d}h\) 积分主路径 |
| 气泡上升速度 **\(u_b\)** | `bubble_rise_velocity` | 与 **(4-4)** 一致：\(u_b = \psi_b\,(u_0-u_{mf}) + u_{b,\mathrm{single}}(d_b)\)，默认 **\(\psi_b=0{,}76\)** |
| 相分率 | `calc_epsilon_b` / `calc_epsilon_d` | \(\epsilon_b=(u_0-u_{mf})/u_b\)，\(\epsilon_d=1-\epsilon_b\) |
| 备用算法 | `bubble_dynamics.py` 内 | **`darton_bubble_diameter`**、**`integrate_bubble_diameter`（Hilligardt ODE）**、`bubble_diameter_along_bed(method=...)` 均已实现，**但 `Cell` 未调用** |

### F.2 与 Hamel 论文/规格的差异（缺口）

| 项目 | 现状 | 若要对标论文 |
|------|------|----------------|
| **\(d_b\) 模型** | 全床统一用 Mori–Wen 单值（随 **`h_center`** 变化，见 `CellGeometry.h_center`） | 原文若用 **(4-6) Hilligardt** 沿 \(h\) 积分，需在 `Cell` 或 reactor 层改为 ODE 路径或分段标定 |
| **\(\psi_b\), \(\xi_b\)** | `bubble_rise_velocity` / ODE 中有默认；**当前 `Cell` 走 Mori–Wen + 默认 \(\psi_b\)**，**\(\xi_b\) 不参与**主路径 | **(4-6)** 标定时 \(\xi_b\)、\(\lambda_b\) 等需与 **Tabelle 2 / 实验炉** 一致 |
| **`N_or`（分布板孔数）** | `mori_wen_bubble_diameter` 等默认 **`N_or=100`**，进入 **`d_b0`、有效高度 `h_eff`** | 应改为 **真实分布板设计** 或按炉型校核，否则 **\(d_b\)、\(K_{bd}\)** 系统偏差 |
| **\(A_{\mathrm{Bett}}\)** | \(\pi D_{\mathrm{Bett}}^2/4\)，来自 `CellGeometry.D_bed` | 与反应器一致即可；多孔/非圆截面需工程修正时**未建模** |
| **分段床 / 锥段** | 单一直筒 `D_bed` | `validation_cases.json` 中部分工况为锥形入口；**当前 1D cell 无变径** |

### F.3 建议动作

1. **确认论文/验证工况**采用的 \(d_b\) 公式：若坚持 Mori–Wen，记录与 **Kapitel 4** 中 Mori–Wen 式的对应关系；若需 Hilligardt，则在 `calc_hydrodynamics` 中增加 `method` 开关并回归 `sanity_checks` 中 **\(K_{bd}\)** 区间。  
2. 将 **`N_or`**（及必要时 **\(\psi_b\)**）提升为 **`ReactorConfig` / 几何配置** 的可调字段，避免写死 100。  
3. **`physics/constants.py` 与 `core/constants.py` 均含 `n_b=2{,}7`**：传质用 `mass_transfer` 已 `from src.core.constants import n_b`；建议 **删重复或显式 re-export**，避免日后两处不一致。

---

## G. 传质与自由板（`mass_transfer.py` / `freeboard.py`）

### G.0 **Gleichung 3.50**（p. 35）：总传质系数 \(K_{bd}\)（气泡 ↔ 悬浮相）

**Equation 3.50**（**page 35**）定义气泡相与悬浮（乳化）相之间的**总传质系数** **\(K_{bd}\)**，是**流体力学子模型与化学摩尔衡算的主连接**。Hamel 采用**对流 + 扩散**的混合形式，并针对**高压气化**作了衔接（与常数交换率的经典模型如 **Kunii & Levenspiel** 不同，后者在原文中常被对比）。

**数学形式（与论文一致，两种写法等价）：**
$$K_{bd} = \frac{\dot{V}_{b,r}}{V_b} + \sqrt{\frac{144 \cdot D_g \cdot \epsilon_{mf} \cdot u_b}{\pi \cdot d_b^3}} = \frac{3 \cdot u_{b,r}}{2 \cdot d_b} + \sqrt{\frac{144 \cdot D_g \cdot \epsilon_{mf} \cdot u_b}{\pi \cdot d_b^3}}$$

1. **对流项 \(\dfrac{3\,u_{b,r}}{2\,d_b}\)**（与 \(\dot{V}_{b,r}/V_b\) 一致）：描述气体**穿流气泡**的体积通量相对气泡体积的贡献；**\(u_{b,r}\)** 为 **Equation 3.44** 中的相对气泡速度（含 **\((P/P_0)^{-0{,}15}\)** 压力修正）。  
2. **扩散项 \(\sqrt{\dfrac{144\,D_g\,\varepsilon_{mf}\,u_b}{\pi d_b^3}}\)**：气泡边界处的**分子扩散**交换，基于 **Sit & Grace (1981)** 的渗透理论（penetration theory）形式。

| 符号 | 含义 | 单位 |
|------|------|------|
| **\(K_{bd}\)** | 总传质系数 | **1/s** |
| **\(d_b\)** | 局部气泡直径 | m |
| **\(D_g\)** | 气体扩散系数（**Gl. 5.20**） | m²/s |
| **\(\varepsilon_{mf}\)** | 最小流化空隙率 | — |
| **\(u_b\)** | 气泡上升速度 | m/s |
| **\(P\)** | 反应器绝对压力（经 **Gl. 3.44** 进入 \(u_{b,r}\)） | Pa |

**物理要点**：\(P\) 升高时 **\(u_{b,r}\)** 减小 → **对流项减弱** → **高压下传质“回落”**，与 **HTW 等加压气化炉**观测一致；**\(K_{bd}\)** 进而影响相间交换 **\(\dot N_{ex}\)** 与组分预测。

### G.1 床内相间传质（已接入 `Cell`）

| 项目 | 位置 | 原著锚点（检索） |
|------|------|------------------|
| **\(u_{b,r}\) ≡ `u_br`** | `calc_u_br(...)` | **Gl. 3.44**（p. 32）；**\(n_b=2{,}7\)**（**Gl. 3.23**，p. 27）；**\(P_0\)**：论文常 **101 300 Pa** / 代码 **101 325 Pa** |
| **\(D_g\)** | `species.gas_diffusivity_correlation(T, P)` | **Gl. 5.20**；与 **\(K_{bd}\)** 共用 |
| **\(K_{bd}\)** | `calc_kbd(u_br, d_b, D_g, eps_mf, u_b)` | **Gl. 3.50**（p. 35），**Sit & Grace (1981)**；排版中亦见 **(4-12)**；实现即上式右端（`u_br` 即 \(u_{b,r}\)） |
| **敏感性** | — | **\(K_{bd}\)** 随 **\(d_b\)**、**\(u_b\)**、**\(P\)**（经 **\(u_{b,r}\)**）变化；故 §F 中 **\(N_{or}\)、\(d_b\)** 影响传质与 **\(\dot N_{ex}\)** |

### G.2 仍缺或弱耦合的部分

| 项目 | 说明 |
|------|------|
| **\(n_b\) 压力指数** | 当前固定 **\(-0{,}15\)**（**Gl. 3.44**）；若高压工况需与文献不同指数，应单列为标定量 |
| **自由板区几何** | `ReactorConfig` 有 **`H_freeboard`、`D_bed`**，但 **`reactor.py` 求解循环仅覆盖床层 `n_cells`，未在自由板增加 cell 或解析段** |
| **自由板子程序** | `freeboard.py` 提供 **`calc_u_gb`（默认 `1.53*u_b`）**、**`calc_cd_haider`**、**`calc_beta_a`**，用于颗粒速度/阻力/轴向衰减的 **工程近似** | **无任何函数被 `Reactor`/`Cell` 调用** |
| **夹带与出口** | 未实现 **颗粒夹带流率 → 返回床层或损失**；未实现 **自由板内二次反应/换热**（仅床内动力学） |
| **测试占位** | `tests/test_table2_LU.py` 中 **`H_freeboard=0.0` + TODO** 与上述一致 |

### G.3 建议动作

1. 明确产品需求：仅床层预测 vs **含自由板出口温度/组成**。后者需 **轴向延伸网格** 或 **半经验出口模块**。  
2. 若先做最小集成：在 **最顶 cell 出口** 调用 `calc_u_gb` / `calc_beta_a` 等，把 **颗粒相稀释或碳流损失** 以源项或边界反馈回固体衡算（需新增方程）。  
3. 校核 **`n_b`**：与加压实验（2.5 MPa）对比 **K_bd 或示踪响应** 后再锁定 2.7。

---

## G2. Gibbs 与燃料 S/N（`thermodynamics/`、`SolidProps`）

> **原著锚点**：**TechSpec §5.5**；Hamel **Anhang** / **Gl. A.28** 类元素衡算（与 `get_element_release` 注释一致；**以 PDF 中附录编号为准**）

### G2.1 固体侧元素输入（缺省与来源）

| 字段 | 默认 | 作用 |
|------|------|------|
| **`sulfur_fraction`** | `0` | 燃料硫 **质量分数**（干基或项目约定需与 `moisture_wt` 一致）；用于 **`get_element_release("S")`** 与 Gibbs 元素衡算 |
| **`sulfur_volatile_frac`** | `0.5` | 硫在热解中挥发比例 vs 留在炭中；影响 **挥发硫 vs 炭硫** 的分流（见 `Cell.get_element_release`） |
| **`nitrogen_fraction`** | `0` | 燃料氮 **质量分数**；**挥发 N** 由 `m_solid_zu * nitrogen_fraction` 进入元素衡算 |

**缺口**：若不做工业/元素分析标定，**S/N 总进入量为 0**，则 **`use_gibbs_minor=True` 时硫、氮系平衡仍退化**（无元素可分）。应将 **Table 2 / CASE_LU** 中的 **`S_dry`、`N_dry`** 等换算进 `SolidProps`（`reactor._build_cells` 当前 **未**从 `ReactorConfig` 传入硫、氮字段）。

### G2.2 `use_gibbs_minor` 双路径

| 模式 | 行为 | 依赖 |
|------|------|------|
| **`False`（默认）** | 硫氧化走 **`rate_R9`**（简化动力学占位参数，见 §B） | 仅需气相 **H2S、O2** 等浓度 |
| **`True`** | **R9 关闭**，在 `calc_reactions` 中由 **`calc_minor_species_gibbs`** 向 **`R_gas_d` 加松弛源项**，趋向 Gibbs 平衡摩尔流 | `thermodynamics/minor_species.py`：**硫** `H2S/SO2/COS`，**氮** `NH3/HCN/NO`；**`gibbs_molar` / 元素衡算** |

**缺口与注意**：

1. **氮元素**：候选含 **NH3、HCN、NO**，**不含 N₂**；大气化剂中的 **N₂** 仍由主气相守恒处理，**Gibbs 块不负责把 N₂ 分解**。  
2. **`get_element_release`**：为 Gibbs 提供 **S/N 可用量**（挥发固体 + 气相）；若 **硫、氮分数仍为 0**，则 **微量组分求解输入元素量为 0**，结果无意义。  
3. **数值**：Gibbs 最小化对 **`gibbs_molar`**、初值敏感；**`k_relax`**（cell 内取 1.0）过大可能振荡，过小收敛慢——属 **实现调参**，非文献常数。  
4. **与 R1–R4 耦合**：炭消耗会改变 **炭中硫** 释放路径；当前框架有该项，**细粒度含硫矿物形态未区分**。

### G2.3 建议动作

1. 在 **`ReactorConfig`** 增加 **`S_dry_wt` / `N_dry_wt`（或与 C/H/O 一致的 wt%）**，构建 **`SolidProps.sulfur_fraction` / `nitrogen_fraction`**（注意干湿基换算）。  
2. 运行 **Table 2** 时显式选择 **`use_gibbs_minor`** 并对比 **R9 路径**，记录 **出口 H2S、NH₃** 量级。  
3. 若论文给出 **挥发硫比例**，用 **`sulfur_volatile_frac`** 标定，而非长期保留 0.5 默认。

---

## H. Reactor / Cell 集成与边界

> 优先级：🟢～🟡

| 缺失项 | 当前状态 |
|--------|----------|
| `ReactorConfig.O2_feed` 等 | **已实现**：设置 **`ER`**（及 **`primary_agent`**、`S_dry`、`steam_to_o2_molar`）后由 **`core/feed_inlet.py`** 自动计算；**`ER=None`** 时保持手动三股流率（兼容基准测试） |
| `mu_g,20` | Cell 内 **1.8e-5 常数**；论文建议按混合气或 Wilke 改进 |
| 干燥 / 热解 | **已进 `calc_reactions` 源项链**（先干燥后热解）；仍需参数标定与工况验证 |
| 粒径迁移 | **`calc_solid_balance` 注释已说明未实现** |
| **TAR 摩尔焓** | 未 `register_tar_component_properties` 时 **能量项中 tar 焓为 0**（`NotImplementedError` 捕获） |

---

## H2. 脚本 TODO 对照清单（按源码注释）

> 以下条目直接对应 `src/` 与 `tests/` 中仍存在的 `TODO` 注释，优先级按对模型主结果影响排序。

| 优先级 | 文件 | TODO 摘要 | 影响 | 建议处理 |
|---|---|---|---|---|
| ✅ | `src/core/reactor.py` + `feed_inlet.py` | （已关闭）`O2_feed` 由 ER 自动计算 | — | `ReactorConfig(ER=..., primary_agent=...)` |
| 🔴 高 | `src/core/cell.py` | `mu_g,20=1.8e-5` 常数占位 | 影响 `u_mf`、`K_bd`、交换与反应速率 | 引入混合气粘度（如 Wilke）替代常数 |
| 🔴 高 | `src/core/cell.py` | 粒径迁移 `m_left/m_right` 未实现 | 固相守恒简化，影响炭反应面与转化 | 按守恒方程 **Gl. 2-3** / **Eq. 2-3**（`specs/01_conservation_equations.md`）增加粒径类迁移项 |
| 🔴 高 | `src/kinetics/gas_reactions.py` | R9 仍为简化一级动力学参数 | 硫组分预测偏差风险高 | 用文献参数/实验标定替换占位 |
| 🟡 中 | `src/thermal/drying.py` | 干燥物性参数仍为典型值 | 影响干燥速率与热解起始时机 | 按煤种/生物质分组给参数集 |
| 🟡 中 | `src/kinetics/char_reactions.py` | R4 严格 SCM+L-H 联立仍为近似 | 高压/高 CO₂ 条件下精度受限 | 逐步切换为联立求解版本 |
| 🟡 中 | `src/physics/bubble_dynamics.py` | `N_or`/几何默认值仍占位 | `d_b` 与 `K_bd` 对装置敏感 | 将 `N_or` 参数化到 `ReactorConfig` |
| 🟢 低 | `tests/test_table2_LU.py` | 蒸汽/氧比占位、自由板未集成 | 验证场景与论文工况不完全一致 | 完成自由板耦合后再收紧该测试 |

---

## 优先级排序与建议检索路径（修订）

```
第 1 步：统一论文/代码口径
   - R8：Gl. 5.45 驱动力；WGSR 平衡 — Benson 型 vs Hamel/Chen 型 exp(−36,893+4019/T) 二选一并回归测试
   - 焦油：Gl. 5.59 + Tabelle 5.4；Serio / Corella 分相 + TAR1/TAR2 浓度定义
   - R1：Tabelle 6.3；Hausmüll 分支如需与 Braunkohle/Holz 并列

第 2 步：边界与进料
   - ReactorConfig 由 CASE_LU（ER、剂种）自动生成

第 3 步：子模型耦合
   - Gl. 4.4–4.12（干燥/DAEM）→ cell 源项
   - TAR NASA 或等价焓模型，避免焓衡偏差

第 4 步：收尾
   - mu_g 混合规则、粒径迁移、自由板集成；Gl. 3.44/3.50 与 (4-4)/(4-12) 一致性质检
```

---

## 快速填写模板

```text
[R1 / Tabelle 6.3] Braunkohle / Holz / Hausmüll / Hobbs = ______

[R8 / Gl. 5.45] WGSR 平衡：Benson exp(4577,8/T−4,33) vs Hamel/Chen exp(−36,893+4019/T) = ______

[TAR / Gl. 5.59 + Tabelle 5.4] 气泡 Serio / 悬浮 Corella×rho_Kat = Y/N
      C_tar = C_TAR1 + C_TAR2 = Y/N

[Reactor 进料] O2, H2O, N2 由 ER + primary_agent + 元素分析 = 已自动（feed_inlet）

[Gl. 4.4–4.12 干燥+DAEM] 接入 Cell = 已完成（待标定）

[TAR 热力学] register_tar_component_properties = 待填 / 已填

[检索用字符串] Gleichung 4.2; Gleichung 4.4 h_v'; Gleichung 4.6 T_r; Kapitel 4 p.52; Gleichung 3.44 p.32; Gleichung 3.50 p.35; Tabelle 5.4; Tabelle 6.3
```
