# BFB 模型人工输入需求清单

> 目的：汇总当前模型中所有需要人工提供或确认的信息（方程、常数、参数、经验系数、边界条件、工况假设）。
>
> 范围：`src/` 全模块 + 执行计划中的 `PAUSE-1 ~ PAUSE-11`。
>
> 更新日期：2026-03-18

---

## 1. 高优先级（缺失将直接影响模型可信度）

### 1.1 Phase 2 流体力学

#### A) 最小流化参数（PAUSE-2）

- 模块：`src/physics/minimum_fluidization.py`
- 需要人工输入：
  - `eps_mf`（最小流化空隙率）
  - `phi_s`（颗粒球形度）
- 关键方程（Ergun）：

```text
Ar = [150(1-eps_mf)/(phi_s*eps_mf^3)]*Re_mf + [1.75/(phi_s*eps_mf^3)]*Re_mf^2
u_mf = Re_mf * mu_g / (rho_g * d_p)
```

- 当前默认：`eps_mf=0.45`, `phi_s=0.86`（代码中）/ 计划建议 `phi_s=0.8`
- 风险：`u_mf` 对这两个参数高度敏感，直接影响后续 `u_b`、`eps_b`、`K_bd`。

#### B) 气泡动力学关键参数（PAUSE-3）

- 模块：`src/physics/bubble_dynamics.py`
- 需要人工输入：
  - 初始气泡直径公式 `d_b0` 的最终选型（Darton/Mori-Wen/设备标定）
  - 气泡寿命 `lambda_b` 的加压修正公式（指数与系数）
  - `xi_b`（Hilligardt ODE 形状参数）
  - 分布板信息：`A_bed`、`N_or`
- 关键方程：

```text
u_b = psi_b*(u0-u_mf) + u_b,single
dd_b/dh = [2/(9*pi) * eps_b^(1/3)/(1 - xi_b*(6/pi)^(1/3)*eps_b^(1/3))]
         - d_b/(3*lambda_b*u_b)
```

- 当前默认：
  - `psi_b=0.76`
  - `A_bed=0.283 m^2`, `N_or=100`
  - `lambda_b ~ (P/P0)^(-0.2)`（经验写法）
  - `xi_b=0.35`
- 风险：`d_b(h)`、慢泡/快泡判别、`K_bd` 将显著偏移。

#### C) Richardson-Zaki 分段（PAUSE-4）

- 模块：`src/physics/phase_fractions.py`
- 需要人工输入：
  - `n_RZ(Re_s)` 分段公式与分界点（项目最终采用哪套文献）
- 当前实现：有分段，但仍需你确认是否采用该版本。
- 风险：影响相分率与固相速度估计。

#### D) 传质/自由板区参数（PAUSE-5, PAUSE-6）

- 模块：`src/physics/mass_transfer.py`, `src/physics/freeboard.py`
- 需要人工输入：
  - `D_g` 混合扩散系数公式或工况标定口径
  - `beta_A` 衰减模型最终形式
  - Haider-Levenspiel 相关式系数是否按当前版本固定
- 关键方程：

```text
u_br = n_b * u_d * (P/P0)^(-0.15)
K_bd = 3*u_br/(2*d_b) + sqrt(144*D_g*eps_mf*u_b/(pi*d_b^3))
```

- 当前状态：已实现上述形式，但 `D_g` 口径需最终确认（现有 `species.md` 给出 Eq.5.20 形式）。

---

### 1.2 Phase 3 动力学参数（PAUSE-7~9）

#### E) R1-R4 炭反应全参数（PAUSE-7）

- 模块：`src/kinetics/char_reactions.py`
- 需要人工输入：
  - R1-R4 的完整 Arrhenius 参数（`k0`, `E`）
  - R4（Weeda L-H）吸附参数（`K1/K2/K3` 对应 `A0/E`）
  - `d_core` 演化关系、灰层有效扩散 `D_eff` 修正
- 当前状态：多组参数为“占位/典型值”，标注待 Table 5.2 确认。
- 关键提醒：R4 吸附相关指数符号必须按文献定义，不能拍脑袋改号。

#### F) R5-R9 均相反应参数（PAUSE-8）

- 模块：`src/kinetics/gas_reactions.py`
- 需要人工输入：
  - R5-R9 的 `A0/E` 完整参数（Table 6.3 口径）
  - R8 的平衡常数表达式（Chen 1987 具体形式）
  - R8 煤灰催化因子模型
- 当前状态：参数多为典型值或估值，待文献回填。
- 关键方程：

```text
R8: r = k_f * (p_CO*p_H2O - p_CO2*p_H2/K_eq), 其中 p 使用 atm
```

#### G) R10-R11 焦油反应参数（PAUSE-9）

- 模块：`src/kinetics/tar_reactions.py`
- 需要人工输入：
  - R10/R11 的 `k0`, `E`
  - 浓度指数（当前暂按一级）
- 当前状态：参数为估值，且明确 TODO。

---

### 1.3 Phase 4 干燥/热解（PAUSE-10, PAUSE-11）

#### H) 干燥模型原始方程与物性（PAUSE-10）

- 模块：`src/thermal/drying.py`
- 需要人工输入：
  - Agarwal 原始 PDE、边界条件、适用假设
  - 关键物性：`lambda`, `rho`, `cp`, `H_evap` 的燃料相关口径
  - 蒸发前沿处理是否保留当前简化
- 当前状态：Crank-Nicolson 替代模型可运行，但文献参数未锁定。

#### I) DAEM 参数与产物分配（PAUSE-11）

- 模块：`src/thermal/devolatilization.py`
- 需要人工输入：
  - DAEM 参数：`A`, `E0`, `sigma`
  - 热解产物分配：`tar/gas/char` 比例
- 当前状态：默认值为褐煤典型量级，待工况拟合。

---

## 2. 中优先级（影响精度与工程可用性）

### 2.1 物化性质口径统一

- 模块：`src/core/species.py`, `src/core/cell.py`
- 需要人工确认：
  - 混合气黏度基准 `mu_g_20` 是否按组分动态计算（当前 `cell.py` 用常数 `1.8e-5`）
  - tar 代理组分映射策略是否固定为 `TAR1/TAR2`
  - 极端工况下是否引入真实气体修正（当前理想气体）

### 2.2 反应器边界/进料策略

- 模块：`src/core/reactor.py`, `tests/test_table2_LU.py`
- 需要人工确认：
  - `O2_feed`、`H2O_feed` 的由 ER 反推规则
  - 蒸汽/氧比（`tests/test_table2_LU.py` 当前 `H2O = 0.8*O2`）
  - 自由板区是否并入主求解（当前有模块但集成仍简化）

### 2.3 固相粒径迁移

- 模块：`src/core/cell.py`
- 需要人工输入：
  - `m_left/m_right` 粒径类间迁移闭合关系（收缩颗粒）
- 当前状态：守恒式已搭好，迁移项未完整实现。

---

## 3. 低优先级（可后置完善）

- `app.py` 中展示参数与核心求解参数完全一致性复核
- `tests/test_table2_LU.py` 的严格门控（将 NOTE 项改为强断言前需先锁定上文参数）

---

## 4. 常数与可配置项总表（需人工“确认是否固定”）

以下并非全部“缺失”，但建议由你明确“固定值/可调范围/数据来源”：

| 名称 | 当前值 | 模块 | 建议动作 |
|---|---:|---|---|
| `Rg` | 8.314 | `core/constants.py` | 固定，确认即可 |
| `g` | 9.81 | `core/constants.py` | 固定，确认即可 |
| `P0` | 101325 Pa | `core/constants.py` | 固定，确认即可 |
| `n_b` | 2.7 | `core/constants.py`/传质 | 确认是否工况相关 |
| `psi_b` | 0.76 | 气泡动力学 | 确认分布板类型是否一致 |
| `eps_mf` | 0.45 | 流化 | 与燃料类型绑定 |
| `phi_s` | 0.86（代码） | 流化 | 与计划建议值统一 |
| `mu_g_20` | 1.8e-5 | `cell.py` | 建议改为随组分计算 |

---

## 5. 建议的人机交互补数顺序

1. 先锁定 Phase 2：`eps_mf/phi_s`, `d_b0/lambda_b/xi_b`, `D_g/beta_A`
2. 再锁定 Phase 3：R1-R11 全动力学参数（Table 5.2/6.3 + R8 催化/平衡）
3. 再锁定 Phase 4：干燥 PDE 与 DAEM 参数
4. 最后做 Phase 5/6 的严格门控（把 NOTE 改为 assert）

---

## 6. 可直接填写模板（复制即用）

```text
[Phase 2]
eps_mf = 
phi_s = 
d_b0 formula = 
lambda_b pressure correction = 
xi_b = 
D_g formula = 
beta_A formula = 

[Phase 3]
R1: k0= , E=
R2: k0= , E=
R3: k0= , E=
R4: kf(k0,E)= ; kb(k0,E)= ; kc(k0,E)=
R5: k0= , E=
R6: k0= , E=
R7: A= , E_T=
R8: k0= , E= , K_eq(T)= , ash_catalyst=
R9: k0= , E=
R10: k0= , E= , order=
R11: k0= , E= , order=

[Phase 4]
Drying PDE/BC:
lambda_w=
rho_w=
cp_w=
H_evap=
DAEM: A= , E0= , sigma=
tar/gas/char split=
```

