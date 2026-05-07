# Phase 1 审计矩阵：核心代码 vs specs/docs

本文档是调试与验证计划的**执行起点**。
目标不是立刻改代码，而是先固定：
1. 哪些文件必须逐行对照；
2. 它们应该实现哪些方程/算法；
3. 当前已确认的规范-实现偏差；
4. 后续每一轮 coding / debugging / testing 的检查顺序。

---

## 1. 本轮审计范围

优先覆盖以下核心路径：

- `src/core/cell.py`
- `src/core/reactor.py`
- `src/solvers/cell_solver.py`
- `src/solvers/global_nr_solver.py`
- `src/core/feed_inlet.py`
- `src/physics/bubble_dynamics.py`
- `src/physics/phase_fractions.py`
- `src/physics/mass_transfer.py`

对应规范与说明文档：

- `docs/CLAUDE.md`
- `specs/00_implementation_scope.md`
- `specs/01_conservation_equations.md`
- `specs/02_hydrodynamics.md`
- `specs/03_drying_devolatilization.md`
- `specs/04_kinetics.md`
- `specs/species.md`
- `docs/validation_gap_analysis.md`
- `docs/solver_audit_and_stability_report.md`
- `docs/cursor_skills_playbook.md`

---

## 2. 工作流约束（已确认）

来自 `docs/CLAUDE.md` 与 `.cursor/rules/*.mdc`：

- 全程使用 SI 单位：Pa, K, m, mol, kg, s
- Arrhenius 速率只能走 `src/kinetics/arrhenius.py`
- 禁止把 `n_RZ` 固定成常数 4.65，必须按 `Re_s` 分段
- 修改 `src/` / `tests/` / 关键 `docs` 后，必须运行 `tests/sanity_checks.py`
- 优先遵循分阶段调试：先子模型，再单 cell，再 reactor coupling，再文献验证

推荐调试技能：

- `docs/cursor_skills_playbook.md` 中的 `bfb-scientific-audit`
- `bfb-numerical-guard`
- `gstack-investigate`

---

## 3. 文件到规范的主映射

| 代码文件 | 主要规范来源 | 应实现内容 | 当前状态 |
|---|---|---|---|
| `src/core/cell.py` | `specs/01_conservation_equations.md`, `specs/04_kinetics.md`, `specs/03_drying_devolatilization.md`, `specs/02_hydrodynamics.md` | 单 cell 的流体力学、相间交换、R1-R11、干燥/热解源项、气/固/能量残差 | 已读，需逐段核对 |
| `src/core/reactor.py` | `specs/00_implementation_scope.md`, `docs/validation_gap_analysis.md`, `docs/CLAUDE.md` | 多 cell 耦合、边界传播、返料、GS / global NR 总控 | 已读，存在结构性偏差 |
| `src/solvers/cell_solver.py` | `docs/CLAUDE.md`, `docs/solver_audit_and_stability_report.md` | 单 cell 非线性求解策略与数值保护 | 已读，文档表述与实现有偏差 |
| `src/solvers/global_nr_solver.py` | `docs/solver_audit_and_stability_report.md`, `docs/validation_gap_analysis.md` | 全局阻尼 NR、缩放、线搜索、FD Jacobian | 已读，需和文档声明对齐 |
| `src/core/feed_inlet.py` | `tests/validation_case_utils.py`, `docs/CLAUDE.md`, `specs/species.md` | ER → O2/H2O/N2 进料换算 | 已读，待与 validation case 字段逐项比对 |
| `src/physics/bubble_dynamics.py` | `specs/02_hydrodynamics.md` | `u_b`, `d_b(h)`, 慢泡/快泡, Hilligardt / Mori-Wen | 已读，确认存在默认模型差异 |
| `src/physics/phase_fractions.py` | `specs/02_hydrodynamics.md` | `n_RZ`, `epsilon_b`, `epsilon_d` | 已读，形式基本一致 |
| `src/physics/mass_transfer.py` | `specs/02_hydrodynamics.md` | `u_br`, `K_bd` | 已读，形式基本一致 |

---

## 4. 已确认的高优先级偏差（第一轮）

以下偏差已通过代码与文档直接对照确认，不是猜测。

### A. 求解架构与论文不等价

**规范/文档：**
- `specs/00_implementation_scope.md`
- `docs/validation_gap_analysis.md`

**代码：**
- `src/core/reactor.py`
- `src/solvers/global_nr_solver.py`

**结论：**
- 论文目标架构是全炉联立 NR + 分块 Jacobian + Vorabrechnung + 拓扑连接矩阵
- 当前主实现仍以 `Reactor._solve_gauss_seidel()` 扫描为主
- 即使已有 `global_nr`，也仍不是论文原始离散的完全等价实现

**影响：**
- 很多出口偏差不一定来自单个公式错误，而可能来自整体求解结构差异

### B. `cell_solver` 文档口径与实际实现不一致

**规范/文档：**
- `docs/CLAUDE.md` Phase 5.2 仍写 `scipy.optimize.fsolve`

**代码：**
- `src/solvers/cell_solver.py`

**结论：**
- 实际实现优先使用 `root(method="hybr")`，失败后回退 `least_squares`
- 这不一定是错误，但必须在审计文档里显式记录，避免后续误判

### C. 流体力学默认气泡直径模型与 specs 主叙述不一致

**规范/文档：**
- `specs/02_hydrodynamics.md` 更强调 Hilligardt ODE

**代码：**
- `src/core/cell.py` 中 `calc_hydrodynamics()`
- `src/physics/bubble_dynamics.py`

**结论：**
- 当前默认走 `mori_wen_bubble_diameter()`
- Hilligardt ODE 保留为备用，不是活跃主路径

**影响：**
- `d_b`, `u_b`, `epsilon_b`, `K_bd` 全链路都会受影响
- 后续验证必须明确“是在验证 Mori-Wen 版本还是 Hilligardt 版本”

### D. `u0` 由当前状态反推，强耦合风险很高

**规范/文档：**
- `docs/validation_gap_analysis.md`

**代码：**
- `src/core/cell.py` 的 `calc_hydrodynamics()`

**结论：**
- 当 `u0_target is None` 时，`u0` 由当前 `N_b + N_d`, `T`, `P` 反推
- 这会把化学转化与流体力学强耦合到同一步中

**影响：**
- 反应导致总摩尔流变化，进而反过来改 `epsilon_b`, `V_b`, `K_bd`
- 这很可能是整炉不稳、假稳态、气体组成异常的核心因素之一

### E. Phase 1/验证基线参数与文献工况可能错位

**规范/文档：**
- `tests/validation_case_utils.py`
- `docs/validation_gap_analysis.md`

**代码：**
- `tests/validation_case_utils.py`
- `tests/test_table2_LU_global_nr.py`

**结论：**
- 当前基线里固定使用 `d_p=0.5e-3`, `rho_s=1000` 或 `1400`
- 文献/JSON 工况对粒径、返料、床层几何、自由板等可能并不等价

**影响：**
- 在 solver 还没校稳前，直接拿 LU 案例做全局拟合会把结构误差、参数误差、数值误差混在一起

### F. `Cell.calc_reactions()` 需按 thesis 编号重新解释 active paths

**规范/文档：**
- `specs/04_kinetics.md`
- `docs/reaction_numbering.md`

**代码：**
- `src/core/cell.py`
- `src/kinetics/*.py`

**已确认事实：**
- 当前主路径已实际注入 `rate_R1`, `rate_R2`, `rate_R3`, `rate_R4_effective`, `rate_R5_*`, `rate_R6`, `rate_R7`, `rate_R8`, `rate_R9`, `rate_R10`, `rate_R11_*`, `rate_R12`
- 但这些**实现标签**与 **thesis 编号**并不一一同名：
  - `rate_R3` = thesis **R4** 加氢气化
  - `rate_R4_effective` = thesis **R3** Boudouard
  - `rate_R6` = thesis **R7** CH4 oxidation
  - `rate_R7` = thesis **R9** methane reforming
  - `rate_R12` = thesis **R6** H2 oxidation
  - `rate_R9` = implementation sulfur placeholder，不属于 thesis R1–R11 主骨架

**影响：**
- 若脚本/文档继续沿用 legacy label 解释热平衡，会把 thesis **R6/R7/R9** 的物理意义错挂到错误的配置旋钮上
- 后续验证必须明确：讨论的是 **thesis 编号** 还是 **implementation label**

---

## 5. 第一批执行顺序（固定）

### Step 1 — Core/Solver 逐段对照
目标：
- 给 `src/core/cell.py`
- `src/core/reactor.py`
- `src/solvers/cell_solver.py`
- `src/solvers/global_nr_solver.py`

建立函数级审计表：
- 函数名
- 对应 spec/doc
- 预期方程/算法
- 当前实现摘要
- 差异等级：`exact` / `compatible` / `deviation` / `unknown`

### Step 2 — 建立单 cell 解耦审计脚本
目标：
- no-reaction case
- exchange-only case
- reaction-only case

检查：
- `gas balance`
- `solid balance`
- `energy balance`
- 元素守恒
- 反应计量

### Step 3 — 建立流体力学独立审计脚本
目标：
- 固定 `u0`, `T`, `P`, `d_p`
- 输出 `u_mf`, `d_b`, `u_b`, `epsilon_b`, `K_bd`
- 对照 `specs/02_hydrodynamics.md` 与 `tests/sanity_checks.py`

### Step 4 — 再回到 HTW LU
条件：
- 只有当单 cell 守恒与流体力学审计通过后
- 才继续调 `tests/test_table2_LU.py` / `tests/test_table2_LU_global_nr.py`

---

## 6. 立即行动项

### 已完成
- 读 `docs/CLAUDE.md`
- 读 `.cursor` 规则与技能说明
- 读 `specs/00`, `01`, `02`, `03`, `04`, `species`
- 读 `src/core/reactor.py`, `src/core/cell.py`, `src/solvers/global_nr_solver.py`, `src/solvers/cell_solver.py`, `src/core/feed_inlet.py`
- 运行 `tests/test_table2_LU_global_nr.py`

### 下一步（正在执行）
1. 对 `src/core/cell.py` 做函数级逐段审计
2. 明确 R3 / R4 / R12 在主路径中的状态
3. 起草第一个独立审计脚本：单 cell 守恒检查

---

## 7. 审计记录规范

从本文件开始，后续每个审计结论都应带：

- `Evidence`：代码文件 + 规范文件
- `Impact`：对温度 / 气体组成 / 碳转化 / 收敛性的影响
- `Action`：补测试、改实现、还是先记录为架构偏差

建议差异标签：

- `exact`：与规范一致
- `compatible`：公式不同写法但等价/可接受
- `deviation`：明确偏离规范或默认路径不同
- `unknown`：尚未证实

---

## 8. 当前阶段结论

本项目现在最需要的不是立刻“调参数”，而是先把以下三类问题拆开：

1. **规范偏差**：例如 Mori-Wen 默认替代 Hilligardt、`root`/`least_squares` 替代文档中的 `fsolve`
2. **结构偏差**：例如 GS 扫描并不等价于论文全局 NR
3. **实现错误**：例如某个反应没真正进入源项、某个守恒式或焓口径不闭合

只有把这三类分开，后续验证和修正才不会互相污染。
