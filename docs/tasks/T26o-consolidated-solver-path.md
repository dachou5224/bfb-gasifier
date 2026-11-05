# T26o：实现 ADR 0054（合并求解路径）

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类与运行记录，CHARTER §8） |
| 前置任务 | ADR 0054 已批准（PR #146）；限时诊断 PR #145、#147 已合并 |
| 分支 | `task/T26o-consolidated-solver-path` |
| 时间上限 | 人工工作 2 个工作日；运行墙钟按 ADR 0037。**先删后加**：`hamel/` 净增以 ≤ 0 为目标，超过 +40 行即停下上报（目标所有者 2026-11-02 决定第 3 项）；`tools/` 净增不超过 +150 行 |

## 来由

限时求解路径诊断（`docs/research/solver-path-bounded-diagnostic.md`）表明，现行规则下外层第 1 轮内层走不动。在 `tools/` 中试过三项求解手段，合用后内层可以收敛：

- 气体出流的对数变量；
- O₂ 出流按 R5 级数取幂变量；
- 三层试探规则。

ADR 0054 把这三项与初值构造合并为一条求解路径，取代 ADR 0045 第 3–4 项、0049、0052 及其修订。本卡把它写进 `hamel/`，同时删除被取代的代码。

本卡只实现求解手段。O4（Corella 有效域）与 O7（灰的闭合）不在本卡内，后者正在做闭合复核（ADR 0055，待起草）。

## 执行要求（给 cursor）

- **新会话，在主仓库目录 `bfb-gasifier/` 内工作，不建 worktree**：确认工作树干净后，`git fetch origin && git switch -c task/T26o-consolidated-solver-path origin/main`。
- **不在仓库目录以外创建或写入任何文件**：临时文件放 `.scratch/`，报告放 `reports/t26o-validation/`。
- **开工前必读**：
  - 本卡与 ADR 0054（全文，含审查记录）；
  - ADR 0017 第 6–8 项、ADR 0043、ADR 0046、ADR 0052（被取代，读懂再删）；
  - 研究底稿第 2–4 节；
  - 诊断工具 `tools/solver_path_diagnostic.py` 中的 `transforms`、`Problem.feasible`、`final_violations`、`run_inner`：它们是 ADR 0054 第 1–2 项的诊断写法，**参照其行为，不得把工具代码整块搬进 `hamel/`**；
  - `hamel/solver/newton.py`、`driver.py`（`_feasible`、`_infeasible_clause`、`prepare`）、`initial.py`（`solid_feasible`、`solve_solids`、`construct`）、`precalc.py`（`_drain`、`_drain_terms`）、`gates.py`（`arrays_allowed`）。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26o 一行。

## 范围内

### A. 先删

1. `initial.solid_feasible`、`initial._infeasible_clause` 及其说明。
2. `driver._feasible`、`_infeasible_clause`、`_bound_hit` 中试探点上的存量与惰性物检查；`gates.arrays_allowed` 若只剩此用途，一并删除。
3. `precalc._drain` 的符号保护在初值构造交替中间的用途：符号检查移到第 C 项第 6 条与外层预计算。
4. ADR 0043“固体子系统可解性预检”的残留（若 `hamel/` 中仍有）。
5. 删除列表逐条写进 PR，附行数。

### B. 实现（ADR 0054 第 1–2 项）

1. **变量变换**：放在 `hamel/solver/newton.py` 或一个新的小模块中，以 y ↔ x 一对函数实现，残差函数不知道变换。
   - O₂ 出流取 s = Ṅ^n：床区悬浮相 n = 0,5，床区气泡相与自由空间 n = 0,3。
   - 其余每轮内层起点上 Ṅ > 0 的气体出流取 ln Ṅ。
   - 起点上 Ṅ = 0 的分量（第 1 关为 Ar）不变换。
   - 差分与列缩放的尺度：对数列 1，幂列 ṅ_ref^n，其余按 ADR 0043。
   - 幂变量的允许域与逆映射严格按 ADR 0054 第 1 项：s < 0 拒绝并减半，不评价残差；不得平方或取绝对值折返；s = 0 可评价。
2. **三层试探规则**：内层与初值构造共用一个函数。
   - 第 1 层（可评价）：自由空间有效域在试探点与差分点上改为缩步，不再 DOMAIN_STOP；起点与收敛解超域仍停止。
   - 第 2 层（试探点约束）：不变换的气体 ≥ −τ_n，幂变量 s ≥ 0，温度 > 0；初值构造另要求最底床格灰分 > 0。
   - 第 3 层（收敛解检查）：存量 ≥ −τ，各床格惰性物 ≥ 0，不变换的气体 ≥ −τ_n；初值构造另要求排灰灰分量 ≥ 0。
   - τ = `stock_roundoff_relative`·M_tot，τ_n = `stock_roundoff_relative`·ṅ_ref；舍入值按原值报告，不改写。
3. 每次迭代记录：y 上最大的步长分量（名与值）、雅可比求解是否成功、第 3 层检查结果（逐格惰性物）。

### C. 初值构造（ADR 0054 第 3 项）

4. 固相子系统沿用现有结构，试探规则换成第 B 项第 2 条。
5. 交替中间的预计算：A5 余额为负时取同式的负值（`_drain_terms`），不停止。
6. 交替结束时检查：最后一次预计算的排灰灰分量 ≥ 0，最后一次固相解满足第 3 层；不满足即停止并报告。
7. 外层每轮预计算的 A5 符号检查不变（ADR 0054 第 4 项）。

### D. 测试

8. 按 ADR 0054“实现约束”逐条写测试：
   - √x 型标量问题在 s 上一步到位、在 x 上过冲；
   - 负 s 被拒且不评价残差；s = 0 可评价；跨零步经减半落回 s ≥ 0 后才接受；n = 0,5 时 ±s 不给出相同结果；
   - 零存量负方向在三层规则下可接受；
   - Ar 列在 −τ_n 内不阻挡整步；
   - 分组雅可比在 y 上与逐列差分逐位相同（ADR 0046 第 5 项）；
   - 初值构造在交替中间负 A5 不停止、结束时停止。
9. 删除或改写被取代的测试，逐条写进 PR。

### E. 验证（ADR 0054 预先登记，不得改动判据）

10. 工况 1，外层 1 轮，不外推 Corella：
    - 初值构造 4 次交替结束，A5 排灰灰分量 +0,153 g/s（±1 %）；
    - 第 1 轮内层不出现 O1、O2、O3、O5 类停止，只能停在最底床格 900 °C 处 Corella 拒绝差分点（O4）。
11. 在 `tools/` 中临时外推 Corella（`hamel.kinetics.homogeneous._CORELLA_HIGH_K = +∞`，不进 `hamel/`），第 1 轮内层 CONVERGED：
    - 迭代 ≤ 100，三块归一残差 ≤ 1·10⁻⁶；
    - 各格温度与 `reports/research-solver-path/diagnostic/outer.json` 中 `rounds[0]` 的 `bed_temperature_c`、`freeboard_temperature_c` 相差 < 1 K。
12. 每次迭代报告 y 上最大步长分量与雅可比求解是否成功；出现 Δy ≥ 10⁶ 或雅可比不可解，停止上报。
13. 任一条不满足即停止上报，不调参数。报告写入 `reports/t26o-validation/`。
14. 默认 pytest、`tools/check_citations.py`、`tools/drift_check.py`。

### F. 工具收尾

15. `tools/solver_path_diagnostic.py` 改为调用正式实现，只保留 `corellaX`、`o2proj` 等诊断变体与报告；`tools/alternation_diagnostic.py` 的符号放宽随第 C 项删除。工具的已有运行记录不重生成。

## 范围外（明确禁止）

- 不改方程、残差、闭合、冻结值与收敛门；不加 clip、下限截断或开关。
- 不外推 Corella 进入 `hamel/`（O4 未裁定）；不改 A5 外层符号检查、床灰粒级或排灰（O7，闭合复核中）。
- 不引入 O₂ 投影（ADR 0054 审查记录“待实现时复核”）：第 E 项若在外层后续轮复现痕量 O₂ 越零，停止上报，不自行处理。本卡只跑外层 1 轮。
- 不整块迁移诊断工具代码；不恢复归档代码。

## 审查记录

- PR #149 审查（codex）：通过。
- **2026-11-04：总司令批准**，可派发。
- PR #152 审查（codex）：两轮。第 1 轮 P1：`Transform` 以 ndarray 作 dataclass 默认值，Python ≥ 3.11 导入失败——改用 default_factory；P2：O4 判定只看“雅可比不能求解”前缀——改为同时要求 Corella 范围拒绝原文与最底床格距 900 °C < 0,01 K，冻结组重跑。第 2 轮：通过（全量 451 passed，唯一失败为已披露的 Bild 7.20 补页检查，与本卡无关）。
- **2026-11-05：目标所有者同意接受 `hamel/` 净增 +95**（超出本卡 +40 停止线；新增为 ADR 0054 第 1、2 项要求的变换与三层规则，被取代代码已删；`hamel/` 合计 10 623 行）。
- **2026-11-05：总司令批准交付**。非阻断记录：`tools/t26o_validate.py` 的 `_BOTTOM_BOUNDARY_K = 900.0` 实为 °C，名称后缀有误，下次改该工具时更正。O4、O7 未由本卡处理；O7 由 ADR 0055 处理，实现卡另起草。
