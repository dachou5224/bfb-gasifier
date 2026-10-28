# T26h：自由空间颗粒轨迹改为逐格常系数解析解（ADR 0047）

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，见“执行要求”） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | ADR 0047 已批准（须目标所有者明确批准） |
| 分支 | `task/T26h-cellwise-trajectories` |
| 时间上限 | 3 个工作日；净增代码上限 300 行（不含测试）；`hamel/` 净增须 ≤ 0（预计为负），否则停下上报 |

## 来由

目标所有者 2026-10-28 决定先做两项提高计算效率的改动，以便后续试算。本卡按 ADR 0047 把自由空间轨迹从连续 u_gb(h) 的数值积分改为 Hamel S. 45–49 的逐格常系数解析解。可与 T26i 并行；T26g 在两者之后派发。

## 执行要求（给 cursor）

- **新会话、新 worktree**：`git worktree add ../bfb-T26h -b task/T26h-cellwise-trajectories origin/main`。
- **开工前必读**：`docs/CHARTER.md`（§5、§6.2、§7）、本卡、ADR 0047、ADR 0016（第 6–11 条）、ADR 0026 第 2、3 项、ADR 0028、ADR 0034、ADR 0035 第 4 项、ADR 0036、`docs/research/freeboard-trajectory-cellwise.md` 与 `tools/cellwise_trajectory_study.py`（闭式的推导与核验），以及 `hamel/freeboard/motion.py`、`accounting.py`、`hamel/solver/precalc.py`。
- 涉及 Hamel 式 (3.82)–(3.105)：按 `.claude/skills/hamel-source-check/SKILL.md` 提交原文答卷 `refs/hamel/answers/T26h-cellwise-trajectories.md`（可引用研究底稿的推导，但须自行对照原页）。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26h 一行。

## 范围内

1. 各格气速 ū_g,i = 该格内 u_gb(h) 的高度平均（ADR 0047 第 1 项）。
2. 格内闭式推进与事件（ADR 0047 第 2 项）：Δ > 0、Δ < 0、Δ ≈ 0 三种闭式；到达格界、u_P 穿过 ū_g、转向三类事件；跨界面换格。
3. 归宿与捕获（ADR 0047 第 3 项）：界面稳定点、R* 为两格之并、按 δ_u 与越界距离判定捕获、捕获时刻按最后一次进入 R*；带出、回落；中性格停止。
4. 记账与输运系数、零自由参考存量、各闭合不变（ADR 0047 第 4 项）；自由状态滑移按段解析计算。
5. 删除被撤销的部分（ADR 0047 第 5 项）：数值积分、`freeboard_max_step_s`、格内根搜索、势能井、逐根 δ_h 及其测试；保留的字段不变。
6. 测试：按 ADR 0047“实现约束”各项。
7. 报告 `reports/t26h-trajectories/`：工况 1 初值点第 1 轮预计算的新旧参考系数（a、b、拆分、θ、悬停位置与份额）对照、样本数加倍的加密检查、预计算耗时对照。只报告，不作通过门。

## 范围外（明确禁止）

- 不改 u_gb(h) 的定义（ADR 0016 第 7 条）、阻力式、记账规则、输运系数的定义与各闭合；不加球形度分布（另行裁定）。
- 不运行第 1 关的完整牛顿求解与整批工况。
- 不在残差中加 clip、下限、锚点或开关。

## 交付物

- 代码、测试、原文答卷、对照报告；`docs/STATUS.md` 中 T26h 一行

## 验收证据

- `pytest` 默认集通过且 5 分钟内；`check_citations.py`、`drift_check.py` 通过；`hamel/` 净增 ≤ 0。

## 停止条件

到达时间或行数上限、`hamel/` 净增为正、发现 ADR 0047 未覆盖的情形（例如格内中性、事件无法由闭式确定），或需要越出范围时，停下并在 PR / STATUS 中写明。
