# T26v：实现 ADR 0059（Corella 高温段外推至 950 °C；敏感性 S20）

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | ADR 0059 已接受；T26u 已合并（PR #170） |
| 分支 | `task/T26v-corella-extrapolation` |
| 时间上限 | 人工工作 1 个工作日，另加两次外层运行（各不超过 8 小时墙钟，可并行）。`hamel/` 净增目标 ≤ +30 行，超过 +50 行停下上报（现为 10 746 行，上限 11 000）。`tools/` 与 `reports/` 脚本合计净增不超过 +150 行 |

## 目标

正式路径不再依赖 `tools` 中的临时补丁：Corella 速率在 900–950 °C 按 ADR 0059 外推，高于 950 °C 拒绝。工况 1、φ = 1 在主线默认设置下复现 T26t 的收敛解，并完成 S20。服务第 1 关。

## 执行要求（给 cursor）

- 新会话，在主仓库目录 `bfb-gasifier/` 内工作，不建 worktree：确认工作树干净后，`git fetch origin && git switch -c task/T26v-corella-extrapolation origin/main`。
- 不在仓库目录以外创建或写入任何文件；临时文件放 `.scratch/`。
- 会话结束时工作树干净，开 PR，同一 PR 内更新 STATUS 的 T26v 一行。
- 动手前先读 ADR 0059、ADR 0010 第 9 条、ADR 0022 C4，以及 `hamel-source-check` skill。

## 范围内

### A. 速率式（`hamel/kinetics/homogeneous.py`）

1. 保留 `_CORELLA_LOW_K`（600 °C）、`_CORELLA_SPLIT_K`（740 °C）、`_CORELLA_HIGH_K`（900 °C，原表上界），新增停止线常数 950 + 273,15 K（命名自定，PR 写明）。
2. `r11_corella_rate` 的判定：
   - T < 600 °C 或 T > 950 °C：拒绝，抛 `KineticsInputError`，消息写明“600–900 °C 为原表区间，ADR 0059 外推至 950 °C 为止”。
   - 600 °C ≤ T ≤ 950 °C：按原两段计算，740 °C 归高温段。
   - 900 °C < T ≤ 950 °C 时，按新增的关键字参数处理：外推（默认）照高温段公式；截断（S20）用 900 °C 的速率常数。参数只接受这两个已批准取值，其它值报错。
3. docstring 改写为现行裁定：`解释性选择：ADR 0010 第 9 条；解释性选择（外推）：ADR 0059`，保留 `Source: Hamel Tab.5.4, S. 97`。

### B. 配置与传递

1. `SolverConfig`（`hamel/cases/config.py`）新增一个字段，只接受两个已批准取值（外推为默认，S20 为截断），做法同 `serio_pair`／`approved_serio`。字段名自定，PR 写明。
2. 经预计算或冻结量传到 `residual.py` 中床区悬浮相的 `r11_corella_rate` 调用（第 441 行附近），做法同 Serio 对的传递。
3. 不改其它速率式；不改自由空间（只用 Serio，ADR 0022 C3）。

### C. 敏感性运行与报告（`tools/`）

1. `tools/campaign.py` 增加 S20：`_ROWS` 一行（“900 °C 以上 Corella 截断”，来源 ADR 0059 第 5 项），`run_specs` 的顺序在 S19b 之后，`overlay` 与配置替换按 S15 的做法。现有各行一字不改。
2. 必报诊断（ADR 0059 第 4 项）：在收敛后的报告中加入床区最高温度（K 与 °C，不取整）、所在格（自上而下的序号与衡算下标）、超出 900 °C 的幅度、超出 900 °C 的格数。放在 `tools/` 的报告路径中（如 `tools/driver_report.py` 或 `tools/run_prep.attach_converged`），PR 写明位置。不改 `driver.py` 对 `tools.run_prep.attach_converged` 的既有调用方式。
3. ADR 0030 维护：补列 F32 与 S20，审查记录加一行“按 ADR 0059 补列 F32 与 S20（本清单只汇总，不作新裁定）”，并把运行总数改为 21 次敏感性运行。
4. 现有历史工具（`t26o_validate`、`t26p_validate`、`o7_checks`、`t26t_validate` 等）把 `_CORELLA_HIGH_K` 置为 +∞ 的写法不改；只在它们的测试因本卡失败时做最小适配，并逐条写进 PR。

### D. 验证

1. **单元测试**（`tests/test_kinetics.py` 等）：
   - 599,9 °C 与 950,1 °C 拒绝；600 °C、740 °C、899,9 °C、900 °C、900,1 °C、949,9 °C、950 °C 接受。
   - 外推：920 °C 的值等于高温段公式；900 °C 处两种模式相等。
   - 截断：920 °C、950 °C 的值都等于 900 °C 的值；900 °C 以下两种模式逐位相等。
   - 非法模式值、非法配置值报错。
   - 原“高于 900 °C 拒绝”的测试改为 950 °C，逐条写进 PR。
2. **床温不超过 900 °C 时逐位不变**：用 `tools/t26q_golden.py` 的现有模式（含 `boundary`）在 `origin/main` 与本分支各跑一次，`np.array_equal` 比较。预期全部相等；凡有不等，逐项列出并说明是否来自床温高于 900 °C 的点（基准中原先被 Corella 拒绝的路径现在可能继续算下去，这类差异是预期的，须逐项指明）。
3. **基准外层（主线默认设置，`tools` 中不打任何 Corella 补丁）**：工况 1、φ = 1，`driver.solve_case`，外层上限 26 轮，墙钟不超过 8 小时。
   - 报告：逐轮外层门四项、内层状态与迭代次数；终态；通过外层门后的 2(b) 复核、元素与能量平衡、壁面散热相对差、底部残差、非负、置零分量；第 C.2 项的必报诊断。
   - **与 T26t D.3 对照**：截取外层通过后、负存量置零之前**首次** `consistency=True` 的 2(b) 复核预计算的全部入参（旧工具保存在 `.scratch/t26t/final_args.pkl` 中的正是这一份，不与置零后的重检混用），与该 pickle 逐位比较（数组用 `np.array_equal`，其它值用 `==`），并对照外层轮数（17）与内层迭代总数（181）。
   - 旧 pickle 中 `kwargs['solver']` 是旧 `SolverConfig`，没有第 B.1 项新增的字段。比较前先给旧配置补入该字段的默认值（外推），再比较；新字段必须等于外推，其余既有字段全部逐位相同。PR 中写明补入方式。
   - 预期逐位相同。不同时停下上报，并给出全程试探点中 Corella 被拒（T > 950 °C）的次数与最高温度。
   - 无论是否相同，都记录全程试探点上最底床格的最高温度与 Corella 拒绝次数。
4. **S20 外层**：同第 3 项，配置改为截断。报告终态、轮数、外层门与收敛后各项检验、必报诊断，以及与基准逐格的床区温度差（K）。S20 未收敛的，如实报告，不调参数重跑。
5. 默认 pytest、`tools/check_citations.py`、`tools/drift_check.py`。已知失败 `test_cli_skip_render_does_not_rewrite_png` 照旧报告。

验证工具写成 `tools/t26v_validate.py`（本卡唯一的新工具），可以 `import` 已合并的 `tools/t26t_validate.py`、`tools/t26q_golden.py` 中的函数，不要复制。

## 范围外（明确禁止）

- 计算、查看或报告验收量（出口组成、T̄_FRB、X_C 等七项）。S20 与基准的对照只用第 D.4 项列出的量；验收量的对照留给第 1 关报告。
- 改 950 °C、改 600 °C 下限、改 740 °C 归属、改其它速率式或自由空间的焦油反应。
- 改收敛门、容差、求解手段（ADR 0054–0058）、ADR 0012 的阈值。
- 改 ADR 0010、0022 的正文（ADR 0059 已加注记），修改其它已接受 ADR 或历史任务卡的正文（第 C.3 项的 ADR 0030 维护除外）。
- 迁移或删除 `driver.py` 对 `tools.run_prep.attach_converged` 的调用；`hamel/` 新增对 `tools/` 的导入。
- 其它敏感性运行（S1–S19b）与第 1 关验收。

## 原文答卷

不另交。原文依据已由 ADR 0059 核对（PDF 58 / S. 96–97，表 5.4）；本卡不新增方程或常数，只实现已裁定的适用范围。PR 中给出 `r11_corella_rate` docstring 的引用与 ADR 0059 的对应。

## 交付物

- 代码改动与测试。
- `tools/t26v_validate.py` 与报告 `reports/t26v-validation/report.md`（运行数据不入库）。
- PR 正文：按文件的行数净变化（`hamel/`、`tools/`、`reports/` 分列）；第 D 项各项结果。
- STATUS 的 T26v 一行；ADR 0030 维护行。

## 停止条件

到达时间上限，或需要越出范围时，停下并在 PR／`docs/STATUS.md` 中写明做到哪里、卡在哪里、需要什么决策。以下情况也停下上报：

- 第 D.3 项与 T26t D.3 不逐位相同；
- 基准外层未收敛，或收敛解的床区任一格高于 950 °C；
- `hamel/` 净增超过 +50 行，或 `tools/` 与 `reports/` 脚本合计净增超过 +150 行；
- 第 D.2 项出现无法归因于床温高于 900 °C 的差异。

S20 未收敛不是停止条件，如实报告即可。

## 审查记录

- PR #171 审查（codex）第 1 轮（`84c2119`）两项 P2：(1) ADR 0059 第 3 项与 ADR 0022 注记把 C4 的“原文适用范围”改成 600–950 °C，与第 2 项“950 °C 不是 Corella 适用范围”冲突；已改为“原文温区仍为 600–900 °C，获批计算区间扩展至 600–950 °C，900–950 °C 属解释性外推”，并写明 S20 截断只是敏感性例外。(2) 第 D.3 项逐位对照未规定新增配置字段与截取时点：旧 pickle 的 `SolverConfig` 没有新字段；已规定补入默认外推值后比较，其余字段逐位相同，截取置零前首次 `consistency=True` 的 2(b) 入参。
- PR #171 复审（codex）第 2 轮（`1a7fea0`）通过。总司令批准本卡（2026-11-19），ADR 0059 接受（合并 PR #171）后可派发。
