# T26w：实现 ADR 0059 修订（950 °C 只判最终收敛解）与路径温度报告

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，在主仓库目录内切任务分支，不建 worktree） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | T26v（PR #172，未合并）；ADR 0059 已接受 |
| 分支 | `task/T26w-corella-converged-check`，以 `task/T26v-corella-extrapolation` 的 `7977437` 为基（ADR 0059 尚未进入 `main`） |
| 时间上限 | 人工半个工作日，另加一次外层运行（≤ 8 h 墙钟）。`hamel/` 净增目标 ≤ +40，超过 +80 停下上报。`tools/` +25 左右 |

## 目标

把 950 °C 从“速率式定义域上限”改判为“最终收敛解的检查线”（ADR 0059 第 2 项，2026-11-23 修订），并实现第 2a 项的路径温度必报。工况 1、φ = 1 在主线默认设置下复现取证运行的收敛解。

## 原因（T26v 的停止）

T26v 的基准在外层第 1 轮第 13 次内层迭代停止：雅可比差分点被 Corella 的 950 °C 上限拒绝。取证（撤掉上限）显示**内层收敛解本身**也两次越过 950 °C：第 1 轮 957,44 °C、第 3 轮 951,72 °C，第 4 轮起 ≤ 913,3 °C，最终 908,945 °C；全程试探点最高 972,98 °C。因此“只把停止线从差分点移开”不成立，必须整体改判为“只判终解”。

## 范围内

1. `hamel/kinetics/homogeneous.py`：删除温度上限判定（只保留 T < 600 °C 拒绝）；把 950 °C 常数改为导出的 `CORELLA_CHECK_K`；docstring 写明“上限不在此处判定”。
2. `hamel/solver/driver.py`：
   - 新增 `converged_bed_peak_reason(state)`：终解床区峰值 > 950 °C 返回理由串，否则空串。
   - `solve_case` 在置零后重检通过之后调用；越界返回 `NOT_CONVERGED`（ADR 0022 C4 的收敛解分支）。
   - 逐轮内层收敛解的床区峰值写入 `diagnostics["bed_peak_per_round"]`；全程试探点（含差分点）最底床格最高温度写入 `diagnostics["trial_bottom_peak_k"]`，只报告。
3. `tools/run_prep.py`、`tools/driver_report.py`：`bed_temperature_peak` 增 `per_round`、`trial_bottom_peak_k`、`check_line_c`；报告新增两行（检查线、逐轮超 900 °C 的轮次）。
4. `docs/adr/0059-…` 第 2、2a、3、4 项与风险、ADR 0010 第 9 条注记、ADR 0022 C4 注记、ADR 0030 F32 与 S20：已在本 PR 内改。
5. `docs/STATUS.md`：T26w 一行；跨环境比较规则写入“已决定（目标所有者）”。

## 范围外

- 不改 600 °C 下限、740 °C 归属、`high_temperature` 的两个取值与 S20 的做法。
- 不改残差、收敛门、容差、求解手段（ADR 0054–0058）与其它速率式。
- 不重跑 S20、不跑 S1–S19b、不计算验收量（出口组成、T̄_FRB、X_C）。
- 不引入功能开关或双轨路径；温度越界只有一处判定（终解）。

## 验证

1. `pytest`：受影响用例已改（950 °C 以上不再拒绝；终解检查线两侧 949 过、950 过、957,44 停）。
2. 主线基准（工况 1、φ = 1，外层上限 26，墙钟 8 h）：终态 `CONVERGED`；终态床区温度场与取证 `.scratch/t26v/evidence.json` 第 16 轮**相对差 ≤ 1×10⁻⁶**（跨环境不要求逐位；规则见 STATUS“已决定”）；逐轮峰值表与取证逐轮吻合（同容差）；报告出现第 1 轮 957,44 °C 与第 3 轮 951,72 °C。
3. `tools/check_citations.py`、`tools/drift_check.py`。

## 停止条件

基准未收敛；终态与取证超容差；`hamel/` 净增 > +80；需要改残差／收敛门／其它速率式。
