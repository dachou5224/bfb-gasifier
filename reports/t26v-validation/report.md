# T26v 验证报告：Corella 高温段外推至 950 °C（ADR 0059）与 S20

本报告不记录出口组成、T̄_FRB、X_C 或其它验收量。全部运行用 `data/cases/vtt_pressurized.json` 的工况 1、φ = 1。

## 结论

| 项 | 结果 |
|---|---|
| A 速率式（600–950 °C，950 °C 停止线，外推／截断） | 已实现 |
| B 配置与传递（`SolverConfig.corella_high_temperature`，经冻结量到 R11） | 已实现 |
| C 敏感性 S20、必报诊断、ADR 0030 维护 | 已实现 |
| D.1 单元测试 | 通过（pytest 515 passed、1 既有失败、5 deselected，78 s） |
| D.2 与 `origin/main` 逐位比较（含 boundary） | 全部数值逐位相同；唯一差异是新增的元数据字符串 |
| D.3 基准外层（主线默认设置，不打补丁） | **未收敛**：第 1 轮第 13 次内层迭代因雅可比差分点越过 950 °C 停止线而停止。未到 2(b)，无 17 轮／181 次 |
| D.3 诊断对照（仅诊断：撤掉 950 °C 停止线） | 17 轮、181 次内层迭代收敛；2(b) 入参与 T26t D.3 的 pickle **逐位相同** |
| D.4 S20（900 °C 截断） | **未收敛**：同一停止线，第 1 轮第 14 次内层迭代停止；没有可比的收敛温度场 |

**结论**：ADR 0059 的 950 °C 停止线在本工况的主线路径上会被触发：求解在第 1 轮第 13 次迭代时要评价一个越过 950 °C 的雅可比差分点。撤掉停止线的诊断运行全程用到 973 °C 的取值，且 2(b) 入参与 T26t D.3 逐位相同，说明差异只来自停止线，不是本卡的代码改动。是否让停止线对差分点生效，超出本卡范围，交总司令／目标所有者裁定。

## A、B、C 的实现

- `hamel/kinetics/homogeneous.py`：新增 `_CORELLA_STOP_K = 950 + 273,15 K` 与 `_CORELLA_HIGH_TEMPERATURE = ("extrapolate", "truncate")`。`r11_corella_rate` 新增关键字参数 `high_temperature`（默认 `extrapolate`）；600–950 °C 外接受理，950 °C 以上抛 `KineticsInputError`（消息含“600–900 °C 为原表区间，ADR 0059 外推至 950 °C 为止”）；`truncate` 时 900 °C 以上用 `_CORELLA_HIGH_K` 的速率常数。docstring 含 `解释性选择：ADR 0010 第 9 条；解释性选择（外推）：ADR 0059` 与 `Source: Hamel Tab.5.4, S. 97`。
- `hamel/cases/config.py`：新增 `corella_high_temperature: str = "extrapolate"`，`__post_init__` 拒绝其它取值；`approved_corella_high_temperature` 同 `approved_serio` 的做法。
- `hamel/solver/frozen.py`、`precalc.py`：冻结量新增 `corella_high_temperature`，由预计算按 `approved_corella_high_temperature` 写入。
- `hamel/solver/residual.py`：床区悬浮相的 R11 调用传入 `high_temperature=self.frozen.corella_high_temperature`。自由空间的 R11（Serio）不变。
- `tools/campaign.py`：`_ROWS` 增 `S20`（“900 °C 以上 Corella 截断”，ADR 0059 第 5 项），顺序在 S19b 之后，`overlay` 不改工况副本，`config_for` 取 `truncate`。现有各行未改。
- `tools/run_prep.py`：`bed_temperature_peak`（ADR 0059 第 4 项）写进 `attach_converged` 的 `diagnostics["bed_temperature_peak"]`：峰值 K 与 °C（不取整）、自上而下序号（1 起）、衡算下标（0 为床顶）、超出 900 °C 的幅度与格数。
- `tools/driver_report.py`：Markdown 报告增“床区最高温度（ADR 0059 第 4 项）”一节。
- `docs/adr/0030-stage1-frozen-list.md`：补列 F32 与 S20，运行总数改为 21 次敏感性运行，审查记录加一行（只汇总，不作新裁定）。

## D.1 单元测试

- `tests/test_kinetics.py`：599,9 °C 与 950,1 °C 拒绝；600、740、899,9、900、900,1、920、949,9、950 °C 接受且等于高温段公式；截断在 900,1–950 °C 等于 900 °C 的值，900 °C 以下两种模式逐位相同；非法模式值报错。原“高于 900 °C 拒绝”的断言改为 950 °C。
- `tests/test_sensitivity_params.py`：`config_for("S20")` 为 `truncate`、基准为 `extrapolate`、S20 在运行表末位；`approved_corella_high_temperature` 与 `SolverConfig` 非法取值报错。
- `tests/test_hamel_skeleton.py`：字段登记加 `corella_high_temperature`。
- `tests/test_driver.py`：运行表断言加入 S20，去重后的改动数 21 → 22。

## D.2 与 `origin/main` 的逐位比较

`origin/main`（`1605069`）导出到 `.scratch/t26v/main/` 后运行，两份都用 `tools/t26q_golden.py` 的现有模式：

| 基准 | main 数组数 | 本分支数组数 | missing | extra | different |
|---|---|---|---|---|---|
| `capture`（φ = 0,3 第 1 轮前 3 次内层迭代 + φ = 1 初值构造拒绝终态） | 418 | 419 | 空 | `data.phi03.first.corella_high_temperature` | 空 |
| `boundary`（T26t 第 17 轮与 D.3 的 2(b) 两份预计算入参） | 160 | 160 | 空 | 空 | 空 |

唯一多出的是本卡新增的冻结量元数据字符串，值为 `extrapolate`，不是数值。全部数值数组 `np.array_equal` 相同，说明床温不超过 950 °C 的路径逐位不变。

## D.3 基准外层（主线默认设置）

`python -m tools.t26v_validate baseline --limit-s 28800`，`SolverConfig(outer_round_limit=26)`，不打任何 Corella 补丁。

| 项 | 值 |
|---|---|
| 终态 | `NOT_CONVERGED` |
| 停止位置 | 第 1 轮、第 13 次内层迭代 |
| 原因 | 雅可比不能求解：差分点不可评价（差分步长减半后仍被 Corella 拒绝，切换前） |
| 墙钟 | 2 227,5 s |
| 试探点最底床格最高温度 | 1 238,057 928 464 187 9 K = 964,907 928 464 187 9 °C |
| 全程 Corella 拒绝次数（T > 950 °C） | 82（含初值构造的试探点；被拒点的最高温度同上） |
| 2(b) 复核／元素与能量平衡／壁面散热／底部残差／非负／置零分量 | 未到；无收敛解 |

未生成 `final_args_baseline.pkl`，故第 D.3 项的逐位比较不成立。

## D.3 诊断对照（不是基准，不进入主线）

为把“实现有误”与“停止线生效”分开，只在诊断进程里把 `homogeneous._CORELLA_STOP_K` 置为 `+inf`，其余与基准完全相同：

| 项 | 值 |
|---|---|
| 终态 | `CONVERGED` |
| 轮数／内层迭代 | 17／181 |
| 墙钟 | 8 193,1 s |
| 试探点最底床格最高温度 | 1 246,128 712 991 987 9 K = 972,978 712 991 987 9 °C |
| 全程 Corella 拒绝次数 | 0 |
| 块残差 | 气相 1,884 675 629·10⁻⁷、固相 3,331 787 029·10⁻⁷、能量 2,441 885 623·10⁻⁸ |
| 元素平衡 C／H／O | 1,658 629 964·10⁻⁸／7,053 670 226·10⁻⁹／4,742 105 107·10⁻⁹ |
| 能量平衡 | 8,878 623 911·10⁻¹¹ |
| 壁面散热 | 21 700,0 W，目标 21 700,0 W，相对差 0 |
| 底部残差 | 1,639 087 793·10⁻⁷ |
| 非负 | 通过 |
| 置零分量 | 2 个：`freeboard[b=3,geo=4].char.free.k4` = −1,149 146 378·10⁻¹⁸⁹ kg、`freeboard[b=3,geo=4].ash.free.k4` = −1,963 506 675·10⁻¹⁹³ kg；低于阈值的分量 0 |
| 床区最高温度（ADR 0059 第 4 项） | 1 182,095 235 938 270 7 K = 908,945 235 938 270 7 °C，自上而下第 8 格（衡算下标 7），超出 900 °C 8,945 235 938 270 57 K，共 1 格 |

与 T26t D.3 的对照（`_differences` 递归比较，数组 `np.array_equal`，其它值 `==`）：旧 pickle 的 `SolverConfig` 缺 `corella_high_temperature`，按字段名取 `vars()` 后补入默认 `extrapolate`；比较结果 `equal = true`、差异列表为空。外层轮数 17、内层迭代总数 181 也相同。

## D.4 S20（900 °C 截断）

`python -m tools.t26v_validate s20 --limit-s 28800`，配置只改 `corella_high_temperature="truncate"`。

| 项 | 值 |
|---|---|
| 终态 | `NOT_CONVERGED`（原因同基准：差分点越过 950 °C 停止线，第 1 轮第 14 次内层迭代） |
| 墙钟 | 2 246,3 s |
| 试探点最底床格最高温度 | 1 238,146 089 576 089 K = 964,996 089 576 089 °C |
| 全程 Corella 拒绝次数 | 108 |
| 与基准逐格床区温度差 | 无（两者都没有收敛的床区温度场） |

S20 的“未收敛”由停止线造成，不是截断的导数拐点造成的；两者只要不越过 950 °C，速率值逐位相同（D.1 已固定）。按要求没有调参重跑。

## 行数

| 目录 | 增 | 删 | 净 |
|---|---|---|---|
| `hamel/` | 35 | 7 | **+28**（目标 ≤ +30，停止线 +50） |
| `tools/`（已跟踪文件） | 41 | 4 | +37 |
| `tools/t26v_validate.py`（新） | 181 | 0 | +181 |
| `tools/` 合计 | — | — | **+218**（卡上停止线 +150，**超出 68 行**） |

`tools/` 超限来自验证入口本身：两次外层运行、试探点温度与 Corella 拒绝计数、逐轮外层门、旧 pickle 的字段补全与逐位比较。`tools/t26t_validate.py` 的 `outer()` 会在 `tools` 里把 `_CORELLA_HIGH_K` 置为 `+∞`（本卡禁止的补丁），且无法运行 S20，故不能复用。没有为压行删去卡上要求的报告项。

## 检查

- 默认 pytest：515 passed、1 failed（既有的 `test_cli_skip_render_does_not_rewrite_png`，本地缺 `pdfplumber`）、5 deselected，78 s。
- `python tools/check_citations.py`：通过。
- `python tools/drift_check.py`：提交后只剩既有的“两周无工况指标变化”；新增求解器字段 `corella_high_temperature` 由 ADR 0059 第 5 项批准，字段名登记在 ADR 0030 的 F32 行（本卡授权的维护）。

## 需要裁定

1. 950 °C 停止线是否应约束**雅可比差分点**。ADR 0022 C4 的原文是对“试探点”判定不可行；本卡的实现与既有的差分点处理（T26o 的 O4）都把不可评价的差分点记为求解停止。诊断显示该路径需要 973 °C 的取值，撤掉停止线后与 T26t 逐位相同。
2. 若维持现状，第 1 关的主线运行会在第 1 轮停止，等于把工况 1 判为“超出动力学参数适用范围”（CHARTER §2.2 的验收无法开始）。
