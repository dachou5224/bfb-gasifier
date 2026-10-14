# T10：VTT 加压气化炉输入数据补全

| 字段 | 内容 |
|---|---|
| 执行者 | cursor；图上读数由目标所有者人工完成 |
| 审查方 | 总司令（数据录入类，CHARTER §8） |
| 前置任务 | T2（原文索引）；第二阶段需要目标所有者提供原始文献 |
| 分支 | `task/T10-vtt-pressurized-inputs` |
| 时间上限 | 第一阶段 2 个工作日；第二阶段视文献到位情况另计 |
| 状态 | 待开始 |

## 目标

为章程 v0.7 §3 的第 1、2 关（VTT 加压工况 1–4，莱茵褐煤）建立逐项注明出处的输入数据文件。做法与 T7 相同，不同之处是：Hamel 明确说明这批数据取自公开文献，所以可以分两个阶段，先从 Hamel 原文录入，再用 Hamel 引用的原始文献补缺口。

## 范围内

### 第一阶段：只用 Hamel 原文

1. **输入清单**：与 T7 相同，以表 2.1、表 2.2（PDF 16 / S. 13）为准逐项列出，不从旧代码反推。
2. **来源**：
   - 表 7.2（PDF 76 右 / S. 133）：工况 1–8 的运行参数（压力、燃料量、一次风、二次风、蒸汽、λ、循环、燃料分析）；
   - §7.1.3 正文（S. 132–137）：装置几何、二次风位置等；
   - Bild 7.7（PDF 76 左 / S. 132）：装置示意与**模拟采用的 cell 划分**；注意 Hamel 写明旋风尺寸原文献没有给出，模型中假设为两个 cell，这一假设按原文录入并注明出处。
3. **范围**：工况 1–4（褐煤）为必做；工况 5–8（泥炭）顺带录入，因为与 1–4 同一张表，几乎没有额外工作量，但不必为它们补缺口。
4. **格式**：与 T7 完全相同，输出 `data/cases/vtt_pressurized.json`。
   - 每个字段包含 `value`、`unit`（SI）、`source`、`status`；`source` 写成“表号或图号 + `S. <数字>`”（Hamel 印刷页码），`derived` 字段写 `formula`；
   - 键名沿用 `data/cases/wsv400.json` 中已有的同义键名（例如 `fuel_mass_flow`、`primary_air_volume_flow`），不要另起名字；
   - 以 m³(V_N) 给出的流量先只录入体积流量；换算为摩尔流量须等 ADR 0004 获目标所有者批准后再做。
5. **工况标签核对**：`data/validation_cases.json` 收录的 VTT 加压工况编号为 1、6、13、14，而表 7.2 的模拟编号只有 1–8。查明 13、14 的来历（是否来自原始文献的试验编号，或其它表格），写进缺口清单，不修改旧文件。
6. **缺口清单** `data/cases/vtt_pressurized_gaps.md`：格式同 T7，列出 `missing` 字段、模型中的用途、原文查找过的位置、可能的处理；另列与 `validation_cases.json` 的差异。
7. **图上读数清单**：Bild 7.7 中需要读取的几何量与 cell 高度，交目标所有者用工具读取，不得目测估读。

### 第二阶段：用 Hamel 引用的原始文献补缺口

8. 目标所有者如能取得以下文献（`refs/hamel/index.json` 的 `cited_references` 键名），文献 PDF 只放在本地，不入库（CHARTER §11）：
   - `KURKELA_1992`：Kurkela & Ståhlberg, *Air Gasification of Peat, Wood and Brown Coal in a Pressurized Fluidized-Bed Reactor. 1. Carbon Conversion, Gas Yields and Tar Formation*（Hamel 表 7.2 的数据来源）；
   - `KURKELA_1991`：Kurkela, Ståhlberg, Laatikainen, *Pressurized Fluid-bed Gasification Experiments with Wood, Peat and Coal at VTT in 1991–1992. Part 1. Test Facilities …*（正文 S. 133 引作 Kurkela et al. (1993)，年份写法与书目不同，以书目条目为准核对）。
9. 按 CHARTER §4.2，这两篇文献是 Hamel 明确作为 VTT 装置数据来源引用的，可用于补齐**装置与运行数据**（如热值、进气温度、床料、粒径分布、几何尺寸、散热）。不得从中引入 Hamel 没有采用的模型或关联式。
10. 用文献补齐的字段：`status` 为 `source`；`source` 同时写明 Hamel 的引用位置和文献的页码与表号，例如 `Hamel S. 133 引用；Kurkela & Ståhlberg 1992, Table 2, p. 25`。`S. <数字>` 必须是 Hamel 的印刷页码，因为偏离检查以它判定出处。
11. 文献与 Hamel 表 7.2 不一致时，以 Hamel 为准，把差异写进缺口清单，由总司令裁定。

## 范围外

- 不修改 `data/validation_cases.json`、`data/cases/wsv400.json`。
- 不用出口值反推任何输入（CHARTER §1）。
- 不沿用旧代码的默认值填补缺口。
- 不写 `hamel/` 代码。

## 原文答卷

不涉及方程，不需要答卷；每个字段必须有可核对的页码。

## 交付物

- `data/cases/vtt_pressurized.json`
- `data/cases/vtt_pressurized_gaps.md`
- `data/cases/vtt_pressurized_reading_list.md`（图上读数清单），以及目标所有者读数后的回填结果
- 第二阶段完成后，更新上述三个文件

## 验收证据

- 表 2.1、2.2 的每一项都有对应字段（`source` / `derived` / `missing` 之一）。
- 总司令抽查 10 个字段，与原页（以及第二阶段的原始文献）一致。
- `derived` 字段的公式可以复算；对工况 1，用 λ、燃料分析和一次风、二次风做一次与 ADR 0004 相同的相容性检查，结果写进缺口清单（ADR 0004 批准前，只作检查，不用于换算）。
- `python3 tools/drift_check.py` 通过。

## 停止条件

第一阶段时间上限到期，或表 2.1/2.2 中某个输入在原文里的含义不明确时，停下并写进缺口清单。第二阶段在文献到位后开始。
