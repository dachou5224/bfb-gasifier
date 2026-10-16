# T15：VTT 加压装置数据补全（T10 第三阶段）

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，见“执行要求”）；图上读数由目标所有者用工具完成 |
| 审查方 | 交付物（数据录入）：总司令（CHARTER §8）。本任务卡：由总司令起草，按 CHARTER v0.9 §8 由 codex 或目标所有者另开的独立会话审查，总司令不审 |
| 前置任务 | T10 第一、二阶段（#25、#26、#31）；ADR 0005、0006、0011 已批准；ADR 0009 已批准（#53，2026-10-16 合并） |
| 分支 | `task/T15-vtt-plant-data` |
| 时间上限 | 2 个工作日（不含目标所有者读图时间） |
| 状态 | 草稿 |

## 目标

第 1 关（VTT 加压工况 1，莱茵褐煤，不循环）的原文层面裁定已经齐备，剩下的主要缺口是装置与工况数据。本卡把这些缺口逐项查清：能从 Hamel 原文、Hamel 引用的 VTT 文献或已批准的 VTT 文献中查到的，按出处录入；查不到的，写成缺口清单，交总司令裁定。

完成后，`data/cases/vtt_pressurized.json` 中工况 1 装配所需的每个字段要么有出处，要么列为缺口并说明影响。

## 要补的字段

以 `data/cases/vtt_pressurized.json` 中 `status: missing` 的字段和 ADR 0009、0011 新增的工况输入为准。工况 1 必做；工况 2–4 共用的装置字段顺带录入。

| 字段 | 用途 | 线索 |
|---|---|---|
| `cell_count`、`cell_heights` | cell 划分 | Hamel Bild 7.7（PDF 76 左 / S. 132）是装置示意与模拟采用的 cell 划分，由目标所有者读图；读图边界见下文“Bild 7.7 的读图边界” |
| `agent_feed_height` | 二次风高度 | Hamel §7.1.3 只写二次风喷入床层上方的自由空间（S. 134）；具体高度查 VTT 文献 |
| `bed_height` | 操作床高（膨胀床高），ADR 0009 第 3 项的几何存量需要 | VTT 文献；注意 `bed_zone_height` 1,2 m 是床区筒体高度，不是操作床高（ADR 0005） |
| `bed_inventory_composition`、`bed_inventory_diameter`，新增 `bed_particle_density` | 床料组成、平均粒径、颗粒密度 ρ_P（ADR 0009 第 3 项，ADR 0011 第 8、9 条的 Geldart 分类） | VTT 文献 |
| `ash_withdrawal_count`，新增排灰流率与位置 | ADR 0009 第 6 项：排灰流率与所在 cell 为工况输入 | Hamel 表 2.1 只要求位置与数目；VTT 文献 |
| 新增 `distributor_type` | ADR 0011 第 8、9 条：式 (3.11)、(3.14) 按布风器类型选取 | VTT 文献（技术布风器、喷嘴板还是多孔板） |
| 新增 `solid_feed_temperature` | ADR 0009 第 4 项：惰性补料的加入温度；燃料加入温度 | VTT 文献；查不到列为缺口 |
| `fuel_size_class_mass_fractions`（工况 1） | 燃料初始粒径分布（表 2.2） | VTT 文献的筛分分析；KURKELA_1992 p. 3 只给“筛到 4 mm 以下” |
| `inlet_gas_mole_fractions`（工况 1） | 进气组成 | 原文与文献若只写“空气”，列为缺口，由总司令裁定（不得自行取 21 %） |
| 吹扫氮温度（新增，例如 `nitrogen_purge_temperature`） | ADR 0006 第 3 项：工况 1、2 计入的料斗吹扫氮须有温度，能量衡算需要；ADR 0006 规定装配时为 missing、不得取默认值 | VTT 文献（料斗加压与给料系统的说明）；查不到列为缺口 |
| `wall_heat_transfer_coefficient` | 壁面传热（表 2.1） | VTT 文献若只有热损失（KURKELA_1992 表 6 为能量平衡差额），记录该数并列为缺口 |
| `cyclone_geometry`、`connecting_line_geometry`、`return_line_geometry` | 旋风与连接管 | Hamel S. 133 写明文献未给旋风尺寸，模型中按两个 cell 处理；工况 1 不循环。按原文录入这一点，其余列为缺口 |

另外，若 VTT 文献给出床层压降，一并录入（`status: cited`），注明测点位置。它可以作为床料质量的旁证，但本卡不据此推算床高。

## Bild 7.7 的读图边界

Bild 7.7 的图注是 „Schematische Darstellung … mit der bei der Simulation verwendeten Zelleneinteilung“，即**示意图**：图中有 cell 分界，但没有轴向刻度、尺寸线或比例尺。工具能精确读出像素，不等于图是按真实比例绘制的。因此：

- **可以从图上录入的**：cell 个数，各 cell 分界、燃料口、进气口、二次风口、床面与自由空间、旋风之间的**拓扑与先后关系**（例如某进气口位于第几个 cell）。录为 `source`，出处为 Bild 7.7。
- **不可以仅凭图上像素录入的**：有单位的高度（cell 高度、进料高度、二次风高度、操作床高）。不得用总高 3,5 m、床区筒体 1,2 m 或自由空间 2,3 m 去标定整张示意图，再把像素比例换成高度。
- 有单位的高度只能来自原文标注，或来自同一装置的几何尺寸加上可证明的定位关系（例如某 cell 分界正好与床区和自由空间的交界重合，且原文或图中写明）。无法证明时，这些字段保持 `missing`，cell 高度的划分方式另交总司令作解释性裁定。
- 目标所有者的工具读数照常执行。读图清单与回填结果须把“图上读得的像素与计数”和“若按比例换算所需的额外假定”分开记录；后者只作参考，不进入 `source` 或 `derived` 字段。

## 来源与顺序

按以下顺序查，前一级查到的不必再查后一级；不同来源的数值不一致时，全部记录，交总司令裁定。

1. **Hamel 原文**：§7.1.3（S. 132–137）与 Bild 7.7。cell 划分与图上尺寸由 cursor 准备读图清单（要读什么、在图的什么位置），**由目标所有者用工具读取**，不得目测估读（T7 规则）。
2. **KURKELA_1992**（`docs/kurkela1992.pdf`，Hamel 表 7.2 的数据来源，CHARTER §4.2）：T10 第二阶段已查过；本卡只补查上表中尚未查过的字段，并写明查过的页码。
3. **KURKELA_1991**（VTT Publications 161，Part 1：*Test Facilities and Gasification Experiments with Sawdust*，本地 `docs/VTT161_part1.pdf`，扫描件，无文本层，SHA256 前缀 6cc286c8a9bb）：
   - **先核实 Hamel 在正文何处、为何引用它**：`refs/hamel/index.json` 中该条的引用位置是按作者与年份自动匹配的，可能与 KURKELA_1992 混淆。逐处查看页面图像，写明 Hamel 引用它的页码和用途，据此确定 CHARTER §4.2 的使用范围。
   - **再核实装置是否相同**：这份报告记录的是 1991–1992 年的试验，Hamel 表 7.2 的试验来自 1988–1990 年（KURKELA_1992）。逐项对照两份文献的装置尺寸（床区与自由空间直径和高度、燃料口高度、布风器、旋风），写明一致与不一致之处。只有能确认适用于 1988–1990 年试验装置的数值才可录入；装置改动过的，不录入，写进缺口。
   - 页面图像：用 PyMuPDF 渲染到主工作树的 `refs/vtt161/pages/`（每页 PNG，三位数页号），`.gitignore` 增加一行 `refs/vtt161/pages`。扫描件没有文本层，引用以页面图像为准。
4. **VTT Publications 249 Part 2**（`docs/VTTPubl.249_Part2.pdf`）：目标所有者 2026-10-14 已批准，**范围限定**为：只用于补齐 VTT 加压装置的装置与运行数据，只有逐字段确认与 Hamel 表 7.2 为同一次试验的值才可采用（见 STATUS“已决定”）。这份报告是 1991–1994 年的试验，大概率对不上表 7.2；若对不上，写明“不适用”，不录入。

不得打开 `docs/kurkela1992_part2.pdf` 和 `docs/VTT_part3.pdf`：这两份 Hamel 未引用，目标所有者未批准使用。

## 录入格式

- 沿用 `data/cases/vtt_pressurized.json` 的格式与键名（与 T7、T10 相同）：每个字段含 `value`、`unit`（SI）、`source`、`status`；`cited` 字段写 `cited_source`（文献与页码）；`derived` 字段写 `formula`；`missing` 字段写查找过的位置。
- 新增键名尽量沿用 GLOSSARY 与 ADR 中的名称；新增的物理参数配置字段须符合 CHARTER §9 的出处要求。
- 所有 `cited` 字段在模型使用前须由总司令成组登记 ADR（同 ADR 0005 的做法）；本卡只录入，不作模型输入裁定。
- 更新 `data/cases/vtt_pressurized_gaps.md`：每个仍缺的字段写明模型何处用到、查过哪些来源、可能的处理方式。

## 执行要求（给 cursor）

- **新会话**：开一个新的 cursor 会话，不要沿用做过 T4、T8、T9、T11–T14 或 worktree 清理的会话。
- **新 worktree**：
  ```
  git worktree add ../bfb-T15 -b task/T15-vtt-plant-data origin/main
  ```
  按 `.claude/skills/hamel-source-check/SKILL.md` 第一节第 2 项链接 Hamel 页面图像。
- **开工前必读**：`docs/CHARTER.md`（§4 来源层级、§9 配置字段出处）、本卡、`docs/tasks/T10-vtt-pressurized-inputs.md`、`data/cases/vtt_pressurized_gaps.md`、`docs/adr/0005-vtt-pressurized-cited-fields.md`、`docs/adr/0009-bed-solids-transport-inventory.md`（第 3、4、6 项）、`docs/adr/0011-balances-hydrodynamics-source-gaps.md`（第 8、9 条）、`docs/STATUS.md` 中关于 VTT 文献的已决定条目。
- **读图**：字形不清时，不要把推断的字符当成录入值；先按量纲和物理定义判断并写明理由，仍无法判断的，记下文献、页、位置，交目标所有者判读。
- 会话结束时工作树干净，开 PR，审查方填总司令。

## 范围外（明确禁止）

- 不写 `hamel/` 下的任何代码。
- 不用出口值反推任何输入（CHARTER §1）；不用床层压降推算床高。
- 不沿用旧代码或旧 `validation_cases.json` 的默认值填补缺口。
- 不自行取进气组成（如 21 % O₂）、壁面传热系数或床料物性的“典型值”。
- 不打开未获批准的文献。

## 交付物

- `data/cases/vtt_pressurized.json` 的更新
- `data/cases/vtt_pressurized_gaps.md` 的更新
- Bild 7.7 的读图清单，以及目标所有者读数后的回填
- KURKELA_1991 的引用位置核实与装置一致性对照（写进 gaps 文件的一节）
- `.gitignore` 的一行；`docs/STATUS.md` 中 T15 一行的更新

## 验收证据

- 上表每个字段都有结论：录入（附出处页码）或缺口（附查过的位置）。
- KURKELA_1991 的使用写明 Hamel 的引用位置和装置一致性对照结果。
- 总司令抽查 10 个字段，与原页一致。
- `python3 tools/drift_check.py` 通过。

## 停止条件

到达时间上限，或发现 VTT 文献中的装置与 1988–1990 年试验装置不一致而无法判断某字段是否适用时，停下，写进缺口清单。“查不到”本身是有效结论。
