# VTT 加压输入缺口与编号核对

对照范围：表 2.1、表 2.2（S. 13）的每一项都在 `vtt_pressurized.json` 里有字段。Hamel 给出的数取自表 7.2（PDF 76 右 / S. 133，以及 PDF 77 左 / S. 134 的续表）和 §7.1.3（S. 132–137）。文献补上的数 `status` 为 `cited`，`source` 写 Hamel 在 S. 133 引用 Kurkela & Ståhlberg（1992）的位置。旧文件没有改。这些 `cited` 字段已由 ADR 0005 成组批准。其中五个字段只作记录，T9 不得装配。料斗吹扫氮见 ADR 0006。

## 表 7.2 的范围

S. 133 是工况 1–8：1–4 为莱茵褐煤，5–8 为芬兰泥炭。S. 134 的「Fortsetzung Tabelle 7.2」继续列出工况 9–15：9–13 为芬兰泥炭，14–15 为锯末（表内 Sägemehl，行末换行印成 Säge-mehl）。

任务卡里「表 7.2 的模拟编号只有 1–8」与原页不符。9–15 仍是同一张表，按任务卡对 5–8 的做法随表录入，不为它们单独补缺口。共享的装置缺口对全部工况同样成立。

褐煤的 C、H、O、N、S 与灰分之和为 100,0 Gew.-%。锯末为 100,0 Gew.-%。泥炭为 100,05 Gew.-%（S 印到 0,25，其余多为一位小数）。这 0,05 个百分点留在原值上，没有改动。

工况 5 和 8 的 Dampf 是短横，记为 0。其余工况的蒸汽都有数字。

## 编号 13、14 的来历

`validation_cases.json` 里的 `CASE_VTT_PRESSURIZED_PEAT_13` 和 `CASE_VTT_PRESSURIZED_SAWDUST_14` 就是表 7.2 续表的 Simulation Nr. 13 和 14，不是另一套文献试验号。

| 旧文件编号 | 表 7.2 的列 | 印刷页 | 燃料 | 循环 |
|---|---|---|---|---|
| 1 | 1 | S. 133 | Rhein. Braunk. | Nein |
| 6 | 6 | S. 133 | Finn. Torf | Nein |
| 13 | 13 | S. 134 | Finn. Torf | Ja |
| 14 | 14 | S. 134 | Sägemehl | Nein |

旧文件把这四列的出处都写成 Table 7.2 (p.133)。1 和 6 在 S. 133；13 和 14 在续表 S. 134。Hamel 用的就是 Simulation Nr.，Bild 7.10 的图注也写符号里的号码是表 7.2 的 Versuchsnummer。

第一阶段不能判断 Kurkela 原文的试验编号是否与这列号相同。那是第二阶段、拿到文献之后的事。

## 第二阶段补上的数

Hamel 在 S. 133 写模拟所用的数据取自 Kurkela & Ståhlberg（1992）。这篇就是 *Fuel Processing Technology* 31 (1992) 1–21。表 4、表 5 的压力、燃料量、水分、一次风、二次风、蒸汽、ER 和循环，与表 7.2 的工况 1–15 逐列相同。表 7.2 已有的数没有改。

§7.1.3（S. 134）的三句位置说明现在是 `shared` 里的字段，`status: source`：

| 字段 | 值 |
|---|---|
| `steam_feed_location` | 所有工况的蒸汽都和一次风一起从布风板送入 |
| `secondary_air_location` | 二次风喷入床层上方的自由空间 |
| `fuel_feed_location` | 燃料在紧贴布风板上方加入 |

文献给出的高度和装置尺寸如下。对得上表 7.2 的才写入。

| 字段 | 值 | 文献 |
|---|---|---|
| `fuel_feed_height` | 0,12 m | p. 5–6，分配器以上的进料口。两台螺旋进这一个口 |
| `bed_zone_height` | 1,2 m | 表 3，床区筒体高度。不是操作床高 |
| `freeboard_height` | 2,3 m | 表 3 |
| `distributor_to_outlet_height` | 3,7 m | p. 4，分配器到气体出口管 |
| `pressure_vessel_height` | 5,6 m | p. 4，压力容器外壳 |
| `pressure_vessel_inside_diameter` | 0,6 m | p. 4 |
| `thermal_design_power` | 400 kW | 表 3，用泥炭时的最大热输入 |
| `ash_sample_pipe_diameter` | 35 mm | p. 6，底部床料取样管 |
| 泥炭 `particle_size_range` | 0–4 mm | p. 3 |
| 褐煤 `particle_size_max` | 4 mm | p. 3，只给上限 |
| 三种燃料的固定碳（干基） | 27,4 / 16,8 / 42,7 Gew.-% | 表 1 平均值 |

`distributor_to_outlet_height`、`pressure_vessel_height`、`pressure_vessel_inside_diameter`、`ash_sample_pipe_diameter`、`thermal_design_power` 的 note 写明：只作记录，不作模型输入（ADR 0005）。`reactor_height` 现为 3,5 m，`status: derived`，公式是 bed_zone_height + freeboard_height，与 S. 133 的 Gesamtanlagenhöhe 3,5 m 一致（ADR 0005）。`plant_height` 仍按原词记 3,5 m。

空气和蒸汽温度、收到基低位热值只在表 6 的六个代表点里。其中五列对得上表 7.2，温度都是 185 °C（458,15 K）：

| 工况 | 对上的条件 | 热值（收到基） | 表 6 的 ER |
|---|---|---|---|
| 1 | 褐煤，不循环，蒸汽/无灰干基 0,17，燃料 46,9 kg/h | 21,0 MJ/kg | 0,32；表 5 和表 7.2 是 0,33，λ 不改 |
| 3 | 褐煤，循环，ER 0,44，蒸汽比 0,20，燃料 41,4 kg/h | 20,9 MJ/kg | 0,44 |
| 7 | 泥炭，不循环，ER 0,32，蒸汽比 0,11，燃料 61,6 kg/h | 17,3 MJ/kg | 0,32 |
| 11 | 泥炭，循环，ER 0,40，蒸汽比 0,18，燃料 45,0 kg/h | 17,3 MJ/kg | 0,40 |
| 14 | 锯末，不循环，蒸汽比 0,08，燃料 42,8 kg/h | 16,8 MJ/kg | 0,27；表 5 和表 7.2 是 0,28，λ 不改 |

表 6 还有一列锯末：ER 0,38，蒸汽比 0,38，热值 17,6 MJ/kg。表 7.2 的工况 15 是 ER 0,39、蒸汽比 0,34，对不上，没有采用。其余工况的进气温度和热值仍是 `missing`。

文献与表 7.2 不一致、因而没有改表 7.2 的地方：

- 泥炭氧（差减）表 1 是 33,55 wt%，表 7.2 是 33,6。按 33,55 与 C、H、N、S、灰分相加为 100,00；按 33,6 为 100,05。元素仍用 33,6。
- 标准立方米在表 3、表 4 的注里定义为 273,15 K、101,3 kPa。摩尔流量仍按 ADR 0004 的 101 325 Pa，不改。
- 表 4、表 5 另有石灰石进料和氮气进料，表 7.2 没有这两路。氮气已写入各工况的 `nitrogen_purge_volume_flow`（cited）和按 ADR 0004 换算的 `nitrogen_purge_amount_flow`（derived）。泥炭九列（表 4，p. 8，工况 5–13）是 8,8、3,6、3,7、2,3、4,2、3,9、3,2、2,7、2,4 Nm³/h。褐煤和锯末（表 5，p. 9，工况 1–4、14–15）是 6,2、3,2、3,5、2,5、2,0、3,5 Nm³/h。只有工况 1、2 的 note 写按 ADR 0006 暂定计入模型；其余写模型是否计入待裁定（ADR 0006 第 4 项）。吹扫氮的温度仍缺，没有填。石灰石是 1,8、0、0、1,1、0、0、0、0、0,3 kg/h（泥炭）和 0（褐煤、锯末），不在 ADR 0006 范围内，没有写入。
- 床层温度、自由空间温度、床层压降和热损失是试验结果，不是进料。表 6 的热损失是差额，没有换成传热系数。
- 文献写两台旋风。Hamel 在 S. 133 按一个旋风假设两个 cell，`cyclone_cell_count` 仍是 2。旋风尺寸两边都没有。
- 床料成分没有给出。p. 16 的石英砂是建议另做的试验，不是这批试验的床料。

`bed_height` 仍空着。表 3 的 1,2 m 是床区筒体。各列床层压降从 0 到 4,3 kPa，不是同一个操作床高。

## 清单对应

| 表 2.1 / 2.2 的项 | 字段 |
|---|---|
| 流化床直径 | `bed_diameter` = 0,15 m |
| 自由空间直径 | `freeboard_diameter` = 0,25 m |
| 反应器高度 | `reactor_height` = 3,5 m（derived，ADR 0005）；正文的 3,5 m 仍记在 `plant_height` |
| 连接管、返料管尺寸 | `connecting_line_geometry`、`return_line_geometry` |
| 旋风尺寸 | `cyclone_geometry` 为 missing；两格假设记在 `cyclone_cell_count` = 2 |
| 燃料进口个数与高度 | `fuel_feed_count` = 1；`fuel_feed_location` 为紧贴布风板上方；`fuel_feed_height` = 0,12 m（cited） |
| 气化剂进口个数与高度 | `agent_feed_count` = 2；`agent_feed_height` 为 missing |
| 排灰口 | `ash_withdrawal_count` |
| 返料位置 | `return_location`；各工况开关是 `recirculation` |
| 壁面传热系数 | `wall_heat_transfer_coefficient` |
| 燃料量 | `fuel_mass_flow` |
| 工业分析、元素分析、热值 | `fuels` 与各工况水分；固定碳为 cited。热值只在工况 1、3、7、11、14 |
| 燃料粒径分布 | `fuel_size_class_mass_fractions` |
| 各进气口体积流量 | `primary_air_volume_flow`、`secondary_air_volume_flow`，摩尔流量按 ADR 0004。吹扫氮另记 `nitrogen_purge_volume_flow` |
| 进入气体组成 | `inlet_gas_mole_fractions` |
| 进入气体温度与压力 | 工况 1、3、7、11、14 的 `inlet_temperature` 为 185 °C；其余仍 missing。`pressure` 取表 7.2 的 Druck |
| 床料组成与粒径 | `bed_inventory_composition`、`bed_inventory_diameter` |
| 床高 | `bed_height` |
| 蒸汽 | `steam_mass_flow`（表 7.2 的 Dampf；与一次风一起从 Düsenboden 进入） |

`secondary_air_amount_flow` 与已有的 `primary_air_amount_flow` 同一换算。WSV400 的二次风是短横，所以那个文件没有这个键。

## 缺口

| 字段 | 模型里用在哪里 | 原文查过的位置 | 可能的处理 |
|---|---|---|---|
| `heating_value` | 能量衡算、进料焓 | 表 6 只覆盖工况 1、3、7、11、14 | 其余工况对不上代表点，不外推。不能用旧文件里的数 |
| `inlet_gas_mole_fractions` | 进气组成 | 表 2.2；表 7.2 只写 Luft/Dampf；1992 年文献也没有空气的 y_i | 相容性检查里的 0,21 不是进气组成。吹扫氮已按 ADR 0006 另记，没有并进空气 |
| `inlet_temperature` | 进气焓 | 表 6 的 185 °C 只用于对上的五列 | 其余工况不外推。旧文件的 293 K 没有采用 |
| `bed_inventory_composition`、`bed_inventory_diameter` | 床料衡算与流动 | 表 2.2；§7.1.3 | 原文未给。不能用旧文件的 1–4 mm |
| `bed_height` | 表 2.2 的操作床高 | 表 3 的 1,2 m 是床区筒体，已记入 `bed_zone_height` | 见读数清单。不把 1,2 m 当作 Wirbelbetthöhe |
| `cell_count`、`cell_heights`、`agent_feed_height` | cell 划分与二次风高度 | Bild 7.7；文献只说二次风在床层上方 | 见读数清单。燃料口高度已由文献记为 0,12 m。旧文件的 800 mm 没有采用 |
| `wall_heat_transfer_coefficient` | 壁面散热 | 表 2.1；§7.1.3 | 原文未给系数 |
| `connecting_line_geometry`、`cyclone_geometry`、`return_line_geometry` | 表 2.1 的下游设备 | S. 133 写明旋风实际尺寸没有给出；图上无尺寸线 | 旋风按原文假设为两个 cell，已记入 `cyclone_cell_count`。尺寸仍缺。第二阶段若文献里其实有尺寸，与 S. 133 这句话不一致，把差异写在这里，由总司令裁定 |
| `ash_withdrawal_count` | 排灰位置 | §7.1.3 未写；Bild 7.7 待看 | 见读数清单 |
| `fuel_size_class_mass_fractions` | 燃料粒径分级 | KURKELA_1992, Fig. 1 是累积曲线，没有质量分数表 | 筛分上限已记入燃料字段。曲线不目测 |

一次风、二次风和吹扫氮的摩尔流已按 ADR 0004 换算，V_M = 0,022414 m³/mol。它们不再列入上表。`reactor_height` 已按 ADR 0005 记为 3,5 m，不再列入上表。吹扫氮的温度仍缺。

S. 135–137 查过，没有另外的输入数。S. 135 起是 Bild 7.8，S. 136–137 是 Bild 7.9 和 Bild 7.10，都是结果图。S. 135 写 Kurkela 与 Ståhlberg 没有给出床内和自由空间测温点的位置，这是测点，不是进气温度。

## 工况 1 的空气量相容性

干燃料量与理论需氧量按任务卡的式子，原子量 M_C = 12,011、M_H = 1,008、M_O = 15,999、M_S = 32,06 kg/kmol。元素质量分数用表 7.2 褐煤的 wf 值。λ_O2 = 0,33。空气含氧 0,21 只用于这次核算。

- m_dry = 46,9 × (1 − 0,118) = 41,3658 kg/h
- n_O2,st = 41,3658 × (0,638/12,011 + 0,046/(4×1,008) + 0,003/32,06 − 0,262/(2×15,999)) = 2,334366 kmol/h
- 空气量 = 2,334366 × 0,33 / 0,21 × 22,414 = 82,221 m³(V_N)/h
- 表 7.2 一次风加二次风 = 51,8 + 30,1 = 81,9 m³(V_N)/h
- 相对差 (82,221 − 81,9) / 81,9 = 0,392%

低于任务卡举例的 2%。输入没有改。差值与 λ 只印到 0,01、两股风量各印到 0,1 m³/h 同一量级。

同一节还有三句约数，用来核对读数，不代替表值：

- 一次风份额 51,8/81,9 = 0,632。S. 134 写工况 1 约 60% 的空气从 Düsenboden 进入。
- 床内空气比 0,33 × 0,632 = 0,209。S. 134 写床内约 λ_O2 = 0,2。
- 蒸汽与干燃料之比：工况 1 为 6,8/41,3658 = 0,164，正文写 0,16；工况 6 为 15,5 / (65,5×(1−0,148)) = 15,5/55,806 = 0,278，正文写 0,28。S. 135 写的就是对无水燃料量的这个比。

## 正文与图框不一致的地方

没有取舍，两处都保留。

| 位置 | 原文 | 表 7.2 |
|---|---|---|
| Bild 7.7 图框，S. 132 | Betriebsdruck = 0,5–1,0 MPa | 工况 12、14、15 为 0,4 MPa。正文 S. 132 只写「压力直到 1,0 MPa」，没有下限 |
| Bild 7.7 图框，S. 132 | Durchsatz ≈ 25–75 kg Brennstoff/h | 工况燃料量从 26,5 到 76,3 kg/h。这是图框里的装置范围，各工况仍用表中的 ṁ_B |
| S. 133 | Gesamtanlagenhöhe 3,5 m | 没有另一条叫 Reaktorhöhe 的数。ADR 0005 将 `reactor_height` 取为床区 1,2 m 与自由空间 2,3 m 之和，与这 3,5 m 一致 |

图框里的直径 0,15 m / 0,25 m、气化剂 Luft/Dampf，与 S. 132 正文一致。

## `docs/VTTPubl.249_Part2.pdf` 没有用来填数

这份文件是 VTT Publications 249（Espoo 1995）：*Pressurised fluidised-bed gasification experiments with wood, peat and coal at VTT in 1991–1994. Part 2.* ISBN 951-38-4787-X。Hamel 引用的是 Part 1（1993，VTT Publication 161，ISBN 951-38-4399-8），书目键名 `KURKELA_1991`。Part 2 的前言写明试验装置的详细描述在 Part 1。它自己的数据是 28 个 setpoint，包括硬煤和热气体过滤，不是表 7.2 那一组。对不上的列没有采用。

p. 12 把莱茵褐煤写成筛到 4 mm 以下、泥炭筛到 0–4 mm，与 1992 年论文 p. 3 相同。采用的出处仍是 1992 年论文。

## 与 `validation_cases.json` 的输入差异

旧文件的出口值不在这次对照里。工况 1、6、13、14 的燃料量、一次风、二次风、蒸汽、压力、λ、循环，以及 13、14 的元素分析，与表 7.2 相应列一致。下面是没有采用的差异。

| 项 | 本文件 | 旧文件 |
|---|---|---|
| 13、14 的页码 | 续表 S. 134 | 写成 Table 7.2 (p.133) |
| 工况 13 的水分 | 16,1 Gew.-% | proximate 里写成 16,0。16,0 是工况 12 那一列 |
| 工况 6、13 的水分区间 | 不设区间，各列单独记 | 11,3–18,5，把锯末和其他泥炭列混进一条 |
| 进气温度 | 对上表 6 的五列为 185 °C；其余仍无 | 293 K |
| 床高 | 待读 | 1,2 m，注明来自表 7.5 |
| 二次风高度 | 待读 | 分配器以上 800 mm |
| 粒径 | 原文无 | 1–4 mm，或锯末 0,1–1 mm |
| 旋风效率 | 原文无；只记两个 cell | 工况 13 写约 90–95% |
| 由 λ 回算的空气量 | 不采用。工况 1 的核算见上，表值仍是 81,9 | 工况 6、13、14 另有一组 calculated_flows |
