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
- 表 4、表 5 另有石灰石进料和氮气进料，表 7.2 没有这两路。氮气已写入各工况的 `nitrogen_purge_volume_flow`（cited）和按 ADR 0004 换算的 `nitrogen_purge_amount_flow`（derived）。泥炭九列（表 4，p. 8，工况 5–13）是 8,8、3,6、3,7、2,3、4,2、3,9、3,2、2,7、2,4 Nm³/h。褐煤和锯末（表 5，p. 9，工况 1–4、14–15）是 6,2、3,2、3,5、2,5、2,0、3,5 Nm³/h。只有工况 1、2 的 note 写按 ADR 0006 暂定计入模型；其余写模型是否计入待裁定（ADR 0006 第 4 项）。吹扫氮温度与燃料加入温度按 ADR 0014 C2 记为 298,15 K，见 T29。石灰石是 1,8、0、0、1,1、0、0、0、0、0,3 kg/h（泥炭）和 0（褐煤、锯末），不在 ADR 0006 范围内，没有写入。
- 床层温度、自由空间温度、床层压降和热损失是试验结果，不是进料。表 6 的热损失是差额，没有换成传热系数。
- 文献写两台旋风。Hamel 在 S. 133 按一个旋风假设两个 cell，`cyclone_cell_count` 仍是 2。旋风尺寸两边都没有。
- 床料成分原文没有给出。p. 16 的石英砂是建议另做的试验，不是这批试验的床料。ADR 0014 A1 的种类裁定没有质量分数，字段仍 missing，见 T29。

`bed_height` 按 ADR 0014 A2 记为 1,2 m（derived，字段在 shared）。表 3 的 1,2 m 仍是床区筒体；操作床高取床面标记与细筒段顶边重合，不由各列压降推算。工况 1 的床料质量另记 `bed_mass`。

## 清单对应

| 表 2.1 / 2.2 的项 | 字段 |
|---|---|
| 流化床直径 | `bed_diameter` = 0,15 m |
| 自由空间直径 | `freeboard_diameter` = 0,25 m |
| 反应器高度 | `reactor_height` = 3,5 m（derived，ADR 0005）；正文的 3,5 m 仍记在 `plant_height` |
| 连接管、返料管尺寸 | `connecting_line_geometry`、`return_line_geometry` |
| 旋风尺寸 | `cyclone_geometry` 为 missing；两格假设记在 `cyclone_cell_count` = 2 |
| 燃料进口个数与高度 | `fuel_feed_count` = 1；`fuel_feed_location` 为紧贴布风板上方；`fuel_feed_height` = 0,12 m（cited） |
| 气化剂进口个数与高度 | `agent_feed_count` = 2；`agent_feed_height` = 1,4875 m（derived，ADR 0014 B） |
| 排灰口 | `ash_withdrawal_count` |
| 返料位置 | `return_location`；各工况开关是 `recirculation` |
| 壁面传热系数 | `wall_heat_transfer_coefficient` |
| 燃料量 | `fuel_mass_flow` |
| 工业分析、元素分析、热值 | `fuels` 与各工况水分；固定碳为 cited。热值只在工况 1、3、7、11、14 |
| 燃料粒径分布 | 工况 1 的 `fuel_size_class_mass_fractions` 为筛分级份额；等表面积换算见 T29 |
| 各进气口体积流量 | `primary_air_volume_flow`、`secondary_air_volume_flow`，摩尔流量按 ADR 0004。吹扫氮另记 `nitrogen_purge_volume_flow` |
| 进入气体组成 | 工况 1 的 `inlet_gas_mole_fractions` 为干空气与纯蒸汽（ADR 0014 C1）；其余工况仍 missing |
| 进入气体温度与压力 | 工况 1、3、7、11、14 的 `inlet_temperature` 为 185 °C；其余仍 missing。`pressure` 取表 7.2 的 Druck |
| 床料组成与粒径 | `bed_inventory_composition` 仍 missing（见 T29）；`bed_inventory_diameter` 仍 missing |
| 床高 | `bed_height` = 1,2 m（derived，ADR 0014 A2） |
| 蒸汽 | `steam_mass_flow`（表 7.2 的 Dampf；与一次风一起从 Düsenboden 进入） |

`secondary_air_amount_flow` 与已有的 `primary_air_amount_flow` 同一换算。WSV400 的二次风是短横，所以那个文件没有这个键。

## 缺口

| 字段 | 模型里用在哪里 | 原文查过的位置 | 可能的处理 |
|---|---|---|---|
| `heating_value` | 能量衡算、进料焓 | 表 6 只覆盖工况 1、3、7、11、14 | 其余工况对不上代表点，不外推。不能用旧文件里的数 |
| `inlet_gas_mole_fractions` | 进气组成 | 表 2.2；表 7.2 只写 Luft/Dampf | 工况 1 已按 ADR 0014 C1 录入。工况 2–15 仍缺。吹扫氮不并进空气 |
| `inlet_temperature` | 进气焓 | 表 6 的 185 °C 只用于对上的五列 | 工况 1、3、7、11、14 的 458,15 K 不改。其余工况不外推。燃料与吹扫氮的温度不在这个字段 |
| `bed_inventory_composition` | 床料种类 | ADR 0014 A1 | 未录入。裁定是焦炭与灰、不含砂，没有质量分数；字段单位是 kg/kg |
| `bed_inventory_diameter` | 床料流动、Geldart | ADR 0014 X3；ADR 0015 | 仍缺。不把燃料筛分当作床料粒径 |
| `bed_particle_density` | ADR 0009 的 ρ_P | ADR 0014 A4 | 不录主值。外层由床区积分等于工况 1 的 `bed_mass` 反求 |
| `cell_count` | 反应器分格 | Bild 7.7 | 已录入 16（床 8、扩径 1、自由空间直筒 7），不含旋风。见下文“cell 划分” |
| `distributor_orifice_count`、`distributor_outlet_height` | 式 (3.37)、(3.39) | ADR 0019 | 已录入 139,5 与 0。第二块板 95,6 只在 note |
| `wall_heat_transfer_coefficient` | 壁面散热 | 表 2.1；Table 6, p. 18 | 工况 1 对上的差额是 21,7 kW，不是传热系数 |
| `connecting_line_geometry`、`cyclone_geometry`、`return_line_geometry` | 表 2.1 的下游设备 | S. 133；KURKELA_1992 p. 6；KURKELA_1991 印刷页 7–10 | 两台旋风，没有尺寸。`cyclone_cell_count` 仍是 2 |
| `ash_withdrawal_count`、`ash_withdrawal_mass_flow` | ADR 0009 第 6 项的排灰位置与流率 | KURKELA_1992 p. 6；KURKELA_1991 印刷页 8 | 仍缺。35 mm 与 38 mm 不一致，取样管原值不改，流率没有 |
| `fuel_size_class_mass_fractions` | 燃料粒径分级 | ADR 0014 D2；Fig. 1 | 工况 1 已录入筛分级份额。`<0,1 mm` 的代表粒径 ADR 未写明，等表面积换算未做。工况 2–15 仍缺 |
| `distributor_type` | ADR 0011 第 8、9 条 | KURKELA_1992 p. 5 | 已录入 multiorifice plate distributor。孔数另见 `distributor_orifice_count` |

一次风、二次风和吹扫氮的摩尔流已按 ADR 0004 换算，V_M = 0,022414 m³/mol。它们不再列入上表。`reactor_height` 已按 ADR 0005 记为 3,5 m。`bed_height`、`cell_heights`、`agent_feed_height`、`solid_feed_temperature`、`nitrogen_purge_temperature` 已按 ADR 0014 录入，不再列入上表。吹扫氮温度是 298,15 K。

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

T15 又看了封面、摘要、前言和目录（PDF 页 1–6，印刷页含 3–6）。前言写明这是 1991–1994 年计划的第二部分，装置描述在 VTT Publications 161。摘要的 28 个稳态点、0,4–1,0 MPa，对不上表 7.2 的 1988–1990 年工况。没有逐字段对上同一次试验的值，本轮仍不录入。

## 与 `validation_cases.json` 的输入差异

旧文件的出口值不在这次对照里。工况 1、6、13、14 的燃料量、一次风、二次风、蒸汽、压力、λ、循环，以及 13、14 的元素分析，与表 7.2 相应列一致。下面是没有采用的差异。

| 项 | 本文件 | 旧文件 |
|---|---|---|
| 13、14 的页码 | 续表 S. 134 | 写成 Table 7.2 (p.133) |
| 工况 13 的水分 | 16,1 Gew.-% | proximate 里写成 16,0。16,0 是工况 12 那一列 |
| 工况 6、13 的水分区间 | 不设区间，各列单独记 | 11,3–18,5，把锯末和其他泥炭列混进一条 |
| 进气温度 | 对上表 6 的五列为 185 °C；其余仍无。燃料与吹扫氮另为 298,15 K（ADR 0014 C2） | 293 K |
| 床高 | 1,2 m（ADR 0014 A2）。不是由表 7.5 或压降推得 | 1,2 m，注明来自表 7.5 |
| 二次风高度 | 1,4875 m（ADR 0014 B）。不是 1991 年报告的 1,5 m、1,9 m | 分配器以上 800 mm |
| 粒径 | 工况 1 为 Fig. 1 的筛分级份额。不是旧文件的 1–4 mm | 1–4 mm，或锯末 0,1–1 mm |
| 旋风效率 | 原文无；只记两个 cell | 工况 13 写约 90–95% |
| 由 λ 回算的空气量 | 不采用。工况 1 的核算见上，表值仍是 81,9 | 工况 6、13、14 另有一组 calculated_flows |

## T15：引用核实、装置对照与 cell 拓扑

本段由新的 cursor 会话写入（2026-10-17）。没有沿用做过 T4、T8、T9、T11–T14 的会话。来源顺序是 Hamel，然后只补查 KURKELA_1992 里尚未写入的字段，然后 KURKELA_1991，最后 VTT 249。

### Hamel 在何处引用 KURKELA_1991

`refs/hamel/index.json` 里 `KURKELA_1991` 的引用位置是作者和年份的自动匹配，和 `KURKELA_1992`、Leppälahti 与 Kurkela（1991）混在一起。对照页面图像之后，正文里指向 VTT Publications 161 的只有一处。

PDF 76 右 / S. 133，§7.1.3。Hamel 先写模拟数据取自 Kurkela und Stahlberg（1992），并写旋风实际尺寸没有给出、模型按 Bild 7.7 用两个 cell。紧接着的句子是：„Weitere Einzelheiten zu der Versuchsanlage des VTT sind zum Beispiel in Kurkela et al. (1993) zu finden.“ 书目里的 Kurkela、Ståhlberg 与 Laatikainen（1993）就是 VTT Publications 161（ISBN 951-38-4399-8，封面为 Espoo 1993）。用途是装置的进一步细节，不是另一套工况，也不是动力学。

其余自动匹配都不是这本书：

| 页面 | 图像上的引文 | 实际文献 |
|---|---|---|
| PDF 13 右 / S. 7 | Kurkela und Stahlberg (1992)，白云石降低焦油 | KURKELA_1992 |
| PDF 56 / S. 92–93 | Kurkela und Stahlberg (1992)；Leppälahti und Kurkela (1991) | 1992 年加压试验；1991 年 *Fuel* 上的常压泥炭，不是 VTT 161 |
| PDF 77–78 / S. 134–137 | Kurkela und Stahlberg (1992) | KURKELA_1992 |
| PDF 79–80 / S. 138–141 | Leppälahti und Kurkela (1991) | 常压装置，§7.1.4 |
| PDF 83–84 / S. 147–148，PDF 87 / S. 154–155 | Kurkela und Stahlberg (1992)，以及常压那篇 | 同上 |

因此 CHARTER §4.2 的使用范围就是 S. 133 这句话：补装置细节。不从这本书引入模型或关联式。

### 装置是否与 1988–1990 年试验时相同

KURKELA_1992 的试验是 1988 年 4 月到 1990 年 3 月（p. 6）。KURKELA_1991 印刷页 7 写 1988–1990 年做的是泥炭和褐煤的杂质，结果已在较早的文献里；这本书自己的试验是 1991–1992 年的木材和硬煤。印刷页 5 的前言写这套 PDU 从 1988 年就在用。直径一致，几处竖向尺寸和接管直径不一致。不一致的数不录入。

| 项目 | KURKELA_1992（1988–1990） | KURKELA_1991（装置说明，印刷页） | 处理 |
|---|---|---|---|
| 床径 | 0,15 m，表 3，p. 5 | 15 cm，印刷页 7 与表 1（印刷页 10） | 已有，不改 |
| 自由空间直径 | 0,25 m，表 3 | 25 cm，同上 | 已有，不改 |
| 床区筒体高度 | 1,2 m，表 3 | 未单列 | 不把 4,2 m 写进筒体高度 |
| 自由空间高度 | 2,3 m，表 3 | 未单列 | 不改 |
| 分配器到气体出口 | 3,7 m，p. 4 | 4,2 m，印刷页 7 与表 1 的 Reactor height | 不一致，4,2 m 不录入 |
| 燃料口高度 | 分配器以上 0,12 m，p. 5–6 | 锁斗系统的进料口约在分配器以上 15 cm，印刷页 8 | 不一致，仍用 0,12 m |
| 布风器 | multiorifice plate distributor，p. 5 | 同一短语；另有两块孔板（1 mm 孔、开孔率 0,62 %；2 mm 孔、开孔率 1,7 %），印刷页 8 | 短语已录入。孔径和开孔率 1992 年论文没有，不录入 |
| 二次风 | above the fluidized bed，没有高度，p. 5 | 两根管，各八个喷嘴，高度 1,5 m 和 1,9 m，印刷页 8 | 高度不录入。竖向总尺寸已经对不上，不能确认这两处高度属于 1988–1990 年的装置 |
| 底部管子 | 35 mm 取样管，p. 6 | 分配器中心 38 mm 内径管排底灰，印刷页 8 | 不一致。已录入的 35 mm 不改，38 mm 不录入，个数和流率仍缺 |
| 旋风 | 两台，没有尺寸，p. 6 | 印刷页 10：1988 年只用旋风净化；陶瓷过滤器在 1989 年试验前建成；1991 年春加了气体冷却器 | 没有尺寸。过滤器和冷却器的尺寸不录入 |
| 床料 | p. 16 的石英砂是建议另做的试验 | 锯末试验提到惰性床料，不是表 7.2 的褐煤工况 | 组成、粒径、ρ_P 仍缺 |
| 操作床高 | 没有 | 没有膨胀床高 | 仍缺 |
| 给料温度、吹扫氮温度 | 没有 | 锁斗用氮气加压，没有温度 | 仍缺 |
| 进气组成 | air and steam | Air, steam | 仍缺，不取 0,21 |
| 燃料分级质量分数 | Fig. 1 累积曲线，p. 3 只有筛分上限 | 不是本卡要补的 1988–1990 年褐煤筛分表 | 仍缺 |

KURKELA_1992 本轮补查的页：p. 3（筛分）、p. 5（布风器、二次风、锁斗氮）、p. 6（床料取样、两台旋风）、p. 8–9（表 4、表 5 的床层压降）、p. 16（石英砂是建议）、p. 18（表 6 的热损失差额）。p. 4、表 1、表 3 的几何在 T10 已经录入，本轮只用来和 1991 年报告对照。

### 本轮录入

`distributor_type` = multiorifice plate distributor，`cited`，KURKELA_1992 p. 5。式 (3.11)、(3.14) 要的技术布风器、多孔板或喷嘴板，这篇英文没有对应词，留给总司令。

`bed_pressure_drop` 按表 4、表 5 写入各工况，单位 Pa。工况 1 为 2,2 kPa。锯末最后一列印的是 0,0 kPa，照录。测点位置两篇都没有写明。这些数不用来推床高。

`cell_count` = 16，`source`，只计 Bild 7.7 的反应器。`cell_heights` 仍是 missing。

### cell 划分

S. 16（PDF 18 左）：„die beliebig erfolgen kann“。划分可以任意选取。S. 79（PDF 49 右）的 „eine realistische Gesamtzellenanzahl von 20“ 是在算方程个数时举的规模，不是 VTT 这张图的格子数。

Bild 7.7（PDF 76 左 / S. 132）没有刻度。下面的个数是柱内连续横线的拓扑，没有把像素换成米。床区 8 格的横线曾与放大图上的编号对齐。

自下而上：

| 段 | 格数 | 图上的标记 |
|---|---|---|
| 床区（标 Wirbelschicht 的细柱） | 8 | 最下一格有一次风箭头和燃料箭头。返料管接在从下数第二格的侧面。床面的倒三角在最上一格的顶部分界上 |
| 锥形扩径 | 1 | 在床区与直筒之间，没有再分 |
| 自由空间直筒（标 Freiraum） | 7 | 二次风箭头在最下一格的下沿，紧挨着锥段。最上一格接出水平管 |
| 去旋风的水平管、返料管 | 没有再分格 | 各画成一根管道 |
| 旋风 | 2 | 与已录入的 `cyclone_cell_count` 一致，不计入 16 |

床区一格的轴向间距大约是自由空间直筒一格的一半。这只说明疏密方向，不给出米制高度。

其它图只作对照，同样不量像素：

| 图 | 床区与自由空间 |
|---|---|
| Bild 7.2，HTW Wesseling，PDF 73 左 / S. 126 | 上部自由空间约 6 格，其下还有更密的侧面进气段（约 3 格加约 5 格），细床区 2 格 |
| Bild 7.3，Berrenrath，PDF 73 右 / S. 127 | 直段约 10 格，下部锥段约 6 格，床的标注在锥段 |
| Bild 7.11，VTT 常压，PDF 79 左 / S. 138 | 上部较疏，底部约 8 格较密；旋风另画 |
| Bild 7.14，WSV400，PDF 81 左 / S. 142 | 自由空间直筒 7 格，下部锥段和床更密 |
| Bild 7.8、7.9 左侧小图，PDF 77 右、PDF 78 左 | 自由空间约 8 格，床约 8 格，床的间距约为自由空间的一半；顶部旋风 2 格 |

Bild 7.7 与 7.8、7.9、7.11、7.14 都是床比自由空间密。Bild 7.2 的侧面进气段很密，床本身只有 2 格，和加压 VTT 不是同一套习惯。这些个数给总司令写 cell 高度的 ADR 用。`cell_heights` 不在本卡填。

Bild 7.8、7.9 剖面转折点的 ξ 由目标所有者用工具读，清单在 `vtt_pressurized_reading_list.md`。ξ 按全部部件的总高度归一化，连接管和旋风高度未知，不能换成米。

## T29：已裁定输入的录入（2026-10-18）

只把 ADR 0014、0019 已经写明的数写入 `vtt_pressurized.json`。下表是这次写入的字段。T15 各段里的“仍缺”是 2026-10-17 当时的检索结果。

| 字段 | 值 | status | 位置 |
|---|---|---|---|
| `bed_height` | 1,2 m | derived，ADR 0014 A2；Bild 7.7, S. 132；KURKELA_1992 表 3, p. 5 | shared。15 个工况共用 |
| `cell_heights` | 8 × 0,15 m + 8 × 0,2875 m | derived，ADR 0014 B；S. 132；表 3, p. 5 | shared。自下而上，不含旋风 |
| `agent_feed_height` | 1,4875 m | derived，ADR 0014 B。1,2 m + 2,3 m / 8 | shared。自由空间第 2 格下沿 |
| `solid_feed_temperature` | 298,15 K | derived，ADR 0014 C2；表 2.2, S. 13 | shared。只覆盖燃料 |
| `nitrogen_purge_temperature` | 298,15 K | derived，ADR 0014 C2；表 2.2, S. 13 | shared。只覆盖料斗吹扫氮 |
| `cases.1.inlet_gas_mole_fractions` | 干空气 y_O₂ = 0,21、y_N₂ = 0,79；蒸汽 y_H₂O = 1 | derived，ADR 0014 C1；表 7.2, S. 133 | 只工况 1 |
| `cases.1.bed_mass` | 2200 × 0,017671 / 9,81 = 3,9629153924566762 kg | derived，ADR 0014 A3；表 5, p. 9；表 2.2, S. 13 | 只工况 1 |
| `cases.1.fuel_size_class_mass_fractions` | 0,116、0,086、0,223、0,117、0,143、0,249、0,066 | derived，ADR 0014 D2；Fig. 1, p. 4；S. 133 | 筛分级，只工况 1 |
| `distributor_orifice_count` | 139,5 | derived，ADR 0019；KURKELA_1991 印刷页 8；式 (3.37), S. 29 | shared。n₀ = φ·(D/d_h)² |
| `distributor_outlet_height` | 0 m | derived，ADR 0019 第 1 项；式 (3.39), S. 29 | shared |

工况 1、3、7、11、14 的 `inlet_temperature` 仍是 458,15 K（KURKELA_1992 表 6, p. 18，cited），没有改。其余工况的 `inlet_temperature` 仍是 missing。

`bed_particle_density` 不录主值。note 写明由 ADR 0014 A4 在外层计算。ε_mf、φ_S 不录。

孔数 note 写明：年代一致性未确认，目标所有者 2026-10-18 批准作跨时期外推假设。第二块板 n₀ = 95,6 只写在 note 中。

停下、没有自定的两项：

1. `<0,1 mm` 一级的代表粒径。ADR 0014 D2 只写它是解释性选择、在数据 PR 登记；`kurkela1992_fig1.json` 也写明未给出。Hamel S. 144 是把木片边长矩阵换成等表面积球的直径分级，没有从这条筛分曲线算出该档代表粒径的式子。筛分级质量分数已录入，等表面积换算没有做。
2. `bed_inventory_composition`。ADR 0014 A1 裁定工况 1 的床料是燃料生成的焦炭与灰、不含惰性砂，没有质量分数。字段单位是 kg/kg。两者不是同一个量，值保持 missing。

