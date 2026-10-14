# VTT 加压输入缺口与编号核对

对照范围：表 2.1、表 2.2（S. 13）的每一项都在 `vtt_pressurized.json` 里有字段。数值取自表 7.2（PDF 76 右 / S. 133，以及 PDF 77 左 / S. 134 的续表）和 §7.1.3（S. 132–137）。Bild 7.7（PDF 76 左 / S. 132）上没有尺寸线的量见读数清单。旧文件没有改。本阶段没有打开 Kurkela & Ståhlberg（1992）和 Kurkela et al.（1993）。

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

## 清单对应

| 表 2.1 / 2.2 的项 | 字段 |
|---|---|
| 流化床直径 | `bed_diameter` = 0,15 m |
| 自由空间直径 | `freeboard_diameter` = 0,25 m |
| 反应器高度 | `reactor_height` 为 missing；正文的 3,5 m 是 Gesamtanlagenhöhe，记在 `plant_height` |
| 连接管、返料管尺寸 | `connecting_line_geometry`、`return_line_geometry` |
| 旋风尺寸 | `cyclone_geometry` 为 missing；两格假设记在 `cyclone_cell_count` = 2 |
| 燃料进口个数与高度 | `fuel_feed_count` = 1；`fuel_feed_height` 为 missing |
| 气化剂进口个数与高度 | `agent_feed_count` = 2；`agent_feed_height` 为 missing |
| 排灰口 | `ash_withdrawal_count` |
| 返料位置 | `return_location`；各工况开关是 `recirculation` |
| 壁面传热系数 | `wall_heat_transfer_coefficient` |
| 燃料量 | `fuel_mass_flow` |
| 工业分析、元素分析、热值 | `fuels` 与各工况 `moisture_mass_fraction`；`heating_value` 为 missing |
| 燃料粒径分布 | `fuel_size_class_mass_fractions` |
| 各进气口体积流量 | `primary_air_volume_flow`、`secondary_air_volume_flow`，摩尔流量按 ADR 0004 |
| 进入气体组成 | `inlet_gas_mole_fractions` |
| 进入气体温度与压力 | `inlet_temperature` 为 missing；`pressure` 取表 7.2 的 Druck |
| 床料组成与粒径 | `bed_inventory_composition`、`bed_inventory_diameter` |
| 床高 | `bed_height` |
| 蒸汽 | `steam_mass_flow`（表 7.2 的 Dampf；与一次风一起从 Düsenboden 进入） |

`secondary_air_amount_flow` 与已有的 `primary_air_amount_flow` 同一换算。WSV400 的二次风是短横，所以那个文件没有这个键。

## 缺口

| 字段 | 模型里用在哪里 | 原文查过的位置 | 可能的处理 |
|---|---|---|---|
| `heating_value` | 能量衡算、进料焓 | 表 2.2；表 7.2；§7.1.3，S. 132–137 | 正文的 20–150 MW_el 是电站规模，不是本装置燃料热值。不能用旧文件里的数。第二阶段再查 Hamel 引用的 Kurkela 文献 |
| `inlet_gas_mole_fractions` | 进气组成 | 表 2.2；表 7.2 只写 Luft/Dampf | 空气的 y_i 原文未给。下面相容性检查里的 0,21 只用于核算，不是进气组成 |
| `inlet_temperature` | 进气焓 | 表 2.2；§7.1.3，S. 132–137 | 原文未给温度。旧文件的 293 K 没有出处，没有采用 |
| `bed_inventory_composition`、`bed_inventory_diameter` | 床料衡算与流动 | 表 2.2；§7.1.3 | 原文未给。不能用旧文件的 1–4 mm |
| `bed_height` | 表 2.2 的床高 | §7.1.3 无数值；Bild 7.7 有床层线但无尺寸 | 见读数清单。旧文件的 1,2 m 写明来自表 7.5，表 7.5 不在本阶段来源里，没有采用 |
| `cell_count`、`cell_heights`、`fuel_feed_height`、`agent_feed_height` | cell 划分与进料位置 | Bild 7.7；S. 134 只有区域描述 | 见读数清单。旧文件二次风口「分配器以上 800 mm」没有写在 §7.1.3 里，没有采用 |
| `reactor_height` | 反应器本体高度 | S. 133 的词是 Gesamtanlagenhöhe 3,5 m | 3,5 m 已按原词记入 `plant_height`。是否等于反应器高度，待总司令裁定，没有直接并入 `reactor_height` |
| `wall_heat_transfer_coefficient` | 壁面散热 | 表 2.1；§7.1.3 | 原文未给系数 |
| `connecting_line_geometry`、`cyclone_geometry`、`return_line_geometry` | 表 2.1 的下游设备 | S. 133 写明旋风实际尺寸没有给出；图上无尺寸线 | 旋风按原文假设为两个 cell，已记入 `cyclone_cell_count`。尺寸仍缺。第二阶段若文献里其实有尺寸，与 S. 133 这句话不一致，把差异写在这里，由总司令裁定 |
| `ash_withdrawal_count` | 排灰位置 | §7.1.3 未写；Bild 7.7 待看 | 见读数清单 |
| `fuel_size_class_mass_fractions` | 燃料粒径分级 | 表 7.2 与 §7.1.3 都没有分布图 | 原文未给。旧文件的粒径范围没有采用 |

一次风和二次风的摩尔流已按 ADR 0004 换算，V_M = 0,022414 m³/mol。它们不再列入上表。

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
| S. 133 | Gesamtanlagenhöhe 3,5 m | 没有另一条叫 Reaktorhöhe 的数 |

图框里的直径 0,15 m / 0,25 m、气化剂 Luft/Dampf，与 S. 132 正文一致。

## 与 `validation_cases.json` 的输入差异

旧文件的出口值不在这次对照里。工况 1、6、13、14 的燃料量、一次风、二次风、蒸汽、压力、λ、循环，以及 13、14 的元素分析，与表 7.2 相应列一致。下面是没有采用的差异。

| 项 | 本文件 | 旧文件 |
|---|---|---|
| 13、14 的页码 | 续表 S. 134 | 写成 Table 7.2 (p.133) |
| 工况 13 的水分 | 16,1 Gew.-% | proximate 里写成 16,0。16,0 是工况 12 那一列 |
| 工况 6、13 的水分区间 | 不设区间，各列单独记 | 11,3–18,5，把锯末和其他泥炭列混进一条 |
| 进气温度 | 原文无 | 293 K |
| 床高 | 待读 | 1,2 m，注明来自表 7.5 |
| 二次风高度 | 待读 | 分配器以上 800 mm |
| 粒径 | 原文无 | 1–4 mm，或锯末 0,1–1 mm |
| 旋风效率 | 原文无；只记两个 cell | 工况 13 写约 90–95% |
| 由 λ 回算的空气量 | 不采用。工况 1 的核算见上，表值仍是 81,9 | 工况 6、13、14 另有一组 calculated_flows |
