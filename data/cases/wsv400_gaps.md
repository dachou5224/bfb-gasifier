# WSV400 输入缺口与新旧差异

对照范围：表 2.1、表 2.2（S. 13）的每一项都在 `wsv400.json` 里有字段。下面只列 `missing`，以及和 `data/validation_cases.json` 的输入差异。旧文件没有改。

## 表 7.4 的分组

水分是逐列的：A 为 27,3 Gew.-%，B–E 为 25,82 Gew.-%，F 为 22,8 Gew.-%。

工业分析和元素分析各有两列数字。前一列的横坐标落在 C 列中心，与 A–E 五列的中点一致；后一列落在 F 列。因此 A–E 共用前一组，F 用后一组。F 的 O（wf）按印刷值 4,16 Gew.-% 录入，没有改成别的数。

A–E 的 C、H、O、N、S 与灰分之和为 99,97 Gew.-%。F 的相应之和为 62,54 Gew.-%。F 这一组对不上 100。是否为印刷差错，由总司令裁定；在裁定前不替换。

## 缺口

| 字段 | 模型里用在哪里 | 原文查过的位置 | 可能的处理 |
|---|---|---|---|
| `heating_value` | 能量衡算、进料焓 | 表 2.2；表 7.4；§7.1.5，S. 142–146 | 正文只有装置 Feuerungswärmeleistung 100 kW，不是燃料热值。不能用旧代码里的 22 MJ/kg |
| `primary_air_amount_flow` | 气化剂摩尔流 | 表 7.4 给出 50 m³(V_N)/h | S. 142–146 没有 V_N 的温度和压力，不能换成 mol/s |
| `inlet_gas_mole_fractions` | 进气组成 | 表 2.2；表 7.4 只写 Luft | 原文未给空气的 y_i |
| `inlet_temperature` | 进气焓 | 表 2.2；§7.1.5 | 原文未给。旧文件中的 293 K 没有出处 |
| `bed_inventory_composition` | 床料衡算 | 表 2.2；§7.1.5 | 原文未给。不能用旧代码的床料密度 |
| `bed_inventory_diameter` | 床料流动 | 表 2.2 | 原文未给平均粒径。不能用旧代码的 2.25 mm |
| `bed_height` | 表 2.2 的床高，流动与 cell 划分 | Bild 7.14 有床面记号 | 见读数清单。旧文件的 1.2 m 来自对 Bild 7.16 的 ξ 推断，不是表 7.4 |
| `cell_count`、`cell_heights`、`agent_feed_height` | cell 模型几何 | Bild 7.14 | 见读数清单 |
| `fuel_size_class_mass_fractions` | 燃料粒径分级 | Bild 7.15 | 分级边界已录入；柱高待读。该图只对应 Versuch A。B–F 没有粒径分布 |
| `wall_heat_transfer_coefficient` | 壁面散热 | 表 2.1；§7.1.5 只说矿物保温 | 原文未给系数 |
| `connecting_line_geometry`、`cyclone_geometry`、`return_line_geometry` | 表 2.1 的下游设备 | §7.1.5；Bild 7.14 | 未给尺寸。六个工况 Rezirkulation 均为 Nein |
| `ash_withdrawal_count` | 排灰位置 | 表 2.1；§7.1.5 | 未给出个数和位置 |

## 与 `validation_cases.json` 的输入差异

旧文件的出口值不在这次对照里。

| 项 | 本文件 | 旧文件 |
|---|---|---|
| B–E 水分 | 25,82 Gew.-% | 25.8 |
| B–E 灰分（wf） | 6,34 Gew.-%，与 A 同一组 | 0.23，用了 F 那一列 |
| B–E 的 C | 46,27 Gew.-% | 46.27，但灰分却取了 F 的 0.23，两组被拆开 |
| F 的 O（wf） | 印刷值 4,16 Gew.-% | 41.6，文件内写明是把 4.16 改过 |
| 床高 | 待读 | 1.2 m，由 Bild 7.16 的 ξ 推断 |
| 进气温度 | 原文无 | 293 K |
| 工况 A 的空气量 | 表 7.4 为 50 m³(V_N)/h | 另有一组由 λ 回算的 50.3 m³/h |
| 粒径 | A 的分级边界已录入，柱高待读；B–F 无 | A 写了 2–26 mm 和优势级 14–17 mm |
