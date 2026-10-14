---
status: proposed
---

# 成组批准 VTT 加压装置的文献补充字段（KURKELA_1992）

`data/cases/vtt_pressurized.json` 中有 23 个 `status: cited` 字段，取自 Kurkela & Ståhlberg (1992)，*Fuel Processing Technology* 31, 1–21（索引键 `KURKELA_1992`，本地 `docs/kurkela1992.pdf`）。按 T10 任务卡，`cited` 字段须经成组 ADR 批准后才能被模型装配使用。本 ADR 逐项裁定。

## 来源资格

Hamel 在 §7.1.3 明确写明表 7.2 的运行数据取自 Kurkela und Stahlberg (1992)（PDF 76 右 / S. 133）。按 CHARTER §4.2，该文献在“VTT 加压装置的装置与运行数据”范围内可以使用，不涉及 §4.3 的外部来源批准。不得从中引入 Hamel 没有采用的模型或关联式。

## 同一次试验的依据

- Hamel 表 7.2 是 Kurkela Table 4、5 的转录：泥炭九列的燃料量 68,8 / 65,5 / 61,6 / 63,7 / 76,3 / 62,5 / 45,0 / 26,5 / 50,2 kg/h（Kurkela p. 8，Table 4）与 Hamel 工况 5–13 一一对应。
- Kurkela Table 6（p. 18）的六个代表点中，五个可唯一对应到 Hamel 工况。依据是用 Hamel 表 7.2 的数值独立计算的汽燃比（kg 蒸汽 / kg 燃料 daf）：

| Hamel 工况 | 燃料、是否循环 | 由表 7.2 算得的汽燃比 | Table 6 | Table 6 的 ER（表 7.2 的 λ） |
|---|---|---|---|---|
| 1 | 褐煤，不循环（四个褐煤工况中唯一不循环者） | 0,172 | 0,17 | 0,32（0,33） |
| 3 | 褐煤，循环（循环工况中唯一 λ = 0,44 者） | 0,195 | 0,20 | 0,44（0,44） |
| 7 | 泥炭，不循环 | 0,109 | 0,11 | 0,32（0,32） |
| 11 | 泥炭，循环 | 0,179 | 0,18 | 0,40（0,40） |
| 14 | 锯末，不循环 | 0,077 | 0,08 | 0,27（0,28） |

  Table 6 的第六个点（锯末，ER 0,38，汽燃比 0,38）与 Hamel 工况 15（λ 0,39，汽燃比 0,34）对应不上，**不采用**。工况 1、14 的 ER 与表 7.2 相差 0,01，按表 7.2 取 λ（以 Hamel 为准）。

## 批准为模型输入（T9 可装配）

| 字段 | 值 | 文献出处 | 说明 |
|---|---|---|---|
| `shared.fuel_feed_height` | 0,12 m | p. 5–6 | 燃料口在布风板以上 0,12 m；与 Hamel S. 134“unmittelbar über dem Düsenboden”相符 |
| `shared.bed_zone_height` | 1,2 m | Table 3, p. 5 | 床区筒体高度（直径 0,15 m）；**不是**操作床高，`bed_height` 仍为 missing |
| `shared.freeboard_height` | 2,3 m | Table 3, p. 5 | 自由空间段高度（直径 0,25 m） |
| `fuels.*.fixed_carbon_mass_fraction_wf` | 褐煤 0,427，泥炭 0,274，锯末 0,168 | Table 1, p. 3 | 干基平均值（表头 “wt%, dry basis”） |
| `fuels.rhein_braunkohle.particle_size_max` | 0,004 m | p. 3 | “sieved below 4 mm” |
| `fuels.finn_torf.particle_size_range` | [0, 0,004] m | p. 3 | “sieved to 0–4 mm” |
| `cases.{1,3,7,11,14}.inlet_temperature` | 458,15 K | Table 6, p. 18 | 空气与蒸汽温度同为 185 °C |
| `cases.{1,3,7,11,14}.heating_value` | 1：21,0；3：20,9；7、11：17,3；14：16,8 MJ/kg | Table 6, p. 18 | **收到基低位热值**（Table 6 “Fuel LHV (wet)”），对应该试验点的燃料含水量；模型使用时须按此基准换算 |

## 只作记录，不作模型输入

以下字段属实，但不是 Hamel 模型所需的输入；T9 不得把它们装配为几何或运行参数：

| 字段 | 值 | 说明 |
|---|---|---|
| `shared.distributor_to_outlet_height` | 3,7 m | 布风板到出气管的距离，含出口连接段 |
| `shared.pressure_vessel_height` | 5,6 m | 压力容器外壳 |
| `shared.pressure_vessel_inside_diameter` | 0,6 m | 压力容器内径（耐火衬里以外） |
| `shared.ash_sample_pipe_diameter` | 0,035 m | 床料取样管，不是排灰口的计数 |
| `shared.thermal_design_power` | 400 kW | 装置最大热输入（泥炭），不是各工况的热功率 |

## 解释性选择：反应器高度取 3,5 m

Hamel 写 VTT 装置的“Gesamtanlagenhöhe beträgt 3,5 m”（S. 133）；Kurkela Table 3 的床区 1,2 m 与自由空间 2,3 m 之和也是 3,5 m；布风板到出气管的 3,7 m 多出的 0,2 m 未说明归属。两处相互印证，因此模型的反应器高度（布风板至自由空间顶部）取 **3,5 m = 1,2 m + 2,3 m**。`shared.reactor_height` 据此可记为 `derived`，`formula` 写明“bed_zone_height + freeboard_height，与 Hamel S. 133 的 Gesamtanlagenhöhe 3,5 m 一致（ADR 0005）”。代码中注明 `解释性选择：ADR 0005`。

## 尚未裁定（第 1 关运行前须定）

1. **料斗吹扫氮气与石灰石**：Kurkela Table 4、5 有、Hamel 表 7.2 没有（褐煤工况 1 的氮气为 6,2 Nm³/h，约占总风量的 7,6%）。北极星是复现 Hamel 的模拟结果，须先判断 Hamel 是否计入；拟在 T3 读出 Bild 7.10（必要时 Bild 7.8）后，用 Hamel 自己的模拟结果判断。
2. **挥发分的基准**：Hamel 表 7.2 标为“Flüchtige (waf)”，数值（褐煤 53,0、泥炭 68,4、锯末 83）与 Kurkela Table 1 的**干基**挥发分相同；按干基，挥发分 + 固定碳 + 灰分 = 100。若按无水无灰基理解，褐煤应为 55,4。这影响热解中挥发分与焦炭的分配，须在热解模块（第 4 章）的原文答卷中结合 Hamel 的公式裁定。
3. **热损失**：Table 6 给出工况 1、3 的热损失（按能量平衡差额）为 21,7、25,0 kW；Hamel 表 2.1 要求的是壁面传热系数。
4. **旋风返料**：Kurkela p. 6 写明只有第一级旋风的细粉返回床层，第二级旋风和过滤器的细粉不返回；Hamel 按两个 cell 假设旋风（S. 133）。返料取自哪一个 cell，须在建立旋风模型时裁定。
5. **工况 2、4 的热值与进气温度**：Table 6 未覆盖，属于第 2 关。

## 证据状态

总司令 2026-10-14 在 `docs/kurkela1992.pdf` 的页面图像上核对了 p. 3（Table 1）、p. 4–6、p. 8（Table 4）、p. 18（Table 6），以及 Hamel PDF 76、77 的表 7.2 与续表。

## 批准记录

待 codex 审查后由总司令批准（CHARTER v0.7 §8：ADR 由总司令批准，目标所有者可否决；本 ADR 不引入 Hamel 未引用的来源）。
