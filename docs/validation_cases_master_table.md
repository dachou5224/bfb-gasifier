# Validation Cases 主表（15 Case / 4 套装置）

> 自动生成自 [`data/validation_cases.json`](../data/validation_cases.json)，脚本 [`scripts/build_validation_cases_master_table.py`](../scripts/build_validation_cases_master_table.py)。 更新日期：2026-09-07

15 个 Hamel 验证工况对应 **4 套物理装置**：Wesseling 中试 2 组、Berrenrath 示范 1 组、VTT 实验台 6 组（常压/加压为同一装置的不同压力工况）、WSV 400 研究炉 6 组。

单位：压力 MPa；流量 Nm³/h（STP）；燃料分析 wt%；出口摩尔分数为 **湿基** 或 **干基**（见「气基」列）。

## 0. 装置分组（Validation 用）

| 装置组 ID | 名称 | 规模 | Case 数 | 运行期 | Hamel | 数据出处 | 分组说明 |
|-----------|------|------|---------|--------|-------|----------|----------|
| `PLANT_HTW_WESSELING` | HTW Wesseling 加压中试炉 | 中试 | 2 | 1989–1992 | §7.1.2, Bild 7.2 | Rheinbraun AG & Krupp-Uhde GmbH (1997); Hamel Ch.7 Table 7.1 / Fig 7.4–7.5 | 同一中试炉两组对比工况（Air/Steam vs O₂/Steam）；可联合验收氧化剂分支与二级风 |
| `PLANT_HTW_BERRENRATH` | HTW Berrenrath 示范工厂 | 工业示范 | 1 | 1986–1997 | §7.1.2, Bild 7.3 | Rheinbraun AG & Krupp-Uhde GmbH (1997); Hamel Ch.7 Table 7.1 / Fig 7.6 | 工业尺度 O₂/Steam 单点；与 Wesseling 同分支但更大径、更少分级 |
| `PLANT_VTT` | VTT Espoo 实验台（常压+加压） | 实验台 | 6 | — | §7.1.3 Bild 7.7 (pressurized); §7.1.4 Bild 7.11 (atmospheric) | Pressurized: Kurkela & Ståhlberg (1992), Kurkela et al. (1993); Atmospheric: Leppälahti & Kurkela (1991) | 按燃料(褐煤/泥炭/锯末)与压力(0.1 vs 0.4–1.0 MPa)分子组；Case D/E 专用于 fuel-N |
| `PLANT_WSV400` | 马格德堡 WSV 400 研究炉 | 100 kW 研究中试 | 6 | — | §7.1.5, Bild 7.14 | Prof. Dr.-Ing. habil. P. Käferstein; Hoyer et al. (1999); Hamel Ch.7 Table 7.4 / Fig 7.16–7.17 | 同一炉 ER 扫描(A–F)；Case A 有轴向剖面，B–F 以 parity 出口为主；F 为不同燃料批次 |

### 0.1 Case → 装置组映射

| Case ID | 装置组 | 规模 | 装置实例名 | Hamel Sim |
|---------|--------|------|------------|-----------|
| CASE_HTW_WESSELING_1 | `PLANT_HTW_WESSELING` | 中试 | HTW Wesseling | 1 |
| CASE_HTW_WESSELING_2 | `PLANT_HTW_WESSELING` | 中试 | HTW Wesseling | 2 |
| CASE_HTW_BERRENRATH | `PLANT_HTW_BERRENRATH` | 工业示范 | HTW Berrenrath | 3 |
| CASE_VTT_PRESSURIZED_COAL_1 | `PLANT_VTT` | 实验台 | VTT Pressurized BFB | 1 |
| CASE_VTT_PRESSURIZED_PEAT_6 | `PLANT_VTT` | 实验台 | VTT Pressurized BFB | 6 |
| CASE_VTT_PRESSURIZED_PEAT_13 | `PLANT_VTT` | 实验台 | VTT Pressurized BFB | 13 |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | `PLANT_VTT` | 实验台 | VTT Pressurized BFB | 14 |
| CASE_VTT_ATMOSPHERIC_PEAT_D | `PLANT_VTT` | 实验台 | VTT Atmospheric BFB | D |
| CASE_VTT_ATMOSPHERIC_PEAT_E | `PLANT_VTT` | 实验台 | VTT Atmospheric BFB | E |
| CASE_WSV400_A | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | A |
| CASE_WSV400_B | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | B |
| CASE_WSV400_C | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | C |
| CASE_WSV400_D | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | D |
| CASE_WSV400_E | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | E |
| CASE_WSV400_F | `PLANT_WSV400` | 100 kW 研究中试 | WSV 400 Magdeburg | F |

### 0.2 跨工况实验覆盖范围

以下范围由主表收录的工况汇总而成，不对应某一个单独 Case。出口组成来自论文图的数字化读数，主要用于说明验证窗口，不能替代原始测量表。

| 指标 | 汇总范围 | 说明 |
|------|----------|------|
| 反应器压力 | 0.1–2.5 MPa | 从常压实验台到加压 HTW 装置 |
| 燃料给料量 | 9.8–25,769 kg/h | 按各试验表给出的入炉燃料流量 |
| 化学计量氧比 | 0.238–0.5 | 空气、空气／水蒸气和氧气／水蒸气工况 |
| 燃料含水量 | 8.4–27.3 wt% | 褐煤、泥炭、锯末与木片 |
| 代表性出口／平均自由空间温度 | 920–1220 K（647–947 ℃） | 共 15 组温度数据 |
| 总碳转化率 | 78–97% | 仅统计有实验值的 8 组工况 |

产品气组成必须按湿基和干基分开汇总：

| 气基 | CO | CO₂ | H₂ | CH₄ | H₂O |
|------|----|-----|----|-----|-----|
| 湿基（9 组） | 9.5–33.5 mol% | 10.5–23.5 mol% | 8–25.5 mol% | 2.8–5.2 mol% | 11–19.5 mol% |
| 干基（6 组） | 14.8–18.8 mol% | 12.8–14.2 mol% | 11–19 mol% | 1.8–2.8 mol% | — |

湿基与干基浓度不能直接合并比较；干基组没有水蒸气浓度。各图读数的数字化误差约为 ±5%，范围端点不宜作高精度设计值使用。

## 1. 总览表

| Case ID | 装置组 | Hamel | 装置 | 燃料 | P [MPa] | ER | 循环 | T_out [K] | Xc [%] | CO | CO₂ | H₂ | CH₄ | H₂O | 气基 |
|---------|--------|-------|------|------|---------|-----|------|-----------|--------|-----|-----|-----|-----|-----|------|
| CASE_HTW_WESSELING_1 | PLANT_HTW_WESSELING | 1 | HTW Wesseling | Rhenish Brown Coal (Rhein. Braunkohle) | 2.5 | 0.337 | 是 | 1100 | 95 | 0.13 | 0.11 | 0.12 | 0.028 | 0.175 | wet |
| CASE_HTW_WESSELING_2 | PLANT_HTW_WESSELING | 2 | HTW Wesseling | Rhenish Brown Coal (Rhein. Braunkohle) | 2.5 | 0.268 | 是 | 1220 | 97 | 0.335 | 0.235 | 0.255 | 0.048 | 0.17 | wet |
| CASE_HTW_BERRENRATH | PLANT_HTW_BERRENRATH | 3 | HTW Berrenrath | Rhenish Brown Coal (Rhein. Braunkohle) | 1 | 0.238 | 是 | 1050 | — | 0.285 | 0.205 | 0.255 | 0.048 | 0.195 | wet |
| CASE_VTT_PRESSURIZED_COAL_1 | PLANT_VTT | 1 | VTT Pressurized BFB | Brown Coal (Braunkohle) — same reactor as peat cases | 0.5 | 0.33 | 否 | 980 | 78 | 0.16 | 0.105 | 0.14 | 0.045 | 0.125 | wet |
| CASE_VTT_PRESSURIZED_PEAT_6 | PLANT_VTT | 6 | VTT Pressurized BFB | Finnish Peat (Finnischer Torf) | 0.5 | 0.29 | 否 | 1050 | 86.5 | 0.11 | 0.14 | 0.09 | 0.035 | 0.185 | wet |
| CASE_VTT_PRESSURIZED_PEAT_13 | PLANT_VTT | 13 | VTT Pressurized BFB | Finnish Peat (Finnischer Torf) | 0.5 | 0.41 | 是 | 1175 | 97 | 0.095 | 0.15 | 0.09 | 0.032 | 0.16 | wet |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | PLANT_VTT | 14 | VTT Pressurized BFB | Sawdust (Sägemehl) | 0.4 | 0.28 | 否 | 1045 | 93 | 0.155 | 0.11 | 0.08 | 0.052 | 0.11 | wet |
| CASE_VTT_ATMOSPHERIC_PEAT_D | PLANT_VTT | D | VTT Atmospheric BFB | Finnish Peat (Finnischer Torf) | 0.1 | 0.26 | 是 | 920 | 82 | 0.118 | 0.148 | 0.12 | 0.042 | 0.13 | wet |
| CASE_VTT_ATMOSPHERIC_PEAT_E | PLANT_VTT | E | VTT Atmospheric BFB | Finnish Peat (Finnischer Torf) | 0.1 | 0.29 | 是 | 1100 | 91 | 0.112 | 0.152 | 0.118 | 0.045 | 0.142 | wet |
| CASE_WSV400_A | PLANT_WSV400 | A | WSV 400 Magdeburg | Natural Wood (Holz) | 0.1 | 0.453 | 否 | 980 | — | 0.152 | 0.132 | 0.11 | 0.022 | — | dry |
| CASE_WSV400_B | PLANT_WSV400 | B | WSV 400 Magdeburg | Natural Wood (Holz) | 0.1 | 0.35 | 否 | 1010 | — | 0.155 | 0.128 | 0.16 | 0.028 | — | dry |
| CASE_WSV400_C | PLANT_WSV400 | C | WSV 400 Magdeburg | Natural Wood (Holz) | 0.1 | 0.4 | 否 | 1115 | — | 0.168 | 0.132 | 0.17 | 0.024 | — | dry |
| CASE_WSV400_D | PLANT_WSV400 | D | WSV 400 Magdeburg | Natural Wood (Holz) | 0.1 | 0.45 | 否 | 1120 | — | 0.178 | 0.135 | 0.18 | 0.02 | — | dry |
| CASE_WSV400_E | PLANT_WSV400 | E | WSV 400 Magdeburg | Natural Wood (Holz) | 0.1 | 0.5 | 否 | 1135 | — | 0.188 | 0.142 | 0.19 | 0.018 | — | dry |
| CASE_WSV400_F | PLANT_WSV400 | F | WSV 400 Magdeburg | Natural Wood (Holz) — different batch | 0.1 | 0.45 | 否 | 1105 | — | 0.148 | 0.13 | 0.15 | 0.025 | — | dry |

## 2. 反应器尺寸

| Case ID | D_bed [m] | D_fb [m] | H_total [m] | H_bed [m] | H_exp [m] | 二级风高度 [mm] | 进料高度 [mm] |
|---------|-----------|----------|-------------|-----------|-----------|-----------------|---------------|
| CASE_HTW_WESSELING_1 | 0.6 | 0.6 | 14.5 | 6 | — | — | — |
| CASE_HTW_WESSELING_2 | 0.6 | 0.6 | 14.5 | 6 | — | — | — |
| CASE_HTW_BERRENRATH | 2.75 | 2.75 | 23 | — | — | — | — |
| CASE_VTT_PRESSURIZED_COAL_1 | 0.15 | 0.25 | 3.5 | 1.2 | — | 800 | — |
| CASE_VTT_PRESSURIZED_PEAT_6 | 0.15 | 0.25 | 3.5 | 1.2 | — | 800 | — |
| CASE_VTT_PRESSURIZED_PEAT_13 | 0.15 | 0.25 | 3.5 | 1.2 | — | 800 | — |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | 0.15 | 0.25 | 3.5 | 1.2 | — | 800 | — |
| CASE_VTT_ATMOSPHERIC_PEAT_D | 0.15 | 0.25 | 3.65 | 1.2 | — | 800 | — |
| CASE_VTT_ATMOSPHERIC_PEAT_E | 0.15 | 0.25 | 3.65 | 1.2 | — | 800 | — |
| CASE_WSV400_A | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |
| CASE_WSV400_B | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |
| CASE_WSV400_C | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |
| CASE_WSV400_D | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |
| CASE_WSV400_E | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |
| CASE_WSV400_F | 0.4 | 0.6 | 3 | 1.2 | — | — | 200 |

## 3. 进料与边界

| Case ID | 燃料 [kg/h] | 氧化剂 | 主空气 | 次空气 | 总空气 | 主O₂ | 次O₂ | 总O₂ | 蒸汽 [kg/h] | ξ_sec | T_in [K] |
|---------|-------------|--------|--------|--------|--------|------|------|------|-------------|-------|----------|
| CASE_HTW_WESSELING_1 | 3377.4 | Air/Steam (Luft/Dampf) | — | 30.1 | 5536.9 | — | — | — | — | 0.6 | 293 |
| CASE_HTW_WESSELING_2 | 6858.9 | O2/Steam (O2/Dampf) | — | — | — | — | 37.5 | 1921.6 | — | 0.6 | 293 |
| CASE_HTW_BERRENRATH | 25769 | O2/Steam (O2/Dampf) | — | — | — | — | — | — | — | — | 293 |
| CASE_VTT_PRESSURIZED_COAL_1 | 46.9 | Air/Steam | 51.8 | 30.1 | 81.9 | — | — | — | 6.8 | 0.6 | 293 |
| CASE_VTT_PRESSURIZED_PEAT_6 | 65.5 | Air/Steam | 60.7 | 21.9 | 82.6 | — | — | — | 15.5 | 0.6 | 293 |
| CASE_VTT_PRESSURIZED_PEAT_13 | 50.2 | Air/Steam | 47 | 41.6 | 88.6 | — | — | — | 7.6 | 0.6 | 293 |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | 42.8 | Air/Steam | 39.1 | 11.6 | 50.7 | — | — | — | 2.9 | 0.6 | 293 |
| CASE_VTT_ATMOSPHERIC_PEAT_D | 10.7 | Air (only, no steam) | 12.9 | 0 | 12.9 | — | — | — | 0 | — | 293 |
| CASE_VTT_ATMOSPHERIC_PEAT_E | 9.8 | Air (only, no steam) | 13 | 0 | 13 | — | — | — | 0 | — | 293 |
| CASE_WSV400_A | 35.2 | Air (only, no steam) | 50 | 0 | 50 | — | — | — | 0 | — | 293 |
| CASE_WSV400_B | 44.4 | Air (only) | 50 | — | 50 | — | — | — | — | — | 293 |
| CASE_WSV400_C | 38.8 | Air (only) | 50 | — | 50 | — | — | — | — | — | 293 |
| CASE_WSV400_D | 34.5 | Air (only) | 50 | — | 50 | — | — | — | — | — | 293 |
| CASE_WSV400_E | 31.1 | Air (only) | 50 | — | 50 | — | — | — | — | — | 293 |
| CASE_WSV400_F | 33.5 | Air (only) | 50 | — | 50 | — | — | — | — | — | 293 |

## 4. 燃料工业/元素分析

| Case ID | M [%] | Ash_d [%] | VM_daf [%] | C_d [%] | H_d [%] | O_d [%] | N_d [%] | S_d [%] |
|---------|-------|-----------|------------|---------|---------|---------|---------|---------|
| CASE_HTW_WESSELING_1 | 16.9 | 11.41 | 53.42 | 61.5 | 4.1 | 21.8 | 0.68 | 0.51 |
| CASE_HTW_WESSELING_2 | 15.8 | 8.31 | 52.7 | 62.9 | 4.03 | 23.3 | 0.88 | 0.58 |
| CASE_HTW_BERRENRATH | 12.4 | 12.1 | 53.38 | 61.1 | 4.4 | 20.7 | 0.7 | 1 |
| CASE_VTT_PRESSURIZED_COAL_1 | 11.8 | 4.3 | 53 | 63.8 | 4.6 | 26.2 | 0.8 | 0.3 |
| CASE_VTT_PRESSURIZED_PEAT_6 | 14.8 | 4.3 | 68.4 | 54.5 | 5.6 | 33.6 | 1.8 | 0.25 |
| CASE_VTT_PRESSURIZED_PEAT_13 | 16 | 4.3 | 68.4 | 54.5 | 5.6 | 33.6 | 1.8 | 0.25 |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | 11.3 | 0.2 | 83 | 50.2 | 6.1 | 43.4 | 0.1 | 0 |
| CASE_VTT_ATMOSPHERIC_PEAT_D | 8.4 | 4.3 | 68.4 | 54.5 | 5.6 | 33.6 | 1.5 | 0.2 |
| CASE_VTT_ATMOSPHERIC_PEAT_E | 8.4 | 4.3 | 68.4 | 54.5 | 5.6 | 33.6 | 1.5 | 0.2 |
| CASE_WSV400_A | 27.3 | 6.34 | 93.62 | 46.27 | 5.94 | 41.01 | 0.4 | 0.01 |
| CASE_WSV400_B | 25.8 | 0.23 | 93.05 | 46.27 | 5.94 | 41.01 | 0.4 | 0.01 |
| CASE_WSV400_C | 25.8 | 0.23 | 93.05 | 46.27 | 5.94 | 41.01 | 0.4 | 0.01 |
| CASE_WSV400_D | 25.8 | 0.23 | 93.05 | 46.27 | 5.94 | 41.01 | 0.4 | 0.01 |
| CASE_WSV400_E | 25.8 | 0.23 | 93.05 | 46.27 | 5.94 | 41.01 | 0.4 | 0.01 |
| CASE_WSV400_F | 22.8 | 0.23 | 93.05 | 51.85 | 6.15 | 41.6 | 0.14 | 0.01 |

## 5. 数据来源

| Case ID | 输入来源 | 输出来源 | 轴向剖面 | Parity |
|---------|----------|----------|----------|--------|
| CASE_HTW_WESSELING_1 | Table 7.1 | Fig. 7.4, 7.5 | 有 | 无 |
| CASE_HTW_WESSELING_2 | Table 7.1 | Fig. 7.4, 7.5 | 有 | 无 |
| CASE_HTW_BERRENRATH | Table 7.1 | Fig. 7.6 | 有 | 无 |
| CASE_VTT_PRESSURIZED_COAL_1 | Table 7.2 (p.133) | Fig. 7.8 (p.135) | 有 | 有 |
| CASE_VTT_PRESSURIZED_PEAT_6 | Table 7.2 (p.133) | Fig. 7.8, 7.9, 7.10 | 有 | 有 |
| CASE_VTT_PRESSURIZED_PEAT_13 | Table 7.2 (p.133) | Fig. 7.8, 7.9, 7.10 | 无 | 有 |
| CASE_VTT_PRESSURIZED_SAWDUST_14 | Table 7.2 (p.133) | Fig. 7.8, 7.9, 7.10 | 无 | 有 |
| CASE_VTT_ATMOSPHERIC_PEAT_D | Table 7.3 (p.139) | Fig. 7.12, 7.13 | 有 | 无 |
| CASE_VTT_ATMOSPHERIC_PEAT_E | Table 7.3 (p.139) | Fig. 7.12, 7.13 | 无 | 有 |
| CASE_WSV400_A | Table 7.4 (p.143) | Fig. 7.16, 7.17 | 有 | 有 |
| CASE_WSV400_B | Table 7.4 | Fig. 7.16, 7.17 | 无 | 有 |
| CASE_WSV400_C | Table 7.4 | Fig. 7.16, 7.17 | 无 | 有 |
| CASE_WSV400_D | Table 7.4 | Fig. 7.16, 7.17 | 无 | 有 |
| CASE_WSV400_E | Table 7.4 | Fig. 7.16, 7.17 | 无 | 有 |
| CASE_WSV400_F | Table 7.4 | Fig. 7.16, 7.17 | 无 | 有 |

## 6. JSON SSOT

结构化全字段副本：[`data/validation_cases_master_table.json`](../data/validation_cases_master_table.json)

原始嵌套 SSOT：[`data/validation_cases.json`](../data/validation_cases.json)

装置原始文献 SSOT：[`data/validation_plant_literature.json`](../data/validation_plant_literature.json)


## 7. 原始文献（按装置分组）

> SSOT：[`data/validation_plant_literature.json`](../data/validation_plant_literature.json)

下列为 Hamel Ch.7 各装置实测数据的**原始出处**；Hamel (2001) VDI 469 为 Table 7.1–7.4 的二次汇编。
无公开 URL 的条目见「获取方式」。

### 7.1 `PLANT_HTW_WESSELING` — HTW Wesseling 加压中试炉

**Hamel 引用**：Table 7.1 Sim 1–2；Fig 7.4 (Sim 1 axial)；Fig 7.5 (Sim 2 axial + exit)

| 角色 | ref_id | 文献 | 链接 | 数据用途 |
|------|--------|------|------|----------|
| 实测数据（primary） | `rheinbraun_krupp_uhde_1997` | Rheinbraun AG; Krupp Uhde GmbH. *Interne Mitteilung: Messdaten aus dem Betrieb der HTW-Druckvergasungsanlage und der HTW-Demonstrationsanlage*. Internal report (cited in Hamel 2001 Literaturverzeichnis), 1997 | 非公开内部报告 | HTW Wesseling Sim 1–2、Berrenrath Sim 3：Table 7.1 输入 + Fig 7.4–7.6 出口/轴向 |
| 二次汇编（Hamel） | `hamel_2001_vdi469` | Hamel, S.. *Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe in atmosphärischen und druckaufgeladenen stationären Wirbelschichten*. Fortschr.-Ber. VDI Reihe 6 Nr. 469, VDI Verlag, Düsseldorf, 2001. ISBN 3-18-346906-5 | [链接](https://www.vdi-verlag.de/) | Chapter 7 Table 7.1–7.4、Fig 7.4–7.17 验证工况数字化来源 |
| 装置工程背景 | `adlhoch_1992a_vgb` | Adlhoch, W.; Keller, J.; Herbert, P. K.. *Das Rheinbraun-HTW-Kohlevergasungsverfahren*. VGB Kraftwerkstechnik, 1992. 72(4): 345-348 | VGB 会员/馆藏 | Wesseling Pilot 装置描述、设计意图 |
| 装置工程背景 | `adlhoch_1992b_cit` | Adlhoch, W.; Meyer, B.; Schumacher, H.-J.. *Stand der HTW-Vergasungstechnologie*. Chemie Ingenieur Technik, 1992. 64(5): 476 | [DOI](https://doi.org/10.1002/cite.330640527) | HTW 技术状态综述；Wesseling 6.5 t/h 设计能力 |
| 装置工程背景 | `schumacher_1993_braunkohle` | Schumacher, H.-J.; Meyer, B.; Brüngel, N.; Renzenbrink, W.. *Das Rheinbraun Hoch-Temperatur-Winkler (HTW)-Verfahren zur Kohleveredelung — Betriebsergebnisse der 25 bar HTW-Druckvergasungsanlage*. Braunkohle, 1993. 45(11): 12-17 | [链接](https://www.osti.gov/etdeweb/biblio/5641048) | HTW Wesseling 运行期、Agent 类型、褐煤/烟煤试验背景 |
| 装置工程背景 | `hamel_et_al_1999_iscs` | Hamel, S.; Krumm, W.; Abraham, R.; Mittelstadt, A.; Karkowski, G.. *Modelling of Pressurized Fluidized Bed Gasification in Comparison with Experimental Data from a Commercial Scale and Pilot Scale HTW-Gasification Plant*. Proceedings of 4th International Symposium on Coal Science, Beijing, 1999. pp. 411-420 | 会议论文 | HTW 示范/中试尺度模型-实测对比背景 |
| 装置工程背景 | `netl_htw_overview` | thyssenkrupp Uhde / NETL compilation. *High Temperature Winkler (HTW) Coal Gasification*. NETL Gasification Technologies Conference paper (GTC9808P), 1998 | [链接](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) | HTW 两装置运行史与规模背景（非 Table 7.1 逐点数据） |

### 7.2 `PLANT_HTW_BERRENRATH` — HTW Berrenrath 示范工厂

**Hamel 引用**：Table 7.1 Sim 3；Fig 7.6 (axial + exit)

| 角色 | ref_id | 文献 | 链接 | 数据用途 |
|------|--------|------|------|----------|
| 实测数据（primary） | `rheinbraun_krupp_uhde_1997` | Rheinbraun AG; Krupp Uhde GmbH. *Interne Mitteilung: Messdaten aus dem Betrieb der HTW-Druckvergasungsanlage und der HTW-Demonstrationsanlage*. Internal report (cited in Hamel 2001 Literaturverzeichnis), 1997 | 非公开内部报告 | HTW Wesseling Sim 1–2、Berrenrath Sim 3：Table 7.1 输入 + Fig 7.4–7.6 出口/轴向 |
| 二次汇编（Hamel） | `hamel_2001_vdi469` | Hamel, S.. *Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe in atmosphärischen und druckaufgeladenen stationären Wirbelschichten*. Fortschr.-Ber. VDI Reihe 6 Nr. 469, VDI Verlag, Düsseldorf, 2001. ISBN 3-18-346906-5 | [链接](https://www.vdi-verlag.de/) | Chapter 7 Table 7.1–7.4、Fig 7.4–7.17 验证工况数字化来源 |
| 装置工程背景 | `adlhoch_1992a_vgb` | Adlhoch, W.; Keller, J.; Herbert, P. K.. *Das Rheinbraun-HTW-Kohlevergasungsverfahren*. VGB Kraftwerkstechnik, 1992. 72(4): 345-348 | VGB 会员/馆藏 | Wesseling Pilot 装置描述、设计意图 |
| 装置工程背景 | `schumacher_1993_process` | Schumacher, H.-J.; et al.. *HTW process — secondary air staging in freeboard (cited in Hamel §7.1.2)*. See Schumacher et al. (1993) Braunkohle; also NETL HTW overview, 1993 | [链接](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) | 自由板二级风/温度控制工艺背景 |
| 装置工程背景 | `hamel_et_al_1999_iscs` | Hamel, S.; Krumm, W.; Abraham, R.; Mittelstadt, A.; Karkowski, G.. *Modelling of Pressurized Fluidized Bed Gasification in Comparison with Experimental Data from a Commercial Scale and Pilot Scale HTW-Gasification Plant*. Proceedings of 4th International Symposium on Coal Science, Beijing, 1999. pp. 411-420 | 会议论文 | HTW 示范/中试尺度模型-实测对比背景 |
| 装置工程背景 | `netl_htw_overview` | thyssenkrupp Uhde / NETL compilation. *High Temperature Winkler (HTW) Coal Gasification*. NETL Gasification Technologies Conference paper (GTC9808P), 1998 | [链接](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) | HTW 两装置运行史与规模背景（非 Table 7.1 逐点数据） |

### 7.3 `PLANT_VTT` — VTT Espoo 实验台（常压+加压）

**Hamel 引用**：Table 7.2 (pressurized Sim 1–15)；Table 7.3 (atmospheric Cases A–E)；Fig 7.7–7.10 (pressurized)；Fig 7.11–7.13 (atmospheric)

| 角色 | ref_id | 文献 | 链接 | 数据用途 |
|------|--------|------|------|----------|
| 实测数据（primary） | `kurkela_stahlberg_1992_part1` | Kurkela, E.; Ståhlberg, P.. *Air gasification of peat, wood and brown coal in a pressurized fluidized-bed reactor. I. Carbon conversion, gas yields and tar formation*. Fuel Processing Technology, 1992. 31(1): 1-21 | [DOI](https://doi.org/10.1016/0378-3820(92)90038-R) · [链接](https://www.sciencedirect.com/science/article/pii/037838209290038R) · [链接](https://cris.vtt.fi/en/publications/air-gasification-of-peat-wood-and-brown-coal-in-a-pressurized-flu) | VTT 加压 BFB Table 7.2 工况、出口气、C 转化；Fig 7.8/7.9/7.10 原始实测 |
| 实测数据（primary） | `leppalahti_kurkela_1991_fuel` | Leppälahti, J.; Kurkela, E.. *Behaviour of nitrogen compounds and tars in fluidized bed air gasification of peat*. Fuel, 1991. 70(4): 491-497 | [DOI](https://doi.org/10.1016/0016-2361(91)90026-7) · [链接](https://www.sciencedirect.com/science/article/abs/pii/0016236191900267) | VTT 常压 BFB Table 7.3 Cases A–E；NH₃/HCN、tar；Fig 7.12/7.13 原始实测 |
| 二次汇编（Hamel） | `hamel_2001_vdi469` | Hamel, S.. *Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe in atmosphärischen und druckaufgeladenen stationären Wirbelschichten*. Fortschr.-Ber. VDI Reihe 6 Nr. 469, VDI Verlag, Düsseldorf, 2001. ISBN 3-18-346906-5 | [链接](https://www.vdi-verlag.de/) | Chapter 7 Table 7.1–7.4、Fig 7.4–7.17 验证工况数字化来源 |
| 装置工程背景 | `kurkela_et_al_1993_vtt161` | Kurkela, E.; Ståhlberg, P.; Laatikainen, J.. *Pressurized fluid-bed gasification experiments with wood, peat and coal at VTT in 1991–1992. Part 1. Test facilities and gasification experiments with sawdust*. VTT Publications 161, 1993. ISBN 951-38-4399-8 | [链接](https://cris.vtt.fi/en/publications/pressurized-fluidized-bed-gasification-experiments-with-wood-peat/) · [链接](https://www.osti.gov/etdeweb/servlets/purl/593634) | VTT 加压装置工程描述（Hamel §7.1.3 引用） |
| 装置工程背景 | `kurkela_et_al_1995_vtt249` | Kurkela, E.; Laatikainen-Luntama, J.; Ståhlberg, P.. *Pressurised fluidised-bed gasification experiments with wood, peat and coal at VTT in 1991–1994. Part 2. Experiences from peat and coal gasification and hot gas filtration*. VTT Publications 249, 1995. ISBN 951-38-4949-X | [链接](https://cris.vtt.fi/en/publications/pressurised-fluidised-bed-gasification-experiments-with-wood-peat/) | VTT 加压装置补充运行经验 |
| 装置工程背景 | `kurkela_stahlberg_1992_part2` | Kurkela, E.; Ståhlberg, P.. *Air gasification of peat, wood and brown coal in a pressurized fluidized-bed reactor. II. Formation of nitrogen compounds*. Fuel Processing Technology, 1992. 31(1): 23-32 | [DOI](https://doi.org/10.1016/0378-3820(92)90039-S) · [链接](https://www.sciencedirect.com/science/article/abs/pii/037838209290039S) | VTT 加压 NH₃/HCN（Peat/Coal/Wood）；补充 Table 7.2 氮组分 |
| 装置工程背景 | `leppalahti_kurkela_1992_particle` | Leppälahti, J.; Kurkela, E.. *Particle size distribution for peat feed (cited in Hamel §7.1.4 for Table 7.3)*. See Leppälahti & Kurkela (1991) and related VTT peat gasification reports, 1992 | Hamel 引用 1992 粒径分布 | VTT 常压泥炭颗粒粒径分布 |

### 7.4 `PLANT_WSV400` — 马格德堡 WSV 400 研究炉

**数据提供**：Table 7.4 实验数据由 Prof. Dr.-Ing. habil. P. Käferstein（Universität Magdeburg）提供（Hamel §7.1.5）

**Hamel 引用**：Table 7.4 Cases A–F；Fig 7.14 (geometry)；Fig 7.16 (Case A axial)；Fig 7.17 (Cases A–F parity)

| 角色 | ref_id | 文献 | 链接 | 数据用途 |
|------|--------|------|------|----------|
| 实测数据（primary） | `hoyer_et_al_1999_bmbf` | Hoyer, L.; Härtel, J.; Käferstein, P.; Gohla, M.; Tepper, H.; Jäger, W.; Pfund, M.. *Entwicklung und Projektierung einer Wirbelschicht-Vergasungsanlage zur dezentralen energetischen Verwertung von Abfallfraktionen und biogenen Reststoffen*. Forschungsbericht BMBF (FUEGO 001620167), 1999 | 德国 BMBF 项目报告 | WSV 400 装置设计与运行描述（Hamel §7.1.5） |
| 二次汇编（Hamel） | `hamel_2001_vdi469` | Hamel, S.. *Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe in atmosphärischen und druckaufgeladenen stationären Wirbelschichten*. Fortschr.-Ber. VDI Reihe 6 Nr. 469, VDI Verlag, Düsseldorf, 2001. ISBN 3-18-346906-5 | [链接](https://www.vdi-verlag.de/) | Chapter 7 Table 7.1–7.4、Fig 7.4–7.17 验证工况数字化来源 |
| 装置工程背景 | `hamel_et_al_2000b_dgmk` | Hamel, S.; Krumm, W.; Tepper, H.; Käferstein, P.. *Theoretische und experimentelle Untersuchung des Betriebsverhaltens einer Wirbelschichtanlage zur Holzvergasung*. DGMK-Tagungsbericht 2000-1, 2000. pp. 7-14 | DGMK 会议论文 | WSV 400 木材气化实验与模型对比补充 |

### 7.5 全文献书目索引

| ref_id | 短引 | DOI | 链接 | 获取说明 |
|--------|------|-----|------|----------|
| `hamel_2001_vdi469` | Hamel (2001) | — | [链接](https://www.vdi-verlag.de/) | 本地 PDF：docs/_2001_VDI-Dis._…Brennstoffe.pdf（.gitignore）；Table 7.1–7.4 二次汇编 SSOT |
| `rheinbraun_krupp_uhde_1997` | Rheinbraun & Krupp-Uhde (1997) | — | 非公开内部报告 | 非公开内部报告；Hamel Table 7.1 实测出口/工况的直接数据来源，需通过 RWE/thyssenkrupp 档案或 Hamel 1999 合作方 |
| `hamel_et_al_1999_iscs` | Hamel et al. (1999) | — | 会议论文 | 会议论文；与 Rheinbraun 1997 内部数据及 HTW 运行经验相关，公开获取需检索 ISCS 1999 论文集 |
| `adlhoch_1992a_vgb` | Adlhoch et al. (1992a) | — | VGB 会员/馆藏 | VGB 会员/馆藏；HTW 工艺与中试/示范装置工程背景 |
| `adlhoch_1992b_cit` | Adlhoch et al. (1992b) | [10.1002/cite.330640527](https://doi.org/10.1002/cite.330640527) | [DOI](https://doi.org/10.1002/cite.330640527) | — |
| `schumacher_1993_braunkohle` | Schumacher et al. (1993) | — | [链接](https://www.osti.gov/etdeweb/biblio/5641048) | OSTI/ETDEWEB 书目记录；1989–1992 Wesseling 25 bar 中试运行总结 |
| `schumacher_1993_process` | Schumacher et al. (1993) [process] | — | [链接](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) | NETL 文档为 HTW 工艺英文综述，含 Berrenrath 示范厂运行统计 |
| `kurkela_stahlberg_1992_part1` | Kurkela & Ståhlberg (1992) I | [10.1016/0378-3820(92)90038-R](https://doi.org/10.1016/0378-3820(92)90038-R) | [DOI](https://doi.org/10.1016/0378-3820(92)90038-R) · [链接](https://www.sciencedirect.com/science/article/pii/037838209290038R) · [链接](https://cris.vtt.fi/en/publications/air-gasification-of-peat-wood-and-brown-coal-in-a-pressurized-flu) | — |
| `kurkela_stahlberg_1992_part2` | Kurkela & Ståhlberg (1992) II | [10.1016/0378-3820(92)90039-S](https://doi.org/10.1016/0378-3820(92)90039-S) | [DOI](https://doi.org/10.1016/0378-3820(92)90039-S) · [链接](https://www.sciencedirect.com/science/article/abs/pii/037838209290039S) | — |
| `kurkela_et_al_1993_vtt161` | Kurkela et al. (1993) | — | [链接](https://cris.vtt.fi/en/publications/pressurized-fluidized-bed-gasification-experiments-with-wood-peat/) · [链接](https://www.osti.gov/etdeweb/servlets/purl/593634) | VTT 报告全文可通过 VTT CRIS / OSTI 获取；装置几何与 cyclone 细节 |
| `kurkela_et_al_1995_vtt249` | Kurkela et al. (1995) | — | [链接](https://cris.vtt.fi/en/publications/pressurised-fluidised-bed-gasification-experiments-with-wood-peat/) | VTT Part 2；泥炭/褐煤扩展试验与热气体过滤 |
| `leppalahti_kurkela_1991_fuel` | Leppälahti & Kurkela (1991) | [10.1016/0016-2361(91)90026-7](https://doi.org/10.1016/0016-2361(91)90026-7) | [DOI](https://doi.org/10.1016/0016-2361(91)90026-7) · [链接](https://www.sciencedirect.com/science/article/abs/pii/0016236191900267) | — |
| `leppalahti_kurkela_1992_particle` | Leppälahti & Kurkela (1992) | — | Hamel 引用 1992 粒径分布 | Hamel 引用 1992 粒径分布；常压 Case D/E 进料筛分 |
| `hoyer_et_al_1999_bmbf` | Hoyer et al. (1999) | — | 德国 BMBF 项目报告 | 德国 BMBF 项目报告；需 DTU/机构仓储或联系 Magdeburg 课题组 |
| `hamel_et_al_2000b_dgmk` | Hamel et al. (2000b) | — | DGMK 会议论文 | DGMK 会议论文；与 WSV 400 / Käferstein 合作 |
| `netl_htw_overview` | Uhde/NETL HTW overview | — | [链接](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) | 英文公开综述；Berrenrath 67000 h、Wesseling/Berrenrath 尺度对比 |
