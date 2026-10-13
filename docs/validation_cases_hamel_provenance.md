# Validation Cases — Hamel 原文溯源

> SSOT 数据：`data/validation_cases.json`  
> 表格摘录 SSOT：`data/hamel_table7_extractions.json`  
> 自动合并脚本：`scripts/enrich_validation_cases_from_hamel.py`

## 来源文献

**二次汇编（本项目 Table 7.x 数字化直接来源）**

- Hamel S (2001) VDI Fortschr.-Ber. 469, ISBN 3-18-346906-5 — Chapter 7 Table 7.1–7.4；Chapter 8 parity（Bild 8.2–8.5）
- OCR 底稿：`data/hamel_1999_extracted.txt`

**各装置实测数据的原始文献（含 DOI / 可访问链接）**

完整书目、按装置分组与获取说明见：

- [`docs/validation_cases_master_table.md`](validation_cases_master_table.md) **§7 原始文献**
- SSOT JSON：[`data/validation_plant_literature.json`](../data/validation_plant_literature.json)

| 装置组 | 主要原始数据出处 |
|--------|------------------|
| `PLANT_HTW_WESSELING` / `PLANT_HTW_BERRENRATH` | Rheinbraun & Krupp-Uhde (1997) 内部报告（非公开）；工程背景 Adlhoch (1992)、Schumacher (1993)；[NETL HTW 综述 PDF](https://netl.doe.gov/sites/default/files/netl-file/GTC9808P.pdf) |
| `PLANT_VTT` | Kurkela & Ståhlberg (1992) [Part I DOI](https://doi.org/10.1016/0378-3820(92)90038-R)、[Part II DOI](https://doi.org/10.1016/0378-3820(92)90039-S)；Leppälahti & Kurkela (1991) [Fuel DOI](https://doi.org/10.1016/0016-2361(91)90026-7)；VTT Publ. 161 [CRIS](https://cris.vtt.fi/en/publications/pressurized-fluidized-bed-gasification-experiments-with-wood-peat/) |
| `PLANT_WSV400` | Hoyer et al. (1999) BMBF 报告（需馆藏）；Hamel et al. (2000b) DGMK；数据由 Käferstein 提供 |

## 15 Case 与 Hamel 表格对应

| Case ID | Table | Sim | 装置组 | 装置 |
|---------|-------|-----|--------|------|
| `CASE_HTW_WESSELING_1` | 7.1 | 1 | `PLANT_HTW_WESSELING` | HTW Wesseling Luft/Dampf |
| `CASE_HTW_WESSELING_2` | 7.1 | 2 | `PLANT_HTW_WESSELING` | HTW Wesseling O₂/Dampf |
| `CASE_HTW_BERRENRATH` | 7.1 | 3 | `PLANT_HTW_BERRENRATH` | HTW Berrenrath O₂/Dampf |
| `CASE_VTT_PRESSURIZED_COAL_1` | 7.2 | 1 | `PLANT_VTT` | VTT 加压褐煤 |
| `CASE_VTT_PRESSURIZED_PEAT_6` | 7.2 | 6 | `PLANT_VTT` | VTT 加压泥炭 |
| `CASE_VTT_PRESSURIZED_PEAT_13` | 7.2 | 13 | `PLANT_VTT` | VTT 加压泥炭 + 循环 |
| `CASE_VTT_PRESSURIZED_SAWDUST_14` | 7.2 | 14 | `PLANT_VTT` | VTT 加压锯末 |
| `CASE_VTT_ATMOSPHERIC_PEAT_D` | 7.3 | D | `PLANT_VTT` | VTT 常压泥炭 + 循环 |
| `CASE_VTT_ATMOSPHERIC_PEAT_E` | 7.3 | E | `PLANT_VTT` | VTT 常压泥炭 + 循环 |
| `CASE_WSV400_A`–`F` | 7.4 | A–F | `PLANT_WSV400` | Magdeburg WSV 400 木材 |

**4 套物理装置**：Wesseling 中试 2 case、Berrenrath 示范 1 case、VTT 实验台 7 case（常压/加压同一装置）、WSV 6 case。详见 [`validation_cases_master_table.md`](validation_cases_master_table.md) §0 与 `validation_cases.json` → `validation_plants`。

## 已从原文补齐（inputs）

- **Table 7.1–7.4**：燃料流量、λ_O₂、空气/O₂ 分股、蒸汽 kg/h（VTT/WSV）、燃料工业/元素分析（按 sim 分行）
- **§7.1.2–7.1.5 工程描述**：反应器直径、总高、床高（VTT **1.2 m** 来自 Table 7.5；HTW Wesseling **6 m** 膨胀床；WSV **D_bed=0.4 / D_fb=0.6**）
- **二级风**：HTW ξ≈0.6；VTT 自由板 secondary air port **800 mm**  над distributor
- **T_inlet**：环境 **293 K**（Hamel 正文默认）

每个 case 节点新增：

- `hamel_simulation_id`
- `hamel_table_ref`
- `hamel_plant_engineering`（装置级工程摘要）

## 仍缺 / 仅图形可读（outputs 或未入表）

| 类别 | 说明 |
|------|------|
| HTW O₂/Steam **绝对蒸汽流量** | Table 7.1 仅给 λ_O₂ 与 agent 类型，无 kg/h 蒸汽列 |
| HTW **primary air/O₂ 分股** | 多级分布器，Table 7.1 未分列；仅 secondary 有实验值 |
| **床层惰性材料 / d_p** | 全文未给 HTW/VTT 床料粒径 |
| Berrenrath **bed_height_m** | Table 7.1 无；仅知单点进料、H=23 m |
| VTT Peat 13/14、WSV B–F **出口标量** | **Fig 7.10/7.17 已数字化**（`hamel_fig7_parity_digitized.json`） |
| VTT Atm **Case E** | **Fig 7.13 已数字化**（无 Bild 7.12 轴向剖面） |
| **NH₃/HCN ppm** | Bild 7.13 未含；需 Leppälahti (1991) 原文 |
| **H₂S** | HTW 出口未报告 |

## 维护

1. 更新 `data/hamel_table7_extractions.json`（人工核对 OCR 与 PDF）
2. 运行 `python3 scripts/enrich_validation_cases_from_hamel.py`
3. **Parity 出口**：更新 `data/hamel_fig7_parity_digitized.json` → `python3 scripts/merge_fig7_parity_outputs.py`
4. 渲染底图：`python3 scripts/render_hamel_fig7_parity_pages.py`
5. 同步 [`validation_cases_chemistry_map.md`](validation_cases_chemistry_map.md) §2 表

*2026-06 首版 enrichment：missing 标记 82→~60（inputs 为主）。*
