# T3：第 7、8 章出口值数字化与工况标签核对

| 字段 | 内容 |
|---|---|
| 执行者 | 图清单与标签核对准备：cursor；看图读数：目标所有者（人工，用专用数字化工具） |
| 前置任务 | T0 合并（可与 T2 并行，有 T2 的图号索引更好） |
| 分支 | `task/T3-digitize` |
| 时间上限 | 3 个工作日 |
| 状态 | 进行中（2026-10-13 收窄：只做出口值；同日再收窄：不提取第 8 章汇总误差带图，HTW／Berrenrath 出口值暂缓） |

## 目标

得到可靠的 Hamel **出口模拟值**（出口或平均自由空间温度、出口气体组成、碳转化率），作为章程 §2.2 的主判据，并把模拟值与实验值严格分开。

轴向（沿炉高）剖面暂不数字化，见“暂缓项与触发条件”。

## 范围内

1. **图清单**（cursor）：列出第 7、8 章所有 Bild，标注每张图的类型（轴向剖面 / parity / 其它）、对应工况、PDF 页与印刷页。剖面图只登记，不读数。
2. **出口值读数**（目标所有者）：只读能给出**单个工况**出口模拟值的图：Bild 7.10（VTT 加压）、7.13（VTT 常压）、7.17（WSV400）。用 WebPlotDigitizer 或同类专用工具读数，不使用 LLM 目测估读。
   - **不提取**第 8 章的汇总 parity 图（Bild 8.2–8.5）：其误差带即 Hamel 自称的精度范围，已记录在 `validation_cases.json` 的 `validation_performance_summary`；图上的点与第 7 章各工况 parity 图重复，且多数无法辨认所属工况。
   - 判断是否还需读某张图的标准：某个工况的某个出口量是否已有明确的模拟值来源；已有则不读。
3. **一个剖面点的例外**：在 Bild 7.16（WSV400 A 轴向剖面）上只读出口端的温度和 CO（模拟线和实验点各一个），用于第 4 项裁定 WSV400 A 的数值矛盾。
4. **工况标签核对**（cursor 准备材料，目标所有者裁定）：建立 `data/digitized/case_id_map.json`，把图中的工况标签对应到 `data/validation_cases.json` 中的 `CASE_*`。已发现的疑点：
   - `hamel_fig7_parity_digitized.json` 中 `fig_7_17_wsv400` 的 A 为 T_FRB 1135 K、CO 0.172；而 `validation_cases.json` 的 `CASE_WSV400_A` 为 T_out 980 K、CO 0.152。可能是标签错配，也可能是温度定义不同：980 K 可能取自 Bild 7.16 剖面的出口端，1135 K 则是平均自由空间温度。须对照原图（含第 3 项读数）和表 7.4 裁定。
   - `fig_7_13_vtt_atmospheric` 有 A–E 五个工况，而 `validation_cases.json` 只收录了 D、E。
   - 该文件已有 `case_id_map`（11 条：VTT 加压 1/6/13/14 → Bild 7.10；WSV400 A–F → Bild 7.17；VTT 常压 E → Bild 7.13；VTT 常压 D 未映射）。这些映射尚未经过核实，须逐项复核，不能直接作为已确认的对应关系；上面 WSV400 A 的数值矛盾就出现在已映射的条目上。复核结论写入新的 `data/digitized/case_id_map.json`，旧文件保持只读、原样不动。
5. **新旧比较**：出口值新读数与 `hamel_fig7_parity_digitized.json` 逐点比较，偏差超过估计误差的点复读并说明原因。

每张图输出一个文件 `data/digitized/fig_X_Y.json`，按图的类型使用下面两种格式之一（总司令 2026-10-14 裁定）。两种格式的公共字段：`figure`、`pdf_page`、`printed_page`、`side`、`caption_de`、`read_date`、`tool`、`est_error`，以及可复算的标定记录 `calibration`（见下）。

**parity 图（按工况汇总）**：

```json
{
  "figure": "7.17", "pdf_page": 83, "printed_page": "146", "side": "left",
  "tool": "…", "est_error": {"T_FRB_K": 8, "CO": 0.008, "CH4": 0.004},
  "basis_gas": "dry",
  "units": {"T_FRB_K": "K", "CO": "mol/mol", "CO2": "mol/mol", "H2": "mol/mol", "CH4": "mol/mol"},
  "quantity_definition": {"T_FRB_K": "平均自由空间温度", "CO": "出口"},
  "points_format": "[实验值（x，Messung）, 模拟值（y，Sim.）]",
  "cases": [
    {"case_label_in_figure": "A", "case_id": "CASE_WSV400_A",
     "T_FRB_K": [1066, 1048], "CO": [0.172, 0.169]}
  ],
  "not_in_source": ["H2O", "X_C"],
  "calibration": {"…": "见下"}
}
```

- 每个量都是 `[实验值, 模拟值]` 一对；未读的一侧写 `null`。
- `units` 必须写明每个量的单位，不能让读者从字段名去猜。
- `est_error` 写成按量给出的数值（与 `units` 同单位），供验收时判断“边界待定”（CHARTER §2.2）。

**轴向剖面图（逐序列）**：沿用 `series` 格式，每个序列含 `kind`（`sim` / `meas`）、`variable`、`basis`、`unit`、`quantity_definition`、`points`（`[ξ, 值]`）。

**标定记录 `calibration`（两种格式都必须有）**：写明渲染方式（工具、分辨率或缩放、裁剪原点），用于标定的网格线像素位置及其对应的坐标值，以及每个读数点所取的像素坐标（例如方框四边描边中线的中点）。任何人都能据此复算读数。若用脚本读数，脚本提交到 `tools/`（按 CHARTER §8 由 codex 审查）。

## 暂缓项与触发条件

- 轴向剖面图（如 Bild 7.4–7.6、7.8、7.16 等）的完整数字化暂缓。
- HTW Wesseling 1、2 与 Berrenrath 的出口模拟值暂缓：第 7 章没有它们的 parity 图，只有轴向剖面（Bild 7.4–7.6）。推进到章程 §3 第 5、6 关时，再读剖面模拟曲线末端的出口点（做法同第 3 项的 Bild 7.16 例外读数），届时另开任务卡。
- **触发条件**：某个工况首次通过章程 §2.2 的出口验收后，数字化该工况的剖面模拟曲线，做一次定性复核（峰值位置、趋势、二次风处跳变），防止“沿程错误、出口碰巧吻合”。届时另开任务卡。

## 范围外

- 不修改 `validation_cases.json` 的既有数值；发现错误时写进核对报告，由目标所有者决定是否修正。
- 不对缺失数据做任何推断或填补。

## 交付物

- `data/digitized/*.json`、`case_id_map.json`；
- `data/digitized/REPORT.md`：图清单（含暂缓的剖面图）、新旧差异、标签核对结论、无法读取的项。
- **出口值覆盖表**（写在 `REPORT.md` 中）：行为 15 个 `CASE_*`，列为 T、CO、CO₂、H₂、H₂O、CH₄、X_C；每格填模拟值的来源（图号与 `data/digitized/` 文件名）、“原文无”或“暂缓（HTW 剖面末端）”。

## 验收证据

- 总司令抽查 2 张 parity 图，各 5 个点，与原图一致（在估计误差内）。
- 出口值覆盖表完整：WSV400 A–F、VTT 加压 1/6/13/14、VTT 常压 D/E 的每一格都有来源或“原文无”；HTW 三个工况标“暂缓”。
- VTT 常压 D 在 Bild 7.13 中的标签已核对（旧映射未包含 D）。
- WSV400 A 的 980 K / 1135 K 矛盾有明确裁定。

## 停止条件

时间上限到期，或某张图无法可靠区分模拟值与实验值时，停下说明。
