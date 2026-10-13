# T3：第 7、8 章出口值数字化与工况标签核对

| 字段 | 内容 |
|---|---|
| 执行者 | 图清单与标签核对准备：cursor；看图读数：目标所有者（人工，用专用数字化工具） |
| 前置任务 | T0 合并（可与 T2 并行，有 T2 的图号索引更好） |
| 分支 | `task/T3-digitize` |
| 时间上限 | 3 个工作日 |
| 状态 | 进行中（2026-10-13 收窄范围：先做出口值，轴向剖面暂缓） |

## 目标

得到可靠的 Hamel **出口模拟值**（出口或平均自由空间温度、出口气体组成、碳转化率），作为章程 §2.2 的主判据，并把模拟值与实验值严格分开。

轴向（沿炉高）剖面暂不数字化，见“暂缓项与触发条件”。

## 范围内

1. **图清单**（cursor）：列出第 7、8 章所有 Bild，标注每张图的类型（轴向剖面 / parity / 其它）、对应工况、PDF 页与印刷页。剖面图只登记，不读数。
2. **出口值读数**（目标所有者）：只读 parity 与出口值相关的图，已知包括 Bild 7.10（VTT 加压）、7.13（VTT 常压）、7.17（WSV400）、8.2–8.5（第 8 章汇总 parity，覆盖含 HTW 在内的全部工况），以第 1 项清单为准。用 WebPlotDigitizer 或同类专用工具读数，不使用 LLM 目测估读。
3. **一个剖面点的例外**：在 Bild 7.16（WSV400 A 轴向剖面）上只读出口端的温度和 CO（模拟线和实验点各一个），用于第 4 项裁定 WSV400 A 的数值矛盾。
4. **工况标签核对**（cursor 准备材料，目标所有者裁定）：建立 `data/digitized/case_id_map.json`，把图中的工况标签对应到 `data/validation_cases.json` 中的 `CASE_*`。已发现的疑点：
   - `hamel_fig7_parity_digitized.json` 中 `fig_7_17_wsv400` 的 A 为 T_FRB 1135 K、CO 0.172；而 `validation_cases.json` 的 `CASE_WSV400_A` 为 T_out 980 K、CO 0.152。可能是标签错配，也可能是温度定义不同：980 K 可能取自 Bild 7.16 剖面的出口端，1135 K 则是平均自由空间温度。须对照原图（含第 3 项读数）和表 7.4 裁定。
   - `fig_7_13_vtt_atmospheric` 有 A–E 五个工况，而 `validation_cases.json` 只收录了 D、E。
   - 该文件已有 `case_id_map`（11 条：VTT 加压 1/6/13/14 → Bild 7.10；WSV400 A–F → Bild 7.17；VTT 常压 E → Bild 7.13；VTT 常压 D 未映射）。这些映射尚未经过核实，须逐项复核，不能直接作为已确认的对应关系；上面 WSV400 A 的数值矛盾就出现在已映射的条目上。复核结论写入新的 `data/digitized/case_id_map.json`，旧文件保持只读、原样不动。
5. **新旧比较**：出口值新读数与 `hamel_fig7_parity_digitized.json` 逐点比较，偏差超过估计误差的点复读并说明原因。

每张图输出一个文件 `data/digitized/fig_X_Y.json`，每个序列必须包含：

```json
{
  "figure": "7.17", "pdf_page": 0, "printed_page": "",
  "case_label_in_figure": "", "case_id": "CASE_...",
  "series": [
    {"kind": "sim | meas", "variable": "CO", "basis": "wet | dry", "unit": "mol/mol",
     "quantity_definition": "出口 | 平均自由空间 | …",
     "points": [[0.0, 0.0]], "tool": "WebPlotDigitizer 4.x", "est_error": ""}
  ]
}
```

parity 图的 `points` 为 `[实验值, 模拟值]` 对，并在 `x_axis`/`y_axis` 字段中注明哪一轴是 Messung、哪一轴是 Simulation。

## 暂缓项与触发条件

- 轴向剖面图（如 Bild 7.4–7.6、7.8、7.16 等）的完整数字化暂缓。
- **触发条件**：某个工况首次通过章程 §2.2 的出口验收后，数字化该工况的剖面模拟曲线，做一次定性复核（峰值位置、趋势、二次风处跳变），防止“沿程错误、出口碰巧吻合”。届时另开任务卡。

## 范围外

- 不修改 `validation_cases.json` 的既有数值；发现错误时写进核对报告，由目标所有者决定是否修正。
- 不对缺失数据做任何推断或填补。

## 交付物

- `data/digitized/*.json`、`case_id_map.json`；
- `data/digitized/REPORT.md`：图清单（含暂缓的剖面图）、新旧差异、标签核对结论、无法读取的项。

## 验收证据

- 总司令抽查 2 张 parity 图，各 5 个点，与原图一致（在估计误差内）。
- 每个 `CASE_*` 都有明确的出口模拟值来源，或明确标注“原文无模拟值”。
- WSV400 A 的 980 K / 1135 K 矛盾有明确裁定。

## 停止条件

时间上限到期，或某张图无法可靠区分模拟值与实验值时，停下说明。
