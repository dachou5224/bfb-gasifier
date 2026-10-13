# T3：第 7、8 章图表数字化与工况标签核对

| 字段 | 内容 |
|---|---|
| 执行者 | 辅助执行者（目标所有者指派） |
| 前置任务 | T0 合并（可与 T2 并行，有 T2 的图号索引更好） |
| 分支 | `task/T3-digitize` |
| 时间上限 | 3 个工作日 |
| 状态 | 草稿 |

## 目标

得到可靠的 Hamel **模拟值**序列（出口与轴向剖面），作为章程 §2 的主判据，并把模拟值与实验值严格分开。

## 范围内

1. 先列出第 7、8 章所有 Bild（图）的清单，标注每张图的类型（轴向剖面 / parity / 其它）和对应工况。已知至少包括 Bild 7.4–7.6（HTW 剖面）、7.8（VTT 加压剖面）、7.10（VTT 加压 parity）、7.13（VTT 常压 parity）、7.16–7.17（WSV400）、8.2–8.5（汇总 parity），以清单为准。
2. 用 WebPlotDigitizer 或同类专用工具读数，不使用 LLM 目测估读。
3. 每张图输出一个文件 `data/digitized/fig_7_X.json`，每个序列必须包含：

```json
{
  "figure": "7.4", "pdf_page": 0, "printed_page": "",
  "case_label_in_figure": "", "case_id": "CASE_...",
  "series": [
    {"kind": "sim | meas", "variable": "CO", "basis": "wet | dry", "unit": "mol/mol",
     "x_axis": {"name": "xi", "definition": "", "unit": ""},
     "points": [[0.0, 0.0]], "tool": "WebPlotDigitizer 4.x", "est_error": ""}
  ]
}
```

4. **工况标签核对**：建立 `data/digitized/case_id_map.json`，把图中的工况标签对应到 `data/validation_cases.json` 中的 `CASE_*`。已发现的疑点：
   - `hamel_fig7_parity_digitized.json` 中 `fig_7_17_wsv400` 的 A 为 T_FRB 1135 K、CO 0.172；而 `validation_cases.json` 的 `CASE_WSV400_A` 为 T_out 980 K、CO 0.152。可能是标签错配，也可能是温度定义不同（平均自由空间温度与出口温度），必须对照原图和表 7.4 裁定。
   - `fig_7_13_vtt_atmospheric` 有 A–E 五个工况，而 `validation_cases.json` 只收录了 D、E。
   - 该文件已有 `case_id_map`（11 条：VTT 加压 1/6/13/14 → Bild 7.10；WSV400 A–F → Bild 7.17；VTT 常压 E → Bild 7.13；VTT 常压 D 未映射）。这些映射尚未经过核实，须逐项复核，不能直接作为已确认的对应关系；上面 WSV400 A 的数值矛盾就出现在已映射的条目上。复核结论写入新的 `data/digitized/case_id_map.json`，旧文件保持只读、原样不动。
5. 新读数与旧文件（`hamel_fig7_parity_digitized.json`、`validation_cases.json` 中的 `axial_profiles`）逐点比较，偏差超过估计误差的点复读并说明原因。

## 范围外

- 不修改 `validation_cases.json` 的既有数值；发现错误时写进核对报告，由目标所有者决定是否修正。
- 不对缺失数据做任何推断或填补。

## 交付物

- `data/digitized/*.json`、`case_id_map.json`；
- `data/digitized/REPORT.md`：图清单、新旧差异、标签核对结论、无法读取的项。

## 验收证据

- 总司令抽查 3 张图，各 5 个点，与原图一致（在估计误差内）。
- 每个 `CASE_*` 都有明确的模拟值来源，或明确标注“原文无模拟值”。

## 停止条件

时间上限到期，或某张图无法可靠区分模拟线与实验点时，停下说明。
