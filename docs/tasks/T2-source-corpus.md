# T2：Hamel 原文语料与索引

| 字段 | 内容 |
|---|---|
| 执行者 | 辅助执行者（目标所有者指派） |
| 前置任务 | T0 合并 |
| 分支 | `task/T2-source-corpus` |
| 时间上限 | 2 个工作日 |
| 状态 | 草稿 |

## 目标

让任何 agent 都能在几秒内按式号、表号、图号找到 Hamel 原文的页面图像和文本，从而让“先读原文”在成本上可行。

## 范围内

1. 从本地 PDF（`docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf`）逐页导出：
   - `refs/hamel/pages/NNN.png`（300 dpi，NNN 为 PDF 页序，三位补零）；
   - `refs/hamel/pages/NNN.txt`（逐页文本）。
   这两类文件只在本地保留，列入 `.gitignore`。
2. 建立 PDF 页序与印刷页码的对照表。
3. 建立 `refs/hamel/index.json`，结构如下：

```json
{
  "page_map": {"<pdf_page>": "<printed_page>"},
  "equations": {
    "2.7": {"pdf_page": 0, "printed_page": "", "section": "", "label_de": "", "short_form": ""}
  },
  "tables": {"7.1": {"pdf_page": 0, "printed_page": "", "caption_de": ""}},
  "figures": {"7.4": {"pdf_page": 0, "printed_page": "", "caption_de": ""}},
  "sections": {"2.1": {"title_de": "", "pdf_page_start": 0}},
  "cited_references": {"<Hamel 文献表中的键>": {"full_citation": "", "cited_for": ["<小节或式号>"]}}
}
```

   `short_form` 是简短的公式转录（LaTeX），只作导航，不作为最终依据。
4. 收录 Hamel 的参考文献表，并标出每篇文献在哪些小节被引用（用于执行章程 §4.2 的来源层级）。

## 范围外

- 不解释、不审计公式含义（那是 T4 的工作）。
- 不读取或参照 `docs/hamel_submodels/` 和旧代码。
- 不提交页面图像、全文或 PDF。

## 可复用的旧材料

`data/hamel_pages_png/`、`data/hamel_eq_crops/`、`data/hamel_*_extract.txt` 在归档快照中，可以作为起点，但必须重新核对，不能直接信任。

## 交付物

- 本地语料目录（不入库）；
- `refs/hamel/index.json` 与生成脚本（入库）。

## 验收证据

- 第 2、3、4、7 章的全部带编号方程、表和图都已收录，并提供计数。
- 总司令随机抽查 20 条，页码与内容全部正确。

## 停止条件

时间上限到期，或 PDF 文本层无法提取时，停下说明情况（例如需要 OCR）。
