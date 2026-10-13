---
status: proposed
---

# 反应编号统一采用 Hamel 的 R1–R11

Hamel 第 5 章开头列出 11 个反应：R1 Verbrennung、R2 heterogene Wassergasreaktion、R3 Boudouard-Reaktion、R4 hydrierende Vergasung、R5 CO-Oxidation、R6 H₂-Oxidation、R7 CH₄-Oxidation、R8 homogene Wassergasreaktion、R9 Methan-Reformierung、R10 Teer-Oxidation、R11 Teer-Umsetzung（PDF 45 / S. 71）。新主线的代码、测试、文档一律使用这套编号。

归档旧代码有 6 处与之不符：旧 R3/R4 互换（加氢气化与 Boudouard）、旧 R6 为 CH₄ 氧化（Hamel R7）、旧 R7 为甲烷重整（Hamel R9）、旧 R12 为 H₂ 氧化（Hamel R6），另有一个 Hamel 列表中不存在的旧 R9（H₂S 氧化）。旧代码在注释里用“仓库编号 R12 = thesis R6”维持双重编号，交接时屡次造成误读。完整对照见 `GLOSSARY.md`；旧 R9 的去留见 `CANDIDATES.md` C4。

## 证据状态

依据 OCR 文本核对（归档 `data/hamel_1999_extracted.txt` 第 45 页）。原页目视复核待 T2 页面图像交付。
