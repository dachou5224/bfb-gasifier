---
status: proposed
---

# 式 (2.6) 第二行的 ṁ_auf,j,i+1 按粒级 k 处理

> **状态：提议中。** 本 ADR 裁定 T8 比对（`refs/hamel/transcription/ch2_ch3_1_diff.md`）“原文有歧义”第 2 条。只使用 Hamel 原文。按 CHARTER §8，由 codex 审查后总司令批准，目标所有者可否决。T11（cell 衡算残差）依赖本裁定。

## 原文情况

- 式 (2.6)（PDF 18 右 / S. 17）第二行第一项印作 ṁ_auf,j,i+1：j 与 i+1 之间没有 k。T8 执行方和 codex 审查都读出这一点（T8 冻结稿 `aa856f9` 式 2.6 的 Q2；比对稿“更正”节）。
- 同一式的其余各项都带 k：ṁ_zu,j,k,i、ṁ_rez,j,k,i、ṁ_aus,j,k,i、ṁ_auf,j,k,i、ṁ_ab,j,k,i、ṁ_ab,j,k,i−1、ṁ_r,j,k,i，以及两项粒级转移。
- S. 18（PDF 19 左）的定义列表写成 ṁ_auf,j,k,i+1：“Aufwärtsgerichteter Stoffstrom der Partikelklasse k der Feststoffsorte j im Bilanzraum i+1”。
- Bild 2.4（PDF 19 左 / S. 18）的标题是“一种固体 j 的一个粒级 k 在 cell i 中的物流衡算”，图中该入流也带 k。

## 裁定

ṁ_auf,j,i+1 视为漏印 k，按 ṁ_auf,j,k,i+1 处理：来自下方 cell i+1 的、固体 j 粒级 k 的上行入流。

依据：

1. 同页定义列表和 Bild 2.4 都带 k；
2. 式 (2.6) 是按粒级 k 逐一写出的衡算（Bild 2.4 标题），一个不分粒级的总流不能与其余按 k 写出的项逐项相加；
3. 式 (2.1) 按每个粒级计 1 个变量，衡算须逐粒级闭合。

另一种读法（故意不写 k，cell i+1 的上行流不按粒级拆分）与上述三点都不相容，不采用。

这是对原文的解释性选择（CHARTER §4.5），印刷式不改。代码中注明 `解释性选择：ADR 0008`。

## 证据状态

总司令 2026-10-15 在页面图像上核对了 PDF 18 右（S. 17）的式 (2.6) 与 PDF 19 左（S. 18）的定义列表和 Bild 2.4 标题。印刷行未见 k 一事，T8 执行方与 codex 已在 PR #36 中一致确认，不再列为字形判读事项。
