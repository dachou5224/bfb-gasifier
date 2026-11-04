# A0055 原文答卷：床灰的粒级、排灰与惰性料衡算

| 项 | 内容 |
|---|---|
| 对应 | ADR 0055（床灰按规定的飞灰份额分为细灰与粗灰） |
| 执行者 | 总司令（Claude）；定位与要点引自三份独立文献审计（codex、opus、deepseek，2026-11-04，报告不入库），总司令复核了逐页文本层 |
| 会话 | 非独立转录任务 |
| 日期 | 2026-11-05 |
| 页面图像 | opus 审计逐页打开图像核对了下列 Hamel 各条；总司令只核了文本层。审查方按 skill 第六节对照页面图像复核 |

Q6 共同：未查看归档代码。

## 1. Tab. 2.1、Tab. 2.2（装置与运行参数）

- **Q1**：PDF 16 右 / S. 13，Tab. 2.1、Tab. 2.2。
- **Q2**：存量（kg）、粒径（m）为输入；床高（m）为输入。
- **Q3**：„Position und Anzahl der … Ascheabzugstellen“；„Zusammensetzung und Partikelgrößenverteilung des Feststoffinventars oder Angabe eines mittleren Partikeldurchmessers“；床高 „über Druckdifferenzmessungen konstant gehalten“。
- **Q4**：床高与压差的关系见 ADR 0014 A2–A4。
- **Q5**：排灰流率、排灰组成、底灰与飞灰的分配，原文未说明。存量的粒度分布是输入，第 7 章没有报告 VTT、HTW 工况的取值。ADR 0055 第 1–3 项是解释性选择。

## 2. Bild 2.2（计算流程）

- **Q1**：PDF 17 左 / S. 14。
- **Q3**：初始化框分列 „Anzahl der Partikelgrößenklassen aller bilanzierten aktiven Feststoffe“ 与 „… des inerten Bettmaterials“；预计算框第 3 步 „Abgleich der Inertmaterialbilanz“。
- **Q5**：“惰性料衡算校正”没有方程，也没有说明是否逐粒级。ADR 0055 第 3 项给出它的内容，是解释性选择。

## 3. 式 (2.6) 与 S. 10（固体衡算对象）

- **Q1**：式 (2.6) PDF 18 右 / S. 17，符号说明 PDF 19 左 / S. 18；S. 10 在 PDF 15 左。
- **Q3**：S. 17 „… verschiedener aktiver Feststoffarten …“；S. 10 „Feststoffentnahmeströme, z.B. im Falle eines Ascheabzugs“，以及 „Partikelgrößenverteilung jedes Feststoffs“。
- **Q5**：惰性料的逐粒级衡算原文没有。S. 10 与 S. 17 的张力原文未解决。ADR 0055 不对床灰逐粒级衡算，只规定两箱分配。

## 4. 简化算例中的惰性料（S. 78–79）

- **Q1**：PDF 49 / S. 78 末–S. 79 首，§5.1.1.2；Tab. 5.1。
- **Q3**：„Der Inertgehalt sowie der Abstrom von Inertmasse ist konstant und für alle Zellen gleich“。
- **Q5**：这是简化算例，不能直接外推到主模型；只作为“惰性料作给定量、不逐粒级求解”的旁证。

## 5. 焦炭的口径（式 4.23、4.26）

- **Q1**：PDF 42 左 / S. 64。
- **Q2**：x_C,Koks = kg C im Koks / kg C im Brennstoff (waf)。
- **Q4**：式 (4.26) 的碳平衡。
- **Q5**：式 (4.23) 两符号之间的运算符，OCR 读作一点（deepseek 建议看原页复核）。按式 (4.26) 的量纲，焦炭只计碳，与 ADR 0015 第 3 项的读法一致。ADR 0055 不依赖这一读法：反应颗粒数已按 ADR 0044 只计焦炭质量。

## 6. 燃尽模型（S. 73、S. 85）

- **Q1**：PDF 46 右 / S. 73，Bild 5.1；PDF 52 / S. 84–85。
- **Q3**：S. 73：shrinking particle 时 „… führt zu einer Abnahme des Partikeldurchmessers“；„Unter Berücksichtigung eines weiteren diffusiven Widerstandes in Form der den Partikel umgebenden Ascheschicht, wird für sehr aschereiche Kohlen häufig das ‚Shrinking core‘-Modell angewendet“（总司令核对文本层）。S. 85 继续采用 shrinking particle。
- **Q5**：燃尽后灰的粒径与去向，原文未说明。所选模型里没有包住反应颗粒的灰壳，但它没有规定脱离的灰有多大，不能据此在灰粒径的写法之间作出唯一选择（PR #150 复审）。ADR 0055 的细灰进入 k = 1、粗灰进入 k = K，都是解释性选择。

## 7. 气化床的含碳量（S. 155）

- **Q1**：PDF 87 右 / S. 155，§8。
- **Q3**：„… stellt sich im stationären Betrieb ein erheblich höherer Kohlenstoffgehalt im Wirbelbett ein als es bei Verbrennungswirbelschichten der Fall ist“（总司令核对文本层）。
- **Q5**：没有数值；只作定性背景。ADR 0055 不设焦炭存量下限。

## 8. KURKELA_1992（Hamel 引用的数据来源）

- **Q1**：`docs/kurkela1992.pdf` p. 7；表 5，p. 9。
- **Q3**：„All data are from such set points where the removal rate of bottom ash was insignificant and the carbon losses were completely composed of the elutriated particles.“；表 5 工况 1 的 ash balance 0,94。
- **Q5**：p. 7 说明的是这批设定点（含工况 1）底灰排出不显著，没有给出底灰与飞灰的数值分配；0,94 是总灰闭合误差。目标所有者 2026-11-04 决定以 p. 7 作 φ 主值 1 的依据（ADR 0055 第 1 项）。

## 9. Wirsum（1998）S. 160（只作定义说明，不在 ADR 0002 的采用范围内）

- **Q1**：`refs/wirsum/pages/090.txt`，PDF 90 / S. 160。
- **Q3**：„Die Asche verbleibt … entweder als Hülle … (‚shrinking-core‘), wobei der Außendurchmesser des Partikels konstant bleibt oder wird vollständig abgetragen (‚shrinking-particle‘).“
- **Q5**：语境是单颗粒燃烧动力学的实验评估，不是床灰衡算；不作 ADR 0055 的采用依据。
