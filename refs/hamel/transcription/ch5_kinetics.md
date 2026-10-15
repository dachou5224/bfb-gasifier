# T4 第 5 章反应动力学独立转录

- 任务卡：T4
- 执行者：cursor
- 日期：2026-10-15
- 本答卷在全新会话中完成。冻结前只接触 `docs/CHARTER.md`、`docs/tasks/T4-transcribe-kinetics.md`、`.claude/skills/hamel-source-check/SKILL.md` 与 `refs/hamel/` 页面图像及索引。未打开 `docs/hamel_submodels/`、归档代码、`docs/adr/CANDIDATES.md`、`GLOSSARY.md` 旧编号表，也未打开 `docs/adr/0001-reaction-numbering.md`。
- 反应编号以 Hamel 印刷页 71 的 R1–R11 为准。
- Q6 在本冻结稿中一律留空，比对阶段另写 `ch5_kinetics_diff.md`，不回改本文件。

## 记数与页码

本书德文排印：逗号是小数点，点是千分位。同一页上已由无歧义数字确认，例如式 (5.10) 的 `2.500` 与 `6.240 K`、式 (5.20) 的 `1.500 K` 与 `101.300 Pa`、式 (5.30) 的 `188.000 J/mol`、式 (5.17) 的 `0,6` / `0,5` / `0,33`。下文每个数字同时给出印刷字形和按这一通例读出的值。式 (5.47) 不按通例裁定，见该节 Q5。

页面图像是对开两页。页眉读数：PDF 45 为印刷页 70–71，PDF 52 为 84–85，PDF 54 为 88–89，PDF 58 为 96–97。左半页印刷页号为 \(2\times\)（PDF 页号）\(-\,20\)。下文“PDF”指图像文件号，“S.”指印刷页。

温度写成 °C 时，SI 换算为 \(T/\mathrm{K} = \vartheta/^\circ\mathrm{C} + 273{,}15\)，只在 Q2 中给出，不改公式。

## 范围

转录对象是第 5 章纳入模型的反应 R1–R4 与 R5–R11：计量、焓、速率、前因子、活化能、级数、单位、温度与压力范围、Hamel 写明采用的文献。

§5.1.1 的式 (5.1)–(5.7) 与两张同号 Tabelle 5.1（S.74 的烧损模型比较，S.79 的年龄分布算例边界）是颗粒年龄与孔模型的对照计算，不是 R1–R4 的最终速率常数。本答卷只转录 S.84–85 的采用结论，不把这些对照式写成模型输入。

---

## 1. 反应总表 R1–R11

**Q1** PDF 45 右半页，S.71，无式号。节题 „5 Modellierung der chemischen Reaktionen“。焓的共同条件印在表头。

**Q2** 表头：Reaktionsenthalpie（\(0{,}1013\,\mathrm{MPa}\)；\(298{,}15\,\mathrm{K}\)）。\(0{,}1013\,\mathrm{MPa} = 1{,}013\times 10^{5}\,\mathrm{Pa}\)。焓的印刷单位是 kJ/mol，SI 为 J/mol 时乘以 \(10^{3}\)。

| 编号 | 名称（原文） | 计量 | 焓 |
|---|---|---|---|
| R1 | Verbrennung | \(\mathrm{C}+\dfrac{1}{\Phi_c}\mathrm{O}_2 \rightarrow \left(2-\dfrac{2}{\Phi_c}\right)\mathrm{CO}+\left(\dfrac{2}{\Phi_c}-1\right)\mathrm{CO}_2\) | \(-393{,}5 < \Delta H_R < -110{,}5\,\mathrm{kJ/mol\,C}\) |
| R2 | heterogene Wassergasreaktion | \(\mathrm{C}+\mathrm{H}_2\mathrm{O}\rightarrow\mathrm{CO}+\mathrm{H}_2\) | \(+131{,}4\,\mathrm{kJ/mol\,C}\) |
| R3 | Boudouard-Reaktion | \(\mathrm{C}+\mathrm{CO}_2\rightarrow 2\,\mathrm{CO}\) | \(+172{,}5\,\mathrm{kJ/mol\,C}\) |
| R4 | hydrierende Vergasung | \(\mathrm{C}+2\,\mathrm{H}_2\rightarrow\mathrm{CH}_4\) | \(-74{,}9\,\mathrm{kJ/mol\,C}\) |
| R5 | Kohlenmonoxid-Oxidation | \(\mathrm{CO}+\dfrac{1}{2}\mathrm{O}_2\rightarrow\mathrm{CO}_2\) | \(-283\,\mathrm{kJ/mol\,CO}\) |
| R6 | Wasserstoff-Oxidation | \(\mathrm{H}_2+\dfrac{1}{2}\mathrm{O}_2\rightarrow\mathrm{H}_2\mathrm{O}\) | \(-241{,}9\,\mathrm{kJ/mol\,H}_2\) |
| R7 | Methan-Oxidation | \(\mathrm{CH}_4+\dfrac{3}{2}\mathrm{O}_2\rightarrow\mathrm{CO}+2\,\mathrm{H}_2\mathrm{O}\) | \(-519{,}4\,\mathrm{kJ/mol\,CH}_4\) |
| R8 | homogene Wassergasreaktion | \(\mathrm{CO}+\mathrm{H}_2\mathrm{O}\rightarrow\mathrm{H}_2+\mathrm{CO}_2\) | \(-41{,}1\,\mathrm{kJ/mol\,CO}\) |
| R9 | Methan-Reformierung | \(\mathrm{CH}_4+\mathrm{H}_2\mathrm{O}\rightarrow\mathrm{CO}+3\,\mathrm{H}_2\) | \(+206{,}3\,\mathrm{kJ/mol\,CH}_4\) |
| R10 | Teer-Oxidation | \(\mathrm{Teer}+\mathrm{O}_2\rightarrow\mathrm{CO}+\mathrm{H}_2\) | \(\Delta H_R<0\) |
| R11 | Teer-Umsetzung | \(\mathrm{Teer}+\mathrm{H}_2\mathrm{O}\rightarrow\mathrm{CO}+\mathrm{H}_2+\mathrm{CH}_4\) | \(\Delta H_R>0\) |

S.71 的 R8 是单向箭。可逆箭出现在 S.92 式 (5.42)，见第 12 节。R10、R11 在此处是示意计量，具体计量在式 (5.52)–(5.57)。

**Q3** „Für die Berechnung der Umsatzraten des Kokses und der gasförmigen Substanzen im Wirbelschichtreaktor werden im Modell die folgenden chemischen Reaktionsgleichungen berücksichtigt“。

**Q4** R1 的系数依赖后文 \(\Phi_c\)。其余计量在本页自洽。R10、R11 的具体物质在 §5.2.6。

**Q5** R1 的焓是区间，原文未在本页按 \(\Phi_c\) 给出函数。R10、R11 未给数值焓。各物质的摩尔质量、\(c_p\) 未在本页给出。

**Q6** 比对阶段填写。

---

## 2. 烧损模型采用结论

**Q1** PDF 52 左半页下至右半页上，S.84–85，无式号。对照图为 Bild 5.10。孔参数 \(\Psi=2\) 在 S.84 末句。

**Q2** 对 \(\Psi=2\)（原文称该参数下的烧损约对应于褐煤实测，并指向第 6 章）：无返料时，孔模型相对 shrinking particle 的偏差约从 1,5 % 到 3,8 %。有固体返料时，最下层计算单元约 6,5 %，最上层约 3,5 %。原文把由烧损表述引起的偏差概括为约 10 %。这些百分数是对照计算的结果，不是速率常数。

**Q3** „Aus diesem Grund wird in dieser Arbeit das Abbrandverhalten zunächst weiterhin mit einem „Shrinking particle“-Ansatz beschrieben.“

**Q4** 该句接在孔模型与 shrinking particle 的碳转化对照之后，对照对象包括气化反应。§5.1.1 的年龄分布与孔模型公式不是这句话所采用的速率式。

**Q5** „zunächst weiterhin“ 未说明此后是否改用其它烧损模型；第 5 章后文没有另一句改口。R1 的速率用外表面积 \(A_P\) 与恒定焦炭密度，R2、R4 用 \(A_C\)，R3 用质量。\(A_C\) 是外表面、内表面还是当前 shrinking-particle 外表面，第 5 章速率式旁没有定义。shrinking density 只作为文献中气化常用的写法被提到，不是本句采用的模型。

**Q6** 比对阶段填写。

---

## 3. R1 计量与机理因子，式 (5.8)–(5.10)

**Q1** PDF 53 左半页，S.86，§5.1.2 Verbrennung des Kokses。文献：Hamel 引用的 Field et al. (1967)、Ross und Davidson (1981)、Arthur (1951)（S.86）。

**Q2**

式 (5.8)，与 S.71 的 R1 相同，下标印刷为 \(\Phi_c\)：

\[
\mathrm{C}+\frac{1}{\Phi_c}\mathrm{O}_2 \rightarrow \left(2-\frac{2}{\Phi_c}\right)\mathrm{CO}+\left(\frac{2}{\Phi_c}-1\right)\mathrm{CO}_2
\]

式 (5.9)，等号左侧印刷为 \(\Phi_C\)：

\[
\Phi_C =
\begin{cases}
\dfrac{2p_c+2}{p_c+2} & d_c \le 0{,}05\cdot 10^{-3}\,\mathrm{m} \\[1em]
\dfrac{2p_c+2-\dfrac{p_c}{0{,}95}\left(10^{3}\cdot d_c-0{,}05\right)}{p_c+2} & 0{,}05\cdot 10^{-3}\,\mathrm{m} \le d_c \le 1\cdot 10^{-3}\,\mathrm{m} \\[1em]
1 & d_c > 1\cdot 10^{-3}\,\mathrm{m}
\end{cases}
\]

\(0{,}05\cdot 10^{-3}\,\mathrm{m} = 5\times 10^{-5}\,\mathrm{m}\)，\(1\cdot 10^{-3}\,\mathrm{m} = 10^{-3}\,\mathrm{m}\)。\(0{,}95\) 无量纲。\(d_c\) 的单位在不等式中印为 m。

式 (5.10)：

\[
p_c = 2.500\cdot\exp\left(-\frac{6.240\,\mathrm{K}}{T_c}\right)
\]

`2.500` 读作 2500，无量纲。`6.240 K` 读作 \(6240\,\mathrm{K}\)。指数中的负号印在分子之前。\(T_c\) 的单位由分母与 K 相除而应为 K，原文未在符号旁再写单位。

**Q3** „Nach Field et al. (1967) und unter Einbeziehung der Ergebnisse von Ross und Davidson (1981) wird der Mechanismusfaktor durch den Ansatz“；\(p_c\) 的句子指向 Arthur (1951)。

**Q4** \(\Phi_C\) 进入后文扩散项 (5.12)、(5.15)。\(p_c\) 只在 (5.9) 中定义。

**Q5** (5.8) 的下标是 \(c\)，(5.9) 等号左侧是 \(C\)。\(T_c\) 与后文颗粒温度 \(T_P\) 是否为同一温度，原文未说明。\(d_c\) 与后文 \(d_P\) 是否为同一直径，原文未说明。

**Q6** 比对阶段填写。

---

## 4. R1 速率与传质，式 (5.11)–(5.22)

**Q1** PDF 53 左半页 S.86 至右半页 S.87。式 (5.16) 没有独立的索引框，印在 (5.15) 与 (5.17) 之间。Hamel 引用的 Ross et al. (1981) 给出发射率范围；LaNauze (1985) 是颗粒温度用表面能量平衡的建议来源。

**Q2**

恒定焦炭密度：

\[
R_C = \frac{\rho_C}{M_C}\cdot\frac{\mathrm{d}V_{P,C}}{\mathrm{d}t}
\tag{5.11}
\]

\[
R_C = -\Phi_C\cdot A_P\cdot k_{d,g}\cdot\left(C_{\mathrm{O}_2,\infty}-C_{\mathrm{O}_2,P}\right)
\tag{5.12}
\]

\[
R_C = -A_P\cdot k_{ch,R1}\cdot C_{\mathrm{O}_2,P}^{n}
\tag{5.13}
\]

正文在 (5.13) 与 (5.14) 之间写反应级数 \(n=1\)。

\[
R_C = -A_P\cdot k_1\cdot C_{\mathrm{O}_2,\infty}
\tag{5.14}
\]

\[
k_1 = \left(\frac{1}{\Phi_C\cdot k_{d,g}}+\frac{1}{k_{ch,R1}}\right)^{-1}
\tag{5.15}
\]

\[
k_{d,g} = \frac{Sh\cdot D_g}{d_P}
\tag{5.16}
\]

\[
Sh = 2 + 0{,}6\cdot Re^{0{,}5}\cdot Sc^{0{,}33}
\tag{5.17}
\]

指数是 \(0{,}33\)，不是分数 \(1/3\)。

\[
Re = \frac{u_r\cdot d_P}{v_g}
\tag{5.18}
\]

分母按放大图记为拉丁字母 \(v_g\)，见 Q5。

\[
Sc = \frac{\mu_g}{\rho_g\cdot D_g}
\tag{5.19}
\]

\[
D_g = 3{,}13\cdot 10^{-4}\left(\frac{\mathrm{m}^{2}}{\mathrm{s}}\right)\cdot\left(\frac{T_m}{1.500\,\mathrm{K}}\right)^{1{,}75}\cdot\left(\frac{101.300\,\mathrm{Pa}}{P}\right)
\tag{5.20}
\]

`1.500 K` 读作 \(1500\,\mathrm{K}\)。`101.300 Pa` 读作 \(101300\,\mathrm{Pa}\)。\(1{,}75\) 只作用在温度比上。\(3{,}13\cdot 10^{-4}\,\mathrm{m}^{2}/\mathrm{s}\) 是印刷值。

\[
T_m = \frac{T_\infty+T_P}{2}
\tag{5.21}
\]

\[
R_C\cdot\Delta h_{R_C} = A_P\cdot\alpha_{P,g}(T_P-T_\infty)+A_P\cdot\sigma\cdot\varepsilon_P\left(T_P^{4}-T_\infty^{4}\right)
\tag{5.22}
\]

颗粒内部温度梯度忽略，定常。发射率：\(0{,}79\le\varepsilon_P\le 0{,}91\)，无量纲。原文给范围，未选单值。\(\sigma\)、\(\alpha_{P,g}\) 未给数值。

\(R_C\) 在 (5.11) 中由 \(\rho_C/M_C\) 与体积变化率相乘；若 \(\rho_C\) 为 \(\mathrm{kg/m}^{3}\)、\(M_C\) 为 \(\mathrm{kg/mol}\)，则 \(R_C\) 为 \(\mathrm{mol/s}\)。原文未在 (5.11) 旁印出这一单位句。\(C\) 的单位未印。\(k_{d,g}\) 由 (5.16) 为 m/s。\(n=1\) 时 \(k_{ch,R1}\) 与 \(k_1\) 同为 m/s。

**Q3** „Mit einer Reaktionsordnung von n = 1 und unter der Annahme stationärer Bedingungen kann die unbekannte Sauerstoffkonzentration an der Koksoberfläche eliminiert werden“。发射率句：„Für die Partikelemissivität geben Ross et al. (1981) ausgehend von Messungen verschiedener Autoren einen Wert im Bereich \(0{,}79\le\varepsilon_P\le 0{,}91\) an.“

**Q4** (5.15) 依赖 (5.9) 的 \(\Phi_C\) 与 (5.23) 的 \(k_{ch,R1}\)。\(D_g\) 依赖 (5.21)。\(Sh\) 依赖 (5.18)、(5.19)。

**Q5**

- (5.14)、(5.15) 的 \(k\) 下标是单独一笔，字形在数字 1 与字母 I 之间。本答卷按合并系数记为 \(k_1\)，待人工辨认。
- (5.18) 分母待人工辨认：放大图为拉丁字母 \(v_g\)。若该符号本应是运动粘度 \(\nu_g\)，原文在本页没有定义它，也没有写出单位。
- \(\varepsilon_P\) 未取范围内的单值。\(\sigma\) 未给数值。\(\alpha_{P,g}\) 未给关联式。
- \(M_C\)、\(\rho_C\) 的数值未在本章给出。
- 浓度 \(C\) 的单位未印。

**Q6** 比对阶段填写。

---

## 5. R1 化学速率常数，式 (5.23) 与表 5.2

**Q1** PDF 54 左半页，S.88。Hamel 引用的 Hobbs et al. (1992)（S.88，式 5.23，Tabelle 5.2）。Hannes (1996) 只被写成综述出处，不是采用的参数来源。

**Q2**

\[
k_{ch,R1} = k_{0,R1}\cdot T\cdot\exp\left(-\frac{E}{R_g\cdot T}\right)
\tag{5.23}
\]

Tabelle 5.2：Zusammenstellung reaktionskinetischer Parameter der Restkoksverbrennung für verschiedene Kohlen nach Hobbs et al. (1992)。\(k_{0,R1}\) 的单位印为 \(\mathrm{m/(s\,K)}\)，\(E/R_g\) 的单位印为 K。

| Kohle | \(k_{0,R1}\) 印刷 | SI | \(E/R_g\) 印刷 | SI |
|---|---|---|---|---|
| HVBA | 1,03 | \(1{,}03\,\mathrm{m/(s\cdot K)}\) | 9.010 | \(9010\,\mathrm{K}\) |
| HVBC | 0,5 | \(0{,}5\,\mathrm{m/(s\cdot K)}\) | 6.310 | \(6310\,\mathrm{K}\) |
| SUBC | 10,4 | \(10{,}4\,\mathrm{m/(s\cdot K)}\) | 11.200 | \(11200\,\mathrm{K}\) |
| Lignite A | 1,22 | \(1{,}22\,\mathrm{m/(s\cdot K)}\) | 10.300 | \(10300\,\mathrm{K}\) |
| alle Kohlen | 2,3 | \(2{,}3\,\mathrm{m/(s\cdot K)}\) | 11.100 | \(11100\,\mathrm{K}\) |

因此 \(k_{ch,R1}\) 的单位是 m/s。\(R_g\) 未给数值；表中直接给 \(E/R_g\)。

**Q3** „Die chemische Reaktionsgeschwindigkeitskonstante wird nach Hobbs et al. (1992) berechnet“。„Im Falle der Rheinischen Braunkohle werden die kinetischen Parameter der Kohleart Lignite A verwendet.“

**Q4** 代入 (5.13)、(5.15)。莱茵褐煤采用 Lignite A 这一行。

**Q5** „alle Kohlen“ 一行与分煤种行的使用关系未说明。除莱茵褐煤以外，原文未规定其它燃料选哪一行。温度、压力适用范围未给。\(T\) 是气体温度还是颗粒温度，式旁未说明。

**Q6** 比对阶段填写。

---

## 6. R2，式 (5.24)、(5.25)

**Q1** PDF 54 左半页，S.88，§5.1.3 Heterogene Wassergasreaktion。Hamel 引用的 Hobbs et al. (1992)。计量与焓同 S.71 的 R2。

**Q2**

\[
R_C = A_C\cdot k_{ch,R2}\cdot C_{\mathrm{H}_2\mathrm{O}}
\tag{5.24}
\]

\[
k_{ch,R2} = 3{,}42\cdot T\cdot\exp(-15.600/T)
\tag{5.25}
\]

`3,42` 读作 3,42。`15.600` 读作 15600，与 \(T\) 相除，按 K 理解。对 \(\mathrm{H}_2\mathrm{O}\) 的幂是 1。没有其它反应物的浓度因子。

若 \(A_C\) 为 \(\mathrm{m}^{2}\)、\(C\) 为 \(\mathrm{mol/m}^{3}\)、\(R_C\) 为 \(\mathrm{mol/s}\)，则 \(k_{ch,R2}\) 为 m/s，`3,42` 的单位为 \(\mathrm{m/(s\cdot K)}\)。这是量纲反推，原文未印。

**Q3** „wird nach Hobbs et al. (1992) berechnet durch“。

**Q4** \(A_C\) 未在本章定义，见第 2 节。

**Q5** 前因子单位、\(A_C\) 的几何意义、\(T\) 是哪一个温度、温度与压力适用范围，原文都未写。级数没有用 „Ordnung“ 这个词写出，只体现在 \(C_{\mathrm{H}_2\mathrm{O}}\) 的一次幂上。

**Q6** 比对阶段填写。

---

## 7. R3，式 (5.26)–(5.33)

**Q1** PDF 54 右半页，S.89，§5.1.4 Boudouard-Reaktion。Hamel 引用的 Weeda (1995)。计量与焓同 S.71 的 R3。

**Q2**

\[
R^{*}_{C,0} = \frac{k_a\cdot P_{\mathrm{CO}_2}}{1+k_b\cdot P_{\mathrm{CO}_2}+k_c\cdot P_{\mathrm{CO}}}
\quad\text{in }1/\mathrm{s}
\tag{5.26}
\]

分母第一项是数字 1。

\[
\begin{align*}
k_{a,0} &= 2{,}72\,(\mathrm{Pa\cdot s})^{-1} \tag{5.27}\\
k_{b,0} &= 8{,}31\cdot 10^{-7}\,\mathrm{Pa}^{-1} \tag{5.28}\\
k_{c,0} &= 2{,}088\cdot 10^{-14}\,\mathrm{Pa}^{-1} \tag{5.29}
\end{align*}
\]

\[
\begin{align*}
E_a &= 188.000\,\mathrm{J/mol} = 188000\,\mathrm{J/mol} \tag{5.30}\\
E_b &= -38.900\,\mathrm{J/mol} = -38900\,\mathrm{J/mol} \tag{5.31}\\
E_c &= -24.500\,\mathrm{J/mol} = -24500\,\mathrm{J/mol} \tag{5.32}
\end{align*}
\]

负活化能照印。

\[
R_{C,0} = R^{*}_{C,0}\cdot\frac{m_{j,k}}{M_C}
\quad\text{in mol C/s}
\tag{5.33}
\]

\(j\) 是焦炭种类，\(k\) 是当前计算的粒径级。压力单位由 \(k_{a,0}\) 等的 Pa 确定。这不是幂律级数。

**Q3** „Für die Frequenzfaktoren \(k_a\) bis \(k_c\), deren Temperaturabhängigkeit jeweils nach Arrhenius beschrieben wird, sind nach den experimentellen Befunden, die mittels linearer Regression ermittelt wurden, folgende Werte von Weeda (1995) einzusetzen“。„Die von Weeda (1995) ermittelte Reaktionsgeschwindigkeit ist auf die Anfangsmasse der von ihm untersuchten Partikel bezogen. Der Einfluss der Partikelgröße wird vernachlässigt“.

**Q4** (5.33) 把 (5.26) 的比速率乘上 \(m_{j,k}/M_C\)。\(M_C\) 数值本章未给。

**Q5** 原文写温度关系 „nach Arrhenius“，但没有写出 \(k=k_0\exp(-E/(RT))\)，没有给出通用气体常数的符号和数值，没有说明 \(T\) 是颗粒温度还是气体温度。\(E_b\)、\(E_c\) 为负，照录，原文未另作解释。\(m_{j,k}\) 在公式里是一个符号；正文同时说 Weeda 的比速率相对其所研究颗粒的初始质量，而模型用于当前计算的粒径级。不把 \(m_{j,k}\) 改写成初始质量或当前质量。粒径影响被明确忽略，与 shrinking particle 是否改变 R3 的面积，原文没有再连接。温度适用范围未给。

**Q6** 比对阶段填写。

---

## 8. R4，式 (5.34)、(5.35)

**Q1** PDF 54 右半页，S.89，§5.1.5 Hydrierende Vergasung。Hamel 引用的 Hobbs et al. (1992)。

**Q2** 计量重复为 \(\mathrm{C}+2\,\mathrm{H}_2\rightarrow\mathrm{CH}_4\)，\(-74{,}9\,\mathrm{kJ/mol\,C}\)。

\[
R_C = A_C\cdot k_{ch,R4}\cdot C_{\mathrm{H}_2}
\tag{5.34}
\]

\[
k_{ch,R4} = 3{,}42\cdot 10^{-3}\cdot T\cdot\exp(-15.600/T)
\tag{5.35}
\]

`3,42` 读作 3,42。`15.600` 读作 15600 K。对 \(\mathrm{H}_2\) 的幂是 1。计量系数是 2，速率式不是对 \(\mathrm{H}_2\) 的二次。与 (5.25) 的差别是前因子多 \(10^{-3}\)。单位与 R2 同样只能由量纲反推为：\(k\) 若为 m/s，则前因子为 \(3{,}42\times 10^{-3}\,\mathrm{m/(s\cdot K)}\)。原文未印。

**Q3** „ist die Reaktionsgeschwindigkeitskonstante nach Hobbs et al. (1992) um drei Größenordnungen geringer als für die heterogene Wassergasreaktion.“

**Q4** 与 (5.25) 同形。\(A_C\) 仍无定义。

**Q5** 同 R2：单位、\(A_C\)、温度符号、适用范围未印。正文的 „drei Größenordnungen“ 与 \(10^{-3}\) 数值一致，本答卷仍按式面抄写，不用这句话改数。

**Q6** 比对阶段填写。

---

## 9. 均相反应所在的相

**Q1** PDF 55 左半页，S.90，§5.2 起段。Dennis et al. (1982)、van der Vaart (1988) 是对悬浮相气相燃烧变慢的解释，不是速率参数来源。

**Q2** 无新常数。

**Q3** „Während die heterogenen Umsatzreaktionen des Kokses per Definition nur in der Suspensionsphase berücksichtigt werden, können die Gasphasenreaktionen sowohl im Gasraum der Suspensionsphase als auch in der Blase auftreten.“ „Im Gegensatz dazu katalysiert der Feststoff die homogene Wassergasreaktion und die Methanreformierung innerhalb der Wirbelschicht, deren rein homogene Reaktionsraten dagegen gering sind.“

**Q4** R1–R4 只在悬浮相。R5–R11 可以在悬浮相气体与气泡中发生，后文对 R5、R9 和焦油热解又有更窄的采用句。

**Q5** „gering“ 没有定量。固体催化水煤气与甲烷重整，但本段没有指定用哪一个速率式。

**Q6** 比对阶段填写。

---

## 10. R5，式 (5.36)、(5.37) 与表 5.3

**Q1** PDF 55 左半页 S.90 至右半页 S.91，§5.2.1。Hamel 引用的 Hayhurst und Tucker (1990) 用于悬浮相，Hottel et al. (1965) 用于无固体气泡相。Hannes (1996) 只是综述。

**Q2**

\[
\frac{\mathrm{d}C_{\mathrm{CO}}}{\mathrm{d}t} = k_{ch,R5}\cdot C_{\mathrm{CO}}^{\alpha}\cdot C_{\mathrm{O}_2}^{\beta}\cdot C_{\mathrm{H}_2\mathrm{O}}^{\gamma}
\tag{5.36}
\]

\[
k_{ch,R5} = k_{0,R5}\cdot\exp\left(-\frac{E}{R_g\cdot T}\right)
\tag{5.37}
\]

没有额外的 \(T\) 乘子。Tabelle 5.3 的表题写明这些参数用于 Simulationsmodell。\(k_{0,R5}\) 的单位印为 \((\mathrm{mol}/\mathrm{m}^{3})^{1-\alpha-\beta-\gamma}\cdot 1/\mathrm{s}\)，\(E/R_g\) 的单位为 K。

| 相 | \(k_{0,R5}\) | \(E/R_g\) | \(\alpha\) | \(\beta\) | \(\gamma\) | 来源 |
|---|---|---|---|---|---|---|
| Blase | \(1{,}91\cdot 10^{6}\) | 8.056，即 \(8056\,\mathrm{K}\) | 1 | 0,3 | 0,5 | Hottel et al. 1965 |
| Suspension | \(1{,}8\cdot T^{0{,}5}\) | 短横线 | 1 | 0,5 | 短横线 | Hayhurst und Tucker 1990 |

气泡相的单位指数 \(1-\alpha-\beta-\gamma = 1-1-0{,}3-0{,}5 = -0{,}8\)。悬浮相的 \(k_{0,R5}\) 单元格已含 \(T^{0{,}5}\)，不是纯数字。

**Q3** „In dieser Arbeit findet der Ansatz von Hayhurst und Tucker (1990) für die Oxidation innerhalb der Suspensionsphase und der Ansatz von Hottel et al. (1965) für die Reaktion in der feststofffreien Blasenphase Verwendung.“

**Q4** 表 5.3 代入 (5.36)、(5.37)。

**Q5**

- (5.36) 写成浓度导数等于正的乘积，没有前导负号。原文未说明这是 CO 的净生成还是消耗。
- 悬浮相的 \(E/R_g\) 与 \(\gamma\) 都是短横线，不是数字 0。水的浓度项是省略，还是 \(\gamma=0\)；指数项是省略，还是 \(E=0\)，原文未写。
- 悬浮相的 \(k_{0,R5}\) 已含 \(T^{0{,}5}\)。把它代入 (5.37) 时是否还乘 \(\exp\)，原文未说明。
- 正文指出火焰高温动力学不能随意用于流化床低温，但没有给出数值温度区间。
- \(R_g\) 未给数值。气泡相因为表中已有 \(E/R_g\)，计算时不必单独知道 \(R_g\)。

**Q6** 比对阶段填写。

---

## 11. R6，式 (5.38)、(5.39)

**Q1** PDF 55 右半页，S.91，§5.2.2。Hamel 引用的 Jensen et al. (1995)，原文用 „in Anlehnung an“。

**Q2**

\[
\frac{\mathrm{d}C_{\mathrm{H}_2}}{\mathrm{d}t} = k_{ch,R6}\cdot C_{\mathrm{H}_2}^{3/2}\cdot C_{\mathrm{O}_2}
\tag{5.38}
\]

\(\mathrm{O}_2\) 的幂是 1。

\[
k_{ch,R6} = 1{,}631\cdot 10^{9}\cdot T^{3/2}\cdot\exp\left(-\frac{3.420}{T}\right)
\tag{5.39}
\]

`1,631` 读作 1,631。`3.420` 读作 3420，与 \(T\) 相除。单位原文未给。符号问题与 R5 相同：左端是 \(\mathrm{d}C/\mathrm{d}t\)，没有前导负号。

**Q3** „wird in Anlehnung an Jensen et al. (1995) folgendermaßen berechnet.“

**Q4** 无表。相的分配只有第 9 节的总述，本节没有再分成气泡与悬浮两套常数。

**Q5** 单位、温度与压力范围、\(T\) 是哪一个温度、速率的正负号含义，原文未写。

**Q6** 比对阶段填写。

---

## 12. R7，式 (5.40)、(5.41)

**Q1** PDF 55 右半页，S.91，§5.2.3。Hamel 引用的 de Souza-Santos (1989)。

**Q2**

\[
\frac{\mathrm{d}C_{\mathrm{CH}_4}}{\mathrm{d}t} = k_{ch,R7}\cdot C_{\mathrm{CH}_4}\cdot C_{\mathrm{O}_2}
\tag{5.40}
\]

两级数都是 1。

\[
k_{ch,R7} = \frac{3{,}552\cdot 10^{11}}{T}\cdot\exp\left(-\frac{15.700}{T}\right)
\tag{5.41}
\]

`3,552` 读作 3,552。`15.700` 读作 15700。单位原文未给。没有前导负号。

**Q3** 本节的采用句在式前，指向 de Souza-Santos (1989)。本答卷不另补原文没有的 “wird verwendet” 措辞。

**Q4** 同 R6，没有分相常数。

**Q5** 单位、适用范围、温度符号、正负号含义未写。

**Q6** 比对阶段填写。

---

## 13. R8 平衡与两套速率，式 (5.42)–(5.49)

**Q1** PDF 56–57，S.92–94，§5.2.4。式 (5.42) 在 S.92。式 (5.45)–(5.48) 在 PDF 57 左半页 S.94。式 (5.49) 在同一页下半，紧挨 §5.2.5 标题之前。

**Q2**

式 (5.42)：\(\mathrm{CO}+\mathrm{H}_2\mathrm{O}\rightleftharpoons\mathrm{H}_2+\mathrm{CO}_2\)，\(-41{,}1\,\mathrm{kJ/mol\,CO}\)。此处是可逆箭。

\[
K_{eq}=K_p=\exp\left(-\frac{\sum\nu_i\mu_i^{0}}{RT}\right)
\tag{5.43}
\]

理想气体。\(R\)、\(\mu_i^{0}\) 的数值与单位，以及求和包含哪些物质，式旁未给。

\[
K_{exp}=\frac{P_{\mathrm{CO}_2}\cdot P_{\mathrm{H}_2}}{P_{\mathrm{CO}}\cdot P_{\mathrm{H}_2\mathrm{O}}}
\tag{5.44}
\]

Wen und Chaung (1979) 的速率，单位 mol CO/(s·kg Asche)。他们把指前因子解释为固体催化活性。式 (5.48) 的 \(a_{R8}=0{,}02\) 是他们自己的气流床煤气化模拟所用的值。

式 (5.45)，无前导负号：

\[
\begin{align*}
R_{\mathrm{CO}} ={}& a_{R8}\cdot 2{,}77\cdot 10^{8}\cdot(y_{\mathrm{CO}}-y_{eq,\mathrm{CO}})\\
&\cdot\exp\left(-\frac{13.971}{T}\right)\cdot p^{\left(0{,}5-p/250\right)}\cdot\exp\left(-8{,}91+\frac{5.553}{T}\right)
\end{align*}
\]

`2,77` 读作 2,77。`13.971` 读作 13971。`5.553` 读作 5553。`8,91` 是逗号小数，读作 8,91。`0,5` 读作 0,5。`250` 无单位。\(p\) 的单位未印。

\[
y_{eq,\mathrm{CO}}=\frac{y_{\mathrm{CO}_2}\cdot y_{\mathrm{H}_2}}{k_{eq,R8}\cdot y_{\mathrm{H}_2\mathrm{O}}}
\tag{5.46}
\]

\[
k_{eq,R8}=\exp\left(-36.893+\frac{4.019}{T_g}\right)
\tag{5.47}
\]

点号不是逗号。按千分位读作 \(\exp(-36893+4019/T_g)\)。若把点号改读成小数点，则为 \(\exp(-36{,}893+4{,}019/T_g)\)。两种读法都记在 Q5，本答卷不选取。括号结构是 \(\exp(-36.893+4.019/T_g)\)，\(-36.893\) 不在分母内。

\[
a_{R8}=0{,}02
\tag{5.48}
\]

Chen et al. (1987)，煤灰流化床、大气压，原文写已换算为 SI：

\[
R_{\mathrm{CO}}=-8{,}82\cdot 10^{-7}\cdot\exp\left(\frac{-90.853}{R_g\cdot T_g}\right)\cdot\left[P_{\mathrm{CO}}\cdot P_{\mathrm{H}_2\mathrm{O}}-\frac{P_{\mathrm{H}_2}\cdot P_{\mathrm{CO}_2}}{k_{eq,R8}}\right]
\tag{5.49}
\]

单位 \(\mathrm{mol/(s\cdot kg_{Asche})}\)。`8,82` 读作 8,82。`90.853` 按千分位读作 90853。数字旁没有印 J/mol。前导负号存在。温度符号是 \(T_g\)。

**Q3** (5.48) 的句子是 „wenden Wen und Chaung für die eigene Simulation eines Flugstrom-Kohlevergasers \(a_{R8}=0{,}02\)“。 (5.49) 之前：„schlagen Chen et al. (1987) folgenden reaktionskinetischen Ansatz vor. Umgerechnet auf SI-Einheiten ergibt sich“。

**Q4** (5.46) 使用 (5.47)。(5.49) 的逆反应项也使用 \(k_{eq,R8}\)。(5.45) 使用 \(a_{R8}\) 与 \(y_{eq,\mathrm{CO}}\)。

**Q5**

- Hamel 没有写本模型采用 (5.45) 还是 (5.49)，也没有写 \(a_{R8}\) 取 0,02。0,02 明确属于 Wen und Chaung 自己的气流床模拟。
- (5.45) 在 \(y_{\mathrm{CO}}>y_{eq,\mathrm{CO}}\) 时为正，未说明是 CO 消耗速率还是浓度导数。(5.49) 有前导负号，括号内为正向减逆向。
- \(k_{eq,R8}\) 同时进入摩尔分数 (5.46) 与分压 (5.49)。\(\Delta n=0\) 时 \(K_p\) 与 \(K_y\) 数值相同，但 (5.47) 的两种点号读法都使平衡常数在气化温度下极小。这是原文歧义，不改式。
- \(p/250\) 的 \(p\) 与 250 的单位未说明。
- (5.43) 的 \(R\) 与 (5.49) 的 \(R_g\) 都未给数值。(5.49) 的 90853 若除以 \(R_g T\)，单位取决于所取的 \(R_g\)；原文只说换算到 SI，没有把 J/mol 印在该数字旁。
- Chen 的水煤气式未给温度区间。
- Bild 5.11 比较 \(K_{eq}\) 与 \(K_{exp}\)，不是速率常数，不抄成参数。

**Q6** 比对阶段填写。

---

## 14. R9，式 (5.50)、(5.51)

**Q1** PDF 57，S.94–95，§5.2.5。Hamel 引用的 Chen et al. (1987)。

**Q2**

\[
\mathrm{CH}_4+\mathrm{H}_2\mathrm{O}\rightarrow\mathrm{CO}+3\,\mathrm{H}_2,\qquad +206{,}3\,\mathrm{kJ/mol\,CH}_4
\tag{5.50}
\]

单向箭。大气压流化床，\(770\,^\circ\mathrm{C}\) 至 \(920\,^\circ\mathrm{C}\)，即 \(1043{,}15\,\mathrm{K}\) 至 \(1193{,}15\,\mathrm{K}\)。

\[
R_{\mathrm{CH}_4}=-6{,}113\cdot 10^{-2}\cdot\exp\left(\frac{-137.327}{R_g\cdot T}\right)\cdot P_{\mathrm{CH}_4}
\tag{5.51}
\]

单位 \(\mathrm{mol/(s\cdot kg_{Asche})}\)。`6,113` 读作 6,113，故前因子为 \(6{,}113\times 10^{-2}\)。`137.327` 按千分位读作 137327，单位未印在数字旁。只有 \(P_{\mathrm{CH}_4}\)，对 \(\mathrm{CH}_4\) 为一次，对 \(\mathrm{H}_2\mathrm{O}\) 没有写出级数。前导负号存在。温度符号是 \(T\)，不是 (5.49) 的 \(T_g\)。

**Q3** „Chen et al. (1987) schlagen folgende Reaktionskinetik, die in einer atmosphärischen Wirbelschicht in einem Temperaturbereich von 770 °C bis 920 °C ermittelt wurde, vor“。„wird die Reaktion im Simulationsprogramm nur in der feststoffbeladenen Suspensionsphase berücksichtigt.“

**Q4** 采用范围是悬浮相。动力学动词是 „schlagen … vor“，程序句只限制在哪一相计算，没有再写一遍常数来源。

**Q5** \(R_g\) 未给数值。137327 的能量单位未印在数字旁。(5.51) 没有重复 „SI-Einheiten“ 这句。对水的级数未出现。770 °C 与 920 °C 的边界是开区间还是闭区间，原文写 „von … bis …“，未用不等式。

**Q6** 比对阶段填写。

---

## 15. R10 与 R11 的计量，式 (5.52)–(5.58)

**Q1** PDF 57 右半页 S.95 至 PDF 58 左半页 S.96，§5.2.6。

**Q2** 代表物：莱茵褐煤用苯与萘；生物质用萘与十六烷 \(\mathrm{C}_{16}\mathrm{H}_{34}\)，不用苯。原文给出木材摩尔 H/C 约 1,364，褐煤 0,856。苯加萘配不出 1,364 这一句是原文的理由。

R10：

\[
\begin{align*}
\mathrm{C}_6\mathrm{H}_6+3\,\mathrm{O}_2 &\rightarrow 6\,\mathrm{CO}+3\,\mathrm{H}_2 \tag{5.52}\\
\mathrm{C}_{10}\mathrm{H}_8+5\,\mathrm{O}_2 &\rightarrow 10\,\mathrm{CO}+4\,\mathrm{H}_2 \tag{5.53}\\
\mathrm{C}_{16}\mathrm{H}_{34}+8\,\mathrm{O}_2 &\rightarrow 16\,\mathrm{CO}+17\,\mathrm{H}_2 \tag{5.54}
\end{align*}
\]

R11：

\[
\begin{align*}
\mathrm{C}_6\mathrm{H}_6+\frac{9}{2}\mathrm{H}_2\mathrm{O} &\rightarrow \frac{9}{2}\mathrm{CO}+\frac{9}{2}\mathrm{H}_2+\frac{3}{2}\mathrm{CH}_4 \tag{5.55}\\
\mathrm{C}_{10}\mathrm{H}_8+\frac{23}{3}\mathrm{H}_2\mathrm{O} &\rightarrow \frac{23}{3}\mathrm{CO}+7\,\mathrm{H}_2+\frac{7}{3}\mathrm{CH}_4 \tag{5.56}\\
\mathrm{C}_{16}\mathrm{H}_{34}+\frac{21}{2}\mathrm{H}_2\mathrm{O} &\rightarrow \frac{21}{2}\mathrm{CO}+\frac{33}{2}\mathrm{H}_2+\frac{11}{2}\mathrm{CH}_4 \tag{5.57}
\end{align*}
\]

全局氧化：

\[
\mathrm{C}_m\mathrm{H}_n+\frac{m}{2}\mathrm{O}_2\rightarrow m\cdot\mathrm{CO}+\frac{n}{2}\cdot\mathrm{H}_2
\tag{5.58}
\]

与 (5.52)–(5.54) 的系数一致。两种代表物的摩尔配比没有公式，只给目标 H/C。

**Q3** 生物质不用苯、因为 H/C 约 1,364 而褐煤为 0,856，见 §5.2.6 起段。R10、R11 的引导句分别写 „mit Sauerstoff (Reaktion R10)“ 与 „in Gegenwart von Wasserdampf (Reformierung R11)“。

**Q4** (5.58) 是 (5.52)–(5.54) 的通式。(5.59) 使用这一通式的动力学。

**Q5** 配比方程未给。1,364 带有 „etwa“。正文提到焦油还可与二氧化碳反应，但没有 \(\mathrm{CO}_2\) 的计量式或速率。S.71 的 R10、R11 焓只有符号，这里也没有数值焓。

**Q6** 比对阶段填写。

---

## 16. R10 氧化动力学，式 (5.59)–(5.61)

**Q1** PDF 58 左半页，S.96。Hamel 引用的 Siminski et al. (1972)，zitiert in Smoot und Smith (1985)。引用行的印刷见 Q5。

**Q2**

\[
R_{\mathrm{C}_m\mathrm{H}_n}=\frac{\mathrm{d}C_{\mathrm{C}_m\mathrm{H}_n}}{\mathrm{d}t}=-k_{10}\cdot\exp\left(-\frac{E}{R_g\cdot T}\right)\cdot T\cdot P^{0{,}3}\cdot C_{\mathrm{C}_m\mathrm{H}_n}^{0{,}5}\cdot C_{\mathrm{O}_2}
\tag{5.59}
\]

有前导负号。式外还有一个 \(T\)。\(P\) 的指数是 \(0{,}3\)。焦油级数 \(0{,}5\)，氧级数 1。

烯烃，式 (5.60)：

\[
k_{10}=59{,}8\ \frac{\mathrm{m}^{1{,}5}}{\mathrm{mol}^{0{,}5}\cdot\mathrm{K}\cdot\mathrm{Pa}^{0{,}3}\cdot\mathrm{s}},\qquad E/R_g=12.200\,\mathrm{K}=12200\,\mathrm{K}
\]

芳烃，式 (5.61)：

\[
k_{10}=20.700=20700,\qquad E/R_g=9.650\,\mathrm{K}=9650\,\mathrm{K}
\]

单位与 (5.60) 相同。\(P\) 的单位已由 Pa 写进 \(k_{10}\) 的单位。若 \(C\) 为 \(\mathrm{mol/m}^{3}\)、\(P\) 为 Pa、\(T\) 为 K，则右端量纲为 \(\mathrm{mol/(m}^{3}\cdot\mathrm{s)}\)。这是由印刷单位得到的，不是另补的声明。

**Q3** „Die Kinetik dieser Reaktion berechnet sich nach Siminski et al. (1972) zitiert in Smoot und Smoot und Smith (1985) zu“。图像上 „Smoot und“ 出现两次，见 Q5。

**Q4** 与 (5.58) 及 (5.52)–(5.54) 相连。苯、萘按芳烃一行，十六烷原文没有点名归入烯烃还是芳烃。

**Q5**

- 引用行待人工辨认：放大图读作 „Smoot und Smoot und Smith (1985)“。可能是换行重复，也可能就是这样印的。不擅自删成一次。
- 十六烷用 (5.60) 还是 (5.61)，原文未写。
- 这一套氧化动力学没有按气泡相与悬浮相分开。
- 温度适用范围未给。
- 与第 17 节的一级热解参数是两套过程，原文没有把它们加成一个速率。

**Q6** 比对阶段填写。

---

## 17. 焦油二次反应与表 5.4

**Q1** PDF 58 左半页 S.96 的正文，右半页 S.97 的 Tabelle 5.4 与采用段。

**Q2** 表题：各参数均针对一级单反应。列：Stoff、Reaktor、Temp.-intervall °C、Frequenzfaktor \(1/\mathrm{s}\) *)、Aktivierungsenergie kJ/mol *)、Bettmaterial、Gas、Quelle。

脚注：WS = Wirbelschicht，FB = Festbett，LR = leerer Reaktor。*) 一级，按反应体积，\(R_{Teer}\) 单位 \(\mathrm{mol\cdot(s\cdot m}^{3}\mathrm{)}^{-1}\)。**) 一级，按催化剂质量，\(R_{Teer}\) 单位 \(\mathrm{mol\cdot(s\cdot kg_{Kat})}^{-1}\)；频率因子单位 \(\mathrm{m}^{3}\cdot(\mathrm{s\cdot kg_{Kat}})^{-1}\)。列头仍印 \(1/\mathrm{s}\)，**) 行不是 \(1/\mathrm{s}\)。

| 物质 | 反应器 | 温度 °C | 频率因子 | \(E\) kJ/mol | 床料 | 气体 | 来源 |
|---|---|---|---|---|---|---|---|
| 匹兹堡 #8 热解焦油 | FB | 550–850 | \(8{,}18\cdot 10^{3}\) *) | 86,3 | 石英 | He | Wen und Cain 1984 |
| 生物质热解焦油 | LR | 600–800 | 短横线 | 90,0 | 短横线 | Rohgas/He | Donnot et al. 1985 |
| 同上 | FB | 600–800 | 短横线 | 77,0 | Sand | Rohgas/He | 同上 |
| 同上 | FB | 600–800 | 短横线 | 46,1 | Dolomit | Rohgas/He | 同上 |
| 匹兹堡 #8 热解焦油 | LR | 500–900 | \(5{,}42\cdot 10^{4}\) 与 \(9{,}96\cdot 10^{2}\) *) | 100,5 与 64,5 | 短横线 | He | Serio et al. 1987 |
| 生物质气化焦油 | WS | 600–740 | \(7{,}5\cdot 10^{-4}\) **) | 5,9 | Sand | Rohgas | Corella et al. 1991 |
| 同上 | WS | 740–900 | 0,7 **) | 63,1 | Sand | Rohgas | 同上 |
| 气化焦油 | FB | 550–950 | \(1{,}5\cdot 10^{10}\) **) | 252,0 | Dolomit | Rohgas | Simell 1996 |

`8,18` 是逗号小数，即 \(8{,}18\times 10^{3}\,\mathrm{s}^{-1}\)（*) 行）。Serio 的两对按单元格内上下行对应：\(5{,}42\times 10^{4}\,\mathrm{s}^{-1}\) 对 \(100{,}5\,\mathrm{kJ/mol}\)，\(9{,}96\times 10^{2}\,\mathrm{s}^{-1}\) 对 \(64{,}5\,\mathrm{kJ/mol}\)。Corella 的 \(7{,}5\cdot 10^{-4}\) 指数为负，单位是脚注 ** 的 \(\mathrm{m}^{3}/(\mathrm{s\cdot kg_{Kat}})\)。\(E\) 的 SI 为表值乘 \(10^{3}\,\mathrm{J/mol}\)。

原文没有写出 \(k=k_0\exp(-E/(RT))\)。\(E\) 为 kJ/mol 时 \(R\) 取何值，未写。

**Q3** „wird im Rahmen dieser Arbeit zur Abschätzung der Sekundärreaktionen und der Sekundärpyrolyse des Teers auf die von Corella et al. (1991) und Serio et al. (1987) ermittelten kinetischen Parameter zurückgegriffen. Die Ergebnisse von Serio et al. (1987) zur homogenen Teerpyrolyse werden dabei für die feststofffreie Blasenphase und die in der Wirbelschicht ermittelten Parameter von Corella et al. (1991) in der Suspensionsphase zugrundegelegt.“

同页稍前：„Die Umsetzung des Teeres in Form von Sekundärpyrolyse oder Dampfreformierung wird zusammenfassend durch die stöchiometrischen Reaktionsgleichungen (5.55) und (5.56) beschrieben.“

**Q4** 气泡相用 Serio，悬浮相用 Corella。计量上，二次热解与蒸汽重整被这句指向 (5.55)、(5.56)。固体催化裂解与蒸汽重整。

**Q5**

- 总结句只点 (5.55) 与 (5.56)，没有点 (5.57)。十六烷的蒸汽重整计量是否同样进入模型，原文未写。
- Serio 单元格有两对参数。正文说一级单反应参数显得偏低，三组分模型更高，第三组分在实验条件下视为不转化。表内没有标明气泡相用哪一对，也没有标明两对与两个活性组分的对应。
- Corella 两个温度区间如何选择、740 °C 如何归属，原文未写。
- 速率对 \(\mathrm{H}_2\mathrm{O}\) 或 \(\mathrm{CO}_2\) 的级数没有出现。表是一级单反应。它与 R11 计量式消耗水如何耦合，原文未说明。
- Donnot 的频率因子是短横线。Wen und Cain 的正文另说用活化能的高斯分布；表中仍有一对一级参数。这两处都不是采用句点名的来源。Simell (1996) 在表内，采用句没有点名。
- Arrhenius 形式与 \(R\) 的取值未写。

**Q6** 比对阶段填写。

---

## 待人工辨认

不因此放弃已辨认的其余各项。

1. 式 (5.18)，PDF 53 右，S.87：分母是 \(v_g\) 还是 \(\nu_g\)。本稿按可见字形记拉丁字母 \(v\)。
2. 式 (5.14)、(5.15)，同一页：\(k\) 的下标是 1 还是 I。本稿记 \(k_1\)。
3. 式 (5.47)，PDF 57 左，S.94：`36.893` 与 `4.019` 的点号。按全书通例是千分位；改读成小数点才在数量级上不像是把平衡常数写成了另一个量。不选取。
4. 式 (5.59) 的引用行，S.96：„Smoot und“ 是否印了两次。
