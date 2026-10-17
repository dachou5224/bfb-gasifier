全新会话，转录冻结前未读旧材料

# T16 答卷：Hamel §4.1–4.3 干燥与热解

- 任务：T16
- 执行者：cursor
- 日期：2026-10-17
- 会话：全新会话。冻结前只读 `docs/CHARTER.md`、`docs/tasks/T16-transcribe-pyrolysis.md`、`.claude/skills/hamel-source-check/SKILL.md`。
- 来源：Hamel (2001) 页面图像。OCR 只用于定位，式号、数字、下标以页面图像为准。
- 范围：印刷页 S. 50–71。任务卡点名的式 (4.1)–(4.30)、表 4.1、表 4.2，以及 S. 71 反应焓表。同页印出、否则产率链或 Gibbs 调用会被截断的材料一并转录：表 4.3、表 4.4、无号辅助式、式 (4.31)、图中印出的工况。`index.json` 的方程键在定位时列到 (4.30)；(4.31) 是 S. 70 印面，不在该列表里。
- Q6：冻结时一律不填，比对阶段写入 `ch4_pyrolysis_diff.md`。
- 本文件冻结后不回改。后发现的转录差错只记入比对文件的勘误节。

德文排印：句点是千分位，逗号是小数点。下面保留印面记号，Q2 给出 SI。waf = wasser- und aschefrei；wmf = wasser- und mineralfrei。正文用到这两个缩写时按印面保留。

页面对照（每张 PDF 先左页、后右页）：

| 印刷页 | PDF |
|---|---|
| S. 50–51 | 35 |
| S. 52–53 | 36 |
| S. 54–55 | 37 |
| S. 56–57 | 38 |
| S. 58–59 | 39 |
| S. 60–61 | 40 |
| S. 62–63 | 41 |
| S. 64–65 | 42 |
| S. 66–67 | 43 |
| S. 68–69 | 44 |
| S. 70–71 | 45 |

## 符号表中与产率有关的行

Q1. PDF 8，印刷页 S. XIII 右页符号表。放大图只切到定义句，符号名的字头被裁掉；下列名称按定义句与 OCR 定位对齐，定义句本身来自放大图。

印面定义句（单位均为 kg/kg，除非另注）：

- \(X_{C,\mathrm{Gas}}\)：Masse C im Gas aus der Pyrolyse bezogen auf C im Brennstoff waf
- \(X_{C,\mathrm{Koks}}\)：Masse C im Koks aus der Pyrolyse bezogen auf C im Brennstoff waf
- \(X_{C,\mathrm{Teer}}\)：Masse C im Teer aus der Pyrolyse bezogen auf C im Brennstoff
- \(X_{\mathrm{Fl}}\)：Flüchtigengehalt bezogen auf Brennstoff waf
- \(X_{\mathrm{Fl},0}\)：Flüchtigengehalt nach Standardversuch bezogen auf Brennstoff waf
- \(X_{H,\mathrm{Gas}}\)：Masse H im Gas aus der Pyrolyse bezogen auf Brennstoff waf
- \(X_{\mathrm{Koks},B}\)：Masse des Kokses aus der Pyrolyse bezogen auf Brennstoff waf
- \(X_{\mathrm{Teer},B}\)：Masse des Teeres aus der Pyrolyse bezogen auf Brennstoff waf

同页希腊字母块：\(\varphi_{ms}\) = dimensionslose Position der Trocknungsfront；\(\beta\) = stöchiometrischer Koeffizient。符号表另有 \(w_0\)（初始水 / 湿颗粒）、\(w_{0,\mathrm{tr}}\)（初始水 / 干颗粒）。

Q5. 见式 (4.26) 与式 (4.22)：\(X_{C,\mathrm{Gas}}\) 的定义句是「相对于燃料中的碳、waf」，式 (4.26) 的单位行是「气体中的碳 / 燃料（waf）」。\(X_{C,\mathrm{Teer}}\) 的定义句停在 „bezogen auf C im Brennstoff“，句末没有 waf；式 (4.22) 的单位行写有 (waf)。不把缺字补上。

Q6. 比对阶段填写。

## §4 与 §4.1 的模型陈述

Q1. PDF 35–36，S. 50–52，§4 开篇与 §4.1 假设。无式号。

印面要点，不改写：

- 湿燃料的加热、干燥与热解采用 Agarwal 等 (1986)。只干燥时用准稳态模型；干燥与脱挥发分同时进行时用非稳态模型，由已干外壳的温度分布得到释放速率。Fett 等 (1994)、Dersch (1994) 用于循环流化床中的湿褐煤，Wirsum (1998) 用于污泥团块。
- 动力学类比 Dersch (1994)：大量互相独立的一级平行反应（Hanbaba 等 1968；Anthony 与 Howard 1976）。各反应共用一个频率因子，活化能为高斯分布（Anthony 与 Howard 1976）。
- 产物分成三组：永久气体、油与焦油、残余焦炭。产率依赖燃料组成和热解条件（加热速率、温度、压力）。物种收成仍含全部元素的代表物。产率来自元素平衡、实验和半经验模型。Solomon 等 (1992) 只作为综述被引用。
- §4.1 假设（Agarwal 等 1986）：球形颗粒；无收缩、磨耗、破碎；水分与挥发分均匀；传热与导热为速率控制步骤，传质快；脱挥发分热中性；颗粒内热解产物彼此之间以及与水分之间无二次反应；球形干燥锋面由表及里；湿核处于蒸发温度 \(T_e\)；热解用 Anthony–Howard，活化能正态分布。两段：干燥与热解同时，热解只发生在已干壳层，湿核直径缩小；干燥结束后继续升温与热解。
- Bild 4.1（S. 52）：\(T\) 对 \(r\)。标注有颗粒表面、干燥锋面、湿核，以及 \(T_{WS}\)、\(T_0\)、\(T_e\)、\(r_e\)、\(r_0\)。

Q3. „Die Pyrolyse wird als thermisch neutral angesehen.“（脱挥发分视为热中性；假设清单中的一句。）

Q4. Agarwal 等 (1986)、(1984)；Dersch (1994)；Anthony 与 Howard (1976)。

Q5. 无编号假设，不是方程。热解产物如何分到单元格，本章没有式子，见文末「单元格接口」。

Q6. 比对阶段填写。

## 式 (4.1)

Q1. PDF 36，S. 52 左页，§4.1，式 (4.1)。

Q2. \(T\) 为 K，\(t\) 为 s，\(r\) 为 m，\(a\) 为 m²/s。印面不写单位。

印面，适用于 \(r_e \le r \le r_0\)：

\[
\frac{\partial T}{\partial t} = \frac{a}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial T}{\partial r}\right)
\]

Q3. 已干外壳上的傅里叶方程，物性取常数。

Q4. 球形导热；\(a\) 由表 4.1 或燃料物性给出，原文未在式边写出代入关系。

Q5. 无。

Q6. 比对阶段填写。

## 式 (4.2)

Q1. PDF 36，S. 52 左页，式 (4.2)。

Q2. \(\lambda_s\) 为 W/(m·K)，\(\alpha\) 为 W/(m²·K)，\(\dot q\) 为 W/m²，温度 K。

印面：

\[
\lambda_s\left.\frac{\mathrm{d}T}{\mathrm{d}r}\right|_{r=R_0} = \alpha\,(T_a - T_s) = \dot q(t)
\]

Q3. 表面热流等于对流换热。

Q4. 式 (4.9) 给 \(\alpha\)。

Q5. 极限写 \(R_0\)，图 4.1 与式 (4.1) 写 \(r_0\)。环境温度写 \(T_a\)，图 4.1 的床温写 \(T_{WS}\)。原文未说明 \(R_0\) 是否即 \(r_0\)，也未说明 \(T_a\) 是否即 \(T_{WS}\)。

Q6. 比对阶段填写。

## 式 (4.3)

Q1. PDF 36，S. 52 左页，式 (4.3)。

Q2. \(h'_v\) 为 J/kg，\(w_{0,\mathrm{tr}}\) 无量纲（kg 水 / kg 干颗粒），\(\rho_s\) 为 kg/m³，\(r_e\) 为 m，\(t\) 为 s。左侧为 W/m² 量级的导热。

印面：

\[
\lambda_s\left.\frac{\mathrm{d}T}{\mathrm{d}r}\right|_{r=r_e} = h'_v\, w_{0,\mathrm{tr}}\, \rho_s\, \frac{\mathrm{d}r_e}{\mathrm{d}t}
\]

Q3. 干燥锋面上的热流等于修正蒸发焓与锋面移动所对应的水释放。

Q4. \(h'_v\) 由式 (4.4) 定义。\(w_{0,\mathrm{tr}}\) 见符号表，是相对于干颗粒质量的初始水，不是式中的另一个因子。

Q5. 印面右侧是 \(\mathrm{d}r_e/\mathrm{d}t\)，没有负号。锋面向内时该导数为负，而已干壳层的 \(\partial T/\partial r\) 为正。

- 读法 A（印面）：不添加负号。
- 建议读法：右侧取 \(-\mathrm{d}r_e/\mathrm{d}t\)，使两侧在锋面内移时同号。理由只是热流方向，不是另一个印出的字符。

不自行改式。

Q6. 比对阶段填写。

## 式 (4.4)

Q1. PDF 36，S. 52 左页，式 (4.4)。

Q2. \(h'_v\)、\(h_v\) 为 J/kg；\(c_w\)、\(c_s\) 为 J/(kg·K)；温度为 K。

印面：

\[
h'_v = h_v + \left(c_w + \frac{c_s}{w_{0,\mathrm{tr}}}\right)(T_e - T_0)
\]

Q3. 修正焓 = 水的蒸发焓，加上把湿核中的固体和水从初始温度带到 \(T_e\) 的显热。

Q4. 正文把初始温度称作 \(T_{P,0}\)。

Q5. 式中印的是 \(T_0\)，正文写 \(T_{P,0}\)。原文未声明二者是同一个温度。\(h'_v\) 带撇，\(h_v\) 不带撇；不要把第二项的蒸发焓也加上撇。

Q6. 比对阶段填写。

## 式 (4.5)–(4.8)

Q1. PDF 36，S. 53 右页，式 (4.5)–(4.8)。准稳态降阶，Agarwal 等 (1984)。

Q2. \(\varphi_{ms}\) 无量纲。Bi 无量纲。\(\alpha\) 为 W/(m²·K)，\(R_0\) 为 m，\(\lambda_s\) 为 W/(m·K)。

印面：

\[
\frac{\mathrm{d}}{\mathrm{d}r}\left(r^2\frac{\mathrm{d}T}{\mathrm{d}r}\right) = 0
\tag{4.5}
\]

\[
T|_{r=R_0} = T_s,\qquad T|_{r=r_e} = T_e
\tag{4.6}
\]

\[
\varphi_{ms} = \frac{r_e}{r_0}
\tag{4.7}
\]

\[
\mathrm{Bi} = \frac{\alpha R_0}{\lambda_s}
\tag{4.8}
\]

式 (4.8) 行尾的逗号是句子标点（„…, kann Gleichung (4.5) integriert werden“），不是公式的一部分。

Q3. 准稳态下温度只由锋面位置决定，总干燥时间 \(\tau_s\) 可由边界条件 (4.3) 算出。

Q4. 式 (4.1) 的准稳态简化；符号表 \(\varphi_{ms}\)。

Q5.

- 积分得到的 \(T(\varphi_{ms})\) 以及 \(\tau_s\) 的闭合式没有印出。正文只说可以积分、可以由 \(\theta_s = t/\tau_s\) 得到锋面位置、表面温度和热流。
- 式 (4.7) 的希腊字母在符号表中是 \(\varphi_{ms}\)。式旁字形曾被读成大写 \(\Phi\)，与符号表是同一个锋面位置符号，本转录记 \(\varphi_{ms}\)。
- \(R_0\) 与 \(r_0\) 的关系同式 (4.2)。

Q6. 比对阶段填写。

## 式 (4.9)

Q1. PDF 36，S. 53 右页，式 (4.9)。Kunii 与 Levenspiel (1977)。

Q2. Nu、Re、Pr 无量纲。系数印面为 \(1{,}2\)，不是 1.2 的另一种写法以外的数。

印面：

\[
\mathrm{Nu}_p = 2 + 1{,}2 \cdot \mathrm{Re}^{1/2}\cdot \mathrm{Pr}^{1/3}
\]

Q3. 颗粒努塞尔数。

Q4. 用于式 (4.2) 的 \(\alpha\)。特征长度、特征速度原文未在本章定义。

Q5. 系数是 \(1{,}2\)。不改成 Ranz–Marshall 的 0.6。Re、Pr 的特征量原文未说明。

Q6. 比对阶段填写。

## 式 (4.10)

Q1. PDF 36，S. 53 右页，式 (4.10)。Anthony 与 Howard (1976)。

Q2. 左端无量纲。\(k_0\) 为 1/s，\(E\) 为 J/mol，\(R\) 为 J/(mol·K)，\(T\) 为 K，\(t\) 为 s。

印面（尚未释放的挥发分质量分率）：

\[
\frac{m_{V^*}-m_V}{m_{V^*}}
= \int_0^\infty \exp\left[-k_0\int_0^t \exp\left(-\frac{E}{RT}\right)\,\mathrm{d}t\right] f(E)\,\mathrm{d}E
\]

Q3. „Massenanteil der noch nicht freigesetzten Flüchtigen.“

Q4. 温度历程由干燥模型给出。\(f(E)\) 为式 (4.11)。

Q5. 星号印在 \(V\) 上，记作 \(m_{V^*}\)。最终可释放量是输入，本章没有把它写成某一张表的挥发分。

Q6. 比对阶段填写。

## 式 (4.11)

Q1. PDF 36，S. 53 右页，式 (4.11)。

Q2. \(f(E)\) 为 1/(J/mol)，使 \(\int f(E)\,\mathrm{d}E\) 无量纲。\(E_0\)、\(\sigma_E\) 为 J/mol。

印面：

\[
f(E) = \left[\sigma_E\,(2\pi)^{1/2}\right]^{-1}
\exp\left[\frac{(E-E_0)^2}{2\sigma_E^2}\right]
\]

Q3. 活化能的高斯正态分布。S. 50 已写 „Gaußsche Normalverteilung“。

Q4. 式 (4.10)、(4.12)。表 4.1、表 4.2 给 \(E_0\)、\(\sigma_E\)。

Q5. 指数里没有负号。这不是缺笔猜测：放大后分子 \((E-E_0)^2\) 前没有减号。前置因子则是正态密度的归一化因子。

- 读法 A（印面）：指数为 \(+(E-E_0)^2/(2\sigma_E^2)\)。
- 建议读法：指数为 \(-(E-E_0)^2/(2\sigma_E^2)\)。理由：正文称高斯正态分布，且前置因子已是正态密度的系数。另一种读法仍是印面，不把减号写入印面转录。

Q6. 比对阶段填写。

## 式 (4.12)

Q1. PDF 37，S. 54 左页，式 (4.12)。

Q2. 左端无量纲。\(R_0\)、\(r\) 为 m。

印面：

\[
\frac{m_{V^*}-m_V}{m_{V^*}}
= \frac{3}{R_0}
\int_0^{R_0}
\left\{
\int_0^\infty
\exp\left[-k_0\int_0^t\exp\left(-\frac{E}{RT}\right)\,\mathrm{d}t\right]
f(E)\,\mathrm{d}E
\right\}
r^2\,\mathrm{d}r
\]

分母 \(R_0\) 的像素：分数线上方是 „3“，线下方是带下标 0 的 \(R\)，没有上标 3。

Q3. 对整颗球做体积加权。温度 \(T(t,r)\) 由非稳态干燥模型给定（„vorgegeben“）。

Q4. 式 (4.10)、(4.11)；非稳态温度场。二次温度函数的系数没有印出。

Q5.

- 读法 A（印面）：前置因子 \(3/R_0\)。
- 建议读法：\(3/R_0^3\)。理由：左端无量纲，\(\int_0^{R_0} r^2\,\mathrm{d}r\) 的量纲是长度³，\(3/R_0^3\) 才使右端无量纲。印面没有这个上标，建议读法不替换印面。
- 湿核内是否热解、积分是否仍从 0 积到 \(R_0\)，原文没有讨论。
- 二次逼近的系数原文未印。

Q6. 比对阶段填写。

## 表 4.1 与图 4.2、图 4.3

Q1. PDF 37，S. 54 左页，表 4.1，§4.1.1 湿褐煤。参数取自 Dersch (1994)。

Q2. 印面 / SI：

| 符号 | 印面 | SI |
|---|---|---|
| \(E_0\) | 192.000 J/mol | \(1{,}92\times 10^5\) J/mol |
| \(\sigma_E\) | 40.000 J/mol | \(4{,}0\times 10^4\) J/mol |
| \(k_{0,P}\) | \(1{,}67\cdot 10^{13}\) 1/s | \(1{,}67\times 10^{13}\) s⁻¹ |
| \(\rho_P\) | 1.250 kg/m³ | 1250 kg/m³ |
| \(c_P\) | 1.256 J/(kg K) | 1256 J/(kg·K) |
| \(a\) | \(0{,}1\cdot 10^{-6}\) m²/s | \(1\times 10^{-7}\) m²/s |

图 4.2：褐煤，\(T_{WS}=538\,\mathrm{K}\)，\(u_0=2{,}0\,\mathrm{m/s}\)，\(\bar d_{P,B}=4{,}4\,\mathrm{mm}\)。图 4.3：\(T_{WS}=713\,\mathrm{K}\)，\(u_0=1{,}5\,\mathrm{m/s}\)，粒径同上。曲线上的读数不作为转录数字。

Q3. 床温低时脱挥发分不完全，挥发分含量趋近一个依赖温度的极限，该极限作为模拟输入。残余水约 3%，与床温无关，当作这些温度下不释放的化学结合水，模型输入作了相应修正。

Q4. 式 (4.10)–(4.12)。实验：Agarwal 等 (1986) 对比 Agarwal 等 (1984)，氮气、铝颗粒流化床，煤放在丝网篮中，干冰淬冷。

Q5. 式中的符号是 \(k_0\)、\(\rho_s\)、\(c_s\)，表中是 \(k_{0,P}\)、\(\rho_P\)、\(c_P\)。原文未写它们是同一组量。频率因子的指数 13 与导温系数的指数 −6 来自表的放大图，不是逐像素认位。

Q6. 比对阶段填写。

## §4.1.2、表 4.2、图 4.4–4.6

Q1. PDF 38–39，S. 56–58，§4.1.2 木材，表 4.2。动力学：Pyle 与 Zaror (1984)，转引 Liliedahl 与 Sjöström (1998)。换热仍用 Kunii 与 Levenspiel (1977)。

Q2. 表 4.2 印面 / SI：

| 符号 | 印面 | SI |
|---|---|---|
| \(E_0\) | 67.500 J/mol | 67500 J/mol |
| \(\sigma_E\) | 13.500 J/mol | 13500 J/mol |
| \(k_{0,P}\) | 2.500 1/s | 2500 s⁻¹ |
| \(\rho_P\) | 500 kg/m³ | 500 kg/m³ |
| \(c_P\) | 1.670 J/(kg K) | 1670 J/(kg·K) |
| \(\lambda_P\) | 0,1256 W/(m K) | 0.1256 W/(m·K) |

行名：\(\lambda_P\) 为 „Wärmeleitfähigkeit des Partikels“；\(\rho_P\)、\(c_P\) 写明是干燥颗粒。进入模型的燃料量是干燥燃料的导热系数、比热容和密度。

图中印出的工况（干燥家具废木，纯热解、无干燥）：

| 图 | 试验号 | \(T_{WS}\) | 速度 | 粒径 | \(m_B\) |
|---|---|---|---|---|---|
| 4.4 | 008 | 798,2 K | 0,81 m/s | 2 mm | 2,46 g |
| 4.5 | 028 | 978,8 K | 0,96 m/s | 2 mm | 2,37 g |
| 4.6 | 956 | 1.218,2 K | \(u_0=1{,}04\) m/s | \(d_{P,B}=2\) mm | 1,6 g |

1.218,2 = 1218.2 K。图 4.6 的图题指向表 6.1（第 6.1.2 节），工况印在图内。

Q3. 实验步骤与燃料见第 6 章。

Q4. 式 (4.9)–(4.12)，表 4.2。

Q5. 符号名 \(k_{0,P}\)、\(\rho_P\)、\(c_P\)、\(\lambda_P\) 与式中 \(k_0\)、\(\rho_s\)、\(c_s\)、\(\lambda_s\) 的对应，原文未写明。

Q6. 比对阶段填写。

## §4.2 选用的产率步骤

Q1. PDF 39–40，S. 58–60，§4.2 至 §4.2.1 开头。无新式号。

加热速率类别（Wen 与 Dutta 1979），印在 S. 58：

1. 慢速：Aufheizrate \(\ll 1\,\mathrm{K/s}\)
2. 中速：\(5\)–\(100\,\mathrm{K/s}\)
3. 快速：\(500\)–\(10^5\,\mathrm{K/s}\)
4. 闪速：\(> 10^6\,\mathrm{K/s}\)

常压、慢速到中速时，有机质中的残焦份额主要随热解温度升高而下降。升压使残焦和永久气体增加（Howard 等 1981），原因是挥发分扩散受阻、焦油二次分解。焦炭气化是最慢的一步。网络模型 FG-DVC、FLASHCHAIN、CPD 被提到，本文不用。经验关联被提到（Weimer 与 Clough 1981；Govind 与 Shah 1984）。Merrick (1983) 在流化床燃烧模型中被引用（Jensen 等 1995；Hannes 1996）；正文描述了 Merrick 的元素平衡想法，紧接着写出本文实际采用的步骤。

本文采用的步骤（S. 60）：

- 热解中的燃料当作封闭系统。
- 先定焦炭和焦油产率。
- 在模拟模型中，焦炭当作纯碳。
- 流化床条件下焦油的摩尔 H/C 近似等于原煤（Baumann 与 Möller 1991；Khan 1989；Tyler 1980）。
- 煤的焦油产率用 Ko 等 (1988) 的半经验模型；焦炭产率由工业分析（Kurzanalyse）给定。
- Ko 对 30 种以上的煤整体拟合较好，单种煤可以偏离很大。若已有详细的焦炭或焦油数据（文献或本文第 6 章），用那些数据，不用 Ko。
- 不在焦炭和焦油中的剩余元素构成永久气体。永久气体组成由剩余元素质量分数加上化学平衡估计。

Q3. „Im Simulationsmodell wird der Koks als reiner Kohlenstoff behandelt.“

Q4. 后面的 §4.2.1–§4.3。封闭系统是这一句陈述的，不是从单元格平衡反推的。

Q5. 「工业分析给定焦炭产率」没有写出 Kurzanalyse 的哪一项、哪一种基准（收到基、干基、waf）进入模型。这是挥发分基准问题在本章里的原文边界：本章没有出现表 7.2，也没有出现 Kurkela。不把本章以外的数字拿来裁定。

Q6. 比对阶段填写。

## 式 (4.13)、(4.14)

Q1. PDF 40，S. 60–61，§4.2.1.1，Ko 等 (1987, 1988)。

式 (4.13) 是图式，不是代数式：Kohle —[1] 桥键热解→ Metaplast —[2] 输运→ Teer。步骤 [1] 依赖煤种，步骤 [2] 带压力依赖。交联（Cross-linking）在此等同于再聚合，并等同于成焦。

式 (4.14) 印面是文字比：

\[
x'_{\mathrm{Teer}}
= \frac{\text{Anzahl labiler Brückenbindungen}\cdot\text{Gehalt an verfügbarem Wasserstoff}}{\text{Anzahl der Cross-links}}
\]

C、H、O、S 和 OH 的量都相对于水与矿物无基（wmf）燃料。Ko 的数据范围：加热 \(100\)–\(1000\,\mathrm{K/s}\)，粒径 \(50\)–\(100\,\mu\mathrm{m}\)，从褐煤到无烟煤，压力 \(100\,\mathrm{Pa}\) 到 \(9\,\mathrm{MPa}\)。正文说方程只简写，细节须查 Ko 等 (1988)。

Q2. \(x'_{\mathrm{Teer}}\) 在后文系数分支里与 15、31 比较，那些阈值按质量百分数使用；式 (4.14) 本身没有印单位。

Q3. 最大焦油产率由 C、H、O、S 和反应器压力给出。

Q4. 式 (4.15)–(4.18) 分别定义文字比中的三项。

Q5. 式 (4.14) 没有印成 \(N_{lB}\,x_{vH}/N_{CL}\)。不把文字比改写成符号式当作印面。

Q6. 比对阶段填写。

## 式 (4.15)–(4.18) 与两条无号式

Q1. PDF 40，S. 61。

印面：

\[
N_{lB}
= \left[\frac{(1-f_c)\cdot x_{C,\mathrm{wmf}}\cdot 100}{12}\right]^{1{,}8}
\tag{4.15}
\]

\[
N_{CL}
= \frac{x_{O,\mathrm{wmf}}\cdot 100}{16}
+ \frac{x_{\mathrm{Sorg},\mathrm{wmf}}\cdot 100}{32{,}066}
\qquad
(x_{O,\mathrm{wmf}}\cdot 100 > 3{,}5\ \mathrm{Gew.\%\ wmf})
\tag{4.16}
\]

\[
N_{CL}
= \frac{3{,}5}{16}
+ \frac{x_{\mathrm{Sorg},\mathrm{wmf}}\cdot 100}{32{,}066}
\qquad
(x_{O,\mathrm{wmf}}\cdot 100 \le 3{,}5)
\tag{4.17}
\]

\[
x_{vH} = x_{H,\mathrm{wmf}}\cdot 100 - \frac{x_{\mathrm{OH},\mathrm{wmf}}\cdot 100}{17}
\tag{4.18}
\]

无号：

\[
x_{\mathrm{OH},\mathrm{wmf}} = 33{,}2 - 0{,}32\cdot x_{C,\mathrm{wmf}}\cdot 100
\]

\[
f_c = 0{,}830526 - 2{,}008147\cdot x_{C,\mathrm{wmf}} + 2{,}241218\cdot x_{C,\mathrm{wmf}}^{2}
\]

定义：\(x_{C,\mathrm{wmf}}\) 为 kg C / kg 煤 (wmf)，氧、氢、OH、有机硫同此基准。\(f_c\) 为芳香度。

Q2. \(N_{lB}\)、\(N_{CL}\) 按「个数」使用，式边没有另给单位。\(1{,}8\) 是整个方括号的指数，不是前面的乘数。32,066 为硫的摩尔质量，g/mol，印面不写单位。3,5 为 Gew.-%（wmf）。\(f_c\) 与 \(x_{\mathrm{OH},\mathrm{wmf}}\) 的长数字来自放大图通读。

Q3. „Anzahl labiler Brückenbindungen“ / „Anzahl der Cross-links“ / „Gehalt an verfügbarem Wasserstoff“。

Q4. 式 (4.14)、(4.20)。原子量 12、16、17、32,066 直接印在式中。

Q5.

- 式 (4.15) 的下标字形在 l 与 I 之间不分。建议读法 \(N_{lB}\)（labil）。另一种读法 \(N_{IB}\)。
- 式 (4.16) 是严格大于 3,5，式 (4.17) 是 \(\le 3{,}5\)。边界点 3,5 只属于 (4.17)。
- \(f_c\) 的多项式里 \(x_{C,\mathrm{wmf}}\) 没有乘 100，而 \(x_{\mathrm{OH},\mathrm{wmf}}\) 的式中有 \(\cdot 100\)。照印，不统一。

Q6. 比对阶段填写。

## 式 (4.19)、(4.20) 与 \(L_P\)

Q1. PDF 41，S. 62，式 (4.19)、(4.20)。Parr (1932)，经 Ko 等 (1987, 1988)。

印面：

\[
x_{\mathrm{Mineral}} = 1{,}08\cdot x_{\mathrm{Asche}} + 0{,}55\cdot x_{\mathrm{Schwefel}}
\tag{4.19}
\]

\[
x_{\mathrm{Teer},\mathrm{wmf}} = a(L_P)\cdot x_{\mathrm{Asche}} + b(L_P)\cdot x'_{\mathrm{Teer}}
\tag{4.20}
\]

系数：

- \(x'_{\mathrm{Teer}} < 15\)：\(a=2\)，\(b=0\)
- \(15 \le x'_{\mathrm{Teer}} \le 31\)：

\[
a = \frac{1}{0{,}021533 + 0{,}02865\cdot L_P} - 36
\]

\[
b = 0{,}50803 + 0{,}696487\cdot L_P - 0{,}06959\cdot L_P^{2}
\]

- \(x'_{\mathrm{Teer}} > 31\)：

\[
a = 11{,}24071 + 9{,}743707\cdot L_P - 0{,}91326\cdot L_P^{2},\qquad b=0
\]

\[
L_P = 1 - \log_{10}(P\cdot 10^{-6})
\]

\(P\le 2{,}5\cdot 10^{6}\,\mathrm{Pa}\) 时 \(P=P\)；\(P> 2{,}5\cdot 10^{6}\,\mathrm{Pa}\) 时 \(P=2{,}5\cdot 10^{6}\,\mathrm{Pa}\)。

15 与 31 都属于中间支（下界 \(\le\) 的左端是 15，上界是 \(\le 31\)；第一支是严格小于 15，第三支是严格大于 31）。

Q2. \(P\) 为 Pa。\(L_P\) 无量纲。\(x_{\mathrm{Teer},\mathrm{wmf}}\) 的单位行未在式边重复；名称带 wmf。

Q3. „Der Schwefelgehalt der organischen Substanz \(x_{\mathrm{Sorg},\mathrm{wmf}}\) wird ebenfalls mit Gleichung (4.19) ermittelt.“

Q4. 式 (4.14) 的 \(x'_{\mathrm{Teer}}\)。灰与全硫来自元素分析时，Ko 等建议用 Parr。

Q5.

- 有机硫「也用式 (4.19) 求得」，但没有印出由 \(x_{\mathrm{Mineral}}\)、\(x_{\mathrm{Asche}}\)、\(x_{\mathrm{Schwefel}}\) 得到 \(x_{\mathrm{Sorg},\mathrm{wmf}}\) 的代数。不补 Parr 的移项。
- 式 (4.20) 把 \(a(L_P)\) 乘在 \(x_{\mathrm{Asche}}\) 上。原文未说明 \(x_{\mathrm{Asche}}\) 与 \(x_{\mathrm{Teer},\mathrm{wmf}}\) 是否同为 wmf，也未解释灰分项的含义。
- 硬煤的矿物与灰可以差很多（Brandt 1981），这是正文的提醒，不是换算式。
- 关联来自 50–100 μm。正文写明外推到流化床常用粒径是关键的，图 4.7 已有明显偏差。因此 Ko 只用于没有足够热解数据的煤的焦油估计。

Q6. 比对阶段填写。

## 图 4.7

Q1. PDF 41，S. 63，图 4.7。

纵轴与横轴文字：Gew.-% bez. auf Brennstoff\(_{\mathrm{(waf)}}\)，刻度约 0–40。图例包括：Ko 等 (1988) 之后的文献（星号）、流化床热解 Tyler (1980)、Khan (1989)、Güell 与 Kandiyoti (1993)。脚注：星号数据 „bezogen auf Kohle(wmf)“。图题：计算与测量的焦油比较，计算用 Ko 等 (1988)。

Q5. 坐标轴写 waf，星号脚注写 wmf。原文没有换算。不裁定哪一个是模型输入基准。

Q6. 比对阶段填写。

## 式 (4.21)、(4.22) 与表 4.3

Q1. PDF 41–42，S. 63–64，§4.2.1.2 莱茵褐煤。数据：Hegermann (1990)，加压流化床热解，正文写实验为 \(0{,}2\)–\(2{,}5\,\mathrm{MPa}\)、\(800\)–\(1.000^\circ\mathrm{C}\)。碳分配用经验数方程近似。式中压力为帕斯卡，温度为开尔文。

印面：

\[
x_{C,\mathrm{Koks}}
= \xi_1\cdot P_{WS}^{\xi_2}\cdot \xi_3^{(10^{-6}\cdot P_{WS})}\cdot T_{WS}^{\xi_4}
\tag{4.21}
\]

单位：kg C im Koks / kg C im Brennstoff (waf)。

\[
x_{C,\mathrm{Teer}}
= \xi_1\cdot P_{WS}^{\xi_2}\cdot T_{WS}^{\xi_3} + \xi_4
\tag{4.22}
\]

单位：kg C im Teer / kg C im Brennstoff (waf)。

两式的 \(\xi\) 不通用。(4.21) 里 \(\xi_3\) 是底数、\(\xi_4\) 是温度指数。(4.22) 里 \(\xi_3\) 是温度指数、\(\xi_4\) 是加项。

表 4.3（最小二乘，Hegermann 1990）：

|  | 焦炭中的碳，式 (4.21) | 焦油中的碳，式 (4.22) |
|---|---|---|
| \(\xi_1\) | 0,3027 | −0,2168 |
| \(\xi_2\) | 0,1047 | 0,04148 |
| \(\xi_3\) | 0,9356 | 0,1566 |
| \(\xi_4\) | −0,08578 | 1,333 |

适用范围：\(1.073\,\mathrm{K} \le T_{WS} \le 1.273\,\mathrm{K}\)，\(0{,}1\,\mathrm{MPa} \le P_{WS} \le 2{,}5\,\mathrm{MPa}\)。1.073 = 1073。

图 4.8 图题引用式 (4.21)，焦炭碳份额对压力，三条温度 1.073、1.173、1.273 K。图 4.9 图题引用式 (4.22)，焦油碳份额。曲线读数不转录。正文写：与第 6 章常压流化床自有数据的符合也令人满意。

Q2. \(P_{WS}\) 在式中为 Pa（正文声明）。表的适用范围行写 MPa。\(T_{WS}\) 为 K。\(x\) 为 kg C / kg C（waf）。

Q3. 碳在焦炭与焦油之间的分配由 Hegermann 的元素平衡用数方程近似。

Q4. 式 (4.23)、(4.24) 把碳份额换成相对燃料 waf 的质量。

Q5. 正文实验压力从 0,2 MPa 起，表 4.3 的适用范围从 0,1 MPa 起。两处都照录，不取一个。\(\le\) 来自表的放大图。

Q6. 比对阶段填写。

## 式 (4.23)、(4.24)

Q1. PDF 42，S. 64 左页。

印面：

\[
x_{\mathrm{Koks},B(\mathrm{waf})} = x_{C,B(\mathrm{waf})}\cdot x_{C,\mathrm{Koks}}
\tag{4.23}
\]

单位：kg Koks / kg Brennstoff (waf)。正文：焦炭总质量相对于无水无灰燃料。焦炭在模型中是纯碳，故焦炭质量等于焦炭中的碳质量。

\[
x_{\mathrm{Teer},B(\mathrm{waf})}
= x_{C,B(\mathrm{waf})}\cdot x_{C,\mathrm{Teer}}
\cdot\left(1+\frac{x_{H,B(\mathrm{waf})}}{x_{C,B(\mathrm{waf})}}\right)
\tag{4.24}
\]

单位：kg Teer / kg Brennstoff (waf)。

Q3. „das molare H/C-Verhältnis des Teeres gleich dem der Ausgangskohle ist“。

Q4. 式 (4.21)、(4.22)。燃料的 \(x_{C,B(\mathrm{waf})}\)、\(x_{H,B(\mathrm{waf})}\) 是输入，本章不给数值。

Q5. 式 (4.24) 的代数是焦油质量 = 焦油中的碳 + 按原煤质量 H/C 配上的氢，没有氧、氮、硫项。对褐煤，原文只写摩尔 H/C 等于原煤。「焦油只含 H 和 C」这句话印在 §4.2.2，对象是木材。不把木材那句推广成褐煤焦油不含氧的裁定。摩尔 H/C 相等与质量 H/C 相等在只含 C、H 时是一回事；原文对褐煤没有写「只含 C、H」。

Q6. 比对阶段填写。

## §4.2.1.3 泥炭、表 4.4、图 4.10、图 4.11

Q1. PDF 43，S. 66，§4.2.1.3。Arpiainen 与 Lappi (1989)，芬兰泥炭，流化床。仍用式 (4.21)、(4.22)。

Q3. „entfallen hier die druckabhängigen Terme.“ 实验只在常压下做。

表 4.4：

|  | 焦炭中的碳，式 (4.21) | 焦油中的碳，式 (4.22) |
|---|---|---|
| \(\xi_1\) | 5,499 | \(-6{,}331\cdot 10^{-5}\) |
| \(\xi_2\) | 0 | 0 |
| \(\xi_3\) | 1 | 1,373 |
| \(\xi_4\) | −1,43 | 1,068 |

适用范围：\(773\,\mathrm{K} \le T_{WS} \le 973\,\mathrm{K}\)，\(P_{WS}=0{,}1\,\mathrm{MPa}\)。

图 4.10 纵轴：\(x_{C,\mathrm{Koks}}\)，单位 kg C im Koks / kg C im Brennstoff (waf)。图题：按式 (4.21)、随床温、泥炭焦炭中的碳；Teilbild B 为焦炭碳的测量与计算比较。

图 4.11 纵轴：\(x_{C,\mathrm{Teer}}\)，单位 kg C im Teer / kg C im Brennstoff (waf)，图内文字是 Torf–Teer。图题印面：

> In Abhängigkeit der Wirbelschichttemperatur nach Gleichung (4.21) berechneter Kohlenstoffanteil im Teer … (Teilbild A). Vergleich zwischen Messung und Berechnung des C-Anteils im Koks (Teilbild B).

Q5. 交给目标所有者，不裁定。

- 原文说压力项省略，但没有印出省略之后的代数。
- 读法 A：式 (4.21) 的结构不变，表中 \(\xi_2=0\)、\(\xi_3=1\) 使压力因子在数值上成为 1，于是 \(x_{C,\mathrm{Koks}}=5{,}499\cdot T_{WS}^{-1{,}43}\)。在 773–973 K，这个值的量级远小于图 4.10 纵轴（约 0,2–0,6）。
- 读法 B：压力项从式中删去后，\(\xi\) 对应另一个没有印出的温度式。原文没有那个式子。
- 焦油列代入式 (4.22)（\(\xi_2=0\)）时，数量级与图 4.11 的 0,1–0,5 同量级。这只说明两列不能互换，不说明焦炭列应如何改写。
- 图 4.11 图题写式 (4.21) 和 „C-Anteil im Koks“，图本身是 Teer。照录，不改成 (4.22)。

Q6. 比对阶段填写。

## 式 (4.25) 与 §4.2.2

Q1. PDF 43–44，S. 67–68，§4.2.2 木材，式 (4.25)。

印面（放大到式子本身）：

\[
\frac{x_{\mathrm{Fl}}}{x_{\mathrm{Fl},0}}
= 1 - \exp\left\{\left[\xi_1\cdot(T_{WS}-573)\right]^3 - \xi_2\right\}
\tag{4.25}
\]

参数：\(\xi_1 = -0{,}002437\)，\(\xi_2 = 0{,}3718\)。

关联估计的是流化床挥发分产率相对于标准试验挥发分（DIN 51720）的比值。本章没有再定义 DIN 51720 的基准；符号表写 \(X_{\mathrm{Fl}}\)、\(X_{\mathrm{Fl},0}\) „bezogen auf Brennstoff waf“。图 4.12 纵轴文字：Relative Freisetzung Flüchtiger \(x_{\mathrm{Fl}}/x_{\mathrm{Fl},0}\)，in Gew.-% bez. auf Brennstoff (waf)。图例：Birkenrinde（Arpiainen 与 Lappi 1989）、Holzabfall（Horne 与 Williams 1996）、Birkenspäne（Corella 等 1988）、Holz Möbelproduktion（eigene Messung）。图 4.13：式 (4.25) 对温度，以及奇偶图，各种木材。

Q3. 两句照录：

- „Als Flüchtige sind hier die Permanentgase zu verstehen, die bei Umgebungsbedingungen gasförmig vorliegen.“
- „Koks nur aus C und der Teer aus H und C mit einem H/C-Verhältnis entsprechend der organischen Substanz des Brennstoffs.“

作者若把产物分成气、焦、焦油，关联里只用气体。焦油量可由这条关联加上焦炭量得到。元素分配之后，用 §4.3 计算气体组成。

Q2. 左端是两个质量份额之比，无量纲。573 为 K。\(\xi_1\) 的印面使 \(\xi_1\cdot(T-573)\) 无量纲（K⁻¹ 乘 K）。图轴另写 Gew.-%。

Q4. 标准试验挥发分 \(x_{\mathrm{Fl},0}\) 是输入。焦炭 = 纯 C、焦油 = H 与 C 且 H/C 等于燃料有机质，是这一节写明的假设。然后进入 §4.3。

Q5.

- 建议读法就是上面的印面：立方作用在整个乘积 \(\xi_1\cdot(T_{WS}-573)\) 上。更宽的截图曾把立方只看成 \((T_{WS}-573)\) 的指数。另一种读法：\(1-\exp[(\xi_1\cdot(T_{WS}-573)^3)-\xi_2]\)。以式子所在的放大图为准，建议读法取前一种。用印出的 \(\xi\) 时，前一种在约 600 K 附近给出明显低于 1 的平台，并随温度升向 1，与图 4.13 的形状相符；这是倾向的理由，不是改写印面的理由。括号没有做逐像素认位。
- 纵轴同时写比值和 Gew.-%。原文未说明百分号是比值的单位，还是两个 waf 质量分数各自的单位。
- DIN 51720 的操作定义（干基、waf、收到基）本章没有印。符号表与图轴写 waf。不把标准外部定义引进来。

Q6. 比对阶段填写。

## 式 (4.26)、(4.27)

Q1. PDF 44，S. 69 右页，§4.3 开头。

印面：

\[
x_{C,\mathrm{Gas}}
= x_{C,B(\mathrm{waf})} - x_{\mathrm{Koks},B(\mathrm{waf})} - x_{C,\mathrm{Teer}}\cdot x_{C,B(\mathrm{waf})}
\tag{4.26}
\]

单位行：kg C im Gas / kg Brennstoff (waf)。

\[
x_{H,\mathrm{Gas}}
= x_{H,B(\mathrm{waf})}
- x_{C,B(\mathrm{waf})}\cdot x_{C,\mathrm{Teer}}\cdot\left(\frac{x_{H,B(\mathrm{waf})}}{x_{C,B(\mathrm{waf})}}\right)
\tag{4.27}
\]

单位行：kg H im Gas / kg Brennstoff (waf)。不把第二项约成 \(x_{H,B}\cdot x_{C,\mathrm{Teer}}\)。

其余原在煤中的元素全部进入气相。正文没有为 O、N、S 再写一条与 (4.26) 同形的式子。

Q3. 剩余元素分配给气态热解产物。

Q4. 式 (4.23)、(4.24) 与 §4.2.2 的焦炭、焦油假设。干燥释放的水是否并入热解水，§4.3 没有写。

Q5.

- 符号表 \(X_{C,\mathrm{Gas}}\)：相对于燃料中的碳、waf。式 (4.26) 单位行：相对于燃料（waf）。两处都照录。不选一个作为「正确基准」。
- 符号表 \(X_{H,\mathrm{Gas}}\) 是相对于燃料 waf，与式 (4.27) 的单位行一致。
- 符号表 \(X_{C,\mathrm{Teer}}\) 的定义句没有 waf，式 (4.22) 的单位行有 (waf)。式 (4.26) 使用的 \(x_{C,\mathrm{Teer}}\) 与哪一条对齐，原文没有再声明。
- O、N、S 以及「全部进入气相」如何写成后面的 \(b_i\)，没有式子。

Q6. 比对阶段填写。

## 式 (4.28)–(4.31)：调用附录 A1

Q1. PDF 45，S. 70 左页，§4.3。Holub 与 Vonka (1976)。

气相要考虑的物质举例：CO、CO₂、CH₄、C₂H₄、NH₃、H₂S、H₂O、H₂。原文是 „z.B.“，不是封闭清单。平衡组成由自由焓最小化得到。平衡时 \(G\) 最小。恒压恒温下：

\[
G = \sum_{j=1}^{N_S} \mu_j n_j = \mathrm{Min!}
\tag{4.28}
\]

\(N_S\) = Anzahl der Substanzen。

\[
\mathrm{d}G = \sum_{j=1}^{N_S}\frac{\partial G}{\partial n_j}\,\mathrm{d}n_j = 0
\tag{4.29}
\]

\[
\sum_{j=1}^{N_S} a_{ij} n_j - b_i = 0,
\qquad i=1,\ldots,N_E
\tag{4.30}
\]

\(N_E\) = Anzahl der Elemente。

紧接的正文：用 \(\beta_{i,j}\) 表示元素 \(i\) 在物质 \(j\) 分子式中的化学计量系数，\(b_i\) 为元素 \(i\) 的物质的量，\(N_j\) 为物质 \(j\) 的物质的量。拉格朗日乘子把优化变成非线性方程组。定义

\[
F = \sum_{j=1}^{N_S}\mu_j n_j
+ \sum_{i=1}^{N_E}\lambda_i^{*}\left(\sum_{j=1}^{N_S} a_{ij} n_j - b_i\right)
\tag{4.31}
\]

\(F\) 的最小与式 (4.28) 在约束 (4.30) 下的最小相同。

调用句照录：„Eine detaillierte Herleitung dieser Methode zur Berechnung des chemischen Gleichgewichtes ist im Anhang A1 zu finden.“

本任务不重写附录 A1。§4.3 给出的接口就是：恒定 \(P\)、\(T\) 下最小化自由焓；元素平衡 (4.30) 为约束；拉格朗日函数为 (4.31)；推导在附录 A1；方法来源 Holub 与 Vonka (1976)。

Q2. \(G\)、\(F\) 为自由焓，式边未写 J。\(n_j\)、\(b_i\) 为物质的量。\(a_{ij}\) 无量纲（每分子的原子数）。\(\lambda_i^{*}\) 为乘子。

Q3. 见上面的附录 A1 调用句。

Q4. 式 (4.26)、(4.27) 以及「其余元素全部进气相」提供元素量。附录 A1 提供算法，本章不展开。

Q5.

- 式用 \(a_{ij}\)、\(n_j\)；正文用 \(\beta_{i,j}\)、\(N_j\)。符号表把 \(\beta\) 列为化学计量系数。不把 \(a_{ij}\) 改记成 \(\beta_{i,j}\)。
- „z.B.“ 的物种清单不是封闭集。原文未说明清单是否就这些。
- 这个平衡所用的温度和压力，§4.3 没有给数值，也没有写是床温、颗粒温度还是反应器压力。
- 干燥水分与热解生成的 H₂O 是否相加，原文未写。
- 式 (4.31) 不在定位用的 index 方程列表（该列表止于 4.30）。它印在 S. 70，是指向附录 A1 的那一段，因此收入本答卷。

Q6. 比对阶段填写。

## 单元格接口（S. 12，不在第 4 章）

Q1. PDF 16，S. 12 左页。第 4 章本身只给颗粒时间，不给床高。任务要求记录产物如何进入源项；第 4 章没有源项式。印出的唯一分配句在 S. 12，转录在此，供比对，不把它写成第 4 章的方程。

印面：在流动、分离和气固停留时间已经算出之后，在「加入的燃料参与所算出的固体混合」这一假设下，水和热解产物的释放按反应器高度确定。由这些结果，把每一格分到的热解产物份额，以及每一格的固体进出流，按流动模型赋给单元格。映射到单元格模型上的预计算结果提供平衡所需的物质流和能量流。单元格温度和物质流又通过外层迭代回到「Vorabrechung」，重新算预计算，直到收敛（图 2.2）。

Q5. 没有印出：释放量对高度的公式、颗粒时间如何变成高度、一格的份额如何计算、源项写在哪一条平衡的哪一项。不补。床层固体输运不在本任务内。

Q6. 比对阶段填写。

## S. 71 反应焓表

Q1. PDF 45，S. 71 右页。这是第 5 章开头，不是第 4 章的编号表。任务卡要求收入。表头：Reaktionsenthalpie (0,1013 MPa; 298,15 K)。正文：模型为计算焦炭和气态物质的转化率，计入下列反应。

Q2. 0,1013 MPa = 0.1013 MPa；298,15 K。焓为 kJ/mol，分母上的物种见各行。

| 编号 | 名称 | 印面 | \(\Delta H_R\) |
|---|---|---|---|
| R1 | Verbrennung | \(\mathrm{C}+\dfrac{1}{\Phi_c}\mathrm{O_2}\rightarrow\left(2-\dfrac{2}{\Phi_c}\right)\mathrm{CO}+\left(\dfrac{2}{\Phi_c}-1\right)\mathrm{CO_2}\) | \(-393{,}5 < \Delta H_R < -110{,}5\) kJ/mol C |
| R2 | heterogene Wassergasreaktion | \(\mathrm{C}+\mathrm{H_2O}\rightarrow\mathrm{CO}+\mathrm{H_2}\) | +131,4 kJ/mol C |
| R3 | Boudouard-Reaktion | \(\mathrm{C}+\mathrm{CO_2}\rightarrow 2\,\mathrm{CO}\) | +172,5 kJ/mol C |
| R4 | hydrierende Vergasung | \(\mathrm{C}+2\,\mathrm{H_2}\rightarrow\mathrm{CH_4}\) | −74,9 kJ/mol C |
| R5 | Kohlenmonoxid-Oxidation | \(\mathrm{CO}+\dfrac{1}{2}\mathrm{O_2}\rightarrow\mathrm{CO_2}\) | −283 kJ/mol CO |
| R6 | Wasserstoff-Oxidation | \(\mathrm{H_2}+\dfrac{1}{2}\mathrm{O_2}\rightarrow\mathrm{H_2O}\) | −241,9 kJ/mol H₂ |
| R7 | Methan-Oxidation | \(\mathrm{CH_4}+\dfrac{3}{2}\mathrm{O_2}\rightarrow\mathrm{CO}+2\,\mathrm{H_2O}\) | −519,4 kJ/mol CH₄ |
| R8 | homogene Wassergasreaktion | \(\mathrm{CO}+\mathrm{H_2O}\rightarrow\mathrm{H_2}+\mathrm{CO_2}\) | −41,1 kJ/mol CO |
| R9 | Methan-Reformierung | \(\mathrm{CH_4}+\mathrm{H_2O}\rightarrow\mathrm{CO}+3\,\mathrm{H_2}\) | +206,3 kJ/mol CH₄ |
| R10 | Teer-Oxidation | \(\mathrm{Teer}+\mathrm{O_2}\rightarrow\mathrm{CO}+\mathrm{H_2}\) | \(\Delta H_R<0\) |
| R11 | Teer-Umsetzung | \(\mathrm{Teer}+\mathrm{H_2O}\rightarrow\mathrm{CO}+\mathrm{H_2}+\mathrm{CH_4}\) | \(\Delta H_R>0\) |

R8 的箭头在放大图中是单向 \(\rightarrow\)。R7 是 \(\frac{3}{2}\mathrm{O_2}\)，产物是 CO 与 2 H₂O，不是完全燃烧。R1 的机理因子下标读作 \(\Phi_c\)。

Q3. 表头条件即 0,1013 MPa 与 298,15 K。

Q4. 这些反应不属于 §4.3 的热解平衡。§4.3 的物种举例与 R1–R11 不是同一张清单。不把 R1–R11 说成热解 Gibbs 最小化的物种集。

Q5. \(\Phi_c\) 的下标 c 来自放大图，未做像素认位。R10、R11 只有不等式，没有数值。Teer 的化学式没有印出。

Q6. 比对阶段填写。

## 覆盖

已转录：式 (4.1)–(4.30)，以及印在同页的 (4.31)；表 4.1、4.2、4.3、4.4；S. 71 的 R1–R11；§4.1 假设；§4.2 的产率步骤；§4.2.2 对挥发分的定义；S. 12 的单元格赋源句。

第 4 章没有印出、因而不能在本答卷里补上的内容：准稳态温度的闭合式；非稳态二次温度函数的系数；式 (4.12) 的上标（只作为建议读法）；式 (4.11) 指数的负号（只作为建议读法）；泥炭压力项省略后的代数；有机硫从式 (4.19) 到 \(x_{\mathrm{Sorg},\mathrm{wmf}}\) 的移项；O、N、S 进入 \(b_i\) 的式子；热解平衡的 \(T\) 与 \(P\)；干燥水与热解水的关系；颗粒时间到床层高度、再到各格源项的公式；表 7.2 与 Kurkela 的数字（不在 S. 50–71）。

## 建议读法一览（不替换印面）

| 位置 | 印面 | 建议读法 | 理由 | 另一种读法 |
|---|---|---|---|---|
| 式 (4.3) | \(\mathrm{d}r_e/\mathrm{d}t\) | 右侧取负号 | 锋面内移时导热与印面右侧反号 | 照印面 |
| 式 (4.11) | 指数无负号 | 指数加负号 | 正文写正态分布，前置因子已是正态密度 | 照印面 |
| 式 (4.12) | \(3/R_0\) | \(3/R_0^3\) | 左端无量纲，\(\int r^2\mathrm{d}r\) 为长度³ | 照印面 |
| 式 (4.15) 下标 | l / I 不分 | \(N_{lB}\) | 左侧文字是 labile Brückenbindungen | \(N_{IB}\) |
| 式 (4.25) | 立方在整个乘积上 | 与印面相同 | 式子放大图；与图 4.13 的温度形状相符 | 立方只在 \((T-573)\) 上 |

泥炭表 4.4 代入式 (4.21) 与图 4.10 的量级不合，物理上不能靠改一个符号决定，见该节 Q5，交目标所有者。图 4.11 的图题与图不一致，同样只记录。符号表与式 (4.26) 的碳基准不一致，只记录。
