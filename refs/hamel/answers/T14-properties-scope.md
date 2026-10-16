# T14 原文答卷：物性与气相化学平衡

| 项 | 内容 |
|---|---|
| 任务卡 | `docs/tasks/T14-properties-scope.md` |
| 执行者 | cursor |
| 会话 | 新开的 cursor 会话。未沿用做过 T4、T8、T9、T11、T12、T13 或 worktree 清理的会话 |
| 工作树 | `../bfb-T14`，分支 `task/T14-properties-scope` |
| 日期 | 2026-10-16 |
| 审查方 | codex |

只查证 Hamel 原文（页面图像；逐页文本只用于定位）。未打开 Hamel 引用的文献，未打开 `docs/hamel_submodels/` 和归档代码。不选物性数据来源，不引入外部数据库或关联式。第 6 问一律写“不适用”。

德文排版：点是千分位，逗号是小数点。下文在给出原文写法的同时写出 SI。压力 0,1013 MPa、符号表的 101300 Pa 与式 (5.20) 的 101.300 Pa 是同一个数：101300 Pa。

全书检索参考文献表（PDF 104–111）和正文后：没有 Reid / Prausnitz / Poling、VDI-Wärmeatlas、Landolt、JANAF、Barin、Knacke、Perry、NASA 多项式。PDF 106 上 Holub und Vonka 条目里的 “Reidel” 是出版社名，不是物性手册。PDF 22 的字母序列落在 “dreidimensionalen” 里，不是文献 Reid。

---

## 1. 模型用到的物性量

符号表为 PDF 6–9 / S. VIII–XV。下表每个符号一条六问记录。有数值或关联式的，公式转录在第 2 节，本节 Q5 指向该节。

### 符号 \(h\)（比焓）与 \(G\)（自由焓）

❓ **Q1 出处**：比焓：符号表 PDF 6 右 / S. IX；使用于式 (2.8)，PDF 19 右 / S. 19。自由焓 \(G\)：符号表 PDF 6 右 / S. IX，单位 J；使用于式 (4.28)、(A.1)、(A.5)、(A.8)。
❓ **Q2 变量与单位**：符号表把比焓写成 J/kg 与 J/mol 两种。式 (2.8) 中气体项是 \(\dot N\,h\)，固体项是 \(\dot m\,h\)。若 \(\dot N\) 为 mol/s，则气体比焓为 J/mol；若 \(\dot m\) 为 kg/s，则固体比焓为 J/kg。原文没有写明这两套单位如何换算。\(G\) 的单位符号表为 J。
❓ **Q3 德文原句**：式 (2.8) 下句 „Die Energieströme der aktiven Gase und Feststoffe beinhalten die Standardbildungsenthalpie.“
❓ **Q4 上游依赖**：式 (2.8) 把物流乘比焓后求和，并加上 \(\dot m_s h_s\)。同页写明势能和动能略去。壁面项是 \(\dot Q_{W,i}\) 与 \(\dot Q_{WÜ,i}\)，系数是装置输入，不是由物性关联式算出。
❓ **Q5 原文空白**：没有 \(h=h^\circ+\int c_p\,\mathrm{d}T\) 或任何生成焓数值。标准生成焓的参考温度、参考压力在本页未写。符号表在 S. IX 先把 \(h\) 定义为高度（m），后文又把 \(h\) 定义为比焓；式 (2.8) 里的 \(h\) 是比焓。下标 \(s\) 在本页没有定义，\(h_s\) 的比热和生成焓原文未给出。
❓ **Q6**：不适用。

### \(c_P\)、\(c_{P,s}\)、\(c_{P,w}\)（比热）

❓ **Q1 出处**：符号表 PDF 6 左 / S. VIII。固体与水的数值只出现在干燥模型，见第 2 节式 (4.4)、表 4.1、表 4.2。气体比热没有单独符号行。
❓ **Q2 变量与单位**：三者单位都是 J/(kg·K)。\(c_P\) 未标明是哪一相。
❓ **Q3 德文原句**：„spezifische Wärmekapazität“；„spezifische Wärmekapazität des Feststoffs“；„spezifische Wärmekapazität des Wassers“.
❓ **Q4 上游依赖**：式 (4.4) 的加热项用到水和固体比热。式 (4.9) 的 \(\mathrm{Pr}\) 若要闭合，还需要气体比热，但原文没有写出 \(\mathrm{Pr}\) 的定义。
❓ **Q5 原文空白**：气体比热原文未给出。符号表的水比热是 \(c_{P,w}\)，式 (4.4) 印刷为 \(c_w\)；同页 \(c_W\) 是气泡阻力系数。建议读法：式 (4.4) 的 \(c_w\) 是水的比热。另一种读法：它就是阻力系数 \(c_W\)，量纲放不进焓修正，不采用为建议。式 (2.8) 的固体比焓是否使用表 4.1 的 \(c_P\)，原文未说明（见该表记录）。
❓ **Q6**：不适用。

### \(h_v\)、\(h'_v\)、\(T_e\)（水的蒸发焓与蒸发温度）

❓ **Q1 出处**：符号表 PDF 6 右 / S. IX；定义式 (4.4)，PDF 36 左 / S. 52。\(T_e\)：符号表 PDF 8 左 / S. XII。
❓ **Q2 变量与单位**：\(h_v\)、\(h'_v\) 为 J/kg。\(T_e\) 为 K。式 (4.4) 见第 2 节。
❓ **Q3 德文原句**：„spezifische Verdampfungsenthalpie des Wassers“；„korrigierte Verdampfungsenthalpie“.
❓ **Q4 上游依赖**：式 (4.3) 的干燥前沿能量衡算使用 \(h'_v\)。
❓ **Q5 原文空白**：\(h_v\) 的数值原文未给出。\(T_e\) 的数值原文未给出。正文在 S. 148–149（PDF 84、85，重复扫描）只说燃料水的蒸发焓要由反应器温度支付，仍然没有数值。
❓ **Q6**：不适用。

### \(\mu_j\)、\(\mu_j^0\)（化学势与标准化学势）

❓ **Q1 出处**：\(\mu_j\)：符号表 PDF 9 左 / S. XIV，单位 J/mol。\(\mu_j^0\)：式 (5.43)，PDF 56 左 / S. 92；式 (A.7)、(A.12)，PDF 92 / S. 164–165。符号表没有单独的 \(\mu^0\) 行。
❓ **Q2 变量与单位**：J/mol。式 (A.12) 把 \(\mu_j^0/(RT)+\ln(P/P_0)\) 收成与组成无关的 \(c_j\)。
❓ **Q3 德文原句**：符号表 „chemisches Potential der Substanz j“；式 (A.7) 前 „Mit dem chemischen Potential für ideale Gase“.
❓ **Q4 上游依赖**：式 (5.43) 的 \(K_{eq}\) 与附录 A1 的目标函数都用它。热解气组成见第 4 节。
❓ **Q5 原文空白**：\(\mu_j^0(T)\) 的函数形式和数值原文未给出。式 (5.43) 没有写“与附录 A1 同一套数据”。两处符号同为 \(\mu^0\)。可选解释：(a) 同一套，因为没有第二套符号；(b) 原文未把数据来源连在一起，数值仍空缺。无论 (a) 或 (b)，数值都是原文未给出。式 (5.43) 里的气体常数印刷为 \(R\)，符号表的通用气体常数是 \(R_g\)。建议读作同一个常数。\(R_g\) 的数值原文未给出。
❓ **Q6**：不适用。

### \(P_0\)、反应焓标注的标准状态

❓ **Q1 出处**：\(P_0\)：符号表 PDF 7 右 / S. XI，„Bezugsdruck 101300 Pa“，单位列另写 Pa。反应焓表头：PDF 45 右 / S. 71。式 (5.20) 的参考压力：PDF 53 右 / S. 87。式 (A.7) 的 \(P/P_0\)：PDF 92 左 / S. 164。
❓ **Q2 变量与单位**：101300 Pa = 0,1013 MPa。反应焓的温度是 298,15 K。
❓ **Q3 德文原句**：„Bezugsdruck 101300 Pa“；表头 „Reaktionsenthalpie (0,1013 MPa; 298,15 K)“.
❓ **Q4 上游依赖**：式 (A.7)、(A.12) 的压力项；式 (5.20) 的压力修正；S. 71 的反应焓是生成焓数据的检验量，本身不是 \(\mu^0(T)\)。
❓ **Q5 原文空白**：没有一句说全书只认这一个标准状态。式 (2.8) 的生成焓参考态未写。298,15 K 只标在 S. 71 的反应焓上，没有写成 \(\mu^0\) 或比焓的积分下限。式 (A.7) 的分子是大写 \(P\) 还是小写 \(p\)：建议读作总压 \(P\)，与符号表的 \(P\) 一致；另一种读法是该组分的分压。
❓ **Q6**：不适用。

### \(R_g\)（通用气体常数）与 \(M\)（摩尔质量）

❓ **Q1 出处**：\(R_g\)：符号表 PDF 8 左 / S. XII。\(M\)：符号表 PDF 7 左 / S. X。式 (5.11) 用到碳的摩尔质量符号，PDF 53 右 / S. 87，没有数值。
❓ **Q2 变量与单位**：\(R_g\) 为 J/(mol·K)。\(M\) 为 kg/mol。
❓ **Q3 德文原句**：„allgemeine Gaskonstante“；„Molmasse“.
❓ **Q4 上游依赖**：平衡式 (5.43)、(A.7) 需要气体常数。气体密度若按理想气体计算，需要 \(M\) 与 \(R_g\)。
❓ **Q5 原文空白**：\(R_g\) 的数值原文未给出。各组分摩尔质量原文未给出。没有写出 \(\rho_g=PM/(R_g T)\)。
❓ **Q6**：不适用。

### \(\rho\)、\(\rho_g\)、\(\rho_{g,20}\)、\(\rho_P\)（密度）

❓ **Q1 出处**：符号表 PDF 9 左 / S. XIV。\(\rho_g\) 用于式 (3.20)、(5.19)。干颗粒密度的数值见表 4.1、表 4.2。惰性物料密度的一个算例见表 5.1（第 2 节）。
❓ **Q2 变量与单位**：kg/m³。\(\rho_{g,20}\) 标明 20 °C。
❓ **Q3 德文原句**：„Gasdichte“；„Gasdichte bei 20 °C“；„Dichte des Feststoffs“.
❓ **Q4 上游依赖**：式 (3.20) 把悬浮相动力粘度和运动粘度用 \(\rho_g(1-\varepsilon_d)\) 连起来。式 (5.19) 的 Schmidt 数用 \(\rho_g\)。S. 22 式 (3.3) 前和 S. 28 式 (3.29) 后只说用环境条件下的物性或气、固物性，没有给数值，不当作数据来源。
❓ **Q5 原文空白**：\(\rho_g\)、\(\rho_{g,20}\) 的数值和温度、压力关联式原文未给出。焦炭密度：式 (5.11) 写明假设焦炭密度恒定，但没有给数值。灰的密度原文未给出。表 5.1 的 2000 kg/m³ 只是该算例的惰性物料密度，不是全书的床料密度。
❓ **Q6**：不适用。

### \(\mu_g\)、\(\mu_{g,20}\)、\(\nu_g\)、\(\mu_d\)、\(\nu_d\)（粘度）

❓ **Q1 出处**：符号表 PDF 9 左 / S. XIV。气体运动粘度用于式 (5.18)，PDF 53 右 / S. 87。气体动力粘度用于式 (5.19)。悬浮相粘度的关联式是式 (3.20)–(3.22)，PDF 23 右 / S. 27，见第 2 节。
❓ **Q2 变量与单位**：动力粘度 kg/(m·s)，即 Pa·s。运动粘度 m²/s。\(\mu_{g,20}\) 标明 20 °C。
❓ **Q3 德文原句**：„dynamische Viskosität des Gases“；„kinematische Viskosität des Gases“.
❓ **Q4 上游依赖**：氧气扩散膜的 Reynolds 数与 Schmidt 数。气泡 Reynolds 数式 (3.19) 用悬浮相运动粘度。
❓ **Q5 原文空白**：气体粘度的数值和温度关联式原文未给出，包括 20 °C 的值。悬浮相粘度有关联式，但它不是气体粘度。
❓ **Q6**：不适用。

### \(D_g\)（扩散系数）

❓ **Q1 出处**：符号表 PDF 6 左 / S. VIII，„Diffusivität“，m²/s。氧气的式子是式 (5.20)，PDF 53 右 / S. 87。使用于式 (5.16)、(5.19)。
❓ **Q2 变量与单位**：m²/s。式 (5.20) 只写了氧气。
❓ **Q3 德文原句**：„Die Berechnung der Diffusivität des Sauerstoffs erfolgt nach Field et al. (1967)“.
❓ **Q4 上游依赖**：膜温度式 (5.21)，总压 \(P\)，以及式 (5.16) 的传质系数。
❓ **Q5 原文空白**：其他气体的扩散系数原文未给出。Field et al. (1967) 的书没有打开；Hamel 给出的就是式 (5.20) 这一支。引用范围只到氧气、压力和膜平均温度。
❓ **Q6**：不适用。

### \(\lambda\)、\(\lambda_s\)、\(\lambda_g\)、\(a\)（导热系数与导温系数）

❓ **Q1 出处**：符号表 PDF 8 右 / S. XIII。\(\lambda_s\) 用于式 (4.3)、(4.8) 的 Biot 数。木材颗粒的 \(\lambda_P\) 见表 4.2。褐煤表 4.1 给的是导温系数 \(a\)，不是 \(\lambda\)。
❓ **Q2 变量与单位**：导热系数 W/(m·K)。导温系数 m²/s。
❓ **Q3 德文原句**：„Wärmeleitfähigkeit des Feststoffs“；„Wärmeleitfähigkeit des Gases“；„Temperaturleitfähigkeit“.
❓ **Q4 上游依赖**：干燥壳层导热。式 (4.9) 若要由 Nu 得到 \(\alpha\)，需要气体导热系数，原文没有写出这个换算。
❓ **Q5 原文空白**：气体导热系数的数值原文未给出。褐煤的导热系数原文未给出（只给了 \(a\)、\(\rho_P\)、\(c_P\)，三者能否互相换算，原文没有写）。壁面传热不另要这些量：壁面系数是输入。
❓ **Q6**：不适用。

### \(\alpha\)（传热系数）与 \(\mathrm{Pr}\)、式 (4.9) 的 \(\mathrm{Re}\)

❓ **Q1 出处**：\(\alpha\)：符号表 PDF 8 右 / S. XIII，W/(m²·K)。颗粒传热关联式 (4.9)：PDF 36 右 / S. 53。颗粒表面能量衡算式 (5.22) 使用 \(\alpha_{P,g}\)，PDF 53 右 / S. 87。
❓ **Q2 变量与单位**：W/(m²·K)。Nu、Re、Pr 无量纲。
❓ **Q3 德文原句**：„Der Wärmeübergang von der Wirbelschicht auf den zu betrachtenden Partikel wird entsprechend dem Vorschlag von Kunii und Levenspiel (1977) berechnet“.
❓ **Q4 上游依赖**：式 (4.2) 的表面热流。式 (5.22) 的对流项。
❓ **Q5 原文空白**：式 (4.9) 没有定义 \(\mathrm{Re}\) 和 \(\mathrm{Pr}\)。气体的 \(\mu\)、\(c_P\)、\(\lambda\) 因此仍然原文未给出。不把外部定义写进模型。壁面的 Wärmedurchgangskoeffizient 在表 2.1 / 表 2.2 的输入清单里（PDF 16 / S. 12–13），原文未给出由物性计算壁面系数的关联式。
❓ **Q6**：不适用。

### \(\varepsilon_P\)、\(\sigma_S\)（发射率与斯特藩–玻尔兹曼常数）

❓ **Q1 出处**：符号表 PDF 8 右 / S. XIII（\(\varepsilon_P\)）、PDF 9 左 / S. XIV（\(\sigma_S\)）。数值范围与用法在式 (5.22) 下，PDF 53 右 / S. 87。
❓ **Q2 变量与单位**：\(\varepsilon_P\) 无量纲。\(\sigma_S\) 的单位是 W/(m²·K⁴)。
❓ **Q3 德文原句**：„Für die Partikelemissivität geben Ross et al. (1981) … einen Wert im Bereich \(0{,}79\le\varepsilon_P\le 0{,}91\) an.“
❓ **Q4 上游依赖**：式 (5.22) 的辐射项。颗粒温度方法的句子引用 LaNauze (1985)，范围只是“用表面能量平衡求颗粒温度”，没有另给物性数。
❓ **Q5 原文空白**：区间内取哪一个值，原文未说明。\(\sigma_S\) 的数值原文未给出。符号表写 \(\sigma_S\)，式 (5.22) 的辐射项印刷为 \(\sigma\)（未见下标 S）。建议读作同一个常数。另一种读法：式中的 \(\sigma\) 是另一个量。Ross et al. (1981) 没有打开；Hamel 给出的就是这个闭区间。
❓ **Q6**：不适用。

### 燃料、焦炭、焦油、惰性床料

❓ **Q1 出处**：燃料工业分析、元素分析与热值是输入清单，PDF 16 / S. 13，表 2.2 的 Betriebsparameter 条文。焦油模型物质：PDF 57 右 / S. 95。惰性床料比热：无专节。式 (2.8) 的 \(h_s\) 见上。
❓ **Q2 变量与单位**：热值在第 7 章图上作为产品气热值，单位 MJ/m³(V_N)，不是燃料生成焓。元素分析的质量分数符号见表 4 的 \(x_{C,B}\) 等，单位 kg/kg。
❓ **Q3 德文原句**：S. 95 „Für den Fall der Vergasung Rheinischer Braunkohle wurden die Substanzen Benzol und Naphtalin ausgewählt“；生物质 „anstatt Benzol Hexadekan \(\mathrm{C}_{16}\mathrm{H}_{34}\)“.
❓ **Q4 上游依赖**：热解元素衡算式 (4.26)、(4.27) 先扣掉焦炭和焦油，再对气相做自由焓最小化。R10、R11 的反应焓在 S. 71 只给符号 \(\Delta H_R<0\) 与 \(\Delta H_R>0\)。
❓ **Q5 原文空白**：由热值换算燃料生成焓的式子原文未给出。焦炭、焦油、苯、萘、十六烷的标准生成焓和比热原文未给出。惰性床料的比热原文未给出。Brandt (1981) 在 PDF 41 / S. 62 只被引来说明烟煤的矿物质与灰分可以差很多，不是焓数据；该书没有打开。
❓ **Q6**：不适用。

### 理想气体与压力范围

❓ **Q1 出处**：附录开篇，PDF 91 右 / S. 163。式 (5.43) 的 \(K_{eq}=K_p\)，PDF 56 左 / S. 92，写明 „für ideales Gas“. 装置压力举例：PDF 76 / S. 133 一带写反应器压力直到 1,0 MPa；表 7.2 所在的 PDF 77 列出 0,4 MPa 至 1,0 MPa 的加压工况。
❓ **Q2 变量与单位**：压力单位在这些表里是 MPa。
❓ **Q3 德文原句**：„alle betrachteten Komponenten wie ideale Gase verhalten und nur im gasförmigen Zustand vorliegen“.
❓ **Q4 上游依赖**：化学势式 (A.7) 与 \(K_p\) 的写法依赖这条假设。
❓ **Q5 原文空白**：没有一句说流动模型里的 \(\rho_g\) 也按理想气体计算，也没有讨论 0,4–1,0 MPa 下这条假设是否保持。状态方程原文未给出。
❓ **Q6**：不适用。

---

## 2. 原文给出的数值与关联式

每条单独作答。S. 22 与 S. 28 只提到物性，没有数值，不立为数据来源。

### 表 4.1 褐煤干颗粒物性（Dersch 1994）

❓ **Q1 出处**：PDF 37 左 / S. 54，表 4.1。节 4.1.1。
❓ **Q2 变量与单位**：原文与 SI：

| 行 | 原文 | SI |
|---|---|---|
| \(\rho_P\) | 1.250 kg/m³ | 1250 kg/m³ |
| \(c_P\) | 1.256 J/(kg·K) | 1256 J/(kg·K) |
| \(a\) | \(0{,}1\cdot 10^{-6}\) m²/s | \(1\cdot 10^{-7}\) m²/s |

同表还有动力学三行（\(E_0=192.000\) J/mol，\(\sigma_E=40.000\) J/mol，\(k_{0,P}=1{,}67\cdot 10^{13}\) 1/s）。它们不是物性模块的输入，本卡只记录它们与物性行出自同一张表、同一句采用语。
❓ **Q3 德文原句**：„Die verwendeten Parameter sind für Braunkohle aus der Arbeit von Dersch (1994) übernommen (Tabelle 4.1).“ 表题限定为 „Trocknung und Pyrolyse von feuchter Braunkohle in der Wirbelschicht“.
❓ **Q4 上游依赖**：干燥与热解颗粒模型（导热、热容、密度）。Dersch (1994) 没有打开。文献表：PDF 105，Dersch, J.，1994，Siegen 学位论文。
❓ **Q5 原文空白**：采用范围按表题和上一句，是湿褐煤在流化床中的干燥与热解颗粒模型。原文没有说这些 \(c_P\)、\(\rho_P\) 也用于式 (2.8) 的固体焓。可选解释：(a) 只用于第 4.1.1 节的颗粒模型；(b) 也用于 cell 能量衡算。待总司令裁定。褐煤的导热系数不在表内。
❓ **Q6**：不适用。

### 表 4.2 木材干颗粒物性

❓ **Q1 出处**：PDF 38 左 / S. 56，表 4.2。节 4.1.2。
❓ **Q2 变量与单位**：物性行的原文与 SI：

| 行 | 原文 | SI |
|---|---|---|
| \(\rho_P\) | 500 kg/m³ | 500 kg/m³ |
| \(c_P\) | 1.670 J/(kg·K) | 1670 J/(kg·K) |
| \(\lambda_P\) | 0,1256 W/(m·K) | 0,1256 W/(m·K) |

同表动力学行：\(E_0=67.500\) J/mol，\(\sigma_E=13.500\) J/mol，\(k_{0,P}=2.500\) 1/s（即 2500 s⁻¹）。
❓ **Q3 德文原句**：正文 „Als kinetische Parameter der Holzpyrolyse wurden die von Pyle und Zaror (1984), zitiert in Liliedahl und Sjöström (1998), veröffentlichten Werte verwendet (siehe Tabelle 4.2).“ 下一句把导热系数、比热和干燃料密度称为进入模型的 „brennstoffspezifische Größen“. 表题则写 „Parameter nach Pyle und Zaror (1984)“.
❓ **Q4 上游依赖**：木材热解颗粒模型。两篇文献都没有打开。
❓ **Q5 原文空白**：物性三行的来源有两种读法。(a) 表题把整张表都归于 Pyle und Zaror (1984)；(b) 正文只把动力学三行归于该文献，物性三行没有另外点名。建议读法：(a)，因为表题覆盖全部行；同时保留 (b)。适用范围按节 4.1.2，是木材热解颗粒模型。原文没有说它们用于式 (2.8)。
❓ **Q6**：不适用。

### 式 (4.4) 修正蒸发焓

❓ **Q1 出处**：PDF 36 左 / S. 52，式 (4.4)。
❓ **Q2 变量与单位**：
\[
h'_v=h_v+\left(c_w+\frac{c_s}{w_{0,tr}}\right)(T_e-T_0)
\]
各项为 J/kg，比热为 J/(kg·K)，\(w_{0,tr}\) 为 kg/kg（符号表：干颗粒基准的初始水质量）。\(T_0\) 与正文里的初始温度 \(T_{P,0}\) 是否同一个量，见 Q5。
❓ **Q3 德文原句**：„Die korrigierte Verdampfungsenthalpie \(h'_v\) setzt sich zusammen aus der Verdampfungsenthalpie des Wassers \(h_v\)“.
❓ **Q4 上游依赖**：式 (4.3)。\(c_s\) 与表 4.1 的 \(c_P\) 是否同一数，原文未写。
❓ **Q5 原文空白**：\(h_v\)、\(T_e\) 无数值。正文说从 \(T_{P,0}\) 加热到 \(T_e\)，式中印刷为 \(T_0\)。建议读法：\(T_0\) 就是 \(T_{P,0}\)。另一种读法：\(T_0\) 是另一个参考温度。\(c_w\)、\(c_s\) 与符号表 \(c_{P,w}\)、\(c_{P,s}\) 的对应见第 1 节。
❓ **Q6**：不适用。

### 式 (4.8)、(4.9) 颗粒传热

❓ **Q1 出处**：式 (4.8) Biot 数、式 (4.9) Nu，PDF 36 / S. 52–53。
❓ **Q2 变量与单位**：式 (4.9)
\[
Nu_P=2+1{,}2\cdot Re^{1/2}\cdot Pr^{1/3}
\]
式 (4.8) \(Bi=\alpha R_0/\lambda_s\)（原文把颗粒半径写成 \(R_0\)）。
❓ **Q3 德文原句**：见第 1 节 Kunii und Levenspiel (1977)。
❓ **Q4 上游依赖**：表面热流。Kunii und Levenspiel (1977) 没有打开。Hamel 写出的就是这一支 Nu 关联式。
❓ **Q5 原文空白**：Re、Pr 的定义，以及气体 \(\mu\)、\(c_P\)、\(\lambda\)、\(\rho\)，原文未给出。
❓ **Q6**：不适用。

### 式 (3.20)–(3.22) 悬浮相粘度

❓ **Q1 出处**：PDF 23 右 / S. 27。Hilligardt (1986)。
❓ **Q2 变量与单位**：
\[
\mu_d=\nu_d\cdot\rho_g\cdot(1-\varepsilon_d)
\]
\[
d_P<0{,}2\,\mathrm{mm}:\quad \mu_d=f_d\,d_P^{2/3},\quad f_d=232\,\frac{\mathrm{kg}}{\mathrm{s}\cdot\mathrm{m}^{5/3}}
\]
\[
d_P\ge 0{,}2\,\mathrm{mm}:\quad \mu_d=f_d\,d_P^{3/2},\quad f_d=2{,}5\cdot 10^{5}\,\frac{\mathrm{kg}}{\mathrm{s}\cdot\mathrm{m}^{5/2}}
\]
单位指数是叠分式。放大后分子是 5，不是 1。与粒径指数 2/3、3/2 相乘后得到 Pa·s。若把单位指数读成 1/3 或 1/2，量纲不再是粘度。建议读法：5/3 与 5/2。另一种读法：1/3 与 1/2，列在这里，不采用为建议。
❓ **Q3 德文原句**：„Die dynamische Viskosität der Suspensionsphase wird in Abhängigkeit des Partikeldurchmessers nach Hilligardt (1986) wie folgt berechnet“.
❓ **Q4 上游依赖**：仍要气体密度 \(\rho_g\)。Hilligardt (1986) 没有打开。
❓ **Q5 原文空白**：代入 \(d_P\) 时用 m 还是 mm，原文只把阈值写成 0,2 mm。\(f_d\) 的单位用 m，建议 \(d_P\) 用 m。另一种读法：阈值用 mm，意味着粒径也用 mm，则 232 与 \(2{,}5\cdot 10^5\) 不能直接代入。这是原文空白，待总司令确认。此关联式不是气体粘度。
❓ **Q6**：不适用。

### 式 (5.16)–(5.22) 与氧气扩散系数

❓ **Q1 出处**：PDF 53 右 / S. 87。
❓ **Q2 变量与单位**：
\[
k_{d,g}=\frac{Sh\cdot D_g}{d_P},\quad
Sh=2+0{,}6\,Re^{0{,}5}\,Sc^{0{,}33},\quad
Re=\frac{u_r d_P}{\nu_g},\quad
Sc=\frac{\mu_g}{\rho_g D_g}
\]
\[
D_g=3{,}13\cdot 10^{-4}\,\frac{\mathrm{m}^2}{\mathrm{s}}\left(\frac{T_m}{1.500\,\mathrm{K}}\right)^{1{,}75}\left(\frac{101.300\,\mathrm{Pa}}{P}\right)
\]
\[
T_m=\frac{T_\infty+T_P}{2}
\]
\[
R_C\cdot\Delta h_{R_c}=A_P\cdot\alpha_{P,g}(T_P-T_\infty)+A_P\cdot\sigma\cdot\varepsilon_P(T_P^4-T_\infty^4)
\]
\(0{,}79\le\varepsilon_P\le 0{,}91\)。1.500 K 即 1500 K，101.300 Pa 即 101300 Pa。
❓ **Q3 德文原句**：氧气扩散见第 1 节；发射率见第 1 节；颗粒温度 „Nach einem Vorschlag von LaNauze (1985)“.
❓ **Q4 上游依赖**：\(D_g\) 只用于氧气膜扩散。Field et al. (1967)、Ross et al. (1981)、LaNauze (1985) 都没有打开。
❓ **Q5 原文空白**：\(\nu_g\)、\(\mu_g\)、\(\rho_g\)、\(\alpha_{P,g}\) 的数值原文未给出。\(\sigma\) 的数值原文未给出，与符号表 \(\sigma_S\) 的关系见第 1 节。\(\varepsilon_P\) 在区间内的取值原文未说明。其他气体的 \(D_g\) 原文未给出。
❓ **Q6**：不适用。

### S. 71 反应焓（0,1013 MPa；298,15 K）

❓ **Q1 出处**：PDF 45 右 / S. 71。R1 至 R11 列在第 5 章开头，不是附录。
❓ **Q2 变量与单位**：原文为 kJ/mol，SI 为 J/mol（乘 1000）。

| 反应 | 原文 | 基准 |
|---|---|---|
| R1 | \(-393{,}5<\Delta H_R<-110{,}5\) | kJ/mol C |
| R2 | \(+131{,}4\) | kJ/mol C |
| R3 | \(+172{,}5\) | kJ/mol C |
| R4 | \(-74{,}9\) | kJ/mol C |
| R5 | \(-283\) | kJ/mol CO |
| R6 | \(-241{,}9\) | kJ/mol H₂ |
| R7 | \(-519{,}4\) | kJ/mol CH₄ |
| R8 | \(-41{,}1\) | kJ/mol CO |
| R9 | \(+206{,}3\) | kJ/mol CH₄ |
| R10 | \(\Delta H_R<0\) | 未给数值 |
| R11 | \(\Delta H_R>0\) | 未给数值 |

R1 的不等号是印刷式：下界 −393,5，上界 −110,5，中间随 \(\Phi_c\) 变化。R8 的 −41,1 kJ/mol 与同页反应式 \(\mathrm{CO}+\mathrm{H_2O}\to\mathrm{H_2}+\mathrm{CO_2}\) 写在一起；Bild 5.11 的图内也标了 −41,1 kJ/mol（PDF 56 右 / S. 93）。
❓ **Q3 德文原句**：„Reaktionsenthalpie (0,1013 MPa; 298,15 K)“.
❓ **Q4 上游依赖**：这些数可以直接检验生成焓数据是否自洽。原文没有给出生成焓表，也没有写由这张表反算 \(\mu^0\) 的步骤。
❓ **Q5 原文空白**：R10、R11 没有数值。没有把反应焓声明为 \(\mu^0\) 或 \(K_{eq}\) 的输入。参考态只标在这张表上。
❓ **Q6**：不适用。

### 式 (5.43)

❓ **Q1 出处**：PDF 56 左 / S. 92，式 (5.43)。图在 PDF 56 右 / S. 93，Bild 5.11。
❓ **Q2 变量与单位**：
\[
K_{eq}=K_p=\exp\left(-\frac{\sum\nu_i\mu_i^0}{RT}\right)
\]
\(K_{eq}\) 无量纲。计量系数的希腊字母在此扫描上与 \(\nu\)、\(\upsilon\) 同形；建议读作化学计量系数 \(\nu_i\)。另一种读法：字母 \(\upsilon_i\)。求和上下限原文未写。
❓ **Q3 德文原句**：„die Gleichgewichtskonstante \(K_{eq}\) für ideales Gas“.
❓ **Q4 上游依赖**：ADR 0010 第 1 条要求 \(k_{eq,R8}\) 取这一式，由物性模块计算。该要求是裁定，不是本页原句。本页没有指向附录 A1。
❓ **Q5 原文空白**：\(\mu_i^0\) 的数值和温度函数原文未给出。标准压力没有写进这一式。与附录 A1 是否同一套数据，原文未明说。
❓ **Q6**：不适用。

### 表 5.1 的惰性物料密度

❓ **Q1 出处**：PDF 49 / S. 79 一带，表 5.1，算例边界。
❓ **Q2 变量与单位**：\(\rho_{inert}=2000\) kg/m³。同表还有 \(T\)、粒径等算例数，不是物性关联式。
❓ **Q3 德文原句**：表内 „Dichte des Inertmaterials“.
❓ **Q4 上游依赖**：该表是一组计算边界，不是符号表里的 \(\rho_P\) 定义。
❓ **Q5 原文空白**：没有说 2000 kg/m³ 适用于全部工况。惰性床料比热原文未给出。
❓ **Q6**：不适用。

### 文献表中与物性有关的条目（只记 Hamel 页码和引用范围）

| 作者、年份 | Hamel 文献表 | Hamel 正文里的引用范围 |
|---|---|---|
| Holub, R.; Vonka, P., 1976，*The Chemical Equilibrium of Gaseous Systems*，D. Reidel，Dordrecht | PDF 106 | 热解气自由焓最小化（S. 69）；阻尼参数与收敛判据（S. 167）。不是 \(\mu^0\) 数据表。书未打开 |
| Hoffmann 1971 | PDF 92 左 / S. 164 正文点名；文献表有 Hoffmann 条目 | 只用于 Lagrange 乘子。不是物性 |
| Dersch 1994 | PDF 105 | 表 4.1 的湿褐煤干燥、热解参数。书未打开 |
| Pyle and Zaror 1984；Liliedahl and Sjöström 1998 | 文献表在 PDF 109 一带 | 表 4.2，范围见该表记录。书未打开 |
| Field et al. 1967 | PDF 105 | 式 (5.20) 的氧气扩散系数；同书还被第 5、6 章用于燃烧速率，那不是本卡的物性范围。书未打开 |
| Ross et al. 1981 | PDF 109 一带 | 只给出 \(\varepsilon_P\) 的区间。书未打开 |
| LaNauze 1985 | PDF 105 | 只作为颗粒表面能量平衡的方法。书未打开 |
| Kunii und Levenspiel 1977 | 文献表有 Kunii 条目 | 式 (4.9) 这一支 Nu。书未打开 |
| Hilligardt 1986 | 文献表有 | 式 (3.20)–(3.22)。书未打开 |
| Brandt 1981，*Brennstoffe und Verbrennungsrechnung* | PDF 104 | PDF 41 / S. 62：烟煤矿物质与灰分之差。不是焓或热容。书未打开 |
| Agarwal et al. 1984、1986 | PDF 104 | 干燥模型结构与图 4.2、4.3 的测量，不是表 4.1 的物性数。书未打开 |

没有找到气体比热、生成焓、气体粘度、气体导热系数、气体密度的手册条目。

---

## 3. 附录 A1：式 (A.1)–(A.40)

位置：PDF 91 右至 PDF 94 左 / S. 163–168。A1 在 S. 168 结束；S. 169 起是附录 A2 的试验汇编，不是平衡算法。

假定（PDF 91 右 / S. 163）：所考虑的组分都按理想气体，并且只以气态出现。允许把一部分气体组分当作已知的源汇项，这些组分不处于该平衡中。\(\mu_j^0\) 的形式和数值原文未给出。收敛阈值 \(\varepsilon\) 的数值原文未给出。阻尼系数只给 \(0<\eta<1\)，没有选定值。

未编号、但解方程时要用的定义写在式 (A.26) 与式 (A.27) 之间（PDF 93 左 / S. 166），单独列在 (A.26) 之后。

印刷上反复出现“对已经是某个 \(j\) 的方程再对 \(j\) 求和”。下面各条照印刷转录，并在 Q5 给出建议读法。建议读法不是把印刷式改掉，而是交给总司令确认。

### 式 (A.1)

❓ **Q1 出处**：PDF 91 右 / S. 163，式 (A.1)。索引无坐标框，式在该半页中部。
❓ **Q2 变量与单位**：\(G=\sum_{j=1}^{N_S}\mu_j N_j\)。\(G\) 为 J，\(\mu_j\) 为 J/mol，\(N_j\) 为 mol。\(N_S\) 是气体组分数。
❓ **Q3 德文原句**：„Für ein Gemisch von \(N_S\) Substanzen ist die Gibbs'sche freie Enthalpie durch die Summe der chemischen Potentiale der Komponenten gegeben“.
❓ **Q4 上游依赖**：\(\mu_j\) 来自式 (A.7)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.2)

❓ **Q1 出处**：PDF 91 右 / S. 163，式 (A.2)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}N_j+\sum_{j=N_S-N_Q+1}^{N_S}a_{ij}N_j^+=b_i^0,\quad i=1,N_E
\]
\(a_{ij}\) 是元素 \(i\) 在物质 \(j\) 中的原子数。\(b_i^0\) 是元素 \(i\) 的总摩尔数。\(N_j^+\) 是已知源汇物质的摩尔数。\(N_Q\) 是源汇物质的个数。\(N_E\) 是元素个数。
❓ **Q3 德文原句**：„Massenbilanz unter Berücksichtigung der \(N_Q\) bekannten Quell- oder Senkensubstanzen“.
❓ **Q4 上游依赖**：元素总量和被固定的 \(N_j^+\) 是输入。
❓ **Q5 原文空白**：哪些物质算源汇，本页没有点名。与第 4 章焦炭、焦油的关系见第 4 节。
❓ **Q6**：不适用。

### 式 (A.3)

❓ **Q1 出处**：PDF 91 右 / S. 163，式 (A.3)。
❓ **Q2 变量与单位**：\(\sum_{j=1}^{N_S-N_Q}a_{ij}N_j=b_i^+\)，\(i=1,N_E\)。\(b_i^+\) 是扣掉源汇之后、留在平衡气体里的元素摩尔数。
❓ **Q3 德文原句**：„korrigierte Massenanteil \(b_i^+\) des Elementes \(i\) am verbleibenden Gasgemisch“. 这里的 „Massenanteil“ 与后文摩尔衡算混用；符号 \(b_i^+\) 的单位按上下文是 mol，不是 kg/kg。
❓ **Q4 上游依赖**：式 (A.4)。
❓ **Q5 原文空白**：用词 „Massenanteil“ 与摩尔量不一致。建议按后式把 \(b_i^+\) 当作元素摩尔数。另一种读法：它真是质量分数，则与 \(a_{ij}N_j\) 的单位不合。
❓ **Q6**：不适用。

### 式 (A.4)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.4)。
❓ **Q2 变量与单位**：\(b_i^+=b_i^0-\sum_{j=N_S-N_Q+1}^{N_S}a_{ij}N_j^+\)。
❓ **Q3 德文原句**：式前承接 „verbleibenden Gasgemisch“.
❓ **Q4 上游依赖**：式 (A.2)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.5)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.5)。
❓ **Q2 变量与单位**：
\[
G=\sum_{j=1}^{N_S}\mu_j N_j=\sum_{j=1}^{N_S-N_Q}\mu_j N_j+\sum_{j=N_S-N_Q+1}^{N_S}\mu_j N_j^+
\]
❓ **Q3 德文原句**：„unter Berücksichtigung der konstanten Anteile durch die Quell-/Senken-Substanzen“.
❓ **Q4 上游依赖**：式 (A.1)、(A.6)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.6)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.6)。
❓ **Q2 变量与单位**：\(\sum_{j=N_S-N_Q+1}^{N_S}\mu_j N_j^+=\mathrm{konst.}\)
❓ **Q3 德文原句**：式上的 „konstanten Anteile“.
❓ **Q4 上游依赖**：\(N_j^+\) 已知，故这一项在最小化时当作常数。
❓ **Q5 原文空白**：常数的数值原文未给出，也不需要单独给定，只要 \(N_j^+\) 与 \(\mu_j\) 已知。\(\mu_j\) 仍依赖原文未给出的 \(\mu^0\)。
❓ **Q6**：不适用。

### 式 (A.7)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.7)。
❓ **Q2 变量与单位**：
\[
\mu_j=\mu_j^0+RT\ln\left(\frac{P}{P_0}\right)+RT\ln\left(\frac{N_j}{N}\right)
\]
\(R\) 为 J/(mol·K)，\(T\) 为 K，压力比无量纲。
❓ **Q3 德文原句**：„Mit dem chemischen Potential für ideale Gase“.
❓ **Q4 上游依赖**：\(\mu_j^0(T)\)、\(T\)、\(P\)、\(P_0\)、组成。
❓ **Q5 原文空白**：\(\mu_j^0\) 原文未给出。\(P_0\) 的数值不在本页，符号表为 101300 Pa。压力符号的大小写见第 1 节。
❓ **Q6**：不适用。

### 式 (A.8)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.8)。
❓ **Q2 变量与单位**：在式 (A.5) 上加上 \(\sum_{i=1}^{N_E}\lambda_i^*\left[\sum_{j=1}^{N_S-N_Q}a_{ij}N_j-\left(b_i^0-\sum_{j=N_S-N_Q+1}^{N_S}a_{ij}N_j^+\right)\right]\)。\(\lambda_i^*\) 是 Lagrange 乘子。
❓ **Q3 德文原句**：„Mit Einführung Lagrange'scher Multiplikatoren … (Hoffmann 1971)“.
❓ **Q4 上游依赖**：式 (A.2)、(A.5)。Hoffmann 1971 没有打开，范围只是乘子方法这一句。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.9)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.9)。放大后可见三个求和号。
❓ **Q2 变量与单位**：印刷式
\[
\frac{\partial G}{\partial N_j}=\sum_{j=1}^{N_S-N_Q}\mu_j+\sum_{i=1}^{N_E}\lambda_i^*\sum_{j=1}^{N_S-N_Q}a_{ij}=0
\]
❓ **Q3 德文原句**：„Mit der Ableitung nach den Variablen \(N_j\) und \(\lambda_i^*\)“.
❓ **Q4 上游依赖**：式 (A.8)。
❓ **Q5 原文空白**：印刷式把对单个 \(N_j\) 的导数又对 \(j\) 求和，右端不再随所求的那个 \(j\) 变化，与“对变量 \(N_j\) 求导”不能同时成立。建议读法：对每个 \(j=1,\ldots,N_S-N_Q\)，
\[
\mu_j+\sum_{i=1}^{N_E}\lambda_i^* a_{ij}=0
\]
另一种读法：三个求和都保留。待总司令确认。
❓ **Q6**：不适用。

### 式 (A.10)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.10)。
❓ **Q2 变量与单位**：
\[
\frac{\partial G}{\partial\lambda_i^*}=\sum_{j=1}^{N_S-N_Q}a_{ij}N_j-\left(b_i^0-\sum_{j=N_S-N_Q+1}^{N_S}a_{ij}N_j^+\right)=0
\]
❓ **Q3 德文原句**：与 (A.9) 同一句。
❓ **Q4 上游依赖**：式 (A.2)、(A.4)。它就是元素衡算。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.11)

❓ **Q1 出处**：PDF 92 左 / S. 164，式 (A.11)。索引无坐标框，式在该半页底部。
❓ **Q2 变量与单位**：印刷式把
\[
\frac{\mu_j^0}{RT}+\ln\left(\frac{P}{P_0}\right)+\ln\left(\frac{N_j}{N}\right)
\]
从 \(j=1\) 加到 \(N_S-N_Q\)，再加上 \(\sum_i(\lambda_i^*/RT)\sum_j a_{ij}\)，令其等于 0。
❓ **Q3 德文原句**：„mit dem Produkt aus Gaskonstante und Temperatur dividiert“.
❓ **Q4 上游依赖**：式 (A.7) 代入式 (A.9)。
❓ **Q5 原文空白**：与 (A.9) 同一问题。建议读法：去掉最外层对 \(j\) 的求和，每个平衡组分一条方程。另一种读法：保留印刷的求和。
❓ **Q6**：不适用。

### 式 (A.12)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.12)。
❓ **Q2 变量与单位**：\(c_j=\mu_j^0/(RT)+\ln(P/P_0)\)。无量纲。
❓ **Q3 德文原句**：„von der Zusammensetzung unabhängigen Konstanten \(c_j\)“.
❓ **Q4 上游依赖**：\(\mu_j^0(T)\)、\(T\)、\(P\)、\(P_0\)。
❓ **Q5 原文空白**：\(c_j\) 的数值原文未给出，因为 \(\mu_j^0\) 原文未给出。
❓ **Q6**：不适用。

### 式 (A.13)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.13)。
❓ **Q2 变量与单位**：\(\lambda_i=-\lambda_i^*/(RT)\)。
❓ **Q3 德文原句**：„neuer Lagrange'scher Multiplikatoren \(\lambda_i\)“.
❓ **Q4 上游依赖**：式 (A.11) 的乘子项。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.14)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.14)。
❓ **Q2 变量与单位**：印刷式
\[
\sum_{j=1}^{N_S-N_Q}\left[c_j+\ln\left(\frac{N_j}{N}\right)\right]-\sum_{i=1}^{N_E}\lambda_i\sum_{j=1}^{N_S-N_Q}a_{ij}=0
\]
❓ **Q3 德文原句**：„folgt aus Gleichung (A.11)“.
❓ **Q4 上游依赖**：式 (A.11)–(A.13)。
❓ **Q5 原文空白**：与 (A.9) 同一问题。建议读法：对每个 \(j\)，\(c_j+\ln(N_j/N)-\sum_i\lambda_i a_{ij}=0\)。另一种读法：保留对 \(j\) 的求和。后文式 (A.15) 已经按每个 \(j\) 来写，与建议读法一致，与印刷的求和不一致。
❓ **Q6**：不适用。

### 式 (A.15)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.15)。
❓ **Q2 变量与单位**：
\[
\frac{N_j}{N}=y_j=\exp\left(\sum_{i=1}^{N_E}\lambda_i a_{ij}-c_j\right),\quad j=1,N_S-N_Q
\]
❓ **Q3 德文原句**：„nach den Molanteilen auflösen“.
❓ **Q4 上游依赖**：式 (A.14) 的建议读法。若坚持 (A.14) 的印刷求和，则推不出这一式。
❓ **Q5 原文空白**：见 (A.14)。
❓ **Q6**：不适用。

### 式 (A.16)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.16)。
❓ **Q2 变量与单位**：\(N_j=N\exp\left(\sum_{k=1}^{N_E}\lambda_k a_{kj}-c_j\right)\)，\(j=1,N_S-N_Q\)。正文说把 (A.15) 括号里的 \(i\) 换成 \(k\)，只是为了形式。
❓ **Q3 德文原句**：„Aus formalen Gründen wird \(i\) durch \(k\) … ersetzt“.
❓ **Q4 上游依赖**：式 (A.15)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.17)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.17)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}\left[N\exp\left(\sum_{k=1}^{N_E}\lambda_k a_{kj}-c_j\right)\right]-b_i^+=0,\quad i=1,N_E
\]
❓ **Q3 德文原句**：„Gleichung (A.16) wird in Gleichung (A.10) eingesetzt“.
❓ **Q4 上游依赖**：式 (A.10)、(A.16)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.18)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.18)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}\exp\left(-c_j+\sum_{k=1}^{N_E}\lambda_k a_{kj}\right)-b_i^+\,z=0
\]
正文写明此式由 (A.17) 除以总摩尔数 \(N\) 得到，因此 \(z\) 表示 \(1/N\)。没有单独的定义式 \(z:=1/N\)。
❓ **Q3 德文原句**：„Division der Gleichung (A.17) durch die Gesamtmolzahl \(N\) liefert“.
❓ **Q4 上游依赖**：式 (A.17)。未知量后来定为 \(\lambda_k\) 与 \(z\)。
❓ **Q5 原文空白**：\(z\) 不在符号表里。建议读法：\(z=1/N\)。另一种读法：\(z\) 是独立的新变量，与 \(1/N\) 无关；这与“除以 \(N\)”这句话不合。
❓ **Q6**：不适用。

### 式 (A.19)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.19)。索引无坐标框，式在该半页下部。
❓ **Q2 变量与单位**：\(\sum_{j=1}^{N_S-N_Q}y_j+\sum_{j=N_S-N_Q+1}^{N_S}y_j^+=1\)。
❓ **Q3 德文原句**：„die Summe der Stoffkonzentrationen Eins ergeben muss“. 符号是摩尔分数 \(y\)，用词是浓度。
❓ **Q4 上游依赖**：式 (A.18) 有 \(N_E\) 个方程和 \(N_E+1\) 个未知数，这一式补上闭合。
❓ **Q5 原文空白**：源汇物质的 \(y_j^+\) 如何从已知的 \(N_j^+\) 得到，要等到式 (A.20) 才写成 \(N_j^+/N\)。
❓ **Q6**：不适用。

### 式 (A.20)

❓ **Q1 出处**：PDF 92 右 / S. 165，式 (A.20)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}\exp\left(-c_j+\sum_{i=1}^{N_E}\lambda_i a_{ij}\right)+\frac{1}{N}\sum_{j=N_S-N_Q+1}^{N_S}N_j^+=1
\]
右端条件印刷为 \(i=1,N_E\)，但式中的 \(i\) 已经在内层求和里用掉。
❓ **Q3 德文原句**：„Die Gleichung für die Stoffkonzentration (A.15) in die Bedingung (A.19) eingesetzt“.
❓ **Q4 上游依赖**：式 (A.15)、(A.19)。
❓ **Q5 原文空白**：\(i=1,N_E\) 与本式没有自由指标 \(i\)。建议读法：这是一个标量闭合条件，旁边的 \(i=1,N_E\) 是从元素方程带下来的，不属于这一式。另一种读法：每个元素各有一条，但印刷式没有留下元素下标。
❓ **Q6**：不适用。

### 式 (A.21)

❓ **Q1 出处**：PDF 93 左 / S. 166，式 (A.21)。
❓ **Q2 变量与单位**：
\[
z=\frac{1}{b_1^+}\sum_{j=1}^{N_S-N_Q}a_{1j}\exp\left(-c_j+\sum_{k=1}^{N_E}\lambda_k a_{kj}\right)
\]
用的是第 1 个元素（下标 1）。
❓ **Q3 德文原句**：„Auflösen der 1. Gleichung von (A.18) nach \(z\)“.
❓ **Q4 上游依赖**：式 (A.18) 的 \(i=1\)。
❓ **Q5 原文空白**：哪一种元素必须排在第 1 位，原文未说明。\(b_1^+=0\) 时此式不能用，原文未说明。
❓ **Q6**：不适用。

### 式 (A.22)

❓ **Q1 出处**：PDF 93 左 / S. 166，式 (A.22)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}\exp\left(-c_j+\sum_k\lambda_k a_{kj}\right)-\frac{b_i^+}{b_1^+}\sum_{j=1}^{N_S-N_Q}a_{1j}\exp\left(-c_j+\sum_k\lambda_k a_{kj}\right)=0
\]
❓ **Q3 德文原句**：„\(z\) aus Gleichung (A.21) in die restlichen Gleichungen von (A.18) eingesetzt“.
❓ **Q4 上游依赖**：式 (A.18)、(A.21)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.23)

❓ **Q1 出处**：PDF 93 左 / S. 166，式 (A.23)。\(i=2,N_E\)。
❓ **Q2 变量与单位**：印刷结构是两个方括号相乘后等于 0：前一个括号里是 \(\sum_j a_{ij}-(b_i^+/b_1^+)\sum_j a_{1j}\)，后一个括号里是 \(\exp(-c_j+\sum_k\lambda_k a_{kj})\)，指数仍带自由的 \(j\)。
❓ **Q3 德文原句**：„Aus (A.22) erhält man durch Umformen“.
❓ **Q4 上游依赖**：式 (A.22)。
❓ **Q5 原文空白**：指数依赖 \(j\)，不能从对 \(j\) 的求和里无条件提出。建议读法：指数留在求和内，
\[
\sum_{j=1}^{N_S-N_Q}\left(a_{ij}-\frac{b_i^+}{b_1^+}a_{1j}\right)\exp\left(-c_j+\sum_{k=1}^{N_E}\lambda_k a_{kj}\right)=0
\]
这才是 (A.22) 的改写。另一种读法：照印刷，指数在求和之外。待总司令确认。
❓ **Q6**：不适用。

### 式 (A.24)

❓ **Q1 出处**：PDF 93 左 / S. 166，式 (A.24)，两行。
❓ **Q2 变量与单位**：第一行是 \(\sum_j\exp(-c_j+\sum_k\lambda_k a_{kj})\)。第二行是
\[
+\frac{1}{b_1^+}\sum_{j=1}^{N_S-N_Q}a_{1j}\exp(\cdots)\sum_{l=N_S-N_Q+1}^{N_S}N_l^+=1
\]
这里源汇摩尔数印刷为大写 \(N_l^+\)。
❓ **Q3 德文原句**：„Durch Einsetzen von \(z\) aus Gleichung (A.21) in die zusätzliche Bedingung (A.20)“.
❓ **Q4 上游依赖**：式 (A.20)、(A.21)。
❓ **Q5 原文空白**：与下一式的小写 \(n_l^+\) 是否同一符号，见 (A.25)。
❓ **Q6**：不适用。

### 式 (A.25)

❓ **Q1 出处**：PDF 93 左 / S. 166，式 (A.25)。
❓ **Q2 变量与单位**：印刷式把
\[
\sum_{j=1}^{N_S-N_Q}1+\frac{1}{b_1^+}\sum_{j=1}^{N_S-N_Q}a_{1j}\sum_{l=N_S-N_Q+1}^{N_S}n_l^+
\]
与 \(\exp(-c_j+\sum_k\lambda_k a_{kj})\) 相乘，等于 1。源汇量在此式印刷为小写 \(n_l^+\)。
❓ **Q3 德文原句**：„Gleichung (A.24) lässt sich vereinfachen zu“.
❓ **Q4 上游依赖**：式 (A.24)。
❓ **Q5 原文空白**：两处。(1) 指数又被放到求和外，自由指标 \(j\) 还在。建议读法：指数留在每个 \(j\) 的项内，与 (A.24) 一致。另一种读法：照印刷提出指数。(2) \(n_l^+\) 与 (A.24) 的 \(N_l^+\)、前文的 \(N_j^+\)。建议读法：同一个源汇摩尔数，字母大小写是印刷差异。另一种读法：\(n_l^+\) 是另一个量。待总司令确认。
❓ **Q6**：不适用。

### 式 (A.26) 及未编号的 \(d_i\)、\(q_{ij}\)

❓ **Q1 出处**：PDF 93 左 / S. 166。式 (A.26) 有编号；\(d_i\) 与 \(q_{ij}\) 紧接其下，没有式号。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}q_{ij}\exp\left(-c_j+\sum_{k=1}^{N_E}\lambda_k a_{kj}\right)=d_i,\quad i=1,N_E
\]
\(i=1\)：\(d_1=1\)，\(q_{1j}=1+(1/b_1^+)a_{1j}\sum_l n_l^+\)。
\(i=2,N_E\)：\(d_i=0\)，\(q_{ij}=a_{ij}-(b_i^+/b_1^+)a_{1j}\)。
求和 \(\sum_l n_l^+\) 的上下限与 (A.25) 相同。未知量是 \(\lambda_k\)，\(k=1,N_E\)。
❓ **Q3 德文原句**：„mit \(N_E\) Gleichungen für die \(\lambda_k\) Unbekannten \((k=1,N_E)\) zu lösen“.
❓ **Q4 上游依赖**：式 (A.23)、(A.25)。
❓ **Q5 原文空白**：\(q_{ij}\) 沿用了 (A.25) 的小写 \(n\) 和“指数在求和外”的结构。建议读法同 (A.25)。正文说未知量是 \(\lambda_k\) 共 \(N_E\) 个，而 \(z\) 已用 (A.21) 消去。
❓ **Q6**：不适用。

### 式 (A.27)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.27)。
❓ **Q2 变量与单位**：印刷为
\[
\frac{\partial F(\lambda)}{\partial\lambda}\,\Delta\lambda=-F(\lambda)
\]
没有矩阵括号。下文称之为 \(N_E\) 个线性方程，修正量是 \(\Delta\lambda\)。
❓ **Q3 德文原句**：„Anwendung des Newton-Verfahrens auf das Gleichungssystem (A.26)“.
❓ **Q4 上游依赖**：式 (A.26)。
❓ **Q5 原文空白**：初值不在这一式，而在 (A.35) 之后。步长如何与阻尼 \(\eta\) 配合，见 (A.30)。
❓ **Q6**：不适用。

### 式 (A.28)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.28)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}q_{ij}\left(\sum_{s=1}^{N_E}a_{sj}\right)\exp\left(-c_j+\sum_{k=1}^{N_E}\lambda_k a_{kj}\right)\Delta\lambda_s=-F(\lambda)
\]
右端没有指标，左端 \(\Delta\lambda_s\) 在求和里。
❓ **Q3 德文原句**：„\(N_E\) lineare Gleichungen zur Berechnung der Korrekturen … \(\Delta\lambda\)“.
❓ **Q4 上游依赖**：式 (A.26)、(A.27)。
❓ **Q5 原文空白**：\(\left(\sum_s a_{sj}\right)\) 被提出后不再带 \(\Delta\lambda_s\) 的对应关系，印刷上 \(\Delta\lambda_s\) 乘在整个求和之外或之内，扫描上它紧挨指数之后。建议读法：雅可比元是 \(q_{ij}\exp(\cdots)\,a_{sj}\)，对 \(s\) 求和后等于 \(-F_i\)。另一种读法：照版面，先把 \(\sum_s a_{sj}\) 收成一个标量再乘 \(\Delta\lambda_s\)。待总司令确认。
❓ **Q6**：不适用。

### 式 (A.29)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.29)。索引无坐标框。
❓ **Q2 变量与单位**：
\[
\sum_{s=1}^{N_E}\sum_{j=1}^{N_S-N_Q}q_{ij}\,a_{sj}\,y_j^{(l)}\,\Delta\lambda_s^{(l)}=d_i-\sum_{j=1}^{N_S-N_Q}q_{ij}y_j^{(l)},\quad i=1,N_E
\]
\(y_j^{(l)}\) 在此式表示迭代 \(l\) 上的 \(\exp(-c_j+\sum_k\lambda_k a_{kj})\)，不是已经归一的摩尔分数。原文没有另给定义句。
❓ **Q3 德文原句**：„Umstellen von Gleichung (A.28)“.
❓ **Q4 上游依赖**：式 (A.28)。这里的雅可比已经是 \(q_{ij}a_{sj}y_j^{(l)}\)，与 (A.28) 的建议读法一致。
❓ **Q5 原文空白**：\(y_j^{(l)}\) 的定义原文没有单列。建议读法：它就是 (A.26) 里的指数项在第 \(l\) 次迭代的值。另一种读法：它是摩尔分数 \(y_j\)。
❓ **Q6**：不适用。

### 式 (A.30)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.30)。
❓ **Q2 变量与单位**：\(\lambda_k^{(2)}=\lambda_k^{(l)}+\eta\,\Delta\lambda_k^{(l)}\)。\(\eta\) 无量纲。上标 (2) 是更新后的值，(l) 是当前迭代。
❓ **Q3 德文原句**：„empfehlen Holub und Vonka (1976) die Verwendung eines Dämpfungsparameters \(\eta\)“.
❓ **Q4 上游依赖**：式 (A.29) 的 \(\Delta\lambda\)。Holub und Vonka 1976 没有打开。范围只是这个阻尼参数。
❓ **Q5 原文空白**：只写 \(0<\eta<1\)，没有选定值，也没有说随迭代变化。选定 \(\eta\) 会改变收敛路径，不改变残差定义；本卡不选值。
❓ **Q6**：不适用。

### 式 (A.31)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.31)。
❓ **Q2 变量与单位**：\(\sum_{k=1}^{N_E}a_{kj}\lambda_k-c_j<0\)，印刷的 \(j\) 范围是 \(j=1,N_S\)。
❓ **Q3 德文原句**：„sollte so gewählt werden, daß die Ungleichung … erfüllt ist, was der Bedingung \(0\le y_j\le 1\) entspricht“.
❓ **Q4 上游依赖**：式 (A.15)：\(y_j=\exp(\sum_k\lambda_k a_{kj}-c_j)\)，所以指数为负才能使 \(y_j<1\)。\(y_j\ge 0\) 由指数函数本身保证，原文没有另写。
❓ **Q5 原文空白**：范围印刷到 \(N_S\)，平衡组分只到 \(N_S-N_Q\)。建议读法：\(j=1,\ldots,N_S-N_Q\)。另一种读法：源汇物质也要满足这条不等式。待总司令确认。
❓ **Q6**：不适用。

### 式 (A.32)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.32)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}y_j=\frac{b_i^0}{N}-\frac{1}{N}\sum_{l=N_S-N_Q+1}^{N_S}a_{il}N_l^+,\quad i=1,N_E
\]
❓ **Q3 德文原句**：„indem man die Massenbilanz (A.3) durch die Gesamtmolzahl teilt“.
❓ **Q4 上游依赖**：式 (A.3)。此处源汇量是大写 \(N_l^+\)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.33)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.33)。
❓ **Q2 变量与单位**：
\[
N=\frac{\sum_{i=1}^{N_E}\left[b_i^0-\sum_{l=N_S-N_Q+1}^{N_S}a_{il}N_l^+\right]}{\sum_{i=1}^{N_E}\sum_{j=1}^{N_S-N_Q}a_{ij}y_j}
\]
❓ **Q3 德文原句**：„über alle Elemente summiert und nach \(N\) auflöst“.
❓ **Q4 上游依赖**：式 (A.32)。
❓ **Q5 原文空白**：分母为零时原文未说明。
❓ **Q6**：不适用。

### 式 (A.34)

❓ **Q1 出处**：PDF 93 右 / S. 167，式 (A.34)。
❓ **Q2 变量与单位**：
\[
D=\sum_{i=1}^{N_E}\left|1-\sum_{j=1}^{N_S-N_Q}a_{ij}\frac{N_j}{b_i^+}\right|<\varepsilon
\]
\(D\) 与 \(\varepsilon\) 无量纲。
❓ **Q3 德文原句**：„Als Abbruchkriterium geben Holub und Vonka (1976) folgende Größe an“.
❓ **Q4 上游依赖**：元素衡算 (A.3)。Holub und Vonka 1976 没有打开。范围只是这个判据的形式。
❓ **Q5 原文空白**：\(\varepsilon\) 的数值原文未给出。\(b_i^+=0\) 时分母无定义，原文未说明。式内是 \(N_j\) 不是 \(y_j\)。
❓ **Q6**：不适用。

### 式 (A.35)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.35)。上一页末句 „Aus Gleichung (A.15) folgt“.
❓ **Q2 变量与单位**：\(\ln y_j=\sum_{i=1}^{N_E}\lambda_i a_{ij}-c_j\)，\(j=1,N_S-N_Q\)。
❓ **Q3 德文原句**：„Zur vollständigen Beschreibung der Methode fehlt noch die Bestimmung der Startwerte“.
❓ **Q4 上游依赖**：式 (A.15)。这是初值段的起点，不是另一套平衡条件。
❓ **Q5 原文空白**：初值 \(y_j^{(l)}\) 取什么，下一式才出现，仍然没有给数值。
❓ **Q6**：不适用。

### 式 (A.36)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.36)。
❓ **Q2 变量与单位**：
\[
\sum_{i=1}^{N_E}\lambda_i a_{ij}-c_j=\ln y_j\approx\ln y_j^{(l)}+\frac{1}{y_j^{(l)}}\left(y_j-y_j^{(l)}\right)
\]
❓ **Q3 德文原句**：„Eine Taylorentwicklung der linken Seite von Gleichung (A.35) um \(y_j^{(l)}\)“. 印刷把展开对象说成 (A.35) 的左端，展开变量却是 \(y_j\)。左端 \(\ln y_j\) 才是 \(y_j\) 的函数。
❓ **Q4 上游依赖**：式 (A.35)。
❓ **Q5 原文空白**：建议读法：对 \(\ln y_j\) 在 \(y_j^{(l)}\) 处展开，式子本身如此。另一种读法：严格按“左端”一词，但左端在 (A.35) 里已经等于 \(\ln y_j\)。不影响后式。\(y_j^{(l)}\) 的初值原文未给出。
❓ **Q6**：不适用。

### 式 (A.37)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.37)。
❓ **Q2 变量与单位**：
\[
y_j\approx y_j^{(l)}\left(1-\ln y_j^{(l)}\right)+y_j^{(l)}\left(\sum_{k=1}^{N_E}\lambda_k a_{kj}-c_j\right)
\]
\(j=1,N_S-N_Q\)。与 (A.36) 的代数展开一致。
❓ **Q3 德文原句**：„Daraus ergibt sich“.
❓ **Q4 上游依赖**：式 (A.36)。
❓ **Q5 原文空白**：初值仍未给出。
❓ **Q6**：不适用。

### 式 (A.38)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.38)。
❓ **Q2 变量与单位**：\(\sum_{j=1}^{N_S-N_Q}a_{ij}y_j=b_i^+/N\)，\(i=1,N_E\)。
❓ **Q3 德文原句**：„Die Massenbilanz (A.3) unter Verwendung der Molanteile“.
❓ **Q4 上游依赖**：式 (A.3)。这里没有再写源汇项，因为 \(b_i^+\) 已经扣过。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.39)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.39)。
❓ **Q2 变量与单位**：
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}\left[y_j^{(l)}(1-\ln y_j^{(l)})+y_j^{(l)}\left(\sum_{k=1}^{N_E}\lambda_k a_{kj}-c_j\right)\right]=\frac{b_i^+}{N}
\]
❓ **Q3 德文原句**：„durch einsetzen von Gleichung (A.37) in Gleichung (A.38)“.
❓ **Q4 上游依赖**：式 (A.37)、(A.38)。
❓ **Q5 原文空白**：无。
❓ **Q6**：不适用。

### 式 (A.40)

❓ **Q1 出处**：PDF 94 左 / S. 168，式 (A.40)。附录 A1 的最后一式。
❓ **Q2 变量与单位**：印刷为
\[
\sum_{j=1}^{N_S-N_Q}a_{ij}y_j^{(l)}\sum_{k=1}^{N_E}\lambda_k a_{kj}=b_i^+-\sum_{j=1}^{N_S-N_Q}a_{ij}y_j^{(l)}(1-\ln y_j^{(l)})+\sum_{j=1}^{N_S-N_Q}a_{ij}y_j^{(l)}c_j
\]
放大后右端第一项是 \(b_i^+\)，没有除以 \(N\)。
❓ **Q3 德文原句**：„und umformen der Gleichung (A.39) ergibt sich“. 页末说明 (A.40) 的求解仍用式 (A.27) 的 Newton 法。
❓ **Q4 上游依赖**：式 (A.39) 的右端是 \(b_i^+/N\)。
❓ **Q5 原文空白**：由 (A.39) 移项，右端第一项应是 \(b_i^+/N\)，左端才与摩尔分数同量纲。印刷是 \(b_i^+\)。建议读法：\(b_i^+/N\)。另一种读法：照印刷为 \(b_i^+\)。待总司令确认。初值 \(y_j^{(l)}\) 原文未给出。
❓ **Q6**：不适用。

---

## 4. 附录 A1 在主模型中的用途

### 热解气组成（S. 69–70）

❓ **Q1 出处**：PDF 44 右 / S. 69，式 (4.26)、(4.27)；PDF 45 左 / S. 70，式 (4.28)–(4.31)。指向附录的句子在 S. 70。
❓ **Q2 变量与单位**：式 (4.26)、(4.27) 是 waf 基准的 kg/kg。式 (4.28) \(G=\sum_{j=1}^{N_S}\mu_j n_j=\mathrm{Min!}\)。式 (4.30) \(\sum_j a_{ij}n_j-b_i=0\)。式 (4.31) 是带乘子的 \(F\)。这里的摩尔数符号是 \(n_j\)，附录用 \(N_j\)。
❓ **Q3 德文原句**：„Mit der Vorgabe der in der Gasphase zu berücksichtigenden Substanzen, z.B. CO, CO₂, CH₄, C₂H₄, NH₃, H₂S, H₂O, H₂, erfolgt die Berechnung der Gleichgewichtszusammensetzung der Gasphase auf der Basis der Arbeit von Holub und Vonka (1976) durch Minimierung der Freien Enthalpie.“ S. 70：„Eine detaillierte Herleitung dieser Methode … ist im Anhang A1 zu finden.“
❓ **Q4 上游依赖**：焦炭和焦油的元素量先用式 (4.26)、(4.27) 从燃料元素里扣掉。其余元素 „gehen vollständig in die Gasphase über“.
❓ **Q5 原文空白**：
- „z.B.“ 表示举例，不是封闭的物质清单。
- 最小化时的温度、压力只写了 „bei konstantem Druck und konstanter Temperatur“，没有说等于床层温度、颗粒温度或某一操作压力。
- 第 4 章没有使用附录的符号 \(N_Q\)。焦炭和焦油是在进平衡之前从元素衡算里扣掉的，不是写成气相里的已知源汇摩尔数 \(N_j^+\)。可选解释：(a) 扣掉之后的 \(b_i\) 就是附录的 \(b_i^+\)，不再另设源汇；(b) 焦炭和焦油就是 \(N_Q\) 种源汇物质。原文没有把两者写成同一句话。
- Holub und Vonka 1976 没有打开。范围就是这段所写的气相自由焓最小化。
❓ **Q6**：不适用。

### 其他调用

❓ **Q1 出处**：正文里写出 „Anhang A1“ 的，是 S. 70 这一句。目录在 PDF 5。S. 167 用 Holub und Vonka (1976) 说明阻尼和收敛，没有再写“见附录”。
❓ **Q2 变量与单位**：不适用。
❓ **Q3 德文原句**：见上。
❓ **Q4 上游依赖**：PDF 62、69、96、97 的 „Anhang“ 指的是附录里的浓度图（A2.4 一类），不是 A1 的平衡算法。
❓ **Q5 原文空白**：式 (5.43) 使用 \(\mu_i^0\)，但没有写“见附录 A1”，也没有写与附录共用一套数据。见式 (5.43) 的记录。
❓ **Q6**：不适用。

---

## 5. 论文锚点

只列出可读的量和范围，不评价。

### S. 71 反应焓

❓ **Q1 出处**：PDF 45 右 / S. 71。
❓ **Q2 变量与单位**：第 2 节那张表。状态 0,1013 MPa、298,15 K。R1 是区间，R10 与 R11 只有符号。
❓ **Q3 德文原句**：„Reaktionsenthalpie (0,1013 MPa; 298,15 K)“.
❓ **Q4 上游依赖**：生成焓若一旦选定，这些数是直接的检验量。原文没有给出生成焓。
❓ **Q5 原文空白**：R10、R11 不能作数量锚点。
❓ **Q6**：不适用。

### Bild 5.11

❓ **Q1 出处**：PDF 56 右 / S. 93。数字化文件 `data/digitized/fig_5_11.json`（2026-10-16）。
❓ **Q2 变量与单位**：横轴自由空间温度，K。纵轴是按式 (5.43) 计算的均相水煤气反应平衡常数，无量纲。文件中的点是实线，不是各文献的 \(K_{exp}\) 散点。
❓ **Q3 德文原句**：图题 „Berechnete Gleichgewichtskonstante der homogenen Wassergasreaktion … in Abhängigkeit der Freiraumtemperatur“.
❓ **Q4 上游依赖**：式 (5.43) 的 \(\mu_i^0\)。图内另标反应焓 −41,1 kJ/mol，与 S. 71 的 R8 相同。
❓ **Q5 原文空白**：曲线约从 1010 K 处的 1,4 进入图框。文件取点 1025–1250 K，\(K_{eq}\) 从 1,314 到 0,648；读数误差 \(T\) 为 3 K，\(K_{eq}\) 为 0,01。1025 K 以下不取点。这些是读数文件里的范围，不是 Hamel 另给的表。
❓ **Q6**：不适用。

### 可作对照、但混有动力学的图

❓ **Q1 出处**：Bild 4.2、4.3，PDF 37 右 / S. 55：湿褐煤的水分与挥发分对时间。Bild 4.4–4.6，PDF 38 / S. 56–57：木材挥发分释放。
❓ **Q2 变量与单位**：质量分数对时间（s）。图中的床温、表观速度、粒径是该次比较的工况，不是物性常数。
❓ **Q3 德文原句**：Bild 4.2 图题 „Vergleich des gemessenen und berechneten Wasser- und Flüchtigengehalts“.
❓ **Q4 上游依赖**：同时依赖表 4.1 或表 4.2 的物性和热解动力学。
❓ **Q5 原文空白**：不能单独检验生成焓或 \(K_{eq}\)。数字化不在本卡。
❓ **Q6**：不适用。

---

## 汇总：物性量

“有依据”指原文给出了数值或关联式。“原文未给出”指模型用得到、但没有数值或闭合式。“有歧义”指用法或符号还有待裁定。单位均为 SI。

| 量 | 符号 | 单位 | 使用位置 | 原文数据或来源 |
|---|---|---|---|---|
| 自由焓 | \(G\) | J | (4.28)、(A.1) | 有公式，无数值 |
| 比焓 | \(h\) | J/kg 或 J/mol | (2.8) | 原文未给出。含标准生成焓，参考态未写。与高度符号 \(h\) 重名 |
| 水的蒸发焓 | \(h_v\) | J/kg | (4.4) | 原文未给出 |
| 修正蒸发焓 | \(h'_v\) | J/kg | (4.3)、(4.4) | 有公式，\(h_v\) 与 \(T_e\) 未给出 |
| 比热 | \(c_P\) | J/(kg·K) | 符号表；干燥模型 | 褐煤 1256、木材 1670，只标明用于干燥/热解颗粒。气体比热原文未给出。是否用于 (2.8) 有歧义 |
| 固体比热 | \(c_{P,s}\) | J/(kg·K) | 符号表；式 (4.4) 印刷为 \(c_s\) | 与 \(c_P\) 是否同一数，原文未说明 |
| 水的比热 | \(c_{P,w}\) | J/(kg·K) | 符号表；式 (4.4) 印刷为 \(c_w\) | 数值原文未给出。与阻力系数 \(c_W\) 不同符号 |
| 标准化学势 | \(\mu_j^0\) | J/mol | (5.43)、(A.7)、(A.12) | 原文未给出。与附录是否同一套，原文未明说 |
| 化学势 | \(\mu_j\) | J/mol | (A.1)、(A.7) | 有公式，依赖 \(\mu_j^0\) |
| 参考压力 | \(P_0\) | Pa | 符号表；(A.7)；(5.20) | 有依据：101300 Pa。全书是否统一，见第 1 节 |
| 反应焓参考态 | — | 0,1013 MPa；298,15 K | S. 71 | 有依据，只标在反应焓表 |
| 通用气体常数 | \(R_g\)（式中常写 \(R\)） | J/(mol·K) | (5.43)、(A.7) | 原文未给出数值 |
| 摩尔质量 | \(M\) | kg/mol | 符号表；式 (5.11) 的碳 | 各组分数值原文未给出 |
| 气体密度 | \(\rho_g\)、\(\rho_{g,20}\) | kg/m³ | (3.20)、(5.19) | 原文未给出。20 °C 也没有数值 |
| 固体密度 | \(\rho_P\) | kg/m³ | 表 4.1、4.2；符号表 | 褐煤干颗粒 1250；木材 500。焦炭、灰的通用值原文未给出 |
| 惰性物料密度 | \(\rho_{inert}\) | kg/m³ | 表 5.1 算例 | 2000，只属于该算例 |
| 气体动力粘度 | \(\mu_g\)、\(\mu_{g,20}\) | Pa·s | (5.19) | 原文未给出 |
| 气体运动粘度 | \(\nu_g\) | m²/s | (5.18) | 原文未给出 |
| 悬浮相粘度 | \(\mu_d\)、\(\nu_d\) | Pa·s，m²/s | (3.19)–(3.22) | 有依据：Hilligardt 关联式。单位指数建议 5/3、5/2 |
| 氧气扩散系数 | \(D_g\) | m²/s | (5.16)、(5.19)、(5.20) | 有依据：Field et al. (1967) 的那一支。其他气体原文未给出 |
| 导热系数 | \(\lambda\)、\(\lambda_s\)、\(\lambda_g\) | W/(m·K) | (4.3)、(4.8)；表 4.2 | 木材颗粒 0,1256。气体与褐煤的 \(\lambda\) 原文未给出 |
| 导温系数 | \(a\) | m²/s | 表 4.1 | 褐煤干颗粒 \(1\cdot 10^{-7}\)。只标明用于干燥/热解 |
| 颗粒传热系数 | \(\alpha\) | W/(m²·K) | (4.9)、(5.22) | 有 Nu 关联式。Re、Pr 与气体物性原文未给出 |
| 壁面传热系数 | — | — | 表 2.1 输入；式 (2.7) 的 \(\dot Q_W\) | 装置输入。原文未给出物性关联式。本卡不把它当作待选物性 |
| 颗粒发射率 | \(\varepsilon_P\) | — | (5.22) | 有依据：\(0{,}79\le\varepsilon_P\le 0{,}91\)。区间内取值原文未说明 |
| 斯特藩–玻尔兹曼常数 | \(\sigma_S\)；式 (5.22) 印刷为 \(\sigma\) | W/(m²·K⁴) | (5.22) | 数值原文未给出。下标有歧义 |
| 蒸发温度 | \(T_e\) | K | (4.4)；符号表 | 原文未给出 |
| 燃料热值、元素分析 | — | 热值为输入量 | 表 2.2 清单 | 原文未给出由热值换算生成焓的式子 |
| 焦油模型物质 | 苯、萘；生物质为萘与十六烷 | — | S. 95 | 物质名称有依据。生成焓与比热原文未给出 |
| 焦炭生成焓与比热 | — | — | 式 (2.8) 需要固体焓 | 原文未给出 |
| 惰性床料比热 | — | J/(kg·K) | 若 \(h_s\) 是床料 | 原文未给出。\(s\) 的含义不在本卡裁定 |
| 理想气体 | — | — | 附录 A1；式 (5.43) | 平衡计算有依据。\(\rho_g\) 的状态方程原文未给出 |

符号表中的物性类符号都在上表。阻力系数 \(c_b\)、\(c_W\) 不是物性，不列入。高度符号 \(h\) 与比焓同字，已在比焓一行说明。

## 汇总：论文锚点

| 图或表 | 页码 | 可读的量与范围 |
|---|---|---|
| S. 71 反应焓表 | PDF 45 右 / S. 71 | R1 区间 −393,5 至 −110,5 kJ/mol C；R2 +131,4；R3 +172,5；R4 −74,9 kJ/mol C；R5 −283 kJ/mol CO；R6 −241,9 kJ/mol H₂；R7 −519,4 kJ/mol CH₄；R8 −41,1 kJ/mol CO；R9 +206,3 kJ/mol CH₄。状态 0,1013 MPa、298,15 K。R10、R11 只有正负号 |
| Bild 5.11 实线 | PDF 56 右 / S. 93；`data/digitized/fig_5_11.json` | \(K_{eq}(T)\)，1025–1250 K，由 1,314 降至 0,648。读数误差 3 K 与 0,01。1025 K 以下、图框外不取点 |
| Bild 4.2、4.3 | PDF 37 右 / S. 55 | 水分与挥发分对时间。混有动力学，不是纯物性锚点 |
| Bild 4.4–4.6 | PDF 38 / S. 56–57 | 木材挥发分释放。同上 |

## 交给总司令确认的读法

物理上可以给出建议、但仍保留另一种读法的，有：

1. 式 (A.9)、(A.11)、(A.14)：建议去掉对 \(j\) 的外层求和。
2. 式 (A.23)、(A.25) 及 \(q_{ij}\) 的写法：建议指数留在对 \(j\) 的求和内。
3. 式 (A.25) 的 \(n_l^+\) 与 (A.24) 的 \(N_l^+\)：建议是同一个源汇摩尔数。
4. 式 (A.28) 的雅可比：建议是 \(q_{ij}\exp a_{sj}\)，与 (A.29) 一致。
5. 式 (A.31) 的 \(j\) 上限：建议到 \(N_S-N_Q\)。
6. 式 (A.40) 右端：建议 \(b_i^+/N\)。
7. 式 (3.21)、(3.22) 的单位指数：建议 5/3 与 5/2；\(d_P\) 建议用 m。
8. 表 4.1、表 4.2 的比热和密度是否进入式 (2.8)：原文未说明。
9. 式 (5.43) 与附录 A1 是否同一套 \(\mu^0\)：原文未明说。数值两边都未给出。

没有一处字形是物理上也无法判断、因而必须交给目标所有者肉眼判读的。上面 9 条都还有一种印刷读法，确认前不选数据来源，也不写代码。
