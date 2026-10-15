全新会话，转录冻结前未读旧材料。

# T8 原文答卷：§3.1 床层流体力学

- 任务卡：T8
- 执行者：cursor
- 会话：全新会话，转录冻结前未读旧材料
- 日期：2026-10-15
- 语料：`refs/hamel/` 页面图像。公式以页面图像为准。

范围：§3.1.1–§3.1.4，印刷页 S. 21–35。式号 3.1–3.48，以及 S. 35 上的式 (3.49)、式 (3.50)。式 (3.49)、式 (3.50) 在 S. 36 另有一组，属于 §3.2.1，索引标记 `ambiguous_number`。引用这四个式号必须带印刷页。§3.2 不展开，文末只做式号区分。

§3.1 没有出现第 2 章的 \(\dot m_{auf}\)、\(\dot m_{ab}\)、\(\dot m_{aus}\)。S. 21 说流动模型确定向上和向下的气固流，但本节印出的是速度、含率、气泡尺寸和相间交换系数。与式 (2.6) 的关系见 `ch2_balances.md` 的头号问题，这里不另作选择。

Q6 在冻结时一律写“比对阶段填写”。

速度在原文给出数字时单位为 m/s。\(K_{bd}\) 在 S. 35 以 s\(^{-1}\) 出现。其余单位在各式下写明；原文未写的，记为原文未说明。

## Eq. 3.1（Ergun，最小流化）

❓ **Q1 出处**：PDF 20 右半页 / S. 21，§3.1.1，Eq. 3.1。
❓ **Q2 变量与单位**：

\[
Ar=150\cdot\frac{1-\varepsilon_{mf}}{\phi_S^2\cdot\varepsilon_{mf}^3}\cdot Re_{mf}+\frac{1{,}75}{\phi_S\cdot\varepsilon_{mf}^3}\cdot Re_{mf}^2
\]

\(Ar\)、\(Re_{mf}\) 由式 (3.2) 定义，无量纲。\(\varepsilon_{mf}\) 最小流化空隙率，无量纲。\(\phi_S\) 球形度，无量纲。第一项分母是 \(\phi_S^2\)，第二项分母是 \(\phi_S\)。
❓ **Q3 德文原句**：„Ergun (1952) gibt zur Berechnung der minimalen Fluidisierungsgeschwindigkeit folgende … Gleichung an“。
❓ **Q4 上游依赖**：Hamel 引用的 Ergun（1952）（式 3.1，S. 21）。\(Ar\)、\(Re_{mf}\) 见式 (3.2)。
❓ **Q5 原文空白**：本式把 \(Ar\) 与 \(Re_{mf}\) 连在一起，没有单独解出 \(u_{mf}\)。求解形式见式 (3.4)。\(\varepsilon_{mf}\) 与 \(\phi_S\) 的取值原文在本式处未给。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.2（Ar 与 Re 的定义）

❓ **Q1 出处**：PDF 20 右半页 / S. 21，§3.1.1，Eq. 3.2。
❓ **Q2 变量与单位**：

\[
Ar=\frac{g\cdot\bar d_P^3}{\nu_g^2}\cdot\frac{\rho_P-\rho_g}{\rho_g},\qquad Re_{mf}=\frac{d_P\cdot u_{mf}}{\nu_g}
\]

\(Ar\) 中的粒径是 \(\bar d_P\)（有平均号）。\(Re_{mf}\) 中的粒径是 \(d_P\)（无平均号）。\(g\) 重力加速度，SI 为 m/s\(^2\)，原文未写。\(\nu_g\) 气体运动粘度，SI 为 m\(^2\)/s，原文未写。\(\rho\) 密度，SI 为 kg/m\(^3\)，原文未写。\(u_{mf}\) 最小流化空塔速度，后文数字用 m/s。
❓ **Q3 德文原句**：„mit folgender Definition der Archimedes-Zahl und der Reynolds-Zahl“。
❓ **Q4 上游依赖**：供式 (3.1) 使用。
❓ **Q5 原文空白**：同一式中一处有平均号、一处没有。可选解释：(a) \(Re_{mf}\) 也应用平均粒径，平均号漏印（影响：两种粒径混用会改变 \(Re\)）；(b) 两处故意不同（影响：须另有一个不加平均号的 \(d_P\)）。待总司令裁定。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.3（由环境条件下的最小流化反求球形度）

❓ **Q1 出处**：PDF 21 左半页 / S. 22，§3.1.1，Eq. 3.3。式号在段末，公式在下一段。
❓ **Q2 变量与单位**：

\[
\phi_S\cdot\bar d_P=\frac{1{,}75\cdot\rho_g\cdot u_{mf,0}^2}{2\cdot\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}+\sqrt{\frac{150\cdot(1-\varepsilon_{mf})\cdot\mu_g\cdot u_{mf,0}}{\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}+\left(\frac{1{,}75\cdot\rho_g\cdot u_{mf,0}^2}{2\cdot\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}\right)^2}
\]

\(u_{mf,0}\)：环境条件下测得的最小流化速度。\(\mu_g\) 气体动力粘度。原文未写单位。长度与速度用 SI 时，左端为 m。
❓ **Q3 德文原句**：„kann zunächst die Sphärizität des Wirbelgutes bestimmt werden“。
❓ **Q4 上游依赖**：Hamel 引用的 Werther（1977）（S. 22）：先在环境条件下测定松动速度，再外推到操作条件。本式用测得的松动速度、对应空隙率、平均粒径和环境条件下的物性求球形度。
❓ **Q5 原文空白**：环境条件的温度和压力原文未给数值。“外推”的具体步骤除本式和式 (3.4) 外没有再写。\(\varepsilon_{mf}\) 写的是“所属的空隙率”，未说明是测量值还是固定值。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.4（操作条件下的最小流化速度）

❓ **Q1 出处**：PDF 21 左半页 / S. 22，§3.1.1，Eq. 3.4。索引 `locator: page_image`，无标签框；公式在式 (3.3) 的说明之后。
❓ **Q2 变量与单位**：

\[
u_{mf}=42{,}9\cdot\frac{(1-\varepsilon_{mf})\cdot\nu_g}{\phi_S\cdot\bar d_P}\left(\sqrt{1+\frac{3{,}1\cdot 10^{-4}\cdot\varepsilon_{mf}^3\cdot g\cdot(\phi_S\cdot\bar d_P)^3\cdot(\rho_P-\rho_g)}{(1-\varepsilon_{mf})^2\cdot\nu_g^2\cdot\rho_g}}-1\right)
\]

\(u_{mf}\) 为速度。常数 42,9 与 \(3{,}1\cdot 10^{-4}\) 按原页读出，原文未写它们的单位；与式 (3.1) 联立时二者无量纲。
❓ **Q3 德文原句**：„Mit bekannter Sphärizität des Wirbelgutes erhält man die minimale Fluidisierungsgeschwindigkeit bei Betriebsbedingungen.“
❓ **Q4 上游依赖**：\(\phi_S\) 来自式 (3.3)。物性为操作条件下的气体物性这一句，原文写的是“操作条件下的最小流化速度”，没有逐个声明 \(\nu_g\)、\(\rho_g\) 必须换成操作条件。
❓ **Q5 原文空白**：42,9 与 \(3{,}1\cdot 10^{-4}\) 的有效数字以原页为准，原文未给出更长的小数。原文未说明本式与直接解式 (3.1) 是否应视为同一关系。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.5（Ergun 型关联的一般形式）

❓ **Q1 出处**：PDF 21 左半页 / S. 22，§3.1.1，Eq. 3.5。
❓ **Q2 变量与单位**：

\[
\frac{Ar}{Re_{mf}}=C_1+C_2\cdot Re_{mf}^{m}
\]

\(C_1\)、\(C_2\)、\(m\) 为系数。原文给出文献中 \(C_1\) 约 18 至 34，\(C_2\) 约 0,0356 至 0,067（Verweyen 1994）。\(m\) 没有给数值。
❓ **Q3 德文原句**：„folgende allgemeine Form aufweisen (Lippens und Mulder 1993)“。
❓ **Q4 上游依赖**：Hamel 引用的 Lippens und Mulder（1993）（式 3.5，S. 22）。\(C_1\)、\(C_2\) 的范围：Hamel 引用的 Verweyen（1994）（S. 22）。
❓ **Q5 原文空白**：指数 \(m\) 原文未给数值。Hamel 没有从上述常数对里选定一对。同页写 Ergun 的结果最好，并写 Olowson und Almstedt（1991）、King und Harrison（1982）表明较高压力下仍可用 Ergun。这是对文献的评述，不是一条新的计算式。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.6（气泡特征数）

❓ **Q1 出处**：PDF 21 右半页 / S. 23，§3.1.2，Eq. 3.6。
❓ **Q2 变量与单位**：

\[
\alpha_b=\frac{u_b}{u_d}\cdot\varepsilon_b
\]

\(u_b\) 气泡上升速度，\(u_d\) 悬浮相内气体的实际速度，\(\varepsilon_b\) 气泡相份额。三者中速度的单位原文在后文用 m/s；\(\alpha_b\)、\(\varepsilon_b\) 无量纲。\(\alpha_b<1\) 为慢气泡，\(\alpha_b>1\) 为快气泡。
❓ **Q3 德文原句**：„als das Verhältnis der Blasenaufstiegsgeschwindigkeit \(u_b\) zur tatsächlichen Geschwindigkeit des Gases … \(u_d\) definiert ist“。
❓ **Q4 上游依赖**：Hamel 引用的 Davidson und Harrison（1963）（S. 23）作为两种流型的来源。Geldart 分类：Hamel 引用的 Geldart（1973）。慢/快气泡的文字说明见同页，图见 Bild 3.1（PDF 22 左半页 / S. 24）。
❓ **Q5 原文空白**：定义句只说 \(u_b\) 与 \(u_d\) 之比，印刷式在这个比之外还乘了 \(\varepsilon_b\)。可选解释：(a) 以印刷式为准，定义句省略了 \(\varepsilon_b\)（影响：慢/快的判据含气泡份额）；(b) 以定义句为准，\(\varepsilon_b\) 是误印（影响：判据只是 \(u_b/u_d\)）。待总司令裁定。原文未说明 \(\alpha_b=1\) 归入哪一类。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.7（床层截面气体流量衡算）

❓ **Q1 出处**：PDF 22 左半页 / S. 24，§3.1.2，Eq. 3.7。
❓ **Q2 变量与单位**：

\[
u_0\cdot A=\dot V_b+\dot V_{b,r}+\dot V_d
\]

\(u_0\) 空塔速度，\(A\) 截面积，\(\dot V\) 体积流量。可见气泡流、相对气泡的穿流、悬浮相穿流。原文假定气体不可压缩，并忽略悬浮相气体返混。单位：速度 m/s 时体积流量为 m\(^3\)/s，原文未写出。
❓ **Q3 德文原句**：„folgende Bilanzgleichung (siehe Bild 3.2)“。
❓ **Q4 上游依赖**：Hamel 引用的 Clift und Grace（1985）（S. 24）的三种气体输运。图：Bild 3.2，PDF 22 右半页 / S. 25。
❓ **Q5 原文空白**：不可压缩的适用范围（压力、床高）原文未给数字界限。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.8（可见气泡气体流）

❓ **Q1 出处**：PDF 22 左半页 / S. 24，§3.1.2，Eq. 3.8。
❓ **Q2 变量与单位**：\(\dot V_b=u_b\cdot A\cdot\varepsilon_b\)。体积流量。
❓ **Q3 德文原句**：„Für den sichtbaren Blasengasstrom gilt“。
❓ **Q4 上游依赖**：代入式 (3.7)。
❓ **Q5 原文空白**：无额外空白。\(u_b\) 与 \(\varepsilon_b\) 的计算式在后文。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.9（相对气泡的穿流）

❓ **Q1 出处**：PDF 22 左半页 / S. 24，§3.1.2，Eq. 3.9。
❓ **Q2 变量与单位**：\(\dot V_{b,r}=u_{b,r}\cdot A\cdot\varepsilon_b\)。\(u_{b,r}\) 为气体相对气泡的速度。
❓ **Q3 德文原句**：„für die Durchströmung relativ zur Blase“。
❓ **Q4 上游依赖**：\(u_{b,r}\) 的数值关系见式 (3.23)，加压时见式 (3.44)。
❓ **Q5 原文空白**：无。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.10（悬浮相穿流）

❓ **Q1 出处**：PDF 22 左半页 / S. 24，§3.1.2，Eq. 3.10。
❓ **Q2 变量与单位**：\(\dot V_d=u_d\cdot A\cdot(1-\varepsilon_b)\)。
❓ **Q3 德文原句**：„und für die Durchströmung der Suspensionsphase“。
❓ **Q4 上游依赖**：Bild 3.2 将 \((1-\varepsilon_b)\) 标为密相份额。
❓ **Q5 原文空白**：\(u_d\) 是空塔基准还是密相真实速度，式 (3.6) 的定义句称“悬浮相内气体的实际速度”，式 (3.10) 又把 \(u_d\) 乘 \((1-\varepsilon_b)\) 得到密相体积流量。两种用法是否同一 \(u_d\)，原文没有再声明。可选解释：(a) \(u_d\) 始终是按全截面定义的密相贡献速度，式 (3.6) 的“实际速度”用词较宽（影响：式 (3.6) 的分母不是颗粒间隙速度）；(b) 两处 \(u_d\) 不是同一个量（影响：不能把式 (3.10) 直接代进式 (3.6)）。待总司令裁定。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.11（Hilligardt，悬浮相速度）

❓ **Q1 出处**：PDF 22 右半页 / S. 25，§3.1.2，Eq. 3.11。
❓ **Q2 变量与单位**：

\[
u_d=u_{mf}+\frac13(u_0-u_{mf})
\]

速度，m/s。适用范围写明：Geldart D，工业布风器，三维床，按空截面计的过剩气速至 0,5 m/s。
❓ **Q3 德文原句**：„fand Hilligardt für technische Gasverteiler in dreidimensionalen Wirbelschichten“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.11，S. 25）。Geldart（1973）用于颗粒分类。
❓ **Q5 原文空白**：与式 (3.12) 并列，原文未写本文对给定颗粒选用哪一式。见式 (3.12)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.12（Wein，悬浮相速度）

❓ **Q1 出处**：PDF 22 右半页 / S. 25，§3.1.2，Eq. 3.12。
❓ **Q2 变量与单位**：\(u_d=1{,}45\cdot u_{mf}\)。Geldart B。Wein 的空塔速度至 3 m/s；沙与灰的过剩气速为 2,87 与 2,97 m/s。原文写：悬浮相速度高于最小流化速度，且不随过剩气速再升高。
❓ **Q3 德文原句**：„Seine Ergebnisse korreliert er folgendermaßen“。
❓ **Q4 上游依赖**：Hamel 引用的 Wein（1992）（式 3.12，S. 25）。
❓ **Q5 原文空白**：式 (3.11) 与式 (3.12) 都印出，没有选择句。可选解释：(a) 按 Geldart 类别分用（影响：B 与 D 的 \(u_d\) 不同）；(b) 只采用其中一式（影响：另一式成为未使用的引述）。待总司令裁定。1,45 的有效数字以原页为准。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.13（气泡上升速度）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.13。索引 `locator: page_image`。
❓ **Q2 变量与单位**：

\[
u_b=\psi_b\cdot(u_0-u_{mf})+\upsilon_b\cdot u_{b,i}
\]

\(\psi_b\)：进入可见气泡流的过剩气体份额。\(\upsilon_b\)：气泡群相对单气泡的修正。\(u_{b,i}\) 单气泡速度，见式 (3.17)。速度为 m/s，两个系数无量纲。\(\upsilon_b\) 的字形见文末待肉眼确认（读作与 \(\psi\) 同类的希腊字母，不是拉丁字母 \(v\)）。
❓ **Q3 德文原句**：„Hilligardt und Werther (1985) haben eine Korrelation für die Blasenaufstiegsgeschwindigkeit vorgeschlagen“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt und Werther（1985）（式 3.13，S. 26）。加压修正见式 (3.43)。
❓ **Q5 原文空白**：\(\psi_b\) 有式 (3.14) 与式 (3.15) 两套，原文未选定。系数的字母形状待肉眼确认。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.14（Hilligardt 的 \(\psi_b\)）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.14。
❓ **Q2 变量与单位**：B 类颗粒，Hilligardt（1986）：工业布风器 \(\psi_b=0{,}76\)；多孔板 \(\psi_b=0{,}67\)。无量纲。
❓ **Q3 德文原句**：„beträgt für Gruppe-B-Partikel nach Hilligardt (1986)“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.14，S. 26）。
❓ **Q5 原文空白**：同页接着给出 Wein 的关联并写可以不再区分布风器。本文用哪一个，未写。A、D 类的 \(\psi_b\) 原文未给。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.15（Wein 的 \(\psi_b\)）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.15。
❓ **Q2 变量与单位**：

\[
\psi_b=0{,}17\cdot u_{mf}^{-0{,}33}\qquad\text{für }0{,}014\,\mathrm{m/s}<u_{mf}<0{,}18\,\mathrm{m/s}
\]

\(u_{mf}\) 以 m/s 代入时 \(\psi_b\) 无量纲。原文写该参数这时只还取决于所加入的固体。
❓ **Q3 德文原句**：„Wein (1992) führt folgende Korrelation für den mittleren Parameter \(\psi_b\) … ein“。
❓ **Q4 上游依赖**：Hamel 引用的 Wein（1992）（式 3.15，S. 26）。
❓ **Q5 原文空白**：区间端点是否包含，原文用的是严格不等号。区间外的 \(\psi_b\) 原文未说明。与式 (3.14) 的取舍见该式 Q5。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.16（气泡群修正 \(\upsilon_b\)）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.16。
❓ **Q2 变量与单位**：B 类颗粒，Hilligardt（1986）：

\[
\upsilon_b=\begin{cases}0{,}64&d_B<0{,}1\,\mathrm{m}\\2{,}0\cdot d_B^{0{,}5}&0{,}1\le d_B\le 1{,}0\,\mathrm{m}\\2{,}0&d_B>1{,}0\,\mathrm{m}\end{cases}
\]

\(d_B\) 以 m 计时，中段 \(2{,}0\cdot d_B^{0{,}5}\) 的系数单位为 m\(^{-0{,}5}\)，原文未写出该单位。两端为无量纲。
❓ **Q3 德文原句**：„Für Gruppe-B-Partikel gilt nach Hilligardt (1986)“。
❓ **Q4 上游依赖**：用于式 (3.13)。\(d_B\) 与后文 \(d_b\) 是否同一直径，原文未声明。
❓ **Q5 原文空白**：符号是 \(d_B\) 而不是 \(d_b\)。可选解释：(a) 同一气泡直径的另一种印刷；(b) 特指某一个气泡直径定义。待总司令裁定。A、D 类未给。\(d_B=0{,}1\,\mathrm{m}\) 落在中段（含等号），下段是严格小于。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.17（单气泡上升速度）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.17。
❓ **Q2 变量与单位**：

\[
u_{b,i}=\sqrt{\frac{1}{c_b}\cdot\frac{4}{3}\cdot d_b\cdot g}
\]

\(c_b\) 单气泡阻力系数，见式 (3.18)。\(d_b\) 气泡直径。
❓ **Q3 德文原句**：„Die Steiggeschwindigkeit einer isolierten Einzelblase“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.17，S. 26）。
❓ **Q5 原文空白**：根号内是 \((1/c_b)\cdot(4/3)\cdot d_b\cdot g\)，原文未再写成 \(\sqrt{g d_b}\) 的倍数。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.18（气泡阻力系数）

❓ **Q1 出处**：PDF 23 左半页 / S. 26，§3.1.2，Eq. 3.18。
❓ **Q2 变量与单位**：

\[
c_b=\frac{16}{Re_b}+2{,}64
\]

原文：对整个雷诺数范围有效。\(Re_b\) 见式 (3.19)。
❓ **Q3 德文原句**：„mit der für den gesamten Reynolds-Zahlenbereich gültigen Gleichung“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（同页）。
❓ **Q5 原文空白**：16 与 2,64 的来源页只指向 Hilligardt，没有再推导。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.19（气泡雷诺数）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.19。
❓ **Q2 变量与单位**：

\[
Re_b=\frac{u_b\cdot d_b}{\nu_d}
\]

分母是悬浮相运动粘度 \(\nu_d\)，不是 \(\nu_g\)。无量纲。
❓ **Q3 德文原句**：„Die Reynolds-Zahl der Blase ist definiert als“。
❓ **Q4 上游依赖**：\(\nu_d\) 与 \(\mu_d\) 的关系见式 (3.20)。
❓ **Q5 原文空白**：无。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.20（悬浮相两种粘度）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.20。
❓ **Q2 变量与单位**：

\[
\mu_d=\nu_d\cdot\rho_g\cdot(1-\varepsilon_d)
\]

\(\mu_d\) 动力粘度，\(\nu_d\) 运动粘度，\(\varepsilon_d\) 悬浮相空隙率。由式 (3.21) 的单位可把 \(\mu_d\) 读成 kg/(m·s)，即 Pa·s；原文没有把 Pa·s 印出来。
❓ **Q3 德文原句**：„Die dynamische und die kinematische Viskosität der Suspensionsphase sind dabei durch die Beziehung … miteinander verknüpft (Hilligardt 1986)“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.20，S. 27）。\(\mu_d\) 由式 (3.21) 或 (3.22) 计算，再由本式得 \(\nu_d\)。这一计算顺序是根据“先计算 \(\mu_d\)”的下一句读出的。
❓ **Q5 原文空白**：密度用的是 \(\rho_g\)，再乘 \((1-\varepsilon_d)\)。原文没有写成颗粒密度。不另改写。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.21（\(d_P<0{,}2\,\mathrm{mm}\) 时的 \(\mu_d\)）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.21。
❓ **Q2 变量与单位**：

\[
\mu_d=f_d\cdot d_P^{2/3},\qquad f_d=232\,\frac{\mathrm{kg}}{\mathrm{s}\cdot\mathrm{m}^{5/3}}
\]

条件：\(d_P<0{,}2\,\mathrm{mm}\)。\(f_d\) 的分母是 \(\mathrm{m}^{5/3}\)。\(d_P^{2/3}\) 与该单位相乘得到 kg/(m·s)。
❓ **Q3 德文原句**：„Für \(d_P<0{,}2\,\mathrm{mm}\)“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.21，S. 27）。
❓ **Q5 原文空白**：判据用 mm，\(f_d\) 用 m。代入前是否把 \(d_P\) 换成 m，原文未写。可选解释：(a) \(d_P\) 以 m 代入（影响：与 \(f_d\) 的米制单位一致）；(b) \(d_P\) 以 mm 代入（影响：\(\mu_d\) 差 \(10^{-2}\) 量级）。待总司令裁定。\(d_P=0{,}2\,\mathrm{mm}\) 不属于本式（严格小于）。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.22（\(d_P\ge 0{,}2\,\mathrm{mm}\) 时的 \(\mu_d\)）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.22。
❓ **Q2 变量与单位**：

\[
\mu_d=f_d\cdot d_P^{3/2},\qquad f_d=2{,}5\cdot 10^{5}\,\frac{\mathrm{kg}}{\mathrm{s}\cdot\mathrm{m}^{5/2}}
\]

指数是 \(3/2\)，分母是 \(\mathrm{m}^{5/2}\)。相乘同样得到 kg/(m·s)。
❓ **Q3 德文原句**：„Für \(d_P\ge 0{,}2\,\mathrm{mm}\)“。
❓ **Q4 上游依赖**：同式 (3.21)。
❓ **Q5 原文空白**：同式 (3.21) 的 mm 与 m。等号落在本式。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.23（穿流速度）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.23。
❓ **Q2 变量与单位**：

\[
u_{b,r}=n_b\cdot u_d,\qquad n_b=2{,}7
\]

❓ **Q3 德文原句**：„nach Messungen von Hilligardt (1986)“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.23，S. 27）。加压形式为式 (3.44)。
❓ **Q5 原文空白**：2,7 未给不确定度。\(u_d\) 来自式 (3.11) 或 (3.12) 的哪一个，随那两式的取舍而定。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.24（气泡体积份额）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.24。
❓ **Q2 变量与单位**：

\[
\varepsilon_b=\frac{u_0-u_d}{u_b+n_b u_d-u_d}
\]

无量纲。
❓ **Q3 德文原句**：„nach dem sichtbaren Blasenvolumenanteil auflösen“。
❓ **Q4 上游依赖**：由式 (3.8)–(3.10) 与 (3.23) 代入式 (3.7) 得到。同页写明这一推导。
❓ **Q5 原文空白**：分母为零或 \(\varepsilon_b\) 超出 \([0,1]\) 时如何处理，原文未说明。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.25（悬浮相空隙率）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.25。
❓ **Q2 变量与单位**：

\[
\varepsilon_d=\varepsilon_{mf}\cdot\left(\frac{u_d}{u_{mf}}\right)^{1/n_{RZ}}
\]

指数是 \(1/n_{RZ}\)，不是 \(n_{RZ}\)。
❓ **Q3 德文原句**：„Hilligardt (1986) verwendet den Ansatz von Richardson und Zaki (1954)“。
❓ **Q4 上游依赖**：Hamel 引用的 Richardson und Zaki（1954），经 Hilligardt（1986）使用（式 3.25，S. 27）。\(n_{RZ}\) 见式 (3.26) 或 (3.27)。同页引 Kunii und Levenspiel（1991）说明悬浮相空隙率可以高于最小流化空隙率。
❓ **Q5 原文空白**：\(n_{RZ}\) 用式 (3.26) 还是式 (3.27)，原文在介绍 Rowe 时没有写“本文只采用其中之一”。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.26（Richardson–Zaki 指数）

❓ **Q1 出处**：PDF 23 右半页 / S. 27，§3.1.2，Eq. 3.26。
❓ **Q2 变量与单位**：

\[
n_{RZ}=\begin{cases}
4{,}65 & Re_S\le 0{,}2\\
4{,}4\cdot Re_S^{-0{,}03} & 0{,}2<Re_S\le 1{,}0\\
4{,}4\cdot Re_S^{-0{,}1} & 1{,}0<Re_S\le 500\\
2{,}4 & Re_S>500
\end{cases}
\]

指数是 \(-0{,}03\) 与 \(-0{,}1\)。\(Re_S\) 见式 (3.28)。
❓ **Q3 德文原句**：式前为 „Danach gilt“，承 Richardson–Zaki 句。
❓ **Q4 上游依赖**：同式 (3.25)。
❓ **Q5 原文空白**：\(Re_S=0{,}2\) 落在第一段（含等号）。原文未说明 \(Re_S\le 0\)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.27（Rowe 的指数）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.27。式属于上一页开始的 Rowe 句。
❓ **Q2 变量与单位**：

\[
n_{RZ}=\frac{2{,}35\cdot\left(2+0{,}175\cdot Re_S^{0{,}75}\right)}{1+0{,}175\cdot Re_S^{0{,}75}}
\]

❓ **Q3 德文原句**：„Wein (1992) schlägt die Verwendung des Ansatzes von Rowe (1987) … vor“（句在 S. 27 底部，式在 S. 28）。
❓ **Q4 上游依赖**：Hamel 引用的 Rowe（1987），由 Wein（1992）建议使用（S. 27–28）。
❓ **Q5 原文空白**：这是 Wein 的建议，原文没有写“本文改用 Rowe、停用式 (3.26)”。两式并存。待总司令裁定是否只保留一个。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.28（单颗粒沉降雷诺数）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.28。
❓ **Q2 变量与单位**：

\[
Re=\frac{u_S\cdot d_P\cdot\rho_g}{\mu_g}
\]

印刷符号是 \(Re\)，不是 \(Re_S\)。同页写两种情形下的雷诺数都用单颗粒沉降速度构成。\(\mu_g\) 是气体动力粘度，不是式 (3.20) 的 \(\mu_d\)。
❓ **Q3 德文原句**：„die Reynolds-Zahl mit der Einzelkornsinkgeschwindigkeit gebildet wird“。
❓ **Q4 上游依赖**：\(u_S\) 见式 (3.29)、(3.31)。式 (3.26)、(3.27) 写成 \(Re_S\)。
❓ **Q5 原文空白**：式 (3.28) 的符号是 \(Re\)，前后文是 \(Re_S\)。可选解释：(a) 同一量，符号印得不一致；(b) 不是同一个雷诺数。待总司令裁定。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.29（单颗粒沉降速度）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.29。
❓ **Q2 变量与单位**：

\[
u_S=\sqrt{\frac{4\cdot g\cdot d_P\cdot(\rho_P-\rho_g)}{3\cdot\rho_g\cdot c_w}}
\]

分母密度是 \(\rho_g\)。\(c_w\) 见式 (3.30)。
❓ **Q3 德文原句**：„Lösung der Kräftebilanz für die umströmte Einzelkugel“。
❓ **Q4 上游依赖**：力平衡，本页没有另引文献作者作为本式的来源。
❓ **Q5 原文空白**：颗粒被写成单球。非球形是否只通过 \(c_w\) 进入，原文未说明。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.30（Martin 阻力系数）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.30。
❓ **Q2 变量与单位**：

\[
c_w=\frac13\left(\sqrt{\frac{72}{Re_S}}+1\right)^2
\]

❓ **Q3 德文原句**：„wird zur Berechnung des Widerstandsbeiwertes die Korrelation von Martin (1984) verwendet“。
❓ **Q4 上游依赖**：Hamel 引用的 Martin（1984）（式 3.30，S. 28）。同页写 Bürkle（1991）、Dersch（1994）、Heinbockel（1995）也如此用。
❓ **Q5 原文空白**：式中写 \(Re_S\)，式 (3.28) 定义的是 \(Re\)。见式 (3.28) 的 Q5。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.31（解出的沉降速度）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.31。
❓ **Q2 变量与单位**：

\[
u_S=18\cdot\left(\frac{\nu_g}{d_P}\right)\cdot\left(\sqrt{1+\frac19\cdot\sqrt{Ar}}-1\right)^2
\]

结构是：先对 \(Ar\) 开方，乘 \(1/9\)，加 1，再开方，减 1，然后整个括号平方。外面不再套一层只盖住“减 1”的根号。
❓ **Q3 德文原句**：„kann Gleichung (3.29) nach der gesuchten Partikelsinkgeschwindigkeit aufgelöst werden“。
❓ **Q4 上游依赖**：由式 (3.30) 与 (3.28) 解式 (3.29)。\(Ar\) 用式 (3.32)。
❓ **Q5 原文空白**：原文未给出这一代数变形的中间步骤。18 与 \(1/9\) 以原页为准。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.32（沉降用的 Ar）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.32。
❓ **Q2 变量与单位**：

\[
Ar=\frac{g\cdot d_P^3\cdot(\rho_P-\rho_g)}{\nu_g^2\cdot\rho_g}
\]

粒径无平均号，与式 (3.2) 的 \(\bar d_P\) 不同。
❓ **Q3 德文原句**：„mit der Archimedes-Zahl“。
❓ **Q4 上游依赖**：供式 (3.31)。
❓ **Q5 原文空白**：与式 (3.2) 是否必须用同一个粒径定义，原文未说明。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.33（气泡长大）

❓ **Q1 出处**：PDF 24 左半页 / S. 28，§3.1.2，Eq. 3.33。
❓ **Q2 变量与单位**：

\[
\frac{dd_b}{dh}=\left[\left(\frac{2}{9\pi}\right)^{1/3}\cdot\frac{\varepsilon_b^{1/3}}{1-\xi_b\left(\dfrac{6}{\pi}\right)^{1/3}\varepsilon_b^{1/3}}-\frac{d_b}{3\lambda_b u_b}\right]
\]

\(h\) 为高度。左端无量纲（长度/长度）。S. 29：第一项是聚并引起的直径增加，第二项是分裂引起的减小。\(\xi_b\) 见式 (3.34)，\(\lambda_b\) 见式 (3.35)。
❓ **Q3 德文原句**：„Eine … Berechnungsleichung für das Blasenwachstum wurde von Hilligardt (1986) vorgestellt“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.33，S. 28）。加压形式为式 (3.41)。
❓ **Q5 原文空白**：微分方程的积分方法和床高方向的离散，原文未在本式处给出。初始直径与起始高度见式 (3.36)–(3.40)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.34（\(\xi_b\)）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.34。
❓ **Q2 变量与单位**：

\[
\xi_b=\begin{cases}1-\alpha_b^3&0\le\alpha_b\le 1\\0&\alpha_b>1\end{cases}
\]

指数是 3。
❓ **Q3 德文原句**：„Für \(\xi_b\) gilt nach Hilligardt (1986)“。
❓ **Q4 上游依赖**：\(\alpha_b\) 为式 (3.6)。Hamel 引用的 Hilligardt（1986）（式 3.34，S. 29）。
❓ **Q5 原文空白**：\(\alpha_b<0\) 未写。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.35（气泡平均寿命）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.35。
❓ **Q2 变量与单位**：

\[
\lambda_b=280\,\frac{u_{mf}}{g}
\]

原文称为 mittlere Lebensdauer，故 \(\lambda_b\) 是时间。\(u_{mf}/g\) 的 SI 为 s。280 未标单位，按此式为无量纲。原文写 „näherungsweise“。
❓ **Q3 德文原句**：„Für die mittlere Lebensdauer einer Blase gilt nach Hilligardt (1986) näherungsweise“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.35，S. 29）。加压形式为式 (3.42)。
❓ **Q5 原文空白**：280 的有效数字只有印刷出来的三位。适用的颗粒与床尺寸原文未在本式旁重复。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.36（多孔帽起始气泡直径）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.36。
❓ **Q2 变量与单位**：

\[
d_{b,0}=1{,}3\cdot\left(\frac{\dot V_0^2}{g}\right)^{0{,}2}
\]

工业布风器。\(d_{b,0}\) 为长度。1,3 的单位原文未单列；\(\dot V_0^2/g\) 的 0,2 次方为长度时，1,3 无量纲。
❓ **Q3 德文原句**：„Für technische Gasverteiler verwendet Hilligardt (1986) die Beziehung von Davidson und Schüler (1960)“。
❓ **Q4 上游依赖**：Hamel 引用的 Davidson und Schüler（1960），经 Hilligardt（1986）使用（式 3.36，S. 29）。\(\dot V_0\) 见式 (3.37)。
❓ **Q5 原文空白**：多孔板不用本式，见式 (3.38)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.37（单孔体积流量）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.37。
❓ **Q2 变量与单位**：\(\dot V_0=u_0\cdot A/n_0\)。\(n_0\) 为布风板上出口孔的个数，无量纲。
❓ **Q3 德文原句**：„Dabei ist \(n_0\) die Anzahl der Austrittsöffnungen im Düsenboden“。
❓ **Q4 上游依赖**：供式 (3.36)。
❓ **Q5 原文空白**：孔不完全均分流量时如何处理，原文未说明。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.38（多孔板起始气泡直径）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.38。
❓ **Q2 变量与单位**：

\[
d_{b,0}=\frac{2{,}78}{g}\cdot(u_0-u_{mf})^2
\]

\((u)^2/g\) 为长度，故 2,78 按此式无量纲。
❓ **Q3 德文原句**：„die Blasenstartgröße nach der von Kunii und Levenspiel (1991) angegebenen Beziehung“。
❓ **Q4 上游依赖**：Hamel 引用的 Kunii und Levenspiel（1991）（式 3.38，S. 29）。用于多孔板。
❓ **Q5 原文空白**：与式 (3.36) 的分界是布风器类型，原文如此写。混合布风器未说明。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.39（开始长大的高度，有出口孔）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.39。
❓ **Q2 变量与单位**：

\[
h_0=h_A+\frac{d_{b,0}}{2}
\]

\(h_A\)：出口孔高出布风板的高度。长度。
❓ **Q3 德文原句**：„Die Höhe \(h_0\), bei der das Blasenwachstum beginnt“。
❓ **Q4 上游依赖**：Hamel 引用的 Hilligardt（1986）（式 3.39，S. 29）。
❓ **Q5 原文空白**：\(h_A\) 的数值原文未给，作为装置输入。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.40（多孔板的起始高度）

❓ **Q1 出处**：PDF 24 右半页 / S. 29，§3.1.2，Eq. 3.40。索引 `locator: page_image`。
❓ **Q2 变量与单位**：\(h_0=d_{b,0}/2\)。多孔板。长度。
❓ **Q3 德文原句**：„dass das Blasenwachstum in der Höhe \(h_0=d_{b,0}/2\) beginnt“。
❓ **Q4 上游依赖**：同页，Hilligardt（1986）的处理。
❓ **Q5 原文空白**：无。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.41（加压后的气泡长大）

❓ **Q1 出处**：PDF 25 右半页 / S. 31，§3.1.3，Eq. 3.41。
❓ **Q2 变量与单位**：

\[
\frac{dd_b}{dh}=\left[\left(\frac{2}{9\pi}\right)^{1/3}\cdot\frac{\varepsilon_b^{1/3}\cdot\left(\dfrac{P}{P_0}\right)^{P_0/P}}{1-\xi_b\cdot\left(\dfrac{6}{\pi}\right)^{1/3}\cdot\varepsilon_b^{1/3}}-\frac{d_b}{3\cdot\lambda_b\cdot u_b}\right]
\]

\((P/P_0)\) 的指数读作 \(P_0/P\)，不是一个纯数字。\(P=P_0\) 时该因子为 1，式 (3.33) 的第一项结构得以保留。指数较小，列入文末待肉眼确认。
❓ **Q3 德文原句**：„hat Heinbockel (1995) die Gleichung (3.33) folgendermaßen abgeändert“。
❓ **Q4 上游依赖**：Hamel 引用的 Heinbockel（1995）（式 3.41，S. 31），用来计入随压力增强的气泡聚并。实验比较对象为 Olowson und Almstedt（1990）。同页写：所加的幂次在大气压力下应回到未修改的模型。
❓ **Q5 原文空白**：指数 \(P_0/P\) 的字形待肉眼确认。若确认如此，它与“大气压下因子为 1”一致。\(P_0\) 的数值（哪一个大气压）原文未给。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.42（加压后的气泡寿命）

❓ **Q1 出处**：PDF 25 右半页 / S. 31，§3.1.3，Eq. 3.42。
❓ **Q2 变量与单位**：

\[
\lambda_b=280\cdot\frac{u_{mf}}{g}\cdot\left(\frac{P}{P_0}\right)^{-0{,}7}
\]

指数是 \(-0{,}7\)。\(P=P_0\) 时回到式 (3.35)。
❓ **Q3 德文原句**：„Die Gleichung (3.35) lautet in modifizierter Form“。
❓ **Q4 上游依赖**：Hamel 引用的 Heinbockel（1995）（式 3.42，S. 31）。指数来源同页指向 Olowson und Almstedt（1990）的穿流测量，并说式 (3.42)–(3.44) 的指数由此确定（这句话在 S. 32）。
❓ **Q5 原文空白**：\(-0{,}7\) 未给拟合范围。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.43（加压后的气泡上升速度）

❓ **Q1 出处**：PDF 25 右半页 / S. 31，§3.1.3，Eq. 3.43。
❓ **Q2 变量与单位**：

\[
u_b=\psi_b\cdot(u_0-u_{mf})\cdot\left(\left(\frac{P}{P_0}\right)^{0{,}2}-1\right)+\upsilon_b\cdot u_{b,i}
\]

按页面图像，括号内是 \((P/P_0)^{0{,}2}-1\)。指数 0,2。
❓ **Q3 德文原句**：„durch Einführen eines mit dem Druck ansteigenden Koeffizienten vor der Überschussgeschwindigkeit in Gleichung (3.13)“。
❓ **Q4 上游依赖**：Hamel 引用的 Heinbockel（1995）（式 3.43，S. 31）。测量现象来自 Olowson und Almstedt（1990）：气泡变小的同时上升速度上升。
❓ **Q5 原文空白**：\(P=P_0\) 时第一项为零，不能回到式 (3.13)。同页又写幂次修正应保证大气压力下仍是未修改的模型。两处原文并存。可选解释：(a) 式 (3.43) 的 “\(-1\)” 按印刷保留，大气压下过剩气速项消失（影响：常压气泡速度只剩 \(\upsilon_b u_{b,i}\)）；(b) “\(-1\)” 是误印，系数应为 \((P/P_0)^{0{,}2}\)（影响：常压回到式 (3.13)，且系数随压力上升）。待总司令裁定。请肉眼确认 PDF 25 右半页、式 (3.43) 在 \((P/P_0)^{0{,}2}\) 之后是否有 \(-1\)。读图结果是有。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.44（加压后的穿流速度）

❓ **Q1 出处**：PDF 26 左半页 / S. 32，§3.1.3，Eq. 3.44。
❓ **Q2 变量与单位**：

\[
u_{b,r}=n_b\cdot u_d\cdot\left(\frac{P}{P_0}\right)^{-0{,}15}
\]

指数是 \(-0{,}15\)。\(P=P_0\) 时回到式 (3.23)。
❓ **Q3 德文原句**：„durch Modifikation von Gleichung (3.23)“。
❓ **Q4 上游依赖**：Hamel 引用的 Heinbockel（1995）（式 3.44，S. 32）。穿流随压力下降的测量：Olowson und Almstedt（1990，1992）。Bild 3.3（同页）比较的是按式 (3.44) 的计算结果，图中条件为 \(T=25^\circ\mathrm{C}\)，\(d_P=0{,}7\,\mathrm{mm}\)，\(\rho_P=2600\,\mathrm{kg/m^3}\)，\(h=0{,}55\,\mathrm{m}\)。这些是图上的比较条件，不是本式的适用范围声明。
❓ **Q5 原文空白**：图注条件是否就是式 (3.44) 的适用界限，原文未写。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

Bild 3.4、Bild 3.5 在 PDF 26 右半页 / S. 33，是式 (3.14)–(3.44) 的计算结果与 Olowson und Almstedt（1990）测量的比较图，没有新的编号方程。

## Eq. 3.45（相间交换的定义）

❓ **Q1 出处**：PDF 27 左半页 / S. 34，§3.1.4，Eq. 3.45。
❓ **Q2 变量与单位**：

\[
-\frac{1}{V_b}\cdot\frac{dN_{i,ex}}{dt}=K_{bd}\cdot(C_{i,b}-C_{i,d})
\]

\(V_b\) 气泡体积。\(N_{i,ex}\) 组分 \(i\) 的交换量。\(C_{i,b}\)、\(C_{i,d}\) 为气泡相与密相（悬浮相）中组分 \(i\) 的浓度。\(K_{bd}\) 在 S. 35 以 s\(^{-1}\) 给出。\(N\) 与 \(C\) 是摩尔还是质量，原文未写。左端有负号，浓度差是气泡减去密相。
❓ **Q3 德文原句**：„Der Austausch der Gaskomponente i zwischen Blase und Suspension, bezogen auf das Blasenvolumen, ist definiert durch“。
❓ **Q4 上游依赖**：同页用 \(\alpha_b\) 区分慢气泡与快气泡，并引 Davidson und Harrison（1963）。
❓ **Q5 原文空白**：负号与浓度差同时出现。可选解释：(a) \(N_{i,ex}\) 的增加表示组分进入气泡，因而浓度差要写成密相减气泡，印刷的 \((C_{i,b}-C_{i,d})\) 与负号一起相当于密相减气泡；(b) \(N_{i,ex}\) 的正方向相反。原文没有定义 \(N_{i,ex}\) 的正方向。待总司令裁定。与式 (2.4) 中 \(\dot N_{ex,bd}\) 的对应关系，本章没有写出。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.46（Kunii–Levenspiel 的串联阻力）

❓ **Q1 出处**：PDF 27 左半页 / S. 34，§3.1.4，Eq. 3.46。
❓ **Q2 变量与单位**：

\[
K_{bd}=\frac{1}{\dfrac{1}{K_{bc}}+\dfrac{1}{K_{cd}}}
\]

单位与 \(K_{bc}\)、\(K_{cd}\) 相同。S. 35 把 \(K_{bd}\) 用在 s\(^{-1}\)。
❓ **Q3 德文原句**：„das von Kunii und Levenspiel (1977) vorgeschlagene Modell“。
❓ **Q4 上游依赖**：Hamel 引用的 Kunii und Levenspiel（1977）（式 3.46，S. 34）。本页介绍该模型常用于快气泡，没有写“本文采用式 (3.46) 作为唯一交换式”。本文写出的计算式是式 (3.50)。
❓ **Q5 原文空白**：见式 (3.50)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.47（\(K_{bc}\) 与 \(K_{cd}\)）

❓ **Q1 出处**：PDF 27 左半页 / S. 34，§3.1.4，Eq. 3.47。
❓ **Q2 变量与单位**：

\[
K_{bc}=4{,}5\cdot\frac{u_{mf}}{d_b}+5{,}85\cdot\frac{D_g^{0{,}5}\cdot g^{0{,}25}}{d_b^{1{,}25}},\qquad K_{cd}=6{,}78\cdot\sqrt{\frac{D_g\cdot u_b\cdot\varepsilon_{mf}}{d_b^3}}
\]

\(D_g\) 气体扩散系数。各项在 \(d_b\) 以 m、\(u\) 以 m/s、\(D_g\) 以 m\(^2\)/s、\(g\) 以 m/s\(^2\) 时为 s\(^{-1}\)。原文未逐项写出这些单位。
❓ **Q3 德文原句**：式 (3.46) 的 „mit“。
❓ **Q4 上游依赖**：同属 Kunii und Levenspiel（1977）。同页：\(K_{bc}\) 的第一项和式 (3.48) 的第一项是穿流引起的对流交换，其余项是扩散。
❓ **Q5 原文空白**：5,85 与 6,78 的有效数字以原页为准。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.48（Sit 与 Grace）

❓ **Q1 出处**：PDF 27 左半页 / S. 34，§3.1.4，Eq. 3.48。
❓ **Q2 变量与单位**：

\[
K_{bd}=\frac{2\cdot u_{mf}}{d_b}+\sqrt{\frac{144\cdot D_g\cdot\varepsilon_{mf}\cdot u_b}{\pi\cdot d_b^3}}
\]

根号内系数是 144，分母有 \(\pi\)。用于聚并气泡。
❓ **Q3 德文原句**：„Sit und Grace (1981) schlagen einen ähnlichen Ansatz für koaleszierende Blasen vor“。
❓ **Q4 上游依赖**：Hamel 引用的 Sit und Grace（1981）（式 3.48，S. 34）。式 (3.50) 的扩散项与本式第二项相同。
❓ **Q5 原文空白**：本式与式 (3.50) 的对流项不同。原文没有说用本式代替式 (3.50)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.49（S. 35，Preto 的对流项）

必须带页码。S. 36 的式 (3.49) 是另一个方程，见文末。

❓ **Q1 出处**：PDF 27 右半页 / S. 35，§3.1.4，Eq. 3.49。索引 `ambiguous_number: true`。
❓ **Q2 变量与单位**：

\[
K_{bd}=4{,}5\cdot\chi\cdot\frac{u_{mf}}{d_b}
\]

\(\chi=1\)：慢的、被穿流的气泡。快气泡的 \(\chi<1\)，Preto（1986）给出的最小界限是 \(\chi=0{,}1\)。
❓ **Q3 德文原句**：„Preto (1986), zitiert in Verweyen (1994), schlägt … folgenden Ansatz vor“（句在 S. 34 末，式在 S. 35）。
❓ **Q4 上游依赖**：Hamel 引用的 Preto（1986），转引自 Verweyen（1994）（式 3.49，S. 35）。
❓ **Q5 原文空白**：快气泡的 \(\chi\) 只给了下限 0,1，没有给本文使用的函数。S. 35 写明：本文的对流交换按穿流计算，不把本式当作本文的最终交换式。见式 (3.50)。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## Eq. 3.50（S. 35，本文的相间交换系数）

必须带页码。S. 36 的式 (3.50) 是另一个方程，见文末。

❓ **Q1 出处**：PDF 27 右半页 / S. 35，§3.1.4，Eq. 3.50。索引 `ambiguous_number: true`。
❓ **Q2 变量与单位**：

\[
K_{bd}=\frac{\dot V_{b,r}}{V_b}+\sqrt{\frac{144\cdot D_g\cdot\varepsilon_{mf}\cdot u_b}{\pi\cdot d_b^3}}=\frac{3}{2}\cdot\frac{u_{b,r}}{d_b}+\sqrt{\frac{144\cdot D_g\cdot\varepsilon_{mf}\cdot u_b}{\pi\cdot d_b^3}}
\]

单位 s\(^{-1}\)。扩散项与式 (3.48) 的第二项相同。对流项写成穿流体积流量除以气泡体积，并等于 \(\frac32 u_{b,r}/d_b\)。
❓ **Q3 德文原句**：„In dieser Arbeit wird der konvektive Gasaustausch zwischen Blase und Suspension anhand der in Kapitel 4.1.3 beschriebenen Berechnung der Blasendurchströmung ermittelt.“
❓ **Q4 上游依赖**：对流部分指向原文印刷的 „Kapitel 4.1.3“ 的穿流计算。扩散部分：Hamel 引用的 Sit und Grace（1981）（式 3.50，S. 35）。\(u_{b,r}\) 在本章是式 (3.23)、(3.44)。\(\dot V_{b,r}\) 是式 (3.9)。
❓ **Q5 原文空白**：

- 印刷文字是 Kapitel 4.1.3，不是 §3.1。Bild 2.2 也把流动力学标成 Kap. 4。本印刷本目录里，穿流方程在 §3.1（式 3.23、3.44）。可选解释：(a) 正文沿用另一套章号，4.1.3 指本章的穿流计算（影响：实现时用式 3.23/3.44）；(b) 4.1.3 指目录中的另一处，而该处不在 §3.1（影响：不能把式 3.23 直接当作此处的依据）。待总司令裁定。
- 由 \(\dot V_{b,r}/V_b\) 得到 \(\frac32 u_{b,r}/d_b\) 时，原文没有印出气泡体积和穿流面积的几何假定。可选解释：(a) 球体积 \(\pi d_b^3/6\) 与截面积 \(\pi d_b^2/4\) 相除得到 3/2（影响：气泡必须按球）；(b) 3/2 是独立系数，不依赖该几何（影响：不能用别的气泡形状去改 3/2）。待总司令裁定。
- 同页：Deinert 等（1986）报告高温下 \(K_{bd}>40\,\mathrm{s^{-1}}\)，式 (3.46)–(3.49) 算不出。把式 (3.50) 与 Kapitel 4.1 的流动计算外推到床温约 \(900^\circ\mathrm{C}\)，小气泡得到 \(K_{bd}\approx 33\,\mathrm{s^{-1}}\)。这是原文给出的数量级，不是本式的拟合常数。Bild 3.6 的比较条件（同页）：\(T_{WS}=298{,}15\,\mathrm{K}\)，\(d_P=0{,}7\,\mathrm{mm}\)，\(\rho_P=2600\,\mathrm{kg/m^3}\)，\(\varepsilon_{mf}=0{,}46\)，\(\phi_S=0{,}8\)，\(u_0-u_{mf}=0{,}6\,\mathrm{m/s}\)，\(d_b=0{,}3\,\mathrm{m}\)。这些是图的条件。
❓ **Q6 与归档旧代码的差异**：比对阶段填写。

## 式号区分：S. 36 上的 (3.49) 与 (3.50)

这两式在 §3.2.1，不在本卡的展开范围内。只记录它们与 S. 35 不是同一方程。

- PDF 28 左半页 / S. 36，§3.2.1，Eq. 3.49：\(\delta_e=\bar d_P\)。原文：气泡破裂时，气泡顶盖的固体层厚度等于床料平均粒径。Hamel 引用的 Pemberton und Davidson（1986），并写与 Pyle（1965）一致。
- PDF 28 左半页 / S. 36，§3.2.1，Eq. 3.50：\(F_{0,r}=\bar d_P\cdot\rho_P\cdot(1-\varepsilon_{mf})\cdot\pi\cdot d_b^2/2\)。一次气泡破裂由顶盖机制抛出的固体量。

S. 35 的 (3.49)、(3.50) 是相间交换系数。S. 36 的 (3.49)、(3.50) 是自由空间抛出。引用时写 Hamel Eq. 3.49, S. 35 或 S. 36。

## 待目标所有者肉眼判读

1. PDF 25 右半页，S. 31，式 (3.41)，\((P/P_0)\) 的指数。读图结果是 \(P_0/P\)。
2. PDF 25 右半页，S. 31，式 (3.43)，\((P/P_0)^{0{,}2}\) 之后是否有 \(-1\)。读图结果是有。
3. PDF 23 左半页，S. 26，式 (3.13)，\(u_{b,i}\) 前的系数。读图结果是希腊字母 \(\upsilon\)，不是拉丁字母 \(v\)。
