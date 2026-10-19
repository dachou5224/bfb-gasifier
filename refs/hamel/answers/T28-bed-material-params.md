# T28 原文答卷：床料的 ε_mf、φ_S、u_mf,0 与 Geldart 分界

| 项 | 内容 |
|---|---|
| 任务卡 | `docs/tasks/T28-bed-material-params.md` |
| 执行者 | cursor |
| 会话 | 新开的 cursor 会话 |
| 工作树 | `../bfb-T28`，分支 `task/T28-bed-material-params`，基点 `origin/main`（`ce020ae`） |
| 日期 | 2026-10-18；同日补第 2 项 |
| 审查方 | 总司令 |
| 本次范围 | 任务卡第 1、2、3 项。A/C 分界只有图，未数字化 |

只查证。Hamel 部分以页面图像为准，逐页文本只用于定位。Geldart (1973) 用主工作树本地文件 `refs/external/geldart1973-types-of-gas-fluidization.pdf`（不入库）；该 PDF 的文字层是 OCR，数值与式子以页面图像为准。不选定 \(\varepsilon_{mf}\)、\(\varphi_S\)，不计算 Geldart 类别，不写 `hamel/` 代码。未打开归档代码、`docs/hamel_submodels/` 和旧 parity 读数。第 6 问一律写“未查看归档”。

德文排版：点是千分位，逗号是小数点。

对 VTT 与 HTW 的褐煤焦炭（或焦炭加灰）床料，Hamel 没有给出 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf,0}\)、床料密度和床料粒径。这五项对第 7 章工况都是原文未给出。

书中另有数值，所指都不是这些工况的床料。下表只作检索索引，依据在后文各节。

| 印刷位置 | 量 | 原文数值 | 所指 |
|---|---|---|---|
| Bild 3.6，S. 35 | \(\varepsilon_{mf}\)、\(\varphi_S\) | 0,46；0,8 | 该图的计算条件。物料未命名。\(T_{WS}=298{,}15\,\mathrm{K}\)，\(d_P=0{,}7\,\mathrm{mm}\)，\(\rho_P=2600\,\mathrm{kg/m^3}\)，\(u_0-u_{mf}=0{,}6\,\mathrm{m/s}\) |
| S. 45，§3.2.4 | \(\bar\varphi_S\) | 0,8 | 自由空间轨迹；无球形度数据时的程序输入。分布区间 \(0{,}3\le\varphi_S\le 1{,}0\) |
| S. 44–45 | \(\varphi_S\) | 1 到 0,026 | Haider und Levenspiel（1989）阻力关联的拟合范围 |
| 式 (3.15)，S. 26 | \(u_{mf}\) | 区间 \(0{,}014\)–\(0{,}18\,\mathrm{m/s}\) | Wein（1992）的 \(\psi_b\) 关联的适用区间，不是某床的测量值 |
| Bild 3.4 与邻图，S. 32–33 | \(d_P\)、\(\rho_P\) | \(0{,}7\,\mathrm{mm}\)；\(2600\,\mathrm{kg/m^3}\) | Olowson und Almstedt（1990）的测量条件。\(T=25^\circ\mathrm{C}\)，\(h=0{,}55\,\mathrm{m}\)。图上的 \(0{,}1\)–\(0{,}6\,\mathrm{m/s}\) 是 \(u_0-u_{mf}\) |
| Bild 3.7，S. 42 | \(d_P\)、\(\rho_P\) | \(1{,}0\,\mathrm{mm}\)；\(2200\,\mathrm{kg/m^3}\) | 该图的计算条件。\(T_{WS}=1073\,\mathrm{K}\)，\(p=0{,}1\,\mathrm{MPa}\)。无 \(\varepsilon_{mf}\)、\(\varphi_S\) |
| 表 5.1，S. 79 | \(\rho_{inert}\) | \(2000\,\mathrm{kg/m^3}\) | 第 5 章动力学比较算例的惰性料密度。燃料行写的是焦炭 |
| 表 7.5，S. 147 | \(u_0\) | \(0{,}75\,\mathrm{m/s}\) | 参数研究标准工况的空塔速度，床高 \(1{,}2\,\mathrm{m}\)。不是 \(u_{mf}\) |
| S. 25 | 过剩气速 | \(2{,}87\) 与 \(2{,}97\,\mathrm{m/s}\) | Wein 所用的砂与灰。空塔速度至 \(3\,\mathrm{m/s}\)。这两种物料的 \(u_{mf}\) 原文未给出 |
| 表 7.1–7.5 | \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)、床料密度、床料粒径 | 原文未给出 | HTW、VTT 加压、VTT 常压、WSV400、参数研究标准工况 |

---

## 1. 符号表中的定义（无数值）

❓ **Q1 出处**：\(u_{mf}\)：PDF 8 左 / S. XII。\(\varepsilon_{mf}\)、\(\varphi_S\)、\(\bar\varphi_S\)：PDF 8 右 / S. XIII。\(\rho_P\)、\(\mu_g\)、\(\nu_g\)：PDF 9 左 / S. XIV。\(d_P\)、\(\bar d_P\)、\(d_{P,B}\)、\(\bar d_{P,B}\)：PDF 6 左 / S. VIII。\(u_{mf,0}\) 没有单独的符号行。
❓ **Q2 变量与单位**：\(u_{mf}\) 为 m/s。\(\varepsilon_{mf}\)、\(\varphi_S\)、\(\bar\varphi_S\) 的单位列是“−”。\(\rho_P\) 为 kg/m³。\(\mu_g\) 为 kg/(m·s)（即 Pa·s）。\(\nu_g\) 为 m²/s。\(d_P\) 与 \(\bar d_P\) 为 m；\(d_{P,B}\) 是燃料粒径，不是床料粒径。
❓ **Q3 德文原句**：„Leerrohrgeschwindigkeit bei Minimalfluidisation“；„Lückengrad bei Minimalfluidisierung“；„Sphärizität der Partikel“；„mittlere Sphärizität der Partikel“；„Dichte des Feststoffs“；„mittlerer Partikeldurchmesser“；„Brennstoffpartikeldurchmesser“.
❓ **Q4 上游依赖**：式 (3.1)–(3.4) 使用 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)、\(\bar d_P\)、\(\rho_P\) 以及 \(\mu_g\) 或 \(\nu_g\)。表 2.2 把床料粒径列为输入（第 11 节）。
❓ **Q5 原文空白**：符号表没有数值。\(u_{mf}\) 的德文用 Fluidisation，\(\varepsilon_{mf}\) 的德文用 Fluidisierung，两处都没有给出取值。\(u_{mf,0}\) 只出现在式 (3.3)。同一式 (3.2) 里 Archimedes 数用 \(\bar d_P\)、Reynolds 数用 \(d_P\)，符号表把二者定义为两个符号。气体密度与粘度的数值不在符号表中。
❓ **Q6**：未查看归档。

## 2. 式 (3.1)、式 (3.2)（Ergun 形式）

❓ **Q1 出处**：PDF 20 右 / S. 21，§3.1.1。
❓ **Q2 变量与单位**：

\[
Ar=150\cdot\frac{1-\varepsilon_{mf}}{\varphi_S^2\cdot\varepsilon_{mf}^3}\,Re_{mf}+\frac{1{,}75}{\varphi_S\cdot\varepsilon_{mf}^3}\,Re_{mf}^2
\tag{3.1}
\]

\[
Ar=\frac{g\cdot\bar d_P^3}{\nu_g^2}\cdot\frac{\rho_P-\rho_g}{\rho_g},\qquad Re_{mf}=\frac{d_P\cdot u_{mf}}{\nu_g}
\tag{3.2}
\]

\(Ar\)、\(Re_{mf}\)、\(\varepsilon_{mf}\)、\(\varphi_S\) 无量纲。\(g\) 的符号表单位是 m/s²（S. VIII）。\(\nu_g\) 为 m²/s，\(\rho\) 为 kg/m³，\(u_{mf}\) 为 m/s，粒径为 m。150 与 1,75 按原式书写。
❓ **Q3 德文原句**：„Ergun (1952) gibt zur Berechnung der minimalen Fluidisierungsgeschwindigkeit folgende, auf der Bilanzierung der Trägheits- und Reibungskräfte beruhende, Gleichung an“.
❓ **Q4 上游依赖**：正文写明，床内流动特征量需要操作条件下的空塔速度和最小流化时的空隙率。式 (3.3) 由环境条件下的测量反求 \(\varphi_S\cdot\bar d_P\)，式 (3.4) 再求操作条件下的 \(u_{mf}\)。
❓ **Q5 原文空白**：本页没有 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\) 的数值，也没有点名 VTT 或 HTW 的床料。式 (3.2) 左端粒径带平均号、右端 Reynolds 数的粒径不带平均号，原文未说明二者是否同一个直径。
❓ **Q6**：未查看归档。

## 3. 式 (3.3)、式 (3.4)（阻塞项 X6：\(u_{mf,0}\) 与操作条件下的 \(u_{mf}\)）

❓ **Q1 出处**：PDF 21 左 / S. 22，§3.1.1，紧接在 Werther（1977）那段之后。
❓ **Q2 变量与单位**：

\[
\varphi_S\cdot\bar d_P
=\frac{1{,}75\cdot\rho_g\cdot u_{mf,0}^2}{2\cdot\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}
+\sqrt{
\frac{150\cdot(1-\varepsilon_{mf})\cdot\mu_g\cdot u_{mf,0}}{\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}
+\left(\frac{1{,}75\cdot\rho_g\cdot u_{mf,0}^2}{2\cdot\varepsilon_{mf}^3\cdot(\rho_P-\rho_g)\cdot g}\right)^2
}
\tag{3.3}
\]

\[
u_{mf}=42{,}9\cdot\frac{(1-\varepsilon_{mf})\cdot\nu_g}{\varphi_S\cdot\bar d_P}
\left(
\sqrt{1+\frac{3{,}1\cdot 10^{-4}\cdot\varepsilon_{mf}^3\cdot g\cdot(\varphi_S\cdot\bar d_P)^3\cdot(\rho_P-\rho_g)}{(1-\varepsilon_{mf})^2\cdot\nu_g^2\cdot\rho_g}}-1
\right)
\tag{3.4}
\]

\(u_{mf}\) 与 \(u_{mf,0}\) 为 m/s。式 (3.3) 用动力粘度 \(\mu_g\)（kg/(m·s)），式 (3.4) 用运动粘度 \(\nu_g\)（m²/s）。42,9 与 \(3{,}1\cdot 10^{-4}\) 在式中是无单位系数。左端 \(\varphi_S\cdot\bar d_P\) 的单位是 m。
❓ **Q3 德文原句**：„Werther (1977) vor, die Lockerungsgeschwindigkeit experimentell bei Umgebungsbedingungen zu bestimmen“；„kann zunächst die Sphärizität des Wirbelgutes bestimmt werden“；„Mit bekannter Sphärizität des Wirbelgutes erhält man die minimale Fluidisierungsgeschwindigkeit bei Betriebsbedingungen.“
❓ **Q4 上游依赖**：式 (3.3) 的输入是测得的环境条件 Lockerungsgeschwindigkeit \(u_{mf,0}\)、与之相应的 Lückengrad \(\varepsilon_{mf}\)、平均粒径，以及环境条件下的物性。式 (3.4) 再用已知的球形度求操作条件下的 \(u_{mf}\)。Hamel 引用的是 Werther（1977）。参考文献索引把该条放在 PDF 111 左 / S. 202；本答卷未再打开该页核对书目排版，也未打开 Werther 原文。
❓ **Q5 原文空白**：全书没有哪一个工况的 \(u_{mf,0}\) 测量值。褐煤焦炭、焦炭加灰、VTT、HTW 的 \(\varepsilon_{mf}\) 与 \(\varphi_S\) 原文未给出。正文说式 (3.3) 使用环境条件下的物性，式面上只有 \(u_{mf}\) 带下标 0，\(\rho_g\)、\(\mu_g\)、\(\varepsilon_{mf}\) 没有这个下标。原文未说明操作条件下的 \(\varepsilon_{mf}\) 是否改用另一个数。原文未说明模拟程序对第 7 章工况是走式 (3.3)–(3.4)，还是走下一节的常数 \(C_1\)、\(C_2\)。
❓ **Q6**：未查看归档。

## 4. 式 (3.5) 与文献中的 \(C_1\)、\(C_2\)

❓ **Q1 出处**：PDF 21 左 / S. 22，紧接式 (3.4)。
❓ **Q2 变量与单位**：

\[
\frac{Ar}{Re_{mf}}=C_1+C_2\cdot Re_{mf}^{m}
\tag{3.5}
\]

\(C_1\) 的文献汇编值是 18 到 34，\(C_2\) 是 0,0356 到 0,067（Verweyen 1994）。二者无量纲。指数印为 \(m\)，本页没有给 \(m\) 的数值。
❓ **Q3 德文原句**：„wird dieser Einfluss durch die Wahl von konstantem \(C_1\) und \(C_2\) eliminiert“；„für \(C_1\) Werte von 18 bis 34 und für \(C_2\) Werte von 0,0356 bis 0,067 (Verweyen 1994)“；„Der Ansatz von Ergun (1952) liefert nach Lippens und Mulder (1993) … die besten Ergebnisse.“
❓ **Q4 上游依赖**：Hamel 写明，式 (3.1) 里的 \(C_1\)、\(C_2\) 由球形度和空隙率决定；换成常数后这一影响被拿掉。Lippens und Mulder（1993）比较了 30 个以上的 \(u_{mf}\) 关联和 80 组已发表测量。同页还写 Olowson und Almstedt（1991）与 King und Harrison（1982）表明，更高压力下 Lockerungsgeschwindigkeit 仍可用 Ergun 的式子计算。
❓ **Q5 原文空白**：Hamel 没有从 18–34 和 0,0356–0,067 里选一对，也没有把它们写成 VTT 或 HTW 的床料参数。紧接着的句子写，这些差别表明球形度和最小流化空隙率的影响不宜忽略，并且各式在实验范围之外还要检验适用范围。本页同时写了“常数 \(C_1\)、\(C_2\) 消去这一影响”和“Ergun 的结果最好”。原文未说明程序采用哪一种。指数 \(m\) 原文未给出。
❓ **Q6**：未查看归档。

## 5. Bild 3.6 上的 \(\varepsilon_{mf}=0{,}46\)、\(\varphi_S=0{,}8\)

❓ **Q1 出处**：PDF 27 右 / S. 35，Bild 3.6 图内条件框。图题：不同计算式给出的气泡–悬浮相交换系数。
❓ **Q2 变量与单位**：图内为 \(T_{WS}=298{,}15\,\mathrm{K}\)，\(d_P=0{,}7\,\mathrm{mm}=7\cdot 10^{-4}\,\mathrm{m}\)，\(\rho_P=2600\,\mathrm{kg/m^3}\)，\(\varepsilon_{mf}=0{,}46\)，\(\varphi_S=0{,}8\)，\(u_0-u_{mf}=0{,}6\,\mathrm{m/s}\)，\(d_b=0{,}3\,\mathrm{m}\)。图内球形度符号没有上划线。\(u_{mf}\) 本身没有单独的数，给出的是过剩气速。
❓ **Q3 德文原句**：图内符号行即条件；正文在图前有 „Stoffaustauschkoeffizienten von \(K_{bd}\approx 33\,\mathrm{s}^{-1}\) für kleine Blasendurchmesser“，该句对应的温度是 \(900^\circ\mathrm{C}\)。
❓ **Q4 上游依赖**：此图比较 Kunii und Levenspiel（1977）、Sit und Grace（1981）、Preto（1986）与式 (3.50)。粒径和密度与 S. 32–33 的 Olowson 图相同，图上没有写“同一物料”。
❓ **Q5 原文空白**：物料名称原文未给出。没有一句把 0,46 和 0,8 写成褐煤焦炭、VTT 或 HTW 的床料参数。图内温度是 298,15 K，同页正文的 \(900^\circ\mathrm{C}\) 是另一句关于 \(K_{bd}\approx 33\,\mathrm{s}^{-1}\) 的估算。可选解释：(a) 0,46 与 0,8 只属于这张比较图；(b) 程序在没有别的数据时沿用它们。两种读法都改变 X6。待总司令裁定。本答卷不选取。
❓ **Q6**：未查看归档。

## 6. S. 45 的平均球形度 0,8（自由空间轨迹）

❓ **Q1 出处**：PDF 32 右 / S. 45，§3.2.4，式 (3.86)–(3.87) 与表 3.1 之间。定义句在 PDF 32 左 / S. 44。
❓ **Q2 变量与单位**：球形度无量纲。轨迹计算设 \(0{,}3\le\varphi_S\le 1{,}0\)。无数据时 \(\bar\varphi_S=0{,}8\)，符号带上划线，与符号表的“mittlere Sphärizität”一致。Haider und Levenspiel（1989）的拟合从 \(\varphi_S=1\) 到 0,026。0,026 是该关联覆盖的盘状颗粒，不是床料取值。
❓ **Q3 德文原句**：„Der entsprechende Mittelwert der Sphärizität ist eine Eingabegröße des Simulationsprogramms, sind keine Angaben bekannt, wird den Berechnungen eine mittlere Sphärizität von \(\bar\varphi_S=0{,}8\) zugrundegelegt.“ 定义句：„Das Verhältnis der Oberfläche einer dem Partikel volumengleichen Kugel bezogen auf die tatsächliche Oberfläche des Partikels.“
❓ **Q4 上游依赖**：这段属于自由空间颗粒轨迹的阻力系数，式 (3.82) 起。它不在 §3.1.1 的最小流化推导里。
❓ **Q5 原文空白**：原文未说明这句里的“keine Angaben”是否覆盖第 7 章的床料，也未说明 \(\bar\varphi_S\) 就是式 (3.3)、(3.4) 里的 \(\varphi_S\)。可选解释：(a) 床料没有球形度数据时，式 (3.3) 的 \(\varphi_S\) 也取 0,8；(b) 0,8 只进入自由空间轨迹。待总司令裁定。本答卷不选取。
❓ **Q6**：未查看归档。

## 7. Olowson 比较图上的粒径与密度（无 \(\varepsilon_{mf}\)、无 \(\varphi_S\)）

❓ **Q1 出处**：PDF 26。左页 S. 32 的 \(u_{b,r}\) 图与右页 S. 33 的 Bild 3.4 都有条件框。Bild 3.5（气泡直径，S. 33 下部）的可见图例是各条 \(u_0-u_{mf}\) 曲线。图题写测量来自 Olowson und Almstedt（1990），计算用式 (3.14)–(3.44)。
❓ **Q2 变量与单位**：\(T=25^\circ\mathrm{C}=298{,}15\,\mathrm{K}\)，\(d_P=0{,}7\,\mathrm{mm}\)，\(\rho_P=2600\,\mathrm{kg/m^3}\)，\(h=0{,}55\,\mathrm{m}\)。曲线标的是 \(u_0-u_{mf}=0{,}1\)、0,2、0,3、0,4、0,6 m/s。\(u_{mf}\) 没有单独的数。条件框里没有 \(\varepsilon_{mf}\) 和 \(\varphi_S\)。
❓ **Q3 德文原句**：„Messwerte nach Olowson und Almstedt 1990“；„Berechnung für \(u_0-u_{mf}\)“.
❓ **Q4 上游依赖**：正文（S. 32）说式 (3.42)–(3.44) 的指数用这组气泡穿流测量拟合。物料名称在这两页的图内原文未给出。
❓ **Q5 原文空白**：这是所引实验的颗粒条件，原文未把它写成 Hamel 气化炉的床料。\(u_{mf}\) 原文未给出，只给出过剩气速。
❓ **Q6**：未查看归档。

## 8. Bild 3.7 的粒径与密度

❓ **Q1 出处**：PDF 31 左 / S. 42，Bild 3.7。
❓ **Q2 变量与单位**：\(T_{WS}=1073\,\mathrm{K}\)，\(p=0{,}1\,\mathrm{MPa}\)，\(d_P=1{,}0\,\mathrm{mm}\)，\(\rho_P=2200\,\mathrm{kg/m^3}\)。图内没有 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)。
❓ **Q3 德文原句**：图题 „Blasendurchmesser und -geschwindigkeit sind entsprechend der in Kapitel 3.1 vorgestellten Gleichungen berechnet.“
❓ **Q4 上游依赖**：图题写气泡直径和速度按第 3.1 节的式子计算。物料名称原文未给出。
❓ **Q5 原文空白**：原文未说明 1,0 mm 与 2200 kg/m³ 属于哪一种床料，也未说明它用于第 7 章。
❓ **Q6**：未查看归档。

## 9. Wein 的砂与灰；式 (3.12)

❓ **Q1 出处**：PDF 22 右 / S. 25。Geldart D 与式 (3.11) 在同一半页的上一段。
❓ **Q2 变量与单位**：Wein（1992）的空塔速度至 3 m/s。他用的 Wirbelgut 是砂和灰，过剩气速分别为 2,87 m/s 与 2,97 m/s。关联是 \(u_d=1{,}45\cdot u_{mf}\)，式 (3.12)。1,45 无量纲。砂和灰的 \(u_{mf}\)、粒径、密度原文未给出。Hilligardt 的 Geldart D 实验，按空截面计的过剩气速至 0,5 m/s；式 (3.11) 为 \(u_d=u_{mf}+\frac13(u_0-u_{mf})\)，适用于技术布风器、三维床。
❓ **Q3 德文原句**：„bei den von ihm verwendeten Wirbelgütern Sand und Asche Überschussgeschwindigkeiten von 2,87 bzw. 2,97 m/s“.
❓ **Q4 上游依赖**：Hamel 引用 Wein（1992）与 Hilligardt（1986），以及 Geldart（1973）的分类名称。分类边界不在本页。
❓ **Q5 原文空白**：砂与灰的 \(u_{mf}\) 原文未给出。原文未说明这两种物料是 VTT 或 HTW 的床料。
❓ **Q6**：未查看归档。

## 10. 式 (3.15) 的 \(u_{mf}\) 区间

❓ **Q1 出处**：PDF 23 左 / S. 26。式 (3.14) 在同一半页：Gruppe-B，技术布风器 \(\psi_b=0{,}76\)，多孔板 \(\psi_b=0{,}67\)。
❓ **Q2 变量与单位**：

\[
\psi_b=0{,}17\cdot u_{mf}^{-0{,}33}\qquad\text{für }0{,}014\,\mathrm{m/s}<u_{mf}<0{,}18\,\mathrm{m/s}
\tag{3.15}
\]

区间端点的单位是 m/s。这是关联的适用区间，不是某一床的 \(u_{mf}\)。
❓ **Q3 德文原句**：„die nur noch vom eingesetzten Feststoff abhängig ist“.
❓ **Q4 上游依赖**：Hamel 引用 Wein（1992）。式 (3.14) 写明 Gruppe-B-Partikel。
❓ **Q5 原文空白**：没有把这个区间里的某个数写成 VTT 或 HTW 的 \(u_{mf}\)。Geldart B 的边界原文未印出（见第 16 节）。
❓ **Q6**：未查看归档。

## 11. 表 2.2：床料粒径是输入，表内无数值

❓ **Q1 出处**：PDF 16 右 / S. 13，Tabelle 2.2。
❓ **Q2 变量与单位**：表内没有床料粒径、密度、\(\varepsilon_{mf}\)、\(\varphi_S\) 或 \(u_{mf}\) 的数。床高也没有数。图中的 \(\dot m_B\)、\(x_i(d_{P,B})\) 是燃料流与燃料粒径。
❓ **Q3 德文原句**：„Zusammensetzung und Partikelgrößenverteilung des Feststoffinventars oder Angabe eines mittleren Partikeldurchmessers des Feststoffinventars“；„Wirbelbetthöhe … über Druckdifferenzmessungen konstant gehalten wird“.
❓ **Q4 上游依赖**：床料组成和粒径分布，或者一个平均粒径，是模拟的输入。床高用压差保持恒定，所以床高也是输入。
❓ **Q5 原文空白**：表 2.2 不提供第 7 章任一工况的这些输入的数值。
❓ **Q6**：未查看归档。

## 12. 第 7 章各工况与表 7.1–7.5

❓ **Q1 出处**：HTW 正文 PDF 72–75，表 7.1 在 PDF 74 左 / S. 128（图像已看）。VTT 加压 PDF 76–77，表 7.2 在 PDF 76 右 / S. 133（图像已看），续表在 PDF 77。VTT 常压 PDF 79，表 7.3 表头在 S. 139（图像已看表头）。WSV400 PDF 81–82，表 7.4 在 PDF 81 右 / S. 143。参数研究 PDF 83 右 / S. 147，表 7.5（图像已看）。第 8 章为 PDF 87–90 / S. 154–161。
❓ **Q2 变量与单位**：表 7.1、7.2、7.5 的行是压力、燃料量、气化剂、\(\lambda_{O2}\)、再循环、燃料水/灰/挥发分/元素。没有 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)、床料密度、床料粒径。表 7.5 另有床高 \(1{,}2\,\mathrm{m}\) 和 \(u_0=0{,}75\,\mathrm{m/s}\)；标准工况燃料是泥炭，装置是 Kurkela und Stahlberg（1992）的加压气化炉，\(\lambda_{O2}=0{,}3\)，蒸汽/燃料（无水无灰）\(0{,}15\,\mathrm{kg/kg}\)，有再循环。\(u_0\) 是空塔速度。表 7.3 表体、表 7.4 由逐页文本定位，行的种类与表 7.1、7.2 相同，未见床料参数行。
❓ **Q3 德文原句**：表 7.2 标题 „Zusammenstellung der wichtigsten Parameter der berechneten Betriebsfälle der VTT-Druckvergasungsanlage“. S. 139：„Die im Simulationsprogramm eingesetzte Partikelgrößenverteilung entspricht der von Leppälahti und Kurkela (1992) angegebenen.“ S. 144（PDF 82，逐页文本）：木片按等表面积折成球形颗粒的直径级，见 Bild 7.15；该句的对象是燃料试样。
❓ **Q4 上游依赖**：表 7.2 的数据 Hamel 写明取自 Kurkela und Stahlberg（1992）。S. 139 的粒径分布指向 Leppälahti und Kurkela（1992），本页没有分布的数值，句子也没有写“床料”。
❓ **Q5 原文空白**：第 7 章对 HTW、VTT 加压（含工况 1）、VTT 常压、WSV400 和表 7.5 的标准工况，床料的 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)、\(u_{mf,0}\)、颗粒密度和床料粒径都是原文未给出。S. 139 只给出文献名，没有粒径数字。第 8 章逐页文本没有 Lückengrad、Sphärizität、Minimalfluidisation、Bettmaterial 的数值；该章图是气体组成和温度的 parity。
❓ **Q6**：未查看归档。

## 13. 表 5.1（S. 79）的惰性料密度 2000 kg/m³

❓ **Q1 出处**：PDF 49 右 / S. 79。此表编号与 S. 74 的另一张表 5.1 重复；S. 74 那张是烧损模型汇编，不是本表。引用须带页码。
❓ **Q2 变量与单位**：\(\rho_{inert}=2000\,\mathrm{kg/m^3}\)。同表还有 \(u_0=1\,\mathrm{m/s}\)，\(T_{WS}=1123\,\mathrm{K}\)，压力 0,1 与 1,0 MPa，4 个 cell，6 个粒径级。没有各粒级的直径，没有 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf}\)。
❓ **Q3 德文原句**：„Dichte des Inertmaterial \(\rho_{inert}\)“；表题 „Randbedingungen und Annahmen die den vergleichenden Berechnungen zum Einfluß der Vergasungskinetik zugrundeliegen“. 燃料行是 „Koks“.
❓ **Q4 上游依赖**：这是第 5 章一套比较计算的边界条件。正文写加入的是“新鲜”焦炭颗粒，六个粒级质量份额相等。惰性料的名称原文未给出。
❓ **Q5 原文空白**：原文未说明 2000 kg/m³ 用于 VTT 或 HTW。它不是表 7.2 工况 1 的床料密度。
❓ **Q6**：未查看归档。

## 14. 附录与其它粒径、密度（不是第 7 章床料）

❓ **Q1 出处**：符号表把 \(d_{P,B}\) 定义为燃料粒径（S. VIII）。附录 A2.4 的图在 PDF 98–103（S. 176 起，top/bottom）。干燥模型的颗粒密度在 PDF 37–38 / S. 54–57 的表中（逐页文本：干颗粒密度 1250 kg/m³ 与 500 kg/m³）。PDF 58 的焦油动力学表有一列 Bettmaterial，条目为 Quarz、Sand、Dolomit 等名称（逐页文本）。
❓ **Q2 变量与单位**：附录图上的 \(d_{P,B}\) 是该次单颗粒试验的燃料或其余焦粒径，并标有样品质量 \(m_B\) 和 \(u_{0,0}\)。这些图不提供气化炉床料的 \(\varepsilon_{mf}\)、\(\varphi_S\) 或 \(u_{mf}\)。本答卷不把各次试验的 \(d_{P,B}\) 逐条登录为床料粒径。
❓ **Q3 德文原句**：符号表 „Brennstoffpartikeldurchmesser“.
❓ **Q4 上游依赖**：附录 A2 汇编的是热解、燃烧和气化试验，供第 4、5 章动力学使用。
❓ **Q5 原文空白**：附录没有第 7 章装置的床料 \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf,0}\) 或床料密度。第 5、6 章把“Formfaktor”用作反应级数 \(n\) 的别名（约 0,5–0,6），那个词不是球形度 \(\varphi_S\)。
❓ **Q6**：未查看归档。

## 15. Hamel 使用 Geldart 类别名称、但没有印出边界

❓ **Q1 出处**：PDF 21 右 / S. 23（细颗粒 A、B、C，粗颗粒 D，并写 siehe Geldart 1973）。PDF 22 右 / S. 25（Hilligardt 的 D 类；Wein 的 B 类）。PDF 23 左 / S. 26（式 (3.14) 的 Gruppe B）。PDF 25 / S. 30–31 的逐页文本有“压力升高时，一类颗粒的行为接近更细的下一类”。参考文献条目在 PDF 105。
❓ **Q2 变量与单位**：这些页给出类别名称和上节已录的速度、系数。没有 A/B、B/D、A/C 的分界式或分界图。
❓ **Q3 德文原句**：„feinen Partikeln (Geldart A, B, C; siehe Geldart 1973)“；„Grobkornwirbelschichten (Geldart D)“.
❓ **Q4 上游依赖**：类别名称用来区分式 (3.11)、(3.12)、(3.14) 的适用范围。判据本身指向 Geldart（1973）。
❓ **Q5 原文空白**：Hamel 没有印出分类边界。只凭这些页不能判断 VTT 工况 1 的床料属于哪一类。
❓ **Q6**：未查看归档。

---

## 16. Geldart (1973) 的 A/B、B/D、A/C 分界

文献：D. Geldart, *Types of Gas Fluidization*, Powder Technology, 7 (1973) 285–292。本地文件不入库。下面的印刷页是期刊页码；PDF 第 1 页即 p. 285。只录分界。式 (1)、(2) 是最小鼓泡速度与粒径的关联，式 (3) 是最小流化速度关联，式 (9) 起是与其它判据的比较，都不转录。不计算 VTT 床料的类别。

符号以 p. 291 的符号表为准（页面图像）：

| 符号 | 原文 | 单位 |
|---|---|---|
| \(d_{sv}\) | surface/volume diameter of particle | cm |
| \(U_0\) | superficial velocity of gas at minimum fluidization | cm s\(^{-1}\) |
| \(U_{MB}\) | superficial velocity of gas at minimum bubbling condition | cm s\(^{-1}\) |
| \(K_{MB}\) | constant in eqn. (1) | s\(^{-1}\) |
| \(\rho_s\) | density of particle (including any internal porosity) | g cm\(^{-3}\) |
| \(\rho_f\) | density of fluidizing gas | g cm\(^{-3}\) |
| \(\mu\) | viscosity of fluidizing gas | g cm\(^{-1}\) s\(^{-1}\) |
| \(\varepsilon_0\) | bed voidage at minimum fluidization velocity | — |
| \(d_B\) | frontal diameter of bubble | cm |
| \(g\) | gravitation constant | 981 cm s\(^{-2}\) |
| \(x\) | weight fraction of particles in each size range | — |

Geldart 的 \(U_0\) 是最小流化时的空塔速度。p. 288 把平均粒径写成 \(\bar d_{sv}(=1/\sum x/d_{sv})\)，符号表称 \(d_{sv}\) 为表面/体积直径。行内求和号两边都用了 \(d_{sv}\)，没有另给粒级粒径的符号。

### 16.1 A/B：式 (4)、(5)、(6)

❓ **Q1 出处**：p. 289（PDF 第 5 页）右栏，§4.3。式 (6) 画在 Fig. 3（p. 290）上，标为 XY。p. 288 的正文用了另一个不等号，见 Q5。
❓ **Q2 变量与单位**：属于 A 类的判据印为

\[
\frac{U_{MB}}{U_0}\ge 1
\tag{4}
\]

\[
\frac{8\times 10^{-4}\, g\, d_{sv}(\rho_s-\rho_f)}{K_{MB}\mu}\le 1
\tag{5}
\]

\[
(\rho_s-\rho_f)\,d'\le 225
\tag{6}
\]

式 (5) 中 \(d_{sv}\) 是一次方。式 (6) 把 \(d_{sv}\)（cm）换成 \(d'\)（粒径，μm）。同页写明，空气、室温、常压时 \(K_{MB}=100\)，\(\mu=1.8\times 10^{-4}\) poise。\(K_{MB}=100\) 的单位条件在 p. 288：\(U_{MB}\) 用 cm s\(^{-1}\)，\(d_{sv}\) 用 cm。式 (6) 的 225 旁边没有再印单位；同页例子把密度差写成 g cm\(^{-3}\)（“density difference 1 g cm\(^{-3}\) larger than about 250 μm”）。
❓ **Q3 原句**：„Thus for a powder to belong to group A“；„Now for air at room temperature and pressure \(K_{MB}=100\) and \(\mu=1.8\times 10^{-4}\) poise.“；„Equation (6) has been plotted on Fig. 3 as XY“. p. 288：„(\(U_{MB}/U_0>1\) for group A powders)“.
❓ **Q4 上游依赖**：p. 289 写明式 (5) 由式 (1) 与式 (3) 代入式 (4)。式 (1)、(3) 本身不在本节约录。Fig. 3 的图题：„Powder classification diagram for fluidization by air (ambient conditions).“ 横轴 „Mean particle size \(d_{sv}\) (μm)“，纵轴 „Density difference \((\rho_s-\rho_f)\) (g/cm\(^3\))“。
❓ **Q5 原文空白**：p. 288 正文是 \(U_{MB}/U_0>1\)，式 (4) 是 \(\ge 1\)。两处都照印，这里不取舍。式 (6) 的界是 225，同页散文把密度差 1 g cm\(^{-3}\) 的 B 类写成 “larger than about 250 μm”。式 (6) 的适用句是空气、室温、常压。p. 291 写，压力升高会增大式 (6) 右端的常数，从而把 Fig. 3 上的 XY 向右移；该页没有给出加压后的新常数，并写 „we do not know the effect on \(K_{MB}\) of simultaneously reducing the gas density“. 因此式 (6) 没有给出非空气或非常压的数值形式。
❓ **Q6**：未查看归档。

### 16.2 B/D：式 (7)、(8)，以及喷动床一句

❓ **Q1 出处**：式 (7) 在 p. 289 右栏末，§4.4。式 (8) 与喷动床一句在 p. 290 左栏。
❓ **Q2 变量与单位**：作者写，气泡上升慢于颗粒间隙气速的密度–粒径组合用式 (7)。属于 D 类的数值形式是式 (8)。

\[
\left(g\frac{d_B}{2}\right)^{1/2}
\le
\frac{8\times 10^{-4}\,(\rho_s-\rho_f)\,g\,d_{sv}^{2}}{\mu\varepsilon_0}
\tag{7}
\]

\[
(\rho_s-\rho_f)\,(d')^{2}\ge 10^{6}
\tag{8}
\]

式 (7) 左端的指数是 \(1/2\)；下一页写 „in eqn. (7) we are considering \(\sqrt{d_B}\)“. 式 (7) 右端的 \(d_{sv}\) 是平方。代入式 (8) 时作者取 \(d_B=25\) cm，大颗粒体系 \(\varepsilon_0\approx 0.4\)，空气 \(\mu=1.8\times 10^{-4}\) g cm\(^{-1}\) s\(^{-1}\)，并用 \(d'\)（μm）替换 \(d_{sv}\)（cm）。\(10^{6}\) 旁边没有另印单位。第二条判据没有式子：D 类粉末能在深于 30 cm 的床里维持稳定喷动。已报道为可喷动的密度–粒径组合在 Fig. 3 上用叉号表示。叉号的坐标不在这里读取。
❓ **Q3 原句**：„The distinction between groups B and D is not as clear-cut as that between A and B“；„so let us choose \(d_B=25\) cm“；„For large-particle systems \(\varepsilon_0\approx 0.4\)“；„we obtain, for group D“；„group D powders are capable of maintaining a stable spout in a bed more than 30 cm deep“.
❓ **Q4 上游依赖**：式 (7) 的右端用了式 (3)。作者接着写，对大颗粒这样用式 (3) „is not strictly justified“，因为流态是过渡区，不是层流；他又写 \(d_B\) 的选取是任意的，因此不改用更复杂的 \(U_0\) 式。喷动床一句引用 Baeyens，并写实验仍在进行。
❓ **Q5 原文空白**：B 与 D 被写成连续过渡，不是一条非此即彼的线。式 (8) 的 \(10^{6}\) 依赖 \(d_B=25\) cm 和 \(\varepsilon_0\approx 0.4\) 这两个作者写入的数。式 (8) 与式 (6) 一样，数值形式对应的是空气和上述代入，不是加压气化条件。p. 290 写，Fig. 3 上的叉号落在“the line”附近或右侧，并说这给式 (8) 一些支持。判据以式 (8) 的印刷式为准，不从图上辨认曲线、不读坐标。30 cm 喷动床这一条没有 \((\rho_s-\rho_f)\)–\(d\) 的式子，叉号也不读。
❓ **Q6**：未查看归档。

### 16.3 A/C：只有图

❓ **Q1 出处**：p. 289 右栏，紧接在式 (6) 的数据说明之后；图是 Fig. 3，p. 290。
❓ **Q2 变量与单位**：没有 A/C 的分界式。原文把 C 与 A 的分界写成 Fig. 3 上的阴影区 PQ。坐标轴与 16.1 的 Fig. 3 相同：横轴 \(d_{sv}\)（μm），纵轴 \((\rho_s-\rho_f)\)（g/cm\(^3\)），空气、环境条件。空心点是作者们认为极难流化的粉末（C 类），半实心点是观察到无气泡床层膨胀的粉末（A 类）。
❓ **Q3 原句**：„The division between groups C and A is represented by the shaded area PQ which has been drawn so as to separate empirically open and half-closed points.“ §3.3（p. 287）对 C 类只有定性句：„Powders which are in any way cohesive belong in this category“，并把它联系到粒径很小、静电或湿粘物料。该节没有粒径或密度的数值界。
❓ **Q4 上游依赖**：PQ 是按图上的点经验画出的。没有写出 PQ 的方程。
❓ **Q5 原文空白**：A/C 只有图。本答卷不目测阴影的坐标，不把图上的刻度读成分界值。是否另开数字化，由总司令决定。
❓ **Q6**：未查看归档。

### 16.4 第 3 节的描述性范围（不是 §4 的分界式）

p. 286，§3.1：A 类是平均粒径小，并且/或者颗粒密度低，“less than about 1.4 g/cm\(^3\)”。p. 286，§3.2：B 类的大多数物料落在 \(40\,\mu\mathrm{m}<d_{sv}<500\,\mu\mathrm{m}\)、\(4\,\mathrm{g/cm^3}>\rho_s>1.4\,\mathrm{g/cm^3}\)，砂是最典型的粉末。这两句写的是 “most materials” 和 “about”，不是式 (6) 或式 (8)。p. 287，§3.4：D 类限于大颗粒和/或很密的颗粒；除最大的气泡外，气泡上升慢于间隙气体。数值界是式 (8)，不在这一节。

---

## 17. VTT 工况 1：若做 Geldart 分类，各量从哪里来

床料按 ADR 0014 A1：燃料生成的焦炭与灰，不含惰性砂。分界已在第 16 节转录。本节只列各量的出处，不把工况 1 代入式 (4)–(8)，不读 Fig. 3，不给出 A、B、C、D 中的任何一类。

### 平均粒径 \(\bar d_P\)

- Hamel 表 7.2 工况 1 没有床料粒径。表 2.2 把它列为输入，没有给数。
- ADR 0015 第 11 项：床区 \(\bar d_P=\sum_k m_k d_k/\sum_k m_k\)，是床区各 cell 中焦炭各粒级的质量平均；灰按所在焦炭粒级的同一粒径计；床区各 cell 用这同一个平均值，用于式 (3.1)–(3.4) 和 Geldart 分类。存量 \(m_k\) 来自式 (2.6) 的牛顿解，在外层更新。
- 同一项写明，外层迭代的初值取燃料粒径分布的质量平均。燃料粒径分布的出处是 ADR 0014 D2：KURKELA_1992 Fig. 1 褐煤累积曲线的工具读数。本答卷不打开该文献，不读分布，不算这个平均。

### 颗粒密度 \(\rho_P\)

- Hamel 表 7.2 工况 1 没有床料颗粒密度。符号表只有定义（kg/m³）。
- ADR 0014 A4：不单独选取 \(\rho_P\)。几何存量 \(\int(1-\varepsilon_b)(1-\varepsilon_d)\,\rho_P\,A\,\mathrm{d}h\) 在床区 0–1,2 m 上等于 A3 的床料质量；\(\rho_P\) 由这个条件确定，并随 \(\varepsilon_b\)、\(\varepsilon_d\) 在外层更新。
- A4 自己写明：\(\varepsilon_b\)、\(\varepsilon_d\)、\(u_{mf}\) 和 Geldart 分类还依赖床料粒径；粒径未定之前，A4 和分类都不能完成。这个依赖环原文和 ADR 都没有拆开。本答卷不拆。
- ADR 0014 A3 另给了由床层压降得到的床料质量，并写了一个床层平均堆密度。堆密度是质量除以床层体积，不是 \(\rho_P\)。本答卷不把堆密度当作颗粒密度。

### \(\varepsilon_{mf}\)、\(\varphi_S\)、\(u_{mf,0}\)

- 对工况 1 的焦炭与灰，这三项都是原文未给出（第 3、12 节）。
- 式 (3.3) 把 \(\varphi_S\cdot\bar d_P\) 写成 \(u_{mf,0}\)、\(\varepsilon_{mf}\) 和物性的函数。式 (3.4) 再用 \(\varphi_S\)、\(\varepsilon_{mf}\)、\(\bar d_P\) 和操作条件下的 \(\nu_g\)、\(\rho_g\) 求 \(u_{mf}\)。
- Bild 3.6 的 0,46 和 0,8，以及 S. 45 的 \(\bar\varphi_S=0{,}8\)，出处见第 5、6 节。它们没有写成工况 1 的床料。本答卷不拿来分类。

### 分界式要用、而工况 1 尚未闭合的量

式 (6)、式 (8) 的数值形式用的是 \((\rho_s-\rho_f)\) 与 \(d'\)（μm）。式 (5) 还用 \(d_{sv}\)（cm）、\(\mu\)、\(g\) 和 \(K_{MB}\)。式 (7) 还用 \(\varepsilon_0\) 与 \(d_B\)。这些式子都不用 \(\varphi_S\)，也不用 Hamel 的 \(u_{mf,0}\)。

- 粒径：Geldart 的 \(d_{sv}\) 是表面/体积平均，行内写成 \(1/\sum x/d_{sv}\)，\(x\) 是各粒级的质量分数。ADR 0015 第 11 项的 \(\bar d_P\) 是床区焦炭的质量加权算术平均。两处定义不同。本答卷不把它们当成同一个数，也不把 ADR 的平均代入式 (6) 或式 (8)。
- 颗粒密度：Geldart 的 \(\rho_s\) 写明包含内孔，单位 g cm\(^{-3}\)。ADR 0014 A4 的 \(\rho_P\) 由床层压降和几何存量确定，单位在 Hamel 符号表里是 kg/m³。两处没有写成同一个量。本答卷不换算、不代入。
- 气体密度与粘度：式 (6)、式 (8) 的数值系数对应空气、室温、常压（式 (8) 还代入了 \(\varepsilon_0\approx 0.4\)、\(d_B=25\) cm）。工况 1 是加压气化，不是这个条件。p. 291 只定性地写压力升高会使 XY 右移，没有新的常数。\(\mu\) 与 \(\rho_f\) 的工况数值不在这里补。
- \(\varepsilon_0\)：符号表定义为最小流化时的床层空隙率。它出现在式 (7)。式 (8) 已经用作者写入的 \(\varepsilon_0\approx 0.4\) 把它代掉。Hamel 没有给出工况 1 的 \(\varepsilon_{mf}\)。本答卷不用 0.4 代替那个空白。
- A/C：只有 Fig. 3 的阴影 PQ，没有式子。不读图，因此也不列出一条“从图上分类还要哪些数”的清单。

### 交给总司令的空白（不在本答卷中取值）

1. 工况 1 的 \(\varepsilon_{mf}\) 与 \(\varphi_S\)：原文未给出。Bild 3.6 与 S. 45 的数是否许用，见第 5、6 节的两种读法。
2. \(u_{mf,0}\)：原文未给出测量值。式 (3.3) 因此缺一个输入。
3. 式 (3.2) 中 \(\bar d_P\) 与 \(d_P\) 是否同一个直径。
4. 程序走式 (3.3)–(3.4)，还是走常数 \(C_1\)、\(C_2\)。\(C_1\)、\(C_2\) 与指数 \(m\) 都没有选定的值。
5. 式 (4) 的 \(\ge\) 与 p. 288 的 \(>\)；式 (6) 的 225 与同页 “about 250 μm”。A/C 是否数字化 Fig. 3 的 PQ。这些都不在本答卷中取值。工况 1 的类别不计算。
