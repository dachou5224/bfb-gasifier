# T13 原文答卷：Wirsum 限定查阅——床区固体输运、床料存量与惰性床料

| 项 | 内容 |
|---|---|
| 任务卡 | T13（`docs/tasks/T13-wirsum-solids-scope.md`） |
| 执行者 | cursor |
| 会话 | 全新会话。未沿用做过 T4、T8、T9、T11、T12 或 worktree 清理的会话。工作树 `../bfb-T13`，分支 `task/T13-wirsum-solids-scope`，自 `origin/main`（`b1e8d3f`）建立 |
| 日期 | 2026-10-16 |
| 审查方 | codex（对照 Wirsum 与 Hamel 原页逐条审查） |
| 本卡做什么 | 只查证，不采用。每处与 Hamel 的差别只列可选解释和影响。采用与否由总司令在 ADR 0009 中裁定 |

Hamel 页面图像按 hamel-source-check 第一节第 2 项，链接到主工作树的 `refs/hamel/pages/`。Wirsum 页面用 PyMuPDF 1.26.5 以 200 dpi 渲染到主工作树的 `refs/wirsum/pages/`（115 页，每页 `NNN.png` 与 `NNN.txt`），本 worktree 以符号链接指向该目录。文本层只用于定位。下面引用的式子和德文原句已对照页面图像。式号以图像为准。

Wirsum 的 PDF 是双页拼版。下文同时给出 PDF 页、半页和印刷页。左半页为偶数印刷页，右半页为奇数印刷页。

未打开 Wirsum 引用的其它文献，未打开 `docs/hamel_submodels/` 和归档代码。参考文献表只用来抄作者和年份（PDF 110–115），没有读那些文献。

已读范围：§2.2 与 §2.3（PDF 16–21 / S. 12–23）、§3.2 至式 (3.26) 为止（PDF 23–26 / S. 26–33；§3.3 标题之后未再往下读）、§3.4.1（PDF 29 右至 PDF 30 左 / S. 39–40；§3.4.2 起未读）、§5.6（PDF 48 右至 PDF 53 / S. 77–87）。检索例外另记在问题 3 和问题 6。

符号：Wirsum 的粒级下标是 \(l\)，固体种类下标是 \(j\)。Hamel 的粒级下标是 \(k\)。两边的 cell 编号方向相同：\(i-1\) 在上，\(i+1\) 在下。Wirsum 的式 (2.6) 是热解释放的物质集合，不是固体衡算；与 Hamel 式 (2.6) 同号不同义。

---

## 问题 1：固体衡算变量

❓ **Q1 出处**：

- 式 (2.1)–(2.3)：PDF 18 左半页 / S. 16。
- 变量分组：PDF 19 左半页 / S. 18。
- 式 (2.12) 与“只衡算活性固体”：PDF 19 右半页 / S. 19。
- 式 (2.12) 各项的定义和 Bild 2.3：PDF 20 左半页 / S. 20。
- 对照：Hamel S. 12（PDF 16 左）、表 2.2（PDF 16 右 / S. 13）、Bild 2.2（PDF 17 左 / S. 14）、式 (2.1)（PDF 17 右 / S. 15）、式 (2.6)（PDF 18 右 / S. 17，定义在 PDF 19 左 / S. 18）。

❓ **Q2 变量与单位**：

- 式 (2.1)：

\[
N_{Bil}=N_{Phs}N_{Gas}+\sum_{k=1}^{N_{Fst}}N_{dp,k}+1
\]

  每个衡算空间的变量个数。\(N_{Phs}\) 是相数，\(N_{Gas}\) 是被衡算的气体种数，\(N_{Fst}\) 是被衡算的固体种数，\(N_{dp,k}\) 是第 \(k\) 种固体的粒级数。末项 \(+1\) 是能量流。不衡算能量时用式 (2.2)，去掉这个 \(+1\)。全炉变量数 \(N_{tot}=N_Z N_{Bil}\)（式 (2.3)）。
- 式 (2.12) 的每一项都是质量流，kg/s。符号表（PDF 7 右 / S. XI）：\(m\) 为 Masse，kg；\(\dot m\) 为 Massenstrom，kg/s。
- 式 (2.12) 没有存量，也没有比例系数。

❓ **Q3 德文原句**：

- S. 18：变量分成三组，„1. gasförmige Stoffströme  2. Feststoffströme  3. Energieströme“。
- S. 19：„Bilanziert werden nur sogenannte „aktive“ Feststoffe, also solche Stoffe, welche an Umsatzreaktionen beteiligt sind.“ 举例是 Koks 和 Kalkstein。
- S. 16：细颗粒当作飞灰，只沿烟气方向运动，„nicht dem blaseninduzierten Vermischungsverhalten der Feststoffe unterworfen“；粗颗粒才计内部返混，并分成粒级。这一条从 Dersch（1994）的做法扩到任意种固体、任意个粒级。
- 式 (2.12)，S. 19：

\[
\dot m_{j,l,i}^{zu}+\dot m_{j,l,i}^{r}+\dot m_{j,l,i+1}^{auf}-\dot m_{j,l,i}^{auf}+\dot m_{j,l,i-1}^{ab}-\dot m_{j,l,i}^{ab}+\dot m_{j,l,i}^{\leftarrow}-\dot m_{j,l,i}^{\rightarrow}=0
\]

  S. 20 的定义：\(\dot m_{j,l,i+1}^{auf}\) 来自下方 cell \(i+1\)；\(\dot m_{j,l,i}^{auf}\) 从本 cell 向上；\(\dot m_{j,l,i-1}^{ab}\) 来自上方 cell \(i-1\)；\(\dot m_{j,l,i}^{ab}\) 从本 cell 向下。\(\dot m^{\leftarrow}\) 是邻级因收缩或长大进入本粒级，\(\dot m^{\rightarrow}\) 是离开本粒级。\(\dot m^{r}\) 是化学反应的源或汇。式中没有离开体系的流。

❓ **Q4 上游依赖**：

粒级数进入式 (2.1) 的和。细颗粒不进入气泡引起的返混。式 (2.12) 的上下流在这一节没有计算式；计算式若有，在 §3.2 和 §5.6（问题 2）。

❓ **Q5 原文空白**：

**结论：Wirsum 有依据。** 单元模型里，每种被衡算固体的每个粒级计 1 个变量；正文把这一组变量称为固体流。被衡算的只有活性固体。这和 ADR 0007 的选择不同：ADR 0007 把未知量取成粒级存量。

和 Hamel 的对照：

- Hamel 式 (2.1)（S. 15）是 \(N_{Bil}=N_{Phs}\cdot N_{Gas}+\sum_j N_{dp,j}+1\)。与 Wirsum 式 (2.1) 同形。Hamel 也是每个 cell、每种固体、每个粒级 1 个变量。
- Hamel 式 (2.6)（S. 17–18）同样是稳态质量流衡算，项里没有存量。Hamel 的定义里有 \(\dot m_{aus}\)（离开该 cell 并且离开体系）。Wirsum 式 (2.12) 没有这一项。
- Hamel §2.2.2 也把式 (2.6) 的对象写成活性固体，举例是焦炭和石灰石。这一点与 S. 19 的原句同向。
- Hamel S. 12 说流动模型先算出向上、向下的固体流，再把每个 cell 的固体进出流按流动模型分配。Wirsum §2.2 也是先算子模型、再映射到 cell（问题 5）。两边都没有在计数式旁边写“未知量是存量”。
- Hamel 表 2.2 的 Feststoffinventar 只写组成和粒径（或一个平均粒径），没有写总质量。Wirsum 表 2.1 把总质量列为输入（问题 3）。Bild 2.2 把活性固体粒级数和惰性床料粒级数分成两行；Wirsum 式 (2.1) 的和只计被衡算的固体，惰性料不在式 (2.12) 里。

计数仍填不上式 (2.12) 里的多股流：向上、向下、外加、反应、粒级迁移不能同时都是自由未知量。原文没有写哪一股是那 1 个变量。

可选解释：

- (a) 那 1 个变量就是某一股固体流，其余各股由预计算给定。影响：与 ADR 0007 的存量未知量不同；牛顿法的尺度是 kg/s。
- (b) 标题写的是流，未知量仍是存量，各股流由流动模型乘份额得到。影响：与 ADR 0007 同向，但 Wirsum 没有印出份额式。
- (c) 细颗粒和粗颗粒不是同一套变量。细颗粒随气体走，粗颗粒才有向上、向下两股。影响：不能把每一种固体都放进式 (2.12) 的四股输运。

待总司令裁定。本卡不采用。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 式 (2.1) 与 Hamel 式 (2.1) 同形，Hamel 已有自己的计数式，不必改用 Wirsum 的式号。式 (2.12) 与 Hamel 式 (2.6) 是同一类衡算，但 Wirsum 没有 \(\dot m_{aus}\)，也没有写未知量是存量。Hamel 没有把式 (2.12) 引作自己的固体衡算。不能用来替换 Hamel 式 (2.6)，也不能用来取消 ADR 0007。是否把“变量是流”当作对 ADR 0007 的修订材料，由 ADR 0009 裁定。

---

## 问题 2：床区整体向上、向下固体流

❓ **Q1 出处**：

- 单元衡算里的上下流：式 (2.12)，S. 19–20，见问题 1。这一节没有给出这两股流的计算式。
- §3.2：PDF 23 右 / S. 27 起，至 PDF 26 右 / S. 33 的式 (3.26)。这里计算的是气泡、空隙率和床料总质量，没有相邻 cell 的 \(\dot m_{auf}\)、\(\dot m_{ab}\)。
- §5.6.1：PDF 48 右 / S. 77 至 PDF 49 左 / S. 78。循环率 \(w\) 被写成沿床高常数。
- §5.6.2.1：PDF 50 右 / S. 81。\(W_1\) 只含砂，\(W_2\) 只含污泥。
- §5.6.3：PDF 51 / S. 82–83，PDF 52 / S. 84–85。式 (5.12)–(5.22)、式 (5.29)–(5.34)。
- §5.6.4：PDF 52 右 / S. 85 至 PDF 53 / S. 86–87。离析参数和 \(D_{ax}\) 的拟合。系数不抄。

❓ **Q2 变量与单位**：

- 式 (2.12) 的 \(\dot m^{auf}\)、\(\dot m^{ab}\)：kg/s，按固体种类和粒级。
- \(W_1\)、\(W_2\)：正文称为 Feststoffvolumenstrom / Volumenstrom / Zirkulationsstrom。由式 (5.13)、式 (5.30) 的量纲，是 m³/s，不是 kg/s。全床一个数，不带 cell 下标。
- \(f_w=V_w/V_b\)（式 (5.32)），\(f_{bw}=V_w/(V_w+V_b)\)（式 (5.33)），无量纲。模型取 \(f_{bw}=0{,}2\)。
- \(\bar u_{auf}\)：m/s。式 (5.12) 的 \(\Omega\) 在 S. 85 带单位 \((\mathrm{cm/s})^{1/2}\)。

❓ **Q3 德文原句与式子**：

S. 81：„Eine nach oben gerichtete Bewegung der Klärschlammpartikel im Nachlauf oder wake einer Blase wird aufgrund der außerordentlichen Größe der Klärschlammpartikel ausgeschlossen.“ 同页：污泥颗粒的向上运动„auf den Blasendächern und im Blasendrift“。于是„der Blasennachlauf ausschließlich aus Wirbelsand (jetsam) besteht“，并区分

- \(W_1\)：„Zirkulationsstrom … welcher nur aus jetsam-Material (Sand) besteht“
- \(W_2\)：„Zirkulationsstrom … welcher nur aus flotsam-Material (Klärschlamm) besteht“

S. 82，Bild 5.16：轴向微元上，bulk 相边界上的流是 \((W_1+W_2)c_{j,b}\)，wake 相边界上的流是 \(W_1 c_{j,w}\)。图上两相都有指向上方的箭头。这是浓度剖面用的体积流，不是式 (2.12) 的 \(\dot m_{auf,j,l,i}\)。

S. 82，式 (5.12)，Nienow et al.（1978）：

\[
\bar u_{auf}=\Omega(u_0-u_{mf})^{1/2}
\]

S. 83，式 (5.13)–(5.16)：

\[
W_2=\bar u_{auf}(1-\varepsilon_{bw})(1-\varepsilon_d)A\,c_{f,b},\qquad c_{f,b}=1-c_{j,b}
\]

\[
W_2=K_1(1-c_{j,b}),\qquad K_1=\bar u_{auf}(1-\varepsilon_{bw})(1-\varepsilon_d)A
\]

式 (5.34)（S. 85）与式 (5.13) 是同一关系，再次写明对象是 Klärschlamm。

S. 85，式 (5.30)，Werther（1976），„Feststoffvolumenstrom im Nachlauf der Blasen“。高分辨率裁图在 \((1-\varepsilon_{mf})\) 与 \((u_0-u_{mf})\) 之间又看到一个 \(f_w\)。紧接着的式 (5.31)，Bilbao et al.（1988）：

\[
W_1=0{,}67\frac{f_{bw}}{1-f_{bw}}(1-\varepsilon_{mf})(u_0-u_{mf})A
\]

正文写这两式等价，差别只在 \(f_w\) 与 \(f_{bw}\) 的定义。由式 (5.32)、式 (5.33)，\(f_w=f_{bw}/(1-f_{bw})\)，等价只需要一个 \(f_w\)。

**目标所有者 2026-10-16 判读：** 第二个 \(f_w\) 是笔误，删去。确认后的式 (5.30) 为

\[
W_1=0{,}67\,f_w(1-\varepsilon_{mf})(u_0-u_{mf})A
\]

这是对印刷字形的判读，不是把该式采用进 Hamel 主模型。\(f_{bw}=0{,}2\) 这一句是清楚的，依据是 Rajan und Wen（1980）和 Bilbao et al.（1988）的试验。

S. 84：不考虑 wake 与 bulk 的相交换时，尾涡里向上的颗粒质量流与高度无关；尾涡被看成把砂从床底送到床面的通道（Daw et al. 1988）。流动量用单组分床计算，忽略少量燃料颗粒；气泡体积份额取一个沿高度的平均值。

式 (5.17)、式 (5.21) 是 jetsam 浓度的微分方程，未知量是 \(c_{j,b}(z)\)，不是 cell 质量流。式 (5.19) \(\alpha=W_1/K_1\)，式 (5.20) \(\beta=D_{ax}A/(HK_1)\)。式 (5.29)（S. 84）把 bulk 中的平均 jetsam 浓度和全床平均浓度用空隙率、\(f_{bw}\) 连起来；尾涡项没有再乘一个浓度，与“尾涡只含砂”同向。式 (5.24)–(5.28) 是这条剖面的解析解，本答卷不重抄积分常数。

§5.6.4 的式 (5.35)–(5.40) 是离析参数和 \(D_{ax}\) 的拟合式。S. 87：系数由约 100 次试验、每次 1000 个点的浓度曲线拟合（Wintruff 1996）。**这些是污泥离析的数值关联式，本答卷不抄系数，也不引入。**

§3.2 与固体流直接有关的只有两句定性说明，没有 cell 流量：S. 27 把“悬浮相气体因固体循环而返混”列为气体输运的第 4 条，并写它通常可忽略（Clift und Grace 1985）。式 (3.7)–(3.25) 给出气泡大小、气泡速度、气泡体积份额和乳液空隙率，供式 (3.26) 和上面的 \(u_0-u_{mf}\)、\(\varepsilon\) 使用。它们是物性与气泡关联式，本答卷不逐条重抄，也不引入。

❓ **Q4 上游依赖**：

\(W_1\) 用到 \(f_w\) 或 \(f_{bw}\)、\(\varepsilon_{mf}\)、\(u_0\)、\(u_{mf}\) 和截面积 \(A\)。\(W_2\) 还用到 \(\Omega\)、\(\varepsilon_{bw}\)、\(\varepsilon_d\) 和当地 flotsam 浓度。这些空隙率和速度来自 §3.2 的气泡计算。式 (2.12) 的 \(\dot m^{auf}\)、\(\dot m^{ab}\) 没有被写成 \(W_1\) 或 \(W_2\)。

❓ **Q5 原文空白**：

**结论：Wirsum 有依据，但是两套不同的量。** 单元衡算的上下质量流没有计算式。GR 模型的循环流是全床的体积流，并且按砂 / 污泥拆开。它不是 Hamel 式 (2.6) 里每个 cell、每个粒级的 \(\dot m_{auf}\)、\(\dot m_{ab}\)。

这套循环流不能直接用于惰性床料和焦炭的全体：

- \(W_1\) 的定义句把它限制为 jetsam，也就是 Wirbelsand。向下返回写在 bulk 相的 \((W_1+W_2)c_{j,b}\) 里，用来闭合浓度剖面，没有写成每个 cell 的 \(\dot M_{ab,i}\)。
- \(W_2\) 的定义句把它限制为 flotsam，也就是 Klärschlamm。污泥被明确排除在尾涡之外。焦炭在式 (2.12) 里是活性固体，这里没有被写成 \(W_2\)。
- S. 85 给出的 \(\Omega\) 是别人试验里的数：煤样颗粒 \(1{,}9\,(\mathrm{cm/s})^{1/2}\)（Nienow et al. 1978），湿秸秆 \(3{,}5\,(\mathrm{cm/s})^{1/2}\)（Bilbao et al. 1988）。这是物料常数，不是砂和焦炭的通用常数。
- S. 16 已把细颗粒排除在气泡引起的固体混合之外。即使是活性固体，也不都走 \(W_1\)。

和 Hamel 的对照：

- Hamel S. 12 和 S. 21 点名流动模型确定向上、向下的固体流，没有计算式（T12 问题 2）。Wirsum 的 §3.2 同样没有把这两股 cell 流写出来。
- Hamel 表 2.2、Bild 2.2、式 (2.1)、式 (2.6) 都没有 \(W_1\)、\(W_2\)、flotsam、jetsam。
- Hamel 写明跟随 Wirsum 的位置见文末“与 T12 的出入”。其中没有一处是床区尾涡循环。Hamel §3.2.1 有床面抛出的尾涡份额和 40 % 抛出（T12 所记 \(f_w=0{,}25\)、\(\zeta_w=0{,}4\)），对象是床面，不是床内相邻 cell。

可选解释：

- (a) \(W_1\) 只描述砂在尾涡里的向上体积流，再靠连续性在 bulk 中向下返回。影响：它覆盖不了焦炭，也给不出每个 cell 的两股 kg/s。
- (b) 把 \(W_1\rho\) 当作 ADR 0007 的 \(\dot M_{auf}\)，再令 \(\dot M_{ab}=\dot M_{auf}\)。影响：原文没有这个密度，也没有“两股相等”；\(W_1\) 沿高度是一个常数，Hamel 的气泡是沿高度长大的。
- (c) 焦炭借用 \(W_2\) 和 \(\Omega=1{,}9\,(\mathrm{cm/s})^{1/2}\)。影响：Wirsum 把 \(W_2\) 写成污泥，\(\Omega\) 来自煤样或秸秆试验；这是外推。

待总司令裁定。本卡不采用。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 不能直接用。Hamel 没有引用 Gibilaro–Rowe 模型，也没有引用 \(W_1\)、\(W_2\) 作为床区固体流。式 (5.12)、式 (5.34) 和式 (5.35)–(5.40) 是大颗粒 flotsam（污泥）的关联，任务卡禁止引入。式 (5.30) 经目标所有者删去第二个 \(f_w\) 之后，仍是全床砂的体积流，不是 Hamel 式 (2.6) 的粒级质量流。

本问题引用、未打开的文献：Gibilaro und Rowe（1974），S. 77；Bilbao et al.（1988），S. 81、S. 85；Werther（1976），S. 85；Nienow、Rowe und Chiba（1978），S. 81、S. 85；Rajan und Wen（1980），S. 85；Daw und Frazier（1988），S. 84；Naimer、Chiba und Nienow（1982），S. 81；Masson und Rios（1984），S. 80；Garcia-Ochoa et al.（1989），S. 81；Hilligardt（1986），S. 29；Clift und Grace（1985），S. 27；Davidson und Harrison（1963），S. 28；Wintruff（1996），S. 85、S. 87。

---

## 问题 3：\(M_i\) 与惰性床料

❓ **Q1 出处**：

- 表 2.1：PDF 16 左 / S. 12。总质量是输入。
- 表 2.2：PDF 16 右 / S. 13。流动子模型计算惰性料和活性固体的 holdup。
- 式 (3.26)：PDF 26 右 / S. 33，仍在 §3.2 末、§3.3 标题之前。
- 符号：PDF 7 右 / S. XI，\(m_s\) 为 „Masse des Feststoffinventars“，kg。
- 式 (2.12) 不含惰性料：S. 19。
- 检索例外见 Q5。

❓ **Q2 变量与单位**：

式 (3.26)：

\[
m_s=\int_0^{H_{Ws}}(1-\varepsilon_b)(1-\varepsilon_d)\rho_p A(h)\,\mathrm{d}h
\]

\(m_s\) 的单位是 kg。\(H_{Ws}\) 是膨胀床高，m。\(\varepsilon_b\)、\(\varepsilon_d\) 无量纲。\(\rho_p\) 是一个颗粒密度，kg/m³。\(A(h)\) 是截面积，m²。式子假设全部固体都在床内，只用一个 \(\rho_p\)。原文没有写这个 \(\rho_p\) 是砂、焦炭，还是两者的平均。

❓ **Q3 德文原句**：

- S. 33：„Für die Masse des Feststoffinventars im Reaktor gilt unter der Annahme, daß sich der gesamte Feststoff innerhalb der Wirbelschicht befindet“.
- 同页：给定总质量时，膨胀床高要迭代；„Geeigneter ist die Vorgabe der Höhe der expandierten Wirbelschicht … woraus direkt das Feststoffinventar nach Gleichung (3.26) bestimmt werden kann.“
- 表 2.1 的运行参数含 „Zusammensetzung, Gesamtmasse und Partikelgrößenverteilung des Feststoffinventars“。
- 表 2.2：床区流动子模型的计算结果包括 „Feststoffbewegung, Feststoffholdup und Feststoffausstoß in den Freiraum (Inertmaterial und aktive Feststoffe)“。自由空间里颗粒的非稳态加热注明 „nur Inertmaterial“。

❓ **Q4 上游依赖**：

\(\varepsilon_b\)、\(\varepsilon_d\) 来自 §3.2 的气泡计算（式 (3.7)–(3.25)，不引入）。床高可以是输入，也可以由总质量反推。式 (2.12) 不使用 \(m_s\)。GR 模型的浓度是体积份额，不是 \(m_s\) 的粒级分解。

❓ **Q5 原文空白**：

**床高、空隙率和床料质量：Wirsum 有依据，即式 (3.26)。** 这是 Hamel 没有印出的那条几何关系。Hamel 表 2.2 把床高当作输入，理由是运行中用压差把床高保持恒定；符号表有膨胀床高，但没有 \((1-\varepsilon)\rho A h\)。两边都倾向于给定床高。Wirsum 多写了由床高算出 \(m_s\)，也允许多给定质量再反推床高。

\(m_s\) 含哪些固体，原文没有用一句话写死。

可选解释：

- (a) \(m_s\) 是床内全部固体。依据是“gesamte Feststoff”和表 2.2 把惰性料、活性固体的 holdup 都放在流动子模型里。影响：一个 \(\rho_p\) 必须代表混合物；活性固体的质量被算进这个积分。
- (b) \(m_s\) 主要是床料，活性固体另算。依据是式 (2.12) 只衡算活性固体，以及 S. 84 计算空隙率时忽略少量燃料颗粒。影响：积分不能当作 ADR 0007 的 \(M_i\)，除非另加活性固体的存量。
- (c) 表 2.1 的总质量是独立输入，式 (3.26) 只是由它求床高的校核。影响：与“给定床高更合适”并列，原文两句都在，没有指定工况用哪一句。

**与 Hamel „Abgleich der Inertmaterialbilanz“ 对应的校正：获准范围内未找到。**

检索（例外，只为查明质量的定义、口径和校正；不引入动力学、物性和污泥关联式）：

| 检索词 | 命中 | 读到哪里 | 和校正的关系 |
|---|---|---|---|
| Abgleich | PDF 17 左 / S. 14 | 读完这一句 | „muß zwischen Vorabrechnung und Zellenmodell iterativ ein Abgleich erfolgen.“ 校正的是预计算和单元模型，不是惰性料衡算 |
| Inertmaterial | PDF 16 右 / S. 13 | 表 2.2 该格 | 流动子模型计算惰性料和活性固体的 holdup；自由空间加热只写惰性料。没有校正步骤 |
| Inertmaterial | PDF 78 右 / S. 137 | 这一句及其前一句（燃烧放热的结果叙述） | „durch Umlauf von Inertmaterial in untere Bettregionen gelangen.“ 说的是残焦随惰性料循环到床层下部。没有校正量、目标值或残差。速率和热值未取 |
| Bettmaterial | PDF 25 左 / S. 30 | 式 (3.11) 旁的一句 | \(\psi\) 取决于床料平均粒径和床的尺寸。这是气泡系数的依赖，不是存量校正 |
| Bettmaterial | PDF 29 右 / S. 39 | §3.4.1 的顶盖厚度句 | 顶盖厚度等于床料平均粒径。见问题 4 |
| Bettmaterial | PDF 47（§5.5，约 S. 74–75） | 命中段 | 细床料混合更好，尾涡向上的固体流更大。这是试验叙述，没有质量公式。曲线未取 |
| Bettmaterial | PDF 58（第 6 章） | 命中句 | 活性颗粒以气泡速度向上或以床料平均速度向下。这是单颗粒停留假设，不是 cell 存量校正。传热式未取 |
| Bettmaterial | PDF 60（第 6 章） | 图题 | 干燥时间随床料粒径变化。未取 |
| Bettmaterial | PDF 74 左 / S. 128 | 命中句及下一句 | 污泥灰不在床内积累，随烟气全部带出；石英砂通常留下一年以上。见问题 4。没有排灰流率 |
| Inertanteil、Bettmasse | 文本层无命中 | — | — |

全书检索没有出现 „Inertmaterialbilanz“ 或 „Abgleich der Inertmaterialbilanz“。已查章节是上表加上 §2.2、§2.3、§3.2、§3.4.1、§5.6。

和 Hamel 的对照：Hamel 的这个校正只出现在 Bild 2.2 的预计算框里，S. 12 列举预计算内容时没有写它的方程（T12 问题 1）。Wirsum 的迭代 Abgleich 有对应的模块分工，但对象是温度和子模型，不是惰性料质量。不能把 S. 14 那句当成 Hamel 图框里的校正。

可选解释：

- (a) Wirsum 不用单独的惰性料校正：砂的分布由 GR 剖面给出，总质量由式 (3.26) 或表 2.1 给出。影响：Hamel 图框里的校正仍然没有着落。
- (b) 式 (3.26) 在“给定床高 / 给定质量”之间迭代，就是某种装料校正，只是没有用 Inertmaterial 这个词。影响：校正的是总质量或床高，不是各 cell 的惰性份额；一个 \(\rho_p\) 的口径仍开着。
- (c) 校正在第 6 章以后的装置计算里，而 S. 137 只是结果句。影响：本卡按例外只读了命中句；若要扩大到整章装置结果，须总司令放开范围。本卡记为范围外未核查的是“第 7 章以后有没有另写校正算法”，不是“有没有这个德文专名”——专名在获准检索里没有。

待总司令裁定。本卡不采用。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 式 (3.26) 是 Hamel 没有印出的存量关系。Hamel 没有引用它。表 2.2 把床高定为输入，与 Wirsum“给定床高更合适”同向，但 Hamel 没有授权用一个 \(\rho_p\) 把空隙率乘成质量。能否把式 (3.26) 当作 ADR 0007 所缺的 \(M_i\)，以及 \(\rho_p\) 取砂还是混合物，由 ADR 0009 裁定。本卡不把它写入模型。S. 14 的 Abgleich 不能用来实现 Bild 2.2 的惰性料校正。

---

## 问题 4：两端与外排

❓ **Q1 出处**：

- 床面抛出：§3.4.1，PDF 29 右 / S. 39 至 PDF 30 左 / S. 40，式 (3.56)–(3.60)。
- 单元衡算的两端：式 (2.12) 与 S. 20 的定义。没有最底一层或最上一层的特例。
- 离析剖面的床底：式 (5.22)，PDF 51 右 / S. 83。
- 排灰：式 (2.12) 无对应项。检索命中 PDF 74 左 / S. 128。

❓ **Q2 变量与单位**：

- 一次爆发由顶盖抛出的质量 \(E_{0,r}\)：kg。式 (3.57) 为 \(d_p\rho_p(1-\varepsilon_{mf})\cdot\pi d_b^2/2\)。半球面面积 \(2\pi(d_b/2)^2=\pi d_b^2/2\)，再乘厚度 \(d_p\) 和堆积密度，量纲是 kg。
- 爆发次数 \(\omega_b\)：式 (3.58)，\(\varepsilon_b u_b\cdot 6/(\pi d_b^3)\)，单位 1/(m² s)。
- 面通量 \(\dot E_{0,r}\)、\(\dot E_{0,w}\)：kg/(m² s)。式 (3.59) 等于式 (3.57) 与式 (3.58) 的乘积：

\[
\dot E_{0,r}=3\varepsilon_b u_b\frac{d_p}{d_b}\rho_p(1-\varepsilon_{mf})
\]

\[
\dot E_{0,w}=0{,}4\,\varepsilon_b u_b(1-\varepsilon_d)f_w\rho_p
\]

  0,4 是尾涡质量中被抛出的份额。模型采用式 (3.59)，不采用式 (3.60)。
- 这些量是床面的面通量。式 (2.12) 的 \(\dot m^{auf}\) 是 kg/s。原文没有印出面通量乘面积等于最上一层 cell 的向上流。

❓ **Q3 德文原句**：

- S. 39：抛出发生在 „Wirbelschichtoberfläche“；„Der Großteil des ausgeworfenen Feststoffs fällt nach einer bestimmten Flugzeit in die Wirbelschicht zurück.“
- 同页：颗粒来自 „Blasendach“ 或 „im Nachlauf der Blase emporgetragenen Feststoff“，称为 „Roof-Mechanismus“ 和 „Wake-Mechanismus“。顶盖厚度 „eine dem mittleren Durchmesser des Bettmaterials identische Schichtdicke“，即式 (3.56) \(\delta_e=d_p\)。
- S. 40：„Im vorliegenden Simulationsmodell wird der Feststoffaustrag an der Wirbelschichtoberfläche nach Gleichung (3.59) berechnet.“ 粗颗粒（Geldart B、D）以顶盖为主；Geldart A 或高气速以尾涡为主。本模型仍用式 (3.59)。
- S. 20：向下流是 „aus der Zelle i nach unten“，没有写最底一层改为 0 或改为排灰。
- S. 83：„Verwendung der gemischten Randbedingung am Boden der Wirbelschicht ergibt“式 (5.22)。这是 jetsam 浓度梯度，未知量是 \(c_{j,b}\)。
- S. 128：„reichert sich die mit dem Klärschlamm eingebrachte Asche nicht im Wirbelbett an, sondern wird vollständig durch den Rauchgasstrom ausgetragen. Das Bettmaterial, in der Regel Quarzsand, verbleibt ein Jahr oder mehr im Reaktor“.

❓ **Q4 上游依赖**：

式 (3.59) 用床面的 \(\varepsilon_b\)、\(u_b\)、\(d_b\)、床料粒径和 \(\varepsilon_{mf}\)、\(\rho_p\)。式 (2.12) 的向上流没有被设成这个面通量。\(W_1\) 把砂送到床面，§3.4.1 再单独计算抛出；两节没有写 \(W_1\) 等于 \(\dot E_0 A\)。

❓ **Q5 原文空白**：

**床面抛出：Wirsum 有依据。** 模型用顶盖机制的式 (3.59)。它是床面的面通量，不是最上一层床区 cell 的 \(\dot m^{auf}\)。

和 Hamel 的对照：Hamel §3.2.1 也有顶盖厚度等于床料平均粒径、一次破裂的抛出量和两种机制合成的床面总抛出（T12 所记式 (3.49)、式 (3.50)、式 (3.66)，单位 kg/(m² s)）。Hamel 同样没有写 \(F_0\) 乘面积等于最上层 cell 的向上流（T12 问题 4）。两边相容的是“抛出是床面通量”。Wirsum 选用顶盖式、给出 0,4 的尾涡抛出式但不用，这些选择 Hamel 没有照引。

可选解释：

- (a) 最上层 cell 的向上固体流等于 \(\dot E_{0,r}\) 乘床面积。影响：要补单位换算；床内其它 cell 的向上流仍然没有。
- (b) \(\dot E_{0,r}\) 只进入自由空间，床区 \(\dot m^{auf}\) 是另一股，床面另做匹配。原文没有匹配式。
- (c) “大部分落回”是最上层 cell 的向下入流。原文只说落回床层，没有指定 cell，也没有指定用 \(\dot m^{ab}\) 还是另一项。

**床底 cell 的向下流：获准范围内未找到。** 已查 §2.3.2.2、S. 20 的符号表、Bild 2.3、§5.6 的床底边界。式 (5.22) 约束的是浓度剖面，不是 \(\dot m^{ab}\)。

可选解释与问题 1 的流定义并列，不选取：(a) 下方没有 cell，该股为 0；(b) 该股就是外排；(c) 仍按内部循环计算，边界另给。原文没有这三句。

**排灰流：式 (2.12) 有依据地没有这一项。** 数值和 cell 位置，获准范围内未找到。S. 128 是装置叙述：污泥灰随烟气离开，石英砂留在床内。没有 kg/s，没有“从哪一个 cell 取出”。

和 Hamel 的对照：Hamel 式 (2.6) 有 \(\dot m_{aus}\)，表 2.1 有排灰点的位置和数目，没有流率（T12 问题 3）。Wirsum 的单元固体衡算连符号都没有。不能用 S. 128 填 Hamel 的 \(\dot m_{aus}\)。

可选解释：

- (a) 污泥灰的离开已经含在床面抛出之后的带出里，所以式 (2.12) 不另列。影响：这只覆盖污泥灰；砂不离开。Hamel 的 \(\dot m_{aus}\) 仍要按 Hamel 自己的定义处理。
- (b) 排灰在获准范围以外的装置章里另有流率。S. 128 没有写。若要查第 7 章的输入表，须总司令扩大范围。本卡把“第 7 章有没有排灰 kg/s”记为范围外未核查。

待总司令裁定。本卡不采用。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 式 (3.56)、式 (3.59) 与 Hamel 已有的床面抛出是同一类量。Hamel 已有自己的抛出式，没有写“改用 Wirsum 式 (3.59)”。式 (3.60) Wirsum 自己不用。这两式都不能填床内 \(\dot M_{auf,i}\)。S. 128 不能填 \(\dot m_{aus}\)。

本问题引用、未打开的文献：Pemberton und Davidson（1986），S. 39–40；Pyle（1965），S. 39。

---

## 问题 5：预计算与单元模型的分工

❓ **Q1 出处**：§2.2，PDF 16 左下 / S. 12 起，PDF 17 / S. 14–15；Bild 2.1 在 PDF 17 左 / S. 14。求解过程在 §2.3.3，PDF 20 右 / S. 21 至 PDF 21 / S. 22–23。

❓ **Q2 变量与单位**：牛顿法的未知量是全部衡算变量，个数由式 (2.1) 给出。原文没有为“内层冻结”单列一个符号。

❓ **Q3 德文原句**：

- S. 14–15：流动、离析、干燥与脱挥发分、自由空间与床之间的传热，„lassen sich nur gesamtheitlich und zunächst unabhängig von der abstrahierten Betrachtung des Reaktors … beschreiben。Diese Prozesse werden durch separate Teilmodelle vor der Formulierung von Bilanzgleichungen berechnet。Die Ergebnisse werden anschließend auf das Zellenmodell abgebildet und liefern die zur Bilanzierung notwendigen Stoff- und Energieströme.“
- 同页：因为这些过程要用到温度等量，„muß zwischen Vorabrechnung und Zellenmodell iterativ ein Abgleich erfolgen.“
- Bild 2.1：由内向外是单元模型、预计算、模型配置。
- S. 21：稳态衡算是非线性代数方程组，求解基于 Wozny（1979）的 Newton–Raphson，并对块三对角雅可比做了适配。S. 22–23 是阻尼：新的一步要让 \(\lVert F\rVert_2\) 下降。非相邻 cell 的物流会在雅可比里留下块外元素，例子是循环床和颗粒对流传热。

❓ **Q4 上游依赖**：预计算的输出被说成衡算所需的物流和能流。单元模型的温度再送回预计算。式 (2.12) 是内层所解的固体方程之一。\(W_1\)、\(W_2\) 和 \(m_s\) 写在流动和离析子模型里，不在牛顿校正公式里。

❓ **Q5 原文空白**：

**模块划分：Wirsum 有依据。** 下表只记原文写了的归属。原文没说内层是否冻结的，写未说明。

| 项目 | 原文放在哪里 | 牛顿内层是否不变 |
|---|---|---|
| 气泡、空隙率、床料总质量 \(m_s\) | §3.2；表 2.2 把床区流动列为子模型 | 未说明。S. 15 说结果提供衡算所需的流，又说要和单元模型迭代 |
| 离析剖面、\(W_1\)、\(W_2\) | 预计算所列举的 Segregation；方程在 §5.6 | 未说明 |
| 干燥、脱挥发分及其释放位置 | 预计算；例子写在 S. 15 | 未说明。温度变了就要重做 |
| 每个衡算空间分到的释放份额 | 由子模型结果分配 | 未说明 |
| 式 (2.12) 的各股固体流 | 单元模型的衡算方程；变量分组称为固体流 | 这些方程是内层所解的对象。哪一股是未知量，见问题 1 |
| 化学反应的源项 \(\dot m^{r}\) | 出现在式 (2.12) 里 | 未说明迭代中是否随浓度、温度更新。速率式本卡不引入 |
| 阻尼牛顿法 | §2.3.3 | 内层解的是衡算方程组。每一次迭代是否重算流动，原文未说明 |

和 Hamel 的对照：Hamel S. 12 同样把流动、离析、干燥和热解放进预计算，再把每个 cell 的固体进出流按流动模型分配，并用外层循环使温度、物流和预计算一致（T12 问题 5）。Bild 2.2 另外有 „Abgleich der Inertmaterialbilanz“，Wirsum 的 Abgleich 不是这一项（问题 3）。Hamel 式 (2.1)、式 (2.6) 与 Wirsum 式 (2.1)、式 (2.12) 的分工相同：计数和衡算在单元模型，流量的算法不在这两式里。

可选解释：

- (a) „abgebildet“表示预计算给出 \(W_1\)、\(m_s\) 和释放份额的数值，内层只解式 (2.12)。影响：与 ADR 0007 的“整体流内层冻结”同向，但原文没写冻结，也没写 \(W_1\) 如何变成式 (2.12) 的 \(\dot m^{auf}\)。
- (b) 分配的是连接关系，式 (2.12) 的各股流仍是牛顿未知量。影响：\(W_1\) 不必是内层常数；则 ADR 0007 要冻结的那个整体流仍然没有被 Wirsum 写出来。
- (c) 内层每一步都重算气泡和 \(W_1\)。S. 15 只写了外层的迭代 Abgleich，没有写内层调用。影响：不能从这句推出内层冻结，也不能推出内层重算。

待总司令裁定。本卡不采用。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 阻尼牛顿法 Hamel 已在 S. 20 写明采用，不在本卡的采用范围内重复裁定（见 ADR 0002）。预计算与单元模型的分工两边同向，但不能把 Wirsum 的 \(W_1\)、式 (3.26) 或污泥拟合塞进 Hamel 的预计算。Hamel 没有写“内层冻结”，Wirsum 也没有。

---

## 问题 6：催化剂质量口径

❓ **Q1 出处**：获准的 §2.2、§2.3、§3.2、§3.4.1、§5.6 里没有以催化剂质量为基准的速率。下面是检索例外。

❓ **Q2 变量与单位**：命中段里看到的速率单位是 kmol/(m³ s)，不是 mol/(s·kg)。没有 \(m_{Kat}\) 或 kg_Kat。本卡不抄那些速率式的常数。

❓ **Q3 德文原句**：

- PDF 23 左 / S. 26：„Beispiel für Gruppe A Partikel sind Katalysatoren für katalytische Wirbelschichtprozesse.“ 这是 Geldart 分类的举例。
- PDF 71 右 / S. 123：„Als Katalysatoren sind Koks, Kalziumoxid und Kalziumsulfat wirksam. Die katalytische Wirkung von Quarzsand und Brennstoffaschen wird im Modell von Johnsson nicht berücksichtigt.“
- PDF 72 左 / S. 124：NO 与 CO 的还原写在 „über den Katalysatoren CaO und Koks“ 上。同页速率的单位写成 kmol/m³s。
- PDF 77 右 / S. 135：Philippek 的试验里，床料对 NH₃ 分解和 NO 还原有催化作用。这是试验结论，没有给出用哪一份质量做分母。

❓ **Q4 上游依赖**：这些句子在氮化学和试验讨论里。它们不进入式 (2.12)、式 (3.26) 或 \(W_1\)。

❓ **Q5 原文空白**：

**结论：获准范围内未找到。** 没有找到悬浮相速率以“床料质量 / 惰性料质量 / 灰质量”为基准的写法，也没有 mol/(s·kg_Kat) 或 m³/(s·kg_Kat)。

已查：§2.2、§2.3、§3.2、§3.4.1、§5.6；检索词 Katalysator、Bettmaterial（问题 3 的表）。Katalysator 的命中页是 PDF 23、71、72、77。读到上表所列的句子为止，没有进入速率常数的数值。

和 Hamel 的对照：ADR 0010 第 9 条所缺的是 Corella 行的催化剂质量口径（悬浮相中的全部床料，还是只计惰性床料或灰），单位是 mol/(s·kg_Kat)。Wirsum 这几句既没有这个单位，也没有在“全部床料 / 惰性料 / 灰”之间做选择。Johnsson 方案点名的催化剂是焦炭、CaO 和 CaSO₄，并写明不考虑石英砂和燃料灰；Philippek 的句子又说床料有催化作用。两句对象不同，都不能填 Hamel 的空位。

可选解释：

- (a) 按 Johnsson 方案，催化剂质量只计焦炭和钙化合物，不计砂和灰。影响：这是氮化学文献的取舍，不是 Corella 焦油速率的口径；Hamel 没有引用这句。
- (b) 按 Philippek 的试验句，床料本身就算催化剂。影响：仍然没有 kg 数的定义，不能乘进速率。
- (c) Wirsum 的均相或催化速率以反应器体积为基准，所以不存在 kg_Kat。影响：即便如此，也不能推出 Hamel 表 5.4 的双星号应当改成体积基准。

待总司令裁定。本卡不采用。悬浮相焦油转化仍按 ADR 0010 第 9 条，在口径裁定前不得实现。

❓ **Q6 与归档旧代码的差异**：不适用。

**能否用于 Hamel 主模型：** 不能。Hamel 表 5.4 的 Corella 行没有在这些句子里出现。本卡不把 Johnsson 或 Philippek 的催化说法写入焦油速率。

本问题引用、未打开的文献：Johnsson，S. 123。文本层写 [60]，同页后文另见 [62]；参考文献表里的 [60] 是床层膨胀一文，与氮化学句子对不上，括号号可能认错，本卡不据此打开文献，也不引用年份。Ogada、Philippek 只在 S. 135 的试验句里点名，年份未从正文读出。

---

## 与 T12 的出入

不改 T12 答卷。

T12 问题 2 的可选解释 (d) 写：Hamel 写明采用 Wirsum（1998）的地方是阻尼牛顿法（S. 20）和自由空间阻力式（S. 45）。本卡核对 Hamel 文本层，采用句还有两处，而且都不是床区固体循环：

- PDF 17 / S. 15：Dersch（1994）的 10 个焦炭粒级，以及 Wirsum（1998）把它扩展到任意个粒级，„wird in dieser Arbeit übernommen“。
- PDF 19 / S. 19：能量流的组成写成 „(Wir-sum 1998)“，断词在行末，紧挨式 (2.8)。

这两处 ADR 0002 已经登记。T12 把采用地点收成两处，少了这两处。少掉的两处不提供 \(\dot M_{auf,i}\) 或 \(M_i\)。T12“床区整体上下固体流原文未给出”的结论，本卡没有看到要改的地方。

其余核对：T12 对 Hamel S. 12、表 2.2、Bild 2.2、式 (2.1)、式 (2.6) 的转述，与本卡用来对照的 Hamel 文本层一致。式 (2.6) 含离开体系的 \(\dot m_{aus}\)，T12 问题 3 已这样写。本卡没有发现需要改 T12 的第二条。

---

## 汇总

| 问题 | Wirsum 结论 | Wirsum 页码 | 与 Hamel 是否相容 | 能否用于 Hamel 主模型的依据 |
|---|---|---|---|---|
| 1. 衡算变量是什么 | Wirsum 有依据：每个粒级 1 个变量，分组名称是固体流；只衡算活性固体。哪一股流是那个变量，有歧义 | S. 16 式 (2.1)；S. 18 分组；S. 19–20 式 (2.12) | 计数式与 Hamel 式 (2.1) 同形；固体衡算与 Hamel 式 (2.6) 同类，但 Wirsum 没有 \(\dot m_{aus}\)。Hamel 也没有把未知量印成存量 | 不能替换 Hamel 式 (2.6)，也不能取消 ADR 0007。Hamel 没有引用式 (2.12) |
| 2. 相邻 cell 的上下固体流 | Wirsum 有依据，但是两套量：式 (2.12) 的质量流没有算法；\(W_1\) 只是砂的尾涡体积流，\(W_2\) 只是污泥的体积流。式 (5.30) 的第二个 \(f_w\) 已由目标所有者判为笔误并删去 | S. 81–85；式 (5.30) 在 PDF 52 右 / S. 85 | Hamel S. 12 只点名上下固体流，没有式子。Hamel 没有 \(W_1\)、\(W_2\)。床面的 \(f_w\) 与 40 % 在 Hamel §3.2.1，不在床内 cell | 不能直接用。污泥关联式禁止引入。删去笔误后的 \(W_1\) 仍不是每个 cell 的 kg/s，也不含焦炭 |
| 3. 床高、空隙率、床料质量 | Wirsum 有依据：式 (3.26)。一个 \(\rho_p\) 代表谁，有歧义。总质量也可作为表 2.1 的输入 | S. 33 式 (3.26)；表 2.1 在 S. 12 | Hamel 把床高定为输入，没有这条积分。方向同为“床高可以给定” | Hamel 没有引用式 (3.26)。能否当作 \(M_i\) 由 ADR 0009 裁定。本卡不采用 |
| 3. 惰性料校正 | 获准范围内未找到与 „Abgleich der Inertmaterialbilanz“ 同名的步骤。已查章节和检索词见问题 3 | S. 14 的 Abgleich 是预计算对单元模型 | Hamel 的校正只在 Bild 2.2 的框里，两边都没有方程 | S. 14 那句不能实现 Hamel 的惰性料校正 |
| 4. 床面抛出与最上层向上流 | 抛出有依据：模型用式 (3.59)，单位 kg/(m² s)。它是否等于最上层 \(\dot m^{auf}\)，有歧义 | S. 39–40 | 与 Hamel 的床面通量同类；Hamel 也没有写成 \(F_0\times\) 面积 | 不能填床内 \(\dot M_{auf,i}\)。Hamel 已有自己的抛出式 |
| 4. 床底向下流 | 获准范围内未找到。式 (5.22) 是浓度边界，不是 \(\dot m^{ab}\) | S. 20；S. 83 式 (5.22) | Hamel 同样未给出（T12 问题 4） | 不能从式 (5.22) 推出床底 \(\dot m^{ab}=0\) 或等于排灰 |
| 4. 排灰 | 式 (2.12) 没有这一项，这一点有依据。流率和 cell 位置获准范围内未找到。S. 128：污泥灰随烟气走，砂留下 | S. 19–20；S. 128 | Hamel 式 (2.6) 有 \(\dot m_{aus}\)，但没有流率 | S. 128 不能填 \(\dot m_{aus}\)。第 7 章是否另有流率，范围外未核查 |
| 5. 预计算与单元模型 | Wirsum 有依据。流动、离析、干燥脱挥发分在预计算；衡算在单元模型里由阻尼牛顿法求解。内层是否冻结，有歧义 | S. 14–15；S. 21–23 | 与 Hamel S. 12 的模块划分同向。Hamel 多一个图框里的惰性料校正 | 分工不能当成采用 \(W_1\) 或式 (3.26) 的理由。阻尼牛顿法 Hamel 已另行采用 |
| 6. 催化剂质量 | 获准范围内未找到 kg_Kat 口径。检索见问题 6 | S. 26、S. 123–124、S. 135 | 填不上 ADR 0010 第 9 条的 Corella 口径 | 不能用于悬浮相焦油转化 |

本卡没有选取上表任何一行里的解释，也没有把任何一式写入 Hamel 主模型。
