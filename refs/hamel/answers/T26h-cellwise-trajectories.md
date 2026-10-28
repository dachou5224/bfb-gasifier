# T26h 原文答卷：自由空间颗粒轨迹的逐格解析解（式 3.82–3.105）

| 项 | 内容 |
|---|---|
| 任务卡 | `docs/tasks/T26h-cellwise-trajectories.md` |
| 执行者 | cursor |
| 会话 | 新会话、新 worktree `../bfb-T26h`（本任务不属独立转录，无独立性要求） |
| 日期 | 2026-10-28 |
| 对照 | Hamel 页面图像 `refs/hamel/pages/032.png`–`034.png`（300 dpi）；OCR 文本只用于定位 |
| 相关 ADR | ADR 0016 第 7、11 条，ADR 0026 第 2、3 项，ADR 0028，ADR 0034，ADR 0036 第 1 项，ADR 0047（本卡的裁定） |

**读法说明。** 下列每节都按六个问题作答。Q3 只引关键的一句。凡原文没有说明的地方写在 Q5，并列出可选解释与影响；本卡采用的那一支注明 ADR 0047 的条号，不在这里另作裁定。式 (3.96)、(3.98)、(3.99)、(3.103) 的印刷面照录，并在 Q5 记下与运动方程不相容的核验结论（ADR 0016 第 11 条已裁定不用）。

**范围。** 只核式 (3.82)–(3.105) 与 S. 45–49 的正文；不核球形度分布的取法（ADR 0047 明确不在本卡内），不核 Tabelle 3.1 的五条阻力式（它们不进入轨迹方程）。

## 覆盖

| 式号／位置 | 印刷页 | PDF 页 | 本文件 |
|---|---|---|---|
| (3.82)、(3.83) | S. 44 | 32 左 | 有 |
| (3.84)–(3.87) | S. 45 | 32 右 | 有（只为确认进入轨迹的是 (3.85)、(3.86)） |
| (3.88)–(3.91) | S. 46 | 33 左 | 有 |
| (3.92)、(3.93) | S. 47 | 33 右 | 有 |
| (3.94)–(3.98) | S. 48 | 34 左 | 有 |
| (3.99)–(3.105) | S. 49 | 34 右 | 有 |
| S. 48 的 Zellenmodell 两句 | S. 48 | 34 左 | 有 |
| S. 49 末段的 50 倍一句 | S. 49 | 34 右 | 有 |

---

## Eq. 3.82（颗粒运动方程）

❓ **Q1 出处**：PDF 32 左半页／S. 44，§3.2.4，式 (3.82)。
❓ **Q2 变量与单位**：u_P — 颗粒速度，m/s；t — 时间，s；阻力项按单位质量力写成 m/s²。原文未给单位表，此处按 SI 读。
❓ **Q3 德文原句**：„Die Bewegungsgleichung des Einzelpartikels lautet“（式前一句）。
❓ **Q4 上游依赖**：$u_r$ 见式 (3.83)；$c_W$ 见式 (3.85)；$z$ 见式 (3.86)。$u_g$ 取什么，正文没有在这一页写。
❓ **Q5 原文空白**：$u_g$ 是否等于 §3.2.3 的鬼泡速度 $u_{gb}(h)$，本节未说明（ADR 0016 第 7 条已裁定取 $u_{gb}(h)$）。原文也没有提到“格内取哪个气速”：S. 48 只说逐格应用，没说格内用平均、格中点还是入口值。可选解释与影响见 ADR 0047 第 1 项（本卡取格内高度平均）。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
\frac{du_P}{dt} = -\frac{3}{4}\cdot\frac{c_W}{d_P}\cdot\frac{\rho_g}{\rho_P}\cdot u_r\,|u_r| - \frac{\rho_P-\rho_g}{\rho_P}\cdot g
$$

## Eq. 3.83（相对速度）

❓ **Q1 出处**：PDF 32 左半页／S. 44，§3.2.4，式 (3.83)。
❓ **Q2 变量与单位**：$u_r = u_P - u_g$，m/s，以颗粒为参照的正方向向上。
❓ **Q3 德文原句**：„mit der Relativgeschwindigkeit“。
❓ **Q4 上游依赖**：式 (3.82)。
❓ **Q5 原文空白**：无。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
u_r = u_P - u_g
$$

## Eq. 3.84–3.87（阻力系数与 $z$）

❓ **Q1 出处**：PDF 32 右半页／S. 45，§3.2.4，式 (3.84)–(3.87)。
❓ **Q2 变量与单位**：$Re = |u_r|\,d_P/\nu_g$；$\nu_g$ 为气体运动黏度，m²/s；$c_W$ 无量纲；$z$ 为 $c_W$ 在高 $Re$ 下的渐近值，无量纲。
❓ **Q3 德文原句**：„Dabei ist $z$ der Wert von $c_W$ für sehr große Reynolds-Zahlen $Re_\infty$“（S. 45）。
❓ **Q4 上游依赖**：式 (3.84)（Haider–Levenspiel）给 $z$；式 (3.85) $c_W = 24/Re + z$ 进入式 (3.82)。
❓ **Q5 原文空白**：$Re_\infty$ 没有数值，式 (3.86) 仍保留 $24/Re_\infty$ 项。ADR 0016 第 8 条已裁定取 $Re_\infty\to\infty$，$z = 73{,}69\,\exp(-5{,}0784\,\varphi_S)$。本卡不改这一条。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
z = c_{W}(Re_\infty) \;\xrightarrow{Re_\infty\to\infty}\; 73{,}69\,\exp(-5{,}0784\,\varphi_S)
$$

## Eq. 3.88（运动方程化为 $u_P$ 的二次多项式）

❓ **Q1 出处**：PDF 33 左半页／S. 46，§3.2.4，式 (3.88)。
❓ **Q2 变量与单位**：$A$ — s/m；$B$ — 1/s；$C$ — m/s²（由下三式的量纲读出）。$\Delta$ 见式 (3.92)。
❓ **Q3 德文原句**：„Mit Gleichung (3.85) lässt sich die Differentialgleichung für die Partikelbewegung (3.82) in folgender Form schreiben“。
❓ **Q4 上游依赖**：式 (3.82)、(3.85)。
❓ **Q5 原文空白**：无。符号约定在式 (3.91) 之后的一段，见本答卷“符号约定”一节。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
\frac{du_P}{dt} = A\cdot u_P^{2} + B\cdot u_P + C
$$

## Eq. 3.89–3.91（系数 $A$、$B$、$C$）

❓ **Q1 出处**：PDF 33 左半页／S. 46，§3.2.4，式 (3.89)–(3.91)，紧接 „mit den Koeffizienten“。
❓ **Q2 变量与单位**：$\mu_g$ — 气体动力黏度，Pa·s；$\rho_g$、$\rho_P$ — kg/m³；$d_P$ — m；$g$ — 9,81 m/s²；$z$ 无量纲。
❓ **Q3 德文原句**：„mit den Koeffizienten“。
❓ **Q4 上游依赖**：式 (3.85) 的 $z$ 项与 $u_r$ 的展开，式 (3.86) 的 $z$。
❓ **Q5 原文空白**：无。括号内外两套符号的适用条件见下一节。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
\begin{aligned}
A &= -\;(+)\;\frac{3}{4}\cdot z\cdot\frac{\rho_g}{\rho_P\,d_P} \\
B &= -\frac{18\,\mu_g}{\rho_P\,d_P^{2}} \;+\; (-)\;\frac{3}{2}\cdot z\cdot\frac{\rho_g}{\rho_P\,d_P}\cdot u_g \\
C &= \frac{18\,\mu_g\,u_g}{\rho_P\,d_P^{2}} \;-\;(+)\;\frac{3}{4}\cdot z\cdot\frac{\rho_g}{\rho_P\,d_P}\cdot u_g^{2}
     \;-\;\frac{\rho_P-\rho_g}{\rho_P}\cdot g
\end{aligned}
$$

括号外：$A = -\tfrac34 z\rho_g/(\rho_P d_P)$，$B$ 的 $z$ 项取 $+\tfrac32$，$C$ 的 $z$ 项取 $-\tfrac34$。括号内：$A = +\tfrac34$、$B$ 的 $z$ 项取 $-\tfrac32$、$C$ 的 $z$ 项取 $+\tfrac34$，其余不变。

核对（把式 (3.82) 的阻力按 $u_r = u_P-u_g$ 展开）：颗粒快于气体时 $-q\,(u_P-u_g)|u_P-u_g| = -q\,u_P^2 + 2q\,u_g\,u_P - q\,u_g^2$，其中 $q = \tfrac34 z\rho_g/(\rho_P d_P)$，故 $A = -q$、$B$ 的 $z$ 项为 $+2q = +\tfrac32 z\rho_g/(\rho_P d_P)$、$C$ 的 $z$ 项为 $-q = -\tfrac34$；括号内（颗粒慢于气体）三项依次换号。与式 (3.89)–(3.91) 的括号位置一致。

## 符号约定（PDF 33 右半页／S. 47，式 3.91 之后、式 3.92 之前）

❓ **Q1 出处**：PDF 33 右半页／S. 47。该段横跨 S. 46 与 S. 47 的接排：式 (3.89)–(3.91) 在 S. 46 左栏，段落在 S. 47 上、式 (3.92) 之前。
❓ **Q2 变量与单位**：不涉及新量。
❓ **Q3 德文原句**：„Bewegt sich ein Partikel langsamer als das ihn umgebende Gas oder fällt ein Partikel entgegen dem Gasstrom, so sind zur Berechnung der Koeffizienten $A$, $B$ und $C$ die in Klammern stehenden Vorzeichen gültig. Ist die Partikelgeschwindigkeit größer als die Gasgeschwindigkeit gelten entsprechend die nicht in Klammern stehenden Vorzeichen.“
❓ **Q4 上游依赖**：式 (3.89)–(3.91)。
❓ **Q5 原文空白**：原文以 $u_P$ 与 $u_g$ 的比较分开两支，但**没有说 $u_g$ 在格内取哪个值**（见式 3.82 的 Q5）。逐格应用时，两支的分界就是该格的 $\bar u_{g,i}$；ADR 0047 第 2 项据此定“$u_P < \bar u_{g,i}$ 用括号内，$u_P > \bar u_{g,i}$ 用括号外”。两支在 $u_P = u_g$ 处给出同一个加速度 $-(\rho_P-\rho_g)g/\rho_P$，因此分界处连续，不需要额外的切换规则。
❓ **Q6 与归档旧代码的差异**：未查看归档。

## Eq. 3.92（判别式）

❓ **Q1 出处**：PDF 33 右半页／S. 47，§3.2.4，式 (3.92)。
❓ **Q2 变量与单位**：$\Delta$ 的量纲与 $B^2$ 同，1/s²。
❓ **Q3 德文原句**：„In Abhängigkeit des Ausdrucks … nimmt die Lösung der Differentialgleichung (3.88) zwei unterschiedliche Formen an (Bronstein 1989).“
❓ **Q4 上游依赖**：式 (3.89)–(3.91)。
❓ **Q5 原文空白**：$\Delta = 0$ 的情形没有给出解。ADR 0047 第 2 项取 $\Delta\to0$ 的极限（$|\Delta| \le 10^{-10} B^2$）作为第三种闭式。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
\Delta = 4AC - B^{2}
$$

## Eq. 3.93–3.95（Fall 1，$\Delta > 0$）

❓ **Q1 出处**：式 (3.93) 在 PDF 33 右半页／S. 47；式 (3.94)、(3.95) 在 PDF 34 左半页／S. 48。
❓ **Q2 变量与单位**：$u_{P,0}$ — 该段起始速度，m/s；$h_{P,0}$ — 该段起始高度（床面以上），m；$t_0$ — 该段起始时刻，s；$\omega$ — 相位，rad。
❓ **Q3 德文原句**：„Für den Fall 1: $\Delta > 0$ ergibt die Integration der Differentialgleichung (3.88) für die Partikelgeschwindigkeit als Funktion der Flugzeit“（S. 47）；„Mit Hilfe der Transformation $u_P\,dt = dh$ liefert die weitere Integration von Gleichung (3.93) die Höhe des Partikels über der Wirbelschichtoberfläche als Funktion der Zeit.“（S. 48）。
❓ **Q4 上游依赖**：式 (3.88)、(3.92)、(3.94)。
❓ **Q5 原文空白**：无。两式与式 (3.88) 相容（代入核验见 ADR 0016 第 11 条与 ADR 0047 的依据）；本卡照用。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
\begin{aligned}
u_P(t) &= \frac{\sqrt{\Delta}}{2A}\tan\!\left[\frac{1}{2}\omega\sqrt{\Delta} + \frac{1}{2}(t-t_0)\sqrt{\Delta}\right] - \frac{B}{2A} \\
\omega &= \frac{2}{\sqrt{\Delta}}\arctan\!\left[\frac{2A\,u_{P,0}+B}{\sqrt{\Delta}}\right] \\
h_P(t) &= h_{P,0} - \frac{1}{A}\ln\!\left[\cos\!\left(-\frac{1}{2}\omega\sqrt{\Delta} - \frac{1}{2}(t-t_0)\sqrt{\Delta}\right)\right]
        - \frac{B}{2A}(t-t_0) + \frac{1}{A}\ln\!\left[\cos\!\left(\frac{1}{2}\omega\sqrt{\Delta}\right)\right]
\end{aligned}
$$

**与本卡的关系。** 以 $\theta_0 = \tfrac12\omega\sqrt{\Delta}$ 代入，即 ADR 0047 第 2 项的第一条闭式。$\theta_0 = \arctan[(2Au_0+B)/\sqrt{\Delta}]$ 与式 (3.94) 的 $\omega$ 等价；由于 $\cos(-\theta)=\cos\theta$，式 (3.95) 的四个余弦项合并为 $\ln[\cos\theta(\tau)/\cos\theta_0]$，即 ADR 的写法。

## Eq. 3.96（Fall 1 的 $h_P(u_P)$）

❓ **Q1 出处**：PDF 34 左半页／S. 48，式 (3.96)。
❓ **Q2 变量与单位**：$h_P(u_P)$ — m；其余同式 (3.93)–(3.95)。
❓ **Q3 德文原句**：„Der Zusammenhang zwischen Höhe des Partikels und seiner Geschwindigkeit lautet“。
❓ **Q4 上游依赖**：式 (3.88) 在 Fall 1 下对高度积分。
❓ **Q5 原文空白**：无空白；但该式的对数项系数印为 $-1/(2A)$，代入核验不满足 $dh/du_P = u_P/[Au_P^2+Bu_P+C]$（ADR 0016 第 11 条，T17 的代入核验）。**本卡不用印刷式 (3.96)**，$h(u)$ 由本答卷自行推导。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
h_P(u_P) = h_{P,0} - \frac{1}{2A}\ln\!\left[\frac{Au_P^{2}+Bu_P+C}{Au_{P,0}^{2}+Bu_{P,0}+C}\right]
  - \frac{B}{A\sqrt{\Delta}}\arctan\!\left[\frac{2Au_P+B}{\sqrt{\Delta}}\right]
  + \frac{B}{A\sqrt{\Delta}}\arctan\!\left[\frac{2Au_{P,0}+B}{\sqrt{\Delta}}\right]
$$

**推导的正确形式（本卡采用）。** 由 $dh/du = u/Q$、$Q = Au^2+Bu+C$、$\Psi = 2Au+B$：

$$
h(u) = h_0 + \frac{1}{2A}\ln\left|\frac{Q}{Q_0}\right| - \frac{B}{2A}\int_{u_0}^{u}\frac{\mathrm du}{Q},
\qquad
\int\frac{\mathrm du}{Q} = \frac{2}{\sqrt{\Delta}}\left[\arctan\frac{\Psi}{\sqrt{\Delta}} - \arctan\frac{\Psi_0}{\sqrt{\Delta}}\right]\ (\Delta>0),
$$

即对数项取 $+\tfrac{1}{2A}$，而不是印刷的 $-\tfrac{1}{2A}$。核验：$\tau=0$ 时 $h=h_0$；对 $u$ 求导回到 $u/Q$。

## Eq. 3.97（$d_P^{*}$）

❓ **Q1 出处**：PDF 34 左半页／S. 48，式 (3.97)。
❓ **Q2 变量与单位**：$d_P^{*}$ — m；分子 $108\,\mu_g^2$；分母 $z\,g\,(\rho_P\rho_g-\rho_g^2)$。
❓ **Q3 德文原句**：„Um die Bedingung $\Delta > 0$ zu erfüllen, muss die Geschwindigkeit des Partikels größer sein als die Gasgeschwindigkeit und der Partikeldurchmesser muss größer sein als der Partikelgrenzdurchmesser.“
❓ **Q4 上游依赖**：$\Delta = 4AC-B^2$ 取颗粒快于气体的那一组符号；此时 $\Delta = 4qw - s^2$，与 $u_g$ 无关，故界限可写成只含物性与 $z$ 的粒径（$q = \tfrac34 z\rho_g/(\rho_P d_P)$，$s = 18\mu_g/(\rho_P d_P^2)$，$w = (\rho_P-\rho_g)g/\rho_P$）。
❓ **Q5 原文空白**：原文没有把 $d_P^{*}$ 定义成扬析切割粒径；它只是 Fall 1 的判别界限（ADR 0016 第 11 条已记）。本卡只用它选择分支，不作粒级切割。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
d_P^{*} = \left(\frac{108\,\mu_g^{2}}{z\,g\,(\rho_P\rho_g-\rho_g^{2})}\right)^{1/3}
$$

## S. 48 的 Zellenmodell 两句（逐格应用的依据）

❓ **Q1 出处**：PDF 34 左半页／S. 48，式 (3.97) 与 Fall 2 之间。
❓ **Q2 变量与单位**：不涉及新量。„Bilanzraum“ 即衡算空间（cell）。
❓ **Q3 德文原句**：„Innerhalb des Zellenmodells ist es zweckmäßig zunächst die Geschwindigkeit des Partikels an den Grenzen des betrachteten Bilanzraumes zu bestimmen, da über die Partikelgeschwindigkeit entschieden wird, ob die Lösung der Differentialgleichung für Fall 1 mit $\Delta > 0$ oder für Fall 2 mit $\Delta < 0$ angewendet werden muss.“ 以及 „Bewegt sich ein Partikel langsamer als das ihn umgebende Gas, fällt ein Partikel zurück oder ist der Partikeldurchmesser kleiner als der Partikelgrenzdurchmesser, nimmt der Faktor $\Delta$ einen Wert kleiner Null an.“
❓ **Q4 上游依赖**：式 (3.92)、(3.97)。
❓ **Q5 原文空白（本卡的核心空白）**：
1. 原文只说“先定边界上的速度以选分支”，**没有说格内气速取常数、更没有说取哪个代表值**。可选解释：(a) 格内取 $u_{gb}(h)$ 的高度平均（ADR 0047 第 1 项，本卡采用）；(b) 格中点值；(c) 入口边界值。影响：最下格（含锥段）(a) 与 (b) 相差 0,0065 m/s（0,53 %），其上各格不超过约 0,0008 m/s；最下格直接改变第一个界面夹住窗口的上界（ADR 0047 第 1 项的数值）。
2. 原文没有给“到达格界、$u_P$ 穿过 $u_g$、转向”这三类事件的求法，也没有说跨格时用哪个格的物性。ADR 0047 第 2 项裁定：格内按闭式推进到最早事件，到达格界时以界面上的速度为初值进入相邻格，用该格的系数。
3. 原文没有给格内 $\Delta \approx 0$ 的处理。ADR 0047 第 2 项取极限（$|\Delta|\le10^{-10}B^2$）。
❓ **Q6 与归档旧代码的差异**：未查看归档。

## Eq. 3.98–3.99（Fall 2，$\Delta < 0$ 的 $u_P(t)$）

❓ **Q1 出处**：式 (3.98) 在 PDF 34 左半页／S. 48；式 (3.99) 在 PDF 34 右半页／S. 49。
❓ **Q2 变量与单位**：$\eta$ — 由 $u_{P,0}$ 定的积分常数，s；$\sqrt{-\Delta}$ — 1/s。
❓ **Q3 德文原句**：„Fall 2: $\Delta < 0$ … Die Geschwindigkeit des Partikels als Funktion der Zeit ist durch … gegeben“。
❓ **Q4 上游依赖**：式 (3.88)、(3.99)。
❓ **Q5 原文空白**：无空白；但印刷式代入不满足初条件：$t = t_0$ 时式 (3.98) 不返回 $u_{P,0}$（研究底稿第 2 节的算例：$0{,}5\to-0{,}53$、$-0{,}5\to0$；ADR 0016 第 11 条）。**本卡不用印刷式 (3.98)、(3.99)**，改用本答卷推导的闭式。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
u_P(t) = \frac{(\sqrt{-\Delta}-B)\exp[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}] - \sqrt{-\Delta} - B}
              {2A\exp[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}] - 2A},
\qquad
\eta = \frac{1}{\sqrt{-\Delta}}\ln\!\left(\frac{\sqrt{-\Delta}+(2Au_{P,0}+B)}{\sqrt{-\Delta}-(2Au_{P,0}+B)}\right)
$$

**本卡采用的等价写法。** $s = \sqrt{-\Delta}$，$\Psi_0 = 2Au_0+B$，

$$
E(\tau) = \frac{\Psi_0+s}{\Psi_0-s}\,e^{-s\tau},\qquad
u(\tau) = \frac{s\,(E+1)/(E-1)-B}{2A},\qquad
h(\tau) = h_0 + \frac{(s-B)\tau - 2\ln\left|\dfrac{E_0-e^{s\tau}}{E_0-1}\right|}{2A},
$$

$E_0 = E(0)$。核验：$u(0) = [s(E_0+1)/(E_0-1)-B]/(2A) = u_0$；$h(0)=h_0$；对 $\tau$ 求导回到式 (3.88)（研究底稿与 T26h 的测试各核一次）。转写自式 (3.98)、(3.99) 的指数形式：$\exp[\eta s-\tau s] = E$，因为 $\eta s = \ln[(\Psi_0+s)/(\Psi_0-s)]$。

## Eq. 3.100–3.102（Fall 2 的 $h_P(t)$）

❓ **Q1 出处**：PDF 34 右半页／S. 49，式 (3.100)–(3.102)。
❓ **Q2 变量与单位**：$\gamma$ 无量纲；$\zeta$ 含 $1/A$，量纲 m；其余同前。
❓ **Q3 德文原句**：„Für die Höhe des Partikels über der Wirbelbettoberfläche als Funktion der Zeit gilt“。
❓ **Q4 上游依赖**：式 (3.98)、(3.101)、(3.102)。
❓ **Q5 原文空白**：无空白。该组式与本卡采用的 $h(\tau)$ 相容（把 $\gamma = E$ 代入，$\zeta$ 即 $\tau=0$ 的积分常数）；本卡直接用上节写出的合并形式，不另置 $\gamma$、$\zeta$。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
h_P(t) = h_{P,0} - \frac{1}{A}\ln[\gamma+1] - \frac{-0{,}5\sqrt{-\Delta}-0{,}5B}{A}\left(\eta-(t-t_0)\right) - \zeta,
\quad
\gamma = e^{\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}}
$$

## Eq. 3.103–3.105（Fall 2 的 $h_P(u_P)$ 与 $\Psi$ 的定义）

❓ **Q1 出处**：PDF 34 右半页／S. 49，式 (3.103)–(3.105)。
❓ **Q2 变量与单位**：$\Psi$、$\Psi_0$ 与 $2Au+B$ 同量纲（1/s）；$h_P(u_P)$ — m。
❓ **Q3 德文原句**：„Die Gleichung der Partikelhöhe über der Wirbelschichtoberfläche in Abhängigkeit der Geschwindigkeit lautet“。
❓ **Q4 上游依赖**：式 (3.98)、(3.104)、(3.105)。
❓ **Q5 原文空白**：无空白；但式 (3.103) 的分子分母配对与 $\Psi$、$\Psi_0$ 的下标相互矛盾：式 (3.104) 给 $\Psi = 2Au_{P,0}+B$（初速），式 (3.105) 给 $\Psi_0 = 2Au_P+B$（当前速度），与“下标 0 表示初值”的通常写法相反；按此代入，对数项不满足 $\int u\,\mathrm du/Q$（ADR 0016 第 11 条）。**本卡不用印刷式 (3.103)–(3.105)**，改用下式。
❓ **Q6 与归档旧代码的差异**：未查看归档。

$$
h(u) = h_0 + \frac{1}{2A}\ln\left|\frac{Q}{Q_0}\right| - \frac{B}{2A\sqrt{-\Delta}}
\ln\left|\frac{(\Psi-s)(\Psi_0+s)}{(\Psi+s)(\Psi_0-s)}\right|,\qquad \Psi = 2Au+B,\ \Psi_0 = 2Au_0+B
$$

## S. 49 末段（解析法的计算量）

❓ **Q1 出处**：PDF 34 右半页／S. 49，§3.2.4 末段。
❓ **Q2 变量与单位**：不涉及新量。
❓ **Q3 德文原句**：„… erlaubt eine um den Faktor 50 schnellere Berechnung der Partikelflugbahnen als die Lösung durch numerische Integration mit Hilfe eines Runge-Kutta-Verfahrens.“
❓ **Q4 上游依赖**：式 (3.93)–(3.105)。
❓ **Q5 原文空白**：原文没有给 50 倍的适用条件（单元数、粒径数、球形度类数）。本卡只把它当作量级参考，不当作验收门；实测耗时写入 `reports/t26h-trajectories/`。
❓ **Q6 与归档旧代码的差异**：未查看归档。

---

## 原文空白汇总

| # | 空白 | 可选解释 | 本卡采用 |
|---|---|---|---|
| 1 | 格内 $u_g$ 的取法 | 高度平均／格中点／入口边界值 | 高度平均（ADR 0047 第 1 项，`解释性选择：ADR 0047`） |
| 2 | 格内事件与跨格物性 | 二分求根／解析反解；用新格物性 | 闭式推进 + 对 $h(\tau)$ 二分（$10^{-12}$ m）+ 新格物性（ADR 0047 第 2 项，`求解手段：ADR 0047`） |
| 3 | $\Delta \approx 0$ | 无印刷解 | 取极限，$\lvert\Delta\rvert\le10^{-10}B^2$（ADR 0047 第 2 项） |
| 4 | 稳定点与捕获区 | 原文无悬停 | 只在“下格 $\bar u_g>v_t>$ 上格 $\bar u_g$”的界面；R* 为两格之并（ADR 0047 第 3 项，细化 ADR 0026 第 2 项） |
| 5 | 捕获判据的 $\delta_u$、$\delta_h$ 与越界距离 | 原文无 | $\lvert u_P\rvert\le\delta_u$ 且按闭式 $h(u)$ 的越界距离 $\le\delta_h$；捕获时刻按最后一次进入 R*（ADR 0047 第 3 项） |
| 6 | 中性格 | 原文无 | $\lvert\bar u_g-v_t\rvert\le10^{-12}v_t$ 报告并停止（ADR 0047 第 3 项） |
| 7 | 格内根、势能井、逐根 $\delta_h$、积分步长 | ADR 0034／0036 的旧机制 | 撤销（ADR 0047 第 5 项） |
| 8 | 打印式 (3.96)、(3.98)、(3.99)、(3.103) | 印刷面照录 | 不用；改用本答卷推导并核验的闭式（ADR 0016 第 11 条、ADR 0047 第 2 项） |

## 核验记录（本卡实际做的对照）

1. 式 (3.88)–(3.91) 的括号内外两套符号，在 $u_P = u_g$ 处给出同一个加速度 $-(\rho_P-\rho_g)g/\rho_P$；逐格分界连续。
2. 慢分支（$u_P<u_g$）的 $\Delta = -(s^2+4qw) < 0$ 恒成立；快分支的 $\Delta = 4qw-s^2$ 与 $u_g$ 无关，$d_P^{*}$ 即 $\Delta=0$。
3. 式 (3.93)–(3.95) 与 ADR 0047 第 2 项的第一条闭式等价（$\theta_0 = \tfrac12\omega\sqrt\Delta$）。
4. 自推的 $\Delta<0$ 与 $\Delta=0$ 闭式对高精度 RK4 的相对差 $<10^{-8}$（`tests/test_freeboard.py`；同一逐格常系数）。
5. 式 (3.96)、(3.98)、(3.99)、(3.103) 的代入不相容，与 ADR 0016 第 11 条、研究底稿第 2 节一致。
