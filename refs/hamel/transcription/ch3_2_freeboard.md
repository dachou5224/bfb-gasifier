全新会话，转录冻结前未读旧材料

# §3.2 自由空间流动力学独立转录

- 任务：T17
- 执行者：cursor
- 日期：2026-10-17
- 范围：§3.2.1–§3.2.4，印刷页 S. 36–S. 49。自 S. 36 的式 (3.49)、(3.50) 起，至式 (3.105)。
- 页面：`refs/hamel/pages/` 的 PDF 028–034。公式按页面图像转录，不用 OCR 文本代替公式。
- 本文件各式的 Q6 在冻结时一律写作「比对阶段填写」。

## 页码

| 印刷页 | PDF | 栏 |
|---|---|---|
| S. 36 | 028 | 左 |
| S. 37 | 028 | 右 |
| S. 38 | 029 | 左 |
| S. 39 | 029 | 右 |
| S. 40 | 030 | 左 |
| S. 41 | 030 | 右 |
| S. 42 | 031 | 左 |
| S. 43 | 031 | 右 |
| S. 44 | 032 | 左 |
| S. 45 | 032 | 右 |
| S. 46 | 033 | 左 |
| S. 47 | 033 | 右 |
| S. 48 | 034 | 左 |
| S. 49 | 034 | 右 |

`refs/hamel/index.json` 的式号表每个编号只保留一条。式 (3.49)、式 (3.50) 在 S. 35 与 S. 36 各出现一次，索引保留的是 S. 35。本答卷引用 S. 36 的这两式时必须带页码，不与 S. 35 混用。索引里 §3.2 的式号从 (3.51) 到 (3.105)，本答卷覆盖这些式号，并加上 S. 36 的 (3.49)、(3.50)。

## 符号与单位

颗粒密度、粒径、体积用下标 P：$\rho_P$、$d_P$、$\bar d_P$、$V_P$。速度下标与 P 同形，本答卷写作 $u_P$。气体下标 g，气泡下标 b。德文小数逗号照录。

本节绝大多数符号没有写出单位。Bild 3.7 写出 $\rho_P = 2200\,\mathrm{kg/m^3}$，横坐标 $u_{b,WS}$ 的单位是 m/s。角度 30°、48° 是正文写出的度。除此以外，本节未给出单位，转录不补单位，也不做 SI 换算。

$\varepsilon_b$、$u_b$、$d_b$、$\varepsilon_{mf}$、$\varepsilon_d$、$u_0$、$u_{mf}$、$h_{WS}$、$A_0$ 在本节直接使用。它们的定义在第 3.1 节。冻结前未打开 §3.1 转录，这里不重写那些式。

## 采用与并列

下面这张表只区分正文自己的措辞。式号后的细节见各式 Q3、Q5。

**Hamel 写明采用**

- 式 (3.72)：„Nach dem Vorschlag von Levy (1986) wird … berechnet“，用于模拟模型的竖直抛出速度分布。
- 式 (3.73)：„wird mit $\sigma = 0{,}6$ angenommen“，平均初速 „in Anlehnung an Martens (1984)“。
- 式 (3.74)–式 (3.76)：鬼泡速度。式 (3.76) 的句子是 „In Anlehnung an Hamdullahpur und MacKay (1986) und Sciazko et al. (1988) wird … angenommen“。
- 式 (3.85)–式 (3.105)：„Daher wird hier dem Vorschlag von Wirsum (1998) gefolgt“。式 (3.86) 用 Haider 与 Levenspiel 的式 (3.84) 计算 $z$。积分形式注明 Bronstein (1989) 的两类解。S. 49 称该解析步骤来自 Wirsum (1998)。
- 缺球形度数据时，$\bar\varphi_S = 0{,}8$ „zugrundegelegt“。轨迹计算假定 $0{,}3 \le \varphi_S \le 1{,}0$。

**正文写出的模型式（不是另一篇文献的并列公式）**

- 式 (3.66)、式 (3.67)：屋顶与尾流合成的总抛出。式 (3.67) 前有 „Modellvorstellung“ 与 „unter Berücksichtigung der Wahrscheinlichkeit einer Einzelblaseneruption“。式 (3.67) 与所引式号不完全一致，见该式 Q5。
- 式 (3.77)–式 (3.83)：单颗粒力平衡与由此得到的运动方程。

**文献并列（正文介绍，未写成模型采用）**

- 式 (3.49)–式 (3.65) 的屋顶、尾流、泡间距、聚并概率，以及 Pemberton 与 Davidson、$f_w = 0{,}25$、$\zeta_w = 0{,}4$、Briens 的 $\sigma = 0{,}4$、Clift 与 Grace 的 30°、Briens 的 48°、Milioli 与 Foster 的式 (3.65)、$n_B$ 可限制为 5。
- 式 (3.68)–式 (3.71)，Hatano 与 Ishida 的速度区间，Bild 3.7 上的各条曲线。式 (3.68) 被指为过高。
- §3.2.3 的自由射流类比（Zenz 与 Weil；George 与 Grace；Horio 等）。这一支没有式号。
- Al Taweel 等 (1989)：Fung 与 Hamdullahpur 用来算颗粒团阻力。正文没有说本模型采用它。
- Tabelle 3.1 的五条球形阻力式。正文把它们放在 „in der Literatur vorgestellte“ 的位置。轨迹阻力用的是式 (3.85)。
- 式 (3.84) 同时是式 (3.86) 的输入，又在 Bild 3.8 里与式 (3.85) 对照。它不是轨迹方程里直接代入的 $c_W$。

二次风、分级风、侧向喷入：已读 S. 36–S. 49 的式与所截正文，没有相应的式或句。本节未说明。

阻力式不需要打开 Wirsum 原书。式 (3.85)、式 (3.86) 已印在 S. 45。

## 固体含量、回落与带出：原文写到哪里

S. 45 印有：

> Die Lösung der Partikelbewegungsgleichung (3.82) erfolgt mit dem Ziel, die Aufenthaltszeiten jeder Partikelklasse im Freiraum zu ermitteln, um somit den Feststoffgehalt in den Bilanzräumen berechnen zu können.

求解式 (3.82) 的目的，是得到每个粒级在自由空间的停留时间，从而计算各平衡区的固体含量。由停留时间到固体含量的换算式，本节没有印出。

S. 48 印有：平衡区边界上的颗粒速度用来决定采用 Fall 1（$\Delta > 0$）还是 Fall 2（$\Delta < 0$）。$d_P^{*}$ 是使 $\Delta > 0$ 成立所需的粒径下限，同时还要求颗粒速度大于气体速度。下一句把三种情况用 „oder“ 并列，并说这时 $\Delta < 0$：颗粒慢于周围气体，颗粒回落，或粒径小于 $d_P^{*}$。这句没有给出回落质量或带出质量。

S. 49 末：解析算法约比 Runge–Kutta 快 50 倍，因此可以把粒径、抛出初速和球形度分得更细，以更好地逼近真实固体浓度剖面。这里仍然没有从轨迹到浓度、回落量或出口带出量的式。

各粒级向上通量、向下通量、返回床面的量和出口带出量：本节未给出。

---

## §3.2.1 抛入自由空间的固体

### 式 (3.49)，S. 36

**Q1** PDF 028，S. 36 左栏，式 (3.49)。与 S. 35 的式 (3.49) 同号，不是同一式。

**Q2** $\delta_e$：气泡顶盖固体层厚度。$\bar d_P$：床料平均粒径。本节未给单位。

**Q3** „die Feststoffschichtdicke des Blasendachs dem mittleren Durchmesser des Bettmaterials entspricht“

**Q4** Pemberton 与 Davidson (1986)，并写明与 Pyle (1965) 一致。

**Q5** 无空白。标注：文献并列。

**Q6** 比对阶段填写。

$$
\delta_e = \bar d_P
$$

### 式 (3.50)，S. 36

**Q1** PDF 028，S. 36 左栏，式 (3.50)。与 S. 35 的式 (3.50) 同号，不是同一式。

**Q2** $F_{0,r}$：一次气泡破裂经屋顶机制抛出的固体量。$\rho_P$ 颗粒密度，$\varepsilon_{mf}$ 最小流化空隙率，$d_b$ 气泡直径。本节未给单位。式面是一次破裂的量，不是单位床面积、单位时间的通量。

**Q3** „die durch den Roof-Mechanismus bei einer Blaseneruption ausgeworfene Feststoffmenge“

**Q4** 式 (3.49)，S. 36。顶盖面积按印出的 $\pi d_b^2 / 2$。

**Q5** 无空白。标注：文献并列（由 Pemberton 与 Davidson 的层厚写出）。

**Q6** 比对阶段填写。

$$
F_{0,r} = \bar d_P \cdot \rho_P \cdot (1 - \varepsilon_{mf}) \cdot \frac{\pi \cdot d_b^2}{2}
$$

### 式 (3.51)，S. 36

**Q1** PDF 028，S. 36 左栏，式 (3.51)。

**Q2** $f_b$：单位床面积、单位时间的气泡破裂次数。$\varepsilon_b$、$u_b$、$d_b$ 为床面气泡量。本节未给单位。

**Q3** „Die Anzahl der Blaseneruptionen $f_b$ bezogen auf die Bettoberfläche und die Zeit“

**Q4** 可见气泡体积流与床面平均气泡直径。二者来自 §3.1，本节不重写。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
f_b = \varepsilon_b \cdot u_b \cdot \frac{6}{\pi \cdot d_b^3}
$$

### 式 (3.52)，S. 37

**Q1** PDF 028，S. 37 右栏，式 (3.52)。

**Q2** 同一符号 $F_{0,r}$，此处是单位面积、单位时间的屋顶抛出质量。空隙率用 $\varepsilon_{mf}$。本节未给单位。

**Q3** 正文称此为屋顶机制造成的、按面积计的抛出固体质量。两相理论下顶盖空隙率取 $\varepsilon_{mf}$。

**Q4** 式 (3.50)，S. 36，乘以式 (3.51)。$\dfrac{\pi d_b^2}{2} \cdot f_b = \dfrac{3\varepsilon_b u_b}{d_b}$，与式面系数 3 一致。

**Q5** 与式 (3.50) 共用符号 $F_{0,r}$，含义从「一次破裂的量」换成「面积时间通量」。见式 (3.53)、式 (3.67)。标注：文献并列中的屋顶通量。

**Q6** 比对阶段填写。

$$
F_{0,r} = \bar d_P \cdot \rho_P \cdot (1 - \varepsilon_{mf}) \cdot \frac{3 \cdot \varepsilon_b \cdot u_b}{d_b}
$$

### 式 (3.53)，S. 37

**Q1** PDF 028，S. 37 右栏，式 (3.53)。

**Q2** 符号仍是 $F_{0,r}$。空隙率改为 $\varepsilon_d$。本节未给单位。

**Q3** „der Lückengrad des Blasendaches dem der Suspensionsphase entspricht“，于是式 (3.52) 改写为此式。

**Q4** 式 (3.52) 的结构；条件是悬浮相空隙率高于 $\varepsilon_{mf}$，且顶盖空隙率等于悬浮相。

**Q5** 式 (3.52) 与式 (3.53) 同号左端、不同空隙率。正文没有写模型在二者中取哪一个。可选解释：悬浮相密于最小流化时用式 (3.53)，否则用式 (3.52)。影响：进入自由空间的屋顶固体量，以及式 (3.67) 里 $(1-\varepsilon_d)$ 与所引式 (3.52) 是否同一套空隙率。不选定。标注：条件式，文献并列。

**Q6** 比对阶段填写。

$$
F_{0,r} = \bar d_P \cdot \rho_P \cdot (1 - \varepsilon_d) \cdot \frac{3 \cdot \varepsilon_b \cdot u_b}{d_b}
$$

### 式 (3.54)，S. 37

**Q1** PDF 028，S. 37 右栏，式 (3.54)。

**Q2** $F_{0,w}$：尾流机制的面积时间抛出通量。$\zeta_w$：尾流固体中被抛上的份额。$f_w$：尾流体积占气泡体积的份额。空隙率是 $\varepsilon_d$。本节未给单位。

**Q3** „Der flächen- und zeitbezogene Feststoffaustrag durch den Wake-Mechanismus wird nach Pemberton und Davidson (1986) … berechnet.“

**Q4** Pemberton 与 Davidson (1986)。S. 37 另写：单泡在床面破裂时抛出以屋顶为主；气速较高或颗粒较小时聚并增强，尾流机制占优。这是定性句，没有式号。

**Q5** $\zeta_w$、$f_w$ 的数值不在本式，见 S. 39 与式 (3.65)。标注：文献并列。

**Q6** 比对阶段填写。

$$
F_{0,w} = \zeta_w \cdot f_w \cdot \varepsilon_b \cdot u_b \cdot \rho_P \cdot (1 - \varepsilon_d)
$$

### 式 (3.55)，S. 37

**Q1** PDF 028，S. 37 右栏，式 (3.55)。

**Q2** $f(\Delta l / d_b)$：以泡间距与泡径之比 $\Delta l/d_b$ 为自变量的分布密度。$\Delta l_m$ 为主值。$\sigma$ 在正文中称为 Standardabweichung。本节未给 $\sigma$ 的单位。

**Q3** „Die nach Gauß normalverteilte Dichtefunktion ergibt sich zu“。间距分布来自 Hatano 与 Ishida (1983)、Hatano 等 (1984)；高斯形式经 Briens、Baron 写入。

**Q4** 床面附近气泡中心距的分布。

**Q5**

1. 指数。整段页面把指数印成 $-\tfrac{1}{2}$。减号在指数分数的分子左侧，紧贴数字 1，单独放大时容易漏读。转录取 $-\tfrac{1}{2}$。
2. $\sigma$ 的含义。正文称 $\sigma$ 为中心距的标准差，下一页又写 Briens 等 (1988) 取 $\sigma = 0{,}4$，没有单位。式面分母是 $\sigma \cdot d_b$，自变量是 $\Delta l/d_b$，前置因子是 $1/(\sigma\sqrt{2\pi})$。建议读法：$\sigma = 0{,}4$ 是无量纲比 $\Delta l/d_b$ 的标准差，与前置因子和 $d(\Delta l/d_b)$ 一致。另一读法：$\sigma$ 是长度，则 $0{,}4$ 缺单位，且分母里的 $d_b$ 会把间距无量纲化两次。不改式面。影响：式 (3.56) 的聚并概率，进而式 (3.64)、式 (3.67) 的尾流抛出。

标注：文献并列。

**Q6** 比对阶段填写。

$$
f\left(\frac{\Delta l}{d_b}\right)
  = \frac{1}{\sigma\sqrt{2\pi}}
    \exp\left(
      -\frac{1}{2}
      \left(\frac{\Delta l - \Delta l_m}{\sigma \cdot d_b}\right)^2
    \right)
$$

### 式 (3.56)，S. 38

**Q1** PDF 029，S. 38 左栏，式 (3.56)。

**Q2** $P'$：后一气泡中心距小于 $\Delta l = d_b$ 的概率。积分上限是数字 1，对应 $\Delta l/d_b = 1$。下限 $-\infty$。本节未给单位。

**Q3** „Die Wahrscheinlichkeit $P'$, dass die sich annähernde Blase den Abstand $\Delta l = d_b$ zur vorherigen Blase unterschreitet“

**Q4** 式 (3.55)。同页：„Für die Standardabweichung $\sigma$ setzt Briens et al. (1988) einen Wert von $\sigma = 0{,}4$ an.“ 这是 Briens 的取值，见式 (3.55) 的 Q5。

**Q5** 积分上限的字形是数字 1，与 $\Delta l = d_b$ 的句子一致。标注：文献并列。

**Q6** 比对阶段填写。

$$
P' = \int_{-\infty}^{1} f\left(\frac{\Delta l}{d_b}\right)\, d\left(\frac{\Delta l}{d_b}\right)
$$

### 式 (3.57)，S. 38

**Q1** PDF 029，S. 38 左栏，式 (3.57)。

**Q2** $\Delta l_m$：两气泡中心的平均距离。本节未给单位。

**Q3** Ishida 等 (1980)，由 Milioli 与 Foster (1995) 引用；假定气泡在床截面上均匀分布。

**Q4** $\varepsilon_b$、$d_b$。

**Q5** 指数按式面为 $-1/3$。标注：文献并列。

**Q6** 比对阶段填写。

$$
\Delta l_m = d_b \cdot \left(\frac{6\cdot\varepsilon_b}{\pi}\right)^{-1/3}
$$

### 式 (3.58)，S. 38

**Q1** PDF 029，S. 38 左栏，式 (3.58)。

**Q2** $\alpha_t$：气泡实际穿过的截面份额。本节未给单位。

**Q3** Baron 等 (1990)：气泡沿 bubble tracks 上升，计算平均间距时只计入被气泡穿过的截面。

**Q4** 式 (3.57) 的结构，分母增加 $\alpha_t$。$\alpha_t$ 见式 (3.59)。

**Q5** 无空白。标注：文献并列。与式 (3.57) 并列，正文没有写模型取哪一个。

**Q6** 比对阶段填写。

$$
\Delta l_m = d_b \cdot \left(\frac{6\cdot\varepsilon_b}{\pi\cdot\alpha_t}\right)^{-1/3}
$$

### 式 (3.59)，S. 38

**Q1** PDF 029，S. 38 左栏，式 (3.59)。

**Q2** $h_{WS}$：床高。$A_0$：式面分母是 $\sqrt{A_0}$，下标为数字 0。$\alpha_t$ 下标为 t。系数 $1{,}74$，指数 $0{,}4$。本节未给单位。

**Q3** Baron 用 Darton 等 (1977)，给出随床高求 $\alpha_t$ 的式。

**Q4** 床高与 $A_0$。

**Q5** 分母根号内的符号按放大图读作 $A_0$，不是 $A$。标注：文献并列。

**Q6** 比对阶段填写。

$$
\alpha_t = 1{,}74 \cdot \left(\frac{h_{WS}}{\sqrt{A_0}} + 4\right)^{0{,}4}
$$

### 式 (3.60)，S. 38

**Q1** PDF 029，S. 38 左栏，式 (3.60)。

**Q2** $P''$：气泡聚并概率。$C$：修正因子，无量纲。

**Q3** 聚并还须满足接近角 $\xi < \xi_{\ldots}$。„Damit gilt als Wahrscheinlichkeit der Blasenkoaleszenz $P''$ schließlich“

**Q4** 式 (3.56) 与下面的 $C$。Clift 与 Grace (1971)，经 Baron 等 (1990) 引用。

**Q5** 角的下标见式 (3.61)。标注：文献并列。

**Q6** 比对阶段填写。

$$
P'' = C \cdot P'
$$

### 式 (3.61)，S. 39

**Q1** PDF 029，S. 39 右栏顶部，式 (3.61)。说明句在 S. 38 底部。

**Q2** $\Delta l$：半径。$\xi$ 的下标见 Q5。$C$ 无量纲。式面未把分数约成 $1-\cos\xi$。

**Q3** „Dabei stellt der Korrekturfaktor $C$ das Verhältnis zwischen der Fläche der Halbkugel mit dem Radius $\Delta l$ und der Fläche der Kugelkappe mit dem Halbwinkel $\xi_{\ldots}$, dem Radius $\Delta l$.“

**Q4** 半球面积 $2\pi(\Delta l)^2$，球冠面积 $2\pi(\Delta l)^2(1-\cos\xi)$。

**Q5** 两处不能并成一种读法。

1. 比值方向。上句把 $C$ 说成半球面积与球冠面积之比，即 $1/(1-\cos\xi)$。式 (3.61) 印出的是球冠面积除以半球面积，结果为 $1-\cos\xi$。式 (3.62) 使用 $1-\cos\xi$，与式 (3.61) 的印刷结果一致，与这句的语序不一致。可选解释：按句子，$C = 1/(1-\cos\xi)$，则式 (3.62) 的因子要改；按式面，$C = 1-\cos\xi$，句子的语序与式不符。影响：聚并概率差一个 $1/(1-\cos\xi)^2$ 量级的因子，进入式 (3.64)、式 (3.67) 的尾流抛出。不选定。
2. 下标。公式裁图把下标读成数字 1。正文裁图把同一角读成字母 l，并称为 Halbwinkel / Annäherungswinkel，给出 30° 与 48°。建议读法：下标是字母 l，记 $\xi_l$。理由：它是角度，不是序号；式 (3.60) 写 $\xi < \xi_{\ldots}$。数字 1 与字母 l 在下标位置形近，故保留读法 $\xi_1$。下面公式按式面先写作 $\xi_1$，不表示已经排除 $\xi_l$。

标注：文献并列。

**Q6** 比对阶段填写。

$$
C = \frac{2\cdot\pi\cdot(\Delta l)^2\cdot(1-\cos\xi_1)}{2\cdot\pi\cdot(\Delta l)^2}
$$

### 式 (3.62)，S. 39

**Q1** PDF 029，S. 39 右栏，式 (3.62)。

**Q2** $P''$：同时满足间距与接近角时的两泡聚并概率。

**Q3** „die Bedingung des Annäherungswinkels erfüllt ist“

**Q4** 式 (3.56)、式 (3.61)。S. 39：Clift 与 Grace (1971) 给该角 30°；Briens 等 (1988) 给 48°。

**Q5** 30° 与 48° 并列，正文没有选择。下标读法同式 (3.61)。因子用的是 $1-\cos$，见式 (3.61) 的比值问题。标注：文献并列。

**Q6** 比对阶段填写。

$$
P'' = (1-\cos\xi_1)\cdot\int_{-\infty}^{1} f\left(\frac{\Delta l}{d_b}\right)\, d\left(\frac{\Delta l}{d_b}\right)
$$

### 式 (3.63)，S. 39

**Q1** PDF 029，S. 39 右栏，式 (3.63)。

**Q2** $P''_k$：$k$ 个气泡在床面附近聚并的概率。本节未给单位。

**Q3** „die Wahrscheinlichkeit, dass $k$ Blasen in der Nähe der Wirbelschichtoberfläche koaleszieren“。假定相继气泡间距相互独立。Briens 等 (1988)。

**Q4** 两泡概率 $P''$。

**Q5** 递推的初值没有印出。式 (3.64) 的求和从 $k = 2$ 起。可选解释：$P''_1 = P''$，于是 $P''_k = (1-P'')^{k-1} P''$。另一解释：下标从别的 $k$ 起算，初值是别的数。影响：多重聚并的权重。不补初值。标注：文献并列。

**Q6** 比对阶段填写。

$$
P''_k = (1-P'')\cdot P''_{k-1}
$$

### 式 (3.64)，S. 39

**Q1** PDF 029，S. 39 右栏，式 (3.64)。

**Q2** $F_{0,w}$：计入多重聚并后的尾流抛出通量。求和从 $k = 2$ 到 $n_B$。本节未给单位。

**Q3** Briens 等 (1988) 假定：$k$ 个气泡在床面附近聚并时，有 $k-1$ 个尾流的固体被抛入自由空间。观察来自 Levy 等 (1983)。

**Q4** 式 (3.54) 的系数，乘以 $\sum (k-1) P''_k$。

**Q5** 同页正文：$k$ 取大于或等于 2，否则不计入尾流抛出，那对应单个破裂气泡。$n_B$ „kann zur Berechnung auf $n_B = 5$ beschränkt werden“（Briens 等 1988）。这是可以限制为 5，不是写成必须取 5。

S. 39 继续写：Pemberton 与 Davidson (1986) 按 George 与 Grace (1981)，取尾流体积份额 $f_w = 0{,}25$，并取尾流固体的 40 % 被抛出，$\zeta_w = 0{,}4$。这是文献取值。紧接着 Milioli 与 Foster 用随 $u_0/u_{mf}$ 变化的乘积代替它，即式 (3.65)。

标注：文献并列。

**Q6** 比对阶段填写。

$$
F_{0,w}
  = \zeta_w\cdot f_w\cdot\varepsilon_b\cdot u_b\cdot\rho_P\cdot(1-\varepsilon_d)
    \cdot\sum_{k=2}^{n_B}(k-1)\cdot P''_k
$$

### 式 (3.65)，S. 40

**Q1** PDF 030，S. 40 左栏，式 (3.65)。

**Q2** 左端是乘积 $f_w\cdot\zeta_w$，无量纲。分支点是 $u_0/u_{mf} = 10$。系数 $0{,}1$、$0{,}5$，指数里的 $9{,}5$。本节未给单位。

**Q3** „Milioli und Foster (1995) wählen folgenden Ansatz“。他们按 Pemberton 与 Davidson (1986) 的数据，把乘积的上限估为 $\zeta_w f_w = 0{,}5$。正文写明：还没有针对多次聚并测定 $f_w$ 或 $\zeta_w$ 的研究。

**Q4** $u_0/u_{mf}$。

**Q5** 在 $u_0/u_{mf} = 10$ 处，左支为 $0{,}1$，右支为 $0{,}5\cdot(1-e^{-0{,}5})$。正文没有写间断点取哪一侧。式 (3.67) 使用符号 $\zeta_w\cdot f_w$，但没有写这个乘积取式 (3.65)，还是取 $0{,}25\times 0{,}4$。叙事上式 (3.65) 就在该乘积的讨论之后。可选解释及影响：取式 (3.65) 则抛出量随气速变；取常数 $0{,}1$ 则不随气速变。不选定。标注：文献并列（Milioli 与 Foster 的选择，Hamel 写出并在后文符号中沿用，见式 (3.67) 的 Q5）。

**Q6** 比对阶段填写。

$$
f_w\cdot\zeta_w =
\begin{cases}
0{,}1 & \dfrac{u_0}{u_{mf}} < 10 \\[1em]
0{,}5\cdot\left(1-e^{\left(9{,}5-\dfrac{u_0}{u_{mf}}\right)}\right) & \dfrac{u_0}{u_{mf}} \ge 10
\end{cases}
$$

### 式 (3.66)，S. 40

**Q1** PDF 030，S. 40 左栏，式 (3.66)。

**Q2** $F_0$：总固体抛出。本节未给单位。

**Q3** 模型设想把两个机制合在一起：单个破裂气泡走屋顶，床面附近聚并的气泡走尾流。

**Q4** $F_{0,r}$ 与 $F_{0,w}$。

**Q5** 没有写 $F_{0,r}$ 取式 (3.52) 还是式 (3.53)。标注：正文的模型式。

**Q6** 比对阶段填写。

$$
F_0 = F_{0,r} + F_{0,w}
$$

### 式 (3.67)，S. 40

**Q1** PDF 030，S. 40 左栏，式 (3.67)。

**Q2** $F_0$ 的展开式。空隙率因子是 $(1-\varepsilon_d)$。求和从 $k = 2$ 到 $n_B$。本节未给单位。

**Q3** „Unter Berücksichtigung der Wahrscheinlichkeit einer Einzelblaseneruption ergeben Gleichung (3.52) und (3.64) eingesetzt in (3.66) folgende Korrelation für den gesamten Feststoffausstoß“

**Q4** 正文点名式 (3.52)、式 (3.64)、式 (3.66)。

**Q5** 三处不能按正文那句话直接代换。

1. 句子引用式 (3.52)，式 (3.52) 的空隙率是 $(1-\varepsilon_{mf})$。式 (3.67) 印出的是 $(1-\varepsilon_d)$，与式 (3.53)、式 (3.64) 相同。
2. 若把式 (3.52) 与式 (3.64) 逐项相加，得不到 $(\zeta_w f_w - 3\bar d_P/d_b)$ 这一组合。单泡破裂概率没有单独的式号。句子里的 „Wahrscheinlichkeit einer Einzelblaseneruption“ 没有写成一个符号。
3. 求和的归一与取值范围，除 $k\ge 2$ 和 $n_B$ 可以限制为 5 之外，没有再写。

可选解释：方括号第一项是未聚并部分的屋顶贡献，第二项在聚并发生时用尾流份额替换屋顶份额。这只是对印刷结构的一种读法，正文没有写出这个替换步骤。影响：床面抛入自由空间的总固体量。不改式，不选定空隙率。标注：正文给出的总抛出式；所引式号与式面不完全一致。

**Q6** 比对阶段填写。

$$
F_0
  = \varepsilon_b\cdot u_b\cdot\rho_P\cdot(1-\varepsilon_d)
    \left[
      \frac{3\cdot\bar d_P}{d_b}
      + \left(\zeta_w\cdot f_w - \frac{3\cdot\bar d_P}{d_b}\right)
        \cdot\sum_{k=2}^{n_B}(k-1)\cdot P''_k
    \right]
$$

---

## §3.2.2 颗粒抛出初速

S. 40 在进入编号式之前写：Wen 与 Chen (1982) 观察到有的抛出速度超过气泡上升速度的 8 倍。这是文献陈述，没有式号。

### 式 (3.68)，S. 40

**Q1** PDF 030，S. 40 左栏，式 (3.68)。

**Q2** $\bar u_{P,0}$：屋顶机制的平均抛出初速。分母是 $12\cdot\bar d_P$。本节未给单位。

**Q3** Pemberton 与 Davidson (1986) 给屋顶机制抛出颗粒的平均初速。下一页指出：无损失的动量传递使这个值高得不真实（Bürkle 1991）。

**Q4** $u_b$、$d_b$、$\bar d_P$。

**Q5** 无空白。标注：文献并列，正文随即指出该值过高。后文采用的平均初速不是此式，见式 (3.73)。

**Q6** 比对阶段填写。

$$
\bar u_{P,0} = \frac{u_b\cdot d_b}{12\cdot\bar d_P}
$$

### 式 (3.69)，S. 40

**Q1** PDF 030，S. 40 左栏，式 (3.69)。

**Q2** 尾流机制的平均抛出初速。系数 $2{,}5$。本节未给单位。

**Q3** „Die Partikelstartgeschwindigkeit bedingt durch den Wake-Mechanismus berechnet sich nach Pemberton und Davidson (1986)“

**Q4** $u_b$。

**Q5** 无空白。标注：文献并列。式 (3.73) 用的是 $1{,}53$，不是 $2{,}5$。

**Q6** 比对阶段填写。

$$
\bar u_{P,0} = 2{,}5\cdot u_b
$$

### 式 (3.70)，S. 41

**Q1** PDF 030，S. 41 右栏，式 (3.70)。

**Q2** 系数 $0{,}092$，指数 $1{,}5$ 与 $-1{,}3$。本节未给 $0{,}092$ 的单位。

**Q3** Demmich 与 Bohnet (1978) 由鼓泡床自由空间上端固体流的测量拟合。Pemberton 与 Davidson 认为此式偏高，因为没有自由射流。

**Q4** $u_b$、$d_b$、$\bar d_P$、$\rho_g/\rho_P$。

**Q5** $0{,}092$ 不是无量纲数，原文未说明单位。标注：文献并列。

**Q6** 比对阶段填写。

$$
\bar u_{P,0}
  = 0{,}092\cdot u_b\cdot d_b^{1{,}5}\cdot\bar d_P^{-1{,}3}
    \cdot\sqrt{\frac{\rho_g}{\rho_P}}
$$

S. 41 在式 (3.70) 与式 (3.71) 之间，Hatano 与 Ishida (1981) 给出区间，没有式号：单泡 $0{,}3 < \bar u_{P,0}/u_b < 1{,}1$；聚并 $1{,}1 < \bar u_{P,0}/u_b < 2{,}5$。文献并列。

### 式 (3.71)，S. 41

**Q1** PDF 030，S. 41 右栏，式 (3.71)。

**Q2** 系数 $7{,}92\cdot 10^{-6}$，指数 $-1{,}65$。没有 $d_b$。本节未给该系数的单位。

**Q3** Son 等 (1988, 1989)，引自 Fung 与 Hamdullahpur (1993b)。适用范围写成环境温度、空气。

**Q4** $u_b$、$\bar d_P$。

**Q5** 系数的单位原文未说明。标注：文献并列。图见 Bild 3.7。

**Q6** 比对阶段填写。

$$
\bar u_{P,0} = 7{,}92\cdot 10^{-6}\cdot u_b\cdot\bar d_P^{-1{,}65}
$$

Martens (1984)：累积质量 50 % 处 $u_{P,0,50\%} = 1{,}53\, u_b$，全部速度落在 $0 < u_{P,0}/u_b < 4$。Berkelmann (1991) 与此相符。这是文献陈述，式 (3.73) 再按它取平均初速。

### 式 (3.72)，S. 41

**Q1** PDF 030，S. 41 右栏，式 (3.72)。

**Q2** $dF$：式面左端。$\sigma$ 出现在高斯密度的标准差位置。本节未给单位。

**Q3** 模拟模型只做轴向一维，因此只要竖直分量（Bürkle 1991）。„Nach dem Vorschlag von Levy (1986) wird die vertikale Partikelstartgeschwindigkeit unter Verwendung einer Gaußschen Normalverteilung berechnet.“ 正文接着说，抛出固体中初速落在 $u_{P,0}+du_{P,0}$ 区间的份额 $F_{0,i}$ 由下式给出。

**Q4** $F_0$，平均初速见式 (3.73)。

**Q5**

1. 正文写份额 $F_{0,i}$ 和区间 $u_{P,0}+du_{P,0}$。式面是 $dF$，没有 $du_{P,0}$。可选解释：$dF$ 是密度，区间份额还要乘 $du_{P,0}$；或 $dF$ 已经是该微分区间的量。不把 $du_{P,0}$ 补进式子。影响：各初速粒级的抛出质量是否归一到 $F_0$。
2. S. 42 把 $\sigma$ 称为 Varianz，又令 $\sigma = 0{,}6$，并放在本式标准差的位置。见式 (3.73)。

标注：Hamel 写明采用（Levy 1986，用于模拟模型）。

**Q6** 比对阶段填写。

$$
dF
  = F_0\cdot\frac{1}{\sigma\sqrt{2\pi}}
    \exp\left(
      -\frac{1}{2}
      \left(\frac{u_{P,0}-\bar u_{P,0}}{\sigma}\right)^2
    \right)
$$

### 式 (3.73)，S. 42

**Q1** PDF 031，S. 42 左栏，式 (3.73)。

**Q2** 平均竖直抛出初速。系数 $1{,}53$。$\sigma = 0{,}6$。本节未给单位。

**Q3** „Die Varianz wird mit $\sigma = 0{,}6$ angenommen und die mittlere Partikelstartgeschwindigkeit wird in Anlehnung an Martens (1984) berechnet“

**Q4** 床面气泡上升速度 $u_b$。Martens 的 50 % 累积速度是 $1{,}53\, u_b$，见上文。这里用约等号，并称为 Mittelwert，不是把 50 % 分位数原样改名。

**Q5** 正文称 $\sigma$ 为 Varianz，数值 $0{,}6$ 放在式 (3.72) 的标准差位置。可选解释：把 $0{,}6$ 当作标准差；或把 $0{,}6$ 当作方差，则式 (3.72) 的 $\sigma$ 应是 $\sqrt{0{,}6}$。影响：初速分布的宽窄，从而改变能飞过各层的颗粒份额。不选定。标注：Hamel 写明采用。

**Q6** 比对阶段填写。

$$
\bar u_{P,0} \approx 1{,}53\cdot u_b
$$

### Bild 3.7，S. 42

**Q1** PDF 031，S. 42 左栏，Bild 3.7。不数字化曲线。

**Q2** 图题：抛出初速随床面气泡速度变化。气泡直径与速度按第 3.1 节的式计算。图内可见 $\rho_P = 2200\,\mathrm{kg/m^3}$。横坐标为 $u_{b,WS}$，单位 m/s，刻度画到 $3{,}5$。图中另有 $d_{b,WS}$ 与抛出速度的标尺。图题没有写出温度、压力和粒径。

**Q3** „Partikelstartgeschwindigkeit in Abhängigkeit der Blasengeschwindigkeit an der Wirbelschichtoberfläche.“

**Q4** 式 (3.68)–式 (3.71) 以及 Hatano 与 Ishida 的区间。这些是图上的文献曲线。

**Q5** 曲线不逐点读取。标注：文献并列的图，不是模型采用的式。

**Q6** 比对阶段填写。

---

## §3.2.3 自由空间里的气体运动

S. 42–S. 43，在式 (3.74) 之前。

气泡在床面破裂后，悬浮相气体与气泡气体不立刻混合，逸出的气泡气体先保持封闭运动。Bürkle (1991) 给出试验装置上的综述。

Zenz 与 Weil (1958) 建议用自由射流耗散来类比。George 与 Grace (1978)、Horio 等 (1980) 沿用这一设想。George 与 Grace (1978)：破裂气泡形成非稳态自由射流，随高度耗散。这一支没有式号。标注：文献并列。

另一支：Hamdullahpur 与 MacKay (1986)、Sciazko 等 (1988) 沿用 Pemberton 与 Davidson (1984) 的鬼泡；Bürkle (1991) 也讨论继续存在的气泡或环涡。下面三式属于这一支。

### 式 (3.74)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.74)。

**Q2** $u_{gb}(h)$：高度 $h$ 处的鬼泡速度。$u_0$：来流速度。$u_{gb,0}$：鬼泡初速。$\beta_A$：衰减系数，下标为 A。本节未给单位。

**Q3** 力平衡可以化成这一函数（Sciazko 等 1988，Fung 与 Hamdullahpur 1993a）。

**Q4** 式 (3.75)、式 (3.76)。

**Q5** 颗粒运动方程里的 $u_g$ 是否就等于 $u_{gb}(h)$，正文没有写。可选解释：$u_g = u_{gb}(h)$；或自由空间主体速度仍是 $u_0$，鬼泡速度只描述气泡气体。影响：式 (3.82) 的相对速度，以及细颗粒是被带着走还是按床面气速计算。不选定。标注：与式 (3.76) 的采用句连在一起的鬼泡速度式。射流类比没有对应的编号式。

**Q6** 比对阶段填写。

$$
u_{gb}(h) = u_0 + (u_{gb,0}-u_0)\cdot e^{-\beta_A h}
$$

### 式 (3.75)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.75)。

**Q2** $\beta_A$ 的下标是 A。$d_{b,WS}$：床面气泡直径。本节未给单位。

**Q3** Hamdullahpur 与 MacKay (1986)、Sciazko 等 (1988)，按床面气泡直径的测量写出。

**Q4** $d_{b,WS}$。

**Q5** 下标读作 A，不是 $\lambda$。标注：与式 (3.76) 一起采用的鬼泡衰减。

**Q6** 比对阶段填写。

$$
\beta_A = \frac{1}{3\cdot d_{b,WS}}
$$

### 式 (3.76)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.76)。

**Q2** $u_{gb,0}$ 取床面气泡上升速度 $u_{b,WS}$。本节未给单位。

**Q3** „In Anlehnung an Hamdullahpur und MacKay (1986) und Sciazko et al. (1988) wird für die Startgeschwindigkeit der Ghost-bubbles die Blasensteiggeschwindigkeit … an der Wirbelschichtoberfläche angenommen“

**Q4** §3.1 的床面气泡速度。

**Q5** 采用句写的是鬼泡初速，没有写 $u_g = u_{gb}$。见式 (3.74) 的 Q5。标注：Hamel 写明采用。

**Q6** 比对阶段填写。

$$
u_{gb,0} = u_{b,WS}
$$

---

## §3.2.4 自由空间里的颗粒运动

S. 43：床面抛出颗粒的轨迹是固体运动和自由空间浓度剖面的基础。单颗粒力平衡，忽略颗粒之间的作用力。力有惯性、浮力、重力和阻力。

### 式 (3.77)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.77)。

**Q2** $F_T$ 惯性，$F_A$ 浮力，$F_G$ 重力，$F_W$ 阻力。

**Q3** „Die Kräftebilanz für das Einzelpartikel unter Berücksichtigung der Trägheits-, Auftriebs-, Gewichts- und Widerstandskraft lautet“

**Q4** 下面四式。

**Q5** 无空白。标注：正文力平衡。

**Q6** 比对阶段填写。

$$
F_T + F_A + F_G + F_W = 0
$$

### 式 (3.78)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.78)。

**Q2** $\rho_P$、$V_P$、$u_P$。本节未给单位。

**Q3** 作用在颗粒上的惯性力。

**Q4** 无上游编号式。

**Q5** 负号与 $du_P/dt>0$ 时力为负一致。正文没有另写「向上为正」。标注：正文力平衡。

**Q6** 比对阶段填写。

$$
F_T = -\rho_P\cdot V_P\cdot\frac{du_P}{dt}
$$

### 式 (3.79)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.79)。

**Q2** 浮力。$\rho_g$、$V_P$、$g$。本节未给单位。

**Q3** 浮力项。

**Q4** 无。

**Q5** 无空白。标注：正文力平衡。

**Q6** 比对阶段填写。

$$
F_A = \rho_g\cdot V_P\cdot g
$$

### 式 (3.80)，S. 43

**Q1** PDF 031，S. 43 右栏，式 (3.80)。

**Q2** 重力。本节未给单位。

**Q3** 重力项。

**Q4** 无。

**Q5** 无空白。标注：正文力平衡。

**Q6** 比对阶段填写。

$$
F_G = -\rho_P\cdot V_P\cdot g
$$

### 式 (3.81)，S. 44

**Q1** PDF 032，S. 44 左栏，式 (3.81)。

**Q2** $c_W$ 阻力系数。相对速度写成 $(u_P-u_g)\,|u_P-u_g|/2$。面积 $\pi d_P^2/4$。本节未给单位。

**Q3** 阻力方向取决于相对速度。$u_P > u_g$，或颗粒下落（$u_P < 0$）时，阻力与惯性力 $F_T$ 反向。颗粒上升但慢于气体（$u_P > 0$ 且 $u_P < u_g$）时，阻力是曳力（Schleppkraft）。

**Q4** $c_W$ 在后面才给出算法。

**Q5** 无空白。标注：正文力平衡。

**Q6** 比对阶段填写。

$$
F_W
  = -c_W\cdot\rho_g\cdot\left(\frac{\pi\cdot d_P^2}{4}\right)
    \cdot\left(\frac{(u_P-u_g)\,|u_P-u_g|}{2}\right)
$$

### 式 (3.82)，S. 44

**Q1** PDF 032，S. 44 左栏，式 (3.82)。

**Q2** $u_r$ 见式 (3.83)。本节未给单位。

**Q3** 把各力代入式 (3.77) 得到自由空间颗粒运动的微分方程。

**Q4** 式 (3.77)–式 (3.81)。由力平衡可得到系数 $3/4$：阻力中的 $1/2$ 与球体积 $\pi d^3/6$、投影面积 $\pi d^2/4$ 合在一起。式面已经是合并后的形式。

**Q5** 无空白。标注：正文导出。S. 45 说求解本式是为了求各粒级停留时间，见文首「固体含量、回落与带出」。

**Q6** 比对阶段填写。

$$
\frac{du_P}{dt}
  = \frac{-3\cdot c_W\cdot\rho_g\cdot u_r\cdot|u_r|}{4\cdot\rho_P\cdot d_P}
    - \frac{(\rho_P-\rho_g)\,g}{\rho_P}
$$

### 式 (3.83)，S. 44

**Q1** PDF 032，S. 44 左栏，式 (3.83)。

**Q2** $u_r = u_P - u_g$。本节未给单位。

**Q3** „mit“

**Q4** 式 (3.82)。$u_g$ 取什么，见式 (3.74) 的 Q5。

**Q5** 无空白。标注：正文定义。

**Q6** 比对阶段填写。

$$
u_r = u_P - u_g
$$

S. 44 继续写：通常假定球形，因为 $c_W$ 可以近似。非球形使用 Tabelle 3.1 被称作一种简化。Martens (1984) 建议形状因子。Haider 与 Levenspiel (1989) 把球形与非球形数据按球形度整理：$\varphi_S$ 为等体积球的表面积除以实际表面积。

### 式 (3.84)，S. 44

**Q1** PDF 032，S. 44 左栏，式 (3.84)。

**Q2** $c_W$ 对 $Re$ 与 $\varphi_S$。适用范围写到 $\varphi_S = 1$ 直至盘状 $\varphi_S = 0{,}026$。系数按式面：$8{,}1716$、$4{,}0655$、$0{,}0964$、$0{,}5565$、$73{,}69$、$5{,}0784$、$5{,}378$、$6{,}2122$。

**Q3** „Gleichung (3.84) repräsentiert die von Haider und Levenspiel (1989) vorgenommene Parameteranpassung“

**Q4** $Re$、$\varphi_S$。本式用来计算式 (3.86) 的 $z$，并在 Bild 3.8 中与式 (3.85) 对照。轨迹方程直接代入的是式 (3.85)，不是本式。

**Q5** 经验数字按两次页面阅读一致。标注：文献公式；同时是已采用的 $z$ 的输入。

**Q6** 比对阶段填写。

$$
\begin{aligned}
c_W
  &= \frac{24}{Re}\left[
       1 + \left(8{,}1716\cdot e^{-4{,}0655\,\varphi_S}\right)
         \cdot Re^{(0{,}0964+0{,}5565\,\varphi_S)}
     \right] \\
  &\quad + \frac{73{,}69\cdot Re\cdot e^{-5{,}0784\,\varphi_S}}
             {Re + 5{,}378\cdot e^{6{,}2122\,\varphi_S}}
\end{aligned}
$$

同页：Fung 与 Hamdullahpur (1993) 用 Al Taweel 等 (1989) 描述颗粒团或颗粒球在式 (3.82) 中的阻力，其中包含形状和团内循环气体。没有单独式号。标注：文献并列，正文没有说本模型采用 Al Taweel。

### 式 (3.85)，S. 45

**Q1** PDF 032，S. 45 右栏，式 (3.85)。

**Q2** $c_W = 24/Re + z$。$z$ 是很大雷诺数 $Re_\infty$ 下的 $c_W$。

**Q3** „Daher wird hier dem Vorschlag von Wirsum (1998) gefolgt, der durch geeignete Wahl der Gleichung für den Widerstandsbeiwert die Partikelbewegungsgleichung in eine analytisch lösbare Form überführt.“ Bürkle (1991) 建议的 Runge–Kutta 被写为计算时间过长。

**Q4** 式 (3.82)。$z$ 由式 (3.86) 从式 (3.84) 算出。

**Q5** $Re_\infty$ 没有给出数值。式 (3.86) 保留 $24/Re_\infty$ 项，没有取 $Re\to\infty$ 的极限。标注：Hamel 写明采用 Wirsum (1998)。未打开 Wirsum 原书；阻力形式已印在本页。

**Q6** 比对阶段填写。

$$
c_W = \frac{24}{Re} + z
$$

### 式 (3.86)，S. 45

**Q1** PDF 032，S. 45 右栏，式 (3.86)。

**Q2** 与式 (3.84) 相同的系数，雷诺数全部换成 $Re_\infty$。

**Q3** „Dabei ist $z$ der Wert von $c_W$ für sehr große Reynolds-Zahlen $Re_\infty$, der unter Zuhilfenahme der von Haider und Levenspiel (1989) vorgeschlagenen Gleichung (3.84) ermittelt wird.“

**Q4** 式 (3.84)。

**Q5** $Re_\infty$ 无数值。可选解释：令 $Re_\infty\to\infty$，则 $24/Re_\infty$ 消失，$z \to 73{,}69\, e^{-5{,}0784\,\varphi_S}$。这不是式面。式面是把式 (3.84) 的 $Re$ 换成 $Re_\infty$ 后的未化简式。不取极限。标注：Hamel 写明采用的 $z$，输入来自 Haider 与 Levenspiel。

**Q6** 比对阶段填写。

$$
\begin{aligned}
z
  &= \frac{24}{Re_\infty}\left[
       1 + \left(8{,}1716\cdot e^{-4{,}0655\,\varphi_S}\right)
         \cdot Re_\infty^{(0{,}0964+0{,}5565\,\varphi_S)}
     \right] \\
  &\quad + \frac{73{,}69\cdot Re_\infty\cdot e^{-5{,}0784\,\varphi_S}}
             {Re_\infty + 5{,}378\cdot e^{6{,}2122\,\varphi_S}}
\end{aligned}
$$

### 式 (3.87)，S. 45

**Q1** PDF 032，S. 45 右栏，式 (3.87)。

**Q2** $Re = |u_r|\, d_P / \nu_g$。分母是运动粘度 $\nu_g$。

**Q3** „mit“

**Q4** 式 (3.83)、式 (3.85)。

**Q5** 无空白。标注：随式 (3.85) 采用。

**Q6** 比对阶段填写。

$$
Re = \frac{|u_r|\cdot d_P}{\nu_g}
$$

S. 45 在式 (3.87) 之后：轨迹计算假定球形度分布在 $0{,}3 \le \varphi_S \le 1{,}0$。球形度平均值是程序输入。若无数据，计算采用 $\bar\varphi_S = 0{,}8$。标注：Hamel 写明采用的缺省值。

Bild 3.8 的说明也在这里：Haider 与 Levenspiel (1989) 的阻力系数与式 (3.85) 的近似对照。图本身在 S. 47。

### Tabelle 3.1，S. 45–S. 46

**Q1** PDF 032，S. 45 右栏起，续到 PDF 033，S. 46 左栏。表题：„Berechnungsgleichungen für den Widerstandsbeiwert $c_W$ eines kugelförmigen Einzelpartikels.“

**Q2** 五条球形单颗粒阻力式。系数按表面。

**Q3** 表题如上。引入时放在文献给出的公式之中。

**Q4** 不进入式 (3.85)。式 (3.88) 的句子写明用式 (3.85) 改写式 (3.82)。

**Q5** 无空白。标注：文献并列，不是轨迹采用的阻力。

**Q6** 比对阶段填写。

| 来源 | 式面 |
|---|---|
| Kaskas (1964) | $c_W = 24/Re + 4/\sqrt{Re} + 0{,}4$ |
| Martin (1984) | $c_W = \dfrac{1}{3}\left(\sqrt{72/Re}+1\right)^2$ |
| Do 等 (1972) | $c_W = \dfrac{24}{Re}(1+0{,}15\cdot Re^{0{,}687}) + \dfrac{0{,}42}{1+4{,}25\cdot 10^{4}\cdot Re^{-1{,}16}}$ |
| Turton 与 Levenspiel，引自 Haider 与 Levenspiel (1989) | $c_W = \dfrac{24}{Re}(1+0{,}1806\cdot Re^{0{,}6459}) + \dfrac{0{,}4251}{1+6880{,}95/Re}$ |
| Muschelknautz (1988)，引自 Hannes (1996) | $c_W = 21{,}5/Re + 6{,}5/\sqrt{Re} + 0{,}23$ |

### 式 (3.88)，S. 46

**Q1** PDF 033，S. 46 左栏，式 (3.88)。

**Q2** 把运动方程写成 $u_P$ 的二次多项式。系数 $A$、$B$、$C$ 见下三式。

**Q3** „Mit Gleichung (3.85) lässt sich die Differentialgleichung für die Partikelbewegung (3.82) in folgender Form schreiben“

**Q4** 式 (3.82)、式 (3.85)。

**Q5** 无空白。标注：随 Wirsum 阻力采用后的正文形式。

**Q6** 比对阶段填写。

$$
\frac{du_P}{dt} = A\cdot u_P^2 + B\cdot u_P + C
$$

符号约定在 S. 47，不在这三式旁边。原句：

> Bewegt sich ein Partikel langsamer als das ihn umgebende Gas oder fällt ein Partikel entgegen dem Gasstrom, so sind zur Berechnung der Koeffizienten $A$, $B$ und $C$ die in Klammern stehenden Vorzeichen gültig. Ist die Partikelgeschwindigkeit größer als die Gasgeschwindigkeit gelten entsprechend die nicht in Klammern stehenden Vorzeichen.

括号内的符号：颗粒慢于周围气体，或逆气流下落。括号外的符号：颗粒速度大于气体速度。

### 式 (3.89)，S. 46

**Q1** PDF 033，S. 46 左栏，式 (3.89)。

**Q2** 式面为 $A = -(+)\dfrac{3}{4}\, z\, \dfrac{\rho_g}{\rho_P\cdot d_P}$。不在括号中的符号是减号，括号中的符号是加号。

**Q3** „mit den Koeffizienten“

**Q4** 式 (3.85) 的 $z$ 项与相对速度的符号。

**Q5** 无空白。照录括号。标注：随式 (3.88) 采用。

**Q6** 比对阶段填写。

$$
A = -(+)\,\frac{3}{4}\cdot z\cdot\frac{\rho_g}{\rho_P\cdot d_P}
$$

### 式 (3.90)，S. 46

**Q1** PDF 033，S. 46 左栏，式 (3.90)。

**Q2** 第一项是 Stokes 项，没有括号备选符号，系数 18，粘度是 $\mu_g$。第二项系数是 $3/2$，括号写在加号位置：括号外为加，括号内为减。

**Q3** 系数 $B$。

**Q4** 式 (3.85)，$u_r = u_P-u_g$ 展开后的一次项。

**Q5** 无空白。标注：随式 (3.88) 采用。

**Q6** 比对阶段填写。

$$
B
  = -\frac{18\cdot\mu_g}{\rho_P\cdot d_P^2}
    + (-)\,\frac{3}{2}\cdot z\cdot\frac{\rho_g}{\rho_P\cdot d_P}\cdot u_g
$$

### 式 (3.91)，S. 46

**Q1** PDF 033，S. 46 左栏，式 (3.91)。

**Q2** 第一项没有括号备选符号。$z$ 项前面是 $-(+)$，即括号外为减、括号内为加。重力项为 $-\dfrac{\rho_P-\rho_g}{\rho_P}\, g$，没有括号备选符号。

**Q3** 系数 $C$。

**Q4** 式 (3.82)、式 (3.85)。

**Q5** 无空白。标注：随式 (3.88) 采用。

**Q6** 比对阶段填写。

$$
C
  = \frac{18\cdot\mu_g\cdot u_g}{\rho_P\cdot d_P^2}
    - (+)\,\frac{3}{4}\cdot z\cdot\frac{\rho_g}{\rho_P\cdot d_P}\cdot u_g^2
    - \frac{\rho_P-\rho_g}{\rho_P}\cdot g
$$

### Bild 3.8，S. 47

**Q1** PDF 033，S. 47 右栏，Bild 3.8。不数字化曲线。

**Q2** $c_W$ 对雷诺数，球形度取不同值。说明写明：实线为式 (3.84)，虚线为式 (3.85)。图上标有 $\varphi_S = 0{,}7$、$0{,}85$、$1{,}0$。坐标是对数轴。

**Q3** S. 45 已写本图是 Haider 与 Levenspiel (1989) 与式 (3.85) 的对照。

**Q4** 式 (3.84)、式 (3.85)。

**Q5** 不读曲线数值。标注：对照图。轨迹采用的是式 (3.85)。

**Q6** 比对阶段填写。

### 式 (3.92)，S. 47

**Q1** PDF 033，S. 47 右栏，式 (3.92)。

**Q2** $\Delta = 4AC - B^2$。这是印出的顺序，不是 $B^2-4AC$。

**Q3** 两类积分按 Bronstein (1989)，取决于该表达式。

**Q4** 式 (3.89)–式 (3.91)。

**Q5** $\Delta = 0$ 的情形没有给出。标注：随解析解采用。

**Q6** 比对阶段填写。

$$
\Delta = 4AC - B^2
$$

### Fall 1：$\Delta > 0$，S. 47–S. 48

### 式 (3.93)，S. 47

**Q1** PDF 033，S. 47 右栏，式 (3.93)。

**Q2** $u_P(t)$。$t_0$ 为该段积分的起始时刻。$\omega$ 见式 (3.94)。

**Q3** „Fall 1: $\Delta > 0$“。对式 (3.88) 积分，得到飞行时间上的颗粒速度。

**Q4** 式 (3.88)、式 (3.92)、式 (3.94)。

**Q5** 正切的宗量是 $\tfrac{1}{2}\omega\sqrt{\Delta} + \tfrac{1}{2}(t-t_0)\sqrt{\Delta}$。$\omega$ 项不含 $(t-t_0)$。标注：解析解，Bronstein (1989) 这类积分。

**Q6** 比对阶段填写。

$$
u_P(t)
  = \frac{\sqrt{\Delta}}{2A}
    \tan\left[
      \frac{1}{2}\omega\sqrt{\Delta}
      + \frac{1}{2}(t-t_0)\sqrt{\Delta}
    \right]
    - \frac{B}{2A}
$$

### 式 (3.94)，S. 48

**Q1** PDF 034，S. 48 左栏，式 (3.94)。

**Q2** $\omega$ 由初速 $u_{P,0}$ 确定。

**Q3** „mit“

**Q4** 式 (3.93) 的初条件。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
\omega = \frac{2}{\sqrt{\Delta}}\,\arctan\left[\frac{2A\, u_{P,0}+B}{\sqrt{\Delta}}\right]
$$

### 式 (3.95)，S. 48

**Q1** PDF 034，S. 48 左栏，式 (3.95)，两行。

**Q2** $h_P(t)$：颗粒在床面以上的高度。$h_{P,0}$ 为该段起始高度。

**Q3** 由 $u_P\, dt = dh$ 积分，得到床面以上高度随时间的变化。床面用词是 Wirbelschichtoberfläche。

**Q4** 式 (3.93)、式 (3.94)。

**Q5** 第二行是 $-\dfrac{B}{2A}(t-t_0) + \dfrac{1}{A}\ln\left[\cos\left(\dfrac{1}{2}\omega\sqrt{\Delta}\right)\right]$。标注：解析解。

**Q6** 比对阶段填写。

$$
\begin{aligned}
h_P(t)
  &= h_{P,0}
    - \frac{1}{A}\ln\left[
        \cos\left(
          -\frac{1}{2}\omega\sqrt{\Delta}
          - \frac{1}{2}(t-t_0)\sqrt{\Delta}
        \right)
      \right] \\
  &\quad - \frac{B}{2A}(t-t_0)
    + \frac{1}{A}\ln\left[
        \cos\left(\frac{1}{2}\omega\sqrt{\Delta}\right)
      \right]
\end{aligned}
$$

### 式 (3.96)，S. 48

**Q1** PDF 034，S. 48 左栏，式 (3.96)。

**Q2** $h_P$ 对 $u_P$。对数项系数是 $-1/(2A)$。两个反正切的系数分别是 $-B/(A\sqrt{\Delta})$ 与 $+B/(A\sqrt{\Delta})$。

**Q3** „Der Zusammenhang zwischen Höhe des Partikels und seiner Geschwindigkeit lautet“

**Q4** 式 (3.88) 在 Fall 1 下对高度积分。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
\begin{aligned}
h_P(u_P)
  &= h_{P,0}
    - \frac{1}{2A}\ln\left[
        \frac{A u_P^2 + B u_P + C}{A u_{P,0}^2 + B u_{P,0} + C}
      \right] \\
  &\quad - \frac{B}{A\sqrt{\Delta}}
      \arctan\left[\frac{2A u_P+B}{\sqrt{\Delta}}\right]
    + \frac{B}{A\sqrt{\Delta}}
      \arctan\left[\frac{2A u_{P,0}+B}{\sqrt{\Delta}}\right]
\end{aligned}
$$

### 式 (3.97)，S. 48

**Q1** PDF 034，S. 48 左栏，式 (3.97)。

**Q2** $d_P^{*}$：Partikelgrenzdurchmesser。分子 $108\,\mu_g^2$。分母 $z\, g\,(\rho_P\rho_g - \rho_g^2)$。指数 $1/3$。

**Q3** „Um die Bedingung $\Delta > 0$ zu erfüllen, muss die Geschwindigkeit des Partikels größer sein als die Gasgeschwindigkeit und der Partikeldurchmesser muss größer sein als der Partikelgrenzdurchmesser.“

**Q4** $\Delta = 4AC-B^2 > 0$，并取颗粒快于气体的那一组符号。$\Delta$ 在这组符号下不依赖 $u_g$，故界限可以写成只含物性与 $z$ 的粒径。式面系数是 108，分母是 $\rho_P\rho_g-\rho_g^2$。

**Q5** $d_P^{*}$ 是 Fall 1 的界限直径，不是带出与回落的分界式。正文没有把 $d_P^{*}$ 定义成扬析切割粒径。见下一段。标注：解析解的判别界限。

**Q6** 比对阶段填写。

$$
d_P^{*}
  = \left(
      \frac{108\,\mu_g^2}{z\, g\,(\rho_P\rho_g-\rho_g^2)}
    \right)^{1/3}
$$

S. 48 紧接的两句：

> Innerhalb des Zellenmodells ist es zweckmäßig zunächst die Geschwindigkeit des Partikels an den Grenzen des betrachteten Bilanzraumes zu bestimmen, da über die Partikelgeschwindigkeit entschieden wird, ob die Lösung der Differentialgleichung für Fall 1 mit $\Delta > 0$ oder für Fall 2 mit $\Delta < 0$ angewendet werden muss.

> Bewegt sich ein Partikel langsamer als das ihn umgebende Gas, fällt ein Partikel zurück oder ist der Partikeldurchmesser kleiner als der Partikelgrenzdurchmesser, nimmt der Faktor $\Delta$ einen Wert kleiner Null an.

平衡区边界上的速度用来选择解析分支。三种情况用 „oder“ 连接，对应 $\Delta < 0$。这里没有回落质量、带出质量或各粒级向上、向下通量。

### Fall 2：$\Delta < 0$，S. 48–S. 49

### 式 (3.98)，S. 48

**Q1** PDF 034，S. 48 左栏，式 (3.98)。

**Q2** $u_P(t)$。指数里是 $\eta\sqrt{-\Delta} - (t-t_0)\sqrt{-\Delta}$。分子末项 $-\sqrt{-\Delta}-B$ 在指数因子之外。分母是 $2A$ 乘同一指数，再减 $2A$。

**Q3** „Die Geschwindigkeit des Partikels als Funktion der Zeit ist durch … gegeben“。前文：„Fall 2: $\Delta < 0$“。

**Q4** 式 (3.88)、式 (3.99)。

**Q5** 分子按页面读作 $(\sqrt{-\Delta}-B)\exp[\cdots] - \sqrt{-\Delta} - B$，末两项不乘指数。放大图与整式裁图一致。标注：解析解。

**Q6** 比对阶段填写。

$$
u_P(t)
  = \frac{
      (\sqrt{-\Delta}-B)\,
      \exp\big[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}\big]
      - \sqrt{-\Delta} - B
    }{
      2A\,
      \exp\big[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}\big]
      - 2A
    }
$$

### 式 (3.99)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.99)。

**Q2** $\eta$ 由初速确定。对数的分子是 $\sqrt{-\Delta} + (2A u_{P,0}+B)$，分母是 $\sqrt{-\Delta} - (2A u_{P,0}+B)$。

**Q3** „mit“

**Q4** 式 (3.98) 的初条件。

**Q5** 分母已在放大图上读全，与分子对称，中间是减号。标注：解析解。

**Q6** 比对阶段填写。

$$
\eta
  = \frac{1}{\sqrt{-\Delta}}
    \ln\left(
      \frac{\sqrt{-\Delta}+(2A\cdot u_{P,0}+B)}
           {\sqrt{-\Delta}-(2A\cdot u_{P,0}+B)}
    \right)
$$

### 式 (3.100)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.100)。

**Q2** $h_P(t)$。本句床面用词是 Wirbelbettoberfläche。系数中的 $0{,}5$ 按式面小数逗号。

**Q3** „Für die Höhe des Partikels über der Wirbelbettoberfläche als Funktion der Zeit gilt“

**Q4** 式 (3.98)、式 (3.101)、式 (3.102)。

**Q5** 无空白。用词与前后文的 Wirbelschichtoberfläche 不同，式面照录。

**Q6** 比对阶段填写。

$$
h_P(t)
  = h_{P,0}
    - \frac{1}{A}\ln[\gamma+1]
    - \frac{(-0{,}5\cdot\sqrt{-\Delta}-0{,}5\cdot B)}{A}
      \big(\eta-(t-t_0)\big)
    - \zeta
$$

### 式 (3.101)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.101)。

**Q2** $\gamma = \exp\big[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}\big]$。

**Q3** „mit“

**Q4** 式 (3.100)。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
\gamma = \exp\big[\eta\sqrt{-\Delta}-(t-t_0)\sqrt{-\Delta}\big]
$$

### 式 (3.102)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.102)。

**Q2** $\zeta$ 只含初值量 $\eta$、$\Delta$、$A$、$B$，不含 $(t-t_0)$。对数内是 $\exp(\eta\cdot\sqrt{-\Delta})+1$。

**Q3** „und“

**Q4** 式 (3.100) 的积分常数。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
\zeta
  = -\frac{1}{A}\ln\big[\exp(\eta\cdot\sqrt{-\Delta})+1\big]
    - \frac{(-0{,}5\cdot\sqrt{-\Delta}-0{,}5\cdot B)}{A}\,\eta
$$

### 式 (3.103)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.103)。

**Q2** $h_P(u_P)$。对数项系数是 $+1/(2A)$，与 Fall 1 式 (3.96) 的 $-1/(2A)$ 不同。第二项系数是 $-B/(2A\sqrt{-\Delta})$。分母里的因子顺序是 $(\Psi_0+\sqrt{-\Delta})\,(\Psi_0-\sqrt{-\Delta})$。

**Q3** „Die Gleichung der Partikelhöhe über der Wirbelschichtoberfläche in Abhängigkeit der Geschwindigkeit lautet“。这里床面用词回到 Wirbelschichtoberfläche。

**Q4** Fall 2 的速度–高度关系。$\Psi$、$\Psi_0$ 见下两式。

**Q5** 无空白。

**Q6** 比对阶段填写。

$$
\begin{aligned}
h_P(u_P)
  &= h_{P,0}
    + \frac{1}{2A}\ln\left[
        \frac{A u_P^2+B u_P+C}{A u_{P,0}^2+B u_{P,0}+C}
      \right] \\
  &\quad - \frac{B}{2A\sqrt{-\Delta}}
      \ln\left[
        \frac{(\Psi-\sqrt{-\Delta})\,(\Psi+\sqrt{-\Delta})}
             {(\Psi_0+\sqrt{-\Delta})\,(\Psi_0-\sqrt{-\Delta})}
      \right]
\end{aligned}
$$

### 式 (3.104)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.104)。

**Q2** $\Psi = 2A\cdot u_{P,0} + B$。用的是起始速度 $u_{P,0}$。

**Q3** „mit“

**Q4** 式 (3.103)。

**Q5** 下标见式 (3.105)。

**Q6** 比对阶段填写。

$$
\Psi = 2\cdot A\cdot u_{P,0} + B
$$

### 式 (3.105)，S. 49

**Q1** PDF 034，S. 49 右栏，式 (3.105)。§3.2.4 到此结束。

**Q2** $\Psi_0 = 2A\cdot u_P + B$。用的是当前速度 $u_P$。下标 0 在当前速度这一侧。

**Q3** 与式 (3.104) 并列的定义。

**Q4** 式 (3.103)。

**Q5** 下标看起来与「0 表示初值」相反：$\Psi$ 含 $u_{P,0}$，$\Psi_0$ 含 $u_P$。式面如此，不调换。影响：若实现时按习惯把下标 0 当成初值，式 (3.103) 的分子分母会互换。标注：解析解的定义，照录。

**Q6** 比对阶段填写。

$$
\Psi_0 = 2\cdot A\cdot u_P + B
$$

S. 49 末段：这里给出的、由 Wirsum (1998) 建议的解析算法，使颗粒轨迹计算比 Runge–Kutta 数值积分快一个因子 50。因此可以把粒径、抛出初速和球形度的分级分得更细，从而更好地逼近真实固体浓度剖面。仍然没有浓度、回落量或带出量的式。

---

## 覆盖

| 式号 | 印刷页 | 本文件 |
|---|---|---|
| (3.49)、(3.50) | S. 36（与 S. 35 同号） | 有，带页码 |
| (3.51) | S. 36 | 有 |
| (3.52)– (3.55) | S. 37 | 有 |
| (3.56)– (3.60) | S. 38 | 有 |
| (3.61)– (3.64) | S. 39 | 有 |
| (3.65)– (3.69) | S. 40 | 有 |
| (3.70)– (3.72) | S. 41 | 有 |
| (3.73)，Bild 3.7 | S. 42 | 有 |
| (3.74)– (3.80) | S. 43 | 有 |
| (3.81)– (3.84) | S. 44 | 有 |
| (3.85)– (3.87)，Tabelle 3.1 起 | S. 45 | 有 |
| (3.88)– (3.91)，Tabelle 3.1 续 | S. 46 | 有 |
| (3.92)、 (3.93)，Bild 3.8 | S. 47 | 有 |
| (3.94)– (3.98) | S. 48 | 有 |
| (3.99)– (3.105) | S. 49 | 有 |
