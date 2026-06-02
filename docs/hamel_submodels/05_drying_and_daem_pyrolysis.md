# 干燥与热解（最小注释版）

## 论文编号锚点

- Kapitel 4.1
- `(4.1)`, `(4.2)`, `(4.3)`, `(4.4)`
- `(4-5)`~`(4.9)`
- `(4.10)`, `(4.11)`, `(4-12)`
- `Tabelle 4.1`

## 原文直引（OCR）

- `[可直接采信] ... Agarwal et al. (1986) ... "quasistationären" Modell ...`
- `[可直接采信] ... simultanen Trocknung und Entgasung ... "instationären" Modell ...`
- `[可直接采信] Tabelle 4.1 ... Aktivierungsenergie ... Standardabweichung ... Frequenzfaktor ...`

## 从 PNG 图像精确转录的公式（四次复核升级）

以下公式从 `data/hamel_pages_png/page_36.png`（300 dpi 渲染）视觉辨识转录。

### Fourier 球坐标热传导方程 (4.1)

$$\frac{\partial T}{\partial t} = \frac{a_s}{r^2}\frac{\partial}{\partial r}\!\left(r^2 \frac{\partial T}{\partial r}\right) \qquad \text{für}\; r_e \le r \le R_0 \qquad (4.1)$$

其中 $a_s = \lambda_s / (\rho_s \cdot c_{p,s})$ 为已干燥外壳的热扩散率。

> **来源**: PAGE 36 左栏，"Grundlage des Trocknungsmodells nach Agarwal et al. (1986) bildet die Fourier'sche Differentialgleichung der eindimensionalen Wärmeleitung mit konstanten Stoffwerten."

### 表面边界条件 (4.2)

$$\lambda_s \frac{dT}{dr}\bigg|_{r=R_0} = \alpha \cdot (T_a - T_s) = q(t) \qquad (4.2)$$

其中 $\alpha$ 为颗粒-床层传热系数，$T_a$ 为床层温度，$T_s$ 为颗粒表面温度。

### 蒸发前沿能量守恒 (4.3)

$$\lambda_s \frac{dT}{dr}\bigg|_{r=r_e} = h_v' \cdot w_{H_2O,P} \cdot \rho_s \cdot \frac{dr_e}{dt} \qquad (4.3)$$

其中 $r_e$ 为蒸发前沿半径，$w_{H_2O,P}$ 为含水质量分数。

### 修正蒸发焓 (4.4)

$$h_v' = h_v + \left(c_w + \frac{c_s}{w_{H_2O,P}}\right) \cdot (T_e - T_0) \qquad (4.4)$$

其中 $h_v$ 为水的蒸发焓，$T_e$ 为蒸发温度，$T_0$ 为初始温度。修正项包括将湿芯固体和水从 $T_0$ 加热到 $T_e$ 所需的焓。

### 初始与稳态边界条件 (4.6)

$$T\big|_{r=R_0} = T_s, \qquad T\big|_{r=r_e} = T_e \qquad (4.6)$$

其中 $T_s$ 为颗粒表面温度，$T_e$ 为蒸发前沿温度（取水的沸点）。此式给出稳态下外壳两侧的热边界条件。

> **来源**: OCR PAGE 36："Randbedingungen: $T|_{r=R_0}=T_s,\; T|_{r=r_e}=T_e$"

### Biot 数定义 (4.8)

$$Bi = \frac{\alpha \cdot R_0}{\lambda_s} \qquad (4.8)$$

其中 $\alpha$ 为颗粒-床层传热系数，$R_0$ 为颗粒半径，$\lambda_s$ 为外壳导热系数。Biot 数表征外部对流传热与内部导热的相对重要性；$Bi \gg 1$ 代表内部传热控制，$Bi \ll 1$ 代表外部传热控制。

> **来源**: OCR PAGE 36："$Bi = \alpha \cdot R_0 / \lambda_s$"

> **注**: 公式 (4.5)（无量纲坐标变换）和 (4.7)（无量纲温度变换）在 OCR 中字符断裂，无法完整转录；本文件仅收录可视觉确认的公式。须凭 page_36.png 图像补充。

### 颗粒-床层换热 Nu 关联式 (4.9)

$$Nu_P = 2 + 1{,}2 \cdot Re^{1/2} \cdot Pr^{1/3} \qquad (4.9)$$

> 按 Kunii and Levenspiel (1977) 的建议计算。

### DAEM 积分方程 (4.10)

$$\frac{m_v' - m_{v,\infty}'}{m_{v,0}'} = \int_0^\infty \exp\!\left[-k_0 \int_0^t \exp\!\left(\frac{-E}{R_g \cdot T}\right)dt'\right] f(E)\,dE \qquad (4.10)$$

其中 $m_v'$ 为尚未释放的挥发分质量分数，$k_0$ 为频率因子，$f(E)$ 为活化能分布函数。

### 活化能高斯分布 (4.11)

$$f(E) = \left[\sigma_E \cdot (2\pi)^{1/2}\right]^{-1} \cdot \exp\!\left[-\frac{(E - E_0)^2}{2\,\sigma_E^2}\right] \qquad (4.11)$$

其中 $E_0$ 为平均活化能，$\sigma_E$ 为标准差。

> 按 Anthony 和 Howard (1976) 的假设，所有平行一阶反应共用同一频率因子 $k_0$，活化能服从高斯分布。

### 颗粒平均挥发分含量（体积积分）(4-12)

$$\frac{m_v' - m_{v,\infty}'}{m_{v,0}'} = \frac{3}{R_0^3} \int_0^{R_0} \left\{\int_0^\infty \exp\!\left[-k_0 \int_0^t \exp\!\left(\frac{-E}{R_g \cdot T}\right)dt'\right] f(E)\,dE\right\} \cdot r^2\,dr \qquad (4\text{-}12)$$

> **可直接采信**（从 page_37.png 清晰读出）。
> 原文："Da die Temperatur T eine Funktion der Zeit t und des Radius r ist, gilt für den mittleren Gehalt an flüchtigen Bestandteilen des Partikels:"
> 将 DAEM 积分 (4.10) 进行球坐标体积平均，考虑颗粒内部径向温度分布 T = T(r, t)。

### Tabelle 4.1 褐煤热解参数（Dersch 1994）

| 参数 | 符号 | 值 | 单位 |
|------|------|------|------|
| 平均活化能 | $E_0$ | 192 000 | J/mol |
| 活化能标准差 | $\sigma_E$ | 40 000 | J/mol |
| 频率因子 | $k_{0,P}$ | 1,67·10¹³ | 1/s |
| 干燥颗粒密度 | $\rho_T$ | 1 250 | kg/m³ |
| 干燥颗粒比热 | $c_T$ | 1 256 | J/(kg·K) |
| 热扩散率 | $a$ | 0,1·10⁻⁶ | m²/s |

> 从 page_37.png Tabelle 4.1 精确转录。

## 参数选择与来源口径（补充）

- 干燥模型：采用 Agarwal et al. 路线，`(4.1)~(4.4)` + `(4.9)` 为控制方程与换热关联。
- 热解动力学：采用 Anthony & Howard 的 DAEM 口径，`(4.10)` 与 `(4.11)` 构成速率积分与活化能分布。
- 径向非均匀温度修正：按 `(4-12)` 对粒内体积分平均，不简化为论文外“均温颗粒”表达。
- `Tabelle 4.1` 参数来源：论文正文明确写明“参数取自 Dersch (1994)”，本分册仅复述该表，不新增外部拟合值。

## 一致性约束

- 不写论文外离散矩阵/积分模板。
- 参数口径以 `Tabelle 4.1` 原文为准。

## PDF 页锚点

- `===== PAGE 35 =====`~`===== PAGE 37 =====`（Agarwal 模型全文，(4.1)~(4.9)）
- `===== PAGE 36 =====`（(4.1)~(4.4) 与 (4.9)~(4.11) 所在段）
- `===== PAGE 37 =====`（`(4-12)` 与 `Tabelle 4.1` 所在段）
- 公式 (4.1)~(4.4) -> PAGE 36 左栏（page_36.png）
- 公式 (4.9)~(4.11) -> PAGE 36 右栏（page_36.png）

## 修改标注（2026-05-07）

- `[M05-1]` 新增“参数选择与来源口径”段，明确干燥/DAEM/体积分平均三层算法关系。
- `[M05-2]` 在页锚点中补充 `PAGE 37` 与 `Tabelle 4.1` 的对应关系。

## 修改标注（2026-05-07，第二轮审计）

- `[M05-3]` 新增公式 `(4.6)`（初始/稳态边界条件）和 `(4.8)`（Biot 数定义），两式均已在 OCR PAGE 36 中确认。同时增加说明：(4.5) 和 (4.7) 因 OCR 字符断裂未能完整转录，须凭 page_36.png 图像补充。
