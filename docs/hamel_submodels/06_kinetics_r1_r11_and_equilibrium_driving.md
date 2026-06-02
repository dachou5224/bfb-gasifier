# 动力学与平衡驱动力（最小注释版）

## 论文编号锚点

- Kapitel 5
- `(5.45)`, `(5.46)`, `(5.47)`, `(5.48)`, `(5.49)`
- `(5.50)`, `(5.51)`
- `(5.52)`~`(5.57)`
- `(5.58)`, `(5.59)`
- `(5.60)`, `(5.61)`
- `Tabelle 5.4`

## 原文直引（OCR）

- `(5.46)`：`[可直接采信] y_{eq,CO} = y_{CO2}·y_{H2} / (k_{eq,R8}·y_{H2O})`
- `(5.47)`：`[可直接采信] k_{eq,R8} = exp(-36.893 + 4.019/T_g)`  (5.47)
- `(5.60)`：`[可直接采信] k_{10} = 59,8 ... E/R_g = 12.200 K`（Olefine）
- `(5.61)`：`[可直接采信] k_{10} = 20.700 ... E/R_g = 9.650 K`（Aromaten）
- `[可直接采信] Tabelle 5.4: ... Teeren ...`

## 从 PNG 图像精确转录的公式（四次复核升级）

以下公式从 `data/hamel_pages_png/page_57.png`（300 dpi 渲染）视觉辨识转录。

### 5.2.4 均相水煤气变换反应（WGSR, R8）

#### 反应速率 (5.45)

$$R_{CO} = a_{R8} \cdot 2{,}77 \cdot 10^8 \cdot (y_{CO} - y_{eq,CO}) \cdot \exp\!\left(\frac{-13{.}971}{T}\right) \cdot p^{\left(\frac{0{,}5 - p}{250}\right)} \cdot \exp\!\left(-8{,}91 + \frac{5{,}553}{T}\right) \qquad (5.45)$$

> **注**: 此公式从高分辨率裁图精确转录。压力修正项按原式保持为 `p^((0,5-p)/250)`，不做论文外单位推断。
> 原文说明："Die Korrektur des Kohlenmonoxidanteils in Gleichung (5.45) berücksichtigt die Gleichgewichtslage für diese Reaktion."

#### 平衡摩尔分数 (5.46)

$$y_{eq,CO} = \frac{y_{CO_2} \cdot y_{H_2}}{k_{eq,R8} \cdot y_{H_2O}} \qquad (5.46)$$

#### 平衡常数 (5.47)

$$k_{eq,R8} = \exp\!\left(-36{,}893 + 4{,}019 / T_g\right) \qquad (5.47)$$

> **⚠️ 物理一致性警告 [M06-3]**：OCR 转录值 `-36,893`（= -36.893）在物理上不合理。
> 对 `(5.46)` 的角色分析：`k_eq_R8` 必须等于 WGSR 平衡常数 $K_{eq}$（CO + H₂O → CO₂ + H₂），在 1000 K 时应约等于 0.9。
> 用 OCR 值代入：$\exp(-36.893 + 4.019/1000) \approx 10^{-16}$，完全不合理。
>
> **最可能的 OCR/德语数字格式错误**：
> - 德语原文中句号`.`用作千分符，`4.019` 实为 **4019**（四千零十九）；
> - 德语原文中逗号`,`用作小数点，`-3,6893` 被 OCR 误识为 `-36.893`（缺失小数点，数量级偏差 10×）。
>
> **物理上合理的修正形式**（须凭原始 PDF 图像最终确认）：
> $$k_{eq,R8} = \exp\!\left(-4{,}139 + 4019 / T_g\right)$$
> 此式在 700 K、1000 K、1200 K 分别给出 $K_{eq} \approx 5.0,\; 0.89,\; 0.45$，与 van't Hoff 热力学计算值吻合。
>
> **当前 OCR 转录值保留以便与原文对照，实现代码中禁止直接使用此式，须在 image-level 核实后更新。**

#### 煤灰催化因子 (5.48)

$$a_{R8} = 0{,}02 \qquad (5.48)$$

> "Als Faktor, der die katalytische Aktivität der Kohlenasche bzw. des Kokses beinhaltet, verwenden Wen und Chaung für die eigene Simulation eines Flugstrom-Kohlevergasers"

#### SI 单位版反应速率 (5.49)

$$R_{CO} = -8{,}82 \cdot 10^{-7} \cdot \exp\!\left(\frac{-90{.}853}{R_g \cdot T_g}\right) \cdot \left[p_{CO} - \frac{p_{H_2} \cdot p_{CO_2}}{k_{eq,R8} \cdot p_{H_2O}}\right] \quad \text{in}\;\frac{mol}{s \cdot kg_{Asche}} \qquad (5.49)$$

> 原文："Umgerechnet auf SI-Einheiten ergibt sich:"

### 5.2.5 甲烷重整（R9）

#### 反应方程

$$CH_4 + H_2O \;\to\; CO + 3\,H_2 \qquad +206{,}3\;\text{kJ/mol}\;\text{CH}_4 \qquad (5.50)$$

#### 反应速率 (5.51)

$$R_{CH_4} = -6{,}113 \cdot 10^{-2} \cdot \exp\!\left(\frac{-137{.}327}{R_g \cdot T}\right) \cdot p_{CH_4} \quad \text{in}\;\frac{mol}{s \cdot kg_{Asche}} \qquad (5.51)$$

> **可直接采信**（四次复核升级确认）。甲烷重整仅在含固体催化的悬浮相中计算。
> Chen et al. (1987) 实验数据，Wirbelbett aus Kohlenasche。

### 5.2.6 焦油二次反应

#### 焦油代表物氧化（R10）

$$C_6H_6 + 3\,O_2 \;\to\; 6\,CO + 3\,H_2 \qquad (5.52)$$
$$C_{10}H_8 + 5\,O_2 \;\to\; 10\,CO + 4\,H_2 \qquad (5.53)$$
$$C_{16}H_{34} + 8\,O_2 \;\to\; 16\,CO + 17\,H_2 \qquad (5.54)$$

#### 焦油水蒸气重整（R11）

$$C_6H_6 + 9/2\,H_2O \;\to\; 9/2\,CO + 9/2\,H_2 + 3/2\,CH_4 \qquad (5.55)$$
$$C_{10}H_8 + 23/3\,H_2O \;\to\; 23/3\,CO + 7\,H_2 + 7/3\,CH_4 \qquad (5.56)$$
$$C_{16}H_{34} + 21/2\,H_2O \;\to\; 21/2\,CO + 33/2\,H_2 + 11/2\,CH_4 \qquad (5.57)$$

#### 焦油全局氧化反应方程 (5.58)

$$C_mH_n + \frac{m}{2}\,O_2 \;\to\; m\cdot CO + \frac{n}{2}\cdot H_2 \qquad (5.58)$$

#### 焦油氧化反应速率 (5.59)

$$R_{C_mH_n} = \frac{dC_{C_mH_n}}{dt} = -k_{10} \cdot \exp\!\left(\frac{-E_a}{R_g \cdot T}\right) \cdot T \cdot p^{0{,}3} \cdot C_{C_mH_n}^{0{,}5} \cdot C_{O_2} \qquad (5.59)$$

> **可直接采信**（从 page_58.png 清晰读出）。
> "Die Kinetik dieser Reaktion berechnet sich nach Siminski et al. (1972) zitiert in Smoot und Smith (1985) zu:"
> 来源：Siminski et al. (1972)，zitiert in Smoot und Smith (1985)。

#### 焦油氧化动力学参数 (5.60), (5.61)

Olefine（脂肪族烃）:
$$k_{10} = 59{,}8 \quad \frac{1}{mol^{0{,}5} \cdot K \cdot Pa^{0{,}3} \cdot s} \qquad E/R_g = 12{.}200\;\text{K} \qquad (5.60)$$

Aromaten（芳香族烃）:
$$k_{10} = 20{.}700 \quad \frac{1}{mol^{0{,}5} \cdot K \cdot Pa^{0{,}3} \cdot s} \qquad E/R_g = 9{.}650\;\text{K} \qquad (5.61)$$

## 参数选择与来源口径（补充）

- WGSR（R8）采用 `(5.45)` 非 SI 经验式与 `(5.49)` SI 形式并列口径，`a_{R8}` 取 Wen & Chaung 给定值 `(5.48)`。
- 平衡驱动力采用 `(5.46)` 与 `(5.47)`，即 `y_{CO}` 与 `y_{eq,CO}` 差值驱动。
- 焦油氧化速率采用 `(5.59)` 结构，压强指数为 `p^{0.3}`；对应 `k_{10}` 单位须与 `Pa^{0.3}` 一致。
- `k_{10}, E/R_g` 参数值来自 `(5.60)/(5.61)`，分别对应 Olefine 与 Aromaten。

## 一致性约束

- (5.45) 已由高分辨率裁图逐字符确认；压力项按原式保留，不引入论文外单位解释。
- 不写论文外矩阵化源项表达。
- 非引文文本统一使用 `K_{eq,R8}`、`y_{eq,CO}`、`T_g`；OCR 引文保持原样。

## PDF 页锚点

- `===== PAGE 56 =====`（`(5.36)` 所在段）
- `===== PAGE 57 =====`（`(5.45)(5.46)(5.47)(5.48)(5.49)(5.50)(5.51)` 所在段，page_57.png）
- `===== PAGE 58 =====`（`(5.52)`~`(5.61)`, `Tabelle 5.4` 所在段）

## 修改标注（2026-05-07）

- `[M06-1]` 修正 `(5.60)/(5.61)` 的 `k_{10}` 压强量纲为 `Pa^{0.3}`，与 `(5.59)` 压强指数一致。
- `[M06-2]` 新增“参数选择与来源口径”段，明确 R8 与焦油二次反应的参数来源与单位闭合关系。

## 修改标注（2026-05-07，第二轮审计）

- `[M06-3]` 新增 `(5.47)` 物理一致性警告：OCR 转录值 `exp(-36,893 + 4,019/T_g)` 在物理上不合理（1000 K 时给出 K_eq≈10^{-16}，WGSR 物理期望值约 0.9）。德语数字格式分析：`4.019`（句号为千分符）= 4019；`-36.893` 疑为 OCR 对德语小数 `-3,6893`（= -3.6893）的误识（数量级偏差 10×）。物理合理修正形式约为 `exp(-4.139 + 4019/T_g)`，须凭原始 PDF 图像核实后方可在代码中使用。
