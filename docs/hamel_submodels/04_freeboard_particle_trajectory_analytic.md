# Freeboard 颗粒轨迹解析（最小注释版）

## 论文编号锚点

- Kapitel 3.2.2~3.2.4
- `(3.71)`, `(3.74)`, `(3.76)`
- `(3.85)`, `(3.88)`, `(3.89)`~`(3.92)`
- `(3.93)`~`(3.105)`

## 原文直引（OCR）

- `(3.71)`：`[可直接采信] ū_{P,0} = 7{,}92·10^{-6} · u_b · d_b^{1{,}5}`  (3.71)
- `[可直接采信]（Martens 1984）50% 质量分位对应 `u_{P,0,50} = 1{,}53 · u_b`
- `(3.74)`：`[可直接采信] u_{gb}(h)=u_0 +(u_{gb,0}-u_0)·e^{-β_A h}`  (3.74)
- `(3.76)`：`[可直接采信] u_{gb,0} = u_{b,WS}`
- `(3.85)`：`[可直接采信] ... Differentialgleichung ... Partikelbewegung ...`
- `(3.88)`：`[可直接采信] du_P/dt = A·u_P² + B·u_P + C`
- `(3.92)`：`[可直接采信] Δ = 4AC - B²`

## 从 PNG 图像精确转录的公式（五次复核升级）

以下公式从 `data/hamel_pages_png/page_33.png` 和 `page_34.png`（300 dpi 渲染）视觉辨识转录。


### Freiraum 气泡与颗粒起始条件 (3.70)~(3.76)

以下公式从 `data/hamel_pages_png/page_30.png` 与 `page_31.png` 视觉辨识转录。

$$u_{P,m}=0{,}092\cdot u_b\cdot d_b^{1{,}5}\cdot \left(\frac{3\,d_P}{d_b}\right)^{0{,}3}\cdot\sqrt{\frac{\rho_g}{\rho_P}} \qquad (3.70)$$

$$\bar{u}_{P,0}=7{,}92\cdot10^{-6}\cdot u_b\cdot d_b^{1{,}5} \qquad (3.71)$$

$$u_{gb}(h)=u_0+(u_{gb,0}-u_0)\cdot e^{-\beta_A h} \qquad (3.74)$$

$$\beta_A=\frac{1}{3\cdot d_{b,WS}} \qquad (3.75)$$

$$u_{gb,0}=u_{b,WS} \qquad (3.76)$$

### Riccati 型颗粒运动方程 (3.88)

$$\frac{du_P}{dt} = A \cdot u_P^2 + B \cdot u_P + C \qquad (3.88)$$

### 系数 A, B, C (3.89)~(3.91)

$$A = -(+)\frac{3}{4} \cdot z \cdot \frac{\rho_g}{\rho_P \cdot d_P} \qquad (3.89)$$

$$B = -\frac{18 \cdot \mu_g}{\rho_P \cdot d_P^2} + (-)\frac{3}{2} \cdot z \cdot \frac{\rho_g}{\rho_P \cdot d_P} \cdot u_g \qquad (3.90)$$

$$C = \frac{18 \cdot \mu_g}{\rho_P \cdot d_P^2} \cdot u_g - (+)\frac{3}{4} \cdot z \cdot \frac{\rho_g}{\rho_P \cdot d_P} \cdot u_g^2 - \frac{(\rho_P - \rho_g)}{\rho_P} \cdot g \qquad (3.91)$$

> **注**：括号内的 (+)/(-) 符号取决于颗粒相对于气体的运动方向：
> - 颗粒速度 **小于** 气体速度或颗粒下落 → 取 **括号内** 符号
> - 颗粒速度 **大于** 气体速度 → 取 **非括号内** 符号

### 判别式 (3.92)

$$\Delta = 4AC - B^2 \qquad (3.92)$$

### Fall 1: Δ > 0

#### 颗粒速度 (3.93)

$$u_P(t) = \frac{\sqrt{\Delta}}{2A} \cdot \tan\!\left[\frac{1}{2}\,\omega\sqrt{\Delta} + \frac{1}{2}(t - t_0)\sqrt{\Delta}\right] - \frac{B}{2A} \qquad (3.93)$$

#### 辅助量 ω (3.94)

$$\omega = \frac{2}{\sqrt{\Delta}} \cdot \arctan\!\left[\frac{2A\,u_{P,0} + B}{\sqrt{\Delta}}\right] \qquad (3.94)$$

#### 颗粒高度（时间函数）(3.95)

$$h_P(t) = h_{P,0} - \frac{1}{A}\ln\!\left[\cos\!\left(-\frac{1}{2}\omega\sqrt{\Delta} - \frac{1}{2}(t-t_0)\sqrt{\Delta}\right)\right] - \frac{B}{2A}(t-t_0) + \frac{1}{A}\ln\!\left[\cos\!\left(\frac{1}{2}\omega\sqrt{\Delta}\right)\right] \qquad (3.95)$$

#### 颗粒高度（速度函数）(3.96)

$$h_P(u_P) = h_{P,0} - \frac{1}{2A}\ln\!\left[\frac{Au_P^2 + Bu_P + C}{Au_{P,0}^2 + Bu_{P,0} + C}\right] - \frac{B}{A\sqrt{\Delta}}\arctan\!\left[\frac{2Au_P + B}{\sqrt{\Delta}}\right] + \frac{B}{A\sqrt{\Delta}}\arctan\!\left[\frac{2Au_{P,0} + B}{\sqrt{\Delta}}\right] \qquad (3.96)$$

#### 颗粒极限直径 (3.97)

$$d_P^{*} = \left(\frac{108\,\mu_g^2}{c_W \cdot \rho_P \cdot \rho_g - \rho_g^2}\right)^{1/3} \qquad (3.97)$$

> Δ > 0 要求颗粒速度大于气体速度且颗粒直径大于极限直径 d_P*。

### Fall 2: Δ < 0

#### 颗粒速度 (3.98)

$$u_P(t) = \frac{(\sqrt{-\Delta} - B)\exp\!\left[\eta\sqrt{-\Delta} - (t-t_0)\sqrt{-\Delta}\right] - \sqrt{-\Delta} - B}{2A\exp\!\left[\eta\sqrt{-\Delta} - (t-t_0)\sqrt{-\Delta}\right] - 2A} \qquad (3.98)$$

#### 辅助量 η (3.99)

$$\eta = \frac{1}{\sqrt{-\Delta}} \ln\!\left(\frac{\sqrt{-\Delta} + (2A\,u_{P,0} + B)}{\sqrt{-\Delta} - (2A\,u_{P,0} + B)}\right) \qquad (3.99)$$

#### 颗粒高度（时间函数）(3.100)

$$h_P(t) = h_{P,0} - \frac{1}{A}\ln|\gamma + 1| - \frac{-0{,}5\cdot\sqrt{-\Delta}\cdot 0{,}5\cdot B}{A}\,(\eta - (t - t_0)) - \varsigma \qquad (3.100)$$

其中：

$$\gamma = \exp\!\left[\eta\sqrt{-\Delta} - (t - t_0)\sqrt{-\Delta}\right] \qquad (3.101)$$

$$\varsigma = -\frac{1}{A}\ln\!\left[\exp(\eta\cdot\sqrt{-\Delta}) + 1\right] - \frac{-0{,}5\cdot\sqrt{-\Delta} - 0{,}5\cdot B}{A}\,\eta \qquad (3.102)$$

#### 颗粒高度（速度函数）(3.103)

$$h_P(u_P) = h_{P,0} + \frac{1}{2A}\ln\!\left[\frac{Au_P^2 + Bu_P + C}{Au_{P,0}^2 + Bu_{P,0} + C}\right] - \frac{B}{2A\sqrt{-\Delta}}\ln\!\left[\frac{(\Psi - \sqrt{-\Delta})\,(\Psi_0 + \sqrt{-\Delta})}{(\Psi + \sqrt{-\Delta})\,(\Psi_0 - \sqrt{-\Delta})}\right] \qquad (3.103)$$

其中：

$$\Psi = 2 \cdot A \cdot u_{P,0} + B \qquad (3.104)$$

$$\Psi_0 = 2 \cdot A \cdot u_P + B \qquad (3.105)$$

## 参数选择与来源口径（补充）

- `u_{P,0}` 起始速度：`(3.70)`（Demnich & Bohnet）与 `(3.71)`（Son et al.）均在原文给出；论文同时引用 Martens (1984) 的 `u_{P,0,50}=1.53u_b` 作为分布口径对照。
- `u_{gb}(h)` 与 `\beta_A`：采用 Ghost-bubble 路线 `(3.74)~(3.76)`，其中 `u_{gb,0}=u_{b,WS}`。
- 阻力相关参数 `z`：来自 Haider & Levenspiel 关联链（通过 `(3.85)` 近似后进入 `(3.89)~(3.91)`），本文档不引入论文外常数替代。
- 分支判据：严格使用 `(3.92)` 的 `\Delta=4AC-B^2` 进行 `\Delta>0` / `\Delta<0` 两分支求解，不扩展论文外统一模板。

## 一致性约束

- 主判据按原文 `Δ > 0 / Δ < 0` 两类（Fall 1 / Fall 2）。
- 不写论文外"Riccati 三分支统一模板"。
- 非引文文本统一使用 `u_{P,0}` 记号；OCR 引文保持原样。
- 系数 A, B, C 中的 (+)/(-) 括号符号为原文表述方式，表示运动方向判定。

## PDF 页锚点

- `===== PAGE 30 =====`~`===== PAGE 34 =====`
- 公式 (3.70)~(3.76) -> PAGE 30~31（page_30.png, page_31.png）
- 公式 (3.88)~(3.92) -> PAGE 33（page_33.png）
- 公式 (3.93)~(3.97) -> PAGE 33~34（page_33.png, page_34.png）
- 公式 (3.98)~(3.105) -> PAGE 34（page_34.png）

## 修改标注（2026-05-07）

- `[M04-1]` 新增“参数选择与来源口径”段，补齐 `u_{P,0}`、`u_{gb}`、`z` 与分支判据的来源链。
- `[M04-2]` 明确 freeboard 轨迹模型的分支策略仅依赖 `(3.92)`，避免论文外求解分支扩展。
