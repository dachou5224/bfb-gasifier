# 床层水力学主链（最小注释版）

## 论文编号锚点

- Kapitel 3.1.2~3.1.4
- `(3.11)`, `(3.12)`, `(3.23)`, `(3.24)`, `(3.33)`
- `(3.41)`, `(3.42)`, `(3.43)`, `(3.44)`
- `(3.46)`~`(3-50)`

## 原文直引（OCR）

- `[可直接采信] ... ub,r = nb·ud ... (3.23)`（含 `nb = 2,7`）

## 从 PNG 图像精确转录的公式（四次复核升级）

以下公式从 `data/hamel_pages_png/page_25.png` 和 `page_27.png`（300 dpi 渲染）视觉辨识转录。

### 3.1.3 高压修正（Heinbockel 链）

原文段落（PAGE 25 右栏）：
> "Zur Berechnung der Strömungsstruktur in der blasenbildenden Wirbelschicht einer unter Druck arbeitenden Wirbelschichtfeuerung oder eines Wirbelschichtvergasers reichen die in Kapitel 3.1 angegebenen Gleichungen nicht aus, da hier der Druck lediglich über die Gasdichte und damit über die minimale Fluidisierungsgeschwindigkeit einfließt."

#### 气泡直径增长率（压力修正）(3.41)

$$\frac{dd_b}{dh} = \frac{\left(\frac{2}{9\pi}\right)^{1/3} \varepsilon_b^{1/3} \left(\frac{P}{P_0}\right)^{P_0/p}}{\left[1 - \varepsilon_b \cdot \left(\frac{P}{P_0}\right)^{1/3} \cdot \varepsilon_b^{1/3}\right]} \cdot \frac{d_b}{3 \cdot \lambda_b \cdot u_b} \qquad (3.41)$$

> **注**: 此公式由 Heinbockel (1995) 对 Hilligardt (1986) 的 (3.33) 进行压力修正得到。
> 五次高分辨率裁图复核确认：指数为 `(P/P_0)^(P_0/p)`。

#### 气泡生存时间修正（压力修正）(3.42)

$$\lambda_b = 280 \cdot \frac{u_{mf}}{g} \cdot \left(\frac{P}{P_0}\right)^{-0{,}7} \qquad (3.42)$$

> **可直接采信**。原文："Die Gleichung (3.35) lautet in modifizierter Form:"
> 来源：Heinbockel (1995) 对 Olowson 和 Almstedt (1990) 实验数据的拟合。

#### 气泡上升速度（压力修正）(3.43)

$$u_b = \psi_b \cdot (u_0 - u_{mf}) \left[\left(\frac{P}{P_0}\right)^{0{,}2} - 1\right] + v_b \cdot u_{b,l} \qquad (3.43)$$

> 原文将压力修正通过 Überschussgeschwindigkeit（过剩速度）前的系数引入 (3.13)。
> 其中 $\psi_b$、$v_b$ 为模型参数。

#### 气泡通过速度（压力修正）(3.44)

$$u_{b,r} = n_b \cdot u_d \cdot \left(\frac{P}{P_0}\right)^{-0{,}15} \qquad (3.44)$$

> **可直接采信**（从 page_26.png 清晰读出）。
> Heinbockel (1995) 对 (3.23) 的压力修正：气泡相对穿透速度随压力升高而降低（指数 -0.15）。

### 3.1.4 气泡-悬浮相传质系数

#### Kunii & Levenspiel 模型 (3.46)~(3.47)

$$K_{bd} = \frac{1}{\frac{1}{K_{bc}} + \frac{1}{K_{cd}}} \qquad (3.46)$$

其中：

$$K_{bc} = 4{,}5 \cdot \frac{u_{mf}}{d_b} + 5{,}85 \cdot \frac{D_g^{0{,}5} \cdot g^{0{,}25}}{d_b^{1{,}25}} \qquad K_{cd} = 6{,}78 \cdot \frac{\sqrt{D_g \cdot u_b \cdot \varepsilon_{mf}}}{d_b^{3/2}} \qquad (3.47)$$

> Kunii und Levenspiel (1977)：分别为 cloud-bubble 和 cloud-suspension 传质阻力的串联模型。

#### Sit & Grace 合并公式 (3.48)

$$K_{bd} = \frac{2 \cdot u_{mf}}{d_b} + \sqrt{\frac{144 \cdot D_g \cdot \varepsilon_{mf} \cdot u_b}{\pi \cdot d_b^3}} \qquad (3.48)$$

> Sit und Grace (1981)，适用于合并气泡（koaleszierende Blasen）。

#### Preto 修正公式 (3.49)

$$K_{bd} = 4{,}5 \cdot \chi \cdot \frac{u_{mf}}{d_b} \qquad (3.49)$$

> 参数 χ 考虑气泡周围流场：慢速穿透气泡 χ = 1；快速气泡最小约 χ = 0,1（Preto 1986）。

#### Sit & Grace (1981) 含 Blasendurchströmung 修正 (3.50)

$$K_{bd} = \frac{u_{b,r}}{V_b} + \sqrt{\frac{144 \cdot D_g \cdot \varepsilon_{mf} \cdot u_b}{\pi \cdot d_b^3}} + \frac{3}{2} \cdot \frac{u_b}{d_b} \cdot \sqrt{\frac{144 \cdot D_g \cdot \varepsilon_{mf} \cdot u_b}{d_b^3}} \qquad (3.50)$$

> **注**: 此公式从 page_27.png 右上角视觉辨识。图像显示三项之和，分别对应：
> 1. 气泡穿透（对流）分量 u_{b,r}/V_b
> 2. 扩散分量（含 π）
> 3. 高 K_bd 修正项（含 3/2 系数）
>
> 原文文段："Mit der Abschätzung des diffusiven Anteils am Stoffaustausch nach Sit und Grace (1981) ergibt sich für den Stoffaustauschkoeffizienten:"
>
> Bild 3.6 对比了四种 K_bd 计算方法在不同压力下的预测值，(3.50) 在 900°C 下预测 K_bd ≈ 33 s⁻¹（小气泡）。

> **⚠️ 量纲一致性警告 [M03-4]**：第一项 `u_{b,r}/V_b` 的量纲为 `[m/s / m³] = [1/(m²·s)]`，与 K_bd 目标量纲 `[1/s]` 不一致。正确的对流分量形式应为 `u_{b,r}/d_b` 量级的表达式（量纲 `[1/s]`）。此问题可能来自图像转录时将气泡直径 `d_b` 误识为气泡体积 `V_b`，或 `V_b` 在原文语境中代表某个无量纲量。须凭 page_27.png 原始图像逐字符核实。

## 参数选择与来源口径（补充）

- 压力修正链 `(3.41)~(3.44)`：来源为 Heinbockel (1995) 对 Hilligardt (1986) 模型的扩展，目的为复现 Olowson & Almstedt (1990) 的高压实验趋势（见 `PAGE 25~26` 叙述）。
- `(3.43)` 中 `\psi_b`、`v_b`：论文文本仅给出“作为压力修正系数引入”的角色，本分册不外推固定数值；具体取值必须由论文对应校准口径或实验拟合给出。
- `(3.49)` 中 `\chi`：按 Preto (1986) 口径，慢速穿透泡取 `\chi=1`，快速泡 `\chi<1` 且下限约 `0.1`。
- `(3.50)`：在 `fitz/layout OCR` 中公式字符存在断裂；当前版本以 `page_27.png` 人工逐字符转录为主，并与该页后续文字说明（Kbd 高温外推）一致。

## 一致性约束

- (3.44) 已从 page_26.png 清晰确认：指数为 -0.15。
- (3.41) 已由高分辨率裁图确认：指数为 `(P_0/p)`。
- 不写论文外 `K_bd` 简化显式式。
- 非引文文本统一使用 `u_{b,r}`、`K_{bd}` 记号；OCR 引文保持原样。

## PDF 页锚点

- `===== PAGE 22 =====`（(3.23) 所在段）
- `===== PAGE 25 =====`（(3.41)(3.42)(3.43) 所在段，page_25.png）
- `===== PAGE 26 =====`（(3.44) 及 Bild 3.3~3.5）
- `===== PAGE 27 =====`（(3.45)~(3.50) 所在段，page_27.png）

## 修改标注（2026-05-07）

- `[M03-1]` 新增“参数选择与来源口径”段，明确 `(3.41)~(3.44)` 的文献归属与用途。
- `[M03-2]` 明确 `\psi_b`、`v_b`、`\chi` 的可用口径，避免论文外固定参数硬编码叙述。
- `[M03-3]` 补充 `(3.50)` 在 OCR 断裂条件下的转录依据（page_27 图像主证据）。

## 修改标注（2026-05-07，第二轮审计）

- `[M03-4]` 新增 `(3.50)` 量纲一致性警告：第一项 `u_{b,r}/V_b` 量纲为 `[1/(m²·s)]`，与 K_bd 目标量纲 `[1/s]` 不符。疑为图像转录时将 `d_b`（直径）误识为 `V_b`（体积），须凭 page_27.png 核实。
