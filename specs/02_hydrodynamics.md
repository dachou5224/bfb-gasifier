# 02 流体力学

**来源**：Hamel (1999)；Hilligardt 模型；Sit & Grace 混合模型

## 1. 最小流化速度 u_mf

Ergun 方程（Ar–Re_mf 关系）：
$$Ar = \frac{150(1-\epsilon_{mf})}{\phi_s \epsilon_{mf}^3} Re_{mf} + \frac{1.75}{\phi_s \epsilon_{mf}^3} Re_{mf}^2$$

实现模块：`src/physics/minimum_fluidization.py`

## 2. 气泡动力学

**气泡上升速度（Hilligardt）：**
$$u_b = \psi_b (u_0 - u_{mf}) + u_{b,i}, \quad \psi_b = 0.76 \text{（工业分布板默认）}$$

**气泡直径 ODE：**
$$\frac{dd_b}{dh} = \left[\frac{2}{9\pi} \cdot \frac{\epsilon_b^{1/3}}{1 - \xi_b(6/\pi)^{1/3}\epsilon_b^{1/3}}\right] - \frac{d_b}{3\lambda_b u_b}$$

**慢泡/快泡判别：**
$$\alpha_b = U_b / U_{mf}; \quad \alpha_b < 1 \Rightarrow \text{慢泡};\quad \alpha_b > 1 \Rightarrow \text{快泡}$$

实现模块：`src/physics/bubble_dynamics.py`

## 3. 相体积分率

n_RZ 必须按 Re_s 分段计算，禁止使用常数 4.65。

实现模块：`src/physics/phase_fractions.py`

## 4. 相间传质系数 K_bd（Gleichung 3.50，p. 35；Sit & Grace）

**Equation 3.50**（Hamel, **page 35**）定义气泡相与悬浮相之间的**总传质系数** \(K_{bd}\) [1/s]，为流体力学与化学摩尔衡算的**主接口**。形式为**对流**与**扩散**之和，适用于高压气化；与 **Kunii–Levenspiel** 等常数交换率模型不同，Hamel 强调 \(P\) 通过 **\(u_{b,r}\)**（**Gl. 3.44**）改变对流项，从而再现高压下传质减弱（如 **HTW** 加压炉）。

**等价写法（与论文一致）：**
$$K_{bd} = \frac{\dot{V}_{b,r}}{V_b} + \sqrt{\frac{144\, D_g\, \epsilon_{mf}\, u_b}{\pi\, d_b^3}} = \frac{3\, u_{b,r}}{2\, d_b} + \sqrt{\frac{144\, D_g\, \epsilon_{mf}\, u_b}{\pi\, d_b^3}}$$

- **对流项** \(\dfrac{3 u_{b,r}}{2 d_b}\)：穿流贡献；\(u_{b,r}\) 与 \(\dot{V}_{b,r}/V_b\) 一致。
- **扩散项** \(\sqrt{\cdots}\)：气泡界面分子扩散，基于 **Sit & Grace (1981)** 渗透理论。

**Gleichung 3.44（\(u_{b,r}\)，代码 `u_br`）：**
$$u_{b,r} = n_b \cdot u_d \cdot (P/P_0)^{-0.15}, \quad n_b = 2.7\ \text{（Gl. 3.23）},\; P_0 \approx 101{,}3\times10^3\,\text{Pa}$$

实现模块：`src/physics/mass_transfer.py`（`calc_u_br`、`calc_kbd`）

## 5. 自由板区

参数：`beta_A`、`u_gb`、`C_D_haider`

实现模块：`src/physics/freeboard.py`
