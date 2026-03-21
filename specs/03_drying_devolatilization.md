# 03 干燥与热解

## 1. 干燥（Hamel Kapitel 4，pp. 52–53）

**实现模块**：`src/thermal/drying.py`（Crank–Nicolson 径向离散；`corrected_evaporation_enthalpy`、`nusselt_particle` 等）

以下为论文中干燥子模型核心式（与 **Dissertation Chapter 4** 一致，便于检索）。

### Gleichung 4.2：外表面热流边界（\(r=R_0\)）

$$\lambda_s \left. \frac{\mathrm{d}T}{\mathrm{d}r} \right|_{r=R_0} = \alpha \cdot (T_{ws} - T_s) = \dot{q}(t)$$

- \(\lambda_s\)：干壳导热系数；\(\alpha\)：对流换热系数（**Gl. 4.9** Nu）；\(T_{ws}\)：床层温度；\(T_s\)：颗粒表面温度。

### Gleichung 4.4：修正蒸发焓

$$h_v' = h_v + \left( c_w + \frac{c_s}{w_{0,tr}} \right) \cdot (T_e - T_0)$$

- \(h_v\)：蒸发潜热；\(c_w\)、\(c_s\)：水与固体比热；\(w_{0,tr}\)：干基初始含水率；\(T_0\)→\(T_e\)：自初温至蒸发温度。

### Gleichung 4.6：干壳区温度边界

$$T|_{r=R_0} = T_s \quad \text{und} \quad T|_{r=r_e} = T_e$$

- \(R_0\)：颗粒外半径；\(r_e\)：蒸发前沿半径；\(T_e\)：蒸发（沸点）温度。

### Gleichung 4.9：颗粒 Nu（用于 \(\alpha\)）

见 `drying.py` 中 `nusselt_particle`；与 **Gl. 4.2** 中 \(\alpha\) 衔接。

NOTE: 解析闭式若不可用，采用 Crank–Nicolson FD 替代 Agarwal 纯解析路径。

---

## 2. 热解（DAEM）

采用 Gauss-Hermite 数值积分求解分布活化能模型。

**验证指标：** 升温至 900°C，积分 100 s，挥发分释放率 > 80%

实现模块：`src/thermal/devolatilization.py`（docs/CLAUDE.md Phase 4.2）
