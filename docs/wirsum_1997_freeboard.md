# Wirsum 1997 exact 自由板核

实现：`src/physics/wirsum_1997_freeboard.py`。生产 `simulate_freeboard` 仍走
`analytical_wirsum`（1.53 `u_b`、Roof+Wake、renormalize、Haider/Hamel 阻力混合），
本核尚未接入全炉。

## 模型名

| `freeboard_model` | 平均抛射 | 床面抛射 | 阻力 |
|---|---|---|---|
| `wirsum_1997_exact` | `1.2 u_b` | 仅 Roof Eq. 3.59，不 renormalize | `c_w = 24/Re + z`，`z ∈ [0.43, 1.0]` |
| `hamel_extended` | `1.53 u_b`（Martens / Hamel 3.73） | Roof + Wake，不 renormalize | Hamel `z(φ)`（未在本核展开） |
| `engineering_freeboard` | `1.53 u_b` | Roof + Wake，再缩放到 available feed | Haider 混合 |
| `analytical_wirsum` | 别名 → `engineering_freeboard` | 同上 | 同上 |

`exact_hamel` / Haider `z(φ)` 不是 `wirsum_1997_exact`。

## 方程

抛射速度（Wirsum §3.4）：

`u_{p,0} = 1.2 u_b`，竖直速度高斯 `σ = 0.6 m/s`。

Roof（Wirsum Eq. 3.59，与 Hamel Eq. 3.53 同构）[kg/(m² s)]：

`F_{0,r} = d_p ρ_p (1 − ε_d) (3 ε_b u_b / d_b)`

轨迹：

`c_w = 24/Re + z`，`du_p/dt = A u_p² + B u_p + C`，`Δ = 4AC − B²`。

`z` 取论文高 Re 极限常数，默认 0.44；不用 `73.69 exp(-5.0784 φ)`。
