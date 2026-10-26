# T26d 验证

失败点取 `reports/t26c-diagnosis/baseline/failure.npz` 与同目录 `frozen.pkl`。雅可比与牛顿步按 ADR 0041 重算。方向检查的 h 与 T26c 第 6 项相同。
小步长指 h ≤ 1e-08。这些 h 里至少有一个可行且 g(h) < g(0)，才记为小步长上的下降方向。

重算 F 与保存 F 的归一差 = 0.000000e+00。残差评价 499 次，用时 62.6 s。
新牛顿步范数 9.445573e+03，保存的牛顿步范数 1.008433e+04。

## 方向

g(0) = 1.211167e+01。g'(0) = -1.211167e+01。

| h | 状态 | g(h) | ‖JΔx − [F(x+hΔx)−F(x)]/h‖ / ‖JΔx‖ | 原因 |
|---|---|---|---|---|
| 1e-01 | reject | — | — | m < 0（下标 0） |
| 1e-02 | reject | — | — | m < 0（下标 0） |
| 1e-03 | reject | — | — | m < 0（下标 45） |
| 1e-04 | reject | — | — | m < 0（下标 52） |
| 1e-05 | reject | — | — | m < 0（下标 44） |
| 1e-06 | ok | 1.211166e+01 | 1.529111e-02 |  |
| 1e-07 | ok | 1.211166e+01 | 1.529198e-02 |  |
| 1e-08 | ok | 1.211167e+01 | 1.529208e-02 |  |
| 1e-09 | ok | 1.211167e+01 | 1.529322e-02 |  |
| 1e-10 | ok | 1.211167e+01 | 1.528616e-02 |  |
| 1e-11 | ok | 1.211167e+01 | 1.571495e-02 |  |
| 1e-12 | ok | 1.211167e+01 | 2.874147e-02 |  |

全部已评价 h 中 g(h) < g(0) 的步长：1e-06、1e-07、1e-08、1e-09、1e-10、1e-11、1e-12。相对差小于 1e-02 的步长：无。
线性化斜率是负的，与 −g(0) 相同到表中的位数，所以 Δx 是线性模型的下降方向。有可行 h 使 g(h) 低于 g(0)。差分与 J 在所测的 h 上都没有达到相对差 1e-2。
小步长上状态为 ok 的 h：1e-08、1e-09、1e-10、1e-11、1e-12。其中 g(h) < g(0) 的 h：1e-08、1e-09、1e-10、1e-11、1e-12。
小步长上是下降方向：True。

## 悬停焦炭与悬停灰

列是失败点上实际存在的悬停未知量。新步长为 √ε·|x|。旧步长为 √ε·max(|x|, s)。
1×、10²×、10⁴× 是新步长的倍数，相对差以新步长的差分列为分母。完整列在 `reports/t26d-validation/point.npz` 的 `hover_columns`。

| 列 | x | 新步长 | 旧步长 | 新列范数 | 与保存列的相对差 | 1× | 10²× | 10⁴× |
|---|---|---|---|---|---|---|---|---|
| freeboard[b=7,geo=0].char.hover.k3 | 1.412145e-11 | 2.104260e-19 | 5.905204e-08 | 7.132233e+07 | 4.181727e+03 | 0.000000e+00 | 1.460417e-06 | 1.489595e-04 |
| freeboard[b=7,geo=0].ash.hover.k3 | 3.689982e-06 | 5.498502e-14 | 5.905204e-08 | 5.458979e+02 | 7.938122e-03 | 0.000000e+00 | 9.846588e-07 | 7.449380e-05 |
| freeboard[b=6,geo=1].char.hover.k3 | 8.774153e-08 | 1.307451e-15 | 5.905204e-08 | 2.285006e+00 | 6.712782e-01 | 0.000000e+00 | 8.943131e-03 | 8.985995e-03 |
| freeboard[b=6,geo=1].ash.hover.k3 | 4.101862e-06 | 6.112250e-14 | 5.905204e-08 | 1.711865e-01 | 4.130708e-03 | 0.000000e+00 | 5.391635e-04 | 5.169157e-04 |

- freeboard[b=7,geo=0].char.hover.k3 绝对值最大的分量：freeboard[b=7,geo=0].ash.hover.k3 = 5.252817e+07；freeboard[b=7,geo=0].ash.free.k1 = -4.803633e+07；freeboard[b=7,geo=0].ash.free.k2 = -4.491833e+06；bed[b=3,geo=4].gas_outflow.suspension.Ar = 0.000000e+00；bed[b=3,geo=4].gas_outflow.suspension.H2S = 0.000000e+00
- freeboard[b=7,geo=0].ash.hover.k3 绝对值最大的分量：freeboard[b=7,geo=0].ash.hover.k3 = -4.020482e+02；freeboard[b=7,geo=0].ash.free.k1 = 3.676679e+02；freeboard[b=7,geo=0].ash.free.k2 = 3.438029e+01；freeboard[b=7,geo=0].gas_outflow.gas.CO = 1.231674e-01；freeboard[b=7,geo=0].gas_outflow.gas.CO2 = -6.486479e-02
- freeboard[b=6,geo=1].char.hover.k3 绝对值最大的分量：freeboard[b=6,geo=1].ash.hover.k3 = 1.681389e+00；freeboard[b=6,geo=1].ash.free.k1 = -1.537609e+00；freeboard[b=6,geo=1].ash.free.k2 = -1.437804e-01；freeboard[b=6,geo=1].gas_outflow.gas.H2 = 6.368632e-02；freeboard[b=6,geo=1].gas_outflow.gas.H2O = -5.837913e-02
- freeboard[b=6,geo=1].ash.hover.k3 绝对值最大的分量：freeboard[b=6,geo=1].gas_outflow.gas.CO = 1.235145e-01；freeboard[b=6,geo=1].ash.hover.k3 = -7.272582e-02；freeboard[b=6,geo=1].ash.free.k1 = 6.650682e-02；freeboard[b=6,geo=1].gas_outflow.gas.CO2 = -6.493594e-02；freeboard[b=6,geo=1].gas_outflow.gas.H2O = 6.357365e-03

## 外层第 1 轮

`outer_round_limit = 1` 只用于取出第 1 轮，不改默认 100。诊断记录在 `reports/t26d-validation/round1/record.json`。
求解状态 NOT_CONVERGED。原因：雅可比不能求解。内层状态 NOT_CONVERGED。外层 1 轮，内层迭代 1 次，墙钟 241.193253167 s。
停止，不重跑。第 1 轮终止状态为 NOT_CONVERGED：雅可比不能求解。

## 初值点的雅可比

外层第 1 轮在内层第 1 次迭代停止，阻尼试探次数为 0，牛顿步没有形成。同一初值上重算：未知量 498 个，其中 24 个为 0；新雅可比最大绝对值 1.04e+18。有 9 列恰好为 0，因为 F(x+√ε|x|) 与 F(x) 相同。

| 列 | x | 新步长 | 旧步长 |
|---|---|---|---|
| 自由空间 8 格的 C10H8 | 1.70e-20 到 2.46e-20 | 2.53e-28 到 3.67e-28 | 1.78e-08 |
| bed[b=7,geo=0].gas_outflow.suspension.O2 | 1.64e-23 | 2.45e-31 | 1.78e-08 |

对这张新雅可比，不缩放的 `solve(J, −F)` 与行列缩放后的求解都返回 Singular matrix。自由空间前两格 C10H8 用旧步长算出的列范数是 3.53e+05 和 4.99e+05。

失败点上的方向检查没有因小步长停止。第 1 轮也不是 DAMPING_FAILED。第 1 轮不能形成牛顿步，所以没有进入基准与 S1–S16。求解规则没有再改。
