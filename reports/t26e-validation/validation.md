# T26e 验证

求解手段：ADR 0042。失败点取 `reports/t26c-diagnosis/baseline/`。探针为 y + h·Δy，F 在 x = exp(y) 处评价。

重算 F 与保存 F 的归一差 = 0.000000e+00。雅可比残差评价 562 次，连同方向探针用时 68.5 s。对数变量 458 个。Δy 范数 6.371225e+06，保存的 x 空间步范数 1.008433e+04。被放大 55 列，到上限仍不可分辨 0 列。

## 失败点方向

g(0) = 1.211167e+01。g'(0) = -1.211167e+01。可行的最大 h = 1.000000e-05。该步 g = 1.211152e+01，下降量 g(0)−g(h) = 1.442579e-04。

| h | 状态 | g(h) | ‖JΔy − [F(y+hΔy)−F]/h‖ / ‖JΔy‖ | 原因 |
|---|---|---|---|---|
| 1e-01 | reject | — | — | x 非有限 |
| 1e-02 | reject | — | — | x 非有限 |
| 1e-03 | reject | — | — | x 非有限 |
| 1e-04 | domain | — | — | 试探点落在自由空间有效域外 |
| 1e-05 | ok | 1.211152e+01 | 4.229326e+00 |  |
| 1e-06 | ok | 1.211165e+01 | 1.144833e-01 |  |
| 1e-07 | ok | 1.211166e+01 | 1.093856e-02 |  |
| 1e-08 | ok | 1.211167e+01 | 4.146147e-03 |  |
| 1e-09 | ok | 1.211167e+01 | 4.263048e-03 |  |
| 1e-10 | ok | 1.211167e+01 | 4.297248e-03 |  |
| 1e-11 | ok | 1.211167e+01 | 4.110715e-03 |  |
| 1e-12 | ok | 1.211167e+01 | 3.156568e-02 |  |

g(h) < g(0) 的步长：1e-05、1e-06、1e-07、1e-08、1e-09、1e-10、1e-11、1e-12。

## 初值点

未知量 498 个，对数变量 458 个，值为 0 而保持线性的流量或存量 24 个。归一残差 1.211191e+01。雅可比 184.7 s，残差评价 685 次。
最小奇异值 3.001512e-34，最大奇异值 3.407447e+05，条件数 1.135243e+39。被放大 58 列，到上限仍不可分辨 24 列。

| 列 | x | 起始步长 | 最终步长 | 上限 |
|---|---|---|---|---|
| freeboard[b=0,geo=7].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=0,geo=7].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=1,geo=6].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=1,geo=6].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=2,geo=5].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=2,geo=5].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=3,geo=4].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=3,geo=4].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=4,geo=3].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=4,geo=3].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=5,geo=2].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=5,geo=2].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=6,geo=1].gas_outflow.gas.C6H6 | 4.002555e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=6,geo=1].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=6,geo=1].ash.hover.k3 | 4.201339e-07 | 1.490116e-08 | 1.490116e-06 | 1.000000e+00 |
| freeboard[b=7,geo=0].gas_outflow.gas.C6H6 | 2.754871e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| freeboard[b=7,geo=0].gas_outflow.gas.C10H8 | 1.695511e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| freeboard[b=7,geo=0].ash.hover.k3 | 4.201339e-07 | 1.490116e-08 | 1.490116e-06 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.suspension.O2 | 1.642539e-23 | 1.490116e-08 | 1.490116e-04 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.suspension.C6H6 | 1.644547e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.bubble.C6H6 | 1.110324e-13 | 1.490116e-08 | 1.490116e-02 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.490116e-08 | 1.000000e+00 | 1.000000e+00 |

| 不可分辨列 | x | 最终步长 |
|---|---|---|
| freeboard[b=0,geo=7].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=1,geo=6].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=2,geo=5].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=3,geo=4].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=4,geo=3].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=5,geo=2].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=6,geo=1].gas_outflow.gas.C10H8 | 2.463410e-20 | 1.000000e+00 |
| freeboard[b=7,geo=0].gas_outflow.gas.C10H8 | 1.695511e-20 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=0,geo=7].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=1,geo=6].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=2,geo=5].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=3,geo=4].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=4,geo=3].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=5,geo=2].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=6,geo=1].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.suspension.C10H8 | 1.012152e-20 | 1.000000e+00 |
| bed[b=7,geo=0].gas_outflow.bubble.C10H8 | 6.833595e-21 | 1.000000e+00 |

## 外层第 1 轮

诊断记录在 `reports/t26e-validation/round1/record.json`。求解状态 NOT_CONVERGED。原因：雅可比不能求解。内层状态 NOT_CONVERGED。外层 1 轮，内层迭代 1 次，墙钟 259.595919125 s。
未通过。停止，不重跑，不另试其它手段。
