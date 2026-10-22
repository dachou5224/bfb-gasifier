# 玩具工况

这是玩具工况，用来核对驱动、收敛门和报告。数字不是第 1 关的解。
S1、S2、S6、S7、S8 已能通过现有接口表达。样例不启动工况进程，表中不列这几次的装置结果。

单次运行墙钟上限 90000.0 s（ADR 0037）。

## 验收

### baseline：CONVERGED

改动：不改冻结项。来源：ADR 0030。
验收状态：不通过。前提：满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | 1100 | 44 |
| X_C | 不通过 | 50 | 25.1 |
| CO | 不通过 | 0.1 | 0.0476 |
| CO2 | 通过 | 0.1 | 0.006 |
| H2 | 通过 | 0.1 | 0.006 |
| CH4 | 边界待定 | 0.01 | 0.002 |
| H2O | 边界待定 | 0.1 | 0.007 |

外部平衡：{'C': 0.0, 'H': 0.0, 'O': 0.0, 'energy': 0.0}。非负：True。残差比：{'gas': 0.0, 'solid': 0.0, 'energy': 0.0}。

### S3：NOT_APPLIED

改动：ε_P = 0.79。来源：ADR 0017 第 17 项。
说明：ε_P 写在残差模块的常数里，evaluate_residual 没有这个参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S4：NOT_APPLIED

改动：ε_P = 0.91。来源：ADR 0017 第 17 项。
说明：ε_P 写在残差模块的常数里，evaluate_residual 没有这个参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S5：NOT_APPLIED

改动：D_g 改为 O2–N2 的 Fuller 二元值。来源：ADR 0021 第 4 项。
说明：预计算只调用式 (5.20) 的相间扩散，没有 Fuller 二元扩散的参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S9：NOT_APPLIED

改动：ε_mf = 0.40，φ_S = 0.6。来源：ADR 0024。
说明：ε_mf 与 φ_S 是预计算模块的常数，不作为参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S10：NOT_APPLIED

改动：ε_mf = 0.40，φ_S = 0.9。来源：ADR 0024。
说明：ε_mf 与 φ_S 是预计算模块的常数，不作为参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S11：NOT_APPLIED

改动：ε_mf = 0.55，φ_S = 0.6。来源：ADR 0024。
说明：ε_mf 与 φ_S 是预计算模块的常数，不作为参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S12：NOT_APPLIED

改动：ε_mf = 0.55，φ_S = 0.9。来源：ADR 0024。
说明：ε_mf 与 φ_S 是预计算模块的常数，不作为参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S13：NOT_APPLIED

改动：悬停样本计为床面回落。来源：ADR 0026 第 6 项。
说明：悬停样本改记床面回落要改自由空间记账，预计算没有这个参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S14：NOT_APPLIED

改动：悬停样本计为出口带出。来源：ADR 0026 第 6 项。
说明：悬停样本改记出口带出要改自由空间记账，预计算没有这个参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S15：NOT_APPLIED

改动：Serio 第二对参数。来源：ADR 0013 第 4 项。
说明：Serio 第二对没有实现，r11_serio_rate 不接收另一组参数。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

### S16：NOT_APPLIED

改动：cell 数加倍。来源：ADR 0013 第 7 项。
说明：二次风固定在自由空间几何下标 1。格子加倍后 1.4875 m 不再落在该格，要改 T25a 的下标。
验收状态：不通过。前提：不满足。

| 量 | 状态 | 值 | 到边界的距离 |
|---|---|---|---|
| T_FRB_K | 不通过 | None | None |
| X_C | 不通过 | None | None |
| CO | 不通过 | None | None |
| CO2 | 不通过 | None | None |
| H2 | 不通过 | None | None |
| CH4 | 不通过 | None | None |
| H2O | 不通过 | None | None |

外部平衡：None。非负：None。残差比：{}。

## 诊断

每次运行的 JSON `diagnostics` 写入 T25b 公开结果里已有的量：Geldart 记录与式 (8) 距离、T_ad、T_init、y_ad、ρ_P、排灰与灰分量、抛出与回落带出系数、悬停标记与 h* 所在的逐根 δ、ū_r/v_t、产率越界与延拓、床层子步 n、k_W。
气泡相总流、d_b/D、与 KURKELA 表 5 的气速对照、自由空间元素与质量闭合、数值积分和闭式的比对、悬停焦炭消耗占比，不在这份公开结果上。本卡不改 T25b 的接口，这些量留空，不另算。
S9–S12 未接线，因此 ρ_P 可行区间没有对应的运行。

## 敏感性差值

| 运行 | 状态 | 差值 |
|---|---|---|
| S3 | NOT_APPLIED | ε_P 写在残差模块的常数里，evaluate_residual 没有这个参数 |
| S4 | NOT_APPLIED | ε_P 写在残差模块的常数里，evaluate_residual 没有这个参数 |
| S5 | NOT_APPLIED | 预计算只调用式 (5.20) 的相间扩散，没有 Fuller 二元扩散的参数 |
| S9 | NOT_APPLIED | ε_mf 与 φ_S 是预计算模块的常数，不作为参数 |
| S10 | NOT_APPLIED | ε_mf 与 φ_S 是预计算模块的常数，不作为参数 |
| S11 | NOT_APPLIED | ε_mf 与 φ_S 是预计算模块的常数，不作为参数 |
| S12 | NOT_APPLIED | ε_mf 与 φ_S 是预计算模块的常数，不作为参数 |
| S13 | NOT_APPLIED | 悬停样本改记床面回落要改自由空间记账，预计算没有这个参数 |
| S14 | NOT_APPLIED | 悬停样本改记出口带出要改自由空间记账，预计算没有这个参数 |
| S15 | NOT_APPLIED | Serio 第二对没有实现，r11_serio_rate 不接收另一组参数 |
| S16 | NOT_APPLIED | 二次风固定在自由空间几何下标 1。格子加倍后 1.4875 m 不再落在该格，要改 T25a 的下标 |

## 数值加密

- pyrolysis_radial_intervals：24
- pyrolysis_time_rtol：1e-06
- freeboard_ejection_samples：48
- freeboard_max_step_s：0.0002
- freeboard_delta_h_m：0.001
- freeboard_delta_u_m_s：0.01
- particle_tolerance_w：0.0001
- particle_max_iterations：25
- coupling_tolerance_k：0.001
- coupling_max_iterations：16
- bed_integration_intervals_per_cell：8
- bed_integration_max_intervals_per_cell：256
- note：逐根 δ 与床层子步 n 来自该次运行的预计算。玩具工况没有这些量。

## 冻结选择

| 编号 | 选择 | 来源 |
|---|---|---|
| F1 | 式 (3.43) 的 −1 为误印；α_b = u_b·ε_d/u_d | ADR 0011 第 7、10 条；ADR 0022 A1 |
| F2 | R8 用式 (5.45)，a_R8 = 0.02 | ADR 0013 第 3 项 |
| F3 | Serio 第一对 | ADR 0013 第 4 项 |
| F4 | 料斗吹扫氮计入，298.15 K | ADR 0013 第 5 项；ADR 0014 C2 |
| F5 | 8×0.15 m + 8×0.2875 m | ADR 0014 B |
| F6 | 床高 1.2 m；床料 3.963 kg；ρ_P 外层反求 | ADR 0014 A2–A4 |
| F7 | 壁面散热 21.7 kW，单一 k_W，T_amb = 298.15 K | ADR 0014 X1；ADR 0029 |
| F8 | 旋风纯分离；连接管不计反应与传热 | ADR 0014 X2、X2a |
| F9 | 第一块布风板，n_0 = 139.5，h_A = 0 | ADR 0019 |
| F10 | ε_mf = 0.46，φ_S = 0.8 | ADR 0024 |
| F11 | 粒级代表径为 Sauter 平均，最细档下限 0.0291 mm | ADR 0023 |
| F12 | ε_P = 0.85 | ADR 0017 第 17 项 |
| F13 | D_g 用 Hamel 式 (5.20)，cell 气体温度 | ADR 0021 第 4 项 |
| F14 | 灰的比热 950 J/(kg·K) | ADR 0018 第 6 项 |
| F15 | ρ_C = ρ_P | ADR 0017 第 19 项 |
| F16 | 自由空间焦油反应只用 Serio | ADR 0022 C3 |
| F17 | 悬停闭合为反应限制的滞留 | ADR 0026 |
| F18 | 热解时间用平均释放时间 | ADR 0015 第 10 项；ADR 0027 |
| F19 | 基准初值扣除 21.7 kW；只有 S1 不扣除 | ADR 0017 第 9 项；ADR 0029 |
| F20 | 内层 1e-6；外层 (a)(b)；元素 0.1 %；能量 0.5 % | ADR 0012 |
| F21 | 复合颗粒 shrinking particle | ADR 0031 |
| F22 | 燃料用同粒级焦炭的系数 | ADR 0032 |
| F23 | 平均泡间距用式 (3.57) | ADR 0028 |
| F24 | 床区 u_r = u_d/ε_d；自由状态为轨迹加权滑移；悬停为 v_t | ADR 0035 第 4 项 |
| F25 | 未收敛时非 B 类用 B 类系数延拓；收敛后按真实类别复核 | ADR 0035 第 2 项 |
