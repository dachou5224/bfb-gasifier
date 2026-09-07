# 固相稳态质量平衡：流变量 vs 库存变量

对照 Hamel (1999) Eq. (2-6) 与 `holdup_transport` 实现。约化问题在 `src/core/solid_balance_reduced.py`，三格、单物种、单粒径，不经过 `Reactor`。数值见 `tests/test_solid_balance_reduced.py`。

Hamel Eq. (2-6) 的未知量是格间质量流：给料 `m_zu`、再循环 `m_rez`、上行 `m_auf`、下行 `m_ab`、系统流出 `m_aus`、反应 `r_m`、粒径迁移。各项单位均为 [kg/s]。式中没有库存 `M` [kg]。

当前生产路径把出流写成 `m_auf = K_auf M`、`m_ab = K_ab M`，从而把未知量换成 `M`。这是一条本构闭合，不是 Eq. (2-6) 本身。`K` 的单位必须是 [1/s]。

## 约化方程

开链边界：床底无下向出流，床顶出流离开系统。给料只进第 0 格，`ṁ_feed,0 = 1 kg/s`。

**Model A（流平衡）** 未知量 `m_out,i` [kg/s]，分流比 `α_i = m_up,i / m_out,i`：

$$
F_i = \dot m_{\mathrm{feed},i} + m_{\mathrm{up},i-1} + m_{\mathrm{down},i+1} - m_{\mathrm{out},i} + \dot r_i
$$

残差单位 [kg/s]，雅可比无量纲。

**Model B（库存平衡）** 未知量 `M_i` [kg]，`m_{\mathrm{up},i}=K_{\mathrm{up},i} M_i`，`m_{\mathrm{down},i}=K_{\mathrm{down},i} M_i`：

$$
F_i = \dot m_{\mathrm{feed},i} + K_{\mathrm{up},i-1} M_{i-1} + K_{\mathrm{down},i+1} M_{i+1}
- (K_{\mathrm{up},i}+K_{\mathrm{down},i}) M_i + \dot r_i - k_{\mathrm{rxn},i} M_i
$$

残差单位 [kg/s]，雅可比单位 [1/s]。`k_rxn` 仅在需要一阶反应时间尺度时出现。

基准 `K`（[1/s]）：`K_up = (0.4, 0.3, 0.5)`，`K_down = (0, 0.2, 0)`，对应 `α = (1, 0.6, 1)`。

## 秩与零空间

| 提法 | 未知量 | `K` | 反应时间尺度 | rank | 零行 / 零列 | cond | 稳态 |
|---|---|---|---|---|---|---|---|
| A 流平衡 | `m_out` [kg/s] | （由 `α` 分流） | 无 | 3 | 0 / 0 | 5.87 | 唯一：`m_out = (5/3, 5/3, 1)` kg/s |
| B 库存 | `M` [kg] | `K>0` 如上 | 无 | 3 | 0 / 0 | 5.71 | 唯一：`M* = (4.167, 3.333, 2.000)` kg；`K_sum M*` 与 A 的 `m_out` 相同 |
| B | `M` | `K=0` | 无，给料 = 0 | 0 | 3 / 3 | ∞ | 任意 `M` 都是解 |
| B | `M` | `K=0` | 无，给料 = 1 kg/s | 0 | 3 / 3 | ∞ | 无有限解；`F` 与 `M` 无关，`F=(1,0,0)` kg/s |
| B + engineering 锚 | `M` | `K=0` | 无，给料 = 1 kg/s | 3 | 0 / 0 | 1 | 人为唯一：`F = source - M` |
| B + 一阶反应 | `M` | `K=0` | `k_rxn = 0.2 s⁻¹` | 3 | 0 / 0 | 1 | 唯一：`M=(5,0,0)` kg |

`K>0` 时 A 与 B 代数等价：同一组流场。B 多出来的信息全部来自本构 `m = K M`。`K` 不能由 Eq. (2-6) 反推。

## `K → 0`

对 Model B，令 `K_up = K_down = 0` 且 `ṙ` 不依赖 `M`：

- 残差退化为 `F = ṁ_feed + ṙ`，不再含 `M`。
- 雅可比为零矩阵：三行三列全零。
- `ṁ_feed + ṙ = 0`：库存不定，任何 `M` 都满足稳态。
- `ṁ_feed + ṙ ≠ 0`：质量持续堆积或抽空，不存在有限库存稳态。

这不是数值病态，是数学问题没闭合。Newton 出现零行/零列是正确信号。

合法处理只有三类：

1. **删自由度**：该桶不再作为 NR 未知量（例如无输运的 VM/moisture 已用 `holdup_transport_char_ash_mask` 挡在 `K·m` 之外）。
2. **补本构**：给出单位为 [1/s] 的 `K` 或一阶 `k_rxn`，使 `m_out = K M` 或 `ṙ = -k_rxn M`。
3. **承认不定**：`K=0` 且净源为 0 时声明 `M` 不唯一，不把残差改写成对 `M` 的罚项。

禁止第四条：把残差改成 `source - M`。`source` 是 [kg/s]，`M` 是 [kg]，相减等价于暗设 `τ = 1 s` 的出流 `m_out = M / 1 s`。Wirsum (1997) 与 Hamel (1999) Eq. (2-6) 都没有这条本构。

## 生产代码里的锚

`src/core/cell_balances.py` 的 `_apply_holdup_transport_numerical_guardrails` 在 `K_sum ≤ 1e-14 s⁻¹` 时：

- 净源 ≈ 0：`F = -M`（disconnected）
- 净源 ≠ 0：`F = net_source - M`（source without transport）

`src/solvers/jacobian_analytic.py` 的 `_holdup_solid_self_derivative` 对这两支都返回 `-1.0`。

单桶对照（char，`K=0`，入流 1 kg/s，`M=0.25 kg`）：生产残差与约化 `k0_anchor='engineering'` 均为 `0.75`，即 `1 kg/s - 0.25 kg`。物理残差应是 `1 kg/s`，与库存无关。

`legacy_stream` 分支同样执行 `… - m_solid`，把 [kg/s] 的源项与 [kg] 的库存相减。那是另一条工程路径，不在本约化问题里展开。

P0 冻结的 Phase2 3+2 smoke 已出现 4 个零行、2 个零列。现有锚没有、也不应再靠加锚去抹平这些奇异方向。

## 建议

- thesis / `wirsum_1997_strict`：残差保持 `net_source - (K_auf+K_ab) M`，`K=0` 时允许零行。不要在 residual 内 clip、floor、anchor。配置 `holdup_k0_anchor='none'`。
- 生产 LU 路径：默认 `holdup_k0_anchor='engineering'`，保留现有 guardrail，但必须继续叫 engineering，不得当作 Eq. (2-6)。
- 开关在 `calc_solid_balance_residual(..., k0_anchor=)`，经 `ReactorConfig.holdup_k0_anchor` 拷到 `Cell.holdup_k0_anchor`。解析雅可比对角与残差同一开关。
