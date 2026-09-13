# Wirsum 1997 Vorabrechnung outer 诊断

实现：`src/solvers/wirsum_1997_outer_diagnostics.py`。
生产 `outer_loop.py` 的工程 Abgleich **未改**；本核只读场、算范数。

Wirsum (1997) §2.2 允许

```text
Vorabrechnung  →  Zellenmodell NR  →  用新状态重算水力  →  再进 NR
```

要看的是固定点

```text
h_{k+1} = H(x*(h_k))
```

是否收缩，而不是 inner residual 是否单独下降。

## 记录量

相邻 outer 的欧氏增量：

`||ε_b^{k+1}−ε_b^k||₂`，`||K_bd||`，`||K_s||`，`||u_b||`，`||T||`

以及当步 inner residual。`K_s` 取 `mean(K_solid_auf)+mean(K_solid_ab)` [1/s]，仅作输运强度诊断。

收缩比 `ρ ≈ ||Δh_{k+1}|| / ||Δh_k||`。逐场分别算；拼接向量只在温度除以 100 K 后用于单一 ρ。

`verdict.non_contractive`：尾部水力场 ρ 均值 > 1。**只下诊断，不加松弛、Anderson 或 inner clip。**

## A/B（三格仿射 map）

机器表：`data/outer_iteration_history.json`。

| map | ρ | inner residual | `non_contractive` |
|---|---|---|---|
| `h* + 0.4 (h−h*)` | 0.4 | 每步 ×0.5 | false |
| `h* + 1.3 (h−h*)` | 1.3 | 每步 ×0.5 | true |

两边 inner residual 同步下降。扩张 map 被标成 non-contractive，说明 **inner 收敛不能代替 outer 收缩**。
