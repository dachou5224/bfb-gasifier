---
status: accepted
---

# 床区两个标量子问题的搜索区间、轮数与原式残差

> **状态：已批准（2026-10-23 补登；随 PR #78 复审通过）。** 本 ADR 登记 PR #78 审查退回的 P2。Hamel 没有给出这些数。它们是 CHARTER §5 允许的求解路径参数，不改原方程。判定仍是求出的点使未修改的原式满足下面的残差。未经批准前，代码不得把它们写成默认值。

## 原文情况

式 (3.13)、(3.19) 把气泡上升速度写成关于 u_b 的隐式代数式。式 (3.10) 的几何存量随 ρ_P 变化，ADR 0014 A4 要求由 M_bed 反求 ρ_P。原文没有搜索区间、轮数或残差容差。ADR 0009、0014 只批准质量闭合的物理关系。ADR 0017 第 6–10 项是外层与牛顿法，不覆盖这两个标量子问题。

## 裁定

调用方把下列参数作为 `HydroSearch` 传入 `bubble_rise_velocity` 与 `particle_density_from_bed_mass`。函数内不设默认值，也不写入 `PhysicalConfig`。T25 装配求解器时再把同一组名字放入 `SolverConfig`。

| 参数 | 登记值 | 作用 |
|---|---|---|
| `velocity_floor` | 1·10⁻¹² m/s | 气泡速度搜索下界 |
| `velocity_margin` | 1 m/s | 上界 = \|excess\| + υ_b·u_b,i 的孤立尺度 + 此余量 |
| `rise_bisection_steps` | 80 | 变号区间上的二分轮数，只负责缩小区间 |
| `rise_residual_tolerance` | 1·10⁻⁸ | 接受条件：\|u_b − excess − υ_b u_b,i\| / max(\|u_b\|, `velocity_floor`) ≤ 此值 |
| `density_start_factor` | 2 | 扫描起点 ρ_P = `density_start_factor`·ρ_g |
| `density_growth_factor` | 1.35 | 下一点 = 当前密度 × 此因子，必须大于 1 |
| `density_scan_limit` | 36 | 扫描次数上限 |
| `density_cap` | 8000 kg/m³ | 超过此密度不再计算。这不是颗粒密度的物理上限 |
| `density_bisection_steps` | 60 | 质量跨过 M_bed 的括号上的二分轮数 |
| `density_mass_tolerance` | 1·10⁻⁶ | 接受条件：\|Σ M_i − M_bed\| / M_bed ≤ 此值 |

不能鼓泡的密度点跳过，继续扫描。跳过不是解。

## 失败时怎么说

轮数用尽仍不满足残差，或登记区间内没有找到变号、没有找到跨过 M_bed 的括号，都报告**数值搜索失败**。不得据此断言原式没有根，也不得断言质量闭合没有物理解。`density_cap` 外面没有被搜索过，更不能当作物理上无解。

返回前必须用原式再算一次残差。只做完登记轮数、不核验残差，不得当作解返回。

## 审查记录

- **2026-10-23：总司令批准（补登）**（CHARTER v0.10 §8，目标所有者可否决）。本 ADR 随 T19 由审查方复审通过（PR #78）并已合并，批准记录当时漏登，由 M2 审计后续（ADR 0040 第 3 项、PR #107）补登。数值不改；实现约束补充：允许以 `precalc.py` 的模块常量实现，注明 `求解器参数：ADR 0020`。
