# 交接：自由空间固相轨迹与炭库存（2026-09-24）

**给审查 agent。** 先读 `docs/CLAUDE.md`，再读本文件。相关计划在 `docs/audit_reports/named_f_e8_review_2026-09-16.md` §5 P3，以及 `docs/handoff_named_f_p2_hydro_2026-09-24.md`。

**日期：** 2026-09-24  
**范围：** 只审查自由空间颗粒轨迹、命运权重，以及 E11 炭库存为什么几乎是 0。不要改速率式，不要把这个状态写成 E8 认证根，不要开 P4。用户未要求则不 commit。

---

## 1. 要回答的问题

当前 E11 解里，自由空间细粒炭库存约 \(4\times 10^{-33}\,\mathrm{kg}\)，粗粒约 \(7.3\times 10^{-5}\,\mathrm{kg}\)。细粒几乎是 0，粗粒相对床内粗粒 \(0.017\,\mathrm{kg}\) 也很小。请判断这是轨迹物理的结果，还是库存公式、温度索引或抛射闭合写错了。

第一次实现把进入量写成 `m_leave`。当时 `w_leave = 0`，库存被解成精确 0。这个公式已经改掉。改完之后粗粒库存不再是 0，但数量级仍是 \(10^{-5}\,\mathrm{kg}\)。审查对象是**改完之后的公式**和它依赖的轨迹输出。

---

## 2. 轨迹算法

实现在 `src/physics/wirsum_1997_freeboard.py` 的 `integrate_particle_fate_with_time`，注释写的是 Hamel (1999) Eq. 3.82–3.96。

对给定初速 \(u_{p,0}\) 做时间积分。步长默认 \(2\times 10^{-3}\,\mathrm{s}\)，上限 \(5\,\mathrm{s}\)。每步用 `_step_branching` 更新速度和位移，气速在这一段里冻结，跨过 \(u_p = u_g\) 换阻力分支。阻力指数 `z` 在 `hamel_z_exp` 下走 Eq. 3.86 的 \(\varphi_s\) 极限，不是 Wirsum 默认 0.44。

停止规则：

| 条件 | 返回 |
|---|---|
| \(u_{p,0}\le 0\) | `returns`，时间 0 |
| 高度达到自由空间高度 | `reaches_next_cell` |
| 已经离开床面之后高度回到 \(\le 0\) | `returns` |
| 积分到 \(t_{\max}\) 仍在空中 | `undecided` |

回落判据是高度回到床面，不是 \(u_p=0\)。

分布积分在 `src/solvers/wirsum_1997_named_f_p2.py` 的 `integrate_launch_distribution`。

- 平均抛射速度：\(u_{p,0,\mathrm{mean}} = 1.53\,u_b\)。系数 `MARTENS_U_P0_OVER_U_B = 1.53`，注释为 Hamel Eq. 3.73 / Martens (1984)。
- 初速分布：`eq372_launch_nodes`。高斯标准差 `WIRSUM_1997_VELOCITY_SIGMA = 0.60 m/s`，在 \(u\ge 0\) 上截断 \(\pm 4\sigma\)，11 个节点，梯形权重再归一。注释写这是 Hamel Eq. 3.72，不是工程均匀离散。
- 每个节点的权重 \(w_i\) 乘事件时间 \(t_i\)，加总得到 `tau_s`。
- `reaches_next_cell` 累加到 `w_leave`，`returns` 累加到 `w_return`，其余累加到 `w_undecided`。
- `occupancy_s` **只**累加未判定那一支的 \(w_i t_i\)。全部回落时它是 0，不能拿来当库存。

`finite_hydro_from_x` 对每个粒径档各积分一次，把 `w_leave`、`w_return`、`w_undecided`、`tau_s`、`occupancy_s` 放进 `NamedHydroSnapshot`。抛射系数另算：

\[
K_{\mathrm{eject}} = \frac{\dot m_{\mathrm{eject}}(M=1)}{M_{\mathrm{geom}}}
\]

`K_eject` 不包含命运权重。实际质量流在 `edges_from_x`：

\[
\dot m_{\mathrm{eject}} = K_{\mathrm{eject}} M,\quad
\dot m_{\mathrm{leave}} = w_{\mathrm{leave}}\dot m_{\mathrm{eject}},\quad
\dot m_{\mathrm{return}} = w_{\mathrm{return}}\dot m_{\mathrm{eject}}
\]

床层固体残差已经使用这三股：给料和回落进入，抛出离开。旋风第一格入口使用 `m_leave_tot`。

---

## 3. 当前炭库存公式

E11 在 E10 的 48 维后再加两个未知量：细粒、粗粒自由空间炭库存。`n_unknowns(11) = 50`。温度仍是最后两项：床温、自由空间温度。

`fb_char_balance`（`src/solvers/wirsum_1997_named_f.py`）现在是

\[
F_i = \dot m_{\mathrm{eject},i} - \frac{M_{\mathrm{fb},i}}{\tau_i} - \dot m_{\mathrm{consume},i}
\]

\(\tau_i\) 用轨迹给出的 `hydro.tau_s`。没有水力快照时才退回 \(H_{\mathrm{fb}}/\max(u_g,10^{-3})\)。R2–R4 用自由空间气体浓度、自由空间温度和 \(M_{\mathrm{fb}}\) 调用现有 `char_hetero_source`，不改速率式。消耗的气体源进自由空间气体行。回落按 `w_return` 从床层回料里扣，离开按 `w_leave` 从旋风入口扣。`w_leave=0` 时旋风入口不因自由空间反应而减少。

无消耗时的显式解是 \(M_{\mathrm{fb}}=\dot m_{\mathrm{eject}}\tau\)，不是 \(\dot m_{\mathrm{leave}}\tau\)。

---

## 4. 当前求解结果

状态文件：`docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e11.npz`。这是公式改正后、用 \(M=\dot m_{\mathrm{eject}}\tau\) 做初值、再做一次内层对数求解的结果。**整炉还没过气体门**，不要当成认证根。E8 正式快照没有改。

在 `outer_P_named` 下读到：

| 量 | 细粒 | 粗粒 |
|---|---|---|
| 床内库存 \(M\) / kg | \(8.9\times 10^{-31}\) | \(0.01738\) |
| \(K_{\mathrm{eject}}\) / s⁻¹ | \(0.01733\) | \(0.01757\) |
| \(w_{\mathrm{leave}}\) | 0 | 0 |
| \(w_{\mathrm{return}}\) | 1 | 1 |
| \(w_{\mathrm{undecided}}\) | 0 | 0 |
| \(\dot m_{\mathrm{eject}}\) / (kg/s) | \(1.54\times 10^{-32}\) | \(3.054\times 10^{-4}\) |
| \(\dot m_{\mathrm{leave}}\) / (kg/s) | 0 | 0 |
| \(\dot m_{\mathrm{return}}\) / (kg/s) | \(1.54\times 10^{-32}\) | \(3.054\times 10^{-4}\) |
| `tau_s` / s | 0.261 | 0.255 |
| `occupancy_s` / s | 0 | 0 |
| 无消耗时 \(M_{\mathrm{fb}}=\dot m_{\mathrm{eject}}\tau\) / kg | \(4.0\times 10^{-33}\) | \(7.77\times 10^{-5}\) |
| 一次内解之后的 \(M_{\mathrm{fb}}\) / kg | \(4.0\times 10^{-33}\) | \(7.35\times 10^{-5}\) |

气速 \(u_g = 0.159\,\mathrm{m/s}\)，自由空间高度 \(2.0\,\mathrm{m}\)。床温约 1365.5 K，自由空间温度约 1384.1 K。

这一次内解的复合残差：气体 L2 \(= 9.25\times 10^{-6}\,\mathrm{mol/s}\)，固体 L2 \(= 1.95\times 10^{-8}\,\mathrm{kg/s}\)，焦油 L2 \(= 1.16\times 10^{-8}\,\mathrm{mol/s}\)，两区能量之和 \(= 0.022\,\mathrm{W}\)。库存行本身约 \(-10^{-33}\) 和 \(-8\times 10^{-9}\,\mathrm{kg/s}\)。气体还高于 \(10^{-6}\) 的门。

细粒库存几乎是 0，是因为床内细粒质量已经是 \(10^{-31}\,\mathrm{kg}\)，抛出流跟着是 0。粗粒全部回落，离开流是 0，但驻留 0.255 s，库存约 \(7\times 10^{-5}\,\mathrm{kg}\)。

---

## 5. 请审查的点

1. **全部回落时库存该不该用 \(\dot m_{\mathrm{eject}}\tau\)。** `tau_s` 是所有命运的加权飞行时间，包括回落。`occupancy_s` 在没有未判定颗粒时是 0。若 Hamel 的自由空间炭只统计能到达旋风的那一支，现在的公式偏大；若统计所有在飞颗粒，用 `m_leave` 会把库存打成 0，那是错的。
2. **`finite_hydro_from_x` 内部仍写 `T = z[-1]`。** E9 起最后一维是自由空间温度，床温在倒数第二维。`outer_P_named` 自己的 DAEM 温度已经改成床温，但调用 `finite_hydro_from_x(z, b)` 时，气速、密度和随后的 \(u_b\)、轨迹仍可能按自由空间温度计算。请核对这是不是命运全部落在 `returns` 的原因之一。E8 的 45 维状态没有这个问题，因为那时最后一维就是唯一温度。
3. **细粒床库存为什么是 \(10^{-31}\,\mathrm{kg}\)。** 这来自 E8/E10 的固体解，不是轨迹积分自己生成的。若细粒档粒径或抛射系数使得细粒在更早的阶段就被耗光，自由空间细粒为 0 是上游结果。
4. **停留时间与 \(H/u_g\) 差一个量级。** \(H/u_g = 2/0.159 \approx 12.6\,\mathrm{s}\)，轨迹 `tau_s` 约 0.25 s。旧公式用前者，会把粗粒库存放大约 50 倍。请确认稳态库存应该用轨迹时间，而不是气体置换时间。
5. **消耗项的分流。** 现在按 `w_return` 减少床层回料、按 `w_leave` 减少旋风入口。`w_leave=0` 时，自由空间气化掉的碳只从回落里扣。这个分流是否和“颗粒在飞行中被消耗、因此回不到床面”一致，需要对照原文，而不是只看残差能否到 0。

---

## 6. 复现

从项目根目录：

```bash
PYTHONPATH=. .venv/bin/python - <<'PY'
import numpy as np
from dataclasses import replace
from src.solvers.wirsum_1997_named_f import outer_P_named, edges_from_x, fb_char_balance
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
x=np.load('docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e11.npz')['x']
nq=replace(freeze_e8_q(Q_W=0.0), stage=11)
q=outer_P_named(x, nq)
ed=edges_from_x(x[:39], q.energy.species)
print('w_leave', q.energy.species.w_leave)
print('w_return', q.energy.species.w_return)
print('tau', q.hydro.tau_s)
print('m_eject', ed['m_eject'])
print('m_leave', ed['m_leave'])
print('M_fb', x[46], x[47])
print('m_eject*tau', ed['m_eject']*q.hydro.tau_s)
PY
```

闭合测试：

```bash
PYTHONPATH=. .venv/bin/python -m pytest tests/test_wirsum_1997_named_f_e8.py::test_e11_freeboard_char_holdup_matches_leaving_flux -q
```

该测试只检查无 R2–R4 时 \(M_{\mathrm{fb}}=\dot m_{\mathrm{eject}}\tau\) 使库存残差为 0，以及有限水蒸气和有限库存时 R2–R4 消耗变大。它不证明轨迹命运或温度索引正确。
