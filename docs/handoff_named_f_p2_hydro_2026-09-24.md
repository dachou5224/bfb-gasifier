# 交接：45 维 E8 水力延拓停在尺度 0.96（2026-09-24）

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读本文件，然后读 `docs/audit_reports/p2_tar_start_2026-09-23.md` 和 `docs/handoff_named_f_p2_2026-09-22.md`。

**分支：** `feat/nr-hamel-wirsum-core-startwert`  
**日期：** 2026-09-24  
**计划位置：** 具名核 P2。P0 雅可比、P1 物料账、P2 映射已经闭合。当前是 P2 的求根：把 45 维 \(H(x)=F(x;P(x;b))=0\) 沿水力尺度推到 1。  
**状态：** 命运为 0 时，水力尺度 **0.9606** 过了分块门。尺度 1 没有到。Eq. 3.72 没开。正式快照没改。P3 没开。  
**Git：** 用户未要求则不 commit。

本核不经 `Reactor`。不要混 reduced 3-cell 与 `thesis_mode`。

---

## 1. 卡点

已经在解 **45 维**方程。E6、E7 仍是 43 维。E8 多出来的两维是两个焦油代表物的气泡相摩尔流量，下标 42、43；温度在 44。

命运保持 0，从已复核的尺度 0.5 延拓，最后通过点是：

| 量 | 值 |
|---|---|
| 文件 | `docs/audit_reports/p2_tar_start_2026-09-23/bfb_p2_hydro_march.npz` |
| 水力尺度 | 0.9605625 |
| 命运尺度 | 0 |
| \(T\) | 1384.7996 K |
| 气体 L2 | \(8.18\times 10^{-7}\) |
| 固体 L2 | \(1.88\times 10^{-8}\) |
| 焦油 L2 | \(4.92\times 10^{-7}\) |
| 能量 | \(-6.08\times 10^{-5}\) W |

尺度 0.5 的文件 `bfb_p2_hydro_accepted.npz` 没有被覆盖，2026-09-24 复核仍过门（气体 \(3.18\times 10^{-7}\)，焦油 \(8.01\times 10^{-10}\)，能量 \(5.38\times 10^{-5}\) W）。

**过不去的是自由空间氧这一行。** 物种顺序是 CO、CO2、H2、H2O、CH4、O2、N2。自由空间从下标 16 起，氧是 **21**。两级旋风氧是 28、35。床层密相氧是 5，气泡氧是 12。

自由空间气体行先写成输运差，再减掉自由空间均相源：

```459:459:src/solvers/wirsum_1997_coupled_f_two_bin_pyro.py
F_fb = u["N_fb"] - u["N_d"] - u["N_b"]
```

```607:608:src/solvers/wirsum_1997_named_f.py
    if int(nq.stage) >= 7:
        F[16:23] = F[16:23] - fb_gas_source(sp, qT, rate_scale=rate_scale, enabled=on)
```

`fb_gas_source` 只用 R5、R6、R8、R12，气泡动力学，\(\rho_{\mathrm{ash}}=0\)。氧精确为 0 时这些速率关掉，第 21 行就等于床层氧之和的相反数。尺度 0.94 附近这一行约 \(-9.2\times 10^{-7}\) mol/s，还在门内。尺度 0.96 附近气泡氧大约 \(4\times 10^{-5}\) mol/s，这一行单独就超过 \(10^{-6}\)。

任何一个可表示的正氧流量都会把自由空间反应点着。在尺度 0.95、当时的合格状态上扫过 \(N_{\mathrm{fb,O_2}}\)：

| \(N_{\mathrm{fb,O_2}}\) / (mol/s) | 第 21 行 / (mol/s) |
|---|---|
| 0 | \(-9.17\times 10^{-7}\) |
| \(10^{-40}\) | \(-8.95\times 10^{-7}\) |
| \(10^{-35}\) | \(-1.94\times 10^{-7}\) |
| \(10^{-34}\) | \(+5.26\times 10^{-7}\) |
| \(10^{-30}\) | \(+2.20\times 10^{-5}\) |
| \(10^{-18}\) | \(+9.12\times 10^{-2}\) |
| \(10^{-6}\) | \(+2.26\times 10^{8}\) |

过零点大约在 \(2\times 10^{-35}\) mol/s。把氧配到这个根上，第 21 行可以到 \(10^{-16}\) 量级，但自由空间 CO、CO2（下标 16、17）各移动约 \(1.8\times 10^{-6}\) mol/s。在尺度 0.94 的合格点上，配平前气体 L2 是 \(9.2\times 10^{-7}\)，配平后变成 \(2.6\times 10^{-6}\)，门反而丢了。动的是 R5（CO 氧化），不是 R12：H2、H2O 那两行几乎不动。

尺度约 0.962、冻结 \(q\) 上做完一次对数最小二乘之后，床层氧之和约 \(4.0\times 10^{-5}\) mol/s，自由空间氧被标量配平到约 \(10^{-30}\) mol/s。冻结残差最大的行是：

| 下标 | 含义 | \(F\) / (mol/s) | \(x\) |
|---|---|---|---|
| 17 | 自由空间 CO2 | \(-1.6\times 10^{-5}\) | 0.38 |
| 16 | 自由空间 CO | \(+1.6\times 10^{-5}\) | 0.052 |
| 43 | 气泡焦油 2 | \(+1.2\times 10^{-5}\) | \(4\times 10^{-14}\) |
| 41 | 悬浮焦油 2 | \(+1.0\times 10^{-5}\) | \(3.6\times 10^{-4}\) |
| 42 | 气泡焦油 1 | \(+7.7\times 10^{-6}\) | \(4\times 10^{-14}\) |
| 40 | 悬浮焦油 1 | \(+6.7\times 10^{-6}\) | \(2.3\times 10^{-4}\) |

内层 `least_squares` 在这些行还没到 \(10^{-6}\) 时就因 `ftol` 停下（一次测量：`nfev=38`，`cost≈2.4×10^{-6}`，冻结气体 L2 仍是 \(2.3\times 10^{-5}\)，焦油 L2 仍是 \(1.9\times 10^{-5}\)）。加函数评估次数没用，它已经停在缩放目标的驻点上。外层 Picard 每刷新一次 \(P(x)\)，气体和焦油大约只削掉一成。从 0.9606 再走 \(\Delta s=2\times 10^{-4}\)，8 次外层后气体 L2 仍是 \(1.05\times 10^{-6}\)，焦油 L2 是 \(8.1\times 10^{-7}\)。

密相氧（下标 5）在这个状态上大约是 \(10^{-19}\) mol/s，仍在联合对数步的活跃集里。自由空间氧和两级旋风氧（21、28、35）已经被排除。差商步长在对数坐标上是 \(10^{-4}\)，温度是 \(10^{-7}\)。

只解自由空间 CO、CO2 这两维，3 次评估就能把第 16、17 行打到 \(10^{-16}\)。把这一步嵌进每次内层之后，外层收缩率没有变好，所以没有留在代码里。只解四个焦油未知量可以把焦油行打到 \(10^{-10}\)，气体 L2 会升到约 \(0.009\)。焦油和主气体必须一起动。

---

## 2. 不要再做的事

- 不要把自由空间氧设成床层氧之和。一次这样的赋值把气体 L2 打到 \(10^{8}\)、能量打到约 \(3.5\times 10^{4}\) W。
- 不要用 \(|F_i|\) 当对数最小二乘的初值再把微量氧放进联合步。
- 不要把氧抬到 \(10^{-18}\) mol/s 再做差商。右导数是硬的。
- 不要代数消去气泡氧。床层氧仍是未知量。
- 不要翻转焦油交换符号。审计已核对它和 Hamel (2.4)/(2.5) 一致。
- 不要放宽气体/固体/焦油 L2 \(<10^{-6}\) 或 \(|E|<1\) W。不要再用含 W 的整体 \(\|F\|_2<10^{-6}\) 当通过条件。
- 不要把 \(s<1\) 或命运 0 写成完整 E8。`recertify_named_e8` 评价的是完整 \(P(x)\)，含命运映射。
- 不要改 `rate_R1`–`rate_R12`。不要把温度拟合到 1100 K。
- 不要覆盖 `bfb_p2_hydro_accepted.npz`。失败候选不得覆盖 `bfb_p2_hydro_march.npz` 里的最后通过点。
- 不要更新 `data/named_e8_endpoint_2026-09-17.json`。那里的 E8 仍是 43 维旧快照。
- 不要开 P3、旋风 \(K\)、返料腿、生产 `hamel_outer_refresh`。
- 用户未要求则不 commit。

---

## 3. 已经装上、不要推倒的部分

45 维残差里，两个焦油代表物各有悬浮相和气泡相，热解给料只进悬浮相，交换按相、等大反号。R10 两相都用该相自己的 \(C_{\mathrm{tar}}\) 和 \(C_{\mathrm{O_2}}\)。阶段 8 以下仍是原来的共享焦油槽，E6/E7 保持 43 维。

验收在 `named_endpoint_blocks`：气体、固体、焦油各自的 L2，能量单独用瓦。E8 延拓用这个门，不再用混合 L2。

内层求根是 `_e8_log_scaled_correct`：已为正的未知量走 \(x=x_{\mathrm{ref}}\exp(v)\)，残差按气体/焦油 0.01 mol/s、固体 0.001 kg/s、能量 1 W 缩放。这是数值工具，不是 Hamel 规定的牛顿法。事后逐项裁剪牛顿步会把方向做成上升方向，不要回到那条路。

`_balance_freeboard_o2` 只在自由空间氧不超过 \(10^{-12}\) mol/s 时，用对分把它配到第 21 行的过零点，并且不让 21、28、35 进入联合对数步。它能过 0.95 到 0.96，过不了上面写的 CO/CO2 和焦油地板。

`solve_named_e8_from_e7` 不再把尺度 0 当成已经接受的根。气泡焦油初值是 \(1.2\times 10^{-5}\)、\(1.8\times 10^{-5}\) mol/s，只为避开 \(\sqrt{C_{\mathrm{tar}}}\) 在 0 处没有导数，不是给料，也不是根。

---

## 4. 下一手

从 `bfb_p2_hydro_march.npz` 继续，命运保持 0。每次冻结 \(q\) 的内解之后，用 `residual_e7_plus_scaled_e8` 重算复合残差，再用 `named_endpoint_blocks` 验收。

要回答的是：在床层气泡氧已经是 \(10^{-5}\) mol/s 的时候，怎样让自由空间氧行和它带动的 CO、CO2、焦油行同时进门。密相氧约 \(10^{-19}\) mol/s 仍在联合差商里，这一列是硬的；把它留在活跃集里，全系统步会在 CO/CO2 还差 \(10^{-5}\) 时停住。只动 CO/CO2 能消那两行，只动焦油会把气体打到 \(10^{-2}\)。

水力尺度 1、命运 0 现场过门之后，才打开 Eq. 3.72。到 \((1,1)\) 再用 `recertify_named_e8`。通过之后才写正式快照。

两种焦油仍共用一个总浓度和固定混合计量，主气体与焦油的浓度分母也还没完全统一。这是单独的模型一致性项，不是这次 0.96 停住的直接原因。

---

## 5. 复现

从项目根目录：

```bash
PYTHONPATH=. .venv/bin/python -c "
import numpy as np
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q, residual_e7_plus_scaled_e8
d=np.load('docs/audit_reports/p2_tar_start_2026-09-23/bfb_p2_hydro_march.npz')
h=residual_e7_plus_scaled_e8(d['x'], freeze_e7_q(), freeze_e8_q(), hydro_scale=float(d['hydro_scale']), fate_scale=0.0)
print(float(d['hydro_scale']), named_endpoint_blocks(h))
"
```

尺度 0 和 0.5 的独立复核仍是：

```bash
PYTHONPATH=. .venv/bin/python docs/audit_reports/p2_tar_start_2026-09-23/verify_saved_states.py
```

快测（不含完整求解和旧 43 维快照）：

```bash
PYTHONPATH=. .venv/bin/python -m pytest tests/test_wirsum_1997_named_f_e8.py -q \
  --deselect tests/test_wirsum_1997_named_f_e8.py::test_e8_saved_endpoint_recertifies_live_not_cached_boolean \
  --deselect tests/test_wirsum_1997_named_f_e8.py::test_e8_solve_does_not_change_rate_or_claim_furnace
```

2026-09-24 这条命令 5 通过、2 跳过。
