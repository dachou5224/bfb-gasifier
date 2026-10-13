# 交接：R1–R4（2026-09-30）

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读本文件。任务来自 `docs/audit_reports/review_abcd_2026-09-28.md`。前一轮 A–D 在 `docs/handoff_named_f_review_abcd_2026-09-28.md`。S1–S3 的交接在 `docs/handoff_named_f_s1_s3_2026-09-30.md`。

**分支：** `feat/nr-hamel-wirsum-core-startwert`  
**日期：** 2026-09-30  
**Git：** 用户未要求则不 commit。不要覆盖 `data/named_e8_endpoint_2026-09-17.json`。

本核仍不经 `Reactor`。不要改 `rate_R1`–`rate_R12`，不要放宽分块门，不要把信赖域拒绝或最小二乘停步当成根。不要打开式 (3.74) 当默认气场，不要把抛射模型切换和求解器改动放进同一次试验。

---

## 1. 四个包的验收

| 包 | 停止条件 | 结果 |
|---|---|---|
| R1 轨迹 | 跨内部面继续积分；各轨迹末态可见。未闭合就保留末态，不加时限 | 协议测试通过。Wesseling 细粒出炉，粗粒因步长不一致停在第 1 格 |
| R2 近根 | 温度 +50 K 过门则再测半尺度；拒绝则留下一个失败点 | 18 s 拒绝。冻结残差下降，复合残差上升。半尺度没有跑 |
| R3 活跃集与接口账 | 零氢气能进入联合未知量；上下流不必相等 | 冷初值 13 s 仍拒绝。氢不再精确为 0。接口按同一股流两侧配对 |
| R4 抛射与输入 | 公式和未决项写全，不改现行抛射 | 式 (3.62)–(3.64) 按原页保留。\(P''_1\) 和接近角未定。`hamel_extended` 未改 |

---

## 2. R1

`continue_undecided_event` 返回命运、附加时间、格内高度和速度。内部上边界被接受后，积分进入下一格，不再把这次越格直接记成未判定。

恒定气场下，一整段和两格串联的回落仍然一致。测试：`tests/test_trajectory_continuation.py`。

Wesseling 1 用原有 15 s 物理续积分预算重算，没有加时。`docs/audit_reports/review_abcd_2026-09-28/r1_wesseling_traces.json`：

| 粒级 | 条数 | 命运 | `w_undecided` |
|---|---:|---|---:|
| \(2\times 10^{-4}\) m | 11 | `furnace_exit` | 0 |
| \(3\times 10^{-4}\) m | 11 | `dt_mismatch` | 1 |

粗粒的一条末态：第 1 格，高度 3.18 m，速度 0.29 m/s，发射速度 0.56 m/s，记录时间 20 s。粗步和半步没有给出同一事件，所以没有把这段时间写入 `tau_s`。

---

## 3. R2

同一个 P3 端点，两个温度各加 50 K。预算 120 s，18 s 停止。`r2_tplus50.json` 含 50 维 \(x\) 和 \(F\)。

| 量 | 值 |
|---|---|
| 原因 | `merit_reject` |
| 内层状态 | `least_squares` |
| 检查点 | `least_squares_x`，不是试探点 |
| 冻结残差 | 0.214，比起点下降 |
| 复合残差 | 0.295，刷新 \(P\) 后上升 |
| 方向标记 | `frozen_down_composite_up` |
| 雅可比预算 | 4，没有升到 8 或 16 |
| 床温 / 自由空间 | 1375.2 K / 1404.0 K |
| 气体门倍数 | 2305 |

门没有通过，半尺度没有加跑。

---

## 4. R3

气泡气体里仍为零、且该行残差不为零的未知量，以初值 0 进入同一次最小二乘的线性非负坐标。不先用一个固定小数种子去做全局残差筛选。

不读端点的冷初值，90 s 预算，13 s 停止。`r3_cold_h2.json`：

| 量 | 值 |
|---|---|
| 原因 | `merit_reject` |
| 气泡氢 | \(1.05\times 10^{-16}\) mol/s |
| 床温 / 自由空间 | 609.9 K / 455.9 K |
| 能量两行 | \(1.13\times 10^{4}\) W，\(1.32\times 10^{4}\) W |
| 气体门倍数 | \(4.3\times 10^{8}\) |
| 雅可比预算 | 8 |

氢不再精确锁在 0。停止原因不是变量被排除。这个数量级还不是展开后的氢流量。

接口账：床层出口的向上流和自由空间入口是同一股，各记一次；反应后的回落流同样两侧配对。向上流和向下流不必相等。出炉流单独记录，外部元素差额不拿这个差去当门。

---

## 5. R4

印刷页 39 的字面记录已经写入 `docs/hamel_submodels/04_freeboard_particle_trajectory_analytic.md`。

- 式 (3.62)：\(P''=(1-\cos\xi_j)\) 乘间距分布从 \(-\infty\) 到 1 的积分。
- 式 (3.63)：\(P''_k=(1-P'')\cdot P''_{k-1}\)。没有改成几何分布幂式。这一行没有 \(P''_1\)。
- 式 (3.64)：wake 通量乘 \(\sum_{k=2}^{n_B}(k-1)P''_k\)。\(n_B\) 计算上可以取 5。
- 接近角：Clift 与 Grace 30°，Briens 等 48°，没有选定。
- Pemberton 与 Davidson 写 \(\zeta_w=0.4\)。George 与 Grace 的约 40% 是从 wake 抛出的固体份额。

现行 `hamel_extended` 仍是 roof 加固定 \(0.4\times 0.25\)。式 (3.65) 的低速支是 0.1，高速支在 \(u_0/u_{mf}=10\) 时为 0.1967，并趋近 0.5。清单函数 `ejection_and_wesseling_open_items`：蒸汽缺失，二次风 30.1 Nm³/h 在 \(\xi\approx 0.6\) 但未布置，10 个燃料粒级没有写入两档炭，燃料生成焓不是高位或低位热值。

---

## 6. 卡点

1. **温度 +50 K 的复合残差在刷新 \(P\) 后变大。** 冻结子问题已经改善。下一步是外层耦合，不是把雅可比预算从 4 加到 16，也不是再缩对数步长。
2. **Wesseling 粗粒的步长细化没有闭合。** 末态已经记下。再加 15 s 飞行时间对不上这条 `dt_mismatch`。
3. **冷初值的氢只离开了精确零。** 气体和两行能量仍远在门外。不能把 \(10^{-16}\) mol/s 读成氢衡算已经打开。
4. **式 (3.63) 缺基值，接近角有两个文献值。** 在选定之前不要把求和权重 \(S\) 填进 `hamel_extended`。

---

## 7. 不要先做的事

- 不要把 1375 K、610 K、456 K 或粗粒末态写成新端点或 Wesseling 预测。
- 不要为了让粗粒 `w_undecided` 变成 0 而加大续积分预算或改用 \(H/u_g\)。
- 不要把式 (3.63) 改写成 \((1-P'')^{k-1}P''\)。
- 不要改速率式，不要放宽 \(10^{-6}\) 或 1 W，不要打开式 (3.74)。
- 不要覆盖 E8 快照。用户未要求则不 commit。

---

## 8. 复现

```bash
PYTHONPATH=. .venv/bin/python -m pytest tests/test_trajectory_continuation.py tests/test_external_state_ledger.py tests/test_hamel_eq365_ejection.py tests/test_wpd_records.py::test_r_series_keeps_the_near_root_direction_and_the_trajectory_end_states -q
```

轨迹末态：`PYTHONPATH=. .venv/bin/python docs/audit_reports/review_abcd_2026-09-28/r1_wesseling_traces.py`  
温度 +50 K：同目录 `r2_tplus50.py`  
冷初值：同目录 `r3_cold_h2.py`
