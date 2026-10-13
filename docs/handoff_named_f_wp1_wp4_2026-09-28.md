# 交接：诊断工作包 1–4（2026-09-28）

审阅工作包 A–D 的结果在 `docs/handoff_named_f_review_abcd_2026-09-28.md`。本文件只保留诊断工作包 1–4。

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读本文件。诊断原文是 `docs/audit_reports/plan_ab_diagnosis_2026-09-28.md`。更早的 A/B 执行记录在 `docs/handoff_named_f_plan_ab_2026-09-28.md`，其中「C 没开始」和旧冷启动 219.8 K 已被本轮覆盖。

**分支：** `feat/nr-hamel-wirsum-core-startwert`  
**日期：** 2026-09-28  
**Git：** 用户未要求则不 commit。不要覆盖 `data/named_e8_endpoint_2026-09-17.json`。

本核仍不经 `Reactor`。不要改 `rate_R1`–`rate_R12`，不要放宽分块门，不要把超时、最小二乘停步或信赖域边界上的温度当成根。

分块门不变：气体、固体、焦油 L2 \(< 10^{-6}\)（mol/s、kg/s、mol/s），\(n\ge 46\) 时最后两行能量各自绝对值 \(< 1\,\mathrm{W}\)。`energy_W` 是两行之和，不能当门。文件里的 `passed` 必须能由重算的 \(F\) 复核。

---

## 1. 这一轮改了什么

诊断四件工作包都做完了。求解入口仍是具名核，没有切到生产 `Cell`。

| 工作包 | 代码 | 验收 |
|---|---|---|
| 1 统一控制器 | `src/solvers/wirsum_1997_named_f_control.py`：`reactivate_zero_unknowns`、`solve_named_controlled`、`advance_parameter_factor` | `docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp1_three_cases.json` |
| 2 已定位的成本 | `coupled_residual_named` 在 stage \(\ge 11\) 时每支残差只算一次 `fb_char_balance` | `wp2_profile.json`。基线 \(F\) 不变 |
| 3 输入合同与独立账 | `src/solvers/wesseling_case_contract.py` | `tests/test_wesseling_case_contract.py` |
| 4 受控重跑与条件性预测 | 同上，外加 `src/solvers/wesseling_two_zone.py` | `wp4_parameter_branches_*.json`、`wp4_wesseling1.json`、`tests/test_wp4_controlled_reruns.py` |

控制器的接受规则：刷新 \(P(x)\) 之后，分块门通过，或复合残差下降，才接受。温度须在 \([250, 2500]\,\mathrm{K}\)，非温度未知量非负。连续两次复合残差拒绝就停。内层 `ftol=1e-8`、`max_nfev=800`、`temperature_factor=1.5`。`F_calls` 含差分雅可比里的 \(F\)，不要和 scipy 的 `nfev` 比较。

零变量重激活：整组气泡或自由空间焦油全为零、且对应进料为正时，先试整组；其余零用该行 \(|F_i|\) 下降来选试探值。全局复合残差在 \(10^{16}\) 量级时看不见单行变化，所以不用它判断单个零。

工作包 2 没有做解析雅可比，也没有雅可比缓存。10 秒剖面：4578 次 \(F\)，458 次/秒。差分雅可比仍是主成本。

---

## 2. 重跑结果

P3 端点仍是 `docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz`。床温 1363.441 K，自由空间 1384.093 K。下面的温度都相对这个端点。

### 2.1 同一参数的初值鲁棒性

记录：`wp1_three_cases.json`。每例 120 s。控制器已接入，本轮工作包 4 没有再跑。

| 案例 | 结果 | 床温 K | 自由空间 K | 主导残差 |
|---|---|---:|---:|---|
| 温度 +50 K | 过门，108 s，3 次外层 | 1363.441 | 1384.093 | 气体 L2 \(5.37\times10^{-7}\) mol/s |
| 正流量/库存 ×0.5 | 超时，1 次外层 | 1342.36 | 1420.26 | 气体 L2 \(2.58\times10^{-3}\) mol/s（门的 2579 倍）。最坏行下标 48 是 5.49 W，气体先出门 |
| 冷初值 | 超时，4 次外层 | 1000.00 | 972.77 | 气体 L2 58.5 mol/s，能量和 \(-2.82\times10^{4}\) W，最坏行下标 49 为 \(-4.82\times10^{4}\) W |

冷初值打开的下标：44、45、7、8、10。自由空间不再掉到 220 K。1000 K 是求解域里的停点，不是过门的解。旧的无条件循环记录 `p3_cold_start.json`（219.8 K）不要再当当前冷启动。

### 2.2 改变参数的物理分支

空气改 `N_in_b`，燃料改 `FuelAtomFeed` 后用 `dry_mass_from_cho_atoms` 重算干质量。起点是 P3 端点，必须走到目标因子并由分块门认证。

| 分支 | 认证到的因子 | 床温 K | 自由空间 K | 记录 |
|---|---:|---:|---:|---|
| 燃料 ×1.05 | 1.05，过门 | 1356.452 | 1354.771 | `wp4_parameter_branches_fuel_rest.json` |
| 空气 ×0.95 | 停在 0.975，过门 | 1354.434 | 1360.628 | `wp4_parameter_branches_air_rest.json` |

燃料 +5% 的路径是先过 +2.5%（138 s），再从该状态过 +5%（154 s）。气体门比值 0.034。

空气 −5%：从端点直接跳到 0.95，120 s 内层不返回。0.975 可以过门。再从 0.975 的解降到 0.95，160 s 内层仍不返回，温度保持 1354.434 / 1360.628 K。该状态下气体 L2 为 0.0388 mol/s（门的 \(3.88\times10^{4}\) 倍），能量第 48 行 −2.61 W。这不是 −5% 的解。

空气 +5% 和燃料 −5% 只有旧的无条件循环记录（`p3_input_pm5_benchmark.json`），没有在控制器上重跑。

### 2.3 E7 链和冷启动

`solve_named_e8_from_e7` 没有改成调用控制器。先前 618 s 停在水力 0.925、命运 0，`endpoint_certified=false`。本轮不重跑。冷启动的控制器结果就是上表，120 s 未过门。

### 2.4 Wesseling 1 两区条件性预测

`wp4_wesseling1.json`。`passed=false`，`prediction=null`。原因 `merit_reject`，151 s。床温 1650 K、自由空间 733.33 K 是 1100 K 乘除 `temperature_factor=1.5` 的信赖域边界，不是预测，也不和论文模拟 1100 K 或测量 1120 K 比较。气体门比值 \(2.02\times10^{8}\)，能量和约 \(-3.50\times10^{5}\) W。

初值上两档炭 `w_undecided=[1,1]`，`tau_s=[0,0]`。自由空间高 8.5 m，轨迹窗口仍是 5 s，未判定不是逃逸。

---

## 3. 堵点

1. **空气 −5% 的第二步内层校正不返回。** 已认证邻居是 −2.5%。再降 2.5% 时，160 s 内接受状态不动，气体残差 0.039 mol/s。下一步若再跑，只给这一步完整时限，不要把整段 5% 从头重跑。没有新的步长或雅可比改动时，不要再加长到 10 分钟。

2. **半尺度和冷初值仍不过门。** 半尺度 120 s 只有 1 次外层，气体、焦油、固体都在门外，能量第 48 行 5.49 W。冷初值能量在 \(10^{4}\) W，气体 58.5 mol/s。零变量 7、8、10、44、45 已打开；诊断里曾列出的气泡下标 9 和 11 不在这次 `opened_zeros` 里。

3. **差分雅可比仍占内层时间。** 一次 50 列雅可比约 100 次 \(F\)。工作包 2 只去掉了残差内部重复的炭衡算。活跃集或轨迹命运变化时必须重算雅可比。

4. **E7 完整链和控制器是两条路径。** 把控制器接进水力延拓才值得重跑 E7。现在重跑会重复水力 0.925、命运 0 的超时。

5. **Wesseling 1 不能当成已验证。** 缺蒸汽、一次风份额、固定床料质量和粒径、燃料热值、循环质量流量。本轮的显式假设是：蒸汽 0，全部表载空气进气泡入口，二次风不放在 \(\xi\approx 0.6\)，热值用具名核示意值 \(-3.96\times10^{6}\) J/kg，壁面热损 0，炭粒径保持 \(2\times10^{-4}\) m 和 \(3\times10^{-4}\) m，床体积按截面积和床高相对具名核默认几何放大，\(V_b:V_d\) 保持 0.10:0.32。两区不能验收二次风局部温度峰。图 7.4 的摩尔分数没有总摩尔流量，单独不能做元素平衡。

6. **Wesseling 2 仍缺蒸汽和一次氧份额。** 不要为它单独调速率。

---

## 4. 不要先做的事

- 不要把 1650 K / 733 K、空气 −5% 未移动的温度，或冷启动的 1000 K / 973 K 写成预测或新端点。
- 不要把 Wesseling 的 14.41 mol/s 氧气塞进具名核默认的 0.40 mol/s 气泡入口，再把它叫做 Wesseling 预测。
- 不要用出口组成反推蒸汽、一次风或床料。不要把数字化 5% 设成模型容差。论文模拟温度和测量温度继续分开。
- 不要加长 5 s 轨迹窗口来消除 `w_undecided=1`，不要用 \(H/u_g\) 当停留时间。
- 不要把燃料热值调到靠近 1100 K。`GEO_FB_RHO_G_T_REF_K` 只是冻结密度参考。
- 不要覆盖 E8 快照，不要重跑 `data/named_e8_recompute_2026-09-17.py` 来「刷新」它。
- 不要改速率式，不要放宽 \(10^{-6}\) 或 1 W，不要打开式 (3.74) 当默认气场。
- 用户未要求则不 commit。

---

## 5. 复现

```bash
PYTHONPATH=. .venv/bin/python -m pytest tests/test_wesseling_case_contract.py tests/test_wp4_controlled_reruns.py -q
```

条件边界和两区 q：`src/solvers/wesseling_case_contract.py`、`src/solvers/wesseling_two_zone.py`。  
参数延拓脚本：`docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp4_parameter_branches.py`。  
Wesseling 求解脚本：同目录 `wp4_wesseling1.py`（180 s，不读 P3 端点）。
