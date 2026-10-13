# WSV400 A 探针

归档上的生产 `Reactor`，论文模式，只补 `probe/` 配置。三组都在 420 s 硬超时前返回。三组都不收敛。默认初值下牛顿步没有离开 Vorabrechnung。床顶气相碳原子是进料的三倍，氢不到进料的三分之一，干基 H₂ 为 0。

数字来自 `probe/out/thesis_default.json`、`probe/out/freeboard_n4.json`、`probe/out/seed_budget_progress.json` 和 `probe/out/run.log`。

## 复现

```text
cd /Users/liuzhen/AI-projects/bfb-T1
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python probe/run_wsv400a.py --setting thesis_default
/Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python probe/run_wsv400a.py --setting freeboard_n4
/Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python probe/run_wsv400a.py --setting seed_budget_progress
```

`Reactor.solve` 的参数三组相同：`max_global_iter=20`，`tol_global=1e-4`，`solver=global_nr`。进程内 `SIGALRM` 420 s。实际耗时：

| 设置 | 开始 (UTC) | 结束 (UTC) | 耗时 | 进程退出 |
|---|---|---|---|---|
| thesis_default | 2026-10-13T05:54:48Z | 2026-10-13T05:59:04Z | 255.3 s | 0，`status=SOLVED` |
| freeboard_n4 | 2026-10-13T05:59:04Z | 2026-10-13T06:04:41Z | 335.7 s | 0，`status=SOLVED` |
| seed_budget_progress | 2026-10-13T06:04:41Z | 2026-10-13T06:09:50Z | 307.7 s | 0，`status=SOLVED` |

方程数：默认与 `seed_budget_progress` 为 12 个求解 cell、452 条方程；`freeboard_n4` 为 16 个 cell、496 条方程。12 个 cell 里有 10 个床层格，外加旋风和返料腿。

## 1. 输入

工况 `CASE_WSV400_A`，共享块 `_WSV400_SHARED`。写入字段见各 JSON 的 `meta.mapping`。

### 直接写入

| 配置字段 | 值 | 来源 |
|---|---|---|
| `fuel_feed` | 35.2/3600 = 0.009777… kg/s | `inputs.fuel.feed_rate_kg_h` |
| `moisture_wt` | 27.3 | 工业分析，湿基 wt% |
| `ash_dry_wt` | 6.34 | 工业分析，干基 wt% |
| `VM_daf` | 93.62 | 工业分析 |
| `C_dry` / `H_dry` / `O_dry` | 46.27 / 5.94 / 41.01 | 干基元素，与共享块 A–E 相同 |
| `S_dry` | 0.01 | 干基 wt%。只进 ER 化学计量 |
| `nitrogen_fraction` | 0.4 | 干基氮 wt%。账本按百分数除以 100 |
| `D_bed` | 0.4 m | 工况与共享块一致 |
| `P` | 1.0×10⁵ Pa | `pressure_MPa=0.1` |
| `T_inlet` | 293 K | `T_inlet_K` |
| `ER` | 0.4532 | 操作条件 |
| `primary_agent` | `air` | “Air (only, no steam)” |
| `steam_to_o2_molar` | 0 | `steam_feed_kg_h=0`。类默认 0.8 会加蒸汽 |
| `recirculation_frac` | 0 | `recirculation=false`。类默认 0.1 |
| `recycle_gas` | false | 同上。类默认 true |
| `fuel_type` | `biomass` | 见下方推断 |
| `thesis_mode` | true | 任务卡要求的论文模式。其余 `*_thesis` 开关保持类默认 |

`__post_init__` 由 ER 算出的进料（`meta.computed_feeds_mol_s`）：

- O₂ = 0.1302853409639094 mol/s
- H₂O = 0
- N₂ = 0.49012104457851635 mol/s

折成干空气 50.06083941197255 Nm³/h（22.414 Nm³/kmol，O₂ 体积分数 0.21）。表上一次风和总空气都是 50.0 Nm³/h，`calculated_flows.air_actual_Nm3_h` 是 50.3。三份数相差约 0.3 Nm³/h，探针没有另写一套摩尔流去贴 50.0。

### 推断

| 字段 | 值 | 说明 |
|---|---|---|
| `H_bed` | 1.2 m | 工况 JSON 自带说明：图 7.16 的 ξ≈0.4，约 40%×3 m。共享块 `bed_height_m` 仍是 `__MISSING__` |
| `H_freeboard` | 1.8 m | 总高 3.0 m 减床高 1.2 m。总高在工况和共享块里都有 |
| `fuel_type=biomass` | — | 燃料类型是 Natural Wood。`TarFuelType` 只有 `coal` 和 `biomass` |

### 沿用类默认，因为 JSON 没有可写的数

| 字段 | 默认 | 原因 |
|---|---|---|
| `n_cells` | 10 | 没有轴向格子数 |
| `n_freeboard_cells` | 0 | 默认不离散自由空间。只在第二组改为 4 |
| `d_p` | 0.00225 m | 燃料给的是 2–26 mm 范围和 14–17 mm 主导档（质量分数 0.22），没有单一粒径。床料 `__MISSING__`。默认 `coarsest_class` 会把质量放到最粗档，所以没有把范围写进 `d_p_min`/`d_p_max` |
| `rho_s` | 1400 kg/m³ | 床料缺失 |
| `phi_s` | 0.86 | 无形状系数 |
| `eps_mf` | 0.45 | 无最小流化空隙率 |
| `HHV_dry` | 22 MJ/kg | 输入里没有热值。`thermal_input_kW=100` 没有换算成热值 |
| `heat_loss_frac` | 0.1 | 无热损失分数 |

### 没有字段可写

- 自由空间直径 0.6 m。`ReactorConfig` 只有 `D_bed`，自由空间沿用床径。
- 燃料入口在喷嘴上方 200 mm。没有轴向位置字段，没有改 `solid_lower_inlet_frac`。
- `S_dry=0.01` 没有进入 `SolidProps.sulfur_fraction`。构造后该字段仍是 0（`balance.solid_sulfur_fraction_on_cell`）。`Reactor.__init__` 组装 `SolidProps` 时只传了 `nitrogen_fraction`。

## 2. 收敛与平衡

三组 `converged`、`converged_outer`、`converged_inner_nr`、`converged_fully` 都是 false。

分块是缩放残差的 RMS 和最大绝对值，尺度与求解器相同：气相参考量为 O₂+H₂O+N₂，固相参考量为燃料质量流率，能量参考量为 `|fuel_feed|×20 MJ/kg`。

| 量 | thesis_default | freeboard_n4 | seed_budget_progress |
|---|---|---|---|
| 气相 RMS | 0.03944022663801242 | 0.03884429228554592 | 0.03044939253984756 |
| 固相 RMS | 0.0840184002907177 | 0.05962558626968479 | 0.0015846301858126713 |
| 能量 RMS | 0.01862825537087718 | 0.016005226725581183 | 0.016825174246219075 |
| 气相 max\|·\| | 0.418630843052881 | 0.41826116093951277 | 0.18849519712681173 |
| 固相 max\|·\| | 0.4638095253763663 | 0.3245507258489592 | 0.018910424288107935 |
| 能量 max\|·\| | 0.05491976063430633 | 0.05765472632664689 | 0.04323100512584842 |
| 外层次数 | 20 | 20 | 20 |
| 内层预算 | 40 / 100 | 45 / 100 | 44 / 100 |
| `norm_history` 条数 | 20 | 23 | 24 |

元素账是进口原子减去床顶 `N_d+N_b` 的原子和床顶向上固体流。氮与 `2×N2_feed` 对齐，说明这两相摩尔流率可以当作离开床顶的气相流，不是把同一股气加了两次。旋风气相原子与床顶相同，返料腿为 0。碳多出来发生在床顶，不是旋风另造的。

进口原子（mol/s，`balance.inlet_mol_atom_s`）：C 0.2738387515148151，H 0.715256443614865，O 0.5909509026926438，N 0.9822721017779769，S 2.2168858395273487×10⁻⁵。

| 缺口（进口−床顶出口，mol/s） | thesis_default | freeboard_n4 | seed_budget_progress |
|---|---|---|---|
| C | −0.5602918220004856 | −0.5603951331170052 | −0.9327588093140462 |
| H | +0.5213074888179297 | +0.521283184576825 | +0.5214213749258831 |
| O | −0.46049098162270785 | −0.45921209889418924 | −0.9127504911255628 |
| N | +0.002030012620944155 | +0.002029906522422631 | +0.0020557574212131113 |
| S | +2.2168858395273487×10⁻⁵ | +2.2168858395273487×10⁻⁵ | +2.2168858395273487×10⁻⁵ |

负号表示床顶气相里的该元素多于进料。默认组床顶气相 C = 0.833952794302232 mol/s，固体向上只带走 C 0.0001777792130686537 mol/s。H 在气相里是 0.1939489547969353 mol/s。

能量残差是各求解 cell 的 `calc_energy_balance()`，单位 W。

| | thesis_default | freeboard_n4 | seed_budget_progress |
|---|---|---|---|
| 代数和 | 59781.88524789665 | 56147.25885033497 | 124730.22919444181 |
| 最大绝对值 | 54919.76063430633 | 57654.726326646894 | 43231.005125848416 |
| 最大所在格 | bed 4，1106.58 K | bed 4，956.58 K | bed 4，822.92 K |

求解器自报碳转化率 0.999497（默认）、1.0（自由空间）、0.997471（预算种子）。工况字段 `carbon_conversion_pct` 是 `__MISSING__`。这几个转化率与上面的碳缺口矛盾：气相碳已经多于进料，不能把它们当成外部碳平衡。

## 3. 最大的 10 条残差，以及迭代趋势

生产求解器不保存逐方程历史。能用的趋势是 `solver.norm_history`、`solver.outer_trend` 里的 `dT_outer` 和 `inner_n_iter`。范数一旦逐位不变，这 10 条方程也不再变。

### thesis_default

20 个外层的 `inner_rms_scaled_final` 都是 0.07600857210458056，`dT_outer` 都是 0，每次内层迭代数都是 1。残差向量停在初值上。

| 名次 | cell | 方程 | 缩放残差 | 原始残差 |
|---|---|---|---|---|
| 1 | bed 1 | solid_class9_ash | −0.4638095253763663 | −0.0034892643502935764 |
| 2 | bed 1 | solid_class8_ash | +0.4638095233863303 | +0.003489264335322427 |
| 3 | bed 2 | solid_class8_ash | +0.42934041128338896 | +0.003492312857888647 |
| 4 | bed 2 | solid_class9_ash | −0.42933496054077347 | −0.0034922685208124697 |
| 5 | bed 4 | gas_dense_CO | −0.418630843052881 | −0.418630843052881 |
| 6 | bed 3 | solid_class9_ash | −0.344495462305092 | −0.0029958989977503453 |
| 7 | bed 3 | solid_class8_ash | +0.3444922637121993 | +0.002995871181241015 |
| 8 | bed 2 | solid_class9_char | −0.31102223395355205 | −0.002529896832861363 |
| 9 | bed 2 | solid_class8_char | +0.31087012073581277 | +0.0025286595234159798 |
| 10 | bed 3 | solid_class9_char | −0.307428906482108 | −0.0026735503180403382 |

灰的第 8 档和第 9 档成对、反号、绝对值相同到 1e-9 kg/s 以下。床温是 1180.0 K 到 1014.8 K 的等差列，没有被牛顿改过。

### freeboard_n4

`norm_history` 从 0.07260112458605766 降到第 5 项的 0.061912761750292665，其后各项保持这个数。外层温差只有第一次是 294.2736564988886 K，其后为 0。第一次内层走了 4 步，之后每次 1 步。

最大残差换成 bed 4 密相 CO（缩放 −0.418261）。第 2 到第 5 名仍是 bed 5 上 class 8/9 的灰和炭，成对反号。bed 5、bed 6 的密相 CO 排在第 6。前 10 名见该 JSON 的 `residuals.top10`。

床温在 bed 4 为 956.5777777777778 K，bed 5 跳到 1164.9532005084013 K。自由空间四格温度 1012.3993745095815、947.4536679774676、882.6630674876196、816.8146309857781 K，平均 914.8326852401117 K。

### seed_budget_progress

`norm_history` 从 0.03612589345250573 降到 0.03554870705335245，从第 12 项起保持该值。外层温差只有第一次是 150.0 K，其后为 0。固相 RMS 落到 0.0015846301858126713，灰的粒径对不再占据前 10。前 10 变成密相 CO（bed 5 到 bed 9）、bed 1 上密相/气泡 N₂ 的一对反号，以及 bed 0 密相 O₂。bed 5 密相 CO 的缩放残差是 −0.18849519712681173。

床 0 温度 600.0 K，向上升到床顶 920.61 K。

## 4. 出口对照

对应关系采用 `data/hamel_fig7_parity_digitized.json` 的 `case_id_map`：`CASE_WSV400_A` → 图 7.17 标签 A。同一工况 JSON 的轴向出口写的是 T=980 K、CO=0.152，与图 7.17 A 的 T_FRB 模拟值 1135 K、CO 模拟值 0.168 不是同一个数。T3 还没裁定。下表的 Hamel 模拟值和实验值都取图 7.17 A，干基摩尔分数。温度对照用结果里的 `reactor_exit_T`。没有自由空间格子时，这个温度是床顶/旋风温度，不是自由空间平均温度。

图 7.17 A：T 模拟 = 实验 = 1135 K；CO 0.168 / 0.172；CO₂ 0.128 / 0.132；H₂ 0.165 / 0.175；CH₄ 0.020 / 0.022（斜线前为模拟，后为实验）。

| | thesis_default | freeboard_n4 | seed_budget_progress |
|---|---|---|---|
| `reactor_exit_T` / K | 1014.8 | 1038.3728081672841 | 920.6070933550121 |
| ΔT 相对模拟 / K | −120.2 | −96.627 | −214.393 |
| 干基 CO | 0.4728409476167274 | 0.47380684153297015 | 0.5085908945523996 |
| 干基 CO₂ | 0.14050274002802332 | 0.13953428538457804 | 0.18144618689012362 |
| 干基 H₂ | 0 | 4.791220692093686×10⁻⁷ | 1.1584292719962747×10⁻⁷ |
| 干基 CH₄ | 0.016494898251017617 | 0.016486351095639326 | 0.02104427931054223 |
| 自由空间平均温度 / K | 无自由空间格 | 914.83 | 无自由空间格 |

ΔT 都大于 20 K。CO 相对模拟值的差是 +0.3048409476167274、+0.3058068415329701、+0.34059089455239955。H₂ 的差是 −0.165、−0.1649995208779308、−0.1649998841570728。CH₄ 的绝对差小于 0.004，CO₂ 在默认组和自由空间组的差是 +0.012502740028023318 和 +0.011534285384578041。状态不收敛，H₂ 和 CO 也不在图上，所以 CH₄、CO₂ 不能单独算通过。`freeboard_n4` 的自由空间平均温度 914.8326852401117 K 离 1135 K 更远。

## 5. 障碍

1. 初值/求解路径。`src/core/reactor.py` 的 `_solve_global_nr` 在论文模式单次预计算下把外层上限设成 `max(10, max_iter)=20`，内层预算 `max_iter×5=100`。`src/solvers/outer_loop.py` 记下的默认组 20 步范数相同、`dT_outer=0`、每步内层 1 次。`src/workflow/steps/nr_inner_step.py` 的 `build_global_nr_inner_solve_fn` 没有把状态推离 Vorabrechnung。后两组只在第一次外层移动温度，然后同样停住。

2. 初值/求解路径。干基 H₂ 为 0（后两组小于 5×10⁻⁷）。默认组的状态就是初值。快氧化初值在 `src/solvers/vorabrechnung/vorab_x0.py` 的 `project_bed_holdup_fast_oxidation_closure`，由类默认 `vorab_transport_x0_fast_oxidation_closure_thesis=True` 打开。探针没有关这个开关。

3. 方程组装。床顶气相碳是进料的三倍，氢只剩不到三分之一，氧也偏多。气相残差在 `src/core/cell_balances.py` 的 `calc_gas_balance_residual`。氮账是闭合的，所以这不是出口流率取错相。

4. 方程组装。默认组前四名是相邻粒径档的灰，原始残差成对反号（bed 1 为 −0.0034892643502935764 与 +0.003489264335322427 kg/s）。固相与粒径迁移在 `calc_solid_balance_residual` 和 `calc_size_migration`（同文件）。类默认 10 个粒径档，而木片粒径分布没有写进配置，各档直径都是 2.25 mm。

5. 方程组装。`recirculation=false` 且 `n_freeboard_cells=0` 时，论文模式仍建旋风和返料腿。`Reactor._use_explicit_side_block_cells` 只看 `thesis_mode` 和类默认 `explicit_side_blocks_enabled=True`。返料腿气相原子为 0。旋风把床顶已经不平衡的气体原样送出去。

6. 初值/求解路径。能量残差最大格稳定在 bed 4，默认组 54919.76063430633 W，自由空间组 57654.726326646894 W，种子组 43231.005125848416 W。函数是 `calc_energy_balance_residual`。默认组的温度是 1180.0 K 起的等差列，牛顿没有改它。燃料入口在喷嘴上方 200 mm 没有位置字段，新鲜燃料进的是底部格。

7. 输入数据。共享块床高缺失，1.2 m 来自图 7.16 的 ξ 比例。床料、单一粒径、热值、热损失都缺失，用的是类默认（密度 1400 kg/m³、粒径 2.25 mm、热值 22 MJ/kg、热损失 0.1）。自由空间直径 0.6 m 无法写入。这些值没有按出口调整，它们仍会进入水力学和生成焓。

8. 输入数据。`S_dry` 与 `sulfur_fraction` 没有接上，硫缺口等于进料硫 2.2168858395273487×10⁻⁵ mol/s。相对碳、氢缺口可以忽略。一次风 50.0、计算值 50.3 和 ER 折算 50.06083941197255 Nm³/h 的差别也解释不了 CO 从 0.168 到 0.47。

## 6. 判断

**部分。**

全局 NR 驱动能在论文模式默认开关下，于 420 s 内返回有限的分块残差和一份出口组成。这一层可以当作归档求解器的地图来读：外层循环、残差分组、出口干基都在。

cell 衡算不能作为 `hamel/` 的迁移底稿。默认初值上牛顿一步都不走；碳、氢、氧在床顶就不闭合；三组初值或离散预算的干基 H₂ 是 0、4.791220692093686×10⁻⁷ 和 1.1584292719962747×10⁻⁷，干基 CO 是 0.4728409476167274、0.47380684153297015 和 0.5085908945523996，而图 7.17 A 的模拟值是 H₂ 0.165、CO 0.168。换种子模式或加上 4 个自由空间格，只改变停住之前的那一次温度移动，不改变这个出口。
