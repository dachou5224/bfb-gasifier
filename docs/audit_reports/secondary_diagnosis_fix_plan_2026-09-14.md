# 二次诊断实施计划（2026-09-14）

对照 `secondary_diagnosis_2026-09-14.md` §9。已做 **阶段 A–F5**。未宣称全炉收敛。Lapple 旋风不是主线。

不以当前 x0 能量残差升降作为物性验收。不用实验 CO。

## 阶段 A — 两条基线契约

| 路径 | 入口 | 用途 |
|---|---|---|
| `reduced_3cell_x0` | `build_initialized_lu_reactor` | 历史 P0–P6 诊断；`thesis_mode=False`，只 `generate_initial_x0` |
| `phase2_3bed_2fb_init` | `build_initialized_thesis_freeboard_reactor` | 生产初始化；`thesis_mode=True`，走 `init_precalc` |

核：`src/solvers/wirsum_1997_baseline_contract.py`  
证据：`data/wirsum_alignment_stage_a_baselines.json`

验收：每条快照标明 path / init_kind / 未知量单位 / 气固集合 / cache / 源项冻结 / 自由板归属 / 工程开关。禁止把 reduced 测试通过写成生产已修。

## 阶段 B — 物性与 R12 转录

1. **生成焓**：纯碳 char = 0；液态水分 = HHV 液态水；VM = `h_f,dry / w_VM`。撤回 2026-07-06 把液态水焓放到 char 列的结论。
2. **R12**：`r = k C_H2^{3/2} C_O2`，`k = 1.631e9 T^{3/2} exp(-3420/T)`。测试期望值按扫描式手写，不把生产参考函数当原文。
3. 固体瞬时 `c_p(T)(T-T_ref)` **本轮不改**（§2.4）。外层源项截断已具名，默认仍不打开 `hamel_outer_refresh`。

验收：标准反应焓与 Hamel 第 5 章同量级；1 kg 干料分解前后焓相等；汽化焓为正；R12 对微量氧的导数有限。

## 阶段 C — 床料 / 库存 / 输运

约化核：`src/core/solid_transport_reduced.py`（等温、无反应，不经 Reactor）。

| 模式 | 行为 |
|---|---|
| `engineering_bottom_threshold`（**生产默认**） | 上层 `M_geo`，床底活性库存 < 15% 几何库存时改用 char+ash；`K_ab` 除以总闭合库存 |
| `species_own_inventory` | 所有格 `m = K M_j`；无 15% 跳变 |
| `solid_transport_bottom_no_penetration` | 床底 `K_ab=0`，不把排灰藏进下行系数 |

`nr_solver_mode=wirsum_1997_strict` 打开物种库存 + 床底无穿透。现有 LU 默认不变。

验收：物种自分母下 `K M` 还原重建物流；工程混合分母有可测 gap；15% 阈值处库存跳变可复现。粒径迁移尚未接入本核。

## 阶段 D — 预计算上下文与内外层契约

核：`src/solvers/wirsum_1997_vorab_context.py`  
政策具名：`src/solvers/vorabrechnung/vorab_policy.py`  
证据：`data/wirsum_alignment_stage_d_vorab_context.json`

| 配置 | 默认 | 行为 |
|---|---|---|
| `vorab_outer_source_policy` | `auto` | thesis + single-shot → `engineering_single_shot`；否则 → `hamel_outer_refresh` |
| `engineering_single_shot` | 现行 thesis 默认 | init 冻源项；外层只刷水力。**不是** Hamel Bild 2.2 |
| `hamel_outer_refresh` | opt-in | 外层可重算干燥/热解。`wirsum_1997_strict` **不**自动打开 |
| `vorabrechnung_tau_for_cell` | 公式未改 | 具名为 `engineering_H_over_umf`，不是颗粒停留 |
| 自由板轨迹 | 外层不更新 | 显式图跳过 sync；非显式图 refresh 立即返回 |

验收：holdup Startwert 后未冻结 `F(x),F(x)` 可历史依赖（复现约 3.6×10²）；同一冻结上下文下 `F(x),F(x)` 与 `F(x),F(y),F(x)` 为 0。不把默认 `thesis_vorab_sources_single_shot` 改成 False。

## 阶段 E — 最小联立刚性（前三阶）

核：`src/solvers/wirsum_1997_stiffness_ladder.py`  
证据：`data/wirsum_alignment_stage_e_stiffness_ladder.json`

| 阶 | 状态 | 要点 |
|---|---|---|
| 等温固体输运 | 已实现 | 线性 holdup；后向误差 ~1e-17 |
| 等温双相交换 | 已实现 | Hamel Eq. 2.2 示踪；闭式解；raw L2 与 scaled RMS 分名 |
| 交换 + 生产 R12 | 已实现 | `rate_scale<1` 是延拓；只有 `scale=1` 可标物理收敛。直接入口初值线搜索失败 |
| 全气相 R5–R8+R12 | 已实现 | 7 物种 × 两相；无 tar/H2S。元素缺口在 Newton 收敛时检查 |
| 炭 R1–R4 | 已实现 | 库存 [kg] + `m=KM`；表面积 `6M/(ρ d_p)` |
| 能量 | 已实现（未物理收敛） | 完整速率停在 `singular_jacobian`，已分名；不含固体瞬时 cp |

禁止 fastox。小步长 / 线搜索失败单独报告。本核不经 `Reactor`，不改生产默认。`R9/R10/R11` 与固体 cp 积分仍未接入。

## 阶段 F — 床–自由空间界面（关回流）

核：`src/solvers/wirsum_1997_bed_freeboard.py`  
证据：`data/wirsum_alignment_stage_f_bed_freeboard.json`  
测试：`tests/test_wirsum_1997_bed_freeboard.py`

| 项 | 行为 |
|---|---|
| 气界面 | `N_fb_in = N_d + N_b`，相分裂不进自由空间 |
| 抛射 | `wirsum_1997_exact` 只用 Roof 3.59；`hamel_extended` 加 Wake 3.54。本几何 Wake 使抛射多约 1.01 kg/s |
| 命运 | 粗粒 2.25 mm → `returns`；细粉 0.20 mm → `reaches_next_cell`。阻力仍走 Wirsum `c_w`，不是 Hamel `z_exp` |
| 关回流 | `m_recycle=0`；抛射 = 回落 + 离开；离开 = 捕集 + 飞灰 |
| 开回流 | `m_recycle = η m_leave`，`η=0.95`。常数效率对照，不是旋风几何 |

生产 `simulate_freeboard` 未改。外层仍不刷自由板轨迹。

## 阶段 F 后半 — 气/T 内层与 Lapple 旋风

核：`src/physics/wirsum_1997_cyclone.py`，`src/solvers/wirsum_1997_bed_freeboard_coord.py`  
证据：`data/wirsum_alignment_stage_f2_bed_freeboard_coord.json`  
测试：`tests/test_wirsum_1997_bed_freeboard_coord.py`

| 项 | 行为 |
|---|---|
| 内层气 | 示踪交换闭式解（Hamel Eq. 2.2）；`N_fb = N_d + N_b`；不重算抛射 |
| 内层 T | `N cp_N2(T_in) (T_in − T) − Q_W = 0`。`Q_W=60 kW` 时 T：1100 → 1008 K |
| 旋风（对照，非主线） | `engineering_lapple_1951`。**不是** Hamel / 诊断口径 |
| T→η | 热损使 `d50` 从 0.736 μm 升到 0.769 μm；原切割粒径处 η 从 0.5 降到 0.478 |
| 反应量 | F2 仍为 `named_absent`。F3 接下。 |

Hamel 的 Zyklon 是 Bild 7.7 的 1–2 cell 拓扑，VTT 算例没有旋风几何。Lapple 式不得写成 Hamel 原文。生产捕集分数未改。

## 阶段 F3 — 离开炭与内层 R1 源汇

核：`src/solvers/wirsum_1997_reaction_coord.py`  
证据：`data/wirsum_alignment_stage_f3_reaction_coord.json`  
测试：`tests/test_wirsum_1997_reaction_coord.py`

| 项 | 行为 |
|---|---|
| 外层 | 细粉 Hamel 扩展抛射，`m_leave,char = 0.830 kg/s` |
| 库存 | `M` 与 `K_fb M + MW_C r` 闭合。`K_fb=1 1/s` 具名 `engineering_K_fb_leave` |
| 反应 | 只 `rate_R1`。面积 `6M/(ρ d_p)`。气源与阶段 E 的 R1 列相同 |
| 协调 | `r` 只来自 `char_r1_extent`。入口 `C_O2` 上 e1=8.39 mol/s，根上 0.476 mol/s |
| Newton | 完整速率 3 步，`‖F‖₂=3.5×10⁻¹¹`，可标物理收敛。气相 C 与 e1 差 ~10⁻¹¹ |

未做 R2–R4、均相。旋风 cell 见 F4。生产路径未改。

## 阶段 F4 — Hamel 旋风 cell（Bild 7.7）

核：`src/solvers/wirsum_1997_cyclone_cells.py`  
证据：`data/wirsum_alignment_stage_f_cyclone_cells.json`  
测试：`tests/test_wirsum_1997_cyclone_cells.py`

| 项 | 行为 |
|---|---|
| 拓扑 | `fb → cyclone(n=1\|2) → return → bed`。`n` 自选 |
| 气 | Eq. 2.4，单相，`ex=0`，`r=0`。`N` 逐格等于 `N_fb` |
| 固 | Eq. 2.6，`m=KM`。离开 = 飞灰 + 下行 |
| 回料 | 关：`m_rez=0`，`∂F_bed/∂M_ret=0`。开：`m_rez` 进床残差，`∂F_bed/∂M_ret=K_return` |
| K | `frozen_context`。论文没有旋风 K 公式；不是 Lapple，不是常数 η |

F2 的 Lapple 只作工程对照，不进本残差，不是主线。

## 阶段 F5 — 界面与旋风 cell 交接

核：`src/solvers/wirsum_1997_bed_fb_cyclone.py`  
证据：`data/wirsum_alignment_stage_f5_bed_fb_cyclone.json`  
测试：`tests/test_wirsum_1997_bed_fb_cyclone.py`

| 项 | 行为 |
|---|---|
| 气 | `N_fb = Σ(N_d+N_b)` 只算一次，旋风格 0 的 `N` 与它相同 |
| 固 | `m_leave` 来自 Hamel 扩展抛射 + 命运，原样进 Eq. 2.6 |
| 粗粒 | 回落，`m_leave=0`，旋风固相入口为 0 |
| 细粉 | `m_leave>0`；关回流 `m_rez=0`；开回流 `m_rez` 进床残差 |
| 不用 | 常数 η、Lapple、`evaluate_interface` 的旋风切分 |

待决问题与合格答案格式：`docs/handoff_wirsum_af_open_questions_2026-09-14.md`。在 Q1/Q2/Q6 未裁定前，不把分核硬耦成一条 \(F\)，不补旋风经验 \(\eta\)。

## 明确不做

- 不叠加 bed0 热解气初值当主线
- 不改 `thesis_vorab_sources_single_shot` 默认（阶段 D）
- 不把默认生产闭合改成 species_own（避免 silently 改 LU）
- 不把 residual clip / fastox 当修正
- 不跑 150-iter LU
