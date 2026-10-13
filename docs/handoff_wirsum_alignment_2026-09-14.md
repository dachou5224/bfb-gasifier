# 交接：Wirsum 对齐（reduced LU x0）与 VM 口径

**日期：** 2026-09-14  
**分支：** `feat/nr-hamel-wirsum-core-startwert`（相对 origin 已 ahead；本链从 P0 起约 22 个 commit）  
**状态：** 未跑全炉 150-iter LU；未宣称收敛。已把 x0 上 \(10^{13}\) 残差拆到能量/热解气缺口。  
**语言：** 接手用简体中文。先读 `docs/CLAUDE.md`。  
**原始交接：** `/Users/liuzhen/Downloads/cursor_handoff_wirsum_alignment.md`

---

## 0. 60 秒摘要

目标不是「调参让 LU 收敛」，而是重建：

**Wirsum 原始数学问题 → 状态变量 → residual → Newton → freeboard closure**

本轮在 **reduced 3-cell LU x0**（`build_initialized_lu_reactor`，只 `generate_initial_x0`，不走 `init_precalc`）上做诊断核 + A/B。

因果链已经闭合到这一步：

1. Vorabrechnung 把热解 H₂ 和剩余 O₂ 叠进同一相 → 裸 thesis R12（仓库编号，Hamel R6）\(\xi\sim 10^{13}\) mol/s，packed \(\|F\|_2\sim 10^{13}\)。
2. 同相 R12 Startwert（只改 holdup）后 \(\|F\|_2\sim 10^{7}\)，linf 是 **bed0 能量**（约 −10 MW）。能量在 R12 之前就在，只是被气相天文量级盖住。
3. −10 MW 的第一层：`holdup_transport` 下 \(H_{\mathrm{aus,solid}}=K\cdot M\)，而 `generate_initial_x0` 时 `m_zu` 仍为 0，故 \(M=0\)、\(H_{\mathrm{sout}}=0\)，给料焓仍进 \(H_{\mathrm{zu}}\)。
4. `M=m_{\mathrm{in}}/K` 只填 **char/ash** 后，bed0 \(E\) 到 −5.3 MW。轴向接力后 bed1 不再空库存。
5. 剩下的主导桶是 **VM 以固体列进 \(H_{\mathrm{zu}}\)（−3.09 MW）且 \(K_{\mathrm{VM}}=0\)**。Hamel 原文里 Flüchtige **不是** Eq.(2-6) 固体种，释放后应并入气相物质/能量衡算。

**下一手：** 先答 `docs/handoff_wirsum_af_open_questions_2026-09-14.md` §4（Q1–Q9）。A–F 约化核未完成诊断验收。禁止给 VM 加 K，禁止 fastox，不跑全炉 LU。

P0–P5 改了什么、哪些进了生产路径：§2A。P6 数字：§3。

---

## 1. 硬约束（接手不得违反）

| 项 | 规定 |
|---|---|
| 单位 | 全程 SI。Arrhenius 只走 `src/kinetics/arrhenius.py` |
| residual | **禁止** clip / projection / floor / anchor / T fence 掩盖奇异或残差 |
| 命名 | `thesis-exact` 与 `engineering` 必须分开。已有 fastox / 线搜索投影 **不要**改名成 Wirsum 原文 |
| 提交 | 小步、可回退、每步独立 commit；**只加本阶段文件**。勿卷 `SCRIPTS_INDEX.md`、未跟踪 `src/core/freeboard/`、`jacobian_analytic.py`、整份 WIP `cell.py` / `outer_loop.py` / `vorab_x0.py` |
| 验收 | **不用实验 CO / 工业测点**给 numerical bridge 背书 |
| VM | **不**给 VM 加 \(K\)；**不**把 VM 补进 Eq.(2.7) `H_aus` 固体载体 |
| 求解 | 默认不跑 `max_global_iter=150` 的 Phase1 HTW；完整 `pytest -q` 会卡在该路径 |
| 默认 fixture | `build_initialized_lu_reactor` **保持**「只 generate_initial_x0、不种 holdup/R12」。诊断核在自己的 payload 里叠加 Startwert，以便 P0/P5 证据仍可复现 |

相关技能：`.cursor/skills/bfb-diagnosing-bugs`、`bfb-hamel-parity-audit`、`bfb-numerical-guard`、`bfb-scientific-audit`。

---

## 2. 本链已完成（独立核，生产路径大多未接）

生产 `solve_global_nr` / `simulate_freeboard` / `outer_loop.py` **未改成**本核默认。`init_precalc` 里已有更大的工程 fastox，reduced fixture 不走那条路。

| 阶段 | 提交（短） | 核 / 证据 | 结论 |
|---|---|---|---|
| P0 | `f31fc9a` `61a9bfd` `3ad1f1f` | baseline JSON、knob plumbing、`global_residual` 对 packed `x` 纯函数 | `snap_trace_o2_holdups` 已移出 residual |
| P1 | `92f410e` `b000a6c` | `src/core/solid_balance_reduced.py`，`docs/solid_balance_formulation_audit.md` | Eq.(2-6) 约束的是质量流。`K=0` 且反应不依赖 `M` 时雅可比为零。`holdup_k0_anchor='engineering'` 是 `source−M`（\(\tau=1\,\mathrm{s}\)），不是 Hamel |
| P2 | `0c5b35b` `477cde2` | `src/solvers/wirsum_1997_strict.py` | \(J\Delta x=-F\)，\(\alpha=1,1/2,\ldots\)，仅 raw \(\|F\|_2\)。生产 NR 仍用缩放 RMS |
| P3 | `f43215d` | `src/physics/wirsum_1997_freeboard.py` | `wirsum_1997_exact` vs `engineering_freeboard` 分名。生产 `simulate_freeboard` 未改 |
| P4 | `9ba52c3` `897a3ab` | `src/solvers/wirsum_1997_outer_diagnostics.py` | inner 收敛 ≠ outer 收缩 |
| P5 | `6347493` `fcd3ed8` | `src/solvers/wirsum_1997_lu_reentry.py` | 8 项清单；x0 \(\|F\|_2=9.29\times10^{13}\)；open = elemental / mass / energy |
| 热解 An=b | `3022c54` `8bb3f99` | `hegermann_yield.py`，`pyrolysis_carbon_split.py` | Hegermann 后先分 Koks/Teer 再 Gibbs；禁止事后把剩余 C 划回 Koks |
| 元素分解 | `64cddf9` `9291a9e` | `wirsum_1997_elemental_decomp.py` | H/O 净口 ~\(10^{-2}\)；R12 元素中性，只解释质量残差 |
| 同相 R12 Startwert | `3e108c4` | `wirsum_1997_startwert_r12.py`，`data/wirsum_alignment_p6_startwert_r12.json` | 只做 `H2+0.5O2→H2O`，不做 R5/R6/R10/R1 |
| 残差族 | `42d1ec5` | `wirsum_1997_residual_decomp.py` | R12 后 linf = `c0.energy`，\(\|F\|_\infty/\|F\|_2=0.993\) |
| Eq.2.7 分项 | `c958831` | `wirsum_1997_energy_decomp.py` | \(H_{\mathrm{sout}}=0\) 因 \(M=0\) |
| holdup Startwert | `49d20a1` | `wirsum_1997_holdup_startwert.py` | \(M=m_{\mathrm{in}}/K\)，仅 \(K>0\)（char/ash） |
| 剩余能量桶 | `9ad1567` | `wirsum_1997_energy_leftover.py` | 桶之和 = live \(E\) |
| 轴向接力 | `098c8cf` | 同 holdup 模块 `apply_km_holdup_startwert_axial_sweep` | bed1/2 不再空库存 |
| VM 口径（对话，无新 commit） | — | Hamel PAGE 18–19、§4.3 PAGE 69–70 | 见 §5 |
| F 界面 | 未 commit | `wirsum_1997_bed_freeboard.py`，`data/wirsum_alignment_stage_f_bed_freeboard.json` | 关回流；`N_d+N_b`；Wirsum Roof ≠ Hamel Roof+Wake；粗粒回落、细粉离床后 `η=0.95` |
| F 后半 | 未 commit | `wirsum_1997_cyclone.py`，`wirsum_1997_bed_freeboard_coord.py` | 内层示踪+焓 T；Lapple `η(Q)`；反应量 `named_absent` |
| F3 | 未 commit | `wirsum_1997_reaction_coord.py` | 离开炭 + `rate_R1`；入口 e1≠根；`‖F‖₂~10⁻¹¹` |
| F4 | 未 commit | `wirsum_1997_cyclone_cells.py` | Bild 7.7；`m_rez` 进床残差；K 冻结；不用 Lapple |
| F5 | 未 commit | `wirsum_1997_bed_fb_cyclone.py` | 同一 `N_fb`/`m_leave` 交接；粗粒不进旋风 |

明细见下节 §2A。P6 以后的 Startwert / 能量桶仍见 §3。

---

## 2A. P0–P5 修改明细（给诊断 agent）

原始交接要求：**先重建数学问题，再谈收敛**。P0–P5 只做这件事。生产 `solve_global_nr` / `simulate_freeboard` / `outer_loop.py` **没有**改成本核默认。真正改到生产路径的，只有 P0.3（residual 纯度）和 P1.2（`holdup_k0_anchor` 开关，默认仍是 `engineering`）。

### 2A.0 生产路径改动一览

| 阶段 | 改了的生产文件 | 行为变化 |
|---|---|---|
| P0.3 | `src/solvers/global_nr_iteration.py`，`src/solvers/global_nr_solver.py` | `snap_trace_o2_holdups` 从 `F(x)` 内挪到 **pre-NR** 与 **accepted-step**。Jacobian / 线搜索评 `F` 时不再偷偷改 holdup |
| P1.2 | `src/core/cell_balances.py`，`src/core/cell.py`，`src/core/reactor.py`，`src/gasifier/models/balances.py` | 新增 `holdup_k0_anchor`。默认 `'engineering'` = 原 `source−M`，LU 数值不变。`'none'` 时 `K=0` 残差保持 [kg/s]，允许零行 |
| P2 | `src/core/reactor.py` | 新增 `nr_solver_mode`。`'wirsum_1997_strict'` 只关一批 config 工程桥，**不改** `solve_global_nr` 的缩放 RMS 接受准则 |
| P3–P5 | 无生产求解路径改动 | 独立核 + 测试 + JSON |

其余 P0–P5 文件都是诊断核、测试、冻结脚本或 docs。不要把它们当成已经接到全炉。

---

### 2A.1 P0 — 冻结现场，不改公式

**目标：** 改公式或求解器之前，先把当前工程栈的数字钉死。不追全炉收敛。

#### P0.1 baseline（`f31fc9a`）

| 新增 | 作用 |
|---|---|
| `scripts/freeze_wirsum_alignment_baseline.py` | 写 `data/wirsum_alignment_p0_baseline.json` |
| `tests/test_wirsum_alignment_baseline.py` | 回归：reduced 3-cell x0 的 packed 范数、Jacobian 结构 |
| `data/wirsum_alignment_p0_baseline.json` | git `fdeee0e` 时刻的工程栈，**不是** thesis-exact 声明 |

reduced 3-cell LU **只 x0、不跑 Newton**（`n_var=105`）：

| 量 | 值 |
|---|---|
| packed \(\|F\|_2\) | \(9.29\times10^{13}\) |
| packed \(\|F\|_\infty\) | \(3.29\times10^{13}\) |
| rms_scaled | \(1.13\times10^{11}\) |
| Jacobian nnz / 零行 / 零列 | 1294 / 0 / 0 |
| \(\kappa_2\)（dense） | \(2.33\times10^{22}\) |
| 轴向 T [K] | 1325, 1232.25, 1139.5 |
| 出口湿摩尔分率 | H₂ 0.085、CO 0.117、CO₂ 0.020、CH₄ 0.012 |

同一 JSON 里还有 `phase1_htw_lu_solve`（工程栈跑了 3 个 inner / 2 个 outer，**未收敛**，\(\|F\|_2\sim10^9\)）和 Phase2 3+2 smoke。Phase2 smoke 已出现 **4 个零行、2 个零列**——后来 P1 证明这是 `K=0` 公式奇异，不是实现 bug。

当时 `tests/sanity_checks.py` 失败项只有 **「10. Pyrolysis elemental allocator」**。P5 热解 An=b 才处理这件事；不要用 P0 的 sanity 失败去改 residual。

#### P0.2 knob plumbing（`61a9bfd`）

只加 `tests/test_wirsum_knob_plumbing.py`。证明 `ReactorConfig` 上若干开关能走到 Cell / 速率 / 残差，避免后面 Wirsum-strict 关桥时改到死配置。P1.2 给同一测试补了 `holdup_k0_anchor` 断言。

#### P0.3 residual 纯度（`3ad1f1f`）——本阶段唯一直接改 NR 循环的提交

**问题：** `snap_trace_o2_holdups` 原先在评 `F(x)` 时改 holdup。Jacobian 差分和线搜索看到的不是同一个纯函数。

**改法：**

- `global_nr_solver.snap_trace_o2_holdups` 仍在，但是 **投影**，不是 residual 的一部分。
- `global_nr_iteration.py`：只在 **NR 开始前** 和 **接受步之后** 调用。
- `tests/test_wirsum_residual_purity.py`：同一 packed `x` 连评两次 `F`，holdup 不变。

**接手约束：** 不要把 snap 塞回 `global_residual`。需要 snap 时走 pre-step / accepted-step。

---

### 2A.2 P1 — 固相未知量是流，不是库存

**目标：** 对照 Hamel Eq.(2-6) 写清：论文约束的是 [kg/s] 物流；库存 `M` [kg] 只在本构 \(m=KM\) 之后才进入方程。

#### P1.1 约化三格问题（`92f410e`）

| 新增 | 作用 |
|---|---|
| `src/core/solid_balance_reduced.py` | 三格、单物种、单粒径，**不经过** `Reactor` |
| `tests/test_solid_balance_reduced.py` | Model A（流）vs Model B（库存）秩、条件数、稳态 |
| `docs/solid_balance_formulation_audit.md` | 公式、零空间、禁止事项 |

约化结论（给料 1 kg/s 进床底）：

| 提法 | rank | 稳态 |
|---|---|---|
| A 流平衡 `m_out` | 3，cond≈5.87 | 唯一：`m_out=(5/3, 5/3, 1)` kg/s |
| B 库存，`K>0` | 3，cond≈5.71 | 唯一：`M*=(4.167, 3.333, 2.000)` kg；`K_sum M*` 与 A 的流相同 |
| B，`K=0`，无反应，给料=0 | **0**，3 零行零列 | 任意 `M` 都是解 |
| B，`K=0`，给料=1 kg/s | **0** | **无有限解**；`F=(1,0,0)` 与 `M` 无关 |
| B + `source−M` 锚 | 3，cond=1 | 人为唯一。单位混用 [kg/s]−[kg]，暗设 \(\tau=1\,\mathrm{s}\) |
| B + 一阶 `k_rxn=0.2 s⁻¹`，`K=0` | 3 | 合法满秩：`M=(5,0,0)` kg |

`K>0` 时 A 与 B **代数等价**。B 多出来的信息全部来自本构 \(m=KM\)，`K` 不能由 Eq.(2-6) 反推。

`K=0` 且 \(\dot r\) 不依赖 `M` 时，零行是 **数学问题没闭合**，不是 Newton 病态。合法处理只有：删自由度、补单位为 [1/s] 的本构、或承认 `M` 不定。**禁止**把残差改成 `source−M`。

生产代码里原先就有这条锚：`_apply_holdup_transport_numerical_guardrails` 在 \(K_{\mathrm{sum}}\le10^{-14}\,\mathrm{s}^{-1}\) 时写 `F=net_source-M`；`jacobian_analytic._holdup_solid_self_derivative` 对角回 `-1.0`。单桶对照：char、`K=0`、入流 1 kg/s、`M=0.25` kg → 生产残差 = 0.75 = `1 kg/s − 0.25 kg`。物理残差应是 1 kg/s，与库存无关。

`legacy_stream` 分支同样执行 `… − m_solid`，也是工程路径。

#### P1.2 显式开关（`b000a6c`）——生产默认不变

`calc_solid_balance_residual(..., k0_anchor=)`，经 `ReactorConfig.holdup_k0_anchor` 拷到 `Cell.holdup_k0_anchor`。

| 值 | 行为 | 用途 |
|---|---|---|
| `'engineering'`（**默认**） | 保持原 `source−M` / disconnected `-M` | 现有 LU 数字不变 |
| `'none'` | `K=0` 时残差 = 净源 [kg/s]，允许零行 | thesis / `wirsum_1997_strict` |

测试：`tests/test_cell_balances.py`，`tests/test_wirsum_knob_plumbing.py`。

**接手约束：** 不要把默认改成 `'none'` 却宣称「LU 行为没变」。不要在 residual 里再加第四种锚。P6 的 \(M=m_{\mathrm{in}}/K\) 只填 **`K>0` 的 char/ash**，VM/moisture 保持 0，就是按本审计办的。

---

### 2A.3 P2 — Wirsum §2.3.3 最小阻尼 Newton

**目标：** 在约化问题上实现原文 Eq.2.21–2.23，并与工程缩放 RMS 对照。生产 `solve_global_nr` **继续**用缩放 RMS。

#### 核（`0c5b35b`）

`src/solvers/wirsum_1997_strict.py`：

1. 组装 \(F(x), J(x)\)
2. 解 \(J\Delta x=-F\)（Eq.2.21）
3. 试探 \(\alpha=1,1/2,1/4,\ldots\)
4. 仅当 merit **严格下降** 时接受。原文 Eq.2.23 的 merit 是 raw \(\|F\|_2\)

**不含** T fence、fastox、syngas floor、split/energy merit、refresh-aware accept、residual clip。

`apply_wirsum_1997_strict_config(cfg)` 关掉一批可由 config 开的工程桥，并强制 `holdup_k0_anchor='none'`。它 **不**把 `solve_global_nr` 的接受准则改成 \(\|F\|_2\)。要走原文线搜索，必须调 `wirsum_1997_strict_newton`。

`ReactorConfig.nr_solver_mode`：`'engineering'` | `'wirsum_1997_strict'`。

说明：`docs/wirsum_1997_strict_newton.md`。测试：`tests/test_wirsum_1997_strict.py`。

#### A/B（`477cde2`）

机器表：`data/wirsum_1997_strict_ab.json`。三格 holdup，`x0=0.05` kg。

| 问题 | raw+fresh | scaled+fresh | raw+frozen | scaled+frozen |
|---|---|---|---|---|
| 线性 `K>0` | 1 步，α=1 | 同左 | 同左（J 常数） | 同左 |
| 非线性 MM 汇 | 5 步，α=1 | 同左 | 40 步未到 `1e-10` | 同 frozen |
| MM 行加权 \((1,10^3,10^{-3})\) | 8 步，α 含 1/8 | 5 步，α=1 | **线搜索失败** | 40 步爬行 |
| `K=0` 无出流本构 | `singular_jacobian`，rank 0，3 零行 | — | — | — |

结论：

- 线性闭合问题不需要缩放或刷新 J。
- 非线性但尺度均匀时，难度在 **frozen J**，不在 raw vs scaled。
- 残差行尺度差 \(10^6\) 时，raw \(\|F\|_2\) 会强迫回溯；scaled RMS 仍满步。
- `K=0` 零行是 P1 的公式奇异，Newton 实现正确地报出来。

**接手约束：** 不要为了让 raw \(\|F\|_2\) 好看去缩放 / clip residual。不要把生产默认 merit 改成 raw 却继续用工程桥。

---

### 2A.4 P3 — 自由板：原文核与工程别名分名

**目标：** 把 Wirsum 原文抛射/阻力与现用工程自由板拆开。**不改** `simulate_freeboard`。

`src/physics/wirsum_1997_freeboard.py`（`f43215d`），说明 `docs/wirsum_1997_freeboard.md`，测试 `tests/test_wirsum_1997_freeboard.py`。

| `freeboard_model` | 平均抛射 | 床面抛射 | 阻力 |
|---|---|---|---|
| `wirsum_1997_exact` | \(1.2\,u_b\) | 仅 Roof Eq.3.59，**不** renormalize | \(c_w=24/\mathrm{Re}+z\)，\(z\in[0.43,1.0]\)，默认 0.44 |
| `hamel_extended` | \(1.53\,u_b\)（Martens / Hamel 3.73） | Roof + Wake，不 renormalize | Hamel \(z(\varphi)\)（本核未展开） |
| `engineering_freeboard` | \(1.53\,u_b\) | Roof + Wake，再缩放到 available feed | Haider 混合 |
| `analytical_wirsum` | **别名 → `engineering_freeboard`** | 同上 | 同上 |

原文轨迹：竖直速度高斯 \(\sigma=0.6\,\mathrm{m/s}\)；Roof \(F_{0,r}=d_p\rho_p(1-\varepsilon_d)(3\varepsilon_b u_b/d_b)\) [kg/(m² s)]；\(\mathrm{d}u_p/\mathrm{d}t=Au_p^2+Bu_p+C\)。

`exact_hamel` / Haider \(z(\varphi)=73.69\exp(-5.0784\varphi)\) **不是** `wirsum_1997_exact`。P5 清单第 6 项确认：生产名仍解析到 `engineering_freeboard`；`wirsum_1997_exact` **未接入全炉**。

**接手约束：** 不要把 `analytical_wirsum` 改名成 exact。不要在未对照 Roof 3.59 的情况下改生产抛射系数。

---

### 2A.5 P4 — 外层 Vorabrechnung：inner 降残差 ≠ outer 收缩

**目标：** Wirsum §2.2 的固定点是 \(h_{k+1}=H(x^*(h_k))\)。只记录场增量，**不加** 松弛 / Anderson / inner clip。生产 `outer_loop.py` 未改。

`src/solvers/wirsum_1997_outer_diagnostics.py`（`9ba52c3`），说明 `docs/wirsum_1997_outer_diagnostics.md`。

记录相邻 outer 的 \(\|\varepsilon_b\|_{2}\)、\(\|K_{bd}\|\)、\(\|K_s\|\)、\(\|u_b\|\)、\(\|T\|\)，以及当步 inner residual。\(K_s\) 取 `mean(K_solid_auf)+mean(K_solid_ab)` [1/s]，只作输运强度诊断。收缩比 \(\rho\approx\|\Delta h_{k+1}\|/\|\Delta h_k\|\)。`verdict.non_contractive`：尾部水力场 \(\rho\) 均值 > 1。

A/B（`897a3ab`，`data/outer_iteration_history.json`）：三格 **仿射 map**，不跑 LU。

| map | \(\rho\) | inner residual | `non_contractive` |
|---|---|---|---|
| \(h^*+0.4(h-h^*)\) | 0.4 | 每步 ×0.5 | false |
| \(h^*+1.3(h-h^*)\) | 1.3 | 每步 ×0.5 | true |

两边 inner residual **同步下降**。扩张 map 仍标成 non-contractive。一次快照（P5 清单第 5 项）只能记下场，**不能**算 \(\rho\)。

**接手约束：** 不要因为 inner \(\|F\|\) 下降就宣布 outer 收敛。不要在本核里加阻尼来「修」\(\rho>1\)。

---

### 2A.6 P5 — 在干净架构上重诊 LU x0（不追收敛）

**目标：** 按固定顺序看闭合，**第 7、8 项只对照 P0 baseline**，不用实验 CO/T。默认只看 Vorabrechnung x0，不跑 `reactor.solve()`。

#### 八项清单（`6347493`，冻结 `fcd3ed8`）

`src/solvers/wirsum_1997_lu_reentry.py`，说明 `docs/wirsum_1997_lu_reentry.md`，证据 `data/wirsum_alignment_p5_lu_reentry.json`。

禁止出现在 composition / temperature 判据里的键：`experiment_CO`、`measured_CO`、`y_CO_target`、`calibration_objective`、`lu_plant_parity`。

reduced 3-cell x0 当时的开口（**尚未**叠加 P6 R12/holdup Startwert）：

| 项 | 状态 | 数字 |
|---|---|---|
| 1 elemental | **open** | C ~ \(10^{-12}\) mol-atom/s；H/O max abs ≈ 0.0078 |
| 2 mass | **open** | 气相 residual L2：bed0 \(6.04\times10^{13}\)，bed1 \(6.10\times10^{13}\)，bed2 \(3.55\times10^{13}\)；固相 bed0 0.50 kg/s，bed1/2 = 0 |
| 3 energy | **open** | bed0 \(E\approx-10.72\) MW（P6 才拆因） |
| 4 inner NR | 记录 | packed \(\|F\|_2=9.29\times10^{13}\)，与 P0 同类 |
| 5 outer | 单帧 | 一次快照，不算 \(\rho\) |
| 6 freeboard | 命名 | 生产仍是 `engineering_freeboard` |
| 7 composition | observation | vs P0 出口摩尔分率，不是测点 |
| 8 temperature | observation | vs P0 轴向 T |

`open_items = [elemental_closure, mass_closure, energy_closure]`。

#### 热解 An=b（`3022c54`，`8bb3f99`）——P0 sanity 第 10 项

| 文件 | 作用 |
|---|---|
| `src/thermal/hegermann_yield.py` | Hegermann Tabelle 4.3 / Eq.4.21–4.22 的 Koks/Teer 碳分率，与生产源项路径对齐 |
| `src/core/pyrolysis_carbon_split.py` | **先**按 Hegermann 切 Koks/Teer，**再**对剩余元素做 Gibbs `An=b` |
| `tests/sanity_checks.py` | 第 10 项改为同一顺序 |

sanity 池 `(10,12,4)` 在 `koks=0` 时 `An=b` 不可行。**禁止**事后把剩余 C 划回 Koks。这是分配器可行性，不是 packed \(10^{13}\) 的原因。

#### 元素分解（`64cddf9`，冻结 `9291a9e`）

`src/solvers/wirsum_1997_elemental_decomp.py`，说明 `docs/wirsum_1997_elemental_decomp.md`，证据 `data/wirsum_alignment_p5_elemental_decomp.json`。

两件必须分开的事：

1. **质量残差 \(10^{13}\)**：裸 thesis R6 / 仓库 R12（Hamel Eq.5.38–5.39），2.5 MPa、~1300 K，\(\xi\sim5.64\times10^{13}\) mol/s。`H2+0.5O2→H2O` **元素中性**，`leftover_over_r12_extent ~ 10⁻¹⁶`。去掉 R12 后 H/O 净口几乎不动。
2. **H/O 净口 ~ \(10^{-2}\)**：bed0 气相原子流 ~ 35 mol-atom H/s，固相 VM+水分约 −35，相对差 ~ \(2\times10^{-4}\)。这条 init 上 `_vm_*_cache` **空**（`pyrolysis_cache_used=false`），活热解写在 `R_solid`。bed1/2 无固相源，净口来自均相反应大数对消。

P6 的同相 R12 Startwert 是接着本条做的：错的是 x0 把 H₂ 与 O₂ 叠在同一相，不是 R12 速率式本身。

---

### 2A.7 P0–P5 复现（不跑 150-iter LU）

```bash
python -m pytest \
  tests/test_wirsum_alignment_baseline.py \
  tests/test_wirsum_knob_plumbing.py \
  tests/test_wirsum_residual_purity.py \
  tests/test_solid_balance_reduced.py \
  tests/test_wirsum_1997_strict.py \
  tests/test_wirsum_1997_freeboard.py \
  tests/test_wirsum_1997_outer_diagnostics.py \
  tests/test_wirsum_1997_lu_reentry.py \
  tests/test_hegermann_yield.py \
  tests/test_pyrolysis_carbon_split.py \
  tests/test_wirsum_1997_elemental_decomp.py -q
```

冻结脚本：`scripts/freeze_wirsum_alignment_baseline.py`，`scripts/freeze_wirsum_1997_{strict_ab,outer_history,lu_reentry,elemental_decomp}.py`。

P0–P5 之后的因果链（R12 → 能量 → holdup → VM）见 §3，不要从本节省略到「直接调 T / 加 \(K_{\mathrm{VM}}\)」。

P6 复现（不跑 Newton；P0–P5 见 §2A.7）：

```bash
python -m pytest tests/test_wirsum_1997_startwert_r12.py \
  tests/test_wirsum_1997_residual_decomp.py \
  tests/test_wirsum_1997_energy_decomp.py \
  tests/test_wirsum_1997_holdup_startwert.py \
  tests/test_wirsum_1997_energy_leftover.py -q
```

冻结脚本：`scripts/freeze_wirsum_1997_{startwert_r12,residual_decomp,energy_decomp,holdup_startwert,energy_leftover,holdup_sweep}.py`

---

## 3. 关键数值（reduced 3-cell LU x0）

工况来自 `tests/global_nr_fixtures.build_initialized_lu_reactor` → `build_phase1_htw_lu_reactor_config`。bed0 \(T=1325\,\mathrm{K}\)，`heat_loss_mode=hamel_eq27_explicit`，\(Q_W\approx 0.688\,\mathrm{MW}\)。

### 3.1 R12 / packed residual

| 量 | 关 R12 Startwert | 开同相 R12 |
|---|---|---|
| 同相 H₂∩O₂ \(\xi\) | 2.36 mol/s | 0 |
| 源项 R12 \(\xi\) | \(5.68\times10^{13}\) mol/s | 7.73 mol/s |
| packed \(\|F\|_2\) | \(9.29\times10^{13}\) | \(1.05\times10^{7}\) |
| 主导方程 | 气相 `c1.N_d.H2O` | **`c0.energy`** |
| 能量 L2 | 已是 \(1.07\times10^{7}\) W | \(1.05\times10^{7}\) W |
| 气相 L2 | \(9.29\times10^{13}\) | \(6.05\times10^{4}\)（bed0）；bed1/2 ~11 |
| R5 \(\xi\) | \(\sim 5.5\times10^{4}\) | \(\sim 2.0\times10^{4}\)，**不是**计量对 |

R12 速率：`rate_R12`，扫描原式 \(r=k C_{\mathrm{H2}}^{3/2} C_{\mathrm{O2}}\)，\(k=1.631\times10^{9}T^{3/2}\exp(-3420/T)\)（Hamel Eq.5.38–5.39；thesis R6 = 仓库 R12）。旧转录 \(C_{\mathrm{H2}}C_{\mathrm{O2}}^{1/2}\) 已于 2026-09-14 撤回。巨大速率仍在，阶数修正不单独保证收敛。化学计量元素中性，不制造 H/O 净口。

### 3.2 为何 x0 会 H₂∩O₂

`generate_initial_x0` 无 `cfg` 时走 legacy：`o2_remaining=O2_feed·exp(-5ξ)`，再把剩余 O₂ **和** 热解 H₂ 叠进同一 holdup。  
生产 `init_precalc` 在其后调用 `project_bed_holdup_fast_oxidation_closure`（R12→R5→R6→R10+可选 R1）。那是 **engineering fastox**，与本核 `same_phase_r12_startwert` 分名。  
Newton 线搜索投影 `nr_line_search_fast_oxidation_projection_thesis` 是第三条路径，默认关。

### 3.3 Eq.2.7（R12 后、holdup 种之前）

`E=H_{\mathrm{in}}-H_{\mathrm{out}}-Q`，与 `calc_energy_balance` 对齐（误差 <1 W）。

bed0：

| 项 | 值 |
|---|---|
| \(H_{\mathrm{in}}\) | −11.20 MW |
| \(H_{\mathrm{zu,solid}}\) | −8.20 MW（char −5.11 + VM −3.09） |
| \(H_{\mathrm{zu,gas}}\) | −2.79 MW |
| \(H_{\mathrm{sout}}\) | **0** |
| \(H_{\mathrm{gout}}\) | −1.16 MW |
| \(Q_W\) | 0.69 MW |
| \(E\) | **−10.72 MW** |
| `m_solid` | **全 0** |
| \(K_{\mathrm{auf,char/ash}}\) | \(\approx 0.079\,\mathrm{s}^{-1}\) |
| \(K_{\mathrm{VM}}\) | **0** |

`generate_initial_x0` 用 `m_zu+m_in+m_rez` 种子；当时三者为 0。BC 之后才有 `m_zu`（char 0.322、VM 0.369、水分 0.159、ash 0.089 kg/s），库存未回填。

### 3.4 `M=m_{\mathrm{in}}/K` 之后

| | 一次扫描 | 轴向接力（种 i → BC → 种 i+1） |
|---|---|---|
| bed0 \(E\) | −5.30 MW | −5.47 MW |
| bed1 \(E\) | −4.42 MW，`holdup_zero` | −2.71 MW，char 4.79 kg |
| bed2 \(E\) | −0.67 MW，`holdup_zero` | −0.59 MW |
| \(M_{\mathrm{VM}}\) | 0 | 0 |

bed0 在一次扫描后的可加总桶（`data/wirsum_alignment_p6_energy_leftover.json`）：

| 桶 | MW | \|份额\| |
|---|---|---|
| VM：\(H_{\mathrm{zu}}-H_{\mathrm{sout}}\) | **−3.09** | 50% |
| 气体：\(H_{\mathrm{in}}-H_{\mathrm{out}}\) | −1.84 | 30% |
| \(-Q_W\) | −0.69 | 11% |
| char 感温/质量差 | +0.45 | 7% |
| ash | −0.13 | 2% |

char 生成焓进出已对上（\(H_{\mathrm{sout,char}}\approx-5.55\) MW）。气体桶主要是冷气化剂进、出口气焓未带上热解原子。

---

## 4. 已定位的问题（按优先级）

### P-A. reduced x0 不走完整 Vorabrechnung

`build_initialized_lu_reactor` 只调用 `generate_initial_x0(...)`，**不传 `cfg` / `cell_budgets`**，不跑 `run_init_and_precalc_for_global_nr`。因此：

- 无 fastox、无 holdup 回填、无热解气 dense seed
- `_vm_*_cache` 在元素分解时未写入（`pyrolysis_cache_used=false`）
- 与生产 `HAMEL_WIRSUM_CORE` 的 x0 **不是同一对象**

诊断时必须写明看的是哪条链。不要用实验 CO 或 Phase2 o8 rms 直接给本链背书。

### P-B. 同相快氧化是 Startwert，不是 F(x)

裸 R12 在物理上「对」：H₂ 与 O₂ 同轨就该炸。错的是 **初值把不相容物种叠在一起**。本核只做同相 R12。CO∩O₂ 仍约 1.04 mol/s（bed0），R5 约 \(2\times10^{4}\) mol/s 且不是计量对——**不是** \(10^{7}\) packed 残差的主因。

### P-C. 固相未知量是流还是库存（P1，已书面化）

Hamel Eq.(2-6) 的未知量是 \(m_{\mathrm{zu}},m_{\mathrm{auf}},m_{\mathrm{ab}},\ldots\) [kg/s]。  
`holdup_transport` 用 \(m=KM\) 把未知量换成 \(M\) [kg]。\(K=0\) 时稳态不约束 \(M\)。  
`holdup_k0_anchor='engineering'` 用 `source−M` 人为满秩，**禁止**再当 thesis 默认。

### P-D. VM 被当成第四固体种记账（Hamel 不一致）

见 §5。能量上表现为：给料 VM 焓进 \(H_{\mathrm{zu}}\)，轴向 \(K_{\mathrm{VM}}=0\) 故不进 \(H_{\mathrm{aus}}\)，x0 气相 holdup 又没接住 §4.3 产物。

### P-E. 热解气 Startwert 生产路径半截、reduced 路径没有

`seed_pyrolysis_gas_dense_holdup_for_nr_x0`（`vorab_x0.py`）只在 `_vm_cache_valid` 且走 `init_precalc` 时运行；注释写明只种 **NH3/TAR** dense holdup，不是 §4.3 全套主气。  
reduced fixture 从未调用它。

### P-F. 工作区污染

下列 **不要整包提交**：`src/core/freeboard/`、`jacobian_analytic.py`、`src/solvers/wirsum_damping.py`、大 WIP `cell.py` / `cell_pyrolysis.py` / `outer_loop.py` / `vorab_x0.py`。  
`tests/global_nr_fixtures.py` 在本机曾是未跟踪文件；本链测试依赖它，接手先确认仓库里是否已跟踪。

---

## 5. Hamel 原文：VM / Flüchtige 口径（已核对，无新 commit）

**证据顺序：** 分册 → PDF 页 → 代码。禁止用代码反向补文档。

| 出处 | 论文在说什么 |
|---|---|
| PAGE 18，§2.2.2 | 固体种 = **活性固体**（Koks、Kalkstein、多种 Koks）× 粒径 \(k\)；另有惰性床料。举例里 **没有 VM** |
| Eq.(4.10)、(4-12)，PAGE 36–37 | \(m_v'\) = 颗粒内 **尚未释放** 的挥发分质量分数（DAEM），不是格间固体物流 |
| §4.3 PAGE 69 | 先算 Koks + Teer；**剩余元素全部进气相**；Gibbs 物种例：CO, CO₂, CH₄, C₂H₄, NH₃, H₂S, H₂O, H₂ |
| Kap.5 开篇 | 释放产物是 **Koks、Teer und Gas**，再与气化剂反应 |
| Eq.(2-4)/(2-5) | 热解气 = 气相源 \(N_{r,d/b,j}\) |
| Eq.(2.7)/(2.8) PAGE 19 | 焓绑在 **气体 + 活性固体** 物流上，含生成焓 |

对应代码：

- `S_VM`：给料/释放 **预算列**，不是 Feststoffsorte。
- `K_VM=0`、`m_auf/ab` 的 VM 置零（`solid_transport.py`）：与「VM 不轴向输运」一致，**保持**。
- `R_gas` ← `calc_drying_pyrolysis_sources` / Hegermann + Gibbs：结构与 §4.3 一致。
- **缺口：** x0 把未释放 VM 记在固体 \(H_{\mathrm{zu}}\)，同时 \(N_d/N_b\) 未带上同一批原子的气相焓。

本地摘录：`data/hamel_1999_extracted.txt`（PAGE 18–19）、`data/hamel_pdf_p43_extract.txt`、`data/hamel_pdf_p44_extract.txt`。分册：`docs/hamel_submodels/02_cell_balances_and_exchange.md`、`05_drying_and_daem_pyrolysis.md`。

---

## 6. 建议下一手（按顺序，仍是诊断）

先建环，再改生产入口。每步独立核 + pytest + `data/wirsum_alignment_*.json`。

### 6.1 热解气 Startwert（优先）

**问题：** 同一批 C/H/O 在 \(H_{\mathrm{zu,VM}}\)（−3.09 MW）里，不在 \(H_{\mathrm{gout}}\)。

**做法：**

1. 在 reduced 3-cell 上，BC + 轴向 holdup sweep 之后，跑 `residuals()`，读 `_vm_gas_source_cache`（若仍空，先弄清为何 cache 不写——`build_initialized_lu_reactor` 历史上不写）。
2. 把 §4.3 主气（CO/CO₂/H₂/H₂O/CH₄，以及 TAR/NH₃ 若在未知量内）按 dense 源项写入 `N_d`，使 \(H_{\mathrm{gout}}\) 出现对应生成焓。
3. A/B：VM 固体桶应被气体桶对消到感温/组成差量级；packed 能量残差应明显下降。
4. **命名：** 例如 `pyrolysis_gas_holdup_startwert`。不要叫 fastox，不要改 residual。
5. 对照已有 `seed_pyrolysis_gas_dense_holdup_for_nr_x0`：它目前只种 NH3/TAR，且依赖 `_vm_cache_valid`。扩展或并行走独立核，避免一次改完 WIP `vorab_x0.py`。

**禁止：** \(K_{\mathrm{VM}}>0\)；VM 列写入 `_solid_energy_outlet_rates`；把剩余 C 划回 Koks。

### 6.2 气体焓桶 −1.84 MW

冷气化剂 \(H_{\mathrm{zu,gas}}\approx-2.79\) MW vs \(H_{\mathrm{gout}}\approx-1.16\) MW。热解气进 holdup 后应重测：可能大部分被 6.1 吃掉，剩下才是 T 或组成 Startwert 问题。**先做 6.1 再拆。**

### 6.3 \(Q_W\approx0.69\) MW

`hamel_eq27_explicit` 显式壁面热。不是 bug。不要为了对齐残差而关热损。

### 6.4 bed1 接力后仍有 ~−2.7 MW

桶上 `other_solid_in` 与 `char_H_zu_minus_H_sout` 正负大数对消（上行焓 vs 当地 T 出流）。先确认 6.1 是否改变 bed0 出口气/固焓再决定是否值得单独拆 T。

### 6.5 接到 `generate_initial_x0` 的时机

诊断核稳定（R12 + holdup sweep + 热解气）后再考虑 **opt-in** 接到 fixture。默认 fixture 先保持「裸 x0」，否则 P0/P5 JSON 不可复现。

### 6.6 不要现在做的事

- 150-iter `solve_global_nr` / Phase1 HTW ladder
- 恢复 `extent_limiters_enabled`
- 工程 fastox 当 thesis、线搜索投影当 Startwert
- 用实验 CO 或 Fig.7 判断本链 x0 是否「对」
- 改 Eq.(2.7) 公式本身
- 给 VM 加输运系数或 `source−M` 锚

---

## 7. 文件地图

### 本链独立核（优先读）

| 路径 | 阶段 | 作用 |
|---|---|---|
| `src/core/solid_balance_reduced.py` | P1 | 流 vs 库存约化问题 |
| `src/solvers/wirsum_1997_strict.py` | P2 | Eq.2.21–2.23；`apply_wirsum_1997_strict_config` |
| `src/physics/wirsum_1997_freeboard.py` | P3 | `wirsum_1997_exact` vs `engineering_freeboard` |
| `src/solvers/wirsum_1997_outer_diagnostics.py` | P4 | outer 场增量 / \(\rho\)，不加松弛 |
| `src/solvers/wirsum_1997_lu_reentry.py` | P5 | 八项清单 |
| `src/thermal/hegermann_yield.py` | P5 | Tabelle 4.3 Koks/Teer |
| `src/core/pyrolysis_carbon_split.py` | P5 | Hegermann 后再 Gibbs |
| `src/solvers/wirsum_1997_elemental_decomp.py` | P5 | 源项元素 / R12 对消 |
| `src/solvers/wirsum_1997_startwert_r12.py` | P6 | 同相 R12 Startwert |
| `src/solvers/wirsum_1997_residual_decomp.py` | P6 | packed F 按气/固/能切片 |
| `src/solvers/wirsum_1997_energy_decomp.py` | P6 | Eq.2.7 分项 |
| `src/solvers/wirsum_1997_energy_leftover.py` | P6 | 剩余桶（必须加总到 live \(E\)） |
| `src/solvers/wirsum_1997_holdup_startwert.py` | P6 | \(M=m_{\mathrm{in}}/K\)；轴向接力 |
| `src/solvers/wirsum_1997_bed_freeboard.py` | F | 关回流床–自由空间界面；Wirsum exact ≠ Hamel 扩展 |
| `src/physics/wirsum_1997_cyclone.py` | F2 | Lapple 切割粒径；不是 Hamel Bild 7.7 |
| `src/solvers/wirsum_1997_bed_freeboard_coord.py` | F2 | 内层气/T + 外层轨迹固体 |
| `src/solvers/wirsum_1997_reaction_coord.py` | F3 | 离开炭与内层 R1 同一套 r |
| `src/solvers/wirsum_1997_cyclone_cells.py` | F4 | Hamel 旋风 1–2 cell + 返料进残差 |
| `src/solvers/wirsum_1997_bed_fb_cyclone.py` | F5 | 床–FB 界面与旋风 cell 交接 |

### P0–P5 审计文档

| 路径 | 内容 |
|---|---|
| `docs/solid_balance_formulation_audit.md` | P1：Eq.(2-6) vs \(m=KM\) |
| `docs/wirsum_1997_strict_newton.md` | P2：merit / frozen J / `K=0` |
| `docs/wirsum_1997_freeboard.md` | P3：模型名与 Roof 3.59 |
| `docs/wirsum_1997_outer_diagnostics.md` | P4：inner ≠ outer |
| `docs/wirsum_1997_lu_reentry.md` | P5：清单顺序 |
| `docs/wirsum_1997_elemental_decomp.md` | P5：H/O 净口 vs R12 |

### 生产对照（P0.3 / P1.2 已改过；其余只读，勿整文件提交）

| 路径 | 注意 |
|---|---|
| `src/solvers/global_nr_iteration.py` | P0.3：snap 只在 pre-NR / accepted-step |
| `src/solvers/global_nr_solver.py` | P0.3：`snap_trace_o2_holdups` 是投影，不是 `F(x)` |
| `src/core/cell_balances.py` | P1.2：`holdup_k0_anchor`；能量 \(H_{\mathrm{in}}-H_{\mathrm{out}}-Q\) |
| `src/core/cell.py` / `src/core/reactor.py` | P1.2 / P2：`holdup_k0_anchor`、`nr_solver_mode` |
| `src/solvers/vorabrechnung/vorab_x0.py` | `generate_initial_x0`、`project_bed_holdup_fast_oxidation_closure`、`seed_pyrolysis_gas_dense_holdup_for_nr_x0` |
| `src/workflow/steps/init_precalc_step.py` | 生产 x0 后 fastox / seed 顺序 |
| `src/core/connectivity/solid_transport.py` | \(K\) 只给 char/ash；VM 轴向置零 |
| `src/kinetics/gas_reactions.py` | `rate_R12` |
| `tests/global_nr_fixtures.py` | reduced 3-cell 入口 |

### 证据 JSON

`data/wirsum_alignment_p0_baseline.json`  
`data/wirsum_1997_strict_ab.json`  
`data/outer_iteration_history.json`  
`data/wirsum_alignment_p5_lu_reentry.json`  
`data/wirsum_alignment_p5_elemental_decomp.json`  
`data/wirsum_alignment_p6_startwert_r12.json`  
`data/wirsum_alignment_p6_residual_decomp.json`  
`data/wirsum_alignment_p6_energy_decomp.json`  
`data/wirsum_alignment_p6_holdup_startwert.json`  
`data/wirsum_alignment_p6_energy_leftover.json`  
`data/wirsum_alignment_p6_holdup_sweep.json`  
`data/wirsum_alignment_stage_a_baselines.json`  
`data/wirsum_alignment_stage_c_solid_transport.json`  
`data/wirsum_alignment_stage_d_vorab_context.json`  
`data/wirsum_alignment_stage_e_stiffness_ladder.json`  
`data/wirsum_alignment_stage_f_bed_freeboard.json`  
`data/wirsum_alignment_stage_f2_bed_freeboard_coord.json`  
`data/wirsum_alignment_stage_f3_reaction_coord.json`  
`data/wirsum_alignment_stage_f_cyclone_cells.json`  
`data/wirsum_alignment_stage_f5_bed_fb_cyclone.json`

---

## 8. 给诊断 agent 的最小反馈环

症状（用户原话级）：**「Hamel 的 VM 应并入气相物质/能量衡算；bed0 在 holdup Startwert 后仍有 −3 MW 固体 VM 焓。」**

建议环（由快到慢）：

```bash
# 1. 现有桶必须仍加总
python -m pytest tests/test_wirsum_1997_energy_leftover.py tests/test_wirsum_1997_holdup_startwert.py -q

# 2. 看热解源项是否进了 R_gas / cache（自己写 20 行 probe，勿先改生产）
#    build_initialized_lu_reactor → BC → axial sweep → residuals
#    打印：sum(R_gas) 的 H2/CO/H2O/CH4、_vm_cache_valid、H_gout、H_zu_vm

# 3. 有环之后再写独立核 + A/B JSON
```

判定成功（仅诊断阶段）：

- 新桶分解仍 `sum_matches_residual`
- 热解原子从固体 VM 列转到气相 \(N\) 后，`vm_H_zu_minus_H_sout` 被气体项对消，而不是被 \(K_{\mathrm{VM}}\) 或 T clip 抹掉
- 不出现实验 CO / `y_CO_target` 作为判据

---

## 9. 与 2026-08 大交接的关系

`docs/handoff_nr_convergence_hamel_isomorph_2026-08.md` 的 CORE 栈（fastox、T 帽、fence、820 K 种子等）是 **工程收敛路径**。本链是 **数学问题重建**，刻意不接那些开关。

该文档 §5.61 已写过「bed0 焓洞是 VM zone=3 + `K_VM=0`」和「勿把 VM 补进 Eq.2.7」。本轮用 reduced x0 把同一洞量化到 −3.09 MW，并补上 Hamel §2.2.2 / §4.3 的物种口径：正确闭合是 **气相 \(N_r\)**，不是固体 VM 出流。
