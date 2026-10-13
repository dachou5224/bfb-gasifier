# 具名核速率延拓（parameter homotopy）

本路径在明确的具名约化残差上使用速率参数延拓和自选初值，以寻找完整速率端点；延拓不属于 Hamel/Wirsum 原文步骤。端点按未投影的目标残差独立验收，求解器停止原因、参数推进是否逐档收敛以及原文之外的修正措施分别报告。

本文记录具名 \(F\) 谱系里，对 **R12 / R5 / R6 / R10** 这类快速耗氧氧化接入完整速率时用的 **参数延拓**。它是求解器技术，**不是** Hamel 方程，也不是生产默认。审查：`docs/audit_reports/rate_continuation_review_2026-09-15.md`。

对照：

- 方程 SoT：`docs/hamel_submodels/01`（Newton 框架）、`06`（R1–R11 速率式）
- Wirsum 阻尼 Newton：`docs/wirsum_1997_nr_solver_notes.md`、`docs/wirsum_1997_strict_newton.md`
- 一般延拓家族：`docs/siegen_cell_model_lineage.md` §13.3
- Startwert / 工程桥：`docs/vorab_init_abgleich_vs_outer_preinner.md`、`docs/hamel_engineering_bridge_notes.md`
- 刚性诊断：`docs/audit_reports/secondary_diagnosis_2026-09-14.md` 阶段 E

---

## 1. 定位：延拓改的是路径，不是模型

目标问题始终是完整速率的代数根

\[
F(x;q)=0,\qquad s=1.
\]

延拓构造一族更容易的问题 \(F(x;s)=0\)，\(s\) 从 \(0\) 走到 \(1\)。`solver_stop_reason` / `physical_convergence` 只描述 Newton 停机（merit 目标 \(10^{-10}\)，报告阈 \(10^{-6}\)）；二者可以并存，`line_search` 本身不是“反应太快”的证明。完整速率端点另由 `endpoint_certified` 验收：\(s=1\)、未投影、分块残差、定义域。中间档的收敛只证明「这一档的 \(F(\cdot;s)\) 有根」，不得写成炉温根、不得切进生产 `Reactor`。

不改 `rate_R12` / `rate_R5_*` / `rate_R6` / `rate_R10` 的阶数、温度幂或前因子。不恢复 fastox、O₂ extent limiter。

---

## 2. 原因：同相燃料∩氧是硬源

Hamel 第 5 章把均相氧化写成浓度幂律，例如

- R12（论文 R6）：\(r_{12}=k C_{\mathrm{H_2}}^{3/2}C_{\mathrm{O_2}}\)（Eq. 5.38–5.39）。氧一阶：\(C_{\mathrm{O_2}}=0\) 且 \(C_{\mathrm{H_2}}>0\) 时 \(r=0\)，但右导数 \(\partial r/\partial C_{\mathrm{O_2}}|_{0^+}=k C_{\mathrm{H_2}}^{3/2}\neq 0\)。
- R5 气泡：\(r_5\propto C_{\mathrm{CO}}C_{\mathrm{O_2}}^{0.3}C_{\mathrm{H_2O}}^{0.5}\)（Eq. 5.36）；悬浮相 Table 5.3 含氧 \(0.5\) 次。\(0<\beta<1\) 时分数幂在零处导数发散，**不能**套用 R12 的有限右导数。CO 一阶在 \(C_{\mathrm{CO}}=0\) 时可保留有限右导数。
- R6（论文 R7）：\(r_6\propto C_{\mathrm{CH_4}}C_{\mathrm{O_2}}\)（Eq. 5.40–5.41）
- R10：\(r_{10}=k_{10}e^{-E/RT}\,T\,P^{0.3}C_{\mathrm{tar}}^{0.5}C_{\mathrm{O_2}}\)（Eq. 5.59）。氧一阶保留右导数；\(C_{\mathrm{tar}}^{0.5}\) 在 0 处发散。

\(C_{\mathrm{O_2}}=0\) 时这些速率归零。稳态限氧来自物种衡算，不是 extent limiter。若某相 holdup 里燃料与 \(\mathrm{O_2}\) 同时为正，源项 \(\xi=rV\) 可达 \(10^5\)–\(10^{13}\,\mathrm{mol/s}\)，相对格内摩尔流 \(10^{-2}\) 是硬源：说明该初值远离守恒解，不等于原模型允许稳态无限耗氧。

Wirsum §2.3.3 的阻尼 Newton 只接受 \(\|F\|_2\) 严格下降（Eq. 2.23）。在硬源点上：

1. \(\|F\|_2\) 被源项行撑到天文数字；
2. 线性化 \(J\Delta x=-F\) 给出的步会把组分打到负值，或点燃下一条氧化（R10 产 \(\mathrm{H_2}/\mathrm{CO}\) 后再碰到 R12/R5）；
3. 线搜索所有 \(\alpha=2^{-j}\) 都被拒，求解器报 `line_search`，不是「式子写错」。

因此需要一条 **从已有数学根出发、逐步打开这一条反应** 的路径，而不是在完整速率上盲目换初值。

---

## 3. 与 Hamel / Wirsum 原文的关系

延拓 **没有** 对应的 Hamel 方程号。论文给的是模型与 Newton 框架；本核把「远初值」拆成可跟踪的参数路径。

| 原文 | 论文实际写了什么 | 本核怎么用 | 不是什么 |
|------|------------------|------------|----------|
| Hamel §2.1 Startwertwahl | 配置阶段给 **全部** 守恒变量初值；可沿用上一轮解（Bild 2.2 Konfiguration） | 同相重叠 \(\xi=\min(N_{\mathrm{fuel}},N_{\mathrm{O_2}}/\nu)\) 只改 \(x_0\)，**不改** \(F\) | 不是把 \(r\) 乘 \(s\)；不是 fastox |
| Hamel §2.3 / Eq. 2.9 | Wozny NR + Wirsum 阻尼；块三对角 | 每一档 \(s\) 上跑 `wirsum_1997_strict_newton` | 延拓参数 \(s\) 不是 Eq. 2.9 的一部分 |
| Wirsum §2.3.3 Eq. 2.21–2.23 | \(J\Delta x=-F\)，仅当 \(\|F\|_2\) 下降才接受 | 每一档仍用该阻尼 Newton；补对角 / 复合 \(P(x)\) 须另行声明 | 原文没有 \(s\) 网格 |
| Hamel Ch.5 速率式 | 完整 \(k(T)\)、浓度阶、计量 | \(s<1\) 只缩放 **源项增量**；\(s=1\) 回到原式 | 不得把 \(s<1\) 的根当成物理模型 |
| 生产 `r5_scale` / `nr_reaction_continuation_*` | 工程旋钮 | 生产默认应保持 1；具名核的 \(s\) **不写入** `ReactorConfig` | 两套 scale 不得混比 |

Hamel 处理「初值离根太远」的原文手段是 **Startwertwahl + 阻尼 Newton**，不是参数同伦。本核两者都用：Startwert 把 \(x\) 搬到快流形附近；延拓是原文之外的求解辅助。内层并不总是「纯 Wirsum \(J=F'\)」：亏秩时给微量组分补对角是 **修正 Newton**，不是等价缩放；E4 起默认残差是复合 \(H(x)=F(x;P(x))\)（`composite_P`），与冻结 q 的 `frozen_q` 分名，不得用 `hamel_outer_refresh=False` 暗示内层 q 全冻。

数学上这是标准 **convex homotopy**（见谱系文档 §13.3），不是伪瞬态 \(\Psi tc\)，也不是伪弧长。当前实现是 **离散 \(s\) 网格上的自然参数延拓**：上一档的 \(x^*(s_k)\) 作为 \(s_{k+1}\) 的初值。解支若回折，这种网格会失败；具名核迄今走的是单调打开一条反应，尚未需要弧长。

---

## 4. 算法

### 4.1 残差族

已有根 \(F_A(x_A)=0\)（例如 F_E0+R12）。要打开下一条反应 \(B\)（例如 R5）。令 \(F_{A\cup B}\) 为两条反应都开、**完整速率** 的残差。定义

\[
F(x;s)=F_A(x)+s\bigl(F_{A\cup B}(x)-F_A(x)\bigr),\qquad s\in[0,1].
\]

性质：

- \(s=0\)：\(F(\cdot;0)=F_A\)，旧根仍是根；
- \(s=1\)：\(F(\cdot;1)=F_{A\cup B}\)，目标模型；
- \(F_{A\cup B}-F_A\) 在**固定 \(x\)** 下只出现在 \(B\) 的物种源项。当前焓流能量方程已含生成焓；固定 \(x\) 时能量行不应凭空再加一份反应热。组成随**不同根**变化后，焓流自然改变。二者不要混成「源项再加一份反应热」。

实现：`residual_*_plus_scaled_*`。R12 第一刀没有「旧反应 \(A\)」时，退化为 `coupled_residual_named(..., rate_scale=s)`，即所有已启用反应的源项同乘 \(s\)。

### 4.2 \(s\) 网格

硬源量级要求从极小 \(s\) 起步。典型无 Startwert 网格：

```
0, 1e-16, 1e-12, 1e-8, 1e-6, 1e-4, 1e-2, 1.0
```

做了同相 Startwert 之后，\(s=0\) 会按旧模型把刚消掉的重叠 **再混回去**。这时从有限速率起步：

```
1e-4, 1e-2, 1e-1, 1.0
```

R10 的网格更稀、且 **默认不做** holdup \(\mathrm{O_2}=0\) Startwert（见 §6.4）。档点是数值参数，不依据实验 CO 或目标温度选取。

两种推进模式（`advance_mode`）：

- `root_correcting`：该档 \(\|F\|_2<10^{-6}\) 且在定义域才升 \(s\)；失败则插入中点、不跳档。R5 有 Startwert 时默认此模式。
- `heuristic_search`：只承诺提供候选初值；有接受的线搜索步即可推进。R10 默认此模式，不得称为逐根跟踪。

### 4.3 内层 Newton

每一档调用 `wirsum_1997_strict_newton`：

1. 组装 \(F(x;s)\)、\(J=\partial F/\partial x\)（旧块 + \(s\times\) 新反应解析块；\(T\) 列常对 \(F(x;P(x))\) 差商，含 DAEM）；
2. 线性步 `named_linear_step`：微量组分尺度 + 行列均衡。诊断区分原 \(J\) 与补对角 \(J\)（`used_trace_diag`、`linear_residual_raw_l2` / `linear_residual_modified_l2`、方向导数 \(\phi'\)）。补对角 `_TRACE_DIAG=1` 是修正 Newton，不是解析导数，也不是等价缩放；
3. \(\alpha=1,1/2,1/4,\ldots\) 到 `lambda_min=1/4096`。失败不是“所有正 \(\alpha\) 均无效”。域外拒绝、非有限值、无下降与线性步质量分开记；
4. `classify_newton_stop` 给出 `solver_stop_reason`。`small_step` ≠ 物理根；`physical_convergence` 要求完整速率且 \(\|F\|_2<10^{-6}\)，仍不是分块端点认证。

试探越界（\(N<0\) 或 \(T\le 0\)）缩短 \(\alpha\)，不改残差定义。细炭 \(M=0\) 时只修正搜索方向（`bound_feasible_direction`），**残差本身不 clip**。投影会把所评价的函数改成 \(\widetilde G(x,s)=G(\Pi_+x,s)\)，且返回的 \(J\) 不再是它的导数。

### 4.4 Startwert（只改 \(x_0\)）

同相计量投影，Hamel §2.1 的 *Startwertwahl* 在本核的最小实现：

| 反应 | \(\xi\) | 产物 |
|------|---------|------|
| R12 | \(\min(N_{\mathrm{H_2}},2N_{\mathrm{O_2}})\) | \(\mathrm{H_2O}\) |
| R5 | \(\min(N_{\mathrm{CO}},2N_{\mathrm{O_2}})\) | \(\mathrm{CO_2}\) |
| R6 | \(\min(N_{\mathrm{CH_4}},N_{\mathrm{O_2}}/1.5)\) | \(\mathrm{CO}+2\mathrm{H_2O}\) |
| R10 | 共享焦油池，\(\min(N_{\mathrm{tar}},N_{\mathrm{O_2}}/\nu_{\mathrm{O_2}})\) | \(\mathrm{CO}+\mathrm{H_2}\) |
| R11 | 共享焦油池，\(\min(N_{\mathrm{tar}},N_{\mathrm{H_2O}}/\nu_{\mathrm{H_2O}})\) | \(\mathrm{CO}+\mathrm{H_2}+\mathrm{CH_4}\) |

规则：

- 按相独立（R10 的焦油是全格槽，先密相后气泡）；
- 同步 FB / 旋风气接口拷贝，避免接口行被旧组成拉回去；
- **不进** \(F(x)\)，不是生产 `project_bed_holdup_fast_oxidation_closure`；
- 流形上必须保留 **右导数**（例如 \(C_{\mathrm{CO}}=0\) 但 \(C_{\mathrm{O_2}}>0\) 时 \(\partial r_5/\partial C_{\mathrm{CO}}\neq 0\)），否则雅可比亏秩。

Startwert 只消 **holdup 重叠**。入口持续供氧时，\(N_{\mathrm{O_2}}=0\) 且产物仍在，Newton 一旦把微量氧加回来，下一条氧化会爆。R10 因此默认 **不** 用 \(\mathrm{O_2}=0\) 当 \(s=1\) 初值。R11 焦油耗尽后入口水蒸气不会再生成 tar，默认可用 \(\mathrm{tar}\cap\mathrm{H_2O}\) Startwert。

### 4.5 一次接入一条反应

顺序固定，与具名谱系一致（仓库编号；括号为 Hamel 论文号）：

\[
\mathrm{R12\,(R6)}\ \to\ \mathrm{R5}\ \to\ \mathrm{R8}\ \to\ \mathrm{R6\,(R7)}\ \to\ \mathrm{R7\,(R9)}\ \to\ \mathrm{R2\text{–}R4}\ \to\ \mathrm{F_{E3}\text{–}E6}\ \to\ \mathrm{R10}\ \to\ \mathrm{R11}\ \to\ \mathrm{F_{E7}}\ \to\ \mathrm{F_{E8}}.
\]

F_E7 不是又一条床内速率式，而是取消登记项 `freeboard_reaction`：在已有 43 维核上给自由板混合气加均相 CSTR 源。再下一步是 F_E8（`K_eject_in_this_P`）。旋风 K 与生产 `hamel_outer_refresh` 永不取消。

每一刀的 \(A\) 是上一刀的数学根。禁止把生产 `r5_scale` 与具名 \(s\) 混在一张表里比 \(T\) 或 CO。

---

## 5. 和仓库里其它 “scale / continuation” 的区别

| 名字 | 文件 / 旋钮 | 作用 | 物理根？ |
|------|-------------|------|----------|
| 具名 \(s\) | `wirsum_1997_named_f_r*.py` 的 `r*_scale` | 打开下一条反应的同伦参数 | 仅 \(s=1\) |
| `coupled_residual_named(..., rate_scale=)` | `wirsum_1997_named_f.py` | 已启用反应源项同乘；R12 第一刀用 | 仅 `=1` |
| 生产 `ReactorConfig.r5_scale` 等 | `src/core/reactor.py` | 工程缩放，留下就改模型 | 默认必须是 1 |
| `nr_reaction_continuation_stages_thesis` | 生产 NR | 按高度/阶段降反应强度的工程桥 | 默认 `(1.0,)`，不开放 |
| Check1 continuation | `nr_outer_preinner_on_check1_continuation_thesis` | 外层 Check1 失败后是否继续 NR | 与速率同伦无关 |

---

## 6. 对比算例（具名核，不经 `Reactor`）

数字来自 `tests/test_wirsum_1997_named_f_r*.py` 锁住的行为，以及 R10 接入时的记录。\(T\) 是 **这一格数学根**，不是 HTW 炉温。

### 6.1 R12：延拓足以走到 \(\mathrm{H_2}\) 耗尽流形

| 状态 | \(\|F\|_2\) | 备注 |
|------|-------------|------|
| 默认 \(x_0\)（同相 \(\mathrm{H_2}\cap\mathrm{O_2}\)） | \(>10^{12}\) | 硬源，不是公式写错 |
| 同相 R12 Startwert 之后 | \(<2\times10^{5}\) | 重叠消掉，能量行仍在 |
| 从 F_E0 根延拓到 \(s=1\) | \(<10^{-6}\) | 秩 40；\(N_{\mathrm{H_2}}<10^{-8}\)；\(T<1000\,\mathrm{K}\) |

测试：`test_r12_solve_does_not_change_rate_or_claim_furnace`。入口：`solve_named_r12_only`。

### 6.2 R5：只有 Startwert + 延拓才到根

在 F_E0+R12 根上，气泡 \(\mathrm{CO}\cap\mathrm{O_2}\) 仍是硬源（\(|dF_{\mathrm{CO},b}|\) 比密相大约 \(10^{2}\)）。

| 路径 | \(s\) 网格 | `physical_convergence` | 观察 |
|------|------------|------------------------|------|
| 无 Startwert，从 \(s=0\) 延拓 | \(0,\ldots,1\) | **否** | \(s=10^{-6}\) 仍可局部收敛；\(s=1\) 接不住 |
| 同相 R5 Startwert 后从 \(s=10^{-4}\) | \(10^{-4}\to1\) | **是** | \(T\) 相对 R12 根升高 \(>200\,\mathrm{K}\)，落在 \(1000\)–\(1200\,\mathrm{K}\)；仍声明非炉温 |

恒等式：`r5_scale=0` 的残差与 F_E0+R12 逐元相同（`test_r5_scale_zero_matches_r12_residual`）。

### 6.3 R6：与 R5 同构

无 Startwert：\(s=10^{-8}\) 可局部收敛，\(s=1\) 失败。  
有 Startwert（再补一刀 R5，清 R6 产的 \(\mathrm{CO}\cap\mathrm{O_2}\)）：\(s=1\) 收敛，\(T\) 从约 \(1098\,\mathrm{K}\) 升到约 \(1314\,\mathrm{K}\)，两相 \(\mathrm{CH_4}<10^{-4}\,\mathrm{mol/s}\)。

### 6.4 R10：错误 Startwert 比没有更糟

F_E6 根上 \(\xi_{\mathrm{R10}}\sim10^{5}\,\mathrm{mol/s}\)。holdup 上把两相 \(\mathrm{O_2}\) 消到 0、留下焦油与 \(\mathrm{H_2}/\mathrm{CO}\) 之后：

- 气泡入口仍供氧，\(F_{\mathrm{O_2},b}\sim 0.4\,\mathrm{mol/s}\)；
- Newton 加回微量氧 → R12/R5 爆炸；
- 能量行约 \(+4.6\,\mathrm{kW}\)，单独抬 \(T\) 约 \(+30\,\mathrm{K}\) 也救不了物种步。

从 F_E6 根做延拓（不走 \(\mathrm{O_2}=0\)）。默认 `advance_mode=heuristic_search`，`x_nonneg_projection=False`。去掉残差 clip 后，旧表中 \(s=1\) 的 \(\|F\|_2\sim 10^{-9}\) 数字来自投影路径，不得再锁成正确行为。端点须用 `certify_named_full_rate_endpoint` 对未投影 \(F_{\mathrm{E6+R10}}\) 独立验收。

入口：`solve_named_e6_then_r10`。测试：`tests/test_wirsum_1997_named_f_r10.py`、`tests/test_wirsum_1997_named_f_continuation.py`。

### 6.5 R11：焦油转化，不升 F_E7

在 F_E6+R10 根上打开 Table 5.4：气泡 Serio 一级 \(k(T)C_{\mathrm{tar}}\)，悬浮 Corella \(k(T)\rho_{\mathrm{cat}}C_{\mathrm{tar}}\)，\(\rho_{\mathrm{cat}}=M_{\mathrm{ash}}/V_d\)。不改 `rate_R11_*`。未知量仍是 43；`freeboard_reaction` 仍省略。\(C_{\mathrm{tar}}=0\) 时保留一级右导数。默认 Startwert + `root_correcting`。

入口：`solve_named_r10_then_r11`。测试：`tests/test_wirsum_1997_named_f_r11.py`。

### 6.6 F_E7：自由板均相 CSTR，不升 F_E8

在 F_E6+R10+R11 根上打开 `fb_gas_source`：体积 `A H_{\mathrm{fb}}`，气泡相均相 R5/R6/R7/R8/R12，`ρ_ash=0` 故 R8 关闭。浓度与 Arrhenius 用状态 \(T=x[-1]\)（`q_at_T`），不是冻结 geo 1100 K。不含床内焦油 R10/R11，也不把生产 `_freeboard_reaction_step` / 轨迹积分搬进具名核。

未知量仍是 43。残差族 \(F_{\mathrm{E6+R10+R11}}+s(F_{\mathrm{E7}}-F_{\mathrm{E6+R10+R11}})\)。两相在床内可以各自无氧，混合成 \(N_{\mathrm{fb}}=N_d+N_b\) 后仍可能出现燃料 ∩ O₂（R11 根上微量氧即可让 R12 源到 \(10^4\,\mathrm{mol/s}\)）。默认只对 `N_fb` 做 R12→R5→R6→R5 Startwert，再同步旋风气拷贝，然后直接 \(s=1\)：小 \(s\) 会把混合约束里的氧加回 `N_fb`，右导数再次硬源。Startwert 的计量焓差大约 \(10^{-4}\,\mathrm{W}\)；混合 \(\|F\|_2\) 把瓦特和 mol/s 加在一起，不得拿来代替 `endpoint_certified`。默认 `root_correcting`，\(s=1\) 以分块认证为准。`K_eject` 仍省略。

入口：`solve_named_e7_from_r11`。测试：`tests/test_wirsum_1997_named_f_e7.py`。

### 6.7 F_E8：抛射 \(K\) 进 \(P(x)\)，不更新旋风 K

在 F_E7 根上取消 `K_eject_in_this_P`。`outer_P_named` 用冻结 `geo_fb.rho_g` 相对参考温度 \(T_{\mathrm{ref}}=1100\,\mathrm{K}\) 做理想气体缩放，再算 \(K_{\mathrm{eject}}=m(M=1)/M_{\mathrm{geom}}\) 与平均抛射命运 \(w_{\mathrm{leave}}\)。Roof/Wake 流量本身不含 \(T\)；命运随气密可能翻转。平均初速下，本几何 800 K 粗档是轨迹**超时未判定**，不得记成逃逸。旋风 `k` 仍冻。这不是生产 `hamel_outer_refresh`。

E7 根约 \(1420\,\mathrm{K}\)：细档仍离开、粗档仍回落，\(K\) 与冻结值相同，故 \(F_{\mathrm{E8}}=F_{\mathrm{E7}}\)，根不移动。未知量仍是 43。

入口：`solve_named_e8_from_e7`。测试：`tests/test_wirsum_1997_named_f_e8.py`。P0 后完整链写入 `data/named_e8_endpoint_2026-09-17.json`；短检把其中的 \(x\) 当初值并现场 `certify_named_full_rate_endpoint`，不得用文件里的布尔字段当通过证据。重算：`data/named_e8_recompute_2026-09-17.py`。

---

## 7. 验收与禁止事项

**算通过（完整速率端点）：**

- \(s=1\)，所有已声明启用反应恢复完整；状态在定义域 \(x[:-1]\ge 0\)、\(T>0\)；
- 不经延拓差值与 clip 重新装配目标 \(F\)；`endpoint_certified`；
- 分块：气体 mol/s、固体 kg/s、焦油 mol/s、能量 W（能量绝对误差 \(<1\,\mathrm{W}\)）；并核验 q 一致性（`q_consistency_l2`）；
- `solver_stop_reason`、`advance_mode`、`branch_tracked`、`coupling_mode` 分别报告；
- `rate_formula_changed is False` 只表示未改 Arrhenius 公式，**不能**代替源项/方程恒等性；
- `not_a_furnace_convergence_claim is True`。抵达不同根时报告分支，不用「温度偏离 1100 K」作正确性判据；
- 生产 `hamel_outer_refresh` 仍关；`thesis_vorab_sources_single_shot` 仍为生产默认。

**不算通过：**

- 中间 \(s<1\) 的残差小；
- 把 `line_search` / `small_step` / `max_iter` 说成物理根；也不把所有 `line_search` 一律排除在诊断之外；
- 用补 \(Q_W\)、残差 clip 或补对角冒充解析导数，把根按到 \(1100\,\mathrm{K}\)；
- 改速率阶数、恢复 fastox、开放生产 continuation 旋钮；
- 混用 `composite_P` 与 `frozen_q` 却当成同一固定点。

原文 raw L2 可保留为算法对照，但它混合不同单位；`\|F\|<10^{-6}` 和满秩也不是完整物理认证。

---

## 8. 代码入口

| 反应 | 残差 / 求解 | Startwert |
|------|-------------|-----------|
| 共享推进 / 认证 | `wirsum_1997_named_f_continuation.py` | — |
| R12 | `wirsum_1997_named_f_r12.py` | `wirsum_1997_startwert_r12.py` |
| R5 | `wirsum_1997_named_f_r12_r5.py` | `wirsum_1997_startwert_r5.py` |
| R8 | `wirsum_1997_named_f_r12_r5_r8.py` | （无重叠硬源，只延拓） |
| R6 | `wirsum_1997_named_f_r12_r5_r8_r6.py` | `wirsum_1997_startwert_r6.py` |
| R7 | `wirsum_1997_named_f_r12_r5_r8_r6_r7.py` | （弱源，只延拓） |
| R2–R4 | `wirsum_1997_named_f_r234.py` | （固体块，只延拓） |
| R10 | `wirsum_1997_named_f_r10.py` | `wirsum_1997_startwert_r10.py`（可选，默认不用） |
| R11 | `wirsum_1997_named_f_r11.py` | `wirsum_1997_startwert_r11.py`（默认用） |
| F_E7 | `wirsum_1997_named_f_e7.py` | `wirsum_1997_startwert_e7.py`（只改 `N_fb`） |
| F_E8 | `wirsum_1997_named_f_e8.py` | （E7 根上命运不翻转，无 Startwert） |

登记：`src/solvers/wirsum_1997_named_f_registry.py`。内层 Newton：`src/solvers/wirsum_1997_strict.py`。停机分类：`classify_newton_stop`。端点认证：`certify_named_full_rate_endpoint`。E4 起耦合：`composite_P`（默认）/ `frozen_q`。
