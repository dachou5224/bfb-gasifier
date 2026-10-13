# 交接：具名 \(F\) 反应接入到 F_E8（2026-09-16）

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读：

- 本文件（进展、未完成项、未启动原因、待决问题）
- 延拓算法：`docs/wirsum_1997_named_f_rate_continuation.md`
- 登记：`src/solvers/wirsum_1997_named_f_registry.py`
- Q 裁定：`docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md`（尤其 Q4、Q6、Q7）
- 更早长交接：`docs/handoff_wirsum_alignment_2026-09-14.md`

**分支：** `feat/nr-hamel-wirsum-core-startwert`  
**日期：** 2026-09-16  
**Git：** 具名 \(F\) 谱系（`wirsum_1997_named_f*.py`、对应 tests、延拓文档）**整包未 commit**。已提交的 HEAD 仍停在 holdup Startwert 轴向接力（`098c8cf`）。`wirsum_1997_startwert_r12.py` 更早已进仓。  
**状态：** 具名核上「一次开一条反应 / 取消一项登记省略」已走到 F_E8。生产 `Reactor` 默认未切。未跑全炉 LU。未宣称炉温根，未宣称完整 \(q=P(x)\)。  
**语言：** 简体中文。

目标仍是重建 **Wirsum/Hamel 原始代数问题** \(F(x;q)=0\)，不是调参让生产 LU 收敛。本核不经 `Reactor`。两条路径不得混比：reduced 3-cell vs `thesis_mode`。

---

## 1. 一句话

反应接入已经走完；F_E8 是登记谱系最后一项会取消的省略。剩下的不是下一条 `rate_*`，而是外层 \(P(x)\)、旋风 \(K\)、返料腿、自由板反应集扩维，以及要不要把具名核接到生产路径。这些都还没启动，原因见 §4。下一手必须先选定问题，不要默认续写 F_E9。

---

## 2. 硬约束（接手不得违反）

| 项 | 规定 |
|---|---|
| 单位 | 全程 SI。Arrhenius 只走 `src/kinetics/arrhenius.py` |
| 速率式 | 不改 `rate_R1`–`rate_R12` 的阶数、温度幂、前因子。不恢复 fastox / O₂ extent limiter |
| 生产默认 | 不切 `hamel_outer_refresh`。`thesis_vorab_sources_single_shot` 保持生产默认 |
| 炉温 | 不把根按到 \(1100\,\mathrm{K}\)。`GEO_FB_RHO_G_T_REF_K=1100` 只是冻结 `geo_fb.rho_g` 的参考温度 |
| \(s<1\) | 不是物理模型。完整速率端点用 `endpoint_certified`，不用混合 \(\\|F\\|_2\) 当门 |
| 路径 | 不混 reduced 3-cell 与 `thesis_mode`；不拿生产 `r5_scale` 与具名 \(s\) 比 \(T\) 或 CO |
| 验收 | 气体/固体/焦油 L2 \(<10^{-6}\)（对应单位 mol/s 或 kg/s），能量 \(\lvert E\rvert<1\,\mathrm{W}\) |
| 提交 | 用户未要求则不 commit |

相关技能：`.cursor/skills/bfb-diagnosing-bugs`、`bfb-hamel-parity-audit`、`bfb-numerical-guard`、`bfb-scientific-audit`。

---

## 3. 已完成

### 3.1 登记阶段 F_E0–F_E8

未知量：E0–E2 = 40；E3–E5 = 41（+灰）；E6–E8 = 43（+TAR1/TAR2）。温度始终在最后。

| 阶段 | 取消的省略 | 入口 | 测试 |
|------|------------|------|------|
| F_E0 | —（仅 R1 在物种块） | `coupled_residual_named` / `wirsum_1997_coupled_f_energy.py` | `tests/test_wirsum_1997_named_f_stages.py` |
| F_E1 | R2–R4 | `wirsum_1997_named_f_r234.py` | `test_wirsum_1997_named_f_r234.py` |
| F_E2 | 均相 R5/R6/R7/R8/R12 | 见下表分刀 | 各 `r12_r5_*.py` |
| F_E3 | `ash_inventory` | `wirsum_1997_named_f_e3.py` | `test_wirsum_1997_named_f_e3.py` |
| F_E4 | `DAEM_X_from_T` | `wirsum_1997_named_f_e4.py` | `test_wirsum_1997_named_f_e4.py` |
| F_E5 | `energy_wall_heat`（\(Q_W\) 进 q，默认 0） | `wirsum_1997_named_f_e5.py` | `test_wirsum_1997_named_f_e5.py` |
| F_E6 | `tar_transport`（焦油出口升为 \(x\)，尚未开 R10） | `wirsum_1997_named_f_e6.py` | `test_wirsum_1997_named_f_e6.py` |
| F_E7 | `freeboard_reaction` | `solve_named_e7_from_r11` | `test_wirsum_1997_named_f_e7.py` |
| F_E8 | `K_eject_in_this_P` | `solve_named_e8_from_e7` | `test_wirsum_1997_named_f_e8.py` |

E4 起默认耦合 `composite_P`：\(H(x)=F(x;P(x))\)。对照接口 `frozen_q` 仍在，不得与 `composite_P` 当成同一固定点。

### 3.2 反应接入顺序（仓库编号；括号为 Hamel 论文号）

\[
\mathrm{R1}\ \text{（F\_E0）}\ \to\ \mathrm{R12\,(R6)}\ \to\ \mathrm{R5}\ \to\ \mathrm{R8}\ \to\ \mathrm{R6\,(R7)}\ \to\ \mathrm{R7\,(R9)}\ \to\ \mathrm{R2\text{–}R4}\ \to\ \mathrm{E3\text{–}E6}\ \to\ \mathrm{R10}\ \to\ \mathrm{R11}\ \to\ \mathrm{E7}\ \to\ \mathrm{E8}.
\]

| 反应 | 求解入口 | Startwert |
|------|----------|-----------|
| R12 | `solve_named_r12_only` | `wirsum_1997_startwert_r12.py` |
| R5 | `wirsum_1997_named_f_r12_r5.py` | `wirsum_1997_startwert_r5.py` |
| R8 | `wirsum_1997_named_f_r12_r5_r8.py` | 无重叠硬源，只延拓 |
| R6 | `wirsum_1997_named_f_r12_r5_r8_r6.py` | `wirsum_1997_startwert_r6.py` |
| R7 | `wirsum_1997_named_f_r12_r5_r8_r6_r7.py` | 弱源，只延拓 |
| R2–R4 | `wirsum_1997_named_f_r234.py` | 固体块，只延拓 |
| R10 | `solve_named_e6_then_r10` | `startwert_r10.py` 可选，**默认不用**（\(\mathrm{O_2}=0\) 比没有更糟） |
| R11 | `solve_named_r10_then_r11` | `wirsum_1997_startwert_r11.py`（默认用） |
| F_E7 | `solve_named_e7_from_r11` | `wirsum_1997_startwert_e7.py`（只改 `N_fb`，然后直接 \(s=1\)） |
| F_E8 | `solve_named_e8_from_e7` | 本根上命运不翻转，无 Startwert |

共享推进 / 认证：`wirsum_1997_named_f_continuation.py`（`certify_named_full_rate_endpoint`）。  
E7/E8 停机后若 Newton 把已认证的根走坏，回退到上一阶段 \(x\)。混合 \(\\|F\\|_2\) 把 W 和 mol/s 加在一起，不得当门。

### 3.3 F_E8 当前事实

- `outer_P_named` 在 `stage>=8` 用冻结 `geo_fb.rho_g` 相对 \(T_{\mathrm{ref}}=1100\,\mathrm{K}\) 做理想气体缩放，再算 \(K_{\mathrm{eject}}=m(M=1)/M_{\mathrm{geom}}\) 与平均抛射命运 \(w_{\mathrm{leave}}\)。
- Roof/Wake 流量本身不含 \(T\)。旋风 `k` 不更新。
- E7 根约 \(1420\,\mathrm{K}\)：细档离开、粗档回落，\(K\) 与冻结值相同，故 \(F_{E8}=F_{E7}\)，根不移动。
- 廉价探针：\(T=800\,\mathrm{K}\) 时粗档也可改为离开（\(w=[1,1]\)），\(K\) 仍不变。这不是本根上的翻转。
- 测试：`PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_wirsum_1997_named_f_e8.py --tb=short` → **4 passed**（含 E7 链，约 8 min，2026-09-16）。

### 3.4 当前核刻意保持的窄边界

这些是已实现模型的边界，不是漏掉的 `rate_*`：

- 自由板：单槽均相 CSTR，体积 \(A H_{\mathrm{fb}}\)，气泡动力学，\(\rho_{\mathrm{ash}}=0\)（R8 关）。只用 R5/R6/R7/R8/R12。不含 R10/R11，不含炭 R1–R4。不是生产 `_freeboard_reaction_step` / 轨迹积分。
- 床内气泡 R8 同样因 \(\rho_{\mathrm{ash}}=0\) 关闭。
- 焦油无相间交换；\(C_{\mathrm{tar}}\) 取全格摩尔分率。
- 默认路径 \(Q_W=0\)：E5 在登记上取消了壁面热省略，能量行数值上仍无壁面热。
- 两级分队列在 F_E0 物种核里；生产 `calc_size_migration` 未移植。
- 旋风气接口是拷贝约束，不是独立反应格。

---

## 4. 未完成项，以及暂时未启动的原因

下面每一项都写了「是什么 / 为什么现在不动 / 动之前要先定什么」。没有登记成 F_E9。

### 4.1 登记上仍省略、且本谱系声明永不取消

| 项 | 现状 | 未启动原因 |
|---|---|---|
| `cyclone_K_update_P` | `omitted_named(8)` 仍在；E8 断言 `q.k == nq.k` | 登记 `_CANCEL` 明确「旋风 K 永不取消」。E8 计划只把抛射 \(K\) 放进 \(P(x)\)。旋风 \(K\) 进外层会改回流分流与灰库存，属于另一套 \(P\)，不是反应延拓的下一刀。Q2 对捕集/停留的工程常数尚未取得可进主线的原文闭合。 |
| `production_hamel_outer_refresh` | 生产 `vorab_outer_source_policy` 不是该值；具名 `composite_P` 只覆盖 DAEM + \(K_{\mathrm{eject}}\) | Q4 列了开放前提 ①–⑦（元素守恒、焓参考态、\(F\) 无历史、内层 q 冻结、外层不重写 N/T、\(q_{k+1}=P(x_k)\) 与最终 \(F\) 双验收、网格可解释）。这些没有齐。打开它会把生产默认从 single-shot 换成整包刷新，与「不切生产」直接冲突。 |
| `return_leg_holdup` | 在 `OMITTED_BASE`，**不在** `_CANCEL` | Q6 第一张联立 F 规定本轮无返料腿库存、无额外未知量。加上它要扩维（新 \(M\)、新方程、与旋风回流边配对）。在旋风 \(K\) 仍冻结时，新库存方程容易变成「冻结系数 × 新未知量」的病态闭合。 |

### 4.2 反应相关、但故意没开的扩维

| 项 | 现状 | 未启动原因 |
|---|---|---|
| 自由板 R10/R11 | `fb_gas_source` 明确不含焦油 | E7 的取消项是 `freeboard_reaction`，规格写成均相 CSTR。床内焦油已在 R10/R11 根上耗尽或受出口方程约束；把同一 lumped tar 再送进自由板，要先定义焦油如何离开床（现在没有相间 tar 交换，也没有 FB tar 未知量）。 |
| 自由板炭 R1–R4 | 未进 `N_fb` 源 | 夹带炭走轨迹派生占据，不在 43 维 \(x\) 里。Q6：FB 内固体若反应，要么把 \(M_{\mathrm{fb}}\) 纳入未知量，要么做含反应的轨迹消元。冻结旧 \(M\) 却用新 \(N/T\) 消耗，不守恒。E7/E8 都没有选这两条。 |
| 仓库 R9（H₂S） | `rate_R9` 是占位；MAJOR 无硫 | 具名核 7 气不含 H₂S。开它要改物种维数，不是延拓参数 \(s\)。 |
| 生产 `calc_size_migration` | 未移植；两级分用 cohort 代表 | Q7：代码前系数 3 的半径/体积定义未核验，禁止不经推导改成 6，也禁止把生产迁移默认塞进具名核。 |

### 4.3 外层 \(P(x)\) 仍只覆盖两块

当前 `outer_P_named`：`stage>=4` 更新 DAEM \(n_{\mathrm{py}}(T)\)；`stage>=8` 更新 \(\rho_g(T)\)、\(K_{\mathrm{eject}}\)、\(w_{\mathrm{leave}}\)。下面这些 **随状态变、但现在冻在 q 里**：

| 冻住的量 | 未启动原因 |
|---|---|
| 水力：\(u_b\)、\(\varepsilon\)、\(K_{\mathrm{bd}}\)、气泡尺寸 | 它们进交换与浓度。一开就会移动已认证的反应根，无法归因「只开了外层」。需要单独的冻结对照 + 分块验收，不能挂在反应延拓后面当 F_E9。 |
| 占据年龄 / 热解释放权重 \(p_i(t)\) | Q4：这是建议的占据积分，不是论文闭合式。现在 DAEM 用单一 \(\tau\)。 |
| 干燥剩余水、未转化 VM 作为 \(x\) | 现为给料映射进 q。升为未知量要配对输运方程，否则会偷偷增加出口。 |
| 完整 Bild 2.2 外层刷新 | 见 §4.1 的 `hamel_outer_refresh`。 |

### 4.4 求解器与认证上已知、但未另开任务的问题

| 问题 | 证据 | 未作为下一刀的原因 |
|---|---|---|
| \(w_{\mathrm{leave}}\) 是离散 0/1 | \(T=800\,\mathrm{K}\) 粗档可翻；E7 根不翻 | 本根上 \(F_{E8}=F_{E7}\)，FD \(1\,\mathrm{K}\) 的 T 列碰不到跳跃。要处理跳跃，须先决定：事件牛顿、光滑化，还是把命运留在外层网格。未裁定前改雅可比没有验收。 |
| T 列对复合 \(H(x)\) 用 \(1\,\mathrm{K}\) 差商 | E4/E6/E8 同构 | 反应根已用分块认证锁住。加密差商是数值卫生，不改变模型。 |
| E7 小 \(s\) 会把混合约束里的氧加回 `N_fb` | 曾导致无限对分，\(\\|F\\|_2\sim10^{-6}\) | 已改为 Startwert 后直接 \(s=1\)。不要「补全 \(s\) 网格」当回归。 |
| 亏秩时微量组分补对角 | `named_linear_step` | 这是修正 Newton，须保持声明，不得写成 Wirsum 原文 \(J=F'\)。 |
| R8 \(k_{\mathrm{eq}}\) OCR | `docs/hamel_submodels/06_*.md`；代码用 `exp(-3.6893+4019/T)` | 动力学债，不是延拓债。改它会移动所有均相根，须单独相对 PDF 图像验收。 |

### 4.5 明确不做、除非用户改口

| 项 | 原因 |
|---|---|
| 把具名核切进生产 `solve_global_nr` / `outer_loop` | 本链约束：不改生产默认；Q4 未齐；reduced vs thesis_mode 禁止混比。 |
| 用 \(Q_W\) 把 \(T\) 拟合到 \(1100\,\mathrm{K}\) | 禁止。E5 的 \(Q_W\) 是显式 q。 |
| 对 Hamel 图 7 / Table 7 宣称 parity | 当前 \(T\sim1420\,\mathrm{K}\) 是单格数学根，几何/给料是具名冻结 q，不是 HTW 全炉。 |
| 跑 `max_global_iter=150` 的 Phase1 HTW | 旧交接：完整 pytest 会卡在该路径；与具名核无关。 |
| 给 VM 加 \(K\)、恢复 fastox、残差 clip 当根 | 贯穿本链。 |

---

## 5. 待决问题（动手写代码前需要裁定）

下列问题没有默认答案。选错会把已认证的 43 维根作废，或把工程开关标成 thesis。

**Q-A. 下一问题是哪一类？**  
（1）更完整的 \(P(x)\)（水力 / 占据 / 旋风 \(K\)）；（2）自由板焦油或夹带炭；（3）返料腿库存；（4）具名核 vs 生产路径对照，仍不切默认；（5）停在 F_E8，先 commit 未跟踪文件。  
未选定则不要新增 `wirsum_1997_named_f_e9.py`。

**Q-B. 若选旋风 \(K\) 进 \(P(x)\)：用哪条分离律？**  
Hamel Bild 7.7 的 1–2 cell + 冻结 \(K\) 已在 `wirsum_1997_cyclone_cells.py`。Lapple / 常数 \(\eta=0.95\) 已标对照、不得回主线。生产 `freeboard_cyclone_capture_*_frac` 是 engineering。没有第四条「随 \(\rho_g(T)\) 更新的原文 \(K(x)\)」。

**Q-C. 若选自由板固体反应：未知量还是轨迹消元？**  
Q6 只给了这两条合法路线。现在细粉已按命运离开床顶，但 \(M_{\mathrm{fb}}\) 不在 \(x\)。选定之前不要把 `rate_R1` 乘进自由板气体行。

**Q-D. \(w_{\mathrm{leave}}\) 跳跃是否视为本核必须处理的非光滑点？**  
E8 根上不触发。若要把低温命运翻转当成物理分支，需要独立的事件/延拓设计，不能靠现有 \(1\,\mathrm{K}\) FD。

**Q-E. `composite_P` 与生产 refresh 的关系如何命名下去？**  
现在：具名核 `coupling_mode="composite_P"`，生产 refresh 仍关。禁止用 `hamel_outer_refresh=False` 暗示内层 q 全冻（E4 起 DAEM 已随 \(T\) 变，E8 起 \(\rho_g\) 也变）。若扩展 \(P\)，继续用具名字段，不要复用生产旋钮名。

**Q-F. 未跟踪文件何时入库？**  
整包 named-F 仍是 `??`。用户此前要求未明确时不 commit。入库时应与 09-14 的 bed/freeboard/cyclone 约化核一起考虑，避免只交延拓求解器、缺底层 \(F\)。

---

## 6. 验收命令

从仓库根：

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q \
  tests/test_wirsum_1997_named_f_stages.py \
  tests/test_wirsum_1997_named_f_continuation.py \
  tests/test_wirsum_1997_named_f_e8.py \
  --tb=short
```

E8 文件会重算 E6→R10→R11→E7 链，约 8 min。不要用 `pytest -k` 扫多文件：`-k` 会作用到整个收集集。

端点字段：`endpoint_certified`、`endpoint_cert.gas_mol_s_l2`、`energy_W_abs`、`Q_W_is_fitted_to_1100K is False`、`production_hamel_outer_refresh is False`、`not_a_furnace_convergence_claim is True`。

---

## 7. 代码地图（未跟踪为主）

```
src/solvers/wirsum_1997_named_f_registry.py   # 阶段名、省略、n_unknowns
src/solvers/wirsum_1997_named_f.py            # 残差装配、outer_P_named、fb_gas_source
src/solvers/wirsum_1997_named_f_continuation.py
src/solvers/wirsum_1997_named_f_r12.py … r234.py、r10.py、r11.py
src/solvers/wirsum_1997_named_f_e3.py … e8.py
src/solvers/wirsum_1997_startwert_{r5,r6,r10,r11,e7}.py
docs/wirsum_1997_named_f_rate_continuation.md
tests/test_wirsum_1997_named_f_*.py
```

底层物种/能量核（本谱系之前已有，多数亦未整包 commit）：`wirsum_1997_coupled_f_two_bin_pyro.py`、`wirsum_1997_coupled_f_energy.py`、`wirsum_1997_cyclone_cells.py`、`wirsum_1997_bed_freeboard.py`。

生产路径：`src/core/reactor.py`、`cell.py`、`outer_loop` **不是**本核入口。

---

## 8. 建议的下一手（仅在 Q-A 选定之后）

1. 若停在 F_E8：按用户指示 commit 未跟踪具名核，提交说明写清「登记终点，不是全模型」。
2. 若对照生产：只做只读残差/源项对比，不改 `ReactorConfig` 默认。
3. 若开更完整 \(P(x)\)：先冻结对照表（水力、旋风 \(K\)、占据），再一次改一块，每块独立 `endpoint_certified`。
4. 不要从 F_E8 自动开始旋风 \(K\) 或自由板焦油。
