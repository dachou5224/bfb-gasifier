# 交接：A–F 进展、阻碍、待决问题（2026-09-14）

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读：

- 诊断：`docs/audit_reports/secondary_diagnosis_2026-09-14.md`（尤其 §1、§9）
- 实施记录：`docs/audit_reports/secondary_diagnosis_fix_plan_2026-09-14.md`
- 本链长交接：`docs/handoff_wirsum_alignment_2026-09-14.md`
- 原文：同目录 Hamel 1999/2001 VDI PDF、Wirsum 1997 PDF（诊断 §1 已给页码）

**分支：** `feat/nr-hamel-wirsum-core-startwert`  
**状态：** A–F 约化核未整包 commit。未跑全炉 150-iter LU。未宣称 \(q=\mathcal P(x)\) 已协调。  
**语言：** 简体中文。答案必须带方程号 / PDF 页 / 印刷页，或明确写「原文无此式」。

目标仍是重建 **Wirsum/Hamel 原始数学问题**，不是调参让 LU 收敛。两条路径不得混比：

| 路径 | 入口 | 特征 |
|---|---|---|
| reduced 3-cell | `tests.global_nr_fixtures.build_initialized_lu_reactor` | 只 `generate_initial_x0`，`thesis_mode=False` |
| Phase2 init | `build_initialized_thesis_freeboard_reactor` | `thesis_mode=True`，`init_precalc`，源项可 single-shot 冻结 |

生产默认几乎未改。`thesis-exact` 与 `engineering` 必须分名。禁止 clip / fastox / 实验 CO 当数值桥。

---

## 1. 一句话结论

诊断 §9 的 A–F **没有改完**。已拆开问题并做出可复现的约化核；**同一套 \(F(x;q)=0\)、\(q=\mathcal P(x)\) 尚未重建**。下面 §3 是阻碍，§4 是必须由你给出具体答案后我才能继续写代码的问题。§5 是我自己能做、不必等你的活。

---

## 2. 已完成（约化核，除非另注）

生产路径真正改过的，仍主要是更早的 P0.3（`F(x)` 纯度）和 P1.2（`holdup_k0_anchor`，默认仍 `engineering`）。本轮 A–F 多数是独立核 + 具名开关。

| 阶段 | 核 / 证据 | 已对齐诊断的部分 |
|---|---|---|
| A | `src/solvers/wirsum_1997_baseline_contract.py`，`data/wirsum_alignment_stage_a_baselines.json` | 两条基线的 path / 未知量 / cache / 源项冻结 / 工程开关已快照 |
| B | `energy_balance_defaults.py`，`rate_R12`，`tests/test_hamel_r12_pdf_literal.py` | char hf=0；水分=液态水；VM=`h_f,dry/w_VM`；R12 为 \(k C_{H_2}^{3/2}C_{O_2}\) |
| C | `src/core/solid_transport_reduced.py` | 等温无反应；`species_own` vs 15% 工程分母；床底无穿透 opt-in |
| D | `wirsum_1997_vorab_context.py`，`vorab_policy.py` | 冻结后 `F(x),F(x)` 与 `F(x),F(y),F(x)` 为 0；single-shot **具名为工程**，不是 Bild 2.2 |
| E | `wirsum_1997_stiffness_ladder.py`，`wirsum_1997_stiffness_full_gas.py` | 输运→交换→R12→全气→炭；`physical_convergence` 仅 `model_is_full_rate` 且 \(\|F\|_2<10^{-6}\) |
| F1 | `wirsum_1997_bed_freeboard.py` | 关回流；`N_fb=N_d+N_b`；Roof 3.59 vs Roof+Wake 3.54 |
| F3 | `wirsum_1997_reaction_coord.py` | 离开炭与内层 R1 共用 `char_r1_extent`；入口 e1≠根上 e1 |
| F4–F5 | `wirsum_1997_cyclone_cells.py`，`wirsum_1997_bed_fb_cyclone.py` | Bild 7.7 的 1–2 cell；`m_rez` 进床残差；同一 `N_fb`/`m_leave` 交接 |

**已撤回 / 已标对照、不得当主线：**

- Lapple / 常数 \(\eta=0.95\)（`src/physics/wirsum_1997_cyclone.py`，F2）。那是本链自作主张，不是诊断也不是 Hamel。
- 热解气 Startwert 不当本轮主线（诊断 §8）。
- 默认 `thesis_vorab_sources_single_shot`、默认 15% 库存未改。

相关测试（本机）：F / F3 / F4 / F5 约化核已绿。完整 `pytest -q` 会卡 Phase1 HTW `max_global_iter=150`，不要当默认环。

---

## 3. 阻碍项（不是「再写一个核」就能过）

这些卡住的是 **原文缺口、两套方程尚未接到同一 \(F\)、或生产默认不能 silently 改**。

1. **没有接到生产，就不能宣称诊断完成。** 约化核绿 ≠ LU 已修。诊断验收多次写「不能把 toy 通过写成生产已修」。
2. **外层 \(q=\mathcal P(x)\) 未实现。** 生产外层仍不刷干燥/热解/自由板轨迹。D 只证明：冻结则 `F` 纯，不冻结则历史依赖。
3. **F 不是一条联立残差。** 界面、示踪/T、R1、旋风 cell 是四套核。诊断要的是「轨迹更新固体，内层解气/T，同一物流与反应量」。
4. **能量未过关。** B 的固体瞬时 \(c_p\) 未改；E 能量阶 `singular_jacobian`，不含固体 cp 积分。不能用 x0 能量残差升降验收物性。
5. **旋风 \(K\) 原文没有计算公式。** F4 标 `frozen_context`。补经验 \(\eta\) 会再偏出口径。
6. **轨迹阻力不是 Hamel 扩展的全文。** 抛射用了 3.53/3.54/3.73；命运 ODE 仍用 Wirsum \(c_w=24/\mathrm{Re}+z\)。
7. **粒径迁移、热解一次映射、R9–R11、全炉元素/能量** 均未进本链核。
8. **提交策略：** 用户未要求则不 commit。工作区相对 origin 已 ahead，另有大量未跟踪 WIP。接手不要整包提交。

---

## 4. 需要你给出具体答案的问题

请按题作答。每题的「合格答案」写在题末。不要给「可以先用工程近似」除非你同时给出 **原文无此式** 的书面结论和具名开关名。

查原文时用诊断 §1 的页码表（PDF 页 vs 印刷页）。图像放大优先于 OCR。禁止用实验 CO / Fig.7 / 工业测点决定开关。

### Q1. 自由空间轨迹：Hamel 扩展的阻力式到底是哪一条？

**现状。** `src/physics/wirsum_1997_freeboard.py` 里：

- `wirsum_1997_exact`：\(u_{p,0}=1.2 u_b\)，只 Roof 3.59，\(c_w=24/\mathrm{Re}+z\)，\(z\in[0.43,1]\)
- `hamel_extended`：\(u_{p,0}=1.53 u_b\)（Martens / Hamel 3.73），Roof+Wake 3.54，**阻力仍走 Wirsum ABC**
- 生产 `simulate_freeboard` 的 `analytical_wirsum` 是工程混合，已声明不是 Wirsum 原文
- 生产另有 `calc_cd_hamel_eq384`（`src/physics/freeboard.py`）：Haider 型 \(c_W(\mathrm{Re},\varphi_s)\)，含 \(73.69\exp(-5.0784\varphi_s)\)

诊断原文：「轨迹闭合采用 Hamel 的扩展形式，Wirsum 可作为来源/对照，不应把 Wirsum 的较早模型等同于完整 Hamel。」

**我不确定。** 「Hamel 扩展」是 (a) 只加 Wake + Martens 抛射，阻力仍用 Wirsum §3.4，还是 (b) 必须换成 Eq. 3.84 / 某 \(z(\varphi)\)，还是 (c) 生产 `global_projected_hamel` 才是 thesis？

**请给出：**

1. 逐项表：平均抛射 \(u_{p,0}/u_b\)、Roof、Wake、阻力 ODE、速度高斯 \(\sigma\)，每一项的 **Hamel 方程号 + PDF 页 + 印刷页**，或「Wirsum only」或「原文无」。
2. 一句裁定：`hamel_extended` 的命运分类现在该不该改阻力。若该改，写出 \(du_p/dt=\ldots\) 的应采用式（含 \(c_W\) 或 \(A,B,C\)），以及 \(\varphi_s\) 从哪来。
3. 若 Eq. 3.84 与 Wirsum ABC 在 Hamel 正文中的关系（替代 / 并列 / 仅评价用）。

**我拿到答案后会做：** 只改 `wirsum_1997_freeboard.py` 的 `hamel_extended` 规格 + 命运分类；生产 `simulate_freeboard` 默认不动。

**不合格答案：** 继续用 Lapple；用「和实验夹带更像」选阻力；把 `analytical_wirsum` 改名为 Hamel 原文。

---

### Q2. 旋风 cell 的 \(K_\mathrm{auf}\)、\(K_\mathrm{ab}\) 原文有没有式？

**现状。** Hamel §7.1 / Bild 7.7：Zyklon 用 1–2 个 cell，与床同一套 Eq. 2.4–2.6。VTT 明确 **没有旋风几何**，仍假定两格。论文没有分级效率曲线。

F4 用冻结 \(K_\mathrm{auf}=0.25\,\mathrm{s}^{-1}\)、\(K_\mathrm{ab}=0.75\,\mathrm{s}^{-1}\)、\(K_\mathrm{return}=1\,\mathrm{s}^{-1}\)，标 `frozen_context`。生产用停留时间 \(K=1/\Delta t\) 再乘常数捕集分数，那是工程。

**我不确定。** 旋风分离在 thesis 里是 (a) 仅靠 cell 衡算 + 外层给的 \(K\)（\(K\) 属 \(q\)），(b) 有未找到的 \(K\) 或 \(\eta(d_p)\) 式，还是 (c) 退化水力（`eps_b=0`）自动决定固体上下行？

**请给出：**

1. 在 Hamel / Wirsum 全文检索 Zyklon / Abscheide / \(K\) / Feststoff 后的结论：有式则抄原式（符号、单位、页码）；无式则写 **「原文无旋风 \(K\) 或 \(\eta\) 公式」**。
2. \(K\) 是否属于外层 \(q=\mathcal P(x)\)。若是，\(\mathcal P\) 用什么输入（\(N,T,P\)、几何、粒径）更新 \(K\)。
3. 返料腿是否必须单独成 cell（Bild 7.7 的 Rückführung），还是旋风下行直接 `m_rez` 进床即可。
4. 生产 `freeboard_cyclone_capture_*_frac` 与 `_estimate_cyclone_residence_time_s` 分别是不是 thesis。各给 `thesis` / `engineering` / `unknown`。

**我拿到答案后会做：** 有式就实现具名核；无式就保持 `frozen_context`，并在文档写死「不得用 Lapple/常数 \(\eta\) 填坑」。可能加 opt-in 把生产捕集分数标成 engineering。

**不合格答案：** 再引入 Lapple / Barth / Leith-Licht 并称 Hamel；用 HTW 捕集率反推 \(K\)。

---

### Q3. 固体焓：\(c_p(T)(T-T_\mathrm{ref})\) 要不要改成积分？

**现状。** 诊断 §2.4：`cp_char`、`cp_ash` 随 \(T\) 变，当前用 \(c_p(T)(T-T_\mathrm{ref})\)。若是瞬时比热应积分；若是平均比热须给出参考温度和来源。本轮 **故意没改**。E 能量阶也不含该积分，且停在 `singular_jacobian`。

**我不确定。** Hamel Eq. 2.8 / 第 5 章固体焓到底是瞬时积分、平均 \(c_p\)，还是多项式焓？与气体 `enthalpy_molar`（Gordon–McBride，含 hf）如何对齐参考态？

**请给出：**

1. 原文对固体 \(h_s(T)\) 的定义（方程号、页码）。写成
   \[
   h = h_f + \int_{T_\mathrm{ref}}^{T} c_p(T')\,dT'
   \]
   或 \(h=h_f+\bar c_p(T-T_\mathrm{ref})\)，并指明 \(T_\mathrm{ref}\)。
2. `cp_char` / `cp_ash` 的文献来源：瞬时还是平均。若代码与原文不一致，给出应改的公式，不要只说「最好积分」。
3. 改完之后，E 能量阶的未知量集合：是否仍是「全气 + 炭库存 + 一个 \(T\)」？奇异常不常被固体 cp 解释？若你判断奇异是另一原因，指出哪一列/哪一行线性相关（单位、变量名）。

**我拿到答案后会做：** 只在约化能量核改固体焓；生产默认先不动，除非你明确说生产 `solid_enthalpy_mode` 已错且给出替换名。

**不合格答案：** 用「改完 x0 的 \(E\) 下降了」当物性验收（诊断 §2.4 禁止）。

---

### Q4. 干燥/热解停留：能否继续用 \(H/u_{mf}\)？

**现状。** `vorabrechnung_tau_for_cell` 具名为 `engineering_H_over_umf`。诊断 D：须用有来源的颗粒停留/释放分布，不能为了接近完全而取 \(H/u_{mf}\)。Wirsum §2.2：释放量和释放位置映射到各格，预计算与单元模型迭代协调。默认外层 **不** 重算干燥/热解（`engineering_single_shot`）。

**我不确定。** 原文停留/轴向分配的具体式是什么？Hegermann Tabelle 4.3 + Gibbs 已在仓库里，缺的是 **时间/位置映射** 还是 **载体（VM 固体列 vs 气相 \(N_r\)）**？

**请给出：**

1. Hamel / Wirsum 中 \(\tau\) 或释放分配 \(\chi_i\) 的原式（页码）。若只有定性描述，逐句引用并说明代码里不能合法写成 \(H/u_{mf}\) 的原因。
2. \(q\) 中热解部分的最小变量清单：各格 \(N_{r,j}\)、焦/焦油质量、剩余元素、焓载体。指出哪些进 \(x\)、哪些进 \(q\)。
3. 现有 `hegermann_yield.py` + `pyrolysis_carbon_split.py` 是否已覆盖「先分 Koks/Teer 再 Gibbs」；缺的是否只是 **写入气相 \(N\) 且从固体 VM 列撤走同一批 C/H/O/焓**。
4. `hamel_outer_refresh` 何时允许打开：必须先满足哪几条验收，才不会变成「为了降残差而刷源项」。

**我拿到答案后会做：** 先写约化核：一份 \(q\) 含 \(\chi_i\) 与气相源，内层只读；**不** 给 VM 加 \(K\)，**不** 把 VM 写入 Eq. 2.7 固体出口。默认 single-shot 先不改，除非你列出打开条件且已满足。

**不合格答案：** 用热解气 Startwert 当主线修 VM 焓桶；把剩余 C 划回 Koks。

---

### Q5. 阶段 E 能量阶为何奇异？下一步改公式还是改未知量？

**现状。** `wirsum_1997_stiffness_full_gas.py` 能量阶：完整速率，`singular_jacobian`。气体用 `enthalpy_molar`；炭 hf=0；无固体瞬时 cp 积分。诊断 E 要求能量残差达标，且小步长不得标物理收敛（已分名）。

**我不确定。** 奇异来自 (a) 能量行与某物种行不独立，(b) 缺固体焓使 \(\partial E/\partial T\) 或 \(\partial E/\partial M\) 退化，(c) 约化状态（无 tar、无壁面、等温几何）本身少一个独立能量方程，还是 (d) FD/缩放问题（诊断 §7.1）？

**请给出：**

1. 在该约化问题里，能量方程对哪些未知量的导数必须非零（列出 \(\partial E/\partial T\)、\(\partial E/\partial N_j\)、\(\partial E/\partial M\) 的应有单位和主导项）。
2. 奇异是 **结构秩缺陷** 还是 **尺度造成的数值秩不足**。给一个可复现的判别（例如列缩放后再看 SVD）。
3. 若结构缺陷：缺的是哪一项原文（壁面 \(Q_W\)、固体 cp、相变、还是不该把 \(T\) 放进这个核）。
4. 在缺陷消除之前，是否禁止把该阶标 `implemented` 为可验收。我倾向保持「已实现、未物理收敛」。

**我拿到答案后会做：** 按你指的项改约化能量核；仍然分名 `small_step` / `singular_jacobian` / `physical_convergence`。

**不合格答案：** 加 T fence、能量 clip、或恢复 fastox 让 Jacobian「看起来满秩」。

---

### Q6. F 联立残差的最小合法未知量集合是什么？

**现状。** 诊断 F：轨迹更新固体；内层解气/T；同一物流与反应量；外层判据含 T、两相、水力/交换、固体双向与库存、热解、床面回落/夹带、全炉元素与能量。

现有分核：

- 界面：`N_d,N_b\to N_fb`，抛射/回落/离开
- 示踪交换 + 具名 \(Q_W\) 定 \(T\)（F2，工程焓身份）
- R1 + \(M=m_\mathrm{leave}/K_{fb}\)（F3；\(K_{fb}\) 是 `engineering_K_fb_leave`）
- 旋风 2 cell + 返料（F4/F5；\(K\) 冻结）

**我不确定。** 下一步是把这些 **塞进一个 \(F(x)\)**，还是先保持分核、只统一物流句柄？自由空间在 Hamel 里是单相还是仍保留床的两相？\(K_{fb}=1\,\mathrm{s}^{-1}\) 有没有原文，还是必须删掉、改用轨迹停留？

**请给出：**

1. 一张表：每个未知量（名称、单位、所在 cell 类型：bed / fb / cyc / return）、对应方程（2.4/2.5/2.6/2.7 或轨迹，不是残差）。标 `x` 或 `q`。
2. 自由空间气相：单相总量（现 F1）还是两相。给方程号。
3. 离开炭库存：由轨迹通量直接进旋风，还是 FB 格也有 \(M\) 与 \(K\)。若有 \(K\)，式从哪来。
4. 第一张「联立 \(F\)」允许省略什么（写 `omitted, named`）：例如无热解、无 R2–R4、无能量、旋风 \(K\) 冻结。
5. 外层判据 10 项里，这一张 \(F\) 必须覆盖哪几项才算 F 的第一张合格联立，哪几项可以留到下一张。

**我拿到答案后会做：** 只实现你列出的 \(x\) 与方程；省略项写进 payload。不把 F2 的 \(Q_W\) 焓身份冒充 Eq. 2.7。

**不合格答案：** 「先全部耦合并跑 LU」；用 \(\eta\) 代替旋风 cell。

---

### Q7. 粒径迁移要不要在接到 F 之前进入阶段 C？

**现状。** 诊断 C：先无反应等温双向输运，再加入反应与粒径迁移。约化 C **没有** 粒径迁移。生产有 `calc_size_migration`，默认 `per_class_holdup`。

**请给出：**

1. Hamel Eq. 2.6 的 \(m_\mathrm{in/out}\)（级分转移）原式或算法步骤（页码）。
2. 生产 `calc_size_migration` 与原文是 `thesis` / `engineering` / `unknown`。若 unknown，缺哪一段。
3. 本链是否必须在 F 联立之前把迁移纳入 `solid_transport_reduced`。是/否 + 一句理由（按诊断顺序，不是按「功能完整」）。

**我拿到答案后会做：** 要则只进约化 C，默认生产开关不动；不要则在 C 文档写死「迁移不在本核，F 亦不假装有多粒径闭合」。

---

### Q8. 何时允许改生产默认？

**现状。** `species_own_inventory`、`hamel_outer_refresh`、R12 新式、生成焓修正：有的已进生产（B 的焓和 R12），有的仅 opt-in。诊断怕 silently 改 LU。

**请给出一张表**，每行一个开关/公式：

| 对象 | 现在 | 是否已是 thesis | 允许改生产默认的前提（可检验） | 现在改？ |
|---|---|---|---|---|

至少包括：`rate_R12`、固体 hf 向量、`thesis_vorab_sources_single_shot`、`solid_transport_inventory_mode`、`holdup_k0_anchor`、自由板 `trajectory_model`、旋风捕集分数。

**我拿到答案后会做：** 只改你标「现在改？= 是」且前提已满足的项；其余保持分名开关。

---

### Q9. 诊断 §7：生产 residual 的投影如何处理才算「原式」？

**现状。** 通用 `unpack_cell()` 仍将负物流裁零、温度裁到上下限。诊断：合法状态域宜由试探步保证，越界拒绝，不要改传入的 \(x\)。约化 strict 核本身不投影；把生产 \(F\) 交给它 ≠ 无投影原式。

**请给出：**

1. 哪些投影是「真实物理分段」（例如 \(N\ge 0\) 作为互补条件）必须保留在模型里，哪些是纯数值、必须移出 \(F\)。
2. 移出之后，越界步的合法行为：拒绝 \(\alpha\)、还是扩域后让 \(F\) 在域外有定义。
3. 是否允许在 **约化核验收完成之前** 动生产 `unpack_cell`。我倾向否。

---

## 5. 不必等你、我可以继续做的

仅当 §4 未回答时，下一手仍限制在这些：

- 保持两条基线快照可复现；新核继续写 `data/wirsum_alignment_*.json`。
- 文档上继续把 Lapple / \(\eta=0.95\) / \(H/u_{mf}\) / single-shot 标成 engineering。
- 不跑 150-iter LU；不给 VM 加 \(K\)；不恢复 fastox。
- 若你只先答 Q1+Q2+Q6，优先做 **F 联立的最小 \(F\)**（按你的未知量表），而不是再开新的对照经验式。

---

## 6. 文件地图（本轮）

| 路径 | 阶段 | 注意 |
|---|---|---|
| `src/solvers/wirsum_1997_baseline_contract.py` | A | 两条路径契约 |
| `src/core/energy_balance_defaults.py` | B | 生成焓已改生产 |
| `src/kinetics/gas_reactions.py` `rate_R12` | B | 已改生产为扫描式 |
| `src/core/solid_transport_reduced.py` | C | 默认生产仍 15% |
| `src/solvers/wirsum_1997_vorab_context.py` | D | 不改默认 single-shot |
| `src/solvers/wirsum_1997_stiffness_ladder.py` | E | 能量未物理收敛 |
| `src/physics/wirsum_1997_freeboard.py` | F | exact vs extended 分名 |
| `src/physics/wirsum_1997_cyclone.py` | 对照 | **不是主线** |
| `src/solvers/wirsum_1997_cyclone_cells.py` | F4 | Bild 7.7 |
| `src/solvers/wirsum_1997_bed_fb_cyclone.py` | F5 | 同一 `N_fb`/`m_leave` |
| `src/solvers/wirsum_1997_reaction_coord.py` | F3 | 未接到 F5 |

证据 JSON：`data/wirsum_alignment_stage_{a,c,d,e,f,f2,f3,f5}_*.json` 以及 `data/wirsum_alignment_stage_f_cyclone_cells.json`。

---

## 7. 你回交接时请用的格式

对 §4 每一题：

```text
### Qn 裁定
- 结论：（不超过五句）
- 原文：（方程号，PDF 页，印刷页；或「原文无」）
- 代码应改：（文件 / 函数 / 具名模式）
- 现在不要改：（明确列出）
- 验收命令：（一条 pytest 或一条 20 行探针，不要全炉 LU）
```

答不全就标 `unknown`，不要用工程式填成 thesis。
