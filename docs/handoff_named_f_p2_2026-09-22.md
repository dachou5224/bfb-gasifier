# 交接：具名 \(F\) P2 有限水力（2026-09-22）

**给下一任 agent。** 先读 `docs/CLAUDE.md`，再读：

1. 本文件（当前工作面、已闭合、阻塞、下一手）
2. 审查与阶段计划：`docs/audit_reports/named_f_e8_review_2026-09-16.md`（尤其 §2.3–2.4、§5 P2）
3. 上一份长交接：`docs/handoff_named_f_e8_2026-09-16.md`（E0–E8 登记谱系、硬约束）
4. 延拓算法：`docs/wirsum_1997_named_f_rate_continuation.md`
5. 上一会话：[具名核 P2 延拓](8b68c0b2-eafa-4d1b-9109-aaae7c282a1c)

**分支：** `feat/nr-hamel-wirsum-core-startwert`（相对 `origin` ahead 22；工作树另有大量未提交改动）  
**日期：** 2026-09-22  
**Git：** 具名 \(F\) 整包仍未 commit。已提交 HEAD 仍是 `098c8cf`（holdup Startwert 轴向接力）。用户未要求则不 commit。  
**状态：** P0 / P1 / **P2 映射验收已闭合**。E8 在完整有限水力 \(P(x;b)\) 下**没有**认证根；快照未按 P2 刷新。P3–P4 未启动。  
**语言：** 简体中文。

目标仍是重建 Wirsum/Hamel 原始代数问题 \(F(x;q)=0\)（SI），不是调参让生产 LU 收敛。本核不经 `Reactor`。不要混 reduced 3-cell 与 `thesis_mode`。

---

## 1. 一句话

P2 要的 \(q=P(x;b)\) 已经装上：理想气体 \(u_g,\rho_g\)、Hilligardt \(u_b\)、Hamel \(K_{bd}\)、Eq. 3.72 初速分布命运、未判定不是逃逸、旋风 \(k\) 冻结。验收测试（无求根）已过。  
**卡住的是**：从已认证 E7 根走到 \(s=1\) 的 \(H(x)=F(x;P(x;b))=0\)。\(q\) 插值 + 冻结 \(q\) Picard 可以走到 \(s=0.6\)；再往完整水力走，Newton 方向对非线性残差不再下降。不要把 E7 根写成 E8 认证，不要放宽 `_PHYSICAL_L2` / `ENDPOINT_L2`。

下一手是 **E8 水力延拓到 \(s=1\) 并现场重认证、重写快照**。不要开 P3，不要开旋风 \(K\) / 返料腿 / 生产 `hamel_outer_refresh`。

---

## 2. 硬约束（接手不得违反）

| 项 | 规定 |
|---|---|
| 单位 | 全程 SI。Arrhenius 只走 `src/kinetics/arrhenius.py` |
| 速率式 | 不改 `rate_R1`–`rate_R12` 的阶数、温度幂、前因子。审查里「不得改 rate」= 不得无依据调参；原文单位/相别错了才许修 |
| 生产默认 | 不切 `hamel_outer_refresh`。`thesis_vorab_sources_single_shot` 保持生产默认 |
| 炉温 | 不把根按到 \(1100\,\mathrm{K}\)。`GEO_FB_RHO_G_T_REF_K=1100` 只是冻结参考温度，不是拟合炉温 |
| \(s<1\) | 不是物理模型。完整速率端点用 `certify_named_full_rate_endpoint`，不用混合 \(\\|F\\|_2\) 当门 |
| 路径 | 不混 reduced 3-cell 与 `thesis_mode` |
| 验收阈值 | 气体/固体/焦油 L2 \(<10^{-6}\)（mol/s 或 kg/s），能量 \(\lvert E\rvert<1\,\mathrm{W}\)。**不要放宽** |
| 未开项 | 旋风 \(K\) 更新、返料腿库存、生产 `hamel_outer_refresh`、P3 分区能量、P4 对照生产 |
| 提交 | 用户未要求则不 commit |

相关技能：`.cursor/skills/bfb-diagnosing-bugs`、`bfb-hamel-parity-audit`、`bfb-numerical-guard`、`bfb-scientific-audit`。  
`bfb-phase-manager` 的 Phase 1–6 是生产路径门禁；本工作是具名核旁路，不要用它来改 `reactor.py`。

---

## 3. 已闭合（按审查 §5）

### 3.1 登记链 E3–E8（P2 之前）

未知量 E6–E8 = 43。E4 起 `composite_P`：\(H(x)=F(x;P(x))\)。  
R11 用缩放气泡预测子 \(\xi=\Delta s\cdot N_{\mathrm{tar}}\) + 气泡内反向 WGSR，**不要**把完整 holdup SW 倒进 \(N_d\)。  
E7 默认 Startwert 后直接 \(s=1\)。  
现场认证入口：`certify_named_full_rate_endpoint`（`src/solvers/wirsum_1997_named_f_continuation.py`）。

P2 改 \(P\) 之前，E7 根仍现场认证（2026-09-21 探针）：

- 气体 L2 \(\approx 1.25\times 10^{-8}\)
- 能量 \(\approx 2.4\times 10^{-9}\,\mathrm{W}\)
- \(T\approx 1411.63\,\mathrm{K}\)

快照 `data/named_e8_endpoint_2026-09-17.json` 里的 **E7 \(x\)** 仍可用作初值（加载后必须 `recertify_named_e7`）。  
其中的 **E8 \(x\)** 是 P2 之前的根（当时 \(F_{E8}\approx F_{E7}\)），**不能**当 P2 后的认证证据。

### 3.2 P0（雅可比 / 证据块）

`jacobian_named_e8`：N/M 列在冻结 \(q=P(x)\) 下对 \(F(x;q)\) 差商；T 列对 \(H(x)\) 差商（\(\Delta T=1\,\mathrm{K}\)）。  
P0 测试要求 T 列与 \(1\,\mathrm{K}\) 差商一致，不要改这个装配来迁就延拓。  
线性步四块证据：`raw_J` / `modified_matrix` / `feasible_direction` / `rollback`。  
`apply_e8_certified_rollback`：仅当**新点未认证且旧点对 F_E8 已认证**才退回。P2 后 E7 根对 F_E8 **不认证**，回退不会触发。

测试：`tests/test_wirsum_1997_named_f_p0.py`。

### 3.3 P1（物料 / 能量账）

- 燃料碳用 Hegermann / 煤 waf C 基准（`src/solvers/wirsum_1997_named_f_p1_ledger.py`，账本 `wirsum_1997_pyrolysis_ledger.py`）
- DAEM 未释放 VM 走 `unconverted_solid_outlet_throughput`，不再从出口消失
- 焦油稀组分槽与 EOS 分开；C/H/O/灰/焓独立账
- P1 测试：`tests/test_wirsum_1997_named_f_p1.py`（6 项，此前已过）

### 3.4 P2 映射（本轮已验收）

实现：`src/solvers/wirsum_1997_named_f_p2.py`，E8 的 `outer_P_named`（`named_f.py` stage≥8）调用它。

| 件 | 做法 |
|---|---|
| 不可变 \(b\) | `NamedReferenceInputs`，`freeze_stage_q` 时 `capture_named_reference`。P 只读 \(b\)，不得把更新后的 `geo_fb.rho_g` 再当 \(\rho_{\mathrm{ref}}\) |
| \(u_g,\rho_g\) | 理想气体：\(u_g=\dot n R_g T/(PA)\)，\(\rho=(P\,\mathrm{MW})/(R_g T)\)。床层气 = \(N_d+N_b\) |
| \(u_b\) | 冻结 \(d_b\) 上 Hilligardt Eq. 3.13；策略 `LEGACY_HYDRO_BUBBLE_VELOCITY_STRATEGY=hilligardt_eq313` |
| \(K_{bd}\) | Hamel Eq. 3.50 |
| 命运 | Eq. 3.72 高斯梯形（**不是** `truncated_pm_2sigma_uniform_xi`）。默认 11 bins，\(\mathrm{d}t=2\times10^{-3}\,\mathrm{s}\)，\(t_{\max}=5\,\mathrm{s}\)，\(\sigma=0.60\,\mathrm{m/s}\)，Martens \(u_{p0}=1.53 u_b\) |
| 未判定 | `w_undecided` 是占据，**不是**逃逸。\(w_{\mathrm{leave}}+w_{\mathrm{return}}+w_{\mathrm{undecided}}=1\)。P3 之前没有 \(M_{fb}\) 未知量 |
| 冻结 | \(V_d,V_b,\varepsilon_b,d_b\) 在 \(b\) 里；旋风 \(k\) 冻结；`psi_b=0.76` 工业分布板 |
| 密度缩放 | 已去掉孤立的 \(\rho_{\mathrm{ref}}\cdot 1100/T\)。`rho_g_fb_ref` 保持 \(b.\rho_{g,\mathrm{ref}}\) |

轨迹带时间：`src/physics/wirsum_1997_freeboard.py` 的 `integrate_particle_fate_with_time` → `(fate, t_event, z_end)`。

测试：`tests/test_wirsum_1997_named_f_p2.py`（无历史 / \(u_g\) 随总气量 / Eq. 3.72 节点 / bins·dt·\(t_{\max}\) 加密 / F 与 J 同 \(q\)）。  
J 对照只查固体行 `[14,15,37,38,39]` 与列 `M_coarse`（rtol=0.05），避开 R12 大数行。

**最后一次映射套件（2026-09-21，P2 延拓接线后仍过）：**

```bash
PYTHONPATH=. .venv/bin/python -m pytest \
  tests/test_wirsum_1997_named_f_p2.py \
  tests/test_wirsum_1997_named_f_p0.py \
  tests/test_wirsum_1997_named_f_e8.py \
  -k "not test_e8_saved_endpoint and not test_e8_solve"
```

**21 passed, 2 deselected。** 这两个 deselected 在 P2 的 \(P\) 改变后**预期失败**，直到 E8 根重算并刷新快照。

---

## 4. 未完成：E8 在完整有限水力下的根

### 4.1 为什么 E8 不再等于 E7

在快照 E7 根上，`outer_P_named(..., freeze_e8_q())` 给出（约数）：

| 量 | E7 冻结 q | P2 有限 P |
|---|---:|---:|
| \(u_g\) | \(1.20\,\mathrm{m/s}\) | \(\approx 0.16\,\mathrm{m/s}\) |
| \(\rho_g\) | \(0.80\,\mathrm{kg/m^3}\) | \(\approx 6.07\,\mathrm{kg/m^3}\) |
| \(K_{bd}\) | \(6.57\,\mathrm{s}^{-1}\) | \(\approx 0.65\,\mathrm{s}^{-1}\) |
| \(w_{\mathrm{leave}}\) | \([1,0]\)（平均初速） | \([0,0]\)（分布积分，全回流） |

冻结 \(1.20\,\mathrm{m/s}\) / \(0.80\,\mathrm{kg/m^3}\) 是常压烟气量级，与 \(P=2.5\,\mathrm{MPa}\) 的摩尔流不一致。这是 P2 要暴露的事实，不是实现错误。  
此时 \(\|F_{E8}(x_{E7})\|\approx 7.19\)，能量行 \(\approx -7.14\,\mathrm{W}\)，气体 L2 \(\approx 0.87\,\mathrm{mol/s}\)。直接 Newton 4 步线搜索失败。

### 4.2 已试过、不要重走的死路

1. **残差混合** \(F_7+s(F_8-F_7)\)：任意 \(s>0\) 都在 \(F_8\) 里代入**完整**全回流 \(P\)。\(s=10^{-4}\) 就 singular / 线搜索失败。
2. **从 E7 直接 \(s=1\) Newton**（含 certified rollback）：E7 对 F_E8 不认证，回退不触发；\(\|F\|\approx 7\)。
3. **整步冻 T**：能量行几乎只能靠 \(T\) 消；冻 T 后第一步就不下降。
4. **能量专用 \(T\) + 冻 T 的 N/M**：能量可到 \(\sim 10^{-11}\,\mathrm{W}\)；床层相间交换残差留下，且 N/M 方向对非线性 \(F\) 不下降。
5. **把贴零物种抬到 \(10^{-18}\)**：毁掉已走到的 \(s=0.6\) 根。
6. **放宽 `ENDPOINT_L2`**：禁止。

能量行对 \(T\) 近似线性：\(\partial E/\partial T\approx -158\,\mathrm{W/K}\)。\(s=10^{-4}\) 时能量专用 \(\Delta T\approx -4.5\times10^{-6}\,\mathrm{K}\) 即可把 \(E\) 打到 \(10^{-9}\)。  
完整 Newton 的 \(\Delta T\) 大约大 \(25\times\)（耦合进 R12），会把自由板化学打飞。

### 4.3 当前仍在代码里的延拓

`src/solvers/wirsum_1997_named_f_e8.py`：

- `blend_hydro_named_q(q7, q8, hydro_scale, fate_scale=None)`  
  气速/密度/\(K_{bd}\)/\(D_g\) 走 `hydro_scale`；\(w\)/\(K_{\mathrm{eject}}\) 走 `fate_scale`（默认与 hydro 同步）。旋风 \(k\) 不混。
- `residual_e7_plus_scaled_e8` = \(F(x;q(s))\)，不是 \(F_7+s(F_8-F_7)\)。
- `_picard_at_hydro_scale`：冻结 \(q(s)\) 用活跃集 Newton，再按 \(P(x)\) 刷新；复合 \(H\) 写进 `last.F`。
- `_e8_active_step`：\(x_i<10^{-14}\) 的未知量从线性步拿掉（E7 根上细炭 \(M\) 与多数 \(\mathrm{O_2}\) 已经贴零）。
- 计划尺度：`(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)`，失败则 `shrink_scale`。

**已证实：** 这条 Picard 路可以从 E7 根走到 **\(s=0.6\)**，复合 \(\|H\|<10^{-6}\)。  
检查点（若机器上还在）：`/tmp/named_chain/e8_hydro_x.npy`（\(s=0.6\)，\(T\approx 1411.47\,\mathrm{K}\)）。  
`/tmp/named_chain/{e3,e4,e5,e6,r10,r11,e7,e8}_x.npy` 是 P2 之前的链；其中 `e7_x.npy` / 快照 E7 \(x\) 仍应能 `recertify_named_e7`。

**墙在 \(s=0.6\to 1\)：** 即使 \(\Delta s=2\times10^{-4}\)，复合 \(\|H\|\) 按 \(\sim 7\cdot\Delta s\) 跳起（能量 + 成对床层交换）。Picard 第一步线搜索后停住。  
贴零指标（\(s=0.6\) 检查点）：`[5, 14, 21, 28, 35, 37, 38]`（密相/自由板/旋风 \(\mathrm{O_2}\)，细炭 \(M\)，旋风固体库存）。  
留下的气体残差是 **成对的** \(F_d\approx -F_b\)（CO / CO2 / H2 / H2O / N2），即 \(K_{bd}\) 从 \(6.57\to 0.65\) 的交换再平衡，不是自由板源。

线性化 \(F\cdot Jp = -\|F\|^2\) 成立，但沿该方向的真实 \(F\) 上升（哪怕 \(\alpha=10^{-3}\)）。冻 \(q\) 的 J 与非线性 \(H\) 在贴零 \(\mathrm{O_2}\) + R12 处对不齐。

### 4.4 快照

`data/named_e8_endpoint_2026-09-17.json`

- `date` 字段已是 `2026-09-21`，note 仍写 P1/R11，**不是** P2 有限水力根
- `source_hash` 当时 `2efd5f0e…`；P2 改了 `named_f.py` / `e8.py` / 新增 `named_f_p2.py`，hash 已过期
- `endpoint_certified_at_write` **不是**验收；必须 `recertify_named_e7` / `recertify_named_e8`
- 重写入口：`data/named_e8_recompute_2026-09-17.py` 的 `summarize` + `q_eject_snapshot`；或按上一会话脚本从快照 E7 出发只重算 E8

**在 \(s=1\) 现场认证之前不要写快照。**

---

## 5. 建议的下一手（只做这件事）

优先顺序：

1. **相间交换 Startwert（只改 \(x_0\)，不改 \(F\)）**  
   在新 \(K_{bd}\) 下调整 \(N_d\leftrightarrow N_b\) 的分配，保持 \(N_d+N_b\)（从而 \(u_g\)）不变，把成对 \(F_d=-F_b\) 压下去，再 Picard。上一轮手写 1D 转移 FD 把残差弄大了，需要按 `species_residual` 的交换式重做。
2. **伪弧长 / 以 \(\|x\|\) 为参数**  
   \(s=0.6\) 像同伦折点：再增 \(s\) 无下降方向。不要为了过门放松阈值。
3. **分字段延拓**  
   `fate_scale` 与 `hydro_scale` 已经分开。先只开 \(u_g/\rho_g/K_{bd}\)（`fate_scale=0`，保持 \(w=[1,0]\)），再开命运。代码已预留，**Picard 默认仍同步**。
4. 走到 \(s=1\) 且 `recertify_named_e8` 为真之后：  
   重写快照（含 `finite_hydro_P: True`、新 `source_hash`、`q_P`、原始 \(F\)）；跑完整  
   `tests/test_wirsum_1997_named_f_e8.py`（含 `test_e8_saved_endpoint_*` 与 `test_e8_solve_*`）。

不要用旧根 \(T\) 不变当验收。不要开 P3（分区能量 / \(M_{fb}\)）来绕过这个根。

---

## 6. 关键文件

| 路径 | 角色 |
|---|---|
| `src/solvers/wirsum_1997_named_f_p2.py` | \(b\)、水力快照、\(u_g\)、Eq. 3.72、`finite_hydro_from_x` |
| `src/solvers/wirsum_1997_named_f.py` | `NamedStageQ.reference/hydro`、`outer_P_named` stage≥8、`overlay_eject_from_p` |
| `src/solvers/wirsum_1997_named_f_e8.py` | F_E8、q 插值、Picard、活跃集、`solve_named_e8_from_e7` |
| `src/solvers/wirsum_1997_named_f_e7.py` | F_E7；E7 根仍是本阶段入口 |
| `src/solvers/wirsum_1997_named_f_continuation.py` | `ENDPOINT_L2=1e-6`、`certify_named_full_rate_endpoint` |
| `src/solvers/wirsum_1997_named_f_endpoint_io.py` | 快照 IO、现场重认证、源码 hash |
| `src/physics/wirsum_1997_freeboard.py` | `integrate_particle_fate_with_time`（已改，**已跟踪** `M`） |
| `src/solvers/wirsum_1997_named_f_p1_ledger.py` | P1 账（已闭合，勿回退） |
| `tests/test_wirsum_1997_named_f_{p0,p1,p2,e8}.py` | 映射 / 认证门 |
| `data/named_e8_endpoint_2026-09-17.json` | 旧链；E7 可用，E8 过期 |
| `data/named_e8_recompute_2026-09-17.py` | 全链重算脚本（约数分钟；P2 后不要当已通过） |

具名核整包在 git 里是 `??`（未跟踪），外加已跟踪的 `src/physics/wirsum_1997_freeboard.py` 有修改。工作树里还有大量**无关**的生产/脚本脏文件，不要把它们卷进这次提交。

---

## 7. 开工命令

```bash
# 映射仍应通过（约 12 s）
PYTHONPATH=. .venv/bin/python -m pytest \
  tests/test_wirsum_1997_named_f_p2.py \
  tests/test_wirsum_1997_named_f_p0.py \
  tests/test_wirsum_1997_named_f_e8.py \
  -k "not test_e8_saved_endpoint and not test_e8_solve" -q

# E7 是否仍认证
PYTHONPATH=. .venv/bin/python - <<'PY'
import numpy as np
from src.solvers.wirsum_1997_named_f_endpoint_io import load_named_e8_snapshot, recertify_named_e7
x = np.asarray(load_named_e8_snapshot()["e7"]["x"])
print(recertify_named_e7(x))
PY
```

求根时用快照 E7 \(x\)，不要从 R11 重走八分钟链，除非 E7 现场认证失败。

---

## 8. 明确不要做的事

- 不要启动 P3 / P4
- 不要更新旋风 \(K\)、返料腿、生产 `hamel_outer_refresh`
- 不要改 `rate_*` 来让根移动
- 不要把 \(s<1\) 的混合残差当成物理模型或认证
- 不要把未认证的 E7 \(x\) 标成 `endpoint_certified` 的 E8
- 不要 commit（除非用户明确要求）；若提交，只加具名核文件，不要扫整个脏工作树
- 不要宣称炉温 / HTW / 全炉 LU
