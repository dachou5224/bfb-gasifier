# 交接文档：反应源项 Hamel 同构 × Phase2 NR 收敛平台期

**日期：** 2026-08-06  
**状态：** 进行中 / 未 fully converge  
**接手对象：** 下一任 agent（数值收敛 + 可选文献求证）  
**语言：** 交接说明用中文；代码/配置键名保持英文。

---

## 0. 30 秒结论

床层 extent 限幅已关，快氧化 Startwert + LS 投影 + 窄 bed1 Pre-Inner 已把残差从天文量级压到 **split≈0.037 平台**，但 `converged_fully=False`。  
继续扫 Wirsum 阻尼旋钮**无效**。平台期根因是线搜索 `merit=max(rms,split,energy)` 下的 **energy↓ / split↑ 对抗**。  
Hamel 正文**无**该对抗的处方；算法细则锁在 **Wirsum (1998) / Wozny (1983) / 1999 Final Report** 等**尚未到手全文**的文献里。  
下一刀应优先：**平台期 cell/物种级 split 主导项**，或 **refresh-aware / split 约束的接受规则 A/B**——不要扩大 LSQ，不要恢复床层限幅当物理修复。

---

## 1. 任务边界与硬约束

### 1.1 主线目标

- 推进「反应源项 Hamel 同构（去 extent 限幅）」下的 HTW LU Phase2 全局 NR 收敛。
- **勿改**用户计划 `.plan.md`（若存在）。
- 床层代数：**不恢复** extent 限幅；可用 Wirsum 阻尼 / 初值 / 自由板 dt 界；自由板与床层政策可分离。

### 1.2 项目硬约束（必读）

| 文档 | 用途 |
|------|------|
| `docs/CLAUDE.md` | Phase 顺序、SI、Arrhenius 工厂、禁止事项 |
| `.cursor/rules/claude-agent-instructions.mdc` | 同上快速参考 |
| `docs/vorab_init_abgleich_vs_outer_preinner.md` | Init / Outer Pre-Inner / Check1–2 语义 |
| `docs/hamel_engineering_bridge_notes.md` | 已知工程桥与错误路线 |
| `docs/hamel_isomorphic_gap_implementation_checklist.md` | 同构台账 |
| `docs/siegen_cell_model_lineage.md` | Siegen 谱系 + 已读/未读文献 |

技能（按需）：

- `.cursor/skills/bfb-diagnosing-bugs/SKILL.md` — 先建反馈环再假设  
- `.cursor/skills/bfb-numerical-guard/SKILL.md` — 求解器/Jacobian  
- `.cursor/skills/bfb-hamel-parity-audit/SKILL.md` — 方程同构争议  

### 1.3 复现入口（Phase2 LU）

```bash
cd /Users/liuzhen/AI-projects/bfb-gasifier
# 典型 harness：scripts/_lu_harness.py
# ladder_level='p2'，enable_vorab_march_staged=True
# 常用 cfg（同构路径）：
#   extent_limiters_enabled_thesis=False
#   freeboard_extent_dt_bounds_thesis=True
#   vorab_transport_x0_fast_oxidation_closure_thesis=True
#   nr_line_search_fast_oxidation_projection_thesis=True
#   enable_r12=True
#   bed0 T seed ≈ 1250 K（cold-anchor 关）
```

验收指标（至少记录）：

- `rms_scaled_final` / `rms_scaled_gas_phase_split_final` / `rms_scaled_energy_final`
- `converged_fully`、`nr_outer_iters`、`nr_accepted_lambda_history`
- **注意：** `energy_final` 常是 best-so-far；晚期活状态 energy 可能 ~0.040，与报告的 ~0.013 不一致。交接后请同时看 `nr_outer_history` 末几项。

---

## 2. 已落地变更（勿回退除非证据推翻）

| 项 | 位置 / 开关 | 作用 |
|----|-------------|------|
| 床层限幅关 | `ReactorConfig.extent_limiters_enabled_thesis=False` | 同构：裸动力学进 \(F(x)\) |
| 自由板 dt 界 | `freeboard_extent_dt_bounds_thesis=True` | 与床层解耦，防自由板炸温 |
| 快氧化 Startwert | `vorab_x0.py`；`vorab_transport_x0_fast_oxidation_closure_thesis` | R12→R5→R6→R10；合计 holdup 氧化再按水力重分；在 init abgleich **之后**调用 |
| NR LS 快氧化投影 | `nr_line_search_fast_oxidation_projection_thesis` | 防 Newton 试探把 H₂+O₂ 叠回 → 裸 R12 天文残差 |
| bed1 Pre-Inner | Phase2：`nr_outer_preinner_bed1_phase_split_thesis`，`max_bed_index≈2` | 压 split；**非 Hamel 处方** |
| 全局 split 守卫 | `preprojection.py` | 单格 LSQ 若抬高全局 split → 回滚 |

相关测试：`tests/test_vorabrechnung_hamel_x0_transport.py`（含 R10）。

---

## 3. 收敛轨迹（摘要）

| 阶段 | o8 量级 | 备注 |
|------|---------|------|
| 仅去限幅 / 裸 R12 | 天文 / stall | — |
| +快氧化 x₀ | rms~0.72→0.14 | — |
| +R10+重分 | split~0.07，energy~0.076 | — |
| +LS fastox | split~0.05–0.07，energy~0.027 | — |
| +bed1 preinner（现默认） | **split~0.037**，energy(best)~**0.013** | **平台；未 fully converge** |

---

## 4. 已做且结论明确的数值实验

### 4.1 阻尼作为「替代 LSQ」主路径 — 失败

证据：

- `data/phase2_extent_raw_damping_vs_lsq_ab_o8.json`
- `data/phase2_extent_raw_damping_step_model_ab_o8.json`

结论：关 LSQ 后 `prefer_full_step=False`、`λ₀=0.25/0.125`、更多 halvings、关 split-merit、`lm`/`ptc` → **末态锁在同一点**（split≈0.0376，energy≈0.070）。线搜索一直在跑，但**从未作为替代 LSQ 系统有效过**；现已证实无效。

### 4.2 线搜索接受轨迹 — 平台期几乎全失败

证据：

- `data/phase2_extent_raw_ls_accept_probe_o8.json`
- `data/phase2_extent_raw_ls_fail_root_o8.json`

| 配置 | LS 失败率 | 含义 |
|------|-----------|------|
| baseline LSQ on | ~20/32 | 晚期 outer 连续失败 |
| LSQ off | **20/20** | 零接受步 |
| LSQ on + soft damp | ~20/32 | 略降 split、抬高 energy |

代码路径：`src/solvers/global_nr_iteration.py`（线搜索循环）、`src/solvers/global_nr_solver.py`（`_line_search_merit` / `_line_search_trial_acceptable`）。

### 4.3 merit 分量拆解 — energy↔split 对抗

证据：`data/phase2_extent_raw_ls_merit_decomp_o8.json`  
插桩方式：monkey-patch `_line_search_trial_acceptable`（勿当生产代码）。

平台期典型（准确钩子）：

| | merit | split | energy | 主导 |
|--|------:|------:|-------:|------|
| 当前 | ~0.040 | ~0.037 | **~0.040** | energy |
| 最佳试探 | **~0.045** | **~0.045** | ~0.022 | split |
| Δ | + | **+** | − | — |

拒绝原因：**100% `merit_not_decreased`**（不是 split-inflation 守卫）。  
试探总 rms≈0.017（更好）仍被拒。

### 4.4 「energy 降即可接受」— 有害

证据：`data/phase2_extent_raw_energy_desc_accept_ab_o8.json`  
`nr_refresh_hydrodynamics_on_accepted_step_thesis=True`（水力刷新已开）。

| 策略 | split | energy |
|------|------:|-------:|
| baseline | **0.037** | **0.013** |
| energy 降就收 | 0.064 | 0.027 |
| energy 降 + split 软帽 | 0.047 | 0.039 |

结论：用户直觉（energy 更重要、下一步会重算 hydro）**机制上部分成立**，但裸/软门控接受会复现「energy 小降、split 爬升」病理；水力刷新**不能**自动修好 split。

---

## 5. 建议的下一刀（接手优先序）

### P0 — 证据收紧（推荐先做）

1. 平台期失败试探：按 **cell × 物种** 拆 split 残差主导项（bed0/1/2、主组分）。  
2. 固定报告模板：每个 outer 同时记 `outer_merit_composite`、`outer_energy_rms_final`、`gas_phase_split_rms`、`outer_merit_best_so_far`，避免 best-so-far 误导。

### P1 — 数值策略 A/B（小步、有守卫）

1. **Refresh-aware 联合接受**：energy 主导且 energy↓ 时，仅当「试探 + 水力刷新」后的 split/merit 不劣于当前才接受。  
2. **Energy 主导时的 split 约束步**：冻结/限制 dense–bubble holdup 扰动，只动能量相关 DOF。  
3. **分阶段 merit**（先 energy 后 split）— 明确标为工程桥，非 Hamel。  
4. 可选：有条件 `allow_gd_fallback` 对照（默认 Hamel 路径 LS 失败后常直接停）。

### P2 — 明确不要做的事

- ❌ 恢复床层 `extent_limiters` 当「物理修复」  
- ❌ 再扫 `nr_lambda_init` / halvings / prefer_full_step 当主线  
- ❌ 扩大 `max_bed_index` LSQ（曾把 energy 抬到 ~0.068）  
- ❌ 在缺 Wirsum/Wozny 全文前，把 max-merit / energy-only accept 写成「Hamel 同构」

### 同构 vs 工程桥（措辞纪律）

| 较稳同构 | 工程桥（须标注） |
|----------|------------------|
| 守恒式、水力链、动力学源项结构 | phase-split LSQ Pre-Inner |
| 床层去限幅 | `merit=max(rms,split,energy)` |
| 快氧化 Startwert 思路（计量投影） | split-inflation 守卫、energy-descent 门控 |
| Eq.2.9 阻尼框架（有旋钮） | GD fallback、PTC/LM 步型 |

---

## 6. 文献求证状态（给接手 agent）

### 6.1 Hamel 正文已确认的内容

- §2.3：Wozny NR、块三对角 + Nebenelemente、Eq. 2.9、远初值用 Wirsum **gedämpfte** Newton；**细节见 Wirsum (1998)**。  
- **没有** energy vs phase-split、分量 merit、单支接受规则。  
- 本地 Diss. PDF：`docs/_2001_VDI-Dis._…Brennstoffe.pdf`（见 `.gitignore`）。

### 6.2 已读全文（E1）

| 文献 | 对 NR 的用途 |
|------|----------------|
| Rajan & Wen, *AIChE J.* 26(4), 1980, 642–655 | 细胞/两相**骨架**；数值是分室+elutriation 迭代，**非**全局 NR |
| Wozny et al., *CIT* 59(2), 1987, 170–171 | 2 页短讯；NR+Jacobian 示意；细则回指 **1983 Habilitation** |
| Hamel & Krumm, *Powder Technology* 120, 2001, 105–112 | 验证短文；求解过程指 **1999 Final Report [10]** |
| Weil/Hasselbach/Hamel/Krumm, ISCC Nanjing, 2003 | 平衡 vs 细胞对比；细节回指 Hamel Diss. |
| 锡根官网 Publikationen 清单 PDF | 检索索引；无早期 Diss.；无 NR 细则 |

用户本地路径（仓库外，示例）：

- `…/Documents/AIChE Journal - July 1980 - Rajan - ….pdf`
- `…/Documents/Chemie Ingenieur Technik - February 1987 - Wozny - ….pdf`
- `…/BFB-cell-model/reference/S. Hamel & W. Krumm_2001_….pdf`
- `…/文献/德国siegon大学/1_ConferencePaper_2003_….pdf`
- `…/文献/德国siegon大学/1_mb.uni-siegen.de-Publikationen….pdf`

### 6.3 关键但尚未拿到全文（阻塞「求解器同构」宣称）

| 文献 | 为何关键 | 阻碍 |
|------|----------|------|
| **Wirsum (1998)** VDI-FB Reihe 6 Nr. 383 | Hamel **明文**阻尼 Newton 细节 | Diss./VDI，难公开获取 |
| **Wozny (1983)** Habilitation, Siegen | 块三对角 NR 框架源头（主题是**板式塔**） | Habilitation 流通差 |
| **Hartleben (1983)** Siegen Diss. | 细胞模型 + Wozny 解法接到鼓泡床 | Diss. |
| **Hamel/Fett/Krumm (1999) Final Report** Project IV A4-9021.1 | PT 文称含 mathematical solution procedure | **州科研结题报告**（非企业商务服务，但是非期刊、难检索） |
| Weiß / Bürkle / Dersch / Heinbockel 等 | 同系可能有实现叙述 | 多为 Diss./VDI-FB |

谱系与读后摘要：`docs/siegen_cell_model_lineage.md`（持续更新）。

### 6.4 文献优先级（若继续求证）

1. Wirsum 1998  
2. Wozny 1983（只筛 NR/阻尼/接受规则，忽略精馏物性）  
3. 1999 Final Report IV A4-9021.1  
4. Hartleben 1983（桥接叙述）  
5. 清单中的 INFUB 2002 *Simulation of Fluidised Bed Gasifiers*（实现叙述概率中等）

---

## 7. 关键代码锚点

| 主题 | 路径 |
|------|------|
| 全局 NR 迭代 / LS | `src/solvers/global_nr_iteration.py` |
| merit / accept | `src/solvers/global_nr_solver.py` (`_line_search_merit`, `_line_search_trial_acceptable`) |
| Wirsum 旋钮 | `src/solvers/wirsum_damping.py` |
| 快氧化 x₀ | `src/solvers/vorabrechnung/vorab_x0.py` |
| Outer / Pre-Inner | `src/solvers/outer_loop.py`, `src/core/connectivity/joint_gate.py` |
| Config 默认 | `src/core/reactor.py` (`nr_*_thesis`, `extent_limiters_*`) |
| Phase2 配置 | `tests/validation/config_phase2_freeboard.py`, `scripts/_lu_harness.py` |

---

## 8. 接手检查清单

- [ ] 读 `docs/CLAUDE.md` + 本交接文档 §0–§5  
- [ ] 跑一条 baseline o8，确认 split≈0.037、未 fully converge（与 §3 同量级）  
- [ ] 浏览 `data/phase2_extent_raw_ls_merit_decomp_o8.json` 与 `*_energy_desc_accept_ab_o8.json`  
- [ ] 选定 P0（物种主导）或 P1（refresh-aware accept）之一，先建可重复 JSON 证据再改生产默认  
- [ ] 任何「Hamel 同构」措辞前核对 §5 / §6.3  

---

## 9. 一句话交接

**同构源项与快氧化初值已到位；收敛卡在 max-merit 下的 energy–split 交易；阻尼与裸 energy 接受已证伪；细则文献未到手——用窄 A/B + 物种级证据推进，勿回退限幅、勿扩 LSQ。**
