# Wirsum (1997/1998) — NR / 阻尼 Newton 摘录笔记

**文献：** Manfred Wirsum, *Theoretische Beschreibung und experimentelle Untersuchung der Verbrennung von Klärschlamm in blasenbildenden Wirbelschichtfeuerungen*, Dissertation Universität-GH Siegen, 1997；Fortschr.-Ber. VDI Reihe 6 Nr. **383**, VDI Verlag, Düsseldorf, 1998.

**本地 PDF（gitignore `*.pdf`）：**  
`docs/_1997_VDI-Dis._Wirsum_Klaerschlamm_blasenbildende_Wirbelschichtfeuerungen.pdf`

**OCR 原文摘录：** `data/wirsum_1997_extract/`（`page_020.txt`–`page_021.txt` 等）

**证据等级：** §2.3.3 求解器段落 **E1**（已通读本地 PDF 对应页）；公式个别符号 OCR 模糊处单独标注。

---

## 1. 与 Hamel 的引用关系

Hamel (2001) §2.3 写明：远初值用 **Wirsum 的 gedämpfte Newton**；更新框架见 Hamel Eq. (2.9)。  
本书目录 **§2.3.3「Berechnungsverfahren zur Lösung der Zellbilanzgleichungen」**（印页约 21）即该细节来源。

本书主题是 **鼓泡床污泥燃烧**（非气化），但 **cell 联立 + 块三对角 NR + 阻尼** 与 Hamel 共用 Siegen 求解语法。

---

## 2. §2.3.3 算法要点（已核实）

### 2.1 问题形式

- 稳态 cell 物质/能量积分平衡 → **非线性代数方程组** \(F(x)=0\)。
- 求解基于 **Wozny [149]** 的 Newton–Raphson，并适配 **块三对角（blocktridiagonal）** 雅可比结构。
- 单 cell 方程 \(F_i\) 通常只依赖本格与邻格未知量；例外：非邻格进料/能量（如 CFB、非邻格颗粒对流传热）。

### 2.2 为何要阻尼

原文大意（PDF p.21）：

- 初值 \(x_0\) 靠近解 \(x^*\) 时，Newton 收敛快；
- **远离解时常发散**；
- 因此采用 **「gedämpfte」Newton 变体**。

### 2.3 接受准则（核心，OCR 清晰）

新迭代 **仅当** 残差向量的 **欧几里得范数严格下降** 才接受：

\[
\|F(x^{(\nu+1)})\|_2 < \|F(x^{(\nu)})\|_2
\qquad \text{(Wirsum Eq. 2.23)}
\]

### 2.3b Newton 线性步（OCR 已清晰）

\[
F'(x^{(0)})\,\Delta x = -F(x^{(0)})
\qquad \text{(Wirsum Eq. 2.21)}
\]

（与标准 Newton / Hamel 更新框架同族。）

**没有**出现：

- `merit = max(rms, split_rms, energy_rms)`
- 单独的 energy-descent / split 门控
- CO/syngas inventory collapse guard
- fastox 投影作为线搜索步骤

→ 本仓库上述机制一律保持 **工程桥** 标注（与谱系文档既有结论一致，现有原文背书）。

### 2.4 迭代流程（原文步骤 1–6）

1. 选全体未知量初值 \(x_0\)，迭代计数 \(\nu=0\)。  
2. 组装 \(F(x^{(\nu)})\) 与雅可比 \(F'(x^{(\nu)})\)。  
3. 解线性系 (2.21) 求 \(\Delta x^{(\nu+1)}\)。  
4. 求阻尼指标 \(j\)，使试探点满足 **(2.23)**。OCR 将步长因子印成残缺符号，结合「gedämpft / 二分回溯」语境与本仓库既有实现，**最合理解读**为对满步 \(\Delta x\) 做 **\(2^{-i}\) 回溯**（\(i=0,1,2,\ldots\)），取最小使 \(\|F\|_2\) 下降的 \(i\) 作为 \(j\)。  
5. \(x^{(\nu+1)} = x^{(\nu)} + 2^{-j}\Delta x\)（因子依解读；接受准则以 (2.23) 为准）。  
6. 回到步骤 2。

### 2.5 终止条件（原文列举）

以 \(F(x)=0\) 的实用替代：

1. 达到最大迭代次数；  
2. \(\|x^{(\nu+1)}-x^{(\nu)}\| < \delta\)；  
3. \(\|F(x^{(\nu+1)})\| < \varepsilon\)。

---

## 3. 与本仓库实现的对照

| 项目 | Wirsum §2.3.3 | 本仓库现状 | 判定 |
|------|---------------|------------|------|
| 全局 Newton + 块结构 | 有（Wozny 系） | `solve_global_nr` + 结构化 J | **同族** |
| 阻尼 / 回溯 | \(\|F\|_2\) 下降 +（解读为）\(2^{-i}\) | `wirsum_damping.py`：`n_damp_halvings`、`lambda` 二分 | **骨架同构** |
| 线搜索标量 | **仅** \(\|F\|_2\)（或等价未加权欧氏范数） | **默认**缩放残差 RMS（‖F̂‖₂ 族）；`max(rms,split,energy)` 仅 opt-in | **默认同构**；max-merit = 工程桥 |
| fastox / syngas guard / LSQ preinner | 文中未见 | 多项 thesis 开关（**默认关**；`HAMEL_WIRSUM_CORE`） | **工程桥** |
| Startwert / Vorabrechnung | 文中强调远初值需阻尼；细则在燃烧章 | `vorabrechnung/` | 需另页对照，勿混称「Wirsum 线搜索」 |

代码锚点：

- `src/solvers/wirsum_1997_strict.py` — **Eq. 2.21–2.23 最小实现**（`raw_l2` + `α=1,1/2,…`；无 fence/fastox/split merit）
- `ReactorConfig.nr_solver_mode='wirsum_1997_strict'` — 关掉可由 config 开关的工程桥，并强制 `holdup_k0_anchor='none'`
- `src/solvers/wirsum_damping.py` — 阻尼旋钮（应对齐「二分回溯 + 残差范数下降」叙事）  
- `src/solvers/global_nr_solver.py` — `_line_search_merit`（默认 RMS；max-merit = opt-in 工程扩展）  
- `src/solvers/global_nr_iteration.py` — LS 试探 / fastox / syngas guard

生产 `solve_global_nr` 在 strict 模式下仍用缩放 RMS 接受步；原文 `||F||₂` 只保证在 `wirsum_1997_strict_newton` 上。

具名核打开 R12/R5/R6/R10 时，每一档 \(s\) 的内层仍是上述 Eq. 2.21–2.23；\(s\) 网格本身不是 Wirsum 原文，见 `docs/wirsum_1997_named_f_rate_continuation.md`。

---

## 4. 对当前 NR 主线的含义（H₂ / 平台）

1. **不能**再把 `merit=max(...)`、energy/split 对抗、syngas guard 说成「Wirsum/Hamel 同构」。  
2. **文献对齐已落地为默认**：`nr_line_search_gas_phase_split_merit_thesis=False`、`nr_line_search_energy_merit_thesis=False`；接受 = **缩放后残差 \(\ell_2\)/RMS 下降**（原文未规定无量纲化细节，缩放保留）。  
3. 这 **不自动修复** bed0 O₂–H₂ 刚性（§5.41）；阻尼准则只约束「何时接受 Newton 试探」，不改反应刚度。  
4. 上段 H₂「可强制闭合但 NR 未走到」仍可能与旧工程 merit 景观有关——在新默认下做平台短对照（CO/rms/y_H₂），再决定是否动 fastox/syngas。

---

## 5. 尚缺 / 下一步文献动作

| 项 | 状态 |
|----|------|
| Eq. (2.21) \(F'\Delta x=-F\) | **已确认**（OCR 清晰） |
| \(2^{-i}\) 步长是否原文显式写出 | **解读**；接受准则 (2.23) 已硬核实 |
| Wozny (1983) Habilitation | 仍缺全文；块消元细节仍阻塞 |
| 本书水力 / 颗粒轨迹章 | 与 freeboard Wirsum 轨迹模型可能相关，未在本笔记展开 |
| handoff 全表核查 | 见 `docs/handoff_nr_convergence_hamel_isomorph_2026-08.md` **§6.6** |

**建议下一刀（文献侧）：** 人工目视校对 PDF 印页 §2.3.3 公式 (2.21)–(2.23)，把步长因子从「解读」升为 E1；可选跑「仅 \(\|F\|_2\) merit」短 A/B。
