---
name: bfb-hamel-parity-audit
description: Hamel 同构 parity 审计。按方程号逐条对照 docs/hamel_submodels、Hamel 原论文与 src/ 实现，产出差异表与证据链。在重构子模型、宣称「与 Hamel 一致」、或 gap checklist 更新前使用。
---

# bfb-hamel-parity-audit（Hamel 同构差异审计）

## 与 bfb-scientific-audit 的分工

| 技能 | 粒度 | 典型时机 |
|------|------|----------|
| **bfb-scientific-audit** | 改动级门禁：SI、Arrhenius 工厂、单位注释、引用是否存在 | 合并前 / `/review` 前 |
| **本技能** | 方程级 parity：论文式 ↔ 转录 SoT ↔ 代码式，逐条判定 | 重构子模型、闭合 gap、争议裁决 |

两者**都遵守同一证据链**，本技能负责**深度、可归档**的差异表，不替代 `tests/sanity_checks.py` 或 Fig.7 数值验证。

## 何时启用

- 重构 **流体力学 / 动力学 / 守恒 / NR 结构** 后，需证明与 Hamel 对齐或标注已知偏离。
- 更新 **`docs/hamel_isomorphic_gap_implementation_checklist.md`** 或 `specs/hamel_checked_algorithm_spec.md` §7 前。
- 某条方程（如 Eq.3.50、Eq.5.36）出现 **代码 vs 文档 vs 论文** 三方不一致。
- **禁止**在未完成本审计（或等价人工表格）前，以模拟结果为由进入调参（见 `docs/CLAUDE.md` §一致性优先）。

## 证据链（不可违反）

顺序固定，**禁止反向**：

1. **转录 SoT**：`docs/hamel_submodels/00–07`（页锚点见 `00_readme_and_citation_rules.md`）
2. **已核验摘要**（若存在）：`specs/hamel_checked_algorithm_spec.md` 方程清单与代码映射
3. **Hamel 原论文 PDF**（本地-only）：仅当转录与代码冲突，或转录本身待核时
4. **代码**：`src/` 对应实现与 docstring
5. **已知偏离台账**：`docs/hamel_isomorphic_gap_implementation_checklist.md` §0.5

**禁止**用当前实现反向补文档或把运行结果当文献依据。证据方向只能是 **论文 → 文档 → 代码**。

## 审计范围选择

启动时先与用户确认 **scope**（一次只审一个 scope，避免泛泛而谈）：

| Scope 代号 | 转录分册 | 主要 `src/` 模块 |
|------------|----------|------------------|
| `hydro-core` | `03_hydrodynamics_core_chain.md` | `physics/minimum_fluidization.py`, `bubble_dynamics.py`, `phase_fractions.py`, `mass_transfer.py` |
| `freeboard` | `04_freeboard_particle_trajectory_analytic.md` | `physics/freeboard.py` |
| `exchange-cell` | `02_cell_balances_and_exchange.md` | `core/cell.py`, `physics/exchange*.py`（若存在） |
| `global-nr` | `01_global_nr_and_connectivity.md` | `solvers/global_nr_solver.py`, `core/reactor.py` |
| `drying-pyro` | `05_drying_and_daem_pyrolysis.md` | `physics/drying.py`, `devolatilization.py` |
| `kinetics` | `06_kinetics_r1_r11_and_equilibrium_driving.md` | `kinetics/gas_reactions.py`, `char_reactions.py`, `tar_reactions.py`, `arrhenius.py`, `thermodynamics/equilibrium.py` |
| `a1-gibbs` | `07_appendix_a1_gibbs_solver.md` | `thermodynamics/gibbs*.py`（若存在） |
| `single-eq` | 用户指定方程号 | 由分册索引定位 |

若 scope 未指定，**先问用户**审哪一分册或哪条方程，再开始读文件。

## 执行步骤

### Step 0 — 跑自动化锚点（若 scope 相关）

在对应 scope 下，先执行（通过则记录，不通过则优先修）：

```bash
cd bfb-gasifier && python3 scripts/audit_hamel_consistency.py
```

该脚本检查：关键文件是否引用 `docs/hamel_submodels`、是否残留 stale `specs/0x`、以及 Eq.3.44/3.50、R5/R6/R8 等**数值锚点**。它**不是**全方程覆盖，只是快速 sanity。

### Step 1 — 列出待审方程清单

从选定分册中提取：

- 方程编号（如 Eq.3.41–3.50 链）
- 表格编号（如 Table 5.3、5.4）
- 分册内已标注的 **「部分一致 / 开放 / B 档」** 条目

输出简短清单（编号 + 分册内小节标题），供 Step 2 逐条勾选。

### Step 2 — 逐条 parity 对照（核心）

对清单中**每一条**，按相同模板填一行（不得跳过「无实现」项——须标为 `缺失`）：

1. 读分册中该方程的**完整 LaTeX/文字转录**（含符号定义与单位）
2. 读 `specs/hamel_checked_algorithm_spec.md`（若本地存在）中对应映射行
3. 在 `src/` 定位实现函数；读 docstring 中的 `Source:` / `Hamel (1999) Eq.xxx`
4. 对照：**符号**、**指数**、**单位**、**相（bubble/suspension/dense）**、**压力项（Pa vs atm，R8）**
5. 若分册与代码不一致 → 查 PDF 页锚（`00_readme_and_citation_rules.md`）裁决
6. 若与 gap checklist §0.5 已有记录一致 → 引用该行，不重复发明新结论

**Arrhenius 专项**：速率常数必须经 `src/kinetics/arrhenius.py` 工厂；手写 `exp(-E/RT)` 直接判 `不一致`。

**R8 专项**：压力项 atm 换算须在函数内部完成，外部接口 Pa（见 `docs/CLAUDE.md`）。

### Step 3 — 歧义处理（限定式「grilling」）

仅当 **论文/转录对同一方程存在多义解读** 且无法单靠读文裁决时：

- **一次只问一个**澄清问题（符号含义、相归属、Table 5.3 行对应关系等）
- 每个答案必须指向 **PDF 页码 / 表格行 / 分册段落**，不得「凭实现推测」
- 裁决结果写入差异表「备注」，必要时建议新增 ADR 或更新 gap checklist

**不要**用开放式 grilling 替代 Step 2 的逐条对照。

### Step 4 — 数值 spot-check（可选但推荐）

对 Step 2 中判为 `一致` 或 `部分一致` 的关键速率/系数：

- 优先复用 `scripts/audit_hamel_consistency.py` 已有用例
- 或写**最小** Python 片段：固定 `T, P, C_*`，对比 `model` vs `hamel_*_reference` / 分册给定数值
- 数值不一致 → 降级为 `不一致（数值）` 并引用计算输出

### Step 5 — 汇总与台账更新建议

交付物必须包含：

1. **Parity 差异表**（见下方模板）
2. **统计**：一致 / 部分一致 / 不一致 / 缺失 / 有意偏离 各几条
3. **建议动作**：修代码、修分册、更新 gap checklist、或进入 B 档工程补缝（`hamel_engineering_bridge_notes.md`）
4. 若存在 `不一致` 且未入台账 → 明确建议更新 `hamel_isomorphic_gap_implementation_checklist.md` §0.5

## 输出模板（必须使用）

```markdown
## Hamel Parity 审计 — {scope} — {日期}

**证据链**：hamel_submodels → checked_spec → PDF（若用）→ src/

| 方程/表 | 分册位置 | 论文式（摘要） | 代码位置 | 判定 | 证据 |
|---------|----------|----------------|----------|------|------|
| Eq.3.50 | 03 §… | K_bd = … | mass_transfer.py:calc_kbd | 一致 | audit_hamel_consistency PASS |
| Eq.5.36 R5d | 06 §… | …·T^-0.5·C_H2O^0.5 | gas_reactions.py:rate_R5_suspension | 部分一致 | T 指数一致；… |
| … | | | | | |

**判定枚举**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离（见 gap §x）`

**自动化**：audit_hamel_consistency.py → {PASS/FAIL/未跑}

**开放项**：…

**建议下一步**：…
```

## 完成标准（completion criterion）

本技能结束当且仅当：

- [ ] scope 内清单**每一条**都有表格行（含 `缺失`）
- [ ] 所有 `不一致` / `部分一致` 均有**文件:行或分册段落**证据
- [ ] 未用代码反向改分册；文档变更建议单独列出并标注需 PDF 依据
- [ ] 用户已确认 scope 内无遗漏方程（或明确记录「本次跳过」项及原因）

## 相关文件

- `docs/README.md` — 文档分层索引
- `docs/hamel_submodels/00_readme_and_citation_rules.md` — 页锚点与引用规则
- `docs/hamel_isomorphic_gap_implementation_checklist.md` — 同构台账（本地-only）
- `specs/hamel_checked_algorithm_spec.md` — 已核验映射（本地-only）
- `docs/hamel_engineering_bridge_notes.md` — 论文未给细节的 B 档补缝
- `scripts/audit_hamel_consistency.py` — 自动化锚点
- `.cursor/skills/bfb-scientific-audit/SKILL.md` — 改动级门禁（审计后合并前仍须跑）
- `.cursor/skills/bfb-phase-manager/SKILL.md` — Phase 顺序；parity 未闭合不得跳至调参
- `.cursor/skills/bfb-diagnosing-bugs/SKILL.md` — 数值/预算症状调试；parity 与调参门禁见该 skill §调参门禁
