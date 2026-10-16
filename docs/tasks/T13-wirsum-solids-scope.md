# T13：Wirsum（1998）限定查阅——床区固体输运、床料存量与惰性床料校正

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，见“执行要求”） |
| 审查方 | codex（原文答卷，对照原页逐条审查） |
| 前置任务 | T12（答卷 `refs/hamel/answers/T12-solids-flow-scope.md`）；ADR 0007（已批准）；目标所有者 2026-10-16 批准限定查阅 Wirsum（见 `docs/STATUS.md`“已决定”） |
| 分支 | `task/T13-wirsum-solids-scope` |
| 时间上限 | 2 个工作日 |
| 状态 | 已批准 |

## 背景

T12 确认：Hamel 原文没有给出床区整体向上、向下固体流的计算式，也没有说明 M_i 包含哪些固体、“Abgleich der Inertmaterialbilanz”校正什么；这几处原文都没有引用文献。ADR 0007 的闭合式因此缺少 Ṁ_auf,i、Ṁ_ab,i 和 M_i。

Wirsum（1998）是 Hamel 的直接前身。Hamel 明确沿用了他的阻尼牛顿法（S. 20）、式 (2.8) 的能量流写法（S. 19）、粒级数的扩展（S. 15）和自由空间阻力式（S. 45），但**没有**在床区固体输运上引用他。按 CHARTER §4.3，目标所有者批准在下述限定范围内查阅 Wirsum。

## 来源与范围限定

- 文件：本地 `docs/_1997_VDI-Dis._Wirsum_Klaerschlamm_blasenbildende_Wirbelschichtfeuerungen.pdf`（Fortschr.-Ber. VDI Reihe 6 Nr. 383，Hamel 引作 Wirsum 1998；SHA256 前缀 e57844f81f1b；115 个 PDF 页，多为双页拼版）。**不入库**（CHARTER §11）。
- **只查以下内容**：
  1. 单元模型中固体的衡算变量（§2.3.1 „Wahl der bilanzierten Variablen“）与固体流衡算（§2.3.2.2），以及模型结构中预计算与单元模型的分工（§2.2）；
  2. 床区相邻 cell 之间向上、向下固体流的计算，包括 §5.6 以 GR 模型为基础的分离与固体循环（尤其 §5.6.2.1 „Aufspaltung der Feststoff-Zirkulationsströme“、§5.6.3 模型方程），以及 §3.2 中与之直接相关的流动量；
  3. 床料存量：床高、空隙率、床料质量之间的关系；惰性床料如何衡算或校正（全书检索 Inertmaterial、Bettmaterial、Abgleich）；
  4. 床面向自由空间的固体抛出与床区最上层 cell 向上流的关系（§3.4.1），床底 cell 的边界，以及排灰流的处理；
  5. 悬浮相中“催化剂质量”或“床料质量”作为反应基准时的口径（若 Wirsum 有，ADR 0010 第 9 条待定项需要）。
- **不查**：污泥的干燥、热解、燃烧、污染物生成和实验部分（第 1 章、第 4 章、§5.1–5.5、第 6 章及以后），除非上述五项的式子明确引用它们；不引入 Wirsum 的任何反应动力学或物性关联式。
- 页码线索（按文本层检索所得，仅供定位，以页面图像为准）：§2.3 约在 PDF 16–20；“Abgleich”见 PDF 17；“Inertmaterial”见 PDF 16、78；§3.4.1 约在 PDF 29；§5.6 约在 PDF 48–53。

## 要回答的问题

每个问题都要回答：Wirsum 怎么写（PDF 页、半页、印刷页码、式号、德文原句）；没有写什么；**与 Hamel 原文是否相容**，逐条对照 Hamel 的 S. 12、表 2.2（S. 13）、Bild 2.2（S. 14）、式 (2.1)（S. 15）、式 (2.6)（S. 17–18）。

1. **固体衡算变量**：Wirsum 单元模型中固体的未知量是粒级存量、质量流，还是其它量？若与 ADR 0007 的选择（粒级存量）不同，写明差别及影响。
2. **床区整体上下固体流**：Wirsum 如何计算相邻 cell 之间的向上、向下固体流？是否就是 §5.6 GR 模型中的固体循环流？计算这些流需要哪些输入（气泡量、尾涡份额、分离系数等），这些输入在 Hamel 中是否都有对应量。若 GR 模型的循环流只用于污泥（flotsam/jetsam 分离），写明它能否用于惰性床料与焦炭。
3. **M_i 与惰性床料**：Wirsum 的 cell 固体存量包含哪些固体？床高、空隙率与床料质量的关系式是什么？Wirsum 是否有与 Hamel“Abgleich der Inertmaterialbilanz”对应的步骤，校正的量与条件是什么？
4. **两端与外排**：床区最上层 cell 的向上流与床面抛出（§3.4.1）的关系；床底 cell 的向下流；排灰流的位置与数值来源。
5. **预计算与单元模型的分工**：哪些固体量在预计算中给定、牛顿法内层不变，哪些由单元模型求解。
6. **催化剂质量口径**（若有）：悬浮相反应以床料质量为基准时，Wirsum 用哪一部分固体的质量。

对每一个找到的式子，另写一行“能否用于 Hamel 主模型”的判断依据：Hamel 是否在任何地方提到同一机制或同名量；不能用时写明原因。**是否采用由总司令在 ADR 0009 中裁定**，本卡不作采用结论。

## 执行要求（给 cursor）

- **新会话**：开一个新的 cursor 会话，不要沿用做过 T4、T8、T9、T11、T12 或 worktree 清理的会话。
- **新 worktree**：
  ```
  git worktree add ../bfb-T13 -b task/T13-wirsum-solids-scope origin/main
  ```
  按 `.claude/skills/hamel-source-check/SKILL.md` 第一节第 2 项链接 Hamel 页面图像（对照 Hamel 原文时用）。
- **Wirsum 页面图像**：用 PyMuPDF 把上述 PDF 渲染到主工作树的 `refs/wirsum/pages/`（每页一张 PNG 和一份文本层 TXT，三位数页号），在 worktree 中按 Hamel 页面的方式建符号链接。本卡在 `.gitignore` 中增加一行 `refs/wirsum/pages`（不带斜杠，同 Hamel 的写法），渲染脚本若写成文件，放在本卡 PR 中但不得处理其它 PDF。文本层只用于定位，引用的式子和德文原句必须对照页面图像。
- **开工前必读**：`docs/CHARTER.md`（§4 来源层级）、本卡、`refs/hamel/answers/T12-solids-flow-scope.md`、`docs/adr/0007-solid-balance-variable.md`、`docs/adr/0010-kinetics-source-gaps.md` 第 9 条、上述 skill。
- **来源限制**：只读 Wirsum 的上述范围和 Hamel 原文。不打开 Wirsum 引用的其它文献（例如 GR 模型的原始论文），只列出作者、年份和 Wirsum 的页码。不打开 `docs/hamel_submodels/` 和归档代码。
- **读图**：字形不清时，先按量纲和物理定义判断并写明理由；仍无法判断的，记下 PDF 页、半页和位置，交目标所有者判读。
- 会话结束时工作树干净，开 PR，审查方填 codex。

## 范围内

- 输出 `refs/hamel/answers/T13-wirsum-solids-scope.md`，按上面六个问题分节，每节按 hamel-source-check 的六问格式（第 6 问写“不适用”）；文末汇总表：问题、Wirsum 结论、Wirsum 页码、与 Hamel 是否相容、能否用于 Hamel 主模型的依据。
- `.gitignore` 增加 `refs/wirsum/pages`。
- 在 `docs/STATUS.md` 更新 T13 一行。

## 范围外（明确禁止）

- 不写 `hamel/` 下的任何代码。
- 不作采用结论，不裁定；对每处与 Hamel 的不一致只列可选解释和影响。
- 不读 Wirsum 的范围外章节，不引入 Wirsum 的动力学、物性或污泥专用关联式。
- 不提交 Wirsum 的 PDF、页面图像或大段原文（只引用必要的式子和短句）。
- 不改 T8 冻结稿和 T12 答卷；发现 T12 有误时，在本答卷中另写一节“与 T12 的出入”。

## 交付物

- `refs/hamel/answers/T13-wirsum-solids-scope.md`
- `.gitignore` 的一行
- `docs/STATUS.md` 中 T13 一行的更新

## 验收证据

- 六个问题都有结论（Wirsum 有依据 / Wirsum 未给出 / 有歧义），附 Wirsum 页码、德文原句，以及与 Hamel 的逐条对照。
- 每个式子都有“能否用于 Hamel 主模型”的依据行。
- codex 对照 Wirsum 和 Hamel 页面图像逐条审查。

## 停止条件

到达时间上限，或需要读范围外章节才能回答时，停下，提交已完成的部分，写明卡在哪里。“Wirsum 也没有给出”本身是有效结论。
