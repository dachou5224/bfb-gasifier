# T11：cell 衡算残差（式 2.4–2.8）

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话，见“执行要求”） |
| 审查方 | codex（代码类，CHARTER §8） |
| 前置任务 | T8（转录冻结稿 `aa856f9`）；ADR 0007（已批准）；ADR 0008（须批准后才开工） |
| 分支 | `task/T11-cell-residuals` |
| 时间上限 | 2 个工作日 |
| 状态 | 草稿 |

## 目标

在 `hamel/balances/` 中写出式 (2.4)–(2.8) 的残差函数：给定一列 cell 的各股物流和能量流，返回每个 cell、每个衡算变量的残差。这是第 1 关 Newton 求解的方程本体。完成后，残差的结构（符号、方向、相邻 cell 的耦合、单位）有测试锁定，后续模块只需提供各股流的数值。

## 执行要求（给 cursor）

- **新会话**：开一个新的 cursor 会话执行本卡，不要沿用做过 T4、T8、T9 或 worktree 清理的会话。本卡不是独立转录任务，但需要干净的上下文，避免带入旧代码的写法。
- **新 worktree**：
  ```
  git worktree add ../bfb-T11 -b task/T11-cell-residuals origin/main
  ```
  按 `.claude/skills/hamel-source-check/SKILL.md` 第一节第 2 项链接页面图像。
- **开工前必读**：`docs/CHARTER.md`（尤其 §5、§6、§7）、本卡、`refs/hamel/transcription/ch2_balances.md`、`refs/hamel/transcription/ch2_ch3_1_diff.md`、`docs/adr/0007-solid-balance-variable.md`、`docs/adr/0008-eq26-size-class-subscript.md`、`hamel/cases/balance.py`（看现有代码风格）。
- **不得打开**归档代码中的 cell 衡算实现（tag `archive-2026-10-12` 的 `src/core/cell_balances.py`、`src/core/solid_balance_reduced.py` 等）。T1 已证明旧衡算会凭空产生碳；本卡只依据冻结稿和 ADR 从头写。
- 会话结束时工作树干净，开 PR，审查方填 codex。

## 编号与方向（按原文，不得改）

cell 自上而下编号：i−1 在上，i+1 在下（ADR 0007“编号方向”；Bild 2.4，S. 18）。

- 气体只向上流：Ṅ_d,j,i 和 Ṅ_b,j,i 是从 cell i 流入 i−1 的流；cell i 的入流是 Ṅ_d,j,i+1 和 Ṅ_b,j,i+1。
- 固体：ṁ_auf,j,k,i 从 i 流入 i−1；ṁ_ab,j,k,i 从 i 流入 i+1。
- 能量：Ḣ_auf,i 从 i 流入 i−1；Ḣ_ab,i 从 i 流入 i+1。

## 范围内

函数放在 `hamel/balances/`，按需分文件（例如 `gas.py`、`solid.py`、`energy.py`）。全部用 numpy 数组，第 0 轴是 cell（下标 0 为最上层 cell），SI 单位。

1. **式 (2.4)、(2.5) 气体物流残差**（S. 16）。输入每个 cell、每种气体的 Ṅ_zu、Ṅ_rez、Ṅ_r（各相分开）、出流 Ṅ_d,j,i、Ṅ_b,j,i，以及交换流 Ṅ_ex,bd,j,i，单位 mol/s。
   - 残差严格照印刷式：悬浮相 0 = Ṅ_zu,d + Ṅ_rez,d + Ṅ_d,i+1 + Ṅ_r,d − Ṅ_d,i − Ṅ_ex；气泡相 0 = Ṅ_zu,b + Ṅ_rez,b + Ṅ_b,i+1 + Ṅ_r,b − Ṅ_b,i + Ṅ_ex。
   - 最下层 cell 没有来自 i+1 的入流；最上层 cell 的出流离开床层。函数另外返回最上层 cell 的两相出流，供上游耦合。
   - Ṅ_ex 只按印刷式的正负号进入残差。它的物理正方向（由哪一相进入哪一相）尚未裁定（T8 歧义第 3 条），本卡不写浓度差、不写交换速率，docstring 注明这一点。
2. **式 (2.6) 固体残差，按 ADR 0007 的存量形式**（S. 17）。输入：
   - 存量 m_j,k,i（kg），形状 (N_z, N_class)，多种固体时按固体种类分别调用或在第 1 轴展开，在 PR 中说明取舍；
   - K_auf,i、K_ab,i（s⁻¹），形状 (N_z,)，**作为输入**，函数内不计算 K（ADR 0007 尚未裁定第 2、3 项）；
   - ṁ_zu、ṁ_rez、ṁ_aus、ṁ_r、ṁ_←、ṁ_→（kg/s），均作为输入。

   残差按 ADR 0007 裁定第 3 条的供体下标写出，来自 i+1 的上行入流按粒级 k（ADR 0008）。每个 cell 都减去自身的 K_auf·m 与 K_ab·m；最上层 cell 的上行出流和最下层 cell 的下行出流不进入任何 cell，函数把这两股流作为边界流返回，不自行判定它们去向（ADR 0007 尚未裁定第 4 项）。
3. **式 (2.7) 能量残差**（S. 19）。输入每个 cell 的 Ḣ_zu、Ḣ_rez、Ḣ_auf、Ḣ_ab、Ḣ_aus、Q̇_W、Q̇_WÜ（W），残差 = Ḣ_zu,i + Ḣ_rez,i + Ḣ_auf,i+1 + Ḣ_ab,i−1 − Ḣ_auf,i − Ḣ_ab,i − Ḣ_aus,i − Q̇_W,i − Q̇_WÜ,i。两端 cell 的处理同第 2 项：缺失的邻居入流为零，越出床层的出流作为边界流返回。
4. **式 (2.8) 物流的能量流**（S. 19）。给定一股物流的各相气体摩尔流与摩尔焓、各粒级固体质量流与比焓，返回能量流（W）：Σ_l Σ_j Ṅ_l,j·h_l,j + Σ_j Σ_k ṁ_j,k·h_j,k + ṁ_s·h_s。
   - ṁ_s·h_s 是**必填参数**（W），不给默认值。下标 s 的含义尚未裁定（T8 歧义第 4 条），由调用方显式给出；docstring 注明。
   - 焓本身（含标准生成焓）由调用方给出，本卡不实现物性。
5. **测试** `tests/test_balances_residuals.py`，期望值按 CHARTER §7 的优先级，只用手算和性质测试：
   - 手算样例：2–3 个 cell，每个式子至少一组手算数值（在测试注释里写出算式）；
   - 守恒：把所有 cell 的残差相加，内部的相邻流逐对抵消，只剩外部进出、反应源项、交换项（气体两相相加时也抵消）和边界流；
   - 方向：只在最下层 cell 加入一股外部流，检查它只出现在该 cell 的残差中，并经由出流进入上一层；
   - 单位：式 (2.6) 残差全部是 kg/s 量级的项，构造一个 m 很大、K 为 0 的样例，确认残差中没有 m 直接出现（禁止 kg/s 与 kg 相减）；
   - 输入检查：形状不一致、出现 NaN 时报错。负值不截断、不报错，按原值进入残差（CHARTER §5）。
6. **docstring 引用**：每个公开函数写 `Source: Hamel Eq.2.x, S. N`，页码见 `refs/hamel/index.json`（式 2.4、2.5 为 S. 16；式 2.6 为 S. 17；式 2.7、2.8 为 S. 19）。用到 ADR 的地方注明 `解释性选择：ADR 0007` / `ADR 0008`。运行 `python3 tools/check_citations.py` 必须通过。

## 范围外（明确禁止）

- 不计算任何一股流的数值：K_auf、K_ab、交换流、反应源项、粒级转移、焓都是输入。流动模型、动力学、物性、热解不在本卡。
- 不写 Newton 求解器、不组装雅可比、不定义未知向量的排列（留给 `hamel/solver/` 的任务卡）。
- 不裁定 T8 的其它歧义：交换项正方向、式 (2.8) 的 s、式 (2.9) 的逆号、两端边界流的去向。遇到需要裁定时停下，写进 PR。
- 不对残差加任何修改：不 clip、不设下限、不加锚点项、不按条件关闭某一项（CHARTER §5）。
- 不新增配置字段（`PhysicalConfig`、`SolverConfig` 不改）。
- 不从归档代码迁移任何函数。

## 原文答卷

本卡的原文依据是 T8 冻结稿 `refs/hamel/transcription/ch2_balances.md` 中式 2.4–2.8 的六问，以及 ADR 0007、0008，已经审查，不另交答卷。若发现冻结稿与原页不符，停下，在 PR 中写明页码和位置，不自行改动冻结稿。

## 交付物

- `hamel/balances/` 下的残差函数（预计合计不超过 400 行，不含测试）
- `tests/test_balances_residuals.py`
- `docs/STATUS.md` 中 T11 一行的状态更新

## 验收证据

- `pytest` 全部通过，默认测试集仍在 5 分钟内。
- `python3 tools/check_citations.py` 与 `python3 tools/drift_check.py` 通过。
- codex 逐式对照 S. 16、S. 17、S. 19 的印刷式和 ADR 0007、0008，核对每一项的正负号和下标。
- 守恒测试、方向测试、单位测试各至少一个，测试注释写出手算过程。

## 停止条件

到达时间上限，或需要越出范围（例如某一项必须先裁定才能写）时，停下并在 PR / `docs/STATUS.md` 中写明：做到哪里、卡在哪里、需要什么决策。
