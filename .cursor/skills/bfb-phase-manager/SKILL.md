---
name: bfb-phase-manager
description: BFB 项目开发顺序门禁。强制遵循 docs/CLAUDE.md 中 Phase 1→6，启动新任务前检查前置 Phase 是否完成并运行 tests/sanity_checks.py。
---

# bfb-phase-manager（开发进度管理）

## 何时启用

- **每次新任务 / 新会话开始**、或准备实现**跨 Phase** 功能（例如直接改 `reactor.py` 却未确认 Phase 3 动力学）时。
- 用户要求「按 CLAUDE 开发顺序」或涉及 **Phase 4（干燥/热解）**、**Phase 5（求解器）**、**Phase 6（验证/UI）** 时，必须先执行本技能。

## Phase 顺序（权威来源）

严格遵循 `docs/CLAUDE.md` 中 **「开发顺序（严格按此顺序，不跳步）」**：

| Phase | 内容摘要 |
|-------|----------|
| 1 | `species.py`、`arrhenius.py`、`validation_cases.json` |
| 2 | 流体力学（u_mf、气泡、相含率、自由板） |
| 3 | 动力学（炭/气/tar 反应） |
| 4 | `drying.py`、`devolatilization.py` |
| 5 | `cell.py`、`cell_solver.py`、`reactor.py`、`sanity_checks.py` |
| 6 | `test_table2_LU.py`、`app.py` |

## 任务启动检查清单

1. **明确当前任务对应的 Phase 与文件**（只改允许范围内的文件；若必须跨 Phase，先说明依赖风险）。
2. **前置 Phase 是否已具备可运行基础**（例如改 Phase 5 前，Phase 3–4 相关模块是否已存在且可被导入）。
3. **运行 sanity**（对 `src/`、`tests/`、影响数量级的 `specs/`/`docs/` 有实质修改时，收尾前必须执行）：
   ```bash
   cd bfb-gasifier && python3 tests/sanity_checks.py
   ```
   - 期望：`ALL 11 SANITY CHECKS PASSED`（以脚本当前输出为准）。
   - **未通过则不得宣称该阶段「已完成」**；应修复或记录阻塞原因。

## 与复杂实现前的 gstack 技能配合

在进入 **Phase 4 或 Phase 5** 的大型改动前，建议先 **@** 使用 `.cursor/skills` 中的 **`gstack-plan-eng-review`**（见 `docs/cursor_skills_playbook.md`），锁定架构与数据流，避免与单位制、Arrhenius 约定冲突。

## 相关文件

- `docs/CLAUDE.md`（Phase 列表与修订后工作流）
- `tests/sanity_checks.py`
- `docs/validation_gap_analysis.md`（算法层差异）
