# T0：快照归档与主线重置

| 字段 | 内容 |
|---|---|
| 执行者 | 总司令（Claude），每个外部动作前由目标所有者确认 |
| 前置任务 | 章程 v0.1 获批 |
| 分支 | `archive/2026-10-12-snapshot`、`chore/charter-reset` |
| 时间上限 | 半天 |
| 状态 | 进行中 |

## 目标

保留现有全部工作（一行不丢），并得到一个只含章程、文档与数据的干净主线起点。

## 范围内

1. **确认停工**：目标所有者确认 codex、cursor 已停止 named_f 及相关工作，没有 agent 正在写工作树。
2. **备份仅本地文件**：`.gitignore` 中的 `*.pdf`、`docs/hamel_submodels/`、`docs/hamel_isomorphic_gap_implementation_checklist.md`、`specs/hamel_checked_algorithm_spec.md` 不会进入 git。将它们打包备份到仓库外（如 `~/bfb-gasifier-local-backup-2026-10-12.tar.gz`），记录文件清单与校验和。
3. **快照提交**：从当前分支 `feat/nr-hamel-wirsum-core-startwert` 建 `archive/2026-10-12-snapshot`，`git add -A`（遵守 `.gitignore`）后提交，打 tag `archive-2026-10-12`。提交前列出大于 5 MB 的文件，请目标所有者确认。
4. **推送远端**：经目标所有者确认后，推送快照分支和 tag。
5. **主线重置**：从 `main` 建 `chore/charter-reset`：
   - 删除旧代码与旧文档：`src/`、`tests/`、`scripts/`、`archive/`、`specs/`、`app.py`、`plan.md`、`status_report.md`、`GEMINI.md`、`docs/CLAUDE.md`、各 handoff、`docs/audit_reports/`、`.cursor/skills/bfb-*`、`.cursor/rules/*`，以及 `data/` 中除保留项以外的文件；
   - 保留：`docs/blog/`、`docs/siegen_cell_model_lineage.md`、`data/validation_cases.json`、`data/hamel_fig7_parity_digitized.json`、`pyproject.toml`、`requirements.txt`、`.gitignore`；
   - 放入章程草稿的最终版本：`docs/CHARTER.md`、`docs/STATUS.md`、`docs/tasks/`，以及重写后的 `AGENTS.md`（简短入口＋原博客写作规则）和 `README.md`（说明正在重建，旧版见归档 tag）；
   - `.gitignore` 增加 `refs/hamel/pages/`。
6. 提 PR，由目标所有者合并。
7. **本地清理**：经确认后删除 `venv/`（546 MB），只保留 `.venv`。

## 范围外

- 不修改任何旧代码内容。
- 不对快照做任何“顺手整理”。
- 未经确认不推送、不合并、不删除。

## 交付物

- 快照分支与 tag；仅本地文件的备份包及清单。
- `chore/charter-reset` 的 PR。

## 验收证据

- `git ls-tree -r --name-only archive-2026-10-12 | wc -l` 与快照前工作树文件数（不含被忽略文件）一致。
- 备份包清单与校验和。
- PR diff 中只有删除和新增章程/文档文件，没有修改过的旧代码。

## 停止条件

任何一步与预期不符（例如仍有 agent 在写工作树、出现异常大文件），立即停下，请目标所有者决定。
