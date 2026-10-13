---
name: bfb-scripts-hygiene
description: BFB scripts 目录瘦身与维护。盘点 audit/benchmark 冗余、登记 SCRIPTS_INDEX、迁移 _deprecated、合并 _lu_harness 重复逻辑。瘦身前后须跑 validate_slimming_plan 与 slimming_audit 测试。
---

# bfb-scripts-hygiene（Scripts 瘦身维护）

管理 `scripts/` 膨胀：区分 **gate 主链** vs **一次性实验**，登记索引、安全归档、合并重复。  
**不是**通用架构重构（勿用 mattpocock `improve-codebase-architecture`）。

## 与兄弟技能的分工

| 技能 | 何时用 |
|------|--------|
| **本技能** | 脚本太多、重复 audit、未登记脚本、想迁 `_deprecated/` |
| `bfb-diagnosing-bugs` | 用脚本**查 bug**；本技能管脚本**存废** |
| `gitnexus-impact-analysis` | 删/移前查引用链 |
| `gitnexus-refactoring` | 安全 extract 公共模块 |

## 何时启用

- 根目录 `scripts/*.py` 增多或 **未出现在 `SCRIPTS_INDEX.md`**
- 多个 `audit_*` / `benchmark_phase2_*` / `calibrate_*` 逻辑重复
- 准备删除 `run_phase_gate_XX.py` 等薄包装
- 用户说「脚本冗余 / 瘦身 / 清理 scripts / deprecated」

## 政策（权威来源）

`scripts/SCRIPTS_INDEX.md` §Policy：

1. 根目录只留 **主链入口 + 高频 audit**
2. 历史/一次性 → `scripts/_deprecated/`（见 `_deprecated/README.md`）
3. 新增脚本前 **先查是否可复用** `_lu_harness`、`run_module_audits.py`

分级台账见 [`SCRIPT-TIERS.md`](SCRIPT-TIERS.md)（2026-07-03 盘点：34 未登记，11 零引用候选归档）。

---

## Phase 1 — 盘点（只读）

### 1.1 数量与索引差

```bash
cd bfb-gasifier
# 根目录脚本数
find scripts -maxdepth 1 -name '*.py' | wc -l
# 对比 SCRIPTS_INDEX Active 列表
grep '^- `' scripts/SCRIPTS_INDEX.md | wc -l
```

### 1.2 引用扫描（删前必做）

对每个候选脚本 `foo.py`：

```bash
rg -l 'scripts/foo\.py|foo\.py' --glob '*.py' --glob '*.md' .
```

或用 **GitNexus**：`impact({target: "scripts/foo.py", direction: "upstream"})`

### 1.3 Gate 依赖（不可删）

以下编排器定义的脚本为 **Tier S**：

| 编排器 | 路径 |
|--------|------|
| 分阶段测试 | `scripts/run_test_sequence.py` → `STAGES` |
| 模块审计 | `scripts/run_module_audits.py` → `MODULE_AUDITS` |
| Phase gate | `scripts/run_phase_gate.py` |
| 瘦身回归 | `tests/test_codebase_slimming_plan.py` |
| 诊断地图 | `.cursor/skills/bfb-diagnosing-bugs/SYMPTOM-MAP.md` |

`tests/slimming_plan_utils.py` 可编程收集 gate pytest 节点与 module audit 路径。

### 1.4 重复度启发式

同一脚本内若同时出现：

- `from _lu_harness import build_phase2_reactor, run_phase2_init, write_json`
- 自写 `argparse` + `json.dump` + init/solve 循环

→ 标记为 **合并候选**（见 SCRIPT-TIERS §Tier M）。

---

## Phase 2 — 分级决策

按 [`SCRIPT-TIERS.md`](SCRIPT-TIERS.md) 打标签：

| Tier | 动作 |
|------|------|
| **S** Gate 主链 | 登记 `SCRIPTS_INDEX.md`；禁止删 |
| **A** 文档/验证工具 | 登记；可重构不可归档 |
| **B** 零引用实验 | → `_deprecated/<category>/` 或 `scripts/experiments/` |
| **C** 薄包装 | 删；改文档为统一入口 |
| **M** 合并簇 | extract → `_audit_common.py` 或扩展 `_lu_harness.py` |

**一次 PR 只处理一个 Tier**（例如仅 Tier C 删包装，或仅 Tier B 迁 3 个 sweep）。

---

## Phase 3 — 执行（小步可回滚）

### 3.1 归档（Tier B）

1. 移到 `scripts/_deprecated/{nr|cell|calibration|...}/`
2. 在 `_deprecated/README.md` 表格追加一行
3. 从根目录删除原文件
4. **不**改 `src/` 除非 gate 仍引用旧路径

### 3.2 删包装（Tier C）

例：`run_phase_gate_20.py` → 文档改为 `python3 scripts/run_phase_gate.py 20`  
更新 `SCRIPTS_INDEX.md`；跑 `tests/test_run_phase_gate_cli.py`。

### 3.3 合并（Tier M）

1. Extract 纯函数到 `scripts/_audit_common.py`（或 `_lu_harness.py` 若全 Phase2 通用）
2. 原脚本改为薄 CLI 或 deprecated
3. 禁止新增 logging；保持现有 argparse 风格

### 3.4 登记 INDEX（Tier S/A）

在 `SCRIPTS_INDEX.md` §Active Scripts 追加 `- \`name.py\`` 与一行用途。

同步更新 `tests/test_codebase_slimming_plan.py`：

- `INDEX_ACTIVE_SCRIPTS` — 新登记的主链脚本
- `MOVED_TO_DEPRECATED` — 新迁移映射

---

## Phase 4 — 验收（完成标准）

全部通过才可宣称瘦身完成：

```bash
cd bfb-gasifier && pytest tests/test_codebase_slimming_plan.py -m slimming_audit -q --tb=short
cd bfb-gasifier && python3 scripts/validate_slimming_plan.py
cd bfb-gasifier && python3 tests/sanity_checks.py
```

- [ ] 无 active 代码 `import` `scripts/_deprecated`（slimming 测试已覆盖）
- [ ] `run_test_sequence` / `run_module_audits` 引用的脚本均存在
- [ ] `SCRIPTS_INDEX.md` 与根目录 **Tier S+A** 脚本一致
- [ ] 迁移脚本在 `_deprecated/README.md` 有记录

`validate_slimming_plan.py` 返回非 0 → **禁止**继续删脚本；先修 Tier A/B 失败项。

---

## Phase 5 — 输出模板

```markdown
## Scripts 瘦身 — {范围} — {日期}

**Tier 处理**: B×3 归档, C×5 删包装, S×8 登记 INDEX

| 脚本 | 动作 | 引用复核 |
|------|------|----------|
| `sweep_phase2_*.py` | → `_deprecated/calibration/` | rg: 0 hits |
| … | | |

**validate_slimming_plan**: PASS/FAIL  
**sanity_checks**: PASS/FAIL  
**待办**: 合并 audit_init_* 簇（Tier M，另开 PR）
```

---

## 禁止事项

- ❌ 未跑 `slimming_audit` 就批量删根目录脚本
- ❌ 把 gate 脚本直接删而不改 `run_test_sequence` / `run_module_audits`
- ❌ 把 `_deprecated/` 脚本移回主链而不更新 slimming 测试
- ❌ 为瘦身引入 logging、抽象基类或新设计模式
- ❌ 用模拟结果「证明」可删脚本（删依据是 **引用 + gate**，不是数值）

---

## 相关文件

- [`SCRIPT-TIERS.md`](SCRIPT-TIERS.md) — 分级台账与候选列表
- `scripts/SCRIPTS_INDEX.md` — Active / Deprecated 政策
- `scripts/_lu_harness.py` — Phase2 harness SSOT
- `scripts/validate_slimming_plan.py` — 瘦身验收 runner
- `tests/test_codebase_slimming_plan.py` — 瘦身假设编码
- `tests/slimming_plan_utils.py` — gate 引用收集
- `scripts/_deprecated/README.md` — 迁移记录
