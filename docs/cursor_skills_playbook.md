# Cursor 技能手册（BFB Gasifier）

本文档说明 **本仓库自定义 Agent Skills** 的用途，以及何时 **@** 引用工作区根目录下的 **gstack** 技能副本。

---

## 一、本仓库自定义技能（`bfb-gasifier/.cursor/skills/`）

在 Cursor 对话中用 **`@`** 引用对应目录下的 `SKILL.md` 即可激活（例如 `@bfb-gasifier/.cursor/skills/bfb-scientific-audit/SKILL.md`）。

| 技能目录 | 职责 | 建议触发时机 |
|----------|------|----------------|
| **`bfb-scientific-audit`** | 科学建模审计：SI 单位、`arrhenius.py` 工厂函数、单位注释、Hamel (1999) 方程溯源 | Validate 阶段；使用 `/review` 类审查**前** |
| **`bfb-phase-manager`** | 强制 `docs/CLAUDE.md` Phase 1→6 顺序；任务启动时检查前置 Phase 与 `sanity_checks.py` | **新任务/新会话开始**；跨 Phase 改动前 |
| **`bfb-numerical-guard`** | 数值稳定性：`fsolve` / `global_nr` 残差缩放、`np.exp` 保护、稀疏 Jacobian 结构 | 修改求解器、残差或收敛逻辑后；NR 不收敛时 |

与已有技能并列参考：

- `hydrodynamics-skill`、`kinetics-balance-skill`（同目录）

---

## 二、建议激活的 gstack 技能（副本路径）

gstack 技能已复制到 **仓库根目录** `AI-projects/.cursor/skills/`，命名形如 `gstack-<命令名>`。在 BFB 场景下的**建议用法**如下。

### 1. `gstack-plan-eng-review`（对应 `/plan-eng-review`）

- **何时用**：准备实现 **Phase 4（干燥/热解）** 或 **Phase 5（全局求解器 / 多 cell 耦合）** 等复杂逻辑**之前**。
- **目的**：锁定架构、数据流、状态边界与测试矩阵，减少单位混用、硬编码常数、隐式耦合。
- **引用示例**：`@.cursor/skills/gstack-plan-eng-review/SKILL.md`（从 `AI-projects` 根打开时）

### 2. `gstack-investigate`（对应 `/investigate`）

- **何时用**：`tests/sanity_checks.py` **失败**；或 **Newton–Raphson / fsolve 不收敛**、残差异常但原因不明时。
- **目的**：按「先根因、后补丁」系统化排查，避免盲目改阻尼或容差。
- **引用示例**：`@.cursor/skills/gstack-investigate/SKILL.md`

### 3. `gstack-qa`（对应 `/qa`）

- **何时用**：验证 **`app.py`** Streamlit 界面；检查 **UI 展示的物理量**（温度、组成、流量等）与 **后端计算** 是否一致。
- **说明**：gstack 原文依赖浏览器自动化工具链；在 Cursor 中主要采用其 **测试思路与清单**，工具可用本地浏览器 + 手工对照或项目已有 E2E 方式替代。
- **引用示例**：`@.cursor/skills/gstack-qa/SKILL.md`

---

## 三、推荐组合（工作流）

1. **任务开始** → `bfb-phase-manager` → 确认 Phase 与 `sanity_checks.py`。
2. **大改 Phase 4/5 前** → `gstack-plan-eng-review`。
3. **提交 / Validate 前** → `bfb-scientific-audit` →（可选）gstack `review` 流程。
4. **动求解器或残差** → `bfb-numerical-guard`；仍失败 → `gstack-investigate`。
5. **动 UI** → `gstack-qa` 思路验证 `app.py`。

---

## 四、维护说明

- **CLAUDE.md** 中 Phase 与 `sanity_checks` 要求变更时，请同步更新 **`bfb-phase-manager`** 技能与本文档。
- gstack 副本更新：见仓库根 `.cursor/skills/README.md`（若存在）。
