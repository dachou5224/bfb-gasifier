# BFB Gasifier 项目核心准则

作为本项目的人工智能工程师，你必须始终遵循以下准则，这些准则源自 `docs/CLAUDE.md` 和科学建模最佳实践。

## 1. 强制性开发顺序
必须严格按 Phase 1→6 的顺序推进任务。在开始新 Phase 前，确保前置 Phase 的所有函数已实现并通过单元测试。

## 2. 科学建模审计 (Critical)
- **单位制**：全程 SI (Pa, K, m, mol, kg, s)。R8 方程内部换算 atm 是唯一允许的例外。
- **动力学**：禁止在 `arrhenius.py` 之外直接编写指数速率方程。必须使用 `k_hobbs`, `k_standard`, `k_jensen_r7`。
- **常数**：禁止硬编码。必须从 `src/core/constants.py` 导入 `Rg`, `g`, `P0` 等。
- **稳定性**：所有的 `np.exp` 必须有溢出保护（使用工厂函数已自带此功能），所有的关键计算函数必须有 `assert` 校验维度。

## 3. 验证门控
每次修改代码后，**必须**运行：
`python3 tests/sanity_checks.py`
只有当输出为 `ALL 11 SANITY CHECKS PASSED` 时，任务才算完成。

## 4. 引用要求
每个物理/化学方程的实现函数，必须在 docstring 中标注 Hamel (1999) 的方程编号。
