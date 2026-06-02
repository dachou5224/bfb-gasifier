# Hamel 同构级未闭合项落地清单

目标：将“数学同构级复现”剩余差距拆为可执行任务。  
依据：`docs/hamel_original_algorithm_reconstruction.md`（Final v1）与 `docs/hamel_evidence_matrix.md`。

## 0. 验收口径（固定）

- 目标级别：数学同构级（连接矩阵 + 侧块 Jacobian + Fortran 结构尽量一致）。
- R10 口径：按原文单位记号执行（当前核验口径为 Pa 指数项），统一实现与文档。

---

## 1. P0：连接矩阵同构（Verbindungsmatrix）

### 1.1 任务
- 在求解主链中引入可配置拓扑矩阵，而非仅线性 bed/freeboard 序列。
- 显式支持旋风/连接管等拓扑 cell 的入出边定义。

### 1.2 目标文件
- `src/core/reactor.py`
- `src/core/connectivity.py`
- `src/core/connectivity_graph.py`
- `src/core/freeboard_bridge.py`（如涉及桥接同步）

### 1.3 验收标准
- 在同一算例下，拓扑图可表达 bed -> freeboard -> cyclone -> return_leg -> bed 回路。
- 连接关系进入残差组装，而不是仅后处理回写。

---

## 2. P0：Jacobian 侧块同构（主带 + Nebenelemente）

### 2.1 任务
- 将返料/侧链耦合显式注入 Jacobian 结构，不再依赖“外循环显式刷新”替代。
- 固化 unknown ordering 与 block layout，形成可复核说明。

### 2.2 目标文件
- `src/solvers/global_nr_solver.py`
- `src/solvers/structured_jacobian.py`

### 2.3 验收标准
- Jacobian 结构化输出可区分：主三对角块、侧块、边界块。
- 回路耦合变量对残差的偏导不再为“隐式/外迭代近似”。

---

## 3. P0：Vorabrechnung 同构化（single-shot 链）

### 3.1 任务
- 锁定 thesis strict 预算路径为主路径（single-shot 源项冻结）。
- 明确与工程 outer-refresh 分支的隔离关系，避免双口径混淆。

### 3.2 目标文件
- `src/solvers/vorabrechnung.py`
- `src/core/reactor.py`
- `docs/hamel_evidence_matrix.md`（状态收口）

### 3.3 验收标准
- 预算模块在同构模式下只执行一次核心源项生成。
- 内层 NR 不重复重算 DAEM/预算源项。

---

## 4. P1：R10 单位口径统一（论文记号口径）

### 4.1 任务
- 清理 `tar_reactions.py` 中“待核验”遗留注释，替换为固定口径说明。
- 在规格文档写死论文单位口径，禁止多口径并存。

### 4.2 目标文件
- `src/kinetics/tar_reactions.py`
- `specs/hamel_checked_algorithm_spec.md`
- `docs/hamel_evidence_matrix.md`

### 4.3 验收标准
- 所有 R10 计算路径只接受 Pa 输入并按论文口径计算。
- 文档与代码不再出现“bar/Pa 待定”表述。

---

## 5. P1：A1 同构性审计闭合

### 5.1 任务
- 对 `gibbs_minimizer` 与 `gibbs_hamel_reduced` 给出“同构项/近似项”逐项对照。
- 输出一份最终 A1 同构审计结论，替代当前争议状态。

### 5.2 目标文件
- `src/thermodynamics/gibbs_minimizer.py`
- `src/thermodynamics/gibbs_hamel_reduced.py`
- `docs/hamel_evidence_matrix.md`

### 5.3 验收标准
- 每个 A1 方程项都能在实现中定位到对应代码段或明确标记“近似替代”。
- 审计结论从“争议”变为“同构/非同构（有证据）”。

---

## 6. P2：自由板 exact_hamel 主线化判定

### 6.1 任务
- 在同构模式下增加 `exact_hamel` 的固定回归集，不以短预算结果判断。
- 若未达标，保持实验分支地位并明确阻塞项。

### 6.2 目标文件
- `src/core/freeboard_segment.py`
- `docs/hamel_evidence_matrix.md`

### 6.3 验收标准
- 相同算例下，`exact_hamel` 在收敛后可稳定复现实验/文献趋势。
- 默认路径切换必须有对照数据支撑。

---

## 7. 文档收口顺序（建议）

1. 更新 `docs/hamel_evidence_matrix.md`（同步一致性结论）。  
2. 更新 `docs/hamel_original_algorithm_reconstruction.md`（Final v2）。  
3. 更新 `specs/hamel_checked_algorithm_spec.md`（同步实现细节口径）。  
4. 对外只引用上述三份文档。

---

## 8. 最小测试门禁（每完成一项必须执行）

- `python3 tests/sanity_checks.py`
- 针对改动模块的最小回归脚本（如 R10、freeboard、global NR 专项脚本）
- 同构模式下至少 1 个固定案例做“结构化 Jacobian + 拓扑回路”验证
