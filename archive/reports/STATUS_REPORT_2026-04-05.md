# BFB 气化炉状态报告

## 当前状态

### 1. 已确认稳定的修正
- **species index / thermo mapping 审计已通过**，确认 `GAS_SPECIES`、`GAS_SPECIES_INDEX`、NASA 7 系数映射没有错位。
- **phase exchange conservation 审计已通过**，bubble / dense 相间交换守恒仍成立。
- **solver / conservation sequence 风险已修正**：
  - `homotopy rate_multiplier` 不再被 `residuals()` 覆盖。
  - cell 接受态、GS 外层 `state_relax` 之后会重新评估 residual / metric。
  - drying/pyrolysis 的 `R_solid` 不再累积污染。
- **reactor-level drying / pyrolysis 传播物理已修正**：
  - recycle solid stream 只回流 `char/ash`
  - axial propagated solid stream 只传播 `char/ash`
  - LU 工况下 drying / pyrolysis 已收缩到底部低层，`fresh-feed` 的 moisture / VM 基本在 `cell 1` 前后完成释放。

### 2. 本轮新完成的改进
- **O2 限速逻辑已从“整池统一缩放”改为“bubble / dense 分相限速”**。
- **`N_ex(O2)` 已计入分相可用氧供给**：
  - dense 相得到正向 `N_ex(O2)` 时，可用于 dense 相耗氧反应；
  - bubble 相得到负向 `N_ex(O2)` 时，可用于 bubble 相耗氧反应。
- `scripts/audit_oxygen_reaction_trace.py` 已同步到分相 O2 供给 / limiter 口径。
- 新增并固定回归测试：
  - `tests/test_cell_kinetics.py`
  - 覆盖分相 O2 limiter
  - 覆盖 `N_ex(O2)` 对 dense 相可用氧的贡献

### 3. 当前验证状态（以当前稳定实现为准）
- **sanity gate**：`ALL 11 SANITY CHECKS PASSED`
- **聚焦回归**：`16 passed`
- **LU 基线结果**（当前稳定状态）：
  - `T_exit ≈ 1139.7 K`，出口温度表现良好
  - `carbon_conv ≈ 90.6%`
  - `exit_gas_dry` 约为：
    - `CO ≈ 0.219`
    - `CO2 ≈ 0.0737`
    - `H2 ≈ 0.145`
    - `CH4 ≈ 0.0118`
    - `O2 ≈ 0.0017`
    - `N2 ≈ 0.547`

### 4. 当前核心判断
- 当前主问题**不再是**：
  - species index 错位
  - thermo 数据错绑
  - drying / pyrolysis 全床传播
  - 相间交换守恒错误
  - solver 顺序覆盖 earlier conservation
- 当前主问题**集中在后段气组成仍未对上**：
  - `CO` 仍偏高
  - `CO2` 仍偏低
  - `H2` 仍偏高一点
  - `CH4` 仍偏低
- `N2` 不是主矛盾。按干基平衡估算，目标 `N2` 大约在 `~0.53`，当前 `0.547` 只高约 `0.016`。

### 5. 本轮诊断结论
- `audit_oxygen_reaction_trace.py` 明确显示：
  - **bubble 相耗氧需求，特别是 `R5b`，仍异常主导 O2 预算**
  - 在很多 cell 中，`R5b` 远大于 dense 相对应耗氧路径
  - 分相 O2 limiter 已经把问题显性化，但尚未根治 gas composition 偏差
- 我尝试过一条 **R5b 公式修改分支**（给 bubble R5 补 `H2O^0.5` / 压力基准口径），但会把 LU validation 拉坏，因此**已完整回退**，未保留。

---

## 下一步计划

### 主线目标
继续缩小 LU validation 的剩余偏差，优先解决：
1. `CO` 偏高
2. `CO2` 偏低
3. `H2` / `CH4` 分布仍不对

### 下一步具体动作
1. **继续拆解 bubble 相耗氧路径**
   - 重点检查 `R5b` 与 `R6b` 的量级来源
   - 区分是动力学表达式本身过猛，还是 bubble 相组分 / 停留时间 / 体积使用口径偏激进

2. **检查 bubble 相氧化是否缺少物理约束**
   - 是否应受额外 mixing / residence-time / transfer window 约束
   - 是否存在“按整相体积一次性完全快反”导致的放大

3. **继续对照 dense 相后段二次转化**
   - 保持现有分相 O2 limiter + `N_ex(O2)` 供给逻辑不动
   - 继续检查 `R8`、dense 相 `R5d` 以及 H2O 供给是否足以把 CO 往 CO2 拉

4. **谨慎避免无效分支**
   - 不再重复使用已验证会恶化 LU 的 `R5b` 公式改法
   - 后续改动优先采用“小步、可回退、先审计后实现”的方式

### 当前推荐工作顺序
1. 先保留当前稳定实现作为基线
2. 进一步做 `R5b/R6b` 来源分解审计
3. 仅在有明确证据时再改 bubble 氧化链路
4. 每次改动后继续以：
   - `python3 tests/sanity_checks.py`
   - LU validation 审计
   - oxygen trace 审计
   作为门控
