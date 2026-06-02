# Hamel 结构桥接工程备忘录

更新时间：2026-05-29

目的：记录 bed/freeboard/side-elements 同构化过程中，Hamel 原文有结构证据但没有逐项实现细节的工程补缝方案。本文不是新的物理模型来源；它说明当前 Python 实现如何把 Hamel 的连接矩阵、Nebenelemente、freeboard trajectory 思路落成一个变量数与方程数一致、Jacobian 可组装、残差可审计的 NR 系统。

## 1. 论文证据与工程缺口

Hamel 原文给出了三类关键证据：

- bed 单元使用 dense/bubble 两相气体守恒。
- freeboard 使用 ghost-bubble/颗粒轨迹闭合，不应继续作为 bubbling-bed cell 求解。
- 连接矩阵与 Nebenelemente 表明 cyclone/return-leg 这类侧元素要进入整体拓扑和 Jacobian，而不是后处理。

原文没有直接给出这些工程细节：

- bed 顶层的 `N_d + N_b` 如何投影为 freeboard 的单相 gas inlet。
- freeboard/cyclone/return_leg 的 NR unknown/residual 宽度如何与 bed 不同而仍能拼成全局向量。
- side element 对 bed[0] 的返料偏导如何放入 structured Jacobian。
- freeboard solid trajectory closure 何时刷新，如何避免覆盖 NR 已接受的 gas/T 状态。
- result finalizer 是否允许为了“好看”而重跑 closure 并突变已求解状态。

当前实现把这些缺口明确为工程桥接层，而不是把所有 cell 伪装成 bubbling-bed cell。

## 2. Cell-type-specific NR layout

当前 NR layout 的原则是：每类 cell 只求解它真正拥有的动态/代数 unknown。

| cell type | gas unknown | solid unknown | temperature unknown | 说明 |
|---|---:|---:|---:|---|
| bed | dense + bubble gas | char/ash holdup | yes | Hamel bed 两相主模型 |
| freeboard | total single gas | no NR solid holdup | yes | 固体由 trajectory/closure 负责 |
| cyclone | no gas NR unknown | char/ash holdup | no | solid separation side element |
| return_leg | no gas NR unknown | char/ash holdup | no | solid return side element |

工程规则：

- `n_var(cell)` 不再由统一 bubbling-bed 模板隐式推导，而是由 cell type 决定。
- residual projection 必须满足 `len(projected_residual(cell)) == n_var(cell)`。
- freeboard 不再保留 `N_b = 0` bubble-anchor row；它只有 total gas residual 和 energy residual。
- cyclone/return_leg 不求解反应、两相交换或 gas residual；它们只作为 solid transport side blocks。
- `gas_phase_split` 诊断只统计真实 bed cells，避免 freeboard/cyclone 的退化行污染两相指标。

这个设计是对 Hamel 结构的工程同构化：论文支持 cell 类型不同，但没有给出 Python 向量 layout；这里把 layout 明文化。

## 3. Bed -> freeboard 桥接

bed 顶层到 freeboard 第一格采用显式通量投影：

```text
N_g,fb,in = N_d,bed,out + N_b,bed,out
N_b,fb,in = 0
```

工程含义：

- bubble split residual 只存在于 bed cell 内，不跨入 freeboard。
- freeboard 对 bed 内 dense/bubble 的纯重分配不敏感，只对 total gas 变化敏感。
- freeboard 第一格 energy inlet 使用 bed 顶层温度；gas composition 使用合并后的 total gas。
- freeboard solid 不是 NR unknown；bed top entrainment 通过 trajectory/closure 转成 freeboard solid holdup/outflow coefficients。

已验证/观察：

- 单相入口合并能减少退化段 split 残差，但不是单独的收敛突破。
- 后续 residual hotspot 仍可能出现在 bed 内两相气体或 energy 方程；这时不应把 freeboard 恢复成两相来“凑方阵”。

## 4. Freeboard solid trajectory 与 age/size closure

Hamel §5.1.1 指出 particle age/size classes 会影响 freeboard trajectory 和反应历史。当前实现采用 closure-only quadrature：

- 在 bed 顶层生成 age/size launch samples。
- trajectory 使用扩展后的粒径/年龄样本。
- 输出再聚合回原有 size class 数量。
- 不增加 NR unknown 数量。

这是一个有文献方向、但具体离散方式属于工程实现的桥接方法。它解决的是“均一粒径/密度导致 freeboard 逃逸固体过少”的结构问题，但不把 Hamel 的完整 size-class state 直接塞进全局 Jacobian。

当前 Phase2 LU 默认：

```text
freeboard_age_quadrature_bins = 8
freeboard_age_quadrature_max_age = 0.98
```

## 5. Side elements 与 structured Jacobian

side elements 的工程原则：

- `bed -> freeboard -> cyclone -> return_leg -> bed` 必须作为显式拓扑存在。
- cyclone/return_leg 不伪装成 bubbling-bed cells。
- cyclone/return_leg 只承担固体分离、返料和边界传递。
- return-leg 热固体返料通过 solid stream 温度进入 bed[0] energy balance；return_leg 自身不引入独立 temperature unknown。

Jacobian 规则：

- 主 bed/freeboard 链仍按 band/block structure 装配。
- last freeboard -> cyclone、cyclone -> return_leg、return_leg -> bed[0] 作为 side blocks 进入结构化 Jacobian。
- 对 `band_plus_side_elements_structured`，禁止 local-BC/affected-row pruning，使用 full-BC 组装，避免 recycle/side coupling 被剪掉。
- 小网格 dense-vs-structured 已作为门禁；标准 Phase2 仍需要长算例持续验证。

这部分是典型“论文没给代码细节”的工程补缝：Hamel 说明有 Nebenelemente/连接矩阵，但没有给当前 Python residual/Jacobian assembly 的剪枝边界。

## 6. Outer refresh 与 closure 不突变原则

一个容易误踩的坑是：为了让 freeboard closure 与当前 bed-top state 一致，在 outer refresh 或 result finalizer 中重跑 closure 并回写 explicit freeboard/side state。

当前结论：

- solver path 中不默认强刷 freeboard solid closure。
- result finalizer 只报告 closure freshness/gap，不突变 NR 已接受状态。
- `preserve_gas_state=True` 的 solid-only closure refresh 保留为诊断/审计工具，不作为默认求解路径。

原因：

- 强刷 closure 可以让 bed-top/freeboard gap 在报告层面变小，但会改变 side recycle 与 freeboard solid state，反而使 residual/RMS 变差。
- finalizer 后验突变会造成“返回的 reactor state 与最后 residual 不一致”，应禁止。

因此当前结果中保留 gap 诊断字段，例如：

```text
freeboard_pre_result_refresh_gap_current_bed_top_kg_s
```

它是未闭合耦合的可见诊断，不是 finalizer 应偷偷修掉的问题。

## 7. 已排除的错误路线

- 用 `N_b = 0` anchor 让所有 cell 看起来同宽：可以凑方阵，但破坏 freeboard/side-element 物理语义。
- 把 freeboard/cyclone 继续纳入 phase-split residual 统计：会制造假 hotspot。
- 每个 outer 都强制刷新 freeboard closure：gap 可被压小，但 side recycle 与 residual 质量会变差。
- line-search merit 直接改成 max-abs 主导：短 probe 中更保守但收敛更差。
- 把底部主气化剂 dense/bubble split 固定为调参常数：Hamel 只把
  `Startwertwahl` 放在 configuration，但把 gas/solid streams、porosity、
  residence times 和 phase exchange 放在 `Vorabrechnung`。因此 Phase1/Phase2
  global-NR 主线改为从预计算水力学通量解析入口相分配。

## 8. Bottom gas inlet split

Hamel 没有直接给出 distributor plate 后的 `N_zu,d / N_zu,b` 固定比例。当前实现的桥接口径是：

```text
dense_frac = u_d * (1 - eps_b) / u0
bubble_frac = 1 - dense_frac
```

证据边界：

- 论文 Kapitel 2.1 明确 `Vorabrechnung` 负责上/下行气固流、孔隙率、相间交换和停留时间。
- Kapitel 3.1.2 Eq.3.7-3.10 把总 superficial gas flux 分解为 visible bubble gas、bubble throughflow 和 dense/suspension flow。
- 因此 bottom inlet split 属于 pre-calculation 派生边界条件，而不是独立 kinetic tuning knob。
- `gas_inlet_dense_frac` 仍保留为 `fixed` 策略和预计算尚不可用时的 fallback，不能再解释为 Hamel 证据支持的目标参数。

## 9. 当前收敛状态与后续方向

最新 Phase2 LU 默认采用：

```text
nr_jacobian_lag_steps = 4
freeboard_age_quadrature_bins = 8
```

10-outer probe 的代表性结果：

```text
rms_scaled_final ≈ 0.00182
max_abs_scaled_final ≈ 0.01383
```

解释：

- RMS 和 component RMS 已基本收敛。
- 严格失败点主要是少数 max residual rows 尚高于 0.01。
- 下一步应定位 bed 内 gas/energy 局部残差，而不是回退 freeboard/side layout。

## 10. 维护规则

后续修改遵循以下规则：

- 若新增 cell type，必须先定义 NR layout，再定义 residual projection。
- 若新增 side topology，必须进入 connectivity graph 和 Jacobian side block，而不是只做 result summary。
- 若要刷新 closure，必须明确是 solver-state mutation 还是 diagnostic-only。
- 若与 Hamel 原文不完全相同，文档必须标注“论文证据”和“工程补缝”边界。
