# Hamel 验证工况化学地图（Validation Cases Chemistry Map）

> **SSOT 数据**：[`data/validation_cases.json`](../data/validation_cases.json)（Hamel 1999 Ch.7–8，2026-03 数字化）  
> **用途**：多工况验证、调参前的物理化学参考；**不得**作为单 case `r*_scale` 凑数依据  
> **关联**：[`optimization_handoff_summary.md`](optimization_handoff_summary.md) §9、[`CLAUDE.md`](CLAUDE.md)「一致性优先」

---

## 1. 文档目的与使用边界

Hamel 论文共提取 **15 个 CASE**（HTW / VTT / WSV400）。本文件归纳：

1. 进料、ER、氧化剂、压力、出口组成与温度的 **跨工况规律**  
2. HTW 炉 **密相 / 气泡相 / 自由板** 的主导反应（定性，来自轴向剖面 `_observations`）  
3. 后续模型验收与调参的 **物理优先级**（边界 → 温度 → char 链 → WGSR，最后才是 diagnostic scale）

Hamel 自述模型对 CO/CO₂/H₂ 的验证带宽约 **±10–15%**（`validation_performance_summary`）。单点零误差不是论文标准。

**禁止用法**：

- ❌ 在 `CASE_HTW_WESSELING_1`（LU）上标定 `{r2,r4,r5,r8}` 后外推到 O₂/Steam 或泥炭/木材工况  
- ❌ 用出口 CO/CO₂ 数值大小直接反推「应把 r5 拧到 0.25」而无轴向/进风/温度证据  

**推荐用法**：

- ✅ 新 case 偏离时，先查 **§3 化学分支** 与 **§4 轴向分工**，判断缺的是哪类物理闭合  
- ✅ 至少 **2 个不同分支** case 联合验收趋势，再讨论弱标定（见 §6）

---

## 2. 全工况一览

| Case ID | 装置 | 规模 | 燃料 | 氧化剂 | ER | P [MPa] | T_exit [K] | CO | CO₂ | H₂ | CH₄ | 二级风 [Nm³/h] | 轴向剖面 |
|---------|------|------|------|--------|-----|---------|------------|-----|-----|-----|-----|----------------|----------|
| `HTW_WESSELING_1` | HTW Wesseling | 中试 | 莱茵褐煤 | **Air+Steam** | 0.337 | 2.5 | 1100 | 0.130 | 0.110 | 0.120 | 0.028 | 30.1 (air) | ✓ |
| `HTW_WESSELING_2` | HTW Wesseling | 中试 | 莱茵褐煤 | **O₂+Steam** | 0.268 | 2.5 | 1220 | 0.404 | 0.283 | 0.307 | 0.058 | 37.5 (O₂) | ✓ |
| `HTW_BERRENRATH` | HTW Berrenrath | 工业 | 莱茵褐煤 | **O₂+Steam** | 0.238 | 1.0 | 1050 | 0.354 | 0.255 | 0.317 | 0.060 | 0 | ✓ |
| `VTT_PRESSURIZED_COAL_1` | VTT 加压 | 小试 | 褐煤 | Air+Steam | 缺 | 0.5 | 980 | 0.183 | 0.120 | 0.160 | 0.051 | 0 | ✓ |
| `VTT_PRESSURIZED_PEAT_6` | VTT 加压 | 小试 | 芬兰泥炭 | Air+Steam | 0.29 | 0.5 | 1050 | 0.135 | 0.172 | 0.110 | 0.043 | 21.9 | ✓ |
| `VTT_PRESSURIZED_PEAT_13` | VTT 加压 | 小试 | 芬兰泥炭 | Air+Steam | 0.41 | 0.5 | — | — | — | — | — | 41.6 | ✓ |
| `VTT_PRESSURIZED_SAWDUST_14` | VTT 加压 | 小试 | 锯末 | Air+Steam | 0.28 | 0.4 | — | — | — | — | — | 11.6 | ✓ |
| `VTT_ATMOSPHERIC_PEAT_D` | VTT 常压 | 小试 | 芬兰泥炭 | **Air only** | 0.26 | 0.1 | 920 | 0.136 | 0.170 | 0.138 | 0.048 | 0 | ✓ |
| `VTT_ATMOSPHERIC_PEAT_E` | VTT 常压 | 小试 | 芬兰泥炭 | Air only | 0.29 | 0.1 | — | — | — | — | — | 0 | ✓ |
| `WSV400_A` | WSV Magdeburg | 中试 | 天然木材 | Air only | 0.453 | 0.1 | 980 | 0.152 | 0.132 | 0.110 | 0.022 | 0 | ✓ |
| `WSV400_B`–`F` | WSV Magdeburg | 中试 | 天然木材 | Air only | 0.35–0.50 | 0.1 | 多数缺 | 多数缺 | 多数缺 | 多数缺 | 多数缺 | 0 | 部分 |

**干基 CO 份额**（CO/(CO+CO₂)，仅有完整出口数据的 case）：

| Case | CO/(CO+CO₂) | 备注 |
|------|-------------|------|
| HTW_WESSELING_1 | **0.54** | Phase2 主验证 case（LU） |
| HTW_WESSELING_2 | **0.59** | O₂/Steam，绝对 CO≈3× LU |
| HTW_BERRENRATH | **0.58** | 工业尺度 O₂/Steam |
| VTT_COAL_1 | **0.60** | 加压褐煤，无二级风 |
| VTT_PEAT_6 | **0.44** | 高水分泥炭，CO₂ 偏高 |
| VTT_PEAT_D | **0.44** | 常压 Air，T_exit 低 |
| WSV400_A | **0.54** | 木材，高 ER |

---

## 3. 化学分支分区（第一性分类）

比 ER 更根本的是 **氧化剂化学分支** 与 **燃料含水/VM**。

### 3.1 Air + Steam（中等气化 — LU 所在类）

- **代表**：Wesseling 1、VTT 褐煤/泥炭（加压）  
- **出口量级**：CO **0.13–0.18**，CO₂ **0.11–0.17**，H₂ **0.11–0.16**  
- **N₂ 稀释大** → 绝对合成气中等；**CO/(CO+CO₂) 稳定在 ~0.54–0.60**  
- **调参含义**：应复现 **床层 char 产 CO + 适度 WGSR**，而非把出口拉到 O₂/Steam 水平

### 3.2 O₂ + Steam（富合成气、低 ER）

- **代表**：Wesseling 2、Berrenrath  
- **出口量级**：CO **0.35–0.40**，CO₂ **0.25–0.28**，H₂ **0.31–0.32**（约为 LU 的 **2.5–3×**）  
- **关键**：CO/(CO+CO₂) **仍 ~0.58–0.59**，与 LU 接近 — 换氧化剂主要动 **绝对产量**，不是把比例翻到另一套  
- **调参含义**：**禁止**沿用 LU 标定的 global `r*_scale`；须检查 O₂/蒸汽比、ER、无 N₂ 稀释下的 R2/R8 环境

### 3.3 Air only（常压、无蒸汽）

- **代表**：VTT 常压泥炭、WSV 木材  
- **出口量级**：CO **0.14–0.15**，CO₂ **0.13–0.17**；T_exit **920–980 K**（低于 HTW 加压）  
- **CO/(CO+CO₂)** 常 **≤0.54**，泥炭更低（~0.44）— 高 O、高 H₂O 抑制 CO、抬 CO₂  
- **调参含义**：自由板冷却、水分蒸发吸热、无 steam 补氢 — **与 HTW 不是同一模板**

### 3.4 燃料效应（同氧化剂类内）

| 燃料 | 典型 VM / 水分 | Hamel 观察 |
|------|----------------|------------|
| 褐煤 | VM ~53% daf，水分中等 | 验证 **最好**（Ch.8）；床内梯度陡、峰温高 |
| 泥炭 | VM ~68%，**高水分** | 床温 **低 ~100 K**；出口 **CO₂ 偏高**；CH₄ 床内升、二级风处 dip |
| 木材 | 大颗粒、高 VM | CH₄ **低**；干燥区最长；Fig. 7.16 为 **干基** 绘图 |

---

## 4. HTW 轴向：密相 / 气泡相 / 自由板分工

来源：`CASE_HTW_WESSELING_1/2`、`CASE_HTW_BERRENRATH` 的 `outputs.axial_profiles._observations` 与 Fig. 7.4–7.6。

```
ξ ≈ 0.00–0.15   床底 / 锥段
  · O₂ 在气泡相快速消耗 → R1（炭燃烧）+ R5 泡相（CO 氧化）
  · 强放热 → T ↑ 至 ~1020–1100 K
  · 热解 + char 表面 R2/R4 启动；H₂O 在前 30% 消耗最剧烈

ξ ≈ 0.15–0.45   床层主体（密相 / 炭表面主导）
  · 吸热气化：R2 (C+H₂O)、R4 (C+CO₂) → T 由 ~1100 降至 ~1060 K
  · CO、H₂ 在浓相积累；O₂ 已耗尽
  · 出口 CO/CO₂ 基调在此奠定

ξ ≈ 0.45–0.55   床顶 / 稀相过渡
  · 组成接近床层出口

ξ ≈ 0.60        二级风（HTW LU / 部分 VTT）
  · 少量 O₂ → 局部 R5 再氧化 CO/H₂
  · T 峰：Air **1240 K**；O₂/Steam **1380 K**
  · CO、H₂ 瞬态 dip；CO₂ spike
  · 注意：LU 二级 air 仅 ~30 Nm³/h（~0.54% 总 O₂），造峰靠 **局部** 燃烧，非全局 ER

ξ ≈ 0.60–1.0    自由板
  · R8 (WGSR) 再平衡 CO/CO₂/H₂
  · CH₄ 基本平坦（二次焦油/裂解弱）
  · T 缓降 → 出口（LU ~1100 K）
```

### 4.1 相区 ↔ 反应 ↔ 出口控制

| 相区 | Hamel 反应（仓库编号） | 主要控制 |
|------|------------------------|----------|
| 床底气泡相 | R1、R5（有 O₂ 时） | 温升、O₂ 耗尽、部分 CO₂ |
| 床层密相 / 炭表面 | R2、R4、热解气 | **CO、H₂ 主产量**；CO/(CO+CO₂) 基调 |
| ξ≈0.6 二级注入 | 局部 R5 | **T 峰**；CO/H₂ 瞬态 dip |
| 自由板悬浮相 | R8、少量 R5/R6 | **CO ↔ CO₂ 再分配**；H₂ 微调 |

### 4.2 与当前 Phase2 LU 偏差的对应

| 现象 | 物理解读（优先排查） |
|------|----------------------|
| CO 偏低、CO₂ 偏高 | 床层 R2 产 CO 不足；密相/悬浮相 R5 过强；R8 把 CO₂ 拉太高 |
| 缺 ξ≈0.6 T 峰 | secondary 注气与能量；非单纯增大 global r5 |
| bed0 T 钉死 / 剖面 RMSE 大 | 温度 seed、围栏、进风分配 — **先于** r scale |

---

## 5. 跨工况规律（调参物理参考）

1. **氧化剂类型是第一维** — Air+Steam vs O₂+Steam 决定合成气 **绝对水平**（差 2–3×），不是微调。  
2. **CO/(CO+CO₂) 在同类分支内对 ER 不敏感** — LU (0.34) 与 Wesseling 2 (0.27) 比值同属 ~0.55–0.59。  
3. **ER 升高** → 氧化更完全 → CO₂ 倾向上升（WSV 木材 ER 0.35–0.50 系列）。  
4. **二级风** — 有则 ξ≈0.6 T 峰 + CO dip；无则自由板单调降温（VTT 褐煤、Berrenrath 主注入在 ξ≈0.5）。  
5. **压力** — 2.5 MPa (LU) vs 0.1 MPa (常压) 改变气速与 R8 窗口；外推须改 **操作边界**。  
6. **Hamel ±15% 带** — LU 目标 CO/CO₂ = 0.13/0.11 → 合理区间约 **CO 0.11–0.15，CO₂ 0.09–0.13**；超出则属 **闭合类型** 问题，非细调。

---

## 6. 调参优先级（相对 `r*_scale`）

当新 validation case 组成偏离时，按序排查：

| 优先级 | 对象 | 示例杠杆 |
|--------|------|----------|
| P0 | 工况是否落在 §3 哪一格 | 氧化剂、ER、P、燃料水分 — 勿混用 LU scale |
| P1 | 边界 / 进风拓扑 | primary/secondary 分配、`gas_inlet_dense_frac`、steam 比 |
| P2 | 温度链 | `budget_progress` seed、围栏、床层/自由板能量、ξ≈0.6 峰 |
| P3 | Char / 热解时空 | VM 分区、R2/R4 环境（比表面积、T_solid） |
| P4 | 气相再平衡 | R8 驱动力（T、分压）；**诊断** 用 `r8_scale`，非默认 |
| P5 | CO 氧化强度 | R5 泡相/浓相相对供给；**诊断** 用 `r5_scale` |
| ✗ | 单 case 四旋钮 `{r2,r4,r5,r8}` | 仅 diagnostic JSON，**不得**作 Hamel 复现默认 |

**若做弱标定**：至少 2 case 联合、≤3 个有物理对象的因子（如 char 活性、WGSR 催化因子、O₂→dense 分配），held-out case 只验趋势。

---

## 7. 推荐多 case 验收矩阵

| 格子 | 建议 case | 验收重点 |
|------|-----------|----------|
| HTW Air+Steam 基准 | `HTW_WESSELING_1` | CO/CO₂/H₂ ±15%；轴向 T、ξ≈0.6 峰 |
| HTW O₂+Steam 对照 | `HTW_WESSELING_2` | 绝对 CO/H₂ ~3× LU；**CO/(CO+CO₂) 趋势** 不下降 |
| 工业尺度 O₂ | `HTW_BERRENRATH` | 平滑剖面、WGS 接近平衡 |
| 高水分燃料 | `VTT_PRESSURIZED_PEAT_6` | 低床温、高 CO₂ |
| 常压 / 无 steam | `VTT_ATMOSPHERIC_PEAT_D` | T_exit 低、自由板 CO 缓降 |
| 生物质 | `WSV400_A` | 低 CH₄、干燥区 |

**下一步 baseline（零 diagnostic scale）**：`python3 scripts/benchmark_phase2_multi_case_matrix.py` → `data/phase2_multi_case_matrix_latest.json`（LU + Wesseling 2 + VTT_PEAT_6 三分支趋势）。

---

## 8. 数据完整度（JSON 缺口）

| 类别 | 状态 |
|------|------|
| HTW 三案 | **97%+**；轴向 + 出口齐全 |
| VTT 13/14、WSV B–F | 出口干基 **多数缺失** |
| VTT_COAL_1 ER | 缺失 |
| O₂/Steam 蒸汽流量 | Wesseling 2、Berrenrath 蒸汽 **待补**（见 `_missing_data_summary`） |

更新 JSON 后须同步修订 **§2 表格**。

---

## 9. 与本仓库 Phase2 开发的关系

| 项 | 当前口径 | 本文定位 |
|----|----------|----------|
| 主验证 case | `CASE_HTW_WESSELING_1` | §3.1 / §4 LU 行 |
| Phase2 默认 | `budget_progress` seed、fixed50、vm_zone=0、`r5=0.25` | 数值/收敛诊断，**非** Hamel Table 5.x |
| `calibrate_phase2_co_co2_window.py` | diagnostic `r*_scale` 扫描 | 仅 **identifiability**；结果存 `data/phase2_co_co2_*`，不写入本 doc 为 SoT |
| 同构差距 | `hamel_isomorphic_gap_implementation_checklist.md` | 结构闭合 **先于** 本 doc §6 的 P4–P5 |

---

## 10. 维护

- **数据变更**：只改 `data/validation_cases.json`，再更新 §2、§8。  
- **新增 case**：补一行表格 + 判定 §3 分支 + 是否需扩展 §4 轴向模板。  
- **禁止**：将某次 `data/phase2_co_co2_calibration_*.json` 的 optimal scale 写进本文件作为「Hamel 参数」。

## 11. Multi-case 矩阵 baseline（2026-06-29）

脚本：`scripts/benchmark_phase2_multi_case_matrix.py`  
汇总：`data/phase2_multi_case_matrix_latest.json`

| Case | 收敛 | T_exit | CO / CO₂ (干基) | ref CO/CO₂ | 备注 |
|------|------|--------|-----------------|------------|------|
| LU | ✓ | 1112 K | 0.088 / 0.166 | 0.13 / 0.11 | T 在 ±15% 带内；组成仍偏 |
| W2 | ✗ | 1955 K | 0.00 / 0.653 | 0.404 / 0.283 | **禁止** LU 默认外推 |
| Peat6 | ✗ | 1054 K | 0.246 / 0.107 | 0.135 / 0.172 | T 接近；CO₂ 趋势反号 |

趋势门控（未收敛时仅诊断）：W2 CO 倍数、Peat CO₂ 倍数见 JSON `trends`。

---
