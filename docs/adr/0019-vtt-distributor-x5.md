---
status: proposed
---

# VTT 布风器的出口孔数 n_0 与出口高度 h_A（阻塞项 X5）

> **状态：提议中。** 本 ADR 裁定 PR #71 审查发现的阻塞项 X5：式 (3.37) 的出口孔数 n_0 与式 (3.39) 的出口高度 h_A，VTT 加压装置在 KURKELA_1992 中没有给出。依据为 KURKELA_1991（T15 已按 CHARTER §4.2 用于补装置细节）与物理定义。按 CHARTER v0.10 §8，审查方审查后批准，目标所有者可否决。

## 背景

ADR 0014 认定 VTT 的 “multiorifice plate distributor” 属于技术布风器，起始气泡直径取式 (3.36)：d_b,0 = 1,3·(V̇_0²/g)^0,2，其中 V̇_0 = u_0·A/n_0（式 3.37）；开始长大的高度取式 (3.39)：h_0 = h_A + d_b,0/2（S. 29）。T15 记录：KURKELA_1992 只有布风器类型；KURKELA_1991 印刷页 8 另有孔板参数，但该报告的竖向尺寸与 1992 年论文不一致，当时未录入（`data/cases/vtt_pressurized_gaps.md`）。

KURKELA_1991 印刷页 8（总司令 2026-10-18 在页面图像上核对）：

> Primary air and a small amount of steam is introduced into the reactor through a multiorifice plate distributor. Two different distributor plates can be used depending on the required range of fluidizing air flow rate. The first plate is a 10 mm thick slightly conical plate with 1 mm holes and an open area of 0.62 % of the reactor cross-sectional area. The other air distributor is a 10 mm thick horizontal plate with 2 mm holes and a total open area of 1.7 % of the reactor area.

## 能否用于 1988–1990 年的装置

- T15 不录入 KURKELA_1991 的数，理由是竖向尺寸（反应器高、进料口高、二次风高）与 1992 年论文不一致。布风板是水平面上的部件，两篇的水平尺寸（床径 15 cm、自由空间径 25 cm）一致，布风器名称（multiorifice plate distributor）一字不差；KURKELA_1991 印刷页 5 写明这套 PDU 自 1988 年起使用。
- 因此本 ADR 采用 KURKELA_1991 的孔板参数（解释性选择，中等依据）；竖向尺寸仍不采用，T15 的结论不变。

## 裁定

### 1. h_A = 0（物理定义）

式 (3.39) 的 h_A 是出口孔高出布风板的高度，对应带风帽或喷嘴的布风器。两块孔板的孔都开在 10 mm 厚的板上，出口就在板面，h_A = 0，h_0 = d_b,0/2。

### 2. n_0 由开孔率与孔径换算

开孔率按原文是相对反应器截面（直径 D = 0,15 m）：n_0 = φ·(D/d_h)²。

| 孔板 | d_h | φ | n_0 |
|---|---|---|---|
| 第一块（微锥形） | 1 mm | 0,62 % | 139,5 |
| 第二块（水平） | 2 mm | 1,7 % | 95,6 |

φ 只有两位有效数字，n_0 按换算值使用，不取整。

### 3. 用哪一块板

原文只写两块板“按所需的流化风量范围”选用，没有写褐煤试验用了哪一块，也没有给出流量分界。这一选择无法由原文或物理定义确定。

**裁定（解释性选择，首次运行前冻结）**：基准运行取第一块板（d_h = 1 mm，φ = 0,62 %，n_0 = 139,5）。理由仅是原文的列举次序；没有更强的独立证据，不以出口结果选择（ADR 0013）。

敏感性：d_b,0 ∝ n_0^−0,4，两块板的 d_b,0 相差 (139,5/95,6)^0,4 = 1,16 倍，只改变布风板附近的起始气泡直径与 h_0；气泡沿高度按式 (3.33) 长大，起始值的影响向上减弱。首次运行报告须在相同冻结设置下给出第二块板（n_0 = 95,6）的 7 个验收量及差值，只报告，不据此改选；若某一验收量的差值超过其容差的一半，报告中单列，交目标所有者知悉。

### 4. 第一块板为微锥形

式 (3.36)–(3.39) 不含板的倾角；板厚 10 mm 与微锥形不进入模型。中心 38 mm 排灰管（KURKELA_1991）不计入开孔率（原文开孔率按反应器截面给出），也不改变 n_0。

## 实现约束

- `data/cases/vtt_pressurized.json` 增加 `distributor_orifice_count`（139,5，`derived`，公式 φ·(D/d_h)²，来源 KURKELA_1991 印刷页 8）与 `distributor_outlet_height`（0，`derived`，依据本 ADR 第 1 项），作为 T19 函数的输入。录入由数据卡或 T25 完成，本 ADR 不改数据文件。
- 第二块板的 n_0 = 95,6 只用于敏感性运行，不进入工况文件的主值。
- 代码与数据注明 `解释性选择：ADR 0019`（第 2、3 项）与 `ADR 0019`（第 1 项）。

## 审查记录
