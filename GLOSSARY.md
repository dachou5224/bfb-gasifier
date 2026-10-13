# Hamel 鼓泡流化床气化模型

本项目复现 Hamel（2001，VDI Fortschritt-Berichte Reihe 6 Nr. 469）的 cell 模型。术语以 Hamel 原文为准。

页码记作“PDF n / S. m”：n 为本地 PDF 的页序（从 1 起），m 为印刷页码。本 PDF 每页是两张印刷页的对开扫描，且有重复扫描（PDF 84、85 均为 S. 148–149）和缺页（S. 152–153），不能按固定偏移换算；页码一律以 `refs/hamel/index.json` 为准。本表页码已按该索引复核（2026-10-13）。

## 模型结构

**Zellenmodell（cell 模型）**：
把反应器沿高度划分为若干个串联的全混 Bilanzraum，对每个 cell 列写气体、固体与能量衡算，并联立求解的模型。PDF 15–19 / S. 10–19。
代码名：`cell_model`
_Avoid_：zone 模型、两区核、named_f 核

**Bilanzraum / Zelle i（cell）**：
cell 模型中编号为 i 的一个衡算空间。所有 cell 的衡算变量相同，衡算变量总数 = cell 数 × 每个 cell 的变量数（Eq. 2.3）；每个 cell 的变量数见 Eq. 2.1（PDF 17 / S. 15）。PDF 18 / S. 16。
代码名：`cell`
_Avoid_：格、zone、segment（用于床层 cell 时）

**Blasenphase（气泡相）**：
cell 内以无固体气泡形式上升的气相，有自己的气体衡算（Eq. 2.5）。PDF 18 / S. 16。
代码名：`bubble`
_Avoid_：b 相

**Suspensionsphase（悬浮相）**：
cell 内气体与床料、燃料颗粒组成的致密相，有自己的气体衡算（Eq. 2.4）；固体只存在于悬浮相（§2.2.2，气泡视为无固体）。PDF 18 / S. 16–17。
代码名：`suspension`
_Avoid_：乳化相、dense phase、d 相、emulsion

**Wirbelbett（床层）**：
流化床中自布风板到床面的致密区域。
代码名：`bed`

**Freiraum（自由空间）**：
床面以上、固体浓度低几个数量级的区域，用单独的流动模型描述，与床层通过床面上的物质流和能量流耦合。第 3.2 节，PDF 28 起 / S. 36 起。
代码名：`freeboard`
_Avoid_：自由板、freeboard 段

**Zyklon（旋风分离器）**：
把自由空间带出的固体分离并经回料送回床层的装置。回料是 cell 间的非相邻物流，使 Funktionalmatrix 出现非带状元素。PDF 20 / S. 20。
代码名：`cyclone`

## 衡算量

**Feststoffart j（固体种类）**：
参与化学反应的活性固体种类，例如 Koks、石灰石；每种可再分若干 Partikelgrößenklasse。§2.2.2，PDF 18 / S. 17。
代码名：`solid_species`

**Partikelgrößenklasse k（粒级）**：
固体种类 j 的离散粒径类别，固体衡算对每个 (j, k) 分别列写（Eq. 2.6）。粒级数可任意设定（Hamel 采用 Wirsum 1998 的扩展，§2.2，PDF 17 / S. 15）。Eq. 2.6 见 PDF 18 / S. 17。
代码名：`size_class`
_Avoid_：bin、档（用于正式命名时）

**Feststoffmassenstrombilanz（固体质量流衡算）**：
Eq. 2.6 的名称：对每个固体种类 j、粒级 k、cell i 列写的质量流平衡，各项都是质量流（kg/s）。名称见 Bild 2.2 流程图（PDF 17 / S. 14）；方程见 PDF 18 / S. 17，符号说明见 PDF 19 / S. 18。
代码名：`solid_mass_flow_balance`
_Avoid_：库存方程、holdup 方程

**Feststoffinventar（床料存量）**：
床层中的固体总量。其组成与粒径分布（或平均粒径）作为输入给定；床层高度也作为输入参数，因为实际运行中床高通过压差测量保持恒定。表 2.2 运行参数，PDF 16 / S. 13。
代码名：`bed_inventory`
_Avoid_：holdup（指床料总量时）

**cell 内固体质量（库存）**：
cell i 中固体种类 j、粒级 k 的质量（kg），归档代码记作 M 或 holdup，并假定流出流与之成正比（ṁ = K·M）。**待核实**：原文中尚未找到该量的叫法，也未确认它是否为衡算变量，还是 Eq. 2.6 只以质量流为变量、由流动模型给出各流项；见 `docs/adr/CANDIDATES.md` C2、C3。核实前暂用下列代码名，核实后按原文更名。
代码名（暂定）：`solid_mass_in_cell`
_Avoid_：holdup、M（作为标识符时）

**ṁ_auf / ṁ_ab（上行／下行固体流）**：
cell i 中颗粒类 (j, k) 向上、向下流出的质量流，相邻 cell 之间的固体返混由它们描述（Eq. 2.6）。PDF 19 / S. 18。
代码名：`m_up`、`m_down`
_Avoid_：auf/ab 流、upflow closure

**ṁ_zu / ṁ_aus / ṁ_rez（外部加入／流出系统／再循环流）**：
从外部进入 cell i 的流、从 cell i 离开整个系统的流、从其它 cell 回流进入 cell i 的流；气体对应 Ṅ_zu、Ṅ_rez。Eq. 2.6 符号说明，PDF 19 / S. 18。
代码名：`m_in_ext`、`m_out_sys`、`m_rec`

**Ṅ_ex,bd（相间交换流）**：
同一 cell 内气泡相与悬浮相之间组分 j 的交换摩尔流（Eq. 2.4、2.5）。PDF 18 / S. 16。
代码名：`n_exchange_bd`

**Energiestrombilanz（能量衡算）**：
对整个 cell（不分相）列写的能量平衡（Eq. 2.7）；物流携带的能量含标准生成焓（Eq. 2.8，Hamel 引用 Wirsum 1998）。PDF 19 / S. 19。
代码名：`energy_balance`

## 燃料转化

**Trocknung（干燥）**：
燃料颗粒中水分的蒸发，第 4.1 节。
代码名：`drying`

**Pyrolyse（热解）**：
燃料受热释放挥发分、生成 Koks、Teer 与热解气的过程，第 4 章，PDF 35 起 / S. 50 起。
代码名：`pyrolysis`
_Avoid_：devolatilization（作为模块名时）

**Flüchtige（挥发分）**：
热解中释放的物质总称。其在衡算中的归属见 `docs/adr/CANDIDATES.md` C1（待核实）。
代码名：`volatiles`
_Avoid_：VM（在代码标识符中）

**Koks / Restkoks（焦炭／残焦）**：
热解后残留的固体含碳产物，参与 R1–R4。
代码名：`char`
_Avoid_：炭（作为正式术语时）

**Teer（焦油）**：
热解释放的可凝有机物，参与 R10、R11。
代码名：`tar`

**Pyrolysegas（热解气）**：
热解释放的气体，其组成在第 4.3 节通过 Gibbs 自由焓最小化计算，组分含 CO、CO₂、CH₄、C₂H₄、NH₃、H₂S、H₂O、H₂ 等。PDF 44 / S. 69。
代码名：`pyrolysis_gas`

## 反应

编号以 Hamel 第 5 章开头的反应列表为准（PDF 45 / S. 71，已目视核对原页），见 `docs/adr/0001-reaction-numbering.md`。计量式按原文列表照录；Φ_c 为 R1 中决定 CO/CO₂ 产物比例的机理因子，R8 的平衡处理见 §5.2.4。

| Hamel | 反应 | 德文名 | 归档旧代码编号 |
|---|---|---|---|
| R1 | C + (1/Φ_c) O₂ → (2 − 2/Φ_c) CO + (2/Φ_c − 1) CO₂ | Verbrennung | R1 |
| R2 | C + H₂O → CO + H₂ | heterogene Wassergasreaktion | R2 |
| R3 | C + CO₂ → 2 CO | Boudouard-Reaktion | R4 |
| R4 | C + 2 H₂ → CH₄ | hydrierende Vergasung | R3 |
| R5 | CO + ½O₂ → CO₂ | Kohlenmonoxid-Oxidation | R5 |
| R6 | H₂ + ½O₂ → H₂O | Wasserstoff-Oxidation | R12 |
| R7 | CH₄ + 3/2 O₂ → CO + 2 H₂O | Methan-Oxidation | R6 |
| R8 | CO + H₂O → H₂ + CO₂ | homogene Wassergasreaktion | R8 |
| R9 | CH₄ + H₂O → CO + 3 H₂ | Methan-Reformierung | R7 |
| R10 | Teer + O₂ → CO + H₂ | Teer-Oxidation | R10 |
| R11 | Teer + H₂O → CO + H₂ + CH₄ | Teer-Umsetzung | R11 |
| — | H₂S 氧化 | （Hamel 反应列表中没有） | R9 |

_Avoid_：在新代码和文档中使用“仓库编号”“thesis R6”这类双重编号；旧编号只出现在本表和迁移记录中。

## 求解

**Vorabrechnung（预计算）**：
程序三大模块之一（另两个是 Konfiguration 与 Zellenmodell），计算那些只能整体描述、与 cell 划分无关的过程（如流动），在建立 cell 衡算方程之前完成；预计算与 cell 模型反复对齐，直到满足收敛判据。§2.1，PDF 15–16 / S. 11–12；流程见 Bild 2.2（PDF 17 / S. 14）。
代码名：`precalc`
_Avoid_：Vorab、init_precalc（作为概念名时）

**Abgleich（对齐）**：
Vorabrechnung 与 Zellenmodell 之间的外层迭代。PDF 16 / S. 12；Bild 2.2（PDF 17 / S. 14）。
代码名：`outer_iteration`
_Avoid_：outer refresh、Check1/Check2（作为正式术语时）

**Startwert（初值）**：
Newton 迭代所需的全部衡算变量初值；可按用户设定新生成，也可取自之前的模拟结果。PDF 16 / S. 12。
代码名：`initial_guess`
_Avoid_：x0 seed、Startwert stack

**Funktionalmatrix（雅可比矩阵）**：
衡算方程对衡算变量的导数矩阵。cell 一维串联使其呈块三对角结构，非相邻 cell 之间的物流（如旋风回料）会产生块外元素。PDF 20 / S. 20。
代码名：`jacobian`

**gedämpftes Newton-Verfahren（阻尼牛顿法）**：
Hamel 采用的求解方法：基于 Wozny（1983）的 Newton–Raphson 算法，并采用 Wirsum（1998）的阻尼变体。PDF 20 / S. 20。见 `docs/adr/0002-damped-newton-solver.md`。
代码名：`damped_newton`

## 已废弃的概念

以下叫法来自归档代码的双轨或分阶段开发，新主线中不再使用：

- thesis 模式 / CORE / engineering 模式 / 工程桥：新主线只有一种模型，即 Hamel 模型（CHARTER §6.2）。
- named_f、stage N、s 推进：属于归档的两区核同伦实验。
- fastox、温度栏（fence）、残氧种子：属于归档的求解修补；新主线只允许 CHARTER §5 所列、并登记 ADR 的求解路径辅助手段。
