# Gibbs–动力学耦合框架

> 来源：BFB_TechSpec_v11 §5、§5.6；Hamel (1999) 附录 A1

## 概述

将 Gibbs 自由焓平衡与动力学模型耦合，实现：
1. 可逆反应连续驱动力 (1 − Q_p/K_eq)
2. 微量组分（H2S、SO2、COS、NH3、HCN、NO）Gibbs 最小化
3. 求解器层次中每次迭代调用 Gibbs

## 模块映射

| TechSpec §5.5 | 实现文件 |
|---------------|----------|
| equilibrium.py | `src/thermodynamics/equilibrium.py` |
| gibbs_minimizer.py | `src/thermodynamics/gibbs_minimizer.py` |
| minor_species.py | `src/thermodynamics/minor_species.py` |
| calc_gibbs_correction | `Cell.calc_gibbs_correction()` |
| get_element_release | `Cell.get_element_release()` |
| calc_minor_species_gibbs | `Cell.calc_minor_species_gibbs()` |

## 驱动力 (1 − Q_p/K_eq)

| 反应 | Q_p 形式 | 截断 |
|------|----------|------|
| R5 CO 氧化 | P_CO2²/(P_CO²·P_O2) | max(0, ·) |
| R7 CH4 重整 | (P_CO·P_H2³)/(P_CH4·P_H2O) | 不截断（可逆） |
| R8 WGSR | (P_CO2·P_H2)/(P_CO·P_H2O) | 不截断（可逆） |

## 微量组分 Gibbs

- **硫系**：H2S、SO2、COS
- **氮系**：NH3、HCN、NO
- **开关**：`Cell.use_gibbs_minor=True` 时替代 R9

## 数据流

```
residuals() → calc_hydrodynamics() → calc_exchange() → calc_reactions()
                                                              ↓
                                              [动力学 R1–R8, R10, R11]
                                                              ↓
                                              [若 use_gibbs_minor: calc_minor_species_gibbs]
                                                              ↓
                                              [R_gas_d 含 Gibbs 松弛源项]
```

## 元素来源

- `get_element_release('S')`：固体进料 sulfur_fraction + 气相 H2S/SO2/COS
- `get_element_release('N')`：固体进料 nitrogen_fraction + 气相 NH3/HCN/NO/N2

---

## Hamel 硫元素两步守恒（Hamel 论文 §4.2、Eq. A.28）

硫按**两步**从固相进入气相：

### Step 1：热解阶段分配

燃料总硫 $x_{S,wf}$ 分为：
- **挥发硫**：热解释放时进入气相
- **炭硫**：留在炭中，随 R1/R2/R3 炭消耗逐步释放

**守恒方程**：
$$x_{S,wf} \cdot \dot{m}_{B,wf} = \sum_{j} a_{S,j} \cdot \dot{N}_j \cdot M_S + \sum_{k} \dot{m}_{S,Koks,k}$$

### Step 2：Gibbs 平衡求 H2S/SO2/COS

以释放到气相的总硫为约束，用附录 A1 最小化 Gibbs 自由焓，得到 H2S、SO2、COS 分布。

### 分区效应

| 区域 | 气氛 | 主要硫产物 |
|------|------|------------|
| 燃烧区（床底） | 高 O2 | SO2 |
| 气化/自由板区（床上部） | 高 H2、低 O2 | H2S |

出口 H2S 典型 100–500 ppm（取决于燃料 S 含量，Table 7.1 如 Coal Case 1 约 0.51 wt-%）。

### 当前实现与 Hamel 的差异

| Hamel 模型 | 当前实现 |
|------------|----------|
| 挥发硫 + 炭硫 分步释放 | `get_element_release` 含挥发硫 + 炭硫 |
| 炭硫随 R1/R2/R3 消耗逐步释放 | 已计入 `R_solid` 对应的炭消耗 |
| 燃料 S 来自 Table 7.1 | `SolidProps.sulfur_fraction`、`sulfur_volatile_frac` |

### Hamel 附录说明（§8.1）

平衡法假设硫反应完美混合、无限停留时间，自由板区较冷时可能**低估** COS 或 SO2。
