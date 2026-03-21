# Gas Species 与 Tar 定义（species）

## 1. 气相物种数量与清单

当前模型固定采用 **11 个气相物种**：

1. `CO`
2. `CO2`
3. `H2`
4. `H2O`
5. `CH4`
6. `O2`
7. `N2`
8. `H2S`
9. `NH3`
10. `TAR1`
11. `TAR2`

说明：
- `NH3` **始终保留**，用于氮平衡与后续 `NOx` 前驱分析。
- `C2H4` 在实验章节中有测量，但在当前反应网络中通常并入“高碳烃/tar 组”而不单独建方程。

## 2. Tar（Teer）定义与代理组分

tar 定义为燃料热解释放的高碳烃复杂混合物。为便于计算，采用**两组分代理模型**。

### 2.1 燃料相关的代理分子

- **煤（Rhenish Brown Coal）**  
  `TAR1 = C6H6`（Benzene），`TAR2 = C10H8`（Naphthalene），目标 H/C 约 `0.856`。
- **生物质（Wood/Sawdust）**  
  `TAR1 = C10H8`（Naphthalene），`TAR2 = C16H34`（Hexadecane），目标 H/C 约 `1.364`。

结论：tar 代理分子由 feedstock 唯一决定；`NH3` 不受 feedstock 切换影响。

### 2.2 H/C 保真约束

默认假设：释放 tar 的摩尔 H/C 与原始有机燃料（waf）保持一致（无额外实验数据时）。

## 3. Tar 反应化学计量（R10/R11）

### R10（Tar Oxidation）

`C6H6 + 3 O2 -> 6 CO + 3 H2`  
`C10H8 + 5 O2 -> 10 CO + 4 H2`  
`C16H34 + 8 O2 -> 16 CO + 17 H2`

### R11（Tar Steam Reforming）

`C6H6 + 9/2 H2O -> 9/2 CO + 9/2 H2 + 3/2 CH4`  
`C10H8 + 23/3 H2O -> 23/3 CO + 7 H2 + 7/3 CH4`  
`C16H34 + 21/2 H2O -> 21/2 CO + 33/2 H2 + 11/2 CH4`

## 4. Gibbs 微量组分（thermodynamics 模块用）

当 `Cell.use_gibbs_minor=True` 时，H2S、SO2、COS、NH3、HCN、NO 由 Gibbs 自由焓最小化求解。`species.py` 中已补充 SO2、COS、HCN、NO 的 NASA 系数及 `get_atom_count(species, element)`，供 `thermodynamics.gibbs_minimizer` 使用。这些物种的热力学数据不改变主守恒方程中的 11 物种状态向量。

**硫元素两步守恒（Hamel §4.2）**：燃料硫分为挥发硫（热解释放）与炭硫（随 R1/R2/R3 炭消耗释放）。`SolidProps.sulfur_volatile_frac` 控制挥发硫比例（煤典型 0.3–0.7）。

## 5. 与代码实现的对应关系

- `src/core/species.py`
  - 定义 `GAS_SPECIES`（11 物种）
  - 定义 `NH3`、SO2、COS、HCN、NO 热力学参数（Gibbs 用）
  - 定义 `get_atom_count()`（Gibbs 原子数矩阵）
  - 定义 `TAR1/TAR2` 与代理分子映射、H/C 配比求解
- `src/kinetics/tar_reactions.py`
  - 定义 R10/R11 的代理分子计量
  - 提供 `TAR1/TAR2` 显式消耗的等效计量
- `src/thermodynamics/`
  - `equilibrium.py`：K_eq、Q_p、驱动力 (1−Q_p/K_eq)
  - `gibbs_minimizer.py`：附录 A1 Gibbs 最小化
  - `minor_species.py`：硫系/氮系微量组分分布
- `src/core/cell.py`
  - 将 tar 反应速率映射为与 `GAS_SPECIES` 一致的源项向量
  - `calc_minor_species_gibbs()`、`get_element_release()`（Gibbs 集成）

## 6. 物化性质方法口径（与论文一致）

论文没有逐项给出完整 NASA/JANAF 系数表，而是给出函数依赖关系并默认采用标准热力学数据源。当前实现采用“**NASA 多项式 + 经典工程关联**”的组合，与该口径兼容：

- **焓/热容（`h_j`, `c_p`）**
  - 用 NASA 7 系数计算 `Cp(T), H(T), S(T)`，且 `H(T)`包含标准生成焓。
  - 这与“温度相关热力学函数 + 标准生成焓”的论文写法一致。
- **气体密度 `rho_g`**
  - 理想气体：`rho_g = P * M_g / (Rg * T)`。
- **气体粘度 `mu_g`**
  - 幂律：`mu_g = mu_g,20 * (T/293.15)^n`，默认 `n=0.7`。
- **扩散系数 `D_g`**
  - 按 Eq.5.20：`D_g = 3.13e-4 * (T_m/1500)^1.75 * (101300/P)`。

对应代码函数位于 `src/core/species.py`：
- `mean_molar_mass(...)`
- `gas_density_ideal(...)`
- `gas_viscosity_power_law(...)`
- `gas_diffusivity_correlation(...)`

