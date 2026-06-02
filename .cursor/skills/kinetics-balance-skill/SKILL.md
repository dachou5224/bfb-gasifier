---
name: kinetics-balance-skill
description: Implements kinetics and conservation equations for the BFB 1D Cell model, including gas-phase molar balances, solid mass balances with size classes, and global energy balance. Use when implementing or modifying calc_gas_balance, calc_solid_balance, calc_energy_balance, or reaction source terms.
---

# Kinetics & Balance Skill (动力学与守恒方程)

## 使用场景

当需要在一维 BFB 模型中实现或修改以下内容时，应主动应用本 Skill：

- `Cell.calc_gas_balance()`：气泡相/悬浮相气体摩尔守恒
- `Cell.calc_solid_balance()`：固相质量守恒（含粒径类迁移）
- `Cell.calc_energy_balance()`：全局能量守恒方程
- 反应源项 `dotN_r`（气相）、`dotm_r`（固相）与守恒方程的耦合

本 Skill 的执行顺序固定为：

1. 先核对代码与 `docs/hamel_submodels/02_cell_balances_and_exchange.md`、`05_drying_and_daem_pyrolysis.md`、`06_kinetics_r1_r11_and_equilibrium_driving.md` 是否一致。
2. 若不一致，再回到 Hamel 原论文 Kapitel 5 的原始公式、方程号、表格裁决。
3. 只有在代码与文档一致后，才允许根据模拟偏差讨论调参。

禁止依据当前实现去反向补全文档；证据链必须保持 **论文 → 文档 → 代码**。

## 1. 气体摩尔守恒：双相一致性约束

### 1.1 方程结构约束

1. 对每个 cell、每个气体组分 `j`，必须分别在气泡相 (b) 与悬浮相 (d) 建立稳态守恒方程：
   - 悬浮相：对应 `techspec.md` (2-1)
   - 气泡相：对应 `techspec.md` (2-2)
2. 两相方程的结构应显式包含：
   - 入口/再循环项：`dotN_zu`, `dotN_rez`
   - 轴向传输项：`dotN_{*, i+1}` 与 `dotN_{*, i}`（上游/下游 cell）
   - 反应源项：`dotN_{r,*,j,i}`
   - 相间交换：`dotN_ex,bd,j,i`

### 1.2 质量交换项 N_dot_ex 的符号一致性

这是本 Skill 的 **硬约束**，实现 `calc_gas_balance()` 时必须满足：

1. 对任意 cell `i`、组分 `j`：
   - 在悬浮相方程中，`dotN_ex,bd,j,i` 作为 **流出项** 出现（负号）
   - 在气泡相方程中，`dotN_ex,bd,j,i` 作为 **流入项** 出现（正号）
2. 在代码层面，应避免独立构造两套相间交换量，而是通过同一数组/张量引用，以保证数值相等：

```python
dotN_ex_bd: npt.NDArray[np.float64]  # shape: (n_cells, n_species)

# 悬浮相：... - dotN_ex_bd[i, j]
# 气泡相：... + dotN_ex_bd[i, j]
```

3. 建议在单元测试中加入以下不变性检查：
   - 对每个 `(i, j)`，验证两相守恒方程中相间交换项系数符号相反、数值相等。

## 2. 固相质量守恒与粒径类向量化

### 2.1 基本方程要求

1. 固相质量守恒必须遵循 `techspec.md` 方程 (2-3)，对每个燃料类型 `j`、粒径类 `k` 考虑：
   - 入口与再循环：`dotm_zu`, `dotm_rez`
   - 轴向传输：`dotm_auf`, `dotm_ab`
   - 反应源项：`dotm_r`
   - 粒径类迁移：`dotm_left`, `dotm_right`（由颗粒收缩或生长导致）

### 2.2 向量化硬约束

实现 `calc_solid_balance()` 时，必须满足以下硬约束：

1. 对固相质量流 `m_solid[j][k]` 以及相关质量通量，**禁止**使用 Python 原生的多层 `for` 循环（如 `for j in fuels: for k in size_classes:`）在 Python 级别进行逐元素更新。
2. 必须使用 **NumPy 向量化** 或数组运算实现粒径类循环，例如：
   - 使用形状为 `(n_cells, n_fuels, n_sizes)` 的张量
   - 粒径类迁移通过切片与加减操作实现：

```python
# m_solid: (n_cells, n_fuels, n_sizes)
# dotm_left/right: 同形状或可 broadcast 的张量

m_solid[:, :, 1:] += dotm_left[:, :, 1:]    # 质量从更大粒径类迁入
m_solid[:, :, :-1] += dotm_right[:, :, :-1] # 质量从更小粒径类迁入
```

3. 若确需对粒径方向进行卷积/差分操作，应使用：
   - `np.roll`、切片算子或矩阵乘法，而非在 Python 层手写循环。

## 3. 全局能量守恒与焓的一致性

1. 能量守恒方程必须遵循 `techspec.md` (2-4)，并基于总焓流定义 (2-5)：
   - `dotH_i` 同时包含所有气相组分、固相燃料及惰性颗粒、以及可能的外加固相流。
2. 在实现 `calc_energy_balance()` 时：
   - 必须使用与质量守恒一致的摩尔/质量流率张量
   - 焓 `h` 必须包含 **标准生成焓**，以便反应热通过焓差自然体现
3. 建议将焓计算与守恒求和解耦：
   - 物性层：`species.get_enthalpy(T, phase)` 负责给出每个组分的焓
   - 守恒层：仅做 `sum(flow * h)`、`sum(mass_flow * h)` 的张量运算

## 4. 动力学源项与守恒方程耦合

1. 对非均相反应（炭反应）：
    - 速率常数采用 Arrhenius 形式，并考虑 SPM/SCM 模型
    - 生成的 `dotm_r` 必须在固相质量守恒和能量守恒中一致使用
    - 若 `R1–R4` 在代码与文档不一致，必须回到 Hamel 原论文 `5.1.2`~`5.1.5` 与对应表格核对后再改
2. 对均相反应（R5–R11）：
    - 使用动力学速率表达式，不可简单用平衡替代
    - 对 `dotN_r` 的实现，必须保证在多相（b/d）中的分配与反应机理一致
    - 若 `R5–R11` 在代码与文档不一致，必须回到 Hamel 原论文 `5.2.x` 的方程号与表格核对后再改
3. 若存在 Gibbs 自由焓平衡校验：
    - 不允许直接覆盖动力学解
    - 应通过限制接近平衡的方向或缩放速率的方式保持数值稳定与物理可接受性

## 5. 接口与测试建议

建议在 `core/cell.py` 中组织如下接口：

```python
class Cell:
    def calc_gas_balance(self) -> None:
        """更新气泡相与悬浮相的组分摩尔流率。

        要求：
        - 使用同一 dotN_ex_bd 张量在两相方程中以相反符号出现
        - 保证每个 cell、每个组分的守恒方程满足数值平衡（残差用于迭代）
        """

    def calc_solid_balance(self) -> None:
        """以向量化方式更新各燃料、各粒径类的质量流率。

        禁止使用 Python 原生嵌套循环遍历 size classes。
        """

    def calc_energy_balance(self) -> None:
        """基于质量流率与焓计算全局能量守恒残差。"""
```

在测试层面，建议增加：

- 对任意静态状态：两相中相间交换项在整体系统中对每个组分的净贡献为零
- 对固相：总质量（轴向积分 + 所有粒径类求和）仅由进出流与反应决定，与离散细节无关
- 对能量：在关闭化学反应与外壁换热时，验证数值解满足能量平衡
