---
name: hydrodynamics-skill
description: Implements hydrodynamics for the bubbling fluidized bed gasifier, including minimum fluidization velocity, bubble dynamics, and K_bd inter-phase mass transfer. Use when implementing or modifying calc_hydrodynamics, K_bd, bubble velocity, or U_mf-related logic in the BFB 1D model.
---

# Hydrodynamics Skill (BFB 流体力学子模型)

## 使用场景

当需要在 BFB 一维模型中实现或修改以下内容时，应主动应用本 Skill：

- `calc_hydrodynamics()`、`calc_exchange()` 等与两相流体力学相关函数
- 最小流化速度 `u_mf`、气泡速度 `u_b`、气泡直径 `d_b(h)` 等气泡动力学
- 相间传质系数 `K_bd`、气泡体积分数、慢泡/快泡判别逻辑

本 Skill 的执行顺序固定为：

1. 先核对代码与 `docs/hamel_submodels/03_hydrodynamics_core_chain.md`、`00_readme_and_citation_rules.md` 是否一致。
2. 若不一致，再回到 Hamel 原论文对应页码、方程号、表格裁决。
3. 只有在代码与文档一致后，才允许根据模拟偏差讨论调参。

禁止把当前实现当作文档来源反向补写公式；证据链必须保持 **论文 → 文档 → 代码**。

## 1. 最小流化速度与慢泡/快泡判别

1. 必须根据 Ergun 方程计算最小流化条件下的 `Re_mf` 和 `u_mf`，使用 `Ar`–`Re_mf` 的关系（参见 `techspec.md` 中 4.1 节的方程）。
2. 所有后续气泡相关量（如 `u_b`、`d_b`、`ε_b`）都必须显式依赖 `u_mf`，禁止写死常数或与 `u_mf` 脱钩的经验式。
3. 慢泡/快泡状态的判别必须严格使用：
   - \(\alpha_b = U_b / U_{mf}\)
   - `alpha_b < 1` → 慢泡（slow bubbles）
   - `alpha_b > 1` → 快泡（fast bubbles）
4. 在代码结构上，应把判别逻辑封装为辅助函数，避免在多个地方重复硬编码阈值：

```python
def classify_bubble_regime(u_b: float, u_mf: float) -> str:
    """根据 U_b / U_mf 判别慢泡或快泡状态。"""
    alpha_b = u_b / u_mf
    return "slow" if alpha_b < 1.0 else "fast"
```

## 2. 气泡动力学实现要求

1. 气泡上升速度必须采用 Hilligardt 型表达式：
   - \(u_b = \psi_b (u_0 - u_{mf}) + u_{b,i}\)
   - `ψ_b` 默认为约 0.76（工业分布板），但应保留为可配置参数。
2. 气泡直径 \(d_b(h)\) 必须通过微分方程（`dd_b/dh`）数值积分得到，不允许直接写死经验极限值：
   - 实现 ODE 积分时，应保证轴向离散步长与 `Cell` 网格一致。
3. 若实现包含气泡寿命 `λ_b`、气泡相体积分数 `ε_b` 等量，需保证：
   - 所有量的单位与 `techspec.md` 中保持一致
   - 通过 assert 或类型注解校验 Shape，与轴向 cell 数量严格一致

## 3. 相间传质系数 K_bd 的实现规范

### 3.1 方程与物理映射

1. `K_bd` 实现必须遵循 `techspec.md` 4.3 节的 Sit & Grace 混合模型：
   - \(K_{bd} = \frac{3 u_{br}}{2 d_b} + \sqrt{\frac{144 D_g \epsilon_{mf} u_b}{\pi d_b^3}}\)
2. 代码必须完整反映上述两项物理含义：
   - 第一项：对流换气，由修正速度 `u_br` 和气泡直径 `d_b` 决定
   - 第二项：扩散贡献，由气体扩散系数 `D_g`、`ε_mf`、`u_b` 和 `d_b` 决定
3. 对流项中的 `u_br` 必须包含压力修正：
   - \(u_{br} = n_b \cdot u_d \cdot (P / P_0)^{-0.15}\)
   - 禁止忽略 `(P / P0)^(-0.15)` 项或将其硬编码为 1。

### 3.2 与两相理论的一致性

1. 相间传质只发生在气泡相与悬浮相之间，`K_bd` 仅用于构造：
   - \(\dot{N}_{ex,bd,j,i} = K_{bd,i} V_{b,i} (C_{j,b,i} - C_{j,d,i})\)
2. 所有粒子相关参数（如 `ε_mf`、颗粒直径、气体扩散系数 `D_g`）应来自统一的物性/建模模块，而不是在 hydrodynamics 模块中重复定义。
3. 若模型存在不同高度区段或设备段（炉膛、自由板区、连接管道等），`K_bd` 的计算必须与对应段的局部条件（压力、温度、气速）一致。

## 4. 接口设计与示例

建议将流体力学相关逻辑集中在 `physics/hydrodynamics.py` 中，并通过清晰的函数接口暴露给 `Cell`：

```python
from typing import Protocol
import numpy as np
import numpy.typing as npt

class HydrodynamicsModel(Protocol):
    def calc_u_mf(self, rho_g: float, rho_s: float, dp: float, mu_g: float) -> float: ...
    def calc_bubble_velocity(self, u_0: float, u_mf: float, u_b0: float, psi_b: float = 0.76) -> float: ...
    def calc_kbd(
        self,
        u_b: float,
        d_b: float,
        D_g: float,
        eps_mf: float,
        u_d: float,
        n_b: float,
        P: float,
        P0: float = 1.0e5,
    ) -> float: ...
```

实现时应：

- 使用 SI 单位并在函数签名附近注明单位
- 使用断言检查输入参数范围（如 `d_b > 0`, `P > 0`）
- 在 docstring 中注明公式来源，例如：
  - `Ref: Hamel (1999) Gleichung 3.24, 3.44, 3.50`
