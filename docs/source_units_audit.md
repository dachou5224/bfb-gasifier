# `src/` 单位与量纲自检说明

> 自动审查日期：与代码库同步维护；发现新混用请更新本文并修正源码注释。

## 1. 全局约定

| 量 | 约定单位 | 主要使用位置 |
|----|-----------|--------------|
| 温度 `T` | K | 全库 |
| 压力 `P` | Pa | `cell.P`、`rate_R*`、`equilibrium` |
| 气相摩尔流率 `N_*` | mol/s | `Cell.N_b`, `N_d`, 进料 |
| 固相质量流率 `m_solid*` | kg/s | `Cell.m_solid`, `m_solid_in`, `R_solid` |
| 气相反应源项体密度 `R_*`（动力学） | mol/(m³·s) | `calc_reactions` 内 `R_b`, `R_d` |
| 积分到 cell 的源项 `R_gas_*` | mol/s | `R_gas_b = R_b * V_b` |
| 炭异相速率 `r1`…`r4` | mol/(m²·s) | `char_reactions`；× 外表面积 → mol/s |
| 固相反应源项 `R_solid` | kg/s | 与固相守恒 `0 = zu + in + R - out` 一致 |

参考：`core/constants.py` 中 `T_REF = 298.15 K`（热力学参考）。

## 2. 已识别并已修正的问题

### 2.1 炭外表面积（`cell.py`）

- **错误**：用 `m_solid` [kg/s] 当颗粒质量算球数。
- **修正**：用乳化相持料量 `M_inventory ≈ rho_s * (1 - eps_mf) * V_d` [kg] 估算外表面积；`R_solid` 用 `-(Σ r_i)·M_C·A_k` [kg/s]。

### 2.2 能量衡算进流温度（`cell.py` + `reactor.py`）

- **错误**：所有流股用同一 `T`。
- **修正**：`T_zu_gas`、`T_zu_solid`、`T_in_gas`、`T_in_solid` 分设；`reactor` 传播上游温度。

### 2.3 固相显焓（`cell.py`）

- **错误**：`cp * T`。
- **修正**：`cp * (T - T_REF)`，与 `T_REF` 一致。

### 2.4 干燥/热解固相源（`cell.py`，本次审查）

- **错误**：`m_water_release`、`m_vm_release` 已由入口质量流率 [kg/s] 与无量纲进度相乘得到 **[kg/s]**，再除以 `tau_cell` 会变成 **[kg/s²]**，与 `R_solid` [kg/s] 不一致。
- **修正**：固相扣减项 **不再** 除以 `tau_cell`。`tau_cell` 仍仅用于 Agarwal/DAEM 子模型的时间尺度。

## 3. 需持续注意、当前视为一致的约定

| 主题 | 说明 |
|------|------|
| R8 WGSR `rate_R8` | `P_bar = P/1e5` 为 **bar**；文献若用 atm 需单独核对。 |
| `equilibrium.calc_reaction_quotient` | 分压由 `y * P` [Pa] 构成，与 `K_eq` 无量纲形式配套。 |
| `species.gas_diffusivity_correlation` | `101300/P` 为无量纲压力修正，**P 为 Pa**。 |
| `mass_transfer.calc_u_br` | `(P/p_ref)^(-0.15)`，**P 与 p_ref 同单位**。 |
| Gibbs `gibbs_minimizer` | `log(P/P0)`，**P0=101325 Pa**。 |

## 4. 辅助模块（非主回路或不同定义）

| 模块 | 说明 |
|------|------|
| `thermal/drying_rate_for_cell`、`devolatilization_rate_for_cell` | 返回 **基于 kg 燃料/kgdaf 的速率型指标** 或 `/tau` 的简化指标；**未**接入 `Cell.calc_reactions` 主路径时，不与 `R_solid` [kg/s] 混用。 |
| `physics/freeboard.py` | 工程近似，调用前确认速度与阻力公式中长度/速度单位。 |

## 5. 回归建议

```bash
cd bfb-gasifier
python3 tests/sanity_checks.py
pytest tests/ -q -m "not slow"
# 与 validation_cases.json 端到端对比（默认严格；未标定可设 BFB_RELAX_VALIDATION=1）：
pytest tests/test_table2_LU.py -q
python3 tests/test_table2_LU.py
```

---

*本文档由源码审查生成；若修改守恒方程或源项定义，请同步更新 §1–§2。*

**交叉引用**：`docs/missing_parameters_summary.md` 文首「相关文档」指向本文。
