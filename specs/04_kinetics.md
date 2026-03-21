# 04 化学动力学

所有速率常数必须调用 `src/kinetics/arrhenius.py` 中的工厂函数。

## R1–R4 异相反应（炭）

| 反应 | 模型 | 速率函数 | 来源 |
|------|------|----------|------|
| R1 炭燃烧 | SPM（k_l 串联） | `k_hobbs` | Hobbs et al. (1992) |
| R2 炭+H2O | SCM | `k_hobbs` | — |
| R3 炭+H2 | SCM | `k_hobbs` | — |
| R4 Boudouard | Weeda / Langmuir-Hinshelwood | `k_hobbs` | Weeda (1995) |

> k_b、k_c（Boudouard）指数均为**正**（吸附热）。

## R5–R9 均相气相反应

| 反应 | 说明 | 驱动力 | 速率函数 |
|------|------|--------|----------|
| R5 CO 氧化 | 气泡相/悬浮相表达式不同；悬浮相含 `C_H2O^0.5` 项 | max(0, 1−Q_p/K_eq) | `rate_R5_bubble`, `rate_R5_suspension` |
| R6 CH4 氧化 | — | — | `rate_R6` |
| R7 水蒸气重整 | — | (1−Q_p/K_eq)，可逆 | `rate_R7` |
| R8 WGSR | Chen (1987)，含煤灰催化修正 | (1−Q_p/K_eq)，可逆 | `rate_R8` |
| R9 H2S 氧化 | 当 `use_gibbs_minor=False` 时使用；否则由 Gibbs 替代 | — | `rate_R9` |

**Gibbs 微量组分**：当 `Cell.use_gibbs_minor=True` 时，H2S/SO2/COS、NH3/HCN/NO 由 Gibbs 自由焓最小化求解（`thermodynamics.minor_species`），替代 R9。

**K_eq 与驱动力**：R5/R7/R8 的平衡常数与反应商由 `src/thermodynamics/equilibrium.py` 提供（`get_K_eq`、`calc_reaction_quotient`、`calc_gibbs_driving_force`）。R5 用 max(0, 1−Q_p/K_eq)；R7、R8 用 (1−Q_p/K_eq)，允许负值表示逆反应。

## R10–R11 焦油反应（两组分 tar 代理）

来源：Hamel (1999) §5.2.6，Eq.5.59，Table 5.4（pp.95–98）；组分与化学计量见 `specs/species.md`。

| 反应 | 相别 | 模型 | 速率函数 | 备注 |
|------|------|------|----------|------|
| R10 tar 氧化 | 气泡 / 悬浮 | 均相气相，**两相同一热力学** | `rate_R10(T, C_tar, C_O2, P, fuel_type)` | Eq.5.59：\(k_{10}\exp(-E/R_gT)\,T\,P^{0.3}\,C_{\mathrm{tar}}^{0.5}\,C_{O_2}\)；芳香 vs 烯烃/烷烃按代理分率加权 |
| R11 裂解/重整 | 气泡 | Serio et al. (1987) 均相热裂解 | `rate_R11_bubble(T, C_tar)` | 一级对 \(C_{\mathrm{tar}}\)；\(k_0=5.42\times10^4\,\mathrm{s^{-1}}\)，\(E=100.5\,\mathrm{kJ/mol}\) |
| R11 裂解/重整 | 悬浮 | Corella et al. (1991) 催化 | `rate_R11_suspension(T, C_tar, rho_cat)` | \(k_0=0.7\,\mathrm{m^3/(s\cdot kg_{Kat})}\)，\(E=63.1\,\mathrm{kJ/mol}\)，**× 局部催化剂质量密度** \(\rho_{\mathrm{cat}}\)（`Cell._catalyst_bulk_density()`） |

- **浓度**：\(C_{\mathrm{tar}} = C_{\mathrm{TAR1}} + C_{\mathrm{TAR2}}\)（mol/m³）。
- **兼容**：`rate_R11` 等价于 `rate_R11_bubble`（旧 `C_H2O` 参数已忽略）。
