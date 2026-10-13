# Hamel Parity 审计 — hydro 非 thesis 水力链 — 2026-07-06

**证据链**：`docs/hamel_submodels/03_hydrodynamics_core_chain.md` → `src/core/hydrodynamics_defaults.py` → `ReactorConfig` / `Cell` / `calc_cell_hydrodynamics`

**判定口径**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离`

---

## 1. 结论摘要

| 项 | 改前 | 改后 |
|----|------|------|
| `Reactor` 非 thesis 主链 | 已通过 `ReactorConfig` 默认走 Hamel 链 | 不变；常量集中到 `hydrodynamics_defaults.py` |
| 独立 `Cell` / `calc_cell_hydrodynamics` 默认 | **legacy**：`mori_wen` + `u_d=current` | **Hamel 主链**（与 thesis 相同策略名） |
| `thesis_mode` 水力覆盖 | 硬编码字符串重复 | 调用 `apply_hamel_hydrodynamics_to_cell`（安全锁） |
| Mori-Wen / `current` u_d | 隐式默认 | **显式 audit 分支**（须传参 opt-in） |

**总体**：非 thesis 与 thesis 在**床层水力链**上现已同构；`thesis_mode` 额外差异仅为 `holdup_transport` 固相、`exact_hamel` 自由板等（非本 scope）。

---

## 2. Hamel 主链对照表

| Hamel 方程 | 含义 | 策略键 | 实现 | 判定 |
|------------|------|--------|------|------|
| Eq. (3.11)–(3.12) | 悬浮相速度 u_d | `wein_1992_eq312` | `cell_hydrodynamics._resolve_u_d_closure` | **一致** |
| Eq. (3.13)–(3.15) | ψ_b 气泡相互作用 | `wein_1992` | `bubble_interaction_factor` | **一致** |
| Eq. (3.41) | dd_b/dh ODE | `heinbockel_eq341` | `integrate_bubble_diameter` | **一致** |
| Eq. (3.43) | u_b 加压修正 | `heinbockel_eq343` + `v_b=1` | `bubble_rise_velocity` | **一致** |
| Eq. (3.44) | ε_b 可见气泡率 | — | `calc_visible_bubble_fraction(u0,u_b,u_d)` | **一致** |
| Eq. (3.50) | K_bd | — | `calc_kbd` | **一致** |
| Mori-Wen d_b(h) | 常压经验式 | `mori_wen` | `mori_wen_bubble_diameter` | **有意偏离**（legacy audit） |
| `current` u_d | 旧工程闭合 | `current` | `_resolve_u_d_closure` | **有意偏离**（legacy audit） |

---

## 3. 默认配置锚点

```python
# src/core/hydrodynamics_defaults.py
HAMEL_HYDRO_U_D_CLOSURE = "wein_1992_eq312"
HAMEL_HYDRO_BUBBLE_DIAMETER_MODEL = "hilligardt_ode"  # 内部 Heinbockel Eq.3.41 ODE
HAMEL_HYDRO_BUBBLE_VELOCITY_STRATEGY = "heinbockel_eq343"
HAMEL_HYDRO_BUBBLE_ODE_STRATEGY = "heinbockel_eq341"
```

---

## 4. 仍开放（非本项）

- holdup `transport_mask` 等 guardrails（默认已是 `holdup_transport`，见 `hamel_parity_defaults.py`）
- 自由板 `exact_hamel` 已是 Reactor / `simulate_freeboard` 默认（2026-07-07 P1）
- `integrate_bubble_diameter(..., method="mori_wen")` 低层 API 默认仍为 Mori-Wen（仅 mori_wen 分支调用）
