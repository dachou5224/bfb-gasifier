# Hamel 能量守恒 parity 审计（2026-07-06）

## 范围

- Hamel (1999) Eq. (2.7) 全床能量平衡
- Hamel (1999) Eq. (2.8) 气/固焓流求和项
- 实现：`src/core/cell_balances.py`、`src/core/energy_balance_defaults.py`

## Eq. (2.8) 固体焓流

| 论文项 | Hamel 含义 | 实现 | 判定 |
|---|---|---|---|
| 气相 Σ N_j h_j | 物种摩尔流 × 比焓 | `calc_gas_enthalpy_flow` | 一致 |
| 固相 Σ m_jk (h_f + cp·ΔT) | 按组分/粒径 | `calc_solid_enthalpy_flow(..., solid_enthalpy_mode=hamel_eq28_component)` | **一致（默认）** |
| 混合比热 0.3/0.3/0.4 | 非论文式 | `legacy_cp_mix_030304` | audit 分支 |

### 修正项（2026-07-06）——其中第 1 条已于 2026-09-14 撤回

1. **生成焓索引（已撤回）**：当时把 `-15.866e6` 从 moisture 挪到 char。该数是液态水生成焓，不是石墨碳。二次诊断 `secondary_diagnosis_2026-09-14.md` §2 已用独立反应焓和 HHV 恒等式证伪。现行：char = 0，moisture = 液态水，VM 由干料 HHV 反算。
2. **比热**：默认按 char / VM(char cp) / moisture(4180) / ash 分别取值，不再用总量混合 cp。

## Eq. (2.7) 能量残差

| 论文项 | 实现 | 判定 |
|---|---|---|
| H_zu + H_rez + H_auf + H_ab − H_aus | 进出焓流组装 | 一致 |
| 工程热损比例 | `heat_loss_mode=inlet_fraction`（**legacy opt-in**） | 近似（`f·max(H_in,0)`）；H_in≤0 时 Q≈0 |
| 显式 Q_W | `heat_loss_mode=hamel_eq27_explicit`（**生产默认**） | `Q_loss = Q_wall_W + Q_heat_exchanger_W`；未设 `Q_wall_W_total` 时用燃料 HHV 代理 |

残差统一为：`H_in − H_out − Q_loss`（与 `H_in(1−f)−H_out` 在 inlet_fraction 下等价）。

## 测试

- `tests/test_cell_balances.py::test_energy_balance_zero_for_matched_inlet_and_outlet_streams`
- `tests/test_energy_balance_parity.py`（组分 cp、显式热损、热损比例线性）

## 2026-07-06 续：自由板轨迹

| 配置项 | 落点 |
|---|---|
| `resolve_freeboard_heat_loss_bundle(cfg)` | 自由板 ``heat_loss_frac`` + 显式 ``Q_wall_W`` / ``Q_heat_exchanger_W`` |
| `apply_freeboard_substep_gas_enthalpy_heat_loss` | `reaction.py` 背景剖面、`trajectory_hamel.py`、`segment.py` 子步焓扣除 |
| `freeboard_bridge` / `result_exit_summary` | 调用 `simulate_freeboard(..., heat_loss_mode=..., Q_wall_W=...)` |
| 显式 Q 分配 | 自由板 cell 均分 ``fb_q_wall`` / ``fb_q_hex``（``Reactor._apply_heat_loss_distribution``） |

## 2026-07-06 续：Phase2 探针（``scripts/audit_energy_heat_loss_mode_compare.py``）

| 项 | 说明 |
|---|---|
| 入口 | ``--fast --init-only``：3 bed + 2 fb，无 side blocks / staged Vorab |
| 校准 | Vorab 后按 cell ``f · H_in`` 求 implied ``|Q|``（``H_in`` 可为负，见 ``calc_energy_balance_H_in``） |
| 输出 | ``data/phase2_energy_heat_loss_mode_compare_latest.json`` |

**快探针结果（init-only，2026-07-06，修正前）**：signed ``Q ≈ -1217 kW`` 来自旧式 ``f·H_in``；已改为 ``f·max(H_in,0)``，负 ``H_in`` 时 implied Q → 0。

## 仍开放

- 化学反应热仍经源项 `R_gas`/`R_solid` 进入物种焓，未在 Eq.(2.7) 单独列出（与 Hamel 一致）。
- **Phase2 全焓流为负时**，默认 `inlet_fraction` 实际 **`Q_loss=0`**（`f·max(H_in,0)`）；探针用 `resolve_inlet_fraction_calibration_W()` 的 **`abs_H_in_audit_bridge`** 映射显式 `Q_W` 仅用于对比，**非生产默认**。

## 2026-07-07 续：探针修复与结论

| 项 | 说明 |
|---|---|
| 校准时机 | 先在 `inlet_fraction` 模式跑 NR，**求解后**再汇总 implied Q（修复“求解前 H_in 为负 → Q=0”） |
| 桥接口径 | `positive_H_in` 优先；若 `total_W≤1 W` 则 `abs_H_in_audit_bridge`（`f·|H_in|`） |
| 命令 | `python3 scripts/audit_energy_heat_loss_mode_compare.py --fast --max-global-iter 8` |
| 输出 | `data/phase2_energy_heat_loss_mode_compare_latest.json`（含 `calibration_mode`） |
| 判定 | **生产默认已切 `hamel_eq27_explicit`**；出口 T 正式 validation **待反应/能量校正后** |

## 2026-07-06 续：Reactor → Cell 透传

| 配置项 | 落点 |
|---|---|
| `energy_balance_heat_loss_mode` | `Cell.heat_loss_mode` |
| `solid_enthalpy_mode` | `Cell.solid_enthalpy_mode` |
| `Q_wall_W_total` / `Q_heat_exchanger_W_total` | 按 `H_bed`/`H_freeboard` 线性分配至 bed cell 的 `Q_wall_W` / `Q_heat_exchanger_W` |

实现：``Reactor._apply_heat_loss_distribution``、``_resolve_axial_explicit_heat_loss_watts``。
