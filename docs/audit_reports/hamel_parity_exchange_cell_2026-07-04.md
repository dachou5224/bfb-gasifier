## Hamel Parity 审计 — exchange-cell — 2026-07-04（Eq.(2.7)/(2.8) 于 2026-07-06 更新）

**证据链**：hamel_submodels/02 → checked_spec → src/  
**能量专项 SoT**：`docs/audit_reports/hamel_parity_energy_2026-07-06.md`

| 方程/表 | 分册位置 | 论文式（摘要） | 代码位置 | 判定 | 证据 |
|---------|----------|----------------|----------|------|------|
| Eq. (2-4) 悬浮相气体摩尔守恒 | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L18–28 | `0 = N_zu,d,j,i + N_rez,d,j,i + N_d,j,i+1 + N_r,d,j,i - N_d,j,i - ex_bd,j,i` | `src/core/cell_balances.py:calc_gas_balance_residual` L27–64；`src/core/cell.py:calc_gas_balance` L509–526 | 一致 | 代码残差 `res_d = N_zu_d + N_rez_d + N_d_in + R_gas_d - N_d + N_ex`，其中 `N_ex = -ex_bd`（docstring L45–47）；形状断言 `(N_GAS,)`；单位 mol/s。 |
| Eq. (2-5) 气泡相气体摩尔守恒 | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L24–28 | `0 = N_zu,b,j,i + N_rez,b,j,i + N_b,j,i+1 + N_r,b,j,i - N_b,j,i + ex_bd,j,i` | `src/core/cell_balances.py:calc_gas_balance_residual` L27–64；`src/core/cell.py:calc_gas_balance` L509–526 | 一致 | 代码残差 `res_b = N_zu_b + N_rez_b + N_b_in + R_gas_b - N_b - N_ex`（L63），`N_ex = -ex_bd`，符号与论文一致；两相耦合守恒。 |
| 相间交换项 `ex_bd` | 嵌入 Eq. (2-4)/(2-5) | `ex_bd,j,i` 为正表示从悬浮相到气泡相 | `src/core/cell_hydrodynamics.py:calc_phase_exchange` L375–390；`src/core/cell.py:calc_exchange` L275–282 | 一致 | 实现为 `N_ex = K_bd * V_b * (C_b - C_d)`（L390），`K_bd` 来自 Eq. 3.50，`V_b = eps_b * V_cell`；在气体守恒中映射为 `N_ex = -ex_bd`，符号正确。 |
| Eq. (2-6) 固相质量守恒 | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L30–40 | 稳态固相质量衡算含 auf/ab、反应、迁移 | `cell_balances.py:calc_solid_balance_residual`；`calc_size_migration` | **部分一致（P1 默认 holdup_transport + per_class_holdup）** | 生产默认 `holdup_transport` + `per_class_holdup`；legacy opt-in。holdup NR guardrails 为 B-tier accepted gap。 |
| Eq. (2-7) 整单元能量衡算 | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L42–46 | 整单元能流平衡，含外部输入、再循环、相邻单元通量、系统流出、壁面换热与内置换热面换热 | `src/core/cell_balances.py:calc_energy_balance_residual` L398–514；`calc_energy_balance_H_in` L341–395；`src/core/energy_balance_defaults.py`；`src/core/cell.py:calc_energy_balance`；`src/core/reactor.py:_apply_heat_loss_distribution` | 部分一致 | **进出焓流结构与论文一致**（`H_in − H_out − Q_loss`）。默认 `heat_loss_mode=inlet_fraction`：`Q_loss = heat_loss_frac · max(H_in, 0)`（负 `H_in` 时不增焓）。可选 `hamel_eq27_explicit`：`Q_wall_W + Q_heat_exchanger_W`（Reactor 按床高/自由板高度分配）。旋风/返料腿仍用 `_calc_side_temperature_closure_residual` 温度贯通闭合。详见 energy 审计。 |
| Eq. (2.8) 焓流求和表达 | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L42–46 | 与物流绑定的焓流求和：气体项 + 固体项 | `src/core/cell_balances.py:calc_gas_enthalpy_flow` L255–274；`calc_solid_enthalpy_flow` L277–338 | **一致（默认）** | 气相 `N·h_vec`（生成焓 + 显焓）一致。固相默认 `solid_enthalpy_mode=hamel_eq28_component`：char / VM / moisture / ash 分别取 cp 与生成焓（char 生成焓索引已修正）。审计分支 `legacy_cp_mix_030304` 保留 0.3/0.3/0.4 混合比热。docstring 已对齐 Eq. (2.8)。 |
| 连接矩阵/流分类 (`zu`/`rez`/`auf`/`ab`) | `docs/hamel_submodels/02_cell_balances_and_exchange.md` L49–52；gap §0.5 | Verbindungsmatrix 定义 `zu`（给料）、`rez`（再循环）、`auf`（上行）、`ab`（下行） | `src/core/connectivity_graph.py` L16–136；`src/core/connectivity_graph.py:build_connectivity_graph_from_reactor_params` L208–516；`src/core/connectivity/routing.py:_apply_bed_upstream_inlet` L56–95；`src/core/connectivity/routing.py:apply_recirculation_from_connectivity_graph` L167–197 | 一致 | 图边显式标注 `stream_class`（L39）；`auf/ab` 路由床层轴向气固耦合（routing.py L56–95），`rez` 路由返料腿回床层（L167–197）。gap checklist §0.5 判定 VTT 单回路「一致：unexpected_blocks==0」。 |

**判定枚举**：`一致` | `部分一致` | `不一致` | `缺失` | `有意偏离（见 gap §x）`

**自动化**：audit_hamel_consistency.py → PASS

**开放项**：
1. holdup NR guardrails 与 char/ash transport mask 已文档化为 **B-tier accepted gap**（`holdup_transport_char_ash_mask` / `_apply_holdup_transport_numerical_guardrails`）。
2. ~~粒径迁移 `M_inv_total` 代理~~ → 生产默认 **`per_class_holdup`**（2026-07-07）；`total_bed_proxy` opt-in。
3. ~~默认 `inlet_fraction`~~ → 生产默认 **`hamel_eq27_explicit`**（2026-07-07）；Fig.7 出口 T validation **待反应/能量校正**。
4. 旋风/返料腿 side-block 仍用温度贯通闭合，非完整 Eq.(2.7) 能量方程。
5. 化学反应热经 `R_gas`/`R_solid` 进入物种焓，未在 Eq.(2.7) 单独列出（与 Hamel 一致，非缺陷）。

**P1 已闭合（2026-07-07）**：bed 默认 `holdup_transport`（`src/core/hamel_parity_defaults.py`）；自由板默认 `exact_hamel`（`ReactorConfig` + `simulate_freeboard`）。

**建议下一步**：
1. Phase2 全 NR 对比 `inlet_fraction` vs `hamel_eq27_explicit`（`scripts/audit_energy_heat_loss_mode_compare.py`）。
2. 出口 T validation（待反应/能量校正后对标 Fig.7）。
3. 任何改动后重新运行 `python3 scripts/audit_hamel_consistency.py` 与 `python3 tests/sanity_checks.py`。

**测试锚点（Eq.(2.7)/(2.8)）**：`tests/test_energy_balance_parity.py`、`tests/test_freeboard_heat_loss_parity.py`、`tests/test_energy_balance_H_in.py`

---

**计数**：一致 4 / 部分一致 2 / 不一致 0 / 缺失 0 / 有意偏离 0
