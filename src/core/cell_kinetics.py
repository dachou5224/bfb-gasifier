"""单 cell 反应源项组装。

Ref: Hamel (1999) Eq. 2.1 - 2.4
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.kinetics.char_reactions import (
    R4_M_C,
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R7, rate_R8, rate_R9, rate_R12
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble, rate_R11_suspension


# 诊断钩子：脚本可临时设为 "r6_first" / "r10_first"；默认 None → 比例限幅。
# 非 Hamel 同构，禁止作生产默认。
_O2_ALLOC_MODE_OVERRIDE: str | None = None

# 诊断钩子：覆盖 ``extent_limiters_enabled`` 参数。
# None → 用函数参数（生产默认 False = Hamel 同构裸动力学）。
# True/False → 强制开/关工程限幅（燃料预封顶 + 比例 O₂/H₂O）。
_EXTENT_LIMITERS_ENABLED_OVERRIDE: bool | None = None


def _sequential_o2_scales(
    budget: float,
    ordered_demands: list[tuple[str, float]],
) -> dict[str, float]:
    """按优先级顺序分配有限 O₂ 预算，返回各反应 extent 缩放因子 ∈ [0,1]。"""
    scales: dict[str, float] = {}
    remaining = max(float(budget), 0.0)
    for name, demand in ordered_demands:
        dem = max(float(demand), 0.0)
        if dem <= 1e-30:
            scales[name] = 1.0
            continue
        take = min(dem, remaining)
        scales[name] = take / dem
        remaining = max(remaining - take, 0.0)
    return scales


@dataclass
class ReactionSourceBundle:
    """反应源项与 limiter 诊断。"""

    R_gas_b: npt.NDArray[np.float64]
    R_gas_d: npt.NDArray[np.float64]
    R_solid: npt.NDArray[np.float64]
    limit_factor_o2: float
    limit_factor_o2_bubble: float
    limit_factor_o2_dense: float
    limit_factor_co_bubble: float
    limit_factor_ch4_bubble: float
    limit_factor_ch4_dense: float
    limit_factor_h2o: float
    co_supply_bubble: float
    ch4_supply_bubble: float
    ch4_supply_dense: float
    o2_supply_bubble: float
    o2_supply_dense: float
    o2_budget_bubble: float
    o2_budget_dense: float
    o2_demand_bubble_raw: float
    o2_demand_dense_raw: float
    o2_demand_bubble_limited: float
    o2_demand_dense_limited: float
    o2_slack_bubble: float
    o2_slack_dense: float
    o2_shortfall_bubble: float
    o2_shortfall_dense: float
    o2_transfer_potential_bd: float
    h2o_supply: float
    h2o_demand_raw: float
    h2o_demand_limited: float
    net_molar_gas_source_total: float
    net_molar_gas_source_vm: float
    net_molar_gas_source_r5: float
    net_molar_gas_source_r6: float
    net_molar_gas_source_r7: float
    net_molar_gas_source_r8: float
    net_molar_gas_source_r9: float
    net_molar_gas_source_r10: float
    net_molar_gas_source_r11: float
    net_molar_gas_source_r12: float
    net_molar_gas_source_char: float
    extent_r1: float
    extent_r2: float
    extent_r3: float
    extent_r4: float
    # 限幅后均相/焦油 extent [mol/s]（审计用；气泡+悬浮已分列）
    extent_r5_bubble: float = 0.0
    extent_r5_dense: float = 0.0
    extent_r6_bubble: float = 0.0
    extent_r6_dense: float = 0.0
    extent_r7: float = 0.0
    extent_r8: float = 0.0
    extent_r9: float = 0.0
    extent_r10_bubble: float = 0.0
    extent_r10_dense: float = 0.0
    extent_r11_bubble: float = 0.0
    extent_r11_dense: float = 0.0
    extent_r12_bubble: float = 0.0
    extent_r12_dense: float = 0.0
    extent_r1_o2: float = 0.0  # R1 消耗的 O₂ [mol/s] = α·extent_r1


def build_reaction_sources(
    *,
    T: float,
    P: float,
    fuel_type: str,
    target_tar_hc_ratio: float | None = None,
    V_b: float,
    V_d: float,
    C_b: npt.ArrayLike,
    C_d: npt.ArrayLike,
    y_b: npt.ArrayLike,
    y_d: npt.ArrayLike,
    gas_src_vm: npt.ArrayLike,
    solid_sink_vm: npt.ArrayLike,
    areas: npt.ArrayLike,
    solid_d_p: float,
    D_g: float,
    char_conversion: float,
    rho_cat: float,
    rho_ash_kg_m3: float,
    enable_r12: bool,
    use_gibbs_minor: bool,
    gibbs_minor_sources: dict[str, float] | None,
    r4_scale: float,
    r2_scale: float,
    r5_scale: float,
    r5_bubble_scale: float = 1.0,
    r6_scale: float,
    r7_scale: float,
    r8_scale: float,
    rate_multiplier: float,
    N_zu_d: npt.ArrayLike,
    N_d_in: npt.ArrayLike,
    N_zu_b: npt.ArrayLike,
    N_b_in: npt.ArrayLike,
    N_rez_d: npt.ArrayLike,
    N_rez_b: npt.ArrayLike,
    N_ex: npt.ArrayLike,
    solid_shape: tuple[int, int],
    char_index: int,
    solid_d_p_classes: npt.ArrayLike | None = None,
    m_char_classes: npt.ArrayLike | None = None,
    u_slip_r1: float | None = None,
    nu_g_r1: float | None = None,
    extent_limiters_enabled: bool = False,
) -> ReactionSourceBundle:
    """组装气相/固相反应源项。

    Ref: Hamel (1999) Eq. 2.1 - 2.4

    ``extent_limiters_enabled``：False（Hamel 默认）时原始动力学直通 F(x)；
    True 时启用燃料预封顶 + 比例 O₂/H₂O/CH₄ 限幅（工程路径）。
    模块 override ``_EXTENT_LIMITERS_ENABLED_OVERRIDE`` 优先于本参数。

    ``u_slip_r1`` [m/s]、``nu_g_r1`` [m²/s]：供 R1 的 Eq.5.17 Sherwood；
    缺省则 Sh→2（Re→0 极限）。
    """
    idx = GAS_SPECIES_INDEX
    C_b_arr = np.asarray(C_b, dtype=np.float64)
    C_d_arr = np.asarray(C_d, dtype=np.float64)
    y_b_arr = np.asarray(y_b, dtype=np.float64)
    y_d_arr = np.asarray(y_d, dtype=np.float64)
    gas_src_vm_arr = np.asarray(gas_src_vm, dtype=np.float64)
    solid_sink_vm_arr = np.asarray(solid_sink_vm, dtype=np.float64)
    areas_arr = np.asarray(areas, dtype=np.float64)
    if solid_d_p_classes is None:
        d_p_classes_arr = np.full(solid_shape[0], float(solid_d_p), dtype=np.float64)
    else:
        d_p_classes_arr = np.asarray(solid_d_p_classes, dtype=np.float64)
    if m_char_classes is None:
        m_char_arr = np.zeros(solid_shape[0], dtype=np.float64)
    else:
        m_char_arr = np.maximum(np.asarray(m_char_classes, dtype=np.float64), 0.0)
    assert C_b_arr.shape == (N_GAS,)
    assert C_d_arr.shape == (N_GAS,)
    assert y_b_arr.shape == (N_GAS,)
    assert y_d_arr.shape == (N_GAS,)
    assert gas_src_vm_arr.shape == (N_GAS,)
    assert solid_sink_vm_arr.shape == solid_shape
    assert areas_arr.shape == (solid_shape[0],)
    assert d_p_classes_arr.shape == (solid_shape[0],)
    assert m_char_arr.shape == (solid_shape[0],)
    n_ex_arr = np.asarray(N_ex, dtype=np.float64)
    assert n_ex_arr.shape == (N_GAS,)

    R_gas_b = np.zeros(N_GAS, dtype=np.float64)
    R_gas_d = gas_src_vm_arr.copy()
    R_solid = solid_sink_vm_arr.copy()

    r5b = float(r5_scale) * float(r5_bubble_scale) * rate_R5_bubble(
        T, C_b_arr[idx["CO"]], C_b_arr[idx["O2"]], C_b_arr[idx["H2O"]], P, y_b_arr, idx
    ) * float(V_b)
    r6b = float(r6_scale) * rate_R6(T, C_b_arr[idx["CH4"]], C_b_arr[idx["O2"]]) * float(V_b)
    r12b = rate_R12(T, C_b_arr[idx["H2"]], C_b_arr[idx["O2"]], P, y_b_arr, idx) * float(V_b) if enable_r12 else 0.0

    r5d = float(r5_scale) * rate_R5_suspension(T, C_d_arr[idx["CO"]], C_d_arr[idx["O2"]], C_d_arr[idx["H2O"]], P, y_d_arr, idx) * float(V_d)
    r6d = float(r6_scale) * rate_R6(T, C_d_arr[idx["CH4"]], C_d_arr[idx["O2"]]) * float(V_d)
    r12d = rate_R12(T, C_d_arr[idx["H2"]], C_d_arr[idx["O2"]], P, y_d_arr, idx) * float(V_d) if enable_r12 else 0.0

    ext7 = float(r7_scale) * rate_R7(
        T,
        C_d_arr[idx["CH4"]],
        C_d_arr[idx["H2O"]],
        C_d_arr[idx["CO"]],
        C_d_arr[idx["H2"]],
        rho_ash_kg_m3=rho_ash_kg_m3,
    ) * float(V_d)
    ext8 = float(r8_scale) * rate_R8(
        T,
        P,
        float(y_d_arr[idx["CO"]]),
        float(y_d_arr[idx["H2O"]]),
        float(y_d_arr[idx["CO2"]]),
        float(y_d_arr[idx["H2"]]),
        rho_ash_kg_m3=rho_ash_kg_m3,
    ) * float(V_d)
    ext9 = 0.0 if use_gibbs_minor else rate_R9(T, C_d_arr[idx["H2S"]], C_d_arr[idx["O2"]]) * float(V_d)

    C_tar_b = max(float(C_b_arr[idx["TAR1"]] + C_b_arr[idx["TAR2"]]), 0.0)
    C_tar_d = max(float(C_d_arr[idx["TAR1"]] + C_d_arr[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(
        T, C_tar_b, C_b_arr[idx["O2"]], P, fuel_type, target_hc_ratio=target_tar_hc_ratio
    ) * float(V_b)
    ext10d = rate_R10(
        T, C_tar_d, C_d_arr[idx["O2"]], P, fuel_type, target_hc_ratio=target_tar_hc_ratio
    ) * float(V_d)
    ext11b = rate_R11_bubble(T, C_tar_b) * float(V_b)
    ext11d = rate_R11_suspension(T, C_tar_d, rho_cat) * float(V_d)

    stoich_r10 = get_lumped_tar_stoichiometry("R10", fuel_type, target_hc_ratio=target_tar_hc_ratio)
    stoich_r11 = get_lumped_tar_stoichiometry("R11", fuel_type, target_hc_ratio=target_tar_hc_ratio)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))
    nu_h2o_r11 = max(0.0, -float(stoich_r11.get("H2O", 0.0)))

    total_area = float(np.sum(areas_arr))
    total_m_char = float(np.sum(m_char_arr))  # [kg]
    r1_by_area = np.zeros_like(areas_arr)
    r2_by_area = np.zeros_like(areas_arr)
    r3_by_area = np.zeros_like(areas_arr)
    # R4：Weeda/Hamel Eq.5.26 给出 R* [1/s]；Eq.5.33 用质量而非面积
    r4_specific = np.zeros_like(areas_arr)  # [1/s]
    alpha_by_area = np.ones_like(areas_arr)
    p_co2 = max(float(C_d_arr[idx["CO2"]] * Rg * T), 0.0)  # [Pa]
    p_co = max(float(C_d_arr[idx["CO"]] * Rg * T), 0.0)  # [Pa]
    if total_area > 0.0:
        for k, (area, d_p_k) in enumerate(zip(areas_arr, d_p_classes_arr)):
            if area <= 0.0 or d_p_k <= 0.0:
                continue
            dc_k = d_core_from_spm_char_conversion(float(char_conversion), float(d_p_k))
            # R1 模拟主链：§5.1.2 SPM + n=1；d_core 仅保留给 R2–R3 签名兼容
            r1_k, alpha_k = rate_R1(
                T,
                C_d_arr[idx["O2"]],
                float(d_p_k),
                D_g,
                None,
                fuel_type,
                u_slip=u_slip_r1,
                nu_g=nu_g_r1,
                geometry="spm",
            )
            r1_by_area[k] = r1_k
            alpha_by_area[k] = alpha_k
            r2_by_area[k] = float(r2_scale) * rate_R2(T, C_d_arr[idx["H2O"]], float(d_p_k), D_g, dc_k)
            r3_by_area[k] = rate_R3(T, C_d_arr[idx["H2"]], float(d_p_k), D_g, dc_k)
    if total_m_char > 0.0:
        for k, (m_k, d_p_k) in enumerate(zip(m_char_arr, d_p_classes_arr)):
            if m_k <= 0.0:
                continue
            dc_k = d_core_from_spm_char_conversion(float(char_conversion), float(d_p_k))
            r4_specific[k] = float(r4_scale) * rate_R4_effective(
                T, p_co2, p_co, float(d_p_k), D_g, dc_k
            )
    r2_extent_raw = float(np.sum(r2_by_area * areas_arr))
    r1_o2_extent_raw = float(np.sum(alpha_by_area * r1_by_area * areas_arr))

    o2_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["O2"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["O2"]])
        + max(-float(n_ex_arr[idx["O2"]]), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["O2"]]),
        1e-6,
    )
    o2_supply_dense = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["O2"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["O2"]])
        + max(float(n_ex_arr[idx["O2"]]), 0.0)
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["O2"]]),
        1e-6,
    )
    # 燃料进料/源项（mol/s）：须在组装 O₂ demand 之前封顶 extent。
    # 否则 raw R10/R6（高 C_tar·C_O2 或高 C_CH4·C_O2）可把 O₂ 需求抬到 10^5–10^6，
    # 使 lf_O2→0，连真实点火氧化也被同因子掐死。自由板侧已对 R10 做 tar bound。
    # Ref: Hamel 计量守恒；freeboard ``_bound_extent_forward`` tar/O2
    tar_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_zu_b, dtype=np.float64)[idx["TAR2"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["TAR2"]])
        + max(-(float(n_ex_arr[idx["TAR1"]]) + float(n_ex_arr[idx["TAR2"]])), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["TAR2"]]),
        1e-6,
    )
    tar_supply_dense = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_zu_d, dtype=np.float64)[idx["TAR2"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["TAR2"]])
        + max(float(n_ex_arr[idx["TAR1"]]) + float(n_ex_arr[idx["TAR2"]]), 0.0)
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["TAR1"]])
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["TAR2"]])
        + float(gas_src_vm_arr[idx["TAR1"]])
        + float(gas_src_vm_arr[idx["TAR2"]]),
        1e-6,
    )
    co_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["CO"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["CO"]])
        + max(-float(n_ex_arr[idx["CO"]]), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["CO"]]),
        1e-6,
    )
    co_supply_dense = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["CO"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["CO"]])
        + max(float(n_ex_arr[idx["CO"]]), 0.0)
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["CO"]])
        + float(gas_src_vm_arr[idx["CO"]]),
        1e-6,
    )
    ch4_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["CH4"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["CH4"]])
        + max(-float(n_ex_arr[idx["CH4"]]), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["CH4"]]),
        1e-6,
    )
    ch4_supply_dense = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["CH4"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["CH4"]])
        + max(float(n_ex_arr[idx["CH4"]]), 0.0)
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["CH4"]])
        + float(gas_src_vm_arr[idx["CH4"]])
        + float(np.sum(r3_by_area * areas_arr)),
        1e-6,
    )
    h2_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["H2"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["H2"]])
        + max(-float(n_ex_arr[idx["H2"]]), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["H2"]]),
        1e-6,
    )
    h2_supply_dense = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["H2"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["H2"]])
        + max(float(n_ex_arr[idx["H2"]]), 0.0)
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["H2"]])
        + float(gas_src_vm_arr[idx["H2"]]),
        1e-6,
    )

    # 优先级：模块 override > 函数参数（生产默认 False = Hamel 同构）
    if _EXTENT_LIMITERS_ENABLED_OVERRIDE is not None:
        limiters_on = bool(_EXTENT_LIMITERS_ENABLED_OVERRIDE)
    else:
        limiters_on = bool(extent_limiters_enabled)

    # 工程路径：燃料计量预封顶（写入 extent，再进 O₂ demand）。Hamel 同构跳过。
    if limiters_on:
        ext10b = min(float(ext10b), 0.995 * tar_supply_bubble)
        ext10d = min(float(ext10d), 0.995 * tar_supply_dense)
        r6b = min(float(r6b), 0.995 * ch4_supply_bubble)
        r6d = min(float(r6d), 0.995 * ch4_supply_dense)
        # R5 Hamel Eq.5.36: CO + 1/2 O2 → CO2；extent = mol CO/s
        r5b = min(float(r5b), 0.995 * co_supply_bubble)
        r5d = min(float(r5d), 0.995 * co_supply_dense)
        if enable_r12:
            r12b = min(float(r12b), 0.995 * h2_supply_bubble)
            r12d = min(float(r12d), 0.995 * h2_supply_dense)

    bubble_o2_demand = (
        0.5 * r5b
        + 1.5 * r6b
        + 0.5 * r12b
        + nu_o2_r10 * ext10b
    )
    dense_o2_demand = (
        0.5 * r5d
        + 1.5 * r6d
        + 0.5 * r12d
        + r1_o2_extent_raw
        + 1.5 * ext9
        + nu_o2_r10 * ext10d
    )
    o2_budget_bubble = 0.995 * o2_supply_bubble
    o2_budget_dense = 0.995 * o2_supply_dense
    # CO supply limiter for the bubble phase（工程路径护栏）。
    bubble_co_demand_r5 = r5b
    h2o_supply = max(
        float(np.asarray(N_zu_d, dtype=np.float64)[idx["H2O"]])
        + float(np.asarray(N_d_in, dtype=np.float64)[idx["H2O"]])
        + float(np.asarray(N_zu_b, dtype=np.float64)[idx["H2O"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["H2O"]])
        + float(np.asarray(N_rez_d, dtype=np.float64)[idx["H2O"]])
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["H2O"]])
        + float(gas_src_vm_arr[idx["H2O"]]),
        1e-6,
    )
    h2o_consume = r2_extent_raw + ext7 + max(0.0, ext8) + nu_h2o_r11 * (ext11b + ext11d)

    if not limiters_on:
        # Hamel 同构：原始动力学直通（无燃料预封顶 / 比例 O₂/H₂O/CH₄ 限幅）
        limit_factor_o2_bubble = 1.0
        limit_factor_o2_dense = 1.0
        limit_factor_o2 = 1.0
        limit_factor_co_bubble = 1.0
        limit_factor_h2o = 1.0
        limit_factor_ch4_bubble = 1.0
        limit_factor_ch4_dense = 1.0
        r5b_lim = r5b
        r6b_lim = r6b
        r12b_lim = r12b
        ext10b_lim = ext10b
        r5d_lim = r5d
        r6d_lim = r6d
        r12d_lim = r12d
        r1_by_area_lim = r1_by_area.copy()
        ext9_lim = ext9
        ext10d_lim = ext10d
        r2_by_area_lim = r2_by_area.copy()
        ext7_lim = ext7
        ext11b_lim = ext11b
        ext11d_lim = ext11d
        ext8_lim = ext8
        r3_by_area_lim = r3_by_area.copy()
        r4_specific_lim = r4_specific.copy()
    else:
        limit_factor_o2_bubble = min(1.0, o2_budget_bubble / max(bubble_o2_demand, 1e-9))
        limit_factor_o2_dense = min(1.0, o2_budget_dense / max(dense_o2_demand, 1e-9))
        limit_factor_o2 = min(limit_factor_o2_bubble, limit_factor_o2_dense)
        limit_factor_co_bubble = min(1.0, 0.995 * co_supply_bubble / max(bubble_co_demand_r5, 1e-9))

        o2_mode = str(_O2_ALLOC_MODE_OVERRIDE or "proportional").strip().lower()
        if o2_mode in ("r6_first", "r10_first"):
            # 诊断：顺序分配有限 O₂（非 Hamel）。小需求先满足，再按模式优先 R6 或 R10。
            dem_b_base = [
                ("R5", 0.5 * float(r5b)),
                ("R12", 0.5 * float(r12b)),
                ("R6", 1.5 * float(r6b)),
                ("R10", float(nu_o2_r10) * float(ext10b)),
            ]
            dem_d_base = [
                ("R5", 0.5 * float(r5d)),
                ("R12", 0.5 * float(r12d)),
                ("R1", float(r1_o2_extent_raw)),
                ("R9", 1.5 * float(ext9)),
                ("R6", 1.5 * float(r6d)),
                ("R10", float(nu_o2_r10) * float(ext10d)),
            ]
            if o2_mode == "r6_first":
                order_b = ["R5", "R12", "R6", "R10"]
                order_d = ["R5", "R12", "R1", "R9", "R6", "R10"]
            else:
                order_b = ["R10", "R5", "R12", "R6"]
                order_d = ["R10", "R5", "R12", "R1", "R9", "R6"]
            dem_b_map = dict(dem_b_base)
            dem_d_map = dict(dem_d_base)
            sc_b = _sequential_o2_scales(
                o2_budget_bubble, [(n, dem_b_map[n]) for n in order_b]
            )
            sc_d = _sequential_o2_scales(
                o2_budget_dense, [(n, dem_d_map[n]) for n in order_d]
            )
            r5b_lim = r5b * min(sc_b["R5"], limit_factor_co_bubble)
            r6b_lim = r6b * sc_b["R6"]
            r12b_lim = r12b * sc_b["R12"]
            ext10b_lim = ext10b * sc_b["R10"]
            r5d_lim = r5d * sc_d["R5"]
            r6d_lim = r6d * sc_d["R6"]
            r12d_lim = r12d * sc_d["R12"]
            r1_by_area_lim = r1_by_area * sc_d["R1"]
            ext9_lim = ext9 * sc_d["R9"]
            ext10d_lim = ext10d * sc_d["R10"]
            # 诊断用：把“有效”限幅记为预算/总需求（与比例限幅同口径）
            limit_factor_o2_bubble = min(1.0, o2_budget_bubble / max(bubble_o2_demand, 1e-9))
            limit_factor_o2_dense = min(1.0, o2_budget_dense / max(dense_o2_demand, 1e-9))
            limit_factor_o2 = min(limit_factor_o2_bubble, limit_factor_o2_dense)
        else:
            r5b_lim = r5b * min(limit_factor_o2_bubble, limit_factor_co_bubble)
            r6b_lim = r6b * limit_factor_o2_bubble
            r12b_lim = r12b * limit_factor_o2_bubble
            ext10b_lim = ext10b * limit_factor_o2_bubble

            r5d_lim = r5d * limit_factor_o2_dense
            r6d_lim = r6d * limit_factor_o2_dense
            r12d_lim = r12d * limit_factor_o2_dense
            r1_by_area_lim = r1_by_area * limit_factor_o2_dense
            ext9_lim = ext9 * limit_factor_o2_dense
            ext10d_lim = ext10d * limit_factor_o2_dense

        limit_factor_h2o = min(1.0, (0.995 * h2o_supply) / max(h2o_consume, 1e-9))

        r2_by_area_lim = r2_by_area * limit_factor_h2o
        ext7_lim = ext7 * limit_factor_h2o
        ext11b_lim = ext11b * limit_factor_h2o
        ext11d_lim = ext11d * limit_factor_h2o
        ext8_lim = ext8 * limit_factor_h2o if ext8 > 0.0 else ext8
        r3_by_area_lim = r3_by_area.copy()
        r4_specific_lim = r4_specific.copy()

        # CH4 supply after H2O-limited R3 contribution to dense CH4 inventory
        ch4_supply_dense = max(
            float(np.asarray(N_zu_d, dtype=np.float64)[idx["CH4"]])
            + float(np.asarray(N_d_in, dtype=np.float64)[idx["CH4"]])
            + max(float(n_ex_arr[idx["CH4"]]), 0.0)
            + float(np.asarray(N_rez_d, dtype=np.float64)[idx["CH4"]])
            + float(gas_src_vm_arr[idx["CH4"]])
            + float(np.sum(r3_by_area_lim * areas_arr)),
            1e-6,
        )
        limit_factor_ch4_bubble = min(1.0, 0.995 * ch4_supply_bubble / max(r6b_lim, 1e-9))
        dense_ch4_demand = r6d_lim + ext7_lim
        limit_factor_ch4_dense = min(1.0, 0.995 * ch4_supply_dense / max(dense_ch4_demand, 1e-9))
        r6b_lim *= limit_factor_ch4_bubble
        r6d_lim *= limit_factor_ch4_dense
        ext7_lim *= limit_factor_ch4_dense

    if rate_multiplier != 1.0:
        rm = float(np.clip(rate_multiplier, 0.0, 1.0))
        r5b_lim *= rm
        r6b_lim *= rm
        r12b_lim *= rm
        r5d_lim *= rm
        r6d_lim *= rm
        r12d_lim *= rm
        r1_by_area_lim *= rm
        r2_by_area_lim *= rm
        r3_by_area_lim *= rm
        r4_specific_lim *= rm
        ext7_lim *= rm
        ext8_lim *= rm
        ext9_lim *= rm
        ext10b_lim *= rm
        ext10d_lim *= rm
        ext11b_lim *= rm
        ext11d_lim *= rm

    # R5: CO + 1/2 O2 → CO2（Hamel §5.2.1）；R6 产 CO 冲抵部分 CO 消耗
    R_gas_b[idx["CO"]] -= r5b_lim - r6b_lim
    R_gas_b[idx["O2"]] -= 0.5 * r5b_lim + 1.5 * r6b_lim + 0.5 * r12b_lim
    R_gas_b[idx["CO2"]] += r5b_lim
    R_gas_b[idx["H2O"]] += 2.0 * r6b_lim + r12b_lim
    R_gas_b[idx["CH4"]] -= r6b_lim
    R_gas_b[idx["H2"]] -= r12b_lim

    R_gas_d[idx["CO"]] -= r5d_lim - r6d_lim
    R_gas_d[idx["O2"]] -= 0.5 * r5d_lim + 1.5 * r6d_lim + 0.5 * r12d_lim
    R_gas_d[idx["CO2"]] += r5d_lim
    R_gas_d[idx["H2O"]] += 2.0 * r6d_lim + r12d_lim
    R_gas_d[idx["CH4"]] -= r6d_lim
    R_gas_d[idx["H2"]] -= r12d_lim

    R_gas_d[idx["CH4"]] -= ext7_lim
    R_gas_d[idx["H2O"]] -= ext7_lim
    R_gas_d[idx["CO"]] += ext7_lim
    R_gas_d[idx["H2"]] += 3.0 * ext7_lim

    R_gas_d[idx["CO"]] -= ext8_lim
    R_gas_d[idx["H2O"]] -= ext8_lim
    R_gas_d[idx["CO2"]] += ext8_lim
    R_gas_d[idx["H2"]] += ext8_lim

    if ext9_lim > 0.0:
        R_gas_d[idx["H2S"]] -= ext9_lim
        R_gas_d[idx["O2"]] -= 1.5 * ext9_lim
        R_gas_d[idx["H2O"]] += ext9_lim

    for sp, nu in stoich_r10.items():
        if sp in idx:
            R_gas_b[idx[sp]] += ext10b_lim * float(nu)
            R_gas_d[idx[sp]] += ext10d_lim * float(nu)
    for sp, nu in stoich_r11.items():
        if sp in idx:
            R_gas_b[idx[sp]] += ext11b_lim * float(nu)
            R_gas_d[idx[sp]] += ext11d_lim * float(nu)

    r1_extent = r2_extent = r3_extent = r4_extent = r1_o2_extent = 0.0
    if total_area > 0.0:
        r1_extent = float(np.sum(r1_by_area_lim * areas_arr))
        r2_extent = float(np.sum(r2_by_area_lim * areas_arr))
        r3_extent = float(np.sum(r3_by_area_lim * areas_arr))
        r1_o2_extent = float(np.sum(alpha_by_area * r1_by_area_lim * areas_arr))
        R_gas_d[idx["O2"]] -= r1_o2_extent
        R_gas_d[idx["CO"]] += 2.0 * (r1_extent - r1_o2_extent) + r2_extent
        R_gas_d[idx["CO2"]] += 2.0 * r1_o2_extent - r1_extent
        R_gas_d[idx["H2O"]] -= r2_extent
        R_gas_d[idx["H2"]] += r2_extent
        R_gas_d[idx["H2"]] -= 2.0 * r3_extent
        R_gas_d[idx["CH4"]] += r3_extent
        # R1–R3：表面通量 [mol/(m²·s)] × 面积 → 固相消耗 [kg/s]
        R_solid[:, char_index] -= (
            (r1_by_area_lim + r2_by_area_lim + r3_by_area_lim)
            * ATOMIC_MASS_KG_PER_MOL["C"]
            * areas_arr
        )
    if total_m_char > 0.0:
        # Hamel Eq.5.33：R [mol C/s] = R* [1/s] · m / M_C；固相消耗 = R* · m [kg/s]
        r4_extent = float(np.sum(r4_specific_lim * m_char_arr / R4_M_C))
        R_gas_d[idx["CO2"]] -= r4_extent
        R_gas_d[idx["CO"]] += 2.0 * r4_extent
        R_solid[:, char_index] -= r4_specific_lim * m_char_arr

    if use_gibbs_minor and gibbs_minor_sources:
        for sp, ddot in gibbs_minor_sources.items():
            if sp in idx:
                R_gas_d[idx[sp]] += float(ddot)

    if np.any(np.isnan(R_gas_b)) or np.any(np.isnan(R_gas_d)) or np.any(np.isnan(R_solid)):
        R_gas_b = np.nan_to_num(R_gas_b)
        R_gas_d = np.nan_to_num(R_gas_d)
        R_solid = np.nan_to_num(R_solid)

    o2_demand_bubble_limited = (
        0.5 * r5b_lim + 1.5 * r6b_lim + 0.5 * r12b_lim + nu_o2_r10 * ext10b_lim
    )
    o2_demand_dense_limited = (
        0.5 * r5d_lim + 1.5 * r6d_lim + 0.5 * r12d_lim + r1_o2_extent
        + 1.5 * ext9_lim + nu_o2_r10 * ext10d_lim
    )
    o2_slack_bubble = max(o2_budget_bubble - o2_demand_bubble_limited, 0.0)
    o2_slack_dense = max(o2_budget_dense - o2_demand_dense_limited, 0.0)
    o2_shortfall_bubble = max(bubble_o2_demand - o2_budget_bubble, 0.0)
    o2_shortfall_dense = max(dense_o2_demand - o2_budget_dense, 0.0)
    o2_transfer_potential_bd = min(o2_slack_bubble, o2_shortfall_dense)
    net_molar_gas_source_vm = float(np.sum(gas_src_vm_arr))
    # CO + 1/2 O2 → CO2 ⇒ Δn = −0.5 mol gas / mol CO
    net_molar_gas_source_r5 = float(-0.5 * (r5b_lim + r5d_lim))
    net_molar_gas_source_r6 = float(0.5 * (r6b_lim + r6d_lim))
    net_molar_gas_source_r7 = float(2.0 * ext7_lim)
    net_molar_gas_source_r8 = 0.0
    net_molar_gas_source_r9 = float(-1.5 * ext9_lim)
    net_molar_gas_source_r10 = float(
        sum(float(nu) for sp, nu in stoich_r10.items() if sp in idx) * (ext10b_lim + ext10d_lim)
    )
    net_molar_gas_source_r11 = float(
        sum(float(nu) for sp, nu in stoich_r11.items() if sp in idx) * (ext11b_lim + ext11d_lim)
    )
    net_molar_gas_source_r12 = float(-0.5 * (r12b_lim + r12d_lim))
    net_molar_gas_source_char = float(
        (r1_extent - r1_o2_extent) + r2_extent - r3_extent + r4_extent
    )
    net_molar_gas_source_total = float(np.sum(R_gas_b) + np.sum(R_gas_d))

    return ReactionSourceBundle(
        R_gas_b=R_gas_b,
        R_gas_d=R_gas_d,
        R_solid=R_solid,
        limit_factor_o2=float(limit_factor_o2),
        limit_factor_o2_bubble=float(limit_factor_o2_bubble),
        limit_factor_o2_dense=float(limit_factor_o2_dense),
        limit_factor_co_bubble=float(limit_factor_co_bubble),
        limit_factor_ch4_bubble=float(limit_factor_ch4_bubble),
        limit_factor_ch4_dense=float(limit_factor_ch4_dense),
        limit_factor_h2o=float(limit_factor_h2o),
        co_supply_bubble=float(co_supply_bubble),
        ch4_supply_bubble=float(ch4_supply_bubble),
        ch4_supply_dense=float(ch4_supply_dense),
        o2_supply_bubble=float(o2_supply_bubble),
        o2_supply_dense=float(o2_supply_dense),
        o2_budget_bubble=float(o2_budget_bubble),
        o2_budget_dense=float(o2_budget_dense),
        o2_demand_bubble_raw=float(bubble_o2_demand),
        o2_demand_dense_raw=float(dense_o2_demand),
        o2_demand_bubble_limited=float(o2_demand_bubble_limited),
        o2_demand_dense_limited=float(o2_demand_dense_limited),
        o2_slack_bubble=float(o2_slack_bubble),
        o2_slack_dense=float(o2_slack_dense),
        o2_shortfall_bubble=float(o2_shortfall_bubble),
        o2_shortfall_dense=float(o2_shortfall_dense),
        o2_transfer_potential_bd=float(o2_transfer_potential_bd),
        h2o_supply=float(h2o_supply),
        h2o_demand_raw=float(h2o_consume),
        h2o_demand_limited=float(
            r2_extent + ext7_lim + max(0.0, ext8_lim) + nu_h2o_r11 * (ext11b_lim + ext11d_lim)
        ),
        net_molar_gas_source_total=net_molar_gas_source_total,
        net_molar_gas_source_vm=net_molar_gas_source_vm,
        net_molar_gas_source_r5=net_molar_gas_source_r5,
        net_molar_gas_source_r6=net_molar_gas_source_r6,
        net_molar_gas_source_r7=net_molar_gas_source_r7,
        net_molar_gas_source_r8=net_molar_gas_source_r8,
        net_molar_gas_source_r9=net_molar_gas_source_r9,
        net_molar_gas_source_r10=net_molar_gas_source_r10,
        net_molar_gas_source_r11=net_molar_gas_source_r11,
        net_molar_gas_source_r12=net_molar_gas_source_r12,
        net_molar_gas_source_char=net_molar_gas_source_char,
        extent_r1=float(r1_extent),
        extent_r2=float(r2_extent),
        extent_r3=float(r3_extent),
        extent_r4=float(r4_extent),
        extent_r5_bubble=float(r5b_lim),
        extent_r5_dense=float(r5d_lim),
        extent_r6_bubble=float(r6b_lim),
        extent_r6_dense=float(r6d_lim),
        extent_r7=float(ext7_lim),
        extent_r8=float(ext8_lim),
        extent_r9=float(ext9_lim),
        extent_r10_bubble=float(ext10b_lim),
        extent_r10_dense=float(ext10d_lim),
        extent_r11_bubble=float(ext11b_lim),
        extent_r11_dense=float(ext11d_lim),
        extent_r12_bubble=float(r12b_lim),
        extent_r12_dense=float(r12d_lim),
        extent_r1_o2=float(r1_o2_extent),
    )
