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
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R7, rate_R8, rate_R9, rate_R12
from src.kinetics.tar_reactions import get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble, rate_R11_suspension


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


def build_reaction_sources(
    *,
    T: float,
    P: float,
    fuel_type: str,
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
    enable_r12: bool,
    use_gibbs_minor: bool,
    gibbs_minor_sources: dict[str, float] | None,
    r4_scale: float,
    r5_scale: float,
    r6_scale: float,
    r7_scale: float,
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
) -> ReactionSourceBundle:
    """组装气相/固相反应源项。

    Ref: Hamel (1999) Eq. 2.1 - 2.4
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
    assert C_b_arr.shape == (N_GAS,)
    assert C_d_arr.shape == (N_GAS,)
    assert y_b_arr.shape == (N_GAS,)
    assert y_d_arr.shape == (N_GAS,)
    assert gas_src_vm_arr.shape == (N_GAS,)
    assert solid_sink_vm_arr.shape == solid_shape
    assert areas_arr.shape == (solid_shape[0],)
    assert d_p_classes_arr.shape == (solid_shape[0],)
    n_ex_arr = np.asarray(N_ex, dtype=np.float64)
    assert n_ex_arr.shape == (N_GAS,)

    R_gas_b = np.zeros(N_GAS, dtype=np.float64)
    R_gas_d = gas_src_vm_arr.copy()
    R_solid = solid_sink_vm_arr.copy()

    r5b = float(r5_scale) * rate_R5_bubble(T, C_b_arr[idx["CO"]], C_b_arr[idx["O2"]], P, y_b_arr, idx) * float(V_b)
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
    ) * float(V_d)
    ext8 = rate_R8(
        T,
        P,
        float(y_d_arr[idx["CO"]]),
        float(y_d_arr[idx["H2O"]]),
        float(y_d_arr[idx["CO2"]]),
        float(y_d_arr[idx["H2"]]),
    ) * float(V_d)
    ext9 = 0.0 if use_gibbs_minor else rate_R9(T, C_d_arr[idx["H2S"]], C_d_arr[idx["O2"]]) * float(V_d)

    C_tar_b = max(float(C_b_arr[idx["TAR1"]] + C_b_arr[idx["TAR2"]]), 0.0)
    C_tar_d = max(float(C_d_arr[idx["TAR1"]] + C_d_arr[idx["TAR2"]]), 0.0)
    ext10b = rate_R10(T, C_tar_b, C_b_arr[idx["O2"]], P, fuel_type) * float(V_b)
    ext10d = rate_R10(T, C_tar_d, C_d_arr[idx["O2"]], P, fuel_type) * float(V_d)
    ext11b = rate_R11_bubble(T, C_tar_b) * float(V_b)
    ext11d = rate_R11_suspension(T, C_tar_d, rho_cat) * float(V_d)

    stoich_r10 = get_lumped_tar_stoichiometry("R10", fuel_type)
    stoich_r11 = get_lumped_tar_stoichiometry("R11", fuel_type)
    nu_o2_r10 = max(0.0, -float(stoich_r10.get("O2", 0.0)))
    nu_h2o_r11 = max(0.0, -float(stoich_r11.get("H2O", 0.0)))

    total_area = float(np.sum(areas_arr))
    r1_by_area = np.zeros_like(areas_arr)
    r2_by_area = np.zeros_like(areas_arr)
    r3_by_area = np.zeros_like(areas_arr)
    r4_by_area = np.zeros_like(areas_arr)
    alpha_by_area = np.ones_like(areas_arr)
    if total_area > 0.0:
        p_co2 = max(float(C_d_arr[idx["CO2"]] * Rg * T), 0.0)
        p_co = max(float(C_d_arr[idx["CO"]] * Rg * T), 0.0)
        for k, (area, d_p_k) in enumerate(zip(areas_arr, d_p_classes_arr)):
            if area <= 0.0 or d_p_k <= 0.0:
                continue
            dc_k = d_core_from_spm_char_conversion(float(char_conversion), float(d_p_k))
            r1_k, alpha_k = rate_R1(T, C_d_arr[idx["O2"]], float(d_p_k), D_g, dc_k, fuel_type)
            r1_by_area[k] = r1_k
            alpha_by_area[k] = alpha_k
            r2_by_area[k] = rate_R2(T, C_d_arr[idx["H2O"]], float(d_p_k), D_g, dc_k)
            r3_by_area[k] = rate_R3(T, C_d_arr[idx["H2"]], float(d_p_k), D_g, dc_k)
            r4_by_area[k] = float(r4_scale) * rate_R4_effective(T, p_co2, p_co, float(d_p_k), D_g, dc_k)
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
    bubble_o2_demand = (
        r5b
        + 1.5 * r6b
        + 0.5 * r12b
        + nu_o2_r10 * ext10b
    )
    dense_o2_demand = (
        r5d
        + 1.5 * r6d
        + 0.5 * r12d
        + r1_o2_extent_raw
        + 1.5 * ext9
        + nu_o2_r10 * ext10d
    )
    o2_budget_bubble = 0.995 * o2_supply_bubble
    o2_budget_dense = 0.995 * o2_supply_dense
    limit_factor_o2_bubble = min(1.0, o2_budget_bubble / max(bubble_o2_demand, 1e-9))
    limit_factor_o2_dense = min(1.0, o2_budget_dense / max(dense_o2_demand, 1e-9))
    limit_factor_o2 = min(limit_factor_o2_bubble, limit_factor_o2_dense)

    # CO supply limiter for the bubble phase.
    # R5 stoichiometry is 2 CO + O2 → 2 CO2 (coefficient 2), so when the system is
    # O2-limited the CO demand can far exceed the CO actually available to the bubble,
    # producing an irreducible residual that prevents NR convergence.
    # Cap r5b by CO inflow to the bubble (same accounting structure as the O2 limiter).
    co_supply_bubble = max(
        float(np.asarray(N_zu_b, dtype=np.float64)[idx["CO"]])
        + float(np.asarray(N_b_in, dtype=np.float64)[idx["CO"]])
        + max(-float(n_ex_arr[idx["CO"]]), 0.0)
        + float(np.asarray(N_rez_b, dtype=np.float64)[idx["CO"]]),
        1e-6,
    )
    bubble_co_demand_r5 = 2.0 * r5b
    limit_factor_co_bubble = min(1.0, 0.995 * co_supply_bubble / max(bubble_co_demand_r5, 1e-9))

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
    limit_factor_h2o = min(1.0, (0.995 * h2o_supply) / max(h2o_consume, 1e-9))

    r2_by_area_lim = r2_by_area * limit_factor_h2o
    ext7_lim = ext7 * limit_factor_h2o
    ext11b_lim = ext11b * limit_factor_h2o
    ext11d_lim = ext11d * limit_factor_h2o
    ext8_lim = ext8 * limit_factor_h2o if ext8 > 0.0 else ext8
    r3_by_area_lim = r3_by_area.copy()
    r4_by_area_lim = r4_by_area.copy()

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
        r4_by_area_lim *= rm
        ext7_lim *= rm
        ext8_lim *= rm
        ext9_lim *= rm
        ext10b_lim *= rm
        ext10d_lim *= rm
        ext11b_lim *= rm
        ext11d_lim *= rm

    R_gas_b[idx["CO"]] -= 2.0 * r5b_lim - r6b_lim
    R_gas_b[idx["O2"]] -= r5b_lim + 1.5 * r6b_lim + 0.5 * r12b_lim
    R_gas_b[idx["CO2"]] += 2.0 * r5b_lim
    R_gas_b[idx["H2O"]] += 2.0 * r6b_lim + r12b_lim
    R_gas_b[idx["CH4"]] -= r6b_lim
    R_gas_b[idx["H2"]] -= r12b_lim

    R_gas_d[idx["CO"]] -= 2.0 * r5d_lim - r6d_lim
    R_gas_d[idx["O2"]] -= r5d_lim + 1.5 * r6d_lim + 0.5 * r12d_lim
    R_gas_d[idx["CO2"]] += 2.0 * r5d_lim
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

    if total_area > 0.0:
        r1_extent = float(np.sum(r1_by_area_lim * areas_arr))
        r2_extent = float(np.sum(r2_by_area_lim * areas_arr))
        r3_extent = float(np.sum(r3_by_area_lim * areas_arr))
        r4_extent = float(np.sum(r4_by_area_lim * areas_arr))
        r1_o2_extent = float(np.sum(alpha_by_area * r1_by_area_lim * areas_arr))
        R_gas_d[idx["O2"]] -= r1_o2_extent
        R_gas_d[idx["CO"]] += 2.0 * (r1_extent - r1_o2_extent) + r2_extent
        R_gas_d[idx["CO2"]] += 2.0 * r1_o2_extent - r1_extent
        R_gas_d[idx["H2O"]] -= r2_extent
        R_gas_d[idx["H2"]] += r2_extent
        R_gas_d[idx["H2"]] -= 2.0 * r3_extent
        R_gas_d[idx["CH4"]] += r3_extent
        R_gas_d[idx["CO2"]] -= r4_extent
        R_gas_d[idx["CO"]] += 2.0 * r4_extent
        R_solid[:, char_index] -= (
            (r1_by_area_lim + r2_by_area_lim + r3_by_area_lim + r4_by_area_lim)
            * ATOMIC_MASS_KG_PER_MOL["C"]
            * areas_arr
        )
    else:
        r1_extent = r2_extent = r3_extent = r4_extent = r1_o2_extent = 0.0

    if use_gibbs_minor and gibbs_minor_sources:
        for sp, ddot in gibbs_minor_sources.items():
            if sp in idx:
                R_gas_d[idx[sp]] += float(ddot)

    if np.any(np.isnan(R_gas_b)) or np.any(np.isnan(R_gas_d)) or np.any(np.isnan(R_solid)):
        R_gas_b = np.nan_to_num(R_gas_b)
        R_gas_d = np.nan_to_num(R_gas_d)
        R_solid = np.nan_to_num(R_solid)

    o2_demand_bubble_limited = (
        r5b_lim + 1.5 * r6b_lim + 0.5 * r12b_lim + nu_o2_r10 * ext10b_lim
    )
    o2_demand_dense_limited = (
        r5d_lim + 1.5 * r6d_lim + 0.5 * r12d_lim + r1_o2_extent
        + 1.5 * ext9_lim + nu_o2_r10 * ext10d_lim
    )
    o2_slack_bubble = max(o2_budget_bubble - o2_demand_bubble_limited, 0.0)
    o2_slack_dense = max(o2_budget_dense - o2_demand_dense_limited, 0.0)
    o2_shortfall_bubble = max(bubble_o2_demand - o2_budget_bubble, 0.0)
    o2_shortfall_dense = max(dense_o2_demand - o2_budget_dense, 0.0)
    o2_transfer_potential_bd = min(o2_slack_bubble, o2_shortfall_dense)
    net_molar_gas_source_vm = float(np.sum(gas_src_vm_arr))
    net_molar_gas_source_r5 = float(-(r5b_lim + r5d_lim))
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
    )
