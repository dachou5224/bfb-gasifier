"""Freeboard gas reaction step and background profile.

Ref: Hamel (1999) Ch.3 freeboard particle trajectory; OPT-003 split from freeboard_segment.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.constants import Rg, g
from src.core.energy_balance_defaults import (
    DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    apply_freeboard_substep_gas_enthalpy_heat_loss,
)
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    N_GAS,
    gas_density_ideal,
    gas_diffusivity_correlation,
    gas_viscosity_power_law,
)
from src.kinetics.char_reactions import (
    R4_M_C,
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R6, rate_R7, rate_R8, rate_R12
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from src.physics.freeboard import (
    calc_beta_a,
    calc_cd_hamel_eq384,
    calc_cd_haider,
    calc_terminal_velocity_haider,
    calc_u_gb,
)
from src.core.freeboard.transport import (
    _bound_extent_forward,
    _bound_extent_signed,
    _gas_concentrations,
    _solve_T_from_enthalpy,
)
from src.core.freeboard.types import FreeboardBackground

def _freeboard_reaction_step(
    *,
    N: np.ndarray,
    T: float,
    P: float,
    V_seg: float,
    dt: float,
    fuel_type: str,
    rho_cat: float,
    char_area_total: float = 0.0,
    m_char_total: float = 0.0,
    solid_d_p: float | None = None,
    D_g: float | None = None,
    char_conversion: float = 0.0,
    r4_scale: float = 1.0,
    r8_scale: float = 1.0,
    enabled_reactions: tuple[str, ...] | None = None,
    extent_limiters_enabled: bool = True,
) -> tuple[np.ndarray, dict[str, float]]:
    """轴向自由板子步反应。

    ``extent_limiters_enabled`` 此处指 **轴向 dt–库存界**（显式轨迹数值壳），
    与床层代数 ``extent_limiters_enabled_thesis``（改写 F(x) 动力学）解耦。
    生产默认由 ``freeboard_extent_dt_bounds_thesis=True`` 传入。
    """
    from src.core import cell_kinetics as _ck

    if _ck._EXTENT_LIMITERS_ENABLED_OVERRIDE is not None:
        limiters_on = bool(_ck._EXTENT_LIMITERS_ENABLED_OVERRIDE)
    else:
        limiters_on = bool(extent_limiters_enabled)
    enabled = set(enabled_reactions or ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R10", "R11", "R12"))
    y, c = _gas_concentrations(N, T, P)
    idx = GAS_SPECIES_INDEX
    c_tar = max(float(c[idx["TAR1"]] + c[idx["TAR2"]]), 0.0)

    r5 = rate_R5_bubble(T, c[idx["CO"]], c[idx["O2"]], c[idx["H2O"]], P, y, idx) * V_seg if "R5" in enabled else 0.0
    r6 = rate_R6(T, c[idx["CH4"]], c[idx["O2"]]) * V_seg if "R6" in enabled else 0.0
    r12 = rate_R12(T, c[idx["H2"]], c[idx["O2"]], P, y, idx) * V_seg if "R12" in enabled else 0.0
    r7 = rate_R7(
        T, c[idx["CH4"]], c[idx["H2O"]], c[idx["CO"]], c[idx["H2"]],
        rho_ash_kg_m3=max(float(rho_cat), 0.0),
    ) * V_seg if "R7" in enabled else 0.0
    r8 = float(r8_scale) * rate_R8(
        T,
        P,
        float(y[idx["CO"]]),
        float(y[idx["H2O"]]),
        float(y[idx["CO2"]]),
        float(y[idx["H2"]]),
        rho_ash_kg_m3=max(float(rho_cat), 0.0),
    ) * V_seg if "R8" in enabled else 0.0
    r10 = rate_R10(T, c_tar, c[idx["O2"]], P, fuel_type) * V_seg if ("R10" in enabled and c_tar > 0.0) else 0.0
    r11b = rate_R11_bubble(T, c_tar) * V_seg if ("R11" in enabled and c_tar > 0.0) else 0.0
    r11d = (
        rate_R11_suspension(T, c_tar, rho_cat) * V_seg
        if ("R11" in enabled and c_tar > 0.0 and rho_cat > 0.0)
        else 0.0
    )
    total_area = max(float(char_area_total), 0.0)
    m_char = max(float(m_char_total), 0.0)  # [kg]
    r1 = r2 = r3 = 0.0
    r4_specific = 0.0  # [1/s] Hamel Eq.5.26
    alpha = 1.0
    if total_area > 0.0 and solid_d_p is not None and D_g is not None:
        d_core = d_core_from_spm_char_conversion(float(char_conversion), float(solid_d_p))
        if "R1" in enabled:
            # 自由板：无完整轨迹滑移时取 Eq.5.17 的 Re→0 极限 Sh=2
            r1, alpha = rate_R1(
                T,
                c[idx["O2"]],
                float(solid_d_p),
                float(D_g),
                None,
                fuel_type,
                geometry="spm",
            )
        if "R2" in enabled:
            r2 = rate_R2(T, c[idx["H2O"]], float(solid_d_p), float(D_g), d_core)
        if "R3" in enabled:
            r3 = rate_R3(T, c[idx["H2"]], float(solid_d_p), float(D_g), d_core)
    if m_char > 0.0 and "R4" in enabled and solid_d_p is not None and D_g is not None:
        d_core = d_core_from_spm_char_conversion(float(char_conversion), float(solid_d_p))
        p_co2 = max(float(c[idx["CO2"]] * Rg * T), 0.0)  # [Pa]
        p_co = max(float(c[idx["CO"]] * Rg * T), 0.0)  # [Pa]
        r4_specific = float(r4_scale) * rate_R4_effective(
            T, p_co2, p_co, float(solid_d_p), float(D_g), d_core
        )

    sto10 = get_lumped_tar_stoichiometry("R10", fuel_type)
    sto11 = get_lumped_tar_stoichiometry("R11", fuel_type)

    # R5 Hamel: CO + 1/2 O2 → CO2；extent = mol CO/s
    o2_demand = (
        0.5 * r5
        + 1.5 * r6
        + 0.5 * r12
        + max(0.0, -float(sto10.get("O2", 0.0))) * r10
        + alpha * r1 * total_area
    )
    h2o_demand = max(r7, 0.0) + max(r8, 0.0) + max(0.0, -float(sto11.get("H2O", 0.0))) * (r11b + r11d) + r2 * total_area
    if limiters_on:
        o2_limit = min(1.0, 0.995 * max(float(N[idx["O2"]]), 0.0) / max(o2_demand * dt, 1e-12))
        h2o_limit = min(1.0, 0.995 * max(float(N[idx["H2O"]]), 0.0) / max(h2o_demand * dt, 1e-12))
        r5 *= o2_limit
        r6 *= o2_limit
        r10 *= o2_limit
        r12 *= o2_limit
        r11b *= h2o_limit
        r11d *= h2o_limit
        r7 *= h2o_limit
        r1 *= o2_limit
        r2 *= h2o_limit
        if r8 > 0.0:
            r8 *= h2o_limit

    def _extent_fwd(rate_mol_s: float, reactants: dict[str, float]) -> float:
        if limiters_on:
            return _bound_extent_forward(rate_mol_s, dt, N, reactants)
        return max(float(rate_mol_s) * float(dt), 0.0)

    def _extent_signed(
        rate_mol_s: float,
        reactants_fwd: dict[str, float],
        reactants_rev: dict[str, float],
    ) -> float:
        if limiters_on:
            return _bound_extent_signed(rate_mol_s, dt, N, reactants_fwd, reactants_rev)
        return float(rate_mol_s) * float(dt)

    e1 = _extent_fwd(r1 * total_area, {"O2": max(alpha, 1e-12)}) if r1 > 0.0 and total_area > 0.0 else 0.0
    e2 = _extent_fwd(r2 * total_area, {"H2O": 1.0}) if r2 > 0.0 and total_area > 0.0 else 0.0
    e3 = _extent_fwd(r3 * total_area, {"H2": 2.0}) if r3 > 0.0 and total_area > 0.0 else 0.0
    # Hamel Eq.5.33：R [mol C/s] = R* [1/s] · m / M_C
    e4 = (
        _extent_fwd(r4_specific * m_char / R4_M_C, {"CO2": 1.0})
        if r4_specific > 0.0 and m_char > 0.0
        else 0.0
    )
    e5 = _extent_fwd(r5, {"CO": 1.0, "O2": 0.5})
    e6 = _extent_fwd(r6, {"CH4": 1.0, "O2": 1.5})
    e12 = _extent_fwd(r12, {"H2": 1.0, "O2": 0.5})
    e7 = _extent_signed(r7, {"CH4": 1.0, "H2O": 1.0}, {"CO": 1.0, "H2": 3.0})
    e8 = _extent_signed(r8, {"CO": 1.0, "H2O": 1.0}, {"CO2": 1.0, "H2": 1.0})
    e10 = (
        _extent_fwd(
            r10,
            {
                "TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0,
                "O2": max(0.0, -float(sto10.get("O2", 0.0))),
            },
        )
        if r10 > 0.0
        else 0.0
    )
    e11b = (
        _extent_fwd(
            r11b,
            {"TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0},
        )
        if r11b > 0.0
        else 0.0
    )
    e11d = (
        _extent_fwd(
            r11d,
            {"TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0},
        )
        if r11d > 0.0
        else 0.0
    )
    e11 = e11b + e11d

    out = np.array(N, dtype=np.float64, copy=True)

    out[idx["CO"]] += -e5 + e6 + e7 - e8
    out[idx["O2"]] += -0.5 * e5 - 1.5 * e6 - 0.5 * e12
    out[idx["CO2"]] += e5 + e8
    out[idx["H2O"]] += 2.0 * e6 + e12 - e7 - e8
    out[idx["CH4"]] += -e6 - e7
    out[idx["H2"]] += -e12 + 3.0 * e7 + e8
    if total_area > 0.0:
        out[idx["O2"]] += -alpha * e1
        out[idx["CO"]] += (2.0 * (1.0 - alpha) * e1 + e2 + 2.0 * e4)
        out[idx["CO2"]] += (2.0 * alpha - 1.0) * e1 - e4
        out[idx["H2O"]] += -e2
        out[idx["H2"]] += e2 - 2.0 * e3
        out[idx["CH4"]] += e3

    if c_tar > 0.0:
        tar_total = max(float(out[idx["TAR1"]] + out[idx["TAR2"]]), 0.0)
        if tar_total > 1e-12:
            tar1_frac = float(out[idx["TAR1"]]) / tar_total
        else:
            tar1_frac = 0.5
        out[idx["TAR1"]] = max(float(out[idx["TAR1"]]) - tar1_frac * (e10 + e11), 0.0)
        out[idx["TAR2"]] = max(float(out[idx["TAR2"]]) - (1.0 - tar1_frac) * (e10 + e11), 0.0)
        for sp, nu in sto10.items():
            if sp not in idx or sp in {"TAR1", "TAR2"}:
                continue
            out[idx[sp]] += e10 * float(nu)
        for sp, nu in sto11.items():
            if sp not in idx or sp in {"TAR1", "TAR2"}:
                continue
            out[idx[sp]] += e11 * float(nu)

    np.maximum(out, 0.0, out=out)
    return out, {
        "R1": float(e1),
        "R2": float(e2),
        "R3": float(e3),
        "R4": float(e4),
        "R5": float(e5),
        "R6": float(e6),
        "R7": float(e7),
        "R8": float(e8),
        "R10": float(e10),
        "R11": float(e11),
        "R11b": float(e11b),
        "R11d": float(e11d),
        "R12": float(e12),
    }


def _build_freeboard_background_profile(
    *,
    N_in: np.ndarray,
    T_in: float,
    P: float,
    A: float,
    H_freeboard: float,
    n_cells: int,
    u_gb0: float,
    beta_a: float,
    heat_loss_frac: float,
    fuel_type: str,
    z_base_m: float,
    total_height_m: float | None,
    secondary_injection_xi: float | None,
    secondary_N: np.ndarray,
    secondary_T_K: float,
    secondary_injection_mode: str,
    enabled_reactions: tuple[str, ...] | None,
    rho_cat_profile: np.ndarray | None = None,
    heat_loss_mode: str = DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE,
    Q_wall_W: float = 0.0,
    Q_heat_exchanger_W: float = 0.0,
    extent_limiters_enabled: bool = False,
) -> list[FreeboardBackground]:
    """Build a segment-wise frozen gas background profile for the whole freeboard."""
    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    dh = H_freeboard / max(int(n_cells), 1)
    h_total = None if total_height_m is None else float(max(total_height_m, 1e-12))
    xi_inj = None if secondary_injection_xi is None else float(np.clip(secondary_injection_xi, 0.0, 1.0))
    secondary_injection_applied = False
    u_gb_prev = float(u_gb0)
    backgrounds: list[FreeboardBackground] = []

    for i in range(int(n_cells)):
        n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
        u0 = (n_tot * Rg * T / P) / max(A, 1e-12)
        u_gb_in = float(u_gb_prev if i > 0 else u_gb0)
        exp_arg_out = float(np.clip(-beta_a * dh, -200.0, 200.0))
        exp_arg_center = float(np.clip(-beta_a * 0.5 * dh, -200.0, 200.0))
        u_gb_out = u0 + (u_gb_in - u0) * math.exp(exp_arg_out)
        u_gb_center = u0 + (u_gb_in - u0) * math.exp(exp_arg_center)
        beta_dh = float(beta_a * dh)
        if abs(beta_dh) <= 1e-12:
            u_gb = float(u_gb_in)
        else:
            decay_mean = -math.expm1(-beta_dh) / beta_dh
            u_gb = float(u0 + (u_gb_in - u0) * decay_mean)
        z_lower = float(i * dh)
        z_upper = float((i + 1) * dh)
        z_center = float((i + 0.5) * dh)
        xi_seg = (
            float((float(z_base_m) + z_center) / h_total)
            if h_total is not None
            else float(z_center / max(H_freeboard, 1e-12))
        )
        tau = dh / max(0.5 * (u_gb_in + u_gb_out), 1e-6)
        V_seg = A * dh
        y_map = {sp: float(max(N[GAS_SPECIES_INDEX[sp]], 0.0) / n_tot) for sp in GAS_SPECIES}
        rho_g = gas_density_ideal(P, T, y_map)
        mu_g = gas_viscosity_power_law(T, 1.8e-5)
        backgrounds.append(
            FreeboardBackground(
                segment_index=i,
                z_lower_m=z_lower,
                z_upper_m=z_upper,
                z_center_m=z_center,
                u0=float(u0),
                u_gb_in=float(u_gb_in),
                u_gb_out=float(u_gb_out),
                u_gb=float(u_gb),
                u_gb_center=float(u_gb_center),
                tau=float(tau),
                V_seg=float(V_seg),
                rho_g=float(rho_g),
                mu_g=float(mu_g),
            )
        )

        if (
            not secondary_injection_applied
            and xi_inj is not None
            and float(np.sum(secondary_N)) > 0.0
            and xi_seg >= xi_inj
        ):
            H_mix = (
                calc_gas_enthalpy_flow(N, T, h_cache={})
                + calc_gas_enthalpy_flow(secondary_N, float(secondary_T_K), h_cache={})
            )
            N = N + secondary_N
            T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
            secondary_injection_applied = True

        dt_slice = tau / 12.0
        for _ in range(12):
            H_before = calc_gas_enthalpy_flow(N, T, h_cache={})
            N, _ = _freeboard_reaction_step(
                N=N,
                T=T,
                P=P,
                V_seg=V_seg,
                dt=dt_slice,
                fuel_type=fuel_type,
                rho_cat=float(rho_cat_profile[i]) if rho_cat_profile is not None and i < len(rho_cat_profile) else 0.0,
                D_g=gas_diffusivity_correlation(T, P),
                enabled_reactions=enabled_reactions,
                extent_limiters_enabled=extent_limiters_enabled,
            )
            H_after = apply_freeboard_substep_gas_enthalpy_heat_loss(
                H_before,
                heat_loss_mode=heat_loss_mode,
                heat_loss_frac_total=heat_loss_frac,
                n_cells=n_cells,
                n_substeps=12,
                n_slices=1,
                Q_wall_W=Q_wall_W,
                Q_heat_exchanger_W=Q_heat_exchanger_W,
            )
            T = _solve_T_from_enthalpy(N, H_after, T)

        u_gb_prev = float(u_gb_out)

    return backgrounds


def _estimate_rho_cat_profile_from_projected_transport(
    *,
    class_transport_by_seg: list[dict[str, np.ndarray]],
    backgrounds: list[FreeboardBackground],
) -> np.ndarray:
    rho_cat_profile = np.zeros(len(class_transport_by_seg), dtype=np.float64)
    for i, class_transport in enumerate(class_transport_by_seg):
        if i >= len(backgrounds):
            break
        hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
        hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)
        up_char = np.asarray(class_transport["m_dot_auf_char_classes"], dtype=np.float64)
        up_ash = np.asarray(class_transport["m_dot_auf_ash_classes"], dtype=np.float64)
        m_hold_seg = float(np.sum(hold_char + hold_ash))
        m_dot_char = float(np.sum(np.maximum(up_char, 0.0)))
        m_dot_ash = float(np.sum(np.maximum(up_ash, 0.0)))
        m_dot_s_avg = max(m_dot_char + m_dot_ash, 1e-12)
        ash_hold_frac = m_dot_ash / m_dot_s_avg
        rho_cat_profile[i] = (m_hold_seg * ash_hold_frac) / max(float(backgrounds[i].V_seg), 1e-12)
    return rho_cat_profile

