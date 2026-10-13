"""Gas profiles, size-class transport aggregation.

Ref: Hamel (1999) Ch.3 freeboard particle trajectory; OPT-003 split from freeboard_segment.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.constants import Rg, g
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    N_GAS,
    gas_density_ideal,
    gas_diffusivity_correlation,
    gas_viscosity_power_law,
)
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
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

def _gas_concentrations(N: np.ndarray, T: float, P: float) -> tuple[np.ndarray, np.ndarray]:
    n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
    y = np.maximum(N, 0.0) / n_tot
    c = y * (P / max(Rg * T, 1e-12))
    return y, c


def _wet_gas_profile(N: np.ndarray) -> dict[str, float]:
    n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
    y = np.maximum(N, 0.0) / n_tot
    return {sp: float(y[j]) for j, sp in enumerate(GAS_SPECIES)}


def _char_area_and_representative_diameter(
    hold_char_classes: np.ndarray,
    d_p_classes: np.ndarray,
    rho_p: float,
) -> tuple[float, float | None]:
    hold_char = np.asarray(hold_char_classes, dtype=np.float64)
    d_p = np.asarray(d_p_classes, dtype=np.float64)
    if hold_char.size == 0 or d_p.size == 0:
        return 0.0, None
    areas = np.maximum(hold_char, 0.0) * 6.0 / (max(float(rho_p), 1e-12) * np.maximum(d_p, 1e-12))
    total_area = float(np.sum(areas))
    if total_area <= 0.0:
        return 0.0, None
    weights = areas / total_area
    d_p_eff = float(np.sum(weights * d_p))
    return total_area, d_p_eff


def _bound_extent_forward(rate_mol_s: float, dt: float, N: np.ndarray, reactants: dict[str, float]) -> float:
    extent = max(rate_mol_s * dt, 0.0)
    if extent <= 0.0:
        return 0.0
    max_extent = extent
    for sp, nu in reactants.items():
        j = GAS_SPECIES_INDEX[sp]
        max_extent = min(max_extent, max(float(N[j]), 0.0) / max(float(nu), 1e-12))
    return max(max_extent, 0.0)


def _bound_extent_signed(
    rate_mol_s: float,
    dt: float,
    N: np.ndarray,
    reactants_fwd: dict[str, float],
    reactants_rev: dict[str, float],
) -> float:
    extent = rate_mol_s * dt
    if extent >= 0.0:
        return _bound_extent_forward(rate_mol_s, dt, N, reactants_fwd)

    max_extent = abs(extent)
    for sp, nu in reactants_rev.items():
        j = GAS_SPECIES_INDEX[sp]
        max_extent = min(max_extent, max(float(N[j]), 0.0) / max(float(nu), 1e-12))
    return -max(max_extent, 0.0)


def _solve_T_from_enthalpy(N: np.ndarray, H_target: float, T_guess: float) -> float:
    t_lo = 300.0
    t_hi = max(2500.0, float(T_guess) + 300.0)
    h_lo = calc_gas_enthalpy_flow(N, t_lo, h_cache={}) - H_target
    h_hi = calc_gas_enthalpy_flow(N, t_hi, h_cache={}) - H_target
    while h_lo * h_hi > 0.0 and t_hi < 5000.0:
        t_hi *= 1.2
        h_hi = calc_gas_enthalpy_flow(N, t_hi, h_cache={}) - H_target

    if h_lo * h_hi > 0.0:
        return float(np.clip(T_guess, t_lo, t_hi))

    for _ in range(60):
        t_mid = 0.5 * (t_lo + t_hi)
        h_mid = calc_gas_enthalpy_flow(N, t_mid, h_cache={}) - H_target
        if abs(h_mid) < 1e-6:
            return t_mid
        if h_lo * h_mid <= 0.0:
            t_hi = t_mid
            h_hi = h_mid
        else:
            t_lo = t_mid
            h_lo = h_mid
    return 0.5 * (t_lo + t_hi)


def _build_velocity_samples(u_p0: float, sigma: float, n_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """构造 freeboard 初始颗粒速度分布的离散样本。

    Hamel Eq. 3.72/3.73：床层表面抛射速度服从 Gaussian 分布，
    均值为 ``u_p0 = 1.53 * u_b``，``sigma`` 为绝对标准差 [m/s]。
    这里用对称离散样本做最小近似，只用于 solids trajectory / hold-up surrogate。
    """
    n = max(int(n_bins), 1)
    if n == 1:
        return np.array([max(float(u_p0), 1e-6)], dtype=np.float64), np.array([1.0], dtype=np.float64)

    xi = np.linspace(-2.0, 2.0, n, dtype=np.float64)
    exp_arg = np.clip(-0.5 * xi * xi, -200.0, 200.0)
    weights = np.exp(exp_arg)
    weights /= np.sum(weights)
    velocities = np.maximum(float(u_p0) + float(sigma) * xi, 1e-6)
    return velocities.astype(np.float64), weights.astype(np.float64)


def _build_size_class_bundles(
    *,
    d_p_classes: np.ndarray,
    m_char_classes: np.ndarray,
    m_ash_classes: np.ndarray,
    rho_p: float,
    eps_b: float,
    eps_d_void: float,
    u_b: float,
    d_b: float,
    area: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """透传床顶已经给出的固体流量 [kg/s]。

    ``m_char_classes`` 与 ``m_ash_classes`` 是质量流，不是库存。绝对抛射由
    床顶上行流决定，本函数不再用 roof/wake 生成第二份总量。``u_b``、``d_b``、
    ``area`` 只保留调用签名，不改变返回的质量流。
    """
    _ = (rho_p, eps_b, eps_d_void, u_b, d_b, area)
    m_char = np.maximum(np.asarray(m_char_classes, dtype=np.float64), 0.0)
    m_ash = np.maximum(np.asarray(m_ash_classes, dtype=np.float64), 0.0)
    d_p = np.maximum(np.asarray(d_p_classes, dtype=np.float64), 1e-9)
    solids = m_char + m_ash
    if solids.size == 0 or float(np.sum(solids)) <= 1e-12:
        return (
            np.array([float(np.mean(d_p)) if d_p.size else 1e-6], dtype=np.float64),
            np.array([0.5], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
            0.0,
        )
    char_frac = np.divide(m_char, solids, out=np.full_like(solids, 0.5), where=solids > 1e-12)
    return d_p, char_frac, solids, float(np.sum(solids))


def bed_face_class_flows_kg_s(cell_class_kg_s: np.ndarray) -> np.ndarray:
    """穿过床面的回落 [kg/s]。只取贴床面的第 0 格。

    ``cell_class_kg_s`` 形状为 ``(n_cells,)`` 或 ``(n_cells, n_classes)``。
    上面几格的向下流进入相邻下一格，不再加进床层。
    """
    flows = np.asarray(cell_class_kg_s, dtype=np.float64)
    if flows.size == 0:
        return np.zeros(0, dtype=np.float64)
    if flows.ndim == 1:
        return np.array([float(flows[0])], dtype=np.float64)
    return np.array(flows[0], dtype=np.float64, copy=True)


def _aggregate_class_transport(
    *,
    n_classes: int,
    class_idx_samples: np.ndarray,
    hold_samples: np.ndarray,
    upflow_samples: np.ndarray,
    downflow_samples: np.ndarray,
    char_frac_samples: np.ndarray,
    char_conversion_samples: np.ndarray,
) -> dict[str, np.ndarray]:
    hold_char = np.bincount(
        class_idx_samples,
        weights=np.maximum(hold_samples, 0.0) * np.clip(char_frac_samples, 0.0, 1.0),
        minlength=n_classes,
    ).astype(np.float64)
    hold_ash = np.bincount(
        class_idx_samples,
        weights=np.maximum(hold_samples, 0.0) * (1.0 - np.clip(char_frac_samples, 0.0, 1.0)),
        minlength=n_classes,
    ).astype(np.float64)
    up_char = np.bincount(
        class_idx_samples,
        weights=np.maximum(upflow_samples, 0.0) * np.clip(char_frac_samples, 0.0, 1.0),
        minlength=n_classes,
    ).astype(np.float64)
    up_ash = np.bincount(
        class_idx_samples,
        weights=np.maximum(upflow_samples, 0.0) * (1.0 - np.clip(char_frac_samples, 0.0, 1.0)),
        minlength=n_classes,
    ).astype(np.float64)
    down_char = np.bincount(
        class_idx_samples,
        weights=np.maximum(downflow_samples, 0.0) * np.clip(char_frac_samples, 0.0, 1.0),
        minlength=n_classes,
    ).astype(np.float64)
    down_ash = np.bincount(
        class_idx_samples,
        weights=np.maximum(downflow_samples, 0.0) * (1.0 - np.clip(char_frac_samples, 0.0, 1.0)),
        minlength=n_classes,
    ).astype(np.float64)
    hold_char_weighted_x = np.bincount(
        class_idx_samples,
        weights=np.maximum(hold_samples, 0.0) * np.clip(char_frac_samples, 0.0, 1.0) * np.clip(char_conversion_samples, 0.0, 1.0 - 1e-9),
        minlength=n_classes,
    ).astype(np.float64)
    hold_total = hold_char + hold_ash
    up_total = up_char + up_ash
    down_total = down_char + down_ash
    char_conversion_classes = np.divide(
        hold_char_weighted_x,
        np.maximum(hold_char, 0.0),
        out=np.zeros_like(hold_char_weighted_x),
        where=np.maximum(hold_char, 0.0) > 1e-12,
    )
    k_auf = np.divide(up_total, hold_total, out=np.zeros_like(up_total), where=hold_total > 1e-12)
    k_ab = np.divide(down_total, hold_total, out=np.zeros_like(down_total), where=hold_total > 1e-12)
    return {
        "m_hold_char_classes": hold_char,
        "m_hold_ash_classes": hold_ash,
        "char_conversion_classes": np.clip(char_conversion_classes, 0.0, 1.0 - 1e-9),
        "m_dot_auf_char_classes": up_char,
        "m_dot_auf_ash_classes": up_ash,
        "m_dot_ab_char_classes": down_char,
        "m_dot_ab_ash_classes": down_ash,
        "K_auf_classes": k_auf,
        "K_ab_classes": k_ab,
    }


def _resolve_transport_class_mapping(
    *,
    d_p_class: np.ndarray,
    m_dot_class: np.ndarray,
    class_aggregate_indices: np.ndarray | None,
    n_output_size_classes: int | None,
) -> tuple[np.ndarray, int, np.ndarray]:
    """Return aggregate class ids, output width, and representative diameters."""
    n_launch = int(np.asarray(d_p_class).size)
    if class_aggregate_indices is None:
        aggregate_idx = np.arange(n_launch, dtype=np.int64)
    else:
        aggregate_idx = np.asarray(class_aggregate_indices, dtype=np.int64)
        if aggregate_idx.shape != (n_launch,):
            raise ValueError(
                "class_aggregate_indices_bed_top must match expanded launch class count"
            )
        if np.any(aggregate_idx < 0):
            raise ValueError("class_aggregate_indices_bed_top must be non-negative")
    n_output = (
        int(n_output_size_classes)
        if n_output_size_classes is not None
        else (int(np.max(aggregate_idx)) + 1 if aggregate_idx.size else 0)
    )
    if n_output < 0 or (aggregate_idx.size and int(np.max(aggregate_idx)) >= n_output):
        raise ValueError("n_output_size_classes must cover all aggregate class indices")
    weights = np.maximum(np.asarray(m_dot_class, dtype=np.float64), 0.0)
    d_p = np.maximum(np.asarray(d_p_class, dtype=np.float64), 1e-9)
    denom = np.bincount(aggregate_idx, weights=weights, minlength=n_output).astype(np.float64)
    numer = np.bincount(aggregate_idx, weights=weights * d_p, minlength=n_output).astype(np.float64)
    fallback = float(np.mean(d_p)) if d_p.size else 1e-6
    d_p_out = np.divide(numer, denom, out=np.full(n_output, fallback, dtype=np.float64), where=denom > 1e-12)
    return aggregate_idx, n_output, np.maximum(d_p_out, 1e-9)


def _empty_class_transport(n_classes: int) -> dict[str, np.ndarray]:
    zeros = np.zeros(int(n_classes), dtype=np.float64)
    return {
        "m_hold_char_classes": zeros.copy(),
        "m_hold_ash_classes": zeros.copy(),
        "char_conversion_classes": zeros.copy(),
        "m_char_sink_applied_classes": zeros.copy(),
        "m_dot_auf_char_classes": zeros.copy(),
        "m_dot_auf_ash_classes": zeros.copy(),
        "m_dot_ab_char_classes": zeros.copy(),
        "m_dot_ab_ash_classes": zeros.copy(),
        "K_auf_classes": zeros.copy(),
        "K_ab_classes": zeros.copy(),
    }


def _representative_char_conversion(class_transport: dict[str, np.ndarray]) -> float:
    hold_char = np.asarray(class_transport.get("m_hold_char_classes", []), dtype=np.float64)
    x_classes = np.asarray(class_transport.get("char_conversion_classes", np.zeros_like(hold_char)), dtype=np.float64)
    total_hold_char = float(np.sum(np.maximum(hold_char, 0.0)))
    if total_hold_char <= 1e-12 or hold_char.size == 0:
        return 0.0
    return float(
        np.clip(
            np.sum(np.maximum(hold_char, 0.0) * np.clip(x_classes, 0.0, 1.0 - 1e-9)) / total_hold_char,
            0.0,
            1.0 - 1e-9,
        )
    )


def _apply_char_sink_to_class_transport(
    class_transport: dict[str, np.ndarray],
    *,
    diag_acc: dict[str, float],
) -> dict[str, np.ndarray]:
    hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
    total_hold_char = float(np.sum(np.maximum(hold_char, 0.0)))
    sink_total_kg = (
        max(float(diag_acc.get("R1", 0.0)), 0.0)
        + max(float(diag_acc.get("R2", 0.0)), 0.0)
        + max(float(diag_acc.get("R3", 0.0)), 0.0)
        + max(float(diag_acc.get("R4", 0.0)), 0.0)
    ) * ATOMIC_MASS_KG_PER_MOL["C"]
    if sink_total_kg <= 0.0 or total_hold_char <= 1e-12:
        class_transport["m_char_sink_applied_classes"] = np.zeros_like(hold_char)
        return class_transport

    weights = np.divide(
        np.maximum(hold_char, 0.0),
        total_hold_char,
        out=np.zeros_like(hold_char),
        where=np.maximum(hold_char, 0.0) > 0.0,
    )
    sink_by_class = np.minimum(np.maximum(hold_char, 0.0), sink_total_kg * weights)
    x_old = np.asarray(class_transport.get("char_conversion_classes", np.zeros_like(hold_char)), dtype=np.float64)
    frac_burn = np.divide(
        sink_by_class,
        np.maximum(hold_char, 0.0),
        out=np.zeros_like(sink_by_class),
        where=np.maximum(hold_char, 0.0) > 1e-12,
    )
    x_new = np.clip(x_old + frac_burn * (1.0 - x_old), 0.0, 1.0 - 1e-9)
    class_transport["m_hold_char_classes"] = np.maximum(hold_char - sink_by_class, 0.0)
    class_transport["char_conversion_classes"] = x_new
    class_transport["m_char_sink_applied_classes"] = sink_by_class

    hold_char_new = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
    hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)
    up_char = np.asarray(class_transport["m_dot_auf_char_classes"], dtype=np.float64)
    up_ash = np.asarray(class_transport["m_dot_auf_ash_classes"], dtype=np.float64)
    down_char = np.asarray(class_transport["m_dot_ab_char_classes"], dtype=np.float64)
    down_ash = np.asarray(class_transport["m_dot_ab_ash_classes"], dtype=np.float64)
    hold_total = hold_char_new + hold_ash
    up_total = up_char + up_ash
    down_total = down_char + down_ash
    class_transport["K_auf_classes"] = np.divide(up_total, hold_total, out=np.zeros_like(up_total), where=hold_total > 1e-12)
    class_transport["K_ab_classes"] = np.divide(down_total, hold_total, out=np.zeros_like(down_total), where=hold_total > 1e-12)
    return class_transport


def _freeboard_inventory_char_sink_diag(
    *,
    N: np.ndarray,
    T: float,
    P: float,
    V_seg: float,
    dt: float,
    fuel_type: str,
    char_area_total: float,
    solid_d_p: float | None,
    D_g: float | None,
    char_conversion: float,
    m_char_total: float = 0.0,
    r4_scale: float = 1.0,
) -> dict[str, float]:
    """Compute closure-owned char depletion extents without altering closure gas state.

    The freeboard closure owns ``m_hold_char_classes``. Even when closure-layer
    hetero chemistry is disabled for gas/energy coupling, we still need a
    thesis-consistent inventory-depletion pass at the same architectural layer.
    Reuse the exact same R1-R4 kinetics kernel, but keep its effect local to
    hold-up depletion instead of feeding closure gas/energy balances.
    """
    if (
        (char_area_total <= 0.0 and float(m_char_total) <= 0.0)
        or solid_d_p is None
        or D_g is None
        or dt <= 0.0
    ):
        return {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 0.0}
    from src.core.freeboard.reaction import _freeboard_reaction_step

    _, diag = _freeboard_reaction_step(
        N=np.array(N, dtype=np.float64, copy=True),
        T=T,
        P=P,
        V_seg=V_seg,
        dt=dt,
        fuel_type=fuel_type,
        rho_cat=0.0,
        char_area_total=char_area_total,
        m_char_total=float(m_char_total),
        solid_d_p=solid_d_p,
        D_g=D_g,
        char_conversion=char_conversion,
        r4_scale=r4_scale,
        enabled_reactions=("R1", "R2", "R3", "R4"),
    )
    return {name: float(diag.get(name, 0.0)) for name in ("R1", "R2", "R3", "R4")}

