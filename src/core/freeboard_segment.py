"""床层顶部之后的 freeboard 轴向段近似。

当前实现目标：
- 在 bed solve 完成后，追加气相主导的 freeboard 段
- 区分 `bed exit` 与 `reactor exit`
- 用 Hamel Chapter 3 的 `u_gb` / `beta_A` 口径给出轴向停留时间

说明：
- 这是对当前仓库“无 freeboard 拓扑”的最小集成，不等价于论文完整的
  Verbindungsmatrix + cyclone / duct cells。
- 当前版本已加入一个最小的 entrained-solids 骨架：
  用 bed-top 抛射通量估算 freeboard 内 `char/ash` 稀相持量，并据此启用
  `R11` 的催化路径。
- 默认反应侧仍以 homogeneous 路径为主：
  `R5, R6, R7, R10, R11b/s, R12`。
- `R8` 当前实现对应煤灰/焦催化 WGSR。基础 ``ReactorConfig`` 默认反应集不含 `R8`，
  但 thesis-mode / HTW freeboard builder 会显式把 `R8` 加回，用于与当前 Hamel
  对齐审计保持一致。
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src.core.cell_balances import calc_gas_enthalpy_flow
from src.core.constants import g
from src.core.constants import Rg
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, N_GAS, gas_density_ideal, gas_diffusivity_correlation, gas_viscosity_power_law
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R6, rate_R7, rate_R8, rate_R12
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from src.physics.freeboard import calc_beta_a, calc_cd_hamel_eq384, calc_cd_haider, calc_terminal_velocity_haider, calc_u_gb

_F_W = 0.25
_ZETA_W = 0.40


@dataclass
class FreeboardState:
    segment_index: int
    z_center_m: float
    T: float
    N: np.ndarray
    u0: float
    u_gb: float
    u_p_mean: float
    u_t_mean: float
    carry_ratio: float
    tau: float
    m_dot_char: float
    m_dot_ash: float
    m_dot_return_char: float
    m_dot_return_ash: float
    m_hold_char_before_classes: np.ndarray
    m_hold_ash_before_classes: np.ndarray
    char_conversion_before_classes: np.ndarray
    char_conversion_classes: np.ndarray
    m_hold_char_classes: np.ndarray
    m_hold_ash_classes: np.ndarray
    m_char_sink_applied_classes: np.ndarray
    m_dot_auf_char_classes: np.ndarray
    m_dot_auf_ash_classes: np.ndarray
    m_dot_ab_char_classes: np.ndarray
    m_dot_ab_ash_classes: np.ndarray
    K_auf_classes: np.ndarray
    K_ab_classes: np.ndarray
    rho_solid: float
    rho_cat: float
    reaction_diag: dict[str, float]


@dataclass
class FreeboardBackground:
    segment_index: int
    z_lower_m: float
    z_upper_m: float
    z_center_m: float
    u0: float
    u_gb_in: float
    u_gb_out: float
    u_gb: float
    u_gb_center: float
    tau: float
    V_seg: float
    rho_g: float
    mu_g: float


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

    solid_frac = solids / np.sum(solids)
    d_p_mean = float(np.sum(d_p * solid_frac))
    roof = np.full_like(
        d_p,
        d_p_mean * rho_p * max(1.0 - eps_d_void, 0.0) * (3.0 * eps_b * u_b / max(d_b, 1e-12)),
    )
    wake = np.full_like(roof, _ZETA_W * _F_W * eps_b * u_b * rho_p * max(1.0 - eps_d_void, 0.0))
    flux_raw = (roof + wake) * solid_frac
    m_dot_class_raw = flux_raw * area
    total_feed = float(np.sum(solids))
    total_raw = float(np.sum(m_dot_class_raw))
    if total_feed > 0.0 and total_raw > 1e-12:
        # Freeboard closure must remain conservative with respect to the bed-top
        # solids actually available to be ejected. Hamel Eq. 3.67 provides the
        # roof/wake structure for the ejection split, but under the current
        # reactor architecture the absolute freeboard source term must close to
        # the real upward solids feed from the top bed cell, otherwise the
        # closure can create more solids than the bed supplies.
        m_dot_class = m_dot_class_raw * (total_feed / total_raw)
    else:
        m_dot_class = np.zeros_like(m_dot_class_raw)
    char_frac = np.divide(m_char, solids, out=np.full_like(solids, 0.5), where=solids > 1e-12)
    return d_p, char_frac, m_dot_class, float(np.sum(m_dot_class))


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
    r4_scale: float = 1.0,
) -> dict[str, float]:
    """Compute closure-owned char depletion extents without altering closure gas state.

    The freeboard closure owns ``m_hold_char_classes``. Even when closure-layer
    hetero chemistry is disabled for gas/energy coupling, we still need a
    thesis-consistent inventory-depletion pass at the same architectural layer.
    Reuse the exact same R1-R4 kinetics kernel, but keep its effect local to
    hold-up depletion instead of feeding closure gas/energy balances.
    """
    if char_area_total <= 0.0 or solid_d_p is None or D_g is None or dt <= 0.0:
        return {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 0.0}
    _, diag = _freeboard_reaction_step(
        N=np.array(N, dtype=np.float64, copy=True),
        T=T,
        P=P,
        V_seg=V_seg,
        dt=dt,
        fuel_type=fuel_type,
        rho_cat=0.0,
        char_area_total=char_area_total,
        solid_d_p=solid_d_p,
        D_g=D_g,
        char_conversion=char_conversion,
        r4_scale=r4_scale,
        enabled_reactions=("R1", "R2", "R3", "R4"),
    )
    return {name: float(diag.get(name, 0.0)) for name in ("R1", "R2", "R3", "R4")}


def _force_balance_particle_accel(
    *,
    u_p: float,
    u_g: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    phi_s: float,
    mu_g: float,
) -> float:
    rel = float(u_g) - float(u_p)
    re_p = max(rho_g * abs(rel) * max(d_p, 1e-9) / max(mu_g, 1e-12), 1e-9)
    cd = calc_cd_haider(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
    drag = 3.0 * cd * rho_g * abs(rel) * rel / max(4.0 * rho_p * max(d_p, 1e-9), 1e-12)
    return -g * (1.0 - rho_g / max(rho_p, 1e-12)) + drag


def _advance_particle_samples_surrogate(
    *,
    dh: float,
    beta_a: float,
    u_g_asym: float,
    u_g_prev: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    exp_arg_u = float(np.clip(-beta_a * dh, -200.0, 200.0))
    u_p_samples = u_g_asym + (u_p_samples_prev - u_g_asym) * math.exp(exp_arg_u)
    decay_ratio = np.clip(u_g_prev / np.maximum(u_p_samples_prev, 1e-6), 0.05, 10.0)
    exp_arg = np.clip(-beta_a * dh * decay_ratio, -200.0, 200.0)
    m_dot_samples = m_dot_samples_prev * np.exp(exp_arg)
    tau_particles = dh / np.maximum(0.5 * (u_p_samples_prev + u_p_samples), 1e-6)
    hold_samples = 0.5 * (m_dot_samples_prev + m_dot_samples) * tau_particles
    return u_p_samples, m_dot_samples, hold_samples, np.zeros_like(m_dot_samples)


def _advance_particle_samples_force_balance(
    *,
    dh: float,
    u_g_prev: float,
    u_g_next: float,
    rho_g: float,
    rho_p: float,
    d_p: np.ndarray,
    phi_s: float,
    mu_g: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
    char_frac_samples: np.ndarray,
    n_substeps: int = 32,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_steps = max(int(n_substeps), 4)
    u_p_samples = np.array(u_p_samples_prev, dtype=np.float64, copy=True)
    m_dot_samples = np.array(m_dot_samples_prev, dtype=np.float64, copy=True)
    hold_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    return_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    dt_cap = 0.25 * dh / max(max(abs(u_g_prev), abs(u_g_next), float(np.max(np.abs(u_p_samples_prev)))), 1e-3)
    dt_floor = 1e-5

    for j in range(len(u_p_samples)):
        if m_dot_samples_prev[j] <= 0.0:
            u_p_samples[j] = max(u_p_samples_prev[j], 1e-6)
            m_dot_samples[j] = 0.0
            continue

        z = 0.0
        t_seg = 0.0
        u_p = max(float(u_p_samples_prev[j]), 1e-6)
        survived = False

        for _ in range(n_steps * 12):
            frac = float(np.clip(z / max(dh, 1e-12), 0.0, 1.0))
            u_g = u_g_prev + (u_g_next - u_g_prev) * frac
            scale_v = max(abs(u_p), abs(u_g), 1e-3)
            dt = float(np.clip((dh / n_steps) / scale_v, dt_floor, dt_cap))
            accel = _force_balance_particle_accel(
                u_p=u_p,
                u_g=u_g,
                rho_g=rho_g,
                rho_p=rho_p,
                d_p=float(d_p[j]),
                phi_s=phi_s,
                mu_g=mu_g,
            )
            u_new = u_p + accel * dt
            u_avg = 0.5 * (u_p + u_new)
            dz = max(u_avg, 0.0) * dt

            if z + dz >= dh:
                dz_remain = dh - z
                dt_cross = dz_remain / max(u_avg, 1e-6)
                t_seg += dt_cross
                z = dh
                u_p = max(u_new, 1e-6)
                survived = True
                break

            z += dz
            t_seg += dt
            u_p = u_new
            if u_p <= 0.0 and z < dh:
                u_p = 0.0
                break

        hold_samples[j] = float(m_dot_samples_prev[j]) * max(t_seg, 0.0)
        if survived:
            u_p_samples[j] = max(u_p, 1e-6)
        else:
            return_samples[j] = float(m_dot_samples_prev[j])
            u_p_samples[j] = 0.0
            m_dot_samples[j] = 0.0

    return u_p_samples, m_dot_samples, hold_samples, return_samples


def _haider_high_re_drag_constant(phi_s: float) -> float:
    """High-Re constant term z in Wirsum-style c_w ~= 24/Re + z."""
    phi = float(np.clip(phi_s, 1e-3, 1.0))
    return math.exp(4.9050 - 13.8944 * phi + 18.4222 * (phi**2) - 10.2599 * (phi**3))


def _hamel_exact_high_re_drag_constant(phi_s: float) -> float:
    """Hamel thesis high-Re plateau ``z`` from the simplified Haider form.

    Hamel Eq. 3.85/3.86 uses ``c_w ~= 24/Re + z`` and defines ``z`` as the
    very-high-Re asymptotic drag level. For the thesis-aligned exact audit
    path we therefore use the explicit plateau term cited from the dissertation:

        z = 73.69 * exp(-5.0784 * phi_s)
    """
    phi = float(np.clip(phi_s, 1e-3, 1.0))
    return 73.69 * math.exp(-5.0784 * phi)


def _wirsum_abc(
    *,
    u_g: float,
    u_p: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    phi_s: float,
    branch: str | None = None,
) -> tuple[float, float, float]:
    """Construct A/B/C for du_p/dt = A*u_p^2 + B*u_p + C."""
    d_eff = max(float(d_p), 1e-9)
    rho_p_eff = max(float(rho_p), 1e-9)
    z = _haider_high_re_drag_constant(float(phi_s))
    branch_mode = branch
    if branch_mode is None:
        branch_mode = "slower" if float(u_p) < float(u_g) else "faster"
    s = -1.0 if branch_mode == "slower" else 1.0
    k_z = 3.0 * z * float(rho_g) / (4.0 * rho_p_eff * d_eff)
    k_mu = 18.0 * float(mu_g) / (rho_p_eff * d_eff * d_eff)
    if branch_mode == "slower":
        a = -s * k_z
        b = -k_mu + s * 2.0 * k_z * float(u_g)
        c = k_mu * float(u_g) - s * k_z * float(u_g) * float(u_g) - g * (1.0 - float(rho_g) / rho_p_eff)
    elif branch_mode == "faster":
        a = -s * k_z
        b = -k_mu + s * 2.0 * k_z * float(u_g)
        c = k_mu * float(u_g) - s * k_z * float(u_g) * float(u_g) - g * (1.0 - float(rho_g) / rho_p_eff)
    else:
        raise ValueError(f"Unknown Wirsum branch mode: {branch_mode}")
    return float(a), float(b), float(c)


def _wirsum_abc_exact_audit(
    *,
    u_g: float,
    u_p: float,
    rho_g: float,
    rho_p: float,
    d_p: float,
    mu_g: float,
    phi_s: float,
    branch: str | None = None,
) -> tuple[float, float, float, float, float, float]:
    """Return Hamel-style exact audit coefficients plus ``(Re_p, c_w, z_drag)``.

    This function is diagnostic-only for now. The production trajectory solver
    still uses ``_wirsum_abc()`` until the exact path is validated against the
    current regression suite.

    Important Hamel nuance:
    the sign set for ``A/B/C`` is selected from the instantaneous velocity
    branch ``u_p < u_g`` vs. ``u_p > u_g``. The resulting solution family
    (``Delta > 0`` trig branch vs. ``Delta < 0`` hyperbolic branch) is then
    determined only after evaluating ``Delta = 4AC - B^2``. In particular,
    ``u_p > u_g`` does not guarantee ``Delta > 0`` for small particles below
    the gravity-limited critical size.
    """
    d_eff = max(float(d_p), 1e-9)
    rho_p_eff = max(float(rho_p), 1e-9)
    branch_mode = branch
    if branch_mode is None:
        branch_mode = "slower" if float(u_p) < float(u_g) else "faster"

    rel = abs(float(u_g) - float(u_p))
    re_p = max(float(rho_g) * rel * d_eff / max(float(mu_g), 1e-12), 1e-9)
    c_w = calc_cd_hamel_eq384(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
    z_drag = _hamel_exact_high_re_drag_constant(float(phi_s))
    k_z = 3.0 * z_drag * float(rho_g) / (4.0 * rho_p_eff * d_eff)
    k_mu = 18.0 * float(mu_g) / (rho_p_eff * d_eff * d_eff)
    grav = g * (1.0 - float(rho_g) / rho_p_eff)

    # Hamel Eq. 3.85 splits c_w into a viscous Stokes contribution and the
    # high-Re plateau z. The Riccati quadratic term therefore comes only from
    # z * u_rel^2, while the 24/Re part reduces to a linear drag term.
    if branch_mode == "slower":
        a = k_z
        b = -2.0 * k_z * float(u_g) - k_mu
        c = k_z * float(u_g) * float(u_g) + k_mu * float(u_g) - grav
    elif branch_mode == "faster":
        a = -k_z
        b = 2.0 * k_z * float(u_g) - k_mu
        c = -k_z * float(u_g) * float(u_g) + k_mu * float(u_g) - grav
    else:
        raise ValueError(f"Unknown Wirsum branch mode: {branch_mode}")

    return float(a), float(b), float(c), float(re_p), float(c_w), float(z_drag)


def _estimate_exact_hamel_critical_diameter(
    *,
    u_g: float,
    u_p: float,
    rho_g: float,
    rho_p: float,
    mu_g: float,
    phi_s: float,
    branch: str,
    d_lo: float = 1e-5,
    d_hi: float = 5e-3,
    n_scan: int = 48,
) -> float | None:
    """Estimate the critical diameter where the exact Hamel audit has ``Delta = 0``."""
    if branch != "faster":
        return None

    def delta_at(d_p: float) -> float:
        a, b, c, _, _, _ = _wirsum_abc_exact_audit(
            u_g=u_g,
            u_p=u_p,
            rho_g=rho_g,
            rho_p=rho_p,
            d_p=d_p,
            mu_g=mu_g,
            phi_s=phi_s,
            branch=branch,
        )
        return float(4.0 * a * c - b * b)

    d_prev = float(d_lo)
    f_prev = delta_at(d_prev)
    if abs(f_prev) <= 1e-12:
        return d_prev
    for i in range(1, n_scan + 1):
        frac = i / float(n_scan)
        d_cur = float(d_lo * ((d_hi / d_lo) ** frac))
        f_cur = delta_at(d_cur)
        if abs(f_cur) <= 1e-12:
            return d_cur
        if f_prev * f_cur < 0.0:
            lo, hi = d_prev, d_cur
            flo, fhi = f_prev, f_cur
            for _ in range(60):
                mid = math.sqrt(lo * hi)
                fmid = delta_at(mid)
                if abs(fmid) <= 1e-12 or abs(hi - lo) <= 1e-10 * max(1.0, hi):
                    return float(mid)
                if flo * fmid <= 0.0:
                    hi, fhi = mid, fmid
                else:
                    lo, flo = mid, fmid
            return float(math.sqrt(lo * hi))
        d_prev, f_prev = d_cur, f_cur
    return None


def _quadratic_ode_step(u0: float, a: float, b: float, c: float, dt: float) -> tuple[float, float]:
    """Analytical step for du/dt = a*u^2 + b*u + c.

    Returns ``(u1, delta_thesis)`` where ``delta_thesis = 4*a*c - b^2``.
    """
    delta = 4.0 * a * c - b * b
    if dt <= 0.0:
        return float(u0), float(delta)

    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            return float(u0 + c * dt), float(delta)
        exp_bt = math.exp(float(np.clip(b * dt, -80.0, 80.0)))
        u1 = u0 * exp_bt + (c / b) * (exp_bt - 1.0)
        return float(u1), float(delta)

    if delta < -1e-14:
        # Numerically stable root/exponential form; corresponds to thesis Delta < 0.
        k = math.sqrt(-delta)
        r1 = (-b - k) / (2.0 * a)
        r2 = (-b + k) / (2.0 * a)
        den0 = (u0 - r2)
        if abs(den0) < 1e-14:
            den0 = 1e-14 if den0 >= 0.0 else -1e-14
        kappa = (u0 - r1) / den0
        ratio = kappa * math.exp(float(np.clip(-k * dt, -80.0, 80.0)))
        den = 1.0 - ratio
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = (r1 - ratio * r2) / den
        return float(u1), float(delta)

    if abs(delta) <= 1e-14:
        r = -b / (2.0 * a)
        y0 = u0 - r
        den = 1.0 - a * y0 * dt
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = r + y0 / den
        return float(u1), float(delta)

    # Numerically stable trigonometric form; corresponds to thesis Delta > 0.
    omega = math.sqrt(delta)
    y0 = (2.0 * a * u0 + b) / omega
    theta = math.atan(y0) + 0.5 * omega * dt
    half_pi = 0.5 * math.pi
    while theta > half_pi:
        theta -= math.pi
    while theta < -half_pi:
        theta += math.pi
    theta = float(np.clip(theta, -half_pi + 1e-6, half_pi - 1e-6))
    u1 = (-b + omega * math.tan(theta)) / (2.0 * a)
    return float(u1), float(delta)


def _quadratic_ode_state(u0: float, a: float, b: float, c: float, dt: float) -> tuple[float, float, float]:
    """Return ``(u(t), z_delta(t), delta_thesis)`` for ``du/dt = a*u^2 + b*u + c``.

    The height increment is the exact integral of u over [0, dt].
    """
    delta = 4.0 * a * c - b * b
    if dt <= 0.0:
        return float(u0), 0.0, float(delta)

    if abs(a) < 1e-14:
        if abs(b) < 1e-14:
            u1 = float(u0 + c * dt)
            z1 = float(u0 * dt + 0.5 * c * dt * dt)
            return u1, z1, float(delta)
        exp_bt = math.exp(float(np.clip(b * dt, -80.0, 80.0)))
        u1 = float(u0 * exp_bt + (c / b) * (exp_bt - 1.0))
        z1 = float((u0 + c / b) * (exp_bt - 1.0) / b - (c / b) * dt)
        return u1, z1, float(delta)

    if delta < -1e-14:
        k = math.sqrt(-delta)
        r1 = (-b - k) / (2.0 * a)
        r2 = (-b + k) / (2.0 * a)
        den0 = u0 - r2
        if abs(den0) < 1e-14:
            den0 = 1e-14 if den0 >= 0.0 else -1e-14
        q0 = (u0 - r1) / den0
        q = q0 * math.exp(float(np.clip(-k * dt, -80.0, 80.0)))
        den = 1.0 - q
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = float((r1 - q * r2) / den)
        den_q0 = 1.0 - q0
        den_q = 1.0 - q
        if abs(den_q0) < 1e-14:
            den_q0 = 1e-14 if den_q0 >= 0.0 else -1e-14
        if abs(den_q) < 1e-14:
            den_q = 1e-14 if den_q >= 0.0 else -1e-14
        z1 = float(r1 * dt + ((r1 - r2) / k) * math.log(abs(den_q / den_q0)))
        return u1, z1, float(delta)

    if abs(delta) <= 1e-14:
        r = -b / (2.0 * a)
        y0 = u0 - r
        den = 1.0 - a * y0 * dt
        if abs(den) < 1e-14:
            den = 1e-14 if den >= 0.0 else -1e-14
        u1 = float(r + y0 / den)
        z1 = float(r * dt - math.log(abs(den)) / a)
        return u1, z1, float(delta)

    omega = math.sqrt(delta)
    theta0 = math.atan((2.0 * a * u0 + b) / omega)
    theta = theta0 + 0.5 * omega * dt
    half_pi = 0.5 * math.pi
    while theta > half_pi:
        theta -= math.pi
    while theta < -half_pi:
        theta += math.pi
    theta = float(np.clip(theta, -half_pi + 1e-6, half_pi - 1e-6))
    u1 = float((-b + omega * math.tan(theta)) / (2.0 * a))
    z1 = float((-b / (2.0 * a)) * dt + (1.0 / a) * math.log(abs(math.cos(theta0) / max(abs(math.cos(theta)), 1e-14))))
    return u1, z1, float(delta)


def _hamel_delta_neg_time_between_velocities(u_start: float, u_end: float, a: float, b: float, c: float) -> float | None:
    """Hamel Fall 2: exact time difference for ``Delta < 0`` between two signed velocities."""
    delta = 4.0 * a * c - b * b
    if delta >= -1e-14:
        return None
    k = math.sqrt(-delta)

    def primitive(u: float) -> float:
        arg = (2.0 * a * float(u) + b) / max(k, 1e-14)
        arg = float(np.clip(arg, -1.0 + 1e-12, 1.0 - 1e-12))
        return float((-2.0 / k) * np.arctanh(arg))

    dt = primitive(float(u_end)) - primitive(float(u_start))
    if not math.isfinite(dt):
        return None
    return float(abs(dt))


def _hamel_delta_neg_height_from_velocity(
    *,
    u: float,
    u0: float,
    h0: float,
    a: float,
    b: float,
    c: float,
) -> float | None:
    """Hamel Eq. 3.103/3.104/3.105: height as a function of signed velocity for ``Delta < 0``."""
    delta = 4.0 * a * c - b * b
    if delta >= -1e-14 or abs(a) < 1e-14:
        return None
    k = math.sqrt(-delta)
    psi = 2.0 * a * float(u) + b
    psi0 = 2.0 * a * float(u0) + b
    q = a * float(u) * float(u) + b * float(u) + c
    q0 = a * float(u0) * float(u0) + b * float(u0) + c
    if abs(q0) < 1e-14:
        return None

    term1 = (1.0 / (2.0 * a)) * math.log(max(abs(q / q0), 1e-300))
    num = (psi - k) * (psi0 + k)
    den = (psi + k) * (psi0 - k)
    if abs(den) < 1e-14:
        return None
    term2 = -(b / (2.0 * a * k)) * math.log(max(abs(num / den), 1e-300))
    h = float(h0 + term1 + term2)
    if not math.isfinite(h):
        return None
    return h


def _hamel_delta_neg_max_height(
    *,
    u0: float,
    h0: float,
    a: float,
    b: float,
    c: float,
) -> float | None:
    """Maximum reachable height in Hamel Fall 2, obtained from Eq. 3.103 at ``u = 0``."""
    return _hamel_delta_neg_height_from_velocity(
        u=0.0,
        u0=u0,
        h0=h0,
        a=a,
        b=b,
        c=c,
    )


def _hamel_delta_neg_real_roots(a: float, b: float, c: float) -> tuple[float, float] | None:
    """Return ordered real roots of ``A u^2 + B u + C = 0`` for Fall 2."""
    delta = 4.0 * a * c - b * b
    if delta >= -1e-14 or abs(a) < 1e-14:
        return None
    k = math.sqrt(-delta)
    r1 = float((-b - k) / (2.0 * a))
    r2 = float((-b + k) / (2.0 * a))
    return (min(r1, r2), max(r1, r2))


def _hamel_delta_neg_zero_reachable(a: float, b: float, c: float) -> bool:
    """Whether Fall 2 can physically reach ``u_p = 0`` in real time."""
    roots = _hamel_delta_neg_real_roots(a, b, c)
    if roots is None:
        return False
    lo, hi = roots
    return lo <= 0.0 <= hi


def _hamel_delta_neg_boundary_velocity(
    *,
    h_target: float,
    u0: float,
    h0: float,
    a: float,
    b: float,
    c: float,
    u_lo: float = -25.0,
    u_hi: float = -1e-9,
    n_scan: int = 96,
) -> float | None:
    """Solve Hamel Eq. 3.103 for a boundary velocity inside a signed bracket."""
    delta = 4.0 * a * c - b * b
    if delta >= -1e-14:
        return None

    def f_vel(u: float) -> float | None:
        h = _hamel_delta_neg_height_from_velocity(u=u, u0=u0, h0=h0, a=a, b=b, c=c)
        if h is None:
            return None
        return float(h - h_target)

    u_lo = float(u_lo)
    u_hi = float(u_hi)
    if abs(u_hi - u_lo) <= 1e-12:
        return None

    if u_lo > u_hi:
        u_lo, u_hi = u_hi, u_lo

    if u_lo > 0.0:
        grid = np.geomspace(max(u_lo, 1e-9), max(u_hi, u_lo * (1.0 + 1e-9)), int(max(n_scan, 8)))
        u_candidates = [float(v) for v in grid[1:]]
    elif u_hi < 0.0:
        grid = np.geomspace(max(abs(u_hi), 1e-9), max(abs(u_lo), 1e-9), int(max(n_scan, 8)))
        u_candidates = [-float(v) for v in grid[1:]]
    else:
        u_candidates = np.linspace(u_lo, u_hi, int(max(n_scan, 8)))[1:].astype(float).tolist()

    u_prev = float(u_lo)
    f_prev = f_vel(u_prev)
    if f_prev is not None and abs(f_prev) <= 1e-10:
        return u_prev

    for u_cur in u_candidates:
        f_cur = f_vel(float(u_cur))
        if f_cur is None:
            continue
        if abs(f_cur) <= 1e-10:
            return float(u_cur)
        if f_prev is not None and f_prev * f_cur < 0.0:
            lo = float(u_prev)
            hi = float(u_cur)
            flo = float(f_prev)
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                fmid = f_vel(mid)
                if fmid is None:
                    return None
                if abs(fmid) <= 1e-10 or abs(hi - lo) <= 1e-10 * max(1.0, abs(hi), abs(lo)):
                    return float(mid)
                if flo * fmid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    flo = float(fmid)
            return float(0.5 * (lo + hi))
        u_prev = float(u_cur)
        f_prev = f_cur
    return None


def _bisect_event_time(
    fn,
    *,
    t_lo: float,
    t_hi: float,
    target: float = 0.0,
    max_iter: int = 80,
) -> float | None:
    f_lo = float(fn(t_lo) - target)
    f_hi = float(fn(t_hi) - target)
    if abs(f_lo) <= 1e-10:
        return float(t_lo)
    if abs(f_hi) <= 1e-10:
        return float(t_hi)
    if f_lo * f_hi > 0.0:
        return None
    lo = float(t_lo)
    hi = float(t_hi)
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = float(fn(mid) - target)
        if abs(f_mid) <= 1e-10 or abs(hi - lo) <= 1e-10 * max(1.0, hi):
            return float(mid)
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return float(0.5 * (lo + hi))


def _find_first_event_time(
    fn,
    *,
    target: float,
    t_max: float,
    n_scan: int = 64,
) -> float | None:
    t_prev = 0.0
    f_prev = float(fn(t_prev) - target)
    if abs(f_prev) <= 1e-10:
        return 0.0
    for i in range(1, n_scan + 1):
        t_cur = t_max * i / float(n_scan)
        f_cur = float(fn(t_cur) - target)
        if abs(f_cur) <= 1e-10:
            return float(t_cur)
        if f_prev * f_cur < 0.0:
            return _bisect_event_time(fn, t_lo=t_prev, t_hi=t_cur, target=target)
        t_prev = t_cur
        f_prev = f_cur
    return None


def _find_first_event_time_with_growth(
    fn,
    *,
    target: float,
    t_init: float,
    t_cap: float,
    growth: float = 2.0,
    n_scan: int = 64,
) -> float | None:
    """Search for the first event while gradually expanding the time horizon."""
    t_max = max(float(t_init), 1e-9)
    t_cap = max(float(t_cap), t_max)
    while True:
        t_evt = _find_first_event_time(fn, target=target, t_max=t_max, n_scan=n_scan)
        if t_evt is not None:
            return float(t_evt)
        if t_max >= t_cap - 1e-12:
            return None
        t_max = min(float(t_cap), max(t_max * max(float(growth), 1.1), t_max + 1e-6))


def _advance_particle_samples_analytical_wirsum(
    *,
    dh: float,
    z_lower_abs: float = 0.0,
    u_g_prev: float,
    u_g_next: float,
    rho_g: float,
    rho_p: float,
    d_p: np.ndarray,
    phi_s: float,
    mu_g: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
    char_frac_samples: np.ndarray,
    slower_ref_u_prev: np.ndarray | None = None,
    slower_ref_z_prev: np.ndarray | None = None,
    coeff_model: str = "stable_mixed_drag_split",
    n_subsegments: int = 6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Advance solids with boundary-event analytical trajectory updates.

    ``n_subsegments`` is intentionally ignored now. Hamel/Wirsum uses
    boundary-to-boundary analytical residence times rather than internal
    pseudo-spatial stepping within a freeboard cell.
    """
    u_p_samples = np.array(u_p_samples_prev, dtype=np.float64, copy=True)
    m_dot_samples = np.array(m_dot_samples_prev, dtype=np.float64, copy=True)
    hold_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    return_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    slower_ref_u_samples = np.array(
        u_p_samples_prev if slower_ref_u_prev is None else slower_ref_u_prev,
        dtype=np.float64,
        copy=True,
    )
    slower_ref_z_samples = np.array(
        np.zeros_like(u_p_samples_prev) if slower_ref_z_prev is None else slower_ref_z_prev,
        dtype=np.float64,
        copy=True,
    )
    diag = {
        "delta_pos": 0,
        "delta_zero": 0,
        "delta_neg": 0,
        "disc_pos": 0,
        "disc_zero": 0,
        "disc_neg": 0,
        "sign_switch": 0,
        "fallbacks": 0,
        "returns": 0,
        "fallback_no_event": 0,
        "fallback_boundary_velocity": 0,
        "fallback_delta_neg_time": 0,
        "re_p_min": float("inf"),
        "re_p_max": 0.0,
        "cd_min": float("inf"),
        "cd_max": 0.0,
        "re_p_active_min": float("inf"),
        "re_p_active_max": 0.0,
        "cd_active_min": float("inf"),
        "cd_active_max": 0.0,
        "z_drag": 0.0,
        "exact_delta_pos": 0,
        "exact_delta_zero": 0,
        "exact_delta_neg": 0,
        "exact_a_rel_max": 0.0,
        "exact_b_rel_max": 0.0,
        "exact_c_rel_max": 0.0,
        "exact_dp_grav_min": float("inf"),
        "exact_dp_grav_max": 0.0,
        "exact_dp_ratio_min": float("inf"),
        "exact_dp_ratio_max": 0.0,
    }
    u_g = float(0.5 * (u_g_prev + u_g_next))
    z_drag = (
        _hamel_exact_high_re_drag_constant(float(phi_s))
        if coeff_model == "exact_hamel"
        else _haider_high_re_drag_constant(float(phi_s))
    )
    diag["z_drag"] = float(z_drag)

    for j in range(len(u_p_samples)):
        if m_dot_samples_prev[j] <= 0.0:
            u_p_samples[j] = max(u_p_samples_prev[j], 1e-6)
            m_dot_samples[j] = 0.0
            continue

        z = float(z_lower_abs)
        t_seg = 0.0
        u_p = max(float(u_p_samples_prev[j]), 1e-6)
        survived = False
        branch_mode = "slower" if u_p < u_g else "faster"
        slower_ref_u = float(slower_ref_u_samples[j]) if branch_mode == "slower" else None
        slower_ref_z = float(slower_ref_z_samples[j]) if branch_mode == "slower" else None
        for _ in range(6):
            if z >= float(z_lower_abs + dh) - 1e-10:
                survived = True
                break

            rel = abs(float(u_g) - float(u_p))
            re_p = max(float(rho_g) * rel * max(float(d_p[j]), 1e-9) / max(float(mu_g), 1e-12), 1e-9)
            cd = calc_cd_haider(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
            diag["re_p_min"] = min(float(diag["re_p_min"]), float(re_p))
            diag["re_p_max"] = max(float(diag["re_p_max"]), float(re_p))
            diag["cd_min"] = min(float(diag["cd_min"]), float(cd))
            diag["cd_max"] = max(float(diag["cd_max"]), float(cd))
            if rel >= 1e-3:
                diag["re_p_active_min"] = min(float(diag["re_p_active_min"]), float(re_p))
                diag["re_p_active_max"] = max(float(diag["re_p_active_max"]), float(re_p))
                diag["cd_active_min"] = min(float(diag["cd_active_min"]), float(cd))
                diag["cd_active_max"] = max(float(diag["cd_active_max"]), float(cd))

            a_ex, b_ex, c_ex, _, _, _ = _wirsum_abc_exact_audit(
                u_g=u_g,
                u_p=u_p,
                rho_g=rho_g,
                rho_p=rho_p,
                d_p=float(d_p[j]),
                mu_g=mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if coeff_model == "exact_hamel":
                a, b, c = a_ex, b_ex, c_ex
            elif coeff_model == "stable_mixed_drag_split":
                a, b, c = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch=branch_mode,
                )
            else:
                raise ValueError(f"Unknown analytical_wirsum coeff_model: {coeff_model}")
            delta_exact = 4.0 * a_ex * c_ex - b_ex * b_ex
            if delta_exact > 1e-14:
                diag["exact_delta_pos"] += 1
            elif delta_exact < -1e-14:
                diag["exact_delta_neg"] += 1
            else:
                diag["exact_delta_zero"] += 1
            d_p_grav = _estimate_exact_hamel_critical_diameter(
                u_g=u_g,
                u_p=u_p,
                rho_g=rho_g,
                rho_p=rho_p,
                mu_g=mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if d_p_grav is not None:
                diag["exact_dp_grav_min"] = min(float(diag["exact_dp_grav_min"]), float(d_p_grav))
                diag["exact_dp_grav_max"] = max(float(diag["exact_dp_grav_max"]), float(d_p_grav))
                ratio = float(d_p[j]) / max(float(d_p_grav), 1e-12)
                diag["exact_dp_ratio_min"] = min(float(diag["exact_dp_ratio_min"]), ratio)
                diag["exact_dp_ratio_max"] = max(float(diag["exact_dp_ratio_max"]), ratio)
            diag["exact_a_rel_max"] = max(
                float(diag["exact_a_rel_max"]),
                float(abs(a_ex - a) / max(abs(a), abs(a_ex), 1e-12)),
            )
            diag["exact_b_rel_max"] = max(
                float(diag["exact_b_rel_max"]),
                float(abs(b_ex - b) / max(abs(b), abs(b_ex), 1e-12)),
            )
            diag["exact_c_rel_max"] = max(
                float(diag["exact_c_rel_max"]),
                float(abs(c_ex - c) / max(abs(c), abs(c_ex), 1e-12)),
            )
            _, _, delta = _quadratic_ode_state(u_p, a, b, c, 0.0)
            if delta > 1e-14:
                diag["delta_pos"] += 1
                diag["disc_neg"] += 1
            elif delta < -1e-14:
                diag["delta_neg"] += 1
                diag["disc_pos"] += 1
            else:
                diag["delta_zero"] += 1
                diag["disc_zero"] += 1

            z_upper_abs = float(z_lower_abs + dh)

            t_scale = dh / max(max(abs(u_p), abs(u_g), 1e-3), 1e-6)
            t_max = max(20.0 * t_scale, 1.0)

            def vel_at(dt: float) -> float:
                u_dt, _, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(u_dt)

            def z_at(dt: float) -> float:
                _, dz_dt, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(z + dz_dt)

            t_switch = None
            if branch_mode == "faster":
                t_switch = _find_first_event_time(vel_at, target=u_g, t_max=t_max)
            t_turn = None
            if u_p > 0.0:
                t_turn = _find_first_event_time(vel_at, target=0.0, t_max=t_max)
            t_top = _find_first_event_time(z_at, target=z_upper_abs, t_max=t_max)

            candidates = [t for t in (t_switch, t_turn, t_top) if t is not None and t > 1e-10]
            if not candidates:
                if branch_mode == "slower" and u_p > 0.0:
                    ref_u = float(slower_ref_u if slower_ref_u is not None else u_p)
                    ref_z = float(slower_ref_z if slower_ref_z is not None else z)
                    t_to_current = _hamel_delta_neg_time_between_velocities(
                        u_start=ref_u,
                        u_end=float(u_p),
                        a=a,
                        b=b,
                        c=c,
                    )
                    zero_reachable = _hamel_delta_neg_zero_reachable(a, b, c)
                    roots = _hamel_delta_neg_real_roots(a, b, c)

                    if not zero_reachable:
                        u_lo = 1e-9
                        if roots is not None:
                            u_lo = max(float(roots[0]) + 1e-9, 1e-9)
                        # Hamel Fall 2 with two positive roots corresponds to
                        # a particle that never reaches u_p = 0. In that case
                        # Eq. 3.103 / h_max is not the right discriminant;
                        # the particle remains elutriated and we must recover
                        # the top-boundary crossing from the time domain.
                        if coeff_model == "exact_hamel":
                            t_cap_long = max(120.0, 8.0 * t_max)
                            if roots is not None:
                                u_lim = max(float(roots[0]), 1e-6)
                                t_cap_long = max(t_cap_long, 6.0 * dh / u_lim)
                            t_top_long = _find_first_event_time_with_growth(
                                z_at,
                                target=z_upper_abs,
                                t_init=t_max,
                                t_cap=t_cap_long,
                            )
                            if t_top_long is not None and t_top_long > 1e-10:
                                u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_top_long)
                                t_seg += float(t_top_long)
                                u_p = max(float(u_evt), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=ref_u,
                            h0=ref_z,
                            a=a,
                            b=b,
                            c=c,
                            u_lo=u_lo,
                            u_hi=max(float(u_p), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_total = _hamel_delta_neg_time_between_velocities(
                                u_start=ref_u,
                                u_end=float(u_bound_up),
                                a=a,
                                b=b,
                                c=c,
                            )
                            if (
                                t_total is not None
                                and t_to_current is not None
                                and math.isfinite(t_total)
                                and math.isfinite(t_to_current)
                            ):
                                t_cross = max(float(t_total - t_to_current), 0.0)
                                t_seg += float(t_cross)
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                    h_max = _hamel_delta_neg_max_height(
                        u0=ref_u,
                        h0=ref_z,
                        a=a,
                        b=b,
                        c=c,
                    )
                    if h_max is not None and h_max < float(z_upper_abs) - 1e-10:
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=ref_u,
                            u_end=0.0,
                            a=a,
                            b=b,
                            c=c,
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a,
                            b=b,
                            c=c,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a,
                            b=b,
                            c=c,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and t_to_current is not None
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                            and math.isfinite(t_to_current)
                        ):
                            t_back = max(float(t_to_apex_total - t_to_current), 0.0) + max(float(t_down), 0.0)
                            t_seg += float(t_back)
                            diag["returns"] += 1
                            return_samples[j] = float(m_dot_samples_prev[j])
                            u_p_samples[j] = 0.0
                            m_dot_samples[j] = 0.0
                            break
                    elif h_max is not None and h_max >= float(z_upper_abs) - 1e-10:
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=ref_u,
                            h0=ref_z,
                            a=a,
                            b=b,
                            c=c,
                            u_lo=1e-9,
                            u_hi=max(float(u_p), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_total = _hamel_delta_neg_time_between_velocities(
                                u_start=ref_u,
                                u_end=float(u_bound_up),
                                a=a,
                                b=b,
                                c=c,
                            )
                            if (
                                t_total is not None
                                and t_to_current is not None
                                and math.isfinite(t_total)
                                and math.isfinite(t_to_current)
                            ):
                                t_cross = max(float(t_total - t_to_current), 0.0)
                                t_seg += float(t_cross)
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                survived = True
                                break
                diag["fallbacks"] += 1
                diag["fallback_no_event"] += 1
                return_samples[j] = float(m_dot_samples_prev[j])
                u_p_samples[j] = 0.0
                m_dot_samples[j] = 0.0
                break

            t_evt = min(candidates)
            u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_evt)
            z = float(z + dz_evt)
            t_seg += float(t_evt)

            if t_top is not None and abs(t_evt - t_top) <= 1e-8 * max(1.0, t_top):
                u_p = max(float(u_evt), 1e-6)
                z = float(z_upper_abs)
                survived = True
                break

            if t_switch is not None and abs(t_evt - t_switch) <= 1e-8 * max(1.0, t_switch):
                diag["sign_switch"] += 1
                # Preserve the analytical event state continuously at the
                # branch switch instead of hard-resetting the particle speed.
                u_p = max(float(u_g), 1e-9)
                slower_ref_u = float(u_g)
                slower_ref_z = float(z)
                slower_ref_u_samples[j] = float(u_g)
                slower_ref_z_samples[j] = float(z)
                branch_mode = "slower"
                continue

            if t_turn is not None and abs(t_evt - t_turn) <= 1e-8 * max(1.0, t_turn):
                # Hamel Fall 2 uses the signed particle velocity directly.
                # After the apex, solve the lower-boundary return using
                # Eq. 3.103 (height as a function of signed velocity) and
                # the exact Delta<0 time primitive instead of restarting
                # a pseudo-time march from u=0 and relying on fallback.
                u_p = min(float(u_evt), 0.0)
                a2, b2, c2 = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch="slower",
                ) if coeff_model == "stable_mixed_drag_split" else _wirsum_abc_exact_audit(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=mu_g,
                    phi_s=phi_s,
                    branch="slower",
                )[:3]
                apex_z = float(z)
                u_bound = _hamel_delta_neg_boundary_velocity(
                    h_target=float(z_lower_abs),
                    u0=0.0,
                    h0=apex_z,
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if u_bound is None:
                    diag["fallbacks"] += 1
                    diag["fallback_boundary_velocity"] += 1
                    return_samples[j] = float(m_dot_samples_prev[j])
                    u_p_samples[j] = 0.0
                    m_dot_samples[j] = 0.0
                    break
                t_back = _hamel_delta_neg_time_between_velocities(
                    u_start=0.0,
                    u_end=float(u_bound),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if (
                    t_back is None
                    or not math.isfinite(t_back)
                ):
                    diag["fallbacks"] += 1
                    diag["fallback_delta_neg_time"] += 1
                    return_samples[j] = float(m_dot_samples_prev[j])
                    u_p_samples[j] = 0.0
                    m_dot_samples[j] = 0.0
                    break
                t_seg += max(float(t_back), 0.0)
                diag["returns"] += 1
                return_samples[j] = float(m_dot_samples_prev[j])
                u_p_samples[j] = 0.0
                m_dot_samples[j] = 0.0
                break

        hold_samples[j] = float(m_dot_samples_prev[j]) * max(t_seg, 0.0)
        if survived:
            u_p_samples[j] = max(u_p, 1e-6)
            if branch_mode == "slower" and slower_ref_u is not None and slower_ref_z is not None:
                # Crossing into the next cell changes the local background-gas
                # state and therefore the frozen Fall-2 coefficients. Hamel's
                # boundary-to-boundary construction carries the *boundary*
                # state forward as the next cell's initial condition.
                slower_ref_u_samples[j] = max(float(u_p), 1e-6)
                slower_ref_z_samples[j] = float(z)
            else:
                slower_ref_u_samples[j] = max(float(u_p), 1e-6)
                slower_ref_z_samples[j] = 0.0
        else:
            return_samples[j] = float(m_dot_samples_prev[j])
            u_p_samples[j] = 0.0
            m_dot_samples[j] = 0.0
            slower_ref_u_samples[j] = max(float(u_p_samples_prev[j]), 1e-6)
            slower_ref_z_samples[j] = 0.0

    return (
        u_p_samples,
        m_dot_samples,
        hold_samples,
        return_samples,
        slower_ref_u_samples,
        slower_ref_z_samples,
        diag,
    )


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
    solid_d_p: float | None = None,
    D_g: float | None = None,
    char_conversion: float = 0.0,
    r4_scale: float = 1.0,
    enabled_reactions: tuple[str, ...] | None = None,
) -> tuple[np.ndarray, dict[str, float]]:
    enabled = set(enabled_reactions or ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R10", "R11", "R12"))
    y, c = _gas_concentrations(N, T, P)
    idx = GAS_SPECIES_INDEX
    c_tar = max(float(c[idx["TAR1"]] + c[idx["TAR2"]]), 0.0)

    r5 = rate_R5_bubble(T, c[idx["CO"]], c[idx["O2"]], P, y, idx) * V_seg if "R5" in enabled else 0.0
    r6 = rate_R6(T, c[idx["CH4"]], c[idx["O2"]]) * V_seg if "R6" in enabled else 0.0
    r12 = rate_R12(T, c[idx["H2"]], c[idx["O2"]], P, y, idx) * V_seg if "R12" in enabled else 0.0
    r7 = rate_R7(T, c[idx["CH4"]], c[idx["H2O"]], c[idx["CO"]], c[idx["H2"]]) * V_seg if "R7" in enabled else 0.0
    r8 = rate_R8(T, P, float(y[idx["CO"]]), float(y[idx["H2O"]]), float(y[idx["CO2"]]), float(y[idx["H2"]])) * V_seg if "R8" in enabled else 0.0
    r10 = rate_R10(T, c_tar, c[idx["O2"]], P, fuel_type) * V_seg if ("R10" in enabled and c_tar > 0.0) else 0.0
    r11b = rate_R11_bubble(T, c_tar) * V_seg if ("R11" in enabled and c_tar > 0.0) else 0.0
    r11d = (
        rate_R11_suspension(T, c_tar, rho_cat) * V_seg
        if ("R11" in enabled and c_tar > 0.0 and rho_cat > 0.0)
        else 0.0
    )
    total_area = max(float(char_area_total), 0.0)
    r1 = r2 = r3 = r4 = 0.0
    alpha = 1.0
    if total_area > 0.0 and solid_d_p is not None and D_g is not None:
        d_core = d_core_from_spm_char_conversion(float(char_conversion), float(solid_d_p))
        if "R1" in enabled:
            r1, alpha = rate_R1(T, c[idx["O2"]], float(solid_d_p), float(D_g), d_core, fuel_type)
        if "R2" in enabled:
            r2 = rate_R2(T, c[idx["H2O"]], float(solid_d_p), float(D_g), d_core)
        if "R3" in enabled:
            r3 = rate_R3(T, c[idx["H2"]], float(solid_d_p), float(D_g), d_core)
        if "R4" in enabled:
            p_co2 = max(float(c[idx["CO2"]] * Rg * T), 0.0)
            p_co = max(float(c[idx["CO"]] * Rg * T), 0.0)
            r4 = float(r4_scale) * rate_R4_effective(T, p_co2, p_co, float(solid_d_p), float(D_g), d_core)

    sto10 = get_lumped_tar_stoichiometry("R10", fuel_type)
    sto11 = get_lumped_tar_stoichiometry("R11", fuel_type)

    o2_demand = r5 + 1.5 * r6 + 0.5 * r12 + max(0.0, -float(sto10.get("O2", 0.0))) * r10 + alpha * r1 * total_area
    h2o_demand = max(r7, 0.0) + max(r8, 0.0) + max(0.0, -float(sto11.get("H2O", 0.0))) * (r11b + r11d) + r2 * total_area
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

    e1 = _bound_extent_forward(r1 * total_area, dt, N, {"O2": max(alpha, 1e-12)}) if r1 > 0.0 and total_area > 0.0 else 0.0
    e2 = _bound_extent_forward(r2 * total_area, dt, N, {"H2O": 1.0}) if r2 > 0.0 and total_area > 0.0 else 0.0
    e3 = _bound_extent_forward(r3 * total_area, dt, N, {"H2": 2.0}) if r3 > 0.0 and total_area > 0.0 else 0.0
    e4 = _bound_extent_forward(r4 * total_area, dt, N, {"CO2": 1.0}) if r4 > 0.0 and total_area > 0.0 else 0.0
    e5 = _bound_extent_forward(r5, dt, N, {"CO": 2.0, "O2": 1.0})
    e6 = _bound_extent_forward(r6, dt, N, {"CH4": 1.0, "O2": 1.5})
    e12 = _bound_extent_forward(r12, dt, N, {"H2": 1.0, "O2": 0.5})
    e7 = _bound_extent_signed(r7, dt, N, {"CH4": 1.0, "H2O": 1.0}, {"CO": 1.0, "H2": 3.0})
    e8 = _bound_extent_signed(r8, dt, N, {"CO": 1.0, "H2O": 1.0}, {"CO2": 1.0, "H2": 1.0})
    e10 = _bound_extent_forward(r10, dt, N, {"TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0, "O2": max(0.0, -float(sto10.get("O2", 0.0)))}) if r10 > 0.0 else 0.0
    e11b = _bound_extent_forward(r11b, dt, N, {"TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0}) if r11b > 0.0 else 0.0
    e11d = _bound_extent_forward(r11d, dt, N, {"TAR1": 1.0 if float(N[idx["TAR1"]]) >= float(N[idx["TAR2"]]) else 0.0}) if r11d > 0.0 else 0.0
    e11 = e11b + e11d

    out = np.array(N, dtype=np.float64, copy=True)

    out[idx["CO"]] += -2.0 * e5 + e6 + e7 - e8
    out[idx["O2"]] += -e5 - 1.5 * e6 - 0.5 * e12
    out[idx["CO2"]] += 2.0 * e5 + e8
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
) -> list[FreeboardBackground]:
    """Build a segment-wise frozen gas background profile for the whole freeboard."""
    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    dh = H_freeboard / max(int(n_cells), 1)
    seg_loss = 1.0 - math.pow(max(1.0 - heat_loss_frac, 0.0), 1.0 / max(n_cells, 1))
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
        substep_loss = 1.0 - math.pow(max(1.0 - seg_loss, 0.0), 1.0 / 12.0)
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
            )
            T = _solve_T_from_enthalpy(N, H_before * (1.0 - substep_loss), T)

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


def _advance_particle_samples_global_projected_hamel(
    *,
    backgrounds: list[FreeboardBackground],
    rho_p: float,
    d_p: np.ndarray,
    phi_s: float,
    u_p_samples_prev: np.ndarray,
    m_dot_samples_prev: np.ndarray,
    char_frac_samples: np.ndarray,
    char_conversion_samples: np.ndarray,
    class_idx_samples: np.ndarray,
    coeff_model: str,
) -> dict[str, object]:
    """Track each particle sample continuously across the whole freeboard and project back to cells."""
    n_cells = len(backgrounds)
    n_classes = int(np.max(class_idx_samples)) + 1 if len(class_idx_samples) else 0
    hold_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    upflow_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    downflow_samples_by_seg = [np.zeros_like(m_dot_samples_prev, dtype=np.float64) for _ in range(n_cells)]
    exit_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    exit_u_samples = np.zeros_like(u_p_samples_prev, dtype=np.float64)
    total_return_samples = np.zeros_like(m_dot_samples_prev, dtype=np.float64)
    diag = {
        "delta_pos": 0,
        "delta_zero": 0,
        "delta_neg": 0,
        "disc_pos": 0,
        "disc_zero": 0,
        "disc_neg": 0,
        "sign_switch": 0,
        "fallbacks": 0,
        "returns": 0,
        "fallback_no_event": 0,
        "fallback_boundary_velocity": 0,
        "fallback_delta_neg_time": 0,
        "projection_crossings": 0,
    }
    coeff_diag = {
        "re_p_min": float("inf"),
        "re_p_max": 0.0,
        "cd_min": float("inf"),
        "cd_max": 0.0,
        "re_p_active_min": float("inf"),
        "re_p_active_max": 0.0,
        "cd_active_min": float("inf"),
        "cd_active_max": 0.0,
        "z_drag": 0.0,
        "exact_delta_pos": 0,
        "exact_delta_zero": 0,
        "exact_delta_neg": 0,
        "exact_a_rel_max": 0.0,
        "exact_b_rel_max": 0.0,
        "exact_c_rel_max": 0.0,
        "exact_dp_grav_min": float("inf"),
        "exact_dp_grav_max": 0.0,
        "exact_dp_ratio_min": float("inf"),
        "exact_dp_ratio_max": 0.0,
    }

    for j in range(len(u_p_samples_prev)):
        if m_dot_samples_prev[j] <= 0.0:
            continue
        cell_idx = 0
        z = 0.0
        u_p = max(float(u_p_samples_prev[j]), 1e-6)
        ref_u = float(u_p)
        ref_z = 0.0
        steps = 0
        while 0 <= cell_idx < n_cells and steps < max(12 * n_cells, 32):
            steps += 1
            bg = backgrounds[cell_idx]
            z = min(max(float(z), float(bg.z_lower_m)), float(bg.z_upper_m))
            u_g = float(bg.u_gb)
            branch_mode = "slower" if u_p < u_g else "faster"
            rel = abs(float(u_g) - float(u_p))
            re_p = max(float(bg.rho_g) * rel * max(float(d_p[j]), 1e-9) / max(float(bg.mu_g), 1e-12), 1e-9)
            cd = calc_cd_haider(re_p, float(np.clip(phi_s, 1e-3, 1.0)))
            coeff_diag["re_p_min"] = min(float(coeff_diag["re_p_min"]), float(re_p))
            coeff_diag["re_p_max"] = max(float(coeff_diag["re_p_max"]), float(re_p))
            coeff_diag["cd_min"] = min(float(coeff_diag["cd_min"]), float(cd))
            coeff_diag["cd_max"] = max(float(coeff_diag["cd_max"]), float(cd))
            if rel >= 1e-3:
                coeff_diag["re_p_active_min"] = min(float(coeff_diag["re_p_active_min"]), float(re_p))
                coeff_diag["re_p_active_max"] = max(float(coeff_diag["re_p_active_max"]), float(re_p))
                coeff_diag["cd_active_min"] = min(float(coeff_diag["cd_active_min"]), float(cd))
                coeff_diag["cd_active_max"] = max(float(coeff_diag["cd_active_max"]), float(cd))

            a_ex, b_ex, c_ex, _, _, z_drag = _wirsum_abc_exact_audit(
                u_g=u_g,
                u_p=u_p,
                rho_g=bg.rho_g,
                rho_p=rho_p,
                d_p=float(d_p[j]),
                mu_g=bg.mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            coeff_diag["z_drag"] = float(max(float(coeff_diag["z_drag"]), float(z_drag)))
            if coeff_model == "exact_hamel":
                a, b, c = a_ex, b_ex, c_ex
            else:
                a, b, c = _wirsum_abc(
                    u_g=u_g,
                    u_p=u_p,
                    rho_g=bg.rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=bg.mu_g,
                    phi_s=phi_s,
                    branch=branch_mode,
                )
            delta_exact = 4.0 * a_ex * c_ex - b_ex * b_ex
            if delta_exact > 1e-14:
                coeff_diag["exact_delta_pos"] += 1
            elif delta_exact < -1e-14:
                coeff_diag["exact_delta_neg"] += 1
            else:
                coeff_diag["exact_delta_zero"] += 1
            d_p_grav = _estimate_exact_hamel_critical_diameter(
                u_g=u_g,
                u_p=u_p,
                rho_g=bg.rho_g,
                rho_p=rho_p,
                mu_g=bg.mu_g,
                phi_s=phi_s,
                branch=branch_mode,
            )
            if d_p_grav is not None:
                coeff_diag["exact_dp_grav_min"] = min(float(coeff_diag["exact_dp_grav_min"]), float(d_p_grav))
                coeff_diag["exact_dp_grav_max"] = max(float(coeff_diag["exact_dp_grav_max"]), float(d_p_grav))
                ratio = float(d_p[j]) / max(float(d_p_grav), 1e-12)
                coeff_diag["exact_dp_ratio_min"] = min(float(coeff_diag["exact_dp_ratio_min"]), ratio)
                coeff_diag["exact_dp_ratio_max"] = max(float(coeff_diag["exact_dp_ratio_max"]), ratio)

            _, _, delta = _quadratic_ode_state(u_p, a, b, c, 0.0)
            if delta > 1e-14:
                diag["delta_pos"] += 1
                diag["disc_neg"] += 1
            elif delta < -1e-14:
                diag["delta_neg"] += 1
                diag["disc_pos"] += 1
            else:
                diag["delta_zero"] += 1
                diag["disc_zero"] += 1

            z_lower_abs = float(bg.z_lower_m)
            z_upper_abs = float(bg.z_upper_m)
            t_scale = (z_upper_abs - z) / max(max(abs(u_p), abs(u_g), 1e-3), 1e-6)
            t_max = max(20.0 * max(float(t_scale), 1e-6), 1.0)

            def vel_at(dt: float) -> float:
                u_dt, _, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(u_dt)

            def z_at(dt: float) -> float:
                _, dz_dt, _ = _quadratic_ode_state(u_p, a, b, c, dt)
                return float(z + dz_dt)

            zero_reachable = _hamel_delta_neg_zero_reachable(a, b, c) if branch_mode == "slower" and u_p > 0.0 else True
            t_switch = _find_first_event_time(vel_at, target=u_g, t_max=t_max) if branch_mode == "faster" else None
            t_turn = _find_first_event_time(vel_at, target=0.0, t_max=t_max) if u_p > 0.0 and zero_reachable else None
            t_top = _find_first_event_time(z_at, target=z_upper_abs, t_max=t_max)
            candidates = [t for t in (t_switch, t_turn, t_top) if t is not None and t > 1e-10]

            if not candidates and branch_mode == "slower" and u_p > 0.0:
                roots = _hamel_delta_neg_real_roots(a, b, c)
                if not zero_reachable:
                    u_lim = max(float(roots[0]), 1e-6) if roots is not None else 1e-6
                    t_cap = max(120.0, 6.0 * max(z_upper_abs - z, 1e-9) / u_lim)
                    t_top = _find_first_event_time_with_growth(
                        z_at,
                        target=z_upper_abs,
                        t_init=t_max,
                        t_cap=t_cap,
                    )
                    if t_top is not None and t_top > 1e-10:
                        candidates = [float(t_top)]

            if not candidates:
                if branch_mode == "slower" and u_p > 0.0:
                    h_max = _hamel_delta_neg_max_height(
                        u0=float(ref_u),
                        h0=float(ref_z),
                        a=a,
                        b=b,
                        c=c,
                    )
                    if h_max is not None and h_max < z_upper_abs - 1e-10:
                        t_to_current = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=float(u_p),
                            a=a,
                            b=b,
                            c=c,
                        )
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=0.0,
                            a=a,
                            b=b,
                            c=c,
                        )
                        a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                            u_g=u_g,
                            u_p=0.0,
                            rho_g=bg.rho_g,
                            rho_p=rho_p,
                            d_p=float(d_p[j]),
                            mu_g=bg.mu_g,
                            phi_s=phi_s,
                            branch="slower",
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_current is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and math.isfinite(t_to_current)
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                        ):
                            t_cell = max(float(t_to_apex_total - t_to_current), 0.0) + max(float(t_down), 0.0)
                            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(t_cell, 0.0)
                            downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                            total_return_samples[j] += float(m_dot_samples_prev[j])
                            diag["returns"] += 1
                            break
                diag["fallbacks"] += 1
                diag["fallback_no_event"] += 1
                break

            t_evt = min(candidates)
            u_evt, dz_evt, _ = _quadratic_ode_state(u_p, a, b, c, t_evt)
            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_evt), 0.0)
            z = float(z + dz_evt)

            if t_top is not None and abs(t_evt - t_top) <= 1e-8 * max(1.0, t_top):
                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                u_p = max(float(u_evt), 1e-9)
                z = float(z_upper_abs)
                if branch_mode == "slower":
                    ref_u = max(float(u_p), 1e-9)
                    ref_z = float(z)
                diag["projection_crossings"] += 1
                cell_idx += 1
                if cell_idx >= n_cells:
                    exit_samples[j] = float(m_dot_samples_prev[j])
                    exit_u_samples[j] = float(u_p)
                    break
                continue

            if t_switch is not None and abs(t_evt - t_switch) <= 1e-8 * max(1.0, t_switch):
                diag["sign_switch"] += 1
                ref_u = float(u_g)
                ref_z = float(z)
                # In the thesis-aligned exact path, ``u_p = u_g`` is the
                # analytical branch handoff point. Resolve the remaining part
                # of the current cell with Fall 2 boundary formulas instead of
                # re-entering the time-domain marcher from a numerically
                # nudged velocity. This preserves per-sample boundary states
                # and avoids collapsing multiple trajectories onto a single
                # lower-root attractor inside the same cell.
                if coeff_model == "exact_hamel":
                    a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                        u_g=u_g,
                        u_p=float(u_g),
                        rho_g=bg.rho_g,
                        rho_p=rho_p,
                        d_p=float(d_p[j]),
                        mu_g=bg.mu_g,
                        phi_s=phi_s,
                        branch="slower",
                    )
                    zero_reachable_2 = _hamel_delta_neg_zero_reachable(a2, b2, c2)
                    roots_2 = _hamel_delta_neg_real_roots(a2, b2, c2)
                    if not zero_reachable_2:
                        u_lo = 1e-9
                        if roots_2 is not None:
                            u_lo = max(float(roots_2[0]) + 1e-9, 1e-9)
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=float(ref_u),
                            h0=float(ref_z),
                            a=a2,
                            b=b2,
                            c=c2,
                            u_lo=u_lo,
                            u_hi=max(float(ref_u), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_cross = _hamel_delta_neg_time_between_velocities(
                                u_start=float(ref_u),
                                u_end=float(u_bound_up),
                                a=a2,
                                b=b2,
                                c=c2,
                            )
                            if t_cross is not None and math.isfinite(t_cross):
                                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cross), 0.0)
                                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                ref_u = float(u_p)
                                ref_z = float(z)
                                diag["projection_crossings"] += 1
                                cell_idx += 1
                                if cell_idx >= n_cells:
                                    exit_samples[j] = float(m_dot_samples_prev[j])
                                    exit_u_samples[j] = float(u_p)
                                    break
                                continue
                    h_max = _hamel_delta_neg_max_height(
                        u0=float(ref_u),
                        h0=float(ref_z),
                        a=a2,
                        b=b2,
                        c=c2,
                    )
                    if h_max is not None and h_max < float(z_upper_abs) - 1e-10:
                        t_to_apex_total = _hamel_delta_neg_time_between_velocities(
                            u_start=float(ref_u),
                            u_end=0.0,
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        u_bound_down = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_lower_abs),
                            u0=0.0,
                            h0=float(h_max),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        t_down = _hamel_delta_neg_time_between_velocities(
                            u_start=0.0,
                            u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                            a=a2,
                            b=b2,
                            c=c2,
                        )
                        if (
                            u_bound_down is not None
                            and t_to_apex_total is not None
                            and t_down is not None
                            and math.isfinite(t_to_apex_total)
                            and math.isfinite(t_down)
                        ):
                            t_cell = max(float(t_to_apex_total), 0.0) + max(float(t_down), 0.0)
                            hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cell), 0.0)
                            downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                            total_return_samples[j] += float(m_dot_samples_prev[j])
                            diag["returns"] += 1
                            break
                    elif h_max is not None and h_max >= float(z_upper_abs) - 1e-10:
                        u_bound_up = _hamel_delta_neg_boundary_velocity(
                            h_target=float(z_upper_abs),
                            u0=float(ref_u),
                            h0=float(ref_z),
                            a=a2,
                            b=b2,
                            c=c2,
                            u_lo=1e-9,
                            u_hi=max(float(ref_u), 1e-9),
                        )
                        if u_bound_up is not None:
                            t_cross = _hamel_delta_neg_time_between_velocities(
                                u_start=float(ref_u),
                                u_end=float(u_bound_up),
                                a=a2,
                                b=b2,
                                c=c2,
                            )
                            if t_cross is not None and math.isfinite(t_cross):
                                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_cross), 0.0)
                                upflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                                u_p = max(float(u_bound_up), 1e-9)
                                z = float(z_upper_abs)
                                ref_u = float(u_p)
                                ref_z = float(z)
                                diag["projection_crossings"] += 1
                                cell_idx += 1
                                if cell_idx >= n_cells:
                                    exit_samples[j] = float(m_dot_samples_prev[j])
                                    exit_u_samples[j] = float(u_p)
                                    break
                                continue
                if float(u_g) > 1e-9:
                    u_p = max(float(np.nextafter(float(u_g), -math.inf)), 1e-9)
                else:
                    u_p = 1e-9
                continue

            if t_turn is not None and abs(t_evt - t_turn) <= 1e-8 * max(1.0, t_turn):
                a2, b2, c2, _, _, _ = _wirsum_abc_exact_audit(
                    u_g=u_g,
                    u_p=0.0,
                    rho_g=bg.rho_g,
                    rho_p=rho_p,
                    d_p=float(d_p[j]),
                    mu_g=bg.mu_g,
                    phi_s=phi_s,
                    branch="slower",
                )
                u_bound_down = _hamel_delta_neg_boundary_velocity(
                    h_target=float(z_lower_abs),
                    u0=0.0,
                    h0=float(z),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                t_back = _hamel_delta_neg_time_between_velocities(
                    u_start=0.0,
                    u_end=float(u_bound_down) if u_bound_down is not None else float("nan"),
                    a=a2,
                    b=b2,
                    c=c2,
                )
                if u_bound_down is None or t_back is None or not math.isfinite(t_back):
                    diag["fallbacks"] += 1
                    diag["fallback_delta_neg_time"] += 1
                    break
                hold_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j]) * max(float(t_back), 0.0)
                downflow_samples_by_seg[cell_idx][j] += float(m_dot_samples_prev[j])
                total_return_samples[j] += float(m_dot_samples_prev[j])
                diag["returns"] += 1
                break

    class_transport_by_seg: list[dict[str, np.ndarray]] = []
    hold_time_mean_by_seg: list[float] = []
    hold_time_max_by_seg: list[float] = []
    for i in range(n_cells):
        if n_classes <= 0:
            class_transport_by_seg.append(_empty_class_transport(0))
            hold_time_mean_by_seg.append(0.0)
            hold_time_max_by_seg.append(0.0)
            continue
        positive_mask = np.maximum(m_dot_samples_prev, 0.0) > 1e-12
        hold_time_samples = np.divide(
            hold_samples_by_seg[i],
            np.maximum(m_dot_samples_prev, 1e-12),
            out=np.zeros_like(hold_samples_by_seg[i]),
            where=positive_mask,
        )
        active_times = hold_time_samples[(hold_time_samples > 1e-12) & positive_mask]
        hold_time_mean_by_seg.append(float(np.mean(active_times)) if active_times.size else 0.0)
        hold_time_max_by_seg.append(float(np.max(active_times)) if active_times.size else 0.0)
        class_transport_by_seg.append(
            _aggregate_class_transport(
                n_classes=n_classes,
                class_idx_samples=class_idx_samples,
                hold_samples=hold_samples_by_seg[i],
                upflow_samples=upflow_samples_by_seg[i],
                downflow_samples=downflow_samples_by_seg[i],
                char_frac_samples=char_frac_samples,
                char_conversion_samples=char_conversion_samples,
            )
        )
    return {
        "class_transport_by_seg": class_transport_by_seg,
        "diag": diag,
        "coeff_diag": coeff_diag,
        "hold_time_mean_by_seg": hold_time_mean_by_seg,
        "hold_time_max_by_seg": hold_time_max_by_seg,
        "exit_samples": exit_samples,
        "exit_u_samples": exit_u_samples,
        "return_samples": total_return_samples,
    }


def _simulate_freeboard_global_projected_hamel(
    *,
    N_in: np.ndarray,
    T_in: float,
    P: float,
    D_bed: float,
    H_freeboard: float,
    n_cells: int,
    u_b_bed_top: float,
    d_b_bed_top: float,
    eps_b_bed_top: float,
    eps_d_void_bed_top: float,
    rho_solid_bed_top: float,
    d_p_classes_bed_top: np.ndarray,
    m_char_classes_bed_top: np.ndarray,
    m_ash_classes_bed_top: np.ndarray,
    char_conversion_bed_top: float,
    class_aggregate_indices_bed_top: np.ndarray | None,
    n_output_size_classes: int | None,
    fuel_type: str,
    heat_loss_frac: float,
    trajectory_coeff_model: str,
    u_gb_scale: float,
    beta_a_scale: float,
    velocity_sigma: float,
    velocity_bins: int,
    phi_s_bed_top: float,
    z_base_m: float,
    total_height_m: float | None,
    secondary_injection_xi: float | None,
    secondary_O2_mol_s: float,
    secondary_N2_mol_s: float,
    secondary_H2O_mol_s: float,
    secondary_T_K: float,
    secondary_injection_mode: str,
    enabled_reactions: tuple[str, ...] | None,
    inventory_char_sink_enabled: bool,
) -> dict:
    A = math.pi * (D_bed**2) / 4.0
    dh = H_freeboard / n_cells
    beta_a = calc_beta_a(max(d_b_bed_top, 1e-6)) * float(max(beta_a_scale, 0.0))
    u_gb0 = calc_u_gb(max(u_b_bed_top, 1e-6), scale=float(max(u_gb_scale, 1e-6)))
    u_p0_mean = 1.53 * max(u_b_bed_top, 1e-6)
    u_p_seed, velocity_weights = _build_velocity_samples(
        u_p0_mean,
        sigma=float(max(velocity_sigma, 0.0)),
        n_bins=int(max(velocity_bins, 1)),
    )
    d_p_class, char_frac_class, m_dot_class0, m_dot_eject0 = _build_size_class_bundles(
        d_p_classes=np.asarray(d_p_classes_bed_top, dtype=np.float64),
        m_char_classes=np.asarray(m_char_classes_bed_top, dtype=np.float64),
        m_ash_classes=np.asarray(m_ash_classes_bed_top, dtype=np.float64),
        rho_p=float(rho_solid_bed_top),
        eps_b=float(eps_b_bed_top),
        eps_d_void=float(eps_d_void_bed_top),
        u_b=float(u_b_bed_top),
        d_b=float(d_b_bed_top),
        area=A,
    )
    aggregate_idx, n_transport_classes, d_p_transport_class = _resolve_transport_class_mapping(
        d_p_class=d_p_class,
        m_dot_class=m_dot_class0,
        class_aggregate_indices=class_aggregate_indices_bed_top,
        n_output_size_classes=n_output_size_classes,
    )
    m_dot_samples0 = np.repeat(m_dot_class0, len(velocity_weights)) * np.tile(velocity_weights, len(m_dot_class0))
    char_frac_samples = np.repeat(char_frac_class, len(velocity_weights))
    char_conversion_samples = np.full_like(m_dot_samples0, float(np.clip(char_conversion_bed_top, 0.0, 1.0 - 1e-9)), dtype=np.float64)
    class_idx_samples = np.repeat(aggregate_idx, len(velocity_weights))
    d_p_samples = np.repeat(d_p_class, len(velocity_weights))
    u_p_samples0 = np.tile(u_p_seed, len(m_dot_class0))
    secondary_N = np.zeros(N_GAS, dtype=np.float64)
    secondary_N[GAS_SPECIES_INDEX["O2"]] = float(max(secondary_O2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["N2"]] = float(max(secondary_N2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["H2O"]] = float(max(secondary_H2O_mol_s, 0.0))

    backgrounds = _build_freeboard_background_profile(
        N_in=N_in,
        T_in=T_in,
        P=P,
        A=A,
        H_freeboard=H_freeboard,
        n_cells=n_cells,
        u_gb0=u_gb0,
        beta_a=beta_a,
        heat_loss_frac=heat_loss_frac,
        fuel_type=fuel_type,
        z_base_m=z_base_m,
        total_height_m=total_height_m,
        secondary_injection_xi=secondary_injection_xi,
        secondary_N=secondary_N,
        secondary_T_K=secondary_T_K,
        secondary_injection_mode=secondary_injection_mode,
        enabled_reactions=enabled_reactions,
    )
    traj = _advance_particle_samples_global_projected_hamel(
        backgrounds=backgrounds,
        rho_p=max(float(rho_solid_bed_top), 1e-6),
        d_p=d_p_samples,
        phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
        u_p_samples_prev=u_p_samples0,
        m_dot_samples_prev=m_dot_samples0,
        char_frac_samples=char_frac_samples,
        char_conversion_samples=char_conversion_samples,
        class_idx_samples=class_idx_samples,
        coeff_model=str(trajectory_coeff_model),
    )
    rho_cat_profile = _estimate_rho_cat_profile_from_projected_transport(
        class_transport_by_seg=traj["class_transport_by_seg"],
        backgrounds=backgrounds,
    )
    if np.any(rho_cat_profile > 0.0):
        backgrounds = _build_freeboard_background_profile(
            N_in=N_in,
            T_in=T_in,
            P=P,
            A=A,
            H_freeboard=H_freeboard,
            n_cells=n_cells,
            u_gb0=u_gb0,
            beta_a=beta_a,
            heat_loss_frac=heat_loss_frac,
            fuel_type=fuel_type,
            z_base_m=z_base_m,
            total_height_m=total_height_m,
            secondary_injection_xi=secondary_injection_xi,
            secondary_N=secondary_N,
            secondary_T_K=secondary_T_K,
            secondary_injection_mode=secondary_injection_mode,
            enabled_reactions=enabled_reactions,
            rho_cat_profile=rho_cat_profile,
        )
        traj = _advance_particle_samples_global_projected_hamel(
            backgrounds=backgrounds,
            rho_p=max(float(rho_solid_bed_top), 1e-6),
            d_p=d_p_samples,
            phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
            u_p_samples_prev=u_p_samples0,
            m_dot_samples_prev=m_dot_samples0,
            char_frac_samples=char_frac_samples,
            char_conversion_samples=char_conversion_samples,
            class_idx_samples=class_idx_samples,
            coeff_model=str(trajectory_coeff_model),
        )

    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    states: list[FreeboardState] = []
    total_return_char = 0.0
    total_return_ash = 0.0
    secondary_injection_applied = False
    secondary_injection_segment: int | None = None
    h_total = None if total_height_m is None else float(max(total_height_m, 1e-12))
    xi_inj = None if secondary_injection_xi is None else float(np.clip(secondary_injection_xi, 0.0, 1.0))
    seg_loss = 1.0 - math.pow(max(1.0 - heat_loss_frac, 0.0), 1.0 / max(n_cells, 1))
    secondary_observation = {
        "z_m": [],
        "xi_reactor": [],
        "slice_global_index": [],
        "segment_index": [],
        "slice_index_within_segment": [],
        "pre_mix_T": [],
        "post_mix_pre_rxn_T": [],
        "post_rxn_T": [],
        "pre_mix_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_mix_pre_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
    }

    for i, bg in enumerate(backgrounds):
        class_transport = traj["class_transport_by_seg"][i]
        hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
        hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)
        x_before = np.asarray(class_transport["char_conversion_classes"], dtype=np.float64).copy()
        up_char = np.asarray(class_transport["m_dot_auf_char_classes"], dtype=np.float64)
        up_ash = np.asarray(class_transport["m_dot_auf_ash_classes"], dtype=np.float64)
        down_char = np.asarray(class_transport["m_dot_ab_char_classes"], dtype=np.float64)
        down_ash = np.asarray(class_transport["m_dot_ab_ash_classes"], dtype=np.float64)
        char_area_total, solid_d_p_eff = _char_area_and_representative_diameter(
            hold_char,
            d_p_transport_class,
            max(float(rho_solid_bed_top), 1e-6),
        )
        char_conversion_seg = _representative_char_conversion(class_transport)
        m_hold_seg = float(np.sum(hold_char + hold_ash))
        return_char_seg = float(np.sum(np.maximum(down_char, 0.0)))
        return_ash_seg = float(np.sum(np.maximum(down_ash, 0.0)))
        total_return_char += return_char_seg
        total_return_ash += return_ash_seg
        m_dot_char = float(np.sum(np.maximum(up_char, 0.0)))
        m_dot_ash = float(np.sum(np.maximum(up_ash, 0.0)))
        m_dot_s_avg = max(m_dot_char + m_dot_ash, 1e-12)
        ash_hold_frac = m_dot_ash / m_dot_s_avg
        rho_cat = (m_hold_seg * ash_hold_frac) / max(bg.V_seg, 1e-12)
        u_t_samples = np.array(
            [
                calc_terminal_velocity_haider(
                    rho_g=bg.rho_g,
                    rho_p=max(float(rho_solid_bed_top), 1e-6),
                    d_p=max(float(dp), 1e-9),
                    phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                    mu_g=bg.mu_g,
                    g_acc=g,
                )
                for dp in d_p_samples
            ],
            dtype=np.float64,
        )
        flow_avg = 0.5 * (np.maximum(up_char + up_ash, 0.0) + np.maximum(down_char + down_ash, 0.0))
        u_p_mean = float(bg.u_gb)
        u_t_mean = float(np.mean(u_t_samples)) if len(u_t_samples) else 0.0
        carry_ratio = float(bg.u_gb / max(u_t_mean, 1e-12))
        xi_seg = (
            float((float(z_base_m) + bg.z_center_m) / h_total)
            if h_total is not None
            else float(bg.z_center_m / max(H_freeboard, 1e-12))
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
            secondary_injection_segment = i

        pre_mix_T = float(T)
        pre_mix_wet = _wet_gas_profile(N)
        post_mix_pre_rxn_T = float(T)
        post_mix_pre_rxn_wet = _wet_gas_profile(N)
        dt_slice = bg.tau / 12.0
        substep_loss = 1.0 - math.pow(max(1.0 - seg_loss, 0.0), 1.0 / 12.0)
        diag_acc: dict[str, float] = {}
        for k in range(12):
            H_before = calc_gas_enthalpy_flow(N, T, h_cache={})
            N, diag_step = _freeboard_reaction_step(
                N=N,
                T=T,
                P=P,
                V_seg=bg.V_seg,
                dt=dt_slice,
                fuel_type=fuel_type,
                rho_cat=rho_cat,
                char_area_total=char_area_total,
                solid_d_p=solid_d_p_eff,
                D_g=gas_diffusivity_correlation(T, P),
                char_conversion=char_conversion_seg,
                enabled_reactions=enabled_reactions,
            )
            T = _solve_T_from_enthalpy(N, H_before * (1.0 - substep_loss), T)
            if k == 0:
                diag_acc = {name: 0.0 for name in diag_step}
            for name, value in diag_step.items():
                diag_acc[name] += float(value)

        hold_char_before = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64).copy()
        hold_ash_before = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64).copy()
        inventory_sink_diag = (
            _freeboard_inventory_char_sink_diag(
                N=N,
                T=T,
                P=P,
                V_seg=bg.V_seg,
                dt=bg.tau,
                fuel_type=fuel_type,
                char_area_total=char_area_total,
                solid_d_p=solid_d_p_eff,
                D_g=gas_diffusivity_correlation(T, P),
                char_conversion=char_conversion_seg,
            )
            if inventory_char_sink_enabled
            else {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 0.0}
        )
        class_transport = _apply_char_sink_to_class_transport(class_transport, diag_acc=inventory_sink_diag)
        hold_char = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64)
        hold_ash = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64)

        secondary_observation["z_m"].append(float(bg.z_center_m))
        secondary_observation["xi_reactor"].append(float(xi_seg))
        secondary_observation["slice_global_index"].append(len(states))
        secondary_observation["segment_index"].append(int(i))
        secondary_observation["slice_index_within_segment"].append(0)
        secondary_observation["pre_mix_T"].append(float(pre_mix_T))
        secondary_observation["post_mix_pre_rxn_T"].append(float(post_mix_pre_rxn_T))
        secondary_observation["post_rxn_T"].append(float(T))
        post_rxn_wet = _wet_gas_profile(N)
        for sp in GAS_SPECIES:
            secondary_observation["pre_mix_wet_gas"][sp].append(float(pre_mix_wet[sp]))
            secondary_observation["post_mix_pre_rxn_wet_gas"][sp].append(float(post_mix_pre_rxn_wet[sp]))
            secondary_observation["post_rxn_wet_gas"][sp].append(float(post_rxn_wet[sp]))

        states.append(
            FreeboardState(
                segment_index=int(i),
                z_center_m=float(bg.z_center_m),
                T=float(T),
                N=N.copy(),
                u0=float(bg.u0),
                u_gb=float(bg.u_gb),
                u_p_mean=float(u_p_mean),
                u_t_mean=float(u_t_mean),
                carry_ratio=float(carry_ratio),
                tau=float(bg.tau),
                m_dot_char=float(m_dot_char),
                m_dot_ash=float(m_dot_ash),
                m_dot_return_char=float(return_char_seg),
                m_dot_return_ash=float(return_ash_seg),
                m_hold_char_before_classes=np.array(hold_char_before, copy=True),
                m_hold_ash_before_classes=np.array(hold_ash_before, copy=True),
                char_conversion_before_classes=np.array(x_before, copy=True),
                char_conversion_classes=np.array(class_transport["char_conversion_classes"], copy=True),
                m_hold_char_classes=np.array(hold_char, copy=True),
                m_hold_ash_classes=np.array(hold_ash, copy=True),
                m_char_sink_applied_classes=np.array(class_transport["m_char_sink_applied_classes"], copy=True),
                m_dot_auf_char_classes=np.array(up_char, copy=True),
                m_dot_auf_ash_classes=np.array(up_ash, copy=True),
                m_dot_ab_char_classes=np.array(down_char, copy=True),
                m_dot_ab_ash_classes=np.array(down_ash, copy=True),
                K_auf_classes=np.array(class_transport["K_auf_classes"], copy=True),
                K_ab_classes=np.array(class_transport["K_ab_classes"], copy=True),
                rho_solid=float(m_hold_seg / max(bg.V_seg, 1e-12)),
                rho_cat=float(rho_cat),
                reaction_diag=diag_acc,
            )
        )

    exit_samples = np.asarray(traj["exit_samples"], dtype=np.float64)
    entrained_exit_char = float(np.sum(exit_samples * char_frac_samples))
    entrained_exit_ash = float(np.sum(exit_samples * (1.0 - char_frac_samples)))
    return {
        "states": states,
        "exit_T": float(T),
        "exit_N": N,
        "exit_gas": _wet_gas_profile(N),
        "beta_A": float(beta_a),
        "trajectory_model": "global_projected_hamel",
        "trajectory_solver": "global_projected_hamel",
        "trajectory_coeff_model": str(trajectory_coeff_model),
        "trajectory_diag": {k: int(v) for k, v in traj["diag"].items()},
        "trajectory_coeff_diag": {
            "re_p_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["re_p_min"]) else traj["coeff_diag"]["re_p_min"]),
            "re_p_max": float(traj["coeff_diag"]["re_p_max"]),
            "cd_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["cd_min"]) else traj["coeff_diag"]["cd_min"]),
            "cd_max": float(traj["coeff_diag"]["cd_max"]),
            "re_p_active_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["re_p_active_min"]) else traj["coeff_diag"]["re_p_active_min"]),
            "re_p_active_max": float(traj["coeff_diag"]["re_p_active_max"]),
            "cd_active_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["cd_active_min"]) else traj["coeff_diag"]["cd_active_min"]),
            "cd_active_max": float(traj["coeff_diag"]["cd_active_max"]),
            "z_drag": float(traj["coeff_diag"]["z_drag"]),
            "exact_delta_pos": int(traj["coeff_diag"]["exact_delta_pos"]),
            "exact_delta_zero": int(traj["coeff_diag"]["exact_delta_zero"]),
            "exact_delta_neg": int(traj["coeff_diag"]["exact_delta_neg"]),
            "exact_a_rel_max": float(traj["coeff_diag"]["exact_a_rel_max"]),
            "exact_b_rel_max": float(traj["coeff_diag"]["exact_b_rel_max"]),
            "exact_c_rel_max": float(traj["coeff_diag"]["exact_c_rel_max"]),
            "exact_dp_grav_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["exact_dp_grav_min"]) else traj["coeff_diag"]["exact_dp_grav_min"]),
            "exact_dp_grav_max": float(traj["coeff_diag"]["exact_dp_grav_max"]),
            "exact_dp_ratio_min": float(0.0 if not math.isfinite(traj["coeff_diag"]["exact_dp_ratio_min"]) else traj["coeff_diag"]["exact_dp_ratio_min"]),
            "exact_dp_ratio_max": float(traj["coeff_diag"]["exact_dp_ratio_max"]),
        },
        "secondary_injection_applied": bool(secondary_injection_applied),
        "secondary_injection_segment": secondary_injection_segment,
        "secondary_local_refine": 1,
        "secondary_injection_mode": str(secondary_injection_mode),
        "entrained_eject_flux_kg_m2_s": float(m_dot_eject0 / max(A, 1e-12)),
        "entrained_eject_char_ash_kg_s": float(m_dot_eject0),
        "entrained_exit_char_kg_s": entrained_exit_char,
        "entrained_exit_ash_kg_s": entrained_exit_ash,
        "entrained_return_char_kg_s": float(total_return_char),
        "entrained_return_ash_kg_s": float(total_return_ash),
        "secondary_observation": secondary_observation,
        "profiles": {
            "T": [float(st.T) for st in states],
            "z_m": [float(st.z_center_m) for st in states],
            "tau_s": [float(st.tau) for st in states],
            "u0_m_s": [float(st.u0) for st in states],
            "u_gb_m_s": [float(st.u_gb) for st in states],
            "u_p_mean_m_s": [float(st.u_p_mean) for st in states],
            "u_t_mean_m_s": [float(st.u_t_mean) for st in states],
            "carry_ratio": [float(st.carry_ratio) for st in states],
            "hold_time_mean_s": [float(v) for v in traj["hold_time_mean_by_seg"]],
            "hold_time_max_s": [float(v) for v in traj["hold_time_max_by_seg"]],
            "reaction_diag": {
                name: [float(st.reaction_diag.get(name, 0.0)) for st in states]
                for name in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R11b", "R11d", "R12")
            },
            "entrained_char_kg_s": [float(st.m_dot_char) for st in states],
            "entrained_ash_kg_s": [float(st.m_dot_ash) for st in states],
            "entrained_return_char_kg_s": [float(st.m_dot_return_char) for st in states],
            "entrained_return_ash_kg_s": [float(st.m_dot_return_ash) for st in states],
            "entrained_solid_density_kg_m3": [float(st.rho_solid) for st in states],
            "entrained_catalyst_density_kg_m3": [float(st.rho_cat) for st in states],
            "solid_holdup_char_before_classes_kg": [np.asarray(st.m_hold_char_before_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_before_classes_kg": [np.asarray(st.m_hold_ash_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_before_classes": [np.asarray(st.char_conversion_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_classes": [np.asarray(st.char_conversion_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_char_classes_kg": [np.asarray(st.m_hold_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_classes_kg": [np.asarray(st.m_hold_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_char_sink_applied_classes_kg": [np.asarray(st.m_char_sink_applied_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_char_classes_kg_s": [np.asarray(st.m_dot_auf_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_ash_classes_kg_s": [np.asarray(st.m_dot_auf_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_char_classes_kg_s": [np.asarray(st.m_dot_ab_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_ash_classes_kg_s": [np.asarray(st.m_dot_ab_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_auf_classes_1_s": [np.asarray(st.K_auf_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_ab_classes_1_s": [np.asarray(st.K_ab_classes, dtype=np.float64).tolist() for st in states],
            "wet_gas": {
                sp: [
                    float(max(st.N[GAS_SPECIES_INDEX[sp]], 0.0) / max(np.sum(np.maximum(st.N, 0.0)), 1e-12))
                    for st in states
                ]
                for sp in GAS_SPECIES
            },
        },
    }


def simulate_freeboard(
    *,
    N_in: np.ndarray,
    T_in: float,
    P: float,
    D_bed: float,
    H_freeboard: float,
    n_cells: int,
    u_b_bed_top: float,
    d_b_bed_top: float,
    eps_b_bed_top: float,
    eps_d_void_bed_top: float,
    rho_solid_bed_top: float,
    d_p_classes_bed_top: np.ndarray,
    m_char_classes_bed_top: np.ndarray,
    m_ash_classes_bed_top: np.ndarray,
    char_conversion_bed_top: float = 0.0,
    class_aggregate_indices_bed_top: np.ndarray | None = None,
    n_output_size_classes: int | None = None,
    fuel_type: str,
    heat_loss_frac: float = 0.0,
    trajectory_model: str = "analytical_wirsum",
    trajectory_coeff_model: str = "stable_mixed_drag_split",
    u_gb_scale: float = 1.0,
    beta_a_scale: float = 1.0,
    velocity_sigma: float = 0.60,
    velocity_bins: int = 5,
    phi_s_bed_top: float = 0.86,
    z_base_m: float = 0.0,
    total_height_m: float | None = None,
    secondary_injection_xi: float | None = None,
    secondary_O2_mol_s: float = 0.0,
    secondary_N2_mol_s: float = 0.0,
    secondary_H2O_mol_s: float = 0.0,
    secondary_T_K: float = 293.15,
    secondary_local_refine: int = 1,
    secondary_injection_mode: str = "lumped",
    enabled_reactions: tuple[str, ...] | None = None,
    inventory_char_sink_enabled: bool = False,
) -> dict:
    if H_freeboard <= 0.0 or n_cells <= 0:
        return {
            "states": [],
            "exit_T": float(T_in),
            "exit_N": np.array(N_in, dtype=np.float64, copy=True),
            "exit_gas": {sp: float(v) for sp, v in zip(GAS_SPECIES, np.maximum(N_in, 0.0) / max(np.sum(np.maximum(N_in, 0.0)), 1e-12))},
        }

    if trajectory_model not in {"surrogate_exp", "force_balance", "analytical_wirsum", "global_projected_hamel"}:
        raise ValueError(f"未知的 freeboard trajectory_model: {trajectory_model}")
    if trajectory_coeff_model not in {"stable_mixed_drag_split", "exact_hamel"}:
        raise ValueError(f"未知的 analytical_wirsum coeff_model: {trajectory_coeff_model}")
    if secondary_injection_mode not in {"lumped", "distributed_uniform"}:
        raise ValueError(f"未知的 secondary_injection_mode: {secondary_injection_mode}")

    if trajectory_model == "global_projected_hamel":
        return _simulate_freeboard_global_projected_hamel(
            N_in=N_in,
            T_in=T_in,
            P=P,
            D_bed=D_bed,
            H_freeboard=H_freeboard,
            n_cells=n_cells,
            u_b_bed_top=u_b_bed_top,
            d_b_bed_top=d_b_bed_top,
            eps_b_bed_top=eps_b_bed_top,
            eps_d_void_bed_top=eps_d_void_bed_top,
            rho_solid_bed_top=rho_solid_bed_top,
            d_p_classes_bed_top=d_p_classes_bed_top,
            m_char_classes_bed_top=m_char_classes_bed_top,
            m_ash_classes_bed_top=m_ash_classes_bed_top,
            char_conversion_bed_top=char_conversion_bed_top,
            class_aggregate_indices_bed_top=class_aggregate_indices_bed_top,
            n_output_size_classes=n_output_size_classes,
            fuel_type=fuel_type,
            heat_loss_frac=heat_loss_frac,
            trajectory_coeff_model=trajectory_coeff_model,
            u_gb_scale=u_gb_scale,
            beta_a_scale=beta_a_scale,
            velocity_sigma=velocity_sigma,
            velocity_bins=velocity_bins,
            phi_s_bed_top=phi_s_bed_top,
            z_base_m=z_base_m,
            total_height_m=total_height_m,
            secondary_injection_xi=secondary_injection_xi,
            secondary_O2_mol_s=secondary_O2_mol_s,
            secondary_N2_mol_s=secondary_N2_mol_s,
            secondary_H2O_mol_s=secondary_H2O_mol_s,
            secondary_T_K=secondary_T_K,
            secondary_injection_mode=secondary_injection_mode,
            enabled_reactions=enabled_reactions,
            inventory_char_sink_enabled=inventory_char_sink_enabled,
        )

    A = math.pi * (D_bed**2) / 4.0
    dh = H_freeboard / n_cells
    beta_a = calc_beta_a(max(d_b_bed_top, 1e-6)) * float(max(beta_a_scale, 0.0))
    u_gb0 = calc_u_gb(max(u_b_bed_top, 1e-6), scale=float(max(u_gb_scale, 1e-6)))
    u_gb_prev = float(u_gb0)
    # Hamel Eq. 3.73/3.76: ghost-bubble gas speed and particle launch speed
    # are distinct. u_gb,0 = u_b,ws while mean particle launch speed is
    # u_p,0 ~= 1.53 * u_b,ws.
    u_p0_mean = 1.53 * max(u_b_bed_top, 1e-6)
    u_p_samples_prev, velocity_weights = _build_velocity_samples(
        u_p0_mean,
        sigma=float(max(velocity_sigma, 0.0)),
        n_bins=int(max(velocity_bins, 1)),
    )
    N = np.array(N_in, dtype=np.float64, copy=True)
    T = float(T_in)
    states: list[FreeboardState] = []
    seg_loss = 1.0 - math.pow(max(1.0 - heat_loss_frac, 0.0), 1.0 / max(n_cells, 1))
    d_p_classes = np.asarray(d_p_classes_bed_top, dtype=np.float64)
    m_char_classes = np.asarray(m_char_classes_bed_top, dtype=np.float64)
    m_ash_classes = np.asarray(m_ash_classes_bed_top, dtype=np.float64)
    d_p_class, char_frac_class, m_dot_class0, m_dot_eject0 = _build_size_class_bundles(
        d_p_classes=d_p_classes,
        m_char_classes=m_char_classes,
        m_ash_classes=m_ash_classes,
        rho_p=float(rho_solid_bed_top),
        eps_b=float(eps_b_bed_top),
        eps_d_void=float(eps_d_void_bed_top),
        u_b=float(u_b_bed_top),
        d_b=float(d_b_bed_top),
        area=A,
    )
    aggregate_idx, n_transport_classes, d_p_transport_class = _resolve_transport_class_mapping(
        d_p_class=d_p_class,
        m_dot_class=m_dot_class0,
        class_aggregate_indices=class_aggregate_indices_bed_top,
        n_output_size_classes=n_output_size_classes,
    )
    m_dot_samples0 = np.repeat(m_dot_class0, len(velocity_weights)) * np.tile(velocity_weights, len(m_dot_class0))
    char_frac_samples = np.repeat(char_frac_class, len(velocity_weights))
    char_conversion_samples = np.full_like(m_dot_samples0, float(np.clip(char_conversion_bed_top, 0.0, 1.0 - 1e-9)), dtype=np.float64)
    class_idx_samples = np.repeat(aggregate_idx, len(velocity_weights))
    d_p_samples = np.repeat(d_p_class, len(velocity_weights))
    u_p_samples_prev = np.tile(u_p_samples_prev, len(m_dot_class0))
    slower_ref_u_prev = np.array(u_p_samples_prev, dtype=np.float64, copy=True)
    slower_ref_z_prev = np.zeros_like(u_p_samples_prev, dtype=np.float64)
    f0 = m_dot_eject0 / max(A, 1e-12)
    total_return_char = 0.0
    total_return_ash = 0.0
    secondary_injection_applied = False
    secondary_injection_segment: int | None = None
    secondary_N = np.zeros(N_GAS, dtype=np.float64)
    secondary_N[GAS_SPECIES_INDEX["O2"]] = float(max(secondary_O2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["N2"]] = float(max(secondary_N2_mol_s, 0.0))
    secondary_N[GAS_SPECIES_INDEX["H2O"]] = float(max(secondary_H2O_mol_s, 0.0))
    xi_inj = None if secondary_injection_xi is None else float(np.clip(secondary_injection_xi, 0.0, 1.0))
    h_total = None if total_height_m is None else float(max(total_height_m, 1e-12))
    refine_local = int(max(secondary_local_refine, 1))
    secondary_observation = {
        "z_m": [],
        "xi_reactor": [],
        "slice_global_index": [],
        "segment_index": [],
        "slice_index_within_segment": [],
        "pre_mix_T": [],
        "post_mix_pre_rxn_T": [],
        "post_rxn_T": [],
        "pre_mix_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_mix_pre_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
        "post_rxn_wet_gas": {sp: [] for sp in GAS_SPECIES},
    }
    trajectory_diag = {
        "delta_pos": 0,
        "delta_zero": 0,
        "delta_neg": 0,
        "disc_pos": 0,
        "disc_zero": 0,
        "disc_neg": 0,
        "sign_switch": 0,
        "fallbacks": 0,
        "returns": 0,
        "fallback_no_event": 0,
        "fallback_boundary_velocity": 0,
        "fallback_delta_neg_time": 0,
    }
    trajectory_coeff_diag = {
        "re_p_min": float("inf"),
        "re_p_max": 0.0,
        "cd_min": float("inf"),
        "cd_max": 0.0,
        "re_p_active_min": float("inf"),
        "re_p_active_max": 0.0,
        "cd_active_min": float("inf"),
        "cd_active_max": 0.0,
        "z_drag": 0.0,
        "exact_delta_pos": 0,
        "exact_delta_zero": 0,
        "exact_delta_neg": 0,
        "exact_a_rel_max": 0.0,
        "exact_b_rel_max": 0.0,
        "exact_c_rel_max": 0.0,
        "exact_dp_grav_min": float("inf"),
        "exact_dp_grav_max": 0.0,
        "exact_dp_ratio_min": float("inf"),
        "exact_dp_ratio_max": 0.0,
    }

    for i in range(n_cells):
        n_tot = max(float(np.sum(np.maximum(N, 0.0))), 1e-12)
        u0 = (n_tot * Rg * T / P) / max(A, 1e-12)
        u_gb_in = float(u_gb_prev if i > 0 else u_gb0)
        exp_arg_out = float(np.clip(-beta_a * dh, -200.0, 200.0))
        exp_arg_center = float(np.clip(-beta_a * 0.5 * dh, -200.0, 200.0))
        u_gb_out = u0 + (u_gb_in - u0) * math.exp(exp_arg_out)
        u_gb = u0 + (u_gb_in - u0) * math.exp(exp_arg_center)
        z_center = (i + 0.5) * dh
        xi_local = z_center / max(H_freeboard, 1e-12)
        xi_seg = xi_local if h_total is None else float((float(z_base_m) + z_center) / h_total)
        tau = dh / max(0.5 * (u_gb_in + u_gb_out), 1e-6)
        V_seg = A * dh
        dt = tau / 12.0

        injection_applied_here = False
        if (
            not secondary_injection_applied
            and xi_inj is not None
            and float(np.sum(secondary_N)) > 0.0
            and xi_seg >= xi_inj
        ):
            secondary_injection_applied = True
            injection_applied_here = True
            secondary_injection_segment = i

        y_map = {sp: float(max(N[GAS_SPECIES_INDEX[sp]], 0.0) / n_tot) for sp in GAS_SPECIES}
        rho_g = gas_density_ideal(P, T, y_map)
        mu_g = gas_viscosity_power_law(T, 1.8e-5)
        u_t_samples = np.array(
            [
                calc_terminal_velocity_haider(
                    rho_g=rho_g,
                    rho_p=max(float(rho_solid_bed_top), 1e-6),
                    d_p=max(float(dp), 1e-9),
                    phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                    mu_g=mu_g,
                    g_acc=g,
                )
                for dp in d_p_samples
            ],
            dtype=np.float64,
        )
        if trajectory_model == "analytical_wirsum":
            u_p_samples, m_dot_samples, hold_samples, return_samples, slower_ref_u, slower_ref_z, diag_seg = _advance_particle_samples_analytical_wirsum(
                dh=dh,
                z_lower_abs=float(i * dh),
                u_g_prev=u_gb_in,
                u_g_next=u_gb_out,
                rho_g=rho_g,
                rho_p=max(float(rho_solid_bed_top), 1e-6),
                phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                mu_g=mu_g,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
                d_p=d_p_samples,
                char_frac_samples=char_frac_samples,
                slower_ref_u_prev=slower_ref_u_prev,
                slower_ref_z_prev=slower_ref_z_prev,
                coeff_model=str(trajectory_coeff_model),
            )
            for key, val in diag_seg.items():
                if key in trajectory_diag:
                    trajectory_diag[key] += int(val)
                elif key in trajectory_coeff_diag:
                    if key.endswith("_min"):
                        trajectory_coeff_diag[key] = min(float(trajectory_coeff_diag[key]), float(val))
                    elif key.endswith("_max"):
                        trajectory_coeff_diag[key] = max(float(trajectory_coeff_diag[key]), float(val))
                    elif key.startswith("exact_delta_"):
                        trajectory_coeff_diag[key] += int(val)
                    else:
                        trajectory_coeff_diag[key] = float(val)
        elif trajectory_model == "force_balance":
            u_p_samples, m_dot_samples, hold_samples, return_samples = _advance_particle_samples_force_balance(
                dh=dh,
                u_g_prev=u_gb_in,
                u_g_next=u_gb_out,
                rho_g=rho_g,
                rho_p=max(float(rho_solid_bed_top), 1e-6),
                phi_s=float(np.clip(phi_s_bed_top, 1e-3, 1.0)),
                mu_g=mu_g,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
                d_p=d_p_samples,
                char_frac_samples=char_frac_samples,
            )
            trajectory_diag["fallbacks"] += 0
        else:
            u_p_samples, m_dot_samples, hold_samples, return_samples = _advance_particle_samples_surrogate(
                dh=dh,
                beta_a=beta_a,
                u_g_asym=u0,
                u_g_prev=u_gb_in,
                u_p_samples_prev=u_p_samples_prev,
                m_dot_samples_prev=m_dot_samples0,
            )
        m_hold_seg = float(np.sum(hold_samples))
        return_char_seg = float(np.sum(np.maximum(return_samples, 0.0) * char_frac_samples))
        return_ash_seg = float(np.sum(np.maximum(return_samples, 0.0) * (1.0 - char_frac_samples)))
        class_transport = _aggregate_class_transport(
            n_classes=n_transport_classes,
            class_idx_samples=class_idx_samples,
            hold_samples=hold_samples,
            upflow_samples=m_dot_samples,
            downflow_samples=return_samples,
            char_frac_samples=char_frac_samples,
            char_conversion_samples=char_conversion_samples,
        )
        char_area_total, solid_d_p_eff = _char_area_and_representative_diameter(
            np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64),
            d_p_transport_class,
            max(float(rho_solid_bed_top), 1e-6),
        )
        char_conversion_seg = _representative_char_conversion(class_transport)
        flow_avg = 0.5 * (m_dot_samples0 + m_dot_samples)
        m_dot_s_avg = float(np.sum(flow_avg))
        rho_solid_seg = m_hold_seg / max(V_seg, 1e-12)
        m_dot_char = float(np.sum(flow_avg * char_frac_samples))
        m_dot_ash = float(np.sum(flow_avg * (1.0 - char_frac_samples)))
        ash_hold_frac = m_dot_ash / max(m_dot_s_avg, 1e-12)
        rho_cat = (m_hold_seg * ash_hold_frac) / max(V_seg, 1e-12)
        u_p_mean = float(np.sum(flow_avg * np.maximum(u_p_samples, 0.0)) / max(m_dot_s_avg, 1e-12))
        u_t_mean = float(np.sum(flow_avg * u_t_samples) / max(m_dot_s_avg, 1e-12))
        carry_ratio = float(u_gb / max(u_t_mean, 1e-12))
        total_return_char += float(return_char_seg)
        total_return_ash += float(return_ash_seg)

        n_slices = refine_local if injection_applied_here else 1
        dh_slice = dh / float(n_slices)
        tau_slice = tau / float(n_slices)
        v_slice = A * dh_slice
        seg_loss_slice = 1.0 - math.pow(max(1.0 - seg_loss, 0.0), 1.0 / float(max(n_slices, 1)))
        substep_loss = 1.0 - math.pow(max(1.0 - seg_loss_slice, 0.0), 1.0 / 12.0)
        secondary_slice_N = secondary_N / float(n_slices) if (injection_applied_here and secondary_injection_mode == "distributed_uniform") else secondary_N

        for j in range(n_slices):
            slice_z_center = i * dh + (j + 0.5) * dh_slice
            slice_xi = (
                float((float(z_base_m) + slice_z_center) / h_total)
                if h_total is not None
                else float(slice_z_center / max(H_freeboard, 1e-12))
            )
            pre_mix_T = float(T)
            pre_mix_wet = _wet_gas_profile(N)
            mixed_here = False
            if injection_applied_here:
                if secondary_injection_mode == "lumped":
                    if j == 0:
                        H_mix = (
                            calc_gas_enthalpy_flow(N, T, h_cache={})
                            + calc_gas_enthalpy_flow(secondary_N, float(secondary_T_K), h_cache={})
                        )
                        N = N + secondary_N
                        T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
                        mixed_here = True
                else:
                    H_mix = (
                        calc_gas_enthalpy_flow(N, T, h_cache={})
                        + calc_gas_enthalpy_flow(secondary_slice_N, float(secondary_T_K), h_cache={})
                    )
                    N = N + secondary_slice_N
                    T = _solve_T_from_enthalpy(N, H_mix, max(T, float(secondary_T_K)))
                    mixed_here = True
            post_mix_pre_rxn_T = float(T)
            post_mix_pre_rxn_wet = _wet_gas_profile(N)
            dt_slice = tau_slice / 12.0
            diag_acc = {}
            for k in range(12):
                H_before = calc_gas_enthalpy_flow(N, T, h_cache={})
                N, diag_step = _freeboard_reaction_step(
                    N=N,
                    T=T,
                    P=P,
                    V_seg=v_slice,
                    dt=dt_slice,
                    fuel_type=fuel_type,
                    rho_cat=rho_cat,
                    char_area_total=char_area_total,
                    solid_d_p=solid_d_p_eff,
                    D_g=gas_diffusivity_correlation(T, P),
                    char_conversion=char_conversion_seg,
                    enabled_reactions=enabled_reactions,
                )
                T = _solve_T_from_enthalpy(N, H_before * (1.0 - substep_loss), T)
                if k == 0:
                    diag_acc = {name: 0.0 for name in diag_step}
                for name, value in diag_step.items():
                    diag_acc[name] += float(value)

            hold_char_before = np.asarray(class_transport["m_hold_char_classes"], dtype=np.float64).copy()
            hold_ash_before = np.asarray(class_transport["m_hold_ash_classes"], dtype=np.float64).copy()
            x_before = np.asarray(class_transport["char_conversion_classes"], dtype=np.float64).copy()
            inventory_sink_diag = (
                _freeboard_inventory_char_sink_diag(
                    N=N,
                    T=T,
                    P=P,
                    V_seg=v_slice,
                    dt=tau_slice,
                    fuel_type=fuel_type,
                    char_area_total=char_area_total,
                    solid_d_p=solid_d_p_eff,
                    D_g=gas_diffusivity_correlation(T, P),
                    char_conversion=char_conversion_seg,
                )
                if inventory_char_sink_enabled
                else {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 0.0}
            )
            class_transport = _apply_char_sink_to_class_transport(class_transport, diag_acc=inventory_sink_diag)

            if mixed_here:
                secondary_observation["z_m"].append(float(slice_z_center))
                secondary_observation["xi_reactor"].append(float(slice_xi))
                secondary_observation["slice_global_index"].append(len(states))
                secondary_observation["segment_index"].append(int(i))
                secondary_observation["slice_index_within_segment"].append(int(j))
                secondary_observation["pre_mix_T"].append(float(pre_mix_T))
                secondary_observation["post_mix_pre_rxn_T"].append(float(post_mix_pre_rxn_T))
                secondary_observation["post_rxn_T"].append(float(T))
                post_rxn_wet = _wet_gas_profile(N)
                for sp in GAS_SPECIES:
                    secondary_observation["pre_mix_wet_gas"][sp].append(float(pre_mix_wet[sp]))
                    secondary_observation["post_mix_pre_rxn_wet_gas"][sp].append(float(post_mix_pre_rxn_wet[sp]))
                    secondary_observation["post_rxn_wet_gas"][sp].append(float(post_rxn_wet[sp]))

            states.append(
                FreeboardState(
                    segment_index=int(i),
                    z_center_m=slice_z_center,
                    T=float(T),
                    N=N.copy(),
                    u0=float(u0),
                    u_gb=float(u_gb),
                    u_p_mean=float(u_p_mean),
                    u_t_mean=float(u_t_mean),
                    carry_ratio=float(carry_ratio),
                    tau=float(tau_slice),
                    m_dot_char=float(m_dot_char),
                    m_dot_ash=float(m_dot_ash),
                    m_dot_return_char=float(return_char_seg / float(n_slices)),
                    m_dot_return_ash=float(return_ash_seg / float(n_slices)),
                    m_hold_char_before_classes=np.array(hold_char_before, copy=True),
                    m_hold_ash_before_classes=np.array(hold_ash_before, copy=True),
                    char_conversion_before_classes=np.array(x_before, copy=True),
                    char_conversion_classes=np.array(class_transport["char_conversion_classes"], copy=True),
                    m_hold_char_classes=np.array(class_transport["m_hold_char_classes"], copy=True),
                    m_hold_ash_classes=np.array(class_transport["m_hold_ash_classes"], copy=True),
                    m_char_sink_applied_classes=np.array(class_transport["m_char_sink_applied_classes"], copy=True),
                    m_dot_auf_char_classes=np.array(class_transport["m_dot_auf_char_classes"], copy=True),
                    m_dot_auf_ash_classes=np.array(class_transport["m_dot_auf_ash_classes"], copy=True),
                    m_dot_ab_char_classes=np.array(class_transport["m_dot_ab_char_classes"], copy=True),
                    m_dot_ab_ash_classes=np.array(class_transport["m_dot_ab_ash_classes"], copy=True),
                    K_auf_classes=np.array(class_transport["K_auf_classes"], copy=True),
                    K_ab_classes=np.array(class_transport["K_ab_classes"], copy=True),
                    rho_solid=float(rho_solid_seg),
                    rho_cat=float(rho_cat),
                    reaction_diag=diag_acc,
                )
            )
        u_p_samples_prev = u_p_samples
        m_dot_samples0 = m_dot_samples
        if trajectory_model == "analytical_wirsum":
            slower_ref_u_prev = slower_ref_u
            slower_ref_z_prev = slower_ref_z
        u_gb_prev = float(u_gb_out)

    exit_gas = _wet_gas_profile(N)
    entrained_exit_total = float(np.sum(m_dot_samples0))
    entrained_exit_char = float(np.sum(m_dot_samples0 * char_frac_samples))
    entrained_exit_ash = float(np.sum(m_dot_samples0 * (1.0 - char_frac_samples)))
    return {
        "states": states,
        "exit_T": float(T),
        "exit_N": N,
        "exit_gas": exit_gas,
        "beta_A": float(beta_a),
        "trajectory_model": str(trajectory_model),
        "trajectory_solver": (
            ("wirsum_analytical_exact_hamel_coeffs" if trajectory_coeff_model == "exact_hamel" else "wirsum_analytical")
            if trajectory_model == "analytical_wirsum"
            else ("force_balance_numeric" if trajectory_model == "force_balance" else "surrogate_exp")
        ),
        "trajectory_coeff_model": (
            str(trajectory_coeff_model)
            if trajectory_model == "analytical_wirsum"
            else ("full_force_balance" if trajectory_model == "force_balance" else "surrogate_decay")
        ),
        "trajectory_diag": {k: int(v) for k, v in trajectory_diag.items()},
        "trajectory_coeff_diag": {
            "re_p_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["re_p_min"]) else trajectory_coeff_diag["re_p_min"]),
            "re_p_max": float(trajectory_coeff_diag["re_p_max"]),
            "cd_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["cd_min"]) else trajectory_coeff_diag["cd_min"]),
            "cd_max": float(trajectory_coeff_diag["cd_max"]),
            "re_p_active_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["re_p_active_min"]) else trajectory_coeff_diag["re_p_active_min"]),
            "re_p_active_max": float(trajectory_coeff_diag["re_p_active_max"]),
            "cd_active_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["cd_active_min"]) else trajectory_coeff_diag["cd_active_min"]),
            "cd_active_max": float(trajectory_coeff_diag["cd_active_max"]),
            "z_drag": float(trajectory_coeff_diag["z_drag"]),
            "exact_delta_pos": int(trajectory_coeff_diag["exact_delta_pos"]),
            "exact_delta_zero": int(trajectory_coeff_diag["exact_delta_zero"]),
            "exact_delta_neg": int(trajectory_coeff_diag["exact_delta_neg"]),
            "exact_a_rel_max": float(trajectory_coeff_diag["exact_a_rel_max"]),
            "exact_b_rel_max": float(trajectory_coeff_diag["exact_b_rel_max"]),
            "exact_c_rel_max": float(trajectory_coeff_diag["exact_c_rel_max"]),
            "exact_dp_grav_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["exact_dp_grav_min"]) else trajectory_coeff_diag["exact_dp_grav_min"]),
            "exact_dp_grav_max": float(trajectory_coeff_diag["exact_dp_grav_max"]),
            "exact_dp_ratio_min": float(0.0 if not math.isfinite(trajectory_coeff_diag["exact_dp_ratio_min"]) else trajectory_coeff_diag["exact_dp_ratio_min"]),
            "exact_dp_ratio_max": float(trajectory_coeff_diag["exact_dp_ratio_max"]),
        },
        "secondary_injection_applied": bool(secondary_injection_applied),
        "secondary_injection_segment": secondary_injection_segment,
        "secondary_local_refine": int(refine_local),
        "secondary_injection_mode": str(secondary_injection_mode),
        "entrained_eject_flux_kg_m2_s": float(f0),
        "entrained_eject_char_ash_kg_s": float(m_dot_eject0),
        "entrained_exit_char_kg_s": entrained_exit_char,
        "entrained_exit_ash_kg_s": entrained_exit_ash,
        "entrained_return_char_kg_s": float(total_return_char),
        "entrained_return_ash_kg_s": float(total_return_ash),
        "secondary_observation": secondary_observation,
        "profiles": {
            "T": [float(st.T) for st in states],
            "z_m": [float(st.z_center_m) for st in states],
            "tau_s": [float(st.tau) for st in states],
            "u0_m_s": [float(st.u0) for st in states],
            "u_gb_m_s": [float(st.u_gb) for st in states],
            "u_p_mean_m_s": [float(st.u_p_mean) for st in states],
            "u_t_mean_m_s": [float(st.u_t_mean) for st in states],
            "carry_ratio": [float(st.carry_ratio) for st in states],
            "reaction_diag": {
                name: [float(st.reaction_diag.get(name, 0.0)) for st in states]
                for name in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R11b", "R11d", "R12")
            },
            "entrained_char_kg_s": [float(st.m_dot_char) for st in states],
            "entrained_ash_kg_s": [float(st.m_dot_ash) for st in states],
            "entrained_return_char_kg_s": [float(st.m_dot_return_char) for st in states],
            "entrained_return_ash_kg_s": [float(st.m_dot_return_ash) for st in states],
            "entrained_solid_density_kg_m3": [float(st.rho_solid) for st in states],
            "entrained_catalyst_density_kg_m3": [float(st.rho_cat) for st in states],
            "solid_holdup_char_before_classes_kg": [np.asarray(st.m_hold_char_before_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_before_classes_kg": [np.asarray(st.m_hold_ash_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_before_classes": [np.asarray(st.char_conversion_before_classes, dtype=np.float64).tolist() for st in states],
            "char_conversion_classes": [np.asarray(st.char_conversion_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_char_classes_kg": [np.asarray(st.m_hold_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_holdup_ash_classes_kg": [np.asarray(st.m_hold_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_char_sink_applied_classes_kg": [np.asarray(st.m_char_sink_applied_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_char_classes_kg_s": [np.asarray(st.m_dot_auf_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_upflow_ash_classes_kg_s": [np.asarray(st.m_dot_auf_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_char_classes_kg_s": [np.asarray(st.m_dot_ab_char_classes, dtype=np.float64).tolist() for st in states],
            "solid_downflow_ash_classes_kg_s": [np.asarray(st.m_dot_ab_ash_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_auf_classes_1_s": [np.asarray(st.K_auf_classes, dtype=np.float64).tolist() for st in states],
            "solid_K_ab_classes_1_s": [np.asarray(st.K_ab_classes, dtype=np.float64).tolist() for st in states],
            "wet_gas": {
                sp: [
                    float(max(st.N[GAS_SPECIES_INDEX[sp]], 0.0) / max(np.sum(np.maximum(st.N, 0.0)), 1e-12))
                    for st in states
                ]
                for sp in GAS_SPECIES
            },
        },
    }
