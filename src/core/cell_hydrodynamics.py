"""单 cell 流体力学与相间交换。

Ref: Hamel (1999) Eq. 3.44, 3.50, 3.52, 3.75
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import GAS_SPECIES, N_GAS, gas_density_ideal, gas_diffusivity_correlation, gas_viscosity_power_law
from src.physics.freeboard import calc_beta_a, calc_u_gb
from src.physics.bubble_dynamics import (
    bubble_interaction_factor,
    bubble_rise_velocity,
    integrate_bubble_diameter,
    mori_wen_bubble_diameter,
)
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.physics.minimum_fluidization import compute_u_mf
from src.physics.phase_fractions import calc_emulsion_porosity, calc_epsilon_b, calc_n_rz, calc_visible_bubble_fraction


@dataclass
class HydrodynamicsBundle:
    """流体力学与交换所需状态。"""

    u_mf: float
    u_b: float
    d_b: float
    eps_b: float
    eps_d: float
    eps_d_voidage: float
    K_bd: float
    V_b: float
    V_d: float
    u0: float
    u_d: float
    n_rz: float


def calc_cell_hydrodynamics(
    *,
    T: float,
    P: float,
    N_b: npt.ArrayLike,
    N_d: npt.ArrayLike,
    D_bed: float,
    dh: float,
    h_center: float,
    N_or: int,
    rho_s: float,
    d_p: float,
    eps_mf: float,
    phi_s: float,
    u0_target: float | None = None,
    u_d_closure: str = "current",
    bubble_diameter_model: str = "mori_wen",
    psi_b_strategy: str = "wein_1992",
    lambda_strategy: str = "hamel_280",
    xi_strategy: str = "hamel_regime",
    bubble_velocity_strategy: str = "hilligardt_eq313",
    bubble_ode_strategy: str = "hilligardt_eq333",
    cell_type: str = "bed",
    freeboard_u_bed_top: float | None = None,
    freeboard_d_b_bed_top: float | None = None,
    freeboard_height_from_bed: float | None = None,
    freeboard_u_gb_scale: float = 1.0,
    freeboard_eps_d_voidage: float | None = None,
    dense_segment_eps_d_voidage: float | None = None,
) -> HydrodynamicsBundle:
    """计算单格流体力学状态。

    Ref: Hamel (1999) Eq. 3.44, 3.50, 3.52, 3.75
    """
    N_b_arr = np.asarray(N_b, dtype=np.float64)
    N_d_arr = np.asarray(N_d, dtype=np.float64)
    assert N_b_arr.shape == (N_GAS,)
    assert N_d_arr.shape == (N_GAS,)

    A_bed = np.pi / 4.0 * float(D_bed) ** 2
    V_cell = A_bed * float(dh)
    if u0_target is None:
        u0 = ((np.sum(np.maximum(N_b_arr, 0.0)) + np.sum(np.maximum(N_d_arr, 0.0))) * Rg * T / P) / A_bed
    else:
        u0 = float(u0_target)
    segment = str(cell_type).strip().lower()

    if segment == "freeboard":
        assert freeboard_u_bed_top is not None, "freeboard_u_bed_top is required for freeboard hydrodynamics"
        assert freeboard_d_b_bed_top is not None, "freeboard_d_b_bed_top is required for freeboard hydrodynamics"
        assert freeboard_eps_d_voidage is not None, "freeboard_eps_d_voidage is required for freeboard hydrodynamics"
        h_rel = max(
            float(h_center) if freeboard_height_from_bed is None else float(freeboard_height_from_bed),
            0.0,
        )
        u_gb0 = calc_u_gb(float(max(freeboard_u_bed_top, 0.0)), scale=float(max(freeboard_u_gb_scale, 1e-12)))
        beta_a = calc_beta_a(float(max(freeboard_d_b_bed_top, 1e-12)))
        exp_arg_u = float(np.clip(-beta_a * h_rel, -200.0, 200.0))
        u_g = float(u0 + (u_gb0 - u0) * math.exp(exp_arg_u))
        eps_d_voidage = float(np.clip(freeboard_eps_d_voidage, 0.0, 1.0))
        return HydrodynamicsBundle(
            u_mf=0.0,
            u_b=u_g,
            d_b=0.0,
            eps_b=0.0,
            eps_d=1.0,
            eps_d_voidage=eps_d_voidage,
            K_bd=0.0,
            V_b=0.0,
            V_d=float(V_cell),
            u0=float(u0),
            u_d=u_g,
            n_rz=0.0,
        )

    if segment == "cyclone":
        eps_d_voidage = 1.0 if dense_segment_eps_d_voidage is None else float(np.clip(dense_segment_eps_d_voidage, 0.0, 1.0))
        return HydrodynamicsBundle(
            u_mf=0.0,
            u_b=float(u0),
            d_b=0.0,
            eps_b=0.0,
            eps_d=1.0,
            eps_d_voidage=eps_d_voidage,
            K_bd=0.0,
            V_b=0.0,
            V_d=float(V_cell),
            u0=float(u0),
            u_d=float(u0),
            n_rz=0.0,
        )

    if segment == "return_leg":
        eps_d_voidage = float(eps_mf if dense_segment_eps_d_voidage is None else np.clip(dense_segment_eps_d_voidage, 0.0, 1.0))
        return HydrodynamicsBundle(
            u_mf=0.0,
            u_b=0.0,
            d_b=0.0,
            eps_b=0.0,
            eps_d=1.0,
            eps_d_voidage=eps_d_voidage,
            K_bd=0.0,
            V_b=0.0,
            V_d=float(V_cell),
            u0=0.0,
            u_d=0.0,
            n_rz=0.0,
        )

    assert segment == "bed", f"Unsupported cell_type={cell_type!r}"
    y_d = np.maximum(N_d_arr, 0.0)
    y_d_sum = float(np.sum(y_d))
    if y_d_sum < 1e-12:
        y_map = {sp: 1.0 / len(GAS_SPECIES) for sp in GAS_SPECIES}
    else:
        y_map = {sp: float(y) for sp, y in zip(GAS_SPECIES, y_d / y_d_sum) if y > 0.0}

    rho_g = gas_density_ideal(P, T, y_map)
    mu_g = gas_viscosity_power_law(T, 1.8e-5)
    D_g = gas_diffusivity_correlation(T, P)
    u_mf = compute_u_mf(rho_g, rho_s, d_p, mu_g, eps_mf, phi_s)
    psi_b = bubble_interaction_factor(u_mf, strategy=psi_b_strategy)
    d_b, u_b, u_d = _resolve_bubble_chain(
        strategy=bubble_diameter_model,
        u0=u0,
        u_mf=u_mf,
        P=P,
        D_bed=D_bed,
        dh=dh,
        h_center=h_center,
        N_or=N_or,
        A_bed=A_bed,
        eps_mf=eps_mf,
        psi_b=float(psi_b),
        u_d_closure=u_d_closure,
        lambda_strategy=lambda_strategy,
        xi_strategy=xi_strategy,
        bubble_velocity_strategy=bubble_velocity_strategy,
        bubble_ode_strategy=bubble_ode_strategy,
    )
    re_s = max(rho_g * max(u_mf, 1e-9) * max(d_p, 1e-9) / max(mu_g, 1e-12), 0.0)
    n_rz = calc_n_rz(re_s)
    # Hamel Eq. 3.24 explicitly defines the visible bubble fraction via ``u_d``
    # and the bubble through-flow relation Eq. 3.23; using the excess-gas proxy
    # ``(u0-u_mf)/u_b`` here overstates bed bubble hold-up for the LU thesis path.
    eps_b = float(np.clip(calc_visible_bubble_fraction(u0, u_b, u_d), 0.01, 0.7))
    eps_d = 1.0 - eps_b
    eps_d_voidage = calc_emulsion_porosity(u_d, u_mf, eps_mf, n_rz)
    V_b = eps_b * V_cell
    V_d = eps_d * V_cell
    u_br = calc_u_br(u_d, P)
    K_bd = calc_kbd(u_br, d_b, D_g, eps_mf, u_b)

    return HydrodynamicsBundle(
        u_mf=float(u_mf),
        u_b=float(u_b),
        d_b=float(d_b),
        eps_b=eps_b,
        eps_d=float(eps_d),
        eps_d_voidage=float(eps_d_voidage),
        K_bd=float(K_bd),
        V_b=float(V_b),
        V_d=float(V_d),
        u0=float(u0),
        u_d=float(u_d),
        n_rz=float(n_rz),
    )


def _resolve_bubble_chain(
    *,
    strategy: str,
    u0: float,
    u_mf: float,
    P: float,
    D_bed: float,
    dh: float,
    h_center: float,
    N_or: int,
    A_bed: float,
    eps_mf: float,
    psi_b: float,
    u_d_closure: str,
    lambda_strategy: str,
    xi_strategy: str,
    bubble_velocity_strategy: str,
    bubble_ode_strategy: str,
) -> tuple[float, float, float]:
    """Resolve local ``d_b / u_b / u_d`` with either legacy MW or audit ODE path.

    Notes
    -----
    ``hilligardt_ode`` is exposed as an audit branch. It uses a short fixed-point
    iteration because ``u_d`` depends on ``u_b`` while the ODE itself also uses
    ``u_d`` through ``alpha_b = u_b / u_d``. Default main-path behavior remains
    ``mori_wen``.
    """
    name = str(strategy).strip().lower()

    if name == "mori_wen":
        d_b = mori_wen_bubble_diameter(h_center, u0, u_mf, D_bed, A_bed=A_bed, N_or=N_or)
        u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b, P=P, strategy=bubble_velocity_strategy)
        u_d = _resolve_u_d_closure(
            strategy=u_d_closure,
            u0=u0,
            u_mf=u_mf,
            u_b=u_b,
            eps_mf=eps_mf,
        )
        return float(d_b), float(u_b), float(u_d)

    if name == "hilligardt_ode":
        d_b = mori_wen_bubble_diameter(h_center, u0, u_mf, D_bed, A_bed=A_bed, N_or=N_or)
        u_b = bubble_rise_velocity(u0, u_mf, d_b, psi_b, P=P, strategy=bubble_velocity_strategy)
        h_span = max(float(h_center) + 0.5 * float(dh), float(dh))
        n_points = max(8, int(np.ceil(h_span / max(float(dh), 1e-9))) + 4)

        for _ in range(5):
            u_d = _resolve_u_d_closure(
                strategy=u_d_closure,
                u0=u0,
                u_mf=u_mf,
                u_b=u_b,
                eps_mf=eps_mf,
            )
            h_arr, d_arr = integrate_bubble_diameter(
                u0=u0,
                u_mf=u_mf,
                P=P,
                H_bed=h_span,
                u_d=u_d,
                D_bed=D_bed,
                n_points=n_points,
                method="hilligardt_ode",
                lambda_strategy=lambda_strategy,
                xi_strategy=xi_strategy,
                velocity_strategy=bubble_velocity_strategy,
                ode_strategy=bubble_ode_strategy,
                psi_b=psi_b,
            )
            d_new = float(np.interp(float(h_center), h_arr, d_arr))
            u_b_new = float(
                bubble_rise_velocity(u0, u_mf, d_new, psi_b, P=P, strategy=bubble_velocity_strategy)
            )
            if abs(d_new - d_b) <= max(1e-6, 1e-3 * max(d_b, 1e-6)):
                d_b = d_new
                u_b = u_b_new
                break
            d_b = d_new
            u_b = u_b_new

        u_d = _resolve_u_d_closure(
            strategy=u_d_closure,
            u0=u0,
            u_mf=u_mf,
            u_b=u_b,
            eps_mf=eps_mf,
        )
        return float(d_b), float(u_b), float(u_d)

    raise ValueError(
        f"Unsupported bubble_diameter_model={strategy!r}; "
        "expected 'mori_wen' or 'hilligardt_ode'"
    )


def _resolve_u_d_closure(
    *,
    strategy: str,
    u0: float,
    u_mf: float,
    u_b: float,
    eps_mf: float,
    n_b: float = 2.7,
) -> float:
    """Resolve the suspension/emulsion gas velocity closure `u_d`.

    Notes
    -----
    The current main path historically used:
        u_d = max(u_mf / eps_mf, min(u0, u_b))

    which collapses to `u_d ~= u0` for the solved LU state because `u_b > u0`.
    Alternative closures are exposed for audit without changing the default.

    Literature-oriented audit closures:
    - ``hilligardt_eq311``:
        u_d = u_mf + (u0 - u_mf) / 3
    - ``wein_1992_eq312``:
        u_d = 1.45 * u_mf
    """
    name = str(strategy).strip().lower()
    u0 = float(max(u0, 0.0))
    u_mf = float(max(u_mf, 0.0))
    u_b = float(max(u_b, 1e-12))
    eps_mf = float(max(eps_mf, 1e-12))

    if name == "current":
        return float(max(u_mf / eps_mf, min(u0, u_b)))

    if name == "umf_over_epsmf":
        return float(u_mf / eps_mf)

    if name == "hilligardt_eq311":
        return float(u_mf + (u0 - u_mf) / 3.0)

    if name in {"wein_1992_eq312", "wein_1992", "wein_eq312"}:
        u_d_potential = float(1.45 * u_mf)
        if u0 > u_d_potential:
            return u_d_potential
        # Hamel thesis logic for Geldart-B / Wein closure: when the fixed
        # dense-phase carrying capacity would exceed the actual superficial
        # velocity, collapse toward the micro-bubbling boundary instead of
        # allowing u_d >= u0 to force a negative Eq.3.24 numerator.
        return float(u_mf + (u0 - u_mf) * 0.999)

    if name == "backsolve_visible_epsb":
        eps_b_target = float(calc_epsilon_b(u0, u_mf, u_b))
        denom = max(1.0 + eps_b_target * (n_b - 1.0), 1e-12)
        u_d = (u0 - eps_b_target * u_b) / denom
        return float(min(max(u_d, 0.0), min(u0, u_b)))

    raise ValueError(
        f"Unsupported u_d_closure={strategy!r}; "
        "expected 'current', 'umf_over_epsmf', 'hilligardt_eq311', "
        "'wein_1992_eq312', or 'backsolve_visible_epsb'"
    )


def calc_phase_exchange(
    *,
    K_bd: float,
    V_b: float,
    C_b: npt.ArrayLike,
    C_d: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """计算气泡相到悬浮相的组分交换。

    Ref: Hamel (1999) Eq. 2.2
    """
    C_b_arr = np.asarray(C_b, dtype=np.float64)
    C_d_arr = np.asarray(C_d, dtype=np.float64)
    assert C_b_arr.shape == (N_GAS,)
    assert C_d_arr.shape == (N_GAS,)
    return float(K_bd) * float(V_b) * (C_b_arr - C_d_arr)
