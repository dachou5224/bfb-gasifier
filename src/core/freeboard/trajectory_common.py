"""Wirsum/Hamel ODE helpers, bisect/event search.

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


def hamel_delta_pos_height_from_velocity(
    *,
    u: float,
    u0: float,
    h0: float,
    a: float,
    b: float,
    c: float,
) -> float | None:
    """Hamel Eq. (3.96)：Δ>0 分支高度–速度反函数。

    Source: ``docs/hamel_submodels/04_freeboard_particle_trajectory_analytic.md``; Hamel Eq. 3.96
    """
    delta = 4.0 * a * c - b * b
    if delta <= 1e-14 or abs(a) < 1e-14:
        return None
    omega = math.sqrt(delta)
    q = a * float(u) * float(u) + b * float(u) + c
    q0 = a * float(u0) * float(u0) + b * float(u0) + c
    if q <= 0.0 or q0 <= 0.0:
        return None
    term1 = -(1.0 / (2.0 * a)) * math.log(max(q / q0, 1e-300))
    arg_u = (2.0 * a * float(u) + b) / omega
    arg_u0 = (2.0 * a * float(u0) + b) / omega
    term2 = -(b / (a * omega)) * (math.atan(arg_u) - math.atan(arg_u0))
    h = float(h0 + term1 + term2)
    if not math.isfinite(h):
        return None
    return h


def calc_critical_diameter_eq397(
    *,
    mu_g: float,
    rho_g: float,
    rho_p: float,
    c_w: float,
) -> float | None:
    """Hamel Eq. (3.97) 颗粒极限直径 d_P* [m]（Δ>0 分支存在条件）。

    d_P* = (108 μ_g² / (c_W ρ_P ρ_g - ρ_g²))^(1/3)

    Source: ``docs/hamel_submodels/04_freeboard_particle_trajectory_analytic.md``; Hamel Eq. 3.97
    """
    denom = float(c_w) * float(rho_p) * float(rho_g) - float(rho_g) ** 2
    if denom <= 0.0:
        return None
    radicand = 108.0 * float(mu_g) ** 2 / denom
    if radicand <= 0.0:
        return None
    return float(radicand ** (1.0 / 3.0))


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

