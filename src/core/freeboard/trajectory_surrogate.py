"""Surrogate / force-balance particle advances.

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

