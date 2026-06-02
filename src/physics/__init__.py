"""物理子模型：流体力学、相分率、自由板区。"""

from .bubble_dynamics import (
    bubble_interaction_factor,
    bubble_diameter_ode,
    bubble_rise_velocity,
    classify_bubble_regime,
    initial_bubble_diameter,
    integrate_bubble_diameter,
    single_bubble_velocity,
    bubble_lifetime,
    resolve_xi_b,
)
from .freeboard import calc_beta_a, calc_cd_haider, calc_terminal_velocity_haider, calc_u_gb
from .mass_transfer import calc_kbd, calc_u_br
from .minimum_fluidization import calc_archimedes, calc_re_mf, calc_u_mf, compute_u_mf
from .phase_fractions import calc_bulk_solid_holdup, calc_epsilon_b, calc_epsilon_d, calc_n_rz

__all__ = [
    "bubble_diameter_ode",
    "bubble_interaction_factor",
    "bubble_rise_velocity",
    "classify_bubble_regime",
    "initial_bubble_diameter",
    "integrate_bubble_diameter",
    "single_bubble_velocity",
    "bubble_lifetime",
    "resolve_xi_b",
    "calc_beta_a",
    "calc_cd_haider",
    "calc_terminal_velocity_haider",
    "calc_u_gb",
    "calc_kbd",
    "calc_u_br",
    "calc_archimedes",
    "calc_re_mf",
    "calc_u_mf",
    "compute_u_mf",
    "calc_bulk_solid_holdup",
    "calc_epsilon_b",
    "calc_epsilon_d",
    "calc_n_rz",
]
