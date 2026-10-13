"""FreeboardState / FreeboardBackground dataclasses.

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

