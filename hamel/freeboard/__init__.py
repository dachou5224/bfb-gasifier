"""自由空间：床面抛出、气速、单颗粒轨迹与参考记账。"""

from hamel.freeboard._common import FreeboardError
from hamel.freeboard.accounting import (
    CellTransport,
    ReferenceAccounting,
    StableEquilibrium,
    TrajectoryFate,
    account_freeboard,
    stable_equilibria,
    trace_until_fate,
)
from hamel.freeboard.ejection import (
    approach_factor,
    coalescence_weight,
    ejection_flux,
    ejection_mass_flow,
    ejection_velocity_samples,
    pair_coalescence_probability,
    spacing_ratio,
)
from hamel.freeboard.motion import (
    closed_form_state,
    discriminant,
    drag_coefficient,
    drag_limit_z,
    integrate_freeboard,
    integrate_motion,
    motion_coefficients,
    particle_acceleration,
    terminal_settling_velocity,
)
from hamel.freeboard.transport import cross_section_m2, ghost_bubble_velocity, superficial_velocity

__all__ = [
    "CellTransport",
    "FreeboardError",
    "ReferenceAccounting",
    "StableEquilibrium",
    "TrajectoryFate",
    "account_freeboard",
    "approach_factor",
    "closed_form_state",
    "coalescence_weight",
    "cross_section_m2",
    "discriminant",
    "drag_coefficient",
    "drag_limit_z",
    "ejection_flux",
    "ejection_mass_flow",
    "ejection_velocity_samples",
    "ghost_bubble_velocity",
    "integrate_freeboard",
    "integrate_motion",
    "motion_coefficients",
    "pair_coalescence_probability",
    "particle_acceleration",
    "spacing_ratio",
    "stable_equilibria",
    "superficial_velocity",
    "terminal_settling_velocity",
    "trace_until_fate",
]
