"""Freeboard axial segment package (OPT-003 split).

Public API preserved via ``src.core.freeboard_segment`` re-export shim.
"""

from src.core.freeboard.types import FreeboardBackground, FreeboardState
from src.core.freeboard.trajectory_common import (
    _bisect_event_time,
    _find_first_event_time,
    _find_first_event_time_with_growth,
    _quadratic_ode_step,
    _quadratic_ode_state,
    _wirsum_abc,
    _wirsum_abc_exact_audit,
)
from src.core.freeboard.transport import (
    _aggregate_class_transport,
    _apply_char_sink_to_class_transport,
    _build_size_class_bundles,
    _build_velocity_samples,
    _char_area_and_representative_diameter,
    _empty_class_transport,
    _freeboard_inventory_char_sink_diag,
    _gas_concentrations,
    _representative_char_conversion,
    _resolve_transport_class_mapping,
    _solve_T_from_enthalpy,
    _wet_gas_profile,
)
from src.core.freeboard.reaction import (
    _build_freeboard_background_profile,
    _estimate_rho_cat_profile_from_projected_transport,
    _freeboard_reaction_step,
)
from src.core.freeboard.trajectory_surrogate import (
    _advance_particle_samples_force_balance,
    _advance_particle_samples_surrogate,
)
from src.core.freeboard.trajectory_wirsum import _advance_particle_samples_analytical_wirsum
from src.core.freeboard.trajectory_hamel import (
    _advance_particle_samples_global_projected_hamel,
    _simulate_freeboard_global_projected_hamel,
)
from src.core.freeboard.segment import simulate_freeboard

__all__ = [
    "FreeboardBackground",
    "FreeboardState",
    "simulate_freeboard",
    "_advance_particle_samples_analytical_wirsum",
    "_advance_particle_samples_force_balance",
    "_advance_particle_samples_global_projected_hamel",
    "_advance_particle_samples_surrogate",
    "_aggregate_class_transport",
    "_apply_char_sink_to_class_transport",
    "_bisect_event_time",
    "_build_freeboard_background_profile",
    "_build_size_class_bundles",
    "_build_velocity_samples",
    "_char_area_and_representative_diameter",
    "_empty_class_transport",
    "_estimate_rho_cat_profile_from_projected_transport",
    "_find_first_event_time",
    "_find_first_event_time_with_growth",
    "_freeboard_inventory_char_sink_diag",
    "_freeboard_reaction_step",
    "_gas_concentrations",
    "_quadratic_ode_step",
    "_quadratic_ode_state",
    "_representative_char_conversion",
    "_resolve_transport_class_mapping",
    "_simulate_freeboard_global_projected_hamel",
    "_solve_T_from_enthalpy",
    "_wet_gas_profile",
    "_wirsum_abc",
    "_wirsum_abc_exact_audit",
]
