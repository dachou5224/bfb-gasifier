"""物性、组分焓与附录 A1 的气相平衡。"""

from hamel.properties.data import (
    GAS_SPECIES, EquilibriumError, PropertyError, PropertyRangeError,
    gas_constant_j_mol_k, molar_mass_kg_mol, stefan_boltzmann_w_m2_k4,
)
from hamel.properties.equilibrium import (
    adiabatic_gas_equilibrium, balance_deviation, gibbs_minimization,
    initial_lambdas, is_converged,
)
from hamel.properties.thermo import (
    ash_specific_heat_j_per_kg_k, char_specific_enthalpy_j_per_kg,
    chemical_potential_j_mol, density_kg_m3, enthalpy_j_mol, entropy_j_mol_k,
    formation_enthalpy_j_mol, fuel_formation_enthalpy_j_per_kg,
    heat_capacity_j_mol_k, vaporization_enthalpy_j_mol, volume_flow_m3_s,
    water_free_to_ash_free,
)
from hamel.properties.transport import (
    binary_diffusivity_m2_s, chung_conductivity_w_m_k, collision_integral,
    fuller_binary_diffusivity_m2_s, fuller_volume, interphase_diffusivity_m2_s,
    lennard_jones, mixture_thermal_conductivity_w_m_k, mixture_viscosity_pa_s,
    thermal_conductivity_w_m_k, viscosity_pa_s,
)
