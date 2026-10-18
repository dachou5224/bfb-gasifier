"""第 5 章反应速率。给定调用方传入的体积或质量，返回有符号进度速率。"""

from hamel.kinetics.arrhenius import arrhenius_rate_constant
from hamel.kinetics.common import KineticsInputError, ReactionRate
from hamel.kinetics.heterogeneous import (
    mechanism_factor, particle_temperature, r1_combustion_rate,
    r2_steam_gasification_rate, r3_boudouard_rate, r4_hydrogasification_rate,
)
from hamel.kinetics.homogeneous import (
    r5_co_oxidation_rate, r6_hydrogen_oxidation_rate, r7_methane_oxidation_rate,
    r8_equilibrium_constant, r8_water_gas_shift_rate, r9_methane_reforming_rate,
    r10_tar_oxidation_rate, r11_corella_rate, r11_serio_rate,
)
from hamel.kinetics.stoichiometry import reaction_stoichiometry

__all__ = [
    "KineticsInputError", "ReactionRate", "arrhenius_rate_constant",
    "mechanism_factor", "particle_temperature", "reaction_stoichiometry",
    "r1_combustion_rate", "r2_steam_gasification_rate", "r3_boudouard_rate",
    "r4_hydrogasification_rate", "r5_co_oxidation_rate", "r6_hydrogen_oxidation_rate",
    "r7_methane_oxidation_rate", "r8_equilibrium_constant", "r8_water_gas_shift_rate",
    "r9_methane_reforming_rate", "r10_tar_oxidation_rate", "r11_corella_rate",
    "r11_serio_rate",
]
