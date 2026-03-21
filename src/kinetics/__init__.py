"""动力学子模型：Arrhenius 工厂函数、炭反应、气相反应、焦油反应。"""

from .arrhenius import k_hobbs, k_jensen_r7, k_standard
from .char_reactions import (
    D_d_A,
    D_A_ref,
    d_core_from_spm_char_conversion,
    phi_c,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4,
    rate_R4_effective,
    r1_hamel_kinetic_constants,
)
from .gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R9,
    wgsr_equilibrium_constant,
)
from .tar_reactions import rate_R10, rate_R11, rate_R11_bubble, rate_R11_suspension
