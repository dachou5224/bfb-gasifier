"""Q5 16 维能量核诊断。输出不是物理收敛通过标记。"""
from __future__ import annotations

import json

import numpy as np

from src.solvers.wirsum_1997_stiffness_full_gas import (
    default_full_gas_inlets,
    energy_jacobian_quality_at,
)

N_d, N_b = default_full_gas_inlets()
M = 0.010 / 0.080  # [kg]
T = 1100.0  # [K]
x = np.concatenate([N_d, N_b, np.array([M, T], dtype=np.float64)])
payload = {
    "note": "diagnostic_not_physical_convergence",
    "T_retained": True,
    "fictitious_Q_W": False,
    "rate_scale_0": energy_jacobian_quality_at(x, rate_scale=0.0),
    "rate_scale_1": energy_jacobian_quality_at(x, rate_scale=1.0),
}
print(json.dumps(payload, indent=2))
