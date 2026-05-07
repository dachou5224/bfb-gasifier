from __future__ import annotations

import numpy as np

from src.core.cell_balances import calc_gas_balance_residual
from src.core.phase_gate_checks import gate_20_elemental_closure


def test_conservation_elemental_closure_is_hard_gate() -> None:
    """无源码对照时，元素守恒是第一硬约束。"""
    result = gate_20_elemental_closure()
    assert result.passed, result.hard_failures
    assert result.metrics["pyrolysis"]["max_abs_closure"] <= 1e-8
    assert result.metrics["reaction"]["max_abs_closure"] <= 1e-7


def test_conservation_gas_exchange_is_strictly_internal() -> None:
    """相间交换在总气相上应内部抵消。"""
    n_ex = np.array([1.0, -2.0, 0.5] + [0.0] * 8, dtype=np.float64)
    zero = np.zeros_like(n_ex)
    res = calc_gas_balance_residual(
        N_zu_d=zero,
        N_rez_d=zero,
        N_d_in=zero,
        R_gas_d=zero,
        N_d=zero,
        N_ex=n_ex,
        N_zu_b=zero,
        N_rez_b=zero,
        N_b_in=zero,
        R_gas_b=zero,
        N_b=zero,
    )
    np.testing.assert_allclose(res[:11] + res[11:], 0.0)
