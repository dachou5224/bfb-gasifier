"""F_E5：在 F_E4 数学根上打开壁面热 ``Q_W``。不经 ``Reactor``，不切生产。

``Q_W`` 是具名 ``q`` [W]，进入 Hamel Eq. 2.7 能量行 ``H_in − H_out − Q_W``。
禁止用 ``Q_W`` 反推 1100 K。未知量仍是 41；DAEM 外层保持 ``outer_P_named``。

``Q_W = 0`` 时残差与 F_E4 同定义，E4 根即 E5 根。
``Q_W > 0`` 是离开本格的壁面热流，能量行平移 ``−Q_W``，温度下降。
这不是 HTW 全炉 ``Q_wall``，也不是拟合炉温。

雅可比：``Q_W`` 不依赖 ``x``，复用 F_E4 块。

Source: Hamel (1999) Eq. 2.7–2.8; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_e4 import (
    freeze_e4_q,
    jacobian_named_e4,
    residual_named_e4,
    solve_named_e4_from_e3,
)
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import named_r12_r5_r8_r6_newton_step
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_named_f_registry import f_name
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON

# 具名杠杆示例 [W]。与 stages 测试同一量级，不是拟合 1100 K 的热损。
NAMED_Q_W_EXAMPLE_W: float = 1000.0


def freeze_e5_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(5, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def residual_named_e5(
    x: npt.NDArray[np.float64],
    nq5: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """``F(x; P_DAEM(x); Q_W)``。``Q_W`` 只进能量行。"""
    return residual_named_e4(x, nq5)


def jacobian_named_e5(
    x: npt.NDArray[np.float64],
    nq5: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """``Q_W`` 对 ``x`` 的导数为零，T 列仍含 DAEM。"""
    return jacobian_named_e4(x, nq5)


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def e5_source_snapshot(
    x: npt.NDArray[np.float64],
    nq4: NamedStageQ,
    nq5: NamedStageQ,
) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    f4 = residual_named_e4(z, nq4)
    f5 = residual_named_e5(z, nq5)
    dF = f5 - f4
    return {
        "F4_l2": float(np.linalg.norm(f4)),
        "F5_l2": float(np.linalg.norm(f5)),
        "dF_l2": float(np.linalg.norm(dF)),
        "dF_energy": float(dF[-1]),
        "Q_W": float(nq5.Q_W),
        "T_K": float(z[-1]),
        "energy_e4_W": float(f4[-1]),
        "energy_e5_W": float(f5[-1]),
    }


def solve_named_e5_from_e4(
    *,
    Q_W: float = 0.0,
    max_iter: int = 20,
    e4: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E4 根做 F_E5。``Q_W`` 是具名 q，默认 0。"""
    nq4 = freeze_e4_q()
    nq5 = freeze_e5_q(Q_W=Q_W)
    if e4 is None:
        e4 = solve_named_e4_from_e3(max_iter=max_iter)
    x0 = np.asarray(e4["x"], dtype=np.float64).copy()
    snap0 = e5_source_snapshot(x0, nq4, nq5)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_named_e5(z, nq5)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_named_e5(z, nq5)

    last = wirsum_1997_strict_newton(
        res, jac, x0, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_r8_r6_newton_step,
    )
    stop = classify_newton_stop(last, model_is_full_rate=True)
    info = named_r12_scaled_rank(jacobian_named_e5(last.x, nq5), last.F, last.x)
    nqP = outer_P_named(last.x, nq5)
    u = unpack_named(last.x, 5)
    t = float(last.x[-1])
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": t,
        "T_e4_K": float(e4["T_K"]),
        "M_ash": float(u["M_ash"]),
        "Q_W": float(nq5.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": f_name(5),
        "n_unknowns": n_unknowns(5),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "e4_physical_convergence": bool(e4["physical_convergence"]),
        "at_e4_root": snap0,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
    }
