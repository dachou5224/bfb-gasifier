"""F_E6：在 F_E5 数学根上把焦油出口升为未知量。不经 ``Reactor``，不切生产。

未知量 43：物种 39 + ``M_ash`` + ``N_TAR1`` + ``N_TAR2`` + ``T``。
焦油行是出口闭合 ``N_tar = n_py,tar(T)``，不是 R10 氧化。
能量行 Hamel Eq. 2.8 的焦油焓改用 ``x`` 中的 ``N_tar``。

从 F_E5 根嵌入：用 ``outer_P_named`` 的 DAEM 焦油作 Startwert。
``Q_W=0`` 时嵌入点已是根（焦油行 0，能量与 E5 重合）。

雅可比：N/M 用 F_E5 块；焦油列解析；T 列对 ``F(x;P(x))`` 差商。

Source: Hamel (1999) Eq. 2.7–2.8; Wirsum (1997) §2.3.3
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.species import enthalpy_molar
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    coupled_residual_named,
    energy_row_named,
    freeze_stage_q,
    n_unknowns,
    outer_P_named,
    pack_named,
    unpack_named,
)
from src.solvers.wirsum_1997_named_f_e4 import residual_named_e4
from src.solvers.wirsum_1997_named_f_e5 import freeze_e5_q, jacobian_named_e5, solve_named_e5_from_e4
from src.solvers.wirsum_1997_named_f_r12 import named_r12_scaled_rank
from src.solvers.wirsum_1997_named_f_r12_r5_r8_r6 import named_r12_r5_r8_r6_newton_step
from src.solvers.wirsum_1997_named_f_r234 import HOMOG_ON, R234_ON
from src.solvers.wirsum_1997_named_f_registry import f_name
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

_ENABLED: frozenset[str] = HOMOG_ON | R234_ON
_IDX_TAR1: int = 40
_IDX_TAR2: int = 41
_T_FD: float = 1.0  # [K]


def freeze_e6_q(*, Q_W: float = 0.0) -> NamedStageQ:
    return freeze_stage_q(6, Q_W=float(Q_W), enabled_rxns=_ENABLED)


def residual_named_e6(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """``F(x; P_DAEM(x); Q_W)``，含焦油出口行。"""
    return residual_named_e4(x, nq6)


def _x41_from_e6(x43: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], float, float]:
    z = np.asarray(x43, dtype=np.float64)
    return np.concatenate([z[:40], z[-1:]]), float(z[_IDX_TAR1]), float(z[_IDX_TAR2])


def _lift_j41(J41: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    j = np.asarray(J41, dtype=np.float64)
    J = np.zeros((43, 43), dtype=np.float64)
    J[:40, :40] = j[:40, :40]
    J[:40, -1] = j[:40, -1]
    J[-1, :40] = j[-1, :40]
    J[-1, -1] = j[-1, -1]
    return J


def embed_e5_root_in_e6(
    x41: npt.NDArray[np.float64],
    nq6: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """41 维 F_E5 根嵌入 F_E6。焦油用当前 T 的 DAEM ``n_py,tar``。"""
    z = np.asarray(x41, dtype=np.float64)
    nqP = outer_P_named(z, freeze_e5_q(Q_W=nq6.Q_W))
    return pack_named(
        z[:39],
        float(z[-1]),
        6,
        M_ash=float(z[39]),
        N_tar1=float(nqP.energy.n_tar.get("TAR1", 0.0)),
        N_tar2=float(nqP.energy.n_tar.get("TAR2", 0.0)),
    )


def jacobian_named_e6(
    x: npt.NDArray[np.float64],
    nq6: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """N/M 用 E5 块；焦油列 ``I`` 与 ``-h_tar``；T 列含 DAEM。"""
    z = np.asarray(x, dtype=np.float64)
    x41, _n1, _n2 = _x41_from_e6(z)
    nq5 = freeze_e5_q(Q_W=nq6.Q_W)
    J = _lift_j41(jacobian_named_e5(x41, nq5))
    T = float(z[-1])
    J[_IDX_TAR1, _IDX_TAR1] = 1.0
    J[_IDX_TAR2, _IDX_TAR2] = 1.0
    J[-1, _IDX_TAR1] = -float(enthalpy_molar("TAR1", T))
    J[-1, _IDX_TAR2] = -float(enthalpy_molar("TAR2", T))
    xp = z.copy()
    xp[-1] = T + _T_FD
    dF = (residual_named_e6(xp, nq6) - residual_named_e6(z, nq6)) / _T_FD
    J[:, -1] = dF
    return J


def _in_domain(z: npt.NDArray[np.float64]) -> bool:
    return bool(np.min(z[:-1]) >= 0.0 and float(z[-1]) > 0.0)


def e6_source_snapshot(
    x41: npt.NDArray[np.float64],
    nq5: NamedStageQ,
    nq6: NamedStageQ,
) -> dict[str, Any]:
    z41 = np.asarray(x41, dtype=np.float64)
    x43 = embed_e5_root_in_e6(z41, nq6)
    f5 = coupled_residual_named(z41, outer_P_named(z41, nq5), rate_scale=1.0)
    f6 = residual_named_e6(x43, nq6)
    u = unpack_named(x43, 6)
    return {
        "F5_l2": float(np.linalg.norm(f5)),
        "F6_l2": float(np.linalg.norm(f6)),
        "F6_tar1": float(f6[_IDX_TAR1]),
        "F6_tar2": float(f6[_IDX_TAR2]),
        "F6_energy": float(f6[-1]),
        "d_energy": float(f6[-1] - f5[-1]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "T_K": float(x43[-1]),
    }


def solve_named_e6_from_e5(
    *,
    Q_W: float = 0.0,
    max_iter: int = 20,
    e5: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """从 F_E5 根做 F_E6。默认 ``Q_W=0``。"""
    nq5 = freeze_e5_q(Q_W=Q_W)
    nq6 = freeze_e6_q(Q_W=Q_W)
    if e5 is None:
        e5 = solve_named_e5_from_e4(Q_W=Q_W, max_iter=max_iter)
    x0 = embed_e5_root_in_e6(np.asarray(e5["x"], dtype=np.float64), nq6)
    snap0 = e6_source_snapshot(np.asarray(e5["x"], dtype=np.float64), nq5, nq6)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return residual_named_e6(z, nq6)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return jacobian_named_e6(z, nq6)

    last = wirsum_1997_strict_newton(
        res, jac, x0, max_iter=max_iter, trial_in_domain=_in_domain,
        linear_solve_fn=named_r12_r5_r8_r6_newton_step,
    )
    stop = classify_newton_stop(last, model_is_full_rate=True)
    info = named_r12_scaled_rank(jacobian_named_e6(last.x, nq6), last.F, last.x)
    nqP = outer_P_named(last.x, nq6)
    u = unpack_named(last.x, 6)
    t = float(last.x[-1])
    return {
        **stop,
        "x": np.asarray(last.x, dtype=np.float64),
        "F": np.asarray(last.F, dtype=np.float64),
        "rank_last": int(last.rank_last),
        "raw_rank": int(info["raw_rank"]),
        "equilibrated_rank": int(info["equilibrated_rank"]),
        "T_K": t,
        "T_e5_K": float(e5["T_K"]),
        "M_ash": float(u["M_ash"]),
        "N_tar1": float(u["N_tar1"]),
        "N_tar2": float(u["N_tar2"]),
        "Q_W": float(nq6.Q_W),
        "Q_W_is_fitted_to_1100K": False,
        "energy_W": float(energy_row_named(last.x, nqP)),
        "f_name": f_name(6),
        "n_unknowns": n_unknowns(6),
        "enabled_rxns": tuple(sorted(_ENABLED)),
        "e5_physical_convergence": bool(e5["physical_convergence"]),
        "at_e5_root": snap0,
        "rate_formula_changed": False,
        "scaled_linear_solve": True,
        "not_a_furnace_convergence_claim": True,
        "rewrites_unknown_N": False,
        "production_hamel_outer_refresh": False,
        "r10_tar_oxidation": False,
    }
