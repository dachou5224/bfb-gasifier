"""Phase 1/2 LU 验证用例 shared ``Reactor.solve`` kwargs。"""

from __future__ import annotations

from typing import Any

PHASE1_HTW_LU_SOLVE_KWARGS: dict[str, Any] = {
    "max_global_iter": 150,
    "tol_global": 1.0,
    "solver": "global_nr",
    "nr_init_strategy": "vorabrechnung",
    "nr_jacobian_strategy": None,
}

PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS: dict[str, Any] = dict(PHASE1_HTW_LU_SOLVE_KWARGS)

PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS: dict[str, Any] = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS["nr_jacobian_strategy"] = None
