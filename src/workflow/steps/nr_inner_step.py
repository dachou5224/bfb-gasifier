"""流程图：Zellenmodell / Newton–Raphson（Check1 内层）。

封装 ``solve_global_nr`` 调用与 Vorabrechnung 冻结开关生命周期。
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from src.solvers.convergence import InnerConvergence

if TYPE_CHECKING:
    from src.core.cell import Cell
    from src.core.reactor import Reactor, ReactorConfig


def build_global_nr_inner_solve_fn(
    reactor: "Reactor",
    solver_cells: list["Cell"],
    cfg: "ReactorConfig",
    tol_rms: float,
    inner_check1: InnerConvergence,
    verbose: bool,
) -> Callable[[int, float, str, int], dict]:
    """构造 ``run_global_nr_outer_abgleich(..., solve_inner_fn=...)`` 所需的内层求解闭包。"""

    def _solve_inner(inner_iter_cap: int, lambda_seed_outer: float, j_mode: str, j_lag: int) -> dict:
        for cell in solver_cells:
            if not bool(getattr(cell, "_vorab_hydro_cache_valid", False)):
                raise RuntimeError("Vorabrechnung hydrodynamics freeze cache is not initialized for inner NR.")
        for cell in solver_cells:
            cell.enable_inner_nr_vorabrechnung_freeze()
        try:
            # 运行时从模块取符号，便于测试 monkeypatch ``src.solvers.global_nr_solver.solve_global_nr``
            from src.solvers import global_nr_solver as _gnr

            return _gnr.solve_global_nr(
                cells=solver_cells,
                apply_bc_fn=reactor._apply_all_bc_for_nr,
                apply_local_bc_fn=reactor._apply_local_bc_for_nr,
                affected_residual_cells_fn=reactor._affected_nr_residual_cells,
                ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
                ref_solid_kg_s=cfg.fuel_feed,
                ref_energy_W=cfg.fuel_feed * 20e6,
                max_iter=inner_iter_cap,
                tol_rms=tol_rms,
                inner_convergence=inner_check1,
                lambda_init=float(np.clip(lambda_seed_outer, 1.0 / 1024.0, 1.0)),
                jacobian_strategy=j_mode,
                linear_solver_backend=(
                    "structured_direct"
                    if str(j_mode) in {"block_tridiag_structured", "band_plus_side_elements_structured"}
                    else "sparse_direct_fallback"
                ),
                jacobian_lag_steps=j_lag,
                verbose=verbose,
            )
        finally:
            for cell in solver_cells:
                cell.disable_inner_nr_vorabrechnung_freeze()

    return _solve_inner
