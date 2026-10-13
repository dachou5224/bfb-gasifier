"""Cell residual assembly pipeline.

Keep the residual stage ordering outside ``Cell`` so the Hamel-style cell model
sequence can be audited and evolved without growing ``cell.py`` further.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from src.core.cell_balances import assemble_cell_residual_vector

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt

    from src.core.cell import Cell

_FULL_RESIDUAL_GROUPS = frozenset({"gas", "solid", "energy"})


def _run_pre_balance_stages(cell: "Cell", *, rate_multiplier: float) -> None:
    if cell._freeze_vorabrechnung_inner_nr and cell._vorab_hydro_cache_valid:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    else:
        cell.calc_hydrodynamics()
    cell.calc_exchange()
    cell.calc_reactions(rate_multiplier=rate_multiplier)


def execute_cell_residual_pipeline(
    cell: "Cell",
    *,
    rate_multiplier: float = 1.0,
    groups: frozenset[str] | None = None,
    copy: bool = True,
) -> "npt.NDArray[np.float64]":
    """Run the cell residual pipeline in the established physical order.

    Order:
    1. Hydrodynamics or frozen Vorabrechnung hydrodynamics
    2. Inter-phase exchange
    3. Reaction/source update
    4. Gas / solid / energy balances (subset when ``groups`` is set)
    5. Residual vector assembly (full vector only for the default groups)

    When ``groups`` is a strict subset, skip unused balance evaluations.
    Ref: OPT-012 preprojection trial loops (connectivity).
    """
    selected = _FULL_RESIDUAL_GROUPS if groups is None else frozenset(groups)
    if not selected:
        raise ValueError("groups must be non-empty")
    unknown = selected - _FULL_RESIDUAL_GROUPS
    if unknown:
        raise ValueError(f"unknown residual groups: {sorted(unknown)}")

    _run_pre_balance_stages(cell, rate_multiplier=float(rate_multiplier))

    need_gas = "gas" in selected
    need_solid = "solid" in selected
    need_energy = "energy" in selected

    if selected == frozenset({"solid"}):
        res_solid = cell.calc_solid_balance().flatten()
        if copy:
            import numpy as np

            return np.asarray(res_solid, dtype=np.float64).copy()
        n_gas = int(cell.N_d.shape[0])
        solid_off = 2 * n_gas
        cell._work_res[solid_off : solid_off + res_solid.size] = res_solid
        return cell._work_res[solid_off : solid_off + res_solid.size]

    if selected == frozenset({"energy"}):
        res_energy = float(cell.calc_energy_balance())
        import numpy as np

        if copy:
            return np.array([res_energy], dtype=np.float64)
        cell._work_res[-1] = res_energy
        return cell._work_res[-1:]

    res_gas = cell.calc_gas_balance() if need_gas else None
    res_solid = cell.calc_solid_balance().flatten() if need_solid else None
    res_energy = float(cell.calc_energy_balance()) if need_energy else 0.0

    if selected == _FULL_RESIDUAL_GROUPS:
        assemble_cell_residual_vector(
            res_gas=res_gas,
            res_solid=res_solid,
            res_energy=res_energy,
            out=cell._work_res,
        )
        if copy:
            return cell._work_res.copy()
        return cell._work_res

    import numpy as np

    if need_gas:
        gas_arr = np.asarray(res_gas, dtype=np.float64)
        cell._work_res[: gas_arr.size] = gas_arr
    if need_solid:
        n_gas = int(cell.N_d.shape[0])
        solid_off = 2 * n_gas
        solid_arr = np.asarray(res_solid, dtype=np.float64).reshape(-1)
        cell._work_res[solid_off : solid_off + solid_arr.size] = solid_arr
    if need_energy:
        cell._work_res[-1] = res_energy

    if copy:
        if selected == frozenset({"gas"}):
            return np.asarray(res_gas, dtype=np.float64).copy()
        if selected == frozenset({"gas", "energy"}):
            gas_arr = np.asarray(res_gas, dtype=np.float64)
            return np.concatenate((gas_arr, np.array([res_energy], dtype=np.float64)))
        parts: list[np.ndarray] = []
        if need_gas:
            parts.append(np.asarray(res_gas, dtype=np.float64))
        if need_solid:
            parts.append(np.asarray(res_solid, dtype=np.float64).reshape(-1))
        if need_energy:
            parts.append(np.array([res_energy], dtype=np.float64))
        return np.concatenate(parts)

    if selected == frozenset({"gas"}):
        n_gas = int(cell.N_d.shape[0])
        return cell._work_res[: 2 * n_gas]
    if selected == frozenset({"gas", "energy"}):
        return cell._work_res
    return cell._work_res
