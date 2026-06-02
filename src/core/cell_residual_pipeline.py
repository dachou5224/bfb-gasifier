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


def execute_cell_residual_pipeline(
    cell: "Cell",
    *,
    rate_multiplier: float = 1.0,
) -> "npt.NDArray[np.float64]":
    """Run the cell residual pipeline in the established physical order.

    Order:
    1. Hydrodynamics or frozen Vorabrechnung hydrodynamics
    2. Inter-phase exchange
    3. Reaction/source update
    4. Gas / solid / energy balances
    5. Residual vector assembly
    """
    if cell._freeze_vorabrechnung_inner_nr and cell._vorab_hydro_cache_valid:
        cell._apply_frozen_vorabrechnung_hydrodynamics()
    else:
        cell.calc_hydrodynamics()

    cell.calc_exchange()
    cell.calc_reactions(rate_multiplier=rate_multiplier)

    res_gas = cell.calc_gas_balance()
    res_solid = cell.calc_solid_balance().flatten()
    res_energy = cell.calc_energy_balance()
    assemble_cell_residual_vector(
        res_gas=res_gas,
        res_solid=res_solid,
        res_energy=res_energy,
        out=cell._work_res,
    )
    return cell._work_res.copy()

