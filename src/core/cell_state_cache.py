"""Cell cache helpers.

Isolate cache match/store/invalidate logic from ``Cell`` to keep the class
focused on state and thin orchestration.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from src.core.constants import Rg

if TYPE_CHECKING:
    import numpy.typing as npt

    from src.core.cell import Cell


def thermo_cache_matches(cell: "Cell") -> bool:
    if not cell._thermo_cache_valid:
        return False
    return bool(
        float(cell.T) == float(cell._thermo_cache_T)
        and float(cell.P) == float(cell._thermo_cache_P)
        and np.array_equal(cell.N_b, cell._thermo_cache_N_b)
        and np.array_equal(cell.N_d, cell._thermo_cache_N_d)
    )


def get_local_thermo_bundle(cell: "Cell") -> dict[str, "npt.NDArray[np.float64] | float"]:
    if not thermo_cache_matches(cell):
        # 运行时经由 src.core.cell 取符号，保留对 tests 中 monkeypatch
        # ``src.core.cell.gas_diffusivity_correlation`` 的兼容性。
        from src.core import cell as cell_mod

        y_b = cell._mole_fractions("b")
        y_d = cell._mole_fractions("d")
        c_fac = cell.P / (Rg * cell.T)
        cell._thermo_cache_valid = True
        cell._thermo_cache_T = float(cell.T)
        cell._thermo_cache_P = float(cell.P)
        cell._thermo_cache_N_b[:] = cell.N_b
        cell._thermo_cache_N_d[:] = cell.N_d
        cell._thermo_cache_y_b[:] = y_b
        cell._thermo_cache_y_d[:] = y_d
        cell._thermo_cache_C_b[:] = y_b * c_fac
        cell._thermo_cache_C_d[:] = y_d * c_fac
        cell._thermo_cache_D_g = float(cell_mod.gas_diffusivity_correlation(cell.T, cell.P))
    return {
        "y_b": cell._thermo_cache_y_b,
        "y_d": cell._thermo_cache_y_d,
        "C_b": cell._thermo_cache_C_b,
        "C_d": cell._thermo_cache_C_d,
        "D_g": cell._thermo_cache_D_g,
    }


def hydrodynamics_cache_matches(cell: "Cell") -> bool:
    """Whether the cached hydrodynamics still matches the current local state."""
    if not cell._hydro_cache_valid:
        return False
    current_u0_target = np.nan if cell.u0_target is None else float(cell.u0_target)
    current_options = (
        str(cell.cell_type),
        str(cell.u_d_closure),
        str(cell.bubble_diameter_model),
        str(cell.psi_b_strategy),
        str(cell.lambda_strategy),
        str(cell.xi_strategy),
        str(cell.bubble_velocity_strategy),
        str(cell.bubble_ode_strategy),
        repr(cell.freeboard_u_bed_top),
        repr(cell.freeboard_d_b_bed_top),
        repr(cell.freeboard_height_from_bed),
        repr(cell.freeboard_u_gb_scale),
        repr(cell.dense_segment_eps_d_voidage if cell.cell_type != "freeboard" else cell.freeboard_eps_d_voidage),
    )
    same_u0_target = (
        (np.isnan(current_u0_target) and np.isnan(float(cell._hydro_cache_u0_target)))
        or current_u0_target == float(cell._hydro_cache_u0_target)
    )
    return bool(
        float(cell.T) == float(cell._hydro_cache_T)
        and float(cell.P) == float(cell._hydro_cache_P)
        and same_u0_target
        and current_options == cell._hydro_cache_options
        and np.array_equal(cell.N_b, cell._hydro_cache_N_b)
        and np.array_equal(cell.N_d, cell._hydro_cache_N_d)
    )


def store_hydrodynamics_cache(cell: "Cell") -> None:
    current_u0_target = np.nan if cell.u0_target is None else float(cell.u0_target)
    cell._hydro_cache_valid = True
    cell._hydro_cache_T = float(cell.T)
    cell._hydro_cache_P = float(cell.P)
    cell._hydro_cache_u0_target = current_u0_target
    cell._hydro_cache_N_b[:] = cell.N_b
    cell._hydro_cache_N_d[:] = cell.N_d
    cell._hydro_cache_options = (
        str(cell.cell_type),
        str(cell.u_d_closure),
        str(cell.bubble_diameter_model),
        str(cell.psi_b_strategy),
        str(cell.lambda_strategy),
        str(cell.xi_strategy),
        str(cell.bubble_velocity_strategy),
        str(cell.bubble_ode_strategy),
        repr(cell.freeboard_u_bed_top),
        repr(cell.freeboard_d_b_bed_top),
        repr(cell.freeboard_height_from_bed),
        repr(cell.freeboard_u_gb_scale),
        repr(cell.dense_segment_eps_d_voidage if cell.cell_type != "freeboard" else cell.freeboard_eps_d_voidage),
    )


def snapshot_vorabrechnung_hydrodynamics(cell: "Cell") -> None:
    """Freeze current hydrodynamics for use inside the inner global-NR loop."""
    cell._vorab_hydro_cache["u_mf"] = float(cell.u_mf)
    cell._vorab_hydro_cache["u_b"] = float(cell.u_b)
    cell._vorab_hydro_cache["d_b"] = float(cell.d_b)
    cell._vorab_hydro_cache["eps_b"] = float(cell.eps_b)
    cell._vorab_hydro_cache["eps_d"] = float(cell.eps_d)
    cell._vorab_hydro_cache["eps_d_voidage"] = float(cell.eps_d_voidage)
    cell._vorab_hydro_cache["K_bd"] = float(cell.K_bd)
    cell._vorab_hydro_cache["V_b"] = float(cell.V_b)
    cell._vorab_hydro_cache["V_d"] = float(cell.V_d)
    cell._vorab_hydro_cache["u0"] = float(cell.u0)
    cell._vorab_hydro_cache["u_d"] = float(cell.u_d)
    cell._vorab_hydro_cache["n_rz"] = float(cell.n_rz)
    cell._vorab_hydro_cache_valid = True
    cell._vorab_K_solid_auf[:, :] = cell.K_solid_auf
    cell._vorab_K_solid_ab[:, :] = cell.K_solid_ab
    cell._vorab_solid_transport_cache_valid = True


def apply_frozen_vorabrechnung_hydrodynamics(cell: "Cell") -> None:
    assert cell._vorab_hydro_cache_valid, "Frozen Vorabrechnung hydrodynamics must be prepared before inner NR"
    cache = cell._vorab_hydro_cache
    cell.u_mf = float(cache["u_mf"])
    cell.u_b = float(cache["u_b"])
    cell.d_b = float(cache["d_b"])
    cell.eps_b = float(cache["eps_b"])
    cell.eps_d = float(cache["eps_d"])
    cell.eps_d_voidage = float(cache["eps_d_voidage"])
    cell.K_bd = float(cache["K_bd"])
    cell.V_b = float(cache["V_b"])
    cell.V_d = float(cache["V_d"])
    cell.u0 = float(cache["u0"])
    cell.u_d = float(cache["u_d"])
    cell.n_rz = float(cache["n_rz"])
    if cell._vorab_solid_transport_cache_valid:
        cell.K_solid_auf[:, :] = cell._vorab_K_solid_auf
        cell.K_solid_ab[:, :] = cell._vorab_K_solid_ab
    cell._hydro_cache_valid = False


def enable_inner_nr_vorabrechnung_freeze(cell: "Cell") -> None:
    cell._freeze_vorabrechnung_inner_nr = True
    cell._hydro_cache_valid = False


def disable_inner_nr_vorabrechnung_freeze(cell: "Cell") -> None:
    cell._freeze_vorabrechnung_inner_nr = False
    cell._hydro_cache_valid = False


def invalidate_vorabrechnung_cache(cell: "Cell") -> None:
    """Explicitly invalidate Vorabrechnung caches for outer-loop rebuild."""
    cell._vm_cache_valid = False
    cell._vorab_hydro_cache_valid = False
    cell._vorab_solid_transport_cache_valid = False
    cell._freeze_vorabrechnung_inner_nr = False
    cell._vm_cache_T = np.nan
    cell._vm_cache_tau = np.nan
    cell._vm_cache_T_init = np.nan
    cell._vm_cache_m_vm_in = np.nan
    cell._vm_cache_m_moist_in = np.nan

