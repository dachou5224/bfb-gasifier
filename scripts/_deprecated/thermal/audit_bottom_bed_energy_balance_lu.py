"""Audit bottom-bed energy balance terms for the LU Phase2 reactor.

The script is intentionally diagnostic only: it builds the standard Phase2 LU
state, applies NR boundary data, and decomposes the lower bed energy residual
into gas/solid enthalpy streams.  It also reports the temperature root for the
current outlet streams and for a "gas-closed" outlet where the current gas
balance residual is added to the outlet gas flow.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.core.cell import Cell
from src.core.cell_balances import calc_gas_enthalpy_flow, calc_solid_enthalpy_flow
from src.core.reactor import Reactor
from src.core.species import cp_ash, cp_char, cp_sand
from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr
from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config


@dataclass(frozen=True)
class EnergyTerms:
    hin: float
    hout: float
    residual: float
    terms: tuple[tuple[str, float, float, float], ...]


def _h_gas(cell: Cell, flow: np.ndarray, temperature: float) -> float:
    return calc_gas_enthalpy_flow(np.asarray(flow, dtype=np.float64), float(temperature), cell._h_cache)


def _h_solid(cell: Cell, flow: np.ndarray, temperature: float) -> float:
    return calc_solid_enthalpy_flow(
        np.asarray(flow, dtype=np.float64),
        float(temperature),
        ash_dry_wt=cell.solid.ash_dry_wt,
        VM_daf=cell.solid.VM_daf,
        h_f_dry=cell.solid.h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
    )


def _inlet_terms(cell: Cell) -> tuple[tuple[str, float, float, float], ...]:
    terms = [
        ("gas_axial_in", _h_gas(cell, cell.N_b_in + cell.N_d_in, cell.T_in_gas), float(np.sum(cell.N_b_in + cell.N_d_in)), float(cell.T_in_gas)),
        ("gas_fresh_zu", _h_gas(cell, cell.N_zu_b + cell.N_zu_d, cell.T_zu_gas), float(np.sum(cell.N_zu_b + cell.N_zu_d)), float(cell.T_zu_gas)),
        ("gas_recycle", _h_gas(cell, cell.N_rez_b + cell.N_rez_d, cell.T_rez_gas), float(np.sum(cell.N_rez_b + cell.N_rez_d)), float(cell.T_rez_gas)),
        ("solid_recycle", _h_solid(cell, cell.m_solid_rez, cell.T_rez_solid), float(np.sum(cell.m_solid_rez)), float(cell.T_rez_solid)),
        ("solid_axial_in", _h_solid(cell, cell.m_solid_in, cell.T_in_solid), float(np.sum(cell.m_solid_in)), float(cell.T_in_solid)),
        ("solid_up_in", _h_solid(cell, cell.m_solid_auf_in, cell.T_solid_auf_in), float(np.sum(cell.m_solid_auf_in)), float(cell.T_solid_auf_in)),
        ("solid_down_in", _h_solid(cell, cell.m_solid_ab_in, cell.T_solid_ab_in), float(np.sum(cell.m_solid_ab_in)), float(cell.T_solid_ab_in)),
        ("solid_fresh_zu", _h_solid(cell, cell.m_solid_zu, cell.T_zu_solid), float(np.sum(cell.m_solid_zu)), float(cell.T_zu_solid)),
    ]
    return tuple(terms)


def energy_terms(cell: Cell, *, gas_out: np.ndarray | None = None, solid_out: np.ndarray | None = None, temperature: float | None = None) -> EnergyTerms:
    T = float(cell.T if temperature is None else temperature)
    gas = np.asarray(cell.N_b + cell.N_d if gas_out is None else gas_out, dtype=np.float64)
    solid = np.asarray(cell._solid_outflow_rates() if solid_out is None else solid_out, dtype=np.float64)
    terms = _inlet_terms(cell)
    hin = float(sum(term[1] for term in terms))
    hout = float(_h_gas(cell, gas, T) + _h_solid(cell, solid, T))
    residual = float(hin * (1.0 - float(cell.heat_loss_frac)) - hout)
    return EnergyTerms(hin=hin, hout=hout, residual=residual, terms=terms)


def temperature_root(cell: Cell, *, gas_out: np.ndarray, solid_out: np.ndarray, lo: float = 250.0, hi: float = 2500.0) -> float | None:
    base = energy_terms(cell)
    target = float(base.hin * (1.0 - float(cell.heat_loss_frac)))

    def residual_at(T: float) -> float:
        return float(target - (_h_gas(cell, gas_out, T) + _h_solid(cell, solid_out, T)))

    f_lo = residual_at(lo)
    f_hi = residual_at(hi)
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if f_lo * f_hi > 0.0:
        return None
    return float(brentq(residual_at, lo, hi))


def _gas_closed_outlet(cell: Cell) -> np.ndarray:
    full = cell.calc_gas_balance()
    n = cell.N_d.size
    return np.maximum(cell.N_d + cell.N_b + full[:n] + full[n : 2 * n], 0.0)


def _print_cell(cell: Cell, index: int) -> None:
    current = energy_terms(cell)
    gas_out = np.maximum(cell.N_d + cell.N_b, 0.0)
    solid_out = np.maximum(cell._solid_outflow_rates(), 0.0)
    gas_closed = _gas_closed_outlet(cell)
    root_current = temperature_root(cell, gas_out=gas_out, solid_out=solid_out)
    root_gas_closed = temperature_root(cell, gas_out=gas_closed, solid_out=solid_out)

    print(
        f"cell={index} type={cell.cell_type} T={cell.T:.6g} "
        f"residual_W={current.residual:.6g} Hin_W={current.hin:.6g} Hout_W={current.hout:.6g} "
        f"root_current_K={root_current if root_current is not None else 'none'} "
        f"root_gas_closed_K={root_gas_closed if root_gas_closed is not None else 'none'}"
    )
    for name, enthalpy, flow, temperature in current.terms:
        print(f"  {name:16s} H_W={enthalpy:.6g} flow={flow:.6g} T_K={temperature:.6g}")
    print(f"  out_gas          H_W={_h_gas(cell, gas_out, cell.T):.6g} flow={float(np.sum(gas_out)):.6g}")
    print(f"  out_solid        H_W={_h_solid(cell, solid_out, cell.T):.6g} flow={float(np.sum(solid_out)):.6g}")


def main(indices: Iterable[int] = (0, 1)) -> None:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    run_init_and_precalc_for_global_nr(reactor, init_strategy="vorabrechnung", gs_warmup_steps=None)
    reactor._apply_all_bc_for_nr()
    for i in indices:
        _print_cell(reactor.cells[int(i)], int(i))


if __name__ == "__main__":
    main()
