"""Cell 解耦子模块审计。

用法：
    python3 scripts/audit_cell_submodels.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import Cell, CellGeometry, SolidProps, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_balance_residual
from src.core.cell_kinetics import build_reaction_sources
from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, configure_tar_components_by_fuel, gas_diffusivity_correlation, get_atom_count


def build_audit_cell() -> Cell:
    configure_tar_components_by_fuel("coal")
    c = Cell(
        geo=CellGeometry(D_bed=0.6, dh=1.0, h_center=0.5),
        solid=SolidProps(n_size_classes=1, d_p=0.5e-3),
        fuel_type="coal",
    )
    c.T = 1150.0
    c.P = 2.5e6
    idx = GAS_SPECIES_INDEX
    c.N_d[idx["O2"]] = 1.0
    c.N_d[idx["H2O"]] = 2.0
    c.N_d[idx["CO"]] = 0.2
    c.N_d[idx["H2"]] = 0.3
    c.N_d[idx["N2"]] = 20.0
    c.N_b[idx["O2"]] = 0.3
    c.N_b[idx["H2O"]] = 0.8
    c.N_b[idx["CO"]] = 0.1
    c.N_b[idx["N2"]] = 8.0
    c.m_solid[0, S_CHAR] = 0.05
    c.m_solid[0, S_VM] = 0.02
    c.m_solid[0, S_MOISTURE] = 0.01
    c.m_solid_zu[0, S_VM] = 0.02
    c.m_solid_zu[0, S_MOISTURE] = 0.01
    c.calc_hydrodynamics()
    return c


def audit_balance_exchange_only() -> None:
    n_ex = np.array([0.2, -0.1, 0.03] + [0.0] * 8, dtype=np.float64)
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
    print("[balance]")
    print(f"  exchange closure max = {np.max(np.abs(res[:11] + res[11:])):.3e}")


def audit_pyrolysis_module(cell: Cell) -> None:
    bundle = calc_drying_pyrolysis_sources(
        tau=cell.geo.dh / max(cell.u_mf, 1e-3),
        T=cell.T,
        P=cell.P,
        d_p=cell.solid.d_p,
        moisture_wt=cell.solid.moisture_wt,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        pyrolysis_tar_carbon_frac=cell.solid.pyrolysis_tar_carbon_frac,
        fuel_type=cell.fuel_type,
        m_vm_in=float(np.sum(cell.m_solid_zu[:, S_VM] + cell.m_solid_in[:, S_VM])),
        m_moist_in=float(np.sum(cell.m_solid_zu[:, S_MOISTURE] + cell.m_solid_in[:, S_MOISTURE])),
        solid_shape=cell.R_solid.shape,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )
    print("[pyrolysis]")
    print(f"  x_dry={bundle.x_dry:.4f}, x_vm={bundle.x_vm:.4f}")
    print(f"  gas_source_sum={bundle.gas_source.sum():.3e} mol/s")
    print(f"  solid_sink_sum={bundle.solid_sink.sum():.3e} kg/s")


def audit_kinetics_module(cell: Cell) -> None:
    tau = cell.geo.dh / max(cell.u_mf, 1e-3)
    cell.compute_vorabrechnung(tau)
    bundle = build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=cell._concentrations("b"),
        C_d=cell._concentrations("d"),
        y_b=cell._mole_fractions("b"),
        y_d=cell._mole_fractions("d"),
        gas_src_vm=cell._vm_gas_source_cache,
        solid_sink_vm=cell._vm_solid_sink_cache,
        areas=cell._calc_char_surface_area_per_class(),
        solid_d_p=cell.solid.d_p,
        D_g=gas_diffusivity_correlation(cell.T, cell.P),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=None,
        r4_scale=cell.r4_scale,
        r5_scale=cell.r5_scale,
        r6_scale=cell.r6_scale,
        r7_scale=cell.r7_scale,
        rate_multiplier=1.0,
        N_zu_d=cell.N_zu_d,
        N_d_in=cell.N_d_in,
        N_zu_b=cell.N_zu_b,
        N_b_in=cell.N_b_in,
        N_rez_d=cell.N_rez_d,
        N_rez_b=cell.N_rez_b,
        solid_shape=cell.R_solid.shape,
        char_index=S_CHAR,
    )
    gas_total = bundle.R_gas_b + bundle.R_gas_d
    elem = {}
    for el in ("C", "H", "O", "N", "S"):
        gas_elem = 0.0
        for i, sp in enumerate(GAS_SPECIES):
            gas_elem += float(gas_total[i]) * float(get_atom_count(sp, el))
        solid_elem = 0.0
        if el == "C":
            solid_elem += float(np.sum(bundle.R_solid[:, S_CHAR])) / 0.012011
        elif el == "H":
            solid_elem += 2.0 * float(np.sum(bundle.R_solid[:, S_MOISTURE])) / 0.018015
        elif el == "O":
            solid_elem += 1.0 * float(np.sum(bundle.R_solid[:, S_MOISTURE])) / 0.018015
        elem[el] = gas_elem + solid_elem
    print("[kinetics]")
    print(f"  limit_factor_o2={bundle.limit_factor_o2:.4f}")
    print(f"  limit_factor_h2o={bundle.limit_factor_h2o:.4f}")
    print(f"  element_closure={elem}")


def main() -> int:
    cell = build_audit_cell()
    audit_balance_exchange_only()
    audit_pyrolysis_module(cell)
    audit_kinetics_module(cell)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
