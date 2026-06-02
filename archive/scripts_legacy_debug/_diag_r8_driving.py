"""诊断R8 WGSR在零浓度时的driving force和reaction rate。"""
import sys
sys.path.insert(0, ".")

import numpy as np
from src.kinetics.gas_reactions import rate_R8
from src.thermodynamics.equilibrium import (
    calc_gibbs_driving_force, calc_reaction_quotient
)
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from src.kinetics.gas_reactions import wgsr_equilibrium_constant

T = 1200.0; P = 2.5e6

# Rate R8 with zero mole fracs
r8 = rate_R8(T, P, y_CO=0.0, y_H2O=0.0, y_CO2=0.0, y_H2=0.0)
print(f"R8 with all-zero y: {r8:.4f} mol/(m3*s)")

y_zero = {"CO": 0.0, "H2O": 0.0, "CO2": 0.0, "H2": 0.0}
df = calc_gibbs_driving_force("R8", T, P, y_zero, clamp_irreversible=False)
qp = calc_reaction_quotient("R8", y_zero, P)
keq = wgsr_equilibrium_constant(T)
print(f"R8 Q_p={qp:.6e}, K_eq={keq:.6f}, driving={df:.6f}")
print()

cfg = ReactorConfig(
    n_cells=3, H_bed=6.0, D_bed=0.6, P=2_500_000.0,
    T_inlet=293.0, fuel_type="coal",
    fuel_feed=0.938, moisture_wt=16.9, ash_dry_wt=11.41, VM_daf=53.42,
    C_dry=61.5, H_dry=4.1, O_dry=21.8, HHV_dry=22.0,
    ER=0.337, primary_agent="air_steam", steam_to_o2_molar=0.8,
    recirculation_frac=0.0, heat_loss_frac=0.0,
)
reactor = Reactor(cfg)
reactor._set_bottom_cell_feeds()
c0 = reactor.cells[0]
c0.T = 1200.0
c0.calc_hydrodynamics()
c0.compute_vorabrechnung(c0.geo.dh / max(c0.u_mf, 1e-3))

y_d = c0._mole_fractions("d")
print("y_d before reactions:", dict(zip(GAS_SPECIES, [f"{v:.4e}" for v in y_d])))
print()

c0.calc_reactions()
print("R_gas_d all entries:")
for i, sp in enumerate(GAS_SPECIES):
    v = c0.R_gas_d[i]
    print(f"  {sp:6s}: {v:+.4f}")
print("Sum R_gas_d:", c0.R_gas_d.sum())
print()

# Also check ext8 directly
from src.kinetics.gas_reactions import rate_R8
idx = GAS_SPECIES_INDEX
C_d = c0._concentrations("d")
y_d2 = c0._mole_fractions("d")
ext8_raw = rate_R8(T, P, float(y_d2[idx["CO"]]), float(y_d2[idx["H2O"]]),
                   float(y_d2[idx["CO2"]]), float(y_d2[idx["H2"]])) * c0.V_d
print(f"ext8_raw = {ext8_raw:.4f} mol/s (V_d={c0.V_d:.4f} m3)")
