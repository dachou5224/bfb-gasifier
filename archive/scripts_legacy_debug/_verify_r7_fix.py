"""Verify rate_R7 fix."""
import sys; sys.path.insert(0, ".")
from src.kinetics.gas_reactions import rate_R7
from src.core.constants import Rg, P0

T = 1200.0
P = 2.5e6
C_fac = P / (Rg * T)
# Cell 1 inlet concentrations
y_CO = 2.19/61.8; y_H2 = 2.24/61.8; y_CH4 = 0.003/61.8; y_H2O = 6.10/61.8
C_CO  = y_CO * C_fac
C_H2  = y_H2 * C_fac
C_CH4 = y_CH4 * C_fac
C_H2O = y_H2O * C_fac
print(f"C_CH4={C_CH4:.4f}, C_H2O={C_H2O:.4f}, C_CO={C_CO:.4f}, C_H2={C_H2:.4f} mol/m3")

r = rate_R7(T, C_CH4, C_H2O, C_CO, C_H2)
print(f"rate_R7 at cell1 conditions: {r:.6f} mol/(m3*s)  (was ~-596000)")

r_zero = rate_R7(T, 0, 0, 0, 0)
print(f"rate_R7 at zero: {r_zero:.6f}  (should be 0)")

r_fwd = rate_R7(T, 1.0, 10.0, 0.01, 0.01)
print(f"rate_R7 forward (CH4+H2O): {r_fwd:.4f} mol/(m3*s)  (should be positive)")

r_rev = rate_R7(T, 0.001, 0.001, 5.0, 8.0)
print(f"rate_R7 reverse (CO+H2): {r_rev:.4f} mol/(m3*s)  (should be negative, bounded)")
print(f"  (was: -510000 mol/(m3*s) previously)")
