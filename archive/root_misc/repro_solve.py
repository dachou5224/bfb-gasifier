
import sys
import os
sys.path.insert(0, os.path.abspath("."))

from src.core.reactor import Reactor
from tests.validation_case_utils import (
    build_phase1_htw_lu_reactor_config,
    load_case_LU
)

print("Loading case...")
case = load_case_LU()
print("Building config...")
cfg = build_phase1_htw_lu_reactor_config(case)
cfg.n_cells = 3 # Use fewer cells for quick check
print(f"Reactor config: {cfg.n_cells} cells, P={cfg.P} Pa, fuel_feed={cfg.fuel_feed} kg/s")

reactor = Reactor(cfg)
print("Starting solve with verbose=True...")
result = reactor.solve(
    max_global_iter=5,
    tol_global=1.0,
    solver="gauss_seidel",
    verbose=True
)

print("Solve complete.")
T_exit = result["T_profile"][-1]
print(f"Exit Temperature: {T_exit:.2f} K")
