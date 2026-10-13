"""Wesseling 1 两区条件性预测。不读具名核端点，不把 1100 K 当目标。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_control import solve_named_controlled
from src.solvers.wesseling_case_contract import wesseling1_conditional_boundary
from src.solvers.wesseling_two_zone import (
    build_wesseling1_two_zone_q,
    freeboard_wet_mole_fraction,
    gas_phase_element_gap_kmol_h,
    wesseling1_physical_x0,
)

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "docs/audit_reports/plan_ab_diagnosis_2026-09-28/wp4_wesseling1.json"
_BUDGET_S = 180.0
_OBSERVED = {"CO": 0.13, "CO2": 0.11, "H2": 0.12, "CH4": 0.028, "H2O": 0.175}


def main() -> None:
    bound = wesseling1_conditional_boundary()
    nq = build_wesseling1_two_zone_q()
    x0 = wesseling1_physical_x0(nq)
    q0 = outer_P_named(x0, nq)
    out = solve_named_controlled(nq, x0, deadline_s=_BUDGET_S)
    x = np.asarray(out["x"], dtype=np.float64)
    certified = bool(out["passed"])
    report = {
        "conditional": True,
        "passed": certified,
        "reason": out["reason"],
        "elapsed_s": out["elapsed_s"],
        "merit": out["merit"],
        "gate_ratio": out["gate_ratio"],
        "T_bed": float(x[-2]),
        "T_fb": float(x[-1]),
        "prediction": None,
        "assumptions": {
            "steam_kg_h": bound["steam_kg_h"],
            "steam_status": bound["steam_status"],
            "secondary_air_placed_at_xi": bound["secondary_air_placed_at_xi"],
            "formation_enthalpy_basis": bound["formation_enthalpy_basis"],
            "char_diameter_basis": bound["char_diameter_basis"],
            "Q_W": bound["Q_W"],
            "launch_window_s": bound["launch_window_s"],
            "H_bed_status": bound["H_bed_status"],
        },
        "initial_w_undecided": np.asarray(q0.hydro.w_undecided, dtype=float).tolist(),
        "initial_tau_s": np.asarray(q0.hydro.tau_s, dtype=float).tolist(),
    }
    if certified:
        y = freeboard_wet_mole_fraction(x)
        report["prediction"] = {
            "T_fb_minus_paper_simulation_K": float(x[-1]) - float(bound["paper_simulation_temperature_K"]),
            "T_fb_minus_measured_K": float(x[-1]) - float(bound["measured_temperature_K"]),
            "wet_mole_fraction": y,
            "mole_fraction_minus_digitized": {
                name: float(y[name]) - target for name, target in _OBSERVED.items()
            },
            "digitization_rel_not_a_model_allowance": bound["digitization_rel"],
            "gas_phase_element_gap_kmol_h": gas_phase_element_gap_kmol_h(x),
            "gas_phase_gap_omits_solids_and_tar": True,
        }
    _OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(certified, out["reason"], round(out["elapsed_s"], 1), round(float(x[-1]), 2), flush=True)
    print("wrote", _OUT, flush=True)


if __name__ == "__main__":
    main()
