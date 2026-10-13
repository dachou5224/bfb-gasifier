"""Rebuild the conditional s=1 root from the saved state. No hidden cache and no new solve."""
import json
from pathlib import Path
import numpy as np
from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import outer_P_named, coupled_residual_named, unpack_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state
from src.solvers.external_state_ledger import external_imbalance

_CERT = Path("docs/audit_reports/equation_loop_2026-10-08/certified.json")
_INDEX = {
    "bed_gas_0_7": "CO CO2 H2 H2O CH4 O2 N2 mol/s",
    "bubble_gas_7_14": "same species mol/s",
    "bed_char_14_16": "fine, coarse kg",
    "freeboard_gas_16_23": "same species mol/s",
    "cyclone1_gas_23_30": "same species mol/s",
    "cyclone2_gas_30_37": "same species mol/s",
    "cyclone_char_37_39": "M_c1, M_c2 kg",
    "ash_39": "kg",
    "bed_tar_40_42": "TAR1 TAR2 mol/s",
    "bubble_tar_42_44": "TAR1 TAR2 mol/s",
    "freeboard_tar_44_46": "TAR1 TAR2 mol/s",
    "freeboard_char_46_48": "fine, coarse kg",
    "T_bed_48": "K",
    "T_fb_49": "K",
}


def rebuild():
    saved = json.loads(_CERT.read_text())
    x = np.asarray(saved["x"], dtype=float)
    q = outer_P_named(x, q_on_wesseling1_path(float(saved["s"])))
    residual = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, residual, q.hydro)
    ledger = external_imbalance(x, q)
    state = unpack_named(x, 11)
    return {
        "command": "PYTHONPATH=. .venv/bin/python docs/audit_reports/s1_review_2026-10-12/reproduce_s1.py",
        "saved_s": float(saved["s"]),
        "saved_source_sha256": saved["source_sha256"],
        "source_sha256": named_kernel_source_sha256(),
        "call": {"ds": 0.002, "temperature_radius_scale": 4.0, "settling_margin_eta": None, "energy_guard": False},
        "index": _INDEX,
        "assumptions": {
            "steam_kg_h": 0.0,
            "Q_W": 0.0,
            "h_f_dry_J_kg": -3.96e6,
            "h_f_scope": "schematic dry-fuel formation enthalpy, not HHV or LHV",
            "char_diameter_m": [2.0e-4, 3.0e-4],
            "all_computed_gas_enters_bubble_inlet": True,
            "experimental_validation": False,
        },
        "passed": bool(cert["passed"]),
        "gate_ratio": cert["gate_ratio"],
        "energy_rows_W": cert["energy_rows_W"],
        "cho_gap_mol_s": {key: float(ledger["gap_mol_s"][key]) for key in "CHO"},
        "enthalpy_gap_W": float(ledger["enthalpy_W"]["gap"]),
        "transport_closed": bool(cert["transport_closed"]),
        "w_undecided": np.asarray(q.hydro.w_undecided).tolist(),
        "T_bed_K": float(state["T"]),
        "T_fb_K": float(state["T_fb"]),
        "bypass_mol_s": {
            "NH3": float(getattr(q.energy.species, "n_NH3_bypass_mol_s", 0.0)),
            "H2S": float(getattr(q.energy.species, "n_H2S_bypass_mol_s", 0.0)),
        },
    }


if __name__ == "__main__":
    report = rebuild()
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit("conditional s=1 root did not pass the original gates")
