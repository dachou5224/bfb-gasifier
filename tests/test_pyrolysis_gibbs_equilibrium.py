"""热解 Gibbs 必须同时满足元素守恒和化学势驻点。"""

import json
from pathlib import Path

import numpy as np

import src.solvers.vorabrechnung.vorab_gibbs_seed as gibbs
from src.solvers.external_state_ledger import external_imbalance
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state


def _trial_inputs():
    saved = json.loads(
        Path("docs/audit_reports/equation_loop_2026-10-08/advance_0.746.json").read_text()
    )
    captured: dict = {}
    original = gibbs._solve_pyrolysis_gibbs_an_eq_b

    def record(**kwargs):
        captured["elements"] = kwargs["elements"]
        captured["seed"] = kwargs["mole_seed"]
        captured["T"] = kwargs["T_K"]
        captured["P"] = kwargs["P_Pa"]
        return original(**kwargs)

    gibbs._solve_pyrolysis_gibbs_an_eq_b = record
    try:
        outer_P_named(
            np.asarray(saved["trial_x"], dtype=float),
            q_on_wesseling1_path(float(saved["s"])),
        )
    finally:
        gibbs._solve_pyrolysis_gibbs_an_eq_b = original
    return captured


def _composition(result: dict) -> np.ndarray:
    return np.array(
        [result[name] for name in gibbs._PYROLYSIS_GIBBS_CANDIDATES],
        dtype=float,
    )


def test_saved_trial_is_stationary_from_either_start() -> None:
    case = _trial_inputs()
    warm, warm_diag = gibbs._solve_pyrolysis_gibbs_an_eq_b(
        T_K=case["T"], P_Pa=case["P"], elements=case["elements"], mole_seed=case["seed"],
    )
    cold, cold_diag = gibbs._solve_pyrolysis_gibbs_an_eq_b(
        T_K=case["T"], P_Pa=case["P"], elements=case["elements"], mole_seed=None,
    )
    for diag in (warm_diag, cold_diag):
        assert diag["element_feasible"] is True
        assert diag["stationarity_abs"] <= 1.0e-12
        assert diag["max_abs_element_closure"] <= 1.0e-12
        assert diag["converged"] is True
        assert diag["converged_by"] != "not_equilibrium"
    assert np.linalg.norm(_composition(warm) - _composition(cold)) <= 1.0e-8


def test_fine_temperature_response_stays_on_the_same_slope() -> None:
    case = _trial_inputs()
    base, _diag = gibbs._solve_pyrolysis_gibbs_an_eq_b(
        T_K=case["T"], P_Pa=case["P"], elements=case["elements"], mole_seed=case["seed"],
    )
    origin = _composition(base)
    slopes = []
    for step in (1.0e-6, 1.0e-5, 1.0e-4, 0.01, 0.1):
        moved, _diag = gibbs._solve_pyrolysis_gibbs_an_eq_b(
            T_K=case["T"] + step,
            P_Pa=case["P"],
            elements=case["elements"],
            mole_seed=case["seed"],
        )
        slopes.append(np.linalg.norm(_composition(moved) - origin) / step)
    reference = slopes[-2]
    for slope in slopes:
        assert abs(slope - reference) / reference < 0.01


def test_saved_failure_trial_passes_the_original_gates() -> None:
    saved = json.loads(
        Path("docs/audit_reports/equation_loop_2026-10-08/advance_0.746.json").read_text()
    )
    x = np.asarray(saved["trial_x"], dtype=float)
    q = outer_P_named(x, q_on_wesseling1_path(float(saved["s"])))
    residual = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, residual, q.hydro)
    gap = external_imbalance(x, q)["gap_mol_s"]
    assert cert["passed"] is True
    assert max(abs(float(gap[key])) for key in ("C", "H", "O")) < 1.0e-4
    assert float(np.max(np.abs(q.hydro.w_undecided))) <= 1.0e-12
