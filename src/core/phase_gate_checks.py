"""Phase gate checks for staged structural and elemental revisions.

These checks are intentionally read-only. They freeze baseline snapshots and
report whether the repository still has the minimum artifacts and smoke gates
needed to continue the phased plan.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import contextlib
import importlib
import importlib.util
import io
import json
from pathlib import Path
import time
from typing import Any

import numpy as np

from src.core.cell import Cell, CellGeometry, SolidProps, S_CHAR, S_MOISTURE, S_VM
from src.core.cell_balances import calc_gas_balance_residual
from src.core.cell_hydrodynamics import calc_phase_exchange
from src.core.cell_kinetics import build_reaction_sources
from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources
from src.core.elemental_ledger import (
    TRACKED_ELEMENTS,
    build_gas_species_registry,
    combined_element_molar_rates,
    max_absolute_element_rate,
)
from src.core.reactor import Reactor, _propagated_solid_stream, _recycled_solid_stream
from src.core.reaction_numbering import map_impl_diag_to_thesis, reaction_numbering_metadata
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, N_GAS, configure_tar_components_by_fuel, cp_molar, enthalpy_molar, gas_diffusivity_correlation

_REPO_ROOT = Path(__file__).resolve().parents[2]
_BASELINE_PATH = _REPO_ROOT / "data" / "phase_gate_baselines.json"
# final solved-state 的 char-dominant 尾部仍可能出现 O(1e-4) 量级的残余
# VM/moisture 重分布；Phase 4 关注的是上部 fresh-feed release 是否重新点燃，
# 而不是要求尾部库存严格机器精度单调。
_PHASE4_REBOUND_TOL = 1e-3


@dataclass
class GateResult:
    """Unified phase-gate result payload."""

    gate_id: str
    phase: str
    passed: bool
    hard_failures: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    artifacts: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_phase_gate_baselines() -> dict[str, Any]:
    """Load frozen phase-gate baseline snapshots."""

    with open(_BASELINE_PATH, encoding="utf-8") as f:
        return json.load(f)


def _load_validation_utils():
    return importlib.import_module("tests.validation_case_utils")


def _load_repo_script_module(module_name: str, relative_path: str):
    script_path = _REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"无法从 {script_path} 加载模块 {module_name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_sanity_checks() -> tuple[int, str]:
    sanity_module = importlib.import_module("tests.sanity_checks")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = int(sanity_module.main())
    return code, buf.getvalue()


def _build_phase0_cell_fixture() -> Cell:
    configure_tar_components_by_fuel("coal")
    idx = GAS_SPECIES_INDEX
    cell = Cell(
        geo=CellGeometry(D_bed=0.6, dh=1.0, h_center=0.5),
        solid=SolidProps(n_size_classes=1, d_p=0.5e-3),
        fuel_type="coal",
    )
    cell.T = 1150.0
    cell.P = 2.5e6
    cell.N_d[idx["O2"]] = 1.0
    cell.N_d[idx["H2O"]] = 2.0
    cell.N_d[idx["CO"]] = 0.2
    cell.N_d[idx["H2"]] = 0.3
    cell.N_d[idx["N2"]] = 20.0
    cell.N_b[idx["O2"]] = 0.3
    cell.N_b[idx["H2O"]] = 0.8
    cell.N_b[idx["CO"]] = 0.1
    cell.N_b[idx["N2"]] = 8.0
    cell.m_solid[0, S_CHAR] = 0.05
    cell.m_solid[0, S_VM] = 0.02
    cell.m_solid[0, S_MOISTURE] = 0.01
    cell.m_solid_zu[0, S_VM] = 0.02
    cell.m_solid_zu[0, S_MOISTURE] = 0.01
    return cell


def build_phase0_cell_snapshot() -> dict[str, float]:
    """Freeze a reusable cell-level snapshot for later structural parity gates."""

    idx = GAS_SPECIES_INDEX
    cell = _build_phase0_cell_fixture()
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    cell.calc_reactions()
    residual = cell.residuals()

    assert cell.N_ex.shape == cell.N_d.shape == cell.N_b.shape
    assert cell.R_solid.shape == cell.m_solid.shape

    return {
        "T_K": float(cell.T),
        "P_Pa": float(cell.P),
        "u0_m_s": float(cell.u0),
        "u_mf_m_s": float(cell.u_mf),
        "u_b_m_s": float(cell.u_b),
        "d_b_m": float(cell.d_b),
        "eps_b": float(cell.eps_b),
        "eps_d": float(cell.eps_d),
        "K_bd_1_s": float(cell.K_bd),
        "V_b_m3": float(cell.V_b),
        "V_d_m3": float(cell.V_d),
        "N_ex_l1": float(np.linalg.norm(cell.N_ex, ord=1)),
        "N_ex_sum": float(np.sum(cell.N_ex)),
        "N_ex_O2": float(cell.N_ex[idx["O2"]]),
        "R_gas_b_l1": float(np.linalg.norm(cell.R_gas_b, ord=1)),
        "R_gas_d_l1": float(np.linalg.norm(cell.R_gas_d, ord=1)),
        "R_solid_l1": float(np.linalg.norm(cell.R_solid, ord=1)),
        "residual_l2": float(np.linalg.norm(residual)),
        "residual_linf": float(np.max(np.abs(residual))),
    }


def _compute_validation_metrics(reactor_snapshot: dict[str, Any]) -> dict[str, Any]:
    utils = _load_validation_utils()
    case_key = utils.CASE_LU_VALIDATION_KEY
    node = utils.load_validation_case_node(case_key)
    outs = node["outputs"]

    target_T = float(outs["exit_temperature_K"])
    measured_T = float(outs.get("exit_temperature_measured_K", target_T))
    T_exit = float(reactor_snapshot["T_exit_K"])
    err_T = min(
        abs(T_exit - target_T) / max(target_T, 1.0),
        abs(T_exit - measured_T) / max(measured_T, 1.0),
    )

    dry_ref = outs.get("exit_gas_dry_mol_frac", {})
    species_err: dict[str, float] = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        tgt = utils.json_numeric_or_none(dry_ref.get(sp))
        sim = reactor_snapshot["exit_gas_dry"].get(sp)
        if tgt is None or tgt <= 0.0 or sim is None:
            continue
        species_err[sp] = abs(float(sim) - float(tgt)) / float(tgt)

    carbon_target_pct = utils.json_numeric_or_none(outs.get("carbon_conversion_pct"))
    carbon_err = None
    if carbon_target_pct is not None and carbon_target_pct > 0.0:
        carbon_err = abs(reactor_snapshot["carbon_conv"] * 100.0 - carbon_target_pct) / carbon_target_pct

    return {
        "case_key": case_key,
        "err_T_vs_validation_best": float(err_T),
        "species_relative_error": species_err,
        "carbon_relative_error": None if carbon_err is None else float(carbon_err),
    }


def build_phase0_lu_snapshot() -> dict[str, Any]:
    """Freeze the current LU reactor-level baseline snapshot."""

    utils = _load_validation_utils()
    case = utils.load_case_LU()
    cfg = utils.build_phase1_htw_lu_reactor_config(case)
    solve_kwargs = dict(utils.PHASE1_HTW_LU_SOLVE_KWARGS)

    t0 = time.perf_counter()
    result = Reactor(cfg).solve(**solve_kwargs)
    runtime_s = time.perf_counter() - t0

    snapshot = {
        "case_key": utils.CASE_LU_VALIDATION_KEY,
        "solver": solve_kwargs.get("solver"),
        "solve_kwargs": solve_kwargs,
        "T_exit_K": float(result["T_profile"][-1]),
        "carbon_conv": float(result["carbon_conv"]),
        "n_iter": int(result.get("n_iter", 0) or 0),
        "converged": bool(result.get("converged", False)),
        "residual_gs": float(result.get("residual_gs", 0.0) or 0.0),
        "rms_scaled_gs": float(result.get("rms_scaled_gs", 0.0) or 0.0),
        "runtime_s": float(runtime_s),
        "exit_gas_dry": {
            str(k): float(v) for k, v in result["exit_gas_dry"].items() if isinstance(v, (int, float))
        },
    }
    snapshot["validation"] = _compute_validation_metrics(snapshot)
    return snapshot


def _diff_scalar(
    *,
    label: str,
    actual: float,
    baseline: float,
    rtol: float,
    warnings: list[str],
    drift: dict[str, Any],
    atol: float = 0.0,
) -> None:
    abs_diff = abs(actual - baseline)
    rel_diff = abs_diff / max(abs(baseline), 1e-12)
    drift[label] = {
        "actual": float(actual),
        "baseline": float(baseline),
        "abs_diff": float(abs_diff),
        "rel_diff": float(rel_diff),
    }
    if abs_diff > max(atol, rtol * max(abs(baseline), 1e-12)):
        warnings.append(
            f"{label} 漂移超出容差: actual={actual:.6g}, baseline={baseline:.6g}, "
            f"rel_diff={rel_diff:.3e}, rtol={rtol:.3e}, atol={atol:.3e}"
        )


def compare_phase0_snapshots(
    *,
    current_cell: dict[str, float],
    current_reactor: dict[str, Any],
    baseline: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    """Compare current phase-0 snapshots to the frozen baseline."""

    warnings: list[str] = []
    drift: dict[str, Any] = {"cell": {}, "reactor": {}, "runtime": {}}
    tol = baseline["tolerances"]

    for key, rtol in tol["cell_relative"].items():
        _diff_scalar(
            label=f"cell.{key}",
            actual=float(current_cell[key]),
            baseline=float(baseline["cell_snapshot"][key]),
            rtol=float(rtol),
            atol=float(tol["cell_absolute"].get(key, 0.0)),
            warnings=warnings,
            drift=drift["cell"],
        )

    if "N_ex_sum" in baseline["cell_snapshot"]:
        _diff_scalar(
            label="cell.N_ex_sum",
            actual=float(current_cell["N_ex_sum"]),
            baseline=float(baseline["cell_snapshot"]["N_ex_sum"]),
            rtol=0.0,
            atol=float(tol["cell_absolute"].get("N_ex_sum", 0.0)),
            warnings=warnings,
            drift=drift["cell"],
        )

    for key, rtol in tol["reactor_relative"].items():
        if key == "dry_species":
            for sp, val in baseline["reactor_snapshot"]["exit_gas_dry"].items():
                if sp not in current_reactor["exit_gas_dry"]:
                    warnings.append(f"reactor.exit_gas_dry 缺少物种 {sp}")
                    continue
                _diff_scalar(
                    label=f"reactor.exit_gas_dry.{sp}",
                    actual=float(current_reactor["exit_gas_dry"][sp]),
                    baseline=float(val),
                    rtol=float(rtol),
                    warnings=warnings,
                    drift=drift["reactor"],
                )
            continue
        _diff_scalar(
            label=f"reactor.{key}",
            actual=float(current_reactor[key]),
            baseline=float(baseline["reactor_snapshot"][key]),
            rtol=float(rtol),
            warnings=warnings,
            drift=drift["reactor"],
        )

    if bool(current_reactor["converged"]) != bool(baseline["reactor_snapshot"]["converged"]):
        warnings.append(
            "reactor.converged 状态发生变化: "
            f"actual={current_reactor['converged']} baseline={baseline['reactor_snapshot']['converged']}"
        )
    if int(current_reactor["n_iter"]) != int(baseline["reactor_snapshot"]["n_iter"]):
        warnings.append(
            "reactor.n_iter 发生变化: "
            f"actual={current_reactor['n_iter']} baseline={baseline['reactor_snapshot']['n_iter']}"
        )

    baseline_runtime = float(baseline["reactor_snapshot"]["runtime_s"])
    runtime_factor = float(current_reactor["runtime_s"]) / max(baseline_runtime, 1e-9)
    drift["runtime"] = {
        "actual_s": float(current_reactor["runtime_s"]),
        "baseline_s": baseline_runtime,
        "factor": runtime_factor,
    }
    if runtime_factor > float(tol["runtime_warn_factor"]):
        warnings.append(
            f"reactor runtime 超过基线因子阈值: actual={current_reactor['runtime_s']:.2f}s, "
            f"baseline={baseline_runtime:.2f}s, factor={runtime_factor:.2f}"
        )

    return warnings, drift


_PHASE1_FINAL_OWNER_EXPECTATIONS = {
    "u_mf": "Cell.calc_hydrodynamics",
    "u_b": "Cell.calc_hydrodynamics",
    "d_b": "Cell.calc_hydrodynamics",
    "eps_b": "Cell.calc_hydrodynamics",
    "eps_d": "Cell.calc_hydrodynamics",
    "K_bd": "Cell.calc_hydrodynamics",
    "V_b": "Cell.calc_hydrodynamics",
    "V_d": "Cell.calc_hydrodynamics",
    "u0": "Cell.calc_hydrodynamics",
    "N_ex": "Cell.calc_exchange",
    "R_gas_b": "Cell.calc_reactions",
    "R_gas_d": "Cell.calc_reactions",
    "R_solid": "Cell.calc_reactions",
    "_vm_gas_source_cache": "Cell.compute_vorabrechnung",
    "_vm_solid_sink_cache": "Cell.compute_vorabrechnung",
    "_vm_cache_valid": "Cell.compute_vorabrechnung",
    "_vm_cache_T": "Cell.compute_vorabrechnung",
    "_vm_cache_tau": "Cell.compute_vorabrechnung",
    "_vm_cache_m_vm_in": "Cell.compute_vorabrechnung",
    "_vm_cache_m_moist_in": "Cell.compute_vorabrechnung",
    "_work_res": "Cell.residuals",
    "N_d": "cell_solver._unpack_state",
    "N_b": "cell_solver._unpack_state",
    "m_solid": "cell_solver._unpack_state",
    "T": "cell_solver._unpack_state",
}


def build_phase1_structural_contract() -> dict[str, Any]:
    """Build the phase-1 structural mutation-owner contract."""

    mutators = [
        {
            "name": "Cell.calc_hydrodynamics",
            "category": "derived-state",
            "final_owner": True,
            "writes": ["u_mf", "u_b", "d_b", "eps_b", "eps_d", "K_bd", "V_b", "V_d", "u0"],
            "source": "src/core/cell.py:107-135",
        },
        {
            "name": "Cell.calc_exchange",
            "category": "derived-state",
            "final_owner": True,
            "writes": ["N_ex"],
            "source": "src/core/cell.py:137-143",
        },
        {
            "name": "Cell._calc_drying_pyrolysis_gas_source",
            "category": "temporary-side-effect",
            "final_owner": False,
            "writes": ["R_solid"],
            "source": "src/core/cell.py:463-485",
            "note": "Temporary solid sink staging before Cell.calc_reactions performs final bundle write-back.",
        },
        {
            "name": "Cell.compute_vorabrechnung",
            "category": "cache-owner",
            "final_owner": True,
            "writes": [
                "_vm_gas_source_cache",
                "_vm_solid_sink_cache",
                "_vm_cache_valid",
                "_vm_cache_T",
                "_vm_cache_tau",
                "_vm_cache_m_vm_in",
                "_vm_cache_m_moist_in",
            ],
            "source": "src/core/cell.py:498-517",
        },
        {
            "name": "Cell.calc_reactions",
            "category": "source-owner",
            "final_owner": True,
            "writes": ["R_gas_b", "R_gas_d", "R_solid"],
            "source": "src/core/cell.py:255-318",
        },
        {
            "name": "Cell.residuals",
            "category": "assembled-output",
            "final_owner": True,
            "writes": ["_work_res"],
            "source": "src/core/cell.py:449-461",
        },
        {
            "name": "Reactor._set_bottom_cell_feeds",
            "category": "boundary-owner",
            "final_owner": True,
            "writes": [
                "bottom.N_zu_b",
                "bottom.N_zu_d",
                "bottom.m_solid_zu",
                "bottom.T_zu_gas",
                "bottom.T_zu_solid",
                "bottom.T_in_gas",
                "bottom.T_in_solid",
            ],
            "source": "src/core/reactor.py:405-435",
        },
        {
            "name": "Reactor._apply_bottom_recycle",
            "category": "boundary-owner",
            "final_owner": True,
            "writes": [
                "bottom.N_rez_d",
                "bottom.N_rez_b",
                "bottom.m_solid_rez",
                "bottom.T_rez_gas",
                "bottom.T_rez_solid",
            ],
            "source": "src/core/reactor.py:437-478",
        },
        {
            "name": "Reactor._propagate_upstream",
            "category": "boundary-owner",
            "final_owner": True,
            "writes": [
                "curr.N_b_in",
                "curr.N_d_in",
                "curr.T_in_gas",
                "curr.m_solid_in",
                "curr.T_in_solid",
            ],
            "source": "src/core/reactor.py:480-512",
        },
        {
            "name": "cell_solver._unpack_state",
            "category": "primitive-state-owner",
            "final_owner": True,
            "writes": ["N_d", "N_b", "m_solid", "T"],
            "source": "src/solvers/cell_solver.py:52-63",
        },
    ]

    orchestrators = [
        {
            "name": "Cell.residuals",
            "calls": [
                "Cell.calc_hydrodynamics",
                "Cell.calc_exchange",
                "Cell.calc_reactions",
                "assemble_cell_residual_vector",
            ],
            "source": "src/core/cell.py:449-461",
            "note": "Impure evaluation path: recalculates derived state and sources before writing _work_res.",
        },
        {
            "name": "cell_solver.evaluate_cell_state",
            "calls": ["Cell.residuals"],
            "source": "src/solvers/cell_solver.py:65-79",
        },
        {
            "name": "Reactor._evaluate_current_gs_state",
            "calls": [
                "Reactor._set_bottom_cell_feeds",
                "Reactor._apply_bottom_recycle",
                "Reactor._propagate_upstream",
                "cell_solver.evaluate_cell_state",
            ],
            "source": "src/core/reactor.py:522-534",
        },
    ]

    override_risks = [
        {
            "risk": "later conservation override earlier conservation",
            "guard": "R_solid temporary staging in Cell._calc_drying_pyrolysis_gas_source is not treated as final owner; Cell.calc_reactions is the only final source owner.",
        },
        {
            "risk": "boundary arrays overwritten from multiple paths",
            "guard": "Bottom fresh feed, bottom recycle, and upstream propagation each keep separate ownership domains.",
        },
        {
            "risk": "residual evaluation mutates derived state",
            "guard": "Cell.residuals is explicitly marked impure/orchestrated; parity gates must compare outputs, not assume pure read-only residuals.",
        },
    ]

    return {
        "phase": "phase-1",
        "mutators": mutators,
        "orchestrators": orchestrators,
        "final_owner_expectations": dict(_PHASE1_FINAL_OWNER_EXPECTATIONS),
        "override_risks": override_risks,
    }


def validate_phase1_structural_contract(contract: dict[str, Any]) -> list[str]:
    """Validate unique final-owner expectations for the structural parity phase."""

    failures: list[str] = []
    owners_by_field: dict[str, list[str]] = {}
    for item in contract["mutators"]:
        if not bool(item.get("final_owner", False)):
            continue
        for field_name in item.get("writes", []):
            owners_by_field.setdefault(field_name, []).append(item["name"])

    for field_name, expected_owner in contract["final_owner_expectations"].items():
        owners = owners_by_field.get(field_name, [])
        if len(owners) != 1:
            failures.append(
                f"{field_name} final owner count != 1: owners={owners if owners else '[]'}"
            )
            continue
        if owners[0] != expected_owner:
            failures.append(
                f"{field_name} final owner mismatch: expected={expected_owner}, actual={owners[0]}"
            )

    return failures


def gate_00_baseline_frozen(*, strict: bool = False) -> GateResult:
    """Phase 0 gate: baseline artifacts exist and are still reproducible enough to use."""

    hard_failures: list[str] = []
    warnings: list[str] = []
    metrics: dict[str, Any] = {}
    artifacts = {
        "baseline_json": str(_BASELINE_PATH),
        "sanity_script": str(_REPO_ROOT / "tests" / "sanity_checks.py"),
        "phase1_lu_audit_script": str(_REPO_ROOT / "scripts" / "audit_phase1_htw_lu.py"),
    }

    if not _BASELINE_PATH.exists():
        hard_failures.append(f"缺少 baseline JSON: {_BASELINE_PATH}")
        return GateResult(
            gate_id="gate_00_baseline_frozen",
            phase="phase-0",
            passed=False,
            hard_failures=hard_failures,
            warnings=warnings,
            metrics=metrics,
            artifacts=artifacts,
        )

    try:
        baseline_root = load_phase_gate_baselines()
        baseline = baseline_root["gate_00"]
    except Exception as exc:  # pragma: no cover - exercised by smoke gate
        hard_failures.append(f"baseline JSON 读取失败: {exc}")
        return GateResult(
            gate_id="gate_00_baseline_frozen",
            phase="phase-0",
            passed=False,
            hard_failures=hard_failures,
            warnings=warnings,
            metrics=metrics,
            artifacts=artifacts,
        )

    try:
        sanity_code, sanity_stdout = _run_sanity_checks()
    except Exception as exc:  # pragma: no cover - exercised by smoke gate
        hard_failures.append(f"sanity_checks 运行失败: {exc}")
        sanity_code, sanity_stdout = 1, ""

    metrics["sanity"] = {
        "passed": sanity_code == 0,
        "exit_code": int(sanity_code),
        "stdout": sanity_stdout,
    }
    if sanity_code != 0:
        hard_failures.append("python3 tests/sanity_checks.py 未通过")

    try:
        current_cell = build_phase0_cell_snapshot()
        current_reactor = build_phase0_lu_snapshot()
        metrics["current_cell_snapshot"] = current_cell
        metrics["current_reactor_snapshot"] = current_reactor
    except Exception as exc:  # pragma: no cover - exercised by smoke gate
        hard_failures.append(f"phase-0 baseline snapshot 构建失败: {exc}")
        current_cell = {}
        current_reactor = {}

    if not hard_failures:
        drift_warnings, drift = compare_phase0_snapshots(
            current_cell=current_cell,
            current_reactor=current_reactor,
            baseline=baseline,
        )
        warnings.extend(drift_warnings)
        metrics["baseline_drift"] = drift
        metrics["frozen_baseline"] = baseline

    passed = not hard_failures and (not strict or not warnings)
    return GateResult(
        gate_id="gate_00_baseline_frozen",
        phase="phase-0",
        passed=passed,
        hard_failures=hard_failures,
        warnings=warnings,
        metrics=metrics,
        artifacts=artifacts,
    )


def gate_10_structural_parity() -> GateResult:
    """Phase 1 gate: enforce the structural ownership contract.

    Phase-0 baseline drift is still reported for context, but after later solver
    phases intentionally move reactor-level behavior it should not be escalated
    as a hard structural-contract failure.
    """

    phase0 = gate_00_baseline_frozen(strict=False)
    contract = build_phase1_structural_contract()
    contract_failures = validate_phase1_structural_contract(contract)

    hard_failures = list(phase0.hard_failures)
    hard_failures.extend(contract_failures)

    metrics = dict(phase0.metrics)
    metrics["structural_contract"] = contract

    artifacts = dict(phase0.artifacts)
    artifacts["phase_1_contract_source"] = "src/core/phase_gate_checks.py"

    return GateResult(
        gate_id="gate_10_structural_parity",
        phase="phase-1",
        passed=not hard_failures,
        hard_failures=hard_failures,
        warnings=list(phase0.warnings),
        metrics=metrics,
        artifacts=artifacts,
    )


def _phase2_pyrolysis_closure() -> dict[str, Any]:
    cell = _build_phase0_cell_fixture()
    cell.calc_hydrodynamics()
    tau = cell.geo.dh / max(cell.u_mf, 1e-3)
    bundle = calc_drying_pyrolysis_sources(
        tau=tau,
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
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )
    closure = combined_element_molar_rates(
        gas_rates=bundle.gas_source,
        solid_rates=bundle.solid_sink,
        fuel_type=cell.fuel_type,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )
    return {
        "x_dry": float(bundle.x_dry),
        "x_vm": float(bundle.x_vm),
        "closure": closure,
        "max_abs_closure": max_absolute_element_rate(closure),
    }


def _phase2_reaction_closure() -> dict[str, Any]:
    cell = _build_phase0_cell_fixture()
    cell.calc_hydrodynamics()
    cell.calc_exchange()
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
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=S_CHAR,
    )
    closure = combined_element_molar_rates(
        gas_rates=bundle.R_gas_b + bundle.R_gas_d,
        solid_rates=bundle.R_solid,
        fuel_type=cell.fuel_type,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )
    net_impl = {
        "R5": float(bundle.net_molar_gas_source_r5),
        "R6": float(bundle.net_molar_gas_source_r6),
        "R7": float(bundle.net_molar_gas_source_r7),
        "R8": float(bundle.net_molar_gas_source_r8),
        "R9": float(bundle.net_molar_gas_source_r9),
        "R10": float(bundle.net_molar_gas_source_r10),
        "R11": float(bundle.net_molar_gas_source_r11),
        "R12": float(bundle.net_molar_gas_source_r12),
        "char": float(bundle.net_molar_gas_source_char),
        "vm": float(bundle.net_molar_gas_source_vm),
    }
    return {
        "limit_factor_o2": float(bundle.limit_factor_o2),
        "limit_factor_o2_bubble": float(bundle.limit_factor_o2_bubble),
        "limit_factor_o2_dense": float(bundle.limit_factor_o2_dense),
        "limit_factor_h2o": float(bundle.limit_factor_h2o),
        "net_molar_gas_source_total": float(bundle.net_molar_gas_source_total),
        "net_molar_gas_source_vm": float(bundle.net_molar_gas_source_vm),
        "net_molar_gas_source_r5": float(bundle.net_molar_gas_source_r5),
        "net_molar_gas_source_r6": float(bundle.net_molar_gas_source_r6),
        "net_molar_gas_source_r7": float(bundle.net_molar_gas_source_r7),
        "net_molar_gas_source_r8": float(bundle.net_molar_gas_source_r8),
        "net_molar_gas_source_r9": float(bundle.net_molar_gas_source_r9),
        "net_molar_gas_source_r10": float(bundle.net_molar_gas_source_r10),
        "net_molar_gas_source_r11": float(bundle.net_molar_gas_source_r11),
        "net_molar_gas_source_r12": float(bundle.net_molar_gas_source_r12),
        "net_molar_gas_source_char": float(bundle.net_molar_gas_source_char),
        "net_molar_gas_source_by_impl": net_impl,
        "net_molar_gas_source_by_thesis": map_impl_diag_to_thesis(
            {k: v for k, v in net_impl.items() if k.startswith("R")}
        ),
        "reaction_numbering": reaction_numbering_metadata(),
        "closure": closure,
        "max_abs_closure": max_absolute_element_rate(closure),
    }


def gate_20_elemental_closure() -> GateResult:
    """Phase 2 gate: elemental registry and source bundles must close on C/H/O/N/S."""

    hard_failures: list[str] = []
    warnings: list[str] = []

    registry_metrics: dict[str, Any] = {}
    for fuel_type in ("coal", "biomass"):
        registry = build_gas_species_registry(fuel_type)
        thermo_missing: list[str] = []
        zero_species: list[str] = []
        for species in GAS_SPECIES:
            atoms = registry[species]["atoms"]
            if sum(atoms[element] for element in TRACKED_ELEMENTS) <= 0:
                zero_species.append(species)
            try:
                cp_molar(species, 1200.0)
                enthalpy_molar(species, 1200.0)
            except Exception:
                thermo_missing.append(species)
        registry_metrics[fuel_type] = {
            "species": registry,
            "zero_element_species": zero_species,
            "thermo_missing": thermo_missing,
        }
        if zero_species:
            hard_failures.append(f"{fuel_type} registry 存在零元素物种: {zero_species}")
        if thermo_missing:
            hard_failures.append(f"{fuel_type} registry 缺少 thermo coverage: {thermo_missing}")

    pyrolysis_metrics = _phase2_pyrolysis_closure()
    reaction_metrics = _phase2_reaction_closure()

    if pyrolysis_metrics["max_abs_closure"] > 1e-8:
        hard_failures.append(
            f"pyrolysis elemental closure 超限: {pyrolysis_metrics['max_abs_closure']:.3e}"
        )
    if reaction_metrics["max_abs_closure"] > 1e-7:
        hard_failures.append(
            f"reaction elemental closure 超限: {reaction_metrics['max_abs_closure']:.3e}"
        )

    metrics = {
        "registry": registry_metrics,
        "pyrolysis": pyrolysis_metrics,
        "reaction": reaction_metrics,
    }
    artifacts = {
        "phase_2_registry_source": "src/core/elemental_ledger.py",
        "phase_2_gate_source": "src/core/phase_gate_checks.py",
    }
    return GateResult(
        gate_id="gate_20_elemental_closure",
        phase="phase-2",
        passed=not hard_failures,
        hard_failures=hard_failures,
        warnings=warnings,
        metrics=metrics,
        artifacts=artifacts,
    )


def _phase3_algebraic_exchange_invariants() -> dict[str, float]:
    c_b = np.array([2.0, 1.0, 0.5] + [0.0] * (N_GAS - 3), dtype=np.float64)
    c_d = np.array([1.0, 1.5, 0.5] + [0.0] * (N_GAS - 3), dtype=np.float64)
    n_ex = calc_phase_exchange(K_bd=2.5, V_b=0.2, C_b=c_b, C_d=c_d)
    zero = np.zeros(N_GAS, dtype=np.float64)
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
    res_d = res[:N_GAS]
    res_b = res[N_GAS:]
    return {
        "dense_match_max_abs": float(np.max(np.abs(res_d - n_ex))),
        "bubble_match_max_abs": float(np.max(np.abs(res_b + n_ex))),
        "net_zero_max_abs": float(np.max(np.abs(res_d + res_b))),
        "equal_phase_zero_max_abs": float(
            np.max(np.abs(calc_phase_exchange(K_bd=3.0, V_b=0.1, C_b=c_b, C_d=c_b)))
        ),
    }


def _phase3_partition_ledger() -> dict[str, Any]:
    cell = _build_phase0_cell_fixture()
    cell.calc_hydrodynamics()
    cell.calc_exchange()
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
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=S_CHAR,
    )

    res_with = calc_gas_balance_residual(
        N_zu_d=cell.N_zu_d,
        N_rez_d=cell.N_rez_d,
        N_d_in=cell.N_d_in,
        R_gas_d=bundle.R_gas_d,
        N_d=cell.N_d,
        N_ex=cell.N_ex,
        N_zu_b=cell.N_zu_b,
        N_rez_b=cell.N_rez_b,
        N_b_in=cell.N_b_in,
        R_gas_b=bundle.R_gas_b,
        N_b=cell.N_b,
    )
    res_noex = calc_gas_balance_residual(
        N_zu_d=cell.N_zu_d,
        N_rez_d=cell.N_rez_d,
        N_d_in=cell.N_d_in,
        R_gas_d=bundle.R_gas_d,
        N_d=cell.N_d,
        N_ex=np.zeros_like(cell.N_ex),
        N_zu_b=cell.N_zu_b,
        N_rez_b=cell.N_rez_b,
        N_b_in=cell.N_b_in,
        R_gas_b=bundle.R_gas_b,
        N_b=cell.N_b,
    )
    exch_d = res_with[:N_GAS] - res_noex[:N_GAS]
    exch_b = res_with[N_GAS:] - res_noex[N_GAS:]

    net_impl = {
        "R5": float(bundle.net_molar_gas_source_r5),
        "R6": float(bundle.net_molar_gas_source_r6),
        "R7": float(bundle.net_molar_gas_source_r7),
        "R8": float(bundle.net_molar_gas_source_r8),
        "R9": float(bundle.net_molar_gas_source_r9),
        "R10": float(bundle.net_molar_gas_source_r10),
        "R11": float(bundle.net_molar_gas_source_r11),
        "R12": float(bundle.net_molar_gas_source_r12),
    }

    return {
        "N_ex_sum": float(np.sum(cell.N_ex)),
        "N_ex_l1": float(np.linalg.norm(cell.N_ex, ord=1)),
        "dense_match_max_abs": float(np.max(np.abs(exch_d - cell.N_ex))),
        "bubble_match_max_abs": float(np.max(np.abs(exch_b + cell.N_ex))),
        "net_zero_max_abs": float(np.max(np.abs(exch_d + exch_b))),
        "o2_supply_bubble": float(bundle.o2_supply_bubble),
        "o2_supply_dense": float(bundle.o2_supply_dense),
        "o2_budget_bubble": float(bundle.o2_budget_bubble),
        "o2_budget_dense": float(bundle.o2_budget_dense),
        "o2_demand_bubble_raw": float(bundle.o2_demand_bubble_raw),
        "o2_demand_dense_raw": float(bundle.o2_demand_dense_raw),
        "o2_demand_bubble_limited": float(bundle.o2_demand_bubble_limited),
        "o2_demand_dense_limited": float(bundle.o2_demand_dense_limited),
        "o2_slack_bubble": float(bundle.o2_slack_bubble),
        "o2_slack_dense": float(bundle.o2_slack_dense),
        "o2_shortfall_bubble": float(bundle.o2_shortfall_bubble),
        "o2_shortfall_dense": float(bundle.o2_shortfall_dense),
        "o2_transfer_potential_bd": float(bundle.o2_transfer_potential_bd),
        "h2o_supply": float(bundle.h2o_supply),
        "h2o_demand_raw": float(bundle.h2o_demand_raw),
        "h2o_demand_limited": float(bundle.h2o_demand_limited),
        "limit_factor_o2_bubble": float(bundle.limit_factor_o2_bubble),
        "limit_factor_o2_dense": float(bundle.limit_factor_o2_dense),
        "limit_factor_h2o": float(bundle.limit_factor_h2o),
        "net_molar_gas_source_total": float(bundle.net_molar_gas_source_total),
        "net_molar_gas_source_vm": float(bundle.net_molar_gas_source_vm),
        "net_molar_gas_source_r5": float(bundle.net_molar_gas_source_r5),
        "net_molar_gas_source_r6": float(bundle.net_molar_gas_source_r6),
        "net_molar_gas_source_r7": float(bundle.net_molar_gas_source_r7),
        "net_molar_gas_source_r8": float(bundle.net_molar_gas_source_r8),
        "net_molar_gas_source_r9": float(bundle.net_molar_gas_source_r9),
        "net_molar_gas_source_r10": float(bundle.net_molar_gas_source_r10),
        "net_molar_gas_source_r11": float(bundle.net_molar_gas_source_r11),
        "net_molar_gas_source_r12": float(bundle.net_molar_gas_source_r12),
        "net_molar_gas_source_char": float(bundle.net_molar_gas_source_char),
        "net_molar_gas_source_by_impl": net_impl,
        "net_molar_gas_source_by_thesis": map_impl_diag_to_thesis(net_impl),
        "reaction_numbering": reaction_numbering_metadata(),
    }


def gate_30_phase_partition_closure() -> GateResult:
    """Phase 3 gate: phase-local budgets and exchange terms must be self-consistent."""

    hard_failures: list[str] = []
    warnings: list[str] = []

    algebra = _phase3_algebraic_exchange_invariants()
    ledger = _phase3_partition_ledger()

    for key in ("dense_match_max_abs", "bubble_match_max_abs", "net_zero_max_abs", "equal_phase_zero_max_abs"):
        if algebra[key] > 1e-12:
            hard_failures.append(f"algebraic {key} 超限: {algebra[key]:.3e}")

    for key in ("dense_match_max_abs", "bubble_match_max_abs", "net_zero_max_abs"):
        if ledger[key] > 1e-12:
            hard_failures.append(f"actual-cell {key} 超限: {ledger[key]:.3e}")
    if abs(ledger["N_ex_sum"]) > 1e-12:
        hard_failures.append(f"N_ex 总和非零: {ledger['N_ex_sum']:.3e}")

    if ledger["o2_demand_bubble_limited"] - 0.995 * ledger["o2_supply_bubble"] > 1e-9:
        hard_failures.append("bubble 相 O2 limiter 后仍超支")
    if ledger["o2_demand_dense_limited"] - 0.995 * ledger["o2_supply_dense"] > 1e-9:
        hard_failures.append("dense 相 O2 limiter 后仍超支")
    if ledger["h2o_demand_limited"] - 0.995 * ledger["h2o_supply"] > 1e-9:
        hard_failures.append("H2O limiter 后仍超支")

    if ledger["o2_demand_bubble_limited"] - ledger["o2_demand_bubble_raw"] > 1e-12:
        hard_failures.append("bubble 相 O2 limited demand 大于 raw demand")
    if ledger["o2_demand_dense_limited"] - ledger["o2_demand_dense_raw"] > 1e-12:
        hard_failures.append("dense 相 O2 limited demand 大于 raw demand")
    if ledger["h2o_demand_limited"] - ledger["h2o_demand_raw"] > 1e-12:
        hard_failures.append("H2O limited demand 大于 raw demand")

    metrics = {
        "algebraic_exchange": algebra,
        "phase_partition_ledger": ledger,
    }
    artifacts = {
        "phase_3_gate_source": "src/core/phase_gate_checks.py",
        "phase_3_budget_source": "src/core/cell_kinetics.py",
    }
    return GateResult(
        gate_id="gate_30_phase_partition_closure",
        phase="phase-3",
        passed=not hard_failures,
        hard_failures=hard_failures,
        warnings=warnings,
        metrics=metrics,
        artifacts=artifacts,
    )


def _phase4_solid_stream_filters() -> dict[str, Any]:
    sample = np.array([[1.0, 2.0, 3.0, 4.0]], dtype=np.float64)
    recycled = _recycled_solid_stream(sample, frac=0.25)
    propagated = _propagated_solid_stream(sample, frac=0.5)
    return {
        "recycled": recycled.tolist(),
        "propagated": propagated.tolist(),
        "recycle_vm_zero": float(recycled[0, S_VM]) == 0.0,
        "recycle_moisture_zero": float(recycled[0, S_MOISTURE]) == 0.0,
        "propagated_vm_zero": float(propagated[0, S_VM]) == 0.0,
        "propagated_moisture_zero": float(propagated[0, S_MOISTURE]) == 0.0,
    }


def _first_true_suffix(flags: list[bool]) -> int | None:
    for i in range(len(flags)):
        if all(flags[i:]):
            return i
    return None


def _positive_rebound_penalty(values: np.ndarray, scale: float) -> float:
    if values.size <= 1:
        return 0.0
    rebound = np.maximum(np.diff(values), 0.0)
    return float(np.sum(rebound) / max(float(scale), 1e-12))


def _phase4_release_rebound(reactor: Reactor) -> dict[str, float]:
    """Compute rebound using fresh-feed drying/pyrolysis release (Hamel-style extent signal)."""
    if not reactor.cells:
        return {
            "vm_release_rebound": 0.0,
            "moist_release_rebound": 0.0,
            "vm_release_scale": 1.0,
            "moist_release_scale": 1.0,
        }

    vm_release_profile: list[float] = []
    moist_release_profile: list[float] = []
    for cell in reactor.cells:
        cell.calc_hydrodynamics()
        tau = float(cell.geo.dh / max(cell.u_mf, 1e-3))
        m_vm_in = float(
            np.sum(
                np.maximum(
                    cell.m_solid_zu[:, S_VM] + cell.m_solid_rez[:, S_VM] + cell.m_solid_in[:, S_VM],
                    0.0,
                )
            )
        )
        m_moist_in = float(
            np.sum(
                np.maximum(
                    cell.m_solid_zu[:, S_MOISTURE]
                    + cell.m_solid_rez[:, S_MOISTURE]
                    + cell.m_solid_in[:, S_MOISTURE],
                    0.0,
                )
            )
        )
        bundle = calc_drying_pyrolysis_sources(
            tau=tau,
            T=float(cell.T),
            P=float(cell.P),
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
            m_vm_in=m_vm_in,
            m_moist_in=m_moist_in,
            solid_shape=cell.R_solid.shape,
            vm_index=S_VM,
            moisture_index=S_MOISTURE,
        )
        vm_release_profile.append(float(max(-np.sum(bundle.solid_sink[:, S_VM]), 0.0)))
        moist_release_profile.append(float(max(-np.sum(bundle.solid_sink[:, S_MOISTURE]), 0.0)))

    vm_arr = np.asarray(vm_release_profile, dtype=float)
    moist_arr = np.asarray(moist_release_profile, dtype=float)
    vm_scale = max(float(np.max(vm_arr, initial=0.0)), 1e-12)
    moist_scale = max(float(np.max(moist_arr, initial=0.0)), 1e-12)
    return {
        "vm_release_rebound": _positive_rebound_penalty(vm_arr, vm_scale),
        "moist_release_rebound": _positive_rebound_penalty(moist_arr, moist_scale),
        "vm_release_scale": vm_scale,
        "moist_release_scale": moist_scale,
    }


def _phase4_solved_state_extent() -> dict[str, Any]:
    utils = _load_validation_utils()
    case = utils.load_case_LU()
    cfg = utils.build_phase2_htw_lu_freeboard_reactor_config(case)
    # Phase-4 gate audits axial release/recycle morphology. Keep it independent from
    # strict major-Gibbs x0 feasibility so failures stay attributable to phase-4 metrics.
    cfg.vorab_major_gibbs_x0 = False
    reactor = Reactor(cfg)
    solve_kwargs = dict(
        getattr(
            utils,
            "PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS",
            getattr(utils, "PHASE1_HTW_LU_SOLVE_KWARGS", {}),
        )
    )
    result = reactor.solve(**solve_kwargs)

    bot = reactor.cells[0]
    vm_feed = float(np.sum(np.maximum(bot.m_solid_zu[:, S_VM], 0.0)))
    moist_feed = float(np.sum(np.maximum(bot.m_solid_zu[:, S_MOISTURE], 0.0)))

    rows: list[dict[str, Any]] = []
    for i, cell in enumerate(reactor.cells):
        m_vm = float(np.sum(np.maximum(cell.m_solid[:, S_VM], 0.0)))
        m_moist = float(np.sum(np.maximum(cell.m_solid[:, S_MOISTURE], 0.0)))
        m_char = float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0)))
        reactive_total = m_char + m_vm + m_moist
        char_share = m_char / max(reactive_total, 1e-12)
        vm_done = m_vm <= 0.01 * max(vm_feed, 1e-12)
        moist_done = m_moist <= 0.01 * max(moist_feed, 1e-12)
        char_dominant = bool(char_share >= 0.95 and vm_done and moist_done)
        rows.append(
            {
                "cell": i,
                "xi": float(cell.geo.h_center / cfg.H_bed),
                "T_K": float(cell.T),
                "m_vm_out_kg_s": m_vm,
                "m_moist_out_kg_s": m_moist,
                "m_char_out_kg_s": m_char,
                "char_share_out": char_share,
                "vm_done": vm_done,
                "moist_done": moist_done,
                "char_dominant": char_dominant,
            }
        )

    vm_done_cell = next((row["cell"] for row in rows if row["vm_done"]), None)
    moisture_done_cell = next((row["cell"] for row in rows if row["moist_done"]), None)
    char_suffix = _first_true_suffix([bool(row["char_dominant"]) for row in rows])
    final_profile_metrics = dict(result.get("final_profile_metrics", {}) or {})
    history = list(result.get("history", []))
    latest_history = history[-1] if history else {}
    profile_metrics = final_profile_metrics if final_profile_metrics else latest_history
    release_rebound = _phase4_release_rebound(reactor)

    return {
        "converged": bool(result.get("converged", False)),
        "n_iter": int(result.get("n_iter", 0) or 0),
        "vm_done_cell": vm_done_cell,
        "moisture_done_cell": moisture_done_cell,
        "char_dominant_from_cell": char_suffix,
        # legacy stock-based rebound (kept for diagnostics only)
        "vm_rebound": float(profile_metrics.get("vm_rebound", 0.0) or 0.0),
        "moist_rebound": float(profile_metrics.get("moist_rebound", 0.0) or 0.0),
        # Hamel-style phase-4 gate signal: fresh-feed release rebound
        "vm_release_rebound": float(release_rebound["vm_release_rebound"]),
        "moist_release_rebound": float(release_rebound["moist_release_rebound"]),
        "vm_release_scale": float(release_rebound["vm_release_scale"]),
        "moist_release_scale": float(release_rebound["moist_release_scale"]),
        "rows": rows,
    }


def gate_40_axial_extent_recycle() -> GateResult:
    """Phase 4 gate: fresh-feed release must be lower-bed limited and recycle char/ash only."""

    extent = _phase4_solved_state_extent()
    hard_failures: list[str] = []
    warnings: list[str] = []

    n_cells = len(extent["rows"])
    if extent["vm_done_cell"] is None or extent["vm_done_cell"] >= n_cells - 1:
        hard_failures.append("VM 未在上部 char-only 区域形成前耗尽")
    if extent["moisture_done_cell"] is None or extent["moisture_done_cell"] >= n_cells - 1:
        hard_failures.append("moisture 未在上部 char-only 区域形成前耗尽")
    if extent["char_dominant_from_cell"] is None or extent["char_dominant_from_cell"] >= n_cells - 1:
        hard_failures.append("未形成至少覆盖一个上部 cell 的连续 char-dominant 区域")
    if extent["vm_release_rebound"] > _PHASE4_REBOUND_TOL:
        hard_failures.append(f"vm_release_rebound 超限: {extent['vm_release_rebound']:.3e}")
    if extent["moist_release_rebound"] > _PHASE4_REBOUND_TOL:
        hard_failures.append(f"moist_release_rebound 超限: {extent['moist_release_rebound']:.3e}")
    if extent["vm_rebound"] > _PHASE4_REBOUND_TOL:
        warnings.append(f"vm_rebound(库存口径) 超限: {extent['vm_rebound']:.3e}")
    if extent["moist_rebound"] > _PHASE4_REBOUND_TOL:
        warnings.append(f"moist_rebound(库存口径) 超限: {extent['moist_rebound']:.3e}")

    if extent["char_dominant_from_cell"] is not None:
        for row in extent["rows"][extent["char_dominant_from_cell"]:]:
            if not row["vm_done"] or not row["moist_done"] or not row["char_dominant"]:
                hard_failures.append(
                    f"上部 char-only suffix 不连续: cell={row['cell']} vm_done={row['vm_done']} "
                    f"moist_done={row['moist_done']} char_dominant={row['char_dominant']}"
                )
                break

    stream_filters = _phase4_solid_stream_filters()
    if not stream_filters["recycle_vm_zero"]:
        hard_failures.append("recycle solid stream 仍携带 VM")
    if not stream_filters["recycle_moisture_zero"]:
        hard_failures.append("recycle solid stream 仍携带 moisture")
    if not stream_filters["propagated_vm_zero"]:
        hard_failures.append("propagated solid stream 仍携带 VM")
    if not stream_filters["propagated_moisture_zero"]:
        hard_failures.append("propagated solid stream 仍携带 moisture")

    metrics = {
        "extent_audit": extent,
        "solid_stream_filters": stream_filters,
    }
    artifacts = {
        "phase_4_deep_audit_script": "scripts/audit_drying_pyrolysis_extent_lu.py",
        "phase_4_gate_source": "src/core/phase_gate_checks.py",
    }
    return GateResult(
        gate_id="gate_40_axial_extent_recycle",
        phase="phase-4",
        passed=not hard_failures,
        hard_failures=hard_failures,
        warnings=warnings,
        metrics=metrics,
        artifacts=artifacts,
    )
