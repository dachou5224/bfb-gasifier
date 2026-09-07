"""冻结 Wirsum 对齐审计的 P0 数值 baseline（不改模型公式）。

采集：
- sanity_checks
- 可选 pytest 摘要
- reduced 3-cell init 残差 / Jacobian
- Phase1 HTW LU（现有最稳 harness）solve 指标
- Phase2 reduced freeboard smoke（3 bed + 2 fb，2 inner iter）

输出 ``data/wirsum_alignment_p0_baseline.json``，供后续 commit 做机器对比。

不追 full-LU 收敛。Jacobian 条件数仅在 reduced 3-cell 上用 dense 2-norm 估计。
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any

_SCRIPTS = Path(__file__).resolve().parent
_ROOT = _SCRIPTS.parent
for _p in (_ROOT, _SCRIPTS):
    _ps = str(_p)
    if _ps not in sys.path:
        sys.path.insert(0, _ps)

import numpy as np

from _lu_harness import add_json_output_arg, repo_root, write_json

SCHEMA_VERSION = 1
EXIT_SPECIES = ("H2", "CO", "CO2", "CH4")
KNOB_KEYS = (
    "thesis_mode",
    "extent_limiters_enabled_thesis",
    "r2_scale",
    "r4_scale",
    "r5_scale",
    "r6_scale",
    "r7_scale",
    "r8_scale",
    "gas_inlet_split_strategy",
    "bed_solid_state_model",
    "freeboard_trajectory_model",
    "freeboard_trajectory_coeff_model",
    "freeboard_velocity_sigma",
    "freeboard_velocity_bins",
    "n_cells",
    "n_freeboard_cells",
)


def _git_head() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root(),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    if out.returncode != 0:
        return None
    return str(out.stdout).strip() or None


def _finite_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _species_subset(mapping: Any) -> dict[str, float | None]:
    if not isinstance(mapping, dict):
        return {sp: None for sp in EXIT_SPECIES}
    return {sp: _finite_float(mapping.get(sp)) for sp in EXIT_SPECIES}


def snapshot_config_knobs(cfg: Any) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in KNOB_KEYS:
        if not hasattr(cfg, key):
            continue
        value = getattr(cfg, key)
        if isinstance(value, (bool, str, int)) or value is None:
            out[key] = value
        else:
            out[key] = _finite_float(value)
    return out


def _condition_estimate(J: Any, *, zero_rows: int, zero_cols: int) -> dict[str, Any]:
    if int(zero_rows) > 0 or int(zero_cols) > 0:
        return {
            "condition_estimate": None,
            "condition_method": "skipped_singular_zero_row_or_col",
        }
    n = int(J.shape[0])
    if n <= 0:
        return {"condition_estimate": None, "condition_method": "empty"}
    if n > 500:
        return {
            "condition_estimate": None,
            "condition_method": "skipped_n_gt_500",
            "n_var": n,
        }
    dense = J.toarray()
    with np.errstate(divide="ignore", invalid="ignore"):
        cond = float(np.linalg.cond(dense, p=2))
    if not np.isfinite(cond):
        return {
            "condition_estimate": None,
            "condition_method": "dense_2norm_nonfinite",
            "n_var": n,
        }
    return {
        "condition_estimate": cond,
        "condition_method": "dense_2norm",
        "n_var": n,
    }


def _packed_residual(reactor: Any, *, include_jacobian: bool = True) -> dict[str, Any]:
    from src.solvers.global_nr_solver import (
        build_equation_scales,
        build_jacobian_fd,
        global_residual,
        pack_reactor,
    )
    from src.solvers.nr_indexing import n_var

    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    x = pack_reactor(cells)
    F = global_residual(x, cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=float(cfg.fuel_feed),
        ref_energy_W=float(cfg.fuel_feed) * 20e6,
    )
    F_hat = F / eq_scale
    out = {
        "n_solver_cells": int(len(cells)),
        "n_var": int(sum(n_var(c) for c in cells)),
        "residual_l2": _finite_float(np.linalg.norm(F, ord=2)),
        "residual_linf": _finite_float(np.linalg.norm(F, ord=np.inf)),
        "rms_scaled": _finite_float(np.sqrt(np.mean(np.square(F_hat)))),
    }
    if not include_jacobian:
        out.update(
            {
                "jacobian_nnz": None,
                "jacobian_zero_rows": None,
                "jacobian_zero_cols": None,
                "jacobian_strategy": None,
                "condition_estimate": None,
                "condition_method": "skipped",
            }
        )
        return out
    J, meta = build_jacobian_fd(
        x,
        F,
        cells,
        reactor._apply_all_bc_for_nr,
        eq_scale,
        strategy="block_tridiag_structured",
    )
    zero_rows = int(meta.get("zero_rows", 0))
    zero_cols = int(meta.get("zero_cols", 0))
    cond = _condition_estimate(J, zero_rows=zero_rows, zero_cols=zero_cols)
    out.update(
        {
            "jacobian_nnz": int(meta.get("nnz", J.nnz)),
            "jacobian_zero_rows": zero_rows,
            "jacobian_zero_cols": zero_cols,
            "jacobian_strategy": str(meta.get("strategy", "block_tridiag_structured")),
            **cond,
        }
    )
    return out


def _exit_from_cells(reactor: Any) -> dict[str, Any]:
    from src.core.composition import wet_to_dry_mole_fractions
    from src.core.species import GAS_SPECIES

    top = reactor.cells[-1]
    y_wet = top._mole_fractions("combined")
    wet = {sp: float(y_wet[j]) for j, sp in enumerate(GAS_SPECIES)}
    dry = wet_to_dry_mole_fractions(wet)
    return {
        "exit_T_K": _finite_float(top.T),
        "exit_gas_wet": _species_subset(wet),
        "exit_gas_dry": _species_subset(dry),
        "bed_T_K": [_finite_float(c.T) for c in reactor.cells],
    }


def _exit_from_solve_result(result: dict[str, Any]) -> dict[str, Any]:
    t_profile = result.get("T_profile") or []
    counts = result.get("nr_counts") or {}
    return {
        "converged": bool(result.get("converged")),
        "converged_inner_nr": result.get("converged_inner_nr"),
        "converged_outer": result.get("converged_outer"),
        "n_iter": int(result.get("n_iter") or 0),
        "nr_outer_iters": result.get("nr_outer_iters"),
        "residual_l2": _finite_float(result.get("residual")),
        "rms_scaled": _finite_float(result.get("rms_scaled_final")),
        "jacobian_nnz": int(counts.get("jacobian_nnz_last") or 0),
        "jacobian_zero_rows": int(counts.get("jacobian_zero_rows_last") or 0),
        "jacobian_zero_cols": int(counts.get("jacobian_zero_cols_last") or 0),
        "exit_T_K": _finite_float(t_profile[-1] if t_profile else None),
        "exit_gas_wet": _species_subset(result.get("exit_gas")),
        "exit_gas_dry": _species_subset(result.get("exit_gas_dry")),
    }


def capture_reduced_3cell_init(*, include_jacobian: bool = True) -> dict[str, Any]:
    """现有 Jacobian 测试用的 3-cell LU init（无 NR）。"""
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    t0 = perf_counter()
    reactor = build_initialized_lu_reactor(n_cells=3)
    packed = _packed_residual(reactor, include_jacobian=include_jacobian)
    packed.update(_exit_from_cells(reactor))
    packed["wall_s"] = float(perf_counter() - t0)
    packed["knobs"] = snapshot_config_knobs(reactor.config)
    packed["case"] = "reduced_3cell_lu_init"
    packed["note"] = "Vorabrechnung x0 only; no Newton iterations."
    return packed


def capture_phase1_lu_solve(*, max_global_iter: int = 3) -> dict[str, Any]:
    """Phase1 HTW LU 短跑。P0 只冻指标，不追 150-iter 收敛。"""
    from src.core.reactor import Reactor
    from tests.validation_case_utils import (
        PHASE1_HTW_LU_SOLVE_KWARGS,
        build_phase1_htw_lu_reactor_config,
    )

    t0 = perf_counter()
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    solve_kwargs = dict(PHASE1_HTW_LU_SOLVE_KWARGS)
    solve_kwargs["max_global_iter"] = int(max_global_iter)
    result = reactor.solve(**solve_kwargs)
    out = _exit_from_solve_result(result)
    out["wall_s"] = float(perf_counter() - t0)
    out["knobs"] = snapshot_config_knobs(cfg)
    out["case"] = "phase1_htw_lu_solve"
    out["solve_kwargs"] = solve_kwargs
    out["full_harness_max_global_iter"] = int(PHASE1_HTW_LU_SOLVE_KWARGS["max_global_iter"])
    out["note"] = (
        "Short freeze of the existing LU harness. Full test kwargs keep "
        f"max_global_iter={PHASE1_HTW_LU_SOLVE_KWARGS['max_global_iter']}; "
        "P0 does not chase that budget."
    )
    return out


def capture_phase2_reduced_freeboard_smoke() -> dict[str, Any]:
    """现有 3+2 freeboard smoke（2 个 global iter），不追 full LU。"""
    from tests.global_nr_fixtures import run_freeboard_trajectory_smoke

    t0 = perf_counter()
    result = run_freeboard_trajectory_smoke(
        trajectory_model="analytical_wirsum",
        coeff_model="exact_hamel",
    )
    out = _exit_from_solve_result(result)
    out["wall_s"] = float(perf_counter() - t0)
    out["case"] = "phase2_reduced_freeboard_smoke"
    out["solve_kwargs"] = {
        "max_global_iter": 2,
        "freeboard_trajectory_model": "analytical_wirsum",
        "freeboard_trajectory_coeff_model": "exact_hamel",
    }
    out["note"] = (
        "Reduced Phase2 smoke only. Name analytical_wirsum currently mixes "
        "Wirsum/Hamel/Haider; split in later P3 commits."
    )
    return out


def run_sanity_checks() -> dict[str, Any]:
    t0 = perf_counter()
    proc = subprocess.run(
        [sys.executable, "tests/sanity_checks.py"],
        cwd=repo_root(),
        capture_output=True,
        text=True,
    )
    stdout = proc.stdout or ""
    failed: list[str] = []
    current = None
    for line in stdout.splitlines():
        if line.startswith("[") and line.endswith("]"):
            current = line.strip("[]")
        elif current is not None and (line.strip().startswith("FAIL:") or line.strip().startswith("ERROR:")):
            failed.append(current)
    return {
        "returncode": int(proc.returncode),
        "passed": proc.returncode == 0,
        "n_failed": int(len(failed)),
        "failed_checks": failed,
        "wall_s": float(perf_counter() - t0),
        "stdout_tail": "\n".join(stdout.strip().splitlines()[-12:]),
        "stderr_tail": "\n".join((proc.stderr or "").strip().splitlines()[-8:]),
    }


def run_pytest_summary() -> dict[str, Any]:
    env = dict(os.environ)
    env.setdefault("PYTHONHASHSEED", "0")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    t0 = perf_counter()
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--tb=line"],
        cwd=repo_root(),
        capture_output=True,
        text=True,
        env=env,
    )
    lines = [ln.strip() for ln in (proc.stdout or "").splitlines() if ln.strip()]
    return {
        "returncode": int(proc.returncode),
        "passed": proc.returncode == 0,
        "wall_s": float(perf_counter() - t0),
        "summary_line": lines[-1] if lines else "",
        "stderr_tail": "\n".join((proc.stderr or "").strip().splitlines()[-12:]),
    }


def _log(msg: str) -> None:
    print(msg, flush=True)


def build_baseline_payload(
    *,
    include_phase1_solve: bool = True,
    include_phase2_smoke: bool = True,
    include_pytest: bool = False,
    phase1_max_global_iter: int = 3,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "description": (
            "Wirsum 1997 alignment P0 freeze. Numbers are the current engineering "
            "stack, not a claim of thesis-exact Wirsum Newton."
        ),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": _git_head(),
        "pythonhashseed": os.environ.get("PYTHONHASHSEED"),
        "sanity_checks": {},
        "cases": {},
    }
    _log("P0 freeze: sanity_checks")
    payload["sanity_checks"] = run_sanity_checks()
    _log("P0 freeze: reduced_3cell_lu_init")
    payload["cases"]["reduced_3cell_lu_init"] = capture_reduced_3cell_init()
    if include_phase1_solve:
        _log(f"P0 freeze: phase1_htw_lu_solve max_global_iter={phase1_max_global_iter}")
        payload["cases"]["phase1_htw_lu_solve"] = capture_phase1_lu_solve(
            max_global_iter=phase1_max_global_iter
        )
    if include_phase2_smoke:
        _log("P0 freeze: phase2_reduced_freeboard_smoke")
        payload["cases"]["phase2_reduced_freeboard_smoke"] = capture_phase2_reduced_freeboard_smoke()
    if include_pytest:
        _log("P0 freeze: pytest -q")
        payload["pytest"] = run_pytest_summary()
    else:
        payload["pytest"] = {
            "skipped": True,
            "note": "Run with --run-pytest or merge pytest -q summary later.",
        }
    return payload


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_json_output_arg(parser, "wirsum_alignment_p0_baseline.json")
    parser.add_argument("--skip-phase1-solve", action="store_true")
    parser.add_argument("--skip-phase2-smoke", action="store_true")
    parser.add_argument("--run-pytest", action="store_true")
    parser.add_argument("--phase1-max-global-iter", type=int, default=3)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    payload = build_baseline_payload(
        include_phase1_solve=not bool(args.skip_phase1_solve),
        include_phase2_smoke=not bool(args.skip_phase2_smoke),
        include_pytest=bool(args.run_pytest),
        phase1_max_global_iter=int(args.phase1_max_global_iter),
    )
    out = write_json(args.output, payload)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
