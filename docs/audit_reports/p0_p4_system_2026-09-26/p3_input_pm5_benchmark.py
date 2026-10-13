"""进气或燃料 ±5% 的热启动。从 P3 端点出发，每项 120 s。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q
from src.solvers.wirsum_1997_pyrolysis_ledger import dry_mass_from_cho_atoms

_ROOT = Path(__file__).resolve().parents[3]
_X = np.load(_ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz")["x"]
_OUT = _ROOT / "docs/audit_reports/p0_p4_system_2026-09-26/p3_input_pm5_benchmark.json"
_BUDGET_S = 120.0


def _base():
    return replace(freeze_e8_q(Q_W=0.0), stage=11)


def _air(factor: float):
    nq = _base()
    sp = nq.energy.species
    sp2 = replace(sp, N_in_b=np.asarray(sp.N_in_b, dtype=np.float64) * factor)
    return replace(nq, energy=replace(nq.energy, species=sp2))


def _fuel(factor: float):
    nq = _base()
    feed = nq.energy.feed
    scaled = replace(
        feed,
        nC_fuel=float(feed.nC_fuel) * factor,
        nH_fuel=float(feed.nH_fuel) * factor,
        nO_fuel=float(feed.nO_fuel) * factor,
        nN_fuel=float(feed.nN_fuel) * factor,
        nS_fuel=float(feed.nS_fuel) * factor,
        nC_fixed=float(feed.nC_fixed) * factor,
        nC_vm=float(feed.nC_vm) * factor,
        m_moist_kg_s=float(feed.m_moist_kg_s) * factor,
    )
    scaled = replace(scaled, m_dry_kg_s=dry_mass_from_cho_atoms(scaled))
    return replace(nq, energy=replace(nq.energy, feed=scaled))


def _one(name: str, nq) -> dict:
    x = np.asarray(_X, dtype=np.float64).copy()
    x_out = x
    t0 = time.monotonic()
    status = "timeout"
    blocks = None
    worst = None
    outers = 0
    try:
        while time.monotonic() - t0 < _BUDGET_S:
            remain = _BUDGET_S - (time.monotonic() - t0)
            q = outer_P_named(x, nq)
            x = _e8_log_scaled_correct(q, x, deadline_s=remain)
            outers += 1
            q = outer_P_named(x, nq)
            f = coupled_residual_named(x, q, rate_scale=1.0)
            blocks = named_endpoint_blocks(f)
            i = int(np.argmax(np.abs(f)))
            worst = {"index": i, "value": float(f[i])}
            x_out = x
            if bool(blocks["passed"]):
                status = "passed"
                break
            status = "failed"
    except TimeoutError:
        status = "timeout" if blocks is None else "failed"
    elapsed = time.monotonic() - t0
    if elapsed > _BUDGET_S and status != "passed":
        status = "timeout"
    row = {
        "name": name,
        "status": status,
        "elapsed_s": elapsed,
        "budget_s": _BUDGET_S,
        "outer_refreshes": outers,
        "blocks": None if blocks is None else {
            k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in blocks.items()
        },
        "worst_row": worst,
        "T_bed": float(x_out[-2]),
        "T_fb": float(x_out[-1]),
    }
    print(name, status, round(elapsed, 2), flush=True)
    return row


def main() -> None:
    cases = [
        ("air_plus_5pct", _air(1.05)),
        ("air_minus_5pct", _air(0.95)),
        ("fuel_plus_5pct", _fuel(1.05)),
        ("fuel_minus_5pct", _fuel(0.95)),
    ]
    rows = [_one(name, nq) for name, nq in cases]
    _OUT.write_text(json.dumps({"cases": rows}, indent=2), encoding="utf-8")
    print("wrote", _OUT, flush=True)


if __name__ == "__main__":
    main()
