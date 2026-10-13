"""失败点上比较冻结方向和完整复合方向，并拆开热解与水力。"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_e8 import (
    _e8_log_scaled_correct,
    freeze_e8_q,
    residual_block_merit,
)

_ROOT = Path(__file__).resolve().parents[3]
_X = json.loads(
    (_ROOT / "docs/audit_reports/r1_r4_review_2026-09-30/s1_regression.json").read_text()
)["cases"][0]["x"]
_OUT = Path(__file__).resolve().parent / "s1_direction_split.json"


def _merit(x, q) -> float:
    return residual_block_merit(coupled_residual_named(x, q, rate_scale=1.0))


def _pyro_only(q_old, q_new):
    sp_old = q_old.energy.species
    sp_new = q_new.energy.species
    sp = replace(
        sp_new,
        K_eject=np.asarray(sp_old.K_eject, dtype=np.float64).copy(),
        w_leave=np.asarray(sp_old.w_leave, dtype=np.float64).copy(),
        w_return=np.asarray(sp_old.w_return, dtype=np.float64).copy(),
        geo_fb=sp_old.geo_fb,
        geo_ex=sp_old.geo_ex,
        k=sp_old.k,
        D_g=float(sp_old.D_g),
        rho_g_fb_ref=float(sp_old.rho_g_fb_ref),
    )
    return replace(q_old, energy=replace(q_new.energy, species=sp), hydro=q_old.hydro)


def _hydro_only(q_old, q_new):
    sp_old = q_old.energy.species
    sp_new = q_new.energy.species
    sp = replace(
        sp_old,
        K_eject=np.asarray(sp_new.K_eject, dtype=np.float64).copy(),
        w_leave=np.asarray(sp_new.w_leave, dtype=np.float64).copy(),
        w_return=np.asarray(sp_new.w_return, dtype=np.float64).copy(),
        geo_fb=sp_new.geo_fb,
        geo_ex=sp_new.geo_ex,
        D_g=float(sp_new.D_g),
        rho_g_fb_ref=float(sp_new.rho_g_fb_ref),
    )
    return replace(q_old, energy=replace(q_old.energy, species=sp), hydro=q_new.hydro)


def main() -> None:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
    x = np.asarray(_X, dtype=np.float64)
    q = outer_P_named(x, nq)
    base = _merit(x, q)
    stats: dict = {}
    proposed = _e8_log_scaled_correct(
        q, x, max_nfev=8, max_j=4, ftol=1.0e-4, temperature_factor=1.5,
        log_span=4.0, stats=stats,
    )
    step = proposed - x
    rows = []
    for alpha in (0.0, 1.0e-4, 0.25, 1.0):
        y = x + alpha * step
        qy = outer_P_named(y, nq)
        rows.append({
            "alpha": alpha,
            "frozen": _merit(y, q),
            "composite": _merit(y, qy),
            "pyro_only": _merit(y, _pyro_only(q, qy)),
            "hydro_only": _merit(y, _hydro_only(q, qy)),
        })
    slope = {}
    for key in ("frozen", "composite", "pyro_only", "hydro_only"):
        slope[key] = (rows[1][key] - base) / 1.0e-4
    out = {"merit0": base, "rows": rows, "slope_per_alpha": slope, "inner_status": stats.get("inner_status")}
    _OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
