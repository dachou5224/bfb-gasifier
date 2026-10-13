"""从已复核的水力尺度继续，命运保持 0。失败不覆盖最后通过点。"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import (
    _e8_log_scaled_correct,
    blend_hydro_named_q,
    freeze_e8_q,
    residual_e7_plus_scaled_e8,
)

ART = Path(__file__).resolve().parent
SRC = ART / "bfb_p2_hydro_march.npz"


def attempt(z, s, nq7, nq8, t0, n_outer=12):
    x = np.asarray(z, dtype=float).copy()
    last_h = None
    for p in range(n_outer):
        q = blend_hydro_named_q(
            outer_P_named(x, nq7),
            outer_P_named(x, nq8),
            hydro_scale=s,
            fate_scale=0.0,
        )
        x = _e8_log_scaled_correct(q, x)
        h = residual_e7_plus_scaled_e8(
            x, nq7, nq8, hydro_scale=s, fate_scale=0.0
        )
        b = named_endpoint_blocks(h)
        last_h = h
        print(
            f"  s={s:.4f} p={p} gas={b['gas_l2']:.3e} solid={b['solid_l2']:.3e} "
            f"tar={b['tar_l2']:.3e} E={b['energy_W']:.3e} pass={b['passed']} "
            f"dt={time.perf_counter()-t0:.0f}s",
            flush=True,
        )
        if b["passed"]:
            return x, h, True
    return np.asarray(z, dtype=float), last_h, False


def main() -> None:
    z = np.load(SRC)["x"]
    s_ok = float(np.load(SRC)["hydro_scale"])
    nq7, nq8 = freeze_e7_q(), freeze_e8_q()
    h0 = residual_e7_plus_scaled_e8(z, nq7, nq8, hydro_scale=s_ok, fate_scale=0.0)
    b0 = named_endpoint_blocks(h0)
    print("restart", s_ok, "T", float(z[-1]), "pass", b0["passed"], b0, flush=True)
    if not b0["passed"]:
        raise SystemExit("saved state failed block gates")
    ds = 0.0008
    t0 = time.perf_counter()
    history = [{"hydro_scale": s_ok, "passed": True, "T_K": float(z[-1])}]
    while s_ok < 1.0 - 1e-12:
        target = min(1.0, s_ok + ds)
        x, h, ok = attempt(z, target, nq7, nq8, t0)
        b = named_endpoint_blocks(h)
        if ok:
            z = x
            s_ok = target
            np.savez(
                SRC, x=z, F=h, hydro_scale=s_ok, fate_scale=0.0,
            )
            history.append(
                {
                    "hydro_scale": s_ok,
                    "passed": True,
                    "gas_l2": b["gas_l2"],
                    "solid_l2": b["solid_l2"],
                    "tar_l2": b["tar_l2"],
                    "energy_W": b["energy_W"],
                    "T_K": float(z[-1]),
                }
            )
            print("ACCEPTED", s_ok, flush=True)
            ds = min(0.0008, ds * 1.5)
        else:
            print("FAIL", target, flush=True)
            history.append({"hydro_scale": target, "passed": False})
            ds *= 0.5
            if ds < 1.0e-4:
                print("WALL", s_ok, flush=True)
                break
    (ART / "bfb_p2_hydro_march_history.json").write_text(
        json.dumps(history, indent=2), encoding="utf-8"
    )
    print("landed", s_ok, "T", float(z[-1]), flush=True)


if __name__ == "__main__":
    main()
