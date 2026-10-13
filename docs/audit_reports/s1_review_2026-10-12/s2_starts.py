"""Three 60 s starts from target-case initial values. Does not read the certified state."""
import json
import time
from pathlib import Path

import numpy as np

from src.solvers.named_kernel_sources import named_kernel_source_sha256
from src.solvers.wesseling_two_zone import q_on_wesseling1_path, wesseling1_physical_x0
from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_control import certify_controlled_state, solve_named_controlled

_INVENTORY = (14, 15, 37, 38, 39, 46, 47)


def _blocks(x, nq):
    q = outer_P_named(x, nq)
    residual = coupled_residual_named(x, q)
    cert = certify_controlled_state(x, residual, q.hydro)
    return {
        "gate_ratio": cert["gate_ratio"],
        "energy_rows_W": [float(v) for v in residual[-2:]],
        "max_abs_F": float(np.max(np.abs(residual))),
    }


def _run(tag, x0, nq):
    first = _blocks(x0, nq)
    t0 = time.perf_counter()
    out = solve_named_controlled(nq, x0, deadline_s=60.0, temperature_radius_scale=4.0)
    elapsed = time.perf_counter() - t0
    last = _blocks(out["x"], nq)
    restart = out["stats"].get("oxygen_restart") or {}
    return {
        "tag": tag,
        "elapsed_s": elapsed,
        "passed": bool(out["passed"]),
        "reason": out["reason"],
        "outers": out["stats"].get("outers"),
        "rejects": out["stats"].get("rejects"),
        "oxygen_restart": restart.get("status"),
        "F": out["stats"].get("F_calls"),
        "J": out["stats"].get("J_calls"),
        "P": out["stats"].get("P_calls"),
        "first": first,
        "last": last,
    }


def main():
    nq = q_on_wesseling1_path(1.0)
    base = wesseling1_physical_x0(nq)
    warm = base.copy()
    warm[-2] += 100.0
    warm[-1] += 100.0
    scaled = base.copy()
    scaled[list(_INVENTORY)] *= 2.0
    rows = [
        _run("physical", base, nq),
        _run("temperature_plus_100K", warm, nq),
        _run("inventory_times_2", scaled, nq),
    ]
    payload = {"source_sha256": named_kernel_source_sha256(), "starts": rows}
    Path(__file__).with_name("s2_starts.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps([{k: row[k] for k in row if k != "first"} for row in rows], indent=2))


if __name__ == "__main__":
    main()
