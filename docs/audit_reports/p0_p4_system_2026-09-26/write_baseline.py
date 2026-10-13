"""另存 E8 与 P3 可复算包。不覆盖正式 E8 快照。"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_named_f_endpoint_io import (
    SNAPSHOT_PATH,
    load_named_e8_snapshot,
    named_model_source_files,
    q_eject_snapshot,
)

_ROOT = Path(__file__).resolve().parents[3]
_OUT = _ROOT / "data" / "named_p3_baseline_2026-09-26.json"
_P3 = _ROOT / "docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz"
_EXTRA = (
    "src/physics/wirsum_1997_freeboard.py",
    "src/physics/mass_transfer.py",
    "src/kinetics/char_reactions.py",
    "src/core/freeboard/transport.py",
)


def _source_files() -> list[Path]:
    files = list(named_model_source_files())
    files.extend(_ROOT / rel for rel in _EXTRA)
    return sorted(set(files), key=lambda p: p.as_posix())


def _hash(files: list[Path]) -> str:
    h = hashlib.sha256()
    for path in files:
        rel = path.relative_to(_ROOT).as_posix()
        h.update(rel.encode("utf-8"))
        h.update(b"\0")
        h.update(path.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def _state(name: str, stage: int, x: np.ndarray) -> dict:
    nq = replace(freeze_e8_q(Q_W=0.0), stage=stage)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q, rate_scale=1.0)
    blocks = named_endpoint_blocks(f)
    return {
        "name": name,
        "stage": stage,
        "n_unknowns": int(x.size),
        "rebuild_q": {"factory": "freeze_e8_q", "Q_W": 0.0, "stage": stage},
        "x": np.asarray(x, dtype=np.float64).tolist(),
        "F": np.asarray(f, dtype=np.float64).tolist(),
        "blocks": {k: (bool(v) if isinstance(v, (bool, np.bool_)) else float(v)) for k, v in blocks.items()},
        "energy_rows_W": np.asarray(f[-2:], dtype=np.float64).tolist() if x.size >= 46 else [float(f[-1])],
        "P": q_eject_snapshot(x, nq),
        "columns": {
            "numerical_pass": bool(blocks["passed"]),
            "conservation_pass": "not_run",
            "experimental_error": "not_applicable",
            "model_coverage": "named_two_zone_not_wesseling",
        },
    }


def main() -> None:
    if _OUT.resolve() == SNAPSHOT_PATH.resolve():
        raise SystemExit("refusing to overwrite the E8 snapshot")
    e8 = load_named_e8_snapshot()
    x8 = np.asarray(e8["e8"]["x"], dtype=np.float64)
    x3 = np.load(_P3)["x"]
    files = _source_files()
    pkg = {
        "date": "2026-09-26",
        "role": "recomputable numerical baseline",
        "not_experimental_validation": True,
        "not_full_hamel": True,
        "e8_snapshot_path": SNAPSHOT_PATH.relative_to(_ROOT).as_posix(),
        "e8_snapshot_overwritten": False,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "source_hash": _hash(files),
        "source_files": [p.relative_to(_ROOT).as_posix() for p in files],
        "assumptions": [
            "one bed zone plus one freeboard zone",
            "trajectory split into two equal heights; species and energy are not per cell",
            "superficial gas velocity, not Eq. 3.74",
            "lumped char holdup M = m*tau/(1+k*tau)",
            "hamel_extended roof plus fixed wake, not Eq. 3.67",
            "geometric bed gap is not an inert species; its ejection flux is 0",
            "cyclone K frozen",
            "Q_W = 0",
            "constructed coal feed, not CASE_HTW_WESSELING_1",
        ],
        "validation_contract": {
            "first_case": "CASE_HTW_WESSELING_1",
            "second_case": "CASE_HTW_WESSELING_2",
            "basis": "wet mole fraction unless a column is marked dry",
            "temperature": "report exit_temperature_K and exit_temperature_measured_K separately; do not take the closer one",
            "digitization_rel": 0.05,
            "digitization_source": "validation_cases.json outputs._precision",
            "rate_multipliers": "forbidden",
            "relaxed_or_xfail": "not a pass",
            "wesseling_1_gaps": ["bed_material", "bed_particle_diameter_mm", "primary_air split"],
            "wesseling_2_gaps": ["steam_flow_kg_h", "primary_O2 split", "bed_material"],
        },
        "states": [
            _state("E8", 8, x8),
            _state("P3_E11", 11, np.asarray(x3, dtype=np.float64)),
        ],
    }
    _OUT.write_text(json.dumps(pkg, indent=2), encoding="utf-8")
    print("wrote", _OUT)
    for st in pkg["states"]:
        print(st["name"], st["blocks"])


if __name__ == "__main__":
    main()
