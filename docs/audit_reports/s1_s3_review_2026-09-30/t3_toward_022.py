"""从已联合通过的 s=0.02 走向 0.22。失败即停，不把旧状态当成答案。"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("t3", _HERE / "t3_corrected_path.py")
t3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t3)


def main() -> None:
    doc = json.loads((_HERE / "t3_corrected_path.json").read_text())
    base = doc["s_0_02_from_p3"]
    if not base["joint_pass"]:
        raise SystemExit("s=0.02 is not a joint pass")
    x = np.asarray(base["x"], dtype=np.float64)
    rows = []
    for s in (0.08, 0.16, 0.22):
        row = t3._solve(s, x, 320.0)
        rows.append(row)
        out = dict(doc)
        out["toward_022"] = rows
        out["source"] = t3._source_hash()
        (_HERE / "t3_corrected_path.json").write_text(json.dumps(t3._plain(out)), encoding="utf-8")
        if not row["joint_pass"]:
            break
        x = np.asarray(row["x"], dtype=np.float64)
    print("stopped", rows[-1]["s"], rows[-1]["joint_pass"], flush=True)


if __name__ == "__main__":
    main()
