"""具名 F_E8 端点快照：x / P(x) 字段 / F / 源码 hash。

JSON 里的 ``endpoint_certified`` 只是当时现场记录，不能当通过证据。
加载后必须再跑 ``certify_named_full_rate_endpoint``。

Source: docs/audit_reports/named_f_e8_review_2026-09-16.md §5
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from src.solvers.wirsum_1997_named_f import NamedStageQ, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import certify_named_full_rate_endpoint
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q

_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_PATH = _ROOT / "data" / "named_e8_endpoint_2026-09-17.json"

_HASH_EXTRAS: tuple[str, ...] = (
    "src/kinetics/gas_reactions.py",
    "src/kinetics/tar_reactions.py",
    "src/solvers/wirsum_1997_bed_freeboard.py",
    "src/solvers/wirsum_1997_coupled_f_two_bin_pyro.py",
    "src/solvers/wirsum_1997_pyrolysis_ledger.py",
    "src/thermal/hegermann_yield.py",
    "src/solvers/wirsum_1997_startwert_e7.py",
)


def named_model_source_files() -> list[Path]:
    """影响 F_E8 装配的源文件。不含本 IO 模块。"""
    files = sorted((_ROOT / "src" / "solvers").glob("wirsum_1997_named_f*.py"))
    files = [p for p in files if p.name != "wirsum_1997_named_f_endpoint_io.py"]
    files.extend(_ROOT / rel for rel in _HASH_EXTRAS)
    return files


def named_model_source_hash() -> str:
    h = hashlib.sha256()
    for path in named_model_source_files():
        rel = path.relative_to(_ROOT).as_posix()
        h.update(rel.encode("utf-8"))
        h.update(b"\0")
        h.update(path.read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def q_eject_snapshot(x: npt.NDArray[np.float64], nq: NamedStageQ) -> dict[str, Any]:
    """``P(x)`` 的抛射/气密/旋风 K，不是完整 q 对象。"""
    z = np.asarray(x, dtype=np.float64)
    nqP = outer_P_named(z, nq)
    sp = nqP.energy.species
    k = sp.k
    return {
        "K_eject": np.asarray(sp.K_eject, dtype=np.float64).tolist(),
        "w_leave": np.asarray(sp.w_leave, dtype=np.float64).tolist(),
        "w_return": np.asarray(sp.w_return, dtype=np.float64).tolist(),
        "rho_g": float(sp.geo_fb.rho_g),
        "rho_g_fb_ref": float(sp.rho_g_fb_ref),
        "cyclone_k": {
            "K_auf": float(k.K_auf),
            "K_ab": float(k.K_ab),
            "K_return": float(k.K_return),
            "model": str(k.model),
        },
        "Q_W": float(nqP.Q_W),
        "enabled_rxns": sorted(nqP.enabled_rxns),
    }


def recertify_named_e7(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    cert = certify_named_full_rate_endpoint(z, freeze_e7_q(Q_W=0.0), coupling="composite_P")
    return {k: v for k, v in cert.items() if k != "F"}


def recertify_named_e8(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    nq8 = freeze_e8_q(Q_W=0.0)
    cert = certify_named_full_rate_endpoint(z, nq8, coupling="composite_P")
    out = {k: v for k, v in cert.items() if k != "F"}
    out["F"] = np.asarray(cert["F"], dtype=np.float64)
    return out


def load_named_e8_snapshot(path: Path | None = None) -> dict[str, Any]:
    target = SNAPSHOT_PATH if path is None else Path(path)
    return json.loads(target.read_text(encoding="utf-8"))
