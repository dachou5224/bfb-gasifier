"""气相组成换算（湿基 / 干基）。"""

from __future__ import annotations


def wet_to_dry_mole_fractions(exit_gas: dict) -> dict[str, float]:
    """湿基摩尔分数字典 → 干基（扣除 H2O）。

    与文献/JSON 中常以干基给出的出口组成对齐；模型内部守恒为湿基。
    """
    h2o = float(exit_gas.get("H2O", 0.0))
    den = 1.0 - h2o
    if den <= 1e-12:
        return {}
    out: dict[str, float] = {}
    for k, v in exit_gas.items():
        if k in ("H2O", "TAR1", "TAR2") or str(k).startswith("_"):
            continue
        try:
            out[str(k)] = float(v) / den
        except (TypeError, ValueError):
            pass
    return out
