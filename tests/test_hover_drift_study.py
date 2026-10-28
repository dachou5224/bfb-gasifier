"""研究脚本的连续剖面稳定根：排斥界面不算根，已审查的研究数据可复现。

模型自 ADR 0047 起只用「格内常气速的界面稳定点」；本脚本研究连续剖面，用 ADR 0034
的判据自行求根，不能拿模型的判据代替。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from hamel.freeboard import terminal_settling_velocity
from tools.hover_drift_study import roots

_RECORD = Path("reports/t26c-diagnosis/baseline/record.json")
_DRIFT = Path("reports/research-hover-density/drift.json")
_GAS = 1.2
_MU = 2.0e-5
_DENSITY = 1200.0
_PHI = 0.8
_DIAMETER = 2.0e-4


def _column():
    return np.array([0.0, 0.5, 1.0]), np.full(2, _GAS), np.full(2, _MU)


def test_repelling_interface_is_not_a_stable_root():
    """下侧 u_g < v_t、上侧 u_g > v_t 的跃变是排斥的，不能收成稳定根。"""
    settling = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    edges, rho, mu = _column()
    # 剖面整体上穿 v_t，恰好在中点穿过：du_g/dh > 0，不是稳定根。
    def rising(_cell, height):
        return settling + 2.0 * (float(height) - 0.5)

    assert roots(edges, rising, rho, mu, _DIAMETER, _DENSITY, _PHI) == []


def test_attracting_continuous_crossing_is_a_stable_root():
    """u_g 随高度下降、恰好在中点等于 v_t：这是稳定根。"""
    settling = terminal_settling_velocity(_GAS, _MU, _DIAMETER, _DENSITY, _PHI)
    edges, rho, mu = _column()

    def falling(_cell, height):
        return settling + 0.8 * (0.5 - float(height))

    found = roots(edges, falling, rho, mu, _DIAMETER, _DENSITY, _PHI)
    assert len(found) == 1
    assert found[0] == pytest.approx(0.5, abs=1e-6)


@pytest.mark.skipif(not _RECORD.exists() or not _DRIFT.exists(), reason="需要 T26c 的运行记录与研究数据")
def test_reviewed_study_roots_are_reproduced():
    """工况 1 第 1 轮第 3 粒级：0.2875 m 的二次风入口是排斥跃变，不得进入稳定根。"""
    record = json.loads(_RECORD.read_text(encoding="utf-8"))
    column = record["diagnostics"]["rounds"][0]["info"]["sensitivity"]["trajectory"]
    drift = json.loads(_DRIFT.read_text(encoding="utf-8"))
    from tools.hover_drift_study import _profile

    edges, gas, _at = _profile(column)
    rho = np.asarray(column["rho_g"], dtype=float)
    mu = np.asarray(column["mu_g"], dtype=float)
    diameter = float(np.asarray(column["diameter_m"], dtype=float)[2])
    found = roots(edges, gas, rho, mu, diameter, float(column["rho_p"]), float(column["phi_s"]))
    assert found == pytest.approx(drift["initial_roots_m"])
    assert not any(abs(value - 0.2875) <= 1.0e-6 for value in found)
