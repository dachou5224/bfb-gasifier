"""LU 工况 char-conversion / char-kinetics 审计。

目的
----
1. 对比当前实现使用的局部 `Xchar_local = 1 - m_char,out / m_char,in`
2. 对比自底部 fresh char 参考得到的累计 char conversion
3. 量化两种 `Xchar` 口径对异相炭反应（R1/R2/R3/R4）的影响

这个脚本不修改模型，只输出可复现诊断，便于后续决定是修
`d_core/Xchar` 定义，还是继续沿着 R11 / tar / O2 链路调参。

Usage
-----
    python3 scripts/audit_char_conversion_profile.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import S_CHAR
from src.core.constants import Rg
from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES_INDEX, gas_diffusivity_correlation
from src.kinetics.char_reactions import (
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
from tests.validation_case_utils import (
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_reactor_config,
)


def _char_reaction_bundle(cell, x_char: float) -> dict[str, float]:
    idx = GAS_SPECIES_INDEX
    cell.calc_hydrodynamics()
    c_d = cell._concentrations("d")
    d_g = gas_diffusivity_correlation(cell.T, cell.P)
    total_area = float(np.sum(cell._calc_char_surface_area_per_class()))
    if total_area <= 0.0:
        return {"r1_o2": 0.0, "r2": 0.0, "r3": 0.0, "r4": 0.0}

    d_core = d_core_from_spm_char_conversion(float(x_char), float(cell.solid.d_p))
    r1, alpha = rate_R1(
        cell.T,
        c_d[idx["O2"]],
        float(cell.solid.d_p),
        d_g,
        d_core,
        cell.fuel_type,
    )
    r2 = rate_R2(cell.T, c_d[idx["H2O"]], float(cell.solid.d_p), d_g, d_core)
    r3 = rate_R3(cell.T, c_d[idx["H2"]], float(cell.solid.d_p), d_g, d_core)
    p_co2 = max(float(c_d[idx["CO2"]] * Rg * cell.T), 0.0)
    p_co = max(float(c_d[idx["CO"]] * Rg * cell.T), 0.0)
    r4 = rate_R4_effective(cell.T, p_co2, p_co, float(cell.solid.d_p), d_g, d_core)
    return {
        "r1_o2": float(alpha * r1 * total_area),
        "r2": float(r2 * total_area),
        "r3": float(r3 * total_area),
        "r4": float(r4 * total_area),
    }


def main() -> int:
    cfg = build_phase1_htw_lu_reactor_config()
    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    fresh_char_ref = float(np.sum(np.maximum(reactor.cells[0].m_solid_zu[:, S_CHAR], 0.0)))
    bottom_total_char_in = float(
        np.sum(
            np.maximum(
                reactor.cells[0].m_solid_zu[:, S_CHAR]
                + reactor.cells[0].m_solid_rez[:, S_CHAR]
                + reactor.cells[0].m_solid_in[:, S_CHAR],
                0.0,
            )
        )
    )

    print("=" * 118)
    print("LU char-conversion audit")
    print("=" * 118)
    print(
        f"converged={result.get('converged')} n_iter={result.get('n_iter')} "
        f"rms_scaled_gs={result.get('rms_scaled_gs'):.6f} residual_gs={result.get('residual_gs', float('nan')):.3f}"
    )
    print(
        f"T_exit={result['T_profile'][-1]:.2f} K carbon_conv={result['carbon_conv']:.6f} "
        f"fresh_char_ref={fresh_char_ref:.6f} kg/s bottom_total_char_in={bottom_total_char_in:.6f} kg/s"
    )
    print("-" * 118)
    print(
        f"{'cell':>4s} {'T[K]':>8s} {'char_in':>11s} {'char_out':>11s} "
        f"{'X_local':>9s} {'X_cum':>9s} {'R1_old':>10s} {'R1_cum':>10s} "
        f"{'R2_old':>10s} {'R2_cum':>10s} {'R4_old':>10s} {'R4_cum':>10s}"
    )
    print("-" * 118)

    totals = {
        "r1_old": 0.0,
        "r1_cum": 0.0,
        "r2_old": 0.0,
        "r2_cum": 0.0,
        "r3_old": 0.0,
        "r3_cum": 0.0,
        "r4_old": 0.0,
        "r4_cum": 0.0,
    }

    for i, cell in enumerate(reactor.cells):
        char_in = float(
            np.sum(
                np.maximum(
                    cell.m_solid_zu[:, S_CHAR]
                    + cell.m_solid_in[:, S_CHAR]
                    + cell.m_solid_rez[:, S_CHAR],
                    0.0,
                )
            )
        )
        char_out = float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0)))
        x_local = float(cell._compute_char_conversion())
        x_cum = float(np.clip(1.0 - char_out / max(fresh_char_ref, 1e-12), 0.0, 1.0))

        bundle_old = _char_reaction_bundle(cell, x_local)
        bundle_cum = _char_reaction_bundle(cell, x_cum)

        totals["r1_old"] += bundle_old["r1_o2"]
        totals["r1_cum"] += bundle_cum["r1_o2"]
        totals["r2_old"] += bundle_old["r2"]
        totals["r2_cum"] += bundle_cum["r2"]
        totals["r3_old"] += bundle_old["r3"]
        totals["r3_cum"] += bundle_cum["r3"]
        totals["r4_old"] += bundle_old["r4"]
        totals["r4_cum"] += bundle_cum["r4"]

        print(
            f"{i:4d} {cell.T:8.1f} {char_in:11.6f} {char_out:11.6f} "
            f"{x_local:9.4f} {x_cum:9.4f} {bundle_old['r1_o2']:10.4f} {bundle_cum['r1_o2']:10.4f} "
            f"{bundle_old['r2']:10.4f} {bundle_cum['r2']:10.4f} {bundle_old['r4']:10.4f} {bundle_cum['r4']:10.4f}"
        )

    print("-" * 118)
    print(
        "totals: "
        f"R1_old={totals['r1_old']:.4f} R1_cum={totals['r1_cum']:.4f} "
        f"R2_old={totals['r2_old']:.4f} R2_cum={totals['r2_cum']:.4f} "
        f"R3_old={totals['r3_old']:.6f} R3_cum={totals['r3_cum']:.6f} "
        f"R4_old={totals['r4_old']:.4f} R4_cum={totals['r4_cum']:.4f}"
    )
    print(
        "note: X_cum is defined here as 1 - m_char(cell) / fresh_char_feed(bottom), "
        "so it monotonically tracks axial cumulative char depletion."
    )
    print("=" * 118)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
