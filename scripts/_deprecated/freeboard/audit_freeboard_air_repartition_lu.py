#!/usr/bin/env python3
"""LU 工况下的 total-air repartition 审计。

与 `audit_freeboard_secondary_injection_lu.py` 不同：
- 不额外增加总空气
- 保持 `ER` 给定的总 `O2/N2` 不变
- 只把总空气在 primary / secondary 之间重分配
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES
from _nr_monitor import solve_with_nr_monitor
from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    apply_secondary_air_repartition,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_validation_case_node,
)


TARGET_XI_BED = 0.45
TARGET_XI_SPIKE = 0.60
TARGET_XI_EXIT = 1.00
SECONDARY_AIR_FRACS = (0.00, 0.05, 0.10, 0.15, 0.20, 0.25)


def _ref_at(node: dict, species_key: str, xi_target: float) -> float:
    prof = node["outputs"]["axial_profiles"]
    for x, v in zip(prof["xi"], prof[species_key]):
        if v is not None and abs(float(x) - xi_target) < 1e-12:
            return float(v)
    raise KeyError(f"reference {species_key} missing at xi={xi_target}")


def _interp_profile(x: list[float], y: list[float], x_target: float) -> float:
    if x_target <= x[0]:
        return float(y[0])
    if x_target >= x[-1]:
        return float(y[-1])
    for i in range(1, len(x)):
        if x_target <= x[i]:
            x0, x1 = float(x[i - 1]), float(x[i])
            y0, y1 = float(y[i - 1]), float(y[i])
            if abs(x1 - x0) < 1e-12:
                return y1
            w = (x_target - x0) / (x1 - x0)
            return y0 + w * (y1 - y0)
    return float(y[-1])


def _combined_wet_profiles(reactor: Reactor, result: dict) -> tuple[list[float], dict[str, list[float]], list[float]]:
    x = list(result["axial_xi_reactor"])
    t = list(result["T_profile"])
    wet = {sp: [] for sp in GAS_SPECIES}

    for cell in reactor.cells:
        y = cell._mole_fractions("combined")
        for j, sp in enumerate(GAS_SPECIES):
            wet[sp].append(float(y[j]))

    if result["freeboard_active"]:
        fb_wet = result["freeboard_gas_profiles_wet"]
        for sp in GAS_SPECIES:
            wet[sp].extend(float(v) for v in fb_wet.get(sp, []))
    return x, wet, t


def _score_case(values: dict[str, float], ref: dict[str, float]) -> float:
    return (
        abs(values["O2_060"] - ref["O2_060"]) / max(ref["O2_060"], 1e-3)
        + abs(values["T_060"] - ref["T_060"]) / max(ref["T_060"], 1.0)
        + abs(values["CO_060"] - ref["CO_060"]) / max(ref["CO_060"], 1e-3)
        + abs(values["H2_060"] - ref["H2_060"]) / max(ref["H2_060"], 1e-3)
        + abs(values["CO2_060"] - ref["CO2_060"]) / max(ref["CO2_060"], 1e-3)
        + abs(values["dCO"] - ref["dCO"]) / max(abs(ref["dCO"]), 1e-3)
        + abs(values["dH2"] - ref["dH2"]) / max(abs(ref["dH2"]), 1e-3)
        + abs(values["exit_T"] - ref["exit_T"]) / max(ref["exit_T"], 1.0)
    )


def main() -> int:
    node = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    ref = {
        "O2_060": _ref_at(node, "O2_mol_wet", TARGET_XI_SPIKE),
        "T_060": _ref_at(node, "T_K", TARGET_XI_SPIKE),
        "CO_060": _ref_at(node, "CO_mol_wet", TARGET_XI_SPIKE),
        "H2_060": _ref_at(node, "H2_mol_wet", TARGET_XI_SPIKE),
        "CO2_060": _ref_at(node, "CO2_mol_wet", TARGET_XI_SPIKE),
        "CO_045": _ref_at(node, "CO_mol_wet", TARGET_XI_BED),
        "H2_045": _ref_at(node, "H2_mol_wet", TARGET_XI_BED),
        "exit_T": _ref_at(node, "T_K", TARGET_XI_EXIT),
    }
    ref["dCO"] = ref["CO_060"] - ref["CO_045"]
    ref["dH2"] = ref["H2_060"] - ref["H2_045"]

    print("=" * 188)
    print("LU freeboard total-air repartition audit")
    print("=" * 188)
    print(
        "ref @xi=0.60: "
        f"O2={ref['O2_060']:.3f} T={ref['T_060']:.0f}K CO={ref['CO_060']:.3f} "
        f"CO2={ref['CO2_060']:.3f} H2={ref['H2_060']:.3f} dCO={ref['dCO']:+.3f} dH2={ref['dH2']:+.3f}"
    )
    print("-" * 188)
    print(
        f"{'sec_air_frac':>12} {'O2prim':>8} {'O2sec':>8} {'N2prim':>8} {'N2sec':>8} {'seg':>4} {'xi_seg':>8} "
        f"{'O2@0.60':>8} {'T@0.60':>8} {'CO@0.60':>8} {'CO2@0.60':>9} {'H2@0.60':>8} "
        f"{'dCO':>8} {'dH2':>8} {'Texit':>8} {'score':>8}"
    )
    print("-" * 188)

    for frac in SECONDARY_AIR_FRACS:
        cfg = build_phase2_htw_lu_freeboard_reactor_config()
        total_o2 = float(cfg.O2_feed)
        total_n2 = float(cfg.N2_feed)
        apply_secondary_air_repartition(
            cfg,
            secondary_air_frac=frac,
            injection_xi=cfg.freeboard_secondary_injection_xi,
            secondary_T_K=cfg.T_inlet,
        )

        reactor = Reactor(cfg)
        result, monitor = solve_with_nr_monitor(
            reactor,
            PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
            check_x0=True,
        )
        x, wet, t = _combined_wet_profiles(reactor, result)

        values = {
            "O2_060": _interp_profile(x, wet["O2"], TARGET_XI_SPIKE),
            "T_060": _interp_profile(x, t, TARGET_XI_SPIKE),
            "CO_060": _interp_profile(x, wet["CO"], TARGET_XI_SPIKE),
            "CO2_060": _interp_profile(x, wet["CO2"], TARGET_XI_SPIKE),
            "H2_060": _interp_profile(x, wet["H2"], TARGET_XI_SPIKE),
            "CO_045": _interp_profile(x, wet["CO"], TARGET_XI_BED),
            "H2_045": _interp_profile(x, wet["H2"], TARGET_XI_BED),
            "exit_T": _interp_profile(x, t, TARGET_XI_EXIT),
        }
        values["dCO"] = values["CO_060"] - values["CO_045"]
        values["dH2"] = values["H2_060"] - values["H2_045"]
        score = _score_case(values, ref)
        seg = result["freeboard_secondary_injection_segment"]
        xi_seg = (
            result["freeboard_axial_xi"][seg]
            if isinstance(seg, int) and 0 <= seg < len(result["freeboard_axial_xi"])
            else float("nan")
        )

        print(
            f"{frac:>12.2f} {cfg.O2_feed:>8.3f} {cfg.freeboard_secondary_O2_mol_s:>8.3f} "
            f"{cfg.N2_feed:>8.3f} {cfg.freeboard_secondary_N2_mol_s:>8.3f} {str(seg):>4s} {xi_seg:>8.3f} "
            f"{values['O2_060']:>8.3f} {values['T_060']:>8.1f} {values['CO_060']:>8.3f} {values['CO2_060']:>9.3f} "
            f"{values['H2_060']:>8.3f} {values['dCO']:>8.3f} {values['dH2']:>8.3f} {values['exit_T']:>8.1f} {score:>8.3f} "
            f"wall={float(monitor.get('wall_time_s', float('nan'))):.2f}s x0_ok={bool((monitor.get('x0_sanity') or {}).get('ok', False))}"
        )
        o2_gap = abs((cfg.O2_feed + cfg.freeboard_secondary_O2_mol_s) - total_o2)
        n2_gap = abs((cfg.N2_feed + cfg.freeboard_secondary_N2_mol_s) - total_n2)
        if o2_gap > 1e-9 or n2_gap > 1e-9:
            print(
                f"warning: repartition conservation drift at frac={frac:.2f}: "
                f"dO2={o2_gap:.3e} dN2={n2_gap:.3e}"
            )

    print("-" * 188)
    print("interpretation:")
    print("  - 这条审计保持总空气固定，只改变 primary/secondary split。")
    print("  - 若某个 split 能在不抬高 Texit 的前提下复现 xi≈0.6 的局部 spike/dip，则更符合当前 HTW Sim 1 口径。")
    print("  - 若所有 split 都无法保留 O2 spike，下一步应查同段瞬时完全混合/反应假设，而不是继续加总空气。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
