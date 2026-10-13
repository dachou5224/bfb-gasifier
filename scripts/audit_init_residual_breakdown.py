#!/usr/bin/env python3
"""Init 残差分项诊断：定位 A-tier vs legacy 初值 RMS 偏高的具体方程。"""

from __future__ import annotations

from _lu_harness import (
    LADDER_P2,
    apply_config_overrides,
    build_phase2_config,
    data_path,
    repo_root,
    run_phase2_init,
    write_json,
)

REPO_ROOT = repo_root()

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from src.core.reactor import Reactor
from src.core.species import GAS_SPECIES
from src.solvers.global_nr_solver import (
    _line_search_merit,
    _residual_gas_phase_metrics,
    _residual_group_metrics,
    _rms_norm,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.solvers.nr_indexing import (
    gas_mode_for_nr,
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_solid_var,
    n_var,
    solid_component_indices_for_nr,
    solid_offset,
)
P2_OVERRIDES = {
    "nr_line_search_gas_phase_split_merit_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
    "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
    "nr_enforce_gas_phase_split_convergence_thesis": True,
    "vorab_bottom_temperature_cap_K_thesis": 1180.0,
    "freeboard_trajectory_coeff_model": "exact_hamel",
}

SOLID_COMP_NAMES = ("char", "vm", "moisture", "ash")


def _equation_labels_for_cell(cell, cell_index: int) -> list[str]:
    """Human-readable labels aligned with packed NR residual vector."""
    labels: list[str] = []
    gas_idx = main_nr_gas_idx_for_cell(cell)
    species = [GAS_SPECIES[int(i)] for i in gas_idx]
    mode = gas_mode_for_nr(cell)
    if mode == "two_phase":
        for sp in species:
            labels.append(f"cell{cell_index}:{sp}:dense_balance")
        for sp in species:
            labels.append(f"cell{cell_index}:{sp}:bubble_balance")
    elif mode == "single":
        for sp in species:
            labels.append(f"cell{cell_index}:{sp}:gas_balance")
    comp_idx = solid_component_indices_for_nr(cell)
    nk = int(cell.solid.n_size_classes)
    for comp_i in comp_idx:
        cname = SOLID_COMP_NAMES[int(comp_i)] if int(comp_i) < len(SOLID_COMP_NAMES) else f"comp{comp_i}"
        for k in range(nk):
            labels.append(f"cell{cell_index}:{cname}[k{k}]:solid_balance")
    if has_temperature_var(cell):
        labels.append(f"cell{cell_index}:T:energy_balance")
    return labels


def _collect_init_residual_breakdown(reactor: Reactor) -> dict:
    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    reactor._apply_all_bc_for_nr()

    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)

    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_raw = np.asarray(F, dtype=np.float64)
    F_hat = F_raw / np.asarray(scale, dtype=np.float64)

    labels: list[str] = []
    for i, cell in enumerate(cells):
        labels.extend(_equation_labels_for_cell(cell, i))
    if len(labels) != len(F_hat):
        raise RuntimeError(f"label count {len(labels)} != residual count {len(F_hat)}")

    group = _residual_group_metrics(F_hat, cells)
    gas_phase = _residual_gas_phase_metrics(F_hat, cells)
    merit = _line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=bool(getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)),
    )

    per_cell: list[dict] = []
    offset = 0
    for i, cell in enumerate(cells):
        nv = n_var(cell)
        block = F_hat[offset : offset + nv]
        block_raw = F_raw[offset : offset + nv]
        ng = n_gas_var(cell)
        ns = n_solid_var(cell)
        gas_block = block[:ng]
        solid_block = block[ng : ng + ns] if ns else np.zeros(0)
        energy = float(block[-1]) if has_temperature_var(cell) else 0.0
        per_cell.append(
            {
                "cell_index": i,
                "cell_type": str(getattr(cell, "cell_type", "bed")),
                "T_K": float(cell.T),
                "rms_scaled": float(np.sqrt(np.mean(block**2))) if block.size else 0.0,
                "max_abs_scaled": float(np.max(np.abs(block))) if block.size else 0.0,
                "gas_rms": float(np.sqrt(np.mean(gas_block**2))) if gas_block.size else 0.0,
                "solid_rms": float(np.sqrt(np.mean(solid_block**2))) if solid_block.size else 0.0,
                "energy_scaled": energy,
                "energy_raw_W": float(block_raw[-1]) if has_temperature_var(cell) else 0.0,
            }
        )
        offset += nv

    order = np.argsort(-np.abs(F_hat))
    top_n = min(25, len(F_hat))
    top_equations = []
    for rank, j in enumerate(order[:top_n], start=1):
        top_equations.append(
            {
                "rank": rank,
                "label": labels[int(j)],
                "F_hat": float(F_hat[int(j)]),
                "F_raw": float(F_raw[int(j)]),
                "scale": float(scale[int(j)]),
            }
        )

    # Group by equation family prefix (species/phase/type)
    family_rms: dict[str, list[float]] = {}
    for lab, val in zip(labels, F_hat):
        parts = lab.split(":")
        if len(parts) >= 3:
            family = ":".join(parts[1:])
        else:
            family = lab
        family_rms.setdefault(family, []).append(float(val))
    family_summary = {
        fam: {
            "count": len(vals),
            "rms_scaled": float(np.sqrt(np.mean(np.square(vals)))),
            "max_abs_scaled": float(np.max(np.abs(vals))),
        }
        for fam, vals in family_rms.items()
    }
    top_families = sorted(
        family_summary.items(),
        key=lambda kv: kv[1]["rms_scaled"],
        reverse=True,
    )[:20]

    return {
        "n_dof": int(len(F_hat)),
        "rms_scaled": float(_rms_norm(F_hat)),
        "merit_scaled": float(merit),
        "max_abs_scaled": float(np.max(np.abs(F_hat))) if F_hat.size else 0.0,
        "gas_rms": group["gas_rms"],
        "solid_rms": group["solid_rms"],
        "energy_rms": group["energy_rms"],
        "gas_combined_rms": gas_phase["gas_combined_rms"],
        "gas_phase_split_rms": gas_phase["gas_phase_split_rms"],
        "gas_max_abs": group["gas_max_abs"],
        "solid_max_abs": group["solid_max_abs"],
        "energy_max_abs": group["energy_max_abs"],
        "gas_phase_split_max_abs": gas_phase["gas_phase_split_max_abs"],
        "component_max_rms": float(max(group["gas_rms"], group["solid_rms"], group["energy_rms"])),
        "per_cell": per_cell,
        "top_equations": top_equations,
        "top_families_by_rms": [{"family": k, **v} for k, v in top_families],
    }


def _run_case(*, a_tier: bool) -> dict:
    cfg = build_phase2_config()
    for k, v in P2_OVERRIDES.items():
        setattr(cfg, k, v)
    cfg.vorab_a_tier_startwert_budget_thesis = bool(a_tier)

    reactor = Reactor(cfg)
    run_phase2_init(reactor, init_strategy="vorabrechnung")
    breakdown = _collect_init_residual_breakdown(reactor)
    breakdown["label"] = "a_tier" if a_tier else "legacy"
    return breakdown


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "init_residual_breakdown_a_tier_vs_legacy.json",
    )
    args = parser.parse_args()

    legacy = _run_case(a_tier=False)
    a_tier = _run_case(a_tier=True)

    def top_delta(key: str) -> list[dict]:
        legacy_map = {row["label"]: row for row in legacy["top_equations"]}
        out = []
        for row in a_tier["top_equations"][:15]:
            lab = row["label"]
            out.append(
                {
                    "label": lab,
                    "a_tier_F_hat": row["F_hat"],
                    "legacy_F_hat": float(legacy_map.get(lab, {}).get("F_hat", 0.0)),
                    "delta_a_minus_legacy": row["F_hat"] - float(legacy_map.get(lab, {}).get("F_hat", 0.0)),
                }
            )
        return out

    payload = {
        "legacy": legacy,
        "a_tier": a_tier,
        "summary_delta": {
            "rms_scaled": a_tier["rms_scaled"] - legacy["rms_scaled"],
            "merit_scaled": a_tier["merit_scaled"] - legacy["merit_scaled"],
            "gas_rms": a_tier["gas_rms"] - legacy["gas_rms"],
            "solid_rms": a_tier["solid_rms"] - legacy["solid_rms"],
            "energy_rms": a_tier["energy_rms"] - legacy["energy_rms"],
            "gas_phase_split_rms": a_tier["gas_phase_split_rms"] - legacy["gas_phase_split_rms"],
            "max_abs_scaled": a_tier["max_abs_scaled"] - legacy["max_abs_scaled"],
        },
        "a_tier_top_equations_vs_legacy": top_delta("top_equations"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=== Init 残差分组 (legacy → A-tier) ===")
    for key in (
        "rms_scaled",
        "merit_scaled",
        "gas_rms",
        "solid_rms",
        "energy_rms",
        "gas_phase_split_rms",
        "max_abs_scaled",
    ):
        print(
            f"  {key:24s}  legacy={legacy[key]:.6f}  a_tier={a_tier[key]:.6f}  "
            f"Δ={a_tier[key]-legacy[key]:+.6f}"
        )

    print("\n=== A-tier top-10 |F_hat| 方程 ===")
    for row in a_tier["top_equations"][:10]:
        print(f"  {row['rank']:2d}. {row['label']:40s}  F_hat={row['F_hat']:+.4f}  raw={row['F_raw']:+.3e}")

    print("\n=== 按方程族 RMS (A-tier top-8) ===")
    for row in a_tier["top_families_by_rms"][:8]:
        print(f"  {row['family']:32s}  rms={row['rms_scaled']:.4f}  max={row['max_abs_scaled']:.4f}")

    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
