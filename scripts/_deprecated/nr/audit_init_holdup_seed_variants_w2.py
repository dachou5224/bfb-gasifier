#!/usr/bin/env python3
"""W2 init holdup seed 变体对比：阶段 3 子步追踪 + bed0 跳过 reseed 实验。

Usage:
    python3 scripts/audit_init_holdup_seed_variants_w2.py
    python3 scripts/audit_init_holdup_seed_variants_w2.py --case CASE_HTW_WESSELING_2
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable

import numpy as np

from _lu_harness import build_phase2_reactor, data_path, ensure_repo_path, write_json

ensure_repo_path()

from src.core.cell import S_CHAR, S_ASH, S_VM, S_MOISTURE
from src.core.connectivity.solid_transport import (
    seed_holdup_from_inflows,
    update_bed_solid_transport_inflows_from_above,
)
from src.core.reactor import Reactor
from src.solvers.global_nr_solver import (
    _line_search_merit,
    _residual_group_metrics,
    build_equation_scales,
    global_residual,
    pack_reactor,
)
from src.solvers.vorabrechnung import align_vm_devolatilization_zone_vorab_state
from src.workflow.steps import init_precalc_step as init_step

CASE_DEFAULT = "CASE_HTW_WESSELING_2"

STAGE3_SUBSTEPS = (
    "3a_after_bc",
    "3b_after_first_holdup_chain_seed",
    "3c_after_reactive_profile_and_reseed",
)


def _bed_cell_objs(reactor: Reactor) -> list:
    return [c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed"]


def _init_merit(reactor: Reactor) -> dict[str, float]:
    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    reactor._apply_all_bc_for_nr()
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    F = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    F_hat = np.asarray(F, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    group = _residual_group_metrics(F_hat, cells)
    merit = _line_search_merit(
        F_hat,
        cells,
        use_gas_phase_split=bool(getattr(cfg, "nr_line_search_gas_phase_split_merit_thesis", False)),
    )
    return {
        "merit_scaled": float(merit),
        "gas_rms": float(group["gas_rms"]),
        "solid_rms": float(group["solid_rms"]),
        "energy_rms": float(group["energy_rms"]),
    }


def _bed0_holdup_detail(reactor: Reactor) -> dict[str, Any]:
    beds = _bed_cell_objs(reactor)
    if not beds:
        return {}
    c0 = beds[0]
    comp = ("char", "vm", "moist", "ash")
    idx = (S_CHAR, S_VM, S_MOISTURE, S_ASH)
    by_class = {n: float(c0.m_solid[:, i].sum()) for n, i in zip(comp, idx)}
    zu_by = {n: float(c0.m_solid_zu[:, i].sum()) for n, i in zip(comp, idx)}
    k_char = float(c0.K_solid_auf[0, S_CHAR] + c0.K_solid_ab[0, S_CHAR])
    return {
        "m_solid_total_kg": float(c0.m_solid.sum()),
        "m_solid_by_class_kg": by_class,
        "zu_by_class_kg_s": zu_by,
        "ab_in_kg_s": float(c0.m_solid_ab_in.sum()),
        "auf_in_kg_s": float(c0.m_solid_auf_in.sum()),
        "K_char_s_inv": k_char,
        "energy_scaled": float(c0.calc_energy_balance())
        / max(abs(float(reactor.config.fuel_feed)) * 20.0e6, 1.0e6),
        "eps_b": float(c0.eps_b),
    }


def _snapshot(reactor: Reactor, label: str, *, merit: bool = False) -> dict[str, Any]:
    out: dict[str, Any] = {"label": label, "bed0": _bed0_holdup_detail(reactor)}
    if merit:
        out["merit"] = _init_merit(reactor)
    return out


def _seed_bed_holdup_chain_skip_bed0(
    bed_cells: list,
    *,
    top_above_cell=None,
) -> bool:
    """与 ``seed_bed_holdup_chain_from_transport`` 相同，但 bed0 不 reseed。"""
    if not bed_cells:
        return False
    changed_any = False
    for _ in range(max(3 * len(bed_cells), 1)):
        update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
        changed = False
        for cell in bed_cells[1:]:
            changed |= seed_holdup_from_inflows(cell)
        changed_any |= changed
        if not changed:
            break
    update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
    return changed_any


def _seed_bed_holdup_chain_freeze_bed0(
    bed_cells: list,
    bed0_m_solid_frozen: np.ndarray,
    *,
    top_above_cell=None,
) -> bool:
    """全链 reseed，但每轮后将 bed0.m_solid 恢复为阶段 2 末的 transport x0。"""
    if not bed_cells:
        return False
    frozen = np.array(bed0_m_solid_frozen, dtype=np.float64, copy=True)
    changed_any = False
    for _ in range(max(3 * len(bed_cells), 1)):
        update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
        changed = False
        for cell in bed_cells:
            changed |= seed_holdup_from_inflows(cell)
        bed_cells[0].m_solid[:] = frozen
        changed_any |= changed
        if not changed:
            break
    update_bed_solid_transport_inflows_from_above(bed_cells, top_above_cell)
    bed_cells[0].m_solid[:] = frozen
    return changed_any


def _make_reseed_fn(
    variant: str,
    *,
    bed0_m_solid_frozen: np.ndarray | None = None,
) -> Callable[[Reactor], None]:
    beds_fn = _bed_cell_objs

    def _reseed_baseline(reactor: Reactor) -> None:
        reactor._seed_initialized_holdup_chain_for_nr()

    def _reseed_skip_bed0(reactor: Reactor) -> None:
        top = reactor.freeboard_cells[0] if reactor.freeboard_cells else None
        _seed_bed_holdup_chain_skip_bed0(beds_fn(reactor), top_above_cell=top)

    def _reseed_freeze_bed0(reactor: Reactor) -> None:
        assert bed0_m_solid_frozen is not None
        top = reactor.freeboard_cells[0] if reactor.freeboard_cells else None
        _seed_bed_holdup_chain_freeze_bed0(
            beds_fn(reactor),
            bed0_m_solid_frozen,
            top_above_cell=top,
        )

    def _reseed_noop(reactor: Reactor) -> None:
        return None

    if variant == "skip_bed0_all_reseeds":
        return _reseed_skip_bed0
    if variant == "freeze_bed0_transport_x0":
        return _reseed_freeze_bed0
    if variant == "no_holdup_reseed":
        return _reseed_noop
    return _reseed_baseline


def _run_stage3_with_variant(
    reactor: Reactor,
    state: init_step._InitPrecalcState,
    variant: str,
    *,
    bed0_m_solid_frozen: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    """重放阶段 3，按变体控制 holdup reseed。"""
    cfg = reactor.config
    cell_budgets = state.cell_budgets
    reseed_fn = _make_reseed_fn(variant, bed0_m_solid_frozen=bed0_m_solid_frozen)
    substeps: list[dict[str, Any]] = []

    reactor._apply_all_bc_for_nr()
    substeps.append(_snapshot(reactor, STAGE3_SUBSTEPS[0]))

    if variant != "no_holdup_reseed":
        reseed_fn(reactor)
    substeps.append(_snapshot(reactor, STAGE3_SUBSTEPS[1]))

    use_reactive = (
        state.use_a_tier_budget
        and cell_budgets is not None
        and bool(getattr(cfg, "vorab_a_tier_bed_solid_reactive_profile_thesis", True))
    )
    if use_reactive and int(getattr(cfg, "vorab_vm_devolatilization_zone_cells_thesis", 0)) > 0:
        align_vm_devolatilization_zone_vorab_state(
            reactor.cells,
            cell_budgets,
            cfg,
            apply_bc_fn=reactor._apply_all_bc_for_nr,
            reseed_holdup_chain_fn=(
                None
                if variant in {"seed_chain_no_reactive_reseed", "no_holdup_reseed"}
                else reseed_fn
            ),
            refresh_zone_pyrolysis_sources=False,
        )
    substeps.append(_snapshot(reactor, STAGE3_SUBSTEPS[2]))

    if state.use_a_tier_budget:
        reactor._vorab_a_tier_cell_budgets = cell_budgets
    return substeps


def _run_full_init_from_state(
    reactor: Reactor,
    state: init_step._InitPrecalcState,
    variant: str,
    *,
    bed0_m_solid_frozen: np.ndarray | None = None,
) -> dict[str, Any]:
    """从阶段 2 末状态出发：阶段 3（变体）→ 4 → 5，返回最终 merit。"""
    _run_stage3_with_variant(reactor, state, variant, bed0_m_solid_frozen=bed0_m_solid_frozen)
    init_step._run_bottom_zone_preprojections(reactor, state)
    after_stage4 = _snapshot(reactor, "4_bottom_zone_preprojections", merit=True)
    init_step._run_init_abgleich_and_joint_preprojection(reactor, state)
    after_stage5 = _snapshot(reactor, "5_abgleich_joint", merit=True)
    return {
        "variant": variant,
        "stage4": after_stage4,
        "stage5": after_stage5,
        "final_merit": after_stage5.get("merit"),
        "final_bed0_m_solid_kg": after_stage5["bed0"].get("m_solid_total_kg"),
    }


def _fresh_state_through_stage2(case_key: str) -> tuple[Reactor, init_step._InitPrecalcState]:
    reactor = build_phase2_reactor(
        validation_key=case_key,
        ladder_level="p2",
        apply_kinetics=True,
    )
    state = init_step._resolve_init_strategy_and_seed_T(reactor, init_strategy="vorabrechnung")
    state = init_step._run_vorabrechnung_snapshot(reactor, state)
    return reactor, state


def _clone_reactor_from_stage2(case_key: str) -> tuple[Reactor, init_step._InitPrecalcState]:
    """每次变体从独立 reactor 跑 stage2，保证初态一致。"""
    return _fresh_state_through_stage2(case_key)


VARIANTS = (
    "baseline",
    "freeze_bed0_transport_x0",
    "skip_bed0_all_reseeds",
    "seed_chain_no_reactive_reseed",
)


def run_audit(case_key: str) -> dict[str, Any]:
    reactor0, state0 = _fresh_state_through_stage2(case_key)
    bed0_frozen = np.array(_bed_cell_objs(reactor0)[0].m_solid, dtype=np.float64, copy=True)
    stage3_substeps = _run_stage3_with_variant(reactor0, state0, "baseline")

    variant_results: list[dict[str, Any]] = []
    for variant in VARIANTS:
        reactor, state = _clone_reactor_from_stage2(case_key)
        frozen = None
        if variant == "freeze_bed0_transport_x0":
            frozen = np.array(_bed_cell_objs(reactor)[0].m_solid, dtype=np.float64, copy=True)
        try:
            variant_results.append(
                _run_full_init_from_state(reactor, state, variant, bed0_m_solid_frozen=frozen)
            )
        except Exception as exc:
            variant_results.append({"variant": variant, "error": str(exc)})

    return {
        "case_key": case_key,
        "stage3_substeps_baseline": stage3_substeps,
        "variants": variant_results,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=CASE_DEFAULT)
    parser.add_argument(
        "--output",
        type=Path,
        default=data_path("init_holdup_seed_variants_w2.json"),
    )
    args = parser.parse_args()

    print(f"Holdup seed variants audit: {args.case} ...")
    payload = run_audit(args.case)
    write_json(args.output, payload)

    print("\n=== 阶段 3 子步 (baseline, bed0 m_solid kg) ===")
    for step in payload["stage3_substeps_baseline"]:
        b0 = step["bed0"]
        print(
            f"  {step['label']:36s}  total={b0['m_solid_total_kg']:7.2f}  "
            f"char={b0['m_solid_by_class_kg']['char']:6.1f}  ash={b0['m_solid_by_class_kg']['ash']:5.1f}  "
            f"E_scaled={b0['energy_scaled']:+.3f}"
        )

    print("\n=== 变体完整 init merit (stage 5) ===")
    for row in payload["variants"]:
        if "error" in row:
            print(f"  {row['variant']:32s}  ERROR: {row['error']}")
            continue
        m = row.get("final_merit") or {}
        s4 = row["stage4"].get("merit") or {}
        print(
            f"  {row['variant']:32s}  stage4 merit={s4.get('merit_scaled', float('nan')):7.4f}  "
            f"stage5 merit={m.get('merit_scaled', float('nan')):7.4f}  "
            f"bed0 m={row['final_bed0_m_solid_kg']:.2f} kg  "
            f"gas/solid/energy rms={m.get('gas_rms',0):.3f}/{m.get('solid_rms',0):.3f}/{m.get('energy_rms',0):.3f}"
        )

    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
