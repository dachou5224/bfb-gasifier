"""Init bed 水力学：Vorab 气量 u0_target + resident K 闭合。"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from src.workflow.steps import init_precalc_step as init_step
def test_w2_stage3_holdup_and_eps_b_with_vorab_u0_and_resident_k():
    """W2 路径：阶段 3 后 bed0 holdup 不应爆炸，eps_b 应显著高于 N-only 微泡区。"""
    from _lu_harness import build_phase2_reactor

    reactor = build_phase2_reactor(
        validation_key="CASE_HTW_WESSELING_2",
        ladder_level="p2",
        apply_kinetics=True,
    )
    reactor.config.vorab_init_bed_vorab_gas_u0_target_thesis = True

    state = init_step._resolve_init_strategy_and_seed_T(reactor, init_strategy="vorabrechnung")
    state = init_step._run_vorabrechnung_snapshot(reactor, state)
    bed0 = next(c for c in reactor.cells if str(getattr(c, "cell_type", "bed")) == "bed")
    assert float(bed0.eps_b) < 0.15

    saved = init_step._apply_init_bed_vorab_gas_u0_targets(
        reactor,
        T_est=state.T_est,
        cell_budgets=state.cell_budgets,
    )
    assert float(bed0.eps_b) > 0.20
    try:
        init_step._run_post_vorabrechnung_holdup_and_reactive(reactor, state)
    finally:
        init_step._restore_macro_hydrodynamics_seed(reactor, saved)

    m_hold = float(bed0.m_solid.sum())
    assert m_hold < 15.0
    assert float(bed0.eps_b) > 0.15
