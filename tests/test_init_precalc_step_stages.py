"""OPT-013：init_precalc 五阶段编排契约测试。"""

from __future__ import annotations

import numpy as np
import pytest

from src.workflow.steps import init_precalc_step as step


def test_init_precalc_stage_functions_exist():
    for name in (
        "_resolve_init_strategy_and_seed_T",
        "_run_vorabrechnung_snapshot",
        "_run_post_vorabrechnung_holdup_and_reactive",
        "_run_bottom_zone_preprojections",
        "_run_init_abgleich_and_joint_preprojection",
    ):
        assert callable(getattr(step, name))


def test_run_init_and_precalc_invokes_five_stages_in_order(monkeypatch):
    order: list[str] = []

    def _state(**kwargs):
        return step._InitPrecalcState(
            resolved="vorabrechnung",
            T_est=np.array([900.0]),
            cell_budgets=None,
            use_a_tier_budget=False,
            use_cell_mapping=False,
            use_bed_stream_mapping=False,
            use_march=False,
            t_floor_K=1050.0,
            t_cap_K=1325.0,
            init_started=0.0,
            vorab_started=0.0,
            **kwargs,
        )

    def stage1(reactor, *, init_strategy):
        order.append("stage1")
        return _state()

    def stage2(reactor, state):
        order.append("stage2")
        return state

    def stage3(reactor, state):
        order.append("stage3")

    def stage4(reactor, state):
        order.append("stage4")

    def stage5(reactor, state):
        order.append("stage5")

    monkeypatch.setattr(step, "_resolve_init_strategy_and_seed_T", stage1)
    monkeypatch.setattr(step, "_run_vorabrechnung_snapshot", stage2)
    monkeypatch.setattr(step, "_run_post_vorabrechnung_holdup_and_reactive", stage3)
    monkeypatch.setattr(step, "_run_bottom_zone_preprojections", stage4)
    monkeypatch.setattr(step, "_run_init_abgleich_and_joint_preprojection", stage5)

    class _Reactor:
        cells = []
        config = type("Cfg", (), {"vorab_init_bed_vorab_gas_u0_target_thesis": False})()

    result = step.run_init_and_precalc_for_global_nr(_Reactor(), init_strategy=None, gs_warmup_steps=None)
    assert order == ["stage1", "stage2", "stage3", "stage4", "stage5"]
    assert result.resolved_init_strategy == "vorabrechnung"
    assert result.T_est.shape == (1,)


def test_bed0_fence_lower_margin_only_widens_bottom_cell():
    """仅底格加宽 Tmin；上段仍用全局 150 K margin。"""
    from src.core.reactor import Reactor
    from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed_temperature_march_thesis = False
    cfg.nr_temperature_fence_enabled_thesis = True
    cfg.nr_temperature_fence_lower_margin_K_thesis = 150.0
    cfg.nr_temperature_fence_bed0_lower_margin_K_thesis = 350.0
    reactor = Reactor(cfg)
    t_ref = 1150.0  # [K]
    t_est = np.full(len(reactor.cells), t_ref, dtype=np.float64)
    step._install_temperature_fence_from_vorabrechnung(reactor, t_est)
    assert reactor.cells[0]._nr_temperature_min_K == pytest.approx(t_ref - 350.0)
    if len(reactor.cells) > 1:
        assert reactor.cells[1]._nr_temperature_min_K == pytest.approx(t_ref - 150.0)


def test_upper_bed_fence_lower_margin_does_not_change_bed0():
    """上段加宽 Tmin 不改底格；bed1 用独立 margin（handoff §5.64）。"""
    from src.core.reactor import Reactor
    from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed_temperature_march_thesis = False
    cfg.nr_temperature_fence_enabled_thesis = True
    cfg.nr_temperature_fence_lower_margin_K_thesis = 150.0
    cfg.nr_temperature_fence_bed0_lower_margin_K_thesis = 350.0
    cfg.nr_temperature_fence_upper_bed_lower_margin_K_thesis = 350.0
    reactor = Reactor(cfg)
    t_ref = 1150.0  # [K]
    t_est = np.full(len(reactor.cells), t_ref, dtype=np.float64)
    step._install_temperature_fence_from_vorabrechnung(reactor, t_est)
    assert reactor.cells[0]._nr_temperature_min_K == pytest.approx(t_ref - 350.0)
    if len(reactor.cells) > 1:
        assert reactor.cells[1]._nr_temperature_min_K == pytest.approx(t_ref - 350.0)


def test_bed0_nr_startwert_T_clips_to_fence_and_overrides_vorab_profile():
    """末次 Startwert T 贴围栏；上段 T 不动。"""
    from src.core.reactor import Reactor
    from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config

    cfg = build_phase2_htw_lu_freeboard_reactor_config(enable_vorab_march_staged=True)
    cfg.vorab_bed_temperature_march_thesis = False
    cfg.nr_temperature_fence_enabled_thesis = True
    cfg.nr_temperature_fence_lower_margin_K_thesis = 150.0
    cfg.nr_temperature_fence_bed0_lower_margin_K_thesis = 350.0
    cfg.vorab_bed0_nr_startwert_T_K_thesis = 880.0
    reactor = Reactor(cfg)
    t_ref = 1150.0  # [K]
    t_est = np.full(len(reactor.cells), t_ref, dtype=np.float64)
    step._install_temperature_fence_from_vorabrechnung(reactor, t_est)
    reactor.cells[0].T = t_ref
    if len(reactor.cells) > 1:
        t1 = float(reactor.cells[1].T)
    step._apply_bed0_nr_startwert_T(reactor)
    assert reactor.cells[0].T == pytest.approx(880.0)
    if len(reactor.cells) > 1:
        assert reactor.cells[1].T == pytest.approx(t1)
    cfg.vorab_bed0_nr_startwert_T_K_thesis = 500.0
    step._apply_bed0_nr_startwert_T(reactor)
    assert reactor.cells[0].T == pytest.approx(t_ref - 350.0)

