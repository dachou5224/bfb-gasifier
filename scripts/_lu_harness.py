"""Phase2 HTW LU freeboard 审计/基准脚本公共 harness。

集中 ``sys.path``、Reactor 构建、init/solve 与 JSON 输出，避免各脚本重复样板。
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from time import perf_counter
from typing import Any

# OpenBLAS/OMP 多线程 + hash 随机化会使同构 LU 栈出现 CO≈0.15 / 0.21 双吸引子
# （handoff §5.29）。审计 harness **强制** 1 线程；PYTHONHASHSEED 须在启动前设为 0
# （进程内设置无效）。若需多线程：BFB_ALLOW_THREADED_BLAS=1。
if os.environ.get("BFB_ALLOW_THREADED_BLAS", "").strip() not in {"1", "true", "TRUE", "yes"}:
    for _thr_key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[_thr_key] = "1"
if sys.flags.hash_randomization and os.environ.get("PYTHONHASHSEED") not in {"0", "1"}:
    print(
        "[bfb-harness] 警告: 未固定 PYTHONHASHSEED；LU 同构栈可能落到 CO≈0.21 坏吸引子。"
        " 请用 PYTHONHASHSEED=0 python3 ...",
        file=sys.stderr,
    )

REF_CO = 0.13  # [-] 文献 Table 2 LU 出口干基参考
REF_CO2 = 0.11  # [-]

_REPO_ROOT: Path | None = None

# 与 ``benchmark_phase2_convergence_ladder.py`` 一致（SSOT）
LADDER_LEVELS: dict[str, dict[str, object]] = {
    "baseline": {
        "nr_line_search_gas_phase_split_merit_thesis": False,
        "nr_line_search_energy_merit_thesis": False,
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": False,
        "nr_enforce_gas_phase_split_convergence_thesis": False,
        "vorab_bottom_temperature_cap_K_thesis": 1325.0,
        "freeboard_trajectory_coeff_model": "stable_mixed_drag_split",
        "energy_balance_heat_loss_mode": "inlet_fraction",
    },
    "p0": {
        # Wirsum Eq. 2.23：默认 ‖F‖₂/RMS 接受；split/energy max-merit 为历史工程桥（可 opt-in）
        "nr_line_search_gas_phase_split_merit_thesis": False,
        "nr_line_search_energy_merit_thesis": False,
        # Historical ladder step that first enabled accept-step hydro refresh (engineering).
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1325.0,
        "freeboard_trajectory_coeff_model": "stable_mixed_drag_split",
        "energy_balance_heat_loss_mode": "inlet_fraction",
    },
    "p1": {
        "nr_line_search_gas_phase_split_merit_thesis": False,
        "nr_line_search_energy_merit_thesis": False,
        # Keep accept-refresh for historical ladder A/B vs baseline/p0.
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": True,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "stable_mixed_drag_split",
        "energy_balance_heat_loss_mode": "inlet_fraction",
    },
    "p2": {
        "nr_line_search_gas_phase_split_merit_thesis": False,
        "nr_line_search_energy_merit_thesis": False,
        # P1 Hamel outer_fixed: Inner solves F(x;h) with h fixed during Inner.
        "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
        "nr_bottom_bc_aware_bed0_bed1_joint_x0_thesis": True,
        "nr_enforce_gas_phase_split_convergence_thesis": True,
        "vorab_bottom_temperature_cap_K_thesis": 1180.0,
        "freeboard_trajectory_coeff_model": "exact_hamel",
    },
}

LADDER_P2: dict[str, object] = dict(LADDER_LEVELS["p2"])

# Hamel 方程同构 + Wirsum Eq.2.23，不含平台 CO 狩猎桥。
# 保留：去床层 extent、Vorab 快氧化 x₀（Startwert）、outer-fixed hydro、
# Phase2 builder 的 bed1 preinner / exact_hamel。
# 关闭：合计 LS fastox、syngas guard、local2/absorb、rate-gate、hot restart、床顶 solid LSQ。
# 底格：关 R1 置零；开 Eq.2.4/2.5 两相 Startwert + R1 残氧种子；LS 按相 fastox + 分相切向 clip。
# 能量：每步 |ΔT| 帽 40 K，避免底格一刀砸进围栏（§5.53）。生产默认仍 600 K。
# 底格围栏：CORE 下沿 T_ref−350≈800 K（§5.58）；生产默认仍全局 150 K。
# 上段围栏：CORE 仍 150 K。o8 把 bed7–9 钉在 ~1000 K；放到 350 只换来 rms 小数点后第四位（§5.64）。
# 底格 Startwert T：CORE=820 K（§5.60）。850 会砸进 800 K 栏；820 升回 ~840。围栏中心仍是 Vorab T_est。
# 固相：fastox/R1 种子后再对齐 Eq.2.6，并对 bed1 局部 snap（§5.54–5.55）。
# VM 分区：钉死 zone=3（Phase2 / Hamel Kap.2.1 轴向释放）。zone=1 能填 bed0 焓洞，
# 但 o8 把 H₂/CO Ê 打到 3×、rms 变差（§5.61）。勿改 Eq.2.7 / 勿开 VM passthrough。
# H2O 透传：CORE 关。fastox 后再贴是 opt-in 正确性（§5.62）；打开会把 T0 砸进 800 K 栏。
# 产物透传：CORE 关。fastox 后再贴 CO/CO2/H2/CH4/TAR（§5.63）；勿与 H2O 叠开验收。
HAMEL_WIRSUM_CORE: dict[str, object] = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "vorab_transport_x0_char_oxidation_o2_closure_thesis": False,
    "vorab_bed0_eq24_two_phase_startwert_thesis": True,
    "vorab_init_solve_bottom_cell_thesis": False,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_gas_phase_split_merit_thesis": False,
    "nr_line_search_energy_merit_thesis": False,
    "nr_line_search_fast_oxidation_projection_thesis": False,
    "nr_line_search_per_phase_fastox_thesis": True,
    "nr_clip_same_phase_oxidizer_fuel_step_thesis": True,
    "nr_clip_allow_oxidizer_into_fuel_thesis": False,
    "vorab_bed0_r1_oxidizer_seed_mol_s_thesis": 2.0,
    "nr_inner_t_step_default_cap_K_thesis": 40.0,
    "nr_temperature_fence_bed0_lower_margin_K_thesis": 350.0,
    "nr_temperature_fence_upper_bed_lower_margin_K_thesis": None,
    "vorab_bed0_nr_startwert_T_K_thesis": 820.0,
    "vorab_vm_devolatilization_zone_cells_thesis": 3,
    "nr_line_search_syngas_collapse_guard_thesis": False,
    "nr_outer_preinner_upper_co_local2x2_thesis": False,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": False,
    "nr_trace_o2_rate_gate_thesis": False,
    "nr_outer_hot_restart_after_y_co_thesis": False,
    "nr_outer_preinner_bed_top_solid_holdup_thesis": False,
    "vorab_a_tier_post_staged_h2o_passthrough_thesis": False,
    "vorab_a_tier_post_staged_syngas_passthrough_thesis": False,
}


def repo_root() -> Path:
    global _REPO_ROOT
    if _REPO_ROOT is None:
        _REPO_ROOT = Path(__file__).resolve().parent.parent
    return _REPO_ROOT


def ensure_repo_path() -> Path:
    root = repo_root()
    scripts = root / "scripts"
    for p in (root, scripts):
        ps = str(p)
        if ps not in sys.path:
            sys.path.insert(0, ps)
    return root


# 作为 ``from _lu_harness import ...`` 的副作用：``python scripts/foo.py`` 时自动可 import src/tests
ensure_repo_path()


def data_path(relative: str) -> Path:
    return repo_root() / "data" / relative


def apply_baseline_kinetics(cfg: Any) -> None:
    """Force all ``r*_scale = 1.0``.

    Phase2 LU 默认已是 ``r5_scale=1.0``。历史 ``0.25`` 压低曾用于抑制
    bed-top O₂ breakthrough / freeboard overburn（见 ``data/phase2_default_o2_fb_o8.json``）。
    """
    for s in ("r2", "r4", "r5", "r6", "r7", "r8"):
        setattr(cfg, f"{s}_scale", 1.0)


def apply_config_overrides(cfg: Any, overrides: dict[str, object]) -> None:
    for key, value in overrides.items():
        setattr(cfg, key, value)


def co_co2_score(
    co: float,
    co2: float,
    *,
    ref_co: float = REF_CO,
    ref_co2: float = REF_CO2,
) -> float:
    return abs(co - ref_co) / ref_co + abs(co2 - ref_co2) / ref_co2


def bed_cell_indices(reactor: Any) -> list[int]:
    return [
        i
        for i, c in enumerate(reactor.cells)
        if str(getattr(c, "cell_type", "bed")) == "bed"
    ]


def build_phase2_config(**builder_kwargs: Any) -> Any:
    ensure_repo_path()
    from tests.validation_case_utils import build_phase2_htw_lu_freeboard_reactor_config

    return build_phase2_htw_lu_freeboard_reactor_config(**builder_kwargs)


def build_phase2_reactor(
    *,
    case: dict[str, Any] | None = None,
    validation_key: str | None = None,
    builder_kwargs: dict[str, Any] | None = None,
    cfg_overrides: dict[str, object] | None = None,
    ladder_level: str | None = None,
    apply_kinetics: bool = False,
) -> Any:
    """构建 Phase2 freeboard ``Reactor``（可选 ladder 级别、case 与 cfg 覆盖）。"""
    ensure_repo_path()
    from src.core.reactor import Reactor

    merged_builder = dict(builder_kwargs or {})
    if validation_key is not None:
        from tests.validation.case_lu_loader import load_validation_case_flat

        case = load_validation_case_flat(validation_key)
    if case is not None:
        merged_builder["case"] = case
    cfg = build_phase2_config(**merged_builder)
    if ladder_level is not None:
        if ladder_level not in LADDER_LEVELS:
            raise ValueError(f"Unknown ladder level {ladder_level!r}")
        apply_config_overrides(cfg, LADDER_LEVELS[ladder_level])
    if cfg_overrides:
        apply_config_overrides(cfg, cfg_overrides)
    # cfg_overrides 可能翻转 enable_r12；种子策略须在覆盖之后重对齐
    from tests.validation.config_phase2_freeboard import apply_enable_r12_temperature_seed_policy

    apply_enable_r12_temperature_seed_policy(cfg)
    if apply_kinetics:
        apply_baseline_kinetics(cfg)
    return Reactor(cfg)


def phase2_solve_kwargs(**overrides: Any) -> dict[str, Any]:
    ensure_repo_path()
    from tests.validation_case_utils import PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS

    kw = dict(PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS)
    kw.update(overrides)
    return kw


def run_phase2_init(
    reactor: Any,
    *,
    init_strategy: str = "vorabrechnung",
    gs_warmup_steps: int | None = None,
) -> Any:
    ensure_repo_path()
    from src.workflow.steps.init_precalc_step import run_init_and_precalc_for_global_nr

    return run_init_and_precalc_for_global_nr(
        reactor,
        init_strategy=init_strategy,
        gs_warmup_steps=gs_warmup_steps,
    )


def timed_phase2_solve(reactor: Any, **solve_kwargs: Any) -> tuple[dict[str, Any], float]:
    t0 = perf_counter()
    out = reactor.solve(**solve_kwargs)
    return out, perf_counter() - t0


def write_json(path: Path | str, payload: Any, *, ensure_ascii: bool = False) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=ensure_ascii) + "\n",
        encoding="utf-8",
    )
    return out_path


def add_max_global_iter_arg(parser: argparse.ArgumentParser, *, default: int = 25) -> None:
    parser.add_argument("--max-global-iter", type=int, default=default)


def add_json_output_arg(parser: argparse.ArgumentParser, default_relative: str) -> None:
    parser.add_argument("--output", type=Path, default=data_path(default_relative))


def phase1_global_nr_solve_kwargs(**overrides: Any) -> dict[str, Any]:
    ensure_repo_path()
    from tests.validation_case_utils import PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS

    kw = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
    kw.update(overrides)
    return kw


def load_case_lu() -> dict[str, Any]:
    ensure_repo_path()
    from tests.validation_case_utils import load_case_LU

    return load_case_LU()
