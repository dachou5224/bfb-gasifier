"""自由板显式求解图与床层顶格之间的桥接（closure → Cell 状态）。

将 ``simulate_freeboard`` 的输出写回 ``freeboard_cells``，供 NR 边界与 Vorabrechnung
刷新使用；与 Hamel 流程图中 PRECALC/连接区段一致。

Ref: ``docs/gasifier_model_flowchart_trilingual.mmd`` PRECALC / Zellenmodell 衔接；
     ``src/core/freeboard_segment.py``（解析轨迹与反应）
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from src.core.cell import S_ASH, S_CHAR
from src.core.cell_kinetics import build_reaction_sources
from src.core.freeboard_segment import simulate_freeboard
from src.kinetics.char_reactions import d_core_from_spm_char_conversion
from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL

if TYPE_CHECKING:
    from src.core.reactor import Reactor


def _spm_shrink_proxy_from_reaction_partition(
    total_char_conversion_proxy: float,
    combustion_share_proxy: float,
) -> tuple[float, float, float]:
    """Map bed heterogeneous conversion to a reduced size-age shrink proxy.

    Ref: Hamel (1999) §5.1.1.1 and Fig. 5.3: gasification advances particle
    "age" / carbon-density depletion without changing the outer diameter, while
    combustion follows a shrinking-particle path.  The same combustion mass
    sink therefore removes a larger shell volume after gasification has lowered
    the remaining carbon density.
    """
    x_total = float(np.clip(total_char_conversion_proxy, 0.0, 1.0 - 1e-9))
    f_comb = float(np.clip(combustion_share_proxy, 0.0, 1.0))
    x_comb = float(np.clip(x_total * f_comb, 0.0, 1.0 - 1e-9))
    x_gas = float(np.clip(x_total - x_comb, 0.0, 1.0 - 1e-9))
    x_spm = float(np.clip(x_comb / max(1.0 - x_gas, 1e-9), 0.0, 1.0 - 1e-9))
    return x_spm, x_comb, x_gas


def _max_entropy_age_weights(age_grid: np.ndarray, target_mean: float) -> np.ndarray:
    """Discrete maximum-entropy weights on a fixed Hamel age grid."""
    ages = np.asarray(age_grid, dtype=np.float64)
    if ages.size <= 1:
        return np.ones(max(int(ages.size), 1), dtype=np.float64)
    target = float(np.clip(target_mean, float(np.min(ages)), float(np.max(ages))))
    if target <= float(ages[0]) + 1e-12:
        w = np.zeros_like(ages)
        w[0] = 1.0
        return w
    if target >= float(ages[-1]) - 1e-12:
        w = np.zeros_like(ages)
        w[-1] = 1.0
        return w
    lo, hi = -120.0, 120.0
    for _ in range(96):
        lam = 0.5 * (lo + hi)
        a = lam * ages
        w = np.exp(a - float(np.max(a)))
        w /= max(float(np.sum(w)), 1e-300)
        mean = float(np.sum(w * ages))
        if mean < target:
            lo = lam
        else:
            hi = lam
    a = (0.5 * (lo + hi)) * ages
    w = np.exp(a - float(np.max(a)))
    return w / max(float(np.sum(w)), 1e-300)


def _expand_bed_top_age_size_launch_quadrature(
    *,
    bed_top_proxy: dict[str, Any],
    m_char_classes: np.ndarray,
    m_ash_classes: np.ndarray,
    n_age_bins: int,
    max_age: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float | int]]:
    """Expand bed-top launch classes into closure-only age/size quadrature.

    The expanded samples are only used by the freeboard trajectory. They carry
    an aggregate index so the closure can collapse back to the reactor's
    original particle-size classes before syncing explicit freeboard cells.
    """
    d_p_input = np.asarray(bed_top_proxy["d_p_input"], dtype=np.float64)
    d_p_eff = np.asarray(bed_top_proxy["d_p_eff"], dtype=np.float64)
    m_char = np.maximum(np.asarray(m_char_classes, dtype=np.float64), 0.0)
    m_ash = np.maximum(np.asarray(m_ash_classes, dtype=np.float64), 0.0)
    n_base = int(d_p_eff.size)
    if int(n_age_bins) <= 1 or n_base == 0:
        return (
            d_p_eff,
            m_char,
            m_ash,
            np.arange(n_base, dtype=np.int64),
            {
                "enabled": 0,
                "n_age_bins": 1,
                "age_mean": float(bed_top_proxy.get("bed_gasification_age_proxy", 0.0)),
                "age_max": 0.0,
                "tail_weight_ge_0p9": 0.0,
                "expanded_class_count": n_base,
                "d_p_min_m": float(np.min(d_p_eff)) if d_p_eff.size else 0.0,
                "d_p_max_m": float(np.max(d_p_eff)) if d_p_eff.size else 0.0,
            },
        )

    n_bins = int(max(n_age_bins, 1))
    age_hi = float(np.clip(max_age, 0.0, 1.0 - 1e-9))
    age_mean = float(np.clip(bed_top_proxy.get("bed_gasification_age_proxy", 0.0), 0.0, age_hi))
    x_comb = float(np.clip(bed_top_proxy.get("bed_combustion_conversion_proxy", 0.0), 0.0, 1.0 - 1e-9))
    x_local = float(np.clip(bed_top_proxy.get("top_local_char_conversion", 0.0), 0.0, 1.0 - 1e-9))
    age_grid = np.linspace(0.0, age_hi, n_bins, dtype=np.float64)
    weights = _max_entropy_age_weights(age_grid, age_mean)

    d_parts: list[np.ndarray] = []
    char_parts: list[np.ndarray] = []
    ash_parts: list[np.ndarray] = []
    idx_parts: list[np.ndarray] = []
    for age, weight in zip(age_grid, weights):
        x_spm_age = float(np.clip(x_comb / max(1.0 - float(age), 1e-9), 0.0, 1.0 - 1e-9))
        x_spm = float(np.clip(max(x_local, x_spm_age), 0.0, 1.0 - 1e-9))
        d_age = np.array(
            [d_core_from_spm_char_conversion(x_spm, float(dp)) for dp in d_p_input],
            dtype=np.float64,
        )
        d_parts.append(np.minimum(d_p_input, np.maximum(d_age, 1e-9)))
        char_parts.append(m_char * float(weight))
        ash_parts.append(m_ash * float(weight))
        idx_parts.append(np.arange(n_base, dtype=np.int64))

    d_expanded = np.concatenate(d_parts)
    m_char_expanded = np.concatenate(char_parts)
    m_ash_expanded = np.concatenate(ash_parts)
    aggregate_idx = np.concatenate(idx_parts)
    return (
        d_expanded,
        m_char_expanded,
        m_ash_expanded,
        aggregate_idx,
        {
            "enabled": 1,
            "n_age_bins": n_bins,
            "age_mean": age_mean,
            "age_max": age_hi,
            "tail_weight_ge_0p9": float(np.sum(weights[age_grid >= 0.9])),
            "expanded_class_count": int(d_expanded.size),
            "d_p_min_m": float(np.min(d_expanded)) if d_expanded.size else 0.0,
            "d_p_max_m": float(np.max(d_expanded)) if d_expanded.size else 0.0,
        },
    )


def _bed_reaction_bundle_for_entrained_size_proxy(cell):
    """Rebuild the current bed-cell heterogeneous bundle for entrained-size audit."""
    if getattr(cell, "cell_type", "") != "bed":
        return None
    if float(np.sum(np.maximum(cell.m_solid[:, S_CHAR], 0.0))) <= 0.0:
        return None
    thermo = cell._get_local_thermo_bundle()
    areas = cell._calc_char_surface_area_per_class()
    if float(np.sum(np.maximum(areas, 0.0))) <= 0.0:
        return None
    return build_reaction_sources(
        T=cell.T,
        P=cell.P,
        fuel_type=cell.fuel_type,
        V_b=cell.V_b,
        V_d=cell.V_d,
        C_b=thermo["C_b"],
        C_d=thermo["C_d"],
        y_b=thermo["y_b"],
        y_d=thermo["y_d"],
        gas_src_vm=cell._vm_gas_source_cache,
        solid_sink_vm=cell._vm_solid_sink_cache,
        areas=areas,
        solid_d_p=cell.solid.d_p,
        D_g=float(thermo["D_g"]),
        char_conversion=cell._compute_char_conversion(),
        rho_cat=cell._catalyst_bulk_density(),
        enable_r12=cell.enable_r12,
        use_gibbs_minor=cell.use_gibbs_minor,
        gibbs_minor_sources=None,
        r4_scale=cell.r4_scale,
        r5_scale=cell.r5_scale,
        r6_scale=cell.r6_scale,
        r7_scale=cell.r7_scale,
        rate_multiplier=1.0,
        N_zu_d=cell.N_zu_d,
        N_d_in=cell.N_d_in,
        N_zu_b=cell.N_zu_b,
        N_b_in=cell.N_b_in,
        N_rez_d=cell.N_rez_d,
        N_rez_b=cell.N_rez_b,
        N_ex=cell.N_ex,
        solid_shape=cell.R_solid.shape,
        char_index=S_CHAR,
    )


def bed_top_entrained_size_proxy(reactor: "Reactor") -> dict[str, Any]:
    """Estimate an effective bed-top entrained particle diameter.

    Hamel treats outer-diameter shrink mainly through combustion-driven size
    migration (SPM-like), while gasification contributes primarily to
    age/conversion progression. Under the current ``nk=1`` simplification we do
    not have the full 2D size-age matrix, so this helper builds a conservative
    proxy from the current bed heterogeneous reaction field:

    - total bed-side heterogeneous char consumption vs. top-bed upward char
      transport gives a steady entrained conversion proxy;
    - the combustion share of that heterogeneous sink maps the total proxy to
      an SPM-like outer-shrink proxy for the exported diameter.
    """
    top = reactor.cells[-1]
    d_p_raw = np.array(top.solid.d_p_classes, dtype=np.float64, copy=True)
    if not bool(getattr(reactor.config, "freeboard_use_shrunk_bed_top_d_p", True)):
        return {
            "d_p_input": d_p_raw,
            "d_p_eff": d_p_raw,
            "top_local_char_conversion": float(top._compute_char_conversion()),
            "bed_total_char_conversion_proxy": 0.0,
            "bed_combustion_share_proxy": 0.0,
            "bed_combustion_conversion_proxy": 0.0,
            "bed_gasification_age_proxy": 0.0,
            "bed_spm_shrink_proxy": 0.0,
            "bed_hetero_char_consumption_kg_s": 0.0,
            "bed_r1_char_consumption_kg_s": 0.0,
            "bed_top_up_char_kg_s": float(np.sum(np.maximum(top._solid_upflow_rates()[:, S_CHAR], 0.0))),
        }
    top_local_char_conversion = float(np.clip(top._compute_char_conversion(), 0.0, 1.0 - 1e-9))
    bed_hetero_char_consumption = 0.0
    bed_r1_char_consumption = 0.0
    for cell in reactor.cells:
        bundle = _bed_reaction_bundle_for_entrained_size_proxy(cell)
        if bundle is None:
            continue
        bed_r1_char_consumption += float(max(bundle.extent_r1, 0.0)) * ATOMIC_MASS_KG_PER_MOL["C"]
        bed_hetero_char_consumption += float(
            max(bundle.extent_r1 + bundle.extent_r2 + bundle.extent_r3 + bundle.extent_r4, 0.0)
        ) * ATOMIC_MASS_KG_PER_MOL["C"]
    bed_top_up_char = float(np.sum(np.maximum(top._solid_upflow_rates()[:, S_CHAR], 0.0)))
    bed_total_char_conversion_proxy = float(
        np.clip(
            bed_hetero_char_consumption / max(bed_hetero_char_consumption + bed_top_up_char, 1e-12),
            0.0,
            1.0 - 1e-9,
        )
    )
    bed_combustion_share_proxy = float(
        np.clip(
            bed_r1_char_consumption / max(bed_hetero_char_consumption, 1e-12),
            0.0,
            1.0,
        )
    ) if bed_hetero_char_consumption > 1e-12 else 0.0
    age_adjusted_spm_proxy, bed_combustion_conversion_proxy, bed_gasification_age_proxy = (
        _spm_shrink_proxy_from_reaction_partition(
            bed_total_char_conversion_proxy,
            bed_combustion_share_proxy,
        )
    )
    bed_spm_shrink_proxy = float(
        np.clip(
            max(top_local_char_conversion, age_adjusted_spm_proxy),
            0.0,
            1.0 - 1e-9,
        )
    )
    if bed_spm_shrink_proxy <= 0.0:
        return {
            "d_p_input": d_p_raw,
            "d_p_eff": d_p_raw,
            "top_local_char_conversion": top_local_char_conversion,
            "bed_total_char_conversion_proxy": bed_total_char_conversion_proxy,
            "bed_combustion_share_proxy": bed_combustion_share_proxy,
            "bed_combustion_conversion_proxy": bed_combustion_conversion_proxy,
            "bed_gasification_age_proxy": bed_gasification_age_proxy,
            "bed_spm_shrink_proxy": bed_spm_shrink_proxy,
            "bed_hetero_char_consumption_kg_s": bed_hetero_char_consumption,
            "bed_r1_char_consumption_kg_s": bed_r1_char_consumption,
            "bed_top_up_char_kg_s": bed_top_up_char,
        }
    d_p_eff = np.array(
        [d_core_from_spm_char_conversion(bed_spm_shrink_proxy, float(dp)) for dp in d_p_raw],
        dtype=np.float64,
    )
    d_p_eff = np.minimum(d_p_raw, np.maximum(d_p_eff, 1e-9))
    return {
        "d_p_input": d_p_raw,
        "d_p_eff": d_p_eff,
        "top_local_char_conversion": top_local_char_conversion,
        "bed_total_char_conversion_proxy": bed_total_char_conversion_proxy,
        "bed_combustion_share_proxy": bed_combustion_share_proxy,
        "bed_combustion_conversion_proxy": bed_combustion_conversion_proxy,
        "bed_gasification_age_proxy": bed_gasification_age_proxy,
        "bed_spm_shrink_proxy": bed_spm_shrink_proxy,
        "bed_hetero_char_consumption_kg_s": bed_hetero_char_consumption,
        "bed_r1_char_consumption_kg_s": bed_r1_char_consumption,
        "bed_top_up_char_kg_s": bed_top_up_char,
    }


def effective_bed_top_entrained_d_p_classes(reactor: "Reactor") -> tuple[np.ndarray, float]:
    """Return effective bed-top entrained diameters and SPM shrink proxy."""
    proxy = bed_top_entrained_size_proxy(reactor)
    return np.array(proxy["d_p_eff"], dtype=np.float64, copy=True), float(proxy["bed_spm_shrink_proxy"])


def sync_freeboard_cells_from_closure(
    reactor: "Reactor",
    fb: dict[str, Any],
    *,
    preserve_gas_state: bool = False,
) -> None:
    """将 freeboard 闭包 ``fb`` 中的分段状态同步到 ``reactor.freeboard_cells``。"""
    if not reactor._use_explicit_freeboard_cells():
        return
    states = list(fb.get("states", []))
    if not states:
        return
    bed_top = reactor.cells[-1]
    rho_s = float(max(bed_top.solid.rho_s, 1e-9))
    grouped_states: list[list[Any]]
    if len(reactor.freeboard_cells) == len(states):
        grouped_states = [[st] for st in states]
    else:
        grouped_states = [[] for _ in reactor.freeboard_cells]
        for st in states:
            seg_idx = int(getattr(st, "segment_index", -1))
            if 0 <= seg_idx < len(grouped_states):
                grouped_states[seg_idx].append(st)
        if any(len(group) == 0 for group in grouped_states):
            return

    for i, (cell, st_group) in enumerate(zip(reactor.freeboard_cells, grouped_states)):
        st = st_group[-1]
        hold_char = np.sum(
            [np.asarray(getattr(seg, "m_hold_char_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64) for seg in st_group],
            axis=0,
        )
        hold_ash = np.sum(
            [np.asarray(getattr(seg, "m_hold_ash_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64) for seg in st_group],
            axis=0,
        )
        up_char = np.asarray(getattr(st, "m_dot_auf_char_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64)
        up_ash = np.asarray(getattr(st, "m_dot_auf_ash_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64)
        down_char = np.sum(
            [np.asarray(getattr(seg, "m_dot_ab_char_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64) for seg in st_group],
            axis=0,
        )
        down_ash = np.sum(
            [np.asarray(getattr(seg, "m_dot_ab_ash_classes", np.zeros(cell.solid.n_size_classes)), dtype=np.float64) for seg in st_group],
            axis=0,
        )
        active_class_mask = (
            (np.maximum(hold_char, 0.0) + np.maximum(hold_ash, 0.0) > 1e-12)
            | (np.maximum(up_char, 0.0) + np.maximum(up_ash, 0.0) > 1e-12)
            | (np.maximum(down_char, 0.0) + np.maximum(down_ash, 0.0) > 1e-12)
        )
        hold_total = np.maximum(hold_char + hold_ash, 0.0)
        k_auf_char = np.divide(
            np.maximum(up_char, 0.0),
            np.maximum(hold_char, 0.0),
            out=np.zeros_like(hold_char),
            where=np.maximum(hold_char, 0.0) > 1e-12,
        )
        k_auf_ash = np.divide(
            np.maximum(up_ash, 0.0),
            np.maximum(hold_ash, 0.0),
            out=np.zeros_like(hold_ash),
            where=np.maximum(hold_ash, 0.0) > 1e-12,
        )
        k_ab_char = np.divide(
            np.maximum(down_char, 0.0),
            np.maximum(hold_char, 0.0),
            out=np.zeros_like(hold_char),
            where=np.maximum(hold_char, 0.0) > 1e-12,
        )
        k_ab_ash = np.divide(
            np.maximum(down_ash, 0.0),
            np.maximum(hold_ash, 0.0),
            out=np.zeros_like(hold_ash),
            where=np.maximum(hold_ash, 0.0) > 1e-12,
        )

        if not preserve_gas_state:
            cell.T = float(st.T)
            cell.N_b.fill(0.0)
            cell.N_d[:] = np.maximum(np.asarray(st.N, dtype=np.float64), 0.0)
        cell.solid_state_model = "freeboard_closure"
        cell.freeboard_explicit_char_hetero_enabled = bool(getattr(reactor.config, "freeboard_explicit_enable_char_hetero", True))
        cell.m_solid.fill(0.0)
        cell.K_solid_auf.fill(0.0)
        cell.K_solid_ab.fill(0.0)
        cell.m_solid_in.fill(0.0)
        cell.m_solid_auf_in.fill(0.0)
        cell.m_solid_ab_in.fill(0.0)
        if hold_char.shape[0] == cell.solid.n_size_classes:
            cell.m_solid[:, S_CHAR] = np.maximum(hold_char, 0.0)
            cell.m_solid[:, S_ASH] = np.maximum(hold_ash, 0.0)
            cell.K_solid_auf[:, S_CHAR] = np.maximum(k_auf_char, 0.0)
            cell.K_solid_auf[:, S_ASH] = np.maximum(k_auf_ash, 0.0)
            cell.K_solid_ab[:, S_CHAR] = np.maximum(k_ab_char, 0.0)
            cell.K_solid_ab[:, S_ASH] = np.maximum(k_ab_ash, 0.0)
            cell._freeboard_active_char_ash_mask = np.asarray(active_class_mask, dtype=np.float64)
        cell.freeboard_u_bed_top = float(max(bed_top.u_b, 1e-9))
        cell.freeboard_d_b_bed_top = float(max(bed_top.d_b, 1e-9))
        cell.freeboard_height_from_bed = float(max(cell.geo.h_center - float(reactor.config.H_bed), 0.0))
        solid_density = float(np.sum(hold_total) / max(np.pi * cell.geo.D_bed**2 * cell.geo.dh / 4.0, 1e-12))
        cell.freeboard_eps_d_voidage = float(np.clip(1.0 - solid_density / rho_s, 0.0, 1.0))
        cell.calc_hydrodynamics()
        cell.freeboard_explicit_char_sink_applied_kg[:] = np.asarray(
            getattr(st, "m_char_sink_applied_classes", np.zeros(cell.solid.n_size_classes)),
            dtype=np.float64,
        )


def sync_freeboard_hydrodynamic_bridge_from_bed_top(reactor: "Reactor") -> None:
    """Refresh the freeboard ghost-bubble bridge without overwriting NR state.

    Hamel's freeboard closure takes ``u_gb,0`` and ``d_b,WS`` from the bed
    surface.  In the explicit freeboard solver graph, freeboard gas/energy are
    Newton unknowns, so a full closure sync would erase accepted NR state.  This
    bridge-only refresh updates just the hydrodynamic inputs used by the next
    Vorabrechnung snapshot.
    """
    if not reactor._use_explicit_freeboard_solver_graph() or not reactor.freeboard_cells or not reactor.cells:
        return
    bed_top = reactor.cells[-1]
    rho_s = float(max(bed_top.solid.rho_s, 1e-12))
    for cell in reactor.freeboard_cells:
        cell.freeboard_u_bed_top = float(max(bed_top.u_b, 1e-9))
        cell.freeboard_d_b_bed_top = float(max(bed_top.d_b, 1e-9))
        cell.freeboard_height_from_bed = float(max(cell.geo.h_center - float(reactor.config.H_bed), 0.0))
        hold_total = np.maximum(cell.m_solid[:, S_CHAR] + cell.m_solid[:, S_ASH], 0.0)
        solid_density = float(np.sum(hold_total) / max(np.pi * cell.geo.D_bed**2 * cell.geo.dh / 4.0, 1e-12))
        cell.freeboard_eps_d_voidage = float(np.clip(1.0 - solid_density / rho_s, 0.0, 1.0))
        cell._hydro_cache_valid = False


def refresh_explicit_freeboard_transport_from_closure(
    reactor: "Reactor",
    *,
    preserve_gas_state: bool = False,
) -> dict[str, Any] | None:
    """基于床顶格状态调用解析 freeboard，并同步到显式 freeboard cell 链。"""
    if not reactor._use_explicit_freeboard_solver_graph():
        return None
    # 延迟导入，避免与 reactor 模块循环依赖
    from src.core.reactor import _resolve_axial_heat_loss_distribution

    cfg = reactor.config
    top = reactor.cells[-1]
    bed_top_proxy = bed_top_entrained_size_proxy(reactor)
    top_up_char = np.maximum(top._solid_upflow_rates()[:, S_CHAR], 0.0)
    top_up_ash = np.maximum(top._solid_upflow_rates()[:, S_ASH], 0.0)
    d_p_launch, m_char_launch, m_ash_launch, class_aggregate_idx, age_quad_diag = (
        _expand_bed_top_age_size_launch_quadrature(
            bed_top_proxy=bed_top_proxy,
            m_char_classes=top_up_char,
            m_ash_classes=top_up_ash,
            n_age_bins=int(max(getattr(cfg, "freeboard_age_quadrature_bins", 1), 1)),
            max_age=float(getattr(cfg, "freeboard_age_quadrature_max_age", 0.98)),
        )
    )
    d_p_eff_bed_top = np.array(bed_top_proxy["d_p_eff"], dtype=np.float64, copy=True)
    enabled_reactions = tuple(cfg.freeboard_enabled_reactions)
    if not bool(getattr(cfg, "freeboard_closure_enable_char_hetero", True)):
        enabled_reactions = tuple(name for name in enabled_reactions if name not in {"R1", "R2", "R3", "R4"})
    _, freeboard_loss = _resolve_axial_heat_loss_distribution(cfg)
    total_height = float(cfg.H_bed + max(cfg.H_freeboard, 0.0))
    fb = simulate_freeboard(
        N_in=np.maximum(top.N_d + top.N_b, 0.0),
        T_in=float(top.T),
        P=cfg.P,
        D_bed=cfg.D_bed,
        H_freeboard=float(cfg.H_freeboard),
        n_cells=int(cfg.n_freeboard_cells),
        u_b_bed_top=float(top.u_b),
        d_b_bed_top=float(top.d_b),
        eps_b_bed_top=float(top.eps_b),
        eps_d_void_bed_top=float(top.eps_d_voidage),
        rho_solid_bed_top=float(top.solid.rho_s),
        d_p_classes_bed_top=np.array(d_p_launch, dtype=np.float64, copy=True),
        m_char_classes_bed_top=np.array(m_char_launch, dtype=np.float64, copy=True),
        m_ash_classes_bed_top=np.array(m_ash_launch, dtype=np.float64, copy=True),
        char_conversion_bed_top=float(bed_top_proxy["bed_total_char_conversion_proxy"]),
        class_aggregate_indices_bed_top=np.array(class_aggregate_idx, dtype=np.int64, copy=True),
        n_output_size_classes=int(top.solid.n_size_classes),
        fuel_type=cfg.fuel_type,
        heat_loss_frac=float(freeboard_loss),
        trajectory_model=str(cfg.freeboard_trajectory_model),
        trajectory_coeff_model=str(cfg.freeboard_trajectory_coeff_model),
        u_gb_scale=float(max(cfg.freeboard_u_gb_scale, 1e-6)),
        beta_a_scale=float(max(cfg.freeboard_beta_a_scale, 0.0)),
        velocity_sigma=float(max(cfg.freeboard_velocity_sigma, 0.0)),
        velocity_bins=int(max(cfg.freeboard_velocity_bins, 1)),
        phi_s_bed_top=float(top.solid.phi_s),
        z_base_m=float(cfg.H_bed),
        total_height_m=total_height,
        secondary_injection_xi=cfg.freeboard_secondary_injection_xi,
        secondary_O2_mol_s=float(max(cfg.freeboard_secondary_O2_mol_s, 0.0)),
        secondary_N2_mol_s=float(max(cfg.freeboard_secondary_N2_mol_s, 0.0)),
        secondary_H2O_mol_s=float(max(cfg.freeboard_secondary_H2O_mol_s, 0.0)),
        secondary_T_K=float(max(cfg.freeboard_secondary_T_K, 1.0)),
        secondary_local_refine=int(max(cfg.freeboard_secondary_local_refine, 1)),
        secondary_injection_mode=str(cfg.freeboard_secondary_injection_mode),
        enabled_reactions=enabled_reactions,
        inventory_char_sink_enabled=bool(getattr(cfg, "freeboard_explicit_enable_char_hetero", True)),
    )
    fb["bed_top_d_p_input_m"] = np.array(bed_top_proxy["d_p_input"], dtype=np.float64, copy=True)
    fb["bed_top_d_p_eff_m"] = np.array(bed_top_proxy["d_p_eff"], dtype=np.float64, copy=True)
    fb["bed_top_char_conversion"] = float(bed_top_proxy["bed_spm_shrink_proxy"])
    fb["bed_top_char_conversion_used"] = float(bed_top_proxy["bed_total_char_conversion_proxy"])
    fb["bed_top_char_conversion_local"] = float(bed_top_proxy["top_local_char_conversion"])
    fb["bed_top_char_conversion_proxy_total"] = float(bed_top_proxy["bed_total_char_conversion_proxy"])
    fb["bed_top_combustion_share_proxy"] = float(bed_top_proxy["bed_combustion_share_proxy"])
    fb["bed_top_combustion_conversion_proxy"] = float(bed_top_proxy["bed_combustion_conversion_proxy"])
    fb["bed_top_gasification_age_proxy"] = float(bed_top_proxy["bed_gasification_age_proxy"])
    fb["bed_top_age_quadrature_diag"] = dict(age_quad_diag)
    fb["bed_top_r1_char_consumption_kg_s"] = float(bed_top_proxy["bed_r1_char_consumption_kg_s"])
    fb["bed_top_hetero_char_consumption_kg_s"] = float(bed_top_proxy["bed_hetero_char_consumption_kg_s"])
    fb["bed_top_up_char_kg_s"] = float(bed_top_proxy["bed_top_up_char_kg_s"])
    sync_freeboard_cells_from_closure(reactor, fb, preserve_gas_state=preserve_gas_state)
    reactor._last_explicit_freeboard_closure = fb
    return fb
