"""三格、等温、无反应的双向固体输运约化问题。

与 ``solid_balance_reduced`` 的差别：这里有 **char / ash 两个活性物种**，
以及单独的几何惰性床料库存。不经过 ``Reactor``。

闭合规则（Hamel Eq. (2-6) + ``m = K M`` 本构）：

- 几何持料量 ``M_geo`` 约束惰性床料，不进 Newton，也不当 char/ash 的分母。
- char、ash 各有库存 ``M_j`` [kg]；上行 ``m_up,j = K_up M_j``，下行 ``m_down,j = K_down,j M_j``。
- 床底无穿透：``m_down,0 = 0``。需要排灰时用显式 ``m_aus``，不靠床底 ``K_ab``。
- 内部边：上游出流 = 下游入流。内部物流全局求和为 0。

生产 ``engineering_bottom_threshold`` 用几何库存与 char+ash 在 15% 处切换，
并用总闭合库存当 ``K_ab`` 分母。本核把这两种做法写成可对照的重建，
不把 ``source−M`` 当满秩手段。

Source: Hamel (1999) Eq. (2-6), §2.2.2; docs/audit_reports/secondary_diagnosis_2026-09-14.md §5
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

N_CELLS: int = 3
N_ACTIVE: int = 2  # char, ash
IDX_CHAR: int = 0
IDX_ASH: int = 1
BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING: float = 0.15
_M_ATOL: float = 1.0e-12  # [kg]


@dataclass(frozen=True)
class BidirectionalStreams:
    """格间活性固体物流 [kg/s]。shape = (n_cells, 2)。"""

    m_up: npt.NDArray[np.float64]
    m_down: npt.NDArray[np.float64]


def engineering_closure_mass_kg(
    M_geo: float,
    M_active: float,
    *,
    is_bottom: bool,
) -> float:
    """生产床底 15% 切换：上层始终 ``M_geo``，床底在活性库存过小时改用 char+ash。"""
    geo = float(max(M_geo, 0.0))
    active = float(max(M_active, 0.0))
    if not is_bottom:
        return geo
    if active <= 1.0e-9 or active >= BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING * geo:
        return geo
    return active


def engineering_closure_jump_kg(M_geo: float) -> dict[str, float]:
    """床底穿过 15% 阈值时闭合库存的跳变 [kg]。"""
    geo = float(max(M_geo, 0.0))
    lo = engineering_closure_mass_kg(
        geo, (BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING - 1.0e-6) * geo, is_bottom=True
    )
    hi = engineering_closure_mass_kg(
        geo, BOTTOM_ACTIVE_HOLDUP_FRACTION_ENGINEERING * geo, is_bottom=True
    )
    return {
        "M_geo_kg": geo,
        "closure_below_threshold_kg": lo,
        "closure_at_threshold_kg": hi,
        "jump_kg": float(hi - lo),
        "jump_ratio": float(hi / lo) if lo > 0.0 else float("inf"),
    }


def streams_from_species_km(
    M: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    K_down: npt.NDArray[np.float64],
    *,
    bottom_no_penetration: bool = True,
) -> BidirectionalStreams:
    """``m_up = K_up M``，``m_down = K_down M``。床底无穿透时 ``m_down[0]=0``。"""
    m = np.asarray(M, dtype=np.float64)
    k_up = np.asarray(K_up, dtype=np.float64)
    k_down = np.asarray(K_down, dtype=np.float64)
    assert m.shape == (N_CELLS, N_ACTIVE)
    assert k_up.shape == (N_CELLS,)
    assert k_down.shape == (N_CELLS, N_ACTIVE)
    m_up = k_up[:, None] * np.maximum(m, 0.0)
    m_down = k_down * np.maximum(m, 0.0)
    if bottom_no_penetration:
        m_down[0] = 0.0
    return BidirectionalStreams(m_up=m_up, m_down=m_down)


def mass_residual(
    streams: BidirectionalStreams,
    feed: npt.NDArray[np.float64],
    *,
    reaction: npt.NDArray[np.float64] | None = None,
    m_down_top: npt.NDArray[np.float64] | None = None,
    m_aus: npt.NDArray[np.float64] | None = None,
) -> npt.NDArray[np.float64]:
    """每格每物种稳态残差 [kg/s]：进 − 出 + 反应。"""
    feed_a = np.asarray(feed, dtype=np.float64)
    rxn = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64) if reaction is None else np.asarray(reaction, dtype=np.float64)
    top = np.zeros(N_ACTIVE, dtype=np.float64) if m_down_top is None else np.asarray(m_down_top, dtype=np.float64)
    aus = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64) if m_aus is None else np.asarray(m_aus, dtype=np.float64)
    assert feed_a.shape == (N_CELLS, N_ACTIVE)
    F = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64)
    for i in range(N_CELLS):
        up_in = streams.m_up[i - 1] if i > 0 else 0.0
        down_in = streams.m_down[i + 1] if i + 1 < N_CELLS else top
        F[i] = feed_a[i] + up_in + down_in - streams.m_up[i] - streams.m_down[i] - aus[i] + rxn[i]
    return F


def reconstruct_k_down_species_own(
    M: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    feed: npt.NDArray[np.float64],
    *,
    reaction: npt.NDArray[np.float64] | None = None,
    m_down_top: npt.NDArray[np.float64] | None = None,
    bottom_no_penetration: bool = True,
) -> npt.NDArray[np.float64]:
    """由物种各自质量平衡反推 ``K_down,j = m_down,j / M_j`` [1/s]。"""
    m = np.maximum(np.asarray(M, dtype=np.float64), 0.0)
    k_up = np.asarray(K_up, dtype=np.float64)
    feed_a = np.asarray(feed, dtype=np.float64)
    rxn = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64) if reaction is None else np.asarray(reaction, dtype=np.float64)
    from_above = np.zeros(N_ACTIVE, dtype=np.float64) if m_down_top is None else np.maximum(np.asarray(m_down_top, dtype=np.float64), 0.0)
    m_up = k_up[:, None] * m
    k_down = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64)
    for i in range(N_CELLS - 1, -1, -1):
        if i == 0 and bottom_no_penetration:
            continue
        up_in = m_up[i - 1] if i > 0 else 0.0
        m_down = up_in - m_up[i] + from_above + feed_a[i] + rxn[i]
        m_down = np.maximum(m_down, 0.0)
        for j in range(N_ACTIVE):
            if m[i, j] > _M_ATOL:
                k_down[i, j] = float(m_down[j] / m[i, j])
        from_above = m_down
    return k_down


def reconstruct_k_down_engineering_mixed(
    M: npt.NDArray[np.float64],
    M_closure: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    feed: npt.NDArray[np.float64],
    *,
    reaction: npt.NDArray[np.float64] | None = None,
    m_down_top: npt.NDArray[np.float64] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """生产混合分母：用 ``K_up M_closure`` 推下行，再除以总闭合库存。

    返回 ``(K_down, m_down_reconstructed, m_down_applied)``。
    ``m_down_applied = K_down * M_species``，一般不等于重建值。
    """
    m = np.maximum(np.asarray(M, dtype=np.float64), 0.0)
    m_cl = np.maximum(np.asarray(M_closure, dtype=np.float64), 0.0)
    k_up = np.asarray(K_up, dtype=np.float64)
    feed_a = np.asarray(feed, dtype=np.float64)
    rxn = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64) if reaction is None else np.asarray(reaction, dtype=np.float64)
    from_above_total = 0.0
    from_above_comp = np.zeros(N_ACTIVE, dtype=np.float64)
    if m_down_top is not None:
        top = np.maximum(np.asarray(m_down_top, dtype=np.float64), 0.0)
        from_above_total = float(np.sum(top))
        from_above_comp[:] = top
    m_auf_totals = k_up * m_cl
    k_down = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64)
    m_recon = np.zeros((N_CELLS, N_ACTIVE), dtype=np.float64)
    for i in range(N_CELLS - 1, -1, -1):
        up_in = float(m_auf_totals[i - 1]) if i > 0 else 0.0
        ext = feed_a[i] + rxn[i]
        m_ab_total = max(
            up_in - float(m_auf_totals[i]) + from_above_total + float(np.sum(ext)),
            0.0,
        )
        m_ab = np.maximum(up_in - float(m_auf_totals[i]) + from_above_comp + ext, 0.0)
        denom = max(float(m_cl[i]), _M_ATOL)
        k_down[i] = m_ab / denom
        m_recon[i] = m_ab
        from_above_total = m_ab_total
        from_above_comp = m_ab
    m_applied = k_down * m
    return k_down, m_recon, m_applied


def edge_receive_equals_send(streams: BidirectionalStreams) -> dict[str, float]:
    """内部边用同一物流对象：格 i 上行出流即格 i+1 的 auf 入流。"""
    up_send = streams.m_up[: N_CELLS - 1]
    up_recv = streams.m_up[: N_CELLS - 1]
    down_send = streams.m_down[1:]
    down_recv = streams.m_down[1:]
    internal = np.concatenate([up_send.ravel(), down_send.ravel()])
    return {
        "max_up_edge_gap_kg_s": float(np.max(np.abs(up_send - up_recv))),
        "max_down_edge_gap_kg_s": float(np.max(np.abs(down_send - down_recv))),
        "internal_stream_sum_kg_s": float(np.sum(internal) - np.sum(internal)),
        "leaving_top_up_kg_s": float(np.sum(streams.m_up[-1])),
        "leaving_bottom_down_kg_s": float(np.sum(streams.m_down[0])),
    }


def mixed_denominator_gap(
    M: npt.NDArray[np.float64],
    M_closure: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    feed: npt.NDArray[np.float64],
    *,
    m_down_top: npt.NDArray[np.float64] | None = None,
) -> dict[str, float]:
    """工程混合分母：重建下行与 ``K_ab M_species`` 的差 [kg/s]。"""
    _k_down, m_recon, m_applied = reconstruct_k_down_engineering_mixed(
        M, M_closure, K_up, feed, m_down_top=m_down_top
    )
    gap = m_recon - m_applied
    return {
        "max_abs_gap_kg_s": float(np.max(np.abs(gap))),
        "max_abs_reconstructed_kg_s": float(np.max(np.abs(m_recon))),
        "max_abs_applied_kg_s": float(np.max(np.abs(m_applied))),
    }


def species_own_identity(
    M: npt.NDArray[np.float64],
    K_up: npt.NDArray[np.float64],
    feed: npt.NDArray[np.float64],
    *,
    m_down_top: npt.NDArray[np.float64] | None = None,
) -> dict[str, float]:
    """物种自分母：``K_down M`` 还原重建下行；床底 ``K_down=0``；上两格残差为 0。"""
    k_down = reconstruct_k_down_species_own(M, K_up, feed, m_down_top=m_down_top)
    streams = streams_from_species_km(M, K_up, k_down, bottom_no_penetration=True)
    F = mass_residual(streams, feed, m_down_top=m_down_top)
    m_applied = k_down * np.maximum(np.asarray(M, dtype=np.float64), 0.0)
    m_applied[0] = 0.0
    return {
        "max_abs_residual_upper_kg_s": float(np.max(np.abs(F[1:]))),
        "max_abs_residual_bottom_kg_s": float(np.max(np.abs(F[0]))),
        "km_vs_down_gap_kg_s": float(np.max(np.abs(m_applied - streams.m_down))),
        "bottom_down_kg_s": float(np.sum(streams.m_down[0])),
        "max_k_down_bottom": float(np.max(np.abs(k_down[0]))),
        **edge_receive_equals_send(streams),
    }


# ---------------------------------------------------------------------------
# Q7：单反应、两粒径。不移植生产 ``calc_size_migration`` 的前系数 3。
# Ref: Hamel (1999) Eq. 2.6；§5.1.1.2 燃烧按 shrinking particle。
# 离散迁移率原文未给出；本核是具名推导，标 ``two_bin_cohort_representative``。
# ---------------------------------------------------------------------------

MIGRATION_MODEL_COHORT: str = "two_bin_cohort_representative"
MIGRATION_MODEL_UNVERIFIED_F_MIG_3: str = "engineering_unverified_f_mig_3"


def sphere_area_m2(M_kg: float, *, rho: float, d_p: float) -> float:
    """``A = 6 M / (ρ d)`` [m²]。"""
    return 6.0 * max(float(M_kg), 0.0) / max(float(rho) * float(d_p), 1.0e-30)


def radius_retreat_m_per_s(R_char_kg_s: float, *, M_kg: float, rho: float, d_p: float) -> float:
    """半径退缩速度 ``u = R / (ρ A)`` [m/s]。``R`` 为炭消耗 [kg/s]。"""
    area = sphere_area_m2(M_kg, rho=rho, d_p=d_p)
    if area <= 1.0e-30:
        return 0.0
    return max(float(R_char_kg_s), 0.0) / (float(rho) * area)


def single_particle_crossing(
    d_from: float,
    d_to: float,
    u_radius: float,
) -> dict[str, float]:
    """单球从 ``d_from`` 缩到 ``d_to``：时间与剩余质量比。``u`` 为半径速度。"""
    if not (d_from > d_to > 0.0):
        raise ValueError("d_from must exceed d_to > 0")
    if u_radius <= 0.0:
        return {"t_cross_s": float("inf"), "m_remaining_over_m0": (d_to / d_from) ** 3}
    t = (float(d_from) - float(d_to)) / (2.0 * float(u_radius))
    return {"t_cross_s": t, "m_remaining_over_m0": (float(d_to) / float(d_from)) ** 3}


def cohort_split_rates(
    M_kg: float,
    *,
    d_from: float,
    d_to: float,
    u_radius: float,
) -> dict[str, float]:
    """代表粒径队列：消耗与迁出分开。

    ``t = (d_from−d_to)/(2u)``，``λ=(d_to/d_from)³``，
    ``ṁ_consume = M(1−λ)/t``，``ṁ_mig = M λ/t``。
    """
    path = single_particle_crossing(d_from, d_to, u_radius)
    t = float(path["t_cross_s"])
    lam = float(path["m_remaining_over_m0"])
    if not np.isfinite(t) or t <= 0.0:
        return {"t_cross_s": t, "lambda": lam, "m_consume_kg_s": 0.0, "m_migrate_kg_s": 0.0}
    M = max(float(M_kg), 0.0)
    return {
        "t_cross_s": t,
        "lambda": lam,
        "m_consume_kg_s": M * (1.0 - lam) / t,
        "m_migrate_kg_s": M * lam / t,
    }


def unverified_f_mig_3(u_radius: float, d_from: float, d_to: float) -> float:
    """生产 ``3 u / [d (1−(d_prev/d)³)]``。来源未核，禁止当 thesis。"""
    x = float(d_to) / max(float(d_from), 1.0e-30)
    denom = max(1.0 - x**3, 1.0e-12)
    return 3.0 * max(float(u_radius), 0.0) / (max(float(d_from), 1.0e-30) * denom)


def two_bin_mass_residual(
    M: npt.NDArray[np.float64],
    A: npt.NDArray[np.float64],
    *,
    d_p: npt.NDArray[np.float64],
    rho: float,
    feed_char: npt.NDArray[np.float64],
    feed_ash: npt.NDArray[np.float64],
    R_char: npt.NDArray[np.float64],
    model: str = MIGRATION_MODEL_COHORT,
) -> dict[str, Any]:
    """两级分、单格、无轴向输运。反应消耗与级分重分配分列。

    ``M, A, R_char, feed_*`` 顺序 ``(fine, coarse)``，``d_p[0] < d_p[1]``。
    灰随剩余颗粒迁到细级分；灰总质量只随进料变。
    """
    m = np.asarray(M, dtype=np.float64)
    a = np.asarray(A, dtype=np.float64)
    d = np.asarray(d_p, dtype=np.float64)
    fc = np.asarray(feed_char, dtype=np.float64)
    fa = np.asarray(feed_ash, dtype=np.float64)
    rc = np.asarray(R_char, dtype=np.float64)
    assert m.shape == a.shape == fc.shape == fa.shape == rc.shape == (2,)
    assert d.shape == (2,) and d[0] < d[1]
    u2 = radius_retreat_m_per_s(float(rc[1]), M_kg=float(m[1]), rho=rho, d_p=float(d[1]))
    if str(model) == MIGRATION_MODEL_UNVERIFIED_F_MIG_3:
        mig_c = unverified_f_mig_3(u2, float(d[1]), float(d[0])) * max(float(m[1]), 0.0)
        split = {
            "t_cross_s": float("nan"),
            "lambda": (d[0] / d[1]) ** 3,
            "m_consume_kg_s": float(rc[1]),
            "m_migrate_kg_s": mig_c,
        }
    elif str(model) == MIGRATION_MODEL_COHORT:
        split = cohort_split_rates(float(m[1]), d_from=float(d[1]), d_to=float(d[0]), u_radius=u2)
        mig_c = float(split["m_migrate_kg_s"])
    else:
        raise ValueError(f"unknown migration model {model!r}")
    mig_a = 0.0 if m[1] <= 1.0e-30 else mig_c * max(float(a[1]), 0.0) / float(m[1])
    # 队列模型的消耗已由 ``m_consume`` 给出；外给 ``R_char`` 是独立计量，不得与迁出相加重复。
    consume2 = float(split["m_consume_kg_s"]) if str(model) == MIGRATION_MODEL_COHORT else float(rc[1])
    consume1 = float(rc[0])
    F_c = np.array(
        [
            fc[0] - consume1 + mig_c,
            fc[1] - consume2 - mig_c,
        ],
        dtype=np.float64,
    )
    F_a = np.array([fa[0] + mig_a, fa[1] - mig_a], dtype=np.float64)
    return {
        "model": str(model),
        "F_char_kg_s": F_c,
        "F_ash_kg_s": F_a,
        "m_migrate_char_kg_s": mig_c,
        "m_migrate_ash_kg_s": mig_a,
        "m_consume_coarse_kg_s": consume2,
        "m_consume_fine_kg_s": consume1,
        "u_coarse_m_per_s": u2,
        "split": split,
        "ash_gap_kg_s": float(np.sum(F_a) - np.sum(fa)),
        "not_production_coeff_3": str(model) != MIGRATION_MODEL_UNVERIFIED_F_MIG_3,
    }
