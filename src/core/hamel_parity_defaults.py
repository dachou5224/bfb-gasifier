"""Hamel Eq.(2-6) / 自由板轨迹默认配置（thesis / 非 thesis 共用）。

Source: docs/hamel_submodels/02_cell_balances_and_exchange.md Eq. (2-6)
       docs/hamel_submodels/04_freeboard_particle_trajectory_analytic.md
"""

from __future__ import annotations

from typing import Protocol

# Eq. (2-6) holdup_transport：显式 K_auf/K_ab 与 m 出流项（Hamel 主链）
HAMEL_BED_SOLID_STATE_MODEL: str = "holdup_transport"

# 稳态入流压缩分支（工程/审计对照，须显式 opt-in）
LEGACY_BED_SOLID_STATE_MODEL: str = "legacy_stream"

# Hamel z=73.69·exp(-5.0784·φ_s) 阻力系数（exact_hamel）
HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL: str = "exact_hamel"

# Haider 混合阻力 + 工程稳定分支（审计对照）
LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL: str = "stable_mixed_drag_split"

# Eq. (2.7) 显式 Q_W / Q_WÜ（Hamel 主链，2026-07-07）
HAMEL_ENERGY_HEAT_LOSS_MODE: str = "hamel_eq27_explicit"

# f·max(H_in,0) 工程分支（审计 / LADDER baseline 对照）
LEGACY_ENERGY_HEAT_LOSS_MODE: str = "inlet_fraction"

# Eq. (2.3)–(2.6) 粒径迁移：按 size class 本地库存 m_jk（Hamel 主链）
HAMEL_SIZE_MIGRATION_INVENTORY_MODE: str = "per_class_holdup"

# V_d·ρ_s·(1−ε_mf)·(m/m_out_total) 全局代理（审计对照）
LEGACY_SIZE_MIGRATION_INVENTORY_MODE: str = "total_bed_proxy"

# Anhang A1 (A.35)–(A.40) 初值 + (A.34) Holub 收敛判据（Hamel 主链，2026-07-08）
HAMEL_GIBBS_INITIAL_GUESS_MODE: str = "hamel_a40"
HAMEL_GIBBS_TOL_MODE: str = "holub_vonka_a34"

# Taylor 线性化初值 + L2 残差判据（审计对照）
LEGACY_GIBBS_INITIAL_GUESS_MODE: str = "taylor_linearized"
LEGACY_GIBBS_TOL_MODE: str = "l2_residual"

# 新鲜进料粒径档：粗档单峰（Hamel Startwertwahl / Table A.2.1.1 语义）
HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE: str = "coarsest_class"

# 均匀 1/nk（历史工程简化，审计对照）
LEGACY_FEED_SOLID_SIZE_DISTRIBUTION_MODE: str = "uniform"

# Init/outer solid-holdup 预投影：Eq.(2-6) 闭合时不联动主反应动力学（Hamel Vorab 轻预计算）
PREPROJ_SOLID_HOLDUP_REACTION_RATE_MULTIPLIER: float = 0.0

# 气相/能量预投影 LSQ：同理冻结主反应（Vorab 已给源项预算）
PREPROJ_GAS_REACTION_RATE_MULTIPLIER: float = 0.0

# Major phase-split 预投影须与 inner NR 同 chemistry；rm=0 会闭合无反应 split，inner 一步 reopen
PREPROJ_PHASE_SPLIT_REACTION_RATE_MULTIPLIER: float = 1.0


class _SolidCellLike(Protocol):
    solid_state_model: str


def apply_hamel_bed_solid_state_to_cell(cell: _SolidCellLike) -> None:
    """将 bed cell 固相状态模型设为 Hamel holdup_transport 主链。"""
    cell.solid_state_model = HAMEL_BED_SOLID_STATE_MODEL


def resolve_size_migration_inventory_mode(*, configured: str) -> str:
    """Resolve size-migration inventory basis."""
    name = str(configured).strip()
    if name in {HAMEL_SIZE_MIGRATION_INVENTORY_MODE, LEGACY_SIZE_MIGRATION_INVENTORY_MODE}:
        return name
    raise ValueError(
        f"Unsupported size_migration_inventory_mode={configured!r}; "
        f"expected {HAMEL_SIZE_MIGRATION_INVENTORY_MODE!r} or "
        f"{LEGACY_SIZE_MIGRATION_INVENTORY_MODE!r}"
    )


def resolve_feed_solid_mass_fractions(*, n_size_classes: int, mode: str) -> "np.ndarray":
    """进料固相按粒径档质量分数（Hamel：未反应 fresh fuel 集中在粗档）。

    Ref: Hamel (1999) Kap. 2.1 Startwertwahl; Table A.2.1.1 离散粒径类。
    """
    import numpy as np

    nk = max(int(n_size_classes), 1)
    name = str(mode).strip()
    if name == HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE:
        frac = np.zeros(nk, dtype=np.float64)
        frac[-1] = 1.0
        return frac
    if name == LEGACY_FEED_SOLID_SIZE_DISTRIBUTION_MODE:
        return np.ones(nk, dtype=np.float64) / float(nk)
    raise ValueError(
        f"Unsupported feed_solid_size_distribution_mode={mode!r}; "
        f"expected {HAMEL_FEED_SOLID_SIZE_DISTRIBUTION_MODE!r} or "
        f"{LEGACY_FEED_SOLID_SIZE_DISTRIBUTION_MODE!r}"
    )


def resolve_bed_solid_state_model(*, configured: str, thesis_mode: bool) -> str:
    """Resolve bed solid model; thesis_mode 强制 holdup_transport。"""
    if bool(thesis_mode):
        return HAMEL_BED_SOLID_STATE_MODEL
    name = str(configured).strip()
    if name in {HAMEL_BED_SOLID_STATE_MODEL, LEGACY_BED_SOLID_STATE_MODEL}:
        return name
    raise ValueError(
        f"Unsupported bed_solid_state_model={configured!r}; "
        f"expected {HAMEL_BED_SOLID_STATE_MODEL!r} or {LEGACY_BED_SOLID_STATE_MODEL!r}"
    )
