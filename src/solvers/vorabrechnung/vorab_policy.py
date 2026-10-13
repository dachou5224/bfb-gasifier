"""Vorabrechnung 外层源项策略。

``engineering_single_shot``：init 预算一次干燥/热解并冻结；外层只刷水力。
这是当前 ``thesis_mode`` 默认，**不是** Hamel Bild 2.2 的外层完整重算
``q = P(x)``。整个求解期间固定源项不符合 Hamel §2.1 / Wirsum §2.2。

``hamel_outer_refresh``：外层按当前解重算干燥/热解（与水力、K）。
内层仍只读冻结上下文。opt-in，默认不改生产。

``auto``：``thesis_mode and thesis_vorab_sources_single_shot`` → engineering，
否则 → hamel。``wirsum_1997_strict`` 目前不改本开关。

Source: Hamel (1999) §2.1 Bild 2.2; Wirsum (1997) §2.2;
docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 阶段 D
"""

from __future__ import annotations

from typing import Any

from src.core.cell import Cell

VORAB_OUTER_SOURCE_AUTO: str = "auto"
VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT: str = "engineering_single_shot"
VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH: str = "hamel_outer_refresh"
VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF: str = "engineering_H_over_umf"

_VORAB_OUTER_SOURCE_NAMES = frozenset(
    {
        VORAB_OUTER_SOURCE_AUTO,
        VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT,
        VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH,
    }
)


def resolve_vorab_outer_source_policy_name(configured: str) -> str:
    """校验配置名，不把 ``auto`` 展开。"""
    name = str(configured).strip().lower()
    if name not in _VORAB_OUTER_SOURCE_NAMES:
        raise ValueError(
            f"Unknown vorab_outer_source_policy={configured!r}; "
            f"expected {sorted(_VORAB_OUTER_SOURCE_NAMES)}"
        )
    return name


def thesis_vorab_sources_single_shot_enabled(cfg: Any) -> bool:
    """现行 thesis 开关：``thesis_mode`` 且 ``thesis_vorab_sources_single_shot``。

    这是工程 single-shot，不是已核对的 Hamel 外层。显式
    ``vorab_outer_source_policy`` 覆盖 ``vorab_outer_may_refresh_*``。
    """
    return bool(getattr(cfg, "thesis_mode", False)) and bool(
        getattr(cfg, "thesis_vorab_sources_single_shot", True)
    )


def resolve_vorab_outer_source_policy(cfg: Any) -> str:
    """展开 ``auto`` 后的外层源项策略。"""
    raw = resolve_vorab_outer_source_policy_name(
        getattr(cfg, "vorab_outer_source_policy", VORAB_OUTER_SOURCE_AUTO)
    )
    if raw != VORAB_OUTER_SOURCE_AUTO:
        return raw
    if thesis_vorab_sources_single_shot_enabled(cfg):
        return VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT
    return VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH


def vorab_outer_may_refresh_drying_pyro_sources(cfg: Any) -> bool:
    """outer Abgleich / bridge 是否允许重算干燥/DAEM 源项。"""
    return resolve_vorab_outer_source_policy(cfg) == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH


def vorab_refresh_policy_label(cfg: Any) -> str:
    """结果汇总用策略标签。

    ``auto`` 保持历史字符串，避免现有 NR 摘要断言断裂。
    显式政策用新名；语义以 ``resolve_vorab_outer_source_policy`` 为准。
    """
    configured = resolve_vorab_outer_source_policy_name(
        getattr(cfg, "vorab_outer_source_policy", VORAB_OUTER_SOURCE_AUTO)
    )
    if configured == VORAB_OUTER_SOURCE_HAMEL_OUTER_REFRESH:
        return "hamel_outer_refresh_full_P_of_x"
    if configured == VORAB_OUTER_SOURCE_ENGINEERING_SINGLE_SHOT:
        return "engineering_single_shot_sources_outer_refresh_hydrodynamics_fixed_inner"
    if thesis_vorab_sources_single_shot_enabled(cfg):
        return "single_shot_sources_outer_refresh_hydrodynamics_fixed_inner"
    return "outer_refresh_fixed_inner_sources"


def vorab_tau_model_label() -> str:
    """现行 ``vorabrechnung_tau_for_cell`` 的模型名。不是 Hamel 颗粒停留。"""
    return VORAB_TAU_MODEL_ENGINEERING_H_OVER_UMF


def cell_vorab_drying_pyro_sources_frozen(cell: Cell) -> bool:
    return bool(getattr(cell, "_vorab_drying_pyro_sources_frozen", False))


def freeze_vorab_drying_pyro_sources_for_cells(cells: list[Cell]) -> None:
    """Init single-shot 完成后锁定干燥/热解源项缓存。"""
    for cell in cells:
        if not bool(getattr(cell, "_vm_cache_valid", False)):
            continue
        cell._vorab_drying_pyro_sources_frozen = True


def clear_vorab_drying_pyro_source_freeze_for_cells(cells: list[Cell]) -> None:
    """新一轮 global NR 开始前解除冻结（init 将重建源项）。"""
    for cell in cells:
        cell._vorab_drying_pyro_sources_frozen = False
