"""Vorabrechnung single-shot 与 outer-refresh 隔离策略。

Hamel thesis strict 口径：干燥/DAEM 源项在 init 末段预算一次并冻结；
outer Abgleich 仅刷新水力学快照，内层 NR 固定源项。

Ref: Hamel (1999) Bild 2.2; docs/hamel_isomorphic_gap_implementation_checklist.md §3
"""

from __future__ import annotations

from typing import Any

from src.core.cell import Cell


def thesis_vorab_sources_single_shot_enabled(cfg: Any) -> bool:
    """是否启用 thesis 严格 single-shot 干燥/热解源项冻结。"""
    return bool(getattr(cfg, "thesis_mode", False)) and bool(
        getattr(cfg, "thesis_vorab_sources_single_shot", True)
    )


def vorab_outer_may_refresh_drying_pyro_sources(cfg: Any) -> bool:
    """outer Abgleich / bridge 是否允许重算干燥/DAEM 源项。"""
    return not thesis_vorab_sources_single_shot_enabled(cfg)


def vorab_refresh_policy_label(cfg: Any) -> str:
    """结果汇总用策略标签。"""
    if thesis_vorab_sources_single_shot_enabled(cfg):
        return "single_shot_sources_outer_refresh_hydrodynamics_fixed_inner"
    return "outer_refresh_fixed_inner_sources"


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
