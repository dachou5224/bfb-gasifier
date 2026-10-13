"""gasifier 域内模型层。

**生产路径（canonical）**：``src.core`` + ``src.solvers`` + ``src.workflow``。

**Phase0 演示（已弃用）**：下列 ``*Model`` / ``BalanceBuilder`` 在 ``GlobalState`` 上
薄封装 core 子模块，NR 主路径不引用；新功能不得在此层实现。

活跃子包：``models.correlations``（见 ``tests/test_correlation_priority.py``）。
"""

from __future__ import annotations

import importlib
import warnings

__all__ = [
    "BalanceBuilder",
    "ConnectivityModel",
    "HeatTransferModel",
    "HydrodynamicsModel",
    "PyrolysisModel",
    "ReactionModel",
]

# Phase0 GlobalState 演示类 → canonical 实现位置（OPT-011）
_PHASE0_EXPORTS: dict[str, tuple[str, str]] = {
    "BalanceBuilder": (
        "src.gasifier.models.balances",
        "src.core.cell_balances / Cell.residuals",
    ),
    "ConnectivityModel": (
        "src.gasifier.models.connectivity",
        "src.core.connectivity",
    ),
    "HeatTransferModel": (
        "src.gasifier.models.heat_transfer",
        "src.core.reactor.ReactorConfig / calc_energy_balance_residual",
    ),
    "HydrodynamicsModel": (
        "src.gasifier.models.hydrodynamics",
        "src.core.cell_hydrodynamics.calc_cell_hydrodynamics",
    ),
    "PyrolysisModel": (
        "src.gasifier.models.pyrolysis",
        "src.core.cell_pyrolysis.calc_drying_pyrolysis_sources",
    ),
    "ReactionModel": (
        "src.gasifier.models.reaction_model",
        "src.core.cell_kinetics.build_reaction_sources",
    ),
}


def __getattr__(name: str):
    if name not in _PHASE0_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    mod_path, canonical = _PHASE0_EXPORTS[name]
    warnings.warn(
        f"{__name__}.{name} is Phase0 GlobalState demo only; "
        f"implement and validate against {canonical}.",
        DeprecationWarning,
        stacklevel=2,
    )
    mod = importlib.import_module(mod_path)
    return getattr(mod, name)


def __dir__() -> list[str]:
    return sorted(__all__)
