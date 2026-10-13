"""轴向离散全局状态（Phase0 演示，已弃用）。

Canonical 状态容器：``src.core.reactor.Reactor`` + ``src.core.cell.Cell``。
"""

from __future__ import annotations
import dataclasses
from dataclasses import dataclass, field
import numpy as np
import numpy.typing as npt
from src.core.species import N_GAS
from src.core.cell import N_SOLID_COMP

@dataclass
class GlobalState:
    """Phase0 演示用集中状态；NR 主路径使用 ``Reactor`` / ``Cell``。"""
    T: npt.NDArray[np.float64]        # (n_cells,) [K]
    P: npt.NDArray[np.float64]        # (n_cells,) [Pa]
    N_b: npt.NDArray[np.float64]      # (n_cells, N_GAS) [mol/s]
    N_d: npt.NDArray[np.float64]      # (n_cells, N_GAS) [mol/s]
    m_solid: npt.NDArray[np.float64]  # (n_cells, n_size_classes, N_SOLID_COMP) [kg/s or kg]
    _iteration: int = field(default=0, repr=False)

    @property
    def n_cells(self) -> int:
        return self.T.shape[0]

    @property
    def n_size_classes(self) -> int:
        return self.m_solid.shape[1]

    def copy(self) -> GlobalState:
        return GlobalState(
            T=self.T.copy(),
            P=self.P.copy(),
            N_b=self.N_b.copy(),
            N_d=self.N_d.copy(),
            m_solid=self.m_solid.copy(),
            _iteration=int(self._iteration),
        )

    def updated(self, **kwargs) -> GlobalState:
        """返回新状态对象，并将外迭代版本号 +1。"""
        kwargs["_iteration"] = int(self._iteration) + 1
        return dataclasses.replace(self, **kwargs)
