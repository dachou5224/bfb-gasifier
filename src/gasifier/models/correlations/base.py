from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import warnings

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class CorrelationMeta:
    """关联式元信息（名称、论文出处、适用范围）。"""

    name: str
    thesis_ref: str
    valid_range: dict[str, tuple[float, float]]


class Correlation(ABC):
    """经验关联式抽象基类。"""

    meta: CorrelationMeta

    def __call__(self, **kwargs: float | npt.NDArray[np.float64]) -> float | npt.NDArray[np.float64]:
        self._check_range(kwargs)
        return self.compute(**kwargs)

    @abstractmethod
    def compute(self, **kwargs: float | npt.NDArray[np.float64]) -> float | npt.NDArray[np.float64]:
        """Compute correlation output."""
        raise NotImplementedError

    def _check_range(self, kwargs: dict[str, float | npt.NDArray[np.float64]]) -> None:
        """Warn (not raise) when inputs are outside declared validity range."""
        for key, bounds in self.meta.valid_range.items():
            if key not in kwargs:
                continue
            lo, hi = bounds
            val = np.asarray(kwargs[key], dtype=np.float64)
            if np.any(val < lo) or np.any(val > hi):
                warnings.warn(
                    f"{self.meta.name}: {key}={val} out of range [{lo}, {hi}] "
                    f"(ref: {self.meta.thesis_ref})",
                    stacklevel=2,
                )
