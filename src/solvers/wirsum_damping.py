"""Wirsum (1997/98) gedämpftes Newton — Hamel Eq. 2.9 / Wirsum Eq. 2.23 knobs.

Source: Wirsum Diss. §2.3.3 (VDI-FB 6/383); Hamel (2001) §2.3.
Accept trial only if residual Euclidean norm decreases; binary step halving.
See docs/wirsum_1997_nr_solver_notes.md — do NOT equate max(split/energy) merit to Wirsum.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class WirsumDampedNewtonParams:
    """Thesis NR line-search defaults aligned with Hamel Eq. 2.9."""

    lambda_init: float = 1.0  # [-] first trial step scale when prefer_full_step is off
    lambda_min: float = 1.0 / 4096.0  # [-] minimum accepted damping factor
    n_damp_halvings: int = 14  # [-] binary halving trials before failure
    prefer_full_step: bool = True  # [-] try undamped Newton step first


def wirsum_damped_newton_params_from_config(cfg: Any) -> WirsumDampedNewtonParams:
    """Read thesis NR damping knobs from ``ReactorConfig`` (or compatible mapping)."""
    halvings = int(max(getattr(cfg, "nr_damping_halvings_thesis", 14), 1))
    lambda_min = float(max(getattr(cfg, "nr_lambda_min_thesis", 1.0 / 4096.0), 1e-12))
    lambda_init = float(np.clip(getattr(cfg, "nr_lambda_init_thesis", 1.0), lambda_min, 1.0))
    prefer_full = bool(getattr(cfg, "nr_prefer_full_step_thesis", True))
    return WirsumDampedNewtonParams(
        lambda_init=lambda_init,
        lambda_min=lambda_min,
        n_damp_halvings=halvings,
        prefer_full_step=prefer_full,
    )
