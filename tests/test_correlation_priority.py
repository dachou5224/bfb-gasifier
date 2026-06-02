from __future__ import annotations

import warnings

import numpy as np
import pytest

from src.gasifier.models.correlations.base import Correlation, CorrelationMeta
from src.gasifier.models.correlations.drag_models import estimate_drag_coefficient
from src.gasifier.models.correlations.heat_transfer_coeff import compound_heat_loss_fraction
from src.gasifier.models.correlations.reaction_kinetics import effective_rate_multiplier


class _DummyCorrelation(Correlation):
    meta = CorrelationMeta(
        name="dummy",
        thesis_ref="Eq. 0.0, p.0",
        valid_range={"Re": (0.0, 10.0)},
    )

    def compute(self, **kwargs):
        return float(kwargs["Re"]) * 2.0


def test_correlation_range_check_warns_without_raising() -> None:
    corr = _DummyCorrelation()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = corr(Re=11.0)
    assert out == 22.0
    assert caught
    assert "out of range" in str(caught[-1].message)


def test_drag_correlation_known_reference_points() -> None:
    # Re<1: Cd = 24/Re
    assert estimate_drag_coefficient(0.5) == pytest.approx(48.0, rel=1e-10)
    # Re>1000: asymptote
    assert estimate_drag_coefficient(2000.0) == pytest.approx(0.44, rel=1e-10)


def test_heat_loss_compound_mapping_is_physical() -> None:
    assert compound_heat_loss_fraction(0.0, 1.0) == 0.0
    assert compound_heat_loss_fraction(0.2, 0.0) == 0.0
    val = compound_heat_loss_fraction(0.2, 0.5)
    assert 0.0 < val < 0.2


def test_reaction_multiplier_clamp_limits_extremes() -> None:
    assert effective_rate_multiplier(-1.0) == 0.0
    assert effective_rate_multiplier(2.0) == 2.0
    assert effective_rate_multiplier(1e9) == 1.0e6
