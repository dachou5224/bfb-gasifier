"""P0 Wirsum 对齐：baseline JSON 契约与 reduced 残差可复现性。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from freeze_wirsum_alignment_baseline import (  # noqa: E402
    EXIT_SPECIES,
    SCHEMA_VERSION,
    capture_reduced_3cell_init,
    snapshot_config_knobs,
)
from tests.validation_case_utils import (  # noqa: E402
    build_phase1_htw_lu_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_PATH = REPO_ROOT / "data" / "wirsum_alignment_p0_baseline.json"

REQUIRED_PACKED_KEYS = {
    "residual_l2",
    "residual_linf",
    "rms_scaled",
    "jacobian_nnz",
    "jacobian_zero_rows",
    "jacobian_zero_cols",
    "exit_T_K",
    "exit_gas_dry",
}


def test_phase1_and_phase2_knob_snapshot_includes_handoff_keys() -> None:
    p1 = snapshot_config_knobs(build_phase1_htw_lu_reactor_config())
    p2 = snapshot_config_knobs(build_phase2_htw_lu_freeboard_reactor_config())
    for knobs in (p1, p2):
        assert "extent_limiters_enabled_thesis" in knobs
        assert "r2_scale" in knobs
        assert "r4_scale" in knobs
        assert "r5_scale" in knobs
        assert "r7_scale" in knobs
        assert "r8_scale" in knobs
        assert "gas_inlet_split_strategy" in knobs
        assert "bed_solid_state_model" in knobs
    assert "freeboard_trajectory_model" in p2
    assert "freeboard_trajectory_coeff_model" in p2


def test_reduced_3cell_init_metrics_are_finite() -> None:
    payload = capture_reduced_3cell_init(include_jacobian=False)
    assert payload["n_solver_cells"] == 3
    assert payload["n_var"] > 0
    for key in ("residual_l2", "residual_linf", "rms_scaled", "exit_T_K"):
        assert payload[key] is not None
        assert np.isfinite(float(payload[key]))
    for sp in EXIT_SPECIES:
        assert sp in payload["exit_gas_dry"]


@pytest.mark.skipif(not BASELINE_PATH.is_file(), reason="P0 baseline JSON not frozen yet")
def test_p0_baseline_json_schema() -> None:
    data = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    assert data["schema_version"] == SCHEMA_VERSION
    assert "sanity_checks" in data
    assert "failed_checks" in data["sanity_checks"]
    cases = data["cases"]
    assert "reduced_3cell_lu_init" in cases
    reduced = cases["reduced_3cell_lu_init"]
    assert REQUIRED_PACKED_KEYS.issubset(reduced)
    assert int(reduced["jacobian_nnz"]) > 0
    if "phase1_htw_lu_solve" in cases:
        lu = cases["phase1_htw_lu_solve"]
        assert lu["exit_T_K"] is not None
        assert "H2" in lu["exit_gas_dry"]
        assert "CO" in lu["exit_gas_dry"]
