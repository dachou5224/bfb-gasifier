#!/usr/bin/env python3
"""Audit key Hamel source-of-truth consistency anchors.

This audit is intentionally narrow: it checks whether the current codebase still
matches the locked Hamel-aligned implementation choices documented in
``docs/hamel_submodels/00-07``. It is *not* a fit-quality or validation-case audit.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.constants import P0_HAMEL
from src.core.reactor import Reactor, ReactorConfig
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.thermodynamics.equilibrium import get_K_eq
from tests.validation_case_utils import (
    PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase2_htw_lu_freeboard_reactor_config,
    load_case_LU,
)


def _check(label: str, ok: bool, detail: str) -> bool:
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {label}: {detail}")
    return ok


def _contains(path: Path, needle: str) -> bool:
    return needle in path.read_text(encoding="utf-8")


def main() -> int:
    ok = True

    print("Hamel consistency audit")
    print("=" * 80)

    critical_ref_files = [
        REPO_ROOT / "src" / "physics" / "mass_transfer.py",
        REPO_ROOT / "src" / "physics" / "bubble_dynamics.py",
        REPO_ROOT / "src" / "kinetics" / "arrhenius.py",
        REPO_ROOT / "src" / "kinetics" / "gas_reactions.py",
        REPO_ROOT / "src" / "core" / "feed_inlet.py",
        REPO_ROOT / "tests" / "validation_case_utils.py",
    ]
    for path in critical_ref_files:
        has_ref = _contains(path, "docs/hamel_submodels")
        ok &= _check(path.relative_to(REPO_ROOT).as_posix(), has_ref, "references docs/hamel_submodels")

    u_d = 0.2688
    p = 2.5e6
    u_br_model = calc_u_br(u_d, p)
    u_br_ref = 2.7 * u_d * (p / P0_HAMEL) ** (-0.15)
    ok &= _check(
        "Eq. 3.44 u_br",
        math.isclose(u_br_model, u_br_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={u_br_model:.12g}, ref={u_br_ref:.12g}",
    )

    d_b = 0.2121
    D_g = 4.0e-5
    eps_mf = 0.45
    u_b = 0.55
    kbd_model = calc_kbd(u_br_model, d_b, D_g, eps_mf, u_b)
    kbd_ref = 3.0 * u_br_ref / (2.0 * d_b) + math.sqrt((144.0 * D_g * eps_mf * u_b) / (math.pi * d_b**3))
    ok &= _check(
        "Eq. 3.50 K_bd",
        math.isclose(kbd_model, kbd_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={kbd_model:.12g}, ref={kbd_ref:.12g}",
    )

    T = 1200.0
    keq_model = get_K_eq("R8", T)
    keq_ref = math.exp(-3.6893 + 4019.0 / T)
    ok &= _check(
        "Eq. 5.47 K_eq,R8",
        math.isclose(keq_model, keq_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={keq_model:.12g}, ref={keq_ref:.12g}",
    )

    phase1_cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_LU())
    ok &= _check(
        "Phase1 hydrodynamics chain",
        (
            phase1_cfg.hydrodynamics_u_d_closure == "wein_1992_eq312"
            and phase1_cfg.hydrodynamics_bubble_diameter_model == "hilligardt_ode"
            and phase1_cfg.hydrodynamics_psi_b_strategy == "wein_1992"
            and phase1_cfg.hydrodynamics_lambda_strategy == "hamel_280"
            and phase1_cfg.hydrodynamics_xi_strategy == "hamel_regime"
            and phase1_cfg.hydrodynamics_bubble_velocity_strategy == "heinbockel_eq343"
            and phase1_cfg.hydrodynamics_bubble_ode_strategy == "heinbockel_eq341"
        ),
        "thesis-aligned hydrodynamics strategies are active",
    )

    base_cfg = ReactorConfig()
    thesis_cfg = ReactorConfig(thesis_mode=True)
    ok &= _check(
        "Base vs thesis freeboard R8 policy",
        ("R8" not in base_cfg.freeboard_enabled_reactions) and ("R8" in thesis_cfg.freeboard_enabled_reactions),
        f"base_has_R8={'R8' in base_cfg.freeboard_enabled_reactions}, thesis_has_R8={'R8' in thesis_cfg.freeboard_enabled_reactions}",
    )
    ok &= _check(
        "Thesis major Gibbs mode",
        thesis_cfg.major_gibbs_solver_mode == "hamel_reduced",
        f"mode={thesis_cfg.major_gibbs_solver_mode}",
    )

    phase2_cfg = build_phase2_htw_lu_freeboard_reactor_config(load_case_LU())
    phase2_reactor = Reactor(phase2_cfg)
    ok &= _check(
        "Phase2 Jacobian delegation",
        PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS.get("nr_jacobian_strategy") is None,
        "solve kwargs delegate topology-based default via nr_jacobian_strategy=None",
    )
    ok &= _check(
        "Phase2 resolved Jacobian default",
        phase2_reactor._default_nr_jacobian_strategy() in {"block_tridiag_structured", "band_plus_side_elements_structured"},
        f"default={phase2_reactor._default_nr_jacobian_strategy()}",
    )

    print("=" * 80)
    if ok:
        print("HAMEL CONSISTENCY AUDIT PASSED")
        return 0
    print("HAMEL CONSISTENCY AUDIT FAILED")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
