#!/usr/bin/env python3
"""Audit key Hamel source-of-truth consistency anchors.

This audit is intentionally narrow: it checks whether the current codebase still
matches the locked Hamel-aligned implementation choices documented in
``docs/hamel_submodels/00-07``. It is *not* a fit-quality or validation-case audit.
"""

from __future__ import annotations

from _lu_harness import (
    build_phase2_config,
    load_case_lu,
    phase2_solve_kwargs,
    repo_root,
)

REPO_ROOT = repo_root()

import math
import sys
from pathlib import Path

import numpy as np

from src.core.constants import P0_HAMEL
from src.core.hamel_parity_defaults import (
    HAMEL_BED_SOLID_STATE_MODEL,
    HAMEL_ENERGY_HEAT_LOSS_MODE,
    HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL,
    HAMEL_GIBBS_INITIAL_GUESS_MODE,
    HAMEL_GIBBS_TOL_MODE,
    HAMEL_SIZE_MIGRATION_INVENTORY_MODE,
    LEGACY_BED_SOLID_STATE_MODEL,
    LEGACY_ENERGY_HEAT_LOSS_MODE,
    LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL,
    LEGACY_GIBBS_INITIAL_GUESS_MODE,
    LEGACY_GIBBS_TOL_MODE,
    LEGACY_SIZE_MIGRATION_INVENTORY_MODE,
)
from src.core.reactor import Reactor, ReactorConfig
from src.core.species import GAS_SPECIES_INDEX
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.thermodynamics.equilibrium import get_K_eq
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R12
from src.kinetics.hamel_r5_reference import (
    HAMEL_R5_BUBBLE_E_OVER_RG_K,
    HAMEL_R5_BUBBLE_K0,
    hamel_r5_bubble_rate_mol_m3_s,
    hamel_r5_suspension_rate_mol_m3_s,
)
from src.kinetics.hamel_r6_h2_reference import (
    HAMEL_R6_H2_A,
    HAMEL_R6_H2_E_OVER_RG_K,
    hamel_r6_h2_rate_mol_m3_s,
)
from src.kinetics.hamel_r7_ch4_reference import (
    HAMEL_R7_CH4_A,
    HAMEL_R7_CH4_E_OVER_RG_K,
    hamel_r7_ch4_rate_mol_m3_s,
)
from src.kinetics import gas_reactions as gr
from tests.validation_case_utils import build_phase1_htw_lu_global_nr_reactor_config

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
        REPO_ROOT / "src" / "physics" / "minimum_fluidization.py",
        REPO_ROOT / "src" / "physics" / "phase_fractions.py",
        REPO_ROOT / "src" / "physics" / "freeboard.py",
        REPO_ROOT / "src" / "kinetics" / "arrhenius.py",
        REPO_ROOT / "src" / "kinetics" / "gas_reactions.py",
        REPO_ROOT / "src" / "kinetics" / "tar_reactions.py",
        REPO_ROOT / "src" / "core" / "feed_inlet.py",
        REPO_ROOT / "src" / "core" / "cell.py",
        REPO_ROOT / "src" / "core" / "species.py",
        REPO_ROOT / "tests" / "validation_case_utils.py",
    ]
    stale_spec_needles = ("specs/00_", "specs/01_", "specs/02_", "specs/03_", "specs/04_", "specs/species.md")
    for path in critical_ref_files:
        has_ref = _contains(path, "docs/hamel_submodels")
        ok &= _check(path.relative_to(REPO_ROOT).as_posix(), has_ref, "references docs/hamel_submodels")
        text = path.read_text(encoding="utf-8")
        stale = [n for n in stale_spec_needles if n in text]
        ok &= _check(
            f"{path.relative_to(REPO_ROOT).as_posix()} stale specs",
            not stale,
            "none" if not stale else ", ".join(stale),
        )

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

    ok &= _check(
        "Table 5.3 R5 bubble k0",
        math.isclose(gr.R5_bubble_k0, HAMEL_R5_BUBBLE_K0, rel_tol=0.0, abs_tol=0.0),
        f"model={gr.R5_bubble_k0:.6g}, ref={HAMEL_R5_BUBBLE_K0:.6g}",
    )
    ok &= _check(
        "Table 5.3 R5 bubble E/Rg",
        math.isclose(gr.R5_bubble_E_Rg, HAMEL_R5_BUBBLE_E_OVER_RG_K, rel_tol=0.0, abs_tol=0.0),
        f"model={gr.R5_bubble_E_Rg:.6g}, ref={HAMEL_R5_BUBBLE_E_OVER_RG_K:.6g}",
    )

    T_r5 = 1030.0
    C_co, C_o2, C_h2o = 10.0, 0.5, 15.0
    y = np.zeros(11)
    r5b = rate_R5_bubble(T_r5, C_co, C_o2, C_h2o, 2.5e6, y, GAS_SPECIES_INDEX)
    r5b_ref = hamel_r5_bubble_rate_mol_m3_s(T_r5, C_co, C_o2, C_h2o)
    ok &= _check(
        "Eq. 5.36 R5 bubble rate",
        math.isclose(r5b, r5b_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={r5b:.6g}, ref={r5b_ref:.6g}",
    )

    r5d = rate_R5_suspension(T_r5, C_co, C_o2, C_h2o, 2.5e6, y, GAS_SPECIES_INDEX)
    r5d_hamel = hamel_r5_suspension_rate_mol_m3_s(T_r5, C_co, C_o2, C_h2o, temperature_exponent=-0.5)
    r5d_ratio = r5d / r5d_hamel if r5d_hamel > 0.0 else float("nan")
    ok &= _check(
        "Table 5.3 R5 suspension T exponent",
        math.isclose(r5d_ratio, 1.0, rel_tol=0.02),
        f"model/ref={r5d_ratio:.4g} (Hamel Table 5.3: 1.8·T^-0.5)",
    )

    T_hom = 1325.0
    C_h2, C_o2_h2 = 0.3, 0.01
    y_h = np.zeros(11)
    r12 = rate_R12(T_hom, C_h2, C_o2_h2, 2.5e6, y_h, GAS_SPECIES_INDEX)
    r12_ref = hamel_r6_h2_rate_mol_m3_s(T_hom, C_h2, C_o2_h2)
    ok &= _check(
        "Eq. 5.38 R12 (thesis R6) H2 rate",
        math.isclose(r12, r12_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={r12:.6g}, ref={r12_ref:.6g}",
    )
    ok &= _check(
        "Eq. 5.39 R12 constants",
        math.isclose(gr.R12_A, HAMEL_R6_H2_A, rel_tol=0.0, abs_tol=0.0)
        and math.isclose(gr.R12_E_Rg, HAMEL_R6_H2_E_OVER_RG_K, rel_tol=0.0, abs_tol=0.0),
        f"A={gr.R12_A:.6g}, E/Rg={gr.R12_E_Rg:.6g}",
    )

    C_ch4, C_o2_ch4 = 1.0, 0.02
    r6 = rate_R6(T_hom, C_ch4, C_o2_ch4)
    r6_ref = hamel_r7_ch4_rate_mol_m3_s(T_hom, C_ch4, C_o2_ch4)
    ok &= _check(
        "Eq. 5.40 R6 (thesis R7) CH4 rate",
        math.isclose(r6, r6_ref, rel_tol=1e-12, abs_tol=0.0),
        f"model={r6:.6g}, ref={r6_ref:.6g}",
    )
    ok &= _check(
        "Eq. 5.41 R6 constants",
        math.isclose(gr.R6_k0, HAMEL_R7_CH4_A, rel_tol=0.0, abs_tol=0.0)
        and math.isclose(gr.R6_E_Rg, HAMEL_R7_CH4_E_OVER_RG_K, rel_tol=0.0, abs_tol=0.0),
        f"k0={gr.R6_k0:.6g}, E/Rg={gr.R6_E_Rg:.6g}",
    )

    phase1_cfg = build_phase1_htw_lu_global_nr_reactor_config(load_case_lu())
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

    phase2_cfg = build_phase2_config(case=load_case_lu())
    phase2_reactor = Reactor(phase2_cfg)
    ok &= _check(
        "Phase2 Jacobian delegation",
        phase2_solve_kwargs().get("nr_jacobian_strategy") is None,
        "solve kwargs delegate topology-based default via nr_jacobian_strategy=None",
    )
    ok &= _check(
        "Phase2 resolved Jacobian default",
        phase2_reactor._default_nr_jacobian_strategy()
        in {
            "hybrid_block_tridiag_structured",
            "hybrid_band_plus_side_elements_structured",
            "block_tridiag_structured",
            "band_plus_side_elements_structured",
        },
        f"default={phase2_reactor._default_nr_jacobian_strategy()}",
    )

    default_cfg = ReactorConfig()
    default_reactor = Reactor(default_cfg)
    ok &= _check(
        "Default bed solid state model",
        default_cfg.bed_solid_state_model == HAMEL_BED_SOLID_STATE_MODEL
        and all(c.solid_state_model == HAMEL_BED_SOLID_STATE_MODEL for c in default_reactor.cells),
        f"cfg={default_cfg.bed_solid_state_model}",
    )
    ok &= _check(
        "Default freeboard trajectory coeff",
        default_cfg.freeboard_trajectory_coeff_model == HAMEL_FREEBOARD_TRAJECTORY_COEFF_MODEL,
        f"cfg={default_cfg.freeboard_trajectory_coeff_model}",
    )
    ok &= _check(
        "Default size migration inventory mode",
        default_cfg.size_migration_inventory_mode == HAMEL_SIZE_MIGRATION_INVENTORY_MODE
        and all(
            c.size_migration_inventory_mode == HAMEL_SIZE_MIGRATION_INVENTORY_MODE
            for c in default_reactor.cells
        ),
        f"cfg={default_cfg.size_migration_inventory_mode}",
    )
    ok &= _check(
        "Default Gibbs A1 modes",
        HAMEL_GIBBS_INITIAL_GUESS_MODE == "hamel_a40"
        and HAMEL_GIBBS_TOL_MODE == "holub_vonka_a34",
        f"initial={HAMEL_GIBBS_INITIAL_GUESS_MODE}, tol={HAMEL_GIBBS_TOL_MODE}",
    )
    from src.thermodynamics.gibbs_minimizer import (
        GibbsMinimizer,
        normalize_gibbs_initial_guess_mode,
        normalize_gibbs_tol_mode,
    )

    class _AuditSpeciesDB:
        def get_mu0(self, species: str, T: float) -> float:
            _ = T
            return {"H2": 0.0, "H2O": -220_000.0, "O2": 0.0}[species]

        def get_atom_count(self, species: str, element: str) -> int:
            atoms = {
                "H2": {"H": 2, "O": 0},
                "H2O": {"H": 2, "O": 1},
                "O2": {"H": 0, "O": 2},
            }
            return int(atoms.get(species, {}).get(element, 0))

    _, gibbs_diag = GibbsMinimizer(_AuditSpeciesDB()).solve(
        T=1100.0,
        P=2.0e6,
        elements={"H": 3.0, "O": 1.5},
        candidates=["H2", "H2O", "O2"],
        return_diag=True,
    )
    ok &= _check(
        "GibbsMinimizer production defaults",
        normalize_gibbs_initial_guess_mode(None) == HAMEL_GIBBS_INITIAL_GUESS_MODE
        and normalize_gibbs_tol_mode(None) == HAMEL_GIBBS_TOL_MODE
        and gibbs_diag["initial_guess_mode"] == HAMEL_GIBBS_INITIAL_GUESS_MODE
        and gibbs_diag["tol_mode"] == HAMEL_GIBBS_TOL_MODE
        and gibbs_diag.get("converged") is True,
        f"diag={gibbs_diag.get('converged_by')}",
    )
    ok &= _check(
        "Legacy Gibbs audit branches opt-in",
        normalize_gibbs_initial_guess_mode(LEGACY_GIBBS_INITIAL_GUESS_MODE)
        == LEGACY_GIBBS_INITIAL_GUESS_MODE
        and normalize_gibbs_tol_mode(LEGACY_GIBBS_TOL_MODE) == LEGACY_GIBBS_TOL_MODE,
        f"{LEGACY_GIBBS_INITIAL_GUESS_MODE}+{LEGACY_GIBBS_TOL_MODE}",
    )
    legacy_cfg = ReactorConfig(
        bed_solid_state_model=LEGACY_BED_SOLID_STATE_MODEL,
        freeboard_trajectory_coeff_model=LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL,
        size_migration_inventory_mode=LEGACY_SIZE_MIGRATION_INVENTORY_MODE,
        energy_balance_heat_loss_mode=LEGACY_ENERGY_HEAT_LOSS_MODE,
    )
    legacy_reactor = Reactor(legacy_cfg)
    ok &= _check(
        "Legacy audit branches opt-in",
        all(c.solid_state_model == LEGACY_BED_SOLID_STATE_MODEL for c in legacy_reactor.cells)
        and legacy_cfg.freeboard_trajectory_coeff_model == LEGACY_FREEBOARD_TRAJECTORY_COEFF_MODEL
        and legacy_cfg.size_migration_inventory_mode == LEGACY_SIZE_MIGRATION_INVENTORY_MODE
        and legacy_cfg.energy_balance_heat_loss_mode == LEGACY_ENERGY_HEAT_LOSS_MODE,
        "legacy_stream + stable_mixed_drag_split + total_bed_proxy + inlet_fraction",
    )
    ok &= _check(
        "Default energy heat loss mode",
        default_cfg.energy_balance_heat_loss_mode == HAMEL_ENERGY_HEAT_LOSS_MODE,
        f"cfg={default_cfg.energy_balance_heat_loss_mode}",
    )
    from src.core.energy_balance_defaults import resolve_reactor_Q_wall_W_total

    ok &= _check(
        "Default explicit Q_W from fuel HHV proxy",
        resolve_reactor_Q_wall_W_total(default_cfg) > 0.0,
        f"Q_W={resolve_reactor_Q_wall_W_total(default_cfg)/1e6:.3f} MW",
    )

    print("=" * 80)
    if ok:
        print("HAMEL CONSISTENCY AUDIT PASSED")
        return 0
    print("HAMEL CONSISTENCY AUDIT FAILED")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
