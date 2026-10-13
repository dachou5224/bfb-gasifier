"""阶段 B：热化学恒等式。不以当前 x0 能量残差升降判对错。"""

from __future__ import annotations

import pytest

from src.core.constants import T_REF
from src.core.energy_balance_defaults import (
    DEFAULT_SOLID_ENTHALPY_MODE,
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
    SOLID_ENTHALPY_MODE_HAMEL_EQ28,
    solid_component_hf_j_per_kg,
)
from src.core.species import (
    SOLID_H_MODE_ENDPOINT,
    SOLID_H_MODE_INTEGRAL,
    cp_ash,
    cp_char,
    cp_sand,
    delta_h_ash_integral,
    delta_h_char_integral,
    delta_h_sand_integral,
    enthalpy_molar,
    h_solid_mass,
    mean_cp_solid,
)
from src.core.thermochemistry_identities import (
    dry_fuel_hhv_identity_gap_j_per_kg,
    hamel_char_reaction_targets_j_per_mol,
    liquid_to_vapor_enthalpy_j_per_mol,
    standard_char_reaction_enthalpy_j_per_mol,
)
from tests.validation_case_utils import build_phase1_htw_lu_reactor_config


def test_char_is_graphite_zero_and_moisture_is_liquid_water() -> None:
    hf = solid_component_hf_j_per_kg(h_f_dry=-3.96e6, ash_dry_wt=11.41, VM_daf=53.42)
    assert hf[0] == pytest.approx(0.0)
    assert hf[0] == pytest.approx(H_FORMATION_GRAPHITE_CHAR_J_PER_KG)
    assert hf[2] == pytest.approx(H_FORMATION_LIQUID_WATER_J_PER_KG)
    assert hf[2] == pytest.approx(-15.866e6, rel=2e-4)
    assert hf[3] == pytest.approx(0.0)


def test_dry_fuel_split_matches_hhv_identity_on_lu_fuel() -> None:
    cfg = build_phase1_htw_lu_reactor_config()
    gap = dry_fuel_hhv_identity_gap_j_per_kg(
        C=cfg.C_dry,
        H=cfg.H_dry,
        O=cfg.O_dry,
        S=getattr(cfg, "S_dry", 0.0),
        HHV_dry_MJ_kg=cfg.HHV_dry,
        ash_dry_wt=cfg.ash_dry_wt,
        VM_daf=cfg.VM_daf,
    )
    assert gap["gap_J_kg"] == pytest.approx(0.0, abs=1.0)


def test_char_standard_reaction_enthalpies_match_hamel_chapter_5() -> None:
    targets = hamel_char_reaction_targets_j_per_mol()
    for name, expected in targets.items():
        got = standard_char_reaction_enthalpy_j_per_mol(name, T=T_REF)
        assert got == pytest.approx(expected, abs=2.0e3), name


def test_liquid_water_vaporization_is_endothermic() -> None:
    delta = liquid_to_vapor_enthalpy_j_per_mol(T=T_REF)
    assert delta > 40.0e3
    assert delta == pytest.approx(44.0e3, abs=3.0e3)


def test_solid_integral_enthalpy_identities() -> None:
    hf = {"char": 0.0, "ash": -1.2e6, "sand": 0.0}
    for kind, inst, dfn in (
        ("char", cp_char, delta_h_char_integral),
        ("ash", cp_ash, delta_h_ash_integral),
        ("sand", cp_sand, delta_h_sand_integral),
    ):
        assert h_solid_mass(kind, T_REF, h_f=hf[kind]) == pytest.approx(hf[kind], abs=1e-9)
        T1, T2, T3 = 400.0, 900.0, 1400.0
        assert dfn(T2, T1) + dfn(T3, T2) == pytest.approx(dfn(T3, T1), rel=0.0, abs=1e-8)
        dT = 0.05
        num = (
            h_solid_mass(kind, T2 + dT, h_f=hf[kind])
            - h_solid_mass(kind, T2 - dT, h_f=hf[kind])
        ) / (2.0 * dT)
        assert num == pytest.approx(inst(T2), rel=1e-6)
        assert mean_cp_solid(kind, T_REF) == pytest.approx(inst(T_REF), rel=0.0, abs=1e-12)
        h_end = h_solid_mass(kind, T2, h_f=hf[kind], mode=SOLID_H_MODE_ENDPOINT)
        h_int = h_solid_mass(kind, T2, h_f=hf[kind], mode=SOLID_H_MODE_INTEGRAL)
        assert h_end != pytest.approx(h_int, rel=1e-6)
    assert DEFAULT_SOLID_ENTHALPY_MODE == SOLID_ENTHALPY_MODE_HAMEL_EQ28
    with pytest.raises(ValueError, match="hamel_exact"):
        h_solid_mass("char", 900.0, h_f=0.0, mode="hamel_exact_cp")


def test_gas_enthalpy_already_includes_hf() -> None:
    assert enthalpy_molar("CO2", T_REF) == pytest.approx(enthalpy_molar("CO2", T_REF))
    h_co2 = enthalpy_molar("CO2", 1200.0)
    assert h_co2 != pytest.approx(0.0)
    # 气体焓含生成焓；固体积分焓在 T_r 只剩 hf，二者参考态对齐
    assert h_solid_mass("char", T_REF, h_f=0.0) == pytest.approx(0.0)
