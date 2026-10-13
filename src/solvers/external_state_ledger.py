"""外部元素与焓账。从状态和 q 组装进出流，不把残差相加。

氮、硫没有 NH3、H2S 求解组分，不声称全元素闭合。
干料 ``h_f`` 是生成焓字段，不是高位或低位热值，也不用来抵消 ``Q_W``。

Source: docs/audit_reports/wp1_wp4_review_2026-09-28.md 工作包 C
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.elemental_ledger import ATOMIC_MASS_KG_PER_MOL, H2O_MOLAR_MASS_KG_PER_MOL
from src.core.energy_balance_defaults import (
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
)
from src.core.species import (
    SOLID_H_MODE_INTEGRAL,
    TAR_SURROGATE_FORMULA,
    enthalpy_molar,
    get_tar_component_mapping,
    h_solid_mass,
)
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import edges_from_x
from src.solvers.wirsum_1997_named_f import (
    NamedStageQ,
    _H_major,
    _H_tar,
    ash_feed_kg_s,
    char_face_flows_kg_s,
    unpack_named,
)
from src.solvers.wirsum_1997_pyrolysis_ledger import unreleased_vm_enthalpy_W
from src.solvers.wirsum_1997_stiffness_full_gas import MAJOR
_ELEMENTS = ("C", "H", "O", "N", "S")


def major_atoms_mol_s(n_mol_s: npt.NDArray[np.float64]) -> dict[str, float]:
    """MAJOR 摩尔流的原子流 [mol-atom/s]。顺序 CO, CO2, H2, H2O, CH4, O2, N2。"""
    n = np.asarray(n_mol_s, dtype=np.float64)
    assert n.shape == (7,), f"MAJOR flow expected (7,), got {n.shape}"
    return {
        "C": float(n[0] + n[1] + n[4]),
        "H": float(2.0 * n[2] + 2.0 * n[3] + 4.0 * n[4]),
        "O": float(n[0] + 2.0 * n[1] + n[3] + 2.0 * n[5]),
        "N": float(2.0 * n[6]),
        "S": 0.0,
    }


def tar_atoms_mol_s(n_tar1: float, n_tar2: float, *, fuel_type: str = "coal") -> dict[str, float]:
    """焦油代理分子的 C、H [mol-atom/s]。煤是苯和萘，不含 O、N、S。"""
    mapping = get_tar_component_mapping(fuel_type)
    out = {el: 0.0 for el in _ELEMENTS}
    for label, n in (("TAR1", n_tar1), ("TAR2", n_tar2)):
        c_count, h_count = TAR_SURROGATE_FORMULA[mapping[label]]
        out["C"] += float(c_count) * float(n)
        out["H"] += float(h_count) * float(n)
    return out


def _add(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    return {el: float(a.get(el, 0.0)) + float(b.get(el, 0.0)) for el in _ELEMENTS}


def _sub(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    return {el: float(a.get(el, 0.0)) - float(b.get(el, 0.0)) for el in _ELEMENTS}


def external_imbalance(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
) -> dict[str, Any]:
    """任意试探状态的外部不平衡。不求解，也不读取残差向量。"""
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    sp = nq.energy.species
    feed = nq.energy.feed
    ed = edges_from_x(z[:39], sp)
    n_water = float(feed.m_moist_kg_s) / H2O_MOLAR_MASS_KG_PER_MOL
    fuel = {
        "C": float(feed.nC_fuel),
        "H": float(feed.nH_fuel),
        "O": float(feed.nO_fuel),
        "N": float(feed.nN_fuel),
        "S": float(feed.nS_fuel),
    }
    moisture = {"C": 0.0, "H": 2.0 * n_water, "O": n_water, "N": 0.0, "S": 0.0}
    oxidant = _add(major_atoms_mol_s(sp.N_in_d), major_atoms_mol_s(sp.N_in_b))
    inlet = _add(_add(fuel, moisture), oxidant)
    gas = major_atoms_mol_s(u["N_c2"])
    if int(nq.stage) >= 10:
        tar = tar_atoms_mol_s(float(u["N_tar1_fb"]), float(u["N_tar2_fb"]))
        n_tar = {"TAR1": float(u["N_tar1_fb"]), "TAR2": float(u["N_tar2_fb"])}
    else:
        tar = tar_atoms_mol_s(
            float(u["N_tar1"]) + float(u["N_tar1_b"]),
            float(u["N_tar2"]) + float(u["N_tar2_b"]),
        )
        n_tar = {
            "TAR1": float(u["N_tar1"]) + float(u["N_tar1_b"]),
            "TAR2": float(u["N_tar2"]) + float(u["N_tar2_b"]),
        }
    n_char = float(ed["m_fly"]) / ATOMIC_MASS_KG_PER_MOL["C"]
    char = {"C": n_char, "H": 0.0, "O": 0.0, "N": 0.0, "S": 0.0}
    unreleased = {
        "C": float(sp.nC_unreleased),
        "H": float(sp.nH_unreleased),
        "O": float(sp.nO_unreleased),
        "N": float(sp.nN_unreleased),
        "S": float(sp.nS_unreleased),
    }
    n_nh3 = float(getattr(sp, "n_NH3_bypass_mol_s", 0.0))
    n_h2s = float(getattr(sp, "n_H2S_bypass_mol_s", 0.0))
    bypass = {
        "C": 0.0,
        "H": 3.0 * n_nh3 + 2.0 * n_h2s,
        "O": 0.0,
        "N": n_nh3,
        "S": n_h2s,
    }
    outlet = _add(_add(_add(_add(gas, tar), char), unreleased), bypass)
    gap = _sub(inlet, outlet)
    m_eject, m_return, m_leave = char_face_flows_kg_s(z, nq)
    up_kg_s = float(np.sum(m_eject))
    down_kg_s = float(np.sum(m_return))
    leave_kg_s = float(np.sum(m_leave))
    t_exit = float(z[-1])
    h_char = h_solid_mass(
        "char", t_exit, h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG, mode=SOLID_H_MODE_INTEGRAL,
    )
    h_ash = h_solid_mass("ash", t_exit, h_f=-1.2e6, mode=SOLID_H_MODE_INTEGRAL)
    k = sp.k
    fly_frac = float(k.K_auf) / max(float(k.K_auf) + float(k.K_ab), 1.0e-30)
    m_fly_ash = fly_frac * float(sp.K_eject[-1]) * float(u["M_ash"]) if int(nq.stage) >= 3 else 0.0
    h_fuel = float(feed.m_dry_kg_s) * float(feed.h_f_dry_J_kg)
    h_moist = float(feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    h_gas_in = _H_major(sp.N_in_d, nq.energy.T_in) + _H_major(sp.N_in_b, nq.energy.T_in)
    h_gas_out = _H_major(u["N_c2"], t_exit)
    h_tar_out = _H_tar(n_tar, t_exit)
    h_char_out = float(ed["m_fly"]) * h_char
    h_ash_out = m_fly_ash * h_ash
    h_unreleased = unreleased_vm_enthalpy_W(feed, float(sp.m_unreleased_kg_s), t_exit)
    from src.solvers.wirsum_1997_named_f import minor_bypass_enthalpy_W

    h_bypass = minor_bypass_enthalpy_W(sp, t_exit)
    q_w = float(nq.Q_W) if int(nq.stage) >= 5 else 0.0
    h_in = h_gas_in + h_fuel + h_moist
    h_out = h_gas_out + h_tar_out + h_char_out + h_ash_out + h_unreleased + h_bypass
    return {
        "units": {"elements": "mol-atom/s", "enthalpy": "W", "mass": "kg/s"},
        "inlet_mol_s": {"fuel": fuel, "moisture": moisture, "oxidant": oxidant, "total": inlet},
        "outlet_mol_s": {
            "gas": gas,
            "tar": tar,
            "char_fly": char,
            "unreleased": unreleased,
            "minor_bypass": bypass,
            "total": outlet,
        },
        "ash_out_kg_s": m_fly_ash,
        "ash_feed_kg_s": ash_feed_kg_s(nq) if int(nq.stage) >= 3 else 0.0,
        "gap_mol_s": gap,
        "interface_char": {
            "bed_exit_up_kg_s": up_kg_s,
            "freeboard_inlet_up_kg_s": up_kg_s,
            "up_stream_paired": True,
            "freeboard_return_after_reaction_kg_s": down_kg_s,
            "bed_inlet_down_kg_s": down_kg_s,
            "down_stream_paired": True,
            "leave_after_reaction_kg_s": leave_kg_s,
            "up_must_equal_down": False,
            "excluded_from_external_gap": True,
        },
        "enthalpy_W": {
            "gas_in": h_gas_in,
            "fuel_formation": h_fuel,
            "moisture_formation": h_moist,
            "gas_out": h_gas_out,
            "tar_out": h_tar_out,
            "char_fly": h_char_out,
            "ash_fly": h_ash_out,
            "unreleased": h_unreleased,
            "minor_bypass": h_bypass,
            "Q_W": q_w,
            "gap": h_in - h_out - q_w,
        },
        "minor_bypass": {
            "model": "inert_nh3_h2s_no_reaction",
            "covers_subsequent_reactions": False,
            "full_NS_model": False,
            "n_NH3_mol_s": n_nh3,
            "n_H2S_mol_s": n_h2s,
            "enters_major_density": False,
            "enters_particle_residence": False,
            "molar_flow_treatment": "atom and enthalpy leave with the exit stream; not added to the seven-species molar total used for gas velocity",
        },
        "full_element_certification": False,
        "h_o_gap_not_waived": True,
        "unmodeled": {
            "N": "allocated NH3 is an inert bypass, not an NH3 reaction model",
            "S": "allocated H2S is an inert bypass, not an H2S reaction model",
        },
        "formation_enthalpy_is_not_a_heating_value": True,
        "Q_is_not_compensated_by_formation_enthalpy": True,
    }
