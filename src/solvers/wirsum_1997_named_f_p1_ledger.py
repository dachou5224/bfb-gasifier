"""P1 独立物料/焓账。从 (x, q, feed) 流股加和，不把内部残差 F 当作闭合证据。

Source: docs/audit_reports/named_f_e8_review_2026-09-16.md §2.5–2.7, §5 P1
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.elemental_ledger import H2O_MOLAR_MASS_KG_PER_MOL
from src.core.energy_balance_defaults import (
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
)
from src.core.species import (
    MW_C_KG_PER_MOL,
    SOLID_H_MODE_INTEGRAL,
    TAR_SURROGATE_FORMULA,
    configure_tar_components_by_fuel,
    get_tar_component_mapping,
    h_solid_mass,
)
from src.solvers.wirsum_1997_coupled_f_energy import _H_major, _H_tar, q_at_T
from src.solvers.wirsum_1997_coupled_f_minimal import _elem_C, _elem_H, _elem_O
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import IDX_COARSE, edges_from_x
from src.solvers.wirsum_1997_named_f import (
    R10_FUEL_TYPE,
    NamedStageQ,
    ash_feed_kg_s,
    energy_row_named,
    tar_concentration_diagnostics,
    unpack_named,
)
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    FuelAtomFeed,
    map_one_state,
    unreleased_vm_enthalpy_W,
)

_N_SPECIES: int = 39
_ASH_HF_NAMED_J_PER_KG: float = -1.2e6  # 与 energy_row_named 同一灰生成焓


def _tar_atoms(n_tar1: float, n_tar2: float, *, fuel_type: str = R10_FUEL_TYPE) -> dict[str, float]:
    configure_tar_components_by_fuel(fuel_type)
    tar_map = get_tar_component_mapping(fuel_type)
    out = {"C": 0.0, "H": 0.0, "O": 0.0, "N": 0.0, "S": 0.0}
    for label, n in (("TAR1", n_tar1), ("TAR2", n_tar2)):
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_map[label]]
        out["C"] += float(c_tar) * float(n)
        out["H"] += float(h_tar) * float(n)
    return out


def mapping_ledger(
    feed: FuelAtomFeed,
    *,
    P_Pa: float,
    T_K: float,
    tau_s: float,
    X_end: float | None = None,
    occupancy_source: str = "named_p1_ledger",
) -> dict[str, Any]:
    """热解映射账：原料 / 释放物 / 剩余物 / 焓，同一 T、同一碳基准。"""
    from src.solvers.wirsum_1997_named_f import daem_X_at_T

    X = float(X_end) if X_end is not None else daem_X_at_T(T_K, tau_s)
    mapped = map_one_state(
        feed,
        P_Pa=float(P_Pa),
        T_K=float(T_K),
        X_vm=np.array([0.0, X], dtype=np.float64),
        tau_s=np.array([float(tau_s)], dtype=np.float64),
        occupancy_source=occupancy_source,
        fuel_type="coal",
    )
    nC_tar = 0.0
    for cell in mapped["cells"]:
        n_py = cell["n_py_mol_s"]
        ta = _tar_atoms(float(n_py.get("TAR1", 0.0)), float(n_py.get("TAR2", 0.0)))
        nC_tar += ta["C"]
    split = mapped["split"]
    nC_fuel = float(feed.nC_fuel)
    return {
        "mapped": mapped,
        "X_end": float(mapped["X_end"]),
        "carbon_basis": mapped["carbon_basis"],
        "remainder_destination": mapped["remainder_destination"],
        "nC_fuel": nC_fuel,
        "nC_koks": float(mapped["nC_koks"]),
        "nC_tar_released": nC_tar,
        "nC_gas_plus_tar": float(mapped["nC_gas_plus_tar"]),
        "nC_unreleased": float(mapped["nC_unreleased_vm"]),
        "x_C_koks_fuel_waf": float(split["x_C_koks_fuel_waf"]),
        "x_C_teer_fuel_waf": float(split["x_C_teer_fuel_waf"]),
        "x_C_koks_recovered": float(mapped["nC_koks"]) / max(nC_fuel, 1.0e-30),
        "x_C_teer_recovered": nC_tar / max(nC_fuel, 1.0e-30),
        "atom_C_gap": float(mapped["atom_C_gap"]),
        "atom_H_gap": float(mapped["atom_H_gap"]),
        "atom_O_gap": float(mapped["atom_O_gap"]),
        "hegermann_coke_basis_gap": float(mapped["hegermann_coke_basis_gap"]),
        "hegermann_koks_minus_fixed": float(mapped["hegermann_koks_minus_fixed"]),
        "hegermann": mapped["hegermann"],
        "H_in_W": float(mapped["H_in_W"]),
        "H_prod_gas_tar_W": float(mapped["H_prod_gas_tar_W"]),
        "H_koks_W": float(mapped["H_koks_W"]),
        "H_unreleased_W": float(mapped["H_unreleased_W"]),
        "H_ash_W": float(mapped["H_ash_W"]),
        "H_map_gap_W": float(mapped["H_map_gap_W"]),
        "m_ash_kg_s": float(mapped["m_ash_kg_s"]),
        "m_unreleased_kg_s": float(mapped["m_unreleased_kg_s"]),
        "not_from_F": True,
    }


def named_stream_ledger(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
) -> dict[str, Any]:
    """联立核进出口流股账。用 x 与当前 q，不读 F 向量。"""
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-1])
    qe = nq.energy
    qT = q_at_T(qe.species, T)
    u = unpack_named(z, nq.stage)
    ed = edges_from_x(z[:_N_SPECIES], qT)
    feed = qe.feed
    N_in = np.asarray(qe.species.N_in_d, dtype=np.float64) + np.asarray(
        qe.species.N_in_b, dtype=np.float64
    )
    n_moist = float(feed.m_moist_kg_s) / H2O_MOLAR_MASS_KG_PER_MOL
    if int(nq.stage) >= 6:
        n_tar1 = float(u["N_tar1"])
        n_tar2 = float(u["N_tar2"])
    else:
        n_tar1 = float(qe.n_tar.get("TAR1", 0.0))
        n_tar2 = float(qe.n_tar.get("TAR2", 0.0))
    tar = _tar_atoms(n_tar1, n_tar2)
    C_in = float(feed.nC_fuel) + _elem_C(N_in)
    H_in_el = float(feed.nH_fuel) + _elem_H(N_in) + 2.0 * n_moist
    O_in_el = float(feed.nO_fuel) + _elem_O(N_in) + 1.0 * n_moist
    C_out = (
        _elem_C(u["N_c2"])
        + tar["C"]
        + float(ed["m_fly"]) / MW_C_KG_PER_MOL
        + float(qe.species.nC_unreleased)
    )
    H_out_el = _elem_H(u["N_c2"]) + tar["H"] + float(qe.species.nH_unreleased)
    O_out_el = _elem_O(u["N_c2"]) + tar["O"] + float(qe.species.nO_unreleased)
    ash_in = ash_feed_kg_s(nq)
    k = qe.species.k
    fly_frac = float(k.K_auf) / max(float(k.K_auf) + float(k.K_ab), 1.0e-30)
    K_e = float(qe.species.K_eject[IDX_COARSE])
    ash_out = fly_frac * K_e * float(u["M_ash"]) if int(nq.stage) >= 3 else 0.0
    n_tar = n_tar1 + n_tar2
    tar_eos = tar_concentration_diagnostics(
        u["N_d"], u["N_b"], n_tar, T, float(qe.species.geo_ex.P)
    )
    h_char = h_solid_mass(
        "char", T, h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG, mode=SOLID_H_MODE_INTEGRAL
    )
    H_in_W = (
        _H_major(qe.species.N_in_d, qe.T_in)
        + _H_major(qe.species.N_in_b, qe.T_in)
        + float(feed.m_dry_kg_s) * float(feed.h_f_dry_J_kg)
        + float(feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    )
    n_tar_dict = {"TAR1": n_tar1, "TAR2": n_tar2}
    H_out_W = (
        _H_major(u["N_c2"], T)
        + float(ed["m_fly"]) * h_char
        + _H_tar(n_tar_dict, T)
        + unreleased_vm_enthalpy_W(feed, float(qe.species.m_unreleased_kg_s), T)
    )
    if int(nq.stage) >= 3:
        h_ash = h_solid_mass(
            "ash", T, h_f=_ASH_HF_NAMED_J_PER_KG, mode=SOLID_H_MODE_INTEGRAL
        )
        H_out_W = H_out_W + ash_out * h_ash
    Q_W = float(nq.Q_W) if int(nq.stage) >= 5 else 0.0
    H_gap_W = float(H_in_W - H_out_W - Q_W)
    return {
        "T_K": T,
        "C_in": C_in,
        "C_out": C_out,
        "C_gap": abs(C_out - C_in),
        "H_in": H_in_el,
        "H_out": H_out_el,
        "H_gap": abs(H_out_el - H_in_el),
        "O_in": O_in_el,
        "O_out": O_out_el,
        "O_gap": abs(O_out_el - O_in_el),
        "ash_in_kg_s": ash_in,
        "ash_out_kg_s": ash_out,
        "ash_gap_kg_s": abs(ash_out - ash_in),
        "nC_unreleased": float(qe.species.nC_unreleased),
        "nC_koks": float(qe.species.nC_koks),
        "m_unreleased_kg_s": float(qe.species.m_unreleased_kg_s),
        "carbon_basis": qe.species.carbon_basis,
        "remainder_destination": qe.species.remainder_destination,
        "H_in_W": H_in_W,
        "H_out_W": H_out_W,
        "Q_W": Q_W,
        "H_enthalpy_gap_W": abs(H_gap_W),
        "energy_row_W": float(energy_row_named(z, nq)),
        "tar_eos": tar_eos,
        "not_from_F": True,
    }
