"""单 cell 干燥/热解源项组装。

Ref: Hamel (1999) Eq. 4.4, 4.9, 4.10-4.12
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
import numpy.typing as npt

from src.core.elemental_ledger import (
    ATOMIC_MASS_KG_PER_MOL,
    H2O_MOLAR_MASS_KG_PER_MOL,
    vm_element_molar_rates,
)
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    N_GAS,
    TAR_SURROGATE_FORMULA,
    calc_tar_surrogate_fractions,
    get_tar_component_mapping,
)
from src.thermal.devolatilization import daem_conversion_radial
from src.thermal.drying import solve_drying_CN


@dataclass
class PyrolysisSourceBundle:
    """干燥/热解源项输出。"""

    gas_source: npt.NDArray[np.float64]
    solid_sink: npt.NDArray[np.float64]
    x_dry: float
    x_vm: float


def allocate_pyrolysis_products_elemental(
    *,
    nC: float,
    nH: float,
    nO: float,
    fuel_type: str,
    pyrolysis_tar_carbon_frac: float,
    target_tar_hc_ratio: float | None = None,
    nN: float = 0.0,
    nS: float = 0.0,
) -> Dict[str, float]:
    """按元素守恒分配热解产物。

    Ref: Hamel (1999) Eq. 4.10-4.12
    """
    out = {sp: 0.0 for sp in GAS_SPECIES}
    out["NH3"] = max(float(nN), 0.0)
    out["H2S"] = max(float(nS), 0.0)
    nH_eff = max(float(nH) - 3.0 * float(nN) - 2.0 * float(nS), 0.0)

    f_tar = float(np.clip(pyrolysis_tar_carbon_frac, 0.0, 0.9))
    tar_mapping = get_tar_component_mapping(fuel_type)
    try:
        tar_fractions = calc_tar_surrogate_fractions(fuel_type, target_hc_ratio=target_tar_hc_ratio)
    except ValueError:
        # 当 ultimate 分析给出的 H/C 略超代理可表示区间时，回退到该燃料默认 H/C。
        tar_fractions = calc_tar_surrogate_fractions(fuel_type)
    tar_components = []
    for label in ("TAR1", "TAR2"):
        surrogate = tar_mapping[label]
        frac = float(tar_fractions.get(surrogate, 0.0))
        c_tar_i, h_tar_i = TAR_SURROGATE_FORMULA[surrogate]
        tar_components.append((label, frac, c_tar_i, h_tar_i))

    tar_carbon = f_tar * float(nC)
    mean_c_tar = sum(frac * c_tar_i for _, frac, c_tar_i, _ in tar_components)
    nt_total = tar_carbon / max(mean_c_tar, 1e-12)
    n_tar = {label: frac * nt_total for label, frac, _, _ in tar_components}
    carbon_to_tar = sum(c_tar_i * n_tar[label] for label, _, c_tar_i, _ in tar_components)
    hydrogen_to_tar = sum(h_tar_i * n_tar[label] for label, _, _, h_tar_i in tar_components)

    # 主产物优先：CO / H2 / CH4；随后用 CO2 / H2O 做元素守恒闭合 fallback。
    c_left = max(float(nC) - carbon_to_tar, 0.0)
    h_left = max(nH_eff - hydrogen_to_tar, 0.0)
    o_left = max(float(nO), 0.0)

    nco = min(c_left, o_left)
    c_left -= nco
    o_left -= nco

    nch4 = min(c_left, h_left / 4.0)
    c_left -= nch4
    h_left -= 4.0 * nch4

    nco2 = min(c_left, o_left / 2.0)
    c_left -= nco2
    o_left -= 2.0 * nco2

    # 若仍有碳且有氧，回退到 CO（保证 C 不超支）。
    nco_fb = min(c_left, o_left)
    nco += nco_fb
    c_left -= nco_fb
    o_left -= nco_fb

    nh2o = min(h_left / 2.0, o_left)
    h_left -= 2.0 * nh2o
    o_left -= nh2o

    nh2 = max(h_left, 0.0) / 2.0

    out.update(
        {
            "CO": max(nco, 0.0),
            "CO2": max(nco2, 0.0),
            "CH4": max(nch4, 0.0),
            "H2O": max(nh2o, 0.0),
            "H2": max(nh2, 0.0),
            "TAR1": max(n_tar["TAR1"], 0.0),
            "TAR2": max(n_tar["TAR2"], 0.0),
        }
    )
    return out


def calc_drying_pyrolysis_sources(
    *,
    tau: float,
    T: float,
    P: float,
    T_init: float | None = None,
    d_p: float,
    moisture_wt: float,
    ash_dry_wt: float,
    C_dry: float,
    H_dry: float,
    O_dry: float,
    nitrogen_fraction: float,
    sulfur_fraction: float,
    sulfur_volatile_frac: float,
    pyrolysis_tar_carbon_frac: float,
    fuel_type: str,
    m_vm_in: float,
    m_moist_in: float,
    solid_shape: tuple[int, int],
    char_index: int | None = None,
    vm_index: int,
    moisture_index: int,
) -> PyrolysisSourceBundle:
    """计算干燥/热解气相源项与固相扣减。

    Ref: Hamel (1999) Eq. 4.4, 4.9, 4.10-4.12

    参数 `T_init` 用于承接 Vorabrechnung 的轴向固体温度历史（上一 cell 固体入口温度）；
    若未提供则回退 300 K。
    """
    assert len(solid_shape) == 2, f"solid_shape must be 2D, got {solid_shape}"
    if max(float(m_vm_in), 0.0) <= 1e-12 and max(float(m_moist_in), 0.0) <= 1e-12:
        return PyrolysisSourceBundle(
            gas_source=np.zeros(N_GAS, dtype=np.float64),
            solid_sink=np.zeros(solid_shape, dtype=np.float64),
            x_dry=0.0,
            x_vm=0.0,
        )
    t_init = 300.0 if T_init is None else float(T_init)
    dry = solve_drying_CN(
        d_p,
        T,
        t_init,
        moisture_wt,
        max(float(tau), 0.05),
        Nr=12,
        Nt=80,
        pressure_pa=float(P),
        return_history=True,
    )
    radius = max(float(d_p) * 0.5, 1e-12)
    x_vm = daem_conversion_radial(dry["T_history_rt"], dry["t"], dry["r_nodes"]) * (
        1.0 - (min(float(dry["r_evap"][-1]), radius) / radius) ** 3
    )
    x_vm = float(np.clip(x_vm, 0.0, 1.0))
    m_vm_rel = float(m_vm_in) * x_vm

    gas_source = np.zeros(N_GAS, dtype=np.float64)
    gas_source[GAS_SPECIES_INDEX["H2O"]] = (
        float(m_moist_in) * float(dry["X_dry"][-1])
    ) / H2O_MOLAR_MASS_KG_PER_MOL

    vm_atoms = vm_element_molar_rates(
        m_vm=m_vm_rel,
        ash_dry_wt=ash_dry_wt,
        C_dry=C_dry,
        H_dry=H_dry,
        O_dry=O_dry,
        nitrogen_fraction=nitrogen_fraction,
        sulfur_fraction=sulfur_fraction,
        sulfur_volatile_frac=sulfur_volatile_frac,
    )

    # Hamel tar 代理比例按目标 H/C 在线求解；目标默认取当前 fuel ultimate（dry）的原子比。
    target_tar_hc_ratio = None
    if C_dry > 1e-12 and H_dry >= 0.0:
        target_tar_hc_ratio = float((H_dry / 1.00794) / max(C_dry / 12.011, 1e-12))

    prod = allocate_pyrolysis_products_elemental(
        nC=vm_atoms["C"],
        nH=vm_atoms["H"],
        nO=vm_atoms["O"],
        nN=vm_atoms["N"],
        nS=vm_atoms["S"],
        fuel_type=fuel_type,
        pyrolysis_tar_carbon_frac=pyrolysis_tar_carbon_frac,
        target_tar_hc_ratio=target_tar_hc_ratio,
    )
    for sp, value in prod.items():
        gas_source[GAS_SPECIES_INDEX[sp]] += float(value)

    solid_sink = np.zeros(solid_shape, dtype=np.float64)
    carbon_to_gas = (
        float(prod.get("CO", 0.0))
        + float(prod.get("CO2", 0.0))
        + float(prod.get("CH4", 0.0))
    )
    tar_mapping = get_tar_component_mapping(fuel_type)
    for tar_label in ("TAR1", "TAR2"):
        tar_surr = tar_mapping[tar_label]
        c_tar, _ = TAR_SURROGATE_FORMULA[tar_surr]
        carbon_to_gas += c_tar * float(prod.get(tar_label, 0.0))
    carbon_to_char = max(vm_atoms["C"] - carbon_to_gas, 0.0)
    if char_index is not None:
        solid_sink[:, char_index] += carbon_to_char * ATOMIC_MASS_KG_PER_MOL["C"]
    solid_sink[:, moisture_index] -= float(m_moist_in) * float(dry["X_dry"][-1])
    solid_sink[:, vm_index] -= m_vm_rel

    return PyrolysisSourceBundle(
        gas_source=gas_source,
        solid_sink=solid_sink,
        x_dry=float(np.clip(dry["X_dry"][-1], 0.0, 1.0)),
        x_vm=x_vm,
    )
