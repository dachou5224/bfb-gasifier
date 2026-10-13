"""单 cell 干燥/热解源项组装。

Ref: Hamel (1999) Eq. 4.4, 4.9 (干燥/DAEM); §4.3 (热解气 Gibbs 组成); §4.2 (Teer/N/S)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

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
    fuel_tar_target_hc_ratio,
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


def _allocate_pyrolysis_tar_and_minors(
    *,
    nC: float,
    nH: float,
    nO: float,
    nN: float,
    nS: float,
    fuel_type: str,
    pyrolysis_tar_carbon_frac: float,
    target_tar_hc_ratio: float | None,
    pyrolysis_koks_carbon_frac: float = 0.0,
) -> tuple[Dict[str, float], float, float, float]:
    """分出 Koks(C) / Teer / NH₃ / H₂S，返回 (产物, c_gas, h_gas, o_gas) [mol-atom/s]。

    Hamel §4.2–4.3：分率相对**本函数输入 ``nC``**。Eq.4.21–4.22 字面量是燃料 waf 碳。
    Koks 为纯碳，不带走 H。
    剩余 C/H/O 进气相 Gibbs。
    """
    out = {sp: 0.0 for sp in GAS_SPECIES}
    out["NH3"] = max(float(nN), 0.0)
    out["H2S"] = max(float(nS), 0.0)
    nH_eff = max(float(nH) - 3.0 * float(nN) - 2.0 * float(nS), 0.0)

    f_koks = float(np.clip(pyrolysis_koks_carbon_frac, 0.0, 0.95))
    f_tar = float(np.clip(pyrolysis_tar_carbon_frac, 0.0, 0.9))
    if f_koks + f_tar > 0.95:
        f_tar = max(0.0, 0.95 - f_koks)

    carbon_to_koks = f_koks * float(nC)

    tar_mapping = get_tar_component_mapping(fuel_type)
    tar_fractions = calc_tar_surrogate_fractions(fuel_type, target_hc_ratio=target_tar_hc_ratio)
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

    out["TAR1"] = max(n_tar["TAR1"], 0.0)
    out["TAR2"] = max(n_tar["TAR2"], 0.0)

    c_gas = max(float(nC) - carbon_to_koks - carbon_to_tar, 0.0)
    h_gas = max(nH_eff - hydrogen_to_tar, 0.0)
    o_gas = max(float(nO), 0.0)
    return out, c_gas, h_gas, o_gas


def apply_pyrolysis_co_reduce_h2_inject(
    prod: Dict[str, float],
    inject_frac: float,
    *,
    h2_per_trimmed_co: float = 2.0,
) -> Dict[str, float]:
    """灵敏度：削减 Gibbs 热解 CO，并任意注入 H₂（不增加 CO）。

    不改动 §4.3 Gibbs 求解。相对 CH₄→2H₂+CO 分拆，本映射只降 CO、升 H₂，
    用于检验 bed0 H₂ 争 O₂ 能否在**不额外副产 CO** 时逆转 CO<CO₂。

    对 trim = f×CO：CO -= trim；H₂ += h2_per_trimmed_co×trim；H₂O += trim（O 留在气相）；
    脱除 CO 的 C 视为隐式转 char（仅灵敏度，不回流固相扣减）。
    """
    f = float(np.clip(inject_frac, 0.0, 1.0))
    if f <= 0.0:
        return prod
    co0 = float(max(prod.get("CO", 0.0), 0.0))
    trim = f * co0
    if trim <= 0.0:
        return prod
    h2_gain = float(max(h2_per_trimmed_co, 0.0)) * trim
    out = dict(prod)
    out["CO"] = co0 - trim
    out["H2"] = float(max(out.get("H2", 0.0), 0.0)) + h2_gain
    out["H2O"] = float(max(out.get("H2O", 0.0), 0.0)) + trim
    return out


def apply_pyrolysis_ch4_to_h2_co_split(
    prod: Dict[str, float],
    split_frac: float,
) -> Dict[str, float]:
    """灵敏度：将 Gibbs 热解 CH₄ 的一部分重映射为 2 H₂ + CO（气相 C/H 守恒）。

    不改动 §4.3 Gibbs 求解，仅在分配后做组成扰动，用于检验 bed0 H₂ 与 CO 争 O₂
    （R12 vs R5）对出口 CO/CO₂ 的影响。

    对每 1 mol CH₄ 转化：CH₄ → 2 H₂ + CO。
    """
    f = float(np.clip(split_frac, 0.0, 1.0))
    if f <= 0.0:
        return prod
    ch4 = float(max(prod.get("CH4", 0.0), 0.0))
    move = f * ch4
    if move <= 0.0:
        return prod
    out = dict(prod)
    out["CH4"] = ch4 - move
    out["H2"] = float(max(out.get("H2", 0.0), 0.0)) + 2.0 * move
    out["CO"] = float(max(out.get("CO", 0.0), 0.0)) + move
    return out


def allocate_pyrolysis_products_hamel(
    *,
    T_K: float,
    P_Pa: float,
    nC: float,
    nH: float,
    nO: float,
    fuel_type: str,
    pyrolysis_tar_carbon_frac: float,
    target_tar_hc_ratio: float | None = None,
    nN: float = 0.0,
    nS: float = 0.0,
    pyrolysis_koks_carbon_frac: float = 0.0,
    pyrolysis_gibbs_solver_mode: str = "hamel_reduced",
    pyrolysis_gibbs_temperature_floor_K: float = 750.0,
    pyrolysis_ch4_to_h2_co_split_frac: float = 0.0,
    pyrolysis_co_reduce_h2_inject_frac: float = 0.0,
) -> Dict[str, float]:
    """Hamel §4.3：Koks/Teer/N/S 分出后，气相主组分由封闭 Gibbs（Anhang A1）求组成。

    Ref: Hamel (1999) §4.2 Teer/Koks; §4.3 pyrolysis gas; Anhang A1

    ``nC/nH/nO`` 为 **原子** 摩尔流 [mol/s]；Gibbs 接口使用 ``nH/2``、``nO/2``。
    ``pyrolysis_koks_carbon_frac`` / ``pyrolysis_tar_carbon_frac``：相对本函数
    **输入 ``nC``** 的碳分率。Hegermann 字面 ``x_C,*`` 是燃料 waf 碳基准；
    只有 ``nC`` 也是燃料 waf 碳时二者同义。VM 池须先换算，见
    ``wirsum_1997_pyrolysis_ledger.fractions_on_named_carbon_pool``。
    """
    minors, c_gas, h_gas, o_gas = _allocate_pyrolysis_tar_and_minors(
        nC=nC,
        nH=nH,
        nO=nO,
        nN=nN,
        nS=nS,
        fuel_type=fuel_type,
        pyrolysis_tar_carbon_frac=pyrolysis_tar_carbon_frac,
        target_tar_hc_ratio=target_tar_hc_ratio,
        pyrolysis_koks_carbon_frac=pyrolysis_koks_carbon_frac,
    )
    out = {sp: 0.0 for sp in GAS_SPECIES}
    out.update({sp: float(minors.get(sp, 0.0)) for sp in GAS_SPECIES})

    if c_gas <= 1.0e-20 and h_gas <= 1.0e-20 and o_gas <= 1.0e-20:
        return out

    from src.solvers.vorabrechnung import solve_pyrolysis_gas_gibbs_composition

    gas_major = solve_pyrolysis_gas_gibbs_composition(
        T_K=float(T_K),
        P_Pa=float(P_Pa),
        nC_mol=c_gas,
        nH_mol=h_gas / 2.0,
        nO_mol=o_gas / 2.0,
        solver_mode=str(pyrolysis_gibbs_solver_mode),
    )

    for sp in ("CO", "CO2", "CH4", "H2", "H2O"):
        out[sp] = max(float(gas_major.get(sp, 0.0)), 0.0)
    out = apply_pyrolysis_ch4_to_h2_co_split(out, pyrolysis_ch4_to_h2_co_split_frac)
    return apply_pyrolysis_co_reduce_h2_inject(out, pyrolysis_co_reduce_h2_inject_frac)


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
    T_K: float = 1100.0,
    P_Pa: float = 101325.0,
    pyrolysis_koks_carbon_frac: float = 0.0,
) -> Dict[str, float]:
    """热解产物分配（Hamel §4.3 Gibbs + §4.2 Teer）。"""
    return allocate_pyrolysis_products_hamel(
        T_K=T_K,
        P_Pa=P_Pa,
        nC=nC,
        nH=nH,
        nO=nO,
        nN=nN,
        nS=nS,
        fuel_type=fuel_type,
        pyrolysis_tar_carbon_frac=pyrolysis_tar_carbon_frac,
        pyrolysis_koks_carbon_frac=pyrolysis_koks_carbon_frac,
        target_tar_hc_ratio=target_tar_hc_ratio,
    )


def vorab_budget_atom_mol_s_from_cumulative_vm(
    budget: Any,
    *,
    n_daf_per_c: float,
    s_daf_per_c: float,
    sulfur_volatile_frac: float,
) -> dict[str, float]:
    """Vorab 累积 VM 预算 → 原子摩尔流（与 ``vm_element_molar_rates`` 同口径）。"""
    nC = float(max(getattr(budget, "nC_vm_mol_s", 0.0), 0.0))
    nH = 2.0 * float(max(getattr(budget, "nH_vm_mol_s", 0.0), 0.0))
    nO = 2.0 * float(max(getattr(budget, "nO_vm_mol_s", 0.0), 0.0))
    nN = nC * float(n_daf_per_c) if nC > 0.0 else 0.0
    nS = nC * float(s_daf_per_c) * float(sulfur_volatile_frac) if nC > 0.0 else 0.0
    return {"nC": nC, "nH": nH, "nO": nO, "nN": nN, "nS": nS}


def incremental_pyrolysis_products_from_budgets(
    budget_curr: Any,
    budget_prev: Any,
    *,
    T_K: float,
    P_Pa: float,
    fuel_type: str,
    pyrolysis_tar_carbon_frac: float,
    n_daf_per_c: float,
    s_daf_per_c: float,
    sulfur_volatile_frac: float,
    target_tar_hc_ratio: float | None = None,
    pyrolysis_koks_carbon_frac: float = 0.0,
    pyrolysis_gibbs_solver_mode: str = "hamel_reduced",
    pyrolysis_gibbs_temperature_floor_K: float = 750.0,
    pyrolysis_ch4_to_h2_co_split_frac: float = 0.0,
    pyrolysis_co_reduce_h2_inject_frac: float = 0.0,
) -> Dict[str, float]:
    """相邻 Vorab 预算的累积热解产物差分（避免单格大块 VM 增量 Gibbs 塌缩）。"""
    curr_atoms = vorab_budget_atom_mol_s_from_cumulative_vm(
        budget_curr,
        n_daf_per_c=n_daf_per_c,
        s_daf_per_c=s_daf_per_c,
        sulfur_volatile_frac=sulfur_volatile_frac,
    )
    prev_atoms = vorab_budget_atom_mol_s_from_cumulative_vm(
        budget_prev,
        n_daf_per_c=n_daf_per_c,
        s_daf_per_c=s_daf_per_c,
        sulfur_volatile_frac=sulfur_volatile_frac,
    )

    _gibbs_kw = dict(
        target_tar_hc_ratio=target_tar_hc_ratio,
        pyrolysis_gibbs_solver_mode=pyrolysis_gibbs_solver_mode,
        pyrolysis_gibbs_temperature_floor_K=pyrolysis_gibbs_temperature_floor_K,
        pyrolysis_ch4_to_h2_co_split_frac=pyrolysis_ch4_to_h2_co_split_frac,
        pyrolysis_co_reduce_h2_inject_frac=pyrolysis_co_reduce_h2_inject_frac,
        fuel_type=fuel_type,
        pyrolysis_tar_carbon_frac=pyrolysis_tar_carbon_frac,
        pyrolysis_koks_carbon_frac=pyrolysis_koks_carbon_frac,
    )
    prod_curr = allocate_pyrolysis_products_hamel(T_K=T_K, P_Pa=P_Pa, **curr_atoms, **_gibbs_kw)
    prod_prev = allocate_pyrolysis_products_hamel(T_K=T_K, P_Pa=P_Pa, **prev_atoms, **_gibbs_kw)
    return {sp: max(float(prod_curr.get(sp, 0.0)) - float(prod_prev.get(sp, 0.0)), 0.0) for sp in GAS_SPECIES}


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
    m_vm_in_classes: np.ndarray | None = None,
    m_moist_in_classes: np.ndarray | None = None,
    char_index: int | None = None,
    vm_index: int,
    moisture_index: int,
    vm_release_fraction_cap: float | None = None,
    pyrolysis_ch4_to_h2_co_split_frac: float = 0.0,
    pyrolysis_co_reduce_h2_inject_frac: float = 0.0,
    pyrolysis_tar_yield_mode: str = "hegermann",
) -> PyrolysisSourceBundle:
    """计算干燥/热解气相源项与固相扣减。

    Ref: Hamel (1999) Eq. 4.4, 4.9, 4.10-4.12；焦油碳分率 Eq.4.21–4.22（hegermann）

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
    if vm_release_fraction_cap is not None:
        x_vm = min(x_vm, float(np.clip(vm_release_fraction_cap, 0.0, 1.0)))
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

    # Hamel §4.2 / §5.2.6：Teer H/C ≈ 燃料摩尔 H/C（论文原子量 C=12、H=1）
    target_tar_hc_ratio = None
    if C_dry > 1e-12 and H_dry >= 0.0:
        target_tar_hc_ratio = fuel_tar_target_hc_ratio(C_dry, H_dry)

    # coal + hegermann：Eq.4.21–4.22 绝对碳分率 x_C,Koks / x_C,Teer（相对燃料 waf 碳）
    from src.core.pyrolysis_carbon_split import resolve_pyrolysis_carbon_fractions

    f_koks_eff, f_tar_eff = resolve_pyrolysis_carbon_fractions(
        fuel_type=fuel_type,
        yield_mode=pyrolysis_tar_yield_mode,
        P_Pa=float(P),
        T_K=float(T),
        pyrolysis_tar_carbon_frac=float(pyrolysis_tar_carbon_frac),
    )

    prod = allocate_pyrolysis_products_hamel(
        T_K=float(T),
        P_Pa=float(P),
        nC=vm_atoms["C"],
        nH=vm_atoms["H"],
        nO=vm_atoms["O"],
        nN=vm_atoms["N"],
        nS=vm_atoms["S"],
        fuel_type=fuel_type,
        pyrolysis_tar_carbon_frac=f_tar_eff,
        pyrolysis_koks_carbon_frac=f_koks_eff,
        target_tar_hc_ratio=target_tar_hc_ratio,
        pyrolysis_ch4_to_h2_co_split_frac=float(pyrolysis_ch4_to_h2_co_split_frac),
        pyrolysis_co_reduce_h2_inject_frac=float(pyrolysis_co_reduce_h2_inject_frac),
    )
    for sp, value in prod.items():
        gas_source[GAS_SPECIES_INDEX[sp]] += float(value)

    solid_sink = np.zeros(solid_shape, dtype=np.float64)
    nk = int(solid_shape[0])
    uniform_class_frac = np.full(nk, 1.0 / max(nk, 1), dtype=np.float64)
    if m_vm_in_classes is None:
        vm_class_frac = uniform_class_frac
    else:
        vm_classes = np.asarray(m_vm_in_classes, dtype=np.float64)
        assert vm_classes.shape == (nk,), f"m_vm_in_classes shape mismatch: expected {(nk,)}, got {vm_classes.shape}"
        vm_total = float(np.sum(np.maximum(vm_classes, 0.0)))
        vm_class_frac = np.maximum(vm_classes, 0.0) / max(vm_total, 1e-12) if vm_total > 1e-12 else uniform_class_frac
    if m_moist_in_classes is None:
        moist_class_frac = uniform_class_frac
    else:
        moist_classes = np.asarray(m_moist_in_classes, dtype=np.float64)
        assert moist_classes.shape == (nk,), f"m_moist_in_classes shape mismatch: expected {(nk,)}, got {moist_classes.shape}"
        moist_total = float(np.sum(np.maximum(moist_classes, 0.0)))
        moist_class_frac = np.maximum(moist_classes, 0.0) / max(moist_total, 1e-12) if moist_total > 1e-12 else uniform_class_frac
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
        solid_sink[:, char_index] += carbon_to_char * ATOMIC_MASS_KG_PER_MOL["C"] * vm_class_frac
    solid_sink[:, moisture_index] -= float(m_moist_in) * float(dry["X_dry"][-1]) * moist_class_frac
    solid_sink[:, vm_index] -= m_vm_rel * vm_class_frac

    return PyrolysisSourceBundle(
        gas_source=gas_source,
        solid_sink=solid_sink,
        x_dry=float(np.clip(dry["X_dry"][-1], 0.0, 1.0)),
        x_vm=x_vm,
    )
