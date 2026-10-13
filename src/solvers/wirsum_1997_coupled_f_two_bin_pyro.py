"""Q6 联立 F 接入 C 两级分与 D 热解源。不经 ``Reactor``。

39 未知 / 39 方程：在等温单床 + FB 气接口 + 两格旋风上，
床炭改为 ``(fine, coarse)``，热解气/焦源作为冻结 q，不重写出口 N。

粗粒队列消耗与迁出用 ``two_bin_cohort_representative``；细粒就地 R1，
无更小一级。旋风仍对总质量用冻结 K（原文无分级分离律）。
新鲜焦全部进粗级分。焦油、灰库存、能量、R2–R12 仍 ``omitted, named``。

Source: Hamel (1999) Eq. 2.4–2.6, 5.8–5.23, §2.1, §5.1.1.2;
docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q6–Q7, Q4
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.elemental_ledger import H2O_MOLAR_MASS_KG_PER_MOL
from src.core.solid_transport_reduced import (
    MIGRATION_MODEL_COHORT,
    cohort_split_rates,
    radius_retreat_m_per_s,
)
from src.core.species import MW_C_KG_PER_MOL, gas_diffusivity_correlation
from src.kinetics.char_reactions import rate_R1
from src.physics.wirsum_1997_freeboard import FREEBOARD_MODEL_HAMEL_EXTENDED
from src.solvers.wirsum_1997_bed_freeboard import (
    BedFreeboardGeometry,
    eject_bed_top,
    particle_fate_for_mean_launch,
)
from src.solvers.wirsum_1997_coupled_f_minimal import (
    _elem_C,
    _elem_H,
    _elem_O,
    r1_gas_source,
)
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK, cyclone_k_implied
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    CARBON_BASIS_FUEL_WAF,
    REMAINDER_DESTINATION,
    TAR_CONC_DILUTE_SHARED_SLOT,
    FuelAtomFeed,
    dry_mass_from_cho_atoms,
    map_one_state,
)
from src.solvers.wirsum_1997_stiffness_full_gas import (
    IDX,
    MAJOR,
    N_MAJ,
    char_area_m2,
    fd_jacobian,
)
from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry, classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

N_X: int = 39
N_BIN: int = 2
IDX_FINE: int = 0
IDX_COARSE: int = 1
OMITTED: tuple[str, ...] = (
    "R2_R12",
    "tar_transport",
    "ash_inventory",
    "freeboard_reaction",
    "energy_wall_heat",
    "return_leg_holdup",
    "cyclone_K_update_P",
    "full_outer_P",
    "production_f_mig_3",
    "vm_empirical_K",
)
OCCUPANCY_SOURCE: str = "named_single_bed_plug_flow_age"


@dataclass(frozen=True)
class CoupledTwoBinPyroQ:
    """冻结 q。热解源与 K_eject 不乘旧库存当流量。级分顺序 ``(fine, coarse)``。"""

    geo_ex: ExchangeGeometry
    geo_fb: BedFreeboardGeometry
    N_in_d: npt.NDArray[np.float64]
    N_in_b: npt.NDArray[np.float64]
    n_py_d: npt.NDArray[np.float64]  # [mol/s] MAJOR，进密相
    n_moist: npt.NDArray[np.float64]  # [mol/s] 干燥水，进密相
    nC_tar_omitted: float  # [mol-atom/s]
    m_feed: npt.NDArray[np.float64]  # [kg/s]
    K_eject: npt.NDArray[np.float64]  # [1/s]
    w_leave: npt.NDArray[np.float64]  # [-]
    w_return: npt.NDArray[np.float64]  # [-]
    d_p: npt.NDArray[np.float64]  # [m]
    k: CycloneCellK
    D_g: float  # [m²/s]
    H_bed: float  # [m]
    rho_p: float  # [kg/m³]
    occupancy_source: str
    vm_pool_basis: str
    nC_fuel: float  # [mol-atom/s]
    nC_fixed: float  # [mol-atom/s] 工业分析元数据，不是焦炭给料
    rho_g_fb_ref: float  # [kg/m³] 冻结参考气密；outer_P 只从此值缩放
    nC_koks: float = 0.0  # [mol-atom/s] Hegermann 焦炭碳
    nC_unreleased: float = 0.0  # [mol-atom/s]
    nH_unreleased: float = 0.0  # [mol-atom/s]
    nO_unreleased: float = 0.0  # [mol-atom/s]
    nN_unreleased: float = 0.0  # [mol-atom/s]
    nS_unreleased: float = 0.0  # [mol-atom/s]
    m_unreleased_kg_s: float = 0.0  # [kg/s] 未转化 VM 固体出口
    carbon_basis: str = CARBON_BASIS_FUEL_WAF
    remainder_destination: str = REMAINDER_DESTINATION
    tar_conc_model: str = TAR_CONC_DILUTE_SHARED_SLOT
    hegermann_T_clipped: bool = False
    hegermann_P_clipped: bool = False
    n_NH3_bypass_mol_s: float = 0.0  # [mol/s] 无反应旁路，不进七组分
    n_H2S_bypass_mol_s: float = 0.0  # [mol/s] 无反应旁路，不进七组分
    minor_bypass_in_major_density: bool = False


def _k_eject_weights(geo: BedFreeboardGeometry, H_bed: float) -> tuple[float, float, float, str]:
    """``K_eject = m_eject(M=1)/M_geom``。未判定既不是离开也不是回落。"""
    m_eject = eject_bed_top(
        geo=geo, m_char=1.0, m_ash=0.0, model=FREEBOARD_MODEL_HAMEL_EXTENDED
    )
    M_geom = (
        float(geo.rho_p)
        * max(1.0 - float(geo.eps_d), 0.0)
        * max(1.0 - float(geo.eps_b), 0.0)
        * float(geo.area)
        * float(H_bed)
    )
    K = float(m_eject[0]) / max(M_geom, 1.0e-12)
    fate = particle_fate_for_mean_launch(geo=geo, model=FREEBOARD_MODEL_HAMEL_EXTENDED)
    if fate == "returns":
        return K, 0.0, 1.0, fate
    if fate == "reaches_next_cell":
        return K, 1.0, 0.0, fate
    return K, 0.0, 0.0, fate


def _pyro_major_and_tar_c(cell: dict[str, Any]) -> tuple[npt.NDArray[np.float64], float]:
    n_py = cell["n_py_mol_s"]
    src = np.zeros(N_MAJ, dtype=np.float64)
    for name, i in IDX.items():
        src[i] = float(n_py.get(name, 0.0))
    nC_major = _elem_C(src)
    nC_tar = float(cell["atoms"]["C"]) - nC_major
    return src, nC_tar


def moist_inlet_mol_s(m_moist_kg_s: float) -> npt.NDArray[np.float64]:
    """干燥水进密相 [mol/s]。只由给料水分决定，不沿用上一次冻结的数组。"""
    row = np.zeros(N_MAJ, dtype=np.float64)
    row[IDX["H2O"]] = float(m_moist_kg_s) / H2O_MOLAR_MASS_KG_PER_MOL
    return row


def apply_pyrolysis_map(
    q: CoupledTwoBinPyroQ,
    mapped: dict[str, Any],
) -> CoupledTwoBinPyroQ:
    """把 Hegermann/DAEM 映射写入 q：焦炭给料、热解源、未释放出口。不重写 N。"""
    n_py, nC_tar = _pyro_major_and_tar_c(mapped["cells"][0])
    nC_koks = float(mapped["nC_koks"])
    heg = mapped["hegermann"]
    return replace(
        q,
        n_py_d=n_py,
        nC_tar_omitted=float(nC_tar),
        vm_pool_basis=str(mapped["vm_pool_basis"]),
        m_feed=np.array([0.0, nC_koks * MW_C_KG_PER_MOL], dtype=np.float64),
        nC_koks=nC_koks,
        nC_unreleased=float(mapped["nC_unreleased_vm"]),
        nH_unreleased=float(mapped["nH_unreleased"]),
        nO_unreleased=float(mapped["nO_unreleased"]),
        nN_unreleased=float(mapped["nN_unreleased"]),
        nS_unreleased=float(mapped["nS_unreleased"]),
        m_unreleased_kg_s=float(mapped["m_unreleased_kg_s"]),
        carbon_basis=str(mapped["carbon_basis"]),
        remainder_destination=str(mapped["remainder_destination"]),
        tar_conc_model=TAR_CONC_DILUTE_SHARED_SLOT,
        hegermann_T_clipped=bool(heg["T_clipped"]),
        hegermann_P_clipped=bool(heg["P_clipped"]),
        n_NH3_bypass_mol_s=float(mapped["cells"][0]["n_py_mol_s"].get("NH3", 0.0)),
        n_H2S_bypass_mol_s=float(mapped["cells"][0]["n_py_mol_s"].get("H2S", 0.0)),
        minor_bypass_in_major_density=False,
    )


def default_char_feed(*, m_char_kg_s: float = 0.02) -> FuelAtomFeed:
    """示意原料，使固定碳 ≈ ``m_char``。干料质量由 C/H/O/灰还原，不是实验夹带。"""
    nC_fixed = float(m_char_kg_s) / MW_C_KG_PER_MOL
    nC_vm = nC_fixed * (5.342 / 4.658)
    nC_fuel = nC_fixed + nC_vm
    feed = FuelAtomFeed(
        nC_fuel=nC_fuel,
        nH_fuel=1.2 * nC_fuel,
        nO_fuel=0.4 * nC_fuel,
        nN_fuel=0.0,
        nS_fuel=0.0,
        nC_fixed=nC_fixed,
        nC_vm=nC_vm,
        m_moist_kg_s=0.002,
        m_dry_kg_s=0.0,
        h_f_dry_J_kg=-3.96e6,
        ash_dry_wt=11.41,
        VM_daf=53.42,
    )
    return replace(feed, m_dry_kg_s=dry_mass_from_cho_atoms(feed))


def freeze_q(
    *,
    geo_ex: ExchangeGeometry | None = None,
    geo_fb: BedFreeboardGeometry | None = None,
    N_in_d: npt.NDArray[np.float64] | None = None,
    N_in_b: npt.NDArray[np.float64] | None = None,
    m_char_kg_s: float = 0.002,
    H_bed: float = 1.0,
    k: CycloneCellK | None = None,
    d_fine: float = 2.0e-4,
    d_coarse: float = 3.0e-4,
    feed: FuelAtomFeed | None = None,
) -> CoupledTwoBinPyroQ:
    """水力/轨迹冻结核级 ``K_eject``；热解映射一次写入 q。

    默认 ``d_fine`` 离床、``d_coarse`` 回落（``hamel_extended`` 命运）。
    ``λ=(d_fine/d_coarse)³`` 须使队列氧化量不超过入口 O2；1 mm / 0.2 mm
    的 λ≈0.008 会把几乎全部粗粒碳烧掉，本核不存在正根。
    """
    if not (d_fine < d_coarse):
        raise ValueError("d_fine must be < d_coarse [m]")
    ex = geo_ex or ExchangeGeometry()
    fb0 = geo_fb or BedFreeboardGeometry()
    if N_in_d is None:
        N_in_d = np.array([0.05, 0.02, 0.04, 0.03, 0.01, 0.0, 2.0], dtype=np.float64)
    if N_in_b is None:
        N_in_b = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.40, 1.5], dtype=np.float64)
    ck = k or CycloneCellK()
    fuel = feed or default_char_feed(m_char_kg_s=m_char_kg_s)
    mapped = map_one_state(
        fuel,
        P_Pa=float(ex.P),
        T_K=float(ex.T),
        X_vm=np.array([0.0, 1.0], dtype=np.float64),
        tau_s=np.array([1.0], dtype=np.float64),
        occupancy_source=OCCUPANCY_SOURCE,
        fuel_type="coal",
    )
    n_py_d, nC_tar = _pyro_major_and_tar_c(mapped["cells"][0])
    n_moist = moist_inlet_mol_s(float(fuel.m_moist_kg_s))
    nC_koks = float(mapped["nC_koks"])
    m_feed = np.array([0.0, nC_koks * MW_C_KG_PER_MOL], dtype=np.float64)
    d_p = np.array([float(d_fine), float(d_coarse)], dtype=np.float64)
    K = np.zeros(N_BIN, dtype=np.float64)
    wL = np.zeros(N_BIN, dtype=np.float64)
    wR = np.zeros(N_BIN, dtype=np.float64)
    for i, d in enumerate(d_p):
        geo_i = replace(fb0, d_p=float(d))
        K[i], wL[i], wR[i], _fate = _k_eject_weights(geo_i, H_bed)
    D_g = gas_diffusivity_correlation(float(ex.T), float(ex.P))
    heg = mapped["hegermann"]
    return CoupledTwoBinPyroQ(
        geo_ex=ex,
        geo_fb=replace(fb0, d_p=float(d_fine)),
        N_in_d=np.asarray(N_in_d, dtype=np.float64).copy(),
        N_in_b=np.asarray(N_in_b, dtype=np.float64).copy(),
        n_py_d=n_py_d,
        n_moist=n_moist,
        nC_tar_omitted=float(nC_tar),
        m_feed=m_feed,
        K_eject=K,
        w_leave=wL,
        w_return=wR,
        d_p=d_p,
        k=ck,
        D_g=float(D_g),
        H_bed=float(H_bed),
        rho_p=float(fb0.rho_p),
        occupancy_source=OCCUPANCY_SOURCE,
        vm_pool_basis=str(mapped["vm_pool_basis"]),
        nC_fuel=float(fuel.nC_fuel),
        nC_fixed=float(fuel.nC_fixed),
        rho_g_fb_ref=float(fb0.rho_g),
        nC_koks=nC_koks,
        nC_unreleased=float(mapped["nC_unreleased_vm"]),
        nH_unreleased=float(mapped["nH_unreleased"]),
        nO_unreleased=float(mapped["nO_unreleased"]),
        nN_unreleased=float(mapped["nN_unreleased"]),
        nS_unreleased=float(mapped["nS_unreleased"]),
        m_unreleased_kg_s=float(mapped["m_unreleased_kg_s"]),
        carbon_basis=str(mapped["carbon_basis"]),
        remainder_destination=str(mapped["remainder_destination"]),
        tar_conc_model=TAR_CONC_DILUTE_SHARED_SLOT,
        hegermann_T_clipped=bool(heg["T_clipped"]),
        hegermann_P_clipped=bool(heg["P_clipped"]),
        n_NH3_bypass_mol_s=float(mapped["cells"][0]["n_py_mol_s"].get("NH3", 0.0)),
        n_H2S_bypass_mol_s=float(mapped["cells"][0]["n_py_mol_s"].get("H2S", 0.0)),
        minor_bypass_in_major_density=False,
    )


def unpack(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_X,), f"expected {N_X}, got {z.shape}"
    return {
        "N_d": z[0:7].copy(),
        "N_b": z[7:14].copy(),
        "M_fine": float(z[14]),
        "M_coarse": float(z[15]),
        "N_fb": z[16:23].copy(),
        "N_c1": z[23:30].copy(),
        "N_c2": z[30:37].copy(),
        "M_c1": float(z[37]),
        "M_c2": float(z[38]),
    }


def pack(
    *,
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    M_fine: float,
    M_coarse: float,
    N_fb: npt.NDArray[np.float64],
    N_c1: npt.NDArray[np.float64],
    N_c2: npt.NDArray[np.float64],
    M_c1: float,
    M_c2: float,
) -> npt.NDArray[np.float64]:
    return np.concatenate(
        [
            np.asarray(N_d, dtype=np.float64),
            np.asarray(N_b, dtype=np.float64),
            np.array([M_fine, M_coarse], dtype=np.float64),
            np.asarray(N_fb, dtype=np.float64),
            np.asarray(N_c1, dtype=np.float64),
            np.asarray(N_c2, dtype=np.float64),
            np.array([M_c1, M_c2], dtype=np.float64),
        ]
    )


def _M_bins(u: dict[str, Any]) -> npt.NDArray[np.float64]:
    return np.array([u["M_fine"], u["M_coarse"]], dtype=np.float64)


def edges_from_x(x: npt.NDArray[np.float64], q: CoupledTwoBinPyroQ) -> dict[str, Any]:
    """当前 x 与冻结系数唯一决定边流 [kg/s]。"""
    u = unpack(x)
    M = _M_bins(u)
    m_eject = q.K_eject * M
    m_leave = q.w_leave * m_eject
    m_return = q.w_return * m_eject
    m_leave_tot = float(np.sum(m_leave))
    m_auf0 = float(q.k.K_auf) * float(u["M_c1"])
    m_ab0 = float(q.k.K_ab) * float(u["M_c1"])
    m_auf1 = float(q.k.K_auf) * float(u["M_c2"])
    m_ab1 = float(q.k.K_ab) * float(u["M_c2"])
    m_rez = m_ab0
    if m_leave_tot > 1.0e-30:
        m_rez_bin = m_rez * m_leave / m_leave_tot
    else:
        m_rez_bin = np.zeros(N_BIN, dtype=np.float64)
    return {
        "m_eject": m_eject,
        "m_leave": m_leave,
        "m_return": m_return,
        "m_rez_bin": m_rez_bin,
        "m_auf0": m_auf0,
        "m_ab0": m_ab0,
        "m_auf1": m_auf1,
        "m_ab1": m_ab1,
        "m_rez": m_rez,
        "m_fly": m_auf1,
        "m_leave_tot": m_leave_tot,
    }


def extents_from_x(x: npt.NDArray[np.float64], q: CoupledTwoBinPyroQ) -> dict[str, Any]:
    """每级一次 R1；粗粒队列消耗与迁出分开。气固共用同一消耗。"""
    u = unpack(x)
    N_d = u["N_d"]
    ntot = float(np.sum(N_d))
    C_O2 = 0.0
    if abs(ntot) > 1.0e-30:
        C_O2 = (float(N_d[IDX["O2"]]) / ntot) * float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    M = _M_bins(u)
    e1_kin = np.zeros(N_BIN, dtype=np.float64)
    alpha = np.zeros(N_BIN, dtype=np.float64)
    area = np.zeros(N_BIN, dtype=np.float64)
    for i in range(N_BIN):
        r_surf, a = rate_R1(float(q.geo_ex.T), C_O2, float(q.d_p[i]), float(q.D_g))
        area[i] = char_area_m2(float(M[i]), rho_s=q.rho_p, d_p=float(q.d_p[i]))
        e1_kin[i] = float(r_surf) * float(area[i])
        alpha[i] = float(a)
    R_kin = MW_C_KG_PER_MOL * e1_kin
    u_c = radius_retreat_m_per_s(
        float(R_kin[IDX_COARSE]),
        M_kg=float(M[IDX_COARSE]),
        rho=q.rho_p,
        d_p=float(q.d_p[IDX_COARSE]),
    )
    split = cohort_split_rates(
        float(M[IDX_COARSE]),
        d_from=float(q.d_p[IDX_COARSE]),
        d_to=float(q.d_p[IDX_FINE]),
        u_radius=u_c,
    )
    consume = np.array(
        [float(R_kin[IDX_FINE]), float(split["m_consume_kg_s"])],
        dtype=np.float64,
    )
    e1 = consume / MW_C_KG_PER_MOL
    return {
        "e1": e1,
        "e1_kinetic": e1_kin,
        "alpha": alpha,
        "area": area,
        "C_O2": float(C_O2),
        "consume_kg_s": consume,
        "m_migrate_kg_s": float(split["m_migrate_kg_s"]),
        "u_coarse_m_per_s": float(u_c),
        "split": split,
        "model": MIGRATION_MODEL_COHORT,
    }


def coupled_residual(
    x: npt.NDArray[np.float64],
    q: CoupledTwoBinPyroQ,
    *,
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    """一张 F。边流、两级 R1、迁出各算一次。热解只从 q 进入源项。"""
    u = unpack(x)
    ed = edges_from_x(x, q)
    ext = extents_from_x(x, q)
    consume = float(rate_scale) * ext["consume_kg_s"]
    e1 = consume / MW_C_KG_PER_MOL
    src = r1_gas_source({"e1": float(e1[IDX_FINE]), "alpha": float(ext["alpha"][IDX_FINE])})
    src = src + r1_gas_source(
        {"e1": float(e1[IDX_COARSE]), "alpha": float(ext["alpha"][IDX_COARSE])}
    )
    C_d = (u["N_d"] / max(float(np.sum(u["N_d"])), 1.0e-30)) * (
        float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    )
    C_b = (u["N_b"] / max(float(np.sum(u["N_b"])), 1.0e-30)) * (
        float(q.geo_ex.P) / (Rg * float(q.geo_ex.T))
    )
    N_ex = float(q.geo_ex.K_bd) * float(q.geo_ex.V_b) * (C_b - C_d)
    F_d = q.N_in_d + q.n_py_d + q.n_moist - u["N_d"] + N_ex + src
    F_b = q.N_in_b - u["N_b"] - N_ex
    mig = float(rate_scale) * float(ext["m_migrate_kg_s"])
    F_Mf = (
        float(q.m_feed[IDX_FINE])
        + float(ed["m_return"][IDX_FINE])
        + float(ed["m_rez_bin"][IDX_FINE])
        + mig
        - float(ed["m_eject"][IDX_FINE])
        - float(consume[IDX_FINE])
    )
    F_Mc = (
        float(q.m_feed[IDX_COARSE])
        + float(ed["m_return"][IDX_COARSE])
        + float(ed["m_rez_bin"][IDX_COARSE])
        - mig
        - float(ed["m_eject"][IDX_COARSE])
        - float(consume[IDX_COARSE])
    )
    F_fb = u["N_fb"] - u["N_d"] - u["N_b"]
    F_c1g = u["N_fb"] - u["N_c1"]
    F_c2g = u["N_c1"] - u["N_c2"]
    F_Mc1 = ed["m_leave_tot"] + ed["m_ab1"] - ed["m_auf0"] - ed["m_ab0"]
    F_Mc2 = ed["m_auf0"] - ed["m_auf1"] - ed["m_ab1"]
    return np.concatenate(
        [
            F_d,
            F_b,
            np.array([F_Mf, F_Mc], dtype=np.float64),
            F_fb,
            F_c1g,
            F_c2g,
            np.array([F_Mc1, F_Mc2], dtype=np.float64),
        ]
    )


def atom_gap_from_F(x: npt.NDArray[np.float64], q: CoupledTwoBinPyroQ) -> dict[str, float]:
    """任意 x：残差结构应闭合 F 内出现的进出口。焦油碳在 q 汇，不进本恒等式。"""
    F = coupled_residual(x, q)
    u = unpack(x)
    ed = edges_from_x(x, q)
    ext = extents_from_x(x, q)
    gas_F = F[0:7] + F[7:14] - F[16:23] + F[23:30] + F[30:37]
    solid_F = float(F[14]) + float(F[15]) + float(F[37]) + float(F[38])
    N_in = q.N_in_d + q.N_in_b + q.n_py_d + q.n_moist
    m_feed = float(np.sum(q.m_feed))
    consume = float(np.sum(ext["consume_kg_s"]))
    return {
        "C": abs(
            _elem_C(gas_F)
            + solid_F / MW_C_KG_PER_MOL
            + _elem_C(u["N_c2"])
            + ed["m_fly"] / MW_C_KG_PER_MOL
            - _elem_C(N_in)
            - m_feed / MW_C_KG_PER_MOL
        ),
        "H": abs(_elem_H(gas_F) + _elem_H(u["N_c2"]) - _elem_H(N_in)),
        "O": abs(_elem_O(gas_F) + _elem_O(u["N_c2"]) - _elem_O(N_in)),
        "N": abs(float(gas_F[IDX["N2"]]) + float(u["N_c2"][IDX["N2"]]) - float(N_in[IDX["N2"]])),
        "mass_kg_s": abs(solid_F + ed["m_fly"] + consume - m_feed),
        "nC_tar_omitted": float(q.nC_tar_omitted),
    }


def default_x0(q: CoupledTwoBinPyroQ) -> npt.NDArray[np.float64]:
    N_d = q.N_in_d + q.n_py_d + q.n_moist
    N_b = q.N_in_b.copy()
    # 粗粒回床，R1 只在密相。C_O2=0 时 M_coarse 列为零；过量密相 O2 又使全步越界。
    if float(N_d[IDX["O2"]]) <= 1.0e-12:
        N_d[IDX["O2"]] = 0.05
    N_tot = N_d + N_b
    M_fine = 0.03
    M_coarse = 0.03
    m_leave = float(q.K_eject[IDX_FINE] * q.w_leave[IDX_FINE] * M_fine)
    ksum = q.k.K_auf + q.k.K_ab
    Mc = max(m_leave / max(ksum, 1.0e-12), 0.01)
    return pack(
        N_d=N_d,
        N_b=N_b,
        M_fine=M_fine,
        M_coarse=M_coarse,
        N_fb=N_tot.copy(),
        N_c1=N_tot.copy(),
        N_c2=N_tot.copy(),
        M_c1=Mc,
        M_c2=Mc,
    )


def omitted_named() -> dict[str, str]:
    return {name: "omitted, named" for name in OMITTED}


def q_diagnostics(q: CoupledTwoBinPyroQ) -> dict[str, Any]:
    fates = []
    for d in q.d_p:
        geo_i = replace(q.geo_fb, d_p=float(d))
        fates.append(particle_fate_for_mean_launch(geo=geo_i, model=FREEBOARD_MODEL_HAMEL_EXTENDED))
    return {
        "K_eject_1_per_s": q.K_eject.tolist(),
        "K_eject_source": "frozen_linearization_m_eject_over_M_geom",
        "w_leave": q.w_leave.tolist(),
        "w_return": q.w_return.tolist(),
        "d_p_m": q.d_p.tolist(),
        "fate": fates,
        "fresh_char_to_coarse": True,
        "pyrolysis_to_dense_named": True,
        "occupancy_source": q.occupancy_source,
        "vm_pool_basis": q.vm_pool_basis,
        "expected_vm_pool_basis": CARBON_BASIS_FUEL_WAF,
        "carbon_basis": q.carbon_basis,
        "nC_koks": float(q.nC_koks),
        "nC_unreleased": float(q.nC_unreleased),
        "remainder_destination": q.remainder_destination,
        "tar_conc_model": q.tar_conc_model,
        "lambda_size": float((q.d_p[0] / q.d_p[1]) ** 3),
        "nC_tar_omitted": q.nC_tar_omitted,
        "engineering_K_fb_leave": "not_used",
        "cyclone": cyclone_k_implied(q.k),
        "units": {"gas": "mol/s", "solid": "kg/s", "K": "1/s"},
        "n_unknowns": N_X,
        "omitted": omitted_named(),
        "species": list(MAJOR),
        "hides_H_over_umf": False,
        "rewrites_unknown_N": False,
    }


def solve_coupled(
    q: CoupledTwoBinPyroQ,
    *,
    x0: npt.NDArray[np.float64] | None = None,
    rate_scale: float = 1.0,
    max_iter: int = 40,
) -> dict[str, Any]:
    """同一张 F 上的 Wirsum strict Newton。不是 Q6 根与两级分核拼接。"""
    x_start = default_x0(q) if x0 is None else np.asarray(x0, dtype=np.float64)
    scale = np.maximum(np.abs(x_start), 0.05)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        # 越界试探不评价 F、不改写 x。不是 residual 内 clip。
        if float(np.min(z)) < 0.0:
            return np.full(N_X, 1.0e6, dtype=np.float64)
        return coupled_residual(z, q, rate_scale=rate_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return fd_jacobian(
            lambda w: coupled_residual(w, q, rate_scale=rate_scale),
            z,
            scale=scale,
        )

    sol = wirsum_1997_strict_newton(res, jac, x_start, max_iter=max_iter)
    stop = classify_newton_stop(sol, model_is_full_rate=abs(rate_scale - 1.0) < 1.0e-15)
    F2 = coupled_residual(sol.x, q, rate_scale=rate_scale)
    return {
        **stop,
        "x": np.asarray(sol.x, dtype=np.float64),
        "F": np.asarray(sol.F, dtype=np.float64),
        "rate_scale": float(rate_scale),
        "q_purity_gap": float(np.max(np.abs(F2 - sol.F))),
        "rank_last": int(sol.rank_last),
        "cond_last": float(sol.cond_last),
        "edges": edges_from_x(sol.x, q),
        "extent": extents_from_x(sol.x, q),
        "atom_gap": atom_gap_from_F(sol.x, q),
    }


def recycle_changes_fine_root(
    *,
    k_lo: CycloneCellK,
    k_hi: CycloneCellK,
) -> dict[str, Any]:
    """细粒进旋风；只改 ``K_ab``，细库存根必须变。"""
    a = solve_coupled(freeze_q(k=k_lo))
    b = solve_coupled(freeze_q(k=k_hi))
    return {
        "M_fine_lo_kg": float(a["x"][14]),
        "M_fine_hi_kg": float(b["x"][14]),
        "gap_kg": abs(float(a["x"][14]) - float(b["x"][14])),
        "both_physical": bool(a["physical_convergence"] and b["physical_convergence"]),
    }


def build_stage_cd_coupled_payload() -> dict[str, Any]:
    q = freeze_q()
    out = solve_coupled(q)
    return {
        "source": "Wirsum alignment C+D two-bin pyrolysis coupled F",
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "n_unknowns": N_X,
        "physical_convergence": out["physical_convergence"],
        "residual_l2": out["residual_l2"],
        "q_purity_gap": out["q_purity_gap"],
        "rank_last": out["rank_last"],
        "cond_last": out["cond_last"],
        "atom_gap": out["atom_gap"],
        "q": q_diagnostics(q),
        "omitted": omitted_named(),
        "note": (
            "一张 39×39 F。热解在 q，出口 N 不改写。"
            "生产 Reactor 未改。省略项见 omitted。"
        ),
    }
