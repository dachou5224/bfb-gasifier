"""具名 F 谱系：在 F_E0 上一次取消一项省略。不经 ``Reactor``。

各阶段叠加上一阶段已取消项。旋风 K 与生产 ``hamel_outer_refresh`` 保持冻结。

Source: Hamel (1999) Eq. 2.4–2.8, 5.8–5.41;
docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import MW_C_KG_PER_MOL
from src.kinetics.tar_reactions import (
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from src.solvers.wirsum_1997_coupled_f_energy import (
    CoupledEnergyQ,
    _H_major,
    _H_tar,
    default_x0 as energy_x0,
    freeze_energy_q,
    outer_P as outer_P_e0,
    pack as pack_e0,
    q_at_T,
)
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
    IDX_COARSE,
    IDX_FINE,
    N_X as N_SPECIES,
    _M_bins,
    apply_pyrolysis_map,
    moist_inlet_mol_s,
    coupled_residual as species_residual,
    edges_from_x,
    unpack as unpack_species,
)
from src.solvers.wirsum_1997_named_f_p2 import (
    NamedHydroSnapshot,
    NamedReferenceInputs,
    capture_named_reference,
    finite_hydro_from_x,
    k_eject_from_geo,
)
from src.solvers.wirsum_1997_named_f_registry import (
    cancelled_this_stage,
    f_name,
    n_unknowns,
    omitted_named,
)
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    TAR_CONC_DILUTE_SHARED_SLOT,
    TAR_CONC_STRICT_EOS,
    map_one_state,
    unreleased_vm_enthalpy_W,
)
from src.solvers.wirsum_1997_stiffness_full_gas import (
    HOMOG_RXN_NAMES,
    IDX,
    N_MAJ,
    _y_full,
    char_hetero_source,
    energy_jacobian_quality_at,
    fd_jacobian,
    gas_homogeneous_source,
    default_full_gas_inlets,
)
from src.solvers.wirsum_1997_stiffness_ladder import classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton
from src.thermal.devolatilization import daem_conversion

N_SPECIES_BLOCK: int = N_SPECIES
RHO_ASH_Q_KG_M3: float = 75.0  # [kg/m³] E2 以前写在 q，不是库存反演
DAEM_TAU_S: float = 1.0  # [s] 具名等温停留，不是 H/u_mf
GEO_FB_RHO_G_T_REF_K: float = 1100.0  # [K] 冻结 geo_fb.rho_g 的参考温度；理想气体缩放，不是拟合炉温
R10_FUEL_TYPE: str = "coal"
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]
# F_E1/F_E2 捆包仍保留；探针按反应逐个接入，不改具名阶段编号。
NAMED_PROBE_RXNS: tuple[str, ...] = ("R2", "R3", "R4", "R5", "R6", "R7", "R8", "R12")


def default_enabled_rxns(stage: int) -> frozenset[str]:
    """F_E1 = R2–R4；F_E2 起再开均相。R1 始终在 F_E0 物种块里。"""
    on: set[str] = set()
    if int(stage) >= 1:
        on.update(("R2", "R3", "R4"))
    if int(stage) >= 2:
        on.update(HOMOG_RXN_NAMES)
    return frozenset(on)


def stage_for_enabled(enabled: frozenset[str]) -> int:
    """探针布局：有均相用 40 维 F_E2 槽，否则 F_E1；空集退回 F_E0。"""
    if any(name in enabled for name in HOMOG_RXN_NAMES):
        return 2
    if any(name in enabled for name in ("R2", "R3", "R4")):
        return 1
    return 0


@dataclass(frozen=True)
class NamedStageQ:
    """F_E0 的 q 加上本阶段开关。内层只读。``reference`` 是不可变 b。"""

    energy: CoupledEnergyQ
    stage: int
    Q_W: float = 0.0  # [W] E5 起为显式 q，禁止反推 1100 K
    rho_ash: float = RHO_ASH_Q_KG_M3  # [kg/m³]
    tau_daem_s: float = DAEM_TAU_S  # [s]
    enabled_rxns: frozenset[str] = frozenset()
    reference: NamedReferenceInputs | None = None
    hydro: NamedHydroSnapshot | None = None


def freeze_stage_q(
    stage: int,
    *,
    Q_W: float = 0.0,
    enabled_rxns: frozenset[str] | None = None,
    **kwargs: Any,
) -> NamedStageQ:
    qe = freeze_energy_q(**kwargs)
    on = default_enabled_rxns(stage) if enabled_rxns is None else frozenset(enabled_rxns)
    nq = NamedStageQ(energy=qe, stage=int(stage), Q_W=float(Q_W), enabled_rxns=on)
    return replace(nq, reference=capture_named_reference(nq))


def overlay_eject_from_p(nq: NamedStageQ, nq_p: NamedStageQ) -> NamedStageQ:
    """把 ``P(x)`` 的抛射 K/命运/气密写进另一份 q。不改 stage 与 enabled。"""
    sp = nq.energy.species
    spp = nq_p.energy.species
    sp2 = replace(
        sp,
        K_eject=np.asarray(spp.K_eject, dtype=np.float64).copy(),
        w_leave=np.asarray(spp.w_leave, dtype=np.float64).copy(),
        w_return=np.asarray(spp.w_return, dtype=np.float64).copy(),
        geo_fb=spp.geo_fb,
        geo_ex=spp.geo_ex,
        k=spp.k,
        D_g=float(spp.D_g),
        rho_g_fb_ref=float(spp.rho_g_fb_ref),
    )
    return replace(
        nq,
        energy=replace(nq.energy, species=sp2),
        hydro=nq_p.hydro,
        reference=nq.reference if nq.reference is not None else nq_p.reference,
    )


def _T_index(stage: int) -> int:
    return n_unknowns(stage) - 1


def unpack_named(x: npt.NDArray[np.float64], stage: int) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    n = n_unknowns(stage)
    assert z.shape == (n,), f"{f_name(stage)} expected {n}, got {z.shape}"
    u = unpack_species(z[:N_SPECIES_BLOCK])
    u["T"] = float(z[-2]) if int(stage) >= 9 else float(z[-1])
    u["T_fb"] = float(z[-1]) if int(stage) >= 9 else float(u["T"])
    off = N_SPECIES_BLOCK
    if int(stage) >= 3:
        u["M_ash"] = float(z[off])
        off += 1
    else:
        u["M_ash"] = 0.0
    if int(stage) >= 6:
        u["N_tar1"] = float(z[off])
        u["N_tar2"] = float(z[off + 1])
        off += 2
    else:
        u["N_tar1"] = float(u.get("N_tar1", 0.0))
        u["N_tar2"] = 0.0
    if int(stage) >= 8:
        u["N_tar1_b"] = float(z[off])
        u["N_tar2_b"] = float(z[off + 1])
        off += 2
    else:
        u["N_tar1_b"] = 0.0
        u["N_tar2_b"] = 0.0
    if int(stage) >= 10:
        u["N_tar1_fb"] = float(z[off])
        u["N_tar2_fb"] = float(z[off + 1])
        off += 2
    else:
        u["N_tar1_fb"] = 0.0
        u["N_tar2_fb"] = 0.0
    if int(stage) >= 11:
        u["M_fb_fine"] = float(z[off])
        u["M_fb_coarse"] = float(z[off + 1])
    else:
        u["M_fb_fine"] = 0.0
        u["M_fb_coarse"] = 0.0
    return u


def pack_named(
    species_x: npt.NDArray[np.float64],
    T: float,
    stage: int,
    *,
    M_ash: float = 0.0,
    N_tar1: float = 0.0,
    N_tar2: float = 0.0,
    N_tar1_b: float = 0.0,
    N_tar2_b: float = 0.0,
    N_tar1_fb: float = 0.0,
    N_tar2_fb: float = 0.0,
    M_fb_fine: float = 0.0,
    M_fb_coarse: float = 0.0,
    T_fb: float | None = None,
) -> npt.NDArray[np.float64]:
    parts = [np.asarray(species_x, dtype=np.float64)]
    if int(stage) >= 3:
        parts.append(np.array([M_ash], dtype=np.float64))
    if int(stage) >= 6:
        parts.append(np.array([N_tar1, N_tar2], dtype=np.float64))
    if int(stage) >= 8:
        parts.append(np.array([N_tar1_b, N_tar2_b], dtype=np.float64))
    if int(stage) >= 10:
        parts.append(np.array([N_tar1_fb, N_tar2_fb], dtype=np.float64))
    if int(stage) >= 11:
        parts.append(np.array([M_fb_fine, M_fb_coarse], dtype=np.float64))
    parts.append(np.array([T], dtype=np.float64))
    if int(stage) >= 9:
        parts.append(np.array([float(T if T_fb is None else T_fb)], dtype=np.float64))
    return np.concatenate(parts)


def default_x0_named(nq: NamedStageQ) -> npt.NDArray[np.float64]:
    x0 = energy_x0(nq.energy)
    species = x0[:-1]
    T = float(x0[-1])
    m_ash = ash_feed_kg_s(nq)
    n1 = float(nq.energy.n_tar.get("TAR1", 0.0))
    n2 = float(nq.energy.n_tar.get("TAR2", 0.0))
    return pack_named(species, T, nq.stage, M_ash=m_ash, N_tar1=n1, N_tar2=n2)


def ash_feed_kg_s(nq: NamedStageQ) -> float:
    feed = nq.energy.feed
    frac = float(feed.ash_dry_wt)
    if frac > 1.0:
        frac = frac / 100.0
    return float(feed.m_dry_kg_s) * frac


def _conc_major(N: npt.NDArray[np.float64], T: float, P: float) -> npt.NDArray[np.float64]:
    ntot = float(np.sum(N))
    if abs(ntot) < 1.0e-30:
        return np.zeros(N_MAJ, dtype=np.float64)
    return (np.asarray(N, dtype=np.float64) / ntot) * (float(P) / (Rg * float(T)))


def r234_dense_source(
    species_x: npt.NDArray[np.float64],
    q_sp: Any,
    *,
    rate_scale: float,
    enabled: frozenset[str] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """仅 R2/R3/R4。R1 已在 F_E0 物种残差，这里不重复。"""
    on = frozenset(("R2", "R3", "R4")) if enabled is None else (
        frozenset(enabled) & frozenset(("R2", "R3", "R4"))
    )
    u = unpack_species(species_x)
    C = _conc_major(u["N_d"], float(q_sp.geo_ex.T), float(q_sp.geo_ex.P))
    src = np.zeros(N_MAJ, dtype=np.float64)
    consume = np.zeros(2, dtype=np.float64)
    if not on:
        return src, consume
    M = _M_bins(u)
    for i in range(2):
        all_src, e_all = char_hetero_source(
            C,
            u["N_d"],
            geo=q_sp.geo_ex,
            M_kg=float(M[i]),
            d_p=float(q_sp.d_p[i]),
            D_g=float(q_sp.D_g),
            rho_s=float(q_sp.rho_p),
            rate_scale=rate_scale,
            enabled=on,
        )
        src = src + all_src
        consume[i] = MW_C_KG_PER_MOL * float(e_all)
    return src, consume


def homog_sources(
    species_x: npt.NDArray[np.float64],
    q_sp: Any,
    *,
    rate_scale: float,
    rho_ash: float,
    enabled: frozenset[str] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    on = frozenset(HOMOG_RXN_NAMES) if enabled is None else (
        frozenset(enabled) & frozenset(HOMOG_RXN_NAMES)
    )
    u = unpack_species(species_x)
    C_d = _conc_major(u["N_d"], float(q_sp.geo_ex.T), float(q_sp.geo_ex.P))
    C_b = _conc_major(u["N_b"], float(q_sp.geo_ex.T), float(q_sp.geo_ex.P))
    y_d, y_b = _y_full(u["N_d"]), _y_full(u["N_b"])
    src_d = gas_homogeneous_source(
        C_d, y_d, geo=q_sp.geo_ex, V=float(q_sp.geo_ex.V_d),
        phase="dense", rate_scale=rate_scale, rho_ash=float(rho_ash),
        enabled=on,
    )
    src_b = gas_homogeneous_source(
        C_b, y_b, geo=q_sp.geo_ex, V=float(q_sp.geo_ex.V_b),
        phase="bubble", rate_scale=rate_scale, rho_ash=0.0,
        enabled=on,
    )
    return src_d, src_b


def c_tar_components(
    N_gas: npt.NDArray[np.float64],
    n_tar1: float,
    n_tar2: float,
    T: float,
    P: float,
) -> tuple[float, float, float]:
    """一相两个焦油代表物的浓度 [mol/m³]，以及二者之和。

    分母是该相主气体与焦油之和。主气体自己的归一化仍不含焦油。
    Hamel §2.2：气体物种按相取浓度。不改 ``rate_R10`` / ``rate_R11_*``。
    """
    n_g = float(np.sum(N_gas))
    a = max(float(n_tar1), 0.0)
    b = max(float(n_tar2), 0.0)
    n_tot = n_g + a + b
    if n_tot < _NTOT_FLOOR:
        return 0.0, 0.0, 0.0
    scale = (float(P) / (Rg * float(T))) / n_tot
    return a * scale, b * scale, (a + b) * scale


def tar_exchange_into_dense(
    u: dict[str, Any],
    q_sp: Any,
) -> tuple[float, float]:
    """焦油相间交换 [mol/s]，正值从气泡进入悬浮。

    与主气体相同：``K_bd V_b (C_b - C_d)``。
    Hamel Eq. 2.4–2.5；§2.2.1 物质 ``j`` 的交换流。
    """
    T = float(u["T"])
    P = float(q_sp.geo_ex.P)
    k_vb = float(q_sp.geo_ex.K_bd) * float(q_sp.geo_ex.V_b)
    c1d, c2d, _ = c_tar_components(u["N_d"], u["N_tar1"], u["N_tar2"], T, P)
    c1b, c2b, _ = c_tar_components(
        u["N_b"], u["N_tar1_b"], u["N_tar2_b"], T, P
    )
    return k_vb * (c1b - c1d), k_vb * (c2b - c2d)


def _phase_c_tar_pair(
    u: dict[str, Any],
    T: float,
    P: float,
    *,
    phase_local: bool,
) -> tuple[float, float]:
    if phase_local:
        _, _, c_d = c_tar_components(u["N_d"], u["N_tar1"], u["N_tar2"], T, P)
        _, _, c_b = c_tar_components(u["N_b"], u["N_tar1_b"], u["N_tar2_b"], T, P)
        return c_d, c_b
    c = r10_c_tar(u["N_d"], u["N_b"], float(u["N_tar1"]) + float(u["N_tar2"]), T, P)
    return c, c


def r10_nu_major(fuel_type: str = R10_FUEL_TYPE) -> npt.NDArray[np.float64]:
    """Lumped R10 对 MAJOR 的计量。TAR 槽不进这 7 维。"""
    sto = get_lumped_tar_stoichiometry("R10", fuel_type)
    nu = np.zeros(N_MAJ, dtype=np.float64)
    for name, val in sto.items():
        if name in IDX:
            nu[IDX[name]] = float(val)
    return nu


def r10_c_tar(
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    n_tar: float,
    T: float,
    P: float,
) -> float:
    """稀组分共享槽：全格焦油摩尔分率 × P/RT [mol/m³]。两相同一 C_tar。

    主气体各相仍按 ``sum_j C_{p,j}=P/RT`` 归一；焦油不进入该分母。
    因此若把焦油当成同一物相真实组分，总浓度超过 P/RT，超出量为 y_tar。
    这是 R10/R11 采用的稀组分近似，不是严格 EOS。
    """
    n_gas = float(np.sum(N_d) + np.sum(N_b))
    n_t = max(float(n_tar), 0.0)
    n_tot = n_gas + n_t
    if n_tot < _NTOT_FLOOR:
        return 0.0
    return (n_t / n_tot) * (float(P) / (Rg * float(T)))


def tar_concentration_diagnostics(
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    n_tar: float,
    T: float,
    P: float,
) -> dict[str, Any]:
    """焦油共享槽 vs 严格两相 EOS。不改 R10/R11 速率。

    Source: docs/audit_reports/named_f_e8_review_2026-09-16.md §2.7
    """
    C_star = float(P) / (Rg * float(T))
    n_gas = float(np.sum(N_d) + np.sum(N_b))
    n_t = max(float(n_tar), 0.0)
    n_tot = n_gas + n_t
    y_tar = n_t / n_tot if n_tot >= _NTOT_FLOOR else 0.0
    C_tar = y_tar * C_star
    return {
        "model_used_for_R10_R11": TAR_CONC_DILUTE_SHARED_SLOT,
        "strict_eos_model": TAR_CONC_STRICT_EOS,
        "phase_allocation": "instantaneous_well_mixed_equal_C_dense_and_bubble",
        "shared_slot": "bed_exit_tar_flow_not_third_eos_species",
        "tar_denominator": "n_gas_bed + n_tar",
        "major_denominator": "n_gas_bed",
        "C_star_mol_m3": C_star,
        "y_tar_dilute": y_tar,
        "C_tar_dilute_mol_m3": C_tar,
        "dilute_total_C_over_P_RT": 1.0 + y_tar,
        "strict_eos_total_C_over_P_RT": 1.0,
        "eos_relative_excess": y_tar,
        "not_freeboard_tar": True,
    }


def r10_major_and_tar_source(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float, float]:
    """R10 气源 [mol/s] 与焦油源 [mol/s]。只进密相/气泡 MAJOR。

    Hamel Eq. 5.59：两相同一 ``k T P^{0.3} C_tar^{0.5} C_O2``，各用该相浓度。
    E8 起 C_tar 按相；E6–E7 仍是全格共享槽。不改 ``rate_R10``。
    """
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    T = float(u["T"])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    phase_local = int(nq.stage) >= 8
    C_tar_d, C_tar_b = _phase_c_tar_pair(u, T, P, phase_local=phase_local)
    C_d = _conc_major(u["N_d"], T, P)
    C_b = _conc_major(u["N_b"], T, P)
    s = float(rate_scale)
    xi_d = s * float(rate_R10(T, C_tar_d, float(C_d[IDX["O2"]]), P, R10_FUEL_TYPE)) * V_d
    xi_b = s * float(rate_R10(T, C_tar_b, float(C_b[IDX["O2"]]), P, R10_FUEL_TYPE)) * V_b
    nu = r10_nu_major()
    sto = get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)
    xi = xi_d + xi_b
    return (
        nu * xi_d,
        nu * xi_b,
        float(sto.get("TAR1", 0.0)) * xi,
        float(sto.get("TAR2", 0.0)) * xi,
    )


def r11_nu_major(fuel_type: str = R10_FUEL_TYPE) -> npt.NDArray[np.float64]:
    """Lumped R11 对 MAJOR 的计量。TAR 槽不进这 7 维。"""
    sto = get_lumped_tar_stoichiometry("R11", fuel_type)
    nu = np.zeros(N_MAJ, dtype=np.float64)
    for name, val in sto.items():
        if name in IDX:
            nu[IDX[name]] = float(val)
    return nu


def r11_major_and_tar_source(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float, float]:
    """R11 气源 [mol/s] 与焦油源 [mol/s]。不改 ``rate_R11_*``。

    气泡：Serio 一级 ``k(T) C_tar``。悬浮：Corella ``k(T) rho_cat C_tar``。
    E8 起各用该相 C_tar。速率不含 C_H2O；H2O 只出现在计量。
    """
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    T = float(u["T"])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    phase_local = int(nq.stage) >= 8
    C_tar_d, C_tar_b = _phase_c_tar_pair(u, T, P, phase_local=phase_local)
    rho_cat = rho_ash_used(nq, u)
    s = float(rate_scale)
    xi_d = s * float(rate_R11_suspension(T, C_tar_d, rho_cat)) * V_d
    xi_b = s * float(rate_R11_bubble(T, C_tar_b)) * V_b
    nu = r11_nu_major()
    sto = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
    xi = xi_d + xi_b
    return (
        nu * xi_d,
        nu * xi_b,
        float(sto.get("TAR1", 0.0)) * xi,
        float(sto.get("TAR2", 0.0)) * xi,
    )


def tar_consumption_fractions() -> tuple[float, float]:
    """R10/R11 对 TAR1、TAR2 的消耗比例 ``(a1, a2)``。两者相同且和为 1。

    这是数值行变换用的系数，不是另一套 Hamel 算法。
    """
    out = []
    for rxn in ("R10", "R11"):
        sto = get_lumped_tar_stoichiometry(rxn, R10_FUEL_TYPE)
        a1 = -float(sto["TAR1"])
        a2 = -float(sto["TAR2"])
        out.append((a1, a2))
    if abs(out[0][0] - out[1][0]) > 1.0e-12 or abs(out[0][1] - out[1][1]) > 1.0e-12:
        raise ValueError("R10 与 R11 的焦油消耗比例不一致")
    a1, a2 = out[0]
    if abs(a1 + a2 - 1.0) > 1.0e-12:
        raise ValueError("焦油消耗比例之和不是 1")
    return a1, a2


def tar_sum_and_fraction(n1: float, n2: float, a1: float, a2: float) -> tuple[float, float]:
    """``S=N1+N2``，``f=D/S``，``D=a2*N1-a1*N2``。``S=0`` 时比例记为 0。"""
    total = float(n1) + float(n2)
    if total == 0.0:
        return 0.0, 0.0
    fraction = (float(a2) * float(n1) - float(a1) * float(n2)) / total
    return total, fraction


def inventories_from_sum_fraction(
    total: float,
    fraction: float,
    a1: float,
    a2: float,
) -> tuple[float, float]:
    """由总量和比例回到两种库存。比例夹在 ``[-a1, a2]``，两种库存都非负。"""
    span = float(total)
    if span < 0.0:
        span = 0.0
    frac = min(max(float(fraction), -float(a1)), float(a2))
    return span * (float(a1) + frac), span * (float(a2) - frac)


def linear_tar_row_map(row: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """焦油行的纯线性映射 ``T``。``T(0)=0``，可加，不读取状态。

    ``G_sum=F1+F2``，``G_ratio=a2*F1-a1*F2``。参数增量和矩阵列用它。
    真实状态的比例行仍由 ``apply_tar_reaction_transport_rows`` 按输运组装。
    """
    out = np.asarray(row, dtype=np.float64).copy()
    if int(out.size) < 46:
        return out
    a1, a2 = tar_consumption_fractions()
    for i in (40, 42, 44):
        f1 = float(out[i])
        f2 = float(out[i + 1])
        out[i] = f1 + f2
        out[i + 1] = float(a2) * f1 - float(a1) * f2
    return out


def original_residual_from_tar_rows(
    row: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """``linear_tar_row_map`` 的逆。真实状态比例行另由输运组装，不能用这个函数代替。"""
    out = np.asarray(row, dtype=np.float64).copy()
    if int(out.size) < 46:
        return out
    a1, a2 = tar_consumption_fractions()
    for i in (40, 42, 44):
        total = float(out[i])
        ratio = float(out[i + 1])
        out[i] = a1 * total + ratio
        out[i + 1] = a2 * total - ratio
    return out


def tar_transport_ratio(
    n1: float,
    n2: float,
    inflow1: float,
    inflow2: float,
    *,
    a1: float,
    a2: float,
) -> float:
    """``a2*T1 - a1*T2`` [mol/s]。``T`` 不含反应源，所以反应进度在这一行里精确消失。"""
    t1 = float(n1) - float(inflow1)
    t2 = float(n2) - float(inflow2)
    return float(a2) * t1 - float(a1) * t2


def apply_tar_reaction_transport_rows(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    residual: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """把每相两行焦油残差换成总和与输运比例。只用于内层校正，不改变验收残差。

    ``G_sum = F1+F2`` 保留反应方向。``G_ratio`` 由入口、出口和交换组装。
    原方程 ``F1=F2=0`` 与 ``G_sum=G_ratio=0`` 等价。
    """
    row = np.asarray(residual, dtype=np.float64).copy()
    if int(row.size) < 46 or int(nq.stage) < 10:
        return row
    a1, a2 = tar_consumption_fractions()
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    qT = q_at_T(nq.energy.species, float(u["T"]))
    nex1, nex2 = tar_exchange_into_dense(u, qT)
    feed1 = float(nq.energy.n_tar.get("TAR1", 0.0))
    feed2 = float(nq.energy.n_tar.get("TAR2", 0.0))
    pairs = (
        (40, 41, tar_transport_ratio(
            u["N_tar1"], u["N_tar2"], feed1 + nex1, feed2 + nex2, a1=a1, a2=a2,
        )),
        (42, 43, tar_transport_ratio(
            u["N_tar1_b"], u["N_tar2_b"], -nex1, -nex2, a1=a1, a2=a2,
        )),
        (44, 45, tar_transport_ratio(
            u["N_tar1_fb"], u["N_tar2_fb"],
            float(u["N_tar1"]) + float(u["N_tar1_b"]),
            float(u["N_tar2"]) + float(u["N_tar2_b"]),
            a1=a1, a2=a2,
        )),
    )
    for i1, i2, ratio in pairs:
        total = float(row[i1]) + float(row[i2])
        row[i1] = total
        row[i2] = float(ratio)
    return row


def tar_sources_by_phase(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> tuple[float, float, float, float]:
    """``(TAR1, TAR2)`` 的悬浮源与气泡源 [mol/s]。

    与 ``r10`` / ``r11`` 的 ``xi`` 同一口径。E8 用各相浓度；更早阶段两相共用 C_tar。
    热解给料不在这里，只在悬浮相残差。
    """
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    T = float(u["T"])
    P = float(nq.energy.species.geo_ex.P)
    V_d = float(nq.energy.species.geo_ex.V_d)
    V_b = float(nq.energy.species.geo_ex.V_b)
    phase_local = int(nq.stage) >= 8
    c_d, c_b = _phase_c_tar_pair(u, T, P, phase_local=phase_local)
    C_d = _conc_major(u["N_d"], T, P)
    C_b = _conc_major(u["N_b"], T, P)
    s = float(rate_scale)
    s1d = s2d = s1b = s2b = 0.0
    if "R10" in nq.enabled_rxns:
        xi_d = s * float(rate_R10(T, c_d, float(C_d[IDX["O2"]]), P, R10_FUEL_TYPE)) * V_d
        xi_b = s * float(rate_R10(T, c_b, float(C_b[IDX["O2"]]), P, R10_FUEL_TYPE)) * V_b
        sto = get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)
        s1d += float(sto.get("TAR1", 0.0)) * xi_d
        s2d += float(sto.get("TAR2", 0.0)) * xi_d
        s1b += float(sto.get("TAR1", 0.0)) * xi_b
        s2b += float(sto.get("TAR2", 0.0)) * xi_b
    if "R11" in nq.enabled_rxns:
        rho_cat = rho_ash_used(nq, u)
        xi_d = s * float(rate_R11_suspension(T, c_d, rho_cat)) * V_d
        xi_b = s * float(rate_R11_bubble(T, c_b)) * V_b
        sto = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
        s1d += float(sto.get("TAR1", 0.0)) * xi_d
        s2d += float(sto.get("TAR2", 0.0)) * xi_d
        s1b += float(sto.get("TAR1", 0.0)) * xi_b
        s2b += float(sto.get("TAR2", 0.0)) * xi_b
    return s1d, s2d, s1b, s2b


def fb_tar_sources(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> tuple[float, float, npt.NDArray[np.float64]]:
    """自由空间焦油源 [mol/s] 与对应主气体源。

    输运进来的是床层两相焦油之和。自由空间没有催化剂，R11 只用气泡相的
    Serio 气相裂解；R10 用自由空间自己的 ``C_tar`` 和 ``C_O2``。不改速率式。
    """
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    T = float(u["T_fb"])
    P = float(nq.energy.species.geo_ex.P)
    V_fb = float(nq.energy.species.geo_fb.area) * float(nq.energy.species.geo_fb.H_freeboard)
    c1, c2, c_sum = c_tar_components(
        u["N_fb"], u["N_tar1_fb"], u["N_tar2_fb"], T, P,
    )
    _ = (c1, c2)
    C_fb = _conc_major(u["N_fb"], T, P)
    s = float(rate_scale)
    major = np.zeros(N_MAJ, dtype=np.float64)
    s1 = s2 = 0.0
    if "R10" in nq.enabled_rxns:
        xi = s * float(rate_R10(T, c_sum, float(C_fb[IDX["O2"]]), P, R10_FUEL_TYPE)) * V_fb
        sto = get_lumped_tar_stoichiometry("R10", R10_FUEL_TYPE)
        s1 += float(sto.get("TAR1", 0.0)) * xi
        s2 += float(sto.get("TAR2", 0.0)) * xi
        major = major + r10_nu_major() * xi
    if "R11" in nq.enabled_rxns:
        xi = s * float(rate_R11_bubble(T, c_sum)) * V_fb
        sto = get_lumped_tar_stoichiometry("R11", R10_FUEL_TYPE)
        s1 += float(sto.get("TAR1", 0.0)) * xi
        s2 += float(sto.get("TAR2", 0.0)) * xi
        major = major + r11_nu_major() * xi
    return s1, s2, major


def fb_char_balance(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """自由空间两档炭库存 [kg] 的稳态衡算。选定集总式，不用轨迹剩余质量积分。

    无反应时零点是 ``M = m_eject * tau_s``。``tau_s`` 是已判定轨迹的平均飞行时间，
    包括回落。``w_leave=0`` 时离开流为零，库存仍然可以不为零。

    有反应时零点是 ``M = m_eject * tau / (1 + k*tau)``，其中 ``k*M`` 就是本函数算出的
    ``consume``。这假设各条路径的损失比例相同。轨迹积分
    ``M = Σ m_j ∫ S_j(t) dt`` 不在这里计算，也不和集总式混用。
    出口上的消耗另按飞行时间拆开。未判定份额不进出口。不改速率式。
    细粒床库存在 ``1e-6 kg/s`` 以下不可分辨，不从 ``1e-33 kg`` 读出细粒。
    """
    z = np.asarray(x, dtype=np.float64)
    u = unpack_named(z, nq.stage)
    T = float(u["T_fb"])
    qT = q_at_T(nq.energy.species, T)
    ed = edges_from_x(z[:N_SPECIES_BLOCK], qT)
    if nq.hydro is not None:
        tau_raw = np.asarray(nq.hydro.tau_s, dtype=np.float64)
        und = np.asarray(nq.hydro.w_undecided, dtype=np.float64)
        closed = und <= 1.0e-12
        tau = np.where(closed, np.maximum(tau_raw, 1.0e-6), np.inf)
    else:
        tau = np.full(
            2,
            float(qT.geo_fb.H_freeboard) / max(float(qT.geo_fb.u_g), 1.0e-3),
            dtype=np.float64,
        )
    M = np.array([u["M_fb_fine"], u["M_fb_coarse"]], dtype=np.float64)
    m_eject = np.asarray(ed["m_eject"], dtype=np.float64)
    on = nq.enabled_rxns & frozenset(("R2", "R3", "R4"))
    C = _conc_major(unpack_species(z[:N_SPECIES_BLOCK])["N_fb"], T, float(qT.geo_ex.P))
    src = np.zeros(N_MAJ, dtype=np.float64)
    consume = np.zeros(2, dtype=np.float64)
    if on:
        for i in range(2):
            all_src, e_all = char_hetero_source(
                C,
                unpack_species(z[:N_SPECIES_BLOCK])["N_fb"],
                geo=qT.geo_ex,
                M_kg=float(M[i]),
                d_p=float(qT.d_p[i]),
                D_g=float(qT.D_g),
                rho_s=float(qT.rho_p),
                rate_scale=rate_scale,
                enabled=on,
            )
            src = src + all_src
            consume[i] = MW_C_KG_PER_MOL * float(e_all)
    resid = m_eject - M / tau - consume
    return resid, src, consume


def fb_gas_source(
    species_x: npt.NDArray[np.float64],
    q_sp: Any,
    *,
    rate_scale: float,
    enabled: frozenset[str] | None = None,
) -> npt.NDArray[np.float64]:
    """自由板单槽均相源 [mol/s]。只用 R5/R6/R8/R12，气泡动力学，ρ_ash=0。

    ``q_sp.geo_ex.T`` 须已是状态 T（``q_at_T``），不是冻结的 geo 1100 K。
    不含 R10/R11。R7 因 ``ρ_ash=0`` 归零（原文仅载固悬浮相）。
    不是生产 ``_freeboard_reaction_step``。
    """
    on = frozenset(HOMOG_RXN_NAMES) if enabled is None else (
        frozenset(enabled) & frozenset(HOMOG_RXN_NAMES)
    )
    u = unpack_species(species_x)
    C = _conc_major(u["N_fb"], float(q_sp.geo_ex.T), float(q_sp.geo_ex.P))
    y = _y_full(u["N_fb"])
    V_fb = float(q_sp.geo_fb.area) * float(q_sp.geo_fb.H_freeboard)
    return gas_homogeneous_source(
        C, y, geo=q_sp.geo_ex, V=V_fb, phase="bubble",
        rate_scale=rate_scale, rho_ash=0.0, enabled=on,
    )


def rho_ash_used(nq: NamedStageQ, u: dict[str, Any]) -> float:
    if int(nq.stage) >= 3:
        V_d = float(nq.energy.species.geo_ex.V_d)
        return float(u["M_ash"]) / max(V_d, 1.0e-12)
    return float(nq.rho_ash)


def species_block_residual(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float,
    char_pack: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]] | None = None,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    stage = int(nq.stage)
    T_bed = float(z[-2]) if stage >= 9 else float(z[-1])
    T_fb = float(z[-1]) if stage >= 9 else T_bed
    qT = q_at_T(nq.energy.species, T_bed)
    sp = z[:N_SPECIES_BLOCK]
    F = species_residual(sp, qT, rate_scale=rate_scale)
    u = unpack_named(z, nq.stage)
    on = nq.enabled_rxns
    if on & frozenset(("R2", "R3", "R4")):
        src234, cons234 = r234_dense_source(sp, qT, rate_scale=rate_scale, enabled=on)
        F[0:7] = F[0:7] + src234
        F[14] = F[14] - cons234[IDX_FINE]
        F[15] = F[15] - cons234[IDX_COARSE]
    if on & frozenset(HOMOG_RXN_NAMES):
        rho = rho_ash_used(nq, u)
        sd, sb = homog_sources(sp, qT, rate_scale=rate_scale, rho_ash=rho, enabled=on)
        F[0:7] = F[0:7] + sd
        F[7:14] = F[7:14] + sb
    if "R10" in on and int(nq.stage) >= 6:
        sd10, sb10, _t1, _t2 = r10_major_and_tar_source(z, nq, rate_scale=rate_scale)
        F[0:7] = F[0:7] + sd10
        F[7:14] = F[7:14] + sb10
    if "R11" in on and int(nq.stage) >= 6:
        sd11, sb11, _u1, _u2 = r11_major_and_tar_source(z, nq, rate_scale=rate_scale)
        F[0:7] = F[0:7] + sd11
        F[7:14] = F[7:14] + sb11
    if int(nq.stage) >= 7:
        F[16:23] = F[16:23] - fb_gas_source(sp, qT, rate_scale=rate_scale, enabled=on)
        if int(nq.stage) >= 9 and T_fb != T_bed:
            q_fb = q_at_T(nq.energy.species, T_fb)
            F[16:23] = (
                F[16:23]
                + fb_gas_source(sp, qT, rate_scale=rate_scale, enabled=on)
                - fb_gas_source(sp, q_fb, rate_scale=rate_scale, enabled=on)
            )
    if int(nq.stage) >= 10 and "R10" in on or (int(nq.stage) >= 10 and "R11" in on):
        _s1, _s2, major_fb = fb_tar_sources(z, nq, rate_scale=rate_scale)
        F[16:23] = F[16:23] - major_fb
    if int(nq.stage) >= 11:
        if char_pack is None:
            char_pack = fb_char_balance(z, nq, rate_scale=rate_scale)
        _resid, src_c, _cons = char_pack
        _eject, m_return, m_leave = char_face_flows_kg_s(
            z, nq, rate_scale=rate_scale, char_pack=char_pack,
        )
        ed_raw = edges_from_x(sp, qT)
        on_return = np.asarray(ed_raw["m_return"], dtype=np.float64) - m_return
        on_leave = np.asarray(ed_raw["m_leave"], dtype=np.float64) - m_leave
        rez_post = recycle_bins_from_leave(
            np.asarray(m_leave, dtype=np.float64),
            float(ed_raw["m_rez"]),
            np.asarray(ed_raw["m_rez_bin"], dtype=np.float64),
        )
        delta_rez = rez_post - np.asarray(ed_raw["m_rez_bin"], dtype=np.float64)
        F[14] = F[14] - float(on_return[0]) + float(delta_rez[0])
        F[15] = F[15] - float(on_return[1]) + float(delta_rez[1])
        F[16:23] = F[16:23] - src_c
        F[37] = F[37] - float(np.sum(on_leave))
    return F


def ash_residual(nq: NamedStageQ, u: dict[str, Any], ed: dict[str, Any]) -> float:
    """Eq. 2.6 灰：给料 − 隐含飞灰分流。灰不参与 R1–R4。"""
    k = nq.energy.species.k
    k_sum = float(k.K_auf) + float(k.K_ab)
    fly_frac = float(k.K_auf) / max(k_sum, 1.0e-30)
    K_e = float(nq.energy.species.K_eject[IDX_COARSE])
    return ash_feed_kg_s(nq) - fly_frac * K_e * float(u["M_ash"])


def minor_bypass_enthalpy_W(species: Any, T: float) -> float:
    """NH3/H2S 旁路焓流 [W]。不进入七组分密度，也不代表后续反应。"""
    from src.core.species import enthalpy_molar

    return (
        float(species.n_NH3_bypass_mol_s) * float(enthalpy_molar("NH3", float(T)))
        + float(species.n_H2S_bypass_mol_s) * float(enthalpy_molar("H2S", float(T)))
    )


def energy_row_named(x: npt.NDArray[np.float64], nq: NamedStageQ) -> float:
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-1])
    qe = nq.energy
    qT = q_at_T(qe.species, T)
    sp = z[:N_SPECIES_BLOCK]
    u_sp = unpack_species(sp)
    ed = edges_from_x(sp, qT)
    from src.core.energy_balance_defaults import (
        H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        H_FORMATION_LIQUID_WATER_J_PER_KG,
    )
    from src.core.species import SOLID_H_MODE_INTEGRAL, h_solid_mass

    H_in = (
        _H_major(qe.species.N_in_d, qe.T_in)
        + _H_major(qe.species.N_in_b, qe.T_in)
        + float(qe.feed.m_dry_kg_s) * float(qe.feed.h_f_dry_J_kg)
        + float(qe.feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    )
    h_char = h_solid_mass(
        "char", T, h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG, mode=SOLID_H_MODE_INTEGRAL
    )
    n_tar = dict(qe.n_tar)
    if int(nq.stage) >= 6:
        n_tar = {"TAR1": float(unpack_named(z, nq.stage)["N_tar1"]),
                 "TAR2": float(unpack_named(z, nq.stage)["N_tar2"])}
        if int(nq.stage) >= 8:
            ub = unpack_named(z, nq.stage)
            n_tar["TAR1"] = float(ub["N_tar1"]) + float(ub["N_tar1_b"])
            n_tar["TAR2"] = float(ub["N_tar2"]) + float(ub["N_tar2_b"])
    H_unreleased = unreleased_vm_enthalpy_W(
        qe.feed, float(qe.species.m_unreleased_kg_s), T
    )
    H_out = (
        _H_major(u_sp["N_c2"], T)
        + float(ed["m_fly"]) * h_char
        + _H_tar(n_tar, T)
        + H_unreleased
    )
    if int(nq.stage) >= 3:
        h_ash = h_solid_mass("ash", T, h_f=-1.2e6, mode=SOLID_H_MODE_INTEGRAL)
        k = qe.species.k
        fly_frac = float(k.K_auf) / max(float(k.K_auf) + float(k.K_ab), 1.0e-30)
        K_e = float(qe.species.K_eject[IDX_COARSE])
        m_fly_ash = fly_frac * K_e * float(unpack_named(z, nq.stage)["M_ash"])
        H_out = H_out + m_fly_ash * h_ash
    H_out = H_out + minor_bypass_enthalpy_W(qe.species, T)
    Q_W = float(nq.Q_W) if int(nq.stage) >= 5 else 0.0
    return float(H_in - H_out - Q_W)


def recycle_bins_from_leave(
    m_leave_post: npt.NDArray[np.float64],
    m_rez: float,
    m_rez_bin_inlet: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """旋风回料 [kg/s] 按到达旋风的反应后离开流分到粒级。

    总回料仍是 ``K_ab M_c1``。两档离开流都非负且总量为正时，用反应后比例。
    离开流总量为零时粒级组成未定，保持调用方传入的分配，不把回料抹成零。
    """
    leave = np.asarray(m_leave_post, dtype=np.float64)
    inlet = np.asarray(m_rez_bin_inlet, dtype=np.float64).copy()
    total = float(np.sum(leave))
    if total > 1.0e-30 and float(np.min(leave)) >= -1.0e-15:
        return float(m_rez) * leave / total
    return inlet


def path_char_consumption_kg_s(
    consume_kg_s: npt.NDArray[np.float64],
    tau_return_s: npt.NDArray[np.float64],
    tau_leave_s: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """把一份炭消耗按路径飞行时间分成回落与越格 [kg/s]。

    时间是 ``w·t`` 的和。两条路径时间不同时，不按 ``w_return`` 去拆。
    未判定时间不在这两个数组里，因此不进入出口。
    """
    cons = np.asarray(consume_kg_s, dtype=np.float64)
    t_ret = np.maximum(np.asarray(tau_return_s, dtype=np.float64), 0.0)
    t_leave = np.maximum(np.asarray(tau_leave_s, dtype=np.float64), 0.0)
    total = t_ret + t_leave
    on_return = np.divide(cons * t_ret, total, out=np.zeros_like(cons), where=total > 0.0)
    on_leave = np.divide(cons * t_leave, total, out=np.zeros_like(cons), where=total > 0.0)
    return on_return, on_leave


def char_face_flows_kg_s(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float = 1.0,
    char_pack: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """床面炭的抛出、反应后回落、反应后离开 [kg/s]。

    回落和离开扣的是同一份消耗，按各自飞行时间分配。E11 以下没有自由空间消耗。
    """
    z = np.asarray(x, dtype=np.float64)
    sp = nq.energy.species
    ed = edges_from_x(z[:N_SPECIES_BLOCK], sp)
    m_eject = np.asarray(ed["m_eject"], dtype=np.float64)
    consume = np.zeros(2, dtype=np.float64)
    on_return = np.zeros(2, dtype=np.float64)
    on_leave = np.zeros(2, dtype=np.float64)
    if int(nq.stage) >= 11:
        if char_pack is None:
            char_pack = fb_char_balance(z, nq, rate_scale=rate_scale)
        _resid, _src, consume = char_pack
        hydro = nq.hydro
        if hydro is not None:
            on_return, on_leave = path_char_consumption_kg_s(
                consume, hydro.tau_return_s, hydro.tau_leave_s,
            )
        else:
            on_return = np.asarray(sp.w_return, dtype=np.float64) * consume
            on_leave = np.asarray(sp.w_leave, dtype=np.float64) * consume
    m_return = np.asarray(ed["m_return"], dtype=np.float64) - on_return
    m_leave = np.asarray(ed["m_leave"], dtype=np.float64) - on_leave
    return m_eject, m_return, m_leave


def char_face_enthalpy_w(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    T_bed: float,
    T_fb: float,
    *,
    rate_scale: float = 1.0,
    char_pack: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]] | None = None,
) -> float:
    """跨床面炭焓流 [W]。床层失去它，自由空间得到它。

    向上用床温，向下用自由空间温度。全部回落、无反应、两区同温时净值为 0。
    """
    from src.core.energy_balance_defaults import H_FORMATION_GRAPHITE_CHAR_J_PER_KG
    from src.core.species import SOLID_H_MODE_INTEGRAL, h_solid_mass

    m_eject, m_return, _m_leave = char_face_flows_kg_s(
        x, nq, rate_scale=rate_scale, char_pack=char_pack,
    )
    h_up = h_solid_mass(
        "char", float(T_bed), h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    h_down = h_solid_mass(
        "char", float(T_fb), h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    return float(np.sum(m_eject) * h_up - np.sum(m_return) * h_down)


def energy_zones_named(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    T_bed: float,
    T_fb: float,
    char_pack: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]] | None = None,
) -> tuple[float, float]:
    """床层与自由空间能量 [W]。旋风无热损，透传不另列方程。

    离开床面的气体 ``N_d+N_b`` 与固体 ``m_leave`` 用床温，回落固体 ``m_return``
    用自由空间温度。同一股焓流在两区等大反号。
    ``T_bed=T_fb`` 时两区之和等于 ``energy_row_named``。

    Source: docs/audit_reports/named_f_e8_review_2026-09-16.md §5 P3
    """
    z = np.asarray(x, dtype=np.float64)
    qe = nq.energy
    sp = z[:N_SPECIES_BLOCK]
    u_sp = unpack_species(sp)
    q_bed = q_at_T(qe.species, float(T_bed))
    ed = edges_from_x(sp, q_bed)
    from src.core.energy_balance_defaults import (
        H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        H_FORMATION_LIQUID_WATER_J_PER_KG,
    )
    from src.core.species import SOLID_H_MODE_INTEGRAL, h_solid_mass

    H_in = (
        _H_major(qe.species.N_in_d, qe.T_in)
        + _H_major(qe.species.N_in_b, qe.T_in)
        + float(qe.feed.m_dry_kg_s) * float(qe.feed.h_f_dry_J_kg)
        + float(qe.feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    )
    h_leave = h_solid_mass(
        "char", float(T_bed), h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    h_return = h_solid_mass(
        "char", float(T_fb), h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    h_fly = h_return
    n_tar = dict(qe.n_tar)
    n_tar_cross = {"TAR1": 0.0, "TAR2": 0.0}
    if int(nq.stage) >= 6:
        u = unpack_named(z, nq.stage)
        n_tar = {"TAR1": float(u["N_tar1"]), "TAR2": float(u["N_tar2"])}
        if int(nq.stage) >= 8:
            n_tar["TAR1"] = float(u["N_tar1"]) + float(u["N_tar1_b"])
            n_tar["TAR2"] = float(u["N_tar2"]) + float(u["N_tar2_b"])
        if int(nq.stage) >= 10:
            n_tar_cross = {"TAR1": n_tar["TAR1"], "TAR2": n_tar["TAR2"]}
            n_tar = {"TAR1": float(u["N_tar1_fb"]), "TAR2": float(u["N_tar2_fb"])}
    H_gas_cross = _H_major(u_sp["N_d"] + u_sp["N_b"], float(T_bed))
    H_char_cross = char_face_enthalpy_w(
        z, nq, float(T_bed), float(T_fb), char_pack=char_pack,
    )
    H_ash_cross = 0.0
    m_fly_ash = 0.0
    h_ash_out = 0.0
    if int(nq.stage) >= 3:
        h_ash_bed = h_solid_mass("ash", float(T_bed), h_f=-1.2e6, mode=SOLID_H_MODE_INTEGRAL)
        h_ash_out = h_solid_mass("ash", float(T_fb), h_f=-1.2e6, mode=SOLID_H_MODE_INTEGRAL)
        k = qe.species.k
        fly_frac = float(k.K_auf) / max(float(k.K_auf) + float(k.K_ab), 1.0e-30)
        K_e = float(qe.species.K_eject[IDX_COARSE])
        m_fly_ash = fly_frac * K_e * float(unpack_named(z, nq.stage)["M_ash"])
        H_ash_cross = m_fly_ash * h_ash_bed
    H_cross = (
        H_gas_cross
        + H_char_cross
        + H_ash_cross
        + _H_tar(n_tar_cross, float(T_bed))
        + minor_bypass_enthalpy_W(qe.species, float(T_bed))
    )
    H_out = (
        _H_major(u_sp["N_c2"], float(T_fb))
        + float(ed["m_fly"]) * h_fly
        + _H_tar(n_tar, float(T_fb))
        + unreleased_vm_enthalpy_W(
            qe.feed, float(qe.species.m_unreleased_kg_s), float(T_fb)
        )
    )
    if int(nq.stage) >= 3:
        H_out = H_out + m_fly_ash * h_ash_out
    H_out = H_out + minor_bypass_enthalpy_W(qe.species, float(T_fb))
    Q_W = float(nq.Q_W) if int(nq.stage) >= 5 else 0.0
    e_bed = float(H_in - H_cross - Q_W)
    e_fb = float(H_cross - H_out)
    return e_bed, e_fb


def coupled_residual_named(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
    *,
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (n_unknowns(nq.stage),)
    char_pack = (
        fb_char_balance(z, nq, rate_scale=rate_scale) if int(nq.stage) >= 11 else None
    )
    F_s = species_block_residual(z, nq, rate_scale=rate_scale, char_pack=char_pack)
    T = float(z[-1])
    qT = q_at_T(nq.energy.species, T)
    parts: list[npt.NDArray[np.float64]] = [F_s]
    u = unpack_named(z, nq.stage)
    if int(nq.stage) >= 3:
        ed = edges_from_x(z[:N_SPECIES_BLOCK], qT)
        parts.append(np.array([ash_residual(nq, u, ed)], dtype=np.float64))
    if int(nq.stage) >= 6:
        src_t1 = 0.0
        src_t2 = 0.0
        if "R10" in nq.enabled_rxns:
            _sd, _sb, src_t1, src_t2 = r10_major_and_tar_source(
                z, nq, rate_scale=rate_scale
            )
        if "R11" in nq.enabled_rxns:
            _sd11, _sb11, t1, t2 = r11_major_and_tar_source(
                z, nq, rate_scale=rate_scale
            )
            src_t1 = src_t1 + t1
            src_t2 = src_t2 + t2
        if int(nq.stage) >= 8:
            s1d, s2d, s1b, s2b = tar_sources_by_phase(z, nq, rate_scale=rate_scale)
            nex1, nex2 = tar_exchange_into_dense(u, qT)
            feed1 = float(nq.energy.n_tar.get("TAR1", 0.0))
            feed2 = float(nq.energy.n_tar.get("TAR2", 0.0))
            parts.append(
                np.array(
                    [
                        float(u["N_tar1"]) - feed1 - s1d - nex1,
                        float(u["N_tar2"]) - feed2 - s2d - nex2,
                        float(u["N_tar1_b"]) - s1b + nex1,
                        float(u["N_tar2_b"]) - s2b + nex2,
                    ],
                    dtype=np.float64,
                )
            )
            if int(nq.stage) >= 10:
                s1f, s2f, _maj = fb_tar_sources(z, nq, rate_scale=rate_scale)
                parts.append(
                    np.array(
                        [
                            float(u["N_tar1_fb"]) - float(u["N_tar1"]) - float(u["N_tar1_b"]) - s1f,
                            float(u["N_tar2_fb"]) - float(u["N_tar2"]) - float(u["N_tar2_b"]) - s2f,
                        ],
                        dtype=np.float64,
                    )
                )
            if int(nq.stage) >= 11:
                assert char_pack is not None
                parts.append(np.asarray(char_pack[0], dtype=np.float64))
        else:
            parts.append(
                np.array(
                    [
                        float(u["N_tar1"]) - float(nq.energy.n_tar.get("TAR1", 0.0)) - src_t1,
                        float(u["N_tar2"]) - float(nq.energy.n_tar.get("TAR2", 0.0)) - src_t2,
                    ],
                    dtype=np.float64,
                )
            )
    if int(nq.stage) >= 9:
        T_bed = float(z[-2])
        T_fb = float(z[-1])
        e_bed, e_fb = energy_zones_named(z, nq, T_bed, T_fb, char_pack=char_pack)
        parts.append(np.array([e_bed, e_fb], dtype=np.float64))
    else:
        parts.append(np.array([energy_row_named(z, nq)], dtype=np.float64))
    return np.concatenate(parts)


def daem_X_at_T(T: float, tau_s: float) -> float:
    t = np.array([0.0, float(tau_s)], dtype=np.float64)
    Th = np.array([float(T), float(T)], dtype=np.float64)
    return float(daem_conversion(Th, t))


def rebuild_feed_derived_inlets(nq: NamedStageQ) -> NamedStageQ:
    """干燥水、热解源和未释放物流从当前 feed 重建。

    构造时用冻结核的完全转化映射，与 ``freeze_q`` 相同。求解时 ``outer_P`` 再按床温的 DAEM 覆盖热解。
    灰给料和干料入口焓不另存，残差直接读 feed。
    """
    if int(nq.stage) < 4:
        return nq
    feed = nq.energy.feed
    sp = nq.energy.species
    mapped = map_one_state(
        feed,
        P_Pa=float(sp.geo_ex.P),
        T_K=float(sp.geo_ex.T),
        X_vm=np.array([0.0, 1.0], dtype=np.float64),
        tau_s=np.array([1.0], dtype=np.float64),
        occupancy_source=nq.energy.occupancy_source,
        fuel_type="coal",
    )
    n_tar = {
        "TAR1": float(mapped["cells"][0]["n_py_mol_s"].get("TAR1", 0.0)),
        "TAR2": float(mapped["cells"][0]["n_py_mol_s"].get("TAR2", 0.0)),
    }
    species = apply_pyrolysis_map(sp, mapped)
    species = replace(species, n_moist=moist_inlet_mol_s(float(feed.m_moist_kg_s)))
    return replace(nq, energy=replace(nq.energy, species=species, n_tar=n_tar))


def inner_hydro_frozen_residual(
    x: npt.NDArray[np.float64],
    nq: NamedStageQ,
) -> npt.NDArray[np.float64]:
    """内层目标：水力冻结，热解按该状态的床温更新。不是完整 ``P(x)``。"""
    return coupled_residual_named(x, refresh_pyrolysis_named(x, nq), rate_scale=1.0)


def refresh_pyrolysis_named(x: npt.NDArray[np.float64], nq: NamedStageQ) -> NamedStageQ:
    """只按床温重算 DAEM/Hegermann 热解。水力、抛射和交换保持传入的 q。"""
    z = np.asarray(x, dtype=np.float64)
    if int(nq.stage) < 4:
        return nq
    T = float(z[-2]) if int(nq.stage) >= 9 else float(z[-1])
    b = nq.reference if nq.reference is not None else capture_named_reference(nq)
    X = daem_X_at_T(T, nq.tau_daem_s)
    mapped = map_one_state(
        nq.energy.feed,
        P_Pa=float(b.P),
        T_K=T,
        X_vm=np.array([0.0, X], dtype=np.float64),
        tau_s=np.array([float(nq.tau_daem_s)], dtype=np.float64),
        occupancy_source=nq.energy.occupancy_source,
        fuel_type="coal",
    )
    n_tar = {
        "TAR1": float(mapped["cells"][0]["n_py_mol_s"].get("TAR1", 0.0)),
        "TAR2": float(mapped["cells"][0]["n_py_mol_s"].get("TAR2", 0.0)),
    }
    species = apply_pyrolysis_map(nq.energy.species, mapped)
    species = replace(species, n_moist=moist_inlet_mol_s(float(nq.energy.feed.m_moist_kg_s)))
    return replace(nq, energy=replace(nq.energy, species=species, n_tar=n_tar))


def outer_P_named(x: npt.NDArray[np.float64], nq: NamedStageQ) -> NamedStageQ:
    """外层 P(x;b)。E4 起 DAEM；E8 有限水力 + K_eject。不重写 N/T，不更新旋风 K。"""
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-2]) if int(nq.stage) >= 9 else float(z[-1])
    b = nq.reference if nq.reference is not None else capture_named_reference(nq)
    qe = nq.energy
    hydro = nq.hydro
    if int(nq.stage) >= 4:
        X = daem_X_at_T(T, nq.tau_daem_s)
        mapped = map_one_state(
            qe.feed,
            P_Pa=float(b.P),
            T_K=T,
            X_vm=np.array([0.0, X], dtype=np.float64),
            tau_s=np.array([float(nq.tau_daem_s)], dtype=np.float64),
            occupancy_source=qe.occupancy_source,
            fuel_type="coal",
        )
        n_tar = {
            "TAR1": float(mapped["cells"][0]["n_py_mol_s"].get("TAR1", 0.0)),
            "TAR2": float(mapped["cells"][0]["n_py_mol_s"].get("TAR2", 0.0)),
        }
        q2 = apply_pyrolysis_map(qe.species, mapped)
        q2 = replace(q2, n_moist=moist_inlet_mol_s(float(qe.feed.m_moist_kg_s)))
        qe = replace(qe, species=q2, n_tar=n_tar)
    else:
        qe = outer_P_e0(pack_e0(z[:N_SPECIES_BLOCK], T), qe)
    if int(nq.stage) >= 8:
        T_fb = float(z[-1]) if int(nq.stage) >= 9 else T
        N_fb = z[16:23] if int(nq.stage) >= 9 else None
        snap, geo_fb, K_bd = finite_hydro_from_x(
            z, b, T_bed_K=T, T_fb_K=T_fb, N_fb=N_fb,
        )
        hydro = snap
        sp = qe.species
        K = np.zeros(int(b.d_p.size), dtype=np.float64)
        for i, d in enumerate(b.d_p):
            geo_i = replace(geo_fb, d_p=float(d))
            K[i] = k_eject_from_geo(geo_i, b.H_bed)
        geo_ex = replace(sp.geo_ex, K_bd=float(K_bd), T=T, P=float(b.P), V_d=b.V_d, V_b=b.V_b)
        qe = replace(
            qe,
            species=replace(
                sp,
                K_eject=K,
                w_leave=snap.w_leave,
                w_return=snap.w_return,
                geo_fb=geo_fb,
                geo_ex=geo_ex,
                k=b.cyclone_k,
                D_g=snap.D_g,
                rho_g_fb_ref=float(b.rho_g_ref),
            ),
        )
    return NamedStageQ(
        energy=qe,
        stage=nq.stage,
        Q_W=nq.Q_W,
        rho_ash=nq.rho_ash,
        tau_daem_s=nq.tau_daem_s,
        enabled_rxns=nq.enabled_rxns,
        reference=b,
        hydro=hydro,
    )


def solve_inner_named(
    nq: NamedStageQ,
    *,
    x0: npt.NDArray[np.float64] | None = None,
    rate_scale: float = 1.0,
    max_iter: int = 40,
) -> dict[str, Any]:
    n = n_unknowns(nq.stage)
    x_start = default_x0_named(nq) if x0 is None else np.asarray(x0, dtype=np.float64)
    scale = np.maximum(np.abs(x_start), 0.05)
    scale[-1] = max(float(scale[-1]), 50.0)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        if float(np.min(z[:-1])) < 0.0 or float(z[-1]) <= 0.0:
            return np.full(n, 1.0e6, dtype=np.float64)
        return coupled_residual_named(z, nq, rate_scale=rate_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return fd_jacobian(
            lambda w: coupled_residual_named(w, nq, rate_scale=rate_scale),
            z,
            scale=scale,
        )

    sol = wirsum_1997_strict_newton(res, jac, x_start, max_iter=max_iter)
    stop = classify_newton_stop(sol, model_is_full_rate=abs(rate_scale - 1.0) < 1.0e-15)
    return {
        **stop,
        "x": np.asarray(sol.x, dtype=np.float64),
        "F": np.asarray(sol.F, dtype=np.float64),
        "rank_last": int(sol.rank_last),
        "T_K": float(sol.x[-1]),
        "energy_W": float(energy_row_named(sol.x, nq)),
        "f_name": f_name(nq.stage),
        "enabled_rxns": tuple(sorted(nq.enabled_rxns)),
        "rewrites_unknown_N": False,
    }


def probe_reaction_ladder(
    *,
    mode: str = "both",
    max_iter: int = 8,
) -> dict[str, Any]:
    """从 F_E0 按反应接入。``isolated`` 每次只开一项；``cumulative`` 逐项累加。

    不改 F_E1/F_E2 具名捆包，也不宣称炉温收敛。
    """
    if mode not in {"isolated", "cumulative", "both"}:
        raise ValueError(f"probe mode {mode!r} not in isolated|cumulative|both")

    def _one(enabled: frozenset[str], label: str) -> dict[str, Any]:
        stage = stage_for_enabled(enabled)
        nq = freeze_stage_q(stage, enabled_rxns=enabled)
        x0 = default_x0_named(nq)
        F0 = coupled_residual_named(x0, nq)
        i_max = int(np.argmax(np.abs(F0)))
        sol = solve_inner_named(nq, max_iter=max_iter)
        return {
            "label": label,
            "enabled_rxns": tuple(sorted(enabled)),
            "stage": stage,
            "n_unknowns": n_unknowns(stage),
            "x0_l2": float(np.linalg.norm(F0)),
            "x0_linf": float(np.max(np.abs(F0))),
            "x0_argmax": i_max,
            "physical_convergence": bool(sol["physical_convergence"]),
            "converged_newton": bool(sol["converged_newton"]),
            "failure_reason": sol["failure_reason"],
            "residual_l2": float(sol["residual_l2"]),
            "T_K": float(sol["T_K"]),
            "rank_last": int(sol["rank_last"]),
            "iterations": sol.get("iterations"),
        }

    out: dict[str, Any] = {
        "not_a_furnace_convergence_claim": True,
        "physical_convergence": False,
        "baseline": _one(frozenset(), "F_E0_R1_only"),
    }
    if mode in {"isolated", "both"}:
        out["isolated"] = [
            _one(frozenset({name}), f"only_{name}") for name in NAMED_PROBE_RXNS
        ]
    if mode in {"cumulative", "both"}:
        acc: set[str] = set()
        rows: list[dict[str, Any]] = []
        for name in NAMED_PROBE_RXNS:
            acc.add(name)
            rows.append(_one(frozenset(acc), f"upto_{name}"))
        out["cumulative"] = rows
    return out


def homogeneous_full_rate_gate() -> dict[str, Any]:
    """F_E2 前置：16 维完整速率线性求解质量。不是物理收敛通过标记。"""
    N_d, N_b = default_full_gas_inlets()
    x = np.concatenate([N_d, N_b, np.array([0.010 / 0.080, 1100.0], dtype=np.float64)])
    inert = energy_jacobian_quality_at(x, rate_scale=0.0)
    full = energy_jacobian_quality_at(x, rate_scale=1.0)
    may = bool(
        full["linear_solve_quality_checked"]
        and full["raw_rank"] == 16
        and full["raw_condition"] < 1.0e12
    )
    return {
        "inert_raw_rank": inert["raw_rank"],
        "full_raw_rank": full["raw_rank"],
        "full_raw_condition": full["raw_condition"],
        "full_equilibrated_rank": full["equilibrated_rank"],
        "linear_solve_quality_checked": full["linear_solve_quality_checked"],
        "physical_convergence": False,
        "fe2_may_claim_physical_convergence": may,
        "T_retained": True,
        "fictitious_Q_W": False,
    }


def scan_E_named(
    nq: NamedStageQ,
    T_K: npt.NDArray[np.float64],
) -> dict[str, Any]:
    """每个 T 先用 39 维物种根（含本阶段源）再评能量。"""
    from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import solve_coupled as solve_species

    T_arr = np.asarray(T_K, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    x_sp = None
    for T in T_arr:
        T = float(T)
        qT = q_at_T(nq.energy.species, T)
        spec = solve_species(qT, x0=x_sp, max_iter=40)
        x_sp = spec["x"]
        x = pack_named(
            x_sp, T, nq.stage,
            M_ash=ash_feed_kg_s(nq),
            N_tar1=float(nq.energy.n_tar.get("TAR1", 0.0)),
            N_tar2=float(nq.energy.n_tar.get("TAR2", 0.0)),
        )
        E = energy_row_named(x, nq)
        rows.append(
            {
                "T_K": T,
                "energy_W": float(E),
                "species_physical": bool(spec["physical_convergence"]),
            }
        )
    E = np.array([r["energy_W"] for r in rows], dtype=np.float64)
    return {
        "f_name": f_name(nq.stage),
        "T_K": T_arr.copy(),
        "energy_W": E,
        "rows": rows,
        "not_a_furnace_convergence_claim": True,
        "Q_W": float(nq.Q_W) if nq.stage >= 5 else 0.0,
    }


def stage_payload(stage: int, *, Q_W: float = 0.0) -> dict[str, Any]:
    nq = freeze_stage_q(stage, Q_W=Q_W)
    x = default_x0_named(nq)
    F = coupled_residual_named(x, nq)
    gate = homogeneous_full_rate_gate() if int(stage) >= 2 else None
    phys = False
    if int(stage) == 0:
        from src.solvers.wirsum_1997_coupled_f_energy import solve_inner

        inner = solve_inner(nq.energy)
        phys = bool(inner["physical_convergence"])
    elif int(stage) == 2 and gate is not None and not gate["fe2_may_claim_physical_convergence"]:
        phys = False
    return {
        "f_name": f_name(stage),
        "cancelled": cancelled_this_stage(stage),
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "n_unknowns": n_unknowns(stage),
        "F_size": int(F.size),
        "Q_W": float(nq.Q_W) if stage >= 5 else 0.0,
        "Q_W_is_fitted_to_1100K": False,
        "omitted": omitted_named(stage),
        "physical_convergence": phys,
        "homogeneous_gate": gate,
        "rewrites_unknown_N": False,
        "cyclone_K_updated": False,
        "production_hamel_outer_refresh": False,
    }
