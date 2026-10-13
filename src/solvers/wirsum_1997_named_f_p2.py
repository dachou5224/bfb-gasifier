"""P2：有限水力 ``q=P(x;b)``。不经 ``Reactor``，不更新旋风 K。

不可变参考输入 ``b`` 在 ``freeze_stage_q`` 捕获，``P`` 只读 ``b``，
不得把已更新的 ``geo_fb.rho_g`` 再当参考密度。

气速、密度由当前床层气量、组成、T、P、截面积用理想气体一致算出。
气泡速度用冻结 ``d_b`` 上的 Hilligardt Eq. 3.13；``K_bd`` 用 Hamel Eq. 3.50。
抛射命运对 Hamel Eq. 3.72 初速高斯作梯形积分；未判定不是逃逸。

Source: Hamel (1999) Eq. 3.13, 3.44, 3.50, 3.72–3.73, 3.82–3.96;
docs/audit_reports/named_f_e8_review_2026-09-16.md §2.3–2.4, §5 P2
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg, n_b
from src.core.hydrodynamics_defaults import LEGACY_HYDRO_BUBBLE_VELOCITY_STRATEGY
from src.core.species import MOLECULAR_WEIGHT, gas_diffusivity_correlation
from src.physics.bubble_dynamics import bubble_rise_velocity
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.physics.minimum_fluidization import compute_u_mf
from src.physics.wirsum_1997_freeboard import (
    FREEBOARD_MODEL_HAMEL_EXTENDED,
    MARTENS_U_P0_OVER_U_B,
    WIRSUM_1997_VELOCITY_SIGMA,
    drag_z_for_spec,
    continue_undecided_event,
    integrate_particle_fate_with_time,
    integrate_particle_through_cells,
    resolve_freeboard_model,
    u_gb_hamel_374,
)
from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry, eject_bed_top
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK
from src.solvers.wirsum_1997_stiffness_full_gas import MAJOR, N_MAJ

PSI_B_TECHNICAL_DISTRIBUTOR: float = 0.76  # [-] Hamel 工业分布板
LAUNCH_N_BINS_DEFAULT: int = 11
LAUNCH_DT_S: float = 2.0e-3  # [s]
LAUNCH_T_MAX_S: float = 5.0  # [s]
LAUNCH_N_SIGMA: float = 4.0  # [-] Eq. 3.72 截断
_G_MOL_TO_KG: float = 1.0e-3  # [kg/g]
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]


@dataclass(frozen=True)
class NamedReferenceInputs:
    """不可变参考输入 b。P(x;b) 不得改这些字段。"""

    P: float  # [Pa]
    area: float  # [m²]
    H_bed: float  # [m]
    H_freeboard: float  # [m]
    T_ref: float  # [K] 冻结时几何温度，不是拟合炉温
    rho_g_ref: float  # [kg/m³]
    rho_p: float  # [kg/m³]
    mu_g: float  # [Pa·s]
    phi_s: float  # [-]
    d_p: npt.NDArray[np.float64]  # [m]
    d_b: float  # [m] 单区床面冻结气泡直径，不是顶格中心外推
    eps_d: float  # [-]
    eps_b: float  # [-] 单区床面冻结气泡相含率，不是顶格中心外推
    eps_mf: float  # [-]
    u_mf: float  # [m/s]
    psi_b: float  # [-]
    V_d: float  # [m³]
    V_b: float  # [m³]
    cyclone_k: CycloneCellK
    trajectory_continue_s: float = 0.0  # [s] 5 s 窗口之后的物理续飞；0 表示不续
    trajectory_extra_physical_s: float = 0.0  # [s] 窗口和上一项都用尽后，还能再飞的物理时间
    trajectory_max_work_steps: int = 0  # 粗步数上限，不由秒数除以 dt 换算
    trajectory_wall_s: float = 0.0  # [s] 这段续积分的墙钟，不是飞行时间


@dataclass(frozen=True)
class NamedHydroSnapshot:
    """P2 水力快照。通量权重与驻留，不是 FB 库存未知量。"""

    u_g: float  # [m/s]
    rho_g: float  # [kg/m³]
    MW_mix: float  # [kg/mol]
    n_tot: float  # [mol/s]
    u_b: float  # [m/s] 由床层气速和冻结床面 d_b 算出，不是自由空间温度
    u_mf: float  # [m/s]
    K_bd: float  # [1/s]
    D_g: float  # [m²/s]
    w_leave: npt.NDArray[np.float64]
    w_return: npt.NDArray[np.float64]
    w_undecided: npt.NDArray[np.float64]
    tau_s: npt.NDArray[np.float64]  # [s] 已判定轨迹的 w·t 之和
    tau_return_s: npt.NDArray[np.float64]  # [s] 回落路径
    tau_leave_s: npt.NDArray[np.float64]  # [s] 越格路径
    occupancy_s: npt.NDArray[np.float64]  # [s] 未判定份额 × τ
    n_bins: int
    dt: float  # [s]
    t_max: float  # [s]
    cyclone_k_model: str
    w_to_next: npt.NDArray[np.float64]  # 第一格上边界份额。上面还有格子时不进旋风


def capture_named_reference(nq: Any) -> NamedReferenceInputs:
    """从冻结 q 抽出 b。只在 freeze 时调用。"""
    sp = nq.energy.species
    geo = sp.geo_fb
    ex = sp.geo_ex
    d_p = np.asarray(sp.d_p, dtype=np.float64).copy()
    d_hydro = float(np.exp(np.mean(np.log(np.maximum(d_p, 1.0e-12)))))
    eps_mf = float(geo.eps_d)
    u_mf = float(
        compute_u_mf(
            rho_g=float(geo.rho_g),
            rho_s=float(geo.rho_p),
            d_p=d_hydro,
            mu_g=float(geo.mu_g),
            eps_mf=eps_mf,
            phi_s=float(geo.phi_s),
        )
    )
    return NamedReferenceInputs(
        P=float(ex.P),
        area=float(geo.area),
        H_bed=float(sp.H_bed),
        H_freeboard=float(geo.H_freeboard),
        T_ref=float(ex.T),
        rho_g_ref=float(sp.rho_g_fb_ref),
        rho_p=float(geo.rho_p),
        mu_g=float(geo.mu_g),
        phi_s=float(geo.phi_s),
        d_p=d_p,
        d_b=float(geo.d_b),
        eps_d=float(geo.eps_d),
        eps_b=float(geo.eps_b),
        eps_mf=eps_mf,
        u_mf=u_mf,
        psi_b=PSI_B_TECHNICAL_DISTRIBUTOR,
        V_d=float(ex.V_d),
        V_b=float(ex.V_b),
        cyclone_k=sp.k,
    )


def mix_molar_mass_kg_per_mol(N_d: npt.NDArray[np.float64], N_b: npt.NDArray[np.float64]) -> float:
    """床层总气的平均摩尔质量 [kg/mol]。"""
    n = np.asarray(N_d, dtype=np.float64) + np.asarray(N_b, dtype=np.float64)
    assert n.shape == (N_MAJ,)
    n_tot = float(np.sum(n))
    if n_tot < _NTOT_FLOOR:
        return 0.0
    acc = 0.0
    for name, ni in zip(MAJOR, n):
        acc += float(ni) * float(MOLECULAR_WEIGHT[name]) * _G_MOL_TO_KG
    return acc / n_tot


def gas_velocity_and_density(
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    *,
    T: float,
    P: float,
    area: float,
) -> tuple[float, float, float, float]:
    """理想气体：``(u_g [m/s], rho_g [kg/m³], n_tot [mol/s], MW [kg/mol])``。

    ``u_g = ṅ Rg T / (P A)``，``ρ = (P MW)/(Rg T)``。
    """
    n_tot = float(np.sum(N_d) + np.sum(N_b))
    MW = mix_molar_mass_kg_per_mol(N_d, N_b)
    T_use = max(float(T), 1.0)
    P_use = float(P)
    A = max(float(area), 1.0e-12)
    u_g = (n_tot * Rg * T_use) / (P_use * A) if n_tot > 0.0 else 0.0
    rho_g = (P_use * MW) / (Rg * T_use) if MW > 0.0 else 0.0
    return float(u_g), float(rho_g), n_tot, MW


def eq372_launch_nodes(
    u_mean: float,
    sigma: float,
    n_bins: int,
    *,
    n_sigma: float = LAUNCH_N_SIGMA,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Hamel Eq. 3.72 高斯密度在 ``u≥0`` 上的梯形节点与归一权重。

    不是 ``truncated_pm_2sigma_uniform_xi`` 工程离散。
    """
    u0 = float(u_mean)
    sig = max(float(sigma), 0.0)
    n = int(n_bins)
    if n < 1:
        raise ValueError("n_bins must be >= 1")
    if n == 1 or sig <= 0.0:
        return (
            np.array([max(u0, 0.0)], dtype=np.float64),
            np.array([1.0], dtype=np.float64),
        )
    u_lo = max(0.0, u0 - float(n_sigma) * sig)
    u_hi = u0 + float(n_sigma) * sig
    u = np.linspace(u_lo, u_hi, n, dtype=np.float64)
    xi = (u - u0) / sig
    w = np.exp(-0.5 * xi * xi)
    du = float(u[1] - u[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    w = w * du
    s = float(np.sum(w))
    if s <= 0.0:
        w = np.full(n, 1.0 / n, dtype=np.float64)
    else:
        w = w / s
    return u, w


def integrate_launch_distribution(
    geo: BedFreeboardGeometry,
    *,
    n_bins: int = LAUNCH_N_BINS_DEFAULT,
    dt: float = LAUNCH_DT_S,
    t_max: float = LAUNCH_T_MAX_S,
    gas_profile: str = "superficial",
    n_cells: int = 1,
    continue_budget_s: float = 0.0,
    extra_physical_s: float = 0.0,
    max_work_steps: int = 0,
    wall_s: float = 0.0,
    runtime_cap_s: float | None = None,
    resolve_events: bool = False,
) -> dict[str, float]:
    """初速分布积分：离开 / 回落 / 未判定份额与平均驻留。

    ``superficial`` 是这一格里已经掺混后的表观气速，全程取 ``geo.u_g``。
    具名核库存和缺省抛射走这条。``hamel_374`` 是同一格里破裂气泡的
    ``u_gb(h)``，不是缺省场。它的未判定份额只留在 ``w_undecided`` 和
    ``occupancy_s``，不写入 ``tau_s``，也不当成出炉或无限驻留。
    额外物理时间、粗步数和墙钟分开传。未到事件的时间不进入 ``tau_s``。
    """
    spec = resolve_freeboard_model(FREEBOARD_MODEL_HAMEL_EXTENDED)
    u_mean = MARTENS_U_P0_OVER_U_B * float(geo.u_b)
    u_nodes, w_nodes = eq372_launch_nodes(
        u_mean, WIRSUM_1997_VELOCITY_SIGMA, int(n_bins)
    )
    z_drag = drag_z_for_spec(spec, phi_s=geo.phi_s)
    u_of_h = None
    if gas_profile == "hamel_374":
        u_of_h = lambda h: u_gb_hamel_374(h, float(geo.u_g), float(geo.u_b), float(geo.d_b))
    elif gas_profile != "superficial":
        raise ValueError(f"unknown gas_profile {gas_profile}")
    w_leave = 0.0
    w_ret = 0.0
    w_und = 0.0
    tau_ret = 0.0
    tau_leave = 0.0
    occ = 0.0
    traces: list[dict[str, Any]] = []
    n = max(int(n_cells), 1)
    h_each = float(geo.H_freeboard) / n
    heights = [h_each] * n
    for u_p0, wi in zip(u_nodes, w_nodes):
        if n == 1:
            fate, t_ev, _z, _u = integrate_particle_fate_with_time(
                float(u_p0),
                u_g=float(geo.u_g),
                rho_g=float(geo.rho_g),
                rho_p=float(geo.rho_p),
                d_p=float(geo.d_p),
                mu_g=float(geo.mu_g),
                z=z_drag,
                cell_height=float(geo.H_freeboard),
                dt=float(dt),
                t_max=float(t_max),
                u_g_of_h=u_of_h,
                resolve_events=bool(resolve_events) and u_of_h is None,
            )
            if fate == "undecided" and float(continue_budget_s) > 0.0 and not resolve_events:
                fate2, t_add, _z2, _u2, ok, _kind = continue_undecided_event(
                    _u, _z,
                    u_g=float(geo.u_g), rho_g=float(geo.rho_g), rho_p=float(geo.rho_p),
                    d_p=float(geo.d_p), mu_g=float(geo.mu_g), z=z_drag,
                    cell_height=float(geo.H_freeboard), dt=float(dt),
                    extra_s=float(continue_budget_s), u_g_of_h=u_of_h,
                )
                if ok:
                    fate = fate2
                    t_ev = float(t_ev) + float(t_add)
        else:
            fate, t_ev, _cell = integrate_particle_through_cells(
                float(u_p0),
                heights,
                u_g=float(geo.u_g),
                rho_g=float(geo.rho_g),
                rho_p=float(geo.rho_p),
                d_p=float(geo.d_p),
                mu_g=float(geo.mu_g),
                z=z_drag,
                dt=float(dt),
                t_max=float(t_max),
                u_g_of_h=u_of_h,
                continue_budget_s=float(continue_budget_s),
                trace=traces,
                extra_physical_s=float(extra_physical_s),
                max_work_steps=int(max_work_steps),
                wall_s=float(wall_s),
                runtime_cap_s=runtime_cap_s,
                resolve_events=bool(resolve_events) and u_of_h is None,
            )
        if n > 1 and traces:
            traces[-1]["u_launch_m_s"] = float(u_p0)
            traces[-1]["weight"] = float(wi)
        if fate == "reaches_next_cell":
            w_leave += float(wi)
            tau_leave += float(wi) * float(t_ev)
        elif fate == "returns":
            w_ret += float(wi)
            tau_ret += float(wi) * float(t_ev)
        else:
            w_und += float(wi)
            occ += float(wi) * float(t_ev)
    return {
        "w_leave": float(w_leave),
        "w_return": float(w_ret),
        "w_undecided": float(w_und),
        "tau_s": float(tau_ret + tau_leave),
        "tau_return_s": float(tau_ret),
        "tau_leave_s": float(tau_leave),
        "occupancy_s": float(occ),
        "traces": traces,
    }


def bed_face_char_interface(
    geo: BedFreeboardGeometry,
    m_char_kg: npt.NDArray[np.float64],
    d_p_m: npt.NDArray[np.float64],
    H_bed: float,
    H_cell: float,
    *,
    n_bins: int = LAUNCH_N_BINS_DEFAULT,
    dt: float = LAUNCH_DT_S,
    t_max: float = LAUNCH_T_MAX_S,
    gas_profile: str = "superficial",
    is_roof: bool = True,
) -> dict[str, npt.NDArray[np.float64]]:
    """床面抛出只算一次，同一股流进入床顶出流和自由空间这一格的入流。

    ``H_cell`` 只决定格内轨迹和库存，不改变抛出流。越过该格上边界的流量是
    轨迹结果，不是第二份抛射。``is_roof`` 为假时，这份流量进入上一格，
    出炉量为零。未判定时间不进入库存。几何床料分母把炭通量从总床料通量里分开。
    """
    m_char = np.asarray(m_char_kg, dtype=np.float64).reshape(-1)
    d_p = np.asarray(d_p_m, dtype=np.float64).reshape(-1)
    if int(m_char.size) != int(d_p.size):
        raise ValueError("m_char and d_p length mismatch")
    n = int(m_char.size)
    m_eject = np.zeros(n, dtype=np.float64)
    m_return = np.zeros(n, dtype=np.float64)
    m_cross = np.zeros(n, dtype=np.float64)
    m_cell = np.zeros(n, dtype=np.float64)
    w_und = np.zeros(n, dtype=np.float64)
    for i in range(n):
        geo_i = replace(geo, d_p=float(d_p[i]))
        k = k_eject_from_geo(geo_i, H_bed)
        m_eject[i] = k * float(m_char[i])
        dist = integrate_launch_distribution(
            replace(geo_i, H_freeboard=float(H_cell)),
            n_bins=n_bins,
            dt=dt,
            t_max=t_max,
            gas_profile=gas_profile,
        )
        m_return[i] = float(dist["w_return"]) * m_eject[i]
        m_cross[i] = float(dist["w_leave"]) * m_eject[i]
        m_cell[i] = m_eject[i] * float(dist["tau_s"])
        w_und[i] = float(dist["w_undecided"])
    if is_roof:
        m_to_next = np.zeros(n, dtype=np.float64)
        m_exit = m_cross.copy()
    else:
        m_to_next = m_cross.copy()
        m_exit = np.zeros(n, dtype=np.float64)
    return {
        "m_eject": m_eject,
        "m_return": m_return,
        "m_cross": m_cross,
        "m_to_next": m_to_next,
        "m_furnace_exit": m_exit,
        "m_cell": m_cell,
        "w_undecided": w_und,
    }


def bed_solids_ejection_split(
    geo: BedFreeboardGeometry,
    m_char_kg: float,
    m_ash_kg: float,
    H_bed: float,
) -> dict[str, float]:
    """几何床料分母下的总抛射能力与炭、灰通量 [kg/s]。

    ``M_geom`` 来自 ``ρ_p(1-ε_d)(1-ε_b)A H``，不是未知量里的惰性床料。
    炭加灰之外的差额没有物种、粒径和返回路径。它的抛射通量记为 0，
    不能把它算成炭，也不能把它当成已经建账的惰性床料。
    炭通量是 ``K·M_char``，不是把 ``A F_0`` 全部派给炭。
    """
    k = k_eject_from_geo(geo, H_bed)
    m_geom = (
        float(geo.rho_p)
        * max(1.0 - float(geo.eps_d), 0.0)
        * max(1.0 - float(geo.eps_b), 0.0)
        * float(geo.area)
        * float(H_bed)
    )
    tracked = float(m_char_kg) + float(m_ash_kg)
    gap = float(m_geom) - tracked
    return {
        "M_geom_kg": float(m_geom),
        "M_char_kg": float(m_char_kg),
        "M_ash_kg": float(m_ash_kg),
        "geometric_gap_kg": gap,
        "gap_is_explicit_species": False,
        "inert_kg_s": 0.0,
        "K_eject_per_s": float(k),
        "total_capacity_kg_s": float(k) * float(m_geom),
        "char_kg_s": float(k) * float(m_char_kg),
        "ash_kg_s": float(k) * float(m_ash_kg),
    }


def k_eject_from_geo(geo: BedFreeboardGeometry, H_bed: float) -> float:
    """``K_eject = m_eject(M=1)/M_geom`` [1/s]。命运权重另算。"""
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
    return float(m_eject[0]) / max(M_geom, 1.0e-12)


def bed_face_exchange(
    x: npt.NDArray[np.float64],
    b: NamedReferenceInputs,
    *,
    T_bed_K: float,
) -> tuple[npt.NDArray[np.float64], float, BedFreeboardGeometry, float, float, float, float, float, float]:
    """床面 ``K_eject`` [1/s]、``K_bd`` [1/s] 和床面几何。不积分颗粒轨迹。

    Source: Hamel (1999) Eq. 3.13, Eq. 3.50；与 ``finite_hydro_from_x`` 的床面段相同。
    """
    z = np.asarray(x, dtype=np.float64)
    T = float(T_bed_K)  # [K]
    N_d = z[0:N_MAJ]
    N_b = z[7:14]
    u_g, rho_g, n_tot, MW = gas_velocity_and_density(
        N_d, N_b, T=T, P=b.P, area=b.area
    )
    u_b = float(
        bubble_rise_velocity(
            u_g,
            b.u_mf,
            b.d_b,
            b.psi_b,
            strategy=LEGACY_HYDRO_BUBBLE_VELOCITY_STRATEGY,
        )
    )
    u_d = b.u_mf / max(b.eps_mf, 1.0e-12)  # [m/s]
    u_br = calc_u_br(u_d, b.P, n_b_factor=n_b)
    D_g = float(gas_diffusivity_correlation(T, b.P))  # [m²/s]
    K_bd = float(calc_kbd(u_br, b.d_b, D_g, b.eps_mf, u_b, model="hamel_350"))  # [1/s]
    geo_base = BedFreeboardGeometry(
        area=b.area,
        H_freeboard=b.H_freeboard,
        d_p=float(b.d_p[0]),
        rho_p=b.rho_p,
        eps_d=b.eps_d,
        eps_b=b.eps_b,
        u_b=u_b,
        d_b=b.d_b,
        u_g=u_g,
        rho_g=rho_g,
        mu_g=b.mu_g,
        phi_s=b.phi_s,
    )
    K_eject = np.zeros(int(b.d_p.size), dtype=np.float64)  # [1/s]
    for i, d in enumerate(b.d_p):
        K_eject[i] = k_eject_from_geo(replace(geo_base, d_p=float(d)), b.H_bed)
    return K_eject, K_bd, geo_base, u_g, rho_g, n_tot, MW, u_b, D_g


def finite_hydro_from_x(
    x: npt.NDArray[np.float64],
    b: NamedReferenceInputs,
    *,
    T_bed_K: float | None = None,
    T_fb_K: float | None = None,
    N_fb: npt.NDArray[np.float64] | None = None,
    n_bins: int = LAUNCH_N_BINS_DEFAULT,
    dt: float = LAUNCH_DT_S,
    t_max: float = LAUNCH_T_MAX_S,
    H_cell: float | None = None,
    cells_above: int = 0,
    launch_traces: list | None = None,
    runtime_cap_s: float | None = None,
) -> tuple[NamedHydroSnapshot, BedFreeboardGeometry, float]:
    """床层水力用床温与床层气量。返回的 ``geo_fb`` 仍是床面状态，供抛射使用。

    若给出自由空间温度和 ``N_fb``，颗粒轨迹改用该区的表观气速和密度。
    不把式 (3.74) 的气泡破裂气流剖面放进这条路径。

    ``H_cell`` 是单独指定的一格高度 [m]。缺省时整段自由空间分成两格串联：
    回床只计最底格下边界，出炉只计最顶格上边界。``cells_above > 0`` 时，
    调用者给出的这一格不是炉顶，越格份额记在 ``w_to_next``，``w_leave`` 置零。
    均匀气速的轨迹用解析边界事件，不再把 30/60/120 s 续飞当成物理边界。
    渐近滞留仍留在 ``w_undecided``，不写成停留时间。
    """
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-1] if T_bed_K is None else T_bed_K)
    _K_eject, K_bd, geo_base, u_g, rho_g, n_tot, MW, u_b, D_g = bed_face_exchange(
        z, b, T_bed_K=T,
    )
    n_bin = int(b.d_p.size)
    wL = np.zeros(n_bin, dtype=np.float64)
    wR = np.zeros(n_bin, dtype=np.float64)
    wU = np.zeros(n_bin, dtype=np.float64)
    tau = np.zeros(n_bin, dtype=np.float64)
    tau_ret = np.zeros(n_bin, dtype=np.float64)
    tau_leave = np.zeros(n_bin, dtype=np.float64)
    occ = np.zeros(n_bin, dtype=np.float64)
    if N_fb is not None and T_fb_K is not None:
        u_g_fb, rho_g_fb, _, _ = gas_velocity_and_density(
            np.asarray(N_fb, dtype=np.float64),
            np.zeros(N_MAJ, dtype=np.float64),
            T=float(T_fb_K),
            P=b.P,
            area=b.area,
        )
        geo_traj = replace(geo_base, u_g=u_g_fb, rho_g=rho_g_fb)
    else:
        geo_traj = geo_base
    h_cell = float(b.H_freeboard if H_cell is None else H_cell)
    geo_traj = replace(geo_traj, H_freeboard=h_cell)
    w_next = np.zeros(n_bin, dtype=np.float64)
    roof = int(cells_above) <= 0
    for i, d in enumerate(b.d_p):
        geo_i = replace(geo_traj, d_p=float(d))
        continue_s = float(b.trajectory_continue_s)
        dist = integrate_launch_distribution(
            geo_i,
            n_bins=n_bins,
            dt=dt,
            t_max=t_max,
            n_cells=1 if H_cell is not None else 2,
            continue_budget_s=continue_s,
            extra_physical_s=float(b.trajectory_extra_physical_s),
            max_work_steps=int(b.trajectory_max_work_steps),
            wall_s=float(b.trajectory_wall_s),
            runtime_cap_s=runtime_cap_s,
            resolve_events=True,
        )
        if launch_traces is not None:
            for item in dist["traces"]:
                row = dict(item)
                row["d_p_m"] = float(d)
                launch_traces.append(row)
        crossed = float(dist["w_leave"])
        if roof:
            wL[i] = crossed
        else:
            w_next[i] = crossed
        wR[i] = dist["w_return"]
        wU[i] = dist["w_undecided"]
        tau[i] = dist["tau_s"]
        tau_ret[i] = dist["tau_return_s"]
        tau_leave[i] = dist["tau_leave_s"]
        occ[i] = dist["occupancy_s"]
    snap = NamedHydroSnapshot(
        u_g=u_g,
        rho_g=rho_g,
        MW_mix=MW,
        n_tot=n_tot,
        u_b=u_b,
        u_mf=b.u_mf,
        K_bd=K_bd,
        D_g=D_g,
        w_leave=wL,
        w_return=wR,
        w_undecided=wU,
        tau_s=tau,
        tau_return_s=tau_ret,
        tau_leave_s=tau_leave,
        occupancy_s=occ,
        n_bins=int(n_bins),
        dt=float(dt),
        t_max=float(t_max),
        cyclone_k_model=str(b.cyclone_k.model),
        w_to_next=w_next,
    )
    return snap, geo_base, K_bd
