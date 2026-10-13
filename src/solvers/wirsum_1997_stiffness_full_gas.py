"""阶段 E 后半：等温全气相 → 炭 → 能量。不经 ``Reactor``。

气相用生产 ``rate_R5_*`` / ``R6`` / ``R7`` / ``R8`` / ``R12``。
炭用生产 ``rate_R1``–``R4`` / ``extent_R4``。能量用 ``enthalpy_molar``。
``rate_scale<1`` 是延拓。雅可比按未知量尺度差分，微量氧用绝对下限。

Source: Hamel (1999) Eq. 2.2, 5.8–5.51; docs/audit_reports/secondary_diagnosis_2026-09-14.md §9 E
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg
from src.core.species import (
    GAS_SPECIES_INDEX,
    MW_C_KG_PER_MOL,
    N_GAS,
    SOLID_H_MODE_INTEGRAL,
    cp_molar,
    enthalpy_molar,
    h_solid_mass,
)
from src.kinetics.char_reactions import extent_R4, rate_R1, rate_R2, rate_R3
from src.kinetics.gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R12,
)
from src.solvers.wirsum_1997_stiffness_ladder import (
    ExchangeGeometry,
    classify_newton_stop,
    newton_step_backward_error,
)
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

MAJOR: tuple[str, ...] = ("CO", "CO2", "H2", "H2O", "CH4", "O2", "N2")
N_MAJ: int = 7
IDX: dict[str, int] = {name: i for i, name in enumerate(MAJOR)}
# 仓库 R6 = 论文 R7。Hamel Eq. 5.40：CH4 + 1.5 O2 → CO + 2 H2O
NU_R6: npt.NDArray[np.float64] = np.array(
    [1.0, 0.0, 0.0, 2.0, -1.0, -1.5, 0.0], dtype=np.float64
)
_ELEM = {"C": np.array([1, 1, 0, 0, 1, 0, 0], dtype=np.float64),
         "H": np.array([0, 0, 2, 2, 4, 0, 0], dtype=np.float64),
         "O": np.array([1, 2, 0, 1, 0, 2, 0], dtype=np.float64)}
_NTOT_FLOOR: float = 1.0e-30  # [mol/s]
_FD_REL: float = 1.0e-8
_FD_ABS: float = 1.0e-10  # [与 x 同单位的绝对扰动下限]


def _y_full(N: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    y = np.zeros(N_GAS, dtype=np.float64)
    ntot = float(np.sum(N))
    if abs(ntot) < _NTOT_FLOOR:
        return y
    for name, val in zip(MAJOR, N):
        y[GAS_SPECIES_INDEX[name]] = float(val) / ntot
    return y


def _conc(N: npt.NDArray[np.float64], geo: ExchangeGeometry) -> npt.NDArray[np.float64]:
    ntot = float(np.sum(N))
    if abs(ntot) < _NTOT_FLOOR:
        return np.zeros(N_MAJ, dtype=np.float64)
    return (np.asarray(N, dtype=np.float64) / ntot) * (float(geo.P) / (Rg * float(geo.T)))


def fd_jacobian(
    residual_fn: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    x: npt.NDArray[np.float64],
    *,
    scale: npt.NDArray[np.float64] | None = None,
) -> npt.NDArray[np.float64]:
    """按 ``max(ε_abs, ε_rel |x_i|, ε_rel scale_i)`` 扰动。"""
    x = np.asarray(x, dtype=np.float64)
    sc = np.ones_like(x) if scale is None else np.asarray(scale, dtype=np.float64)
    F0 = residual_fn(x)
    J = np.zeros((x.size, x.size), dtype=np.float64)
    for i in range(x.size):
        eps = max(_FD_ABS, _FD_REL * max(abs(float(x[i])), abs(float(sc[i]))))
        xp = x.copy()
        xp[i] = float(x[i]) + eps
        J[:, i] = (residual_fn(xp) - F0) / eps
    return J


N_ENERGY: int = 2 * N_MAJ + 2  # 14 气体 + M + T
# 明示尺度，不是隐式 clip。气体 1 mol/s，库存 1 kg，温度 1000 K。
ENERGY_X_SCALE: npt.NDArray[np.float64] = np.array(
    [1.0] * (2 * N_MAJ) + [1.0, 1000.0], dtype=np.float64
)  # [mol/s]*14, [kg], [K]
# 残差：气体 mol/s，炭衡算 kg/s，能量 W。
ENERGY_F_SCALE: npt.NDArray[np.float64] = np.array(
    [1.0] * (2 * N_MAJ) + [1.0e-3, 1.0e5], dtype=np.float64
)  # [mol/s]*14, [kg/s], [W]


def row_col_equilibrate(
    A: npt.NDArray[np.float64],
    *,
    n_pass: int = 12,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """行列均衡诊断。不改变线搜索 merit。"""
    a = np.asarray(A, dtype=np.float64).copy()
    n = a.shape[0]
    dr = np.ones(n, dtype=np.float64)
    dc = np.ones(n, dtype=np.float64)
    for _ in range(int(n_pass)):
        r = 1.0 / np.maximum(np.max(np.abs(a), axis=1), 1.0e-300)
        a = r[:, None] * a
        dr *= r
        c = 1.0 / np.maximum(np.max(np.abs(a), axis=0), 1.0e-300)
        a = a * c
        dc *= c
    return a, dr, dc


def scale_linear_system(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    *,
    x_scale: npt.NDArray[np.float64],
    F_scale: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """``Ĵ = D_F^{-1} J D_x``，``F̂ = D_F^{-1} F``。尺度必须具名。"""
    xs = np.asarray(x_scale, dtype=np.float64)
    fs = np.asarray(F_scale, dtype=np.float64)
    J_hat = (np.asarray(J, dtype=np.float64) * xs[None, :]) / fs[:, None]
    F_hat = np.asarray(F, dtype=np.float64) / fs
    return J_hat, F_hat


def solve_scaled_newton_step(
    J: npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    *,
    x_scale: npt.NDArray[np.float64],
    F_scale: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], float]:
    """用明示尺度解 ``J Δx = -F``，返回 ``Δx`` 与原始单位后向误差。"""
    J_hat, F_hat = scale_linear_system(J, F, x_scale=x_scale, F_scale=F_scale)
    z = np.linalg.solve(J_hat, -F_hat)
    dx = np.asarray(x_scale, dtype=np.float64) * z
    return dx, float(newton_step_backward_error(J, F, dx))


def analytic_gas_dE_dT_W_per_K(
    N_d: npt.NDArray[np.float64],
    N_b: npt.NDArray[np.float64],
    T: float,
) -> float:
    """冻结 N 时 ``∂E/∂T = -∑_j (N_{d,j}+N_{b,j}) c_{p,j}(T)`` [W/K]。

    Source: Hamel (1999) Eq. 2.7–2.8 气体焓流；固体显热默认 omitted_named。
    """
    acc = 0.0  # [W/K]
    for name, nd, nb in zip(MAJOR, N_d, N_b):
        acc += (float(nd) + float(nb)) * float(cp_molar(name, T))
    return -acc


def energy_jacobian_quality_at(
    x: npt.NDArray[np.float64],
    *,
    rate_scale: float = 1.0,
    geo: ExchangeGeometry | None = None,
    N_in_d: npt.NDArray[np.float64] | None = None,
    N_in_b: npt.NDArray[np.float64] | None = None,
    T_in: float | None = None,
    m_char_in: float = 0.010,
    K_char: float = 0.080,
) -> dict[str, Any]:
    """16 维能量核的尺度/秩/线性求解诊断。不是物理收敛通过标记。

    不删 T，不补虚构 Q_W。
    Source: docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q5
    """
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_ENERGY,)
    geo = ExchangeGeometry() if geo is None else geo
    if N_in_d is None or N_in_b is None:
        N_in_d, N_in_b = default_full_gas_inlets()
    T_in = float(geo.T) if T_in is None else float(T_in)

    def residual_fn(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return energy_residual(
            v, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, T_in=T_in,
            m_char_in=m_char_in, K_char=K_char, rate_scale=rate_scale,
        )

    F = residual_fn(z)
    J = fd_jacobian(residual_fn, z, scale=ENERGY_X_SCALE)
    raw_rank = int(np.linalg.matrix_rank(J))
    raw_cond = float(np.linalg.cond(J))
    J_eq, dr, dc = row_col_equilibrate(J)
    eq_rank = int(np.linalg.matrix_rank(J_eq))
    s = np.linalg.svd(J_eq, compute_uv=False)
    eq_cond = float(s[0] / s[-1]) if s[-1] > 0.0 else float("inf")
    J_hat, F_hat = scale_linear_system(
        J, F, x_scale=ENERGY_X_SCALE, F_scale=ENERGY_F_SCALE
    )
    scaled_rank = int(np.linalg.matrix_rank(J_hat))
    scaled_cond = float(np.linalg.cond(J_hat))
    lin_ok = False
    berr_raw = float("nan")
    berr_scaled = float("nan")
    dx_raw: npt.NDArray[np.float64] | None
    try:
        dx_raw = np.linalg.solve(J, -F)
        berr_raw = float(newton_step_backward_error(J, F, dx_raw))
        lin_ok = True
    except np.linalg.LinAlgError:
        dx_raw = None
    try:
        _dx_s, berr_scaled = solve_scaled_newton_step(
            J, F, x_scale=ENERGY_X_SCALE, F_scale=ENERGY_F_SCALE
        )
        lin_ok = True
    except np.linalg.LinAlgError:
        pass
    N_d, N_b = z[:N_MAJ], z[N_MAJ:2 * N_MAJ]
    T = float(z[-1])  # [K]
    dE_an = analytic_gas_dE_dT_W_per_K(N_d, N_b, T)
    dE_fd = float(J[-1, -1])
    return {
        "n_unknowns": N_ENERGY,
        "T_retained": True,
        "fictitious_Q_W": False,
        "physical_convergence": False,
        "rate_scale": float(rate_scale),
        "merit_kind": "raw_l2",
        "linear_solve_scaling": "explicit_ENERGY_X_SCALE_ENERGY_F_SCALE",
        "raw_rank": raw_rank,
        "equilibrated_rank": eq_rank,
        "explicit_scale_rank": scaled_rank,
        "raw_condition": raw_cond,
        "equilibrated_condition": eq_cond,
        "explicit_scale_condition": scaled_cond,
        "raw_linear_solve_ok": dx_raw is not None,
        "newton_step_backward_error_raw": berr_raw,
        "newton_step_backward_error_scaled": berr_scaled,
        "linear_solve_quality_checked": lin_ok,
        "analytic_dE_dT_W_per_K": dE_an,
        "fd_dE_dT_W_per_K": dE_fd,
        "energy_row_M_W_per_kg": float(J[-1, 2 * N_MAJ]),
        "residual_l2": float(np.linalg.norm(F)),
        "x_scale_units": ["mol/s"] * 14 + ["kg", "K"],
        "F_scale_units": ["mol/s"] * 14 + ["kg/s", "W"],
        "row_scale": dr.tolist(),
        "col_scale": dc.tolist(),
        "F_hat_l2": float(np.linalg.norm(F_hat)),
    }


def default_full_gas_inlets() -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """dense 燃料气，bubble 空气。单位 [mol/s]。"""
    N_d = np.array([0.40, 0.10, 0.30, 0.20, 0.10, 0.0, 8.0], dtype=np.float64)
    N_b = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.80, 3.0], dtype=np.float64)
    return N_d, N_b


HOMOG_RXN_NAMES: tuple[str, ...] = ("R5", "R6", "R7", "R8", "R12")
CHAR_RXN_NAMES: tuple[str, ...] = ("R1", "R2", "R3", "R4")


def _gas_source(
    C: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    V: float,
    phase: str,
    rate_scale: float,
    rho_ash: float,
    enabled: frozenset[str] | None = None,
) -> npt.NDArray[np.float64]:
    """均相源项 [mol/s]。``phase`` 为 ``dense`` 或 ``bubble``。

    ``enabled`` 缺省为 R5/R6/R7/R8/R12 全开；具名 F 探针可只开其中一项。
    """
    on = frozenset(HOMOG_RXN_NAMES) if enabled is None else frozenset(enabled)
    T, P = float(geo.T), float(geo.P)
    i = IDX
    r5 = 0.0
    if "R5" in on:
        r5 = (
            rate_R5_suspension(T, C[i["CO"]], C[i["O2"]], C[i["H2O"]], P, y, GAS_SPECIES_INDEX)
            if phase == "dense"
            else rate_R5_bubble(T, C[i["CO"]], C[i["O2"]], C[i["H2O"]], P, y, GAS_SPECIES_INDEX)
        )
    r6 = rate_R6(T, C[i["CH4"]], C[i["O2"]]) if "R6" in on else 0.0
    # 气泡无固体：R7/R8 灰催化只进 dense。Ref: Hamel Eq. 5.51、Eq. 2.6
    ash = float(rho_ash) if phase == "dense" else 0.0
    r7 = (
        rate_R7(
            T, C[i["CH4"]], C[i["H2O"]], C[i["CO"]], C[i["H2"]],
            rho_ash_kg_m3=ash,
        )
        if "R7" in on else 0.0
    )
    r8 = 0.0
    if "R8" in on:
        r8 = rate_R8(
            T, P, y[GAS_SPECIES_INDEX["CO"]], y[GAS_SPECIES_INDEX["H2O"]],
            y[GAS_SPECIES_INDEX["CO2"]], y[GAS_SPECIES_INDEX["H2"]], ash,
        )
    r12 = (
        rate_R12(T, C[i["H2"]], C[i["O2"]], P, y, GAS_SPECIES_INDEX)
        if "R12" in on else 0.0
    )
    s = float(rate_scale) * float(V)
    src = np.zeros(N_MAJ, dtype=np.float64)
    src[i["CO"]] += s * (-r5 + r7 - r8)
    src[i["CO2"]] += s * (r5 + r8)
    src[i["H2"]] += s * (3.0 * r7 + r8 - r12)
    src[i["H2O"]] += s * (-r7 - r8 + r12)
    src[i["CH4"]] += s * (-r7)
    src[i["O2"]] += s * (-0.5 * r5 - 0.5 * r12)
    src += s * float(r6) * NU_R6
    return src


gas_homogeneous_source = _gas_source


def full_gas_residual(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
    rate_scale: float = 1.0,
    rho_ash: float = 75.0,
) -> npt.NDArray[np.float64]:
    """14 元摩尔流残差 [mol/s]。``rho_ash`` [kg/m³] 是工程灰密度，不是持料反演。"""
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (2 * N_MAJ,)
    N_d, N_b = z[:N_MAJ], z[N_MAJ:]
    C_d, C_b = _conc(N_d, geo), _conc(N_b, geo)
    y_d, y_b = _y_full(N_d), _y_full(N_b)
    N_ex = float(geo.K_bd) * float(geo.V_b) * (C_b - C_d)
    src_d = _gas_source(C_d, y_d, geo=geo, V=geo.V_d, phase="dense",
                        rate_scale=rate_scale, rho_ash=rho_ash)
    src_b = _gas_source(C_b, y_b, geo=geo, V=geo.V_b, phase="bubble",
                        rate_scale=rate_scale, rho_ash=rho_ash)
    return np.concatenate([N_in_d - N_d + N_ex + src_d, N_in_b - N_b - N_ex + src_b])


def element_atom_flows(N_d: npt.NDArray[np.float64], N_b: npt.NDArray[np.float64]) -> dict[str, float]:
    """气相元素摩尔流 [mol/s]。"""
    tot = np.asarray(N_d, dtype=np.float64) + np.asarray(N_b, dtype=np.float64)
    return {el: float(_ELEM[el] @ tot) for el in _ELEM}


def char_area_m2(M_kg: float, *, rho_s: float, d_p: float) -> float:
    """球颗粒外表面积 [m²]：``6 M / (ρ d_p)``。"""
    return 6.0 * max(float(M_kg), 0.0) / max(float(rho_s) * float(d_p), 1.0e-30)


def _char_gas_source(
    C: npt.NDArray[np.float64],
    N: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    M_kg: float,
    d_p: float,
    D_g: float,
    rho_s: float,
    rate_scale: float,
    enabled: frozenset[str] | None = None,
) -> tuple[npt.NDArray[np.float64], float]:
    """异相源项 [mol/s] 与炭消耗 [mol C/s]。

    ``enabled`` 缺省为 R1–R4 全开。具名 F 的 R2/R3/R4 探针不传 R1
    （R1 已在 F_E0 物种残差里）。
    """
    on = frozenset(CHAR_RXN_NAMES) if enabled is None else frozenset(enabled)
    T, P = float(geo.T), float(geo.P)
    i = IDX
    area = char_area_m2(M_kg, rho_s=rho_s, d_p=d_p)
    e1 = 0.0
    alpha = 1.0
    if "R1" in on:
        r1, alpha = rate_R1(T, float(C[i["O2"]]), d_p, D_g)
        e1 = float(rate_scale) * float(r1) * area
    e2 = 0.0
    if "R2" in on:
        e2 = float(rate_scale) * float(rate_R2(T, float(C[i["H2O"]]), d_p, D_g)) * area
    e3 = 0.0
    if "R3" in on:
        e3 = float(rate_scale) * float(rate_R3(T, float(C[i["H2"]]), d_p, D_g)) * area
    e4 = 0.0
    if "R4" in on:
        ntot = float(np.sum(N))
        p_co2 = (float(N[i["CO2"]]) / ntot) * P if abs(ntot) >= _NTOT_FLOOR else 0.0
        p_co = (float(N[i["CO"]]) / ntot) * P if abs(ntot) >= _NTOT_FLOOR else 0.0
        e4 = float(rate_scale) * float(extent_R4(T, p_co2, p_co, M_kg))
    src = np.zeros(N_MAJ, dtype=np.float64)
    src[i["O2"]] += -alpha * e1
    src[i["CO"]] += 2.0 * (1.0 - alpha) * e1 + e2 + 2.0 * e4
    src[i["H2O"]] += -e2
    src[i["H2"]] += e2 - 2.0 * e3
    src[i["CH4"]] += e3
    src[i["CO2"]] += (2.0 * alpha - 1.0) * e1 - e4
    return src, e1 + e2 + e3 + e4


char_hetero_source = _char_gas_source


def char_residual(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
    m_char_in: float,
    K_char: float,
    rate_scale: float = 1.0,
    rho_ash: float = 75.0,
    d_p: float = 2.25e-3,
    D_g: float = 1.5e-4,
    rho_s: float = 1400.0,
) -> npt.NDArray[np.float64]:
    """14 个气相 [mol/s] + 炭库存 [kg]。``K_char`` [1/s]。"""
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (2 * N_MAJ + 1,)
    N_d, N_b, M = z[:N_MAJ], z[N_MAJ:2 * N_MAJ], float(z[-1])
    C_d, C_b = _conc(N_d, geo), _conc(N_b, geo)
    y_d, y_b = _y_full(N_d), _y_full(N_b)
    N_ex = float(geo.K_bd) * float(geo.V_b) * (C_b - C_d)
    g_d = _gas_source(C_d, y_d, geo=geo, V=geo.V_d, phase="dense",
                      rate_scale=rate_scale, rho_ash=rho_ash)
    g_b = _gas_source(C_b, y_b, geo=geo, V=geo.V_b, phase="bubble",
                      rate_scale=rate_scale, rho_ash=rho_ash)
    c_d, e_d = _char_gas_source(C_d, N_d, geo=geo, M_kg=M, d_p=d_p, D_g=D_g,
                                rho_s=rho_s, rate_scale=rate_scale)
    F_d = N_in_d - N_d + N_ex + g_d + c_d
    F_b = N_in_b - N_b - N_ex + g_b
    F_m = float(m_char_in) - float(K_char) * M - MW_C_KG_PER_MOL * e_d
    return np.concatenate([F_d, F_b, np.array([F_m], dtype=np.float64)])


def _H_gas(N: npt.NDArray[np.float64], T: float) -> float:
    """气相焓流 [W]。"""
    return float(sum(float(n) * enthalpy_molar(name, T) for name, n in zip(MAJOR, N)))


def energy_residual(
    x: npt.NDArray[np.float64],
    *,
    geo: ExchangeGeometry,
    N_in_d: npt.NDArray[np.float64],
    N_in_b: npt.NDArray[np.float64],
    T_in: float,
    m_char_in: float,
    K_char: float,
    rate_scale: float = 1.0,
    rho_ash: float = 75.0,
    d_p: float = 2.25e-3,
    D_g: float = 1.5e-4,
    rho_s: float = 1400.0,
    solid_enthalpy_mode: str = "omitted_named",
) -> npt.NDArray[np.float64]:
    """char 残差 + 气体焓 [W]。默认仍省略固体显热。

    ``thermo_consistent_cp_integral`` 才加入炭入口/出口焓流；不宣称能单独修奇异。
    """
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (2 * N_MAJ + 2,)
    T = float(z[-1])
    geo_T = ExchangeGeometry(T=T, P=geo.P, K_bd=geo.K_bd, V_b=geo.V_b, V_d=geo.V_d)
    F_mc = char_residual(
        z[:-1], geo=geo_T, N_in_d=N_in_d, N_in_b=N_in_b,
        m_char_in=m_char_in, K_char=K_char, rate_scale=rate_scale,
        rho_ash=rho_ash, d_p=d_p, D_g=D_g, rho_s=rho_s,
    )
    N_d, N_b, M = z[:N_MAJ], z[N_MAJ:2 * N_MAJ], float(z[2 * N_MAJ])
    H = (_H_gas(N_in_d, T_in) + _H_gas(N_in_b, T_in) - _H_gas(N_d, T) - _H_gas(N_b, T))
    if str(solid_enthalpy_mode) == SOLID_H_MODE_INTEGRAL:
        h_in = h_solid_mass("char", T_in, h_f=0.0, mode=SOLID_H_MODE_INTEGRAL)
        h_out = h_solid_mass("char", T, h_f=0.0, mode=SOLID_H_MODE_INTEGRAL)
        H += float(m_char_in) * h_in - float(K_char) * M * h_out
    return np.concatenate([F_mc, np.array([H], dtype=np.float64)])


def _solve_scaled(
    residual_factory: Callable[[float], Callable],
    x0: npt.NDArray[np.float64],
    scales: tuple[float, ...],
    *,
    jac_scale: npt.NDArray[np.float64],
    max_iter: int = 40,
) -> dict[str, Any]:
    """延拓：失败不覆盖上次成功 ``x``；完整速率失败不得被小参数收敛掩盖。"""
    x = np.asarray(x0, dtype=np.float64).copy()
    last_ok: float | None = None
    out: dict[str, Any] = {}
    for scale in scales:
        res_fn = residual_factory(scale)

        def jac_fn(z: npt.NDArray[np.float64], _res=res_fn) -> npt.NDArray[np.float64]:
            return fd_jacobian(_res, z, scale=jac_scale)

        sol = wirsum_1997_strict_newton(res_fn, jac_fn, x, max_iter=max_iter)
        key = "scale_0" if scale == 0.0 else f"scale_{scale:.0e}"
        rec = {
            **classify_newton_stop(sol, model_is_full_rate=abs(scale - 1.0) < 1.0e-15),
            "rate_scale": float(scale),
            "started_from_scale": last_ok,
        }
        if sol.converged:
            x = sol.x
            last_ok = float(scale)
            rec["x"] = np.asarray(sol.x, dtype=np.float64).tolist()
        out[key] = rec
    out["_last_successful_scale"] = last_ok
    out["_last_successful_x"] = np.asarray(x, dtype=np.float64).tolist()
    return out


def build_full_gas_char_energy_section() -> dict[str, Any]:
    geo = ExchangeGeometry()
    N_in_d, N_in_b = default_full_gas_inlets()
    m_char_in = 0.010  # [kg/s]
    K_char = 0.080  # [1/s]
    T_in = float(geo.T)
    scales = (0.0, 1.0e-12, 1.0e-8, 1.0e-4, 1.0)
    x_gas0 = np.concatenate([N_in_d, N_in_b])
    gas_scale = np.maximum(x_gas0, 0.1)

    def gas_factory(scale: float):
        def _res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return full_gas_residual(
                z, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, rate_scale=scale
            )
        return _res

    gas_cont = _solve_scaled(gas_factory, x_gas0, scales, jac_scale=gas_scale)
    x_gas = np.array(gas_cont.get("_last_successful_x") or x_gas0, dtype=np.float64)
    el_in = element_atom_flows(N_in_d, N_in_b)
    el_out = element_atom_flows(x_gas[:N_MAJ], x_gas[N_MAJ:])
    el_gap = {el: abs(el_out[el] - el_in[el]) for el in el_in}

    x_char0 = np.concatenate([x_gas, np.array([m_char_in / K_char])])
    char_scale = np.concatenate([gas_scale, np.array([1.0])])

    def char_factory(scale: float):
        def _res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return char_residual(
                z, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b,
                m_char_in=m_char_in, K_char=K_char, rate_scale=scale,
            )
        return _res

    char_cont = _solve_scaled(char_factory, x_char0, scales, jac_scale=char_scale)
    x_char = np.array(char_cont.get("_last_successful_x") or x_char0, dtype=np.float64)

    x_e0 = np.concatenate([x_char, np.array([T_in])])
    e_scale = np.concatenate([char_scale, np.array([1000.0])])

    def energy_factory(scale: float):
        def _res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
            return energy_residual(
                z, geo=geo, N_in_d=N_in_d, N_in_b=N_in_b, T_in=T_in,
                m_char_in=m_char_in, K_char=K_char, rate_scale=scale,
            )
        return _res

    energy_cont = _solve_scaled(energy_factory, x_e0, scales, jac_scale=e_scale, max_iter=50)

    def _public_cont(cont: dict[str, Any]) -> dict[str, Any]:
        scales_out = {
            k: {kk: vv for kk, vv in v.items() if kk != "x"}
            for k, v in cont.items()
            if isinstance(v, dict)
        }
        scales_out["last_successful_scale"] = cont.get("_last_successful_scale")
        return scales_out

    return {
        "full_gas": {
            "rung": "full_gas",
            "reactions": ["R5", "R6", "R7", "R8", "R12"],
            "not": ["R9", "R10", "R11", "fastox"],
            "units": {"state": "mol/s"},
            "species": list(MAJOR),
            "r6_stoich": "CH4+1.5O2->CO+2H2O",
            "continuation": _public_cont(gas_cont),
            "element_gap_mol_s": el_gap,
            "element_in_mol_s": el_in,
            "element_out_mol_s": el_out,
        },
        "char": {
            "rung": "char",
            "reactions": ["R1", "R2", "R3", "R4"],
            "units": {"gas": "mol/s", "holdup": "kg", "K": "1/s"},
            "m_char_in_kg_s": m_char_in,
            "K_char_1_per_s": K_char,
            "continuation": _public_cont(char_cont),
            "M_char_kg": float(x_char[-1]) if x_char.size == 15 else None,
        },
        "energy": {
            "rung": "energy",
            "units": {"energy": "W", "T": "K"},
            "T_in_K": T_in,
            "solid_sensible": "omitted_named_not_hamel_cp_integral",
            "char_hf": 0.0,
            "continuation": _public_cont(energy_cont),
            "T_K": (
                float(energy_cont["scale_1e+00"]["x"][-1])
                if "x" in energy_cont.get("scale_1e+00", {})
                else None
            ),
        },
    }
