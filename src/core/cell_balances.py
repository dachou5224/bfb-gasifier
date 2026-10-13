"""单 cell 守恒组装纯函数。

Ref: Hamel (1999) Eq. 2-4, 2-5, 2-6, 2-7, 2.8
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import numpy.typing as npt

from src.core.constants import T_REF
from src.core.energy_balance_defaults import (
    CP_MOISTURE_J_PER_KG_K,
    HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT,
    HEAT_LOSS_MODE_INLET_FRACTION,
    SOLID_ENTHALPY_MODE_HAMEL_EQ28,
    SOLID_ENTHALPY_MODE_LEGACY_CP_MIX,
    solid_component_hf_j_per_kg,
)

from src.core.species import GAS_SPECIES, N_GAS, cp_ash, cp_char, cp_sand, enthalpy_molar

# Solid component column layout (must stay aligned with src/core/cell.py):
# [char, vm, moisture, ash]
S_CHAR, S_VM, S_MOISTURE, S_ASH = 0, 1, 2, 3

# K=0 时是否用 source−M 工程锚。engineering 为生产默认；none 为 thesis/strict。
HOLDUP_K0_ANCHOR_ENGINEERING: str = "engineering"
HOLDUP_K0_ANCHOR_NONE: str = "none"
_HOLDUP_K0_ANCHORS = frozenset(
    {HOLDUP_K0_ANCHOR_ENGINEERING, HOLDUP_K0_ANCHOR_NONE}
)


def resolve_holdup_k0_anchor(configured: str) -> str:
    """校验 ``holdup_k0_anchor``。"""
    name = str(configured).strip().lower()
    if name not in _HOLDUP_K0_ANCHORS:
        raise ValueError(
            f"Unknown holdup_k0_anchor={configured!r}; "
            f"expected {HOLDUP_K0_ANCHOR_ENGINEERING!r} or {HOLDUP_K0_ANCHOR_NONE!r}"
        )
    return name


def _as_float64_array(name: str, value: npt.ArrayLike, shape: tuple[int, ...]) -> npt.NDArray[np.float64]:
    arr = np.asarray(value, dtype=np.float64)
    assert arr.shape == shape, f"{name} shape mismatch: expected {shape}, got {arr.shape}"
    return arr


def calc_gas_balance_residual(
    *,
    N_zu_d: npt.ArrayLike,
    N_rez_d: npt.ArrayLike,
    N_d_in: npt.ArrayLike,
    R_gas_d: npt.ArrayLike,
    N_d: npt.ArrayLike,
    N_ex: npt.ArrayLike,
    N_zu_b: npt.ArrayLike,
    N_rez_b: npt.ArrayLike,
    N_b_in: npt.ArrayLike,
    R_gas_b: npt.ArrayLike,
    N_b: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """组装双相气体摩尔守恒残差。

    Ref: Hamel (1999) Eq. 2-4, 2-5.

    The thesis term ``ex_bd`` is positive from suspension to bubble. This
    implementation stores ``N_ex`` positive from bubble to suspension, i.e.
    ``N_ex = -ex_bd``.
    """
    shape = (N_GAS,)
    N_zu_d_arr = _as_float64_array("N_zu_d", N_zu_d, shape)
    N_rez_d_arr = _as_float64_array("N_rez_d", N_rez_d, shape)
    N_d_in_arr = _as_float64_array("N_d_in", N_d_in, shape)
    R_gas_d_arr = _as_float64_array("R_gas_d", R_gas_d, shape)
    N_d_arr = _as_float64_array("N_d", N_d, shape)
    N_ex_arr = _as_float64_array("N_ex", N_ex, shape)
    N_zu_b_arr = _as_float64_array("N_zu_b", N_zu_b, shape)
    N_rez_b_arr = _as_float64_array("N_rez_b", N_rez_b, shape)
    N_b_in_arr = _as_float64_array("N_b_in", N_b_in, shape)
    R_gas_b_arr = _as_float64_array("R_gas_b", R_gas_b, shape)
    N_b_arr = _as_float64_array("N_b", N_b, shape)

    res_d = N_zu_d_arr + N_rez_d_arr + N_d_in_arr + R_gas_d_arr - N_d_arr + N_ex_arr
    res_b = N_zu_b_arr + N_rez_b_arr + N_b_in_arr + R_gas_b_arr - N_b_arr - N_ex_arr
    return np.concatenate((res_d, res_b))


def holdup_transport_char_ash_mask(shape: tuple[int, ...]) -> npt.NDArray[np.float64]:
    """Char/ash-only axial transport mask for Eq. (2-6) holdup_transport.

    B-tier accepted gap: VM/moisture 不参与 K-based auf/ab 轴向输运；见
    ``docs/audit_reports/hamel_parity_audit_summary_2026-07-04.md`` §2.2。
    """
    mask = np.zeros(shape, dtype=np.float64)
    mask[:, S_CHAR] = 1.0
    mask[:, S_ASH] = 1.0
    return mask


def _apply_holdup_transport_numerical_guardrails(
    *,
    res: npt.NDArray[np.float64],
    m_solid_arr: npt.NDArray[np.float64],
    net_source_arr: npt.NDArray[np.float64],
    k_sum_arr: npt.NDArray[np.float64],
    transport_mask: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """NR 数值 guardrails（engineering，非 Hamel Eq. (2-6)）。

    * ``disconnected_mask``：无源、无输运、无反应的 holdup 桶写成 ``−m``。
    * ``source_without_transport_mask``：有源但 ``K=0`` 时写成 ``net_source−m``，
      把 [kg/s] 与 [kg] 相减（暗设 ``τ = 1 s``）。

    Ref: docs/solid_balance_formulation_audit.md
    """
    disconnected_mask = (np.abs(k_sum_arr) <= 1e-14) & (np.abs(net_source_arr) <= 1e-14)
    source_without_transport_mask = (
        (transport_mask > 0.0)
        & (np.abs(k_sum_arr) <= 1e-14)
        & (np.abs(net_source_arr) > 1e-14)
    )
    out = np.asarray(res, dtype=np.float64)
    if np.any(disconnected_mask):
        out = out.copy()
        out[disconnected_mask] = -m_solid_arr[disconnected_mask]
    if np.any(source_without_transport_mask):
        out = out.copy() if out is res else out
        out[source_without_transport_mask] = (
            net_source_arr[source_without_transport_mask]
            - m_solid_arr[source_without_transport_mask]
        )
    return out


def calc_size_migration(
    *,
    m_solid: npt.ArrayLike,
    R_solid: npt.ArrayLike,
    areas: npt.ArrayLike,
    d_p_classes: npt.ArrayLike,
    rho_s: float,
    eps_mf: float,
    V_d: float,
    V_cell: float,
    char_index: int,
    inventory_mode: str = "per_class_holdup",
) -> npt.NDArray[np.float64]:
    """计算颗粒缩减导致的跨粒径迁移。

    ``per_class_holdup``：按 size class 本地库存 ``m_jk``（Hamel Eq. 2.3–2.6 主链）。
    ``total_bed_proxy``：``ρ_s(1−ε_mf)V_d`` 全局床库存按 ``m/m_out`` 分配（审计对照）。

    Ref: Hamel (1999) Eq. 2.3 - 2.6
    """
    m_solid_arr = np.asarray(m_solid, dtype=np.float64)
    R_solid_arr = np.asarray(R_solid, dtype=np.float64)
    assert m_solid_arr.shape == R_solid_arr.shape, "m_solid and R_solid must share shape"
    nk, n_comp = m_solid_arr.shape
    if nk <= 1:
        return np.zeros_like(m_solid_arr)

    areas_arr = _as_float64_array("areas", areas, (nk,))
    dp_arr = _as_float64_array("d_p_classes", d_p_classes, (nk,))

    m_out_total = float(np.sum(np.maximum(m_solid_arr, 0.0)))
    if m_out_total <= 1e-12:
        return np.zeros_like(m_solid_arr)

    shrink_num = -R_solid_arr[:, char_index]
    denom_area = np.maximum(areas_arr * float(rho_s), 1e-12)
    r_shrink = np.where(areas_arr > 1e-15, shrink_num / denom_area, 0.0)
    r_shrink = np.where(r_shrink > 0.0, r_shrink, 0.0)

    dp_prev = np.concatenate((np.array([0.0]), dp_arr[:-1]))
    denom = 1.0 - (dp_prev / np.maximum(dp_arr, 1e-12)) ** 3
    denom = np.maximum(denom, 1e-4)
    f_mig = np.where(r_shrink > 0.0, (3.0 * r_shrink) / (dp_arr * denom), 0.0)

    mode = str(inventory_mode).strip()
    if mode == "per_class_holdup":
        m_bed = np.maximum(m_solid_arr, 0.0)
    elif mode == "total_bed_proxy":
        V_d_safe = max(float(V_d), 0.05 * float(V_cell))
        M_inv_total = float(rho_s) * (1.0 - float(eps_mf)) * V_d_safe
        m_bed = M_inv_total * (m_solid_arr / m_out_total)
    else:
        raise ValueError(
            f"Unknown size_migration inventory_mode={inventory_mode!r}; "
            "expected 'per_class_holdup' or 'total_bed_proxy'"
        )
    m_left = f_mig[:, None] * m_bed
    # The modeled size grid has no class below the smallest bin. Char
    # disappearance is already represented by R_solid; migration only
    # redistributes remaining particle mass between size classes.
    m_left[0, :] = 0.0

    res = -m_left
    res[:-1, :] += m_left[1:, :]
    return res


def calc_solid_balance_residual(
    *,
    m_solid_zu: npt.ArrayLike,
    m_solid_rez: npt.ArrayLike,
    m_solid_in: npt.ArrayLike,
    R_solid: npt.ArrayLike,
    size_migration: npt.ArrayLike,
    m_solid: npt.ArrayLike,
    m_solid_auf_in: npt.ArrayLike | None = None,
    m_solid_ab_in: npt.ArrayLike | None = None,
    K_solid_auf: npt.ArrayLike | None = None,
    K_solid_ab: npt.ArrayLike | None = None,
    solid_state_model: str = "legacy_stream",
    k0_anchor: str = HOLDUP_K0_ANCHOR_ENGINEERING,
) -> npt.NDArray[np.float64]:
    """组装固相质量守恒残差。

    ``legacy_stream``：Eq. (2-6) 源项求和；无 ``K·m`` 出流项（稳态 ``m`` 由入流
    + 反应 + 迁移闭合）。当连接矩阵提供 ``m_solid_auf_in`` / ``m_solid_ab_in`` 时
    显式计入，与 ``holdup_transport`` 一致。该分支仍含 ``−m_solid``，属 engineering。

    ``holdup_transport``：``F = net_source − (K_auf+K_ab) m``，单位 [kg/s]。
    ``k0_anchor='engineering'``（默认）在 ``K=0`` 时改写为 ``source−m`` / ``−m``；
    ``k0_anchor='none'`` 保持物理残差，允许零行。

    Ref: Hamel (1999) Eq. 2-6; docs/solid_balance_formulation_audit.md
    """
    m_solid_arr = np.asarray(m_solid, dtype=np.float64)
    shape = m_solid_arr.shape
    m_solid_zu_arr = _as_float64_array("m_solid_zu", m_solid_zu, shape)
    m_solid_rez_arr = _as_float64_array("m_solid_rez", m_solid_rez, shape)
    m_solid_in_arr = _as_float64_array("m_solid_in", m_solid_in, shape)
    R_solid_arr = _as_float64_array("R_solid", R_solid, shape)
    size_migration_arr = _as_float64_array("size_migration", size_migration, shape)
    if solid_state_model == "legacy_stream":
        m_solid_auf_in_arr = _as_float64_array(
            "m_solid_auf_in",
            np.zeros(shape, dtype=np.float64) if m_solid_auf_in is None else m_solid_auf_in,
            shape,
        )
        m_solid_ab_in_arr = _as_float64_array(
            "m_solid_ab_in",
            np.zeros(shape, dtype=np.float64) if m_solid_ab_in is None else m_solid_ab_in,
            shape,
        )
        return (
            m_solid_zu_arr
            + m_solid_rez_arr
            + m_solid_in_arr
            + m_solid_auf_in_arr
            + m_solid_ab_in_arr
            + R_solid_arr
            + size_migration_arr
            - m_solid_arr
        )
    if solid_state_model == "freeboard_closure":
        return np.zeros_like(m_solid_arr, dtype=np.float64)
    if solid_state_model != "holdup_transport":
        raise ValueError(f"Unknown solid_state_model: {solid_state_model}")

    m_solid_auf_in_arr = _as_float64_array(
        "m_solid_auf_in",
        np.zeros(shape, dtype=np.float64) if m_solid_auf_in is None else m_solid_auf_in,
        shape,
    )
    m_solid_ab_in_arr = _as_float64_array(
        "m_solid_ab_in",
        np.zeros(shape, dtype=np.float64) if m_solid_ab_in is None else m_solid_ab_in,
        shape,
    )
    K_solid_auf_arr = _as_float64_array(
        "K_solid_auf",
        np.zeros(shape, dtype=np.float64) if K_solid_auf is None else K_solid_auf,
        shape,
    )
    K_solid_ab_arr = _as_float64_array(
        "K_solid_ab",
        np.zeros(shape, dtype=np.float64) if K_solid_ab is None else K_solid_ab,
        shape,
    )
    transport_mask = holdup_transport_char_ash_mask(shape)
    m_solid_auf_in_arr = m_solid_auf_in_arr * transport_mask
    m_solid_ab_in_arr = m_solid_ab_in_arr * transport_mask
    K_solid_auf_arr = K_solid_auf_arr * transport_mask
    K_solid_ab_arr = K_solid_ab_arr * transport_mask
    m_solid_out_arr = (K_solid_auf_arr + K_solid_ab_arr) * m_solid_arr
    net_source_arr = (
        m_solid_zu_arr
        + m_solid_rez_arr
        + m_solid_in_arr
        + m_solid_auf_in_arr
        + m_solid_ab_in_arr
        + R_solid_arr
        + size_migration_arr
    )
    res = net_source_arr - m_solid_out_arr
    k_sum_arr = K_solid_auf_arr + K_solid_ab_arr
    if resolve_holdup_k0_anchor(k0_anchor) == HOLDUP_K0_ANCHOR_NONE:
        return res
    return _apply_holdup_transport_numerical_guardrails(
        res=res,
        m_solid_arr=m_solid_arr,
        net_source_arr=net_source_arr,
        k_sum_arr=k_sum_arr,
        transport_mask=transport_mask,
    )


def calc_gas_enthalpy_flow(
    N: npt.ArrayLike,
    T: float,
    h_cache: Dict[float, npt.NDArray[np.float64]] | None = None,
) -> float:
    """计算气相总焓流。

    Ref: Hamel (1999) Eq. (2.8)（气相 Σ N_j h_j）
    """
    N_arr = _as_float64_array("N", N, (N_GAS,))
    T_key = round(float(T), 3)
    h_vec = None if h_cache is None else h_cache.get(T_key)
    if h_vec is None:
        h_vec = np.array([enthalpy_molar(sp, float(T)) for sp in GAS_SPECIES], dtype=np.float64)
        if h_cache is not None:
            if len(h_cache) >= 4:
                first_key = next(iter(h_cache.keys()))
                h_cache.pop(first_key, None)
            h_cache[T_key] = h_vec
    return float(np.dot(N_arr, h_vec))


def calc_solid_enthalpy_flow(
    m: npt.ArrayLike,
    T: float,
    *,
    ash_dry_wt: float,
    VM_daf: float,
    h_f_dry: float,
    cp_char_fn,
    cp_ash_fn,
    cp_sand_fn,
    solid_enthalpy_mode: str = SOLID_ENTHALPY_MODE_HAMEL_EQ28,
) -> float:
    """计算固相总焓流。

    Ref: Hamel (1999) Eq. (2.8) 只给出 ``Σ m (h_f + h_sensible)``。
    默认 ``hamel_eq28_component`` 用端点瞬时 ``c_p(T)(T−T_ref)``，
    即 ``legacy_endpoint_cp``，**不是**积分焓，也不是已核验的 Hamel 物性系数。
    积分见 ``species.h_solid_mass(..., mode=thermo_consistent_cp_integral)``。
    ``legacy_cp_mix_030304``：0.3/0.3/0.4 混合比热（审计对照）。
    """
    m_arr = np.asarray(m, dtype=np.float64)
    assert m_arr.ndim == 2 and m_arr.shape[1] == 4, f"m shape mismatch: expected (_, 4), got {m_arr.shape}"
    m_comp = np.sum(m_arr, axis=0)
    m_total = float(np.sum(m_arr))
    dT = float(T) - float(T_REF)
    hf = np.array(
        solid_component_hf_j_per_kg(
            h_f_dry=float(h_f_dry),
            ash_dry_wt=float(ash_dry_wt),
            VM_daf=float(VM_daf),
        ),
        dtype=np.float64,
    )

    mode = str(solid_enthalpy_mode or SOLID_ENTHALPY_MODE_HAMEL_EQ28).strip().lower()
    if mode == SOLID_ENTHALPY_MODE_LEGACY_CP_MIX:
        cp_mix = 0.3 * float(cp_char_fn(T)) + 0.3 * float(cp_ash_fn(T)) + 0.4 * float(cp_sand_fn(T))
        return float(np.dot(m_comp, hf) + m_total * cp_mix * dT)

    if mode != SOLID_ENTHALPY_MODE_HAMEL_EQ28:
        raise ValueError(
            f"Unsupported solid_enthalpy_mode={solid_enthalpy_mode!r}; "
            f"expected {SOLID_ENTHALPY_MODE_HAMEL_EQ28!r} or {SOLID_ENTHALPY_MODE_LEGACY_CP_MIX!r}"
        )
    cp_comp = np.array(
        [
            float(cp_char_fn(T)),
            float(cp_char_fn(T)),
            float(CP_MOISTURE_J_PER_KG_K),
            float(cp_ash_fn(T)),
        ],
        dtype=np.float64,
    )
    h_comp = hf + cp_comp * dT
    return float(np.dot(m_comp, h_comp))


def calc_energy_balance_H_in(
    *,
    N_b_in: npt.ArrayLike,
    N_d_in: npt.ArrayLike,
    T_in_gas: float,
    N_zu_b: npt.ArrayLike,
    N_zu_d: npt.ArrayLike,
    T_zu_gas: float,
    N_rez_b: npt.ArrayLike,
    N_rez_d: npt.ArrayLike,
    T_rez_gas: float,
    m_solid_rez: npt.ArrayLike,
    T_rez_solid: float,
    m_solid_in: npt.ArrayLike,
    T_in_solid: float,
    m_solid_zu: npt.ArrayLike,
    T_zu_solid: float,
    h_cache: Dict[float, npt.NDArray[np.float64]] | None = None,
    m_solid_auf_in: npt.ArrayLike | None = None,
    T_solid_auf_in: float | None = None,
    m_solid_ab_in: npt.ArrayLike | None = None,
    T_solid_ab_in: float | None = None,
    ash_dry_wt: float,
    VM_daf: float,
    h_f_dry: float,
    solid_enthalpy_mode: str = SOLID_ENTHALPY_MODE_HAMEL_EQ28,
) -> float:
    """Eq. (2.7) 入流焓 H_in [W]（不含出口项与热损）。"""
    _solid_kw = dict(
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
        h_f_dry=h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
        solid_enthalpy_mode=solid_enthalpy_mode,
    )
    return (
        calc_gas_enthalpy_flow(np.asarray(N_b_in, dtype=np.float64) + np.asarray(N_d_in, dtype=np.float64), T_in_gas, h_cache)
        + calc_gas_enthalpy_flow(np.asarray(N_zu_b, dtype=np.float64) + np.asarray(N_zu_d, dtype=np.float64), T_zu_gas, h_cache)
        + calc_gas_enthalpy_flow(np.asarray(N_rez_b, dtype=np.float64) + np.asarray(N_rez_d, dtype=np.float64), T_rez_gas, h_cache)
        + calc_solid_enthalpy_flow(m_solid_rez, T_rez_solid, **_solid_kw)
        + calc_solid_enthalpy_flow(m_solid_in, T_in_solid, **_solid_kw)
        + calc_solid_enthalpy_flow(
            np.zeros_like(np.asarray(m_solid_in, dtype=np.float64)) if m_solid_auf_in is None else m_solid_auf_in,
            T_in_solid if T_solid_auf_in is None else float(T_solid_auf_in),
            **_solid_kw,
        )
        + calc_solid_enthalpy_flow(
            np.zeros_like(np.asarray(m_solid_in, dtype=np.float64)) if m_solid_ab_in is None else m_solid_ab_in,
            T_in_solid if T_solid_ab_in is None else float(T_solid_ab_in),
            **_solid_kw,
        )
        + calc_solid_enthalpy_flow(m_solid_zu, T_zu_solid, **_solid_kw)
    )


def calc_energy_balance_residual(
    *,
    N_b_in: npt.ArrayLike,
    N_d_in: npt.ArrayLike,
    T_in_gas: float,
    N_zu_b: npt.ArrayLike,
    N_zu_d: npt.ArrayLike,
    T_zu_gas: float,
    N_rez_b: npt.ArrayLike,
    N_rez_d: npt.ArrayLike,
    T_rez_gas: float,
    m_solid_rez: npt.ArrayLike,
    T_rez_solid: float,
    m_solid_in: npt.ArrayLike,
    T_in_solid: float,
    m_solid_zu: npt.ArrayLike,
    T_zu_solid: float,
    N_b: npt.ArrayLike,
    N_d: npt.ArrayLike,
    m_solid: npt.ArrayLike,
    T: float,
    heat_loss_frac: float,
    ash_dry_wt: float,
    VM_daf: float,
    h_f_dry: float,
    h_cache: Dict[float, npt.NDArray[np.float64]] | None = None,
    m_solid_auf_in: npt.ArrayLike | None = None,
    T_solid_auf_in: float | None = None,
    m_solid_ab_in: npt.ArrayLike | None = None,
    T_solid_ab_in: float | None = None,
    solid_enthalpy_mode: str = SOLID_ENTHALPY_MODE_HAMEL_EQ28,
    heat_loss_mode: str = HEAT_LOSS_MODE_INLET_FRACTION,
    Q_wall_W: float = 0.0,
    Q_heat_exchanger_W: float = 0.0,
) -> float:
    """组装全床能量守恒残差。

    Ref: Hamel (1999) Eq. (2.7)
        H_zu + H_rez + H_auf + H_ab − H_aus − Q_W − Q_WÜ = 0
    （化学反应热已含于气/固源项 R_gas / R_solid，不在此式重复扣除。）

    ``heat_loss_mode``:
    - ``inlet_fraction``（默认）：Q_loss = heat_loss_frac · max(H_in, 0)（仅对正入流焓散热）
    - ``hamel_eq27_explicit``：Q_loss = Q_wall_W + Q_heat_exchanger_W
    """
    _solid_kw = dict(
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
        h_f_dry=h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
        solid_enthalpy_mode=solid_enthalpy_mode,
    )
    H_in = calc_energy_balance_H_in(
        N_b_in=N_b_in,
        N_d_in=N_d_in,
        T_in_gas=T_in_gas,
        N_zu_b=N_zu_b,
        N_zu_d=N_zu_d,
        T_zu_gas=T_zu_gas,
        N_rez_b=N_rez_b,
        N_rez_d=N_rez_d,
        T_rez_gas=T_rez_gas,
        m_solid_rez=m_solid_rez,
        T_rez_solid=T_rez_solid,
        m_solid_in=m_solid_in,
        T_in_solid=T_in_solid,
        m_solid_zu=m_solid_zu,
        T_zu_solid=T_zu_solid,
        h_cache=h_cache,
        m_solid_auf_in=m_solid_auf_in,
        T_solid_auf_in=T_solid_auf_in,
        m_solid_ab_in=m_solid_ab_in,
        T_solid_ab_in=T_solid_ab_in,
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
        h_f_dry=h_f_dry,
        solid_enthalpy_mode=solid_enthalpy_mode,
    )
    H_out = calc_gas_enthalpy_flow(np.asarray(N_b, dtype=np.float64) + np.asarray(N_d, dtype=np.float64), T, h_cache) + calc_solid_enthalpy_flow(
        m_solid,
        T,
        **_solid_kw,
    )
    loss_mode = str(heat_loss_mode or HEAT_LOSS_MODE_INLET_FRACTION).strip().lower()
    if loss_mode == HEAT_LOSS_MODE_INLET_FRACTION:
        q_loss = float(heat_loss_frac) * max(float(H_in), 0.0)
    elif loss_mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT:
        q_loss = float(Q_wall_W) + float(Q_heat_exchanger_W)
    else:
        raise ValueError(
            f"Unsupported heat_loss_mode={heat_loss_mode!r}; "
            f"expected {HEAT_LOSS_MODE_INLET_FRACTION!r} or {HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT!r}"
        )
    return H_in - H_out - q_loss


def assemble_cell_residual_vector(
    *,
    res_gas: npt.ArrayLike,
    res_solid: npt.ArrayLike,
    res_energy: float,
    out: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """拼接单 cell 原始残差向量。"""
    res_gas_arr = np.asarray(res_gas, dtype=np.float64)
    res_solid_arr = np.asarray(res_solid, dtype=np.float64).reshape(-1)
    expected = res_gas_arr.size + res_solid_arr.size + 1
    assert out.shape == (expected,), f"out shape mismatch: expected {(expected,)}, got {out.shape}"
    out[: res_gas_arr.size] = res_gas_arr
    out[res_gas_arr.size : res_gas_arr.size + res_solid_arr.size] = res_solid_arr
    out[-1] = float(res_energy)
    if np.any(np.isnan(out)) or np.any(np.isinf(out)):
        return np.nan_to_num(out, nan=1e6, posinf=1e10, neginf=-1e6)
    return out
