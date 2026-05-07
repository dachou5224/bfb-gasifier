"""单 cell 守恒组装纯函数。

Ref: Hamel (1999) Eq. 2-4, 2-5, 2-6, 2-7, 2.8
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import numpy.typing as npt

from src.core.constants import T_REF

# Solid component column layout (must stay aligned with src/core/cell.py):
# [char, vm, moisture, ash]
S_CHAR, S_VM, S_MOISTURE, S_ASH = 0, 1, 2, 3
from src.core.species import GAS_SPECIES, N_GAS, cp_ash, cp_char, cp_sand, enthalpy_molar


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
) -> npt.NDArray[np.float64]:
    """计算颗粒缩减导致的跨粒径迁移。

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

    V_d_safe = max(float(V_d), 0.05 * float(V_cell))
    M_inv_total = float(rho_s) * (1.0 - float(eps_mf)) * V_d_safe

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

    m_bed = M_inv_total * (m_solid_arr / m_out_total)
    m_left = f_mig[:, None] * m_bed

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
) -> npt.NDArray[np.float64]:
    """组装固相质量守恒残差。

    Ref: Hamel (1999) Eq. 2.3
    """
    m_solid_arr = np.asarray(m_solid, dtype=np.float64)
    shape = m_solid_arr.shape
    m_solid_zu_arr = _as_float64_array("m_solid_zu", m_solid_zu, shape)
    m_solid_rez_arr = _as_float64_array("m_solid_rez", m_solid_rez, shape)
    m_solid_in_arr = _as_float64_array("m_solid_in", m_solid_in, shape)
    R_solid_arr = _as_float64_array("R_solid", R_solid, shape)
    size_migration_arr = _as_float64_array("size_migration", size_migration, shape)
    if solid_state_model == "legacy_stream":
        return m_solid_zu_arr + m_solid_rez_arr + m_solid_in_arr + R_solid_arr + size_migration_arr - m_solid_arr
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
    # Thesis holdup transport guardrail:
    # only char/ash participate in K-based axial transport terms.
    transport_mask = np.zeros(shape, dtype=np.float64)
    transport_mask[:, S_CHAR] = 1.0
    transport_mask[:, S_ASH] = 1.0
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
    res = (
        net_source_arr
        - m_solid_out_arr
    )
    # For disconnected holdup buckets (no source, no in/out transport, no reaction),
    # anchor residual to current inventory to avoid zero-row Jacobian singularity.
    k_sum_arr = K_solid_auf_arr + K_solid_ab_arr
    disconnected_mask = (
        np.abs(k_sum_arr) <= 1e-14
    ) & (
        np.abs(net_source_arr) <= 1e-14
    )
    # If a bucket has finite source but zero transport channel, ``net_source - K*m`` is
    # physically unsatisfiable at steady state (K=0). Add an inventory anchor term so
    # Newton has a well-conditioned direction to absorb the incoming source.
    source_without_transport_mask = (
        (transport_mask > 0.0)
        & (np.abs(k_sum_arr) <= 1e-14)
        & (np.abs(net_source_arr) > 1e-14)
    )
    if np.any(disconnected_mask):
        res = np.asarray(res, dtype=np.float64)
        res[disconnected_mask] = -m_solid_arr[disconnected_mask]
    if np.any(source_without_transport_mask):
        res = np.asarray(res, dtype=np.float64)
        res[source_without_transport_mask] = (
            net_source_arr[source_without_transport_mask] - m_solid_arr[source_without_transport_mask]
        )
    return res


def calc_gas_enthalpy_flow(
    N: npt.ArrayLike,
    T: float,
    h_cache: Dict[float, npt.NDArray[np.float64]] | None = None,
) -> float:
    """计算气相总焓流。

    Ref: Hamel (1999) Eq. 2.5
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
    return float(np.dot(np.nan_to_num(N_arr), h_vec))


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
) -> float:
    """计算固相总焓流。

    Ref: Hamel (1999) Eq. 2.5
    """
    m_arr = np.asarray(m, dtype=np.float64)
    assert m_arr.ndim == 2 and m_arr.shape[1] == 4, f"m shape mismatch: expected (_, 4), got {m_arr.shape}"
    w_ash = float(ash_dry_wt) / 100.0
    w_vm = (float(VM_daf) / 100.0) * (1.0 - w_ash)
    hf = np.array([0.0, float(h_f_dry) / max(w_vm, 1e-9), -15.866e6, 0.0], dtype=np.float64)
    cp_mix = 0.3 * float(cp_char_fn(T)) + 0.3 * float(cp_ash_fn(T)) + 0.4 * float(cp_sand_fn(T))
    return float(np.dot(np.sum(np.nan_to_num(m_arr), axis=0), hf) + np.sum(m_arr) * cp_mix * (float(T) - T_REF))


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
) -> float:
    """组装全床能量守恒残差。

    Ref: Hamel (1999) Eq. 2.7
    """
    H_in = (
        calc_gas_enthalpy_flow(np.asarray(N_b_in, dtype=np.float64) + np.asarray(N_d_in, dtype=np.float64), T_in_gas, h_cache)
        + calc_gas_enthalpy_flow(np.asarray(N_zu_b, dtype=np.float64) + np.asarray(N_zu_d, dtype=np.float64), T_zu_gas, h_cache)
        + calc_gas_enthalpy_flow(np.asarray(N_rez_b, dtype=np.float64) + np.asarray(N_rez_d, dtype=np.float64), T_rez_gas, h_cache)
        + calc_solid_enthalpy_flow(
            m_solid_rez,
            T_rez_solid,
            ash_dry_wt=ash_dry_wt,
            VM_daf=VM_daf,
            h_f_dry=h_f_dry,
            cp_char_fn=cp_char,
            cp_ash_fn=cp_ash,
            cp_sand_fn=cp_sand,
        )
        + calc_solid_enthalpy_flow(
            m_solid_in,
            T_in_solid,
            ash_dry_wt=ash_dry_wt,
            VM_daf=VM_daf,
            h_f_dry=h_f_dry,
            cp_char_fn=cp_char,
            cp_ash_fn=cp_ash,
            cp_sand_fn=cp_sand,
        )
        + calc_solid_enthalpy_flow(
            m_solid_zu,
            T_zu_solid,
            ash_dry_wt=ash_dry_wt,
            VM_daf=VM_daf,
            h_f_dry=h_f_dry,
            cp_char_fn=cp_char,
            cp_ash_fn=cp_ash,
            cp_sand_fn=cp_sand,
        )
    )
    H_out = calc_gas_enthalpy_flow(np.asarray(N_b, dtype=np.float64) + np.asarray(N_d, dtype=np.float64), T, h_cache) + calc_solid_enthalpy_flow(
        m_solid,
        T,
        ash_dry_wt=ash_dry_wt,
        VM_daf=VM_daf,
        h_f_dry=h_f_dry,
        cp_char_fn=cp_char,
        cp_ash_fn=cp_ash,
        cp_sand_fn=cp_sand,
    )
    return H_in * (1.0 - float(heat_loss_frac)) - H_out


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
