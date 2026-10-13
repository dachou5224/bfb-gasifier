"""能量行 + 约化外层 ``q=P(x)``。不经 ``Reactor``，不切生产默认。

在 C+D 的 39 方程上增加一个床温未知量和 Hamel Eq. 2.7 总焓残差。
固体焓用 ``thermo_consistent_cp_integral``。``Q_W=0`` 是省略壁面换热，
不是为修奇异虚构的源。气体 ``enthalpy_molar`` 已含生成焓，不再加一遍 hf。

外层只按已接受 ``T`` 重算热解映射；不重写未知 ``N``/``T``，不更新旋风 K。
``K_eject`` 与占据年龄仍冻结。本轮不开放生产 ``hamel_outer_refresh``。

Source: Hamel (1999) Eq. 2.7–2.8, §2.1 Bild 2.2;
docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q3–Q5, Q4⑦
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
import numpy.typing as npt

from src.core.energy_balance_defaults import (
    H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
    H_FORMATION_LIQUID_WATER_J_PER_KG,
)
from src.core.species import (
    SOLID_H_MODE_INTEGRAL,
    enthalpy_molar,
    gas_diffusivity_correlation,
    h_solid_mass,
)
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import (
    N_X as N_SPECIES,
    OCCUPANCY_SOURCE,
    CoupledTwoBinPyroQ,
    apply_pyrolysis_map,
    coupled_residual as species_residual,
    default_char_feed,
    default_x0 as species_x0,
    edges_from_x,
    freeze_q as freeze_species_q,
    solve_coupled as solve_species,
    unpack as unpack_species,
)
from src.solvers.wirsum_1997_cyclone_cells import CycloneCellK
from src.solvers.wirsum_1997_pyrolysis_ledger import (
    FuelAtomFeed,
    map_one_state,
    unreleased_vm_enthalpy_W,
)
from src.solvers.wirsum_1997_stiffness_full_gas import MAJOR, fd_jacobian
from src.solvers.wirsum_1997_stiffness_ladder import ExchangeGeometry, classify_newton_stop
from src.solvers.wirsum_1997_strict import wirsum_1997_strict_newton

N_X: int = N_SPECIES + 1
F_NAME: str = "F_E0"
Q_WALL_W: float = 0.0  # [W] 省略壁面换热，不是虚构闭合
OMITTED: tuple[str, ...] = (
    "energy_wall_heat",
    "cyclone_K_update_P",
    "K_eject_in_this_P",
    "DAEM_X_from_T",
    "R2_R12",
    "tar_transport",
    "ash_inventory",
    "freeboard_reaction",
    "return_leg_holdup",
    "production_hamel_outer_refresh",
)


@dataclass(frozen=True)
class CoupledEnergyQ:
    """物种 q + 原料账。热解源在外层更新，内层只读。"""

    species: CoupledTwoBinPyroQ
    feed: FuelAtomFeed
    T_in: float  # [K]
    n_tar: dict[str, float]  # [mol/s]
    solid_enthalpy_mode: str
    occupancy_source: str


def _H_major(N: npt.NDArray[np.float64], T: float) -> float:
    """MAJOR 焓流 [W]。``enthalpy_molar`` 已含生成焓。"""
    return float(sum(float(n) * enthalpy_molar(name, T) for name, n in zip(MAJOR, N)))


def _H_tar(n_tar: dict[str, float], T: float) -> float:
    return float(sum(float(n) * enthalpy_molar(name, T) for name, n in n_tar.items()))


def q_at_T(q: CoupledTwoBinPyroQ, T: float) -> CoupledTwoBinPyroQ:
    """内层物性随 T；热解源与 ``K_eject`` 仍是传入的 q。"""
    T = float(T)
    return replace(
        q,
        geo_ex=replace(q.geo_ex, T=T),
        D_g=gas_diffusivity_correlation(T, float(q.geo_ex.P)),
    )


def freeze_energy_q(
    *,
    T_in: float = 300.0,
    T_map: float = 1100.0,
    m_char_kg_s: float = 0.002,
    k: CycloneCellK | None = None,
    feed: FuelAtomFeed | None = None,
) -> CoupledEnergyQ:
    fuel = feed or default_char_feed(m_char_kg_s=m_char_kg_s)
    geo = ExchangeGeometry(T=float(T_map))
    q = freeze_species_q(geo_ex=geo, m_char_kg_s=m_char_kg_s, k=k, feed=fuel)
    mapped = map_one_state(
        fuel,
        P_Pa=float(q.geo_ex.P),
        T_K=float(T_map),
        X_vm=np.array([0.0, 1.0], dtype=np.float64),
        tau_s=np.array([1.0], dtype=np.float64),
        occupancy_source=OCCUPANCY_SOURCE,
        fuel_type="coal",
    )
    n_tar = {
        "TAR1": float(mapped["cells"][0]["n_py_mol_s"].get("TAR1", 0.0)),
        "TAR2": float(mapped["cells"][0]["n_py_mol_s"].get("TAR2", 0.0)),
    }
    return CoupledEnergyQ(
        species=q,
        feed=fuel,
        T_in=float(T_in),
        n_tar=n_tar,
        solid_enthalpy_mode=SOLID_H_MODE_INTEGRAL,
        occupancy_source=OCCUPANCY_SOURCE,
    )


def unpack(x: npt.NDArray[np.float64]) -> dict[str, Any]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_X,), f"expected {N_X}, got {z.shape}"
    u = unpack_species(z[:N_SPECIES])
    u["T"] = float(z[-1])
    return u


def pack(species_x: npt.NDArray[np.float64], T: float) -> npt.NDArray[np.float64]:
    return np.concatenate(
        [np.asarray(species_x, dtype=np.float64), np.array([T], dtype=np.float64)]
    )


def energy_residual_W(x: npt.NDArray[np.float64], qe: CoupledEnergyQ) -> float:
    """总焓残差 [W]：入口燃料+进气 − 出口气+飞灰+焦油+未转化 VM − ``Q_W``。

    热解气不是第二套入口；它们由燃料焓经内部分配进入 ``N``。
    未释放挥发分是固体未转化出口，与气体/焦/灰同一 ``T``。
    """
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-1])
    qT = q_at_T(qe.species, T)
    u = unpack_species(z[:N_SPECIES])
    ed = edges_from_x(z[:N_SPECIES], qT)
    H_in = (
        _H_major(qe.species.N_in_d, qe.T_in)
        + _H_major(qe.species.N_in_b, qe.T_in)
        + float(qe.feed.m_dry_kg_s) * float(qe.feed.h_f_dry_J_kg)
        + float(qe.feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    )
    h_char = h_solid_mass(
        "char",
        T,
        h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    H_unreleased = unreleased_vm_enthalpy_W(
        qe.feed, float(qe.species.m_unreleased_kg_s), T
    )
    H_out = (
        _H_major(u["N_c2"], T)
        + float(ed["m_fly"]) * h_char
        + _H_tar(qe.n_tar, T)
        + H_unreleased
    )
    return float(H_in - H_out - Q_WALL_W)


def enthalpy_buckets_W(x: npt.NDArray[np.float64], qe: CoupledEnergyQ) -> dict[str, float]:
    """总焓残差拆桶 [W]。各桶有号之和等于 ``energy_residual_W``。"""
    z = np.asarray(x, dtype=np.float64)
    T = float(z[-1])
    qT = q_at_T(qe.species, T)
    u = unpack_species(z[:N_SPECIES])
    ed = edges_from_x(z[:N_SPECIES], qT)
    H_gas_in = _H_major(qe.species.N_in_d, qe.T_in) + _H_major(qe.species.N_in_b, qe.T_in)
    H_fuel = float(qe.feed.m_dry_kg_s) * float(qe.feed.h_f_dry_J_kg)
    H_moist = float(qe.feed.m_moist_kg_s) * H_FORMATION_LIQUID_WATER_J_PER_KG
    h_char = h_solid_mass(
        "char",
        T,
        h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,
        mode=SOLID_H_MODE_INTEGRAL,
    )
    H_gas_out = _H_major(u["N_c2"], T)
    H_fly = float(ed["m_fly"]) * h_char
    H_tar = _H_tar(qe.n_tar, T)
    H_unreleased = unreleased_vm_enthalpy_W(
        qe.feed, float(qe.species.m_unreleased_kg_s), T
    )
    buckets = {
        "H_gas_in_W": H_gas_in,
        "H_fuel_dry_W": H_fuel,
        "H_moist_liquid_W": H_moist,
        "H_gas_out_W": -H_gas_out,
        "H_fly_char_W": -H_fly,
        "H_tar_omitted_W": -H_tar,
        "H_unreleased_vm_W": -H_unreleased,
        "Q_W": -Q_WALL_W,
    }
    summed = float(sum(buckets.values()))
    return {
        **buckets,
        "sum_W": summed,
        "energy_W": float(energy_residual_W(z, qe)),
        "T_K": T,
        "T_in_K": float(qe.T_in),
        "m_fly_kg_s": float(ed["m_fly"]),
    }


def scan_E_at_species_roots(
    qe: CoupledEnergyQ,
    T_K: npt.NDArray[np.float64],
    *,
    refresh_pyrolysis: bool = False,
    max_iter: int = 40,
) -> dict[str, Any]:
    """固定或按 T 刷新 q，每个 T 先求 39 维物种根再评 E。不改生产默认。"""
    T_arr = np.asarray(T_K, dtype=np.float64)
    assert T_arr.ndim == 1 and T_arr.size >= 2
    rows: list[dict[str, Any]] = []
    x_spec = None
    q_use = qe
    for T in T_arr:
        T = float(T)
        dummy = pack(species_x0(q_use.species), T)
        if refresh_pyrolysis:
            q_use = outer_P(dummy, qe)
        qT = q_at_T(q_use.species, T)
        spec = solve_species(qT, x0=x_spec, max_iter=max_iter)
        x_spec = spec["x"]
        x = pack(x_spec, T)
        buckets = enthalpy_buckets_W(x, q_use)
        rows.append(
            {
                "T_K": T,
                "energy_W": buckets["energy_W"],
                "species_physical": bool(spec["physical_convergence"]),
                "species_residual_l2": float(spec["residual_l2"]),
                "refresh_pyrolysis": bool(refresh_pyrolysis),
                **{k: buckets[k] for k in buckets if k.endswith("_W") or k == "m_fly_kg_s"},
            }
        )
    E = np.array([r["energy_W"] for r in rows], dtype=np.float64)
    phys = np.array([r["species_physical"] for r in rows], dtype=bool)
    signs = np.sign(E[phys]) if np.any(phys) else np.array([], dtype=np.float64)
    sign_change = bool(signs.size >= 2 and np.min(signs) < 0.0 and np.max(signs) > 0.0)
    return {
        "rows": rows,
        "T_K": T_arr.copy(),
        "energy_W": E,
        "species_physical": phys,
        "sign_change": sign_change,
        "all_same_sign": bool(phys.any() and not sign_change and np.all(np.sign(E[phys]) == np.sign(E[phys][0]))),
        "Q_W": Q_WALL_W,
        "refresh_pyrolysis": bool(refresh_pyrolysis),
        "not_a_furnace_convergence_claim": True,
    }


def coupled_residual(
    x: npt.NDArray[np.float64],
    qe: CoupledEnergyQ,
    *,
    rate_scale: float = 1.0,
) -> npt.NDArray[np.float64]:
    z = np.asarray(x, dtype=np.float64)
    assert z.shape == (N_X,)
    T = float(z[-1])
    qT = q_at_T(qe.species, T)
    F_s = species_residual(z[:N_SPECIES], qT, rate_scale=rate_scale)
    return np.concatenate([F_s, np.array([energy_residual_W(z, qe)], dtype=np.float64)])


def default_x0(qe: CoupledEnergyQ) -> npt.NDArray[np.float64]:
    xs = species_x0(qe.species)
    return pack(xs, float(qe.species.geo_ex.T))


def outer_P(x: npt.NDArray[np.float64], qe: CoupledEnergyQ) -> CoupledEnergyQ:
    """``q = P(x)``：用当前 T 重算热解源。不修改 ``x``，不更新旋风 K / ``K_eject``。"""
    T = float(np.asarray(x, dtype=np.float64)[-1])
    mapped = map_one_state(
        qe.feed,
        P_Pa=float(qe.species.geo_ex.P),
        T_K=T,
        X_vm=np.array([0.0, 1.0], dtype=np.float64),
        tau_s=np.array([1.0], dtype=np.float64),
        occupancy_source=qe.occupancy_source,
        fuel_type="coal",
    )
    n_tar = {
        "TAR1": float(mapped["cells"][0]["n_py_mol_s"].get("TAR1", 0.0)),
        "TAR2": float(mapped["cells"][0]["n_py_mol_s"].get("TAR2", 0.0)),
    }
    q2 = apply_pyrolysis_map(qe.species, mapped)
    return CoupledEnergyQ(
        species=q2,
        feed=qe.feed,
        T_in=qe.T_in,
        n_tar=n_tar,
        solid_enthalpy_mode=qe.solid_enthalpy_mode,
        occupancy_source=qe.occupancy_source,
    )


def dump_e0_E_of_T_json(path: str) -> dict[str, Any]:
    """门禁 0：把 F_E0 的 E(T) 基线写成 JSON。不是炉温目标。"""
    import json
    from pathlib import Path

    qe = freeze_energy_q()
    T = np.array([500.0, 600.0, 1100.0, 1500.0], dtype=np.float64)
    out = scan_E_at_species_roots(qe, T, refresh_pyrolysis=False)
    payload = {
        "f_name": F_NAME,
        "n_unknowns": N_X,
        "Q_W": Q_WALL_W,
        "Q_W_meaning": "omitted_named",
        "X_vm": [0.0, 1.0],
        "outer_P_updates_X_vm": False,
        "not_a_furnace_convergence_claim": True,
        "T_K": [float(v) for v in out["T_K"]],
        "energy_W": [float(v) for v in out["energy_W"]],
        "sign_change": bool(out["sign_change"]),
        "rows": [
            {
                "T_K": float(r["T_K"]),
                "energy_W": float(r["energy_W"]),
                "species_physical": bool(r["species_physical"]),
            }
            for r in out["rows"]
        ],
    }
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def omitted_named() -> dict[str, str]:
    return {name: "omitted, named" for name in OMITTED}


def solve_inner(
    qe: CoupledEnergyQ,
    *,
    x0: npt.NDArray[np.float64] | None = None,
    rate_scale: float = 1.0,
    max_iter: int = 40,
) -> dict[str, Any]:
    """固定 q 的内层。物种在当前 T 上求根，再联立能量。"""
    if x0 is None:
        T0 = float(qe.species.geo_ex.T)
        spec = solve_species(q_at_T(qe.species, T0), rate_scale=rate_scale, max_iter=max_iter)
        x_start = pack(spec["x"], T0)
    else:
        x_start = np.asarray(x0, dtype=np.float64)
        T0 = float(x_start[-1])
        spec = solve_species(
            q_at_T(qe.species, T0),
            x0=x_start[:N_SPECIES],
            rate_scale=rate_scale,
            max_iter=max_iter,
        )
        x_start = pack(spec["x"], T0)
    scale = np.maximum(np.abs(x_start), 0.05)
    scale[-1] = max(float(scale[-1]), 50.0)

    def res(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        if float(np.min(z[:-1])) < 0.0 or float(z[-1]) <= 0.0:
            return np.full(N_X, 1.0e6, dtype=np.float64)
        return coupled_residual(z, qe, rate_scale=rate_scale)

    def jac(z: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return fd_jacobian(
            lambda w: coupled_residual(w, qe, rate_scale=rate_scale),
            z,
            scale=scale,
        )

    sol = wirsum_1997_strict_newton(res, jac, x_start, max_iter=max_iter)
    stop = classify_newton_stop(sol, model_is_full_rate=abs(rate_scale - 1.0) < 1.0e-15)
    F2 = coupled_residual(sol.x, qe, rate_scale=rate_scale)
    return {
        **stop,
        "x": np.asarray(sol.x, dtype=np.float64),
        "F": np.asarray(sol.F, dtype=np.float64),
        "species_physical": bool(spec["physical_convergence"]),
        "species_residual_l2": float(spec["residual_l2"]),
        "q_purity_gap": float(np.max(np.abs(F2 - sol.F))),
        "rank_last": int(sol.rank_last),
        "energy_W": float(energy_residual_W(sol.x, qe)),
        "T_K": float(sol.x[-1]),
        "rewrites_unknown_N": False,
    }


def solve_outer(
    qe: CoupledEnergyQ,
    *,
    n_outer: int = 3,
    max_iter: int = 40,
) -> dict[str, Any]:
    """外层 ``q_{k+1}=P(x_k)``。验收 P 与最终 F，不只看残差变小。"""
    q = qe
    x: npt.NDArray[np.float64] | None = None
    hist: list[dict[str, Any]] = []
    last: dict[str, Any] = {}
    for k in range(int(n_outer)):
        last = solve_inner(q, x0=x, max_iter=max_iter)
        x = last["x"]
        x_before = x.copy()
        q_next = outer_P(x, q)
        F_old = coupled_residual(x, q)
        F_new = coupled_residual(x, q_next)
        hist.append(
            {
                "k": k,
                "T_K": float(x[-1]),
                "energy_W": float(energy_residual_W(x, q)),
                "n_py_sum": float(np.sum(q.species.n_py_d)),
                "n_py_next_sum": float(np.sum(q_next.species.n_py_d)),
                "x_unchanged": bool(np.allclose(x_before, x)),
                "F_after_P_l2": float(np.linalg.norm(F_new)),
                "F_before_P_l2": float(np.linalg.norm(F_old)),
                "inner_physical": bool(last["physical_convergence"]),
            }
        )
        q = q_next
    assert x is not None
    F_final = coupled_residual(x, q)
    return {
        "x": x,
        "q": q,
        "F": F_final,
        "history": hist,
        "rewrites_unknown_N": False,
        "rewrites_unknown_T": False,
        "K_eject_updated": False,
        "production_policy_unchanged": True,
        "inner_last": last,
    }


def build_stage_energy_outer_payload() -> dict[str, Any]:
    qe = freeze_energy_q()
    inner = solve_inner(qe)
    outer = solve_outer(qe, n_outer=2)
    return {
        "source": "Wirsum alignment energy + reduced outer P",
        "f_name": F_NAME,
        "solved": False,
        "not_a_furnace_convergence_claim": True,
        "n_unknowns": N_X,
        "Q_W": Q_WALL_W,
        "Q_W_meaning": "omitted_named",
        "outer_P_X_vm_frozen": True,
        "solid_enthalpy_mode": SOLID_H_MODE_INTEGRAL,
        "gas_hf_counted_once": True,
        "inner_physical_convergence": inner["physical_convergence"],
        "species_physical": inner["species_physical"],
        "inner_residual_l2": inner["residual_l2"],
        "energy_W": inner["energy_W"],
        "T_K": inner["T_K"],
        "rank_last": inner["rank_last"],
        "outer_rewrites_N": outer["rewrites_unknown_N"],
        "outer_rewrites_T": outer["rewrites_unknown_T"],
        "outer_history": outer["history"],
        "omitted": omitted_named(),
        "note": (
            "40×40：39 物种物流 + 1 个 T。干料质量与 C/H/O 原子流对齐后"
            "绝热根约 565 K，不是气化炉工况。外层只刷新热解 q。"
            "生产 single-shot 未改。壁面 Q 省略。"
        ),
    }
