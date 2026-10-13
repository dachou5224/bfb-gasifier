#!/usr/bin/env python3
"""inner#2 总步长 Δx 分组掩码：找触发 bed2 R_gas 跳跃的 DOF（handoff §5.13）。

流程：钩住 outer abgleich，在 inner#2 入口 pack x0、出口 pack x1，
对 Δx=x1-x0 做按床/物种/相/固体/T 掩码回放，测量 bed2 R_H2O。

用法：
  python3 scripts/audit_bed2_flip_inner2_dx_mask.py
"""

from __future__ import annotations

import argparse
from typing import Any

import numpy as np

from _lu_harness import (
    build_phase2_reactor,
    data_path,
    phase2_solve_kwargs,
    run_phase2_init,
    timed_phase2_solve,
    write_json,
)

ensure_repo_path = __import__("_lu_harness").ensure_repo_path
ensure_repo_path()

from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX
from src.solvers.nr_indexing import (
    cell_offsets,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_solid_var,
    n_var,
    pack_reactor,
    temperature_offset,
    unpack_reactor,
)
from src.workflow.steps import outer_abgleich_step as oas

ISOMORPH = {
    "extent_limiters_enabled_thesis": False,
    "vorab_transport_x0_fast_oxidation_closure_thesis": True,
    "nr_line_search_fast_oxidation_projection_thesis": True,
    "nr_outer_preinner_bed1_phase_split_thesis": True,
    "nr_include_pyrolysis_gas_species_in_unknowns_thesis": True,
    "nr_refresh_hydrodynamics_on_accepted_step_thesis": False,
    "nr_line_search_syngas_collapse_guard_thesis": True,
    "nr_outer_preinner_upper_co_local2x2_thesis": True,
    "nr_outer_preinner_upper_dense_syngas_absorb_thesis": True,
    "nr_outer_preinner_upper_dense_absorb_species_thesis": ("CO",),
    "nr_outer_preinner_upper_dense_absorb_fraction_thesis": 0.3,
    "nr_outer_preinner_upper_dense_absorb_mode_thesis": "phase_share",
    "nr_outer_preinner_upper_dense_absorb_post_gate_thesis": True,
}


def _bed2_R(cells: list[Any]) -> dict[str, float]:
    c = cells[2]
    c.residuals()
    out: dict[str, float] = {"T": float(c.T)}
    for sp in ("H2O", "H2", "O2", "CO"):
        j = int(GAS_SPECIES_INDEX[sp])
        out[f"R_{sp}"] = float(c.R_gas_d[j] + c.R_gas_b[j])
    return out


def _invalidate(cells: list[Any]) -> None:
    for c in cells:
        c._thermo_cache_valid = False
        c._hydro_cache_valid = False


def _dof_maps(cells: list[Any]) -> dict[str, Any]:
    offs = cell_offsets(cells)
    bed_slices: list[slice] = []
    gas_idx_all: list[int] = []
    solid_idx_all: list[int] = []
    t_idx_all: list[int] = []
    # species -> packed indices (dense+bubble for two-phase)
    sp_idx: dict[str, list[int]] = {sp: [] for sp in GAS_SPECIES}
    dense_idx: list[int] = []
    bubble_idx: list[int] = []
    for i, cell in enumerate(cells):
        base = int(offs[i])
        nv = int(n_var(cell))
        bed_slices.append(slice(base, base + nv))
        ng = int(n_gas_var(cell))
        ns = int(n_solid_var(cell))
        gas_idx_all.extend(range(base, base + ng))
        solid_idx_all.extend(range(base + ng, base + ng + ns))
        if ns + ng < nv:
            t_idx_all.append(base + temperature_offset(cell))
        gidx = main_nr_gas_idx_for_cell(cell)
        # two-phase: [N_d[gidx], N_b[gidx]]
        if ng == 2 * int(gidx.size):
            for k, gi in enumerate(gidx):
                sp = GAS_SPECIES[int(gi)]
                sp_idx[sp].append(base + k)
                sp_idx[sp].append(base + int(gidx.size) + k)
                dense_idx.append(base + k)
                bubble_idx.append(base + int(gidx.size) + k)
        elif ng == int(gidx.size):
            for k, gi in enumerate(gidx):
                sp = GAS_SPECIES[int(gi)]
                sp_idx[sp].append(base + k)
                dense_idx.append(base + k)
    return {
        "bed_slices": bed_slices,
        "gas": np.asarray(gas_idx_all, dtype=np.intp),
        "solid": np.asarray(solid_idx_all, dtype=np.intp),
        "T": np.asarray(t_idx_all, dtype=np.intp),
        "species": {k: np.asarray(v, dtype=np.intp) for k, v in sp_idx.items() if v},
        "dense": np.asarray(dense_idx, dtype=np.intp),
        "bubble": np.asarray(bubble_idx, dtype=np.intp),
        "n": int(offs[-1]),
    }


def _mask_from_indices(n: int, idxs: np.ndarray) -> np.ndarray:
    m = np.zeros(n, dtype=np.float64)
    if idxs.size:
        m[idxs] = 1.0
    return m


def _eval_masked(
    cells: list[Any],
    x0: np.ndarray,
    dx: np.ndarray,
    mask: np.ndarray,
    apply_bc: Any,
) -> dict[str, float]:
    x = x0 + mask * dx
    unpack_reactor(x, cells)
    _invalidate(cells)
    apply_bc()
    _invalidate(cells)
    return _bed2_R(cells)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-global-iter", type=int, default=3)
    args = ap.parse_args()

    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=dict(ISOMORPH))
    run_phase2_init(reactor)
    apply_bc = reactor._apply_all_bc_for_nr

    capture: dict[str, Any] = {}
    orig = oas.run_outer_abgleich_for_global_nr
    ctr = {"inner": 0}

    def wrapped(**kwargs: Any) -> Any:
        cells_kw = kwargs["cells"]
        solve_inner_fn = kwargs["solve_inner_fn"]

        def wrap_inner(
            max_iter: int, lambda_seed: float, jac_mode: str, jac_lag: int
        ) -> dict:
            ctr["inner"] += 1
            i = ctr["inner"]
            if i == 2:
                _invalidate(cells_kw)
                apply_bc()
                capture["cells"] = cells_kw
                capture["x0"] = pack_reactor(cells_kw).copy()
                capture["R0"] = _bed2_R(cells_kw)
            result = solve_inner_fn(max_iter, lambda_seed, jac_mode, jac_lag)
            if i == 2:
                _invalidate(cells_kw)
                apply_bc()
                capture["x1"] = pack_reactor(cells_kw).copy()
                capture["R1"] = _bed2_R(cells_kw)
                capture["inner_result"] = {
                    "n_newton": int(
                        result.get("n_newton_iters_attempted", result.get("n_iter", -1))
                    ),
                    "rms": float(result.get("rms_scaled_final", float("nan"))),
                    "lambdas": list(result.get("accepted_lambda_history") or []),
                }
            return result

        kw = dict(kwargs)
        kw["solve_inner_fn"] = wrap_inner
        return orig(**kw)

    oas.run_outer_abgleich_for_global_nr = wrapped  # type: ignore[assignment]
    try:
        result, wall = timed_phase2_solve(
            reactor, **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter))
        )
    finally:
        oas.run_outer_abgleich_for_global_nr = orig  # type: ignore[assignment]

    if "x0" not in capture or "x1" not in capture or "cells" not in capture:
        raise RuntimeError("inner#2 未捕获到状态")

    cells = capture["cells"]
    x0 = np.asarray(capture["x0"], dtype=np.float64)
    x1 = np.asarray(capture["x1"], dtype=np.float64)
    dx = x1 - x0
    maps = _dof_maps(cells)
    n = int(maps["n"])
    if x0.size != n:
        raise RuntimeError(f"pack size mismatch x0={x0.size} offs={n}")

    trials: list[dict[str, Any]] = []

    def run_trial(label: str, mask: np.ndarray) -> None:
        R = _eval_masked(cells, x0, dx, mask, apply_bc)
        row = {
            "label": label,
            "mask_nnz": int(np.count_nonzero(mask)),
            "dx_norm_masked": float(np.linalg.norm(mask * dx)),
            **R,
            "flipped": bool(R["R_H2O"] > 20.0),
        }
        trials.append(row)
        flag = "FLIP" if row["flipped"] else "ok"
        print(
            f"[{flag}] {label:28s} R_H2O={R['R_H2O']:+8.2f} "
            f"R_H2={R['R_H2']:+8.2f} R_O2={R['R_O2']:+8.2f} "
            f"|dx_m|={row['dx_norm_masked']:.3g}"
        )

    # 基线
    run_trial("none (x0)", np.zeros(n))
    run_trial("full (x1)", np.ones(n))

    # 按床（含 freeboard/side）
    for i, cell in enumerate(cells):
        kind = str(getattr(cell, "cell_type", "bed"))
        m = np.zeros(n)
        m[maps["bed_slices"][i]] = 1.0
        run_trial(f"cell{i}_{kind}_only", m)

    # 排除关键床
    for i in (1, 2, 6):
        m = np.ones(n)
        m[maps["bed_slices"][i]] = 0.0
        run_trial(f"all_except_cell{i}", m)

    # 类型
    run_trial("gas_only", _mask_from_indices(n, maps["gas"]))
    run_trial("solid_only", _mask_from_indices(n, maps["solid"]))
    run_trial("T_only", _mask_from_indices(n, maps["T"]))
    run_trial("dense_only", _mask_from_indices(n, maps["dense"]))
    run_trial("bubble_only", _mask_from_indices(n, maps["bubble"]))

    # 物种（全床）
    for sp in ("O2", "H2O", "H2", "CO", "CO2", "CH4", "N2"):
        idxs = maps["species"].get(sp)
        if idxs is None or idxs.size == 0:
            continue
        run_trial(f"sp_{sp}_allcells", _mask_from_indices(n, idxs))

    # bed2 物种细分
    b2 = maps["bed_slices"][2]
    gidx = main_nr_gas_idx_for_cell(cells[2])
    base2 = int(b2.start)
    ng2 = int(gidx.size)
    two_phase = int(n_gas_var(cells[2])) == 2 * ng2
    for sp in ("O2", "H2O", "H2", "CO", "CO2"):
        pos = [k for k, gi in enumerate(gidx) if GAS_SPECIES[int(gi)] == sp]
        if not pos:
            continue
        idxs = []
        for k in pos:
            idxs.append(base2 + k)
            if two_phase:
                idxs.append(base2 + ng2 + k)
        run_trial(f"bed2_{sp}", _mask_from_indices(n, np.asarray(idxs, dtype=np.intp)))

    # bed1 物种（上游）
    b1 = maps["bed_slices"][1]
    gidx1 = main_nr_gas_idx_for_cell(cells[1])
    base1 = int(b1.start)
    ng1 = int(gidx1.size)
    two_phase1 = int(n_gas_var(cells[1])) == 2 * ng1
    for sp in ("O2", "H2O", "H2", "CO"):
        pos = [k for k, gi in enumerate(gidx1) if GAS_SPECIES[int(gi)] == sp]
        if not pos:
            continue
        idxs = []
        for k in pos:
            idxs.append(base1 + k)
            if two_phase1:
                idxs.append(base1 + ng1 + k)
        run_trial(f"bed1_{sp}", _mask_from_indices(n, np.asarray(idxs, dtype=np.intp)))

    # 组合：全气但禁 O2；全气但禁 bed2
    m = _mask_from_indices(n, maps["gas"])
    o2 = maps["species"].get("O2")
    if o2 is not None:
        m2 = m.copy()
        m2[o2] = 0.0
        run_trial("gas_except_O2", m2)
    m3 = m.copy()
    m3[maps["bed_slices"][2]] = 0.0
    run_trial("gas_except_bed2", m3)

    # beds 0-2 gas only / beds 6-9 gas only
    m_lo = np.zeros(n)
    for i in (0, 1, 2):
        m_lo[maps["bed_slices"][i]] = 1.0
    m_lo = m_lo * _mask_from_indices(n, maps["gas"])
    run_trial("gas_beds0to2", m_lo)
    m_hi = np.zeros(n)
    for i in (6, 7, 8, 9):
        m_hi[maps["bed_slices"][i]] = 1.0
    m_hi = m_hi * _mask_from_indices(n, maps["gas"])
    run_trial("gas_beds6to9", m_hi)

    # 恢复末态以免误导
    unpack_reactor(x1, cells)
    apply_bc()

    exit_co = float((result.get("bed_exit_gas") or {}).get("CO", float("nan")))
    flipped = [t for t in trials if t["flipped"] and t["label"] not in ("full (x1)",)]
    necessary = [
        t["label"]
        for t in trials
        if t["label"].startswith("all_except_") and not t["flipped"]
    ]

    out = {
        "note": "inner#2 Δx group mask → bed2 R_H2O; flip if R_H2O>20",
        "wall_s": float(wall),
        "exit_CO": exit_co,
        "R0": capture["R0"],
        "R1": capture["R1"],
        "inner_result": capture.get("inner_result"),
        "dx_norm": float(np.linalg.norm(dx)),
        "dx_max_abs": float(np.max(np.abs(dx))),
        "trials": trials,
        "flip_inducing_masks": [t["label"] for t in flipped],
        "necessary_beds_from_except": necessary,
    }
    path = data_path("phase2_bed2_flip_inner2_dx_mask.json")
    write_json(path, out)
    print(
        f"wrote {path} flip_inducing={out['flip_inducing_masks']} "
        f"necessary_except={necessary}"
    )


if __name__ == "__main__":
    main()
