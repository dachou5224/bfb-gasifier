#!/usr/bin/env python3
"""Phase2 LU 出口 CO/CO₂ 结构偏差：多因子根因探针（床内 + 自由板 + 裁定）。

单路径 ``reactor.solve()``（禁止 apply_kinetics）。默认 8-outer 快探针；
差分环每次只动一个因子，结果写入 data/，不作为验收标定。

Usage
-----
    python3 scripts/audit_phase2_coco2_structural_rca.py
    python3 scripts/audit_phase2_coco2_structural_rca.py --max-global-iter 8 --skip-variants
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from _lu_harness import (
    REF_CO,
    REF_CO2,
    build_phase2_reactor,
    co_co2_score,
    data_path,
    phase2_solve_kwargs,
    repo_root,
    timed_phase2_solve,
    write_json,
)
from audit_bed_co_co2_budget_lu import _analyze_cell_co_co2, _sum_terms
from src.core.composition import wet_to_dry_mole_fractions
from src.core.connectivity_graph import build_connectivity_from_reactor_config
from src.core.freeboard_bridge import refresh_explicit_freeboard_transport_from_closure
from src.core.species import GAS_SPECIES_INDEX as idx

REPO_ROOT = repo_root()

MAJOR = ("CO", "CO2", "H2", "O2", "H2O", "N2", "CH4")

# Bild 7.4 wet anchors (validation_cases axial_profiles)
BILD74_WET = {
    0.05: {"CO": 0.02, "CO2": 0.035, "H2": 0.01, "H2O": 0.28, "T_K": None},
    0.15: {"CO": 0.075, "CO2": 0.065, "H2": 0.06, "H2O": 0.245, "T_K": None},
    0.30: {"CO": 0.115, "CO2": 0.095, "H2": 0.105, "H2O": 0.21, "T_K": 1080.0},
    0.45: {"CO": 0.13, "CO2": 0.105, "H2": 0.12, "H2O": 0.195, "T_K": 1060.0},
    0.60: {"CO": 0.11, "CO2": 0.125, "H2": 0.105, "H2O": 0.165, "T_K": 1240.0},
    0.80: {"CO": 0.125, "CO2": 0.115, "H2": 0.115, "H2O": 0.17, "T_K": None},
    1.00: {"CO": 0.13, "CO2": 0.11, "H2": 0.12, "H2O": 0.175, "T_K": 1100.0},
}
BILD74_T = {
    0.00: 750.0,
    0.10: 1020.0,
    0.20: 1100.0,
    0.45: 1060.0,
    0.60: 1240.0,
    0.90: 1130.0,
    1.00: 1100.0,
}
TRUE_DRY = {"CO": 0.157, "CO2": 0.133}


def _y_dict(arr: np.ndarray) -> dict[str, float]:
    return {sp: float(arr[idx[sp]]) for sp in MAJOR}


def _N_dict(arr: np.ndarray) -> dict[str, float]:
    return {sp: float(arr[idx[sp]]) for sp in ("CO", "CO2", "H2", "O2", "H2O", "N2")}


def _collect_bed_axial(reactor: Any) -> list[dict[str, Any]]:
    H_bed = float(reactor.config.H_bed)
    rows: list[dict[str, Any]] = []
    for i, cell in enumerate(reactor.cells):
        if str(getattr(cell, "cell_type", "bed")) != "bed":
            continue
        y_b = cell._mole_fractions("b")
        y_d = cell._mole_fractions("d")
        y_c = cell._mole_fractions("combined")
        N_b = np.asarray(cell.N_b, dtype=float)
        N_d = np.asarray(cell.N_d, dtype=float)
        n_ex = np.asarray(cell.N_ex, dtype=float)
        xi = float(cell.geo.h_center / max(H_bed, 1e-12))
        E = float(getattr(cell, "E", 0.0) or 0.0)
        rows.append(
            {
                "i": i,
                "xi": xi,
                "T_K": float(cell.T),
                "E_MW": E / 1e6,
                "eps_b": float(cell.eps_b),
                "K_bd": float(cell.K_bd),
                "V_b": float(cell.V_b),
                "V_d": float(cell.V_d),
                "y_combined": _y_dict(y_c),
                "y_bubble": _y_dict(y_b),
                "y_dense": _y_dict(y_d),
                "N_bubble": _N_dict(N_b),
                "N_dense": _N_dict(N_d),
                "N_ex": _N_dict(n_ex),
                "N_tot_CO": float(N_b[idx["CO"]] + N_d[idx["CO"]]),
                "N_tot_CO2": float(N_b[idx["CO2"]] + N_d[idx["CO2"]]),
                "N_tot_O2": float(N_b[idx["O2"]] + N_d[idx["O2"]]),
            }
        )
    return rows


def _bild_T_at(xi: float) -> float | None:
    keys = sorted(BILD74_T)
    if xi <= keys[0]:
        return BILD74_T[keys[0]]
    if xi >= keys[-1]:
        return BILD74_T[keys[-1]]
    for a, b in zip(keys[:-1], keys[1:]):
        if a <= xi <= b:
            t = (xi - a) / max(b - a, 1e-12)
            return BILD74_T[a] * (1.0 - t) + BILD74_T[b] * t
    return None


def _source_budget(reactor: Any) -> dict[str, Any]:
    rows = [_analyze_cell_co_co2(cell) for cell in reactor.cells]
    H_bed = float(reactor.config.H_bed)
    for i, row in enumerate(rows):
        row["cell_index"] = i
        row["xi_over_H"] = float(reactor.cells[i].geo.h_center / max(H_bed, 1e-12))
    co_sum = _sum_terms(rows, "co_terms")
    co2_sum = _sum_terms(rows, "co2_terms")
    # segment integrals: bed0, bed1-2, bed3-9
    segs = {
        "bed0": [0],
        "bed1_2": [1, 2],
        "bed3_top": list(range(3, len(rows))),
    }
    seg_out: dict[str, Any] = {}
    for name, idxs in segs.items():
        sel = [rows[j] for j in idxs if j < len(rows)]
        seg_out[name] = {
            "co_terms": _sum_terms(sel, "co_terms"),
            "co2_terms": _sum_terms(sel, "co2_terms"),
            "R5_CO_mol_s": sum(
                float(r["co_terms"]["R5b"] + r["co_terms"]["R5d"]) for r in sel
            ),
            "R2_CO_mol_s": sum(float(r["co_terms"]["R2"]) for r in sel),
            "R4_CO_mol_s": sum(float(r["co_terms"]["R4"]) for r in sel),
            "R8_CO_mol_s": sum(float(r["co_terms"]["R8"]) for r in sel),
        }
    return {
        "cells": rows,
        "bed_integrated_co_terms_mol_s": co_sum,
        "bed_integrated_co2_terms_mol_s": co2_sum,
        "segments": seg_out,
    }


def _collect_fb(reactor: Any, result: dict[str, Any]) -> dict[str, Any]:
    cfg = reactor.config
    H_bed = float(cfg.H_bed)
    H_tot = H_bed + float(cfg.H_freeboard)
    graph = build_connectivity_from_reactor_config(cfg, reactor=reactor)
    sec_mode = str(getattr(cfg, "freeboard_secondary_injection_mode", "lumped"))
    sec_o2 = float(getattr(cfg, "freeboard_secondary_O2_mol_s", 0.0))
    sec_xi = getattr(cfg, "freeboard_secondary_injection_xi", None)

    nr_fb: list[dict[str, Any]] = []
    for j, cell in enumerate(getattr(reactor, "freeboard_cells", []) or []):
        y = cell._mole_fractions("combined")
        z = float(cell.geo.h_center)
        xi_tot = z / max(H_tot, 1e-12)
        nr_fb.append(
            {
                "i": j,
                "z_m": z,
                "xi_reactor": xi_tot,
                "T_K": float(cell.T),
                "y": _y_dict(y),
                "N": _N_dict(np.asarray(cell.N_d + cell.N_b, dtype=float)),
            }
        )

    ana = refresh_explicit_freeboard_transport_from_closure(
        reactor, preserve_gas_state=False
    ) or {}
    prof = (ana.get("profiles") or {}).get("wet_gas") or {}
    ana_rows: list[dict[str, Any]] = []
    n_prof = len(prof.get("CO", [])) if isinstance(prof, dict) else 0
    z_arr = prof.get("z_m") or prof.get("z") or []
    T_arr = prof.get("T") or prof.get("T_K") or []
    for k in range(n_prof):
        z = float(z_arr[k]) if k < len(z_arr) else float("nan")
        # freeboard-local z → reactor xi
        xi_r = (H_bed + z) / max(H_tot, 1e-12) if z == z else float("nan")
        row = {
            "i": k,
            "z_m_fb": z,
            "xi_reactor": xi_r,
            "T_K": float(T_arr[k]) if k < len(T_arr) else None,
            "y_wet": {sp: float(prof[sp][k]) for sp in MAJOR if sp in prof},
        }
        ana_rows.append(row)

    # Bild 7.4 freeboard window ξ=0.6–1.0
    bild_fb = {str(k): v for k, v in BILD74_WET.items() if k >= 0.6}
    # model bed-top / exit dry
    bed_exit = result.get("bed_exit_gas_dry") or {}
    exit_dry = result.get("exit_gas_dry") or {}

    # dip/recover metrics on NR fb CO/CO2
    def _extrema(rows: list[dict], key: str) -> dict[str, Any]:
        if not rows:
            return {}
        ys = [float(r["y"][key]) for r in rows]
        return {
            "min": min(ys),
            "max": max(ys),
            "first": ys[0],
            "last": ys[-1],
            "delta_last_first": ys[-1] - ys[0],
            "has_dip_then_recover": (min(ys) < ys[0] - 1e-4) and (ys[-1] > min(ys) + 1e-4),
        }

    return {
        "config": {
            "secondary_injection_mode": sec_mode,
            "secondary_O2_mol_s": sec_o2,
            "secondary_N2_mol_s": float(getattr(cfg, "freeboard_secondary_N2_mol_s", 0.0)),
            "secondary_injection_xi": sec_xi,
            "graph_secondary_injection_active": bool(graph.secondary_injection_active),
            "lumped_excluded_from_nr_graph": sec_mode == "lumped",
            "n_freeboard_cells": int(getattr(cfg, "n_freeboard_cells", 0)),
            "trajectory_model": str(getattr(cfg, "freeboard_trajectory_model", "")),
        },
        "bed_exit_dry": {k: float(bed_exit.get(k, 0.0) or 0.0) for k in ("CO", "CO2", "H2", "O2")},
        "reactor_exit_dry": {k: float(exit_dry.get(k, 0.0) or 0.0) for k in ("CO", "CO2", "H2", "O2")},
        "delta_fb_dry_CO": float(exit_dry.get("CO", 0.0) or 0.0) - float(bed_exit.get("CO", 0.0) or 0.0),
        "delta_fb_dry_CO2": float(exit_dry.get("CO2", 0.0) or 0.0) - float(bed_exit.get("CO2", 0.0) or 0.0),
        "nr_freeboard_cells": nr_fb,
        "nr_CO_extrema": _extrema(nr_fb, "CO"),
        "nr_CO2_extrema": _extrema(nr_fb, "CO2"),
        "analytical_freeboard_rows": ana_rows,
        "analytical_exit_wet": ana.get("exit_gas"),
        "analytical_secondary_segment": ana.get("secondary_injection_segment"),
        "bild74_fb_wet_anchors": bild_fb,
        "bild74_expected": {
            "at_xi_0.6": "CO/H2 dip, CO2 spike, T~1240 K, O2~0.012",
            "after_0.6": "mild CO/H2 recover, CO2 softens toward exit 0.11 wet",
            "delta_CO2_wet_0.6_to_1.0": 0.11 - 0.125,
        },
    }


def _hypothesis_signals(
    axial: list[dict[str, Any]],
    budget: dict[str, Any],
    fb: dict[str, Any],
) -> dict[str, Any]:
    # H1: lower-bed bubble R5 locks CO2 inventory (segment-resolved)
    seg = budget["segments"]
    r5_0 = float(seg["bed0"]["R5_CO_mol_s"])
    r5_12 = float(seg["bed1_2"]["R5_CO_mol_s"])
    r2_12 = float(seg["bed1_2"]["R2_CO_mol_s"])
    r4_12 = float(seg["bed1_2"]["R4_CO_mol_s"])
    r5_all = float(budget["bed_integrated_co_terms_mol_s"].get("R5b", 0.0)) + float(
        budget["bed_integrated_co_terms_mol_s"].get("R5d", 0.0)
    )
    co2_lock_cell = None
    for row in axial:
        if row["i"] >= 1 and abs(row["y_combined"]["CO2"] - axial[-1]["y_combined"]["CO2"]) < 0.005:
            co2_lock_cell = row["i"]
            break
    y0 = axial[0]["y_combined"] if axial else {}
    y2 = axial[2]["y_combined"] if len(axial) > 2 else {}
    ytop = axial[-1]["y_combined"] if axial else {}

    # H2: exchange of O2/CO in lower bed (0–2)
    nex_rows = []
    for row in axial:
        if row["i"] in (0, 1, 2):
            nex_rows.append(
                {
                    "i": row["i"],
                    "N_ex_O2": row["N_ex"]["O2"],
                    "N_ex_CO": row["N_ex"]["CO"],
                    "y_b_O2": row["y_bubble"]["O2"],
                    "y_d_O2": row["y_dense"]["O2"],
                    "y_b_CO": row["y_bubble"]["CO"],
                    "y_d_CO": row["y_dense"]["CO"],
                    "K_bd": row["K_bd"],
                    "eps_b": row["eps_b"],
                }
            )

    # H3: T vs Bild (bed ξ → reactor ξ ≈ ξ_bed * H_bed/H_tot ~ 0.45)
    H_bed_frac = 0.45
    t_cmp = []
    for row in axial:
        xi_r = float(row["xi"]) * H_bed_frac
        t_lit = _bild_T_at(xi_r)
        t_cmp.append(
            {
                "i": row["i"],
                "xi_bed": row["xi"],
                "xi_reactor_approx": xi_r,
                "T_model": row["T_K"],
                "T_bild": t_lit,
                "dT": None if t_lit is None else row["T_K"] - t_lit,
            }
        )

    co2_rise = float(y2.get("CO2", 0) - y0.get("CO2", 0)) if y0 else 0.0
    h1 = {
        "claim": "床内过氧化：下部气泡相 R5 建立过高 CO₂ 库存，中上部落定后难回落",
        "R5_CO_bed0_mol_s": r5_0,
        "R5_CO_bed1_2_mol_s": r5_12,
        "R2_CO_bed1_2_mol_s": r2_12,
        "R4_CO_bed1_2_mol_s": r4_12,
        "R5_CO_bed_integrated_mol_s": r5_all,
        "R5_bed0_share_of_integrated": abs(r5_0) / max(abs(r5_all), 1e-9),
        "R5_over_gasif_ratio_bed1_2": abs(r5_12) / max(abs(r2_12) + abs(r4_12), 1e-9),
        "y_CO_bed0_to_bed2": {
            "bed0": y0.get("CO"),
            "bed2": y2.get("CO"),
            "delta": (None if not y0 else float(y2.get("CO", 0) - y0.get("CO", 0))),
        },
        "y_CO2_bed0_to_bed2": {
            "bed0": y0.get("CO2"),
            "bed2": y2.get("CO2"),
            "delta": co2_rise,
        },
        "y_top_wet": ytop,
        "co2_near_top_from_cell": co2_lock_cell,
        "supported": abs(r5_all) > 2.0 and (co2_rise > 0.02 or abs(r5_0) > abs(r5_12)),
        "note": (
            "本探针：R5 积分以 bed0 气泡相为主；bed1–2 的 R2 可部分回补 CO，"
            "但 CO₂ 湿基在 bed1 附近已抬升并沿床锁定。"
        ),
    }
    h2 = {
        "claim": "两相分流：K_bd/ε_b 使 O₂/CO 在可氧化通道上错配",
        "lower_bed_exchange": nex_rows,
        "supported_coupling": bool(nex_rows)
        and any(abs(r["N_ex_O2"]) > 0.05 or abs(r["N_ex_CO"]) > 0.5 for r in nex_rows),
        "note": "dense_frac 差分在 8-outer 下易病理塌缩，基线 N_ex(CO) 大通量作耦合证据",
    }
    hot = [r for r in t_cmp if r["i"] in (0, 1, 2)]
    h3 = {
        "claim": "传热/温度场：下部过热强化 R5，吸热气化窗口偏短",
        "T_compare_lower_bed": hot,
        "mean_dT_bed012": float(np.mean([r["dT"] for r in hot if r["dT"] is not None]))
        if hot
        else None,
        "supported": bool(hot)
        and all((r["dT"] or 0) > 80.0 for r in hot if r["i"] in (0, 1)),
    }
    bed_co2 = float(fb["bed_exit_dry"].get("CO2", 0))
    h4 = {
        "claim": "自由板缺位：未再现二级风 dip+恢复；lumped 二次风未进 NR 图",
        "lumped_excluded_from_nr_graph": fb["config"]["lumped_excluded_from_nr_graph"],
        "graph_secondary_active": fb["config"]["graph_secondary_injection_active"],
        "delta_fb_dry_CO2": fb["delta_fb_dry_CO2"],
        "nr_CO_has_dip_recover": fb["nr_CO_extrema"].get("has_dip_then_recover"),
        "nr_CO2_delta": fb["nr_CO2_extrema"].get("delta_last_first"),
        "bed_exit_already_high_CO2": bed_co2 > 0.145,
        "fb_alone_insufficient": (bed_co2 - 0.015) > TRUE_DRY["CO2"],
        "supported": fb["config"]["lumped_excluded_from_nr_graph"]
        and abs(fb["delta_fb_dry_CO2"]) < 0.005
        and not fb["nr_CO_extrema"].get("has_dip_then_recover", False),
        "role": (
            "结构债/文献同构缺口；因床顶 CO₂ 已高，FB 即使完美复现 Bild Δ≈−0.015 "
            "仍不够收到真干基 → 不能单独作为下一刀主路径"
        ),
    }
    return {"H1": h1, "H2": h2, "H3": h3, "H4": h4}


def _run_variant(
    label: str,
    *,
    max_global_iter: int,
    cfg_overrides: dict[str, object] | None = None,
    cell_hook: Any | None = None,
) -> dict[str, Any]:
    reactor = build_phase2_reactor(ladder_level="p2", cfg_overrides=cfg_overrides)
    if cell_hook is not None:
        cell_hook(reactor)
    out, wall = timed_phase2_solve(
        reactor, **phase2_solve_kwargs(max_global_iter=max_global_iter, verbose=False)
    )
    bed = out.get("bed_exit_gas_dry") or {}
    ex = out.get("exit_gas_dry") or {}
    return {
        "label": label,
        "wall_s": float(wall),
        "converged": bool(out.get("converged")),
        "rms_E": out.get("rms_scaled_energy_final"),
        "bed_exit_dry": {k: float(bed.get(k, 0.0) or 0.0) for k in ("CO", "CO2", "H2", "O2")},
        "exit_dry": {k: float(ex.get(k, 0.0) or 0.0) for k in ("CO", "CO2", "H2", "O2")},
        "score_vs_harness": co_co2_score(
            float(ex.get("CO", 0.0) or 0.0), float(ex.get("CO2", 0.0) or 0.0)
        ),
        "score_vs_true_dry": abs(float(ex.get("CO", 0.0) or 0.0) - TRUE_DRY["CO"]) / TRUE_DRY["CO"]
        + abs(float(ex.get("CO2", 0.0) or 0.0) - TRUE_DRY["CO2"]) / TRUE_DRY["CO2"],
        "cfg": dict(cfg_overrides or {}),
    }


def _verdict(
    signals: dict[str, Any],
    fb: dict[str, Any],
    variants: list[dict[str, Any]],
    exit_dry: dict[str, float],
) -> dict[str, Any]:
    h1, h2, h3, h4 = signals["H1"], signals["H2"], signals["H3"], signals["H4"]
    var_by = {v["label"]: v for v in variants}

    # interpret variants
    r5_lo = var_by.get("diag_r5_scale_0.10")
    r5_hi = var_by.get("diag_r5_scale_0.50")
    dense_lo = var_by.get("diag_dense_frac_0.25")
    dense_hi = var_by.get("diag_dense_frac_0.75")
    base = var_by.get("baseline")

    r5_causal = False
    if base and r5_lo:
        r5_causal = float(r5_lo["exit_dry"]["CO2"]) < float(base["exit_dry"]["CO2"]) - 0.005

    exchange_causal = False
    if base and dense_lo and dense_hi:
        exchange_causal = abs(
            float(dense_lo["exit_dry"]["CO2"]) - float(dense_hi["exit_dry"]["CO2"])
        ) > 0.008

    factors = [
        {
            "id": "H1_bed_overox",
            "evidence": h1,
            "fix_surface": "下部床气泡相 R5 / O₂–CO 相区分流（parity 先于全局 r5_scale）",
            "status": "PRIMARY" if h1.get("supported") else "WEAK",
        },
        {
            "id": "H2_phase_split",
            "evidence": {
                "supported_coupling": h2.get("supported_coupling"),
                "exchange_causal_variant": exchange_causal,
                "dense_frac_o8_note": "8-outer dense_frac 变体易病理，不作主证据",
            },
            "fix_surface": "K_bd / ε_b / 入口相分配（两相 O₂/CO）",
            "status": "COUPLING" if h2.get("supported_coupling") else "WEAK",
        },
        {
            "id": "H3_temperature",
            "evidence": h3,
            "fix_surface": "下部 T 场（禁止抬 Tin；先查围栏/seed/热损同构）",
            "status": "COUPLING" if h3.get("supported") else "SECONDARY",
        },
        {
            "id": "H4_freeboard",
            "evidence": h4,
            "fix_surface": "lumped→NR 二次风接线 / Bild 7.4 dip+recover 同构（并行结构债）",
            "status": (
                "PARALLEL_STRUCTURAL_GAP"
                if h4.get("supported") and h4.get("fb_alone_insufficient")
                else ("STRUCTURAL_GAP" if h4.get("supported") else "WEAK")
            ),
        },
        {
            "id": "H5_basis",
            "evidence": {
                "harness_REF": {"CO": REF_CO, "CO2": REF_CO2},
                "true_dry": TRUE_DRY,
                "model_exit": exit_dry,
            },
            "fix_surface": "验收口径改真干基 0.157/0.133；勿用 0.13/0.11 当干基",
            "status": "METROLOGY_GATE",
        },
    ]

    # Capacity gate from plan Phase C: bed-top already high → FB cannot be sole knife
    if h4.get("fb_alone_insufficient") or h1.get("supported"):
        next_knife = {
            "path": "bed_inner_R5_O2_split",
            "rationale": (
                "床顶干基 CO₂ 已≈0.15+，即使 FB 复现 Bild ΔCO₂≈−0.015 仍到不了真干基 0.133；"
                "源汇积分显示气泡 R5（本探针主要在 bed0，bed1–2 仍有可观 R5）建立 CO₂ 库存。"
                "主刀：下部床 R5/O₂ 相区分流 parity。自由板二次风 NR 同构为并行结构债，不能单独背锅。"
            ),
            "do_not": [
                "不要用全局 r5_scale 当验收解（先清 bubble R5 / N_ex parity）",
                "不要抬 T_inlet",
                "不要把 bed0 焓地板当 CO/CO₂ 主因重开",
                "不要假设只修 FB 接线即可收到真干基",
            ],
            "falsifiers": [
                "诊断缩放 bed0–2 R5_bubble 后床顶 CO₂ 不降 → 否证 H1 主因，改查 H2/H3",
                "床顶已贴近真干基而反应器出口仍偏 → 改主刀为 FB 二次风同构",
                "仅修 K_bd/入口相分配即同时改善 CO↑CO₂↓ → 升级 H2 为并列主刀",
            ],
        }
    else:
        next_knife = {
            "path": "freeboard_secondary_nr_isomorphism",
            "rationale": "床顶已近真干基；出口偏差由自由板缺位主导",
            "do_not": ["不要用 r*_scale 掩盖"],
            "falsifiers": ["NR 图接入 secondary 后仍无 dip → 工程桥接 B 档"],
        }

    return {
        "factors": factors,
        "variant_summary": variants,
        "r5_scale_moves_CO2": r5_causal,
        "dense_frac_moves_CO2": exchange_causal,
        "next_knife": next_knife,
        "success_criteria": {
            "where_formed": (
                "下部床气泡相 R5（本探针 bed0 约占积分 70%+，bed1–2 续烧）建立过高 CO₂；"
                "中上床湿基 CO₂ 锁定；床顶≈出口"
            ),
            "fb_inert_is": (
                "文献冲突 + lumped 二次风未进 NR 连通图（解析路径可注入，但 NR 显式链近乎平坦、无 dip+recover）"
            ),
            "unique_next_path": next_knife["path"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-global-iter", type=int, default=8)
    parser.add_argument("--skip-variants", action="store_true")
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=data_path("phase2_coco2_structural"),
    )
    args = parser.parse_args()
    prefix = Path(args.out_prefix)
    out_bed = Path(str(prefix) + "_bed_decomp.json")
    out_fb = Path(str(prefix) + "_fb_link.json")
    out_verdict = Path(str(prefix) + "_verdict.json")

    print("=== baseline p2 solve ===", flush=True)
    reactor = build_phase2_reactor(ladder_level="p2")
    assert float(getattr(reactor.config, "r5_scale", 1.0)) == 1.0
    result, wall = timed_phase2_solve(
        reactor,
        **phase2_solve_kwargs(max_global_iter=int(args.max_global_iter), verbose=False),
    )
    axial = _collect_bed_axial(reactor)
    budget = _source_budget(reactor)
    fb = _collect_fb(reactor, result)
    exit_dry = {
        k: float((result.get("exit_gas_dry") or {}).get(k, 0.0) or 0.0)
        for k in ("CO", "CO2", "H2", "O2")
    }
    bed_exit = {
        k: float((result.get("bed_exit_gas_dry") or {}).get(k, 0.0) or 0.0)
        for k in ("CO", "CO2", "H2", "O2")
    }
    signals = _hypothesis_signals(axial, budget, fb)

    bed_payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "max_global_iter": int(args.max_global_iter),
        "wall_s": float(wall),
        "r5_scale": float(reactor.config.r5_scale),
        "converged": bool(result.get("converged")),
        "rms_E": result.get("rms_scaled_energy_final"),
        "exit_dry": exit_dry,
        "bed_exit_dry": bed_exit,
        "axial_phase": axial,
        "source_budget": {
            "bed_integrated_co_terms_mol_s": budget["bed_integrated_co_terms_mol_s"],
            "bed_integrated_co2_terms_mol_s": budget["bed_integrated_co2_terms_mol_s"],
            "segments": budget["segments"],
            "cells_compact": [
                {
                    "i": r["cell_index"],
                    "xi": r["xi_over_H"],
                    "T_K": r["T_K"],
                    "wet_y": r["wet_y"],
                    "co_terms": r["co_terms"],
                    "co2_terms": r["co2_terms"],
                    "raw_r5": {
                        "r5b_lim": r["raw"]["r5b_lim"],
                        "r5d_lim": r["raw"]["r5d_lim"],
                    },
                }
                for r in budget["cells"]
            ],
        },
        "hypothesis_signals": {k: signals[k] for k in ("H1", "H2", "H3")},
        "bild74_T_compare": signals["H3"]["T_compare_lower_bed"],
    }
    write_json(out_bed, bed_payload)
    print(f"wrote {out_bed}", flush=True)

    fb_payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "wall_s": float(wall),
        "freeboard": fb,
        "hypothesis_H4": signals["H4"],
        "capacity_check": {
            "bed_exit_CO2_dry": bed_exit["CO2"],
            "bild_fb_delta_CO2_wet": -0.015,
            "true_dry_CO2": TRUE_DRY["CO2"],
            "even_full_bild_fb_fix_floor": bed_exit["CO2"] - 0.015,
            "still_above_true_dry": bed_exit["CO2"] - 0.015 > TRUE_DRY["CO2"],
        },
    }
    write_json(out_fb, fb_payload)
    print(f"wrote {out_fb}", flush=True)

    variants: list[dict[str, Any]] = [
        {
            "label": "baseline",
            "wall_s": float(wall),
            "converged": bool(result.get("converged")),
            "rms_E": result.get("rms_scaled_energy_final"),
            "bed_exit_dry": bed_exit,
            "exit_dry": exit_dry,
            "score_vs_harness": co_co2_score(exit_dry["CO"], exit_dry["CO2"]),
            "score_vs_true_dry": abs(exit_dry["CO"] - TRUE_DRY["CO"]) / TRUE_DRY["CO"]
            + abs(exit_dry["CO2"] - TRUE_DRY["CO2"]) / TRUE_DRY["CO2"],
            "cfg": {},
        }
    ]
    if not args.skip_variants:
        for label, ov in (
            ("diag_r5_scale_0.10", {"r5_scale": 0.10}),
            ("diag_r5_scale_0.50", {"r5_scale": 0.50}),
            ("diag_dense_frac_0.25", {"gas_inlet_dense_frac": 0.25}),
            ("diag_dense_frac_0.75", {"gas_inlet_dense_frac": 0.75}),
        ):
            print(f"=== variant {label} ===", flush=True)
            variants.append(
                _run_variant(label, max_global_iter=int(args.max_global_iter), cfg_overrides=ov)
            )

    verdict = _verdict(signals, fb, variants, exit_dry)
    verdict_payload = {
        "case_key": "CASE_HTW_WESSELING_1",
        "basis_ref": "data/phase2_coco2_basis_ladder.json",
        "bed_decomp_ref": str(out_bed),
        "fb_link_ref": str(out_fb),
        "model_exit_dry": exit_dry,
        "model_bed_exit_dry": bed_exit,
        "true_dry_target": TRUE_DRY,
        "harness_REF": {"CO": REF_CO, "CO2": REF_CO2},
        "verdict": verdict,
    }
    write_json(out_verdict, verdict_payload)
    print(f"wrote {out_verdict}", flush=True)
    print(json.dumps(verdict["next_knife"], indent=2, ensure_ascii=False))
    print(json.dumps(verdict["success_criteria"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
