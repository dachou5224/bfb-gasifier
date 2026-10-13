#!/usr/bin/env python3
"""热解 Gibbs（§4.3 An=b min-G）充分计算 + 数值稳定性审计。

覆盖：
1. LU coal + Hegermann Koks/Teer：T×P 全产物表
2. SLSQP 诊断：opt_success / element closure / G
3. 多初值扰动：解唯一性
4. 元素池缩放不变性
5. T 连续性（相邻点最大跃变）

用法：
  python scripts/audit_pyrolysis_gibbs_stability.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from src.core.cell_pyrolysis import (
    _allocate_pyrolysis_tar_and_minors,
    allocate_pyrolysis_products_hamel,
)
from src.core.elemental_ledger import vm_element_molar_rates
from src.core.species import fuel_tar_target_hc_ratio
from src.solvers.vorabrechnung.vorab_gibbs_seed import (
    _PYROLYSIS_GIBBS_CANDIDATES,
    _hamel_teer_gas_pool_strict_an_eq_b,
    _solve_pyrolysis_gibbs_an_eq_b,
)
from src.thermal.hegermann_yield import xc_koks_hegermann, xc_teer_hegermann

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "phase2_pyrolysis_gibbs_stability_audit.json"

# LU / VTT 常用干基
C_DRY = 61.5
H_DRY = 4.1
O_DRY = 21.8
ASH = 11.41
N_FRAC = 0.6
S_FRAC = 0.01
S_VOL = 0.5
M_VM = 1.0  # kg/s 归一

TEMPS_K = np.arange(700.0, 1400.0 + 1.0, 25.0)
PRESSURES_PA = np.array([1.01325e5, 1.0e6, 2.5e6, 4.0e6], dtype=float)
GAS_KEYS = list(_PYROLYSIS_GIBBS_CANDIDATES) + ["TAR1", "TAR2", "NH3", "H2S"]


def _lu_atoms() -> dict[str, float]:
    return vm_element_molar_rates(
        m_vm=M_VM,
        ash_dry_wt=ASH,
        C_dry=C_DRY,
        H_dry=H_DRY,
        O_dry=O_DRY,
        nitrogen_fraction=N_FRAC,
        sulfur_fraction=S_FRAC,
        sulfur_volatile_frac=S_VOL,
    )


def _full_product_row(T: float, P: float, atoms: dict[str, float]) -> dict:
    xk = float(xc_koks_hegermann(P, T))
    xt = float(xc_teer_hegermann(P, T))
    hc = fuel_tar_target_hc_ratio(C_DRY, H_DRY)
    minors, c_g, h_g, o_g = _allocate_pyrolysis_tar_and_minors(
        nC=atoms["C"],
        nH=atoms["H"],
        nO=atoms["O"],
        nN=atoms["N"],
        nS=atoms["S"],
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=xt,
        pyrolysis_koks_carbon_frac=xk,
        target_tar_hc_ratio=hc,
    )
    t0 = time.perf_counter()
    err = None
    diag = {}
    prod = {}
    try:
        elems, seed = _hamel_teer_gas_pool_strict_an_eq_b(
            c_gas=c_g, h_gas=h_g, o_gas=o_g
        )
        gas, diag = _solve_pyrolysis_gibbs_an_eq_b(
            T_K=T, P_Pa=P, elements=elems, mole_seed=seed
        )
        prod = allocate_pyrolysis_products_hamel(
            T_K=T,
            P_Pa=P,
            nC=atoms["C"],
            nH=atoms["H"],
            nO=atoms["O"],
            nN=atoms["N"],
            nS=atoms["S"],
            fuel_type="coal",
            pyrolysis_tar_carbon_frac=xt,
            pyrolysis_koks_carbon_frac=xk,
            target_tar_hc_ratio=hc,
        )
    except Exception as exc:  # noqa: BLE001 — 审计需吞并记录
        err = f"{type(exc).__name__}: {exc}"
        gas = {}
        elems = {"C": c_g, "H": h_g, "O": o_g}
    dt_ms = (time.perf_counter() - t0) * 1e3

    mol_gas = {k: float(prod.get(k, 0.0)) for k in GAS_KEYS}
    tot = sum(max(v, 0.0) for v in mol_gas.values()) or 1.0
    y = {k: mol_gas[k] / tot for k in GAS_KEYS}
    # 气相主组分 C 闭合（相对 pool）
    c_used = (
        mol_gas.get("CO2", 0.0)
        + mol_gas.get("CO", 0.0)
        + mol_gas.get("CH4", 0.0)
    )
    return {
        "T_K": T,
        "P_Pa": P,
        "xc_koks": xk,
        "xc_teer": xt,
        "pool_CHO": {"C": c_g, "H": h_g, "O": o_g},
        "mol_s": mol_gas,
        "y_mol": y,
        "gas_C_used": c_used,
        "gas_C_pool": c_g,
        "C_rel_err": abs(c_used - c_g) / max(c_g, 1e-30),
        "opt_success": bool(diag.get("opt_success", False)) if diag else False,
        "converged": bool(diag.get("converged", False)) if diag else False,
        "max_abs_element_closure": float(diag.get("max_abs_element_closure", np.nan))
        if diag
        else float("nan"),
        "G": float(diag.get("G", np.nan)) if diag else float("nan"),
        "message": str(diag.get("message", "")) if diag else "",
        "wall_ms": dt_ms,
        "error": err,
        "minors_TAR2": float(minors.get("TAR2", 0.0)),
    }


def _multistart_uniqueness(
    T: float, P: float, atoms: dict[str, float], *, n_pert: int = 12
) -> dict:
    """对同一 (T,P,pool) 用随机扰动初值，检查解是否落到同一点。"""
    xk = float(xc_koks_hegermann(P, T))
    xt = float(xc_teer_hegermann(P, T))
    hc = fuel_tar_target_hc_ratio(C_DRY, H_DRY)
    _, c_g, h_g, o_g = _allocate_pyrolysis_tar_and_minors(
        nC=atoms["C"],
        nH=atoms["H"],
        nO=atoms["O"],
        nN=atoms["N"],
        nS=atoms["S"],
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=xt,
        pyrolysis_koks_carbon_frac=xk,
        target_tar_hc_ratio=hc,
    )
    elems, seed = _hamel_teer_gas_pool_strict_an_eq_b(c_gas=c_g, h_gas=h_g, o_gas=o_g)
    base, base_diag = _solve_pyrolysis_gibbs_an_eq_b(
        T_K=T, P_Pa=P, elements=elems, mole_seed=seed
    )
    rng = np.random.default_rng(42)
    species = list(_PYROLYSIS_GIBBS_CANDIDATES)
    vecs = []
    Gs = []
    closures = []
    fails = 0
    for _ in range(n_pert):
        # 对数空间扰动 + 再投影到可行 seed 附近
        scale = rng.uniform(0.05, 5.0, size=len(species))
        pert = {sp: float(max(seed.get(sp, 0.0), 1e-14) * scale[i]) for i, sp in enumerate(species)}
        # 均匀混合 seed
        mix = rng.uniform(0.0, 1.0)
        for i, sp in enumerate(species):
            pert[sp] = (1.0 - mix) * pert[sp] + mix * float(max(seed.get(sp, 0.0), 1e-14))
        out, diag = _solve_pyrolysis_gibbs_an_eq_b(
            T_K=T, P_Pa=P, elements=elems, mole_seed=pert
        )
        if not diag.get("converged", False):
            fails += 1
            continue
        vecs.append(np.array([out[sp] for sp in species], dtype=float))
        Gs.append(float(diag["G"]))
        closures.append(float(diag["max_abs_element_closure"]))

    if not vecs:
        return {
            "T_K": T,
            "P_Pa": P,
            "n_ok": 0,
            "n_fail": fails,
            "max_l2_spread": float("nan"),
            "max_rel_G_spread": float("nan"),
            "max_closure": float("nan"),
            "base_converged": bool(base_diag.get("converged", False)),
        }

    stack = np.vstack(vecs)
    # 相对总量归一后比组成
    norms = np.sum(stack, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-30)
    y = stack / norms
    y_mean = np.mean(y, axis=0)
    l2 = np.sqrt(np.sum((y - y_mean) ** 2, axis=1))
    G_arr = np.asarray(Gs, dtype=float)
    G_ref = float(np.median(G_arr))
    rel_G = np.abs(G_arr - G_ref) / max(abs(G_ref), 1.0)
    base_y = np.array([base[sp] for sp in species], dtype=float)
    base_y = base_y / max(float(np.sum(base_y)), 1e-30)
    return {
        "T_K": T,
        "P_Pa": P,
        "n_ok": len(vecs),
        "n_fail": fails,
        "max_l2_spread": float(np.max(l2)),
        "mean_l2_spread": float(np.mean(l2)),
        "max_rel_G_spread": float(np.max(rel_G)),
        "max_closure": float(np.max(closures)),
        "base_converged": bool(base_diag.get("converged", False)),
        "base_y": {sp: float(base_y[i]) for i, sp in enumerate(species)},
        "mean_y": {sp: float(y_mean[i]) for i, sp in enumerate(species)},
    }


def _scale_invariance(T: float, P: float, atoms: dict[str, float]) -> dict:
    """元素池 ×s 后，摩尔产率应近似线性；组成 y 不变。"""
    xk = float(xc_koks_hegermann(P, T))
    xt = float(xc_teer_hegermann(P, T))
    hc = fuel_tar_target_hc_ratio(C_DRY, H_DRY)
    _, c_g, h_g, o_g = _allocate_pyrolysis_tar_and_minors(
        nC=atoms["C"],
        nH=atoms["H"],
        nO=atoms["O"],
        nN=atoms["N"],
        nS=atoms["S"],
        fuel_type="coal",
        pyrolysis_tar_carbon_frac=xt,
        pyrolysis_koks_carbon_frac=xk,
        target_tar_hc_ratio=hc,
    )
    species = list(_PYROLYSIS_GIBBS_CANDIDATES)
    ref_elems, ref_seed = _hamel_teer_gas_pool_strict_an_eq_b(
        c_gas=c_g, h_gas=h_g, o_gas=o_g
    )
    ref, ref_diag = _solve_pyrolysis_gibbs_an_eq_b(
        T_K=T, P_Pa=P, elements=ref_elems, mole_seed=ref_seed
    )
    ref_v = np.array([ref[sp] for sp in species], dtype=float)
    ref_tot = max(float(np.sum(ref_v)), 1e-30)
    ref_y = ref_v / ref_tot

    rows = []
    for s in (0.01, 0.1, 1.0, 10.0, 100.0):
        elems, seed = _hamel_teer_gas_pool_strict_an_eq_b(
            c_gas=s * c_g, h_gas=s * h_g, o_gas=s * o_g
        )
        out, diag = _solve_pyrolysis_gibbs_an_eq_b(
            T_K=T, P_Pa=P, elements=elems, mole_seed=seed
        )
        v = np.array([out[sp] for sp in species], dtype=float)
        tot = max(float(np.sum(v)), 1e-30)
        y = v / tot
        # 期望：n ≈ s * n_ref
        lin_err = float(np.max(np.abs(v - s * ref_v)) / max(s * ref_tot, 1e-30))
        y_l2 = float(np.sqrt(np.sum((y - ref_y) ** 2)))
        rows.append(
            {
                "scale": s,
                "converged": bool(diag.get("converged", False)),
                "closure": float(diag.get("max_abs_element_closure", np.nan)),
                "linearity_rel_max": lin_err,
                "y_l2_vs_ref": y_l2,
            }
        )
    return {
        "T_K": T,
        "P_Pa": P,
        "ref_converged": bool(ref_diag.get("converged", False)),
        "scales": rows,
    }


def _continuity(rows: list[dict], *, P: float) -> dict:
    """固定 P 下，对 T 扫描的组成跃变。"""
    sub = sorted([r for r in rows if abs(r["P_Pa"] - P) < 1.0 and r["error"] is None], key=lambda r: r["T_K"])
    if len(sub) < 2:
        return {"P_Pa": P, "n": len(sub), "max_dy_l2": float("nan")}
    species = list(_PYROLYSIS_GIBBS_CANDIDATES)
    jumps = []
    for a, b in zip(sub[:-1], sub[1:]):
        ya = np.array([a["y_mol"].get(sp, 0.0) for sp in species])
        yb = np.array([b["y_mol"].get(sp, 0.0) for sp in species])
        # 只看主气相（不含 tar）
        sa = max(float(np.sum(ya)), 1e-30)
        sb = max(float(np.sum(yb)), 1e-30)
        ya, yb = ya / sa, yb / sb
        l2 = float(np.sqrt(np.sum((yb - ya) ** 2)))
        jumps.append(
            {
                "T_lo": a["T_K"],
                "T_hi": b["T_K"],
                "dy_l2": l2,
                "dG": float(b["G"] - a["G"]) if np.isfinite(b["G"]) and np.isfinite(a["G"]) else float("nan"),
            }
        )
    worst = max(jumps, key=lambda j: j["dy_l2"])
    return {
        "P_Pa": P,
        "n": len(sub),
        "max_dy_l2": worst["dy_l2"],
        "worst_jump": worst,
        "mean_dy_l2": float(np.mean([j["dy_l2"] for j in jumps])),
    }


def main() -> None:
    atoms = _lu_atoms()
    print("=== LU atoms (m_vm=1 kg/s) ===")
    print({k: round(atoms[k], 4) for k in ("C", "H", "O", "N", "S")})

    grid_rows: list[dict] = []
    n_err = 0
    n_opt_fail = 0
    n_unconv = 0
    closures = []
    times = []

    print(f"\n=== T×P grid: {len(TEMPS_K)} × {len(PRESSURES_PA)} = {len(TEMPS_K)*len(PRESSURES_PA)} ===")
    for P in PRESSURES_PA:
        for T in TEMPS_K:
            row = _full_product_row(float(T), float(P), atoms)
            grid_rows.append(row)
            if row["error"]:
                n_err += 1
                print(f"  FAIL T={T:.0f} P={P:.3e}: {row['error']}")
            else:
                if not row["opt_success"]:
                    n_opt_fail += 1
                if not row["converged"]:
                    n_unconv += 1
                closures.append(row["max_abs_element_closure"])
                times.append(row["wall_ms"])

    # 多初值：代表性点
    multi_pts = [
        (800.0, 2.5e6),
        (1000.0, 2.5e6),
        (1173.15, 2.5e6),
        (1300.0, 2.5e6),
        (1173.15, 1.01325e5),
        (1173.15, 4.0e6),
        (700.0, 2.5e6),
        (1400.0, 2.5e6),
    ]
    print("\n=== multi-start uniqueness ===")
    multi_rows = []
    for T, P in multi_pts:
        m = _multistart_uniqueness(T, P, atoms, n_pert=16)
        multi_rows.append(m)
        print(
            f"  T={T:.0f} P={P:.2e}: ok={m['n_ok']}/{m['n_ok']+m['n_fail']} "
            f"max_l2={m['max_l2_spread']:.3e} max_dG_rel={m['max_rel_G_spread']:.3e} "
            f"closure={m['max_closure']:.3e}"
        )

    print("\n=== scale invariance @ LU 1173K 25bar ===")
    scale_row = _scale_invariance(1173.15, 2.5e6, atoms)
    for s in scale_row["scales"]:
        print(
            f"  s={s['scale']:>7g}: conv={s['converged']} "
            f"lin_err={s['linearity_rel_max']:.3e} y_l2={s['y_l2_vs_ref']:.3e} "
            f"closure={s['closure']:.3e}"
        )

    print("\n=== T-continuity (fixed P) ===")
    cont_rows = []
    for P in PRESSURES_PA:
        c = _continuity(grid_rows, P=float(P))
        cont_rows.append(c)
        w = c.get("worst_jump", {})
        print(
            f"  P={P:.2e}: max_dy_l2={c['max_dy_l2']:.4f} "
            f"@ {w.get('T_lo','?')}→{w.get('T_hi','?')}K  mean={c['mean_dy_l2']:.4f}"
        )

    # 代表性产物表（25 bar）
    slate_25 = [r for r in grid_rows if abs(r["P_Pa"] - 2.5e6) < 1.0 and r["error"] is None]
    slate_25.sort(key=lambda r: r["T_K"])

    summary = {
        "n_grid": len(grid_rows),
        "n_error": n_err,
        "n_opt_success_false": n_opt_fail,
        "n_unconverged_closure": n_unconv,
        "closure_max": float(np.max(closures)) if closures else float("nan"),
        "closure_p99": float(np.percentile(closures, 99)) if closures else float("nan"),
        "closure_median": float(np.median(closures)) if closures else float("nan"),
        "wall_ms_median": float(np.median(times)) if times else float("nan"),
        "wall_ms_p95": float(np.percentile(times, 95)) if times else float("nan"),
        "multistart_max_l2": float(np.nanmax([m["max_l2_spread"] for m in multi_rows])),
        "multistart_max_fail": int(max(m["n_fail"] for m in multi_rows)),
        "scale_max_y_l2": float(max(s["y_l2_vs_ref"] for s in scale_row["scales"])),
        "scale_max_lin_err": float(max(s["linearity_rel_max"] for s in scale_row["scales"])),
        "continuity_max_dy_l2": float(max(c["max_dy_l2"] for c in cont_rows)),
    }

    # 稳定性判据（报告用）
    verdict_parts = []
    ok = True
    if n_err > 0 or n_unconv > 0:
        ok = False
        verdict_parts.append(f"有失败/未闭合点: err={n_err} unconv={n_unconv}")
    if summary["closure_max"] > 1e-6:
        ok = False
        verdict_parts.append(f"closure_max={summary['closure_max']:.3e} > 1e-6")
    if summary["multistart_max_l2"] > 1e-3:
        ok = False
        verdict_parts.append(f"多初值组成散度 max_l2={summary['multistart_max_l2']:.3e}")
    if summary["scale_max_y_l2"] > 1e-3:
        ok = False
        verdict_parts.append(f"缩放不变性破坏 y_l2={summary['scale_max_y_l2']:.3e}")
    if summary["continuity_max_dy_l2"] > 0.15:
        # 警告而非硬失败：Hegermann f_tar(T) 本身可使产物缓变
        verdict_parts.append(
            f"警告: T 连续性最大跃变 dy_l2={summary['continuity_max_dy_l2']:.3f}（可能含 Hegermann 产量变化）"
        )
    if n_opt_fail > 0:
        verdict_parts.append(
            f"注意: SLSQP opt_success=False 共 {n_opt_fail} 点，但 closure 仍可能达标（现有回退逻辑）"
        )
    if ok and not verdict_parts:
        verdict_parts.append("网格全部闭合；多初值/缩放稳定")

    payload = {
        "meta": {
            "fuel": "LU coal",
            "C_H_O_ash_wt": [C_DRY, H_DRY, O_DRY, ASH],
            "m_vm_kg_s": M_VM,
            "candidates": list(_PYROLYSIS_GIBBS_CANDIDATES),
            "temps_K": TEMPS_K.tolist(),
            "pressures_Pa": PRESSURES_PA.tolist(),
            "path": "Hegermann xc_koks+xc_teer → strict An=b min-G (SLSQP)",
        },
        "summary": summary,
        "verdict_ok": ok,
        "verdict_notes": verdict_parts,
        "grid": grid_rows,
        "multistart": multi_rows,
        "scale_invariance": scale_row,
        "continuity": cont_rows,
        "slate_P25bar": [
            {
                "T_K": r["T_K"],
                "xc_koks": r["xc_koks"],
                "xc_teer": r["xc_teer"],
                "y_CO": r["y_mol"]["CO"],
                "y_CO2": r["y_mol"]["CO2"],
                "y_CH4": r["y_mol"]["CH4"],
                "y_H2": r["y_mol"]["H2"],
                "y_H2O": r["y_mol"]["H2O"],
                "y_TAR2": r["y_mol"]["TAR2"],
                "closure": r["max_abs_element_closure"],
                "opt_success": r["opt_success"],
                "G": r["G"],
            }
            for r in slate_25
        ],
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print(f"verdict_ok={ok}")
    for note in verdict_parts:
        print(f"  - {note}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
