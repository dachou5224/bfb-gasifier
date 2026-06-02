#!/usr/bin/env python3
"""Phase 1 审计：HTW Wesseling Table 2 LU（CASE_HTW_WESSELING_1）出口温度与干基主气相。

对照 ``data/validation_cases.json`` 中 ``outputs``，容差见 ``validation_numeric_tolerances()``。
干基对比**必须**使用 ``result[\"exit_gas_dry\"]``（由 ``exit_gas`` 湿基换算；``exit_gas`` 为顶格**两相摩尔流加权和**组成，见 ``Cell._mole_fractions(\"combined\")``）。

用法（在 ``bfb-gasifier/`` 目录下）::

    python3 scripts/audit_phase1_htw_lu.py
    python3 scripts/audit_phase1_htw_lu.py --json
    python3 scripts/audit_phase1_htw_lu.py --strict

环境变量 ``BFB_RELAX_VALIDATION=1``：仅打印对比，退出码 0（与 test_table2_LU 一致）。

Source: tests/validation_case_utils.py；Phase 1 计划
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path
from typing import Any

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from src.core.reactor import Reactor

from tests.validation_case_utils import (
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
    estimate_gas_feeds,
    json_numeric_or_none,
    load_case_LU,
    load_validation_case_node,
    strict_validation_gate,
    validation_numeric_tolerances,
)


def _load_char_ledger_summary_fn():
    mod_path = _REPO / "scripts" / "audit_char_mass_conservation_lu.py"
    spec = importlib.util.spec_from_file_location("audit_char_mass_conservation_lu", mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod.summarize_reactor_char_ledger


def _dry_main_four_sum(y_dry: dict[str, float]) -> float:
    return sum(float(y_dry.get(sp, 0.0) or 0.0) for sp in ("CO", "CO2", "H2", "CH4"))


def _overall_validation_pass(
    *,
    validation_candidate_ok: bool,
    pass_T: bool,
    pass_species: bool,
    carbon_ok: bool,
    relax: bool,
) -> bool:
    if relax:
        return True
    return bool(validation_candidate_ok and pass_T and pass_species and carbon_ok)


def _char_mass_ledger_reasons(
    ledger: dict[str, float],
    *,
    bed_transport_summary: dict[str, float] | None = None,
    max_rel_residual: float = 0.02,
) -> list[str]:
    fresh = abs(float(ledger.get("fresh_char_kg_s", 0.0) or 0.0))
    residual = abs(float(ledger.get("residual_char_kg_s", 0.0) or 0.0))
    migration = abs(float(ledger.get("size_migration_char_kg_s", 0.0) or 0.0))
    reasons: list[str] = []
    if fresh <= 1e-12:
        reasons.append("char_mass_missing_fresh_feed")
        return reasons
    if bed_transport_summary is not None:
        local_norm = float(bed_transport_summary.get("max_cell_char_residual_norm_local", 0.0) or 0.0)
        auf_gap = abs(float(bed_transport_summary.get("max_abs_auf_in_gap_vs_below_kg_s", 0.0) or 0.0))
        ab_gap = abs(float(bed_transport_summary.get("max_abs_ab_in_gap_vs_above_kg_s", 0.0) or 0.0))
        if local_norm > 0.10:
            reasons.append("bed_char_cell_residual>10.0%_local_transport")
        if auf_gap > 1e-9:
            reasons.append("bed_char_auf_in_gap_nonzero")
        if ab_gap > 1e-9:
            reasons.append("bed_char_ab_in_gap_nonzero")
        if migration / fresh > 1e-9:
            reasons.append("char_size_migration_not_conservative")
        return reasons
    if residual / fresh > float(max_rel_residual):
        reasons.append(f"char_mass_residual>{100.0 * float(max_rel_residual):.1f}%_fresh")
    if migration / fresh > 1e-9:
        reasons.append("char_size_migration_not_conservative")
    return reasons


def run_audit(
    *,
    strict: bool,
) -> tuple[dict[str, Any], bool, list[str]]:
    """返回 (结果字典, 是否通过 validation 容差, 警告列表)。"""
    warnings: list[str] = []
    case = load_case_LU()
    feeds = estimate_gas_feeds(case)
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)

    # 低级错误自检：进料与 ER 链一致
    if abs(cfg.O2_feed - feeds["O2_feed"]) > 1e-6:
        warnings.append(
            f"O2_feed 不一致: cfg={cfg.O2_feed} vs estimate={feeds['O2_feed']}（SI/进料链错误）"
        )

    relax = os.environ.get("BFB_RELAX_VALIDATION", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )

    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)
    char_ledger = _load_char_ledger_summary_fn()(
        reactor,
        result,
        mode="phase1",
        solve_kwargs=dict(PHASE1_HTW_LU_SOLVE_KWARGS),
    )
    validation_candidate_ok, validation_candidate_reasons = strict_validation_gate(result)
    char_ledger_summary = char_ledger["system_char_summary"]
    char_mass_reasons = _char_mass_ledger_reasons(
        char_ledger_summary,
        bed_transport_summary=char_ledger.get("bed_transport_summary"),
    )
    if char_mass_reasons:
        validation_candidate_ok = False
        validation_candidate_reasons = list(validation_candidate_reasons) + char_mass_reasons
    raw = load_validation_case_node(CASE_LU_VALIDATION_KEY)
    outs = raw["outputs"]
    tol = validation_numeric_tolerances()

    T_exit = float(result["T_profile"][-1])
    T_json = float(outs["exit_temperature_K"])
    T_meas = float(outs.get("exit_temperature_measured_K", T_json))
    err_T = min(
        abs(T_exit - T_json) / T_json,
        abs(T_exit - T_meas) / max(T_meas, 1.0),
    )

    y_dry = result["exit_gas_dry"]
    y_wet = result["exit_gas"]
    dry_ref = outs.get("exit_gas_dry_mol_frac", {})

    # 干基四组分归一（粗检；JSON 可能不含全部惰性）
    s4 = _dry_main_four_sum(y_dry)
    if s4 > 1.01 or (s4 < 0.05 and not relax):
        warnings.append(
            f"干基 CO+CO2+H2+CH4 合计={s4:.4f}（预期约在 0–1 之间；若异常请查湿基→干基换算）"
        )

    y_n2_wet = float(y_wet.get("N2", 0.0) or 0.0)
    if case.get("primary_agent") == "air_steam" and y_n2_wet < 1e-3 and not relax:
        warnings.append(
            f"空气工况出口湿基 y_N2={y_n2_wet:.4e} 过低（可能守恒/稀释路径错误，非调参可解释）"
        )
    if not validation_candidate_ok:
        warnings.append(
            "严格 validation 前置条件未满足: " + ", ".join(validation_candidate_reasons)
        )

    species_ok: dict[str, bool] = {}
    species_err: dict[str, float] = {}
    for sp in ("CO", "CO2", "H2", "CH4"):
        tgt = json_numeric_or_none(dry_ref.get(sp))
        sim = y_dry.get(sp)
        if tgt is None or tgt <= 0.0 or sim is None:
            continue
        rtol_sp = tol["rtol_CH4"] if sp == "CH4" else tol["rtol_CO_CO2_H2"]
        rel_sp = abs(float(sim) - tgt) / tgt
        species_err[sp] = rel_sp
        species_ok[sp] = rel_sp <= rtol_sp

    carbon_ok = True
    carbon_rel: float | None = None
    X_json = json_numeric_or_none(outs.get("carbon_conversion_pct"))
    sim_X_pct = float(result["carbon_conv"]) * 100.0
    if X_json is not None and X_json > 0.0:
        carbon_rel = abs(sim_X_pct - X_json) / X_json
        carbon_ok = carbon_rel <= tol["rtol_carbon_conv"]

    pass_T = err_T <= tol["rtol_T"]
    pass_species = len(species_ok) == 0 or all(species_ok.values())
    pass_all = _overall_validation_pass(
        validation_candidate_ok=validation_candidate_ok,
        pass_T=pass_T,
        pass_species=pass_species,
        carbon_ok=carbon_ok,
        relax=relax,
    )

    # --strict：自检 WARNING 时失败；与 BFB_RELAX_VALIDATION 并存时以 relax 为准（不强制非零退出）
    if strict and warnings and not relax:
        pass_all = False

    out: dict[str, Any] = {
        "case_key": CASE_LU_VALIDATION_KEY,
        "P_Pa": cfg.P,
        "T_inlet_K": cfg.T_inlet,
        "solver": PHASE1_HTW_LU_SOLVE_KWARGS.get("solver"),
        "solve_kwargs": dict(PHASE1_HTW_LU_SOLVE_KWARGS),
        "heat_loss_frac": cfg.heat_loss_frac,
        "u0_target": cfg.u0_target,
        "feed_check_O2_ok": abs(cfg.O2_feed - feeds["O2_feed"]) < 1e-6,
        "T_exit_K": T_exit,
        "err_T_vs_json_best": err_T,
        "pass_T": pass_T,
        "exit_gas_dry": {k: float(v) for k, v in y_dry.items() if isinstance(v, (int, float))},
        "species_relative_error": species_err,
        "species_pass": species_ok,
        "pass_species": pass_species,
        "carbon_conv_pct": sim_X_pct,
        "carbon_relative_error": carbon_rel,
        "pass_carbon": carbon_ok,
        "char_mass_ledger": char_ledger_summary,
        "bed_char_transport_summary": char_ledger.get("bed_transport_summary", {}),
        "top_char_residual_rows": char_ledger.get("top_char_residual_rows", []),
        "converged": result.get("converged"),
        "converged_outer": result.get("converged_outer"),
        "converged_fully": result.get("converged_fully"),
        "rms_scaled_final": result.get("rms_scaled_final"),
        "rms_scaled_gas_final": result.get("rms_scaled_gas_final"),
        "rms_scaled_solid_final": result.get("rms_scaled_solid_final"),
        "rms_scaled_energy_final": result.get("rms_scaled_energy_final"),
        "rms_scaled_component_max_final": result.get("rms_scaled_component_max_final"),
        "max_abs_scaled_final": result.get("max_abs_scaled_final"),
        "max_abs_scaled_solid_final": result.get("max_abs_scaled_solid_final"),
        "validation_candidate_ok": validation_candidate_ok,
        "validation_candidate_reasons": validation_candidate_reasons,
        "n_iter": result.get("n_iter"),
        "warnings": warnings,
        "relax_validation": relax,
    }

    return out, pass_all, warnings


def main() -> int:
    ap = argparse.ArgumentParser(description="Phase 1 HTW LU 出口温度与干基气相审计")
    ap.add_argument(
        "--json",
        action="store_true",
        help="仅输出 JSON（stdout）",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="存在自检 WARNING 时以非零退出码结束",
    )
    args = ap.parse_args()

    out, ok, warnings = run_audit(strict=args.strict)

    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
    else:
        print("Phase 1 审计 — CASE_HTW_WESSELING_1")
        print(f"  P = {out['P_Pa']:.0f} Pa  |  solver = {out['solver']}  |  heat_loss_frac = {out['heat_loss_frac']}")
        print(f"  feed O2 自检: {'OK' if out['feed_check_O2_ok'] else 'FAIL'}")
        print(f"  T_exit = {out['T_exit_K']:.2f} K  (err vs JSON best = {out['err_T_vs_json_best']:.4f})  pass_T = {out['pass_T']}")
        print(
            f"  validation_candidate_ok = {out['validation_candidate_ok']}  "
            f"converged={out['converged']} outer={out['converged_outer']} fully={out['converged_fully']}  "
            f"rms={out['rms_scaled_final']}"
        )
        print(
            "  residual groups: "
            f"gas={out.get('rms_scaled_gas_final')} "
            f"solid={out.get('rms_scaled_solid_final')} "
            f"energy={out.get('rms_scaled_energy_final')} "
            f"max_abs={out.get('max_abs_scaled_final')}"
        )
        for sp, e in out["species_relative_error"].items():
            p = out["species_pass"].get(sp, False)
            print(f"  干基 {sp}: rel_err={e:.4f}  pass={p}")
        if out["carbon_relative_error"] is not None:
            print(
                f"  碳转化率: {out['carbon_conv_pct']:.1f}%  rel_err={out['carbon_relative_error']:.4f}  pass={out['pass_carbon']}"
            )
        ledger = out.get("char_mass_ledger", {})
        if ledger:
            print(
                "  char ledger: "
                f"fresh={ledger.get('fresh_char_kg_s', 0.0):.4e} kg/s, "
                f"reaction_net={ledger.get('reaction_char_kg_s', 0.0):.4e} kg/s, "
                f"source={ledger.get('reaction_char_source_positive_kg_s', 0.0):.4e} kg/s, "
                f"sink={ledger.get('reaction_char_sink_negative_kg_s', 0.0):.4e} kg/s, "
                f"migration={ledger.get('size_migration_char_kg_s', 0.0):.4e} kg/s, "
                f"residual={ledger.get('residual_char_kg_s', 0.0):.4e} kg/s"
            )
        if warnings:
            print("  WARNING:")
            for w in warnings:
                print(f"    - {w}")
        print(f"  overall_pass = {ok}")

    if not ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
