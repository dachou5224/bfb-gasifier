"""HTW Wesseling Table 2 LU（``CASE_HTW_WESSELING_1``）验证数据加载与 Phase 1 基线配置。

单一事实来源供 ``test_table2_LU.py``、``scripts/audit_phase1_htw_lu.py`` 等复用，
避免 JSON 展平与 ``ReactorConfig`` 字段在多处漂移。

**Phase 1 基线**（与 ``test_table2_LU`` 中 slow 用例一致）：
- ``n_cells=10``, ``H_bed=6.0`` m（JSON ``bed_height_m``，非全炉 ``height_m``）, ``H_freeboard=0``（自由板未集成）
- ``heat_loss_frac=0.1``（``ReactorConfig`` 默认，中试炉热损失文献带 5–10% 见 validation_cases 注释）
- ``u0_target=None``（表观气速走内部默认/闭环）
- 求解（历史 GS 基线对照，需显式）：``ReactorConfig(..., allow_legacy_gs=True)`` + ``Reactor.solve(..., solver=\"gauss_seidel\", ...)``

**Phase 1 global NR 开发口径**：
- 与同一 LU shared config 复用几何/进料/动力学窗口
- 求解：``Reactor.solve(..., solver=\"global_nr\", max_global_iter=20, tol_global=1.0, nr_init_strategy=\"vorabrechnung\", nr_jacobian_strategy=\"block_tridiag_structured\")``
- 该口径当前用于 NR 审计与开发验证，尚不等同于 GS 基线稳定门控

**Phase 2 freeboard-aware 开发口径**：
- 在 shared ``global_nr`` 配置基础上恢复 ``H_freeboard = reactor_height - H_bed``
- 将床层顶部之后的 freeboard 作为独立 gas-only 轴向段求解
- Jacobian 走 ``nr_jacobian_strategy="auto"``，由显式图拓扑自动选择（side-block 路径避免误用纯 block-tridiag）
- 用于区分 ``bed exit`` 与 ``reactor exit`` 的验证/审计

Source: docs/CLAUDE.md；data/validation_cases.json；Phase 1 计划
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from src.core.feed_inlet import compute_gas_feeds_mol_s
from src.core.reactor import ReactorConfig

# 主数据：data/validation_cases.json → CASE_HTW_WESSELING_1（Table 7.1 Sim Nr.1 = 文献 Table 2 LU）
CASE_LU_VALIDATION_KEY = "CASE_HTW_WESSELING_1"

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _validation_cases_path() -> Path:
    return _REPO_ROOT / "data" / "validation_cases.json"


def _numeric_or_none(val: object) -> float | None:
    if isinstance(val, (int, float)):
        return float(val)
    return None


def nm3_h_to_mol_s_stp(flow_nm3_h: float) -> float:
    """Normvolumenstrom [Nm3/h] -> [mol/s]，按 STP 22.414 Nm3/kmol。"""
    return float(flow_nm3_h) * (1000.0 / 22.414) / 3600.0


def air_nm3_h_to_species_mol_s(flow_nm3_h: float) -> dict[str, float]:
    """空气 Normvolumenstrom 转成 `O2/N2` 摩尔流。"""
    n_tot = nm3_h_to_mol_s_stp(flow_nm3_h)
    return {
        "O2_mol_s": 0.21 * n_tot,
        "N2_mol_s": 0.79 * n_tot,
    }


def o2_nm3_h_to_species_mol_s(flow_nm3_h: float) -> dict[str, float]:
    """纯氧 Normvolumenstrom 转成 `O2/N2` 摩尔流。"""
    n_tot = nm3_h_to_mol_s_stp(flow_nm3_h)
    return {
        "O2_mol_s": n_tot,
        "N2_mol_s": 0.0,
    }


def _resolve_shared_ref(section: dict, root: dict | None, section_name: str | None = None) -> dict:
    if not isinstance(section, dict):
        return {}
    ref = section.get("_ref")
    if not ref:
        return dict(section)
    if root is None:
        root = load_validation_json_root()
    shared = root.get(str(ref), {})
    if not isinstance(shared, dict):
        return dict(section)
    if section_name and isinstance(shared.get(section_name), dict):
        shared = shared[section_name]
    merged = dict(shared)
    merged.update({k: v for k, v in section.items() if k != "_ref"})
    return merged


def extract_case_inlet_topology_hints(node: dict, root: dict | None = None) -> dict[str, object]:
    """从任意 validation case 节点中抽取 inlet / recycle / geometry 提示。

    设计目标：
    - 为 freeboard / secondary-air / recycle 开发提供统一 source-of-truth
    - 只做结构化读取与简单单位换算，不臆造缺失参数
    """
    ins = node.get("inputs", {})
    ga = ins.get("gasification_agent", {})
    react = _resolve_shared_ref(ins.get("reactor", {}), root, section_name="reactor")
    op = ins.get("operating_conditions", {})

    primary_air_nm3_h = _numeric_or_none(ga.get("primary_air_Nm3_h"))
    total_air_nm3_h = _numeric_or_none(ga.get("total_air_Nm3_h"))
    total_o2_nm3_h = _numeric_or_none(ga.get("total_O2_Nm3_h"))
    secondary_air_nm3_h = _numeric_or_none(ga.get("secondary_air_Nm3_h"))
    secondary_o2_nm3_h = _numeric_or_none(ga.get("secondary_O2_Nm3_h"))
    steam_feed_kg_h = _numeric_or_none(ga.get("steam_feed_kg_h"))
    secondary_xi = _numeric_or_none(ga.get("secondary_injection_xi_reactor"))
    secondary_agent = str(ga.get("secondary_injection_agent", "")).strip()

    if primary_air_nm3_h is None and total_air_nm3_h is not None and secondary_air_nm3_h is not None:
        primary_air_nm3_h = max(total_air_nm3_h - secondary_air_nm3_h, 0.0)

    primary_o2_nm3_h = _numeric_or_none(ga.get("primary_O2_Nm3_h"))
    if primary_o2_nm3_h is None and total_o2_nm3_h is not None and secondary_o2_nm3_h is not None:
        primary_o2_nm3_h = max(total_o2_nm3_h - secondary_o2_nm3_h, 0.0)

    primary_air_species = None if primary_air_nm3_h is None else air_nm3_h_to_species_mol_s(primary_air_nm3_h)
    primary_o2_species = None if primary_o2_nm3_h is None else o2_nm3_h_to_species_mol_s(primary_o2_nm3_h)
    secondary_air_species = None if secondary_air_nm3_h is None else air_nm3_h_to_species_mol_s(secondary_air_nm3_h)
    secondary_o2_species = None if secondary_o2_nm3_h is None else o2_nm3_h_to_species_mol_s(secondary_o2_nm3_h)

    secondary_agent_nm3_h = secondary_air_nm3_h if secondary_air_nm3_h is not None else secondary_o2_nm3_h
    secondary_agent_species = secondary_air_species if secondary_air_species is not None else secondary_o2_species

    bed_diameter = _numeric_or_none(react.get("bed_diameter_m"))
    freeboard_diameter = _numeric_or_none(react.get("freeboard_diameter_m"))
    legacy_diameter = _numeric_or_none(react.get("diameter_m"))
    if bed_diameter is None:
        bed_diameter = legacy_diameter
    if freeboard_diameter is None:
        freeboard_diameter = legacy_diameter

    bed_height = _numeric_or_none(react.get("bed_height_m"))
    reactor_height = _numeric_or_none(react.get("height_m")) or _numeric_or_none(react.get("total_height_m"))

    recirculation = bool(op.get("recirculation", False))
    recycle_topology = str(op.get("recycle_topology", "")).strip()
    recirculation_note = str(op.get("recirculation_note", "")).strip()
    cyclone_present = bool(op.get("cyclone_present", recirculation))
    if not recycle_topology and recirculation and "cyclone" in recirculation_note.lower():
        recycle_topology = "cyclone return"

    return {
        "plant": str(node.get("plant", "")),
        "secondary_injection_xi": secondary_xi,
        "secondary_injection_agent": secondary_agent,
        "secondary_injection_note": str(ga.get("secondary_injection_note", "")).strip(),
        "total_air_Nm3_h": total_air_nm3_h,
        "total_O2_Nm3_h": total_o2_nm3_h,
        "secondary_air_Nm3_h": secondary_air_nm3_h,
        "secondary_O2_Nm3_h": secondary_o2_nm3_h,
        "secondary_agent_Nm3_h": secondary_agent_nm3_h,
        "secondary_air_O2_mol_s": None if secondary_air_species is None else secondary_air_species["O2_mol_s"],
        "secondary_air_N2_mol_s": None if secondary_air_species is None else secondary_air_species["N2_mol_s"],
        "secondary_O2_O2_mol_s": None if secondary_o2_species is None else secondary_o2_species["O2_mol_s"],
        "secondary_O2_N2_mol_s": None if secondary_o2_species is None else secondary_o2_species["N2_mol_s"],
        "secondary_agent_O2_mol_s": None if secondary_agent_species is None else secondary_agent_species["O2_mol_s"],
        "secondary_agent_N2_mol_s": None if secondary_agent_species is None else secondary_agent_species["N2_mol_s"],
        "primary_air_Nm3_h": primary_air_nm3_h,
        "primary_O2_Nm3_h": primary_o2_nm3_h,
        "primary_air_O2_mol_s": None if primary_air_species is None else primary_air_species["O2_mol_s"],
        "primary_air_N2_mol_s": None if primary_air_species is None else primary_air_species["N2_mol_s"],
        "primary_O2_O2_mol_s": None if primary_o2_species is None else primary_o2_species["O2_mol_s"],
        "primary_O2_N2_mol_s": None if primary_o2_species is None else primary_o2_species["N2_mol_s"],
        "steam_feed_kg_h": steam_feed_kg_h,
        "recirculation": recirculation,
        "recirculation_note": recirculation_note,
        "recycle_topology": recycle_topology,
        "cyclone_present": cyclone_present,
        "reactor_height_m": reactor_height,
        "bed_height_m": bed_height,
        "bed_diameter_m": bed_diameter,
        "freeboard_diameter_m": freeboard_diameter,
        "secondary_air_port_mm_above_nozzle": _numeric_or_none(react.get("secondary_air_port_mm_above_nozzle")),
        "fuel_inlet_mm_above_nozzle": _numeric_or_none(react.get("fuel_inlet_mm_above_nozzle")),
    }


def apply_secondary_air_repartition(
    cfg: ReactorConfig,
    *,
    secondary_air_frac: float,
    injection_xi: float | None = None,
    secondary_T_K: float | None = None,
) -> ReactorConfig:
    """将总空气在 primary / secondary 之间分配，保持总 `ER` 不变。

    说明：
    - 只重分配 `O2/N2`，不改总蒸汽 `H2O_feed`
    - 适用于当前 `CASE_HTW_WESSELING_1` 这类：总空气由 `ER` 给定，
      但 `secondary air` 位置已知、split 未知的工况
    """
    frac = float(max(0.0, min(1.0, secondary_air_frac)))
    total_o2 = float(cfg.O2_feed + max(cfg.freeboard_secondary_O2_mol_s, 0.0))
    total_n2 = float(cfg.N2_feed + max(cfg.freeboard_secondary_N2_mol_s, 0.0))

    cfg.O2_feed = total_o2 * (1.0 - frac)
    cfg.N2_feed = total_n2 * (1.0 - frac)
    cfg.freeboard_secondary_O2_mol_s = total_o2 * frac
    cfg.freeboard_secondary_N2_mol_s = total_n2 * frac
    cfg.freeboard_secondary_H2O_mol_s = 0.0
    if injection_xi is not None:
        cfg.freeboard_secondary_injection_xi = float(injection_xi)
    if secondary_T_K is not None:
        cfg.freeboard_secondary_T_K = float(secondary_T_K)
    return cfg


def apply_case_secondary_stream(
    cfg: ReactorConfig,
    case: dict,
    *,
    secondary_T_K: float | None = None,
) -> ReactorConfig:
    """按 validation case 的 staged secondary stream 重分配总 oxidant。

    语义：
    - `ER` / `compute_gas_feeds_mol_s` 先给出 total oxidant budget
    - case 中若给了 explicit secondary stream，则从 primary feed 中扣出，
      并挂到 freeboard secondary inlet
    """
    sec_o2 = float(case.get("secondary_agent_O2_mol_s") or 0.0)
    sec_n2 = float(case.get("secondary_agent_N2_mol_s") or 0.0)
    total_o2 = float(cfg.O2_feed + max(cfg.freeboard_secondary_O2_mol_s, 0.0))
    total_n2 = float(cfg.N2_feed + max(cfg.freeboard_secondary_N2_mol_s, 0.0))

    cfg.freeboard_secondary_O2_mol_s = min(sec_o2, total_o2)
    cfg.freeboard_secondary_N2_mol_s = min(sec_n2, total_n2)
    cfg.freeboard_secondary_H2O_mol_s = 0.0
    cfg.O2_feed = max(total_o2 - cfg.freeboard_secondary_O2_mol_s, 0.0)
    cfg.N2_feed = max(total_n2 - cfg.freeboard_secondary_N2_mol_s, 0.0)

    if case.get("secondary_injection_xi") is not None:
        cfg.freeboard_secondary_injection_xi = float(case["secondary_injection_xi"])
    cfg.freeboard_secondary_T_K = float(case.get("T_inlet", cfg.T_inlet) if secondary_T_K is None else secondary_T_K)
    return cfg


def _flatten_htw_wesseling_1_case(node: dict) -> dict:
    """将 validation_cases.json 中嵌套结构展平为 Reactor 测试所用扁平字段。"""
    ins = node["inputs"]
    fuel = ins["fuel"]
    pa = fuel["proximate_analysis"]
    ua = fuel["ultimate_analysis_dry_wt_pct"]
    react = ins["reactor"]
    op = ins["operating_conditions"]
    ga = ins["gasification_agent"]
    hints = extract_case_inlet_topology_hints(node)
    agent_type = str(ga.get("type", ""))
    if "O2" in agent_type and "Steam" in agent_type and "Air" not in agent_type:
        primary_agent = "o2_steam"
    else:
        primary_agent = "air_steam"

    return {
        "_comment": node.get("_comment", ""),
        "_source_case_key": CASE_LU_VALIDATION_KEY,
        "reactor": "HTW_pressurised",
        "fuel": "brown_coal_RB",
        "P": float(op["pressure_MPa"]) * 1e6,
        "T_inlet": float(op["T_inlet_K"]),
        "fuel_feed": float(fuel["feed_rate_kg_h"]),
        "primary_agent": primary_agent,
        "primary_agent_note": agent_type,
        "ER": float(op["ER"]),
        "recirculation": bool(hints["recirculation"]),
        "recirculation_note": str(hints["recirculation_note"]),
        "recycle_topology": str(hints["recycle_topology"]),
        "cyclone_present": bool(hints["cyclone_present"]),
        # 轴向离散与流体力学：必须用**床层高度**（Table 7.1 / Section 7.1.2），勿用全炉总高 height_m
        "H_bed": float(react["bed_height_m"]),
        "D_freeboard": float(react["diameter_m"]),
        "bed_diameter_m": float(hints["bed_diameter_m"]),
        "freeboard_diameter_m": float(hints["freeboard_diameter_m"]),
        "bed_height_m": float(react["bed_height_m"]),
        "reactor_height_m": float(react["height_m"]),
        "primary_air_Nm3_h": hints["primary_air_Nm3_h"],
        "total_air_Nm3_h": hints["total_air_Nm3_h"],
        "secondary_air_Nm3_h": hints["secondary_air_Nm3_h"],
        "primary_O2_Nm3_h": hints["primary_O2_Nm3_h"],
        "total_O2_Nm3_h": hints["total_O2_Nm3_h"],
        "secondary_O2_Nm3_h": hints["secondary_O2_Nm3_h"],
        "secondary_injection_xi": hints["secondary_injection_xi"],
        "secondary_injection_agent": str(hints["secondary_injection_agent"]),
        "secondary_injection_note": str(hints["secondary_injection_note"]),
        "secondary_agent_Nm3_h": hints["secondary_agent_Nm3_h"],
        "secondary_agent_O2_mol_s": hints["secondary_agent_O2_mol_s"],
        "secondary_agent_N2_mol_s": hints["secondary_agent_N2_mol_s"],
        "moisture_wt": float(pa["moisture_wt_pct"]),
        "ash_dry_wt": float(pa["ash_dry_wt_pct"]),
        "VM_daf": float(pa["VM_daf_pct"]),
        "C_dry": float(ua["C"]),
        "H_dry": float(ua["H"]),
        "O_dry": float(ua["O"]),
        "N_dry": float(ua["N"]),
        "S_dry": float(ua["S"]),
    }


def load_case_LU() -> dict:
    data_path = _validation_cases_path()
    with open(data_path, encoding="utf-8") as f:
        data = json.load(f)
    if CASE_LU_VALIDATION_KEY not in data:
        raise KeyError(
            f"{data_path} 中缺少键 {CASE_LU_VALIDATION_KEY!r}（Table 2 LU / Table 7.1 Sim 1）"
        )
    return _flatten_htw_wesseling_1_case(data[CASE_LU_VALIDATION_KEY])


def load_validation_json_root() -> dict:
    """读取 validation_cases.json 根对象（含 validation_performance_summary）。"""
    data_path = _validation_cases_path()
    with open(data_path, encoding="utf-8") as f:
        return json.load(f)


def load_validation_case_node(key: str = CASE_LU_VALIDATION_KEY) -> dict:
    """读取 validation_cases.json 中完整 CASE 节点（含 inputs / outputs）。"""
    data = load_validation_json_root()
    if key not in data:
        raise KeyError(f"{_validation_cases_path()} 中缺少键 {key!r}")
    return data[key]


def json_numeric_or_none(val: object) -> float | None:
    if isinstance(val, (int, float)):
        return float(val)
    return None


def validation_numeric_tolerances() -> dict[str, float]:
    """与根节点 validation_performance_summary 对齐的相对容差（解析为标量）。"""
    return {
        "rtol_T": 0.15,
        "rtol_CO_CO2_H2": 0.15,
        "rtol_CH4": 0.20,
        "rtol_carbon_conv": 0.15,
    }


def strict_validation_gate(
    result: dict[str, Any],
    *,
    rms_max: float = 0.15,
) -> tuple[bool, list[str]]:
    """Return whether a solved state is trustworthy enough for strict validation.

    We intentionally separate:
    - numerical convergence quality
    - physical validation against literature targets

    so tests/audits can report "not converged enough for validation" instead of
    misclassifying that case as a model-vs-literature regression.
    """
    reasons: list[str] = []
    converged = bool(
        result.get(
            "converged_fully",
            result.get("converged_outer", result.get("converged", False)),
        )
    )
    if not converged:
        reasons.append("solver_not_converged")

    rms_raw = result.get("rms_scaled_final")
    if rms_raw is None:
        reasons.append("missing_rms_scaled_final")
    else:
        rms = float(rms_raw)
        if not np.isfinite(rms):
            reasons.append("nonfinite_rms_scaled_final")
        elif rms > float(rms_max):
            reasons.append(f"rms_scaled_final>{float(rms_max):.3f}")

    return len(reasons) == 0, reasons


def estimate_gas_feeds(case: dict) -> dict[str, float]:
    """与 ``ReactorConfig(ER=..., primary_agent=...)`` / ``compute_gas_feeds_mol_s`` 一致。"""
    fuel_kg_s = case["fuel_feed"] / 3600.0
    o2, h2o, n2 = compute_gas_feeds_mol_s(
        fuel_feed_kg_s=fuel_kg_s,
        moisture_wt=case["moisture_wt"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
    )
    return {"O2_feed": o2, "H2O_feed": h2o, "N2_feed": n2}


def build_phase1_htw_lu_reactor_config(case: dict | None = None) -> ReactorConfig:
    """构建 Phase 1（Table 2 LU）与 ``test_table2_LU`` slow 用例一致的 ``ReactorConfig``。

    单位：P [Pa]；fuel_feed [kg/s]；元素分析为 JSON 中 dry wt%；S_dry 与 ``ReactorConfig`` 一致。
    """
    if case is None:
        case = load_case_LU()
    return ReactorConfig(
        n_age_classes=1,
        n_cells=10,
        H_bed=case["H_bed"],
        H_freeboard=0.0,
        D_bed=case["bed_diameter_m"],
        P=case["P"],
        T_inlet=case["T_inlet"],
        fuel_type="coal",
        # 2026-04-14: updated from 0.5mm→1.0mm (char diameter after coal
        # devolatilization shrinkage ~0.65× of feed 1.5–3.0mm);
        # rho_s 1000→1200 (lignite char true density) and phi_s 0.86→0.75
        # (irregular char morphology) → u_mf ≈ 0.145 m/s at 2.5 MPa/1200 K,
        # matching Hamel (1999) Table 7.1 reference u_mf ≈ 0.15 m/s.
        rho_s=1200.0,
        d_p=1.0e-3,
        phi_s=0.75,
        eps_mf=0.45,
        fuel_feed=case["fuel_feed"] / 3600.0,
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
        moisture_wt=case["moisture_wt"],
        ash_dry_wt=case["ash_dry_wt"],
        VM_daf=case["VM_daf"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        nitrogen_fraction=case.get("N_dry", 0.0),
        # Hamel-aligned 主线改走 global-NR + Gibbs-minor。
        use_gibbs_minor=False,
        # 2026-04-05 local tuning after branch-acceptance guard:
        # dense≈0.30 with r7≈2.50 is now the best local LU window on both
        # raw and convergence-aware scores, while staying on the same stable
        # ~1140 K / ~0.90 carbon-conversion branch family as dense≈0.35.
        gas_inlet_dense_frac=0.30,
        recirculation_frac=0.1 if case.get("recirculation", True) else 0.0,
        recycle_gas=bool(case.get("recirculation", True)),
        enable_r12=False,
        top_solid_inlet_frac=0.50,
        r4_scale=0.50,
        r5_scale=0.75,
        r7_scale=2.50,
        # NR-only policy: legacy GS path is disabled at ``Reactor.solve`` entry.
        allow_legacy_gs=False,
    )


def build_phase1_htw_lu_global_nr_reactor_config(case: dict | None = None) -> ReactorConfig:
    """Phase 1 LU shared global NR 开发口径。

    说明
    ----
    与 GS 稳定基线分离，避免把 NR 开发阶段的底部 O2 相分配试验反写进
    `test_table2_LU` 依赖的 GS baseline。
    """
    cfg = build_phase1_htw_lu_reactor_config(case)
    # 2026-04-08 cell0 overlap sensitivity:
    # 在当前 NR 主线下，dense≈0.20 比 dense≈0.30/0.40 更能同时改善
    # CO / CO2 / carbon-conversion，而不诉诸直接修改 R6 kinetics。
    cfg.gas_inlet_dense_frac = 0.20
    # 2026-04-10 hydrodynamics source-of-truth locked to Hamel/Wein main chain:
    # - u_d: Wein (1992) Eq.3.12
    # - psi_b: Wein (1992) Eq.3.15
    # - d_b(h): Hilligardt ODE with Heinbockel pressurized extension
    # - lambda_b: Hamel Eq.3.42
    # - xi_b: Hamel Eq.3.34
    cfg.hydrodynamics_u_d_closure = "wein_1992_eq312"
    cfg.hydrodynamics_bubble_diameter_model = "hilligardt_ode"
    cfg.hydrodynamics_psi_b_strategy = "wein_1992"
    cfg.hydrodynamics_lambda_strategy = "hamel_280"
    cfg.hydrodynamics_xi_strategy = "hamel_regime"
    cfg.hydrodynamics_bubble_velocity_strategy = "heinbockel_eq343"
    cfg.hydrodynamics_bubble_ode_strategy = "heinbockel_eq341"
    cfg.use_gibbs_minor = True
    cfg.thesis_mode = True
    return cfg


def build_phase2_htw_lu_freeboard_reactor_config(
    case: dict | None = None,
    n_freeboard_cells: int = 8,
) -> ReactorConfig:
    """Phase 2 LU freeboard-aware 开发口径。

    在当前 shared global NR 配置上追加 freeboard 段，用于让验证口径从
    ``bed exit`` 过渡到更接近文献的 ``reactor exit``。
    """
    if case is None:
        case = load_case_LU()
    cfg = build_phase1_htw_lu_global_nr_reactor_config(case)
    cfg.H_freeboard = max(float(case["reactor_height_m"]) - float(case["bed_height_m"]), 0.0)
    cfg.n_freeboard_cells = int(max(n_freeboard_cells, 0))
    cfg.D_bed = float(case.get("bed_diameter_m", case["D_freeboard"]))
    cfg.freeboard_heat_loss_frac = 0.0
    cfg.freeboard_enabled_reactions = ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R10", "R11", "R12")
    cfg.freeboard_trajectory_model = "analytical_wirsum"
    cfg.freeboard_u_gb_scale = 1.0
    cfg.freeboard_beta_a_scale = 1.0
    cfg.freeboard_velocity_sigma = 0.60
    cfg.freeboard_velocity_bins = 5
    cfg.freeboard_cyclone_capture_char_frac = 0.90
    cfg.freeboard_cyclone_capture_ash_frac = 0.95
    cfg.freeboard_secondary_injection_xi = case.get("secondary_injection_xi")
    return apply_case_secondary_stream(cfg, case)


def build_phase1_htw_lu_refined_config(
    case: dict | None = None,
    n_fine: int = 3,
    dh_fine: float = 0.15,
) -> ReactorConfig:
    """Phase 1 LU 配置 + 底部轴向加密。

    Parameters
    ----------
    case:
        来自 ``load_case_LU()`` 的平铺字典，默认加载验证 JSON。
    n_fine:
        底部加密 cell 数量（默认 3，对应燃烧区 0–0.45 m）。
    dh_fine:
        每个底部加密 cell 的高度 [m]（默认 0.15 m）。

    Design
    ------
    总床高 H_bed 保持不变（与基线一致，用于流体力学计算）。
    底部 n_fine 个 cell 取高度 dh_fine；
    剩余 (H_bed - n_fine*dh_fine) 高度均匀分配给上方若干 cell，
    以尽量接近基线每 cell 0.6 m。
    总 cell 数 = n_fine + n_upper，其中
    n_upper = round((H_bed - n_fine * dh_fine) / 0.6)（向上取整，至少 1）。
    """
    if case is None:
        case = load_case_LU()

    H_bed = float(case["H_bed"])
    H_remaining = H_bed - n_fine * dh_fine
    if H_remaining <= 0:
        raise ValueError(
            f"n_fine={n_fine} × dh_fine={dh_fine} = {n_fine*dh_fine:.2f} m "
            f"already exceeds H_bed={H_bed} m"
        )
    # Number of upper cells: target ~0.6 m each to match baseline resolution above the combustion zone.
    import math
    n_upper = max(1, math.ceil(H_remaining / 0.6))
    dh_upper = H_remaining / n_upper

    n_cells = n_fine + n_upper

    return ReactorConfig(
        n_age_classes=1,
        n_cells=n_cells,
        H_bed=H_bed,
        H_freeboard=0.0,
        D_bed=case["bed_diameter_m"],
        P=case["P"],
        T_inlet=case["T_inlet"],
        fuel_type="coal",
        rho_s=1200.0,
        d_p=1.0e-3,
        phi_s=0.75,
        eps_mf=0.45,
        fuel_feed=case["fuel_feed"] / 3600.0,
        ER=case["ER"],
        primary_agent=case.get("primary_agent", "air_steam"),
        S_dry=case.get("S_dry", 0.0),
        moisture_wt=case["moisture_wt"],
        ash_dry_wt=case["ash_dry_wt"],
        VM_daf=case["VM_daf"],
        C_dry=case["C_dry"],
        H_dry=case["H_dry"],
        O_dry=case["O_dry"],
        nitrogen_fraction=case.get("N_dry", 0.0),
        gas_inlet_dense_frac=0.40,
        recirculation_frac=0.1 if case.get("recirculation", True) else 0.0,
        recycle_gas=bool(case.get("recirculation", True)),
        enable_r12=False,
        top_solid_inlet_frac=0.50,
        r4_scale=0.75,
        r5_scale=0.50,
        allow_legacy_gs=True,
    )


# Phase 1 shared solve kwargs（NR-only，与 Hamel 主路径一致）
PHASE1_HTW_LU_SOLVE_KWARGS: dict[str, Any] = {
    "max_global_iter": 20,
    "tol_global": 1.0,
    "solver": "global_nr",
    "nr_init_strategy": "vorabrechnung",
    "nr_jacobian_strategy": "block_tridiag_structured",
}


# Phase 1 global NR（与 Hamel 论文一致的主路径；``build_phase1_htw_lu_global_nr_reactor_config``）
PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS: dict[str, Any] = dict(PHASE1_HTW_LU_SOLVE_KWARGS)


# Phase 2 freeboard-aware global NR（显式图含 side-elements，Jacobian 走默认拓扑选择）
PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS: dict[str, Any] = dict(PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)
PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS["nr_jacobian_strategy"] = None
