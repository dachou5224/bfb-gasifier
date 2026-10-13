"""ReactorConfig 二次流/secondary oxidant 重分配（Phase 2 freeboard 开发口径）。"""

from __future__ import annotations

from src.core.reactor import ReactorConfig


def apply_secondary_air_repartition(
    cfg: ReactorConfig,
    *,
    secondary_air_frac: float,
    injection_xi: float | None = None,
    secondary_T_K: float | None = None,
) -> ReactorConfig:
    """将总空气在 primary / secondary 之间分配，保持总 ``ER`` 不变。"""
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
    """按 validation case 的 staged secondary stream 重分配总 oxidant。"""
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
