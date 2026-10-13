"""Hamel (1999) Eq. (2.7)/(2.8) 能量衡算默认与模式常量。

Source: docs/hamel_submodels/02_cell_balances_and_exchange.md
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol

from src.core.species import liquid_water_formation_enthalpy_j_per_kg

# 固体组分列顺序与 ``src/core/cell.py`` 一致：[char, vm, moisture, ash]
SOLID_ENTHALPY_MODE_HAMEL_EQ28 = "hamel_eq28_component"  # 生产默认：端点瞬时 cp，非积分
SOLID_ENTHALPY_MODE_LEGACY_CP_MIX = "legacy_cp_mix_030304"
SOLID_ENTHALPY_MODE_ENDPOINT = "legacy_endpoint_cp"  # 与 hamel_eq28 同一公式的具名
SOLID_ENTHALPY_MODE_INTEGRAL = "thermo_consistent_cp_integral"  # 约化 opt-in，本轮不切默认

HEAT_LOSS_MODE_INLET_FRACTION = "inlet_fraction"
HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT = "hamel_eq27_explicit"

# 水分 sensible heat [J/(kg·K)]（液态水，VM/干燥源项另计潜热）
CP_MOISTURE_J_PER_KG_K: float = 4180.0

# 纯碳 char（石墨）标准生成焓 [J/kg]。反应式按 C 计量时必须为 0。
# 旧值 -15.866e6 是液态水生成焓，曾被 2026-07-06 审计误放到 char 列。
H_FORMATION_GRAPHITE_CHAR_J_PER_KG: float = 0.0
H_FORMATION_CHAR_J_PER_KG: float = H_FORMATION_GRAPHITE_CHAR_J_PER_KG

# 液态水标准生成焓 [J/kg]，与 species.HF_H2O_LIQUID_J_PER_MOL 同一基准。
H_FORMATION_LIQUID_WATER_J_PER_KG: float = float(liquid_water_formation_enthalpy_j_per_kg())

# ReactorConfig / Cell 默认（与 Hamel 主链一致，2026-07-07）
DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE = HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT
DEFAULT_SOLID_ENTHALPY_MODE = SOLID_ENTHALPY_MODE_HAMEL_EQ28

# 工程/审计对照：fraction 热损（H_in≤0 时 Q≈0，见 resolve_inlet_fraction_calibration_W）
LEGACY_ENERGY_BALANCE_HEAT_LOSS_MODE = HEAT_LOSS_MODE_INLET_FRACTION


def solid_component_hf_j_per_kg(
    *,
    h_f_dry: float,
    ash_dry_wt: float,
    VM_daf: float,
) -> tuple[float, float, float, float]:
    """固相四列标准生成焓 [J/kg]：char, VM, moisture, ash。

    纯碳 char：石墨基准 0。液态水分：HHV 液态水。
    VM：由 1 kg 干料 HHV 恒等式反算，``h_f,VM = h_f,dry / w_VM``。
    灰：0。

    Source: Hamel (1999) Eq. (2.8)；docs/audit_reports/secondary_diagnosis_2026-09-14.md §2
    """
    w_ash = float(ash_dry_wt) / 100.0
    w_vm = (float(VM_daf) / 100.0) * (1.0 - w_ash)
    return (
        float(H_FORMATION_GRAPHITE_CHAR_J_PER_KG),
        float(h_f_dry) / max(w_vm, 1e-9),
        float(H_FORMATION_LIQUID_WATER_J_PER_KG),
        0.0,
    )


class _EnergyBalanceCellLike:
    """Cell 能量衡算字段协议（供 apply 辅助函数）。"""

    heat_loss_frac: float
    heat_loss_mode: str
    solid_enthalpy_mode: str
    Q_wall_W: float
    Q_heat_exchanger_W: float


def apply_energy_balance_defaults_to_cell(cell: _EnergyBalanceCellLike) -> None:
    """初始化 cell 能量衡算模式为 Hamel 默认。"""
    cell.heat_loss_mode = DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE
    cell.solid_enthalpy_mode = DEFAULT_SOLID_ENTHALPY_MODE
    cell.Q_wall_W = 0.0
    cell.Q_heat_exchanger_W = 0.0


def compound_segment_heat_loss_fraction(total_frac: float, n_segments: int) -> float:
    """将总热损比例按段复合分配（与 bed/freeboard ``heat_loss_frac`` 轴向分配一致）。"""
    total = float(min(max(total_frac, 0.0), 1.0))
    n = int(max(n_segments, 1))
    if total <= 0.0:
        return 0.0
    return float(1.0 - math.pow(max(1.0 - total, 0.0), 1.0 / n))


def apply_freeboard_substep_gas_enthalpy_heat_loss(
    H_before: float,
    *,
    heat_loss_mode: str,
    heat_loss_frac_total: float,
    n_cells: int,
    n_substeps: int = 12,
    n_slices: int = 1,
    Q_wall_W: float = 0.0,
    Q_heat_exchanger_W: float = 0.0,
) -> float:
    """自由板气相焓步进热损（Eq. 2.7 轨迹侧）。

    ``inlet_fraction``：按段/片/子步复合 ``heat_loss_frac_total``，仅扣除 ``max(H, 0)`` 部分。
    ``hamel_eq27_explicit``：每子步扣除 ``(Q_W + Q_WÜ) / (n_cells · n_slices · n_substeps)`` [W]。
    """
    mode = str(heat_loss_mode or HEAT_LOSS_MODE_INLET_FRACTION).strip().lower()
    if mode == HEAT_LOSS_MODE_INLET_FRACTION:
        seg_loss = compound_segment_heat_loss_fraction(heat_loss_frac_total, n_cells)
        seg_loss_slice = compound_segment_heat_loss_fraction(seg_loss, n_slices)
        substep_loss = compound_segment_heat_loss_fraction(seg_loss_slice, n_substeps)
        h_pos = max(float(H_before), 0.0)
        return float(H_before) - substep_loss * h_pos
    if mode == HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT:
        q_total = max(float(Q_wall_W), 0.0) + max(float(Q_heat_exchanger_W), 0.0)
        denom = max(int(n_cells), 1) * max(int(n_slices), 1) * max(int(n_substeps), 1)
        return float(H_before) - q_total / denom
    raise ValueError(
        f"Unsupported heat_loss_mode={heat_loss_mode!r}; "
        f"expected {HEAT_LOSS_MODE_INLET_FRACTION!r} or {HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT!r}"
    )


def distribute_axial_explicit_heat_loss_watts(
    *,
    H_bed: float,
    H_freeboard: float,
    n_bed_cells: int,
    Q_wall_W_total: float,
    Q_heat_exchanger_W_total: float,
) -> tuple[list[float], list[float], float, float]:
    """按床高/自由板高度比例分配 reactor 级 Q_W / Q_WÜ [W]。"""
    bed_height = float(max(H_bed, 0.0))
    freeboard_height = float(max(H_freeboard, 0.0))
    total_height = bed_height + freeboard_height
    n_bed = int(max(n_bed_cells, 0))
    q_wall_total = float(max(Q_wall_W_total, 0.0))
    q_hex_total = float(max(Q_heat_exchanger_W_total, 0.0))

    if total_height <= 0.0 or n_bed <= 0:
        return [0.0] * n_bed, [0.0] * n_bed, 0.0, 0.0

    bed_seg_weight = (bed_height / total_height) / n_bed
    fb_weight = freeboard_height / total_height
    bed_q_wall = [q_wall_total * bed_seg_weight for _ in range(n_bed)]
    bed_q_hex = [q_hex_total * bed_seg_weight for _ in range(n_bed)]
    return bed_q_wall, bed_q_hex, q_wall_total * fb_weight, q_hex_total * fb_weight


@dataclass(frozen=True)
class FreeboardHeatLossBundle:
    """自由板轨迹求解热损配置（与 ReactorConfig 对齐）。"""

    heat_loss_mode: str
    heat_loss_frac: float
    Q_wall_W: float
    Q_heat_exchanger_W: float


def cell_inlet_fraction_heat_loss_W(cell, *, use_abs_H_in: bool = False) -> float:
    """单 bed cell 在 ``inlet_fraction`` 模式下 implied Q_loss [W]。"""
    from src.core.cell_balances import calc_energy_balance_H_in

    H_in = float(
        calc_energy_balance_H_in(
            N_b_in=cell.N_b_in,
            N_d_in=cell.N_d_in,
            T_in_gas=cell.T_in_gas,
            N_zu_b=cell.N_zu_b,
            N_zu_d=cell.N_zu_d,
            T_zu_gas=cell.T_zu_gas,
            N_rez_b=cell.N_rez_b,
            N_rez_d=cell.N_rez_d,
            T_rez_gas=cell.T_rez_gas,
            m_solid_rez=cell.m_solid_rez,
            T_rez_solid=cell.T_rez_solid,
            m_solid_in=cell.m_solid_in,
            T_in_solid=cell.T_in_solid,
            m_solid_zu=cell.m_solid_zu,
            T_zu_solid=cell.T_zu_solid,
            h_cache=getattr(cell, "_h_cache", None),
            m_solid_auf_in=cell.m_solid_auf_in,
            T_solid_auf_in=cell.T_solid_auf_in,
            m_solid_ab_in=cell.m_solid_ab_in,
            T_solid_ab_in=cell.T_solid_ab_in,
            ash_dry_wt=cell.solid.ash_dry_wt,
            VM_daf=cell.solid.VM_daf,
            h_f_dry=cell.solid.h_f_dry,
            solid_enthalpy_mode=getattr(cell, "solid_enthalpy_mode", DEFAULT_SOLID_ENTHALPY_MODE),
        )
    )
    basis = abs(H_in) if bool(use_abs_H_in) else max(H_in, 0.0)
    return float(cell.heat_loss_frac) * basis


def summarize_reactor_inlet_fraction_heat_loss_W(
    reactor,
    *,
    use_abs_H_in: bool = False,
) -> dict[str, float]:
    """汇总 reactor 在 ``inlet_fraction`` 下的 implied 热损功率 [W]。"""
    from src.core.cell_balances import calc_gas_enthalpy_flow
    from src.core.reactor import resolve_freeboard_heat_loss_bundle

    bed_loss = sum(
        cell_inlet_fraction_heat_loss_W(cell, use_abs_H_in=use_abs_H_in) for cell in reactor.cells
    )
    bundle = resolve_freeboard_heat_loss_bundle(reactor.config)
    top = reactor.cells[-1]
    fb_inlet_H = float(
        calc_gas_enthalpy_flow(
            top.N_b + top.N_d,
            float(top.T),
            getattr(top, "_h_cache", None),
        )
    )
    fb_basis = abs(fb_inlet_H) if bool(use_abs_H_in) else max(fb_inlet_H, 0.0)
    freeboard_loss = float(bundle.heat_loss_frac) * fb_basis
    return {
        "bed_W": float(bed_loss),
        "freeboard_W": float(freeboard_loss),
        "total_W": float(bed_loss + freeboard_loss),
    }


def resolve_reactor_Q_wall_W_total(cfg) -> float:
    """解析 reactor 级 Q_W [W]（Eq. 2.7 显式热损）。

    若 ``Q_wall_W_total`` 已设 (>0) 则直接使用；否则在 ``hamel_eq27_explicit`` 下用
    ``heat_loss_frac × fuel_feed × HHV_dry`` 近似标定（MJ/kg → W）。出口 T 正式 validation
    待反应/能量守恒校正后再做 Hamel 对标。
    """
    explicit = float(max(getattr(cfg, "Q_wall_W_total", 0.0), 0.0))
    if explicit > 0.0:
        return explicit
    mode = str(getattr(cfg, "energy_balance_heat_loss_mode", DEFAULT_ENERGY_BALANCE_HEAT_LOSS_MODE)).strip().lower()
    if mode != HEAT_LOSS_MODE_HAMEL_EQ27_EXPLICIT:
        return 0.0
    frac = float(max(getattr(cfg, "heat_loss_frac", 0.0), 0.0))
    fuel_kg_s = float(max(getattr(cfg, "fuel_feed", 0.0), 0.0))
    hhv_mj_kg = float(max(getattr(cfg, "HHV_dry", 0.0), 0.0))
    return frac * fuel_kg_s * hhv_mj_kg * 1e6


def resolve_inlet_fraction_calibration_W(reactor) -> tuple[dict[str, float], str]:
    """为 ``hamel_eq27_explicit`` 校准选择 implied Q 口径。

    生产默认 ``f·max(H_in,0)``；当 Phase2 全焓流为负时回退 ``f·|H_in|`` 审计桥接。
    """
    positive = summarize_reactor_inlet_fraction_heat_loss_W(reactor, use_abs_H_in=False)
    if float(positive["total_W"]) > 1.0:
        return positive, "positive_H_in"
    abs_bridge = summarize_reactor_inlet_fraction_heat_loss_W(reactor, use_abs_H_in=True)
    return abs_bridge, "abs_H_in_audit_bridge"

