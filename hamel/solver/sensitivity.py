"""敏感性几何、悬停改送和已有中间量的诊断打包。ADR 0039。"""
from __future__ import annotations

import math

import numpy as np

from hamel.properties import (
    ash_specific_heat_j_per_kg_k, char_specific_enthalpy_j_per_kg, gas_constant_j_mol_k,
)
from hamel.solver.frozen import AssemblyError

# ADR 0014 B。锥段高度不随格子数缩放。
CONE_HEIGHT_M = 0.2875
# ADR 0014 E。只作对照，不回写流动。
KURKELA_BED_M_S = 0.74
KURKELA_FREEBOARD_M_S = 0.67
_AMBIENT_K = 298.15
_CLOSED_FORM = "闭式比对在报告工具中按 ADR 0040 第 5 项执行，不参与求解，不设通过门。"
_FREEBOARD_CLOSURE = "全自由空间闭合在报告工具中按 ADR 0040 第 5 项逐股记账，不加反应项，不进入残差。"


def diameter_at(height_m: float, bed_diameter_m: float, freeboard_diameter_m: float) -> float:
    """自床面向上的直径。锥段内线性，锥段以上为自由空间直径。ADR 0014 B。"""
    height = float(height_m)
    if height <= 0.0:
        return float(bed_diameter_m)
    if height >= CONE_HEIGHT_M:
        return float(freeboard_diameter_m)
    fraction = height / CONE_HEIGHT_M
    return float(bed_diameter_m) + (float(freeboard_diameter_m) - float(bed_diameter_m)) * fraction


def edge_diameters(count: int, cell_height_m: float, bed_diameter_m: float, freeboard_diameter_m: float) -> np.ndarray:
    """自由空间格边界上的直径。下标 0 是床面。"""
    return np.array([
        diameter_at(index * float(cell_height_m), bed_diameter_m, freeboard_diameter_m)
        for index in range(int(count) + 1)
    ])


def _pieces(z0: float, z1: float):
    if z1 <= CONE_HEIGHT_M or z0 >= CONE_HEIGHT_M:
        return ((z0, z1),)
    return ((z0, CONE_HEIGHT_M), (CONE_HEIGHT_M, z1))


def section_volume_m3(z0: float, z1: float, bed_diameter_m: float, freeboard_diameter_m: float) -> float:
    """圆台体积。整段落在锥段或直筒上时与圆柱或单段圆台是同一个式子。"""
    total = 0.0
    for start, end in _pieces(float(z0), float(z1)):
        height = end - start
        radius_1 = 0.5 * diameter_at(start, bed_diameter_m, freeboard_diameter_m)
        radius_2 = 0.5 * diameter_at(end, bed_diameter_m, freeboard_diameter_m)
        total += math.pi * height * (radius_1**2 + radius_1 * radius_2 + radius_2**2) / 3.0
    return total


def section_wall_m2(z0: float, z1: float, bed_diameter_m: float, freeboard_diameter_m: float) -> float:
    """圆台侧面积。直筒两端直径相同，退回 π·D·h。"""
    total = 0.0
    for start, end in _pieces(float(z0), float(z1)):
        height = end - start
        radius_1 = 0.5 * diameter_at(start, bed_diameter_m, freeboard_diameter_m)
        radius_2 = 0.5 * diameter_at(end, bed_diameter_m, freeboard_diameter_m)
        if radius_1 == radius_2:
            total += math.pi * (radius_1 + radius_2) * height
            continue
        slant = math.sqrt(height**2 + (radius_2 - radius_1) ** 2)
        total += math.pi * (radius_1 + radius_2) * slant
    return total


def interface_balance_index(height_m, bed_top_m, count, cell_height_m) -> int:
    """二次风是界面入口，必须落在某一格的下边界上。返回衡算下标，0 是最上层。

    不插值。床面以上的相对高度除以格高，偏离最近整数超过 1e-8 即失败。
    """
    cell_height = float(cell_height_m)
    if cell_height <= 0.0 or int(count) < 1:
        raise AssemblyError("自由空间格子高度和格数必须为正")
    quotient = (float(height_m) - float(bed_top_m)) / cell_height
    nearest = round(quotient)
    if abs(quotient - nearest) > 1.0e-8 or nearest < 0 or nearest >= int(count):
        raise AssemblyError(f"二次风高度 {height_m} m 没有落在自由空间格子的下边界上")
    return int(count) - 1 - int(nearest)


def solid_enthalpy_w(char_kg_s, ash_kg_s, fuel_kg_s, temperature_k, frozen) -> float:
    """源状态固体温度下的焦炭、灰和燃料焓。与残差里的比焓是同一组函数。"""
    temperature = float(temperature_k)
    fuel = frozen.fuel_formation_j_per_kg + frozen.fuel_heat_capacity_j_per_kg_k * (temperature - _AMBIENT_K)
    return (
        float(char_kg_s) * char_specific_enthalpy_j_per_kg(temperature)
        + float(ash_kg_s) * ash_specific_heat_j_per_kg_k() * (temperature - _AMBIENT_K)
        + float(fuel_kg_s) * fuel
    )


def _portion(flow, part, total) -> float:
    total = float(total)
    if total == 0.0:
        return 0.0
    return float(flow) * float(part) / total


def _item(kind, klass, char, ash, fuel, temperature, frozen, undo):
    if char == 0.0 and ash == 0.0 and fuel == 0.0:
        return None
    return {
        "kind": kind, "klass": int(klass), "char": float(char), "ash": float(ash), "fuel": float(fuel),
        "temperature": float(temperature), "undo": undo,
        "enthalpy": solid_enthalpy_w(char, ash, fuel, temperature, frozen),
    }


def redirect_ledger(
    frozen, state, migration_free, migration_hover, ash_migration_free, ash_migration_hover,
    bed_ash_share, ash_release_free, ash_release_hover,
) -> tuple:
    """进入悬停状态的四项流。A 返回空。

    模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038。
    顶格 a^t 已经是带出，底格 b^t 已经回到床顶，这两股不改送。
    """
    if frozen.hover_closure == "A":
        return ()
    if frozen.hover_closure not in {"B", "C"}:
        raise AssemblyError(f"不能识别的悬停闭合 {frozen.hover_closure}")
    n_free, n_class = frozen.n_freeboard, frozen.n_class
    items = []
    for klass in range(n_class):
        for cell in range(1, n_free):
            coeff = float(frozen.a_hover_per_s[cell, klass])
            total = coeff + float(frozen.a_free_per_s[cell, klass])
            mobile = _mobile_ash(state, frozen, cell)
            item = _item(
                "a_t", klass, coeff * float(state.freeboard_char_free_kg[cell, klass]),
                coeff * float(mobile[klass]),
                _portion(frozen.fuel_up_free_kg_s[cell, klass], coeff, total),
                state.freeboard_temperature_k[cell], frozen, ("free", cell - 1),
            )
            if item:
                items.append(item)
        for cell in range(n_free - 1):
            coeff = float(frozen.b_hover_per_s[cell, klass])
            total = coeff + float(frozen.b_free_per_s[cell, klass])
            mobile = _mobile_ash(state, frozen, cell)
            item = _item(
                "b_t", klass, coeff * float(state.freeboard_char_free_kg[cell, klass]),
                coeff * float(mobile[klass]),
                _portion(frozen.fuel_down_free_kg_s[cell, klass], coeff, total),
                state.freeboard_temperature_k[cell], frozen, ("free", cell + 1),
            )
            if item:
                items.append(item)
        share = float(frozen.ejection_to_hover[klass])
        thrown = float(frozen.k_auf_per_s[0]) * share
        item = _item(
            "ejection", klass, thrown * float(state.bed_char_kg[0, klass]),
            thrown * float(bed_ash_share[0, klass]),
            share * float(frozen.fuel_up_bed_kg_s[0, klass]),
            state.bed_temperature_k[0], frozen, ("free", n_free - 1),
        )
        if item:
            items.append(item)
        for cell in range(n_free):
            if klass + 1 >= n_class:
                break
            theta = float(frozen.theta[cell, klass])
            char_in = float(migration_free[cell, klass + 1] + migration_hover[cell, klass + 1])
            ash_in = float(ash_migration_free[cell, klass + 1] + ash_migration_hover[cell, klass + 1])
            item = _item(
                "theta", klass, theta * char_in, theta * ash_in, 0.0,
                state.freeboard_temperature_k[cell], frozen, ("free", cell),
            )
            if item:
                items.append(item)
    for cell in range(n_free):
        if not bool(frozen.free_absent[cell, 0]):
            continue
        released = float(np.sum(ash_release_free[cell]) + np.sum(ash_release_hover[cell]))
        item = _item("release", 0, 0.0, released, 0.0, state.freeboard_temperature_k[cell], frozen, ("free", cell))
        if item:
            items.append(item)
    return tuple(items)


def _mobile_ash(state, frozen, cell) -> np.ndarray:
    ash = np.array(state.freeboard_ash_free_kg[cell], dtype=float, copy=True)
    char = state.freeboard_char_free_kg[cell]
    for klass in range(1, frozen.n_class):
        if char[klass] <= 0.0:
            ash[klass] = 0.0
    return ash


def redirect_feedback(closure: str, redirect) -> tuple[float, float, float]:
    """本轮改送进入下一轮预计算的三份质量。

    返回床顶回落增量、出口灰、出口燃料。B 把焦炭、灰和燃料加进床层总流递推；
    C 把灰和燃料记为外部出口。顶格 a^t 与底格 b^t 不在改送量里，这里不再加一次。
    模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038。
    """
    char, ash, fuel = float(redirect[0]), float(redirect[1]), float(redirect[2])
    if closure == "B":
        return char + ash + fuel, 0.0, 0.0
    if closure == "C":
        return 0.0, ash, fuel
    return 0.0, 0.0, 0.0


def visible_bubble_m3_s(epsilon_b, u_b, bed_diameter_m) -> np.ndarray:
    """可见气泡流 ε_b·u_b·床面积。下标 0 是最上层。ADR 0017 第 3 项。"""
    area = math.pi * (0.5 * float(bed_diameter_m)) ** 2
    return np.asarray(epsilon_b, dtype=float) * np.asarray(u_b, dtype=float) * area


def bubble_balance_ratio(bubble_mol_s, temperature_k, pressure_pa, visible_m3_s) -> float | None:
    """衡算气泡相摩尔流除以可见气泡流换成的摩尔流。只做单位换算。"""
    volume = float(visible_m3_s)
    if volume <= 0.0:
        return None
    visible_mol = volume * float(pressure_pa) / (gas_constant_j_mol_k() * float(temperature_k))
    if visible_mol == 0.0:
        return None
    return float(bubble_mol_s) / visible_mol


def pack_reported(profile, flow, case, released) -> dict:
    """只收预计算里已经算好的量。闭式比对保持未计算。"""
    from hamel.properties.thermo import volume_flow_m3_s

    visible = visible_bubble_m3_s(profile.epsilon_b, profile.u_b, case.bed_diameter_m)
    ratio = np.asarray(profile.d_b, dtype=float) / float(case.bed_diameter_m)
    top_mol = float(np.sum(flow.free_mol_s[0]))
    top_area = math.pi * (0.5 * float(case.free_diameter_m)) ** 2
    free_u = volume_flow_m3_s(top_mol, float(flow.free_temperature[0]), float(flow.pressure)) / top_area
    n_free = int(case.n_free)
    return {
        "visible_bubble_m3_s": visible,
        "bubble_diameter_over_bed": ratio,
        "cells_db_exceeds_d": int(np.sum(np.asarray(profile.d_b, dtype=float) > float(case.bed_diameter_m))),
        "bed_superficial_m_s": np.asarray(flow.u0, dtype=float),
        "freeboard_superficial_m_s": float(free_u),
        "kurkela_bed_m_s": KURKELA_BED_M_S,
        "kurkela_freeboard_m_s": KURKELA_FREEBOARD_M_S,
        "fuel_release_mass_gap_kg_s": np.asarray(released.mass_gap_kg_s[:n_free], dtype=float),
        "fuel_release_element_gap_kg_s": {
            name: np.asarray(values[:n_free], dtype=float) for name, values in released.element_gap_kg_s.items()
        },
        "freeboard_mass_element_closure": None,
        "freeboard_mass_element_closure_note": _FREEBOARD_CLOSURE,
        "freeboard_closed_form": None,
        "freeboard_closed_form_note": _CLOSED_FORM,
    }
