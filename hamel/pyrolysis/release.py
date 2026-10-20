"""未干燥与已干燥未热解两个库存的稳态释放份额。"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.properties.data import atomic_weight, atoms, molar_mass_kg_mol
from hamel.pyrolysis.particle import PyrolysisError
from hamel.pyrolysis.yields import LigniteYields, PYROLYSIS_GAS_SPECIES, pyrolysis_gas

# ADR 0012 第 3 项：|进 − 出| / 进 ≤ 0,1 %。元素进料为 0 时分母改用燃料质量。
_CLOSURE = 1.0e-3
_ELEMENTS = ("C", "H", "O", "N", "S")


@dataclass(frozen=True)
class ReleaseDistribution:
    """流率按 cell。库存是当量质量；dried_exit_kg_s 是 (1−水分) 乘当量出流。"""

    water_kg_s: np.ndarray
    gas_mol_s: dict[str, np.ndarray]
    benzene_mol_s: np.ndarray
    naphthalene_mol_s: np.ndarray
    char_kg_s: np.ndarray
    ash_kg_s: np.ndarray
    undried_exit_kg_s: np.ndarray
    dried_exit_kg_s: np.ndarray
    undried_kg: np.ndarray
    dried_kg: np.ndarray
    mass_gap_kg_s: np.ndarray
    element_gap_kg_s: dict[str, np.ndarray]


def feed_cell_index(height_m: float, bottom_m, top_m) -> int:
    """下标 0 是最上层。高度落在交界面上时归入上方的 cell。Source: Hamel Eq.2.6, S. 17"""
    bottom, top = _array("cell_bottom_m", bottom_m, 1), _array("cell_top_m", top_m, 1)
    if top.shape != bottom.shape or bottom.size < 1 or np.any(top <= bottom):
        raise PyrolysisError("每个 cell 的上边界必须高于下边界，且与下边界等长")
    height = float(height_m)
    matches = [index for index in range(bottom.size) if bottom[index] <= height <= top[index]]
    if not matches:
        raise PyrolysisError("fuel_feed_height 不在任何一个 cell 内")
    return max(matches, key=lambda index: (top[index], -index))


def release_distribution(
    feed_kg_s, fuel_feed_height_m: float, cell_bottom_m, cell_top_m, k_up_per_s,
    k_down_per_s, k_exit_per_s, tau_dry_s, tau_pyr_s, moisture_mass_fraction: float,
    ash_mass_fraction_wf: float, yields: LigniteYields, temperature_k: float, pressure_pa: float,
) -> ReleaseDistribution:
    """稳态解 W、P。τ = 0 时该库存不设状态，入流即转化流，不做除法。

    输运系数由调用方传入。解释性选择：ADR 0015 第 10 项。
    Source: Hamel Eq.4.23, S. 64
    Source: Hamel Eq.4.26, S. 69
    """
    feed = _array("feed_kg_s", feed_kg_s, 1)
    upward = _array("k_up_per_s", k_up_per_s, 2)
    cells, classes = upward.shape
    downward = _same("k_down_per_s", k_down_per_s, upward.shape)
    leaving = _same("k_exit_per_s", k_exit_per_s, upward.shape)
    dry, pyro = _same("tau_dry_s", tau_dry_s, (classes,)), _same("tau_pyr_s", tau_pyr_s, (classes,))
    if feed.shape != (classes,) or np.any(feed < 0.0) or np.any(upward < 0.0) or np.any(downward < 0.0) or np.any(leaving < 0.0):
        raise PyrolysisError("进料必须是非负粒级向量，输运系数不能为负")
    if np.any(dry < 0.0) or np.any(pyro < 0.0):
        raise PyrolysisError("时间常数不能为负")
    moisture, ash = float(moisture_mass_fraction), float(ash_mass_fraction_wf)
    if not (0.0 <= moisture < 1.0) or not (0.0 <= ash < 1.0):
        raise PyrolysisError("水分或灰分不在定义域内")
    if np.shape(cell_bottom_m) != (cells,):
        raise PyrolysisError("cell 数必须与输运系数的第 0 轴一致")
    cell = feed_cell_index(fuel_feed_height_m, cell_bottom_m, cell_top_m)
    equilibrium = pyrolysis_gas(yields, temperature_k, pressure_pa)
    undried, dried_stock = np.zeros((cells, classes)), np.zeros((cells, classes))
    dry_rate, pyro_rate = np.zeros((cells, classes)), np.zeros((cells, classes))
    for klass in range(classes):
        source = np.zeros(cells)
        source[cell] = feed[klass]
        undried[:, klass], dry_rate[:, klass] = _inventory(upward[:, klass], downward[:, klass], leaving[:, klass], dry[klass], source)
        dried_stock[:, klass], pyro_rate[:, klass] = _inventory(upward[:, klass], downward[:, klass], leaving[:, klass], pyro[klass], dry_rate[:, klass])
    waf = (1.0 - moisture) * (1.0 - ash)
    water = moisture * np.sum(dry_rate, axis=1)
    char = waf * yields.x_char * pyro_rate
    ash_made = (1.0 - moisture) * ash * pyro_rate
    scale = waf * np.sum(pyro_rate, axis=1)
    gas_flow = {name: equilibrium.amounts_mol.get(name, 0.0) * scale for name in PYROLYSIS_GAS_SPECIES}
    benzene, naphthalene = yields.benzene_mol_per_kg * scale, yields.naphthalene_mol_per_kg * scale
    received, dried_elements, water_elements = _compositions(moisture, ash, yields)
    mass_gap, element_gap = _closure(
        feed, cell, upward, downward, leaving, undried, dried_stock, water, char, ash_made,
        benzene, naphthalene, gas_flow, moisture, received, dried_elements, water_elements,
    )
    return ReleaseDistribution(
        water, gas_flow, benzene, naphthalene, char, ash_made,
        _system_exit(upward, downward, leaving, undried),
        (1.0 - moisture) * _system_exit(upward, downward, leaving, dried_stock),
        undried, dried_stock, mass_gap, element_gap,
    )


def _inventory(up, down, exit_rate, tau, source):
    if float(tau) == 0.0:
        return np.zeros_like(source), np.array(source, dtype=float, copy=True)
    matrix = np.diag(up + down + exit_rate + 1.0 / float(tau))
    count = up.size
    for index in range(count):
        if index + 1 < count:
            matrix[index, index + 1] -= up[index + 1]
        if index > 0:
            matrix[index, index - 1] -= down[index - 1]
    try:
        inventory = np.linalg.solve(matrix, source)
    except np.linalg.LinAlgError as exc:
        raise PyrolysisError("释放份额的线性方程奇异") from exc
    return inventory, inventory / float(tau)


def _system_exit(up, down, exit_rate, inventory):
    outflow = (exit_rate * inventory).copy()
    outflow[0, :] += up[0, :] * inventory[0, :]
    outflow[-1, :] += down[-1, :] * inventory[-1, :]
    return outflow


def _neighbor_in(up, down, inventory, index):
    arrived = np.zeros(inventory.shape[1])
    if index + 1 < inventory.shape[0]:
        arrived += up[index + 1] * inventory[index + 1]
    if index > 0:
        arrived += down[index - 1] * inventory[index - 1]
    return arrived


def _closure(feed, cell, up, down, exit_rate, undried, dried, water, char, ash_made, benzene, naphthalene, gas, moisture, received, dried_elements, water_elements):
    carbon, hydrogen = atomic_weight("C") * 1.0e-3, atomic_weight("H") * 1.0e-3
    tar_c = benzene * 6.0 * carbon + naphthalene * 10.0 * carbon
    tar_h = benzene * 6.0 * hydrogen + naphthalene * 8.0 * hydrogen
    gas_elements = {name: np.zeros(undried.shape[0]) for name in _ELEMENTS}
    for species, flow in gas.items():
        for element, count in atoms(species).items():
            if element in gas_elements:
                gas_elements[element] += flow * count * atomic_weight(element) * 1.0e-3
    mass_gap = np.zeros(undried.shape[0])
    element_gap = {name: np.zeros(undried.shape[0]) for name in _ELEMENTS}
    for index in range(undried.shape[0]):
        feed_here = feed if index == cell else np.zeros_like(feed)
        inflow_w, inflow_p = _neighbor_in(up, down, undried, index), _neighbor_in(up, down, dried, index)
        outflow_w = (up[index] + down[index] + exit_rate[index]) * undried[index]
        outflow_p = (up[index] + down[index] + exit_rate[index]) * dried[index]
        gas_mass = sum(gas_elements[name][index] for name in _ELEMENTS)
        mass_gap[index] = (
            np.sum(feed_here) + np.sum(inflow_w) + (1.0 - moisture) * np.sum(inflow_p)
            - np.sum(outflow_w) - (1.0 - moisture) * np.sum(outflow_p) - water[index]
            - np.sum(char[index]) - np.sum(ash_made[index]) - tar_c[index] - tar_h[index] - gas_mass
        )
        for name in _ELEMENTS:
            entered = np.sum(feed_here) * received[name] + np.sum(inflow_w) * received[name] + np.sum(inflow_p) * dried_elements[name]
            solid = np.sum(char[index]) if name == "C" else 0.0
            tar = tar_c[index] if name == "C" else tar_h[index] if name == "H" else 0.0
            exited = (
                np.sum(outflow_w) * received[name] + np.sum(outflow_p) * dried_elements[name]
                + water[index] * water_elements[name] + solid + tar + gas_elements[name][index]
            )
            element_gap[name][index] = entered - exited
    fuel = float(np.sum(feed))
    if fuel <= 0.0:
        raise PyrolysisError("燃料进料为 0 时没有闭合的分母")
    if np.any(np.abs(mass_gap) > _CLOSURE * fuel):
        raise PyrolysisError("释放份额的质量闭合超过 ADR 0012 的 0,1 %")
    for name, gap in element_gap.items():
        inlet = fuel * received[name]
        if np.any(np.abs(gap) > _CLOSURE * (inlet if inlet > 0.0 else fuel)):
            raise PyrolysisError(f"释放份额的 {name} 闭合超过 ADR 0012 的 0,1 %")
    return mass_gap, element_gap


def _compositions(moisture: float, ash: float, yields: LigniteYields):
    water_mass = molar_mass_kg_mol("H2O")
    water = {
        "C": 0.0, "H": 2.0 * atomic_weight("H") * 1.0e-3 / water_mass,
        "O": atomic_weight("O") * 1.0e-3 / water_mass, "N": 0.0, "S": 0.0,
    }
    waf = {"C": yields.x_c, "H": yields.x_h, "O": yields.x_o, "N": yields.x_n, "S": yields.x_s}
    dried = {name: (1.0 - moisture) * (1.0 - ash) * waf[name] for name in _ELEMENTS}
    return {name: dried[name] + moisture * water[name] for name in _ELEMENTS}, dried, water


def _array(name, value, ndim):
    array = np.asarray(value, dtype=float)
    if array.ndim != ndim or not np.all(np.isfinite(array)):
        raise PyrolysisError(f"{name} 必须是有限的 {ndim} 维数组")
    return array


def _same(name, value, shape):
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise PyrolysisError(f"{name} 的形状应为 {shape}")
    return array
