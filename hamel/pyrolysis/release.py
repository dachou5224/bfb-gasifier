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
class FuelTransport:
    """自由空间按归宿拆开的系数。cell 轴与释放函数一致，0 是最上层。

    数组末尾 bed_count 格是床区，仍用单状态系数。解释性选择：ADR 0032。
    """

    bed_count: int
    hover: np.ndarray
    free_absent: np.ndarray
    a_free_per_s: np.ndarray
    a_hover_per_s: np.ndarray
    b_free_per_s: np.ndarray
    b_hover_per_s: np.ndarray
    ejection_to_free: np.ndarray
    ejection_to_hover: np.ndarray


@dataclass(frozen=True)
class ReleaseDistribution:
    """流率按 cell。库存是当量质量；dried_exit_kg_s 是 (1−水分) 乘当量出流。

    带拆分时，自由与悬停数组在床区行上为 0，床区总量在不带后缀的数组里。
    """

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
    undried_free_kg: np.ndarray | None = None
    undried_hover_kg: np.ndarray | None = None
    dried_free_kg: np.ndarray | None = None
    dried_hover_kg: np.ndarray | None = None
    char_free_kg_s: np.ndarray | None = None
    char_hover_kg_s: np.ndarray | None = None
    ash_free_kg_s: np.ndarray | None = None
    ash_hover_kg_s: np.ndarray | None = None
    fuel_up_kg_s: np.ndarray | None = None
    fuel_down_kg_s: np.ndarray | None = None


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
    transport: FuelTransport | None = None,
) -> ReleaseDistribution:
    """稳态解 W、P。τ = 0 时该库存不设状态，入流即转化流，不做除法。

    不传 transport 时每格一个状态，原有调用不变。传入时自由空间按 ADR 0032 拆成自由与悬停。
    输运系数由调用方传入。解释性选择：ADR 0015 第 10 项。
    Source: Hamel Eq.4.23, S. 64
    Source: Hamel Eq.4.26, S. 69
    """
    if transport is not None:
        return _release_states(
            feed_kg_s, fuel_feed_height_m, cell_bottom_m, cell_top_m, k_up_per_s, k_down_per_s,
            k_exit_per_s, tau_dry_s, tau_pyr_s, moisture_mass_fraction, ash_mass_fraction_wf,
            yields, temperature_k, pressure_pa, transport,
        )
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


def _release_states(
    feed_kg_s, fuel_feed_height_m, cell_bottom_m, cell_top_m, k_up_per_s, k_down_per_s,
    k_exit_per_s, tau_dry_s, tau_pyr_s, moisture_mass_fraction, ash_mass_fraction_wf,
    yields, temperature_k, pressure_pa, transport: FuelTransport,
) -> ReleaseDistribution:
    """自由空间两状态。悬停部分没有出流，原地转化。解释性选择：ADR 0032。"""
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
    n_free = cells - int(transport.bed_count)
    if transport.bed_count < 1 or n_free < 1:
        raise PyrolysisError("拆分输运至少要有一个床区格和一个自由空间格")
    _mask("hover", transport.hover, (n_free, classes))
    _mask("free_absent", transport.free_absent, (n_free, classes))
    for name in ("a_free_per_s", "a_hover_per_s", "b_free_per_s", "b_hover_per_s"):
        _same(name, getattr(transport, name), (n_free, classes))
    _same("ejection_to_free", transport.ejection_to_free, (classes,))
    _same("ejection_to_hover", transport.ejection_to_hover, (classes,))
    if np.any(transport.free_absent & ~transport.hover):
        raise PyrolysisError("自由状态不存在的格子必须设有悬停状态")
    if np.any(np.abs(transport.ejection_to_free + transport.ejection_to_hover - 1.0) > 1.0e-12):
        raise PyrolysisError("抛出进入自由与悬停的份额之和必须为 1")
    cell = feed_cell_index(fuel_feed_height_m, cell_bottom_m, cell_top_m)
    equilibrium = pyrolysis_gas(yields, temperature_k, pressure_pa)
    undried = np.zeros((cells, classes))
    dried = np.zeros((cells, classes))
    undried_free, undried_hover = np.zeros((n_free, classes)), np.zeros((n_free, classes))
    dried_free, dried_hover = np.zeros((n_free, classes)), np.zeros((n_free, classes))
    dry_rate, pyro_rate = np.zeros((cells, classes)), np.zeros((cells, classes))
    dry_free, dry_hover = np.zeros((n_free, classes)), np.zeros((n_free, classes))
    pyro_free, pyro_hover = np.zeros((n_free, classes)), np.zeros((n_free, classes))
    for klass in range(classes):
        stock, free_stock, hover_stock, rate, free_rate, hover_rate = _two_inventories(
            upward[:, klass], downward[:, klass], leaving[:, klass], dry[klass],
            _feed_sources(cell, feed[klass], transport, n_free, klass), transport, n_free, klass,
        )
        undried[:, klass], undried_free[:, klass], undried_hover[:, klass] = stock, free_stock, hover_stock
        dry_rate[:, klass], dry_free[:, klass], dry_hover[:, klass] = rate, free_rate, hover_rate
        stock, free_stock, hover_stock, rate, free_rate, hover_rate = _two_inventories(
            upward[:, klass], downward[:, klass], leaving[:, klass], pyro[klass],
            _conversion_sources(dry_rate[:, klass], dry_free[:, klass], dry_hover[:, klass], n_free),
            transport, n_free, klass,
        )
        dried[:, klass], dried_free[:, klass], dried_hover[:, klass] = stock, free_stock, hover_stock
        pyro_rate[:, klass], pyro_free[:, klass], pyro_hover[:, klass] = rate, free_rate, hover_rate
    waf = (1.0 - moisture) * (1.0 - ash)
    water = moisture * np.sum(dry_rate, axis=1)
    char = waf * yields.x_char * pyro_rate
    ash_made = (1.0 - moisture) * ash * pyro_rate
    char_free = waf * yields.x_char * pyro_free
    char_hover = waf * yields.x_char * pyro_hover
    ash_free = (1.0 - moisture) * ash * pyro_free
    ash_hover = (1.0 - moisture) * ash * pyro_hover
    scale = waf * np.sum(pyro_rate, axis=1)
    gas_flow = {name: equilibrium.amounts_mol.get(name, 0.0) * scale for name in PYROLYSIS_GAS_SPECIES}
    benzene, naphthalene = yields.benzene_mol_per_kg * scale, yields.naphthalene_mol_per_kg * scale
    received, dried_elements, water_elements = _compositions(moisture, ash, yields)
    fuel_up, fuel_down, exit_w, exit_p = _state_flows(upward, downward, leaving, undried, dried, undried_free, dried_free, transport, n_free, moisture)
    mass_gap, element_gap = _state_closure(
        feed, cell, water, char, ash_made, benzene, naphthalene, gas_flow, moisture, received,
        dried_elements, water_elements, undried, dried, dry_rate, pyro_rate, transport, n_free,
        upward, downward, leaving, undried_free, dried_free, undried_hover, dried_hover,
    )
    return ReleaseDistribution(
        water, gas_flow, benzene, naphthalene, char, ash_made, exit_w, (1.0 - moisture) * exit_p,
        undried, dried, mass_gap, element_gap, undried_free, undried_hover, dried_free, dried_hover,
        char_free, char_hover, ash_free, ash_hover, fuel_up, fuel_down,
    )


def _feed_sources(cell, feed, transport, n_free, klass):
    """新鲜燃料只进入进料格的入口状态。床区是唯一状态，自由空间优先自由状态。"""
    if float(feed) == 0.0:
        return {}
    if cell >= n_free:
        return {(cell, "bed"): float(feed)}
    if _has_free(transport, cell, klass):
        return {(cell, "free"): float(feed)}
    if bool(transport.hover[cell, klass]):
        return {(cell, "hover"): float(feed)}
    raise PyrolysisError(f"进料格 {cell} 没有可进入的燃料状态")


def _conversion_sources(bed_rate, free_rate, hover_rate, n_free):
    """干燥产物按产生它的状态进入热解。不把悬停部分并回自由状态。ADR 0032。"""
    sources = {}
    for cell, amount in enumerate(bed_rate):
        if cell >= n_free and amount != 0.0:
            sources[(cell, "bed")] = float(amount)
    for cell, amount in enumerate(free_rate):
        if amount != 0.0:
            sources[(cell, "free")] = float(amount)
    for cell, amount in enumerate(hover_rate):
        if amount != 0.0:
            sources[(cell, "hover")] = float(amount)
    return sources


def _two_inventories(up, down, exit_rate, tau, sources, transport, n_free, klass):
    cells = up.size
    stock = np.zeros(cells)
    rate = np.zeros(cells)
    free_stock, hover_stock = np.zeros(n_free), np.zeros(n_free)
    free_rate, hover_rate = np.zeros(n_free), np.zeros(n_free)
    if float(tau) == 0.0:
        _place_instant(rate, free_rate, hover_rate, sources, transport, n_free, klass)
        return stock, free_stock, hover_stock, rate, free_rate, hover_rate
    nodes = _nodes(transport, cells, n_free, klass)
    index = {node: position for position, node in enumerate(nodes)}
    matrix = np.zeros((len(nodes), len(nodes)))
    for node in nodes:
        position = index[node]
        matrix[position, position] += 1.0 / float(tau)
        for coeff, dest in _outflows(node, up, down, exit_rate, transport, n_free, klass):
            if coeff == 0.0:
                continue
            matrix[position, position] += coeff
            if dest is not None:
                matrix[index[dest], position] -= coeff
    rhs = np.zeros(len(nodes))
    for node, amount in sources.items():
        if node not in index:
            raise PyrolysisError(f"格 {node[0]} 的{node[1]}状态不能接受转化量")
        rhs[index[node]] += float(amount)
    try:
        solved = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise PyrolysisError("释放份额的线性方程奇异") from exc
    for node, value in zip(nodes, solved):
        cell, kind = node
        stock[cell] += value
        rate[cell] += value / float(tau)
        if kind == "free":
            free_stock[cell], free_rate[cell] = value, value / float(tau)
        elif kind == "hover":
            hover_stock[cell], hover_rate[cell] = value, value / float(tau)
    return stock, free_stock, hover_stock, rate, free_rate, hover_rate


def _place_instant(rate, free_rate, hover_rate, sources, transport, n_free, klass):
    """τ = 0 时转化留在来源状态，不重新选择自由或悬停。"""
    for (cell, kind), amount in sources.items():
        if kind == "bed":
            rate[cell] += amount
        elif kind == "free":
            if cell >= n_free or not _has_free(transport, cell, klass):
                raise PyrolysisError(f"格 {cell} 没有自由状态来接受即时转化")
            free_rate[cell] += amount
            rate[cell] += amount
        elif kind == "hover":
            if cell >= n_free or not bool(transport.hover[cell, klass]):
                raise PyrolysisError(f"格 {cell} 没有悬停状态来接受即时转化")
            hover_rate[cell] += amount
            rate[cell] += amount
        else:
            raise PyrolysisError(f"不能识别的燃料状态 {kind}")


def _has_free(transport, cell, klass) -> bool:
    return not bool(transport.free_absent[cell, klass])


def _nodes(transport, cells, n_free, klass):
    nodes = []
    for cell in range(cells):
        if cell >= n_free:
            nodes.append((cell, "bed"))
            continue
        if _has_free(transport, cell, klass):
            nodes.append((cell, "free"))
        if bool(transport.hover[cell, klass]):
            nodes.append((cell, "hover"))
        if not nodes or nodes[-1][0] != cell:
            raise PyrolysisError(f"自由空间格 {cell} 粒级没有燃料状态")
    return nodes


def _outflows(node, up, down, exit_rate, transport, n_free, klass):
    cell, kind = node
    if kind == "bed":
        flows = []
        if cell == n_free:
            rise = float(up[cell])
            share_free = float(transport.ejection_to_free[klass])
            share_hover = float(transport.ejection_to_hover[klass])
            if share_free != 0.0:
                flows.append((rise * share_free, _dest_free(transport, n_free - 1, klass)))
            if share_hover != 0.0:
                flows.append((rise * share_hover, _dest_hover(transport, n_free - 1, klass)))
        elif cell > n_free:
            flows.append((float(up[cell]), (cell - 1, "bed")))
        flows.append((float(down[cell]), (cell + 1, "bed") if cell + 1 < up.size else None))
        flows.append((float(exit_rate[cell]), None))
        return flows
    if kind != "free":
        return []
    flows = []
    _add(flows, float(transport.a_free_per_s[cell, klass]), lambda: None if cell == 0 else _dest_free(transport, cell - 1, klass))
    _add(flows, float(transport.a_hover_per_s[cell, klass]), lambda: None if cell == 0 else _dest_hover(transport, cell - 1, klass))
    if cell + 1 == n_free:
        bed = (n_free, "bed")
        _add(flows, float(transport.b_free_per_s[cell, klass]), lambda: bed)
        _add(flows, float(transport.b_hover_per_s[cell, klass]), lambda: bed)
    else:
        _add(flows, float(transport.b_free_per_s[cell, klass]), lambda: _dest_free(transport, cell + 1, klass))
        _add(flows, float(transport.b_hover_per_s[cell, klass]), lambda: _dest_hover(transport, cell + 1, klass))
    return flows


def _add(flows, coeff, dest) -> None:
    if coeff != 0.0:
        flows.append((coeff, dest()))


def _dest_free(transport, cell, klass):
    if _has_free(transport, cell, klass):
        return (cell, "free")
    return _dest_hover(transport, cell, klass)


def _dest_hover(transport, cell, klass):
    if bool(transport.hover[cell, klass]):
        return (cell, "hover")
    if _has_free(transport, cell, klass):
        return (cell, "free")
    raise PyrolysisError(f"自由空间格 {cell} 没有可接受通量的状态")


def _state_flows(up, down, exit_rate, undried, dried, undried_free, dried_free, transport, n_free, moisture):
    cells, classes = undried.shape
    fuel_up, fuel_down = np.zeros((cells, classes)), np.zeros((cells, classes))
    exit_w, exit_p = np.zeros((cells, classes)), np.zeros((cells, classes))
    for klass in range(classes):
        for cell in range(n_free, cells):
            actual_w, actual_p = undried[cell, klass], (1.0 - moisture) * dried[cell, klass]
            fuel_up[cell, klass] = up[cell, klass] * (actual_w + actual_p)
            fuel_down[cell, klass] = down[cell, klass] * (actual_w + actual_p)
            exit_w[cell, klass] = exit_rate[cell, klass] * undried[cell, klass]
            exit_p[cell, klass] = exit_rate[cell, klass] * dried[cell, klass]
        for cell in range(n_free):
            coeff_up = transport.a_free_per_s[cell, klass] + transport.a_hover_per_s[cell, klass]
            coeff_down = transport.b_free_per_s[cell, klass] + transport.b_hover_per_s[cell, klass]
            actual = undried_free[cell, klass] + (1.0 - moisture) * dried_free[cell, klass]
            fuel_up[cell, klass] = coeff_up * actual
            fuel_down[cell, klass] = coeff_down * actual
            if cell == 0:
                exit_w[cell, klass] = coeff_up * undried_free[cell, klass]
                exit_p[cell, klass] = coeff_up * dried_free[cell, klass]
    exit_w[-1] += down[-1] * undried[-1]
    exit_p[-1] += down[-1] * dried[-1]
    return fuel_up, fuel_down, exit_w, exit_p


def _state_closure(
    feed, cell, water, char, ash_made, benzene, naphthalene, gas, moisture, received, dried_elements,
    water_elements, undried, dried, dry_rate, pyro_rate, transport, n_free, up, down, exit_rate,
    undried_free, dried_free, undried_hover, dried_hover,
):
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
        inflow_w, inflow_p = _state_in(index, undried, dried, undried_free, dried_free, transport, n_free, up, down)
        outflow_w, outflow_p = _state_out(index, undried, dried, undried_free, dried_free, transport, n_free, up, down, exit_rate)
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


def _state_in(index, undried, dried, undried_free, dried_free, transport, n_free, up, down):
    w = np.zeros(undried.shape[1])
    p = np.zeros_like(w)
    if index >= n_free:
        if index + 1 < undried.shape[0]:
            w += up[index + 1] * undried[index + 1]
            p += up[index + 1] * dried[index + 1]
        if index > n_free:
            w += down[index - 1] * undried[index - 1]
            p += down[index - 1] * dried[index - 1]
        elif index == n_free:
            coeff = transport.b_free_per_s[-1] + transport.b_hover_per_s[-1]
            w += coeff * undried_free[-1]
            p += coeff * dried_free[-1]
        return w, p
    if index + 1 < n_free:
        w += transport.a_free_per_s[index + 1] * undried_free[index + 1]
        p += transport.a_free_per_s[index + 1] * dried_free[index + 1]
        w += transport.a_hover_per_s[index + 1] * undried_free[index + 1]
        p += transport.a_hover_per_s[index + 1] * dried_free[index + 1]
    elif index + 1 == n_free:
        share_w = up[n_free] * undried[n_free]
        share_p = up[n_free] * dried[n_free]
        w += share_w
        p += share_p
    if index > 0:
        w += (transport.b_free_per_s[index - 1] + transport.b_hover_per_s[index - 1]) * undried_free[index - 1]
        p += (transport.b_free_per_s[index - 1] + transport.b_hover_per_s[index - 1]) * dried_free[index - 1]
    return w, p


def _state_out(index, undried, dried, undried_free, dried_free, transport, n_free, up, down, exit_rate):
    if index >= n_free:
        coeff = up[index] + down[index] + exit_rate[index]
        return coeff * undried[index], coeff * dried[index]
    coeff = transport.a_free_per_s[index] + transport.a_hover_per_s[index] + transport.b_free_per_s[index] + transport.b_hover_per_s[index]
    return coeff * undried_free[index], coeff * dried_free[index]


def _mask(name, value, shape):
    array = np.asarray(value, dtype=bool)
    if array.shape != shape:
        raise PyrolysisError(f"{name} 的形状应为 {shape}")
    return array


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
