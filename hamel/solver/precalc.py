"""外层预计算与第一次初值。顺序：流动、离析、干燥热解、灰校正。
解释性选择：ADR 0017 第 2、9、10 项。不在这里做牛顿迭代。
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from hamel.cases.assemble import assemble_case, load_case_file, resolve_case_key
from hamel.cases.config import approved_emissivity, approved_serio
from hamel.freeboard import FreeboardError, terminal_settling_velocity
from hamel.freeboard.accounting import account_freeboard
from hamel.freeboard.ejection import ejection_mass_flow, ejection_velocity_samples
from hamel.freeboard.transport import ghost_bubble_velocity
from hamel.hydrodynamics.bed import (
    BedConditions, HydroSearch, HydrodynamicsError, bed_profile, bed_upflow, geldart_group,
    local_bed_state, minimum_fluidization_velocity, particle_density_from_bed_mass,
    state_at_height, suspension_superficial_velocity,
)
from hamel.properties.data import (
    GAS_SPECIES, PropertyError, atomic_weight, atoms, gas_constant_j_mol_k, molar_mass_g_mol, molar_mass_kg_mol,
)
from hamel.properties.equilibrium import adiabatic_gas_equilibrium
from hamel.properties.thermo import (
    density_kg_m3, enthalpy_j_mol, fuel_formation_enthalpy_j_per_kg, mixture_enthalpy_j, volume_flow_m3_s,
)
from hamel.properties.transport import fuller_binary_diffusivity_m2_s, fuller_volume, interphase_diffusivity_m2_s, mixture_viscosity_pa_s
from hamel.pyrolysis.particle import PyrolysisError, lignite_particle_properties, particle_release_times
from hamel.pyrolysis.release import FuelTransport, release_distribution
from hamel.pyrolysis.yields import PYROLYSIS_GAS_SPECIES, lignite_yields, path_yields
from hamel.solver.frozen import AssemblyError, FrozenQuantities
from hamel.solver.residual import wall_area_m2
from hamel.solver.sensitivity import diameter_at, edge_diameters, interface_balance_index, pack_reported, section_volume_m3
from hamel.solver.vector import ReactorState
# ADR 0024 的冻结数值。缺字段不再回退。实现修正：ADR 0024、ADR 0039、ADR 0040。
EPSILON_MF = 0.46
PHI_S = 0.8
# ADR 0029：环境温度。壁面热损本身从工况文件读取。
AMBIENT_K = 298.15
# 求解器参数：ADR 0020。模块常量，不进 SolverConfig。数值与 ADR 0020 相同。
SEARCH = HydroSearch(1e-12, 1.0, 80, 1e-8, 2.0, 1.35, 36, 8000.0, 60, 1e-6)
DISTRIBUTOR = "technical"
SEED_CHAR = 1.0e-6
class PrecalcError(ValueError):
    """预计算的有效域不成立。报告并停止，不截断。"""
@dataclass
class PrecalcResult:
    """冻结量、第一次使用的状态，以及冒烟要列出的量。"""
    frozen: FrozenQuantities
    state: ReactorState
    t_ad_k: float
    y_ad: dict[str, float]
    geldart: str
    d_bar_m: float
    d_sv_m: float
    bottom_residual_kg_s: float
    feed_enthalpy_w: float
    equilibrium_enthalpy_w: float
    f0_a_kg_s: float
    hover_balance: tuple[tuple[int, int], ...]
    release_to_gas: float
    t_init_k: float
    hydro_group: str
    geldart_record: tuple
    yield_notes: tuple
    slip_ratios: tuple
    bed_intervals: int
    substep_history: tuple
    reported: dict
def require_geldart_b(surface_volume_diameter_m: float, rho_s_kg_m3: float, rho_f_kg_m3: float) -> str:
    """收敛后的真实类别必须是 B。求解手段：ADR 0035。"""
    group = geldart_group(surface_volume_diameter_m, rho_s_kg_m3, rho_f_kg_m3)
    if group != "B":
        raise PrecalcError(f"Geldart 判类为 {group}，不是 B 类")
    return group
def downward_mass_flow(up, feed, recycle, leaving, reaction, returned):
    """自上而下确定 Ṁ_ab。最底一格保持 0，残差只返回。
    解释性选择：ADR 0009 第 2 项。容差由调用方按 ADR 0017 第 13 项报告。
    Source: Hamel Eq.2.6, S. 17
    """
    upward, added, recycled, exited, reacted = (_column(name, value) for name, value in (
        ("up", up), ("feed", feed), ("recycle", recycle), ("leaving", leaving), ("reaction", reaction),
    ))
    count = upward.size
    if any(item.size != count for item in (added, recycled, exited, reacted)):
        raise PrecalcError("固体流递推的各数组长度必须相同")
    came_back = float(returned)
    down = np.zeros(count)
    incoming = came_back
    for index in range(count - 1):
        down[index] = added[index] + recycled[index] + upward[index + 1] + incoming - upward[index] - exited[index] + reacted[index]
        incoming = down[index]
    residual = added[-1] + recycled[-1] + incoming - upward[-1] - exited[-1] + reacted[-1]
    return down, float(residual)
def wall_coefficient(bed_temperature_k, freeboard_temperature_k, bed_height_m, freeboard_height_m, bed_diameter_m, freeboard_diameter_m, heat_loss_w, *, adiabatic: bool) -> float:
    """k_W = Q / Σ A·(T−298.15 K)。绝热对照由调用方指定为 0。
    解释性选择：ADR 0029 第 4 项。
    Source: Hamel Eq.2.7, S. 19
    """
    if adiabatic:
        return 0.0
    loss = float(heat_loss_w)
    if loss < 0.0:
        raise PrecalcError("壁面热损不能为负")
    areas = _areas(bed_temperature_k, freeboard_temperature_k, bed_height_m, freeboard_height_m, bed_diameter_m, freeboard_diameter_m)
    if areas <= 0.0:
        raise PrecalcError("Σ A_W·(T−298.15 K) 不是正数，不能求 k_W")
    return loss / areas
def initial_precalculation(path, case_id: str = "1", *, adiabatic: bool = False) -> PrecalcResult:
    """第一次外层：先做绝热平衡初值，再把这份初值当作上一轮。ADR 0017 第 9、10 项。"""
    return precalculate(path, case_id, None, adiabatic=adiabatic)
def initial_temperature_k(t_ad: float) -> float:
    """求解手段：ADR 0035。T_init = min(max(T_ad, 1073), 1273)，只用于构造初值。"""
    temperature = float(t_ad)
    if not np.isfinite(temperature):
        raise PrecalcError("T_ad 不是有限值")
    return min(max(temperature, 1073.0), 1273.0)
def _eq8_distance(d_sv, rho_s, rho_f) -> float:
    return (float(rho_s) - float(rho_f)) / 1000.0 * (float(d_sv) * 1.0e6) ** 2 - 1.0e6
def _halt(iteration, quantity, value, cell=None, klass=None):
    text = f"迭代 {iteration}：{quantity}={value}"
    if cell is not None:
        text += f"，cell={cell}"
    if klass is not None:
        text += f"，粒级={klass}"
    raise PrecalcError(text)
def precalculate(path, case_id, previous: ReactorState | None, *, reaction_solid_kg_s=None, previous_drain_kg_s: float = 0.0, previous_fuel_return_kg_s: float = 0.0, adiabatic: bool = False, iteration: int = 0, consistency: bool = False, solver=None, redirected_return_kg_s: float = 0.0, redirected_ash_out_kg_s: float = 0.0, redirected_fuel_out_kg_s: float = 0.0) -> PrecalcResult:
    """用上一轮单元解填满冻结量。previous 为空时构造初值。
    自由空间燃料回落由调用方传入，第一次外层为 0。一次风温度是混合焓的等效温度。
    consistency 为真时按真实类别重算，非 B 或产率式越界即停止。求解手段：ADR 0035。
    B 的改送回落加入床层总流；C 的改送灰和燃料记入排灰出口。默认 0，与未改送逐位相同。
    模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038。
    """
    try:
        return _precalculate(path, case_id, previous, reaction_solid_kg_s, previous_drain_kg_s, previous_fuel_return_kg_s, adiabatic, iteration, consistency, solver, redirected_return_kg_s, redirected_ash_out_kg_s, redirected_fuel_out_kg_s)
    except PrecalcError:
        raise
    except (HydrodynamicsError, PyrolysisError, PropertyError, FreeboardError) as exc:
        raise PrecalcError(f"迭代 {iteration} 停止：{exc}") from exc
def _precalculate(path, case_id, previous, reaction_solid_kg_s, previous_drain_kg_s, previous_fuel_return_kg_s, adiabatic, iteration, consistency, solver=None, redirected_return_kg_s: float = 0.0, redirected_ash_out_kg_s: float = 0.0, redirected_fuel_out_kg_s: float = 0.0) -> PrecalcResult:
    case = _load(path, case_id, solver)
    equilibrium = _equilibrium(case, adiabatic)
    t_ad = float(equilibrium.temperature)
    t_init = initial_temperature_k(t_ad)
    notes: list = []
    state = previous if previous is not None else _blank_state(case, t_init)
    fractions = _class_fractions(state if previous is not None else None, case)
    d_bar, d_sv = _diameters(case.diameters_m, fractions, case.phi_s)
    flow = _flows(case, state, equilibrium, previous is None, t_init)
    gas = _bed_gas(flow)
    rho_p, profile, conditions, intervals, substeps = _adaptive_bed(case, flow, gas, d_bar)
    profile = _exchange_diffusivity(case, flow, profile, conditions)
    group = geldart_group(d_sv, rho_p, gas.rho)
    record = (int(iteration), float(d_sv), group, _eq8_distance(d_sv, rho_p, gas.rho))
    hydro = group if consistency or group == "B" else "B"
    if consistency and group != "B":
        raise PrecalcError(
            f"迭代 {iteration}：Geldart={group}；记录（迭代号, d_sv, 类别, 式(8)距离）={record}；T_ad={t_ad}；T_init={t_init}"
        )
    share = _feed_share(case, flow, conditions, intervals)
    temperatures = state.bed_temperature_k, state.freeboard_temperature_k
    coefficient = wall_coefficient(*temperatures, case.bed_cell_m, case.free_cell_m, case.bed_diameter_m, case.free_diameter_m, case.heat_loss_w, adiabatic=adiabatic)
    if previous is None:
        state = _initial_inventory(case, flow, profile.solids_mass, share, fractions, t_init, notes, iteration, consistency)
    try:
        ejection, accounting, rho_free, mu_free, trajectory = _segregate(case, flow, profile, rho_p, d_bar)
    except FreeboardError as exc:
        raise PrecalcError(
            f"迭代 {iteration} 停止：{exc}；Geldart 记录（迭代号, d_sv, 类别, 式(8)距离）={record}；T_ad={t_ad}；T_init={t_init}"
        ) from exc
    transport, theta, hover, absent = _transport(accounting)
    if case.solver.hover_closure != "A":
        hover = np.zeros_like(hover)
        transport = (hover, *transport[1:])
    bed_u, free_u, hover_u, ratios = _slip(case, flow, profile, accounting, rho_p, hydro, gas, hover, rho_free, mu_free)
    if previous is None:
        _seed_freeboard(case, state, profile.solids_mass, t_init, hover, absent, notes, iteration, consistency)
    upflow = bed_upflow((1.0 - ejection.direct) * ejection.mass_flow, profile.internal_wake_upflow)
    masses = profile.solids_mass
    k_up = upflow / masses
    reaction = np.zeros(case.n_bed) if reaction_solid_kg_s is None else _column("reaction_solid_kg_s", reaction_solid_kg_s)
    returned = freeboard_return_kg_s(
        transport[2]["b_f"][-1], transport[2]["b_t"][-1], state.freeboard_char_free_kg[-1],
        state.freeboard_ash_free_kg[-1], previous_fuel_return_kg_s,
    ) + float(redirected_return_kg_s)
    down, residual = bed_downflow(upflow, case.fuel_kg_s, reaction, returned, previous_drain_kg_s)
    k_down = down / masses
    k_exit = np.zeros_like(k_down)
    if previous_drain_kg_s:
        k_exit[-1] = float(previous_drain_kg_s) / masses[-1]
    released, pyrolysis = _release(case, state, flow, profile, transport, k_up, k_down, k_exit, hydro, notes, iteration, consistency)
    fuel_bed = released.undried_kg[case.n_free:] + (1.0 - float(case.physical.moisture_mass_fraction)) * released.dried_kg[case.n_free:]
    inert = masses - state.bed_char_kg.sum(axis=1) - fuel_bed.sum(axis=1)
    if np.any(inert < 0.0):
        raise PrecalcError("几何存量扣去焦炭和未转化燃料后灰分为负")
    drain_ash, drain_char = _drain(
        case, state, released, inert, masses, accounting, redirected_ash_out_kg_s, redirected_fuel_out_kg_s,
    )
    frozen = _frozen(case, state, flow, profile, share, coefficient, rho_p, k_up, down, transport, theta, hover, absent, released, drain_ash, drain_char, bed_u, free_u, hover_u)
    gas_out = _gas_release(released)
    reported = pack_reported(profile, flow, case, released)
    reported["bed_gas_temperature_deviation_k"] = float(np.max(np.abs(flow.bed_temperature - gas.temperature)))
    reported["reference_return_share"] = np.array([book.return_share for book in accounting], dtype=float)
    reported["reference_elutriation_share"] = np.array([book.elutriation_share for book in accounting], dtype=float)
    reported["phi_s"] = float(case.phi_s)
    reported["n_0"] = float(case.n_0)
    reported["cell_division"] = f"{case.n_bed}×{case.bed_cell_m:g} m + {case.n_free}×{case.free_cell_m:g} m"
    free = case.n_free
    reported["moisture_mass_fraction"] = float(pyrolysis["moisture_mass_fraction"])
    reported["undried_free_kg"] = np.asarray(released.undried_free_kg, dtype=float)
    reported["dried_free_kg"] = np.asarray(released.dried_free_kg, dtype=float)
    reported["undried_hover_kg"] = np.asarray(released.undried_hover_kg, dtype=float)
    reported["dried_hover_kg"] = np.asarray(released.dried_hover_kg, dtype=float)
    reported["undried_bed_kg"] = np.asarray(released.undried_kg[free:], dtype=float)
    reported["dried_bed_kg"] = np.asarray(released.dried_kg[free:], dtype=float)
    reported["pyrolysis"] = pyrolysis
    reported["trajectory"] = trajectory
    return PrecalcResult(
        frozen, state, t_ad, dict(equilibrium.fractions), group, d_bar, d_sv, residual,
        equilibrium.feed_enthalpy, equilibrium.product_enthalpy, ejection.mass_flow, _hover_pairs(hover),
        gas_out / case.fuel_kg_s, t_init, hydro, record, tuple(notes), ratios, intervals, substeps,
        reported,
    )
@dataclass
class _Case:
    physical: object
    solver: object
    diameters_m: np.ndarray
    class_fraction: np.ndarray
    bed_mass_kg: float
    heat_loss_w: float
    n_0: float
    h_a: float
    solid_temperature_k: float
    purge_temperature_k: float
    n_bed: int
    n_free: int
    bed_cell_m: float
    free_cell_m: float
    epsilon_mf: float
    phi_s: float
    secondary_height_m: float
    @property
    def fuel_kg_s(self) -> float:
        return float(self.physical.fuel_mass_flow)
    @property
    def bed_diameter_m(self) -> float:
        return float(self.physical.bed_diameter)
    @property
    def free_diameter_m(self) -> float:
        return float(self.physical.freeboard_diameter)
    @property
    def pressure_pa(self) -> float:
        return float(self.physical.pressure)
@dataclass
class _Equilibrium:
    temperature: float
    fractions: dict[str, float]
    total_mol_s: float
    primary_mol_s: np.ndarray
    secondary_mol_s: np.ndarray
    primary_total: float
    secondary_total: float
    feed_enthalpy: float
    product_enthalpy: float
    amounts: dict[str, float]
@dataclass
class _Flow:
    u0: np.ndarray
    bed_mol_s: np.ndarray
    free_mol_s: np.ndarray
    bed_temperature: np.ndarray
    free_temperature: np.ndarray
    pressure: float
@dataclass
class _Gas:
    rho: float
    mu: float
    nu: float
    diffusivity: float
    temperature: float
    fractions: dict[str, float]
@dataclass
class _Ejection:
    mass_flow: float
    speeds: np.ndarray
    shares: np.ndarray
    direct: float
def _load(path, case_id: str, solver=None) -> _Case:
    config = assemble_case(path, case_id)
    used = config.solver if solver is None else solver
    doc = load_case_file(path)
    key = resolve_case_key(doc, case_id)
    heights = np.asarray(_value(doc, key, "cell_heights"), dtype=float)
    bed = float(heights[0])
    n_bed = int(np.sum(np.isclose(heights, bed)))
    if n_bed < 1 or not np.all(np.isclose(heights[:n_bed], bed)):
        raise PrecalcError("床区格子高度必须相同，且排在前面")
    free = heights[n_bed:]
    if free.size < 1 or not np.all(np.isclose(free, free[0])):
        raise PrecalcError("自由空间格子高度必须相同")
    return _Case(
        config.physical, used, np.asarray(_value(doc, key, "fuel_size_class_diameters"), dtype=float),
        np.asarray(_value(doc, key, "fuel_size_class_mass_fractions"), dtype=float),
        float(_value(doc, key, "bed_mass")), float(_value(doc, key, "wall_heat_loss")),
        float(_value(doc, key, "distributor_orifice_count")), float(_value(doc, key, "distributor_outlet_height")),
        float(_value(doc, key, "solid_feed_temperature")), float(_value(doc, key, "nitrogen_purge_temperature")),
        n_bed, int(free.size), bed, float(free[0]),
        float(_value(doc, key, "bed_minimum_fluidization_voidage")),
        float(_value(doc, key, "bed_sphericity")),
        float(_value(doc, key, "agent_feed_height")),
    )
def _value(doc, key, name):
    case = doc["cases"][key]
    fuel = doc.get("fuels", {}).get(case.get("fuel_id"), {})
    for container in (case, fuel, doc.get("shared", {})):
        obj = container.get(name) if isinstance(container, dict) else None
        if isinstance(obj, dict) and "status" in obj:
            if obj.get("status") == "missing" or obj.get("value") is None:
                raise PrecalcError(f"{name} 没有数值")
            return obj["value"]
    raise PrecalcError(f"工况文件没有 {name}")
def _equilibrium(case: _Case, adiabatic: bool) -> _Equilibrium:
    """全炉进料的气相绝热平衡。基准运行扣除壁面热损，绝热对照不扣。焦炭不从元素池扣除。
    吹扫氮按 ADR 0014 C2 用 298.15 K。冻结量的一次风温度是使混合进气焓与此相同的等效温度。
    """
    physical = case.physical
    air = np.array([0.21, 0.79]) * float(physical.primary_air_amount_flow)
    steam = float(physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    purge = float(physical.nitrogen_purge_volume_flow)
    primary = {name: 0.0 for name in GAS_SPECIES}
    primary["O2"], primary["N2"], primary["H2O"] = air[0], air[1] + purge, steam
    secondary_n = float(physical.secondary_air_amount_flow)
    secondary = {name: 0.0 for name in GAS_SPECIES}
    secondary["O2"], secondary["N2"] = 0.21 * secondary_n, 0.79 * secondary_n
    elements = _elements(case, primary, secondary)
    inlet = float(physical.inlet_temperature)
    enthalpy = _gas_enthalpy(primary, inlet) + _gas_enthalpy({"N2": purge}, case.purge_temperature_k) - _gas_enthalpy({"N2": purge}, inlet)
    enthalpy += _gas_enthalpy(secondary, inlet) + _fuel_enthalpy(case)
    heat = 0.0 if adiabatic else case.heat_loss_w
    result = adiabatic_gas_equilibrium(enthalpy - heat, case.pressure_pa, elements)
    product = mixture_enthalpy_j(result.temperature_k, result.amounts_mol)
    return _Equilibrium(
        float(result.temperature_k), {name: float(result.mole_fractions.get(name, 0.0)) for name in GAS_SPECIES},
        float(result.total_mol), _vector(primary), _vector(secondary), float(sum(primary.values())), secondary_n,
        enthalpy, product, {name: float(result.amounts_mol.get(name, 0.0)) for name in GAS_SPECIES},
    )
def _elements(case, primary, secondary) -> dict[str, float]:
    total = {name: 0.0 for name in ("C", "H", "O", "N", "S")}
    for mixture in (primary, secondary):
        for name, moles in mixture.items():
            for element, count in atoms(name).items():
                total[element] = total.get(element, 0.0) + count * moles
    physical = case.physical
    dry = case.fuel_kg_s * (1.0 - float(physical.moisture_mass_fraction))
    moisture = case.fuel_kg_s * float(physical.moisture_mass_fraction)
    weights = {element: atomic_weight(element) * 1.0e-3 for element in total}
    total["C"] += dry * physical.carbon_mass_fraction_wf / weights["C"]
    total["H"] += dry * physical.hydrogen_mass_fraction_wf / weights["H"] + 2.0 * moisture / molar_mass_kg_mol("H2O")
    total["O"] += dry * physical.oxygen_mass_fraction_wf / weights["O"] + moisture / molar_mass_kg_mol("H2O")
    total["N"] += dry * physical.nitrogen_mass_fraction_wf / weights["N"]
    total["S"] += dry * physical.sulfur_mass_fraction_wf / weights["S"]
    return {name: value for name, value in total.items() if value > 0.0}
def _fuel_enthalpy(case: _Case) -> float:
    physical = case.physical
    moisture = float(physical.moisture_mass_fraction)
    formation = fuel_formation_enthalpy_j_per_kg(
        physical.heating_value, carbon=physical.carbon_mass_fraction_wf * (1.0 - moisture),
        hydrogen=physical.hydrogen_mass_fraction_wf * (1.0 - moisture),
        nitrogen=physical.nitrogen_mass_fraction_wf * (1.0 - moisture),
        sulfur=physical.sulfur_mass_fraction_wf * (1.0 - moisture), moisture=moisture,
    )
    capacity = lignite_particle_properties().heat_capacity_j_kg_k
    return case.fuel_kg_s * (formation + capacity * (case.solid_temperature_k - AMBIENT_K))
def _gas_enthalpy(flows: dict[str, float], temperature: float) -> float:
    return float(sum(amount * enthalpy_j_mol(name, temperature) for name, amount in flows.items() if amount))
def _blank_state(case: _Case, temperature: float) -> ReactorState:
    species = len(GAS_SPECIES)
    return ReactorState(
        np.zeros((case.n_bed, 2, species)), np.zeros((case.n_free, species)),
        np.full(case.n_bed, temperature), np.full(case.n_free, temperature),
        np.zeros((case.n_bed, case.diameters_m.size)), np.zeros((case.n_free, case.diameters_m.size)),
        np.zeros((case.n_free, case.diameters_m.size)), np.zeros((case.n_free, case.diameters_m.size)),
        np.zeros((case.n_free, case.diameters_m.size)),
    )
def _class_fractions(previous, case: _Case) -> np.ndarray:
    if previous is None:
        total = float(case.class_fraction.sum())
        if total <= 0.0:
            raise PrecalcError("燃料粒级质量分数之和不是正数")
        return case.class_fraction / total
    mass = np.asarray(previous.bed_char_kg, dtype=float).sum(axis=0)
    if float(mass.sum()) <= 0.0:
        raise PrecalcError("上一轮床区焦炭为 0，不能求粒径")
    return mass / float(mass.sum())
def _diameters(diameter_m, fraction, phi_s) -> tuple[float, float]:
    if diameter_m.shape != fraction.shape:
        raise PrecalcError("粒径与质量分数的粒级数不一致")
    d_bar = float(np.sum(fraction * diameter_m))
    d_sv = 1.0 / float(np.sum(fraction / (float(phi_s) * diameter_m)))
    return d_bar, d_sv
def _flows(case, state, equilibrium, initial: bool, t_init: float) -> _Flow:
    pressure = case.pressure_pa
    if initial:
        bed_t = np.full(case.n_bed, t_init)
        free_t = np.full(case.n_free, t_init)
        gas_in = equilibrium.primary_total + equilibrium.secondary_total
        if gas_in <= 0.0:
            raise PrecalcError("进气摩尔流不是正数")
        kept = equilibrium.total_mol_s * equilibrium.primary_total / gas_in
        bed = np.tile(_scaled(equilibrium.fractions, kept), (case.n_bed, 1))
        free = np.zeros((case.n_free, len(GAS_SPECIES)))
        index = _secondary_index(case)
        for cell in range(case.n_free):
            amount = equilibrium.total_mol_s if cell <= index else kept
            free[cell] = _scaled(equilibrium.fractions, amount)
    else:
        bed_t = np.asarray(state.bed_temperature_k, dtype=float)
        free_t = np.asarray(state.freeboard_temperature_k, dtype=float)
        bed = np.asarray(state.bed_gas_mol_s, dtype=float).sum(axis=1)
        free = np.asarray(state.freeboard_gas_mol_s, dtype=float)
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2
    u0 = np.array([
        volume_flow_m3_s(float(bed[cell].sum()), float(bed_t[cell]), pressure) / area for cell in range(case.n_bed)
    ])
    return _Flow(u0, bed, free, bed_t, free_t, pressure)
def _scaled(fractions, amount) -> np.ndarray:
    return np.array([fractions[name] * amount for name in GAS_SPECIES])
def _secondary_index(case) -> int:
    try:
        return interface_balance_index(
            case.secondary_height_m, case.n_bed * case.bed_cell_m, case.n_free, case.free_cell_m,
        )
    except AssemblyError as exc:
        raise PrecalcError(str(exc)) from exc
def _cell_diffusivity(mode: str, temperature: float, pressure: float) -> float:
    """该格气体温度下的 D_g。流体力学物性不在这里换。实现修正：ADR 0021、ADR 0040。"""
    if mode == "hamel_5_20":
        return interphase_diffusivity_m2_s(temperature, pressure)
    if mode == "fuller_o2_n2":
        return fuller_binary_diffusivity_m2_s(
            temperature, pressure, molar_mass_g_mol("O2"), molar_mass_g_mol("N2"),
            fuller_volume("O2"), fuller_volume("N2"),
        )
    raise PrecalcError(f"不能识别的相间扩散 {mode}")
def _exchange_diffusivity(case, flow, profile, conditions):
    """各床格 K_bd 用该格气体温度下的 D_g，写法与 S5 逐格分支相同。
    ρ_g、μ_g、ν_g、气泡直径和 u_mf 仍是床区平均值上的剖面。S5 只改公式。
    实现修正：ADR 0021、ADR 0040。
    """
    mode = case.solver.interphase_diffusivity
    k_bd = np.array(profile.k_bd, dtype=float, copy=True)
    for index in range(k_bd.size):
        diffusivity = _cell_diffusivity(mode, float(flow.bed_temperature[index]), float(flow.pressure))
        local = local_bed_state(
            float(flow.u0[index]), float(profile.d_b[index]), float(profile.u_mf),
            replace(conditions, diffusivity=diffusivity), SEARCH,
        )
        k_bd[index] = local.k_bd
    return replace(profile, k_bd=k_bd)
def _bed_gas(flow: _Flow) -> _Gas:
    weights = flow.bed_mol_s.sum(axis=1)
    if float(weights.sum()) <= 0.0:
        raise PrecalcError("床区气体流量不是正数")
    temperature = float(np.dot(weights, flow.bed_temperature) / weights.sum())
    total = flow.bed_mol_s.sum(axis=0)
    fractions = {name: float(total[index] / total.sum()) for index, name in enumerate(GAS_SPECIES)}
    rho = density_kg_m3(temperature, flow.pressure, fractions)
    mu = mixture_viscosity_pa_s(fractions, temperature, flow.pressure)
    return _Gas(rho, mu, mu / rho, interphase_diffusivity_m2_s(temperature, flow.pressure), temperature, fractions)
def _density(case, flow, gas, d_bar, intervals) -> float:
    """反求 ρ_P 时先用 B 类系数。求解手段：ADR 0035。真实类别另记。ADR 0014 A4。"""
    conditions = _conditions(case, gas, d_bar, gas.rho * SEARCH.density_start_factor, "B")
    edges = np.linspace(0.0, case.n_bed * case.bed_cell_m, case.n_bed + 1)
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2
    return particle_density_from_bed_mass(
        case.bed_mass_kg, edges, flow.u0, area, case.n_0, case.h_a, intervals, case.phi_s, conditions, SEARCH,
    )
def _conditions(case, gas, d_bar, rho_p, group) -> BedConditions:
    return BedConditions(
        d_bar, case.epsilon_mf, rho_p, gas.rho, gas.mu, gas.nu, case.pressure_pa, gas.diffusivity,
        case.bed_diameter_m, DISTRIBUTOR, group, True,
    )
def _profile(case, flow, gas, d_bar, rho_p, group, intervals):
    conditions = _conditions(case, gas, d_bar, rho_p, group)
    minimum = minimum_fluidization_velocity(case.epsilon_mf, gas.nu, case.phi_s, d_bar, rho_p, gas.rho)
    edges = np.linspace(0.0, case.n_bed * case.bed_cell_m, case.n_bed + 1)
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2
    profile = bed_profile(edges, flow.u0, minimum, area, case.n_0, case.h_a, intervals, conditions, SEARCH)
    return profile, conditions
def _feed_share(case, flow, conditions, intervals) -> float:
    """布风器处可见气泡流占式 (3.7) 总流的份额，其余进悬浮相。穿流不另作一股进料。
    Source: Hamel Eq.3.7, S. 24
    Source: Hamel Eq.3.8, S. 24
    """
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2
    minimum = minimum_fluidization_velocity(conditions.epsilon_mf, conditions.nu_g, case.phi_s, conditions.d_p, conditions.rho_p, conditions.rho_g)
    local = state_at_height(0.0, float(flow.u0[-1]), minimum, area, case.n_0, case.h_a, intervals, conditions, SEARCH)
    bubble = local.epsilon_b * local.u_b / float(flow.u0[-1])
    if not 0.0 <= bubble <= 1.0:
        raise PrecalcError("进料气泡份额不在 [0, 1]")
    return 1.0 - bubble
def _fuel(case):
    physical = case.physical
    return (
        physical.carbon_mass_fraction_wf, physical.hydrogen_mass_fraction_wf, physical.oxygen_mass_fraction_wf,
        physical.nitrogen_mass_fraction_wf, physical.sulfur_mass_fraction_wf, physical.ash_mass_fraction_wf,
    )
def _yield_at(case, temperature, cell, klass, notes, iteration, consistency):
    """只检查这次实际代入产率式的温度和压力。求解手段：ADR 0035。"""
    pressure = float(case.pressure_pa)
    value = float(temperature)
    if not (1073.0 <= value <= 1273.0 and 1.0e5 <= pressure <= 2.5e6):
        notes.append((int(iteration), "产率式温度或压力", value, cell, klass))
        if consistency:
            _halt(iteration, "产率式温度或压力", value, cell, klass)
        return path_yields(pressure, value, *_fuel(case))
    return lignite_yields(pressure, value, *_fuel(case))
def _initial_inventory(case, flow, masses, share, fractions, temperature, notes, iteration, consistency) -> ReactorState:
    yields = _yield_at(case, temperature, None, None, notes, iteration, consistency)
    ash = _ash_fraction(case, yields.x_char)
    char = masses[:, None] * (1.0 - ash) * np.asarray(fractions)[None, :]
    state = _blank_state(case, temperature)
    state.bed_char_kg = char
    state.bed_gas_mol_s[:, 0, :] = flow.bed_mol_s * share
    state.bed_gas_mol_s[:, 1, :] = flow.bed_mol_s * (1.0 - share)
    state.freeboard_gas_mol_s = flow.free_mol_s
    return state
def _ash_fraction(case, x_char: float) -> float:
    """w_A,0 = x_A / (x_A + x_Koks,B·(1−x_w−x_A))。ADR 0017 第 9 项。"""
    moisture = float(case.physical.moisture_mass_fraction)
    ash = float(case.physical.ash_mass_fraction_wf) * (1.0 - moisture)
    denominator = ash + x_char * (1.0 - moisture - ash)
    if denominator <= 0.0 or not 0.0 <= ash / denominator < 1.0:
        raise PrecalcError("初值灰分比例没有定义")
    return ash / denominator
def _segregate(case, flow, profile, rho_p, d_bar):
    edges = np.linspace(0.0, case.n_free * case.free_cell_m, case.n_free + 1)
    diameters = edge_diameters(case.n_free, case.free_cell_m, case.bed_diameter_m, case.free_diameter_m)
    free_balance = flow.free_mol_s[::-1]
    temperature = flow.free_temperature[::-1]
    flows = np.array([
        volume_flow_m3_s(float(free_balance[cell].sum()), float(temperature[cell]), flow.pressure) for cell in range(case.n_free)
    ])
    rho = np.array([
        density_kg_m3(float(temperature[cell]), flow.pressure, _fractions(free_balance[cell])) for cell in range(case.n_free)
    ])
    mu = np.array([
        mixture_viscosity_pa_s(_fractions(free_balance[cell]), float(temperature[cell]), flow.pressure) for cell in range(case.n_free)
    ])
    def velocity(cell, height, flows=flows, edges=edges, diameters=diameters):
        return ghost_bubble_velocity(cell, height, edges, diameters, flows, float(profile.u_b[0]), float(profile.d_b[0]))
    speeds, shares, direct = ejection_velocity_samples(float(profile.u_b[0]), case.solver.freeboard_ejection_samples)
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2
    mass = ejection_mass_flow(float(profile.epsilon_b[0]), float(profile.u_b[0]), float(profile.d_b[0]), rho_p, float(profile.epsilon_d[0]), d_bar, area)
    books = []
    for klass, diameter in enumerate(case.diameters_m):
        args = (edges, velocity, rho, mu, float(diameter), rho_p, case.phi_s, speeds, shares, direct, case.solver)
        books.append(account_freeboard(*args))
    column = {
        "edges_m": edges, "volume_flow_m3_s": flows, "edge_diameter_m": diameters,
        "rho_g": rho, "mu_g": mu, "bubble_m_s": float(profile.u_b[0]),
        "bubble_diameter_m": float(profile.d_b[0]), "rho_p": float(rho_p),
        "phi_s": float(case.phi_s), "diameter_m": np.asarray(case.diameters_m, dtype=float),
    }
    return _Ejection(mass, speeds, shares, direct), tuple(books), rho, mu, column
def _fractions(flow) -> dict[str, float]:
    total = float(np.asarray(flow).sum())
    if total <= 0.0:
        raise PrecalcError("自由空间气体流量不是正数")
    return {name: float(flow[index] / total) for index, name in enumerate(GAS_SPECIES)}
def _transport(books):
    n_free, n_class = len(books[0].cells), len(books)
    theta = np.zeros((n_free, n_class))
    hover = np.zeros((n_free, n_class), dtype=bool)
    absent = np.zeros_like(hover)
    arrays = {name: np.zeros((n_free, n_class)) for name in ("a_f", "a_t", "b_f", "b_t")}
    to_free = np.zeros(n_class)
    to_hover = np.zeros(n_class)
    for klass, book in enumerate(books):
        launched = book.ejection_free + book.ejection_hover
        if launched <= 0.0:
            to_free[klass], to_hover[klass] = 1.0, 0.0
        else:
            to_free[klass], to_hover[klass] = book.ejection_free / launched, book.ejection_hover / launched
        for cell, item in enumerate(book.cells):
            balance = n_free - 1 - cell
            missing = item.free_state == "自由状态不存在"
            if missing and item.theta != 1:
                raise PrecalcError(f"粒级 {klass} 格 {balance} 自由状态不存在但 θ 不是 1")
            absent[balance, klass] = missing
            theta[balance, klass] = 1.0 if missing else float(item.theta)
            # θ = 0 仍可能有稳定根。悬停由捕获区决定；θ = 1 时冻结量不允许没有悬停状态。
            hover[balance, klass] = bool(missing or cell in book.capture_cells or item.theta == 1)
            if not missing:
                arrays["a_f"][balance, klass] = item.a_f
                arrays["a_t"][balance, klass] = item.a_t
                arrays["b_f"][balance, klass] = item.b_f
                arrays["b_t"][balance, klass] = item.b_t
        if absent[-1, klass]:
            to_free[klass], to_hover[klass] = 0.0, 1.0
    return (hover, absent, arrays, to_free, to_hover, theta), theta, hover, absent
def freeboard_return_kg_s(b_free_bottom, b_hover_bottom, char_free_kg, ash_free_kg, fuel_return_kg_s: float = 0.0) -> float:
    """床面回落 = 最下自由格 (b^f+b^t) 乘上一轮可输运存量，再加未转化燃料回落。
    直接回落不进入自由空间，无反应参考份额不参与递推。ADR 0016 第 10 条；ADR 0032。
    """
    coeff = np.asarray(b_free_bottom, dtype=float) + np.asarray(b_hover_bottom, dtype=float)
    mobile = np.asarray(char_free_kg, dtype=float) + np.asarray(ash_free_kg, dtype=float)
    if coeff.shape != mobile.shape:
        raise PrecalcError("回落系数与自由空间存量的粒级不一致")
    return float(np.sum(coeff * mobile) + float(fuel_return_kg_s))
def bed_downflow(up, fuel_feed_kg_s, reaction, returned, drain_kg_s):
    """整体固体向下流。最底格排灰进入守恒残差，不改写 Ṁ_ab = 0。ADR 0009 第 2 项。"""
    upward = _column("up", up)
    added = np.zeros(upward.size)
    added[-1] = float(fuel_feed_kg_s)
    leaving = np.zeros(upward.size)
    leaving[-1] = float(drain_kg_s)
    reacted = np.zeros(upward.size) if reaction is None else _column("reaction", reaction)
    return downward_mass_flow(upward, added, np.zeros(upward.size), leaving, reacted, returned)
def _release(case, state, flow, profile, packed, k_up, k_down, k_exit, group, notes, iteration, consistency):
    hover, absent, arrays, to_free, to_hover, _theta = packed
    transport = FuelTransport(
        case.n_bed, hover, absent, arrays["a_f"], arrays["a_t"], arrays["b_f"], arrays["b_t"], to_free, to_hover,
        closure=case.solver.hover_closure,
    )
    temperature = float(np.average(flow.bed_temperature, weights=flow.bed_mol_s.sum(axis=1)))
    yields = _yield_at(case, temperature, None, None, notes, iteration, consistency)
    relative = float(np.mean(_bed_slip(flow, profile, group)))
    tau_dry, tau_pyr = [], []
    for diameter in case.diameters_m:
        times = particle_release_times(
            float(diameter), temperature, flow.pressure, float(case.physical.moisture_mass_fraction),
            case.solid_temperature_k, _bed_gas(flow).fractions, relative,
            float(case.physical.volatiles_mass_fraction_waf), lignite_particle_properties(), case.solver,
        )
        tau_dry.append(times.tau_dry_s)
        tau_pyr.append(times.tau_pyr_s)
    cells = case.n_free + case.n_bed
    zeros = np.zeros((cells, case.diameters_m.size))
    upward, downward, leaving = zeros.copy(), zeros.copy(), zeros.copy()
    upward[case.n_free:] = k_up[:, None]
    downward[case.n_free:] = k_down[:, None]
    leaving[case.n_free:] = k_exit[:, None]
    bottoms, tops = _bounds(case)
    result = release_distribution(
        case.fuel_kg_s * case.class_fraction, float(case.physical.fuel_feed_height), bottoms, tops,
        upward, downward, leaving, np.array(tau_dry), np.array(tau_pyr),
        float(case.physical.moisture_mass_fraction), float(case.physical.ash_mass_fraction_wf),
        yields, temperature, flow.pressure,
        transport=transport,
    )
    inputs = {
        "diameter_m": np.asarray(case.diameters_m, dtype=float),
        "temperature_k": temperature, "pressure_pa": float(flow.pressure),
        "moisture_mass_fraction": float(case.physical.moisture_mass_fraction),
        "feed_temperature_k": float(case.solid_temperature_k),
        "gas": {name: float(value) for name, value in _bed_gas(flow).fractions.items()},
        "relative_velocity_m_s": relative,
        "ultimate_volatile_mass_fraction": float(case.physical.volatiles_mass_fraction_waf),
    }
    return result, inputs
def _bounds(case):
    bed_top = case.n_bed * case.bed_cell_m
    free_dz, bed_dz = case.free_cell_m, case.bed_cell_m
    bottoms = np.empty(case.n_free + case.n_bed)
    tops = np.empty_like(bottoms)
    for cell in range(case.n_free):
        geometric = case.n_free - 1 - cell
        bottoms[cell] = bed_top + geometric * free_dz
        tops[cell] = bottoms[cell] + free_dz
    for cell in range(case.n_bed):
        geometric = case.n_bed - 1 - cell
        bottoms[case.n_free + cell] = geometric * bed_dz
        tops[case.n_free + cell] = bottoms[case.n_free + cell] + bed_dz
    return bottoms, tops
def _drain(case, state, released, inert, masses, books, redirected_ash_out_kg_s: float = 0.0, redirected_fuel_out_kg_s: float = 0.0):
    """排灰灰分量 = 燃料带入灰 − 带出灰。ADR 0014 A5。
    C 的改送灰和改送燃料是外部出口，从余额扣除。两项为 0 时与原式相同。
    模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038。
    """
    moisture = float(case.physical.moisture_mass_fraction)
    ash_in = case.fuel_kg_s * float(case.physical.ash_mass_fraction_wf) * (1.0 - moisture)
    top_free = np.asarray(state.freeboard_ash_free_kg[0], dtype=float)
    coeff = np.array([
        (0.0 if book.cells[-1].a_f is None else book.cells[-1].a_f) + (0.0 if book.cells[-1].a_t is None else book.cells[-1].a_t)
        for book in books
    ])
    # cells[-1] 是自由空间顶（edges 从床面向上，最后一格是顶）。衡算下标 0 是顶。
    carried = float(np.sum(top_free * coeff))
    fuel_top = 0.0 if released.fuel_up_kg_s is None else float(released.fuel_up_kg_s[0].sum())
    fuel_drain = float(released.undried_exit_kg_s[case.n_free:].sum() + released.dried_exit_kg_s[case.n_free:].sum())
    ash_out = carried + ash_in / case.fuel_kg_s * (fuel_top + fuel_drain + float(redirected_fuel_out_kg_s))
    ash_out += float(redirected_ash_out_kg_s)
    if ash_out > ash_in:
        raise PrecalcError("带出灰多于燃料带入灰")
    fraction = float(inert[-1] / masses[-1])
    if fraction <= 0.0:
        raise PrecalcError("最底格灰分质量分数不是正数")
    total = (ash_in - ash_out) / fraction
    char = total * np.asarray(state.bed_char_kg[-1], dtype=float) / float(masses[-1])
    return ash_in - ash_out, char
def _frozen(case, state, flow, profile, share, coefficient, rho_p, k_up, down, packed, theta, hover, absent, released, drain_ash, drain_char, bed_u, free_u, hover_u):
    hover_mask, absent_mask, arrays, to_free, to_hover, theta_used = packed
    bed, free = case.n_bed, case.n_free
    pyro_bed, pyro_free = _pyro_flows(released, bed, free)
    fuel_bed = released.undried_kg[free:] + (1.0 - float(case.physical.moisture_mass_fraction)) * released.dried_kg[free:]
    fuel_up = np.zeros_like(fuel_bed) if released.fuel_up_kg_s is None else released.fuel_up_kg_s[free:]
    fuel_down = np.zeros_like(fuel_bed) if released.fuel_down_kg_s is None else released.fuel_down_kg_s[free:]
    fuel_exit = released.undried_exit_kg_s[free:] + released.dried_exit_kg_s[free:]
    char_free = np.zeros((free, case.diameters_m.size)) if released.char_free_kg_s is None else released.char_free_kg_s
    char_hover = np.zeros_like(char_free) if released.char_hover_kg_s is None else released.char_hover_kg_s
    ash_free = np.zeros_like(char_free) if released.ash_free_kg_s is None else released.ash_free_kg_s
    ash_hover = np.zeros_like(char_free) if released.ash_hover_kg_s is None else released.ash_hover_kg_s
    formation = _formation(case)
    emissivity = _approved(approved_emissivity, case.solver.particle_emissivity_case)
    frequency, activation = _approved(approved_serio, case.solver.serio_pair)
    return FrozenQuantities(
        diameter_m=case.diameters_m, bed_pressure_pa=np.full(bed, flow.pressure), freeboard_pressure_pa=np.full(free, flow.pressure),
        k_auf_per_s=k_up, m_dot_ab_kg_s=down, geometric_mass_kg=profile.solids_mass, k_bd_per_s=profile.k_bd,
        v_bubble_exchange_m3=_bubble_volume(case, profile), v_suspension_m3=_suspension_volume(case, profile),
        v_bubble_m3=_bubble_volume(case, profile), v_freeboard_m3=_free_volume(case),
        primary_feed_mol_s=_primary_vector(case), feed_to_suspension=share,
        primary_temperature_k=_primary_temperature(case), secondary_feed_mol_s=_secondary_vector(case),
        secondary_temperature_k=float(case.physical.inlet_temperature), secondary_balance_index=_secondary_index(case),
        pyro_bed_mol_s=pyro_bed, pyro_freeboard_mol_s=pyro_free,
        char_gen_bed_kg_s=released.char_kg_s[free:], char_gen_free_kg_s=char_free, char_gen_hover_kg_s=char_hover,
        ash_gen_bed_kg_s=released.ash_kg_s[free:], ash_gen_free_kg_s=ash_free, ash_gen_hover_kg_s=ash_hover,
        fuel_bed_kg=fuel_bed, fuel_up_bed_kg_s=fuel_up, fuel_down_bed_kg_s=fuel_down, fuel_exit_bed_kg_s=fuel_exit,
        fuel_up_free_kg_s=released.fuel_up_kg_s[:free], fuel_down_free_kg_s=released.fuel_down_kg_s[:free],
        fuel_feed_kg_s=case.fuel_kg_s * case.class_fraction, fuel_formation_j_per_kg=formation,
        fuel_heat_capacity_j_per_kg_k=lignite_particle_properties().heat_capacity_j_kg_k,
        fuel_feed_temperature_k=case.solid_temperature_k, hover=hover_mask, free_absent=absent_mask,
        a_free_per_s=arrays["a_f"], a_hover_per_s=arrays["a_t"], b_free_per_s=arrays["b_f"], b_hover_per_s=arrays["b_t"],
        theta=theta_used, ejection_to_free=to_free, ejection_to_hover=to_hover,
        u_r_bed_m_s=bed_u, u_r_freeboard_m_s=free_u, u_r_hover_m_s=hover_u,
        rho_p_kg_m3=rho_p, drain_char_kg_s=drain_char, drain_ash_kg_s=drain_ash, k_w_w_m2_k=coefficient,
        bed_height_m=case.bed_cell_m, freeboard_height_m=case.free_cell_m, bed_diameter_m=case.bed_diameter_m,
        freeboard_diameter_m=case.free_diameter_m, epsilon_mf=case.epsilon_mf,
        particle_emissivity=emissivity, hover_closure=case.solver.hover_closure,
        serio_frequency_per_s=frequency, serio_activation_j_mol=activation,
        particle_tolerance_w=case.solver.particle_tolerance_w, particle_max_iterations=case.solver.particle_max_iterations,
        coupling_tolerance_k=case.solver.coupling_tolerance_k, coupling_max_iterations=case.solver.coupling_max_iterations,
    )
def _approved(resolver, choice):
    try:
        return resolver(choice)
    except ValueError as exc:
        raise PrecalcError(str(exc)) from exc
def _formation(case) -> float:
    physical = case.physical
    moisture = float(physical.moisture_mass_fraction)
    return fuel_formation_enthalpy_j_per_kg(
        physical.heating_value, carbon=physical.carbon_mass_fraction_wf * (1.0 - moisture),
        hydrogen=physical.hydrogen_mass_fraction_wf * (1.0 - moisture),
        nitrogen=physical.nitrogen_mass_fraction_wf * (1.0 - moisture),
        sulfur=physical.sulfur_mass_fraction_wf * (1.0 - moisture), moisture=moisture,
    )
def _pyro_flows(released, n_bed, n_free):
    bed = np.zeros((n_bed, len(GAS_SPECIES)))
    free = np.zeros((n_free, len(GAS_SPECIES)))
    index = {name: position for position, name in enumerate(GAS_SPECIES)}
    for name in PYROLYSIS_GAS_SPECIES:
        if name in index:
            free[:, index[name]] += released.gas_mol_s[name][:n_free]
            bed[:, index[name]] += released.gas_mol_s[name][n_free:]
    water = released.water_kg_s / molar_mass_kg_mol("H2O")
    free[:, index["H2O"]] += water[:n_free]
    bed[:, index["H2O"]] += water[n_free:]
    free[:, index["C6H6"]] += released.benzene_mol_s[:n_free]
    bed[:, index["C6H6"]] += released.benzene_mol_s[n_free:]
    free[:, index["C10H8"]] += released.naphthalene_mol_s[:n_free]
    bed[:, index["C10H8"]] += released.naphthalene_mol_s[n_free:]
    return bed, free
def _bubble_volume(case, profile):
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2 * case.bed_cell_m
    return profile.epsilon_b * area
def _suspension_volume(case, profile):
    """悬浮相气体体积 ε_d·(1−ε_b)·V。浓度按气体体积使用。解释性选择，原文没有把这个体积写进冻结量。
    Source: Hamel Eq.3.10, S. 24
    """
    area = np.pi * (0.5 * case.bed_diameter_m) ** 2 * case.bed_cell_m
    return profile.epsilon_d * (1.0 - profile.epsilon_b) * area
def _free_volume(case):
    volumes = np.empty(case.n_free)
    for cell in range(case.n_free):
        geometric = case.n_free - 1 - cell
        z0 = geometric * case.free_cell_m
        d0 = diameter_at(z0, case.bed_diameter_m, case.free_diameter_m)
        d1 = diameter_at(z0 + case.free_cell_m, case.bed_diameter_m, case.free_diameter_m)
        if d0 == d1:
            volumes[cell] = np.pi * (0.5 * d0) ** 2 * case.free_cell_m
        else:
            volumes[cell] = section_volume_m3(z0, z0 + case.free_cell_m, case.bed_diameter_m, case.free_diameter_m)
    return volumes
def _primary_vector(case):
    physical = case.physical
    vector = np.zeros(len(GAS_SPECIES))
    index = {name: position for position, name in enumerate(GAS_SPECIES)}
    air = float(physical.primary_air_amount_flow)
    vector[index["O2"]] = 0.21 * air
    vector[index["N2"]] = 0.79 * air + float(physical.nitrogen_purge_volume_flow)
    vector[index["H2O"]] = float(physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    return vector
def _primary_temperature(case) -> float:
    """一个温度，使一次风、蒸汽与吹扫氮的混合焓等于各自温度下的焓。ADR 0014 C2。
    二分区间 200–2500 K、60 次。求解器参数：ADR 0040。
    """
    physical = case.physical
    inlet = float(physical.inlet_temperature)
    air = float(physical.primary_air_amount_flow)
    purge = float(physical.nitrogen_purge_volume_flow)
    steam = float(physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    target = (
        _gas_enthalpy({"O2": 0.21 * air, "N2": 0.79 * air}, inlet)
        + _gas_enthalpy({"H2O": steam}, inlet)
        + _gas_enthalpy({"N2": purge}, case.purge_temperature_k)
    )
    flows = {name: float(amount) for name, amount in zip(GAS_SPECIES, _primary_vector(case)) if amount}
    lo, hi = 200.0, 2500.0
    if (_gas_enthalpy(flows, lo) - target) * (_gas_enthalpy(flows, hi) - target) > 0.0:
        raise PrecalcError("一次风等效温度不在 200–2500 K")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _gas_enthalpy(flows, mid) > target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)
def _secondary_vector(case):
    vector = np.zeros(len(GAS_SPECIES))
    index = {name: position for position, name in enumerate(GAS_SPECIES)}
    air = float(case.physical.secondary_air_amount_flow)
    vector[index["O2"]] = 0.21 * air
    vector[index["N2"]] = 0.79 * air
    return vector
def _bed_slip(flow, profile, group):
    """床区 u_r = u_d/ε_d。解释性选择：ADR 0035。
    Source: Hamel Eq.3.10, S. 24
    Source: Hamel Eq.3.11, S. 25
    """
    out = np.empty(flow.u0.size)
    for index, speed in enumerate(flow.u0):
        dense = suspension_superficial_velocity(profile.u_mf, float(speed), DISTRIBUTOR, group, True)
        if not 0.0 < profile.epsilon_d[index] < 1.0:
            raise PrecalcError("ε_d 不在 (0, 1)，不能求 u_r")
        out[index] = dense / float(profile.epsilon_d[index])
    return out
def _slip(case, flow, profile, books, rho_p, group, gas, hover, rho_free, mu_free):
    """自由状态取轨迹加权滑移，悬停取 v_t。解释性选择：ADR 0035。"""
    bed = _bed_slip(flow, profile, group)
    free_u = np.zeros((case.n_free, len(books)))
    hover_u = np.zeros_like(free_u)
    ratios = []
    for cell in range(case.n_bed):
        for klass, diameter in enumerate(case.diameters_m):
            terminal = terminal_settling_velocity(gas.rho, gas.mu, float(diameter), rho_p, case.phi_s)
            ratios.append(("bed", cell, klass, float(bed[cell]), terminal, float(bed[cell]) / terminal))
    for klass, book in enumerate(books):
        for cell, item in enumerate(book.cells):
            balance = case.n_free - 1 - cell
            terminal = terminal_settling_velocity(float(rho_free[cell]), float(mu_free[cell]), float(case.diameters_m[klass]), rho_p, case.phi_s)
            if item.mean_slip_m_s is not None:
                free_u[balance, klass] = item.mean_slip_m_s
                ratios.append(("free", balance, klass, float(item.mean_slip_m_s), terminal, float(item.mean_slip_m_s) / terminal))
            if hover[balance, klass]:
                hover_u[balance, klass] = terminal
                ratios.append(("hover", balance, klass, terminal, terminal, 1.0))
    return bed, free_u, hover_u, tuple(ratios)
def _adaptive_bed(case, flow, gas, d_bar, group="B"):
    """求解器参数：ADR 0036。整床同一 n，加倍后重求 ρ_P 与 profile，采用较细的一组。"""
    n = int(case.solver.bed_integration_intervals_per_cell)
    limit = int(case.solver.bed_integration_max_intervals_per_cell)
    def once(intervals):
        rho = _density(case, flow, gas, d_bar, intervals)
        profile, conditions = _profile(case, flow, gas, d_bar, rho, group, intervals)
        return rho, profile, conditions
    def rel(left, right):
        den = max(abs(float(left)), abs(float(right)))
        return 0.0 if den == 0.0 else abs(float(left) - float(right)) / den
    rho, profile, conditions = once(n)
    history = []
    while n * 2 <= limit:
        rho_n, profile_n, conditions_n = once(n * 2)
        masses = tuple(rel(a, b) for a, b in zip(profile.solids_mass, profile_n.solids_mass))
        rho_rel = rel(rho, rho_n)
        history.append((n, n * 2, rho_rel, masses))
        rho, profile, conditions, n = rho_n, profile_n, conditions_n, n * 2
        if rho_rel < 1.0e-6 and max(masses) < 1.0e-6:
            return rho, profile, conditions, n, tuple(history)
    raise PrecalcError(f"床层子步到上限 {limit} 仍未过门：{tuple(history)}")
def _seed_freeboard(case, state, masses, temperature, hover, absent, notes, iteration, consistency) -> None:
    """每个存在的状态焦炭取 1e-6·M_tot，灰按 w_A,0/(1−w_A,0)。ADR 0017 第 9 项。"""
    yields = _yield_at(case, temperature, None, None, notes, iteration, consistency)
    ash = _ash_fraction(case, yields.x_char)
    seed = SEED_CHAR * float(np.sum(masses))
    ash_seed = seed * ash / (1.0 - ash)
    for cell in range(hover.shape[0]):
        for klass in range(hover.shape[1]):
            if not absent[cell, klass]:
                state.freeboard_char_free_kg[cell, klass] = seed
                state.freeboard_ash_free_kg[cell, klass] = ash_seed
            if hover[cell, klass]:
                state.freeboard_char_hover_kg[cell, klass] = seed
                state.freeboard_ash_hover_kg[cell, klass] = ash_seed
def _gas_release(released) -> float:
    gas = sum(float(np.sum(flow) * molar_mass_kg_mol(name)) for name, flow in released.gas_mol_s.items())
    tar = float(np.sum(released.benzene_mol_s) * molar_mass_kg_mol("C6H6") + np.sum(released.naphthalene_mol_s) * molar_mass_kg_mol("C10H8"))
    return float(np.sum(released.water_kg_s) + gas + tar)
def _hover_pairs(hover) -> tuple[tuple[int, int], ...]:
    pairs = np.argwhere(hover)
    return tuple((int(cell), int(klass)) for cell, klass in pairs)
def _areas(bed_t, free_t, bed_h, free_h, bed_d, free_d) -> float:
    bed = np.asarray(bed_t, dtype=float)
    free = np.asarray(free_t, dtype=float)
    total = 0.0
    for index, temperature in enumerate(bed):
        total += wall_area_m2("bed", index, bed.size, bed_h, bed_d, free_d) * (float(temperature) - AMBIENT_K)
    for index, temperature in enumerate(free):
        total += wall_area_m2("freeboard", index, free.size, free_h, bed_d, free_d) * (float(temperature) - AMBIENT_K)
    return total
def _column(name, value) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or not np.isfinite(array).all():
        raise PrecalcError(f"{name} 必须是有限的一维数组")
    return array
def _vector(flows: dict[str, float]) -> np.ndarray:
    return np.array([flows[name] for name in GAS_SPECIES])
