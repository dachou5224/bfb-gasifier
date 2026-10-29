"""完整残差 F(x; 冻结量)。反应、交换和焓在试探点上计算，不缓存。"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from hamel.balances import energy_residual, gas_mole_residual, solid_mass_flow_residual
from hamel.kinetics import (
    mechanism_factor, particle_temperature, r1_combustion_rate, r10_tar_oxidation_rate,
    r11_corella_rate, r11_serio_rate, r2_steam_gasification_rate, r3_boudouard_rate,
    r4_hydrogasification_rate, r5_co_oxidation_rate, r6_hydrogen_oxidation_rate,
    r7_methane_oxidation_rate, r8_water_gas_shift_rate, r9_methane_reforming_rate,
    reaction_stoichiometry,
)
from hamel.properties import (
    GAS_SPECIES, ash_specific_heat_j_per_kg_k, char_specific_enthalpy_j_per_kg,
    chemical_potential_j_mol, density_kg_m3, enthalpy_j_mol, gas_constant_j_mol_k,
    molar_mass_kg_mol, stefan_boltzmann_w_m2_k4,
)
from hamel.properties.data import atomic_weight, atoms
from hamel.properties.transport import mixture_thermal_conductivity_w_m_k, mixture_viscosity_pa_s
from hamel.pyrolysis import heat_transfer_coefficient
from hamel.solver.frozen import AssemblyError, FreeboardDomainError, FrozenQuantities, TrialNotEvaluable
from hamel.solver.sensitivity import diameter_at, redirect_ledger, section_volume_m3, section_wall_m2
from hamel.solver.vector import ReactorState, UnknownSlot, unpack_state, unknown_slots
# 外侧温度：ADR 0029。ε_P 在冻结量上，默认 0.85。ADR 0017 第 17 项；ADR 0039。
_AMBIENT_K = 298.15
_INDEX = {name: index for index, name in enumerate(GAS_SPECIES)}
_TARS = ("C6H6", "C10H8")

@dataclass(frozen=True)
class AssemblyResidual:
    """values 与变量名表对齐。各 gap 是全炉加总后的闭合差，装配正确时应为 0。"""
    values: np.ndarray
    slots: tuple[UnknownSlot, ...]
    m_inert_bed_kg: np.ndarray
    char_consumption_bed_kg_s: np.ndarray
    char_migration_bed_kg_s: np.ndarray
    char_consumption_free_kg_s: np.ndarray
    char_consumption_hover_kg_s: np.ndarray
    char_migration_free_kg_s: np.ndarray
    char_migration_hover_kg_s: np.ndarray
    gas_gap_mol_s: np.ndarray
    char_gap_kg_s: float
    ash_gap_kg_s: float
    energy_gap_w: float
    element_gap_kg_s: dict[str, float]
    redirect_char_kg_s: float = 0.0
    redirect_ash_kg_s: float = 0.0
    redirect_fuel_kg_s: float = 0.0
    redirect_enthalpy_w: float = 0.0
    redirect_char_by_class_kg_s: np.ndarray | None = None
    redirect_ash_by_class_kg_s: np.ndarray | None = None
    redirect_fuel_by_class_kg_s: np.ndarray | None = None

def wall_area_m2(zone: str, balance_index: int, count: int, height_m: float, bed_diameter_m: float, freeboard_diameter_m: float) -> float:
    """侧壁面积。锥段是自由空间最下格。章程 v0.11 §1 例外；ADR 0014 X1；ADR 0029。
    Source: Hamel Eq.2.7, S. 19
    """
    if zone == "bed":
        return math.pi * bed_diameter_m * height_m
    geometric = count - 1 - balance_index
    z0 = geometric * height_m
    d0 = diameter_at(z0, bed_diameter_m, freeboard_diameter_m)
    d1 = diameter_at(z0 + height_m, bed_diameter_m, freeboard_diameter_m)
    if d0 == d1:
        return math.pi * d0 * height_m
    return section_wall_m2(z0, z0 + height_m, bed_diameter_m, freeboard_diameter_m)

def evaluate_residual(unknown: np.ndarray, frozen: FrozenQuantities) -> AssemblyResidual:
    """返回式 (2.4)–(2.7)。调用方先解包。自由空间出口即体系出口。
    模型边界选择：ADR 0014 X2、X2a。模型边界选择：ADR 0026。
    解释性选择：ADR 0031、ADR 0032、ADR 0044。求解手段：ADR 0022。
    焦炭的消耗、粒级迁移、灰随迁移与灰随燃尽释放都按单位焦炭质量的消耗系数 α_s 乘存量写，
    对存量线性，不再出现 m_A/m_C，也没有 m_C = 0 的切换（ADR 0044 第 3、4 项）。
    Source: Hamel Eq.2.4, S. 16
    Source: Hamel Eq.2.5, S. 16
    Source: Hamel Eq.2.6, S. 17
    Source: Hamel Eq.2.7, S. 19
    """
    state = unpack_state(unknown, frozen)
    _require_temperature(state)
    molar = molar_mass_kg_mol("C_gr")
    masses = _particle_masses(frozen)
    m_inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
    context = _Context(frozen, state, molar, masses, m_inert)
    context.heterogeneous()
    context.domain()
    context.homogeneous()
    gas_bed, gas_free, gas_gap = context.gases()
    char_bed, char_free, char_hover, char_gap = context.chars()
    ash_free, ash_hover, ash_gap = context.ashes()
    energy_bed, energy_free, energy_gap = context.energies()
    slots = unknown_slots(frozen)
    values = np.empty(len(slots))
    for slot in slots:
        values[slot.index] = _slot_residual(slot, gas_bed, gas_free, char_bed, char_free, char_hover, ash_free, ash_hover, energy_bed, energy_free)
    return AssemblyResidual(
        values, slots, m_inert,
        context.consumption_bed, context.migration_bed,
        context.consumption_free, context.consumption_hover,
        context.migration_free, context.migration_hover,
        gas_gap, char_gap, ash_gap, energy_gap, context.element_gaps(),
        context.redirect_char, context.redirect_ash, context.redirect_fuel, context.redirect_enthalpy,
        context._redirect_mass("char"), context._redirect_mass("ash"), context._redirect_mass("fuel"),
    )

class _Context:
    def __init__(self, frozen, state, molar, masses, m_inert):
        self.frozen, self.state, self.molar = frozen, state, molar
        self.masses, self.m_inert = masses, m_inert
        n_bed, n_free, n_class = frozen.n_bed, frozen.n_freeboard, frozen.n_class
        n_species = len(GAS_SPECIES)
        self.reaction_bed = np.zeros((n_bed, 2, n_species))
        self.reaction_free = np.zeros((n_free, n_species))
        self.consumption_bed = np.zeros((n_bed, n_class))
        self.consumption_free = np.zeros((n_free, n_class))
        self.consumption_hover = np.zeros((n_free, n_class))
        self.migration_bed = np.zeros((n_bed, n_class))
        self.migration_free = np.zeros((n_free, n_class))
        self.migration_hover = np.zeros((n_free, n_class))
        self.ash_release_free = np.zeros((n_free, n_class))
        self.ash_release_hover = np.zeros((n_free, n_class))
        self.ash_migration_free = np.zeros((n_free, n_class))
        self.ash_migration_hover = np.zeros((n_free, n_class))
        self.silent_hover = np.ones((n_free, n_class), dtype=bool)
        self._bed_ash_share = _bed_ash_by_class(state.bed_char_kg, m_inert)
        self.redirect_char = self.redirect_ash = self.redirect_fuel = self.redirect_enthalpy = 0.0
        self._redirect_cache = None
    def heterogeneous(self) -> None:
        """先算异相。悬停有效域在均相之前判定，避免水蒸气为零时 R8 先抛错。"""
        frozen, state = self.frozen, self.state
        for cell in range(frozen.n_bed):
            y, c, fractions, mu, rho = _phase_point(state.bed_gas_mol_s[cell, 0], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
            for klass in range(frozen.n_class):
                gas, consumption, migration, _, _, _silent = self._char_rate(
                    state.bed_char_kg[cell, klass], self._bed_ash_share[cell, klass], klass,
                    state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell], frozen.u_r_bed_m_s[cell], y, c, fractions, mu, rho,
                )
                self.reaction_bed[cell, 0] += gas
                self.consumption_bed[cell, klass] = consumption
                self.migration_bed[cell, klass] = migration
        for cell in range(frozen.n_freeboard):
            y, c, fractions, mu, rho = _phase_point(state.freeboard_gas_mol_s[cell], state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell])
            for klass in range(frozen.n_class):
                self._one_freeboard_state(cell, klass, False, y, c, fractions, mu, rho)
                self._one_freeboard_state(cell, klass, True, y, c, fractions, mu, rho)
    def homogeneous(self) -> None:
        frozen, state = self.frozen, self.state
        for cell in range(frozen.n_bed):
            y, c, fractions, _, _ = _phase_point(state.bed_gas_mol_s[cell, 0], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
            catalyst = float(state.bed_char_kg[cell].sum() + self.m_inert[cell])
            self._homogeneous(self.reaction_bed[cell, 0], "suspension", y, c, fractions, state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell], frozen.v_suspension_m3[cell], self.m_inert[cell], catalyst, serio=False, corella=True)
            y_b, c_b, fractions_b, _, _ = _phase_point(state.bed_gas_mol_s[cell, 1], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
            self._homogeneous(self.reaction_bed[cell, 1], "bubble", y_b, c_b, fractions_b, state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell], frozen.v_bubble_m3[cell], 0.0, 0.0, serio=True, corella=False)
        for cell in range(frozen.n_freeboard):
            y, c, fractions, _, _ = _phase_point(state.freeboard_gas_mol_s[cell], state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell])
            ash = float(state.freeboard_ash_free_kg[cell].sum() + state.freeboard_ash_hover_kg[cell].sum())
            self._homogeneous(self.reaction_free[cell], "freeboard", y, c, fractions, state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell], frozen.v_freeboard_m3[cell], ash, 0.0, serio=True, corella=False)
    def _one_freeboard_state(self, cell, klass, hover, y, c, fractions, mu, rho) -> None:
        frozen, state = self.frozen, self.state
        if hover and not frozen.hover[cell, klass]:
            return
        if not hover and frozen.free_absent[cell, klass]:
            return
        char = state.freeboard_char_hover_kg if hover else state.freeboard_char_free_kg
        ash_mass = state.freeboard_ash_hover_kg if hover else state.freeboard_ash_free_kg
        gas, consumed, migrated, released, ash_migrated, silent = self._char_rate(
            char[cell, klass], ash_mass[cell, klass], klass, state.freeboard_temperature_k[cell],
            frozen.freeboard_pressure_pa[cell],
            frozen.u_r_hover_m_s[cell, klass] if hover else frozen.u_r_freeboard_m_s[cell, klass],
            y, c, fractions, mu, rho,
        )
        self.reaction_free[cell] += gas
        if hover:
            self.consumption_hover[cell, klass] = consumed
            self.migration_hover[cell, klass] = migrated
            self.ash_release_hover[cell, klass] = released
            self.ash_migration_hover[cell, klass] = ash_migrated
            self.silent_hover[cell, klass] = silent
        else:
            self.consumption_free[cell, klass] = consumed
            self.migration_free[cell, klass] = migrated
            self.ash_release_free[cell, klass] = released
            self.ash_migration_free[cell, klass] = ash_migrated
    def element_gaps(self) -> dict[str, float]:
        """反应计量的元素闭合。焦炭只含碳。热解不在这里，它来自冻结的燃料转化。"""
        reaction = self.reaction_bed.sum(axis=(0, 1)) + self.reaction_free.sum(axis=0)
        consumed = self.consumption_bed.sum() + self.consumption_free.sum() + self.consumption_hover.sum()
        gaps = {}
        for element in ("C", "H", "O", "N", "S"):
            moles = 0.0
            for index, name in enumerate(GAS_SPECIES):
                moles += float(reaction[index]) * atoms(name).get(element, 0)
            if element == "C":
                moles -= consumed / self.molar
            gaps[element] = moles * atomic_weight(element) * 1.0e-3
        return gaps
    def domain(self) -> None:
        """ADR 0026 第 5 项。只报错。"""
        frozen, state = self.frozen, self.state
        for cell in range(frozen.n_freeboard):
            solid = 0.0
            for klass in range(frozen.n_class):
                solid += state.freeboard_char_free_kg[cell, klass] + state.freeboard_char_hover_kg[cell, klass]
                solid += state.freeboard_ash_free_kg[cell, klass] + state.freeboard_ash_hover_kg[cell, klass]
                if not frozen.hover[cell, klass]:
                    continue
                if self.silent_hover[cell, klass]:
                    raise FreeboardDomainError(f"自由空间格 {cell} 粒级 k={klass + 1} 的悬停状态没有反应也没有迁移")
            volume = _cell_volume_m3("freeboard", cell, frozen)
            if solid / frozen.rho_p_kg_m3 > volume * (1.0 - frozen.epsilon_mf):
                raise FreeboardDomainError(f"自由空间格 {cell} 的固体体积超过 (1−ε_mf)·V")
        self._ash_inflow_to_smallest_hover()
    def gases(self):
        frozen, state = self.frozen, self.state
        pyro_bed = np.zeros_like(self.reaction_bed)
        pyro_bed[:, 0, :] = frozen.pyro_bed_mol_s
        feed_s = np.zeros_like(self.reaction_bed[:, 0])
        feed_b = np.zeros_like(feed_s)
        feed_s[-1] = frozen.primary_feed_mol_s * frozen.feed_to_suspension
        feed_b[-1] = frozen.primary_feed_mol_s * (1.0 - frozen.feed_to_suspension)
        exchange = self._exchange()
        zeros_bed = np.zeros_like(feed_s)
        bed = gas_mole_residual(
            feed_s, zeros_bed, self.reaction_bed[:, 0] + pyro_bed[:, 0], state.bed_gas_mol_s[:, 0],
            feed_b, zeros_bed, self.reaction_bed[:, 1], state.bed_gas_mol_s[:, 1], exchange,
        )
        feed_free = np.zeros_like(self.reaction_free)
        feed_free[-1] = state.bed_gas_mol_s[0, 0] + state.bed_gas_mol_s[0, 1]
        feed_free[frozen.secondary_balance_index] += frozen.secondary_feed_mol_s
        zeros = np.zeros_like(self.reaction_free)
        free = gas_mole_residual(feed_free, zeros, self.reaction_free + frozen.pyro_freeboard_mol_s, state.freeboard_gas_mol_s, zeros, zeros, zeros, zeros, zeros)
        gas_sum = bed.suspension_mol_s.sum(axis=0) + bed.bubble_mol_s.sum(axis=0) + free.suspension_mol_s.sum(axis=0)
        outlet = state.freeboard_gas_mol_s[0]
        gap = gas_sum - frozen.primary_feed_mol_s - frozen.secondary_feed_mol_s - frozen.pyro_bed_mol_s.sum(axis=0) - frozen.pyro_freeboard_mol_s.sum(axis=0) - self.reaction_bed.sum(axis=(0, 1)) - self.reaction_free.sum(axis=0) + outlet
        return bed, free.suspension_mol_s, gap
    def chars(self):
        frozen, state = self.frozen, self.state
        returning = (frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]) * state.freeboard_char_free_kg[-1]
        upward = frozen.k_auf_per_s[:, None] * state.bed_char_kg
        downward = frozen.m_dot_ab_kg_s[:, None] * state.bed_char_kg / frozen.geometric_mass_kg[:, None]
        leaving = np.zeros_like(state.bed_char_kg)
        leaving[-1] = frozen.drain_char_kg_s
        recycle = np.zeros_like(state.bed_char_kg)
        recycle[0] = returning
        # 模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038
        redirected = self._redirect_mass("char")
        if frozen.hover_closure == "B":
            recycle[0] = recycle[0] + redirected
        entered = np.zeros_like(self.migration_bed)
        entered[:, :-1] = self.migration_bed[:, 1:]
        bed = solid_mass_flow_residual(np.zeros_like(state.bed_char_kg), recycle, leaving, frozen.char_gen_bed_kg_s - self.consumption_bed, entered, self.migration_bed, upward, downward)
        free_residual = np.zeros_like(state.freeboard_char_free_kg)
        hover_residual = np.zeros_like(state.freeboard_char_hover_kg)
        for cell in range(frozen.n_freeboard):
            for klass in range(frozen.n_class):
                inflow_free, inflow_hover = self._char_inflow(cell, klass, returning_unused=None)
                mig_in = self.migration_free[cell, klass + 1] + self.migration_hover[cell, klass + 1] if klass + 1 < frozen.n_class else 0.0
                if not frozen.free_absent[cell, klass]:
                    outflow = (frozen.a_free_per_s[cell, klass] + frozen.a_hover_per_s[cell, klass] + frozen.b_free_per_s[cell, klass] + frozen.b_hover_per_s[cell, klass]) * state.freeboard_char_free_kg[cell, klass]
                    free_residual[cell, klass] = inflow_free + frozen.char_gen_free_kg_s[cell, klass] - outflow - self.consumption_free[cell, klass] - self.migration_free[cell, klass] + (1.0 - frozen.theta[cell, klass]) * mig_in
                if frozen.hover[cell, klass]:
                    hover_residual[cell, klass] = inflow_hover + frozen.char_gen_hover_kg_s[cell, klass] - self.consumption_hover[cell, klass] - self.migration_hover[cell, klass] + frozen.theta[cell, klass] * mig_in
        elutriation = (frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * state.freeboard_char_free_kg[0]
        if frozen.hover_closure == "C":
            elutriation = elutriation + redirected
        char_sum = bed.residual_kg_s.sum() + free_residual.sum() + hover_residual.sum()
        generated = frozen.char_gen_bed_kg_s.sum() + frozen.char_gen_free_kg_s.sum() + frozen.char_gen_hover_kg_s.sum()
        consumed = self.consumption_bed.sum() + self.consumption_free.sum() + self.consumption_hover.sum()
        gap = float(char_sum - generated + consumed + elutriation.sum() + frozen.drain_char_kg_s.sum())
        return bed.residual_kg_s, free_residual, hover_residual, gap
    def ashes(self):
        frozen, state = self.frozen, self.state
        entering_free, entering_hover, ejection = self._ash_entering()
        free_residual = np.zeros_like(entering_free)
        hover_residual = np.zeros_like(entering_hover)
        for cell in range(frozen.n_freeboard):
            for klass in range(frozen.n_class):
                if not frozen.free_absent[cell, klass]:
                    outflow = (frozen.a_free_per_s[cell, klass] + frozen.a_hover_per_s[cell, klass] + frozen.b_free_per_s[cell, klass] + frozen.b_hover_per_s[cell, klass]) * state.freeboard_ash_free_kg[cell, klass]
                    free_residual[cell, klass] = entering_free[cell, klass] - outflow - self.ash_migration_free[cell, klass] - self.ash_release_free[cell, klass]
                if frozen.hover[cell, klass]:
                    hover_residual[cell, klass] = entering_hover[cell, klass] - self.ash_migration_hover[cell, klass] - self.ash_release_hover[cell, klass]
        returning = (frozen.b_free_per_s[-1] + frozen.b_hover_per_s[-1]) * state.freeboard_ash_free_kg[-1]
        elutriation = (frozen.a_free_per_s[0] + frozen.a_hover_per_s[0]) * state.freeboard_ash_free_kg[0]
        held = entering_hover - self.ash_migration_hover - self.ash_release_hover
        redirected = np.zeros(frozen.n_class) if frozen.hover_closure == "A" else held.sum(axis=0)
        if frozen.hover_closure == "B":
            returning = returning + redirected
        elif frozen.hover_closure == "C":
            elutriation = elutriation + redirected
        bed_gap = float(frozen.ash_gen_bed_kg_s.sum() + returning.sum() - ejection.sum() - frozen.drain_ash_kg_s)
        gap = float(free_residual.sum() + hover_residual.sum() + bed_gap - frozen.ash_gen_free_kg_s.sum() - frozen.ash_gen_hover_kg_s.sum() - frozen.ash_gen_bed_kg_s.sum() + elutriation.sum() + frozen.drain_ash_kg_s)
        return free_residual, hover_residual, gap
    def energies(self):
        frozen, state = self.frozen, self.state
        bed_up, bed_down, bed_out, bed_ext = self._bed_enthalpy()
        free_up, free_down, free_ext = self._free_enthalpy()
        bed_ext[0] += free_down[-1]
        free_ext[-1] += bed_up[0]
        # 模型边界选择（敏感性）：ADR 0026 第 6 项、ADR 0038
        leaving_redirect = self._apply_redirect(bed_ext, free_ext)
        wall_bed, wall_free = self._wall("bed"), self._wall("freeboard")
        bed = energy_residual(bed_ext, np.zeros(frozen.n_bed), bed_up, bed_down, bed_out, wall_bed, np.zeros(frozen.n_bed))
        free = energy_residual(free_ext, np.zeros(frozen.n_freeboard), free_up, free_down, np.zeros(frozen.n_freeboard), wall_free, np.zeros(frozen.n_freeboard))
        # 交叉界面已互相放入 external。床底向下流、顶出流、排灰和壁面离开体系。
        feeds = (
            _gas_enthalpy(frozen.primary_feed_mol_s, frozen.primary_temperature_k)
            + _fuel_enthalpy(float(frozen.fuel_feed_kg_s.sum()), frozen.fuel_feed_temperature_k, frozen)
            + _gas_enthalpy(frozen.secondary_feed_mol_s, frozen.secondary_temperature_k)
        )
        gap = float(bed.residual_w.sum() + free.residual_w.sum() - feeds + free_up[0] + bed_down[-1] + bed_out.sum() + wall_bed.sum() + wall_free.sum() + leaving_redirect)
        return bed.residual_w, free.residual_w, gap
    def _ash_entering(self):
        """灰方程与最小粒级有效域共用的入流。含按状态的热解灰与燃尽释放。"""
        frozen, state = self.frozen, self.state
        entering_free = np.zeros_like(frozen.ash_gen_free_kg_s)
        entering_hover = np.zeros_like(frozen.ash_gen_hover_kg_s)
        for cell in range(frozen.n_freeboard):
            for klass in range(frozen.n_class):
                self._deposit(entering_free, entering_hover, cell, klass, False, float(frozen.ash_gen_free_kg_s[cell, klass]))
                self._deposit(entering_free, entering_hover, cell, klass, True, float(frozen.ash_gen_hover_kg_s[cell, klass]))
        for cell in range(frozen.n_freeboard):
            for klass in range(frozen.n_class):
                self._deposit_release(entering_free, entering_hover, cell, self.ash_release_free[cell, klass] + self.ash_release_hover[cell, klass])
                if klass == 0:
                    continue
                self._deposit_migration(entering_free, entering_hover, cell, klass - 1, self.ash_migration_free[cell, klass], self.ash_migration_hover[cell, klass])
        ejection = frozen.k_auf_per_s[0] * self._bed_ash_share[0]
        for klass in range(frozen.n_class):
            self._deposit(entering_free, entering_hover, frozen.n_freeboard - 1, klass, False, ejection[klass] * frozen.ejection_to_free[klass])
            self._deposit(entering_free, entering_hover, frozen.n_freeboard - 1, klass, True, ejection[klass] * frozen.ejection_to_hover[klass])
        for cell in range(frozen.n_freeboard - 1):
            mobile = state.freeboard_ash_free_kg[cell + 1]
            for klass in range(frozen.n_class):
                self._deposit(entering_free, entering_hover, cell, klass, False, frozen.a_free_per_s[cell + 1, klass] * mobile[klass])
                self._deposit(entering_free, entering_hover, cell, klass, True, frozen.a_hover_per_s[cell + 1, klass] * mobile[klass])
        for cell in range(1, frozen.n_freeboard):
            mobile = state.freeboard_ash_free_kg[cell - 1]
            for klass in range(frozen.n_class):
                self._deposit(entering_free, entering_hover, cell, klass, False, frozen.b_free_per_s[cell - 1, klass] * mobile[klass])
                self._deposit(entering_free, entering_hover, cell, klass, True, frozen.b_hover_per_s[cell - 1, klass] * mobile[klass])
        return entering_free, entering_hover, ejection
    def _char_rate(self, m_char, m_ash, klass, temperature, pressure, velocity, y, c, fractions, mu, rho):
        diameter = float(self.frozen.diameter_m[klass])
        particle_k = _particle_k(diameter, temperature, pressure, velocity, c, fractions, mu, rho, self.frozen)
        film = 0.5 * (temperature + particle_k)
        mu_film = mixture_viscosity_pa_s(fractions, film, pressure)
        rho_film = density_kg_m3(film, pressure, fractions)
        phi = mechanism_factor(diameter, particle_k)
        r1 = r1_combustion_rate(m_char, self.frozen.rho_p_kg_m3, diameter, diameter, particle_k, particle_k, particle_k, temperature, c[_INDEX["O2"]], pressure, velocity, mu_film, rho_film)
        # 式 (5.25)、(5.34) 的 T 不在 ADR 0022 C1 点名的式号里，调用方传入格子气温。R3 的合同同样取 cell 温度。
        r2 = r2_steam_gasification_rate(m_char, self.frozen.rho_p_kg_m3, diameter, c[_INDEX["H2O"]], temperature)
        r3 = r3_boudouard_rate(m_char, self.molar, y[_INDEX["CO2"]] * pressure, y[_INDEX["CO"]] * pressure, temperature, gas_constant_j_mol_k())
        r4 = r4_hydrogasification_rate(m_char, self.frozen.rho_p_kg_m3, diameter, c[_INDEX["H2"]], temperature)
        gas = np.zeros(len(GAS_SPECIES))
        char_mass = _apply(gas, "R1", r1.extent_mol_s, self.molar, phi=phi)
        char_mass += _apply(gas, "R2", r2.extent_mol_s, self.molar)
        char_mass += _apply(gas, "R3", r3.extent_mol_s, self.molar)
        char_mass += _apply(gas, "R4", r4.extent_mol_s, self.molar)
        consumption = -char_mass
        particle_mass = float(self.masses[klass])
        previous = 0.0 if klass == 0 else float(self.masses[klass - 1])
        # 单位焦炭质量的消耗系数 α_s（R1、R2、R4 按颗粒数，R3 按焦炭质量），各状态用本状态的输入。
        # 消耗 = α_s·M_C，迁移与灰项都乘 α_s，对存量线性。解释性选择：ADR 0044
        alpha = (r1.specific + r2.specific + r4.specific) * self.molar / particle_mass + r3.specific
        step = 0.0 if particle_mass == previous else previous / (particle_mass - previous)
        char_mig = alpha * m_char * step
        ash_mig = alpha * m_ash * step
        ash_release = alpha * m_ash
        silent = r1.specific == 0.0 and r2.specific == 0.0 and r3.specific == 0.0 and r4.specific == 0.0
        return gas, consumption, char_mig, ash_release, ash_mig, silent
    def _homogeneous(self, gas, phase, y, c, fractions, temperature, pressure, volume, ash, catalyst, serio, corella):
        r5_phase = "suspension" if phase == "suspension" else "bubble"
        water = None if r5_phase == "suspension" else c[_INDEX["H2O"]]
        _apply(gas, "R5", r5_co_oxidation_rate(r5_phase, c[_INDEX["CO"]], c[_INDEX["O2"]], temperature, volume, water).extent_mol_s, self.molar)
        _apply(gas, "R6", r6_hydrogen_oxidation_rate(c[_INDEX["H2"]], c[_INDEX["O2"]], temperature, volume).extent_mol_s, self.molar)
        _apply(gas, "R7", r7_methane_oxidation_rate(c[_INDEX["CH4"]], c[_INDEX["O2"]], temperature, volume).extent_mol_s, self.molar)
        mus = [chemical_potential_j_mol(name, temperature) for name in ("CO", "H2O", "H2", "CO2")]
        _apply(gas, "R8", r8_water_gas_shift_rate(y[_INDEX["CO"]], y[_INDEX["H2O"]], y[_INDEX["H2"]], y[_INDEX["CO2"]], mus[0], mus[1], mus[2], mus[3], temperature, pressure, gas_constant_j_mol_k(), ash).extent_mol_s, self.molar)
        _apply(gas, "R9", r9_methane_reforming_rate(y[_INDEX["CH4"]] * pressure, temperature, gas_constant_j_mol_k(), ash).extent_mol_s, self.molar)
        for tar in _TARS:
            _apply(gas, "R10", r10_tar_oxidation_rate(c[_INDEX[tar]], c[_INDEX["O2"]], temperature, pressure, volume).extent_mol_s, self.molar, tar=tar)
            if serio:
                rate = r11_serio_rate(
                    c[_INDEX[tar]], temperature, gas_constant_j_mol_k(), volume,
                    frequency_per_s=self.frozen.serio_frequency_per_s,
                    activation_energy_j_mol=self.frozen.serio_activation_j_mol,
                )
                _apply(gas, "R11", rate.extent_mol_s, self.molar, tar=tar)
            if corella:
                _apply(gas, "R11", r11_corella_rate(c[_INDEX[tar]], temperature, gas_constant_j_mol_k(), catalyst).extent_mol_s, self.molar, tar=tar)
        del fractions
    def _exchange(self):
        frozen, state = self.frozen, self.state
        exchange = np.zeros_like(state.bed_gas_mol_s[:, 0])
        for cell in range(frozen.n_bed):
            _, c_s, _, _, _ = _phase_point(state.bed_gas_mol_s[cell, 0], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
            _, c_b, _, _, _ = _phase_point(state.bed_gas_mol_s[cell, 1], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
            exchange[cell] = frozen.k_bd_per_s[cell] * frozen.v_bubble_exchange_m3[cell] * (c_s - c_b)
        return exchange
    def _char_inflow(self, cell, klass, returning_unused):
        frozen, state = self.frozen, self.state
        free_in = hover_in = 0.0
        if cell == frozen.n_freeboard - 1:
            thrown = frozen.k_auf_per_s[0] * state.bed_char_kg[0, klass]
            free_in += thrown * frozen.ejection_to_free[klass]
            hover_in += thrown * frozen.ejection_to_hover[klass]
        if cell + 1 < frozen.n_freeboard:
            below = state.freeboard_char_free_kg[cell + 1, klass]
            free_in += frozen.a_free_per_s[cell + 1, klass] * below
            hover_in += frozen.a_hover_per_s[cell + 1, klass] * below
        if cell > 0:
            above = state.freeboard_char_free_kg[cell - 1, klass]
            free_in += frozen.b_free_per_s[cell - 1, klass] * above
            hover_in += frozen.b_hover_per_s[cell - 1, klass] * above
        self._reject_missing_destination(cell, klass, free_in, hover_in)
        return free_in, hover_in
    def _reject_missing_destination(self, cell, klass, free_in, hover_in) -> None:
        frozen = self.frozen
        if frozen.free_absent[cell, klass] and free_in != 0.0:
            raise AssemblyError(f"自由空间格 {cell} 粒级 k={klass + 1} 没有自由状态，却有进入自由状态的通量")
        if frozen.hover_closure != "A":
            return
        if not frozen.hover[cell, klass] and hover_in != 0.0:
            raise AssemblyError(f"自由空间格 {cell} 粒级 k={klass + 1} 没有悬停状态，却有进入悬停状态的通量")
    def _deposit(self, entering_free, entering_hover, cell, klass, hover, amount) -> None:
        if amount == 0.0:
            return
        if hover:
            entering_hover[cell, klass] += amount
        else:
            entering_free[cell, klass] += amount
    def _deposit_release(self, entering_free, entering_hover, cell, amount) -> None:
        self._deposit(entering_free, entering_hover, cell, 0, bool(self.frozen.free_absent[cell, 0]), amount)
    def _deposit_migration(self, entering_free, entering_hover, cell, klass, from_free, from_hover) -> None:
        total = from_free + from_hover
        theta = self.frozen.theta[cell, klass]
        self._deposit(entering_free, entering_hover, cell, klass, False, (1.0 - theta) * total)
        self._deposit(entering_free, entering_hover, cell, klass, True, theta * total)
    def _ash_inflow_to_smallest_hover(self) -> None:
        """与灰方程同一份入流。热解生成与燃尽释放都算流入。"""
        frozen = self.frozen
        if not np.any(frozen.hover[:, 0]):
            return
        _, entering_hover, _ = self._ash_entering()
        for cell in range(frozen.n_freeboard):
            if frozen.hover[cell, 0] and entering_hover[cell, 0] > 0.0:
                raise FreeboardDomainError(f"最小粒级在自由空间格 {cell} 被捕获，且有灰流入悬停状态")
    def _items(self):
        if self._redirect_cache is not None:
            return self._redirect_cache
        frozen = self.frozen
        if frozen.hover_closure == "A":
            self._redirect_cache = ()
            return self._redirect_cache
        items = redirect_ledger(
            frozen, self.state, self.migration_free, self.migration_hover,
            self.ash_migration_free, self.ash_migration_hover, self._bed_ash_share,
            self.ash_release_free, self.ash_release_hover,
        )
        self._redirect_cache = items
        self.redirect_char = float(sum(item["char"] for item in items))
        self.redirect_ash = float(sum(item["ash"] for item in items))
        self.redirect_fuel = float(sum(item["fuel"] for item in items))
        self.redirect_enthalpy = float(sum(item["enthalpy"] for item in items))
        return items
    def _redirect_mass(self, key: str) -> np.ndarray:
        total = np.zeros(self.frozen.n_class)
        for item in self._items():
            total[item["klass"]] += item[key]
        return total
    def _apply_redirect(self, bed_ext, free_ext) -> float:
        items = self._items()
        if not items:
            return 0.0
        for item in items:
            zone, index = item["undo"]
            if zone == "free":
                free_ext[index] -= item["enthalpy"]
            else:
                bed_ext[index] -= item["enthalpy"]
        if self.frozen.hover_closure == "B":
            bed_ext[0] += self.redirect_enthalpy
            return 0.0
        return self.redirect_enthalpy
    def _wall(self, zone: str) -> np.ndarray:
        frozen, state = self.frozen, self.state
        count = frozen.n_bed if zone == "bed" else frozen.n_freeboard
        height = frozen.bed_height_m if zone == "bed" else frozen.freeboard_height_m
        temperature = state.bed_temperature_k if zone == "bed" else state.freeboard_temperature_k
        areas = np.array([wall_area_m2(zone, index, count, height, frozen.bed_diameter_m, frozen.freeboard_diameter_m) for index in range(count)])
        return frozen.k_w_w_m2_k * areas * (temperature - _AMBIENT_K)
    def _bed_enthalpy(self):
        frozen, state = self.frozen, self.state
        upward = np.zeros(frozen.n_bed)
        downward = np.zeros(frozen.n_bed)
        leaving = np.zeros(frozen.n_bed)
        external = np.zeros(frozen.n_bed)
        for cell in range(frozen.n_bed):
            temperature = float(state.bed_temperature_k[cell])
            gas = state.bed_gas_mol_s[cell, 0] + state.bed_gas_mol_s[cell, 1]
            char_up = float((frozen.k_auf_per_s[cell] * state.bed_char_kg[cell]).sum())
            ash_up = frozen.k_auf_per_s[cell] * self.m_inert[cell]
            fuel_up = float(frozen.fuel_up_bed_kg_s[cell].sum())
            upward[cell] = _stream_enthalpy(gas, char_up, ash_up, fuel_up, temperature, frozen)
            char_down = float((frozen.m_dot_ab_kg_s[cell] * state.bed_char_kg[cell]).sum() / frozen.geometric_mass_kg[cell])
            ash_down = frozen.m_dot_ab_kg_s[cell] * self.m_inert[cell] / frozen.geometric_mass_kg[cell]
            downward[cell] = _stream_enthalpy(np.zeros(len(GAS_SPECIES)), char_down, ash_down, float(frozen.fuel_down_bed_kg_s[cell].sum()), temperature, frozen)
            if cell == frozen.n_bed - 1:
                leaving[cell] = _stream_enthalpy(np.zeros(len(GAS_SPECIES)), float(frozen.drain_char_kg_s.sum()), frozen.drain_ash_kg_s, float(frozen.fuel_exit_bed_kg_s[cell].sum()), temperature, frozen)
                external[cell] += _gas_enthalpy(frozen.primary_feed_mol_s, frozen.primary_temperature_k)
                external[cell] += _fuel_enthalpy(float(frozen.fuel_feed_kg_s.sum()), frozen.fuel_feed_temperature_k, frozen)
            else:
                leaving[cell] = _fuel_enthalpy(float(frozen.fuel_exit_bed_kg_s[cell].sum()), temperature, frozen)
        return upward, downward, leaving, external
    def _free_enthalpy(self):
        frozen, state = self.frozen, self.state
        upward = np.zeros(frozen.n_freeboard)
        downward = np.zeros(frozen.n_freeboard)
        external = np.zeros(frozen.n_freeboard)
        for cell in range(frozen.n_freeboard):
            temperature = float(state.freeboard_temperature_k[cell])
            coefficient_up = frozen.a_free_per_s[cell] + frozen.a_hover_per_s[cell]
            coefficient_down = frozen.b_free_per_s[cell] + frozen.b_hover_per_s[cell]
            char_up = float((coefficient_up * state.freeboard_char_free_kg[cell]).sum())
            ash_up = float((coefficient_up * state.freeboard_ash_free_kg[cell]).sum())
            upward[cell] = _stream_enthalpy(state.freeboard_gas_mol_s[cell], char_up, ash_up, float(frozen.fuel_up_free_kg_s[cell].sum()), temperature, frozen)
            char_down = float((coefficient_down * state.freeboard_char_free_kg[cell]).sum())
            ash_down = float((coefficient_down * state.freeboard_ash_free_kg[cell]).sum())
            downward[cell] = _stream_enthalpy(np.zeros(len(GAS_SPECIES)), char_down, ash_down, float(frozen.fuel_down_free_kg_s[cell].sum()), temperature, frozen)
        external[frozen.secondary_balance_index] += _gas_enthalpy(frozen.secondary_feed_mol_s, frozen.secondary_temperature_k)
        return upward, downward, external

def _apply(gas, reaction, extent, molar, phi=None, tar=None) -> float:
    stoich = reaction_stoichiometry(reaction, phi_c=phi, tar=tar)
    for name, coefficient in stoich.items():
        if name == "C":
            continue
        gas[_INDEX[name]] += coefficient * extent
    return stoich.get("C", 0.0) * extent * molar

def _bed_ash_by_class(char_kg, m_inert) -> np.ndarray:
    total = char_kg.sum(axis=1)
    share = np.zeros_like(char_kg)
    for cell, char_total in enumerate(total):
        if char_total > 0.0:
            share[cell] = m_inert[cell] * char_kg[cell] / char_total
        elif m_inert[cell] != 0.0:
            raise TrialNotEvaluable(f"床区格 {cell} 没有焦炭，不能按焦炭比例分配灰")
    return share

def _phase_point(flow_mol_s, temperature_k, pressure_pa):
    total = float(np.sum(flow_mol_s))
    if total <= 0.0:
        raise TrialNotEvaluable("气相总摩尔流不是正数，浓度没有定义")
    fraction = np.asarray(flow_mol_s, dtype=float) / total
    concentration = fraction * pressure_pa / (gas_constant_j_mol_k() * temperature_k)
    named = {name: float(fraction[index]) for index, name in enumerate(GAS_SPECIES)}
    film_mu = mixture_viscosity_pa_s(named, temperature_k, pressure_pa)
    film_rho = density_kg_m3(temperature_k, pressure_pa, named)
    return fraction, concentration, named, film_mu, film_rho

def particle_temperatures(frozen, state) -> dict:
    """各 (cell, 粒级) 的颗粒温度。与残差里的式 (5.22) 是同一次对齐，只供加密检查读取。"""
    bed = np.full((frozen.n_bed, frozen.n_class), np.nan)
    free = np.full((frozen.n_freeboard, frozen.n_class), np.nan)
    hover = np.full((frozen.n_freeboard, frozen.n_class), np.nan)
    for cell in range(frozen.n_bed):
        _y, concentration, fractions, mu, rho = _phase_point(
            state.bed_gas_mol_s[cell, 0], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell],
        )
        ash = _bed_ash_by_class(state.bed_char_kg, frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1))
        for klass in range(frozen.n_class):
            if state.bed_char_kg[cell, klass] <= 0.0 and ash[cell, klass] > 0.0:
                continue
            bed[cell, klass] = _particle_k(
                float(frozen.diameter_m[klass]), float(state.bed_temperature_k[cell]), float(frozen.bed_pressure_pa[cell]),
                float(frozen.u_r_bed_m_s[cell]), concentration, fractions, mu, rho, frozen,
            )
    for cell in range(frozen.n_freeboard):
        _y, concentration, fractions, mu, rho = _phase_point(
            state.freeboard_gas_mol_s[cell], state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell],
        )
        for klass in range(frozen.n_class):
            for hover_state, target in ((False, free), (True, hover)):
                char = state.freeboard_char_hover_kg if hover_state else state.freeboard_char_free_kg
                ash_row = state.freeboard_ash_hover_kg if hover_state else state.freeboard_ash_free_kg
                if hover_state and not frozen.hover[cell, klass]:
                    continue
                if char[cell, klass] <= 0.0 and ash_row[cell, klass] > 0.0:
                    continue
                if char[cell, klass] <= 0.0 and ash_row[cell, klass] <= 0.0 and not hover_state:
                    continue
                velocity = frozen.u_r_hover_m_s if hover_state else frozen.u_r_freeboard_m_s
                target[cell, klass] = _particle_k(
                    float(frozen.diameter_m[klass]), float(state.freeboard_temperature_k[cell]),
                    float(frozen.freeboard_pressure_pa[cell]), float(velocity[cell, klass]),
                    concentration, fractions, mu, rho, frozen,
                )
    return {"bed": bed, "free": free, "hover": hover}
def _particle_k(diameter, temperature, pressure, velocity, concentration, fractions, mu, rho, frozen) -> float:
    """式 (5.22) 与式 (4.9) 对齐。d_c、T_c 与式 (5.23) 的 T 取颗粒温度。ADR 0022 C1；ADR 0017 第 20 项。"""
    del mu, rho
    guess = temperature
    for _ in range(frozen.coupling_max_iterations):
        alpha, _biot = heat_transfer_coefficient(diameter, temperature, guess, pressure, fractions, velocity, 1.0)
        film = 0.5 * (temperature + guess)
        viscosity = mixture_viscosity_pa_s(fractions, film, pressure)
        density = density_kg_m3(film, pressure, fractions)
        phi = mechanism_factor(diameter, guess)
        enthalpy = _r1_enthalpy(phi, guess)
        updated = particle_temperature(
            diameter, diameter, guess, guess, temperature, concentration[_INDEX["O2"]], pressure, velocity,
            viscosity, density, enthalpy, alpha, frozen.particle_emissivity, stefan_boltzmann_w_m2_k4(), guess,
            frozen.particle_tolerance_w, frozen.particle_max_iterations,
        )
        if abs(updated - guess) <= frozen.coupling_tolerance_k:
            return updated
        guess = updated
    raise TrialNotEvaluable("式 (5.22) 的颗粒温度没有和膜温对齐")

def _r1_enthalpy(phi: float, temperature: float) -> float:
    total = 0.0
    for name, coefficient in reaction_stoichiometry("R1", phi_c=phi).items():
        total += coefficient * enthalpy_j_mol("C_gr" if name == "C" else name, temperature)
    return total

def _particle_masses(frozen: FrozenQuantities) -> np.ndarray:
    return frozen.rho_p_kg_m3 * math.pi * frozen.diameter_m**3 / 6.0

def _cell_volume_m3(zone: str, balance: int, frozen: FrozenQuantities) -> float:
    if zone == "bed":
        return math.pi * (0.5 * frozen.bed_diameter_m) ** 2 * frozen.bed_height_m
    height = frozen.freeboard_height_m
    geometric = frozen.n_freeboard - 1 - balance
    z0 = geometric * height
    d0 = diameter_at(z0, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
    d1 = diameter_at(z0 + height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)
    if d0 == d1:
        return math.pi * (0.5 * d0) ** 2 * height
    return section_volume_m3(z0, z0 + height, frozen.bed_diameter_m, frozen.freeboard_diameter_m)

def _require_temperature(state: ReactorState) -> None:
    if np.any(state.bed_temperature_k <= 0.0) or np.any(state.freeboard_temperature_k <= 0.0):
        raise TrialNotEvaluable("温度必须为正才能计算物性")

def _gas_enthalpy(flow_mol_s, temperature_k: float) -> float:
    return float(sum(flow_mol_s[index] * enthalpy_j_mol(name, temperature_k) for index, name in enumerate(GAS_SPECIES)))

def _fuel_enthalpy(mass_kg_s: float, temperature_k: float, frozen: FrozenQuantities) -> float:
    sensible = frozen.fuel_formation_j_per_kg + frozen.fuel_heat_capacity_j_per_kg_k * (temperature_k - _AMBIENT_K)
    return mass_kg_s * sensible

def _stream_enthalpy(gas_mol_s, char_kg_s, ash_kg_s, fuel_kg_s, temperature_k, frozen) -> float:
    ash = ash_kg_s * ash_specific_heat_j_per_kg_k() * (temperature_k - _AMBIENT_K)
    char = char_kg_s * char_specific_enthalpy_j_per_kg(temperature_k)
    return _gas_enthalpy(gas_mol_s, temperature_k) + char + ash + _fuel_enthalpy(fuel_kg_s, temperature_k, frozen)

def _slot_residual(slot, gas_bed, gas_free, char_bed, char_free, char_hover, ash_free, ash_hover, energy_bed, energy_free) -> float:
    cell = slot.balance_index
    if slot.kind == "gas_outflow":
        species = _INDEX[slot.species]
        if slot.zone == "bed":
            phase = 0 if slot.phase == "suspension" else 1
            block = gas_bed.suspension_mol_s if phase == 0 else gas_bed.bubble_mol_s
            return float(block[cell, species])
        return float(gas_free[cell, species])
    if slot.kind == "temperature":
        block = energy_bed if slot.zone == "bed" else energy_free
        return float(block[cell])
    if slot.kind == "char" and slot.zone == "bed":
        return float(char_bed[cell, slot.size_class])
    if slot.kind == "char" and slot.phase == "free":
        return float(char_free[cell, slot.size_class])
    if slot.kind == "char":
        return float(char_hover[cell, slot.size_class])
    if slot.phase == "free":
        return float(ash_free[cell, slot.size_class])
    return float(ash_hover[cell, slot.size_class])
