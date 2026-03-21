"""Cell 类：单个轴向离散单元，物理一致性重构版。

Source: docs/CLAUDE.md Phase 5.1; Hamel & Krumm (2001)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict
import numpy as np
import numpy.typing as npt
from src.core.constants import Rg, T_REF, g
from src.core.species import (
    GAS_SPECIES, GAS_SPECIES_INDEX, MOLECULAR_WEIGHT, N_GAS, TarFuelType,
    cp_molar, enthalpy_molar, gas_density_ideal, gas_diffusivity_correlation,
    gas_viscosity_power_law, cp_char, cp_ash, cp_sand, TAR_SURROGATE_FORMULA,
    get_tar_component_mapping, formation_enthalpy_dry_fuel, get_atom_count
)
from src.kinetics.arrhenius import k_hobbs, k_standard
from src.kinetics.char_reactions import d_core_from_spm_char_conversion, rate_R1, rate_R2, rate_R3, rate_R4_effective
from src.kinetics.gas_reactions import rate_R5_bubble, rate_R5_suspension, rate_R6, rate_R7, rate_R8, rate_R9
from src.kinetics.tar_reactions import TarReactionId, get_lumped_tar_stoichiometry, rate_R10, rate_R11_bubble, rate_R11_suspension
from src.physics.bubble_dynamics import bubble_rise_velocity, mori_wen_bubble_diameter
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.physics.minimum_fluidization import compute_u_mf
from src.physics.phase_fractions import calc_epsilon_b, calc_epsilon_d

S_CHAR, S_VM, S_MOISTURE, S_ASH, N_SOLID_COMP = 0, 1, 2, 3, 4

@dataclass
class CellGeometry:
    D_bed: float = 0.6; dh: float = 0.5; h_center: float = 0.0

@dataclass
class SolidProps:
    rho_s: float = 1400.0; d_p: float = 0.001; phi_s: float = 0.86; eps_mf: float = 0.45; n_size_classes: int = 1
    d_p_classes: npt.NDArray[np.float64] = field(default_factory=lambda: np.array([0.001]))
    mass_fractions: npt.NDArray[np.float64] = field(default_factory=lambda: np.array([1.0]))
    sulfur_fraction: float = 0.0; sulfur_volatile_frac: float = 0.5; nitrogen_fraction: float = 0.0
    catalyst_solid_fraction: float = 1.0; moisture_wt: float = 16.9; ash_dry_wt: float = 11.41; VM_daf: float = 53.42
    C_dry: float = 61.5; H_dry: float = 4.1; O_dry: float = 21.8; HHV_dry_MJ_kg: float = 20.0
    pyrolysis_tar_carbon_frac: float = 0.2; h_f_dry: float = field(init=False)
    def __post_init__(self) -> None:
        self.h_f_dry = formation_enthalpy_dry_fuel(self.C_dry, self.H_dry, self.O_dry, getattr(self, "S_dry", 0.0), self.HHV_dry_MJ_kg)

class Cell:
    def __init__(self, geo: CellGeometry | None = None, solid: SolidProps | None = None, fuel_type: TarFuelType = "coal") -> None:
        self.geo = geo or CellGeometry(); self.solid = solid or SolidProps(); self.fuel_type = fuel_type
        self.use_gibbs_minor = False; self.u0_target = None; self.heat_loss_frac = 0.0
        self.T, self.P = 1200.0, 2_500_000.0
        self.N_b, self.N_d, self.N_b_in, self.N_d_in = np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS)
        self.N_zu_b, self.N_zu_d, self.N_rez_b, self.N_rez_d = np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS)
        nk = self.solid.n_size_classes
        self.m_solid, self.m_solid_in, self.m_solid_zu = np.zeros((nk, N_SOLID_COMP)), np.zeros((nk, N_SOLID_COMP)), np.zeros((nk, N_SOLID_COMP))
        self.u_mf = self.u_b = self.d_b = self.eps_b = self.eps_d = self.K_bd = self.V_b = self.V_d = self.u0 = 0.0
        self.R_gas_b, self.R_gas_d = np.zeros(N_GAS), np.zeros(N_GAS); self.R_solid = np.zeros((nk, N_SOLID_COMP))
        self._vm_gas_source_cache = np.zeros(N_GAS); self._vm_solid_sink_cache = np.zeros((nk, N_SOLID_COMP)); self._vm_cache_valid = False
        self.N_ex = np.zeros(N_GAS); self._work_res = np.zeros(2 * N_GAS + nk * N_SOLID_COMP + 1)
        self.T_zu_gas = self.T_rez_gas = self.T_zu_solid = self.T_in_gas = self.T_in_solid = 293.15

    def _mole_fractions(self, phase: str) -> npt.NDArray[np.float64]:
        N = self.N_b if phase == "b" else self.N_d
        work = np.maximum(N, 0.0); total = np.sum(work)
        if total < 1e-12: work[GAS_SPECIES_INDEX["N2"]] = 1.0; return work
        return work / total

    def _concentrations(self, phase: str) -> npt.NDArray[np.float64]:
        return self._mole_fractions(phase) * (self.P / (Rg * self.T))

    def calc_hydrodynamics(self) -> None:
        """计算流体力学状态（eps_b, eps_d, u_mf, d_b, u_b, K_bd）。
        
        Ref: Hamel (1999) Eq. 3.44, 3.50, 3.52
        """
        A_bed = np.pi / 4.0 * self.geo.D_bed**2; T, P = self.T, self.P
        rho_g = gas_density_ideal(P, T, {sp: float(y) for sp, y in zip(GAS_SPECIES, self._mole_fractions("d"))})
        mu_g = gas_viscosity_power_law(T, 1.8e-5); D_g = gas_diffusivity_correlation(T, P)
        self.u_mf = compute_u_mf(rho_g, self.solid.rho_s, self.solid.d_p, mu_g, self.solid.eps_mf, self.solid.phi_s)
        self.u0 = ( (np.sum(np.maximum(self.N_b, 0.0)) + np.sum(np.maximum(self.N_d, 0.0))) * Rg * T / P ) / A_bed if self.u0_target is None else self.u0_target
        self.d_b = mori_wen_bubble_diameter(self.geo.h_center, self.u0, self.u_mf, self.geo.D_bed)
        self.u_b = bubble_rise_velocity(self.u0, self.u_mf, self.d_b)
        self.eps_b = np.clip(calc_epsilon_b(self.u0, self.u_mf, self.u_b), 0.01, 0.7) # 限制气泡份额
        self.eps_d = 1.0 - self.eps_b; V_cell = A_bed * self.geo.dh; self.V_b, self.V_d = self.eps_b * V_cell, self.eps_d * V_cell
        u_br = calc_u_br(self.u_mf / max(self.solid.eps_mf, 0.01), P)
        self.K_bd = calc_kbd(u_br, self.d_b, D_g, self.solid.eps_mf, self.u_b)

    def calc_exchange(self) -> None:
        self.N_ex = self.K_bd * self.V_b * (self._concentrations("b") - self._concentrations("d"))

    def calc_reactions(self) -> None:
        idx = GAS_SPECIES_INDEX; T, P = self.T, self.P; C_b, C_d = self._concentrations("b"), self._concentrations("d")
        y_b, y_d = self._mole_fractions("b"), self._mole_fractions("d"); D_g = gas_diffusivity_correlation(T, P)
        self.R_solid.fill(0.0); self.R_gas_b.fill(0.0); self.R_gas_d.fill(0.0)
        
        # 1. Bubble Phase Reactions [mol/s]
        r5b = rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], P, y_b, idx) * self.V_b
        r6b = rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]]) * self.V_b
        self.R_gas_b[idx["CO"]] -= 2*r5b - r6b; self.R_gas_b[idx["O2"]] -= r5b + 1.5*r6b
        self.R_gas_b[idx["CO2"]] += 2*r5b; self.R_gas_b[idx["H2O"]] += 2*r6b; self.R_gas_b[idx["CH4"]] -= r6b

        # 2. Suspension Phase Reactions [mol/s]
        gas_src_vm = self._vm_gas_source_cache if self._vm_cache_valid else self._calc_drying_pyrolysis_gas_source(self.geo.dh/max(self.u_mf, 1e-3))
        self.R_gas_d += gas_src_vm
        
        # --- 动力学速率计算 ---
        r5d = rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], P, y_d, idx) * self.V_d
        r6d = rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]]) * self.V_d
        
        areas = self._calc_char_surface_area_per_class(); total_area = np.sum(areas)
        r1 = r2 = 0.0; alpha = 1.0
        if total_area > 0:
            X = self._compute_char_conversion(); dc = d_core_from_spm_char_conversion(X, self.solid.d_p)
            r1, alpha = rate_R1(T, C_d[idx["O2"]], self.solid.d_p, D_g, dc, self.fuel_type)
            # alpha = max(alpha, 0.7) # 移除强制下限
            r2 = rate_R2(T, C_d[idx["H2O"]], self.solid.d_p, D_g, dc)

        # --- 全局氧气限速机制 (O2 Supply Guard) ---
        o2_supply = max(self.N_zu_d[idx["O2"]] + self.N_d_in[idx["O2"]] + max(self.N_ex[idx["O2"]], 0.0), 1e-6)
        total_o2_demand = r5d + 1.5*r6d + alpha * r1 * total_area
        limit_factor = min(1.0, (0.95 * o2_supply) / max(total_o2_demand, 1e-9))
        
        r5d_lim, r6d_lim, r1_lim = r5d * limit_factor, r6d * limit_factor, r1 * limit_factor
        
        # --- 累加悬浮相源项 ---
        self.R_gas_d[idx["CO"]] -= 2*r5d_lim - r6d_lim; self.R_gas_d[idx["O2"]] -= r5d_lim + 1.5*r6d_lim
        self.R_gas_d[idx["CO2"]] += 2*r5d_lim; self.R_gas_d[idx["H2O"]] += 2*r6d_lim; self.R_gas_d[idx["CH4"]] -= r6d_lim
        
        if total_area > 0:
            self.R_gas_d[idx["O2"]] -= alpha * r1_lim * total_area
            self.R_gas_d[idx["CO"]] += (2*(1-alpha)*r1_lim + r2) * total_area # R2 暂不限速
            self.R_gas_d[idx["CO2"]] += (2*alpha-1)*r1_lim * total_area
            self.R_gas_d[idx["H2O"]] -= r2 * total_area; self.R_gas_d[idx["H2"]] += r2 * total_area
            self.R_solid[:, S_CHAR] -= (r1_lim + r2) * 12.011e-3 * areas

    def _calc_char_surface_area_per_class(self) -> npt.NDArray[np.float64]:
        # 数值防护：V_d 不应过小
        V_cell = (np.pi/4.0 * self.geo.D_bed**2) * self.geo.dh
        V_d_safe = max(self.V_d, 0.05 * V_cell)
        M_inv = self.solid.rho_s * (1.0 - self.solid.eps_mf) * V_d_safe
        
        m_out = np.sum(np.maximum(self.m_solid, 0.0))
        if m_out <= 1e-12: return np.zeros(self.solid.n_size_classes)
        
        # 限制 char_frac 以稳定初值演化
        char_frac = np.clip(np.maximum(self.m_solid[:, S_CHAR], 0.0) / m_out, 0.0, 1.0)
        char_mass = M_inv * char_frac
        return char_mass * 6.0 / (self.solid.rho_s * self.solid.d_p_classes)

    def _compute_char_conversion(self) -> float:
        m_in = np.sum(np.maximum(self.m_solid_zu[:, S_CHAR] + self.m_solid_in[:, S_CHAR], 0.0))
        return float(np.clip(1.0 - np.sum(np.maximum(self.m_solid[:, S_CHAR], 0.0))/max(m_in, 1e-12), 0.0, 1.0))

    def calc_gas_balance(self) -> npt.NDArray[np.float64]:
        """气相摩尔守恒残差向量。
        
        Ref: Hamel (1999) Eq. 2.1, 2.2
        """
        return np.concatenate([
            self.N_zu_d + self.N_rez_d + self.N_d_in + self.R_gas_d - self.N_d + self.N_ex,
            self.N_zu_b + self.N_rez_b + self.N_b_in + self.R_gas_b - self.N_b - self.N_ex
        ])

    def calc_solid_balance(self) -> npt.NDArray[np.float64]:
        return self.m_solid_zu + self.m_solid_in + self.R_solid - self.m_solid

    def calc_energy_balance(self) -> float:
        """全床能量平衡残差 [W]。
        
        Ref: Hamel (1999) Eq. 2.7
        """
        H_in = (self._calc_gas_enthalpy_flow(self.N_b_in + self.N_d_in, self.T_in_gas) + 
                self._calc_gas_enthalpy_flow(self.N_zu_b + self.N_zu_d, self.T_zu_gas) +
                self._calc_gas_enthalpy_flow(self.N_rez_b + self.N_rez_d, self.T_rez_gas) +
                self._calc_solid_enthalpy_flow(self.m_solid_in, self.T_in_solid) +
                self._calc_solid_enthalpy_flow(self.m_solid_zu, self.T_zu_solid))
        H_out = self._calc_gas_enthalpy_flow(self.N_b + self.N_d, self.T) + self._calc_solid_enthalpy_flow(self.m_solid, self.T)
        return H_in * (1.0 - self.heat_loss_frac) - H_out

    def _calc_gas_enthalpy_flow(self, N: npt.NDArray[np.float64], T: float) -> float:
        return float(np.dot(np.nan_to_num(N), np.array([enthalpy_molar(sp, T) for sp in GAS_SPECIES])))

    def _calc_solid_enthalpy_flow(self, m: npt.NDArray[np.float64], T: float) -> float:
        w_ash = self.solid.ash_dry_wt / 100.0; w_vm = (self.solid.VM_daf/100.0)*(1-w_ash)
        hf = np.array([0.0, self.solid.h_f_dry/max(w_vm, 1e-9), -15.866e6, 0.0])
        cp = 0.3*cp_char(T) + 0.3*cp_ash(T) + 0.4*cp_sand(T)
        return float(np.dot(np.sum(np.nan_to_num(m), axis=0), hf) + np.sum(m)*cp*(T-T_REF))

    def residuals(self) -> npt.NDArray[np.float64]:
        try:
            self.calc_hydrodynamics(); self.calc_exchange(); self.calc_reactions()
            rg = np.sum(self.N_zu_b + self.N_zu_d + self.N_b_in + self.N_d_in) + 1.0
            rs = np.sum(self.m_solid_zu + self.m_solid_in) + 0.1
            re = abs(self._calc_gas_enthalpy_flow(self.N_b_in+self.N_d_in, self.T_in_gas)) + 1e6
            res = self.calc_gas_balance(); self._work_res[:2*N_GAS] = res / rg
            self._work_res[2*N_GAS : 2*N_GAS+self.solid.n_size_classes*N_SOLID_COMP] = self.calc_solid_balance().flatten() / rs
            self._work_res[-1] = self.calc_energy_balance() / re
        except: self._work_res.fill(1e5)
        return np.nan_to_num(self._work_res, nan=1e5)

    def _calc_drying_pyrolysis_gas_source(self, tau: float) -> npt.NDArray[np.float64]:
        from src.thermal.drying import solve_drying_CN; from src.thermal.devolatilization import daem_conversion_radial
        m_vm_in = np.sum(self.m_solid_zu[:, S_VM] + self.m_solid_in[:, S_VM])
        m_moist_in = np.sum(self.m_solid_zu[:, S_MOISTURE] + self.m_solid_in[:, S_MOISTURE])
        dry = solve_drying_CN(self.solid.d_p, self.T, 300.0, self.solid.moisture_wt, max(tau, 0.05), Nr=12, Nt=80, return_history=True)
        x_vm = daem_conversion_radial(dry["T_history_rt"], dry["t"], dry["r_nodes"]) * (1.0 - (min(dry["r_evap"][-1], self.solid.d_p*0.5)/(self.solid.d_p*0.5))**3)
        m_vm_rel = m_vm_in * np.clip(x_vm, 0, 1); res = np.zeros(N_GAS); res[GAS_SPECIES_INDEX["H2O"]] = (m_moist_in * dry["X_dry"][-1]) / 0.018015
        to_daf = 1.0/max(1.0-self.solid.ash_dry_wt/100.0, 1e-9)
        nC = m_vm_rel * (self.solid.C_dry/100.0)*to_daf/0.012011; nH = m_vm_rel * (self.solid.H_dry/100.0)*to_daf/0.001008; nO = m_vm_rel * (self.solid.O_dry/100.0)*to_daf/0.016
        prod = self._allocate_pyrolysis_products_elemental(nC, nH, nO, m_vm_rel*(self.solid.nitrogen_fraction/100.0)*to_daf/0.014)
        for sp, v in prod.items(): res[GAS_SPECIES_INDEX[sp]] += v
        self.R_solid[:, S_MOISTURE] -= (m_moist_in * dry["X_dry"][-1]); self.R_solid[:, S_VM] -= m_vm_rel
        return res

    def _allocate_pyrolysis_products_elemental(self, nC, nH, nO, nN=0.0) -> Dict[str, float]:
        out = {sp: 0.0 for sp in GAS_SPECIES}; out["NH3"] = nN; nH = max(nH - 3*nN, 0.0)
        tar_surr = get_tar_component_mapping(self.fuel_type)["TAR1"]; ct, ht = TAR_SURROGATE_FORMULA[tar_surr]; ft = np.clip(self.solid.pyrolysis_tar_carbon_frac, 0, 0.9)
        nt = (ft*nC)/ct; nco = min(nO, max(nC-ct*nt, 0.0)); nch4 = max(nC-ct*nt-nco, 0.0); hu = ht*nt + 4*nch4
        if hu > nH: defic = hu-nH; dch4 = min(nch4, defic/4.0); nch4 -= dch4; nt -= (defic-4*dch4)/ht
        out.update({"CO": nco, "CH4": nch4, "TAR1": nt, "H2": max(nH-ht*nt-4*nch4, 0)/2.0})
        return out

    def compute_vorabrechnung(self, tau: float) -> None:
        self.R_solid.fill(0.0); self._vm_gas_source_cache = self._calc_drying_pyrolysis_gas_source(tau)
        self._vm_solid_sink_cache = self.R_solid.copy(); self._vm_cache_valid = True
