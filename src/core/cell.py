"""Cell 类：单个轴向离散单元，物理一致性重构版。

Source: docs/CLAUDE.md Phase 5.1; Hamel & Krumm (2001)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict
import numpy as np
import numpy.typing as npt
from src.core.species import (
    GAS_SPECIES_INDEX, MOLECULAR_WEIGHT, N_GAS, TarFuelType,
    gas_diffusivity_correlation,
    gas_viscosity_power_law,
    formation_enthalpy_dry_fuel,
)
from src.core.cell_balances import (
    calc_energy_balance_residual,
    calc_gas_balance_residual,
    calc_size_migration,
    calc_solid_balance_residual,
)
from src.core.cell_thermo_helpers import (
    calc_gas_enthalpy_flow,
    gas_concentrations,
    gas_mole_fractions,
    solid_enthalpy_flow_for_fuel_props,
)
from src.core.cell_hydrodynamics import calc_cell_hydrodynamics, calc_phase_exchange
from src.core.cell_kinetics import build_reaction_sources
from src.core.cell_minor_species import (
    calc_minor_species_gibbs as calc_minor_species_gibbs_impl,
    element_moles_in_gas,
    get_element_release as get_element_release_impl,
    tracked_minor_element_moles_in_gas,
)
from src.core.cell_pyrolysis import allocate_pyrolysis_products_elemental, calc_drying_pyrolysis_sources
from src.core.cell_residual_pipeline import execute_cell_residual_pipeline
from src.core.cell_state_cache import (
    apply_frozen_vorabrechnung_hydrodynamics,
    disable_inner_nr_vorabrechnung_freeze,
    enable_inner_nr_vorabrechnung_freeze,
    get_local_thermo_bundle,
    hydrodynamics_cache_matches,
    invalidate_vorabrechnung_cache as invalidate_vorabrechnung_cache_impl,
    snapshot_vorabrechnung_hydrodynamics,
    store_hydrodynamics_cache,
    thermo_cache_matches,
)
from src.kinetics.tar_reactions import TarReactionId
from src.physics.phase_fractions import calc_epsilon_d
from src.thermodynamics.equilibrium import calc_gibbs_driving_force

S_CHAR, S_VM, S_MOISTURE, S_ASH, N_SOLID_COMP = 0, 1, 2, 3, 4

# Gibbs 微量组分松弛系数（越大越快趋向平衡，过大可能振荡）
_GIBBS_RELAX_K: float = 1.0

@dataclass
class CellGeometry:
    D_bed: float = 0.6; dh: float = 0.5; h_center: float = 0.0; N_or: int = 100

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
        # S 干基质量分数 [%]：由 sulfur_fraction（kg S / kg 干燃料）换算
        s_wt_pct = self.sulfur_fraction * 100.0
        self.h_f_dry = formation_enthalpy_dry_fuel(self.C_dry, self.H_dry, self.O_dry, s_wt_pct, self.HHV_dry_MJ_kg)
        # 确保粒径按升序排列
        if len(self.d_p_classes) > 1:
            idx = np.argsort(self.d_p_classes)
            self.d_p_classes = self.d_p_classes[idx]
            self.mass_fractions = self.mass_fractions[idx]

class Cell:
    def __init__(self, geo: CellGeometry | None = None, solid: SolidProps | None = None, fuel_type: TarFuelType = "coal") -> None:
        self.geo = geo or CellGeometry(); self.solid = solid or SolidProps(); self.fuel_type = fuel_type
        self.cell_type = "bed"
        self.solid_state_model = "legacy_stream"
        self.use_gibbs_minor = True; self.u0_target = None; self.heat_loss_frac = 0.0
        self.u_d_closure = "current"
        self.bubble_diameter_model = "mori_wen"
        self.psi_b_strategy = "wein_1992"
        self.lambda_strategy = "hamel_280"
        self.xi_strategy = "hamel_regime"
        self.bubble_velocity_strategy = "hilligardt_eq313"
        self.bubble_ode_strategy = "hilligardt_eq333"
        self.freeboard_u_bed_top = None
        self.freeboard_d_b_bed_top = None
        self.freeboard_height_from_bed = None
        self.freeboard_u_gb_scale = 1.0
        self.freeboard_explicit_char_hetero_enabled = True
        self.freeboard_eps_d_voidage = None
        self.dense_segment_eps_d_voidage = None
        self.T, self.P = 1200.0, 2_500_000.0
        self.enable_r12 = True
        self.r4_scale = 1.0
        self.r5_scale = 1.0
        self.r6_scale = 1.0
        self.r7_scale = 1.0
        self.nr_reaction_rate_multiplier = 1.0
        self.N_b, self.N_d, self.N_b_in, self.N_d_in = np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS)
        self.N_zu_b, self.N_zu_d, self.N_rez_b, self.N_rez_d = np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS), np.zeros(N_GAS)
        nk = self.solid.n_size_classes
        self.freeboard_explicit_char_sink_applied_kg = np.zeros(nk, dtype=np.float64)
        self.m_solid, self.m_solid_in, self.m_solid_zu = np.zeros((nk, N_SOLID_COMP)), np.zeros((nk, N_SOLID_COMP)), np.zeros((nk, N_SOLID_COMP))
        self.m_solid_rez = np.zeros((nk, N_SOLID_COMP))
        self.m_solid_auf_in = np.zeros((nk, N_SOLID_COMP))
        self.m_solid_ab_in = np.zeros((nk, N_SOLID_COMP))
        self.T_solid_auf_in = 293.15
        self.T_solid_ab_in = 293.15
        self.K_solid_auf = np.zeros((nk, N_SOLID_COMP))
        self.K_solid_ab = np.zeros((nk, N_SOLID_COMP))
        self.u_mf = self.u_b = self.d_b = self.eps_b = self.eps_d = self.K_bd = self.V_b = self.V_d = self.u0 = 0.0
        self.u_d = 0.0
        self.n_rz = 0.0
        self.eps_d_voidage = 0.0
        self._vorab_bottom_gas_inlet_dense_frac = np.nan
        self.R_gas_b, self.R_gas_d = np.zeros(N_GAS), np.zeros(N_GAS); self.R_solid = np.zeros((nk, N_SOLID_COMP))
        self._vm_gas_source_cache = np.zeros(N_GAS); self._vm_solid_sink_cache = np.zeros((nk, N_SOLID_COMP)); self._vm_cache_valid = False
        self._vm_cache_T = np.nan; self._vm_cache_tau = np.nan
        self._vm_cache_T_init = np.nan
        self._vm_cache_m_vm_in = np.nan; self._vm_cache_m_moist_in = np.nan
        self._hydro_cache_valid = False
        self._hydro_cache_T = np.nan
        self._hydro_cache_P = np.nan
        self._hydro_cache_u0_target = np.nan
        self._hydro_cache_N_b = np.zeros(N_GAS)
        self._hydro_cache_N_d = np.zeros(N_GAS)
        self._hydro_cache_options = ("", "", "", "", "", "", "", "", "", "", "", "", "")
        self._thermo_cache_valid = False
        self._thermo_cache_T = np.nan
        self._thermo_cache_P = np.nan
        self._thermo_cache_N_b = np.zeros(N_GAS)
        self._thermo_cache_N_d = np.zeros(N_GAS)
        self._thermo_cache_y_b = np.zeros(N_GAS)
        self._thermo_cache_y_d = np.zeros(N_GAS)
        self._thermo_cache_C_b = np.zeros(N_GAS)
        self._thermo_cache_C_d = np.zeros(N_GAS)
        self._thermo_cache_D_g = 0.0
        self._vorab_hydro_cache_valid = False
        self._freeze_vorabrechnung_inner_nr = False
        self._vorab_hydro_cache = {
            "u_mf": 0.0,
            "u_b": 0.0,
            "d_b": 0.0,
            "eps_b": 0.0,
            "eps_d": 0.0,
            "eps_d_voidage": 0.0,
            "K_bd": 0.0,
            "V_b": 0.0,
            "V_d": 0.0,
            "u0": 0.0,
            "u_d": 0.0,
            "n_rz": 0.0,
        }
        self._vorab_solid_transport_cache_valid = False
        self._vorab_K_solid_auf = np.zeros((nk, N_SOLID_COMP))
        self._vorab_K_solid_ab = np.zeros((nk, N_SOLID_COMP))
        self._minor_gibbs_warm_start: dict[str, dict[str, float | np.ndarray]] = {
            "S": {"ln_N": 0.0, "lambda": np.zeros(4, dtype=np.float64)},
            "N": {"ln_N": 0.0, "lambda": np.zeros(4, dtype=np.float64)},
        }
        self._minor_gibbs_cache_valid = False
        self._minor_gibbs_cache_signature = None
        self._minor_gibbs_cache_out: dict[str, float] = {}
        self._h_cache: Dict[float, npt.NDArray[np.float64]] = {}
        self.N_ex = np.zeros(N_GAS); self._work_res = np.zeros(2 * N_GAS + nk * N_SOLID_COMP + 1)
        self.T_zu_gas = self.T_rez_gas = self.T_zu_solid = self.T_rez_solid = self.T_in_gas = self.T_in_solid = 293.15

    def _mole_fractions(self, phase: str) -> npt.NDArray[np.float64]:
        """气相摩尔分数。

        phase
            ``"b"`` / ``"d"``：单相；``"combined"`` / ``"total"``：气泡相与悬浮相摩尔流加权和
            （与实验出口气取样 / 全截面混合更接近，见 ``Reactor.solve`` 的 exit_gas）。
        """
        return gas_mole_fractions(self.N_b, self.N_d, phase)

    def _concentrations(self, phase: str) -> npt.NDArray[np.float64]:
        return gas_concentrations(self.N_b, self.N_d, phase, self.P, self.T)

    def _thermo_cache_matches(self) -> bool:
        return thermo_cache_matches(self)

    def _get_local_thermo_bundle(self) -> dict[str, npt.NDArray[np.float64] | float]:
        return get_local_thermo_bundle(self)

    def _hydrodynamics_cache_matches(self) -> bool:
        return hydrodynamics_cache_matches(self)

    def _store_hydrodynamics_cache(self) -> None:
        store_hydrodynamics_cache(self)

    def _snapshot_vorabrechnung_hydrodynamics(self) -> None:
        snapshot_vorabrechnung_hydrodynamics(self)

    def _apply_frozen_vorabrechnung_hydrodynamics(self) -> None:
        apply_frozen_vorabrechnung_hydrodynamics(self)

    def enable_inner_nr_vorabrechnung_freeze(self) -> None:
        enable_inner_nr_vorabrechnung_freeze(self)

    def disable_inner_nr_vorabrechnung_freeze(self) -> None:
        disable_inner_nr_vorabrechnung_freeze(self)

    def calc_hydrodynamics(self) -> None:
        """计算流体力学状态（eps_b, eps_d, u_mf, d_b, u_b, K_bd）。
        
        Ref: Hamel (1999) Eq. 3.44, 3.50, 3.52
        """
        if self._hydrodynamics_cache_matches():
            return
        bundle = calc_cell_hydrodynamics(
            T=self.T,
            P=self.P,
            N_b=self.N_b,
            N_d=self.N_d,
            D_bed=self.geo.D_bed,
            dh=self.geo.dh,
            h_center=self.geo.h_center,
            N_or=self.geo.N_or,
            rho_s=self.solid.rho_s,
            d_p=self.solid.d_p,
            eps_mf=self.solid.eps_mf,
            phi_s=self.solid.phi_s,
            u0_target=self.u0_target,
            u_d_closure=self.u_d_closure,
            bubble_diameter_model=self.bubble_diameter_model,
            psi_b_strategy=self.psi_b_strategy,
            lambda_strategy=self.lambda_strategy,
            xi_strategy=self.xi_strategy,
            bubble_velocity_strategy=self.bubble_velocity_strategy,
            bubble_ode_strategy=self.bubble_ode_strategy,
            cell_type=self.cell_type,
            freeboard_u_bed_top=self.freeboard_u_bed_top,
            freeboard_d_b_bed_top=self.freeboard_d_b_bed_top,
            freeboard_height_from_bed=self.freeboard_height_from_bed,
            freeboard_u_gb_scale=self.freeboard_u_gb_scale,
            freeboard_eps_d_voidage=self.freeboard_eps_d_voidage,
            dense_segment_eps_d_voidage=self.dense_segment_eps_d_voidage,
        )
        self.u_mf = bundle.u_mf
        self.u_b = bundle.u_b
        self.d_b = bundle.d_b
        self.eps_b = bundle.eps_b
        self.eps_d = bundle.eps_d
        self.eps_d_voidage = bundle.eps_d_voidage
        self.K_bd = bundle.K_bd
        self.V_b = bundle.V_b
        self.V_d = bundle.V_d
        self.u0 = bundle.u0
        self.u_d = bundle.u_d
        self.n_rz = bundle.n_rz
        self._store_hydrodynamics_cache()

    def calc_exchange(self) -> None:
        thermo = self._get_local_thermo_bundle()
        self.N_ex = calc_phase_exchange(
            K_bd=self.K_bd,
            V_b=self.V_b,
            C_b=np.asarray(thermo["C_b"], dtype=np.float64),
            C_d=np.asarray(thermo["C_d"], dtype=np.float64),
        )

    def _catalyst_bulk_density(self) -> float:
        """催化剂质量密度 [kg/m³ 床层]，用于 R11 悬浮相 Corella 速率。

        Source: specs/04_kinetics.md R11 悬浮相
        """
        return self.solid.rho_s * (1.0 - self.solid.eps_mf) * max(self.solid.catalyst_solid_fraction, 0.0)

    def _element_moles_in_gas(self, element: str) -> float:
        """悬浮相气体中某元素的总摩尔流 [mol/s]（仅 GAS_SPECIES 状态向量）。"""
        return element_moles_in_gas(self, element)

    def _tracked_minor_element_moles_in_gas(self, element: str) -> float:
        """11 组分状态向量中可参与微量 Gibbs 子系统的元素库存 [mol/s]。"""
        return tracked_minor_element_moles_in_gas(self, element)

    def get_element_release(self, element: str) -> float:
        """可用于 Gibbs 子集衡算的元素量 [mol/s]（气相 + 挥发分释放近似）。"""
        return get_element_release_impl(self, element)

    def calc_gibbs_correction(
        self,
        reaction_id: str,
        phase: str,
        clamp_irreversible: bool = False,
    ) -> float:
        """热力学驱动力诊断辅助 `(1 − Q_p/K_eq)`。

        当前 thesis-aligned 主动力学里仅 R8 显式使用等价平衡驱动力；
        本辅助函数保留给诊断/测试与旧规格书对照。

        Parameters
        ----------
        reaction_id : str
            R5, R7, R8 之一
        phase : str
            'b' 气泡相 或 'd' 悬浮相
        clamp_irreversible : bool
            R5 用 True；R7/R8 用 False
        """
        y = self._mole_fractions("b" if phase == "b" else "d")
        return float(
            calc_gibbs_driving_force(
                reaction_id,
                self.T,
                self.P,
                y,
                GAS_SPECIES_INDEX,
                clamp_irreversible=clamp_irreversible,
            )
        )

    def calc_minor_species_gibbs(
        self,
        use_sulfur: bool = True,
        use_nitrogen: bool = True,
        k_relax: float | None = None,
    ) -> dict[str, float]:
        """Gibbs 微量组分松弛源项 [mol/s]，写入 R_gas_d。"""
        k = _GIBBS_RELAX_K if k_relax is None else float(k_relax)
        return calc_minor_species_gibbs_impl(
            self,
            use_sulfur=use_sulfur,
            use_nitrogen=use_nitrogen,
            k_relax=k,
        )

    def calc_reactions(self, rate_multiplier: float = 1.0) -> None:
        """计算反应源项。

        Parameters
        ----------
        rate_multiplier:
            所有反应速率的比例因子（0–1）。用于源项同伦延拓：在
            ``solve_cell`` 的 homotopy 循环中从 0 逐步升至 1，避免
            刚性初始条件下 least_squares 跳入错误分支。
            默认 1.0（完整速率，与原始行为完全一致）。
        """
        if self.cell_type in {"cyclone", "return_leg"}:
            self.R_gas_b.fill(0.0)
            self.R_gas_d.fill(0.0)
            self.R_solid.fill(0.0)
            return
        T, P = self.T, self.P
        thermo = self._get_local_thermo_bundle()
        C_b = np.asarray(thermo["C_b"], dtype=np.float64)
        C_d = np.asarray(thermo["C_d"], dtype=np.float64)
        y_b = np.asarray(thermo["y_b"], dtype=np.float64)
        y_d = np.asarray(thermo["y_d"], dtype=np.float64)
        D_g = float(thermo["D_g"])

        tau_val = self.geo.dh / max(self.u_mf, 1e-3)
        if self._vm_cache_valid:
            gas_src_vm = self._vm_gas_source_cache.copy()
            solid_sink_vm = self._vm_solid_sink_cache.copy()
        else:
            gas_src_vm = self._calc_drying_pyrolysis_gas_source(tau_val)
            solid_sink_vm = self.R_solid.copy()
        areas = self._calc_char_surface_area_per_class()
        if self.cell_type == "freeboard" and not bool(getattr(self, "freeboard_explicit_char_hetero_enabled", True)):
            areas = np.zeros_like(areas)
        gibbs_minor = self.calc_minor_species_gibbs() if self.use_gibbs_minor else None

        bundle = build_reaction_sources(
            T=T,
            P=P,
            fuel_type=self.fuel_type,
            V_b=self.V_b,
            V_d=self.V_d,
            C_b=C_b,
            C_d=C_d,
            y_b=y_b,
            y_d=y_d,
            gas_src_vm=gas_src_vm,
            solid_sink_vm=solid_sink_vm,
            areas=areas,
            solid_d_p=self.solid.d_p,
            solid_d_p_classes=self.solid.d_p_classes,
            D_g=D_g,
            char_conversion=self._compute_char_conversion(),
            rho_cat=self._catalyst_bulk_density(),
            enable_r12=self.enable_r12,
            use_gibbs_minor=self.use_gibbs_minor,
            gibbs_minor_sources=gibbs_minor,
            r4_scale=self.r4_scale,
            r5_scale=self.r5_scale,
            r6_scale=self.r6_scale,
            r7_scale=self.r7_scale,
            rate_multiplier=rate_multiplier,
            N_zu_d=self.N_zu_d,
            N_d_in=self.N_d_in,
            N_zu_b=self.N_zu_b,
            N_b_in=self.N_b_in,
            N_rez_d=self.N_rez_d,
            N_rez_b=self.N_rez_b,
            N_ex=self.N_ex,
            solid_shape=self.R_solid.shape,
            char_index=S_CHAR,
        )
        self.R_gas_b[:] = bundle.R_gas_b
        self.R_gas_d[:] = bundle.R_gas_d
        self.R_solid[:, :] = bundle.R_solid

    def _calc_char_surface_area_per_class(self) -> npt.NDArray[np.float64]:
        """计算各粒径类的炭表面积 [m²]。

        Ref: Hamel (1999) §5.1，SPM/SCM 炭颗粒表面积
        legacy_stream 模式下，床层稳态炭存量估计：M_bed * char_mass_fraction，
        其中 M_bed = rho_s * (1 - eps_mf) * V_d（悬浮相固体总量），
        char_mass_fraction = m_char_flow / m_total_flow（流率比 = 存量比，假设稳态）。
        holdup_transport/freeboard_closure 模式下，m_solid 已是 cell 内 fuel-solid
        holdup [kg]，不能再乘以包含 inert bed material 的总床料库存。
        比表面积 A = 6 * M_char_inventory / (rho_s * d_p)
        """
        m_char_state = np.maximum(self.m_solid[:, S_CHAR], 0.0)
        if str(self.solid_state_model) in {"holdup_transport", "freeboard_closure"}:
            return m_char_state * 6.0 / (self.solid.rho_s * self.solid.d_p_classes)

        V_cell = (np.pi / 4.0 * self.geo.D_bed**2) * self.geo.dh
        V_d_safe = max(self.V_d, 0.05 * V_cell)
        # 床层悬浮相固体总存量 [kg]（稳态值）
        M_bed = self.solid.rho_s * (1.0 - self.solid.eps_mf) * V_d_safe

        m_solid_total = np.sum(np.maximum(self.m_solid, 0.0))
        if m_solid_total < 1e-10:
            return np.zeros(self.solid.n_size_classes)

        char_frac = m_char_state / m_solid_total
        # 各粒径类炭床层存量 [kg]
        m_char_inventory = M_bed * char_frac  # [kg]
        areas = m_char_inventory * 6.0 / (self.solid.rho_s * self.solid.d_p_classes)
        return areas

    def _solid_upflow_rates(self) -> npt.NDArray[np.float64]:
        """Current upward solid transport rates [kg/s] under the active state model."""
        if str(self.solid_state_model) in {"holdup_transport", "freeboard_closure"}:
            return np.maximum(self.m_solid, 0.0) * np.maximum(self.K_solid_auf, 0.0)
        return np.maximum(self.m_solid, 0.0)

    def _solid_downflow_rates(self) -> npt.NDArray[np.float64]:
        """Current downward solid transport rates [kg/s] under the active state model."""
        if str(self.solid_state_model) in {"holdup_transport", "freeboard_closure"}:
            return np.maximum(self.m_solid, 0.0) * np.maximum(self.K_solid_ab, 0.0)
        return np.zeros_like(self.m_solid)

    def _solid_outflow_rates(self) -> npt.NDArray[np.float64]:
        """Current total outgoing solid transport rates [kg/s] under the active state model."""
        return self._solid_upflow_rates() + self._solid_downflow_rates()

    def _compute_char_conversion(self) -> float:
        if str(self.solid_state_model) == "freeboard_closure":
            # Freeboard solids are supplied by the trajectory closure rather than
            # Eulerian inlet holdup balances, so there is no consistent local
            # feed/outlet budget for a per-cell SPM conversion estimate here.
            return 0.0
        m_in = np.sum(np.maximum(self.m_solid_zu[:, S_CHAR] + self.m_solid_in[:, S_CHAR] + self.m_solid_rez[:, S_CHAR], 0.0))
        m_char_out = float(np.sum(self._solid_outflow_rates()[:, S_CHAR]))
        x_char = 1.0 - m_char_out / max(float(m_in), 1e-12)
        if x_char < 0.0:
            return 0.0
        if x_char > 1.0:
            return 1.0
        return float(x_char)

    def calc_gas_balance(self) -> npt.NDArray[np.float64]:
        """气相摩尔守恒残差向量。
        
        Ref: Hamel (1999) Eq. 2.1, 2.2
        """
        return calc_gas_balance_residual(
            N_zu_d=self.N_zu_d,
            N_rez_d=self.N_rez_d,
            N_d_in=self.N_d_in,
            R_gas_d=self.R_gas_d,
            N_d=self.N_d,
            N_ex=self.N_ex,
            N_zu_b=self.N_zu_b,
            N_rez_b=self.N_rez_b,
            N_b_in=self.N_b_in,
            R_gas_b=self.R_gas_b,
            N_b=self.N_b,
        )

    def calc_solid_balance(self) -> npt.NDArray[np.float64]:
        if str(self.solid_state_model) == "freeboard_closure":
            return np.zeros_like(self.m_solid)
        mig = self._calc_size_migration()
        return calc_solid_balance_residual(
            m_solid_zu=self.m_solid_zu,
            m_solid_rez=self.m_solid_rez,
            m_solid_in=self.m_solid_in,
            m_solid_auf_in=self.m_solid_auf_in,
            m_solid_ab_in=self.m_solid_ab_in,
            R_solid=self.R_solid,
            size_migration=mig,
            m_solid=self.m_solid,
            K_solid_auf=self.K_solid_auf,
            K_solid_ab=self.K_solid_ab,
            solid_state_model=str(self.solid_state_model),
        )

    def _calc_size_migration(self) -> npt.NDArray[np.float64]:
        """计算颗粒由于缩减导致的跨级迁移量 [kg/s]。
        
        Ref: Hamel (1999) Eq. 2.3 - 2.6
        """
        areas = self._calc_char_surface_area_per_class()
        V_cell = (np.pi/4.0 * self.geo.D_bed**2) * self.geo.dh
        return calc_size_migration(
            m_solid=self.m_solid,
            R_solid=self.R_solid,
            areas=areas,
            d_p_classes=self.solid.d_p_classes,
            rho_s=self.solid.rho_s,
            eps_mf=self.solid.eps_mf,
            V_d=self.V_d,
            V_cell=V_cell,
            char_index=S_CHAR,
        )

    def calc_energy_balance(self) -> float:
        """全床能量平衡残差 [W]。
        
        Ref: Hamel (1999) Eq. 2.7
        """
        if self.cell_type in {"cyclone", "return_leg"}:
            return self._calc_side_temperature_closure_residual()
        return calc_energy_balance_residual(
            N_b_in=self.N_b_in,
            N_d_in=self.N_d_in,
            T_in_gas=self.T_in_gas,
            N_zu_b=self.N_zu_b,
            N_zu_d=self.N_zu_d,
            T_zu_gas=self.T_zu_gas,
            N_rez_b=self.N_rez_b,
            N_rez_d=self.N_rez_d,
            T_rez_gas=self.T_rez_gas,
            m_solid_rez=self.m_solid_rez,
            T_rez_solid=self.T_rez_solid,
            m_solid_in=self.m_solid_in,
            T_in_solid=self.T_in_solid,
            m_solid_auf_in=self.m_solid_auf_in,
            T_solid_auf_in=self.T_solid_auf_in,
            m_solid_ab_in=self.m_solid_ab_in,
            T_solid_ab_in=self.T_solid_ab_in,
            m_solid_zu=self.m_solid_zu,
            T_zu_solid=self.T_zu_solid,
            N_b=self.N_b,
            N_d=self.N_d,
            m_solid=self._solid_outflow_rates(),
            T=self.T,
            heat_loss_frac=self.heat_loss_frac,
            ash_dry_wt=self.solid.ash_dry_wt,
            VM_daf=self.solid.VM_daf,
            h_f_dry=self.solid.h_f_dry,
            h_cache=self._h_cache,
        )

    def _calc_gas_enthalpy_flow(self, N: npt.NDArray[np.float64], T: float) -> float:
        return calc_gas_enthalpy_flow(N, T, self._h_cache)

    def _calc_solid_enthalpy_flow(self, m: npt.NDArray[np.float64], T: float) -> float:
        return solid_enthalpy_flow_for_fuel_props(
            m,
            T,
            ash_dry_wt=self.solid.ash_dry_wt,
            VM_daf=self.solid.VM_daf,
            h_f_dry=self.solid.h_f_dry,
        )

    def _calc_side_temperature_closure_residual(self) -> float:
        """Side-element temperature closure residual [W-equivalent].

        Cyclone / return-leg blocks are separators and recycle topology elements,
        not reactive bubbling-bed cells.  Their NR temperature row therefore
        enforces pass-through thermal closure instead of a full cell energy
        balance with reaction/phase terms.
        """
        if self.cell_type == "cyclone" and float(np.sum(np.maximum(self.N_b_in + self.N_d_in, 0.0))) > 1e-12:
            T_ref = float(self.T_in_gas)
        else:
            solid_ab = float(np.sum(np.maximum(self.m_solid_ab_in, 0.0)))
            solid_auf = float(np.sum(np.maximum(self.m_solid_auf_in, 0.0)))
            if solid_ab > 1e-12:
                T_ref = float(self.T_solid_ab_in)
            elif solid_auf > 1e-12:
                T_ref = float(self.T_solid_auf_in)
            else:
                T_ref = float(self.T_in_solid)

        gas_in = np.maximum(self.N_b_in + self.N_d_in, 0.0)
        solid_in = np.maximum(self.m_solid_in + self.m_solid_auf_in + self.m_solid_ab_in, 0.0)
        dT = 1.0
        cp_dot = 0.0
        if float(np.sum(gas_in)) > 1e-12:
            cp_dot += abs(
                self._calc_gas_enthalpy_flow(gas_in, T_ref + dT)
                - self._calc_gas_enthalpy_flow(gas_in, T_ref)
            ) / dT
        if float(np.sum(solid_in)) > 1e-12:
            cp_dot += abs(
                self._calc_solid_enthalpy_flow(solid_in, T_ref + dT)
                - self._calc_solid_enthalpy_flow(solid_in, T_ref)
            ) / dT
        cp_dot = max(float(cp_dot), 1.0e3)
        return cp_dot * (T_ref - float(self.T))

    def residuals(self, rate_multiplier: float | None = None) -> npt.NDArray[np.float64]:
        rm = self.nr_reaction_rate_multiplier if rate_multiplier is None else float(rate_multiplier)
        return execute_cell_residual_pipeline(self, rate_multiplier=rm)

    def _calc_drying_pyrolysis_gas_source(self, tau: float) -> npt.NDArray[np.float64]:
        bundle = calc_drying_pyrolysis_sources(
            tau=tau,
            T=self.T,
            P=self.P,
            T_init=self.T_in_solid,
            d_p=self.solid.d_p,
            moisture_wt=self.solid.moisture_wt,
            ash_dry_wt=self.solid.ash_dry_wt,
            C_dry=self.solid.C_dry,
            H_dry=self.solid.H_dry,
            O_dry=self.solid.O_dry,
            nitrogen_fraction=self.solid.nitrogen_fraction,
            sulfur_fraction=self.solid.sulfur_fraction,
            sulfur_volatile_frac=self.solid.sulfur_volatile_frac,
            pyrolysis_tar_carbon_frac=self.solid.pyrolysis_tar_carbon_frac,
            fuel_type=self.fuel_type,
            m_vm_in=float(np.sum(self.m_solid_zu[:, S_VM] + self.m_solid_in[:, S_VM])),
            m_moist_in=float(np.sum(self.m_solid_zu[:, S_MOISTURE] + self.m_solid_in[:, S_MOISTURE])),
            solid_shape=self.R_solid.shape,
            m_vm_in_classes=np.maximum(self.m_solid_zu[:, S_VM] + self.m_solid_in[:, S_VM], 0.0),
            m_moist_in_classes=np.maximum(self.m_solid_zu[:, S_MOISTURE] + self.m_solid_in[:, S_MOISTURE], 0.0),
            char_index=S_CHAR,
            vm_index=S_VM,
            moisture_index=S_MOISTURE,
        )
        self.R_solid[:, :] = bundle.solid_sink
        return bundle.gas_source

    def _allocate_pyrolysis_products_elemental(self, nC, nH, nO, nN=0.0, nS=0.0) -> Dict[str, float]:
        target_tar_hc_ratio = None
        if self.solid.C_dry > 1e-12 and self.solid.H_dry >= 0.0:
            target_tar_hc_ratio = float((self.solid.H_dry / 1.00794) / max(self.solid.C_dry / 12.011, 1e-12))
        return allocate_pyrolysis_products_elemental(
            nC=float(nC),
            nH=float(nH),
            nO=float(nO),
            nN=float(nN),
            nS=float(nS),
            fuel_type=self.fuel_type,
            pyrolysis_tar_carbon_frac=self.solid.pyrolysis_tar_carbon_frac,
            target_tar_hc_ratio=target_tar_hc_ratio,
        )

    def compute_vorabrechnung(self, tau: float) -> None:
        """重建干燥/热解预计算源项缓存。

        Ref: Hamel (1999) Kapitel 2 Vorabrechnung; Kapitel 4

        在 global-NR 主线路径中，这些缓存由外层 Abgleich 显式刷新，随后
        内层 Newton 迭代固定使用，不随 Jacobian 扰动重算。这里保留阈值复用
        逻辑，供 GS / 非强制 refresh 路径复用。
        """
        tau_val = max(float(tau), 1e-4)
        t_init = float(self.T_in_solid)
        m_vm_in = float(np.sum(self.m_solid_zu[:, S_VM] + self.m_solid_in[:, S_VM]))
        m_moist_in = float(np.sum(self.m_solid_zu[:, S_MOISTURE] + self.m_solid_in[:, S_MOISTURE]))
        if self._vm_cache_valid:
            if (
                abs(float(self.T) - float(self._vm_cache_T)) <= 10.0
                and abs(tau_val - float(self._vm_cache_tau)) <= 1e-3
                and abs(t_init - float(self._vm_cache_T_init)) <= 1e-3
                and abs(m_vm_in - float(self._vm_cache_m_vm_in)) <= 1e-8
                and abs(m_moist_in - float(self._vm_cache_m_moist_in)) <= 1e-8
            ):
                self._snapshot_vorabrechnung_hydrodynamics()
                return
        self.R_solid.fill(0.0)
        self._vm_gas_source_cache = self._calc_drying_pyrolysis_gas_source(tau_val)
        self._vm_solid_sink_cache = self.R_solid.copy()
        self._vm_cache_valid = True
        self._vm_cache_T = float(self.T)
        self._vm_cache_tau = tau_val
        self._vm_cache_T_init = t_init
        self._vm_cache_m_vm_in = m_vm_in
        self._vm_cache_m_moist_in = m_moist_in
        self._snapshot_vorabrechnung_hydrodynamics()

    def invalidate_vorabrechnung_cache(self) -> None:
        """显式失效 Vorabrechnung 缓存，供外层 Abgleich 强制重建。"""
        invalidate_vorabrechnung_cache_impl(self)
