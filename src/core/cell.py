"""Cell 类：单个轴向离散单元，两相理论（气泡相 + 悬浮相）。

每个 cell 包含：
  - 气泡相 (b) 和悬浮相 (d) 的各组分摩尔流率
  - 固体颗粒质量流率（按粒径类）
  - 温度 T（假设两相局部热平衡）
  - 流体力学参数（K_bd, u_b, eps_b 等）

六个核心方法：
  1. calc_hydrodynamics()  — u_mf, u_b, eps_b, K_bd
  2. calc_exchange()       — 相间摩尔交换 N_dot_ex
  3. calc_reactions()      — 所有反应源项
  4. calc_gas_balance()    — 气相摩尔守恒残差
  5. calc_solid_balance()  — 固相质量守恒残差
  6. calc_energy_balance() — 全局能量守恒残差

Source: docs/CLAUDE.md Phase 5.1; Hamel & Krumm (2001) Section 2.1;
        specs/01_conservation_equations.md
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

import numpy as np
import numpy.typing as npt

from src.core.constants import Rg, T_REF, g
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    MOLECULAR_WEIGHT,
    N_GAS,
    TarFuelType,
    cp_molar,
    enthalpy_molar,
    gas_density_ideal,
    gas_diffusivity_correlation,
    gas_viscosity_power_law,
    cp_char,
    cp_ash,
    cp_sand,
    TAR_SURROGATE_FORMULA,
    get_tar_component_mapping,
)
from src.kinetics.arrhenius import k_hobbs, k_standard
from src.kinetics.char_reactions import (
    d_core_from_spm_char_conversion,
    rate_R1,
    rate_R2,
    rate_R3,
    rate_R4_effective,
)
from src.kinetics.gas_reactions import (
    rate_R5_bubble,
    rate_R5_suspension,
    rate_R6,
    rate_R7,
    rate_R8,
    rate_R9,
)
from src.kinetics.tar_reactions import (
    TarReactionId,
    get_lumped_tar_stoichiometry,
    rate_R10,
    rate_R11_bubble,
    rate_R11_suspension,
)
from src.core.species import get_atom_count
from src.thermodynamics.equilibrium import calc_gibbs_driving_force
from src.thermodynamics.minor_species import solve_minor_species
from src.thermal.drying import solve_drying_CN
from src.thermal.devolatilization import daem_conversion_radial
from src.physics.bubble_dynamics import bubble_rise_velocity, mori_wen_bubble_diameter
from src.physics.mass_transfer import calc_kbd, calc_u_br
from src.physics.minimum_fluidization import compute_u_mf
from src.physics.phase_fractions import calc_epsilon_b, calc_epsilon_d


@dataclass
class CellGeometry:
    """Cell 几何参数。"""
    D_bed: float = 0.6          # [m]    床层直径
    dh: float = 0.5             # [m]    cell 轴向高度
    h_center: float = 0.0      # [m]    cell 中心高度


@dataclass
class SolidProps:
    """固体颗粒属性。"""
    rho_s: float = 1400.0       # [kg/m³] 颗粒密度
    d_p: float = 0.001          # [m]     Sauter 平均直径
    phi_s: float = 0.86         # [-]     球形度
    eps_mf: float = 0.45        # [-]     最小流化空隙率
    n_size_classes: int = 1     # [-]     粒径类数
    d_p_classes: npt.NDArray[np.float64] = field(
        default_factory=lambda: np.array([0.001])
    )  # [m] 各粒径类直径
    mass_fractions: npt.NDArray[np.float64] = field(
        default_factory=lambda: np.array([1.0])
    )  # [-] 各粒径类质量分数
    sulfur_fraction: float = 0.0   # [-] 燃料硫质量分数（Gibbs 微量组分用）
    sulfur_volatile_frac: float = 0.5  # [-] 燃料硫中挥发硫比例（Hamel：热解释放）
    nitrogen_fraction: float = 0.0  # [-] 燃料氮质量分数（Gibbs 微量组分用）
    # R11 Corella：悬浮相催化速率 ∝ 催化剂质量密度；1 = 乳化相固体均为催化/床料
    catalyst_solid_fraction: float = 1.0  # [-]
    # Chapter 4 干燥/热解耦合所需基础分析参数
    moisture_wt: float = 16.9          # [wt%] 湿基水分
    ash_dry_wt: float = 11.41          # [wt%] 干基灰分
    VM_daf: float = 53.42              # [wt%] daf 挥发分
    C_dry: float = 61.5                # [wt%] 干基碳
    H_dry: float = 4.1                 # [wt%] 干基氢
    O_dry: float = 21.8                # [wt%] 干基氧
    pyrolysis_tar_carbon_frac: float = 0.2   # [-] 挥发碳进入焦油(TAR1)比例


class Cell:
    """单个轴向离散单元，包含气泡相与悬浮相。

    Source: docs/CLAUDE.md Phase 5.1; specs/01_conservation_equations.md
    """

    def __init__(
        self,
        geo: CellGeometry | None = None,
        solid: SolidProps | None = None,
        fuel_type: TarFuelType = "coal",
    ) -> None:
        self.geo = geo or CellGeometry()
        self.solid = solid or SolidProps()
        self.fuel_type: TarFuelType = fuel_type
        self.use_gibbs_minor: bool = False  # True 时用 Gibbs 替代 R9 微量组分

        self.T: float = 1200.0          # [K]    温度
        self.P: float = 2_500_000.0     # [Pa]   压力

        self.N_b = np.zeros(N_GAS)      # [mol/s] 气泡相各组分出口摩尔流率
        self.N_d = np.zeros(N_GAS)      # [mol/s] 悬浮相各组分出口摩尔流率

        # 入口摩尔流率（从上游 cell / 进料设定）
        self.N_b_in = np.zeros(N_GAS)
        self.N_d_in = np.zeros(N_GAS)

        # 进料项 (zu = Zuführung)
        self.N_zu_b = np.zeros(N_GAS)   # [mol/s] 气泡相进料
        self.N_zu_d = np.zeros(N_GAS)   # [mol/s] 悬浮相进料
        # 循环项 (rez = Rezirkulation)
        self.N_rez_b = np.zeros(N_GAS)
        self.N_rez_d = np.zeros(N_GAS)

        # 固体
        nk = self.solid.n_size_classes
        self.m_solid = np.zeros(nk)     # [kg/s] 固体出口质量流率
        self.m_solid_in = np.zeros(nk)  # [kg/s] 固体入口质量流率
        self.m_solid_zu = np.zeros(nk)  # [kg/s] 固体进料
        # 炭含量质量分数：由燃料工业分析计算（湿基固定碳分率）
        # FC_dry = 1 - ash_dry - VM_daf*(1-ash_dry)；FC_wet = FC_dry*(1-moisture)
        # Source: Hamel (1999) §5.1；炭分率影响持料量估算与 R1–R4 反应面积
        _ash = float(np.clip(self.solid.ash_dry_wt / 100.0, 0.0, 0.95))
        _vm = float(np.clip(self.solid.VM_daf / 100.0, 0.0, 0.95))
        _moi = float(np.clip(self.solid.moisture_wt / 100.0, 0.0, 0.95))
        _fc_dry = max(1.0 - _ash - _vm * (1.0 - _ash), 0.0)
        _fc_wet = _fc_dry * (1.0 - _moi)
        _fc_wet = max(_fc_wet, 0.05)  # 下限：避免数值上无炭
        self.m_char_frac = np.ones(nk, dtype=np.float64) * _fc_wet

        # 壁面热损
        self.Q_wall: float = 0.0   # [W]

        # ---------- 流体力学缓存 ----------
        self.u_mf: float = 0.0
        self.u_b: float = 0.0
        self.d_b: float = 0.0
        self.eps_b: float = 0.0
        self.eps_d: float = 1.0
        self.K_bd: float = 0.0
        self.V_b: float = 0.0      # [m³] 气泡相体积
        self.V_d: float = 0.0      # [m³] 悬浮相体积
        self.u0: float = 0.0       # [m/s] 表观气速

        # ---------- 反应源项缓存 ----------
        self.R_gas_b = np.zeros(N_GAS)  # [mol/s] 气泡相气体反应源项
        self.R_gas_d = np.zeros(N_GAS)  # [mol/s] 悬浮相气体反应源项
        self.R_solid = np.zeros(nk)     # [kg/s]  固体反应源项

        # Vorabrechnung 缓存：每外迭代算一次，内层 fsolve 残差中固定（Hamel 1999 Module 2）
        self._vm_gas_source_cache = np.zeros(N_GAS)   # [mol/s] 干燥+热解气相源
        self._vm_solid_sink_cache = np.zeros(nk)      # [kg/s] 干燥+热解固相汇（与 R_solid 同号）
        self._vm_cache_valid: bool = False

        # ---------- 相间交换缓存 ----------
        self.N_ex = np.zeros(N_GAS)     # [mol/s] b->d 交换

        # ---------- 预分配工作数组（减少 residuals 链中的分配） ----------
        self._work_y_b = np.zeros(N_GAS)
        self._work_y_d = np.zeros(N_GAS)
        self._work_C_b = np.zeros(N_GAS)
        self._work_C_d = np.zeros(N_GAS)
        self._work_R_b = np.zeros(N_GAS)
        self._work_R_d = np.zeros(N_GAS)
        self._work_tar_src = np.zeros(N_GAS)
        self._work_res = np.zeros(2 * N_GAS + nk + 1)
        self._work_h_array = np.zeros(N_GAS)   # 焓流向量化
        self._h_array_T_cache = -1.0            # 缓存温度，避免重复计算
        self._mu_g_cache = -1.0
        self._D_g_cache = -1.0
        self._hydro_TP_cache = (-1.0, -1.0)    # (T, P) 物性缓存

        # 各流股温度 [K]（能量衡算）；由 Reactor 在每步扫描前更新
        self.T_zu_gas: float = 293.15   # 新鲜进料气温度
        self.T_rez_gas: float = 293.15  # 循环气温度（与顶部出口气一致，由 Reactor 更新）
        self.T_zu_solid: float = 293.15  # 固体进料温度
        self.T_in_gas: float = self.T    # 来自下游 cell 的气体入口
        self.T_in_solid: float = self.T  # 来自上游 cell 的固体入口

    def compute_vorabrechnung(self, tau_cell: float) -> None:
        """在**外迭代**开始时计算干燥+热解源项一次，内层 Newton–Raphson（fsolve）中保持固定。

        将结果写入 `_vm_gas_source_cache`、`_vm_solid_sink_cache`，并置 `_vm_cache_valid`。
        避免在每次残差/Jacobian 扰动中重复 Crank–Nicolson + DAEM（Hamel 1999 Module 2）。

        调用方（Reactor.solve）须先 `calc_hydrodynamics()`，以便 `tau_cell`、`V_d` 等与当前状态一致。

        Parameters
        ----------
        tau_cell : float
            有效停留时间 [s]，通常取 `dh / u_mf`。
        """
        self.R_solid.fill(0.0)
        gas_src = self._calc_drying_pyrolysis_gas_source(tau_cell)
        self._vm_gas_source_cache[:] = gas_src
        self._vm_solid_sink_cache[:] = self.R_solid
        self._vm_cache_valid = True

    # ===================================================================
    # 辅助：浓度计算
    # ===================================================================

    def _total_gas_molar_flow(self, phase: str) -> float:
        N = self.N_b if phase == "b" else self.N_d
        return float(np.sum(np.maximum(N, 0.0))) + 1e-30

    def _mole_fractions(self, phase: str) -> npt.NDArray[np.float64]:
        """摩尔分数（湿基），写入预分配工作数组并返回。"""
        work = self._work_y_b if phase == "b" else self._work_y_d
        N = self.N_b if phase == "b" else self.N_d
        np.maximum(N, 0.0, out=work)
        total = float(np.sum(work))
        if total < 1e-30:
            work.fill(0.0)
            work[GAS_SPECIES_INDEX["N2"]] = 1.0
            return work
        work /= total
        return work

    def _concentrations(self, phase: str) -> npt.NDArray[np.float64]:
        """摩尔浓度 [mol/m³] = y_j * P / (Rg * T)，写入预分配工作数组。"""
        y = self._mole_fractions(phase)
        work = self._work_C_b if phase == "b" else self._work_C_d
        np.multiply(y, self.P / (Rg * self.T), out=work)
        return work

    # ===================================================================
    # 1. calc_hydrodynamics
    # ===================================================================

    def calc_hydrodynamics(self) -> None:
        """计算 u_mf, u_b, d_b, eps_b, K_bd, V_b, V_d。

        Source: specs/02_hydrodynamics.md; Hamel (1999) Eq.4.2-4.12
        """
        A_bed = np.pi / 4.0 * self.geo.D_bed**2
        T, P = self.T, self.P

        # mu_g, D_g 仅依赖 T, P，可缓存
        if abs(T - self._hydro_TP_cache[0]) < 1e-4 and abs(P - self._hydro_TP_cache[1]) < 1.0:
            mu_g = self._mu_g_cache
            D_g = self._D_g_cache
        else:
            mu_g_20 = 1.8e-5  # TODO: 从进料气体组成确定
            mu_g = gas_viscosity_power_law(T, mu_g_20)
            D_g = gas_diffusivity_correlation(T, P)
            self._mu_g_cache = mu_g
            self._D_g_cache = D_g
            self._hydro_TP_cache = (T, P)

        y_d = self._mole_fractions("d")
        mole_frac_dict = {sp: float(y_d[i]) for i, sp in enumerate(GAS_SPECIES)}
        rho_g = gas_density_ideal(P, T, mole_frac_dict)

        self.u_mf = compute_u_mf(
            rho_g, self.solid.rho_s, self.solid.d_p, mu_g,
            self.solid.eps_mf, self.solid.phi_s,
        )

        total_gas_mol_s = self._total_gas_molar_flow("d") + self._total_gas_molar_flow("b")
        Q_gas = total_gas_mol_s * Rg * self.T / self.P  # [m³/s]
        self.u0 = Q_gas / max(A_bed, 1e-10)

        self.d_b = mori_wen_bubble_diameter(
            self.geo.h_center, self.u0, self.u_mf, self.geo.D_bed,
        )
        self.u_b = bubble_rise_velocity(self.u0, self.u_mf, self.d_b)

        self.eps_b = calc_epsilon_b(self.u0, self.u_mf, self.u_b)
        self.eps_d = calc_epsilon_d(self.eps_b)

        V_cell = A_bed * self.geo.dh
        self.V_b = self.eps_b * V_cell
        self.V_d = self.eps_d * V_cell

        u_d = self.u_mf / max(self.solid.eps_mf, 0.01)
        u_br = calc_u_br(u_d, self.P)
        self.K_bd = calc_kbd(u_br, self.d_b, D_g, self.solid.eps_mf, self.u_b)

    # ===================================================================
    # 2. calc_exchange
    # ===================================================================

    def calc_exchange(self) -> None:
        """计算相间气体摩尔交换 N_dot_ex [mol/s]。

        N_ex,j = K_bd * V_b * (C_j,b - C_j,d)

        正值 = b -> d（气泡相流出，悬浮相流入）

        Source: specs/01_conservation_equations.md §1; Hamel (1999) Eq.1-1
        """
        C_b = self._concentrations("b")
        C_d = self._concentrations("d")
        self.N_ex = self.K_bd * self.V_b * (C_b - C_d)

    # ===================================================================
    # 3. calc_reactions
    # ===================================================================

    def calc_reactions(self) -> None:
        """计算所有反应源项。

        气泡相：仅均相 R5b, R6, R7, R8, R9, R10, R11
        悬浮相：均相 + 非均相 R1-R4, R5d, R6-R11

        Source: specs/04_kinetics.md; docs/CLAUDE.md Phase 3
        """
        C_b = self._concentrations("b")
        C_d = self._concentrations("d")
        y_b = self._mole_fractions("b")
        y_d = self._mole_fractions("d")
        T = self.T

        idx = GAS_SPECIES_INDEX
        D_g = gas_diffusivity_correlation(T, self.P)
        self.R_solid.fill(0.0)

        # --- 气泡相源项 [mol/(m³·s)] -> [mol/s] 乘以 V_b ---
        R_b = self._work_R_b
        R_b.fill(0.0)
        r5b = rate_R5_bubble(T, C_b[idx["CO"]], C_b[idx["O2"]], self.P, y_b, idx)
        R_b[idx["CO"]] -= 2.0 * r5b
        R_b[idx["O2"]] -= r5b
        R_b[idx["CO2"]] += 2.0 * r5b

        r6b = rate_R6(T, C_b[idx["CH4"]], C_b[idx["O2"]])
        R_b[idx["CH4"]] -= r6b
        R_b[idx["O2"]] -= 1.5 * r6b
        R_b[idx["CO"]] += r6b
        R_b[idx["H2O"]] += 2.0 * r6b

        r7b = rate_R7(T, C_b[idx["CH4"]], C_b[idx["H2O"]], self.P, y_b, idx)
        R_b[idx["CH4"]] -= r7b
        R_b[idx["H2O"]] -= r7b
        R_b[idx["CO"]] += r7b
        R_b[idx["H2"]] += 3.0 * r7b

        r8b = rate_R8(T, self.P, y_b[idx["CO"]], y_b[idx["H2O"]],
                       y_b[idx["CO2"]], y_b[idx["H2"]])
        R_b[idx["CO"]] -= r8b
        R_b[idx["H2O"]] -= r8b
        R_b[idx["CO2"]] += r8b
        R_b[idx["H2"]] += r8b

        # R10/R11 气泡相：R10 两相同动力学；R11 = Serio 均相热裂解（Eq.5.59 / Table 5.4）
        C_tar_mix_b = C_b[idx["TAR1"]] + C_b[idx["TAR2"]]
        r10b = rate_R10(T, C_tar_mix_b, C_b[idx["O2"]], self.P, self.fuel_type)
        r11b = rate_R11_bubble(T, C_tar_mix_b)
        tar_src_b = self._calc_tar_source({"R10": r10b, "R11": r11b})
        R_b += tar_src_b

        self.R_gas_b = R_b * self.V_b  # [mol/s]

        # --- 悬浮相源项 ---
        R_d = self._work_R_d
        R_d.fill(0.0)

        # Chapter 4：先干燥（Agarwal）再热解（DAEM）
        # Vorabrechnung 缓存：由 Reactor 每外迭代调用 compute_vorabrechnung，残差内固定；
        # 无缓存时回退为即时计算（单 cell 测试、首次调用等）。
        if self._vm_cache_valid:
            R_d += self._vm_gas_source_cache / max(self.V_d, 1e-20)
            self.R_solid += self._vm_solid_sink_cache
        else:
            tau_cell = self.geo.dh / max(self.u_mf, 1e-3)
            dry_py_gas = self._calc_drying_pyrolysis_gas_source(tau_cell)
            R_d += dry_py_gas / max(self.V_d, 1e-20)

        r5d = rate_R5_suspension(T, C_d[idx["CO"]], C_d[idx["O2"]], C_d[idx["H2O"]], self.P, y_d, idx)
        R_d[idx["CO"]] -= 2.0 * r5d
        R_d[idx["O2"]] -= r5d
        R_d[idx["CO2"]] += 2.0 * r5d

        r6d = rate_R6(T, C_d[idx["CH4"]], C_d[idx["O2"]])
        R_d[idx["CH4"]] -= r6d
        R_d[idx["O2"]] -= 1.5 * r6d
        R_d[idx["CO"]] += r6d
        R_d[idx["H2O"]] += 2.0 * r6d

        r7d = rate_R7(T, C_d[idx["CH4"]], C_d[idx["H2O"]], self.P, y_d, idx)
        R_d[idx["CH4"]] -= r7d
        R_d[idx["H2O"]] -= r7d
        R_d[idx["CO"]] += r7d
        R_d[idx["H2"]] += 3.0 * r7d

        r8d = rate_R8(T, self.P, y_d[idx["CO"]], y_d[idx["H2O"]],
                       y_d[idx["CO2"]], y_d[idx["H2"]])
        R_d[idx["CO"]] -= r8d
        R_d[idx["H2O"]] -= r8d
        R_d[idx["CO2"]] += r8d
        R_d[idx["H2"]] += r8d

        # R1-R4 异相反应（先于 Gibbs，使 R_solid 可用于 get_element_release）
        # 外表面积由乳化相床层持料量（inventory）估算，而非 m_solid 流率 [kg/s]
        char_area_per_class = self._calc_char_surface_area_per_class()
        total_char_area = float(np.sum(char_area_per_class))
        if total_char_area > 0:
            P_CO2 = y_d[idx["CO2"]] * self.P
            P_CO = y_d[idx["CO"]] * self.P

            X_char = self._compute_char_conversion()
            d_core = d_core_from_spm_char_conversion(X_char, self.solid.d_p)
            r1, alpha = rate_R1(
                T,
                C_d[idx["O2"]],
                self.solid.d_p,
                D_g,
                d_core,
                fuel=self.fuel_type,
            )
            R_d[idx["O2"]] -= alpha * r1 * total_char_area
            R_d[idx["CO"]] += 2.0 * (1.0 - alpha) * r1 * total_char_area
            R_d[idx["CO2"]] += (2.0 * alpha - 1.0) * r1 * total_char_area

            r2 = rate_R2(T, C_d[idx["H2O"]], self.solid.d_p, D_g, d_core)
            R_d[idx["H2O"]] -= r2 * total_char_area
            R_d[idx["CO"]] += r2 * total_char_area
            R_d[idx["H2"]] += r2 * total_char_area

            r3 = rate_R3(T, C_d[idx["H2"]], self.solid.d_p, D_g, d_core)
            R_d[idx["H2"]] -= 2.0 * r3 * total_char_area
            R_d[idx["CH4"]] += r3 * total_char_area

            r4 = rate_R4_effective(T, P_CO2, P_CO, self.solid.d_p, D_g)
            R_d[idx["CO2"]] -= r4 * total_char_area
            R_d[idx["CO"]] += 2.0 * r4 * total_char_area

            M_C = 12.011e-3  # [kg/mol]
            # total_r_char [kg/(m²·s)]；R_solid[k] [kg/s] = -total_r_char * A_k
            total_r_char = (r1 + r2 + r3 + r4) * M_C
            for k in range(self.solid.n_size_classes):
                self.R_solid[k] += -total_r_char * char_area_per_class[k]

        if not self.use_gibbs_minor:
            r9d = rate_R9(T, C_d[idx["H2S"]], C_d[idx["O2"]])
            R_d[idx["H2S"]] -= r9d
            R_d[idx["O2"]] -= 1.5 * r9d
            R_d[idx["H2O"]] += r9d
        else:
            gibbs_result = self.calc_minor_species_gibbs()
            k_relax = 1.0
            V_d = self.V_d if self.V_d > 1e-20 else 1e-20
            for sp, n_gibbs in gibbs_result.items():
                if sp in idx:
                    R_d[idx[sp]] += k_relax * (n_gibbs - max(self.N_d[idx[sp]], 0.0)) / V_d

        # Tar reactions in suspension：R10 同气泡；R11 = Corella × 局部催化剂密度
        rho_cat = self._catalyst_bulk_density()
        C_tar_mix_d = C_d[idx["TAR1"]] + C_d[idx["TAR2"]]
        r10d = rate_R10(T, C_tar_mix_d, C_d[idx["O2"]], self.P, self.fuel_type)
        r11d = rate_R11_suspension(T, C_tar_mix_d, rho_cat)
        tar_src_d = self._calc_tar_source({"R10": r10d, "R11": r11d})
        R_d += tar_src_d

        self.R_gas_d = R_d * self.V_d  # [mol/s]

    def _allocate_pyrolysis_products_elemental(
        self,
        nC: float,
        nH: float,
        nO: float,
    ) -> Dict[str, float]:
        """按元素守恒分配热解产物（主产物 CO/H2/CH4/Teer，必要时 CO2/H2O 兜底）。

        Chapter 4 建议流程：挥发分先按元素守恒释放，再进入后续反应网络。
        """
        nC = max(nC, 0.0)
        nH = max(nH, 0.0)
        nO = max(nO, 0.0)
        out = {"CO": 0.0, "H2": 0.0, "CH4": 0.0, "TAR1": 0.0, "CO2": 0.0, "H2O": 0.0}
        if nC <= 1e-16 and nH <= 1e-16 and nO <= 1e-16:
            return out

        # TAR1 的代理分子按 fuel_type 映射，读取其 C/H 原子数
        tar_map = get_tar_component_mapping(self.fuel_type)
        tar_surrogate = tar_map["TAR1"]
        c_tar, h_tar = TAR_SURROGATE_FORMULA[tar_surrogate]
        f_tar = float(np.clip(self.solid.pyrolysis_tar_carbon_frac, 0.0, 0.95))

        # 先给定焦油碳份额，再尽量用 O 生成 CO，剩余碳优先给 CH4
        n_tar = (f_tar * nC) / max(c_tar, 1e-12)
        n_co = min(nO, max(nC - c_tar * n_tar, 0.0))
        n_ch4 = max(nC - c_tar * n_tar - n_co, 0.0)

        # 氢守恒修正：若 H 不足，先减 CH4，再减 TAR1
        h_used = h_tar * n_tar + 4.0 * n_ch4
        if h_used > nH:
            deficit = h_used - nH
            d_ch4 = min(n_ch4, deficit / 4.0)
            n_ch4 -= d_ch4
            deficit -= 4.0 * d_ch4
            if deficit > 1e-16 and h_tar > 0:
                d_tar = min(n_tar, deficit / h_tar)
                n_tar -= d_tar

        # 重新核算剩余 C/O/H，若仅靠 4 主产物无法闭合，用 CO2/H2O 兜底
        c_used = c_tar * n_tar + n_co + n_ch4
        c_left = max(nC - c_used, 0.0)
        o_left = max(nO - n_co, 0.0)
        h_left = max(nH - (h_tar * n_tar + 4.0 * n_ch4), 0.0)

        if c_left > 1e-16 and o_left > 1e-16:
            # 优先 CO，再 CO2
            d_co = min(c_left, o_left)
            n_co += d_co
            c_left -= d_co
            o_left -= d_co
            d_co2 = min(c_left, o_left / 2.0)
            out["CO2"] += d_co2
            c_left -= d_co2
            o_left -= 2.0 * d_co2

        if o_left > 1e-16 and h_left > 1e-16:
            d_h2o = min(o_left, h_left / 2.0)
            out["H2O"] += d_h2o
            o_left -= d_h2o
            h_left -= 2.0 * d_h2o

        out["CO"] = max(n_co, 0.0)
        out["CH4"] = max(n_ch4, 0.0)
        out["TAR1"] = max(n_tar, 0.0)
        out["H2"] = max(h_left / 2.0, 0.0)
        return out

    def _calc_drying_pyrolysis_gas_source(
        self, tau_cell: float
    ) -> npt.NDArray[np.float64]:
        """Chapter 4 干燥-热解串联耦合源项（返回 mol/s, 写入悬浮相）。

        方法：
        1) 先用 Agarwal/FD 得到蒸发前沿与径向温度史；
        2) 仅对干壳层体积分数触发 DAEM（Eq. 4.12）；
        3) 将释放的水分与挥发分转为气相摩尔源项，并同步扣减固体质量。
        """
        src = np.zeros(N_GAS)
        if tau_cell <= 1e-6:
            return src

        # 干燥：Agarwal (Eq.4.1-4.4, 4.9)
        dry = solve_drying_CN(
            d_p=self.solid.d_p,
            T_bed=self.T,
            T_init=300.0,
            moisture_wt=self.solid.moisture_wt,
            t_total=max(tau_cell, 0.05),
            Nr=12,
            Nt=80,
            h_conv=None,
            u_rel=max(self.u_b - self.u_mf, 0.1),
            return_history=True,
        )
        x_dry = float(dry["X_dry"][-1])
        r_nodes = dry["r_nodes"]
        t_hist = dry["t"]
        T_rt = dry.get("T_history_rt")
        r_e = float(dry["r_evap"][-1])
        R0 = max(self.solid.d_p * 0.5, 1e-9)
        dry_shell_frac = float(np.clip(1.0 - (min(r_e, R0) / R0) ** 3, 0.0, 1.0))

        # 固体入量：新鲜进料 m_solid_zu 含湿基水；上游流入 m_solid_in 已在上方 cell 干燥过，
        # 不得再按全量 moisture_wt 重复蒸发（否则轴向多格会数倍放大 H2O）。
        m_zu_total = float(np.sum(np.maximum(self.m_solid_zu, 0.0)))
        m_in_total = float(np.sum(np.maximum(self.m_solid_in, 0.0)))
        if m_zu_total <= 1e-12 and m_in_total <= 1e-12:
            return src

        # 水分蒸发 -> H2O；仅对新鲜进料 [kg/s]
        m_water_release = m_zu_total * (self.solid.moisture_wt / 100.0) * x_dry
        src[GAS_SPECIES_INDEX["H2O"]] += m_water_release / (MOLECULAR_WEIGHT["H2O"] * 1e-3)

        # 固相扣减（水分）[kg/s]：水分仅来自 zu 流，按 zu 各粒径类比例扣减
        if m_water_release > 1e-18 and m_zu_total > 1e-12:
            for k in range(self.solid.n_size_classes):
                self.R_solid[k] += -m_water_release * (
                    self.m_solid_zu[k] / m_zu_total
                )

        # 热解：仅干壳层参与 DAEM
        if T_rt is not None:
            x_vm_raw = daem_conversion_radial(T_rt, t_hist, r_nodes)
        else:
            x_vm_raw = 0.0
        x_vm_eff = float(np.clip(x_vm_raw * dry_shell_frac, 0.0, 1.0))

        # daf 入料质量：新鲜 zu 按湿基去水；上游流入视为已无水分（干基固体）
        w_m = np.clip(self.solid.moisture_wt / 100.0, 0.0, 0.95)
        w_ash_dry = np.clip(self.solid.ash_dry_wt / 100.0, 0.0, 0.95)
        m_dry_total = m_zu_total * (1.0 - w_m) + m_in_total
        m_daf_in = m_dry_total * (1.0 - w_ash_dry)
        vm_daf = np.clip(self.solid.VM_daf / 100.0, 0.0, 0.95)
        # m_vm_release [kg/s] daf 挥发分释放率（入口 daf 流 × VM 份额 × DAEM 进度）
        m_vm_release = m_daf_in * vm_daf * x_vm_eff

        if m_vm_release <= 1e-12:
            return src

        # 由 dry -> daf 换算 C/H/O，做元素守恒分配
        c_dry = np.clip(self.solid.C_dry / 100.0, 0.0, 1.0)
        h_dry = np.clip(self.solid.H_dry / 100.0, 0.0, 1.0)
        o_dry = np.clip(self.solid.O_dry / 100.0, 0.0, 1.0)
        to_daf = 1.0 / max(1.0 - w_ash_dry, 1e-9)
        c_daf = c_dry * to_daf
        h_daf = h_dry * to_daf
        o_daf = o_dry * to_daf
        cho_sum = max(c_daf + h_daf + o_daf, 1e-12)
        c_eff = c_daf / cho_sum
        h_eff = h_daf / cho_sum
        o_eff = o_daf / cho_sum

        M_C = 12.011e-3
        M_H = 1.00794e-3
        M_O = 15.999e-3
        nC = m_vm_release * c_eff / M_C
        nH = m_vm_release * h_eff / M_H
        nO = m_vm_release * o_eff / M_O

        prod = self._allocate_pyrolysis_products_elemental(nC=nC, nH=nH, nO=nO)
        for sp, n_val in prod.items():
            if sp in GAS_SPECIES_INDEX and n_val > 0:
                src[GAS_SPECIES_INDEX[sp]] += n_val

        # 固相扣减（挥发分）[kg/s]
        for k in range(self.solid.n_size_classes):
            self.R_solid[k] += -m_vm_release * self.solid.mass_fractions[k]

        return src

    def calc_gibbs_correction(
        self, reaction_id: str, phase: str = "d", clamp_irreversible: bool = False
    ) -> float:
        """计算反应热力学驱动力 (1 − Q_p/K_eq)。

        TechSpec §5.5：供动力学速率修正使用。
        """
        y = self._mole_fractions(phase)
        return calc_gibbs_driving_force(
            reaction_id, self.T, self.P, y, GAS_SPECIES_INDEX, clamp_irreversible
        )

    def get_element_release(self, element: str) -> float:
        """统计该 cell 内可用元素量 [mol/s]（Hamel 两步：挥发 + 炭消耗 + 气相）。

        S：挥发硫（热解）+ 炭硫（R1/R2/R3 炭消耗释放）+ 气相 H2S/SO2/COS
        N：固体进料 + 气相 NH3/HCN/NO/N2

        TechSpec §5.5；Hamel 论文 §4.2、Eq. A.28
        """
        M_S, M_N = 32.065e-3, 14.007e-3  # [kg/mol]
        total = 0.0
        if element == "S":
            m_zu = float(np.sum(self.m_solid_zu))
            # 挥发硫：热解释放（Hamel Step 1）
            total += m_zu * self.solid.sulfur_fraction * self.solid.sulfur_volatile_frac / M_S
            # 炭硫：R1/R2/R3 炭消耗时释放（需 R_solid 已计算）
            char_consumption = max(0.0, -float(np.sum(self.R_solid)))
            total += char_consumption * self.solid.sulfur_fraction * (
                1.0 - self.solid.sulfur_volatile_frac
            ) / M_S
        elif element == "N":
            m_zu = float(np.sum(self.m_solid_zu))
            total += m_zu * self.solid.nitrogen_fraction / M_N
        for N_arr in (self.N_d, self.N_b):
            for j, sp in enumerate(GAS_SPECIES):
                total += max(N_arr[j], 0.0) * get_atom_count(sp, element)
        return max(total, 0.0)

    def calc_minor_species_gibbs(
        self,
        use_sulfur: bool = True,
        use_nitrogen: bool = True,
    ) -> Dict[str, float]:
        """调用 Gibbs 最小化，返回微量组分摩尔流率 {组分名: mol/s}。

        TechSpec §5.5：当 use_gibbs_minor 时替代 R9。
        """
        from src.thermodynamics.minor_species import (
            NITROGEN_CANDIDATES,
            SULFUR_CANDIDATES,
        )

        result: Dict[str, float] = {}
        if use_sulfur:
            b_S = self.get_element_release("S")
            if b_S > 1e-20:
                b_H = self.get_element_release("H")
                b_O = self.get_element_release("O")
                b_C = self.get_element_release("C")
                elements = {"S": b_S, "H": max(b_H, 1e-10), "O": max(b_O, 1e-10), "C": max(b_C, 1e-10)}
                try:
                    for k, v in solve_minor_species(
                        self.T, self.P, elements, SULFUR_CANDIDATES
                    ).items():
                        result[k] = v
                except Exception:
                    pass
        if use_nitrogen:
            b_N = self.get_element_release("N")
            if b_N > 1e-20:
                b_H = self.get_element_release("H")
                b_O = self.get_element_release("O")
                b_C = self.get_element_release("C")
                elements = {"N": b_N, "H": max(b_H, 1e-10), "O": max(b_O, 1e-10), "C": max(b_C, 1e-10)}
                try:
                    for k, v in solve_minor_species(
                        self.T, self.P, elements, NITROGEN_CANDIDATES
                    ).items():
                        result[k] = v
                except Exception:
                    pass
        return result

    def _catalyst_bulk_density(self) -> float:
        """悬浮相乳化团内催化剂/床料质量密度 rho_cat [kg/m³]。

        Corella et al. (1991) 中 R11 前置因子须乘以局部催化剂质量浓度：
        rho_cat ≈ rho_s · (1 - eps_mf) · f_cat。

        Source: Hamel (1999) §5.2.6, Table 5.4
        """
        f_cat = float(getattr(self.solid, "catalyst_solid_fraction", 1.0))
        return self.solid.rho_s * (1.0 - self.solid.eps_mf) * f_cat

    def _compute_char_conversion(self) -> float:
        """炭转化率 X [0,1)，用于 SPM 核径 d_core = d_p·(1−X)^(1/3)。

        稳态 cell：以炭质量流为基准，X = 1 − m_char_out / m_char_in。
        m_char_in 为进料与上游入固体中的炭；m_char_out 为当前出口固体中的炭。

        Source: 收缩颗粒/未反应核模型（SPM）；Hamel (1999) §6.3
        """
        nk = self.solid.n_size_classes
        m_char_in = 0.0
        m_char_out = 0.0
        for k in range(nk):
            m_in = max(self.m_solid_zu[k] + self.m_solid_in[k], 0.0)
            m_out = max(self.m_solid[k], 0.0)
            m_char_in += m_in * self.m_char_frac[k]
            m_char_out += m_out * self.m_char_frac[k]
        if m_char_in <= 1e-30:
            return 0.0
        ratio = m_char_out / m_char_in
        X = 1.0 - ratio
        return float(np.clip(X, 0.0, 1.0 - 1e-9))

    def _calc_char_surface_area_per_class(self) -> npt.NDArray[np.float64]:
        """由乳化相床层持料量估算各类炭粒外表面积 [m²]（非质量流率）。

        M_inventory ≈ rho_s * (1 - eps_mf) * V_d，其中 V_d 为乳化相体积；
        类 k 炭质量 M_k = M_inventory * mass_fractions[k] * m_char_frac[k]。

        Source: Hamel (1999) §5.1；两相床空隙模型
        """
        nk = self.solid.n_size_classes
        areas = np.zeros(nk, dtype=np.float64)
        V_d = max(self.V_d, 1e-30)
        M_inventory = self.solid.rho_s * (1.0 - self.solid.eps_mf) * V_d
        for k in range(nk):
            frac = float(self.solid.mass_fractions[k])
            char_frac = float(self.m_char_frac[k])
            M_k = M_inventory * frac * char_frac
            d_k = float(self.solid.d_p_classes[k])
            rho_s = float(self.solid.rho_s)
            if d_k > 0.0 and rho_s > 0.0 and M_k > 0.0:
                n_p = M_k / (rho_s * np.pi / 6.0 * d_k**3)
                areas[k] = n_p * np.pi * d_k**2
        return areas

    def _calc_char_surface_area(self) -> float:
        """cell 中活性炭总外表面积 [m²]（持料量模型）。"""
        return float(np.sum(self._calc_char_surface_area_per_class()))

    def _calc_tar_source(
        self,
        rates: Dict[TarReactionId, float],
    ) -> npt.NDArray[np.float64]:
        """tar 反应源项向量 [mol/(m³·s)]，写入 _work_tar_src 并返回。"""
        self._work_tar_src.fill(0.0)
        for rid, rate in rates.items():
            if rate <= 0:
                continue
            nu = get_lumped_tar_stoichiometry(rid, self.fuel_type)
            for species, coeff in nu.items():
                if species in GAS_SPECIES_INDEX:
                    self._work_tar_src[GAS_SPECIES_INDEX[species]] += coeff * rate
        return self._work_tar_src

    # ===================================================================
    # 4. calc_gas_balance — 气相摩尔守恒残差
    # ===================================================================

    def calc_gas_balance(self) -> npt.NDArray[np.float64]:
        """气相摩尔守恒残差，shape=(2*N_GAS,)，写入 _work_res 并返回视图。

        悬浮相(d): 0 = N_zu_d + N_rez_d + N_d_in + R_gas_d - N_d - N_ex
        气泡相(b): 0 = N_zu_b + N_rez_b + N_b_in + R_gas_b - N_b + N_ex

        Source: specs/01_conservation_equations.md §2.1; Eq.2-1/2-2
        """
        res_d = (
            self.N_zu_d + self.N_rez_d + self.N_d_in
            + self.R_gas_d - self.N_d - self.N_ex
        )
        res_b = (
            self.N_zu_b + self.N_rez_b + self.N_b_in
            + self.R_gas_b - self.N_b + self.N_ex
        )
        self._work_res[:N_GAS] = res_d
        self._work_res[N_GAS:2 * N_GAS] = res_b
        return self._work_res[:2 * N_GAS]

    # ===================================================================
    # 5. calc_solid_balance — 固相质量守恒残差
    # ===================================================================

    def calc_solid_balance(self) -> npt.NDArray[np.float64]:
        """固相质量守恒残差，shape=(n_size_classes,)。

        0 = m_solid_zu + m_solid_in + R_solid - m_solid

        (简化版：未实现粒径类间迁移 m_left/m_right)
        TODO: 实现完整的粒径类间迁移（收缩颗粒 -> 小粒径类转移）

        Source: specs/01_conservation_equations.md §2.2; Eq.2-3
        """
        return self.m_solid_zu + self.m_solid_in + self.R_solid - self.m_solid

    # ===================================================================
    # 6. calc_energy_balance — 全局能量守恒残差
    # ===================================================================

    def calc_energy_balance(self) -> float:
        """全局能量守恒残差 [W]。

        H_in = H(上游气, T_in_gas) + H(新鲜气 zu, T_zu_gas) + H(循环气 rez, T_rez_gas)
             + H(上游固, T_in_solid) + H(进料固, T_zu_solid)
        H_out = H(出口气, T) + H(出口固, T)；气相 h_j(T) 含生成焓；固相显焓自 T_REF。

        Source: specs/01_conservation_equations.md §2.3; Eq.2-4/2-5
        """
        H_gas_in = self._calc_gas_enthalpy_flow(self.N_b_in + self.N_d_in, self.T_in_gas)
        # 新鲜进料气与循环气分开：N_rez 来自炉顶高温，不可用 T_zu_gas（进料温度）
        H_gas_zu = self._calc_gas_enthalpy_flow(
            self.N_zu_b + self.N_zu_d,
            self.T_zu_gas,
        )
        H_gas_rez = self._calc_gas_enthalpy_flow(
            self.N_rez_b + self.N_rez_d,
            self.T_rez_gas,
        )
        H_gas_out = self._calc_gas_enthalpy_flow(self.N_b + self.N_d, self.T)

        H_solid_in = self._calc_solid_enthalpy_flow(self.m_solid_in, self.T_in_solid)
        H_solid_zu = self._calc_solid_enthalpy_flow(self.m_solid_zu, self.T_zu_solid)
        H_solid_out = self._calc_solid_enthalpy_flow(self.m_solid, self.T)

        return float(
            H_gas_in + H_gas_zu + H_gas_rez + H_solid_in + H_solid_zu
            - H_gas_out - H_solid_out - self.Q_wall
        )

    def _enthalpy_array_at_T(self, T: float) -> npt.NDArray[np.float64]:
        """各组分摩尔焓 h_j(T) [J/mol]，写入 _work_h_array 并返回。同 T 复用缓存。

        注意：每次调用 `_calc_gas_enthalpy_flow` 会立即 `np.dot` 消耗数值；
        勿在两次焓流计算之间长期保存 `_work_h_array` 的引用。
        """
        if abs(T - self._h_array_T_cache) < 1e-6:
            return self._work_h_array
        self._h_array_T_cache = T
        for j, sp in enumerate(GAS_SPECIES):
            try:
                self._work_h_array[j] = enthalpy_molar(sp, T)
            except NotImplementedError:
                self._work_h_array[j] = 0.0
        return self._work_h_array

    def _calc_gas_enthalpy_flow(
        self, N: npt.NDArray[np.float64], T: float,
    ) -> float:
        """气相焓流 [W] = sum_j N_j * h_j(T)，向量化实现。"""
        h_arr = self._enthalpy_array_at_T(T)
        return float(np.dot(N, h_arr))

    def _calc_solid_enthalpy_flow(
        self, m: npt.NDArray[np.float64], T: float,
    ) -> float:
        """固相焓流 [W]（char/灰/砂加权 Cp；自 T_REF 起的显焓，与气相 NASA 参考一致）。"""
        cp_avg = 0.3 * cp_char(T) + 0.3 * cp_ash(T) + 0.4 * cp_sand(T)
        return float(np.sum(m) * cp_avg * (T - T_REF))

    # ===================================================================
    # 完整残差向量
    # ===================================================================

    def residuals(self) -> npt.NDArray[np.float64]:
        """计算所有守恒方程残差并返回拼接向量。

        顺序: [gas_d(N_GAS), gas_b(N_GAS), solid(nk), energy(1)]
        写入预分配 _work_res，避免 concatenate 分配。
        """
        self.calc_hydrodynamics()
        self.calc_exchange()
        self.calc_reactions()

        self.calc_gas_balance()
        nk = self.solid.n_size_classes
        self._work_res[2 * N_GAS:2 * N_GAS + nk] = self.calc_solid_balance()
        self._work_res[-1] = self.calc_energy_balance()

        return self._work_res
