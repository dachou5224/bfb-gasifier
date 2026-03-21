"""Reactor 类：多 cell 串联扫描迭代求解器。

沿轴向从底部到顶部依次求解每个 cell：
  1. 设定 cell 的入口条件（上一 cell 的出口 / 底部进料）
  2. 调用 cell_solver 求解
  3. 将出口传递给下一 cell

支持外循环（recirculation）：将出口固体按比例返回底部 cell。

Source: docs/CLAUDE.md Phase 5.3
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np

from src.core.cell import Cell, CellGeometry, SolidProps
from src.core.composition import wet_to_dry_mole_fractions
from src.core.species import (
    GAS_SPECIES,
    GAS_SPECIES_INDEX,
    N_GAS,
    TarFuelType,
    configure_tar_components_by_fuel,
)
from src.core.constants import Rg
from src.core.feed_inlet import compute_gas_feeds_mol_s
from src.solvers.cell_solver import solve_cell


@dataclass
class ReactorConfig:
    """反应器配置。

    **气化剂进料**：若设置 **ER**（当量比），则根据 `feed_inlet.compute_gas_feeds_mol_s`
    自动计算 **O2_feed / H2O_feed / N2_feed**（覆盖下方默认值），并与燃料元素分析、
    **primary_agent**（空气+蒸汽 / 纯氧+蒸汽等）对齐。若 **ER 为 None**，则使用手动给定的
    三股摩尔流率（兼容旧脚本与基准测试）。
    """
    n_cells: int = 10             # [-]    轴向 cell 数
    H_bed: float = 5.0            # [m]    流化床高度
    H_freeboard: float = 9.5      # [m]    自由板区高度
    D_bed: float = 0.6            # [m]    床层直径

    P: float = 2_500_000.0        # [Pa]   操作压力
    T_inlet: float = 300.0        # [K]    进气温度

    fuel_type: TarFuelType = "coal"

    # 固体属性
    rho_s: float = 1400.0         # [kg/m³]
    d_p: float = 0.001            # [m]
    phi_s: float = 0.86
    eps_mf: float = 0.45
    # 颗粒龄期 / 离散直径类数（与 Hamel 多类离散一致）→ SolidProps.n_size_classes；当前默认 1
    n_age_classes: int = 1

    # 进料
    fuel_feed: float = 0.938      # [kg/s] (3377 kg/h)

    # 燃料分析（干基元素用于化学计量氧与 ER）
    moisture_wt: float = 16.9     # [wt%]
    ash_dry_wt: float = 11.41     # [wt%]
    VM_daf: float = 53.42         # [wt%]
    C_dry: float = 61.5           # [wt%]
    H_dry: float = 4.1            # [wt%]
    O_dry: float = 21.8           # [wt%]
    S_dry: float = 0.0            # [wt%] 干基硫（CASE_LU 等需填入）

    # 当量比与气化剂类型（若 ER 非 None，则自动填充 O2/H2O/N2_feed）
    ER: float | None = None       # [-] 当量比；None 表示使用手动 O2_feed 等
    primary_agent: str = "air_steam"  # air_steam | o2_steam | air（见 feed_inlet）
    steam_to_o2_molar: float = 0.8    # [-] H2O/O2 摩尔比（蒸汽量）
    o2_steam_n2_frac_of_o2: float = 0.01  # [-] o2_steam 时 N2/O2 微量载气比

    # 气化剂摩尔流率 [mol/s]；当 ER 给定时常由 `compute_gas_feeds_mol_s` 写入
    O2_feed: float = 0.5
    H2O_feed: float = 0.3
    N2_feed: float = 0.02

    recirculation_frac: float = 0.1  # [-] 固体循环比

    def __post_init__(self) -> None:
        if self.ER is not None:
            o2, h2o, n2 = compute_gas_feeds_mol_s(
                fuel_feed_kg_s=self.fuel_feed,
                moisture_wt=self.moisture_wt,
                C_dry=self.C_dry,
                H_dry=self.H_dry,
                O_dry=self.O_dry,
                ER=float(self.ER),
                primary_agent=self.primary_agent,
                S_dry=self.S_dry,
                steam_to_o2_molar=self.steam_to_o2_molar,
                o2_steam_n2_frac_of_o2=self.o2_steam_n2_frac_of_o2,
            )
            self.O2_feed = o2
            self.H2O_feed = h2o
            self.N2_feed = n2


class Reactor:
    """BFB 反应器：多 cell 串联求解。

    Source: docs/CLAUDE.md Phase 5.3; Hamel & Krumm (2001)
    """

    def __init__(self, config: ReactorConfig) -> None:
        self.config = config
        self.cells: List[Cell] = []
        self._build_cells()

        configure_tar_components_by_fuel(config.fuel_type)

    def _build_cells(self) -> None:
        """创建并初始化所有 cell。"""
        cfg = self.config
        dh = cfg.H_bed / cfg.n_cells

        nk = max(1, int(cfg.n_age_classes))
        d_p_classes = np.full(nk, float(cfg.d_p), dtype=np.float64)
        mass_fractions = np.ones(nk, dtype=np.float64) / float(nk)

        solid = SolidProps(
            rho_s=cfg.rho_s,
            d_p=cfg.d_p,
            phi_s=cfg.phi_s,
            eps_mf=cfg.eps_mf,
            n_size_classes=nk,
            d_p_classes=d_p_classes,
            mass_fractions=mass_fractions,
            moisture_wt=cfg.moisture_wt,
            ash_dry_wt=cfg.ash_dry_wt,
            VM_daf=cfg.VM_daf,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            O_dry=cfg.O_dry,
        )

        self.cells = []
        for i in range(cfg.n_cells):
            geo = CellGeometry(
                D_bed=cfg.D_bed,
                dh=dh,
                h_center=(i + 0.5) * dh,
            )
            cell = Cell(geo=geo, solid=solid, fuel_type=cfg.fuel_type)
            cell.P = cfg.P
            cell.T = 900.0 + 273.15  # 初始猜测 [K]
            # m_char_frac 由 Cell.__init__ 据 SolidProps 计算；与 ReactorConfig 一致时应数值相同
            expected_fc = max(
                (
                    1.0
                    - cfg.ash_dry_wt / 100.0
                    - cfg.VM_daf / 100.0 * (1.0 - cfg.ash_dry_wt / 100.0)
                )
                * (1.0 - cfg.moisture_wt / 100.0),
                0.05,
            )
            assert abs(float(cell.m_char_frac[0]) - expected_fc) < 0.01, (
                f"m_char_frac 与工业分析不一致: {cell.m_char_frac[0]:.4f} vs 期望 {expected_fc:.4f}"
            )

            # 初始摩尔流率猜测
            total_gas = cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed
            for j, sp in enumerate(GAS_SPECIES):
                frac = 0.01  # 默认小量
                if sp == "O2":
                    frac = cfg.O2_feed / total_gas
                elif sp == "H2O":
                    frac = cfg.H2O_feed / total_gas
                elif sp == "N2":
                    frac = cfg.N2_feed / total_gas
                cell.N_b[j] = frac * total_gas * 0.3
                cell.N_d[j] = frac * total_gas * 0.7

            cell.m_solid[:] = cfg.fuel_feed * 0.5  # 初始猜测
            cell._vm_cache_valid = False  # 由 Reactor.solve 每外迭代 Vorabrechnung 置真
            self.cells.append(cell)

    def _set_bottom_cell_feeds(self) -> None:
        """设定底部 cell (i=0) 的进料和进气。"""
        cfg = self.config
        cell = self.cells[0]

        idx = GAS_SPECIES_INDEX
        cell.N_zu_d[idx["O2"]] = cfg.O2_feed * 0.3
        cell.N_zu_d[idx["H2O"]] = cfg.H2O_feed * 0.3
        cell.N_zu_d[idx["N2"]] = cfg.N2_feed * 0.3
        cell.N_zu_b[idx["O2"]] = cfg.O2_feed * 0.7
        cell.N_zu_b[idx["H2O"]] = cfg.H2O_feed * 0.7
        cell.N_zu_b[idx["N2"]] = cfg.N2_feed * 0.7

        cell.m_solid_zu[:] = cfg.fuel_feed

        # 进料/入口流股温度：气体为工况进气；固体为冷态进料（与 validation T_inlet≈293 K 一致）
        cell.T_zu_gas = cfg.T_inlet
        cell.T_zu_solid = 293.15
        cell.T_in_gas = cfg.T_inlet
        cell.T_in_solid = cfg.T_inlet

    def _propagate_upstream(self, i: int) -> None:
        """将 cell[i-1] 的出口设为 cell[i] 的入口。

        注意：cell 编号从底部 (0) 到顶部 (n-1)，
        气体从底向上流动，固体从上向下流动。
        """
        if i == 0:
            return
        prev = self.cells[i - 1]
        curr = self.cells[i]

        # 气体：从下方 cell 流入
        curr.N_b_in[:] = prev.N_b
        curr.N_d_in[:] = prev.N_d
        curr.T_in_gas = prev.T

        # 固体：从上方 cell 流入（反向）
        cfg = self.config
        if i < len(self.cells) - 1:
            above = self.cells[i + 1]
            curr.m_solid_in[:] = above.m_solid
            curr.T_in_solid = above.T
        else:
            curr.m_solid_in[:] = 0.0
            curr.T_in_solid = cfg.T_inlet

    def _apply_all_bc_for_nr(self) -> None:
        """全局残差用：底格进料 + 轴向传播 + 循环返料（与 Gauss–Seidel 物理一致）。"""
        cfg = self.config
        self._set_bottom_cell_feeds()
        self.cells[-1].m_solid_in[:] = 0.0
        for i in range(len(self.cells)):
            self._propagate_upstream(i)
        if cfg.recirculation_frac > 0:
            top = self.cells[-1]
            bot = self.cells[0]
            bot.N_rez_d[:] = top.N_d * cfg.recirculation_frac
            bot.T_rez_gas = top.T
            bot.m_solid_zu[:] = cfg.fuel_feed + top.m_solid * cfg.recirculation_frac

    def solve(
        self,
        max_global_iter: int = 20,
        tol_global: float = 1e-4,
        solver: str = "gauss_seidel",
    ) -> dict:
        """执行多 cell 求解。

        Parameters
        ----------
        solver : str
            ``"gauss_seidel"`` — 逐格 ``fsolve`` 扫描（默认，兼容旧行为）；
            ``"global_nr"`` — Hamel (1999) §2.3 全局阻尼 Newton–Raphson。

        Returns
        -------
        dict with keys:
          'converged'   : bool  外迭代（温度/Vorabrechnung 固定点）是否满足判据（与历史兼容）
          'converged_outer' : bool  同上，显式命名
          'converged_inner_nr' : bool  最后一次内层 global NR 是否满足 RMS(F̂) 容差
          'rms_scaled_final' : float | None  最后一次内层结束时的 RMS(F̂)
          'n_iter'      : int
          'T_profile'   : list[float]  各 cell 温度
          'exit_gas'     : dict          顶部出口气体组成（湿基）
          'exit_gas_dry' : dict          同上，干基（扣除 H2O）
          'carbon_conv'  : float         碳转化率估算

        Source: docs/CLAUDE.md Phase 5.3
        """
        if solver == "global_nr":
            return self._solve_global_nr(max_iter=max_global_iter, tol=tol_global)
        return self._solve_gauss_seidel(
            max_global_iter=max_global_iter, tol_global=tol_global
        )

    def _solve_gauss_seidel(self, max_global_iter: int, tol_global: float) -> dict:
        """逐格 fsolve + 外迭代（原 ``solve`` 实现）。"""
        cfg = self.config

        for g_iter in range(max_global_iter):
            T_old = [c.T for c in self.cells]

            self._set_bottom_cell_feeds()

            # 循环气焓温：与炉顶出口气一致（上一外迭代末的顶格 T；首迭代用初始 guess）
            if cfg.recirculation_frac > 0:
                self.cells[0].T_rez_gas = self.cells[-1].T

            # 固体入口（顶部 cell 进料）
            self.cells[-1].m_solid_in[:] = 0.0

            # Vorabrechnung：用当前 T、流场估计干燥/热解源项，本外迭代内 fsolve 固定（Hamel Module 2）
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
                c = self.cells[i]
                c.calc_hydrodynamics()
                tau = c.geo.dh / max(c.u_mf, 1e-3)
                c.compute_vorabrechnung(tau)

            # 正扫描（底 -> 顶）
            max_res = 0.0
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
                result = solve_cell(self.cells[i])
                max_res = max(max_res, result.get("residual", 0.0))
                if not result["converged"]:
                    pass  # 允许继续迭代

            # 循环：返回部分固体到底部
            if cfg.recirculation_frac > 0:
                top_solid = self.cells[-1].m_solid.copy()
                self.cells[0].N_rez_d[:] = (
                    self.cells[-1].N_d * cfg.recirculation_frac
                )
                self.cells[0].m_solid_zu[:] = (
                    cfg.fuel_feed + top_solid * cfg.recirculation_frac
                )

            # 检查全局收敛：温度变化或残差足够小
            dT = max(abs(c.T - T_old[j]) for j, c in enumerate(self.cells))
            if dT < tol_global:
                break
            if max_res < tol_global * 50:
                break

        # 提取结果
        top = self.cells[-1]
        y_exit = top._mole_fractions("d")
        exit_gas = {sp: float(y_exit[j]) for j, sp in enumerate(GAS_SPECIES)}

        C_in = cfg.fuel_feed * (cfg.C_dry / 100.0) * (1.0 - cfg.moisture_wt / 100.0)
        M_C = 12.011e-3
        C_out_mol = (
            top.N_b[GAS_SPECIES_INDEX["CO"]]
            + top.N_b[GAS_SPECIES_INDEX["CO2"]]
            + top.N_b[GAS_SPECIES_INDEX["CH4"]]
            + top.N_d[GAS_SPECIES_INDEX["CO"]]
            + top.N_d[GAS_SPECIES_INDEX["CO2"]]
            + top.N_d[GAS_SPECIES_INDEX["CH4"]]
        )
        C_out_kg = C_out_mol * M_C
        carbon_conv = C_out_kg / max(C_in, 1e-30)

        return {
            "converged": dT < tol_global if max_global_iter > 0 else False,
            "n_iter": g_iter + 1,
            "T_profile": [c.T for c in self.cells],
            "exit_gas": exit_gas,
            "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas),
            "carbon_conv": min(carbon_conv, 1.0),
        }

    def _solve_global_nr(self, max_iter: int, tol: float) -> dict:
        """全局阻尼 NR + 双层 Vorabrechnung（Hamel 1999 Bild 2.2）。

        Level 1 — 一次 bootstrap：轴向 T 估计 → 水力学 → 化学一致 x₀。
        Level 2 — 外迭代：按当前 T 重算 DAEM 缓存 → 内层 global NR（固定缓存）→ 外收敛判据。
        """
        from src.solvers.global_nr_solver import solve_global_nr
        from src.solvers.vorabrechnung import (
            estimate_axial_T_profile,
            generate_initial_x0,
        )

        cfg = self.config

        # ── Level 1: Bootstrap x₀ ──
        T_profile_est = estimate_axial_T_profile(
            n_cells=cfg.n_cells,
            T_inlet=cfg.T_inlet,
            O2_feed=cfg.O2_feed,
            H2O_feed=cfg.H2O_feed,
            N2_feed=cfg.N2_feed,
            fuel_feed_kg_s=cfg.fuel_feed,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            moisture_wt=cfg.moisture_wt,
            P=cfg.P,
        )

        for i, cell in enumerate(self.cells):
            cell.T = float(T_profile_est[i])

        self._set_bottom_cell_feeds()
        for i in range(len(self.cells)):
            self._propagate_upstream(i)
        for cell in self.cells:
            cell.calc_hydrodynamics()

        generate_initial_x0(
            cells=self.cells,
            O2_feed=cfg.O2_feed,
            H2O_feed=cfg.H2O_feed,
            N2_feed=cfg.N2_feed,
            fuel_feed_kg_s=cfg.fuel_feed,
            C_dry=cfg.C_dry,
            H_dry=cfg.H_dry,
            O_dry=cfg.O_dry,
            moisture_wt=cfg.moisture_wt,
            ash_dry_wt=cfg.ash_dry_wt,
            VM_daf=cfg.VM_daf,
            T_profile=T_profile_est,
            fuel_type=cfg.fuel_type,
        )

        # ── Level 1.5: Gauss–Seidel 暖启动（逐格 fsolve，替代过粗的 Vorabrechnung x₀）──
        N_GS_WARMUP = 3
        for _gs in range(N_GS_WARMUP):
            self._set_bottom_cell_feeds()
            if cfg.recirculation_frac > 0:
                self.cells[0].T_rez_gas = self.cells[-1].T
            self.cells[-1].m_solid_in[:] = 0.0
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
                c = self.cells[i]
                c.calc_hydrodynamics()
                tau = c.geo.dh / max(c.u_mf, 1e-3)
                c.compute_vorabrechnung(tau)
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
                solve_cell(self.cells[i])
            if cfg.recirculation_frac > 0:
                top = self.cells[-1]
                self.cells[0].N_rez_d[:] = top.N_d * cfg.recirculation_frac
                self.cells[0].m_solid_zu[:] = (
                    cfg.fuel_feed + top.m_solid * cfg.recirculation_frac
                )

        # ── Level 2: 外迭代 + 内层 NR ──
        outer_max = max(3, max_iter // 5)
        tol_outer = max(20.0, float(tol) * 5.0)
        # 与全局 tol_global 对齐的 RMS(F̂) 内层容差（工程量级 ~1–5%）
        tol_rms = max(0.01, min(0.05, 0.01 * float(tol)))

        ref_gas = float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed)
        ref_solid = cfg.fuel_feed
        ref_energy = cfg.fuel_feed * 20e6

        all_nr_history: list = []
        outer_converged = False
        last_nr_result: dict = {}
        last_inner_converged = False

        for _outer_iter in range(outer_max):
            T_vorab = np.array([c.T for c in self.cells], dtype=np.float64)

            self._set_bottom_cell_feeds()
            if cfg.recirculation_frac > 0:
                self.cells[0].T_rez_gas = self.cells[-1].T
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
            for cell in self.cells:
                cell.calc_hydrodynamics()
                tau = cell.geo.dh / max(cell.u_mf, 1e-3)
                cell.compute_vorabrechnung(tau)

            def apply_bc_closure() -> None:
                self._apply_all_bc_for_nr()

            inner_max = max(10, max_iter // max(outer_max, 1))
            nr_result = solve_global_nr(
                cells=self.cells,
                apply_bc_fn=apply_bc_closure,
                ref_gas_mol_s=ref_gas,
                ref_solid_kg_s=ref_solid,
                ref_energy_W=ref_energy,
                max_iter=inner_max,
                tol_rms=tol_rms,
                verbose=False,
            )
            last_nr_result = nr_result
            last_inner_converged = bool(nr_result.get("converged", False))
            all_nr_history.extend(nr_result.get("norm_history", []))

            T_new = np.array([c.T for c in self.cells], dtype=np.float64)
            dT_outer = float(np.max(np.abs(T_new - T_vorab)))

            if dT_outer < tol_outer and nr_result.get("converged", False):
                outer_converged = True
                break
            if dT_outer < tol_outer * 0.5:
                outer_converged = True
                break

        last_nr_result = dict(last_nr_result)
        # converged：外迭代（与旧版兼容）；内层 NR 是否收敛单独给出
        last_nr_result["converged_outer"] = outer_converged
        last_nr_result["converged_inner_nr"] = last_inner_converged
        last_nr_result["converged"] = outer_converged
        last_nr_result["norm_history"] = all_nr_history
        last_nr_result["n_iter"] = len(all_nr_history)
        last_nr_result["residual"] = last_nr_result.get(
            "residual",
            float(
                np.linalg.norm(
                    np.concatenate([c.residuals() for c in self.cells])
                )
            ),
        )
        return self._finalize_global_nr_result(last_nr_result)

    def _finalize_global_nr_result(self, nr_result: dict) -> dict:
        """组装 global_nr 返回 dict（与 Gauss–Seidel 键对齐 + NR 诊断）。"""
        cfg = self.config
        top = self.cells[-1]
        y_exit = top._mole_fractions("d")
        exit_gas = {sp: float(y_exit[j]) for j, sp in enumerate(GAS_SPECIES)}

        C_in = cfg.fuel_feed * (cfg.C_dry / 100.0) * (1.0 - cfg.moisture_wt / 100.0)
        M_C = 12.011e-3
        C_out_mol = (
            top.N_b[GAS_SPECIES_INDEX["CO"]]
            + top.N_b[GAS_SPECIES_INDEX["CO2"]]
            + top.N_b[GAS_SPECIES_INDEX["CH4"]]
            + top.N_d[GAS_SPECIES_INDEX["CO"]]
            + top.N_d[GAS_SPECIES_INDEX["CO2"]]
            + top.N_d[GAS_SPECIES_INDEX["CH4"]]
        )
        C_out_kg = C_out_mol * M_C
        carbon_conv = C_out_kg / max(C_in, 1e-30)

        conv_outer = bool(nr_result.get("converged_outer", nr_result["converged"]))
        conv_inner = bool(nr_result.get("converged_inner_nr", False))
        rms_sf = nr_result.get("rms_scaled_final")

        return {
            "converged": nr_result["converged"],
            "converged_outer": conv_outer,
            "converged_inner_nr": conv_inner,
            "converged_fully": conv_outer and conv_inner,
            "rms_scaled_final": rms_sf,
            "n_iter": nr_result["n_iter"],
            "residual": nr_result["residual"],
            "norm_history": nr_result.get("norm_history", []),
            "T_profile": [c.T for c in self.cells],
            "exit_gas": exit_gas,
            "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas),
            "carbon_conv": min(carbon_conv, 1.0),
        }
