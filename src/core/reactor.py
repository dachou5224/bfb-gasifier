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

from src.core.cell import (
    Cell,
    CellGeometry,
    SolidProps,
    S_CHAR,
    S_VM,
    S_MOISTURE,
    S_ASH,
    N_SOLID_COMP,
)
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
    """反应器配置。"""
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
    n_age_classes: int = 1

    # 进料
    fuel_feed: float = 0.938      # [kg/s] (3377 kg/h)

    # 燃料分析
    moisture_wt: float = 16.9     # [wt%]
    ash_dry_wt: float = 11.41     # [wt%]
    VM_daf: float = 53.42         # [wt%]
    C_dry: float = 61.5           # [wt%]
    H_dry: float = 4.1            # [wt%]
    O_dry: float = 21.8           # [wt%]
    S_dry: float = 0.0            # [wt%]
    HHV_dry: float = 22.0         # [MJ/kg]
    u0_target: float | None = None
    heat_loss_frac: float = 0.1
    nitrogen_fraction: float = 0.0

    # 气化剂
    ER: float | None = None
    primary_agent: str = "air_steam"
    steam_to_o2_molar: float = 0.8
    o2_steam_n2_frac_of_o2: float = 0.01

    O2_feed: float = 0.5
    H2O_feed: float = 0.3
    N2_feed: float = 0.02

    recirculation_frac: float = 0.1

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
    """BFB 反应器：多 cell 串联求解。"""

    def __init__(self, config: ReactorConfig) -> None:
        self.config = config
        self.cells: List[Cell] = []
        self._build_cells()
        configure_tar_components_by_fuel(config.fuel_type)

    def _build_cells(self) -> None:
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
            HHV_dry_MJ_kg=cfg.HHV_dry,
            nitrogen_fraction=cfg.nitrogen_fraction,
        )

        for i in range(cfg.n_cells):
            geo = CellGeometry(D_bed=cfg.D_bed, dh=dh, h_center=(i + 0.5) * dh)
            cell = Cell(geo=geo, solid=solid, fuel_type=cfg.fuel_type)
            cell.P = cfg.P
            cell.heat_loss_frac = cfg.heat_loss_frac
            cell.u0_target = cfg.u0_target
            self.cells.append(cell)

    def _set_bottom_cell_feeds(self) -> None:
        cfg = self.config
        cell = self.cells[0]
        idx = GAS_SPECIES_INDEX
        cell.N_zu_d.fill(0.0)
        cell.N_zu_b.fill(0.0)
        cell.m_solid_zu.fill(0.0)

        cell.N_zu_d[idx["O2"]] = cfg.O2_feed * 0.7  # 提高悬浮相分配
        cell.N_zu_d[idx["H2O"]] = cfg.H2O_feed * 0.7
        cell.N_zu_d[idx["N2"]] = cfg.N2_feed * 0.7
        cell.N_zu_b[idx["O2"]] = cfg.O2_feed * 0.3
        cell.N_zu_b[idx["H2O"]] = cfg.H2O_feed * 0.3
        cell.N_zu_b[idx["N2"]] = cfg.N2_feed * 0.3

        w_m = cfg.moisture_wt / 100.0
        w_ash_dry = cfg.ash_dry_wt / 100.0
        w_vm_dry = (cfg.VM_daf / 100.0) * (1.0 - w_ash_dry)
        w_char_dry = 1.0 - w_ash_dry - w_vm_dry
        nk = cfg.n_age_classes
        cell.m_solid_zu[:, S_CHAR] = cfg.fuel_feed * (1.0 - w_m) * w_char_dry / nk
        cell.m_solid_zu[:, S_VM] = cfg.fuel_feed * (1.0 - w_m) * w_vm_dry / nk
        cell.m_solid_zu[:, S_MOISTURE] = cfg.fuel_feed * w_m / nk
        cell.m_solid_zu[:, S_ASH] = cfg.fuel_feed * (1.0 - w_m) * w_ash_dry / nk

        cell.T_zu_gas = cfg.T_inlet
        cell.T_zu_solid = 293.15
        cell.T_in_gas = cfg.T_inlet
        cell.T_in_solid = cfg.T_inlet

    def _propagate_upstream(self, i: int) -> None:
        if i == 0: return
        prev, curr = self.cells[i - 1], self.cells[i]
        curr.N_b_in[:] = prev.N_b
        curr.N_d_in[:] = prev.N_d
        curr.T_in_gas = prev.T
        if i < len(self.cells) - 1:
            above = self.cells[i + 1]
            curr.m_solid_in[:] = above.m_solid
            curr.T_in_solid = above.T
        else:
            curr.m_solid_in[:] = 0.0
            curr.T_in_solid = self.config.T_inlet

    def _apply_all_bc_for_nr(self) -> None:
        cfg = self.config
        self._set_bottom_cell_feeds()
        self.cells[-1].m_solid_in[:] = 0.0
        for i in range(len(self.cells)):
            self._propagate_upstream(i)
        if cfg.recirculation_frac > 0:
            top, bot = self.cells[-1], self.cells[0]
            bot.N_rez_d[:] = top.N_d * cfg.recirculation_frac
            bot.T_rez_gas = top.T
            bot.m_solid_zu = bot.m_solid_zu + top.m_solid * cfg.recirculation_frac

    def solve(self, max_global_iter: int = 20, tol_global: float = 1e-4, solver: str = "gauss_seidel") -> dict:
        if solver == "global_nr":
            return self._solve_global_nr(max_iter=max_global_iter, tol=tol_global)
        return self._solve_gauss_seidel(max_global_iter=max_global_iter, tol_global=tol_global)

    def _solve_gauss_seidel(self, max_global_iter: int, tol_global: float) -> dict:
        cfg = self.config
        for g_iter in range(max_global_iter):
            T_old = [c.T for c in self.cells]
            max_res = 0.0
            self._set_bottom_cell_feeds()
            if cfg.recirculation_frac > 0: self.cells[0].T_rez_gas = self.cells[-1].T
            
            for i in range(len(self.cells)):
                self._propagate_upstream(i)
                cell = self.cells[i]
                old_T = cell.T
                
                # ── 关键改进：反应项预叠加（引入松弛因子） ──
                cell.calc_hydrodynamics()
                cell.compute_vorabrechnung(cell.geo.dh/max(cell.u_mf, 1e-3))
                cell.calc_reactions()
                # 仅叠加 20% 的变化，防止数值振荡和过早枯竭
                cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + 0.2 * cell.R_gas_d, 0.0)
                cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + 0.2 * cell.R_gas_b, 0.0)
                cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + 0.2 * cell.R_solid, 0.0)

                # 求解当前 cell 的守恒方程
                res_cell = solve_cell(cell)
                cell.T = old_T + 0.4 * (cell.T - old_T)
                max_res = max(max_res, res_cell.get("residual", 0.0))
                
            if cfg.recirculation_frac > 0:
                top, bot = self.cells[-1], self.cells[0]
                bot.N_rez_d[:] = top.N_d * cfg.recirculation_frac
                bot.m_solid_zu = bot.m_solid_zu + top.m_solid * cfg.recirculation_frac
            
            dT = max(abs(self.cells[j].T - T_old[j]) for j in range(len(self.cells)))
            if g_iter >= 10 and dT < tol_global: break
        
        top = self.cells[-1]
        y_exit = top._mole_fractions("d")
        exit_gas = {sp: float(y_exit[j]) for j, sp in enumerate(GAS_SPECIES)}
        
        # 碳转化率：基于总炭收支
        bot = self.cells[0]
        m_char_in = float(np.sum(bot.m_solid_zu[:, S_CHAR] + bot.m_solid_in[:, S_CHAR]))
        m_char_out = float(np.sum(top.m_solid[:, S_CHAR]))
        carbon_conv = 1.0 - (m_char_out / max(m_char_in, 1e-12))
        
        return {"converged": g_iter < max_global_iter - 1, "n_iter": g_iter + 1, "T_profile": [c.T for c in self.cells], "exit_gas": exit_gas, "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas), "carbon_conv": float(np.clip(carbon_conv, 0.0, 1.0))}

    def _solve_global_nr(self, max_iter: int, tol: float) -> dict:
        from src.solvers.global_nr_solver import solve_global_nr
        from src.solvers.vorabrechnung import estimate_axial_T_profile, generate_initial_x0
        cfg = self.config
        T_est = estimate_axial_T_profile(n_cells=cfg.n_cells, T_inlet=cfg.T_inlet, O2_feed=cfg.O2_feed, H2O_feed=cfg.H2O_feed, N2_feed=cfg.N2_feed, fuel_feed_kg_s=cfg.fuel_feed, C_dry=cfg.C_dry, H_dry=cfg.H_dry, moisture_wt=cfg.moisture_wt, P=cfg.P)
        for i, cell in enumerate(self.cells): cell.T = float(T_est[i])
        self._set_bottom_cell_feeds()
        for i in range(len(self.cells)): self._propagate_upstream(i)
        for cell in self.cells: cell.calc_hydrodynamics()
        generate_initial_x0(cells=self.cells, O2_feed=cfg.O2_feed, H2O_feed=cfg.H2O_feed, N2_feed=cfg.N2_feed, fuel_feed_kg_s=cfg.fuel_feed, C_dry=cfg.C_dry, H_dry=cfg.H_dry, O_dry=cfg.O_dry, moisture_wt=cfg.moisture_wt, ash_dry_wt=cfg.ash_dry_wt, VM_daf=cfg.VM_daf, T_profile=T_est, fuel_type=cfg.fuel_type)
        for _ in range(5): self._solve_gauss_seidel(max_global_iter=1, tol_global=100.0)
        outer_max = max(3, max_iter // 5)
        all_nr_history = []
        outer_converged = False
        last_nr_result = {}
        for _outer in range(outer_max):
            T_vorab = np.array([c.T for c in self.cells])
            self._apply_all_bc_for_nr()
            for c in self.cells:
                c.calc_hydrodynamics()
                c.compute_vorabrechnung(c.geo.dh / max(c.u_mf, 1e-3))
            nr_result = solve_global_nr(cells=self.cells, apply_bc_fn=self._apply_all_bc_for_nr, ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), ref_solid_kg_s=cfg.fuel_feed, ref_energy_W=cfg.fuel_feed * 20e6, max_iter=max(10, max_iter // outer_max), tol_rms=0.01)
            last_nr_result = nr_result
            all_nr_history.extend(nr_result.get("norm_history", []))
            dT_outer = np.max(np.abs(np.array([c.T for c in self.cells]) - T_vorab))
            if dT_outer < tol * 5.0 and nr_result.get("converged"):
                outer_converged = True
                break
        return self._finalize_global_nr_result(dict(last_nr_result, converged_outer=outer_converged, converged_inner_nr=bool(last_nr_result.get("converged")), converged=outer_converged, norm_history=all_nr_history, n_iter=len(all_nr_history)))

    def _finalize_global_nr_result(self, nr_result: dict) -> dict:
        cfg, top = self.config, self.cells[-1]
        y_exit = top._mole_fractions("d")
        exit_gas = {sp: float(y_exit[j]) for j, sp in enumerate(GAS_SPECIES)}
        bot = self.cells[0]
        m_char_in = float(np.sum(bot.m_solid_zu[:, S_CHAR] + bot.m_solid_in[:, S_CHAR]))
        m_char_out = float(np.sum(top.m_solid[:, S_CHAR]))
        carbon_conv = 1.0 - (m_char_out / max(m_char_in, 1e-12))
        return {"converged": nr_result["converged"], "converged_outer": nr_result.get("converged_outer"), "converged_inner_nr": nr_result.get("converged_inner_nr"), "converged_fully": nr_result.get("converged_outer") and nr_result.get("converged_inner_nr"), "rms_scaled_final": nr_result.get("rms_scaled_final"), "n_iter": nr_result["n_iter"], "residual": nr_result.get("residual", 0.0), "norm_history": nr_result.get("norm_history", []), "T_profile": [c.T for c in self.cells], "exit_gas": exit_gas, "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas), "carbon_conv": float(np.clip(carbon_conv, 0.0, 1.0))}
