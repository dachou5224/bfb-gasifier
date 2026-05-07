from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from src.core.reactor import Reactor

def solve_gauss_seidel_reactor(
    reactor: Reactor,
    max_iter: int,
    tol: float,
    verbose: bool = False,
) -> dict:
    """Legacy Gauss-Seidel solver for BFB reactor."""
    from src.solvers.cell_solver import solve_cell
    from src.core.species import GAS_SPECIES
    from src.core.composition import wet_to_dry_mole_fractions
    from src.core.cell import S_CHAR
    
    cfg = reactor.config
    cells = reactor.cells
    g_iter = 0
    
    for g_iter in range(max_iter):
        T_old = [c.T for c in cells]
        max_res = 0.0
        
        # 边界条件初始化
        reactor._set_bottom_cell_feeds()
        if cfg.recirculation_frac > 0:
            cells[0].T_rez_gas = cells[-1].T
            
        # 扫描迭代
        for i in range(len(cells)):
            reactor._propagate_upstream(i)
            cell = cells[i]
            old_T = cell.T
            
            # 物理准备与预更新 (符合 Hamel 原始 Abgleich 精神)
            cell.calc_hydrodynamics()
            cell.compute_vorabrechnung(cell.geo.dh / max(cell.u_mf, 1e-3))
            cell.calc_reactions()
            
            # 反应项预叠加（松弛因子 0.2）
            cell.N_d[:] = np.maximum(cell.N_d_in + cell.N_zu_d + 0.2 * cell.R_gas_d, 0.0)
            cell.N_b[:] = np.maximum(cell.N_b_in + cell.N_zu_b + 0.2 * cell.R_gas_b, 0.0)
            cell.m_solid[:] = np.maximum(cell.m_solid_in + cell.m_solid_zu + 0.2 * cell.R_solid, 0.0)

            # 局部求解
            res_cell = solve_cell(cell)
            cell.T = old_T + 0.4 * (cell.T - old_T) # 松弛
            max_res = max(max_res, res_cell.get("residual", 0.0))
            
        if cfg.recirculation_frac > 0:
            top, bot = cells[-1], cells[0]
            bot.N_rez_d[:] = top.N_d * cfg.recirculation_frac
            bot.m_solid_zu = bot.m_solid_zu + top.m_solid * cfg.recirculation_frac
            
        dT = max(abs(cells[j].T - T_old[j]) for j in range(len(cells)))
        if g_iter >= 10 and dT < tol:
            break
            
    # 结果封装
    top = cells[-1]
    y_exit = top._mole_fractions("d")
    exit_gas = {sp: float(y_exit[j]) for j, sp in enumerate(GAS_SPECIES)}
    
    bot = cells[0]
    m_char_in = float(np.sum(bot.m_solid_zu[:, S_CHAR] + bot.m_solid_in[:, S_CHAR]))
    m_char_out = float(np.sum(top.m_solid[:, S_CHAR]))
    carbon_conv = 1.0 - (m_char_out / max(m_char_in, 1e-12))
    
    # 额外评估状态（用于测试对标）
    res_gs, rms_gs = reactor._evaluate_current_gs_state()
    
    return {
        "converged": g_iter < max_iter - 1,
        "n_iter": g_iter + 1,
        "T_profile": [c.T for c in cells],
        "exit_gas": exit_gas,
        "exit_gas_dry": wet_to_dry_mole_fractions(exit_gas),
        "carbon_conv": float(np.clip(carbon_conv, 0.0, 1.0)),
        "residual_gs": res_gs,
        "rms_scaled_gs": rms_gs,
    }
