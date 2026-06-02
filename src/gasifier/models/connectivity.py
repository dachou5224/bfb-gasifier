from __future__ import annotations
import numpy as np
from src.gasifier.domain.state import GlobalState
from src.gasifier.domain.config import PlantData, OperatingCondition
from src.core.species import GAS_SPECIES_INDEX, N_GAS
from src.core.cell import S_CHAR, S_ASH, S_MOISTURE, S_VM, N_SOLID_COMP

class ConnectivityModel:
    """连接模型：负责在 GlobalState 中各单元之间路由流束 (Stream Routing)。"""
    
    def __init__(self, plant: PlantData, op: OperatingCondition):
        self.plant = plant
        self.op = op

    def compute_boundary_conditions(self, state: GlobalState, hydro_bundles) -> list[dict]:
        """计算每个单元的边界入流条件。"""
        n_cells = state.n_cells
        idx = GAS_SPECIES_INDEX
        bc_list = []
        
        for i in range(n_cells):
            bc = {
                "N_zu_d": np.zeros(N_GAS),
                "N_zu_b": np.zeros(N_GAS),
                "N_rez_d": np.zeros(N_GAS),
                "N_rez_b": np.zeros(N_GAS),
                "N_d_in": np.zeros(N_GAS),
                "N_b_in": np.zeros(N_GAS),
                "N_ex": np.zeros(N_GAS), # 需由 hydro 决定
                "m_solid_zu": np.zeros((state.n_size_classes, N_SOLID_COMP)),
                "m_solid_rez": np.zeros((state.n_size_classes, N_SOLID_COMP)),
                "m_solid_in": np.zeros((state.n_size_classes, N_SOLID_COMP)),
                "m_solid_outflow": np.zeros((state.n_size_classes, N_SOLID_COMP)),
                "T_in_gas": 293.15,
                "T_zu_gas": self.op.T_inlet,
                "T_rez_gas": 293.15,
                "T_in_solid": 293.15,
                "T_zu_solid": 293.15,
                "T_rez_solid": 293.15,
                "size_migration": np.zeros((state.n_size_classes, N_SOLID_COMP)), # 简化
            }
            
            # 底部单元
            if i == 0:
                bc["N_zu_d"][idx["O2"]] = self.op.O2_feed
                bc["N_zu_d"][idx["H2O"]] = self.op.H2O_feed
                bc["N_zu_d"][idx["N2"]] = self.op.N2_feed
                bc["m_solid_zu"][:, S_CHAR] = self.op.fuel_feed_kg_s / state.n_cells # 极简假设
            else:
                # 气相从下方单元流入 (i-1 -> i)
                bc["N_d_in"] = state.N_d[i-1].copy()
                bc["N_b_in"] = state.N_b[i-1].copy()
                bc["T_in_gas"] = float(state.T[i-1])
                
            # 固相传输 (i+1 -> i)
            if i < n_cells - 1:
                 bc["m_solid_in"] = state.m_solid[i+1].copy()
                 bc["T_in_solid"] = float(state.T[i+1])
                 
            # 传质由当前单元的 hydro 决定
            h = hydro_bundles[i]
            # 此处应调用 calc_phase_exchange
            # 为演示先占位
            bc["N_ex"] = np.zeros(N_GAS)
            
            # 最终出流
            bc["m_solid_outflow"] = state.m_solid[i].copy()
            
            bc_list.append(bc)
            
        return bc_list
