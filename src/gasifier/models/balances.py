from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt
from src.core.cell_balances import (
    calc_gas_balance_residual,
    calc_solid_balance_residual,
    calc_energy_balance_residual,
    assemble_cell_residual_vector,
)
from src.gasifier.domain.state import GlobalState
from src.gasifier.domain.config import SolidProperties

class BalanceBuilder:
    """平衡构建器：根据各物理模型的计算结果组装全局残差向量。"""
    def __init__(self, solid_props: SolidProperties):
        self.solid_props = solid_props

    def compute_cell_residual(
        self,
        cell_idx: int,
        state: GlobalState,
        hydro_bundle, # HydrodynamicsBundle
        pyro_bundle,  # PyrolysisSourceBundle
        reaction_bundle, # ReactionSourceBundle
        boundary_conditions: dict,
        heat_loss_frac: float
    ) -> npt.NDArray[np.float64]:
        """计算单个单元的残差向量。"""
        
        # 1. 气相平衡
        res_gas = calc_gas_balance_residual(
            N_zu_d=boundary_conditions["N_zu_d"],
            N_rez_d=boundary_conditions["N_rez_d"],
            N_d_in=boundary_conditions["N_d_in"],
            R_gas_d=reaction_bundle.R_gas_d,
            N_d=state.N_d[cell_idx],
            N_ex=boundary_conditions["N_ex"],
            N_zu_b=boundary_conditions["N_zu_b"],
            N_rez_b=boundary_conditions["N_rez_b"],
            N_b_in=boundary_conditions["N_b_in"],
            R_gas_b=reaction_bundle.R_gas_b,
            N_b=state.N_b[cell_idx],
        )

        # 2. 固相平衡
        res_solid = calc_solid_balance_residual(
            m_solid_zu=boundary_conditions["m_solid_zu"],
            m_solid_rez=boundary_conditions["m_solid_rez"],
            m_solid_in=boundary_conditions["m_solid_in"],
            R_solid=reaction_bundle.R_solid,
            size_migration=boundary_conditions["size_migration"],
            m_solid=state.m_solid[cell_idx],
            m_solid_auf_in=boundary_conditions.get("m_solid_auf_in"),
            m_solid_ab_in=boundary_conditions.get("m_solid_ab_in"),
            K_solid_auf=boundary_conditions.get("K_solid_auf"),
            K_solid_ab=boundary_conditions.get("K_solid_ab"),
            solid_state_model=boundary_conditions.get("solid_state_model", "legacy_stream"),
        ).flatten()

        # 3. 能量平衡
        res_energy = calc_energy_balance_residual(
            N_b_in=boundary_conditions["N_b_in"],
            N_d_in=boundary_conditions["N_d_in"],
            T_in_gas=boundary_conditions["T_in_gas"],
            N_zu_b=boundary_conditions["N_zu_b"],
            N_zu_d=boundary_conditions["N_zu_d"],
            T_zu_gas=boundary_conditions["T_zu_gas"],
            N_rez_b=boundary_conditions["N_rez_b"],
            N_rez_d=boundary_conditions["N_rez_d"],
            T_rez_gas=boundary_conditions["T_rez_gas"],
            m_solid_rez=boundary_conditions["m_solid_rez"],
            T_rez_solid=boundary_conditions["T_rez_solid"],
            m_solid_in=boundary_conditions["m_solid_in"],
            T_in_solid=boundary_conditions["T_in_solid"],
            m_solid_zu=boundary_conditions["m_solid_zu"],
            T_zu_solid=boundary_conditions["T_zu_solid"],
            N_b=state.N_b[cell_idx],
            N_d=state.N_d[cell_idx],
            m_solid=boundary_conditions["m_solid_outflow"], # 此处需由 Boundary 控制
            T=state.T[cell_idx],
            heat_loss_frac=heat_loss_frac,
            ash_dry_wt=self.solid_props.ash_dry_wt,
            VM_daf=self.solid_props.VM_daf,
            h_f_dry=self.solid_props.h_f_dry,
        )

        # 组装
        n_vars = res_gas.size + res_solid.size + 1
        out = np.zeros(n_vars, dtype=np.float64)
        return assemble_cell_residual_vector(
            res_gas=res_gas,
            res_solid=res_solid,
            res_energy=res_energy,
            out=out
        )
