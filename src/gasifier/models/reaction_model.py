from __future__ import annotations
import numpy as np
from src.gasifier.domain.config import SolidProperties, OperatingCondition, SimulationConfig
from src.gasifier.models.correlations.reaction_kinetics import effective_rate_multiplier
from src.core.cell_kinetics import build_reaction_sources, ReactionSourceBundle
from src.core.cell import S_CHAR, N_SOLID_COMP

class ReactionModel:
    """动力学模型：计算所有单元的化学反应速率与源项。"""
    def __init__(self, solid: SolidProperties, op: OperatingCondition, sim: SimulationConfig):
        self.solid = solid
        self.op = op
        self.sim = sim

    def compute_cell(
        self, 
        T: float, 
        P: float, 
        V_b: float, 
        V_d: float, 
        C_b: np.ndarray, 
        C_d: np.ndarray, 
        y_b: np.ndarray, 
        y_d: np.ndarray, 
        gas_src_vm: np.ndarray, 
        solid_sink_vm: np.ndarray, 
        areas: np.ndarray, 
        D_g: float, 
        char_conversion: float, 
        rho_cat: float, 
        gibbs_minor: dict | None, 
        rate_multiplier: float, 
        N_zu_d: np.ndarray, 
        N_d_in: np.ndarray, 
        N_zu_b: np.ndarray, 
        N_b_in: np.ndarray, 
        N_rez_d: np.ndarray, 
        N_rez_b: np.ndarray, 
        N_ex: np.ndarray
    ) -> ReactionSourceBundle:
        """计算单个单元的反应源项。"""
        bounded_multiplier = effective_rate_multiplier(rate_multiplier)
        return build_reaction_sources(
            T=T,
            P=P,
            fuel_type=self.op.fuel_type,
            V_b=V_b,
            V_d=V_d,
            C_b=C_b,
            C_d=C_d,
            y_b=y_b,
            y_d=y_d,
            gas_src_vm=gas_src_vm,
            solid_sink_vm=solid_sink_vm,
            areas=areas,
            solid_d_p=self.solid.d_p,
            D_g=D_g,
            char_conversion=char_conversion,
            rho_cat=rho_cat,
            enable_r12=self.sim.enable_r12,
            use_gibbs_minor=self.sim.use_gibbs_minor,
            gibbs_minor_sources=gibbs_minor,
            r4_scale=self.sim.r4_scale,
            r5_scale=self.sim.r5_scale,
            r6_scale=self.sim.r6_scale,
            r7_scale=self.sim.r7_scale,
            rate_multiplier=bounded_multiplier,
            N_zu_d=N_zu_d,
            N_d_in=N_d_in,
            N_zu_b=N_zu_b,
            N_b_in=N_b_in,
            N_rez_d=N_rez_d,
            N_rez_b=N_rez_b,
            N_ex=N_ex,
            solid_shape=(self.solid.n_size_classes, N_SOLID_COMP),
            char_index=S_CHAR,
        )
