from __future__ import annotations
from src.gasifier.domain.state import GlobalState
from src.gasifier.domain.config import PlantData, SolidProperties, OperatingCondition, SimulationConfig
from src.core.cell_hydrodynamics import calc_cell_hydrodynamics, HydrodynamicsBundle

class HydrodynamicsModel:
    """流体力学模型：计算所有单元的流动参数。"""
    def __init__(self, plant: PlantData, solid: SolidProperties, op: OperatingCondition, sim: SimulationConfig):
        self.plant = plant
        self.solid = solid
        self.op = op
        self.sim = sim
        self.dh = plant.H_bed / plant.n_cells

    def compute(self, state: GlobalState) -> list[HydrodynamicsBundle]:
        bundles = []
        for i in range(state.n_cells):
            h_center = (i + 0.5) * self.dh
            bundle = calc_cell_hydrodynamics(
                T=state.T[i],
                P=state.P[i],
                N_b=state.N_b[i],
                N_d=state.N_d[i],
                D_bed=self.plant.D_bed,
                dh=self.dh,
                h_center=h_center,
                N_or=self.plant.N_or,
                rho_s=self.solid.rho_s,
                d_p=self.solid.d_p,
                eps_mf=self.solid.eps_mf,
                phi_s=self.solid.phi_s,
                u0_target=self.op.u0_target,
                u_d_closure=self.sim.hydro_u_d_closure,
                bubble_diameter_model=self.sim.hydro_bubble_diameter_model,
                psi_b_strategy=self.sim.hydro_psi_b_strategy,
                lambda_strategy=self.sim.hydro_lambda_strategy,
                xi_strategy=self.sim.hydro_xi_strategy,
                bubble_velocity_strategy=self.sim.hydro_bubble_velocity_strategy,
                bubble_ode_strategy=self.sim.hydro_bubble_ode_strategy,
                cell_type="bed",
            )
            bundles.append(bundle)
        return bundles
