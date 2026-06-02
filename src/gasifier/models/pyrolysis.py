from __future__ import annotations
from src.gasifier.domain.state import GlobalState
from src.gasifier.domain.config import PlantData, SolidProperties, OperatingCondition
from src.core.cell_pyrolysis import calc_drying_pyrolysis_sources, PyrolysisSourceBundle
from src.core.cell import S_VM, S_MOISTURE, S_CHAR, N_SOLID_COMP

class PyrolysisModel:
    """热解与干燥模型：计算水分蒸发和挥发分析出源项。"""
    def __init__(self, plant: PlantData, solid: SolidProperties, op: OperatingCondition):
        self.plant = plant
        self.solid = solid
        self.op = op

    def compute_cell(
        self, 
        T: float, 
        P: float, 
        tau: float, 
        T_init: float, 
        m_vm_in: float, 
        m_moist_in: float
    ) -> PyrolysisSourceBundle:
        """计算单个单元的热解源项。"""
        return calc_drying_pyrolysis_sources(
            tau=tau,
            T=T,
            P=P,
            T_init=T_init,
            d_p=self.solid.d_p,
            moisture_wt=self.solid.moisture_wt,
            ash_dry_wt=self.solid.ash_dry_wt,
            C_dry=self.solid.C_dry,
            H_dry=self.solid.H_dry,
            O_dry=self.solid.O_dry,
            nitrogen_fraction=self.solid.nitrogen_fraction,
            sulfur_fraction=self.solid.sulfur_fraction,
            sulfur_volatile_frac=self.solid.sulfur_volatile_frac,
            pyrolysis_tar_carbon_frac=0.2,  # 可选：从 config 传入
            fuel_type=self.op.fuel_type,
            m_vm_in=m_vm_in,
            m_moist_in=m_moist_in,
            solid_shape=(self.solid.n_size_classes, N_SOLID_COMP),
            char_index=S_CHAR,
            vm_index=S_VM,
            moisture_index=S_MOISTURE,
        )
