from __future__ import annotations
from src.gasifier.domain.config import PlantData, OperatingCondition
from src.gasifier.models.correlations.heat_transfer_coeff import compound_heat_loss_fraction

class HeatTransferModel:
    """热传导与热损失模型：计算反应器的轴向热损失分布。"""
    def __init__(self, plant: PlantData, op: OperatingCondition):
        self.plant = plant
        self.op = op

    def compute_distribution(self) -> list[float]:
        """计算每个单元的热损失比例。"""
        bed_height = float(max(self.plant.H_bed, 0.0))
        freeboard_height = float(max(self.plant.H_freeboard, 0.0))
        total_height = bed_height + freeboard_height

        if total_height <= 0.0:
            return [0.0] * self.plant.n_cells

        n_bed = self.plant.n_cells
        bed_weight_each = (bed_height / total_height) / max(n_bed, 1)
        bed_losses = [compound_heat_loss_fraction(self.op.heat_loss_frac, bed_weight_each) for _ in range(n_bed)]
        
        # 备注：暂不支持自由板热损失
        return bed_losses
