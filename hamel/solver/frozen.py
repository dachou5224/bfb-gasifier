"""内层冻结量。T25b 填写这一结构；T25a 只读取。测试可以直接构造。"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hamel.properties.data import GAS_SPECIES

class AssemblyError(ValueError):
    """装配输入或冻结量不成立。不改残差。"""

class TrialNotEvaluable(AssemblyError):
    """试探点上浓度、颗粒温度或灰分配没有定义。
    求解手段：ADR 0022 C4。T25c 与速率式越界（KineticsInputError、PropertyError）同样判为不可行。
    """

class FreeboardDomainError(AssemblyError):
    """ADR 0026 第 5 项。报告并停止，不截断，也不靠缩短步长把这一项变成可行。"""

def _array(name: str, value: object, shape: tuple[int, ...]) -> np.ndarray:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise AssemblyError(f"{name} 不能转换成浮点数组") from exc
    if array.shape != shape:
        raise AssemblyError(f"{name} 的形状是 {array.shape}，应为 {shape}")
    if not np.isfinite(array).all():
        raise AssemblyError(f"{name} 含 NaN 或无穷大")
    return array

def _mask(name: str, value: object, shape: tuple[int, ...]) -> np.ndarray:
    array = np.asarray(value, dtype=bool)
    if array.shape != shape:
        raise AssemblyError(f"{name} 的形状是 {array.shape}，应为 {shape}")
    return array

@dataclass(frozen=True)
class FrozenQuantities:
    """一轮外层给出的系数和边界量。数组第 0 轴是衡算编号，0 为该区最上层。
    气体出流、交换流、反应项、焓流和 K·m 不在这里。它们在试探点上计算。
    解释性选择：ADR 0017 第 12 项；T18 答卷第 4 节。
    """
    diameter_m: np.ndarray
    bed_pressure_pa: np.ndarray
    freeboard_pressure_pa: np.ndarray
    k_auf_per_s: np.ndarray
    m_dot_ab_kg_s: np.ndarray
    geometric_mass_kg: np.ndarray
    k_bd_per_s: np.ndarray
    v_bubble_exchange_m3: np.ndarray
    v_suspension_m3: np.ndarray
    v_bubble_m3: np.ndarray
    v_freeboard_m3: np.ndarray
    primary_feed_mol_s: np.ndarray
    feed_to_suspension: float
    primary_temperature_k: float
    secondary_feed_mol_s: np.ndarray
    secondary_temperature_k: float
    secondary_balance_index: int
    pyro_bed_mol_s: np.ndarray
    pyro_freeboard_mol_s: np.ndarray
    char_gen_bed_kg_s: np.ndarray
    char_gen_free_kg_s: np.ndarray
    char_gen_hover_kg_s: np.ndarray
    ash_gen_bed_kg_s: np.ndarray
    ash_gen_free_kg_s: np.ndarray
    ash_gen_hover_kg_s: np.ndarray
    fuel_bed_kg: np.ndarray
    fuel_up_bed_kg_s: np.ndarray
    fuel_down_bed_kg_s: np.ndarray
    fuel_exit_bed_kg_s: np.ndarray
    fuel_up_free_kg_s: np.ndarray
    fuel_down_free_kg_s: np.ndarray
    fuel_feed_kg_s: np.ndarray
    fuel_formation_j_per_kg: float
    fuel_heat_capacity_j_per_kg_k: float
    fuel_feed_temperature_k: float
    hover: np.ndarray
    free_absent: np.ndarray
    a_free_per_s: np.ndarray
    a_hover_per_s: np.ndarray
    b_free_per_s: np.ndarray
    b_hover_per_s: np.ndarray
    theta: np.ndarray
    ejection_to_free: np.ndarray
    ejection_to_hover: np.ndarray
    u_r_bed_m_s: np.ndarray
    u_r_freeboard_m_s: np.ndarray
    u_r_hover_m_s: np.ndarray
    rho_p_kg_m3: float
    drain_char_kg_s: np.ndarray
    drain_ash_kg_s: float
    k_w_w_m2_k: float
    bed_height_m: float
    freeboard_height_m: float
    bed_diameter_m: float
    freeboard_diameter_m: float
    epsilon_mf: float
    particle_tolerance_w: float
    particle_max_iterations: int
    coupling_tolerance_k: float
    coupling_max_iterations: int
    def __post_init__(self) -> None:
        n_bed = len(np.asarray(self.bed_pressure_pa))
        n_free = len(np.asarray(self.freeboard_pressure_pa))
        n_class = len(np.asarray(self.diameter_m))
        n_species = len(GAS_SPECIES)
        if n_bed < 1 or n_free < 1 or n_class < 1:
            raise AssemblyError("床区、自由空间和粒级至少各有一个")
        bed, free, klass = (n_bed,), (n_free,), (n_class,)
        groups = (
            (klass, "diameter_m", "fuel_feed_kg_s", "ejection_to_free", "ejection_to_hover", "drain_char_kg_s"),
            (bed, "bed_pressure_pa", "k_auf_per_s", "m_dot_ab_kg_s", "geometric_mass_kg", "k_bd_per_s", "v_bubble_exchange_m3", "v_suspension_m3", "v_bubble_m3", "u_r_bed_m_s"),
            (free, "freeboard_pressure_pa", "v_freeboard_m3"),
            ((n_species,), "primary_feed_mol_s", "secondary_feed_mol_s"),
            ((n_bed, n_species), "pyro_bed_mol_s"),
            ((n_free, n_species), "pyro_freeboard_mol_s"),
            ((n_bed, n_class), "char_gen_bed_kg_s", "ash_gen_bed_kg_s", "fuel_bed_kg", "fuel_up_bed_kg_s", "fuel_down_bed_kg_s", "fuel_exit_bed_kg_s"),
            ((n_free, n_class), "char_gen_free_kg_s", "char_gen_hover_kg_s", "ash_gen_free_kg_s", "ash_gen_hover_kg_s", "fuel_up_free_kg_s", "fuel_down_free_kg_s", "a_free_per_s", "a_hover_per_s", "b_free_per_s", "b_hover_per_s", "theta", "u_r_freeboard_m_s", "u_r_hover_m_s"),
        )
        for shape, *names in groups:
            for name in names:
                object.__setattr__(self, name, _array(name, getattr(self, name), shape))
        object.__setattr__(self, "hover", _mask("hover", self.hover, (n_free, n_class)))
        object.__setattr__(self, "free_absent", _mask("free_absent", self.free_absent, (n_free, n_class)))
        self._check_scalars(n_free)
        self._check_states()
    def _check_scalars(self, n_free: int) -> None:
        share = float(self.feed_to_suspension)
        if not np.isfinite(share) or not 0.0 <= share <= 1.0:
            raise AssemblyError("进料进入悬浮相的份额必须在 [0, 1]")
        index = self.secondary_balance_index
        if isinstance(index, bool) or not isinstance(index, (int, np.integer)):
            raise AssemblyError("二次风所在格必须是整数衡算下标")
        if not 0 <= int(index) < n_free:
            raise AssemblyError("二次风所在格超出自由空间")
        object.__setattr__(self, "secondary_balance_index", int(index))
        positive = {
            "primary_temperature_k": self.primary_temperature_k,
            "secondary_temperature_k": self.secondary_temperature_k,
            "fuel_feed_temperature_k": self.fuel_feed_temperature_k,
            "rho_p_kg_m3": self.rho_p_kg_m3,
            "bed_height_m": self.bed_height_m,
            "freeboard_height_m": self.freeboard_height_m,
            "bed_diameter_m": self.bed_diameter_m,
            "freeboard_diameter_m": self.freeboard_diameter_m,
            "particle_tolerance_w": self.particle_tolerance_w,
            "coupling_tolerance_k": self.coupling_tolerance_k,
        }
        for name, value in positive.items():
            number = float(value)
            if not np.isfinite(number) or number <= 0.0:
                raise AssemblyError(f"{name} 必须是有限正数")
            object.__setattr__(self, name, number)
        for name in ("fuel_formation_j_per_kg", "fuel_heat_capacity_j_per_kg_k", "drain_ash_kg_s", "k_w_w_m2_k"):
            number = float(getattr(self, name))
            if not np.isfinite(number):
                raise AssemblyError(f"{name} 必须是有限数")
            object.__setattr__(self, name, number)
        epsilon = float(self.epsilon_mf)
        if not np.isfinite(epsilon) or not 0.0 < epsilon < 1.0:
            raise AssemblyError("ε_mf 必须在 (0, 1) 内")
        object.__setattr__(self, "epsilon_mf", epsilon)
        for name in ("particle_max_iterations", "coupling_max_iterations"):
            number = getattr(self, name)
            if isinstance(number, bool) or not isinstance(number, (int, np.integer)) or int(number) < 1:
                raise AssemblyError(f"{name} 必须是不小于 1 的整数")
            object.__setattr__(self, name, int(number))
        if np.any(self.geometric_mass_kg <= 0.0) or np.any(self.diameter_m <= 0.0):
            raise AssemblyError("几何存量和粒径必须为正")
        if np.any(np.diff(self.diameter_m) <= 0.0):
            raise AssemblyError("粒径必须严格从小到大")
        if self.m_dot_ab_kg_s[-1] != 0.0:
            raise AssemblyError("床区最下层没有更下方的格子，向下总流必须为 0")
    def _check_states(self) -> None:
        if np.any(self.free_absent & ~self.hover):
            raise AssemblyError("自由状态不存在的格子必须设有悬停状态")
        if np.any(self.char_gen_free_kg_s[self.free_absent] != 0.0) or np.any(self.ash_gen_free_kg_s[self.free_absent] != 0.0):
            raise AssemblyError("自由状态不存在时不能向该状态给生成量")
        if np.any(self.char_gen_hover_kg_s[~self.hover] != 0.0) or np.any(self.ash_gen_hover_kg_s[~self.hover] != 0.0):
            raise AssemblyError("没有悬停状态时不能向悬停状态给生成量")
        if np.any((self.theta != 0.0) & (self.theta != 1.0)):
            raise AssemblyError("θ 只能是 0 或 1")
        if np.any(self.free_absent & (self.theta != 1.0)):
            raise AssemblyError("自由状态不存在时 θ 必须为 1")
        if np.any(~self.hover & (self.theta != 0.0)):
            raise AssemblyError("没有悬停状态时 θ 必须为 0")
        if np.any(np.abs(self.ejection_to_free + self.ejection_to_hover - 1.0) > 1.0e-12):
            raise AssemblyError("床面抛出进入自由与悬停的份额之和必须为 1")
        bottom = self.free_absent[-1]
        if np.any(bottom & (self.ejection_to_free != 0.0)):
            raise AssemblyError("自由空间最下一格没有自由状态时，抛出不能进入自由状态")
    @property
    def n_bed(self) -> int:
        """床区格子数。Source: Hamel Eq.2.4, S. 16"""
        return int(self.bed_pressure_pa.shape[0])
    @property
    def n_freeboard(self) -> int:
        """自由空间格子数。Source: Hamel Eq.2.4, S. 16"""
        return int(self.freeboard_pressure_pa.shape[0])
    @property
    def n_class(self) -> int:
        """粒级数。k = 1 是最小粒级。Source: Hamel Eq.2.6, S. 17"""
        return int(self.diameter_m.shape[0])
