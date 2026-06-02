from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Tuple
from src.core.species import TarFuelType

@dataclass(frozen=True)
class SimulationConfig:
    max_outer_iterations: int = 20
    max_newton_iterations: int = 50
    tol_outer: float = 1e-4
    tol_inner_rms: float = 0.01
    verbose: bool = False
    thesis_mode: bool = False
    thesis_vorab_sources_single_shot: bool = True
    # 流体力学策略
    hydro_u_d_closure: str = "wein_1992_eq312"
    hydro_bubble_diameter_model: str = "hilligardt_ode"
    hydro_psi_b_strategy: str = "wein_1992"
    hydro_lambda_strategy: str = "hamel_280"
    hydro_xi_strategy: str = "hamel_regime"
    hydro_bubble_velocity_strategy: str = "heinbockel_eq343"
    hydro_bubble_ode_strategy: str = "heinbockel_eq341"
    # 反应比例因子与开关
    r4_scale: float = 1.0
    r5_scale: float = 1.0
    r6_scale: float = 1.0
    r7_scale: float = 1.0
    enable_r12: bool = True
    use_gibbs_minor: bool = True

@dataclass(frozen=True)
class PlantData:
    n_cells: int
    H_bed: float
    D_bed: float
    H_freeboard: float = 0.0
    n_freeboard_cells: int = 0
    N_or: int = 100  # Number of orifices

@dataclass(frozen=True)
class SolidProperties:
    rho_s: float = 1400.0
    d_p: float = 0.00225
    phi_s: float = 0.86
    eps_mf: float = 0.45
    n_size_classes: int = 1
    d_p_classes: Tuple[float, ...] = (0.00225,)
    mass_fractions: Tuple[float, ...] = (1.0,)
    # Fuel analysis
    moisture_wt: float = 16.9
    ash_dry_wt: float = 11.41
    VM_daf: float = 53.42
    C_dry: float = 61.5
    H_dry: float = 4.1
    O_dry: float = 21.8
    S_dry: float = 0.0
    HHV_dry_MJ_kg: float = 22.0
    nitrogen_fraction: float = 0.0
    sulfur_fraction: float = 0.0
    sulfur_volatile_frac: float = 0.5

@dataclass(frozen=True)
class OperatingCondition:
    P: float = 2_500_000.0
    T_inlet: float = 300.0
    fuel_feed_kg_s: float = 0.938
    fuel_type: TarFuelType = "coal"
    ER: Optional[float] = None
    primary_agent: str = "air_steam"
    # Gas feeds if ER is None
    O2_feed: float = 0.0
    H2O_feed: float = 0.0
    N2_feed: float = 0.0
    # Other settings
    heat_loss_frac: float = 0.1
    recirculation_frac: float = 0.1
    u0_target: Optional[float] = None
