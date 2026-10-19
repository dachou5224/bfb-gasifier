"""读取 data/properties/。系数不写在代码里。"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2] / "data" / "properties"
# 第 1 关气相。不含 SO2（ADR 0021 第 6 项）、H2O(L)、C(gr)。
GAS_SPECIES = (
    "CO", "CO2", "H2", "H2O", "CH4", "C2H4", "N2", "O2", "NH3", "H2S", "Ar",
    "C6H6", "C10H8",
)
class PropertyError(ValueError):
    """物性或平衡输入不能计算。"""
class PropertyRangeError(PropertyError):
    """温度或对比压力超出规定区间。"""
class EquilibriumError(PropertyError):
    """自由焓最小化没有得到收敛解。"""
@lru_cache(maxsize=1)
def _bundle() -> dict:
    return {
        "constants": json.loads((_ROOT / "constants.json").read_text(encoding="utf-8")),
        "nasa": json.loads((_ROOT / "nasa_thermo.json").read_text(encoding="utf-8")),
        "transport": json.loads((_ROOT / "transport.json").read_text(encoding="utf-8")),
        "ash": json.loads((_ROOT / "ash.json").read_text(encoding="utf-8")),
    }
def gas_constant_j_mol_k() -> float:
    """CODATA 2018 的 R。ADR 0018 第 3 项。"""
    return float(_bundle()["constants"]["gas_constant_J_mol_K"]["value"])
def stefan_boltzmann_w_m2_k4() -> float:
    """CODATA 2018 的 σ。ADR 0018 第 3 项。"""
    return float(_bundle()["constants"]["stefan_boltzmann_W_m2_K4"]["value"])
def standard_pressure_pa() -> float:
    """Hamel 的 P_0。"""
    return float(_bundle()["constants"]["standard_pressure_Pa"]["value"])
def reference_pressure_pa() -> float:
    """NASA 气体标准压力，1 bar。"""
    return float(_bundle()["constants"]["reference_pressure_bar_Pa"]["value"])
def atomic_weight(element: str) -> float:
    weights = _bundle()["constants"]["atomic_weights"]
    if element not in weights:
        raise PropertyError(f"没有 {element} 的原子量")
    return float(weights[element])
def nasa_species(name: str) -> dict:
    species = _bundle()["nasa"]["species"]
    if name not in species:
        raise PropertyError(f"NASA 数据里没有 {name}")
    return species[name]
def atoms(name: str) -> dict[str, int]:
    raw = nasa_species(name)["atoms"]
    return {element: int(round(count)) for element, count in raw.items()}
def molar_mass_g_mol(name: str) -> float:
    """由原子量求和。不用 NASA 记录里的分子量。"""
    return sum(count * atomic_weight(element) for element, count in atoms(name).items())
def molar_mass_kg_mol(name: str) -> float:
    return molar_mass_g_mol(name) * 1e-3
def transport_species(name: str) -> dict:
    species = _bundle()["transport"]["species"]
    if name not in species:
        raise PropertyError(f"输运数据里没有 {name}")
    return species[name]
def ash_record() -> dict:
    return _bundle()["ash"]
def is_gas(name: str) -> bool:
    return int(nasa_species(name)["phase_code"]) == 0
