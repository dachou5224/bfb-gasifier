"""Phase 1 单元测试：species.py + arrhenius.py

运行方式：cd bfb-gasifier && python -m pytest tests/test_phase1.py -v
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pytest

from src.core.species import (
    cp_molar,
    enthalpy_molar,
    entropy_molar,
    formation_enthalpy_298,
    cp_char,
    cp_ash,
    cp_sand,
    GAS_SPECIES,
    MOLECULAR_WEIGHT,
    TAR_SURROGATE_HC_RATIO,
    calc_tar_surrogate_fractions,
    configure_tar_components_by_fuel,
    mean_molar_mass,
    gas_density_ideal,
    gas_viscosity_power_law,
    gas_diffusivity_correlation,
)
from src.kinetics.arrhenius import k_hobbs, k_standard, k_jensen_r7


# ===== species.py 测试 =====

class TestCpMolar:
    """验证各组分在典型温度下的 Cp 值与文献吻合。"""

    def test_cp_co2_300K(self):
        """Cp_CO2(300K) ≈ 37.1 J/(mol·K)（NIST WebBook 参考值）"""
        cp = cp_molar("CO2", 300.0)
        assert 35.0 < cp < 39.0, f"Cp_CO2(300K) = {cp:.2f}, expected ~37"

    def test_cp_h2o_1000K(self):
        """Cp_H2O(1000K) ≈ 41.2 J/(mol·K)"""
        cp = cp_molar("H2O", 1000.0)
        assert 38.0 < cp < 44.0, f"Cp_H2O(1000K) = {cp:.2f}, expected ~41"

    def test_cp_o2_298K(self):
        """Cp_O2(298K) ≈ 29.4 J/(mol·K)"""
        cp = cp_molar("O2", 298.15)
        assert 28.0 < cp < 31.0, f"Cp_O2(298K) = {cp:.2f}, expected ~29.4"

    def test_cp_h2_298K(self):
        """Cp_H2(298K) ≈ 28.8 J/(mol·K)"""
        cp = cp_molar("H2", 298.15)
        assert 27.0 < cp < 30.0, f"Cp_H2(298K) = {cp:.2f}, expected ~28.8"

    def test_cp_ch4_298K(self):
        """Cp_CH4(298K) ≈ 35.7 J/(mol·K)"""
        cp = cp_molar("CH4", 298.15)
        assert 33.0 < cp < 38.0, f"Cp_CH4(298K) = {cp:.2f}, expected ~35.7"

    def test_cp_h2s_298K(self):
        """Cp_H2S(298K) ≈ 34.2 J/(mol·K)"""
        cp = cp_molar("H2S", 298.15)
        assert 32.0 < cp < 37.0, f"Cp_H2S(298K) = {cp:.2f}, expected ~34.2"

    def test_cp_nh3_298K(self):
        """Cp_NH3(298K) ≈ 35.0 J/(mol·K)"""
        cp = cp_molar("NH3", 298.15)
        assert 33.0 < cp < 38.0, f"Cp_NH3(298K) = {cp:.2f}, expected ~35"

    def test_cp_positive(self):
        """所有组分在 300-2000K 范围内 Cp > 0"""
        for species in GAS_SPECIES:
            if species in {"TAR1", "TAR2"}:
                continue
            for T in [300.0, 500.0, 1000.0, 1500.0, 2000.0]:
                cp = cp_molar(species, T)
                assert cp > 0, f"Cp({species}, {T}K) = {cp} <= 0"


class TestEnthalpy:
    """验证焓值与已知标准生成焓吻合。"""

    def test_h_co2_298(self):
        """h_f(CO2, 298K) ≈ -393.51 kJ/mol"""
        h = formation_enthalpy_298("CO2")
        assert -400e3 < h < -387e3, f"h_f(CO2) = {h/1e3:.1f} kJ/mol"

    def test_h_h2o_298(self):
        """h_f(H2O, 298K) ≈ -241.83 kJ/mol"""
        h = formation_enthalpy_298("H2O")
        assert -248e3 < h < -236e3, f"h_f(H2O) = {h/1e3:.1f} kJ/mol"

    def test_h_co_298(self):
        """h_f(CO, 298K) ≈ -110.53 kJ/mol"""
        h = formation_enthalpy_298("CO")
        assert -117e3 < h < -104e3, f"h_f(CO) = {h/1e3:.1f} kJ/mol"

    def test_h_h2_298(self):
        """h_f(H2, 298K) ≈ 0 kJ/mol（元素基准态）"""
        h = formation_enthalpy_298("H2")
        assert -5e3 < h < 5e3, f"h_f(H2) = {h/1e3:.1f} kJ/mol"

    def test_h_o2_298(self):
        """h_f(O2, 298K) ≈ 0 kJ/mol（元素基准态）"""
        h = formation_enthalpy_298("O2")
        assert -5e3 < h < 5e3, f"h_f(O2) = {h/1e3:.1f} kJ/mol"

    def test_h_monotonic(self):
        """H(T) 在 300-1200K 间单调递增"""
        for species in ["CO2", "H2O", "CO", "N2"]:
            temps = np.linspace(300, 1200, 50)
            enthalpies = [enthalpy_molar(species, T) for T in temps]
            for i in range(1, len(enthalpies)):
                assert enthalpies[i] > enthalpies[i-1], (
                    f"H({species}) not monotonic at T={temps[i]:.0f}K"
                )

    def test_h_ch4_298(self):
        """h_f(CH4, 298K) ≈ -74.87 kJ/mol"""
        h = formation_enthalpy_298("CH4")
        assert -82e3 < h < -68e3, f"h_f(CH4) = {h/1e3:.1f} kJ/mol"


class TestTarComponentPlaceholder:
    """TAR1/TAR2 在未注册时应抛出 NotImplementedError。"""

    def test_tar1_cp_raises(self):
        with pytest.raises(NotImplementedError):
            cp_molar("TAR1", 1000.0)

    def test_tar2_enthalpy_raises(self):
        with pytest.raises(NotImplementedError):
            enthalpy_molar("TAR2", 1000.0)


class TestTarSurrogateModel:
    """Tar 两组分代理模型 H/C 一致性验证。"""

    def test_tar_coal_hc_match(self):
        frac = calc_tar_surrogate_fractions("coal")
        r_mix = sum(frac[s] * TAR_SURROGATE_HC_RATIO[s] for s in frac)
        assert abs(r_mix - 0.856) < 1e-10
        assert abs(sum(frac.values()) - 1.0) < 1e-12

    def test_tar_biomass_hc_match(self):
        frac = calc_tar_surrogate_fractions("biomass")
        r_mix = sum(frac[s] * TAR_SURROGATE_HC_RATIO[s] for s in frac)
        assert abs(r_mix - 1.364) < 1e-10
        assert abs(sum(frac.values()) - 1.0) < 1e-12

    def test_tar_custom_hc_out_of_range_raises(self):
        with pytest.raises(ValueError):
            calc_tar_surrogate_fractions("coal", target_hc_ratio=1.2)

    def test_tar_mapping_changes_with_feedstock(self):
        coal_map = configure_tar_components_by_fuel("coal")
        biomass_map = configure_tar_components_by_fuel("biomass")
        assert coal_map == {"TAR1": "C6H6", "TAR2": "C10H8"}
        assert biomass_map == {"TAR1": "C10H8", "TAR2": "C16H34"}

    def test_nh3_always_in_species_list(self):
        assert "NH3" in GAS_SPECIES


class TestSolidCp:
    """固体组分比热验证。"""

    def test_char_cp_range(self):
        """炭 Cp 在 800-2000 J/(kg·K) 范围内（T=1000K）"""
        cp = cp_char(1000.0)
        assert 800 < cp < 2000, f"Cp_char(1000K) = {cp}"

    def test_ash_cp_range(self):
        """灰 Cp 在 700-1300 J/(kg·K) 范围内（T=1000K）"""
        cp = cp_ash(1000.0)
        assert 700 < cp < 1300, f"Cp_ash(1000K) = {cp}"

    def test_sand_cp_range(self):
        """SiO2 Cp 在 800-1400 J/(kg·K) 范围内（T=1000K）"""
        cp = cp_sand(1000.0)
        assert 800 < cp < 1400, f"Cp_sand(1000K) = {cp}"


class TestMolecularWeight:
    """分子量基础检查。"""

    def test_all_species_have_mw(self):
        for s in GAS_SPECIES:
            assert s in MOLECULAR_WEIGHT, f"{s} missing MW"

    def test_mw_values(self):
        assert abs(MOLECULAR_WEIGHT["CO2"] - 44.01) < 0.1
        assert abs(MOLECULAR_WEIGHT["H2O"] - 18.015) < 0.1
        assert abs(MOLECULAR_WEIGHT["N2"] - 28.013) < 0.1


class TestGasMixtureProperties:
    """混合气物性关联测试。"""

    def test_mean_molar_mass_air_like(self):
        y = {"O2": 0.21, "N2": 0.79}
        m = mean_molar_mass(y)
        assert 28.7 < m < 29.1

    def test_gas_density_ideal_n2(self):
        y = {"N2": 1.0}
        rho = gas_density_ideal(P=101325.0, T=300.0, mole_fractions=y)
        assert 1.10 < rho < 1.20

    def test_gas_viscosity_power_law_monotonic(self):
        mu1 = gas_viscosity_power_law(T=500.0, mu_g_20=1.8e-5, n=0.7)
        mu2 = gas_viscosity_power_law(T=1000.0, mu_g_20=1.8e-5, n=0.7)
        assert mu2 > mu1

    def test_gas_diffusivity_pressure_inverse(self):
        d1 = gas_diffusivity_correlation(T_m=1500.0, P=101300.0)
        d2 = gas_diffusivity_correlation(T_m=1500.0, P=202600.0)
        assert abs(d2 / d1 - 0.5) < 1e-12


# ===== arrhenius.py 测试 =====

class TestArrhenius:
    """验证三种 Arrhenius 工厂函数。"""

    def test_k_hobbs_magnitude(self):
        """k_hobbs(k0=1e6, E=1e5, T=1000) 应为正有限值"""
        k = k_hobbs(1e6, 1e5, 1000.0)
        assert np.isfinite(k), "k_hobbs returned non-finite"
        assert k > 0, "k_hobbs returned non-positive"

    def test_k_hobbs_includes_T_factor(self):
        """k_hobbs 含 T 前因子，在相同 exp 项下，T=2000 的 k 约为 T=1000 的 2 倍"""
        k1 = k_hobbs(1.0, 0.0, 1000.0)  # E=0 时 k = k0 * T
        k2 = k_hobbs(1.0, 0.0, 2000.0)
        assert abs(k2 / k1 - 2.0) < 0.01

    def test_k_standard_limit(self):
        """k_standard 在 T->inf（E/RT->0）时趋近 k0"""
        k = k_standard(42.0, 1e5, 1e10)  # 极高温
        assert abs(k - 42.0) < 0.01, f"k_standard at T=1e10 = {k}, expected ~42"

    def test_k_standard_low_T_no_overflow(self):
        """k_standard 在极低温时不溢出（np.clip 生效）"""
        k = k_standard(1e10, 1e6, 1.0)  # T=1K, E/RT 极大
        assert np.isfinite(k), "k_standard overflow at T=1K"
        assert k >= 0

    def test_k_jensen_r7_no_overflow(self):
        """k_jensen_r7 在极低温时不溢出"""
        k = k_jensen_r7(1e10, 1e6, 1.0)
        assert np.isfinite(k), "k_jensen_r7 overflow at T=1K"

    def test_k_jensen_r7_positive(self):
        """k_jensen_r7 在正常范围内返回正值"""
        k = k_jensen_r7(1e8, 20000.0, 1200.0)
        assert k > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
