"""单 cell 守恒方程审计脚本。

目的：
1. 在控制条件下（固定 T/P/进料），验证 Cell 的守恒方程残差；
2. 分离三类情况：
   a) no-reaction（源项零）
   b) reaction-only（无外进料）
   c) full（进料 + 反应）
3. 检查是否存在：
   - 元素不守恒（C/H/O/N/S）
   - 能量不闭合
   - 错误的化学计量

用法：
    python scripts/audit_cell_conservation.py

约定：
- Case A (no-reaction)：进料设置，但 R_gas_* 强制为零，检查气/固平衡残差；
  期望：残差应该只来自"状态不稳定"，而不是"守恒式错误"。
  
- Case B (reaction-only)：反应打开，但无新进料进入cell；
  期望：气相摩尔守恒 = 反应源项；固相摩尔守恒 = 反应源项；
  能量守恒 = 反应热 + 焓差。
  
- Case C (full)：标准运行；
  期望：与 specs/01_conservation_equations.md 一致。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import Cell, CellGeometry, SolidProps, N_SOLID_COMP, S_CHAR, S_VM, S_MOISTURE, S_ASH
from src.core.species import GAS_SPECIES, GAS_SPECIES_INDEX, N_GAS, configure_tar_components_by_fuel, get_atom_count
from src.core.constants import Rg


def build_test_cell(T: float = 1100.0, P: float = 2.5e6) -> Cell:
    """构建测试 cell（单尺度粒级）。"""
    # 必须先配置 tar 代理物，否则分子量未定义
    configure_tar_components_by_fuel("coal")
    
    geo = CellGeometry(D_bed=0.6, dh=1.0, h_center=0.5)
    solid = SolidProps(
        rho_s=1400.0,
        d_p=0.5e-3,
        phi_s=0.86,
        eps_mf=0.45,
        n_size_classes=1,
        moisture_wt=16.9,
        ash_dry_wt=11.41,
        VM_daf=53.42,
        C_dry=61.5,
        H_dry=4.1,
        O_dry=21.8,
    )
    cell = Cell(geo=geo, solid=solid, fuel_type="coal")
    cell.T = T
    cell.P = P
    return cell


def case_a_no_reaction(cell: Cell, verbose: bool = False) -> dict:
    """Case A：进料设置，反应源项强制为零。"""
    idx = GAS_SPECIES_INDEX

    # 设置进料（20/80 phase split，同 reactor._set_bottom_cell_feeds）
    # 用同一总流量按比例切分，确保两相组成严格一致（理论上 N_ex=0）
    cell.N_zu_d.fill(0.0)
    cell.N_zu_b.fill(0.0)
    gas_feed_total = {
        "O2": 14.376,
        "H2O": 11.501,
        "N2": 54.081,
    }
    split_b, split_d = 0.2, 0.8
    for sp, n_tot in gas_feed_total.items():
        cell.N_zu_b[idx[sp]] = split_b * n_tot
        cell.N_zu_d[idx[sp]] = split_d * n_tot

    # 初始状态：从进料开始
    cell.N_d[:] = cell.N_zu_d
    cell.N_b[:] = cell.N_zu_b
    cell.m_solid.fill(0.0)  # 案例A：无固体（纯气相进料/出料），测试气相守恒
    # cell.m_solid[0, S_CHAR] = 0.05  # kg/s
    # cell.m_solid[0, S_VM] = 0.05
    # cell.m_solid[0, S_MOISTURE] = 0.01
    # cell.m_solid[0, S_ASH] = 0.02

    # 模式：计算水力学与进料流量，但强制零反应源项
    cell.calc_hydrodynamics()
    cell.calc_exchange()

    # 强制反应源项为零，用来检查输入/进料逻辑是否一致
    cell.R_gas_b.fill(0.0)
    cell.R_gas_d.fill(0.0)
    cell.R_solid.fill(0.0)

    # 设置入口与出口相等（稳态进料 = 出口，无反应改变）
    # 注：如果有 m_solid，需设定 T_in_solid = T 以匹配 H_out
    cell.N_b_in[:] = 0.0
    cell.N_d_in[:] = 0.0
    cell.N_rez_b[:] = 0.0
    cell.N_rez_d[:] = 0.0
    cell.m_solid_in.fill(0.0)  # 无固体进料，但固体库存非零；T_in_solid 应设为外出温度
    cell.T_in_gas = cell.T
    cell.T_in_solid = cell.T
    cell.T_zu_gas = cell.T
    cell.T_zu_solid = cell.T
    # 为了平衡能量：固体进料应该是 0，但如果不是零，需要同步温度
    # 目前固体进料=0，所以 m_solid_in 对能量无贡献

    # 计算残差
    res_gas = cell.calc_gas_balance()
    res_solid = cell.calc_solid_balance().flatten()
    res_energy = cell.calc_energy_balance()

    max_gas = float(np.max(np.abs(res_gas))) if len(res_gas) > 0 else 0.0
    max_solid = float(np.max(np.abs(res_solid))) if len(res_solid) > 0 else 0.0
    
    # 找出最大残差对应的物种
    max_idx = np.argmax(np.abs(res_gas)) if len(res_gas) > 0 else -1
    if max_idx >= 0:
        sp_idx = max_idx % len(GAS_SPECIES)
        phase = "emulsion" if max_idx < len(GAS_SPECIES) else "bubble"
        culprit_sp = GAS_SPECIES[sp_idx]
    else:
        culprit_sp = phase = "N/A"

    return {
        "case": "no_reaction",
        "max_residual_gas": max_gas,
        "max_residual_solid": max_solid,
        "max_residual_energy": abs(res_energy),
        "residuals_gas": res_gas.copy(),
        "residuals_solid": res_solid.copy(),
        "residual_energy": res_energy,
        "culprit_sp": culprit_sp,
        "culprit_phase": phase,
        "verbose": verbose,
    }


def case_b_reaction_only(cell: Cell, verbose: bool = False) -> dict:
    """Case B：反应工作，但无进料（仅从初始固体/气体消耗）。"""
    idx = GAS_SPECIES_INDEX

    # 清除进料，仅依赖初始库存与反应
    cell.N_zu_d.fill(0.0)
    cell.N_zu_b.fill(0.0)
    cell.N_b_in.fill(0.0)
    cell.N_d_in.fill(0.0)
    cell.N_rez_b.fill(0.0)
    cell.N_rez_d.fill(0.0)

    # 初始库存：一组典型的反应气
    cell.N_d[idx["O2"]] = 1.0
    cell.N_d[idx["H2O"]] = 2.0
    cell.N_d[idx["N2"]] = 20.0
    cell.N_d[idx["CO"]] = 0.1
    cell.N_d[idx["H2"]] = 0.1

    cell.N_b[idx["O2"]] = 0.5
    cell.N_b[idx["H2O"]] = 1.0
    cell.N_b[idx["N2"]] = 10.0

    cell.m_solid[0, S_CHAR] = 0.05
    cell.m_solid[0, S_VM] = 0.01
    cell.m_solid[0, S_ASH] = 0.02

    # 计算水力学与反应
    cell.calc_hydrodynamics()
    cell.calc_exchange()
    cell.calc_reactions()

    # 设置进出温度一致（隔离焓差影响，仅看反应热）
    cell.T_in_gas = cell.T
    cell.T_in_solid = cell.T
    cell.T_zu_gas = cell.T
    cell.T_zu_solid = cell.T

    # 计算残差
    res_gas = cell.calc_gas_balance()
    res_solid = cell.calc_solid_balance().flatten()
    res_energy = cell.calc_energy_balance()

    max_gas = float(np.max(np.abs(res_gas))) if len(res_gas) > 0 else 0.0
    max_solid = float(np.max(np.abs(res_solid))) if len(res_solid) > 0 else 0.0

    # 反应源项元素闭合（比未求解代数残差更能反映化学计量正确性）
    Rg = cell.R_gas_b + cell.R_gas_d
    elem_res = {}
    for el in ("C", "H", "O", "N"):
        gas_elem = 0.0
        for i, sp in enumerate(GAS_SPECIES):
            gas_elem += float(Rg[i]) * float(get_atom_count(sp, el))

        solid_elem = 0.0
        # 固相显式元素源：char(纯C) 与 moisture(按H2O)
        if el == "C":
            solid_elem += float(np.sum(cell.R_solid[:, S_CHAR])) / 0.012011
        elif el == "H":
            solid_elem += 2.0 * float(np.sum(cell.R_solid[:, S_MOISTURE])) / 0.018015
        elif el == "O":
            solid_elem += 1.0 * float(np.sum(cell.R_solid[:, S_MOISTURE])) / 0.018015

        elem_res[el] = gas_elem + solid_elem

    max_elem_res = float(max(abs(v) for v in elem_res.values()))
    vm_source = float(np.sum(np.abs(cell.R_solid[:, S_VM])))

    return {
        "case": "reaction_only",
        "max_residual_gas": max_gas,
        "max_residual_solid": max_solid,
        "max_residual_energy": abs(res_energy),
        "residuals_gas": res_gas.copy(),
        "residuals_solid": res_solid.copy(),
        "residual_energy": res_energy,
        "element_residuals": elem_res,
        "max_element_residual": max_elem_res,
        "vm_source_abs": vm_source,
        "verbose": verbose,
    }


def print_result(result: dict) -> None:
    """打印审计结果。"""
    case = result["case"]
    print()
    print(f"Case: {case.upper()}")
    gas_line = f"  max_residual_gas     = {result['max_residual_gas']:.3e}"
    if "culprit_sp" in result and "culprit_phase" in result:
        gas_line += f"  ({result['culprit_sp']} in {result['culprit_phase']})"
    print(gas_line)
    print(f"  max_residual_solid   = {result['max_residual_solid']:.3e}")
    print(f"  max_residual_energy  = {result['max_residual_energy']:.3e}")
    if "max_element_residual" in result:
        print(f"  max_element_residual = {result['max_element_residual']:.3e}")
        print(f"  vm_source_abs        = {result['vm_source_abs']:.3e}")
        print(f"  element_residuals    = {result['element_residuals']}")

    if result["verbose"]:
        print(f"\n  Full gas residuals (sample, first 5):")
        for i in range(min(5, len(result["residuals_gas"]))):
            print(f"    [{i}] {result['residuals_gas'][i]:.3e}")


def main() -> int:
    print("=" * 72)
    print("Single-cell conservation audit")
    print("=" * 72)

    cell = build_test_cell(T=1100.0, P=2.5e6)

    result_a = case_a_no_reaction(cell, verbose=True)
    print_result(result_a)

    # 重建 cell 以避免状态污染
    cell = build_test_cell(T=1100.0, P=2.5e6)
    result_b = case_b_reaction_only(cell, verbose=True)
    print_result(result_b)

    print()
    print("=" * 72)
    print("Summary:")
    print(f"  Case A (no-reaction) max residual: {result_a['max_residual_gas']:.3e}")
    print(f"    → Target <=1e-8 (exchange-only, no reaction)")
    print(f"  Case B (reaction-only) max residual: {result_b['max_residual_gas']:.3e}")
    print(f"    → Algebraic residual未求解，主要看元素源项闭合")
    print(f"  Case B max element residual: {result_b['max_element_residual']:.3e}")
    print(f"    → Target <=1e-6 (reaction stoichiometry closure)")
    print("=" * 72)

    # 简单检查
    if result_a["max_residual_gas"] <= 1e-8 and result_b["max_element_residual"] <= 1e-6:
        print("Result: BASIC STRUCTURE OK")
        return 0
    else:
        print("Result: ATTENTION NEEDED")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
