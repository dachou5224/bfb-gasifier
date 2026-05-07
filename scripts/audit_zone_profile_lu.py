"""LU 工况分区诊断：热解区 / 燃烧区 / 气化区。

输出每个 cell 的关键状态，定位：
1) O2 耗尽位置
2) VM 释放完成位置
3) 水分释放完成位置

用于排查 NR 主链中的分区耦合与收敛异常。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase1_htw_lu_global_nr_reactor_config,
)
from src.core.reactor import Reactor
from src.core.cell import S_CHAR, S_VM, S_MOISTURE
from src.core.species import GAS_SPECIES_INDEX as idx
from src.solvers.nr_indexing import solid_component_indices_for_nr


def main() -> int:
    cfg = build_phase1_htw_lu_global_nr_reactor_config()
    r = Reactor(cfg)
    # 与 LU global-NR 主验证保持同口径。
    res = r.solve(**PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS)

    bot = r.cells[0]
    top = r.cells[-1]

    m_vm_feed = float(np.sum(np.maximum(bot.m_solid_zu[:, S_VM], 0.0)))
    m_moist_feed = float(np.sum(np.maximum(bot.m_solid_zu[:, S_MOISTURE], 0.0)))
    m_char_feed = float(np.sum(np.maximum(bot.m_solid_zu[:, S_CHAR], 0.0)))

    o2_in_ref = float(cfg.O2_feed)

    print("=" * 108)
    print("LU zone profile audit")
    print("=" * 108)
    residual = float(res.get("residual", float("nan")))
    comp_names = {S_CHAR: "char", S_VM: "vm", S_MOISTURE: "moisture", 3: "ash"}
    unknown_layout: dict[str, list[str]] = {}
    for cell in r._solver_cells_for_nr():
        ctype = str(getattr(cell, "cell_type", "bed"))
        if ctype in unknown_layout:
            continue
        idxs = solid_component_indices_for_nr(cell).tolist()
        unknown_layout[ctype] = [comp_names.get(int(i), str(i)) for i in idxs]
    vm_moist_in_unknown = any(
        ("vm" in comps) or ("moisture" in comps)
        for comps in unknown_layout.values()
    )
    nr_counts = dict(res.get("nr_counts") or {})
    ls_evals = int(nr_counts.get("line_search_evaluations", 0))
    ls_backtracks = int(nr_counts.get("line_search_backtracks", 0))
    ls_failures = int(nr_counts.get("line_search_failures", 0))
    jac_zero_cols = int(nr_counts.get("jacobian_zero_cols_last", 0))
    jac_zero_rows = int(nr_counts.get("jacobian_zero_rows_last", 0))
    clip_last = (res.get("nr_clip_history") or [{}])[-1]
    print(
        f"converged={res['converged']} n_iter={res['n_iter']} residual={residual:.3e} "
        f"rms_scaled_final={res.get('rms_scaled_final', float('nan')):.3e} "
        f"outer={res.get('nr_outer_iters')} budget={res.get('nr_inner_budget_used')}/{res.get('nr_inner_budget_total')} "
        f"T_exit={res['T_profile'][-1]:.1f}K Xc={res['carbon_conv']:.4f}"
    )
    print(
        f"line-search: evals={ls_evals} backtracks={ls_backtracks} failures={ls_failures} "
        f"last_best_trial_rms={float(clip_last.get('best_trial_rms_scaled', float('nan'))):.3e} "
        f"last_best_trial_lambda={clip_last.get('best_trial_lambda')}"
    )
    print(
        f"jacobian(last): zero_cols={jac_zero_cols} zero_rows={jac_zero_rows} "
        f"nnz={nr_counts.get('jacobian_nnz_last', 'n/a')}"
    )
    print(
        f"nr solid unknown layout: {unknown_layout} "
        f"(vm/moisture_in_unknown={vm_moist_in_unknown})"
    )
    print(
        f"fresh feed: O2={cfg.O2_feed:.3f} mol/s, VM={m_vm_feed:.5f} kg/s, "
        f"Moist={m_moist_feed:.5f} kg/s, Char={m_char_feed:.5f} kg/s"
    )
    print("-" * 108)
    print(
        "cell  h(m)   T(K)   O2_out   O2_cons   m_VM     m_Moist  m_Char   "
        "VM_rel(H2Oeq)  VM_out/in  Moist_out/in"
    )
    print("-" * 108)

    o2_exhaust_cell = None
    vm_done_cell = None
    moist_done_cell = None

    vm_frac_prev = 1.0

    for i, c in enumerate(r.cells):
        h = (i + 0.5) * cfg.H_bed / cfg.n_cells

        c.calc_hydrodynamics()
        tau = c.geo.dh / max(c.u_mf, 1e-3)
        c.compute_vorabrechnung(tau)
        c.calc_exchange()
        c.calc_reactions()

        N_out = np.maximum(c.N_d + c.N_b, 0.0)
        o2_out = float(N_out[idx["O2"]])
        o2_cons = float(max(-(c.R_gas_d[idx["O2"]] + c.R_gas_b[idx["O2"]]), 0.0))

        m_vm = float(np.sum(np.maximum(c.m_solid[:, S_VM], 0.0)))
        m_moist = float(np.sum(np.maximum(c.m_solid[:, S_MOISTURE], 0.0)))
        m_char = float(np.sum(np.maximum(c.m_solid[:, S_CHAR], 0.0)))

        m_vm_in_local = float(np.sum(np.maximum(c.m_solid_zu[:, S_VM] + c.m_solid_rez[:, S_VM] + c.m_solid_in[:, S_VM], 0.0)))
        m_moist_in_local = float(np.sum(np.maximum(c.m_solid_zu[:, S_MOISTURE] + c.m_solid_rez[:, S_MOISTURE] + c.m_solid_in[:, S_MOISTURE], 0.0)))
        vm_frac = m_vm / max(m_vm_in_local, 1e-12)
        moist_frac = m_moist / max(m_moist_in_local, 1e-12)

        # 用 VM 气源中的 H2O 项作为“热解/干燥活跃度”代理（仅诊断）
        vm_rel_proxy = float(max(c._vm_gas_source_cache[idx["H2O"]], 0.0))

        print(
            f"{i:>3d}  {h:>4.1f}  {c.T:>6.1f}  {o2_out:>7.3f}  {o2_cons:>8.3f}  "
            f"{m_vm:>7.5f}  {m_moist:>8.5f}  {m_char:>7.5f}  {vm_rel_proxy:>12.4f}  "
            f"{vm_frac:>7.3f}   {moist_frac:>9.3f}"
        )

        if o2_exhaust_cell is None and o2_out <= 0.01 * max(o2_in_ref, 1e-9):
            o2_exhaust_cell = i

        if vm_done_cell is None and m_vm <= 0.01 * max(m_vm_feed, 1e-12):
            vm_done_cell = i

        if moist_done_cell is None and m_moist <= 0.01 * max(m_moist_feed, 1e-12):
            moist_done_cell = i

        vm_frac_prev = vm_frac

    print("-" * 108)
    print("milestones:")
    print(f"  O2 exhausted cell: {o2_exhaust_cell}")
    print(f"  VM released (<=1%) cell: {vm_done_cell}")
    print(f"  Moisture released (<=1%) cell: {moist_done_cell}")

    # 简单分区建议
    if vm_done_cell is not None:
        pyro_zone = f"cell 0 ~ {vm_done_cell}"
    else:
        pyro_zone = "未在床层内完成"

    if o2_exhaust_cell is not None:
        comb_zone = f"cell 0 ~ {o2_exhaust_cell}"
        gasif_zone = f"cell {o2_exhaust_cell} ~ {len(r.cells)-1}"
    else:
        comb_zone = "未在床层内耗尽"
        gasif_zone = "氧气未耗尽，气化区不明确"

    print("zone inference:")
    print(f"  pyrolysis zone: {pyro_zone}")
    print(f"  combustion zone: {comb_zone}")
    print(f"  gasification zone: {gasif_zone}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
