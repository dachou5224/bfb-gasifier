#!/usr/bin/env python3
"""LU 工况下 freeboard 初始夹带固体依据审计。

目的：
1. 用 Hamel 本地提取中的 bed→freeboard 初始抛射通量公式估算 `F0r/F0w/F0`
2. 检查当前 shared global NR 床顶是否存在不可忽略的 entrained-solids basis
3. 对照当前代码口径：freeboard 已有 entrained-solids surrogate，并提供 impl/thesis 双标签反应诊断

参考提取：
- docs/chat-record.md 问答 48
  F0r = d_p_bar * rho_p * (1 - eps_d) * (3 * eps_b * u_b / d_b)
  F0w = zeta_w * f_w * eps_b * u_b * rho_p * (1 - eps_d)
- 其中 f_w = 0.25, zeta_w = 0.40
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.cell import S_ASH, S_CHAR
from src.core.reactor import Reactor
from _nr_monitor import print_nr_monitor, solve_with_nr_monitor
from tests.validation_case_utils import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    build_phase2_htw_lu_freeboard_reactor_config,
)


F_W = 0.25
ZETA_W = 0.40


def main() -> int:
    cfg = build_phase2_htw_lu_freeboard_reactor_config()
    reactor = Reactor(cfg)
    result, monitor = solve_with_nr_monitor(
        reactor,
        PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
        check_x0=True,
    )

    top = reactor.cells[-1]
    area = math.pi * (cfg.D_bed ** 2) / 4.0
    d_p_bar = float(np.average(top.solid.d_p_classes, weights=top.solid.mass_fractions))
    n_size = int(len(top.solid.d_p_classes))
    rho_p = float(top.solid.rho_s)
    eps_b = float(top.eps_b)
    eps_d_void = float(top.eps_d_voidage)
    u_b = float(top.u_b)
    d_b = max(float(top.d_b), 1e-12)

    f0_r = d_p_bar * rho_p * max(1.0 - eps_d_void, 0.0) * (3.0 * eps_b * u_b / d_b)
    f0_w = ZETA_W * F_W * eps_b * u_b * rho_p * max(1.0 - eps_d_void, 0.0)
    f0 = f0_r + f0_w
    m_dot_eject = f0 * area

    m_top_char = float(np.sum(np.maximum(top.m_solid[:, S_CHAR], 0.0)))
    m_top_ash = float(np.sum(np.maximum(top.m_solid[:, S_ASH], 0.0)))
    m_top_char_ash = m_top_char + m_top_ash

    print("=" * 132)
    print("LU freeboard entrained-solids basis audit")
    print("=" * 132)
    print_nr_monitor(monitor, prefix="NR monitor")
    print(
        f"freeboard_active={result['freeboard_active']} "
        f"bed_exit_T={result['bed_T_profile'][-1]:.1f}K reactor_exit_T={result['reactor_exit_T']:.1f}K "
        f"n_iter={result['n_iter']} rms={result.get('rms_scaled_final', float('nan')):.3e}"
    )
    print("-" * 132)
    print(
        f"d_p_bar={d_p_bar:.4e}m rho_p={rho_p:.1f}kg/m3 eps_b={eps_b:.4f} "
        f"eps_d_void={eps_d_void:.4f} u_b={u_b:.4f}m/s d_b={d_b:.4f}m area={area:.4f}m2 "
        f"n_size_classes={n_size}"
    )
    print(
        f"F0r={f0_r:.4f} kg/m2/s  F0w={f0_w:.4f} kg/m2/s  "
        f"F0={f0:.4f} kg/m2/s  m_dot_eject0={m_dot_eject:.4f} kg/s"
    )
    print("-" * 132)
    print(
        f"top_out_char={m_top_char:.4f} kg/s  top_out_ash={m_top_ash:.4f} kg/s  "
        f"top_out_char_ash={m_top_char_ash:.4f} kg/s"
    )
    if m_top_char_ash > 1e-12:
        print(
            f"eject0 / top_out_char_ash = {m_dot_eject / m_top_char_ash:.3f}"
        )
    if result["freeboard_entrained_catalyst_density_kg_m3"]:
        rho_cat_max = max(result["freeboard_entrained_catalyst_density_kg_m3"])
        rho_cat_avg = sum(result["freeboard_entrained_catalyst_density_kg_m3"]) / len(result["freeboard_entrained_catalyst_density_kg_m3"])
        diag_impl = result.get("freeboard_reaction_diag_impl", result["freeboard_reaction_diag"])
        diag_thesis = result.get("freeboard_reaction_diag_thesis", diag_impl)
        r11d_sum_impl = sum(abs(v) for v in diag_impl.get("R11d", []))
        r11_sum_thesis = sum(abs(v) for v in diag_thesis.get("R11", []))
        carry_min = min(result["freeboard_carry_ratio_profile"])
        carry_max = max(result["freeboard_carry_ratio_profile"])
        u0_avg = sum(result["freeboard_u0_profile_m_s"]) / len(result["freeboard_u0_profile_m_s"])
        ugb_avg = sum(result["freeboard_u_gb_profile_m_s"]) / len(result["freeboard_u_gb_profile_m_s"])
        up_avg = sum(result["freeboard_u_p_mean_profile_m_s"]) / len(result["freeboard_u_p_mean_profile_m_s"])
        ut_avg = sum(result["freeboard_u_t_mean_profile_m_s"]) / len(result["freeboard_u_t_mean_profile_m_s"])
        print(
            f"freeboard rho_cat[max]={rho_cat_max:.4e} kg/m3  "
            f"rho_cat[avg]={rho_cat_avg:.4e} kg/m3  "
            f"sum|R11d_impl|={r11d_sum_impl:.4e}  sum|R11_thesis|={r11_sum_thesis:.4e}"
        )
        print(
            f"u0[avg]={u0_avg:.4f} m/s  u_gb[avg]={ugb_avg:.4f} m/s  "
            f"u_p[avg]={up_avg:.4f} m/s  u_t[avg]={ut_avg:.4f} m/s  "
            f"carry_ratio[min/max]={carry_min:.3f}/{carry_max:.3f}"
        )
    print(
        f"entrained_exit_char={result['freeboard_entrained_exit_char_kg_s']:.4f} kg/s  "
        f"entrained_exit_ash={result['freeboard_entrained_exit_ash_kg_s']:.4f} kg/s"
    )
    print(
        f"return_char={result['freeboard_entrained_return_char_kg_s']:.4f} kg/s  "
        f"return_ash={result['freeboard_entrained_return_ash_kg_s']:.4f} kg/s"
    )
    print(
        f"cyclone_capture_char={result['freeboard_cyclone_capture_char_kg_s']:.4f} kg/s  "
        f"cyclone_capture_ash={result['freeboard_cyclone_capture_ash_kg_s']:.4f} kg/s  "
        f"recycle_candidate={result['freeboard_cyclone_recycle_candidate_char_ash_kg_s']:.4f} kg/s"
    )
    print("-" * 132)
    print("interpretation:")
    print("  - 若 m_dot_eject0 显著非零，则 Hamel 口径下 freeboard 不应先验视为纯气相。")
    print("  - 当前代码已包含 entrained-solids surrogate（轨迹 + 返回 + 旋风捕集候选）。")
    print("  - `R11d_impl` 与 `R11_thesis` 同时输出，用于消除编号歧义。")
    print("  - 仍不是 Hamel 原始 Verbindungsmatrix 的完整拓扑实现。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
