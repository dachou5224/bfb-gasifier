"""固定 s=0.16 上，NH3/H2S 若计入自由空间混合物，气速、密度和轨迹差多少。

不改求解器，不把这个差值写进残差。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from src.core.species import MOLECULAR_WEIGHT
from src.solvers.wirsum_1997_bed_freeboard import BedFreeboardGeometry
from src.solvers.wirsum_1997_named_f import outer_P_named
from src.solvers.wirsum_1997_named_f_p2 import (
    LAUNCH_DT_S,
    LAUNCH_N_BINS_DEFAULT,
    LAUNCH_T_MAX_S,
    gas_velocity_and_density,
    integrate_launch_distribution,
)
from src.solvers.wesseling_case_contract import wesseling12_input_ledger
from src.solvers.wesseling_two_zone import q_on_wesseling1_path

_OLD = Path(__file__).resolve().parents[1] / "s1_s3_review_2026-09-30" / "t3_corrected_path.json"
_OUT = Path(__file__).resolve().parent / "u4_inputs.json"
_G_PER_KG = 1.0e-3  # [kg/g]


def _plain(obj):
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def _augmented(n_major: float, mw_major: float, n_nh3: float, n_h2s: float) -> dict:
    mass = (
        float(n_major) * float(mw_major)
        + float(n_nh3) * float(MOLECULAR_WEIGHT["NH3"]) * _G_PER_KG
        + float(n_h2s) * float(MOLECULAR_WEIGHT["H2S"]) * _G_PER_KG
    )
    n_all = float(n_major) + float(n_nh3) + float(n_h2s)
    mw_all = mass / n_all
    return {
        "n_mol_s": n_all,
        "MW_kg_per_mol": mw_all,
        "n_ratio": n_all / float(n_major),
        "MW_ratio": mw_all / float(mw_major),
        "y_NH3": float(n_nh3) / n_all,
        "y_H2S": float(n_h2s) / n_all,
        "major_mole_fraction_ratio": float(n_major) / n_all,
    }


def main() -> None:
    doc = json.loads(_OLD.read_text())
    saved = [row for row in doc["toward_022"] if abs(row["s"] - 0.16) < 1.0e-12 and not row["joint_pass"]][-1]
    x = np.asarray(saved["x"], dtype=np.float64)
    nq = q_on_wesseling1_path(0.16)
    q = outer_P_named(x, nq)
    sp = q.energy.species
    n_nh3 = float(sp.n_NH3_bypass_mol_s)
    n_h2s = float(sp.n_H2S_bypass_mol_s)
    b = q.reference
    n_fb = x[16:23]
    u7, rho7, n7, mw7 = gas_velocity_and_density(
        n_fb, np.zeros(7), T=float(x[-1]), P=float(b.P), area=float(b.area),
    )
    bed = _augmented(float(q.hydro.n_tot), float(q.hydro.MW_mix), n_nh3, n_h2s)
    free = _augmented(n7, mw7, n_nh3, n_h2s)
    base = []
    alt = []
    for i, d in enumerate(np.asarray(b.d_p, dtype=np.float64)):
        geo = BedFreeboardGeometry(
            area=float(b.area),
            H_freeboard=float(b.H_freeboard),
            d_p=float(d),
            rho_p=float(b.rho_p),
            eps_d=float(b.eps_d),
            eps_b=float(b.eps_b),
            u_b=float(q.hydro.u_b),
            d_b=float(b.d_b),
            u_g=u7 * free["n_ratio"],
            rho_g=rho7 * free["MW_ratio"],
            mu_g=float(b.mu_g),
            phi_s=float(b.phi_s),
        )
        dist = integrate_launch_distribution(
            geo,
            n_bins=LAUNCH_N_BINS_DEFAULT,
            dt=LAUNCH_DT_S,
            t_max=LAUNCH_T_MAX_S,
            n_cells=2,
            continue_budget_s=float(b.trajectory_continue_s),
            extra_physical_s=float(b.trajectory_extra_physical_s),
            max_work_steps=int(b.trajectory_max_work_steps),
            wall_s=float(b.trajectory_wall_s),
        )
        base.append({
            "d_p_m": float(d),
            "w_leave": float(q.hydro.w_leave[i]),
            "w_return": float(q.hydro.w_return[i]),
            "w_undecided": float(q.hydro.w_undecided[i]),
            "tau_return_s": float(q.hydro.tau_return_s[i]),
            "tau_leave_s": float(q.hydro.tau_leave_s[i]),
        })
        alt.append({
            "d_p_m": float(d),
            "w_leave": float(dist["w_leave"]),
            "w_return": float(dist["w_return"]),
            "w_undecided": float(dist["w_undecided"]),
            "tau_return_s": float(dist["tau_return_s"]),
            "tau_leave_s": float(dist["tau_leave_s"]),
        })
        print(i, base[-1]["tau_return_s"], alt[-1]["tau_return_s"], flush=True)
    report = {
        "ledger": wesseling12_input_ledger(),
        "state": "t3 toward_022 s=0.16 merit_reject, not a root",
        "bypass_mol_s": {"NH3": n_nh3, "H2S": n_h2s, "sum": n_nh3 + n_h2s},
        "freeboard_seven_mol_s": n7,
        "bypass_share_of_freeboard_plus_bypass": (n_nh3 + n_h2s) / (n7 + n_nh3 + n_h2s),
        "current": {
            "freeboard_u_g_m_s": u7,
            "freeboard_rho_kg_m3": rho7,
            "freeboard_MW_kg_per_mol": mw7,
            "bed_u_g_m_s": float(q.hydro.u_g),
            "bed_rho_kg_m3": float(q.hydro.rho_g),
        },
        "if_bypass_joined_the_freeboard_mixture": free,
        "if_the_same_moles_joined_the_bed_total": bed,
        "trajectory_current": base,
        "trajectory_if_freeboard_mixture_included_bypass": alt,
        "written_into_solver": False,
    }
    _OUT.write_text(json.dumps(_plain(report)), encoding="utf-8")
    print("share", report["bypass_share_of_freeboard_plus_bypass"], "u_ratio", free["n_ratio"], flush=True)


if __name__ == "__main__":
    main()
