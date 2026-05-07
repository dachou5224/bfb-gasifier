"""LU 工况下 R12 开关影响对比。

运行两组配置（仅 R12 开关不同），输出关键 KPI：
- 出口温度 T_exit
- 碳转化率 carbon_conv
- 干基 CO/CO2/H2/CH4

用法：
    /Users/liuzhen/AI-projects/bfb-gasifier/.venv/bin/python scripts/audit_r12_toggle_lu.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.core.reactor import Reactor
from tests.validation_case_utils import (
    build_phase1_htw_lu_reactor_config,
    load_validation_case_node,
    CASE_LU_VALIDATION_KEY,
    PHASE1_HTW_LU_SOLVE_KWARGS,
)


def run_once(enable_r12: bool) -> dict:
    cfg = build_phase1_htw_lu_reactor_config()
    cfg.enable_r12 = enable_r12

    reactor = Reactor(cfg)
    result = reactor.solve(**PHASE1_HTW_LU_SOLVE_KWARGS)

    y_dry = result["exit_gas_dry"]
    return {
        "enable_r12": enable_r12,
        "converged": result.get("converged"),
        "n_iter": result.get("n_iter"),
        "T_exit": float(result["T_profile"][-1]),
        "carbon_conv": float(result["carbon_conv"]),
        "CO_dry": float(y_dry.get("CO", 0.0)),
        "CO2_dry": float(y_dry.get("CO2", 0.0)),
        "H2_dry": float(y_dry.get("H2", 0.0)),
        "CH4_dry": float(y_dry.get("CH4", 0.0)),
    }


def main() -> int:
    ref = load_validation_case_node(CASE_LU_VALIDATION_KEY)["outputs"]
    ref_dry = ref.get("exit_gas_dry_mol_frac", {})
    t_ref = float(ref.get("exit_temperature_K", 1100.0))
    xc_ref = float(ref.get("carbon_conversion_pct", 95.0)) / 100.0
    print("=" * 72)
    print("R12 toggle audit on HTW LU case")
    print("=" * 72)
    print(
        "Reference: "
        f"T={ref.get('exit_temperature_K')} K, "
        f"Xc={ref.get('carbon_conversion_pct')}%, "
        f"dry={ref.get('exit_gas_dry_mol_frac')}"
    )
    print()

    off = run_once(False)
    on = run_once(True)

    def _row(name: str, a: float, b: float, fmt: str = ".4f") -> None:
        diff = b - a
        print(f"{name:16s} OFF={format(a, fmt):>10s}  ON={format(b, fmt):>10s}  Δ(ON-OFF)={format(diff, fmt):>10s}")

    print("Convergence:")
    print(f"  OFF: converged={off['converged']} n_iter={off['n_iter']}")
    print(f"  ON : converged={on['converged']} n_iter={on['n_iter']}")
    print()

    print("KPIs:")
    _row("T_exit [K]", off["T_exit"], on["T_exit"], ".1f")
    _row("carbon_conv", off["carbon_conv"], on["carbon_conv"], ".4f")
    _row("CO_dry", off["CO_dry"], on["CO_dry"])
    _row("CO2_dry", off["CO2_dry"], on["CO2_dry"])
    _row("H2_dry", off["H2_dry"], on["H2_dry"])
    _row("CH4_dry", off["CH4_dry"], on["CH4_dry"])

    def _rel_err(sim: float, ref_val: float) -> float:
        return abs(sim - ref_val) / max(abs(ref_val), 1e-12)

    print()
    print("Relative error vs validation (OFF / ON):")
    print(f"  T_exit      : {_rel_err(off['T_exit'], t_ref):.2%} / {_rel_err(on['T_exit'], t_ref):.2%}  (ref={t_ref:.1f} K)")
    print(f"  carbon_conv : {_rel_err(off['carbon_conv'], xc_ref):.2%} / {_rel_err(on['carbon_conv'], xc_ref):.2%}  (ref={xc_ref:.4f})")
    for sp in ("CO", "CO2", "H2", "CH4"):
        ref_sp = float(ref_dry.get(sp, 0.0) or 0.0)
        if ref_sp <= 0.0:
            continue
        off_v = off[f"{sp}_dry"]
        on_v = on[f"{sp}_dry"]
        print(f"  {sp:<11s}: {_rel_err(off_v, ref_sp):.2%} / {_rel_err(on_v, ref_sp):.2%}  (ref={ref_sp:.4f})")

    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
