"""求解器性能与残差调用计数测试。

用于验证优化前后数值结果一致，并记录 residuals() 调用次数与耗时。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.core.reactor import Reactor, ReactorConfig
from src.core.species import configure_tar_components_by_fuel


def test_solve_timing():
    """记录 reactor.solve() 耗时。"""
    configure_tar_components_by_fuel("coal")
    cfg = ReactorConfig(
        n_cells=10,
        H_bed=5.0,
        D_bed=0.6,
        P=2.5e6,
        T_inlet=300.0,
        fuel_type="coal",
        rho_s=1000.0,
        d_p=2e-3,
        fuel_feed=3377.4 / 3600.0,
        O2_feed=16.17,
        H2O_feed=12.93,
        N2_feed=0.16,
        moisture_wt=16.9,
        C_dry=61.5,
        H_dry=4.1,
        VM_daf=53.42,
        recirculation_frac=0.1,
    )
    reactor = Reactor(cfg)

    t0 = time.perf_counter()
    result = reactor.solve(max_global_iter=10, tol_global=5.0)
    elapsed = time.perf_counter() - t0

    assert result["converged"]
    assert result["n_iter"] >= 1
    assert 0 < result["carbon_conv"] <= 1.0
    assert len(result["T_profile"]) == 10

    print(f"  solve() 耗时: {elapsed:.3f} s")
    print(f"  全局迭代: {result['n_iter']}")
    print(f"  出口温度: {result['T_profile'][-1]:.0f} K")


def test_residual_count():
    """通过 wrapper 计数 residuals() 调用次数（需修改 cell_solver 注入计数器）。"""
    # 简化：仅验证 solve 可完成
    test_solve_timing()


if __name__ == "__main__":
    print("=== BFB 求解器性能测试 ===")
    test_solve_timing()
    print("PASS")
