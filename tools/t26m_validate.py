"""T26m 验证：ADR 0052 修订后的初值构造与基准外层第 1 轮内层。只报告。

复用 T26l 的验证入口（判据、硬超时与报告字段相同），输出目录与标题不同。
判据预先登记：初值构造正常结束，且基准外层第 1 轮内层以 CONVERGED 终止。
"""

from __future__ import annotations

from pathlib import Path

from tools.t26l_validate import main

_OUT = Path(__file__).resolve().parents[1] / "reports" / "t26m-validation"


if __name__ == "__main__":
    main(str(_OUT), "T26m")
