"""T26m 验证：ADR 0052 修订后的初值构造与基准外层第 1 轮内层。只报告。

复用 T26l 的验证入口（判据、硬超时与报告字段相同），输出目录与标题不同。
判据预先登记：初值构造正常结束，且基准外层第 1 轮内层以 CONVERGED 终止。

``--report`` 不重跑求解，只用已有的 ``validation.json`` 重写 ``validation.md``。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from tools import t26l_validate
from tools.t26l_validate import main

_OUT = Path(__file__).resolve().parents[1] / "reports" / "t26m-validation"


if __name__ == "__main__":
    if "--report" in sys.argv:
        payload = json.loads((_OUT / "validation.json").read_text(encoding="utf-8"))
        (_OUT / "validation.md").write_text(t26l_validate._markdown(payload, "T26m"), encoding="utf-8")
    else:
        main(str(_OUT), "T26m")
