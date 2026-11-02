"""诊断：第 1 轮内层阻尼失败点上，差分雅可比矩阵与残差的局部行为。只报告，不改模型。

输入：tools/alternation_diagnostic.py --artifacts 写出的 failure.npz 与 frozen.pkl（只在本地）。
输出：reports/research-solver-path/o2-jacobian.json，含
- 失败点各块的归一残差；
- 沿牛顿方向 x + 2^-i·Δx 的实际归一残差、线性预测与不可评价原因；
- 按变量名末段分组（各气体组分、温度、固相）只沿该组分量移动时，实际与预测的方向导数之差；
- 一个 O₂ 未知量对 CO 行与 O₂ 行的差商随差分步长的变化，以及 ADR 0043 的差分步长。
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

from hamel.solver.residual import evaluate_residual

_ROOT = Path(__file__).resolve().parents[1]
_OUT = _ROOT / "reports" / "research-solver-path" / "o2-jacobian.json"
_CELL = "freeboard[b=6,geo=1].gas_outflow.gas"
_SQRT_EPS = float(np.sqrt(np.finfo(float).eps))


def _group(name: str) -> str:
    return name.rsplit(".", 1)[-1] if ("gas" in name or name.endswith("temperature")) else "solid"


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts")
    args = parser.parse_args(argv)
    folder = Path(args.artifacts)
    data = np.load(folder / "failure.npz", allow_pickle=True)
    with open(folder / "frozen.pkl", "rb") as handle:
        frozen = pickle.load(handle)
    x, dx, J, names = data["x"], data["dx"], data["J"], [str(n) for n in data["names"]]
    blocks = {name: data[name] for name in ("gas", "solid", "energy")}
    refs = {"gas": float(data["n_dot_ref"]), "solid": float(data["m_dot_ref"]), "energy": float(data["q_dot_ref"])}
    weight = np.zeros(x.size)
    for name, rows in blocks.items():
        weight[rows] = 1.0 / refs[name]

    def residual(v):
        return evaluate_residual(v, frozen).values

    merit = lambda values: float(np.linalg.norm(values * weight))
    base = residual(x)
    predicted = J @ dx
    out = {
        "blocks_max": {name: float(np.max(np.abs(base[rows] * weight[rows]))) for name, rows in blocks.items()},
        "merit": merit(base),
        "linear_residual": float(np.max(np.abs((predicted + base) * weight))),
        "line": [],
        "groups": {},
    }
    for index in range(31):
        alpha = 2.0 ** -index
        try:
            values = residual(x + alpha * dx)
        except Exception as exc:
            out["line"].append({"i": index, "error": f"{type(exc).__name__}: {exc}"[:200]})
            continue
        gap = (values - base - alpha * predicted) * weight
        out["line"].append({
            "i": index, "merit": merit(values), "predicted": merit(base + alpha * predicted),
            "worst_row": names[int(np.argmax(np.abs(gap)))],
        })
    alpha = 2.0 ** -27
    groups: dict = {}
    for column, name in enumerate(names):
        groups.setdefault(_group(name), []).append(column)
    for key, columns in groups.items():
        step = np.zeros_like(dx)
        step[columns] = dx[columns]
        if not np.any(step):
            continue
        actual = (residual(x + alpha * step) - base) / alpha * weight
        expected = (J @ step) * weight
        row = int(np.argmax(np.abs(actual - expected)))
        out["groups"][key] = {"max_gap": float(np.abs(actual - expected)[row]), "row": names[row],
                              "actual": float(actual[row]), "predicted": float(expected[row])}
    column, co, o2 = (names.index(f"{_CELL}.{s}") for s in ("O2", "CO", "O2"))
    scale = float(data["scales"][column])
    out["o2_column"] = {
        "name": names[column], "x": float(x[column]), "scale": scale,
        "adr0043_step": _SQRT_EPS * max(abs(float(x[column])), scale),
        "difference_quotients": [],
        "o2_unknowns": {n: float(x[i]) for i, n in enumerate(names) if n.endswith(".O2")},
    }
    for h in (1e-23, 1e-21, 1e-20, 1e-19, 1e-18, 1e-16, 1e-14, 1e-12, 1e-10, out["o2_column"]["adr0043_step"]):
        trial = x.copy()
        trial[column] += h
        quotient = (residual(trial) - base) / h
        out["o2_column"]["difference_quotients"].append({"h": h, "d_co_row": float(quotient[co]), "d_o2_row": float(quotient[o2])})
    out["o2_column"]["jacobian_entries"] = {"d_co_row": float(J[co, column]), "d_o2_row": float(J[o2, column])}
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    _OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(out["blocks_max"]), out["merit"])


if __name__ == "__main__":
    main()
