"""T31 反应路径账的拆分与聚合。只读 ``_apply`` 的记录，不改模型、不改判定。

C.0 的口径：每条记录存带符号的 q = ν_i × extent（mol/s）。拆分永远在保方向的实际记录层做
（``max(q, 0)``／``max(-q, 0)``），聚合只做同层相加；不得先把同一反应跨格／跨相的 q 合并再拆。
"""
from __future__ import annotations

ACCOUNT_SPECIES = ("C", "O2", "H2O", "H2", "CH4", "CO", "CO2")
TAR_SPECIES = ("C6H6", "C10H8")
PHASES = {"bubble": "气泡相", "suspension": "悬浮相", "free": "自由空间", "hover": "自由空间悬停"}
TINY = 1.0e-12


def account(records) -> dict:
    """按 C.0 给出每个物种的总消耗／总生成／带符号净值，并逐格、逐反应分列。

    返回值分别是 ``species``、``cells``、``reactions``。逐格键为 ``(zone, cell)``；
    逐反应的消耗与生成同样在记录层拆开再加（不是先合并 q 再取符号）。
    """
    species = {name: _zero() for name in ACCOUNT_SPECIES}
    cells: dict = {}
    reactions: dict = {name: {} for name in ACCOUNT_SPECIES}
    for record in records:
        for name, value in record["q"].items():
            q = float(value)
            if name not in species:
                continue
            cell = cells.setdefault((record["zone"], int(record["cell"])), {name: _zero() for name in ACCOUNT_SPECIES})
            reaction = reactions[name].setdefault(record["reaction"], _zero())
            _add(species[name], q)
            _add(cell[name], q)
            _add(reaction, q)
    return {"species": species, "cells": cells, "reactions": reactions}


def tar_totals(records) -> dict:
    """焦油代表物（C6H6、C10H8）的消耗／生成／净值，只用于核对 R10／R11 的分项。"""
    totals = {name: _zero() for name in TAR_SPECIES}
    for record in records:
        for name, value in record["q"].items():
            if name in totals:
                _add(totals[name], float(value))
    return totals


def _zero() -> dict:
    return {"consume": 0.0, "generate": 0.0, "net": 0.0}


def _add(target: dict, q: float) -> None:
    target["generate"] += max(q, 0.0)
    target["consume"] += max(-q, 0.0)
    target["net"] += q


def share(part: float, total: float) -> float | None:
    """同符号份额。分母小于 1e-12 时返回 None，报告记“分母可忽略”。"""
    return None if abs(total) < TINY else part / total


def registry(records) -> dict:
    """C.8 的占比登记表。分母一律取自 ``account`` 的同符号总量。"""
    book = account(records)
    water, oxygen = book["species"]["H2O"], book["species"]["O2"]
    reactions = book["reactions"]

    def consumed(name: str, reaction: str) -> float:
        return float(reactions[name].get(reaction, _zero())["consume"])

    def generated(name: str, reaction: str) -> float:
        return float(reactions[name].get(reaction, _zero())["generate"])

    return {
        "water_consume": water["consume"], "oxygen_consume": oxygen["consume"],
        "water_net": water["net"],
        "r2_water_share": share(consumed("H2O", "R2"), water["consume"]),
        "r5_oxygen_share": share(consumed("O2", "R5"), oxygen["consume"]),
        "tar_water_share": share(consumed("H2O", "R9") + consumed("H2O", "R11"), water["consume"]),
        "r6_water_generate": generated("H2O", "R6"),
    }


def aggregated(records) -> list:
    """按 (zone, cell, phase, reaction, phi, tar) 相加，用于 Markdown 表；q 的零项不列。"""
    merged: dict = {}
    for record in records:
        key = (record["zone"], int(record["cell"]), record["phase"], record["reaction"],
               record["phi"], record["tar"])
        row = merged.setdefault(key, {})
        for name, value in record["q"].items():
            row[name] = row.get(name, 0.0) + float(value)
    return [
        {"zone": key[0], "cell": key[1], "phase": key[2], "reaction": key[3],
         "phi": key[4], "tar": key[5], "q": {name: value for name, value in row.items() if value != 0.0}}
        for key, row in sorted(merged.items(), key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3]))
    ]


def markdown_ledger(ledger) -> list:
    """driver_report 的“反应路径账”一节。完整记录留在 .scratch，这里只给聚合表。"""
    records = (ledger or {}).get("records") or []
    if not records:
        return ["## 反应路径账（T31）", "", "这次运行没有反应路径账（未收敛或未开启记录）。", ""]
    book = account(records)
    reg = registry(records)
    lines = ["## 反应路径账（T31）", "",
             "逐记录 q = ν_i×extent（mol/s，产物为正）。消耗／生成在记录层拆开再加，净值另列。",
             "完整逐记录表在 `.scratch/t31/run.json`（不入库），本节只给聚合表。", "",
             "| 物种 | 消耗 mol/s | 生成 mol/s | 净值 mol/s |", "|---|---|---|---|"]
    for name in ACCOUNT_SPECIES:
        row = book["species"][name]
        lines.append(f"| {name} | {row['consume']:.6e} | {row['generate']:.6e} | {row['net']:.6e} |")
    lines.append("")
    for zone, label in (("bed", "床区"), ("freeboard", "自由空间")):
        if not any(key[0] == zone for key in book["cells"]):
            continue
        lines += [f"### {label}逐格（消耗／生成／净值，mol/s）", "",
                  "| 格 | 物种 | 消耗 | 生成 | 净值 |", "|---|---|---|---|---|"]
        for (cell_zone, cell), row in sorted(book["cells"].items(), key=lambda item: item[0][1]):
            if cell_zone != zone:
                continue
            for name in ACCOUNT_SPECIES:
                body = row[name]
                lines.append(f"| {cell} | {name} | {body['consume']:.6e} | {body['generate']:.6e} | {body['net']:.6e} |")
        lines.append("")
    lines += ["### 按 (区, 格, 相, 反应) 聚合（mol/s）", "",
              "| 区 | 格 | 相 | 反应 | φ | 焦油 | " + " | ".join(ACCOUNT_SPECIES) + " |",
              "|---|---|---|---|---|---|" + "---|" * len(ACCOUNT_SPECIES)]
    for row in aggregated(records):
        cells = [f"{row['q'][name]:.6e}" if name in row["q"] else "" for name in ACCOUNT_SPECIES]
        phi = "" if row["phi"] is None else f"{row['phi']:.9f}"
        lines.append(f"| {row['zone']} | {row['cell']} | {row['phase']} | {row['reaction']} | {phi} | {row['tar'] or ''} | "
                     + " | ".join(cells) + " |")
    lines.append("")
    lines += ["### 占比登记（C.8）", "", "| 项 | 占比 | 绝对量 mol/s |", "|---|---|---|"]
    def shown(value):
        return "分母可忽略" if value is None else f"{value:.6f}"

    lines.append(f"| R2 占总耗水 | {shown(reg['r2_water_share'])} | {reg['water_consume']:.6e} |")
    lines.append(f"| R5 占总耗氧 | {shown(reg['r5_oxygen_share'])} | {reg['oxygen_consume']:.6e} |")
    lines.append(f"| R9+R11 占总耗水 | {shown(reg['tar_water_share'])} | {reg['water_consume']:.6e} |")
    lines += ["", f"R6 产水 {reg['r6_water_generate']:.6e} mol/s；水净值 {reg['water_net']:.6e} mol/s。",
              "本表只说明路径占比，不构成因果结论。", ""]
    return lines
