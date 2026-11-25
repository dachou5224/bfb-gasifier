"""T31 反应路径账的拆分与聚合。只读 ``_apply`` 的记录，不改模型、不改判定。

C.0 的口径：每条记录存带符号的 q = ν_i × extent（mol/s）。拆分永远在保方向的实际记录层做
（``max(q, 0)``／``max(-q, 0)``），聚合只做同层相加；不得先把同一反应跨格／跨相的 q 合并再拆。
"""
from __future__ import annotations

import numpy as np

from hamel.properties.data import GAS_SPECIES

ACCOUNT_SPECIES = ("C", "O2", "H2O", "H2", "CH4", "CO", "CO2")
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
    """C.8 的占比登记表。分母一律取自 ``account`` 的同符号总量。

    ``*_contribution`` 是该路径自己的同符号贡献（mol/s），与分母分开；分母另记在 ``*_total``。
    """
    book = account(records)
    water, oxygen = book["species"]["H2O"], book["species"]["O2"]
    reactions = book["reactions"]

    def consumed(name: str, reaction: str) -> float:
        return float(reactions[name].get(reaction, _zero())["consume"])

    def generated(name: str, reaction: str) -> float:
        return float(reactions[name].get(reaction, _zero())["generate"])

    r2 = consumed("H2O", "R2")
    r5 = consumed("O2", "R5")
    tar = consumed("H2O", "R9") + consumed("H2O", "R11")
    return {
        "water_consume": water["consume"], "oxygen_consume": oxygen["consume"],
        "water_net": water["net"],
        "r2_water_share": share(r2, water["consume"]), "r2_water_contribution": r2,
        "r5_oxygen_share": share(r5, oxygen["consume"]), "r5_oxygen_contribution": r5,
        "tar_water_share": share(tar, water["consume"]), "tar_water_contribution": tar,
        "r6_water_generate": generated("H2O", "R6"),
    }


def gas_residual_rows(frozen, state, result, records) -> list:
    """逐格逐相的气相摩尔残差与尺度，供报告 C.3 与 D.2 共用同一份计算。

    原式：入流 + 回流 + 反应 − 出流 ± 相间交换。交换项的相符号只在这里加一次
    （悬浮相减、气泡相加），报告不再取负。尺度按本卡 D.2 定义逐物种取
    max(|入流|, |出流|, Σ|ν×进度|, |交换|, 1e-12)。
    """
    count = len(GAS_SPECIES)
    magnitudes: dict = {}
    for record in records:
        key = (record["zone"], int(record["cell"]), record["phase"])
        bucket = magnitudes.setdefault(key, {})
        for name, q in record["q"].items():
            bucket[name] = bucket.get(name, 0.0) + abs(float(q))
    exchange = result.gas_exchange_mol_s
    rows = []

    def add(zone, cell, phase, sink_key, inflow, outflow, swap, reaction):
        magnitude = np.array([magnitudes.get(sink_key, {}).get(name, 0.0) for name in GAS_SPECIES])
        scale = np.maximum.reduce([
            np.abs(inflow), np.abs(outflow), magnitude, np.abs(swap), np.full(count, TINY),
        ])
        rows.append({
            "zone": zone, "cell": int(cell), "phase": phase,
            "inflow": inflow, "outflow": outflow, "exchange": swap,
            "reaction": reaction, "residual": inflow + reaction - outflow + swap,
            "scale": scale, "magnitude": magnitude,
        })

    for cell in range(frozen.n_bed):
        for index, phase in ((0, "suspension"), (1, "bubble")):
            inflow = state.bed_gas_mol_s[cell + 1, index] if cell + 1 < frozen.n_bed else np.zeros(count)
            if cell == frozen.n_bed - 1:
                share_in = frozen.feed_to_suspension if index == 0 else 1.0 - frozen.feed_to_suspension
                inflow = inflow + frozen.primary_feed_mol_s * share_in
            reaction = result.reaction_bed[cell, index].copy()
            if index == 0:
                reaction = reaction + frozen.pyro_bed_mol_s[cell]
            swap = -exchange[cell] if index == 0 else exchange[cell]
            add("bed", cell, phase, ("bed", cell, phase), inflow, state.bed_gas_mol_s[cell, index], swap, reaction)
    for cell in range(frozen.n_freeboard):
        inflow = state.freeboard_gas_mol_s[cell + 1] if cell + 1 < frozen.n_freeboard else np.zeros(count)
        if cell == frozen.n_freeboard - 1:
            inflow = inflow + state.bed_gas_mol_s[0, 0] + state.bed_gas_mol_s[0, 1]
        if cell == frozen.secondary_balance_index:
            inflow = inflow + frozen.secondary_feed_mol_s
        reaction = result.reaction_free[cell] + frozen.pyro_freeboard_mol_s[cell]
        add("freeboard", cell, "free", ("freeboard", cell, "free"), inflow,
            state.freeboard_gas_mol_s[cell], np.zeros(count), reaction)
    return rows


def aggregated(records) -> list:
    """按 (zone, cell, phase, reaction, phi, tar) 相加，用于 Markdown 表；q 的零项不列。"""
    merged: dict = {}
    for record in records:
        key = (record["zone"], int(record["cell"]), record["phase"], record["reaction"],
               record["phi"], record["tar"])
        row = merged.setdefault(key, {})
        for name, value in record["q"].items():
            if value != 0.0:
                row[name] = row.get(name, 0.0) + float(value)
    rows = []
    for key, row in sorted(merged.items(), key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3])):
        q = {name: value for name, value in row.items() if value != 0.0}
        if q:  # 结构零的调用（例如 extent = 0）不单列
            rows.append({"zone": key[0], "cell": key[1], "phase": key[2], "reaction": key[3],
                         "phi": key[4], "tar": key[5], "q": q})
    return rows


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
    lines += ["### 按 (区, 格, 相, 反应) 聚合（mol/s）", "",
              "| 区 | 格 | 相 | 反应 | φ | 焦油 | " + " | ".join(ACCOUNT_SPECIES) + " |",
              "|---|---|---|---|---|---|" + "---|" * len(ACCOUNT_SPECIES)]
    for row in aggregated(records):
        cells = [f"{row['q'][name]:.6e}" if name in row["q"] else "" for name in ACCOUNT_SPECIES]
        phi = "" if row["phi"] is None else f"{row['phi']:.9f}"
        lines.append(f"| {row['zone']} | {row['cell']} | {row['phase']} | {row['reaction']} | {phi} | {row['tar'] or ''} | "
                     + " | ".join(cells) + " |")
    lines.append("")
    lines += ["### 占比登记（C.8）", "", "| 项 | 占比 | 该路径贡献 mol/s | 分母 mol/s |", "|---|---|---|---|"]

    def shown(value):
        return "分母可忽略" if value is None else f"{value:.6f}"

    lines.append(f"| R2 占总耗水 | {shown(reg['r2_water_share'])} | {reg['r2_water_contribution']:.6e} | {reg['water_consume']:.6e} |")
    lines.append(f"| R5 占总耗氧 | {shown(reg['r5_oxygen_share'])} | {reg['r5_oxygen_contribution']:.6e} | {reg['oxygen_consume']:.6e} |")
    lines.append(f"| R9+R11 占总耗水 | {shown(reg['tar_water_share'])} | {reg['tar_water_contribution']:.6e} | {reg['water_consume']:.6e} |")
    lines += ["", f"R6 产水 {reg['r6_water_generate']:.6e} mol/s；水净值 {reg['water_net']:.6e} mol/s。",
              "“分母”是同一物种的同类总量（`消耗_H2O`／`消耗_O2`），只作占比分母，不是该路径的量。",
              "本表只说明路径占比，不构成因果结论。", ""]
    return lines
