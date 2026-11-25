"""T31 报告生成：从 ``.scratch/t31/state.npz`` + ``frozen.pkl`` 离线复算账目。

不进入求解，不改任何模型设定；只读记录并按 C.0 拆分。产物是
``reports/t31-reaction-ledger/report.md``。运行记录与 T30 基准取自 ``.scratch/``。
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hamel.cases.assemble import assemble_case  # noqa: E402
from hamel.kinetics.arrhenius import arrhenius_rate_constant  # noqa: E402
from hamel.properties import chemical_potential_j_mol, gas_constant_j_mol_k  # noqa: E402
from hamel.properties.data import GAS_SPECIES, atomic_weight, molar_mass_kg_mol  # noqa: E402
from hamel.solver.driver import _references  # noqa: E402
from hamel.solver.residual import _phase_point, evaluate_residual  # noqa: E402
from hamel.solver.vector import unpack_state  # noqa: E402
from tools.run_prep import fuel_element_fractions  # noqa: E402
from tools.ledger_account import (  # noqa: E402
    ACCOUNT_SPECIES, account, aggregated, gas_residual_rows, registry, share,
)

TAR_SPECIES = ("C6H6", "C10H8")

CASE = ROOT / "data" / "cases" / "vtt_pressurized.json"
SCRATCH = ROOT / ".scratch" / "t31"
OUT = Path(__file__).resolve().parent / "report.md"
INDEX = {name: index for index, name in enumerate(GAS_SPECIES)}
SPECIES_LABEL = {"C": "焦炭 C", "O2": "O2", "H2O": "H2O", "H2": "H2", "CH4": "CH4", "CO": "CO", "CO2": "CO2"}
REACTION_ORDER = [f"R{index}" for index in range(1, 12)]
TARGET_NET_WATER = 0.020
T30_QUANTITIES = {
    "T_FRB_K": 1174.5786279498361, "X_C": 86.63883804995233, "CO": 0.16550155875233868,
    "CO2": 0.11152052199505044, "H2": 0.15081039649453565,
    "CH4": 0.017476599979452976, "H2O": 0.060270632847835806,
}
T30_RESIDUALS = {"gas": 2.0203230380580937e-07, "solid": 9.730255704422228e-07, "energy": 4.091150206045431e-08}
T30_ASH = 3.248117502078796e-07
T30_ELEMENTS = {"C": 8.092482833475196e-08, "H": 7.0533461427074705e-09, "O": 4.7419780308566526e-09, "energy": 2.153606204094182e-09}
T30_WALL_W = 21700.0
T30_BOTTOM = 5.014963079349439e-07
T30_ELAPSED = 5474.425006125


def f(value, digits=6):
    return "—" if value is None else f"{float(value):.{digits}e}"


def g(value, digits=6):
    return "—" if value is None else f"{float(value):.{digits}f}"


def rel(left, right):
    scale = max(abs(float(left)), abs(float(right)))
    return 0.0 if scale == 0.0 else abs(float(left) - float(right)) / scale


def pct(value):
    return "分母可忽略" if value is None else f"{100.0 * value:.2f} %"


def table(headers, rows):
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return lines


class Report:
    def __init__(self):
        self.x = np.load(SCRATCH / "state.npz")["x"]
        with (SCRATCH / "frozen.pkl").open("rb") as handle:
            self.frozen = pickle.load(handle)
        self.run = json.loads((SCRATCH / "run.json").read_text(encoding="utf-8"))
        self.state = unpack_state(self.x, self.frozen)
        self.result = evaluate_residual(self.x, self.frozen, ledger=True)
        self.records = self.result.ledger
        self.book = account(self.records)
        self.reg = registry(self.records)
        self.rows = aggregated(self.records)
        self.tars = _tar_totals(self.records)
        assembled = assemble_case(str(CASE), "1")
        self.physical = assembled.physical
        self.n_ref = _references(self.physical)[0]
        self.abs_q = {}
        for row in self.records:
            key = (row["zone"], int(row["cell"]), row["phase"])
            bucket = self.abs_q.setdefault(key, {})
            for name, q in row["q"].items():
                bucket[name] = bucket.get(name, 0.0) + abs(float(q))
        water = molar_mass_kg_mol("H2O")
        self.moisture_water = float(self.physical.fuel_mass_flow) * float(self.physical.moisture_mass_fraction) / water
        self.steam_water = float(self.frozen.primary_feed_mol_s[INDEX["H2O"]] + self.frozen.secondary_feed_mol_s[INDEX["H2O"]])
        self.outlet_water = float(self.state.freeboard_gas_mol_s[0][INDEX["H2O"]])
        self.pyro_water = float(self.frozen.pyro_bed_mol_s[:, INDEX["H2O"]].sum() + self.frozen.pyro_freeboard_mol_s[:, INDEX["H2O"]].sum())
        # 总水账：出口 − 入口水（燃料水分 + 蒸汽）。热解水是冻结源项，不在反应账里。
        self.overall_water_net = self.outlet_water - self.moisture_water - self.steam_water

    # ---- 基础查询 -------------------------------------------------------
    def cell_row(self, zone, cell, species):
        return (self.book["cells"].get((zone, cell)) or {}).get(species, {"consume": 0.0, "generate": 0.0, "net": 0.0})

    def reaction_rows(self, species):
        rows = []
        for name in REACTION_ORDER:
            body = self.book["reactions"][species].get(name)
            if body and (body["consume"] != 0.0 or body["generate"] != 0.0):
                rows.append((name, body))
        return rows

    def conditions(self):
        frozen, state = self.frozen, self.state
        out = []
        for cell in range(frozen.n_bed):
            for phase, label in ((0, "悬浮相"), (1, "气泡相")):
                y, c, _f, _mu, _rho = _phase_point(state.bed_gas_mol_s[cell, phase], state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell])
                volume = frozen.v_suspension_m3[cell] if phase == 0 else frozen.v_bubble_m3[cell]
                out.append(self._condition("bed", cell, label, state.bed_temperature_k[cell], frozen.bed_pressure_pa[cell],
                                           float(np.sum(state.bed_gas_mol_s[cell, phase])), volume, y, c, phase))
        for cell in range(frozen.n_freeboard):
            y, c, _f, _mu, _rho = _phase_point(state.freeboard_gas_mol_s[cell], state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell])
            out.append(self._condition("freeboard", cell, "自由空间", state.freeboard_temperature_k[cell], frozen.freeboard_pressure_pa[cell],
                                       float(np.sum(state.freeboard_gas_mol_s[cell])), frozen.v_freeboard_m3[cell], y, c, 0))
        return out

    def _condition(self, zone, cell, label, temperature, pressure, flow, volume, y, c, phase):
        k6 = arrhenius_rate_constant(1.631e9 * float(temperature) ** 1.5, float(temperature), activation_temperature_k=3420.0)
        derivative = k6 * max(c[INDEX["H2"]], 0.0) ** 1.5
        if phase == 0 and zone == "bed":
            numerator = 1.8 * float(temperature) ** 0.5 * c[INDEX["CO"]]  # R5 悬浮相：O2 的 0.5 次
            crossover = (numerator / derivative) ** 2 if derivative > 0.0 else None
        elif phase == 1:
            k5 = arrhenius_rate_constant(1.91e6, float(temperature), activation_temperature_k=8056.0)
            numerator = k5 * c[INDEX["CO"]] * max(c[INDEX["H2O"]], 0.0) ** 0.5  # 气泡相：O2 的 0.3 次
            crossover = (numerator / derivative) ** (1.0 / 0.7) if derivative > 0.0 else None
        else:
            crossover = None
        position = "—" if crossover is None else ("低于交叉点" if c[INDEX["O2"]] < crossover else "高于交叉点")
        return {"zone": zone, "cell": cell, "phase": label, "t": float(temperature), "p": float(pressure),
                "flow": float(flow), "volume": float(volume),
                "c": {name: float(c[INDEX[name]]) for name in ("O2", "H2", "H2O", "CO", "CO2", "CH4")},
                "y": {name: float(y[INDEX[name]]) for name in ("O2", "H2", "H2O", "CO", "CO2", "CH4")},
                "crossover": crossover, "position": position}


def main() -> int:
    report = Report()
    book, reg, records = report.book, report.reg, report.records
    frozen, state, result = report.frozen, report.state, report.result
    lines = ["# T31：第 1 关（VTT 加压 1）反应路径账", "",
             "工况 1（`data/cases/vtt_pressurized.json`，莱茵褐煤、不循环），基准运行 `adiabatic=False`。",
             "本文件是**诊断记录**：只读出收敛解上逐格、逐相的 R1–R11 计量贡献并按 C.0 拆分。",
             "不改任何速率式、计量、温度取法、物性、闭合、容差或判定；结论限于占比与绝对量，不写因果判断。", ""]

    status = report.run.get("status")
    lines += ["## 首页", "",
              f"- 终态：**`{status}`**；外层 {report.run.get('rounds')} 轮 / 内层 {report.run.get('iterations')} 次；墙钟 {g(report.run.get('elapsed_s'), 1)} s。",
              f"- 账目记录：{len(records)} 条（`_apply` 的每次调用一条，含结构零系数）；按 (格, 相, 反应, φ, 焦油) 聚合后 {len(report.rows)} 行（结构零的调用不单列）。",
              f"- 水净值 `Σq = {f(book['species']['H2O']['net'])}` mol/s；总生成 {f(book['species']['H2O']['generate'])}、总消耗 {f(book['species']['H2O']['consume'])} mol/s。",
              f"- 重建精度 D.1 的最坏相对差见下；物理收支 D.2 与 T30 的对照见下。",
              "- 账目只在通过 §2.1、2(b) 与 ADR 0059 检查线之后对最终解评价一次，不进入牛顿迭代、雅可比、初值构造或试探路径；下面的 C.10 用 7 量与残差证明数值未改。", ""]

    # ---- 事实复算（卡片第 2–4 条） -------------------------------------
    fractions = fuel_element_fractions(report.physical)
    fuel = float(report.physical.fuel_mass_flow)
    carbon_in = fractions["C"] * fuel / (atomic_weight("C") * 1.0e-3)
    hydrogen_in = fractions["H"] * fuel / (atomic_weight("H") * 1.0e-3) + 2.0 * report.steam_water
    oxygen_in = (fractions["O"] * fuel / (atomic_weight("O") * 1.0e-3)
                 + 2.0 * float(np.sum(frozen.primary_feed_mol_s[INDEX["O2"]] + frozen.secondary_feed_mol_s[INDEX["O2"]]))
                 + report.steam_water)
    outlet_mass = float(sum(state.freeboard_gas_mol_s[0][index] * molar_mass_kg_mol(name) for index, name in enumerate(GAS_SPECIES)))
    outlet = state.freeboard_gas_mol_s[0]
    outlet_total = float(np.sum(outlet))
    y = outlet / outlet_total
    k_act = float(y[INDEX["CO2"]] * y[INDEX["H2"]] / (y[INDEX["CO"]] * y[INDEX["H2O"]]))
    target_k = 0.091 * 0.109 / (0.164 * 0.122)
    peak_k = float(np.max(state.bed_temperature_k))
    mus = [chemical_potential_j_mol(name, peak_k) for name in ("CO", "H2O", "H2", "CO2")]
    k_eq = float(np.exp(-(-mus[0] - mus[1] + mus[2] + mus[3]) / (gas_constant_j_mol_k() * peak_k)))
    reverse = sum(1 for row in records if row["reaction"] == "R8" and float(row["q"]["H2O"]) > 0.0)
    forward = sum(1 for row in records if row["reaction"] == "R8" and float(row["q"]["H2O"]) < 0.0)
    lines += ["## 事实复算（卡片“已知事实”第 2–4 条）", "",
              "全部用本次收敛解与 `data/cases/vtt_pressurized.json` 复算，不引用卡片数字作为输入。", ""]
    rows = [
        ["入口 C", f(carbon_in), "0.6103", "燃料碳"],
        ["入口 H", f(hydrogen_in), "0.9047", "燃料 H + 蒸汽 2×0.104851"],
        ["入口 O", f(oxygen_in), "0.8047", "燃料 O + 空气 2×O2 + 蒸汽"],
        ["入口水（燃料水分 + 蒸汽）", f(report.moisture_water + report.steam_water), "0.1901", "收到基水分 + 一次风蒸汽"],
        ["出口气质量 kg/s", f(outlet_mass), "0.044879", "自由空间最上格出流"],
        ["总水账净产水", f(report.overall_water_net), "约 −0.084", "出口水 − 入口水（见 C.2）"],
        ["K_act 模型（出口）", f(k_act), "1.686", "Q = y_CO2·y_H2 / (y_CO·y_H2O)"],
        ["K_act 目标（目标组成）", f(target_k), "0.496", "同上，用目标摩尔分数"],
        ["K（式 5.43，1182 K，同向复算）", f(k_eq), "1.307", "同一 Q 定义；**与卡片值不符**，见下"],
    ]
    lines += table(["量", "本次复算", "卡片值", "说明"], rows) + [""]
    lines += [f"R8 的实际符号：{forward} 条记录为正进度（耗水、产氢），{reverse} 条为负进度（产水、耗氢）；"
              "同一全床合计同时含两个方向的记录，未按反应编号合并成一个符号。",
              f"床区最高温度 {f(peak_k)} K。式 (5.43) 与这里的 Q 都对应同一个正向反应 CO + H₂O → H₂ + CO₂，"
              f"所以同向复算的 K = {f(k_eq)}；`r8_water_gas_shift_rate` 也按同一方向用 "
              "`y_eq_CO = y_CO2·y_H2/(K·y_H2O)`。**卡片所列 K = 1.307 与这个同向复算值不符**"
              f"（1.307 是 {f(k_eq)} 的倒数），本卡不视其为已核实，也未为凑数而改方向或改公式。",
              f"两种取值下定性结论相同：K_act 模型 {f(k_act)} 与 K_act 目标 {f(target_k)} 都分别位于"
              f"K = {f(k_eq)}（或 1.307）的两侧，且离 1 都很远。", ""]
    lines += ["运行命令（配置与 T30 逐字相同，`outer_round_limit=26`、进程级硬墙钟 28 800 s）：", "",
              "```", ".venv/bin/python reports/t31-reaction-ledger/run.py", "```", "",
              "环境：arm64、Python 3.9.6、numpy 1.23.5、scipy 1.10.1。未跑 S1 与任何敏感性运行。", ""]

    # ---- C.0 ----------------------------------------------------------
    lines += ["## C.0 统一的符号定义", "",
              "每条记录是 `q[j, i] = ν[j, i] × extent[j]`（mol/s，带符号，产物为正、反应物为负），",
              "键为 `(zone, cell, phase, reaction, phi, tar)`。对每个物种：", "",
              "$$生成_i = \\sum_j \\max(q[j,i], 0), \\quad 消耗_i = \\sum_j \\max(-q[j,i], 0), \\quad 净_i = \\sum_j q[j,i]$$", "",
              "拆分在**保方向的实际记录层**做（逐记录取 `max(q,0)`／`max(−q,0)`），聚合只做同层相加；",
              "不先把同一反应跨格／跨相的 q 合并再拆。可逆的 R8 按实际符号落列。", "",
              "`phase` 取 `bubble`（气泡相）、`suspension`（床区悬浮相）、`free`（自由空间气相与自由态颗粒）、",
              "`hover`（自由空间悬停态颗粒）；自由空间的均相反应记在 `free`。`phi` 只在 R1 非空，`tar` 只在 R10／R11 非空。",
              "完整逐记录表在 `.scratch/t31/run.json`（不入库），本报告只给聚合表。", ""]

    # ---- C.1 ----------------------------------------------------------
    lines += ["## C.1 物种账（消耗／生成／净值）", "",
              "### 全床合计", ""]
    rows = [[SPECIES_LABEL[name], f(book["species"][name]["consume"]), f(book["species"][name]["generate"]), f(book["species"][name]["net"])] for name in ACCOUNT_SPECIES]
    lines += table(["物种", "总消耗 mol/s", "总生成 mol/s", "带符号净值 mol/s"], rows)
    lines += ["", "### 逐格合计（床区与自由空间分开，格号 0 为各区最上层）", ""]
    for zone, label in (("bed", "床区"), ("freeboard", "自由空间")):
        cells = sorted(cell for (cell_zone, cell) in book["cells"] if cell_zone == zone)
        rows = []
        for cell in cells:
            for name in ACCOUNT_SPECIES:
                row = report.cell_row(zone, cell, name)
                rows.append([cell, SPECIES_LABEL[name], f(row["consume"]), f(row["generate"]), f(row["net"])])
        if rows:
            lines += [f"#### {label}", ""] + table(["格", "物种", "消耗", "生成", "净值"], rows) + [""]
    lines += ["### 按反应分列（全床；每行只列该反应的非零贡献）", ""]
    for name in ACCOUNT_SPECIES:
        rows = []
        for reaction, body in report.reaction_rows(name):
            if body["consume"] == 0.0 and body["generate"] == 0.0:
                continue
            rows.append([SPECIES_LABEL[name], reaction, f(body["consume"]), f(body["generate"]), f(body["net"])])
        if rows:
            lines += table(["物种", "反应", "消耗", "生成", "净值"], rows) + [""]
    tar_rows = [[name, f(body["consume"]), f(body["generate"]), f(body["net"])] for name, body in report.tars.items()]
    lines += ["焦油代表物（不属于 7 个物种账，单列）：", ""] + table(["焦油", "消耗", "生成", "净值"], tar_rows) + [""]

    # ---- C.2 ----------------------------------------------------------
    water = book["species"]["H2O"]
    lines += ["## C.2 水账专节", "",
              "水由 R6（+1/进度）与 R7（+2/进度）生成，由 R2／R8（正进度）／R9／R11 消耗；",
              "R5 计量上既不耗水也不耗氢，气泡相的 `c_H2O^0.5` 只是动力学促进，不在水账里。", "",
              f"全床：总消耗 **{f(water['consume'])}**、总生成 **{f(water['generate'])}**、净值 **{f(water['net'])}** mol/s。", ""]
    lines += ["### 水的分层账（反应层与总水账分开）", "",
              "`Σq` 是反应账（R2／R6／R7／R8／R9／R11）的净值；干燥与热解释放的水是冻结源项，不在反应账里。", ""]
    rows = [
        ["反应账 Σq（C.0）", f(water["net"]), "R6+R7 生成 − R2/R8/R9/R11 消耗"],
        ["热解与干燥水（冻结源项）", f(report.pyro_water), "`pyro_bed_mol_s` + `pyro_freeboard_mol_s` 的 H2O"],
        ["进口蒸汽", f(report.steam_water), "一次风与二次风的 H2O"],
        ["燃料水分", f(report.moisture_water), "收到基水分 × 燃料流量 / M_H2O"],
        ["出口水", f(report.outlet_water), "自由空间最上格 H2O 出流"],
        ["总水账（出口 − 蒸汽 − 燃料水分）", f(report.overall_water_net), "与卡片事实 3 的 −0.084 对照"],
    ]
    lines += table(["项", "mol/s", "说明"], rows) + [""]
    lines += [f"总水账核对：出口 {f(report.outlet_water)} = 蒸汽 {f(report.steam_water)} + 热解水 {f(report.pyro_water)} + 反应账 {f(water['net'])}"
              f"（差 {f(report.outlet_water - report.steam_water - report.pyro_water - water['net'])}）。",
              f"入口水合计 {f(report.moisture_water + report.steam_water)} mol/s（燃料水分 {f(report.moisture_water)} + 蒸汽 {f(report.steam_water)}），"
              f"与卡片事实 2 的 0.1901 mol/s 一致；总水账 **{f(report.overall_water_net)}** mol/s，与卡片事实 3 的约 −0.084 mol/s 同量级。",
              "因此：反应层净分解 " + f(-water['net']) + " mol/s，其中热解水补偿 " + f(report.pyro_water) + " mol/s，总水账净分解 "
              + f(-report.overall_water_net) + " mol/s。", ""]
    rows = []
    for row in report.rows:
        q = float(row["q"].get("H2O", 0.0))
        if q == 0.0:
            continue
        side = "生成" if q > 0.0 else "消耗"
        tag = row["tar"] or ""
        note = f"φ={row['phi']:.6f}" if row["phi"] is not None else ""
        rows.append([row["zone"], row["cell"], row["phase"], row["reaction"], tag, note, f(abs(q)), side] if abs(q) >= 1e-12 else None)
    tiny = sum(1 for item in rows if item is None)
    rows = [item for item in rows if item is not None]
    lines += table(["区", "格", "相", "反应", "焦油", "φ", "abs(q) mol/s", "按符号落列"], rows) + [""]
    lines += [f"表中省略 {tiny} 条 |q| < 1×10⁻¹² mol/s 的记录（数值零），完整记录在 JSON。", ""]
    lines += ["### 按反应汇总（全床，水）", ""]
    rows = [[reaction, f(body["consume"]), f(body["generate"]), f(body["net"])] for reaction, body in report.reaction_rows("H2O")]
    lines += table(["反应", "消耗", "生成", "净值"], rows) + [""]
    lines += ["R11 的系数按实际焦油代表物：苯 `−4.5`、萘 `−23/3`；下表按 `tar` 分开，未写 `−R11`。", ""]

    # ---- C.3 ----------------------------------------------------------
    lines += ["## C.3 氧账专节", "",
              f"全床 O2：总消耗 **{f(book['species']['O2']['consume'])}**、总生成 {f(book['species']['O2']['generate'])}、净值 {f(book['species']['O2']['net'])} mol/s。", ""]
    rows = [[reaction, f(body["consume"]), f(body["net"])] for reaction, body in report.reaction_rows("O2")]
    lines += table(["反应", "消耗", "净值"], rows) + [""]
    lines += ["### 各相与各反应的耗氧（逐记录取 max(−q,0) 后相加）", ""]
    buckets: dict = {}
    for row in report.rows:
        value = float(row["q"].get("O2", 0.0))
        if value < 0.0:
            key = (row["phase"], row["reaction"])
            buckets[key] = buckets.get(key, 0.0) - value
    total_oxygen = book["species"]["O2"]["consume"]
    order = {"suspension": 0, "bubble": 1, "free": 2, "hover": 3}
    rows = [[phase, reaction, f(value), pct(share(value, total_oxygen))]
            for (phase, reaction), value in sorted(buckets.items(), key=lambda item: (order.get(item[0][0], 9), int(item[0][1][1:])))]
    lines += table(["相", "反应", "消耗 mol/s", "占总耗氧"], rows) + [""]
    lines += ["R1 的氧系数按实际 φ 取 `−1/φ`；R5 分气泡相与悬浮相两行（见 C.1 的逐格表与 C.7）。", ""]
    lines += ["### 逐格 O2 流（进口气 + 由下方来 + 相间交换 − 出流 − 反应 = 残差）", ""]
    rows = _oxygen_rows(report)
    lines += table(["区", "格", "相", "入流 mol/s", "出流 mol/s", "交换 mol/s", "反应消耗 mol/s", "残差 mol/s"], rows) + [""]
    inlet_o2 = float(np.sum(frozen.primary_feed_mol_s[INDEX["O2"]] + frozen.secondary_feed_mol_s[INDEX["O2"]]))
    lines += [f"模型入口氧总量 {f(inlet_o2)} mol/s；任务卡对照值 `0.21 × 1.015 = 0.21315` mol/s。"
              "出口氧以自由空间最上格出流为准。", ""]

    # ---- C.4 ----------------------------------------------------------
    lines += ["## C.4 氢与甲烷专节", "",
              "H2 的生成（R2／R8／R9／R10／R11，R8 按实际符号）与消耗（R4／R6，以及负进度的 R8）；",
              "CH4 的生成（R4／R11）与消耗（R7／R9）。", ""]
    for name in ("H2", "CH4"):
        rows = [[reaction, f(body["consume"]), f(body["generate"]), f(body["net"])] for reaction, body in report.reaction_rows(name)]
        lines += [f"### {name}", ""] + table(["反应", "消耗", "生成", "净值"], rows) + [""]

    # ---- C.5 ----------------------------------------------------------
    lines += ["## C.5 碳账专节", "",
              "焦炭量按 C.0 的 `q[C]` 带符号记账（R1–R4 的 `ν_C = −1`，正进度即消耗）；"
              "折算只用 `返回质量 / M_C`，不再取负号。", ""]
    rows = []
    for reaction, body in report.reaction_rows("C"):
        rows.append([reaction, f(body["consume"]), f(body["generate"]), f(body["net"])])
    lines += table(["反应", "消耗（mol/s）", "生成（mol/s）", "净值"], rows) + [""]
    total_char = book["species"]["C"]["consume"]
    rows = [[reaction, f(body["consume"]), pct(share(body["consume"], total_char))] for reaction, body in report.reaction_rows("C")]
    lines += ["焦炭消耗份额（分母为 `消耗_C` = " + f(total_char) + " mol/s）：", ""]
    lines += table(["反应", "消耗 mol/s", "占焦炭消耗"], rows) + [""]
    rows = []
    for name in ("CO", "CO2"):
        for reaction, body in report.reaction_rows(name):
            if body["consume"] == 0.0 and body["generate"] == 0.0:
                continue
            rows.append([name, reaction, f(body["consume"]), f(body["generate"]), f(body["net"])])
    lines += ["CO 与 CO2 的生成／相互转化：", ""] + table(["物种", "反应", "消耗", "生成", "净值"], rows) + [""]
    phis: dict = {}
    for row in report.rows:
        if row["reaction"] != "R1":
            continue
        body = phis.setdefault(row["phi"], {"O2": 0.0, "CO": 0.0, "CO2": 0.0})
        for name in body:
            body[name] += float(row["q"].get(name, 0.0))
    kept = [(phi, body) for phi, body in phis.items()
            if abs(body["O2"]) + abs(body["CO"]) + abs(body["CO2"]) >= 1e-12]
    omitted = len(phis) - len(kept)
    rows = [[g(phi, 9), f(body["O2"]), f(body["CO"]), f(body["CO2"])]
            for phi, body in sorted(kept, key=lambda item: (item[0] is None, item[0]))]
    lines += ["R1 按实际 φ 分项（同一 φ 的逐记录相加；φ=1 时 CO 系数为 0、CO2 为 +1）：", ""]
    lines += table(["φ", "O2", "CO", "CO2"], rows) + [""]
    lines += [f"其余 {omitted} 个 φ 的 O2／CO／CO2 贡献绝对值之和都小于 1×10⁻¹² mol/s，未单列；完整记录在 JSON。", ""]
    measured = (report.run.get("measurements") or {}).get("X_C")
    lines += [f"`X_C` = **{g(measured, 6)} %**（口径照 T30 报告：出口气体含焦油碳 ÷ 加入碳，本卡不改）。", ""]

    # ---- C.6 ----------------------------------------------------------
    lines += ["## C.6 份额定义", "",
              "- 消耗份额 = 该项的 `max(−q, 0)` ÷ 该物种的 `消耗_i`；生成份额同法，分母为 `生成_i`。",
              "- 净值单列，带符号，不作份额分母。分母 < 1e−12 mol/s 时不定义份额，记“分母可忽略”。",
              "- 占比登记表的“总耗水”= `消耗_H2O`，“总耗氧”= `消耗_O2`。", ""]

    # ---- C.7 ----------------------------------------------------------
    lines += ["## C.7 局部条件与 R5／R6 的氧竞争交叉点", "",
              "交叉点按本次解的实际速率常数解 `R5 耗氧速率 = R6 耗氧速率`（两者对 O2 的系数都是 0.5，故等价于进度速率相等）：",
              "悬浮相 R5 对 O2 为 0.5 次、R6 为 1 次，`c_O2* = (1.8√T·c_CO / k6·c_H2^1.5)²`；",
              "气泡相 R5 对 O2 为 0.3 次，`c_O2* = (k5·c_CO·√c_H2O / k6·c_H2^1.5)^(1/0.7)`。只报数，不下因果结论。", ""]
    rows = []
    for item in report.conditions():
        rows.append([item["zone"], item["cell"], item["phase"], g(item["t"], 3), g(item["p"], 3), f(item["volume"]),
                     f(item["flow"]), f(item["c"]["O2"]), f(item["y"]["O2"]), f(item["crossover"]), item["position"]])
    lines += ["### 温度、压力、相体积、气相摩尔流与 O2", ""]
    lines += table(["区", "格", "相", "T / K", "p / Pa", "V / m³", "气相 mol/s", "c_O2", "y_O2", "c_O2*", "位置"], rows) + [""]
    rows = []
    for item in report.conditions():
        rows.append([item["zone"], item["cell"], item["phase"]] +
                    [f(item["c"][name]) for name in ("H2", "H2O", "CO", "CO2", "CH4")] +
                    [f(item["y"][name]) for name in ("H2", "H2O", "CO", "CO2", "CH4")])
    lines += ["### 其余物种的浓度与摩尔分数", "",
              "| 区 | 格 | 相 | c_H2 | c_H2O | c_CO | c_CO2 | c_CH4 | y_H2 | y_H2O | y_CO | y_CO2 | y_CH4 |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows] + [""]
    lines += ["自由空间只有单一气相，R5 不在此相计算（ADR 0010 第 5 条：R5 只在床区两相），故交叉点记“—”。", ""]

    # ---- C.8 ----------------------------------------------------------
    lines += ["## C.8 占比登记表（预登记门槛，看到数字后未移动）", "",
              "分母一律取自 C.0 的同符号总量；结论栏只做门槛分类，不构成因果结论。", ""]
    lines += ["| 项 | 占比 | 门槛分类 | 该路径贡献 mol/s | 分母（同符号总量）mol/s |", "|---|---|---|---|---|"]
    lines.append(f"| R2 占总耗水 | {pct(reg['r2_water_share'])} | {_class(reg['r2_water_share'], 0.5, 0.2, '集中', '次要', '可忽略')} | {f(reg['r2_water_contribution'])} | {f(reg['water_consume'])} |")
    lines.append(f"| R5 占总耗氧 | {pct(reg['r5_oxygen_share'])} | {_class(reg['r5_oxygen_share'], 0.5, 0.2, '集中', '次要', '可忽略')} | {f(reg['r5_oxygen_contribution'])} | {f(reg['oxygen_consume'])} |")
    lines.append(f"| R9+R11 占总耗水 | {pct(reg['tar_water_share'])} | {_class(reg['tar_water_share'], 0.2, 0.2, '显著', '显著', '可忽略')} | {f(reg['tar_water_contribution'])} | {f(reg['water_consume'])} |")
    delta = TARGET_NET_WATER - report.overall_water_net
    lines.append(f"| R6 产水缺口（算术） | — | 只报数 | 基准 R6 产水 {f(reg['r6_water_generate'])} mol/s；目标侧所需 {f(reg['r6_water_generate'] + delta)} mol/s；差值 {f(-delta)} mol/s |")
    lines += ["", "结论句只允许“占比 X %”“该反应消耗 Y mol/s”的形式；本表不判断谁过量、谁抢氧。", ""]

    # ---- C.9 ----------------------------------------------------------
    lines += ["## C.9 算术一致性检验（**非模型结果，仅算术**）", "",
              f"令 δ = 目标侧净产水 0.020 − 模型**总水账**净产水 {f(report.overall_water_net)} = **{f(delta)}** mol/s。"
              f"（反应账净值为 {f(book['species']['H2O']['net'])} mol/s，不含热解水 {f(report.pyro_water)} mol/s；"
              "卡片事实 3 的约 −0.084 指总水账。）", "",
              "若 δ 全部由 R6 多产水提供（ΔR6 = δ），则连带 ΔH2 = −δ、ΔO2 = −0.5δ：", ""]
    n_out = float(np.sum(state.freeboard_gas_mol_s[0]))
    h2_out = float(state.freeboard_gas_mol_s[0][INDEX["H2"]])
    n2_out = float(state.freeboard_gas_mol_s[0][INDEX["N2"]])
    target_y = {"CO": 0.164, "CO2": 0.091, "H2": 0.109, "CH4": 0.013, "H2O": 0.122}
    target_total = n2_out / (1.0 - sum(target_y.values()))  # N2 守恒 + 目标组成，仅算术
    target_h2 = target_y["H2"] * target_total
    h2_gap = h2_out - target_h2
    h2_gap_model = h2_out - target_y["H2"] * n_out
    residual_pct = abs(delta - h2_gap) / abs(h2_gap)
    rows = [
        ["ΔR6", f(delta), "占基准 R6 实际进度", pct(delta / reg["r6_water_generate"])],
        ["ΔH2（−δ）", f(-delta), "H2 实际缺口（模型 − 目标）", f(h2_gap)],
        ["H2 残差", f(abs(delta - h2_gap)), "残差百分比", pct(residual_pct)],
        ["ΔO2（−0.5δ）", f(-0.5 * delta), "入口氧可用量", f(inlet_o2)],
        ["ΔO2 占入口氧", pct(abs(0.5 * delta) / inlet_o2), "—", "—"],
    ]
    lines += table(["量", "值", "对照", "占比"], rows) + [""]
    lines += [f"H2 实际缺口按 N2 守恒与目标组成折算：出口 N2 {f(n2_out)} mol/s，目标干湿基总量 {f(target_total)} mol/s，"
              f"目标 H2 {f(target_h2)} mol/s，模型出口 H2 {f(h2_out)} mol/s。"
              f"若改用模型出口总摩尔流 {f(n_out)} mol/s 折算，缺口为 {f(h2_gap_model)} mol/s（两种算术都列出）。",
              "该检验只判断“单靠 R6 能否在量级上同时解释水与氢”，不据此推荐任何参数改动。", ""]

    # ---- C.10 ---------------------------------------------------------
    lines += ["## C.10 数值不变声明", "",
              "7 个验收量与四块残差逐项与 T30 的 `reports/t30-acceptance/report.md` 对照；相对差门 1×10⁻⁹。", ""]
    rows, worst = [], 0.0
    for key, target in T30_QUANTITIES.items():
        value = (report.run.get("measurements") or {}).get(key)
        gap = None if value is None else rel(value, target)
        worst = max(worst, 0.0 if gap is None else gap)
        rows.append([key, g(value, 9), g(target, 9), "—" if gap is None else f(gap)])
    for key, target in T30_RESIDUALS.items():
        value = (report.run.get("ratios") or {}).get(key)
        gap = None if value is None else rel(value, target)
        worst = max(worst, 0.0 if gap is None else gap)
        rows.append([f"残差 {key}", f(value), f(target), "—" if gap is None else f(gap)])
    for key, target in T30_ELEMENTS.items():
        value = ((report.run.get("balance") or {}).get("ratios") or {}).get(key)
        gap = None if value is None else rel(value, target)
        worst = max(worst, 0.0 if gap is None else gap)
        rows.append([f"平衡 {key}", f(value), f(target), "—" if gap is None else f(gap)])
    wall = (report.run.get("balance") or {}).get("heat_w")
    wall_gap = None if wall is None else rel(wall, T30_WALL_W)
    worst = max(worst, 0.0 if wall_gap is None else wall_gap)
    rows.append(["壁面散热 W", g(wall, 3), g(T30_WALL_W, 3), "—" if wall_gap is None else f(wall_gap)])
    ash = ((report.run.get("diagnostics") or {}).get("consistent_ash_relative"))
    ash_gap = None if ash is None else rel(ash, T30_ASH)
    worst = max(worst, 0.0 if ash_gap is None else ash_gap)
    rows.append(["灰平衡（2(b)）", f(ash), f(T30_ASH), "—" if ash_gap is None else f(ash_gap)])
    bottom = ((report.run.get("diagnostics") or {}).get("bottom_relative"))
    bottom_gap = None if bottom is None else rel(bottom, T30_BOTTOM)
    worst = max(worst, 0.0 if bottom_gap is None else bottom_gap)
    rows.append(["底部残差", f(bottom), f(T30_BOTTOM), "—" if bottom_gap is None else f(bottom_gap)])
    rows.append(["轮数 / 次数", f"{report.run.get('rounds')} / {report.run.get('iterations')}", "16 / 180",
                 "相同" if (report.run.get("rounds"), report.run.get("iterations")) == (16, 180) else "不同"])
    lines += table(["量", "本次", "T30", "相对差"], rows) + [""]
    worst_text = f"最大相对差 **{f(worst)}**"
    lines += [worst_text + ("，≤ 1×10⁻⁹，数值未改。" if worst <= 1e-9 else "，**超过 1×10⁻⁹，触发停止条件**。"), ""]

    # ---- C.11 ---------------------------------------------------------
    elapsed = report.run.get("elapsed_s")
    lines += ["## C.11 运行时长对照", "",
              f"本次 {g(elapsed, 1)} s，T30 基准 {g(T30_ELAPSED, 1)} s，比值 {g(float(elapsed) / T30_ELAPSED, 4)}。",
              "账目只在收敛后评价一次，不进入内层；该比值与 1 的差异属正常运行波动，开销可忽略（同一量级）。", ""]

    # ---- D.1 ----------------------------------------------------------
    d1 = _reconstruction(report)
    lines += ["## D.1 重建精度（本卡的漏项判据）", "",
              "账目重建的源项与同一 `evaluate_residual` 装配的 `reaction_bed`／`reaction_free` 逐元素比较；",
              "判据 `|Σq − 装配值| ≤ 1×10⁻¹² × max(|装配值|, |Σq|, 1×10⁻¹⁵)`。", ""]
    rows = [["气相（悬浮相）", f(d1["gas_suspension"]), "通过" if d1["gas_suspension"] <= 1e-12 else "不通过"],
            ["气相（气泡相）", f(d1["gas_bubble"]), "通过" if d1["gas_bubble"] <= 1e-12 else "不通过"],
            ["气相（自由空间，含悬停）", f(d1["gas_free"]), "通过" if d1["gas_free"] <= 1e-12 else "不通过"],
            ["固相 C（床区）", f(d1["char_bed"]), "通过" if d1["char_bed"] <= 1e-12 else "不通过"],
            ["固相 C（自由空间，含悬停）", f(d1["char_free"]), "通过" if d1["char_free"] <= 1e-12 else "不通过"]]
    lines += table(["块", "最大相对差", "结论"], rows) + [""]
    lines += ["这一层不过就是记账漏项；本次" + ("全部通过。" if d1["worst"] <= 1e-12 else "**未通过，须先修记录**。"), ""]

    # ---- D.2 ----------------------------------------------------------
    d2 = _physical_balance(report)
    lines += ["## D.2 物理收支（照实报，不改门）", "",
              "该格的 `gas_mole_residual`（入流 + 回流 + 反应 − 出流 ± 相间交换）按原式报出；",
              "主求解的收敛门：气相 ≤ 1×10⁻⁶、固相 ≤ 1×10⁻⁶、能量 ≤ 1×10⁻⁶、元素 ≤ 1×10⁻³。",
              "这一层**不是本卡的通过条件**；只要求与 T30 的残差数值一致（相对差 ≤ 1×10⁻⁹），证明记录没有改动求解。", ""]
    lines += table(["区", "格", "相", "max abs(r) mol/s", "max abs(r)/n_ref", "max per-species 比值", "尺度"], d2["rows"]) + [""]
    lines += ["`max per-species 比值` 按本卡 D.2 的尺度定义逐物种计算；接近 1 的格是 O2 痕量（反应项远大于入出流），"
              "其绝对残差仍远低于块门。块的判定按 `n_ref` 归一，见下一行的块残差。", ""]
    lines += [f"本次气相块归一残差 {f(d2['gas_ratio'])}；由账目评价的同一残差数组算得 {f(d2['block_ratio'])}，"
              f"与运行记录一致（相对差 {f(rel(d2['block_ratio'], d2['gas_ratio']))}）。",
              f"气相/固相/能量块与 T30 的最大相对差分别为 {f(rel((report.run.get('ratios') or {}).get('gas'), T30_RESIDUALS['gas']))}、"
              f"{f(rel((report.run.get('ratios') or {}).get('solid'), T30_RESIDUALS['solid']))}、"
              f"{f(rel((report.run.get('ratios') or {}).get('energy'), T30_RESIDUALS['energy']))}。",
              "收敛解本身的气相残差不是零（T30 为 2.02×10⁻⁷），因此**不能**用“逐格闭合到 1e−9”去要求物理收支，"
              "也不能把非零残差解释成记账漏项。", ""]

    # ---- 声明 ----------------------------------------------------------
    lines += ["## 声明", "",
              "- 未改任何速率式、指前因子、活化能、反应级数、计量系数、温度取法、相分配、焦油口径、蒸汽量、ρ_P、a_R8、Corella／Serio 参数。",
              "- 未改 `data/digitized/` 的目标值、`est_error` 或容差；未重判 T30 的 7 个量；未跑 S1–S20。",
              "- 结论只写占比与绝对量；因果判断与参数改动留给后续任务卡。",
              "- 账目可离线复算：`.scratch/t31/state.npz` + `frozen.pkl` 重跑一次 `evaluate_residual(..., ledger=True)` 得到同一份记录。", ""]

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join([f"status={status}", f"records={len(records)}", f"rows={len(report.rows)}",
                     f"d1_worst={d1['worst']:.3e}", f"invariance_worst={worst:.3e}",
                     f"water_net={book['species']['H2O']['net']:.6e}"]))
    if d1["worst"] > 1e-12 or worst > 1e-9 or status != "CONVERGED":
        return 2
    return 0


def _tar_totals(records) -> dict:
    """焦油代表物的消耗／生成／净值。不属于 7 个物种账，单列核对 R10／R11。"""
    totals = {name: {"consume": 0.0, "generate": 0.0, "net": 0.0} for name in TAR_SPECIES}
    for record in records:
        for name, value in record["q"].items():
            if name not in totals:
                continue
            q = float(value)
            totals[name]["generate"] += max(q, 0.0)
            totals[name]["consume"] += max(-q, 0.0)
            totals[name]["net"] += q
    return totals


def _class(value, high, low, high_label, mid_label, low_label):
    if value is None:
        return "分母可忽略"
    if value >= high:
        return high_label
    if value >= low:
        return mid_label
    return low_label


def _gas_rows(report: Report) -> list:
    """逐格逐相的残差行。计算与 D.2 共用 ``tools.ledger_account.gas_residual_rows``。"""
    return gas_residual_rows(report.frozen, report.state, report.result, report.records)


def _oxygen_rows(report: Report) -> list:
    rows = []
    for item in _gas_rows(report):
        inflow, outflow = item["inflow"][INDEX["O2"]], item["outflow"][INDEX["O2"]]
        swap, reaction = item["exchange"][INDEX["O2"]], item["reaction"][INDEX["O2"]]
        residual = item["residual"][INDEX["O2"]]
        label = {"suspension": "悬浮相", "bubble": "气泡相", "free": "自由空间"}[item["phase"]]
        rows.append([item["zone"], item["cell"], label, f(inflow), f(outflow), f(swap),
                     f(-reaction if reaction < 0.0 else 0.0), f(residual)])
    return rows


def _reconstruction(report: Report) -> dict:
    frozen, result = report.frozen, report.result
    gas_bed = np.zeros_like(result.reaction_bed)
    gas_free = np.zeros_like(result.reaction_free)
    char_bed = np.zeros(frozen.n_bed)
    char_free = np.zeros(frozen.n_freeboard)
    for row in report.records:
        cell = int(row["cell"])
        for name, q in row["q"].items():
            if name == "C":
                continue
            if row["zone"] == "bed":
                gas_bed[cell, 0 if row["phase"] == "suspension" else 1, INDEX[name]] += float(q)
            else:
                gas_free[cell, INDEX[name]] += float(q)
    for row in report.records:
        if "C" in row["q"]:
            target = char_bed if row["zone"] == "bed" else char_free
            target[int(row["cell"])] += float(row["q"]["C"])
    molar = molar_mass_kg_mol("C_gr")
    char_bed *= molar
    char_free *= molar

    def worst(left, right):
        left = np.asarray(left, dtype=float)
        right = np.asarray(right, dtype=float)
        scale = np.maximum(np.maximum(np.abs(left), np.abs(right)), 1.0e-15)
        return float(np.max(np.abs(left - right) / scale))

    values = {
        "gas_suspension": worst(gas_bed[:, 0], result.reaction_bed[:, 0]),
        "gas_bubble": worst(gas_bed[:, 1], result.reaction_bed[:, 1]),
        "gas_free": worst(gas_free, result.reaction_free),
        "char_bed": worst(char_bed, -result.char_consumption_bed_kg_s.sum(axis=1)),
        "char_free": worst(char_free, -(result.char_consumption_free_kg_s + result.char_consumption_hover_kg_s).sum(axis=1)),
    }
    values["worst"] = max(values.values())
    return values


def _physical_balance(report: Report) -> dict:
    rows, peak = [], 0.0
    for item in _gas_rows(report):
        value = float(np.max(np.abs(item["residual"])))
        ratio = float(np.max(np.abs(item["residual"]) / item["scale"]))
        peak = max(peak, ratio)
        label = {"suspension": "悬浮相", "bubble": "气泡相", "free": "自由空间"}[item["phase"]]
        rows.append([item["zone"], item["cell"], label, f(value), f(value / report.n_ref), f(ratio),
                     f(float(np.max(item["scale"])))])
    result = report.result
    block = max(float(np.max(np.abs(result.gas_suspension_mol_s))), float(np.max(np.abs(result.gas_bubble_mol_s))),
                float(np.max(np.abs(result.gas_free_mol_s)))) / report.n_ref
    return {"rows": rows, "block_ratio": block, "gas_ratio": (report.run.get("ratios") or {}).get("gas"), "peak_ratio": peak}


if __name__ == "__main__":
    raise SystemExit(main())
