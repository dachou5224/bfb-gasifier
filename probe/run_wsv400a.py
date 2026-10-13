"""用论文模式默认 Reactor 跑 WSV400 A。单次硬超时 420 s。

超时或异常时仍写出当时的残差快照。不修改 src。
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import traceback
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

from probe.wsv400a_config import SETTINGS, build_reactor_config, load_parity_wsv400_a  # noqa: E402
from src.core.cell import S_ASH, S_CHAR, S_MOISTURE, S_VM  # noqa: E402
from src.core.elemental_ledger import gas_element_molar_rates, solid_element_molar_rates  # noqa: E402
from src.core.reactor import Reactor  # noqa: E402
from src.core.species import GAS_SPECIES  # noqa: E402
from src.solvers.global_nr_solver import build_equation_scales, global_residual  # noqa: E402
from src.solvers.nr_indexing import (  # noqa: E402
    cell_offsets,
    gas_mode_for_nr,
    has_temperature_var,
    main_nr_gas_idx_for_cell,
    n_gas_var,
    n_var,
    pack_reactor,
    solid_flat_indices_for_nr,
)

HARD_TIMEOUT_S = 420
_SOLID_NAME = {S_CHAR: "char", S_VM: "vm", S_MOISTURE: "moisture", S_ASH: "ash"}
_M_KG_PER_MOL = {
    "C": 12.011e-3,
    "H": 1.00794e-3,
    "O": 15.999e-3,
    "N": 14.0067e-3,
    "S": 32.065e-3,
}
_H2O_KG_PER_MOL = 18.01528e-3

_STATE: dict[str, Any] = {}


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_jsonable(v) for v in value.tolist()]
    if isinstance(value, np.floating):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(_jsonable(payload), indent=2, ensure_ascii=False) + "\n"
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def _inlet_atoms(cfg: Any) -> dict[str, float]:
    wet = float(cfg.fuel_feed)
    moisture = float(cfg.moisture_wt) / 100.0
    dry = wet * (1.0 - moisture)
    atoms = {
        "C": dry * float(cfg.C_dry) / 100.0 / _M_KG_PER_MOL["C"],
        "H": dry * float(cfg.H_dry) / 100.0 / _M_KG_PER_MOL["H"],
        "O": dry * float(cfg.O_dry) / 100.0 / _M_KG_PER_MOL["O"],
        "N": dry * float(cfg.nitrogen_fraction) / 100.0 / _M_KG_PER_MOL["N"],
        "S": dry * float(cfg.S_dry) / 100.0 / _M_KG_PER_MOL["S"],
    }
    n_h2o_fuel = wet * moisture / _H2O_KG_PER_MOL
    atoms["H"] += 2.0 * n_h2o_fuel
    atoms["O"] += n_h2o_fuel
    atoms["O"] += 2.0 * float(cfg.O2_feed)
    atoms["H"] += 2.0 * float(cfg.H2O_feed)
    atoms["O"] += float(cfg.H2O_feed)
    atoms["N"] += 2.0 * float(cfg.N2_feed)
    return atoms


def _gas_atoms(cell: Any, fuel_type: str) -> dict[str, float]:
    flow = np.maximum(np.asarray(cell.N_d, dtype=np.float64), 0.0) + np.maximum(
        np.asarray(cell.N_b, dtype=np.float64), 0.0
    )
    return gas_element_molar_rates(flow, fuel_type=fuel_type)


def _solid_atoms(rates: np.ndarray, cell: Any) -> dict[str, float]:
    solid = cell.solid
    return solid_element_molar_rates(
        rates,
        ash_dry_wt=float(solid.ash_dry_wt),
        C_dry=float(solid.C_dry),
        H_dry=float(solid.H_dry),
        O_dry=float(solid.O_dry),
        nitrogen_fraction=float(solid.nitrogen_fraction),
        sulfur_fraction=float(solid.sulfur_fraction),
        sulfur_volatile_frac=float(solid.sulfur_volatile_frac),
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )


def _balance(reactor: Reactor) -> dict[str, Any]:
    """外部账：进口原子减去床顶气相流与床顶向上固体流。

    N_d、N_b 按求解器里的 mol/s 使用。循环固体不从床顶向上流里扣除；
    本工况 recirculation_frac=0。硫的固体侧用 SolidProps.sulfur_fraction，
    构造函数没有把 S_dry 写进去。
    """
    cfg = reactor.config
    inlet = _inlet_atoms(cfg)
    top = reactor.cells[-1]
    gas = _gas_atoms(top, str(cfg.fuel_type))
    solid_up = _solid_atoms(np.asarray(top._solid_upflow_rates(), dtype=np.float64), top)
    outlet = {el: float(gas[el]) + float(solid_up[el]) for el in inlet}
    extra: dict[str, Any] = {}
    for name in ("cyclone_cell", "return_leg_cell"):
        side = getattr(reactor, name, None)
        if side is None:
            continue
        extra[name] = _gas_atoms(side, str(cfg.fuel_type))
    gap = {el: float(inlet[el]) - float(outlet[el]) for el in inlet}
    energy_rows = []
    for idx, cell in enumerate(reactor._solver_cells_for_nr()):
        try:
            residual_w = float(cell.calc_energy_balance())
        except Exception as exc:  # noqa: BLE001 — 诊断，不改方程
            residual_w = None
            energy_rows.append(
                {
                    "index": idx,
                    "cell_type": str(getattr(cell, "cell_type", "")),
                    "residual_W": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        energy_rows.append(
            {
                "index": idx,
                "cell_type": str(getattr(cell, "cell_type", "")),
                "T_K": float(cell.T),
                "residual_W": residual_w,
            }
        )
    finite = [row["residual_W"] for row in energy_rows if isinstance(row.get("residual_W"), float)]
    return {
        "method": "inlet_atoms - (bed_top N_d+N_b atoms + bed_top solid upflow atoms)",
        "inlet_mol_atom_s": inlet,
        "bed_top_gas_mol_atom_s": gas,
        "bed_top_solid_upflow_mol_atom_s": solid_up,
        "side_gas_mol_atom_s": extra,
        "gap_mol_atom_s": gap,
        "energy_residual_sum_W": float(sum(finite)) if finite else None,
        "energy_residual_max_abs_W": float(max(abs(v) for v in finite)) if finite else None,
        "energy_residuals": energy_rows,
        "solid_sulfur_fraction_on_cell": float(top.solid.sulfur_fraction),
    }


def _residual_labels(cells: list[Any]) -> list[dict[str, Any]]:
    labels: list[dict[str, Any]] = []
    for idx, cell in enumerate(cells):
        kind = str(getattr(cell, "cell_type", "bed"))
        mode = gas_mode_for_nr(cell)
        gas_idx = main_nr_gas_idx_for_cell(cell)
        names = [GAS_SPECIES[int(i)] for i in gas_idx]
        if mode == "two_phase":
            for phase in ("dense", "bubble"):
                for name in names:
                    labels.append({"cell_index": idx, "cell_type": kind, "equation": f"gas_{phase}_{name}"})
        elif mode == "single":
            for name in names:
                labels.append({"cell_index": idx, "cell_type": kind, "equation": f"gas_{name}"})
        flats = solid_flat_indices_for_nr(cell)
        n_comp = int(cell.m_solid.shape[1])
        for flat in flats:
            klass = int(flat) // n_comp
            comp = int(flat) % n_comp
            labels.append(
                {
                    "cell_index": idx,
                    "cell_type": kind,
                    "equation": f"solid_class{klass}_{_SOLID_NAME.get(comp, str(comp))}",
                }
            )
        if has_temperature_var(cell):
            labels.append({"cell_index": idx, "cell_type": kind, "equation": "energy"})
    return labels


def _ranked_residuals(reactor: Reactor) -> dict[str, Any]:
    cells = reactor._solver_cells_for_nr()
    cfg = reactor.config
    ref_gas = max(float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed), 1.0)
    ref_energy = max(abs(float(cfg.fuel_feed)) * 20.0e6, 1.0e6)
    ref_solid = max(abs(float(cfg.fuel_feed)), 1.0e-9)
    raw = global_residual(pack_reactor(cells), cells, reactor._apply_all_bc_for_nr)
    scale = build_equation_scales(cells, ref_gas, ref_solid, ref_energy)
    scaled = np.asarray(raw, dtype=np.float64) / np.asarray(scale, dtype=np.float64)
    labels = _residual_labels(cells)
    if len(labels) != int(scaled.size):
        raise RuntimeError(f"残差标签 {len(labels)} 与向量长度 {scaled.size} 不一致")
    order = np.argsort(-np.abs(scaled))
    top = []
    for rank, pos in enumerate(order[:10], start=1):
        item = dict(labels[int(pos)])
        item["rank"] = rank
        item["index"] = int(pos)
        item["residual_raw"] = float(raw[int(pos)])
        item["residual_scaled"] = float(scaled[int(pos)])
        top.append(item)
    offsets = cell_offsets(cells)
    blocks = {"gas": [], "solid": [], "energy": []}
    cursor = 0
    for cell in cells:
        ng = n_gas_var(cell)
        nv = n_var(cell)
        if ng:
            blocks["gas"].append(scaled[cursor : cursor + ng])
        sol = solid_flat_indices_for_nr(cell).size
        if sol:
            blocks["solid"].append(scaled[cursor + ng : cursor + ng + sol])
        if has_temperature_var(cell):
            blocks["energy"].append(scaled[cursor + nv - 1 : cursor + nv])
        cursor += nv
    del offsets

    def _rms(parts: list[np.ndarray]) -> float | None:
        if not parts:
            return None
        joined = np.concatenate(parts)
        if joined.size == 0:
            return None
        return float(np.sqrt(np.mean(np.square(joined))))

    return {
        "n_equations": int(scaled.size),
        "rms_scaled": {
            "gas": _rms(blocks["gas"]),
            "solid": _rms(blocks["solid"]),
            "energy": _rms(blocks["energy"]),
        },
        "max_abs_scaled": {
            "gas": float(np.max(np.abs(np.concatenate(blocks["gas"])))) if blocks["gas"] else None,
            "solid": float(np.max(np.abs(np.concatenate(blocks["solid"])))) if blocks["solid"] else None,
            "energy": float(np.max(np.abs(np.concatenate(blocks["energy"])))) if blocks["energy"] else None,
        },
        "top10": top,
        "scale_refs": {
            "ref_gas_mol_s": ref_gas,
            "ref_solid_kg_s": ref_solid,
            "ref_energy_W": ref_energy,
        },
    }


def _outer_trend(history: Any) -> list[dict[str, Any]]:
    if not isinstance(history, list):
        return []
    rows = []
    for item in history:
        if not isinstance(item, dict):
            continue
        rows.append(
            {
                "outer_iter": item.get("outer_iter"),
                "inner_n_iter": item.get("inner_n_iter"),
                "inner_rms_scaled_final": item.get("inner_rms_scaled_final"),
                "outer_energy_rms_final": item.get("outer_energy_rms_final"),
                "preinner_solid_rms_before": item.get("preinner_solid_rms_before"),
                "preinner_solid_rms_after": item.get("preinner_solid_rms_after"),
                "preinner_energy_rms_before": item.get("preinner_energy_rms_before"),
                "preinner_energy_rms_after": item.get("preinner_energy_rms_after"),
                "preinner_split_rms_before": item.get("preinner_split_rms_before"),
                "preinner_split_rms_after": item.get("preinner_split_rms_after"),
                "inner_converged": item.get("inner_converged"),
                "dT_outer": item.get("dT_outer"),
            }
        )
    return rows


def _exit_view(result: dict[str, Any] | None, reactor: Reactor) -> dict[str, Any]:
    parity = load_parity_wsv400_a()
    bed_T = [float(cell.T) for cell in reactor.cells]
    fb_T = [float(cell.T) for cell in getattr(reactor, "freeboard_cells", []) or []]
    view: dict[str, Any] = {
        "bed_T_K": bed_T,
        "bed_T_mean_K": float(np.mean(bed_T)) if bed_T else None,
        "freeboard_cell_T_K": fb_T,
        "freeboard_cell_T_mean_K": float(np.mean(fb_T)) if fb_T else None,
        "parity_fig_7_17_A": parity,
        "correspondence": (
            "采用 hamel_fig7_parity_digitized.json case_id_map："
            "CASE_WSV400_A → fig 7.17 标签 A。"
            "同一 JSON 里 CASE_WSV400_A.outputs 的轴向出口是 T=980 K、CO=0.152，"
            "与图 7.17 A 的 T_FRB=1135 K、CO sim=0.168 不一致。"
            "T3 尚未裁定。本探针的 Hamel 模拟对照用图 7.17 A，不把轴向出口当作 Sim。"
        ),
    }
    if isinstance(result, dict):
        for key in (
            "exit_gas",
            "exit_gas_dry",
            "bed_exit_gas_dry",
            "reactor_exit_T",
            "carbon_conv",
            "converged",
            "converged_outer",
            "converged_inner_nr",
            "converged_fully",
            "rms_scaled_final",
            "rms_scaled_gas_final",
            "rms_scaled_solid_final",
            "rms_scaled_energy_final",
            "max_abs_scaled_gas_final",
            "max_abs_scaled_solid_final",
            "max_abs_scaled_energy_final",
            "n_iter",
            "residual",
            "nr_outer_iters",
            "nr_inner_budget_total",
            "nr_inner_budget_used",
            "nr_merit_early_stop",
            "nr_merit_early_stop_reason",
            "nr_init_s_total",
            "nr_total_s",
            "freeboard_active",
        ):
            if key in result:
                view[key] = result[key]
        dry = result.get("exit_gas_dry") or {}
        if isinstance(dry, dict) and dry:
            sim = parity
            view["delta_vs_fig717_A_sim"] = {
                "T_K": (
                    None
                    if result.get("reactor_exit_T") is None
                    else float(result["reactor_exit_T"]) - float(sim["T_FRB_K"]["sim"])
                ),
                "CO": float(dry.get("CO", float("nan"))) - float(sim["CO"]["sim"]),
                "CO2": float(dry.get("CO2", float("nan"))) - float(sim["CO2"]["sim"]),
                "H2": float(dry.get("H2", float("nan"))) - float(sim["H2"]["sim"]),
                "CH4": float(dry.get("CH4", float("nan"))) - float(sim["CH4"]["sim"]),
            }
            view["delta_vs_fig717_A_meas"] = {
                "T_K": (
                    None
                    if result.get("reactor_exit_T") is None
                    else float(result["reactor_exit_T"]) - float(sim["T_FRB_K"]["meas"])
                ),
                "CO": float(dry.get("CO", float("nan"))) - float(sim["CO"]["meas"]),
                "CO2": float(dry.get("CO2", float("nan"))) - float(sim["CO2"]["meas"]),
                "H2": float(dry.get("H2", float("nan"))) - float(sim["H2"]["meas"]),
                "CH4": float(dry.get("CH4", float("nan"))) - float(sim["CH4"]["meas"]),
            }
    return view


def _slim_result(result: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(result, dict):
        return None
    kept = _exit_view(result, _STATE["reactor"])
    history = result.get("norm_history")
    if isinstance(history, list):
        kept["norm_history"] = [float(v) for v in history]
    kept["outer_trend"] = _outer_trend(result.get("nr_outer_history"))
    return kept


def _dump(status: str, error: str | None = None) -> None:
    reactor: Reactor | None = _STATE.get("reactor")
    result = _STATE.get("result")
    payload: dict[str, Any] = {
        "status": status,
        "error": error,
        "elapsed_s": perf_counter() - float(_STATE["t0"]),
        "timeout_s": HARD_TIMEOUT_S,
        "setting": _STATE["setting"],
        "meta": _STATE["meta"],
        "stage": _STATE.get("stage"),
    }
    if reactor is not None:
        try:
            payload["balance"] = _balance(reactor)
        except Exception:
            payload["balance_error"] = traceback.format_exc()
        try:
            payload["residuals"] = _ranked_residuals(reactor)
        except Exception:
            payload["residuals_error"] = traceback.format_exc()
        try:
            payload["exit"] = _exit_view(result if isinstance(result, dict) else None, reactor)
        except Exception:
            payload["exit_error"] = traceback.format_exc()
    if isinstance(result, dict):
        payload["solver"] = _slim_result(result)
    path = _STATE["output"]
    _write_json(path, payload)
    print(f"[probe] wrote {path} status={status} elapsed={payload['elapsed_s']:.1f}s", flush=True)


def _on_alarm(signum: int, frame: Any) -> None:
    del signum, frame
    _STATE["stage"] = "timeout"
    try:
        _dump("TIMEOUT", error=f"SIGALRM after {HARD_TIMEOUT_S} s")
    finally:
        os._exit(124)


def run_setting(setting: str, output: Path) -> int:
    cfg, meta = build_reactor_config(setting)
    _STATE["setting"] = setting
    _STATE["meta"] = meta
    _STATE["output"] = output
    _STATE["t0"] = perf_counter()
    _STATE["stage"] = "build"
    print(f"[probe] setting={setting} building Reactor", flush=True)
    reactor = Reactor(cfg)
    _STATE["reactor"] = reactor
    n_eq = sum(n_var(cell) for cell in reactor._solver_cells_for_nr())
    print(
        f"[probe] cells={len(reactor._solver_cells_for_nr())} equations={n_eq} "
        f"O2={cfg.O2_feed:.6e} H2O={cfg.H2O_feed:.6e} N2={cfg.N2_feed:.6e} mol/s",
        flush=True,
    )
    signal.signal(signal.SIGALRM, _on_alarm)
    signal.alarm(HARD_TIMEOUT_S)
    _STATE["stage"] = "solve"
    print("[probe] solve start", flush=True)
    try:
        result = reactor.solve(**meta["solve_kwargs"], verbose=False)
    except Exception:
        _STATE["stage"] = "exception"
        signal.alarm(0)
        _dump("EXCEPTION", error=traceback.format_exc())
        return 2
    signal.alarm(0)
    _STATE["result"] = result
    _STATE["stage"] = "solved"
    _dump("SOLVED", error=None)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="WSV400 A 论文模式探针")
    parser.add_argument("--setting", required=True, choices=sorted(SETTINGS))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or (_ROOT / "probe" / "out" / f"{args.setting}.json")
    return run_setting(args.setting, output)


if __name__ == "__main__":
    raise SystemExit(main())
