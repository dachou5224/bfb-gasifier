"""T26r 验证：工况 1、φ = 1，外层 12 轮。Corella 只在本工具外推。不记录出口组成、T̄_FRB、X_C。"""
from __future__ import annotations

import argparse, json, signal, time
from dataclasses import replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import hamel.kinetics.homogeneous as homogeneous
import hamel.solver.driver as driver
from hamel.cases.assemble import assemble_case
from hamel.cases.config import SolverConfig
from hamel.solver.gates import _abs_max
from hamel.solver.newton import SolveTimeout

_ROOT = Path(__file__).resolve().parents[2]
_CASE, _OUT = _ROOT / "data/cases/vtt_pressurized.json", _ROOT / "reports/t26r-validation"


class _Timeout(BaseException):
    """墙钟到点。CHARTER §7。"""


def run(rounds, limit_s, out):
    original, started = homogeneous._CORELLA_HIGH_K, time.monotonic()
    homogeneous._CORELLA_HIGH_K = float("inf")
    n_ref, m_ref, _q = driver._references(assemble_case(str(_CASE), "1").physical)
    box, pre, mon = {"pre": [], "mon": [], "rows": []}, driver.precalculate, driver._monitor

    def precalculate(*args, **kwargs):
        box["pre"].append(pre(*args, **kwargs))
        return box["pre"][-1]

    def monitor(state, frozen, assembly=None):
        value = mon(state, frozen, assembly)
        box["mon"].append(value)
        inert = frozen.geometric_mass_kg - state.bed_char_kg.sum(axis=1) - frozen.fuel_bed_kg.sum(axis=1)
        ash = (box["pre"][-1].reported or {}).get("ash_balance") or {}
        box["rows"].append(dict(
            phi=float(frozen.fly_ash_fraction), w_f=float(frozen.bed_ash_fine_share),
            E=None if ash.get("E") is None else float(ash["E"]),
            target=float(ash.get("target") or frozen.ash_target_kg_s or 0.0),
            ash_relative=value.ash_relative, bed_char_kg=float(np.sum(state.bed_char_kg)),
            min_inert_kg=float(np.min(inert)), inert_kg=inert.tolist(),
            bed_temperature_k=np.asarray(state.bed_temperature_k).tolist(),
        ))
        (Path(out) / "progress.json").write_text(
            json.dumps(box["rows"], ensure_ascii=False, indent=1, default=str), encoding="utf-8",
        )
        print(f"monitor {len(box['rows'])} phi={frozen.fly_ash_fraction} char={np.sum(state.bed_char_kg):.4f}", flush=True)
        return value

    driver.precalculate, driver._monitor = precalculate, monitor
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(_Timeout()))
    signal.setitimer(signal.ITIMER_REAL, float(limit_s))
    summary = {}
    try:
        result = driver.solve_case(
            str(_CASE), "1", config=replace(SolverConfig(), outer_round_limit=int(rounds)),
            deadline=lambda: time.monotonic() - started >= limit_s, artifact_dir=str(out),
        )
        summary.update(status=result.status, reason=result.reason, rounds=result.rounds, iterations=result.iterations)
        initial = result.diagnostics.get("initial") or {}
        summary["initial"] = {key: initial.get(key) for key in ("phi_init", "phi", "note")}
        summary["initial"]["loops"] = len(initial.get("loops") or [])
        m_tot = float(box["pre"][-1].frozen.geometric_mass_kg.sum()) if box["pre"] else 1.0
        gate_from = 3 if initial.get("phi_init") != initial.get("phi") else 2
        for index, item in enumerate(result.diagnostics.get("rounds") or []):
            newton = item.get("newton") or []
            row = box["rows"][index] if index < len(box["rows"]) else {}
            row.update(
                round=item.get("round"), inner=item.get("inner"), inner_reason=item.get("reason"),
                inner_iterations=len(newton), switched=bool(item.get("switched")),
                projected=int(sum(entry.get("projected", 0) for entry in newton)),
                overshoot=float(max((entry.get("overshoot", 0.0) for entry in newton), default=0.0)),
            )
            if 0 < index < len(box["mon"]):
                cur, prev = box["mon"][index], box["mon"][index - 1]
                row["gate"] = dict(
                    dT=_abs_max(cur.temperature_k, prev.temperature_k),
                    solid=_abs_max(cur.solid_flow_kg_s, prev.solid_flow_kg_s) / m_ref,
                    inventory=_abs_max(cur.inventory_kg, prev.inventory_kg) / m_tot,
                    gas=_abs_max(cur.gas_flow_mol_s, prev.gas_flow_mol_s) / n_ref,
                    checked=item.get("round", 0) >= gate_from,
                )
            if index < len(box["rows"]):
                box["rows"][index] = row
            else:
                box["rows"].append(row)
        path = out / "failure.npz"
        if "雅可比" in str(result.reason) and path.exists():
            data = np.load(path, allow_pickle=True)
            matrix, y, x = data["J"], data["y"], data["x"]
            names = list(map(str, data["names"])) if "names" in data.files else []

            def named(axis):
                idx = [i for i in range(matrix.shape[axis]) if not np.any(matrix[i] if axis == 0 else matrix[:, i])]
                return [{"index": i, "name": names[i] if i < len(names) else str(i),
                         "y": float(y[i]), "x": float(x[i])} for i in idx]

            summary["singular"] = {"zero_columns": named(1), "zero_rows": named(0)}
    except (_Timeout, SolveTimeout) as exc:
        summary.update(status="TIMEOUT", reason=str(exc) or "运行超过墙钟上限", rounds=len(box["rows"]))
    except Exception as exc:  # noqa: BLE001
        summary.update(status=type(exc).__name__, reason=str(exc)[:800], rounds=len(box["rows"]))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        driver.precalculate, driver._monitor = pre, mon
        homogeneous._CORELLA_HIGH_K = original
    summary["elapsed_s"] = time.monotonic() - started
    return {"summary": summary, "rounds": box["rows"]}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit-s", type=float, default=3.0 * 3600.0)
    parser.add_argument("--rounds", type=int, default=12)
    parser.add_argument("--out", default=str(_OUT))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    payload = run(args.rounds, args.limit_s, out)
    (out / "extrapolated.json").write_text(json.dumps(payload, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print(json.dumps(payload["summary"], ensure_ascii=False, default=str), flush=True)


if __name__ == "__main__":
    main()
