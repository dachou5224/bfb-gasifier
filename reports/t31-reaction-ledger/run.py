"""T31 的单次基准运行：工况 1、adiabatic=False，配置与 T30 逐字相同（26 轮 / 28 800 s）。

走 ``tools.campaign.execute`` 的单侧入口（与 T30 相同），进程级硬墙钟由它设置；
``artifact_dir`` 使收敛后落盘 ``state.npz`` 与 ``frozen.pkl``，供离线复算账目。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.campaign import execute  # noqa: E402
from tools.driver_report import write_report  # noqa: E402

OUT = ROOT / ".scratch" / "t31"
CASE = str(ROOT / "data" / "cases" / "vtt_pressurized.json")
ROUNDS, LIMIT = 26, 28_800.0


def _trim_report_json(folder) -> None:
    """入库的 driver_report.json 只留结论与账目摘要；完整记录在 .scratch/t31/run.json。"""
    keep = (
        "negative_floor", "final_check", "bottom_relative", "consistent", "consistent_ash_relative",
        "bed_temperature_peak", "bed_peak_per_round", "trial_bed_peak_k", "geldart", "adiabatic",
        "bed_gas_temperature_deviation_k", "hover_drift", "trajectory", "reaction_ledger",
    )
    path = Path(folder) / "driver_report.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    for run in payload.get("runs", []):
        diagnostics = run.get("diagnostics") or {}
        runner = {key: diagnostics[key] for key in keep if key in diagnostics}
        ledger = runner.get("reaction_ledger")
        if ledger:
            ledger.pop("records", None)
            ledger["records_note"] = "完整逐记录表见 .scratch/t31/run.json（不入库）"
        run["diagnostics"] = runner
        run["diagnostics"]["omitted"] = "逐轮数组与完整反应记录见 .scratch/t31/run.json；报告表见 driver_report.md"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")


def main() -> int:
    folder = ROOT / "reports" / "t31-reaction-ledger"
    if "--render" in sys.argv[1:]:
        # 不重跑求解：用 .scratch/t31/run.json 重新渲染报告（改了报告工具时用）。
        record = json.loads((OUT / "run.json").read_text(encoding="utf-8"))
        write_report(folder, [record], title="第 1 关反应路径账运行", toy=False)
        _trim_report_json(folder)
        print(json.dumps({"status": record.get("status"), "rendered": True}, ensure_ascii=False))
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    record = execute("baseline", CASE, False, outer_round_limit=ROUNDS, wall_timeout_s=LIMIT, artifact_dir=str(OUT))
    (OUT / "run.json").write_text(json.dumps(record, ensure_ascii=False), encoding="utf-8")
    write_report(folder, [record], title="第 1 关反应路径账运行", toy=False)
    _trim_report_json(folder)
    print(json.dumps({key: record.get(key) for key in ("status", "rounds", "iterations", "elapsed_s")}, ensure_ascii=False))
    return 0 if record.get("status") == "CONVERGED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
