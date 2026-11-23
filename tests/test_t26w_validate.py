"""T26w 验证入口：硬墙钟记 TIMEOUT（CHARTER §7）与跨环境轮数判据。不跑长工况。"""
import time

import numpy as np

from hamel.solver.newton import SolveTimeout


def test_blocking_solve_is_recorded_as_timeout(tmp_path, monkeypatch):
    from tools import t26w_validate as mod

    monkeypatch.setattr(mod, "_OUT", tmp_path / "baseline.json")
    monkeypatch.setattr(mod, "_EVIDENCE", tmp_path / "missing.json")

    def blocking(*args, **kwargs):
        time.sleep(2.0)
        raise AssertionError("硬墙钟没有中断求解")

    monkeypatch.setattr(mod.driver, "solve_case", blocking)
    real = mod.driver._monitor
    assert mod.main(["--limit-s", "0.2"]) == 2
    record = __import__("json").loads((tmp_path / "baseline.json").read_text(encoding="utf-8"))
    assert record["status"] == "TIMEOUT" and "墙钟上限" in record["reason"]
    assert mod.driver._monitor is real  # 监测补丁已还原


def test_round_count_may_differ_in_cross_environment_comparison():
    """STATUS 2026-11-23 第 5 条：轮数不同不否决；只看已配对轮次与终态。"""
    from tools import t26w_validate as mod

    evidence = {"rounds_done": 17, "iterations": 181, "trial_bottom_max_c": 972.97,
                "rounds": [{"index": i + 1, "bed_peak_c": p, "bed_temperature_k": [1001.0]}
                           for i, p in enumerate([957.44, 901.63, 908.94])]}
    record = {"status": "CONVERGED", "final_bed_k": [1001.0],
              "bed_peak_per_round": [{"round": 1, "peak_k": 957.44 + 273.15},
                                     {"round": 2, "peak_k": 901.63 + 273.15}],
              "trial_bottom_peak_k": 972.97 + 273.15}
    # 取证侧最后一轮 bed_temperature_k 与 record 的终态对齐（只比容差内的相对差）。
    comparison = mod.compare(record, evidence)
    assert comparison["rounds_covered"] == 2
    assert comparison["rounds_differ"] is True
    assert comparison["passed"] is True


def test_failure_path_record_carries_the_path_and_marks_itself_unconverged():
    """ADR 0059 第 2a 项：未收敛／越线停止也必须有记录，且不冒充终解。"""
    from hamel.solver.driver import _path_peak
    from hamel.solver.vector import PathPeaks

    track = PathPeaks()
    track.bottom_k, track.peak_k = 1246.129, 1250.0
    peaks = [{"round": 1, "peak_k": 1230.59, "peak_c": 957.44, "from_top": 8, "balance_index": 7,
              "above_900c_k": 57.44, "above_900c_cells": 1, "above_950c_k": 7.44}]
    record = _path_peak(peaks, track)
    assert record["converged"] is False and record["peak_c"] == 957.44
    assert record["above_950c_k"] == 7.44
    assert record["trial_bed_peak_k"] == 1250.0 and record["trial_bottom_peak_k"] == 1246.129
    assert record["per_round"] == peaks and "round" not in record
    empty = _path_peak([], PathPeaks())
    assert empty["converged"] is False and empty["peak_c"] is None and empty["above_950c_k"] == 0.0
