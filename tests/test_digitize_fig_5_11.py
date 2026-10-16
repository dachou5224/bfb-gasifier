# -*- coding: utf-8 -*-
"""Bild 5.11 读数的复算与失败路径。改副本，不改 data/digitized 里的正式文件；不读 PDF。"""
import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("digitize_511", ROOT / "tools" / "digitize_fig_5_11.py")
DIGIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIGIT)


def load():
    return json.loads((ROOT / "data" / "digitized" / "fig_5_11.json").read_text(encoding="utf-8"))


class Fig511Test(unittest.TestCase):
    def test_recorded_points_match_calibration(self):
        self.assertEqual(DIGIT.check(load()), [])

    def test_curve_decreases_with_temperature(self):
        values = [k for _, k in load()["points"]]
        self.assertTrue(all(a > b for a, b in zip(values, values[1:])))

    def test_crosses_one_near_1100_k(self):
        # 曲线在 1100 K 处穿过 1,0 网格线（图上可目视核对）
        points = dict((t, k) for t, k in load()["points"])
        self.assertAlmostEqual(points[1100], 1.0, delta=0.01)

    def test_edited_point_is_reported(self):
        doc = copy.deepcopy(load())
        doc["points"][3][1] = 1.02
        self.assertTrue(any("points" in p for p in DIGIT.check(doc)))

    def test_column_off_its_temperature_is_reported(self):
        doc = copy.deepcopy(load())
        doc["calibration"]["centerline"][0]["pixel_x"] += 20
        self.assertTrue(any("T=1025" in p for p in DIGIT.check(doc)))

    def test_error_must_be_numeric(self):
        doc = copy.deepcopy(load())
        doc["est_error"]["K_eq"] = "约 0,01"
        self.assertTrue(any("est_error.K_eq" in p for p in DIGIT.check(doc)))


if __name__ == "__main__":
    unittest.main()
