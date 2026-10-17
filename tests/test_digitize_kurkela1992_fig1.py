# -*- coding: utf-8 -*-
"""KURKELA_1992 Fig. 1 读数的复算与失败路径。改副本，不改 data/digitized 里的正式文件；不读 PDF。"""
import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("digitize_kf1", ROOT / "tools" / "digitize_kurkela1992_fig1.py")
DIGIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIGIT)


def load():
    return json.loads((ROOT / "data" / "digitized" / "kurkela1992_fig1.json").read_text(encoding="utf-8"))


class KurkelaFig1Test(unittest.TestCase):
    def test_recorded_values_match_calibration(self):
        self.assertEqual(DIGIT.check(load()), [])

    def test_class_fractions_sum_to_one_hundred(self):
        self.assertAlmostEqual(sum(load()["class_mass_fractions_pct"]), 100.0, delta=0.2)

    def test_points_sit_on_standard_sieves(self):
        # 圆心横坐标与标准筛孔一致（0,1 mm 的圆压在左边框上，允许 6 %）
        doc = load()
        size_mm, _ = DIGIT.axes(doc["calibration"])
        for center, sieve in zip(doc["calibration"]["circle_centers_px"], DIGIT.SIEVES_MM):
            self.assertLess(abs(size_mm(center[0]) / sieve - 1), 0.06)

    def test_edited_point_is_reported(self):
        doc = copy.deepcopy(load())
        doc["points"][0][1] = 7.0
        self.assertTrue(any("points" in p for p in DIGIT.check(doc)))

    def test_center_off_its_sieve_is_reported(self):
        # 把 1,0 mm 的圆心挪到 1,2 mm 附近：应报筛孔不符
        doc = copy.deepcopy(load())
        doc["calibration"]["circle_centers_px"][4][0] += 76
        self.assertTrue(any("1.0 mm" in p for p in DIGIT.check(doc)))

    def test_non_monotonic_curve_is_reported(self):
        doc = copy.deepcopy(load())
        doc["calibration"]["circle_centers_px"][3][1] = 900.0
        doc["points"], doc["class_mass_fractions_pct"] = DIGIT.recompute(doc)
        self.assertTrue(any("不单调" in p for p in DIGIT.check(doc)))

    def test_degenerate_or_reversed_frame_is_reported_not_raised(self):
        for key, pair in (("x_frame_px", [500.0, 500.0]), ("y_frame_px", [700.0, 700.0]),
                          ("x_frame_px", [1888.0, 263.5]), ("y_frame_px", [72.0, 1313.5]),
                          ("x_frame_px", [263.5, float("nan")])):
            doc = copy.deepcopy(load())
            doc["calibration"][key] = pair
            problems = DIGIT.check(doc)  # 不得抛异常
            self.assertTrue(any(key in p or "坐标框" in p for p in problems), (key, pair, problems))

    def test_invalid_error_values_are_reported(self):
        for bad in (True, -1, float("nan"), float("inf"), "约 0,5", None):
            doc = copy.deepcopy(load())
            doc["est_error"]["cumulative_pct"] = bad
            self.assertTrue(any("est_error.cumulative_pct" in p for p in DIGIT.check(doc)), bad)

    def test_error_must_be_numeric(self):
        doc = copy.deepcopy(load())
        doc["est_error"]["cumulative_pct"] = "约 0,5"
        self.assertTrue(any("est_error.cumulative_pct" in p for p in DIGIT.check(doc)))


if __name__ == "__main__":
    unittest.main()
