# -*- coding: utf-8 -*-
"""Bild 3.3 / 3.5 读数的复算与失败路径。改副本，不改正式 JSON；不读 PDF。"""
import copy
import importlib.util
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("digitize_335", ROOT / "tools" / "digitize_fig_3_5.py")
DIGIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIGIT)

PRESSURES = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]


def load(name):
    return json.loads((ROOT / "data" / "digitized" / name).read_text(encoding="utf-8"))


def points(doc, speed):
    for series in doc["series"]:
        if series["u0_minus_umf_m_s"] == speed:
            return dict(series["points"])
    raise KeyError(speed)


class Fig33And35Test(unittest.TestCase):
    def test_recorded_points_match_calibration(self):
        self.assertEqual(DIGIT.check(load("fig_3_3.json")), [])
        self.assertEqual(DIGIT.check(load("fig_3_5.json")), [])

    def test_each_series_has_the_six_adr_pressures(self):
        for name in ("fig_3_3.json", "fig_3_5.json"):
            for series in load(name)["series"]:
                self.assertEqual([p for p, _v in series["points"]], PRESSURES)

    def test_fig_3_3_order_follows_the_leaders(self):
        # 引线：0,1 m/s 在最上，0,6 m/s 在最下
        doc = load("fig_3_3.json")
        self.assertEqual(doc["series_order_from_top"], [0.1, 0.2, 0.3, 0.4, 0.6])
        at_1 = [points(doc, v)[1.0] for v in (0.1, 0.2, 0.3, 0.4, 0.6)]
        self.assertTrue(all(a > b for a, b in zip(at_1, at_1[1:])))

    def test_fig_3_5_order_is_the_opposite(self):
        doc = load("fig_3_5.json")
        self.assertEqual(doc["series_order_from_top"], [0.6, 0.4, 0.3, 0.2, 0.1])
        at_1 = [points(doc, v)[1.0] for v in (0.1, 0.2, 0.3, 0.4, 0.6)]
        self.assertTrue(all(a < b for a, b in zip(at_1, at_1[1:])))

    def test_fig_3_3_top_line_is_read_through_2_5_mpa(self):
        top = points(load("fig_3_3.json"), 0.1)
        self.assertTrue(all(top[p] is not None for p in PRESSURES))
        self.assertTrue(all(top[a] > top[b] for a, b in zip(PRESSURES, PRESSURES[1:])))

    def test_fig_3_3_merged_lines_are_null_with_a_reason(self):
        doc = load("fig_3_3.json")
        for speed in (0.2, 0.3, 0.4, 0.6):
            for p in (2.0, 2.5):
                self.assertIsNone(points(doc, speed)[p])
        for series in doc["calibration"]["series"]:
            if series["u0_minus_umf_m_s"] == 0.1:
                continue
            for sample in series["samples"]:
                if sample["P_MPa"] in (2.0, 2.5):
                    self.assertIn("分不开", sample["note"])
                    self.assertIsNone(sample["pixel_y_center"])

    def test_fig_3_5_five_lines_stay_separated_at_2_5_mpa(self):
        doc = load("fig_3_5.json")
        vals = [points(doc, v)[2.5] for v in (0.1, 0.2, 0.3, 0.4, 0.6)]
        self.assertTrue(all(v is not None for v in vals))
        self.assertTrue(all(a < b for a, b in zip(vals, vals[1:])))

    def test_every_sample_stores_the_pixel_used(self):
        for name in ("fig_3_3.json", "fig_3_5.json"):
            for series in load(name)["calibration"]["series"]:
                for sample in series["samples"]:
                    self.assertIn("pixel_x", sample)
                    self.assertIn("pixel_y_center", sample)
                    self.assertTrue(math.isfinite(sample["pixel_x"]))

    def test_edited_point_is_reported(self):
        doc = copy.deepcopy(load("fig_3_3.json"))
        doc["series"][0]["points"][2][1] = 0.5
        self.assertTrue(any("复算不一致" in p for p in DIGIT.check(doc)))

    def test_swapped_series_order_is_reported(self):
        doc = copy.deepcopy(load("fig_3_3.json"))
        doc["calibration"]["series"][0], doc["calibration"]["series"][-1] = (
            doc["calibration"]["series"][-1],
            doc["calibration"]["series"][0],
        )
        self.assertTrue(any("系列顺序" in p for p in DIGIT.check(doc)))

    def test_shifted_grid_is_reported_even_if_points_stay(self):
        doc = copy.deepcopy(load("fig_3_5.json"))
        doc["calibration"]["y_grid"][1][0] += 40
        self.assertTrue(any("y_grid" in p for p in DIGIT.check(doc)))

    def test_null_without_a_reason_is_reported(self):
        doc = copy.deepcopy(load("fig_3_3.json"))
        sample = doc["calibration"]["series"][1]["samples"][4]
        self.assertIsNone(sample["pixel_y_center"])
        sample["note"] = ""
        self.assertTrue(any("没有说明原因" in p for p in DIGIT.check(doc)))

    def test_error_must_be_numeric(self):
        for bad in (True, -1, float("nan"), float("inf"), None, "约 0,02"):
            doc = copy.deepcopy(load("fig_3_3.json"))
            doc["est_error"]["u_b_r"] = bad
            self.assertTrue(any("est_error.u_b_r" in p for p in DIGIT.check(doc)), bad)


if __name__ == "__main__":
    unittest.main()
