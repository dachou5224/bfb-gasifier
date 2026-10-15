# -*- coding: utf-8 -*-
"""出口读数复算的失败路径。改副本，不改 data/digitized 里的正式文件。"""
import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "digitize", ROOT / "tools" / "digitize_wsv400_outlets.py")
DIGIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIGIT)


def load_pair():
    fig717 = json.loads((ROOT / "data" / "digitized" / "fig_7_17.json").read_text(encoding="utf-8"))
    fig716 = json.loads((ROOT / "data" / "digitized" / "fig_7_16.json").read_text(encoding="utf-8"))
    return fig717, fig716


def load_710():
    return json.loads((ROOT / "data" / "digitized" / "fig_7_10.json").read_text(encoding="utf-8"))


def case_710(doc, label):
    return next(case for case in doc["cases"] if case["case_label_in_figure"] == label)


class OutletCheckTests(unittest.TestCase):
    def test_recorded_files_pass(self):
        self.assertEqual(DIGIT.check(), [])

    def test_non_outlet_height_fails(self):
        fig717, fig716 = load_pair()
        fig716["series"][0]["points"] = [[0.5, 1000]]
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("ξ = 1", text)

    def test_duplicate_series_fails(self):
        fig717, fig716 = load_pair()
        fig716["series"].append(copy.deepcopy(fig716["series"][0]))
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("重复", text)

    def test_multiple_points_fail(self):
        fig717, fig716 = load_pair()
        fig716["series"][0]["points"] = [[1.0, 1000], [0.5, 900]]
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("ξ = 1", text)

    def test_empty_est_error_fails(self):
        fig717, fig716 = load_pair()
        fig717["est_error"] = {}
        fig716["est_error"] = {}
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("7.17 est_error 缺少", text)
        self.assertIn("7.16 est_error 缺少", text)

    def test_missing_units_fail(self):
        fig717, fig716 = load_pair()
        fig716.pop("units")
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("7.16 缺少 units", text)

    def test_illegal_error_values_fail(self):
        fig717, fig716 = load_pair()
        fig717["est_error"]["T_FRB_K"] = True
        fig716["est_error"]["T"] = float("nan")
        fig716["est_error"]["CO"] = float("inf")
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("7.17 est_error.T_FRB_K", text)
        self.assertIn("7.16 est_error.T", text)
        self.assertIn("7.16 est_error.CO", text)

    def test_negative_error_fails(self):
        fig717, fig716 = load_pair()
        fig717["est_error"]["CH4"] = -0.001
        text = "\n".join(DIGIT.check(fig717, fig716))
        self.assertIn("7.17 est_error.CH4", text)

    def test_710_official_without_center_fails(self):
        fig717, fig716 = load_pair()
        fig710 = load_710()
        del fig710["calibration"]["centers"]["CO"]["1"]
        case_710(fig710, "1")["CO"][1] = 0.299
        text = "\n".join(DIGIT.check(fig717, fig716, fig710))
        self.assertIn("7.10 CO 1: 正式读数没有对应的标定中心", text)

    def test_710_center_without_official_fails(self):
        fig717, fig716 = load_pair()
        fig710 = load_710()
        del case_710(fig710, "1")["CO"]
        text = "\n".join(DIGIT.check(fig717, fig716, fig710))
        self.assertIn("7.10 CO 1: 标定中心没有对应的正式读数", text)

    def test_710_duplicate_label_fails(self):
        fig717, fig716 = load_pair()
        fig710 = load_710()
        extra = copy.deepcopy(case_710(fig710, "1"))
        extra["CO"] = [0.999, 0.999]
        fig710["cases"].insert(0, extra)
        text = "\n".join(DIGIT.check(fig717, fig716, fig710))
        self.assertIn("7.10 case_label_in_figure 1 重复", text)
        self.assertIn("7.10 case_id CASE_VTT_PRESSURIZED_1 重复", text)
        self.assertIn("[0.999, 0.999]", text)

    def test_710_quantity_absent_on_both_sides_passes(self):
        fig717, fig716 = load_pair()
        fig710 = load_710()
        del case_710(fig710, "1")["CO"]
        del fig710["calibration"]["centers"]["CO"]["1"]
        self.assertEqual(DIGIT.check(fig717, fig716, fig710), [])


if __name__ == "__main__":
    unittest.main()
