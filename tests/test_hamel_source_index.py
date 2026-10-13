# -*- coding: utf-8 -*-
"""T2 原文导航测试；锚点取自当前 PDF 的页码、印刷标签和图注，不使用模型输出。"""
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("corpus", ROOT / "tools/build_hamel_corpus.py")
CORPUS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CORPUS)

# 原图抽查锚点：类别、标签、PDF页（1起）、印刷页、所在半页。
SOURCE_ANCHORS = [
    ("equations", "2.3", 18, "16", "left"),
    ("equations", "3.2", 20, "21", "right"),
    ("equations", "3.4", 21, "22", "left"),
    ("equations", "3.13", 23, "26", "left"),
    ("equations", "3.40", 24, "29", "right"),
    ("equations", "3.44", 26, "32", "left"),
    ("equations", "3.45", 27, "34", "left"),
    ("equations", "3.48", 27, "34", "left"),
    ("equations", "3.56", 29, "38", "left"),
    ("equations", "3.63", 29, "39", "right"),
    ("equations", "3.72", 30, "41", "right"),
    ("equations", "3.75", 31, "43", "right"),
    ("equations", "3.82", 32, "44", "left"),
    ("equations", "3.84", 32, "44", "left"),
    ("equations", "3.86", 32, "45", "right"),
    ("equations", "4.13", 40, "60", "left"),
    ("equations", "4.21", 41, "63", "right"),
    ("equations", "4.22", 42, "64", "left"),
    ("equations", "4.25", 43, "67", "right"),
    ("equations", "4.30", 45, "70", "left"),
]


class SourceIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = json.loads((ROOT / "refs/hamel/index.json").read_text(encoding="utf-8"))

    def test_twenty_source_image_anchors(self):
        for kind, label, pdf_page, printed, side in SOURCE_ANCHORS:
            with self.subTest(label=label):
                entry = self.index[kind][label]
                self.assertEqual((entry["pdf_page"], entry["printed_page"], entry["side"]),
                                 (pdf_page, printed, side))

    def test_repeated_equation_numbers_do_not_hide_second_formula(self):
        for label in ("3.49", "3.50"):
            entry = self.index["equations"][label]
            self.assertTrue(entry["ambiguous_number"])
            self.assertEqual(entry["printed_page"], "35")
            self.assertEqual(entry["additional_locations"][0]["printed_page"], "36")
            self.assertNotEqual(entry["section"], entry["additional_locations"][0]["section"])

    def test_continuation_pages_and_absent_figures(self):
        for label, page, printed in [("7.10", 78, "137"), ("7.13", 80, "141")]:
            extra = self.index["figures"][label]["additional_locations"]
            self.assertIn((page, printed), [(v["pdf_page"], v["printed_page"]) for v in extra])
        self.assertEqual(self.index["tables"]["7.2"]["additional_locations"][0]["printed_page"], "134")
        for label in ("7.20", "7.21"):
            entry = self.index["figures"][label]
            self.assertEqual(entry["status"], "source_page_missing")
            self.assertIsNone(entry["pdf_page"])
            self.assertIsNone(entry["caption_de"])

    def test_duplicate_and_rotated_scan_page_map(self):
        self.assertEqual(self.index["page_map"]["84"], self.index["page_map"]["85"])
        self.assertEqual(len(self.index["printed_page_map"]["148"]), 2)
        self.assertNotIn("152", self.index["printed_page_map"])
        self.assertEqual(self.index["page_map"]["98"], {"top": "176", "bottom": "177"})

    def test_author_matching_uses_whole_names(self):
        self.assertEqual(CORPUS.author_positions("Koks sowie So (1995)", ["KO"]), [])
        self.assertEqual(CORPUS.author_positions("sowie", ["SO"]), [])
        self.assertEqual(CORPUS.author_positions("Ko et al. (1988)", ["KO"]), [(0, 2)])
        self.assertEqual(CORPUS.author_positions("Al Taweel (1989)", ["ALTAWEEL"]), [(0, 9)])

    def test_different_pdf_is_rejected_before_export(self):
        with tempfile.TemporaryDirectory() as directory:
            pdf = Path(directory) / "other.pdf"
            pdf.write_bytes(b"different source")
            output = Path(directory) / "output"
            result = subprocess.run([sys.executable, str(ROOT / "tools/build_hamel_corpus.py"),
                                     str(pdf), "--output", str(output)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("SHA256", result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
