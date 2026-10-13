# -*- coding: utf-8 -*-
"""T2 原文导航测试；锚点取自当前 PDF 的页码、印刷标签和图注，不使用模型输出。"""
import importlib.util
import json
import shutil
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

    def test_label_bbox_is_in_page_image_frame(self):
        # Eq. 2.1（PDF 17 右半页）：pdfplumber 原坐标 x0=491.7、top=325.21 pt，页 bbox 原点 (25.1, -28.3)，
        # 平移到图像坐标后为 (466.6, 353.5)；300 dpi 原图上该公式行目视位于约 1440–1520 px。
        x0, y0, _, _ = self.index["equations"]["2.1"]["label_bbox_pt"]
        self.assertAlmostEqual(x0, 466.6, delta=0.5)
        self.assertAlmostEqual(y0, 353.5, delta=0.5)
        self.assertTrue(1440 < y0 * 300 / 72 < 1520)
        for kind in ("equations", "tables", "figures"):
            for label, entry in self.index[kind].items():
                if isinstance(entry, dict) and entry.get("label_bbox_pt"):
                    with self.subTest(kind=kind, label=label):
                        bx0, by0, bx1, by1 = entry["label_bbox_pt"]
                        self.assertTrue(0 <= bx0 < bx1 and 0 <= by0 < by1 <= 560)

    def test_repeated_equation_numbers_do_not_hide_second_formula(self):
        for label in ("3.49", "3.50"):
            entry = self.index["equations"][label]
            self.assertTrue(entry["ambiguous_number"])
            self.assertEqual(entry["printed_page"], "35")
            self.assertEqual(entry["additional_locations"][0]["printed_page"], "36")
            self.assertNotEqual(entry["section"], entry["additional_locations"][0]["section"])

    def test_continuation_pages_and_supplement_figures(self):
        for label, page, printed in [("7.10", 78, "137"), ("7.13", 80, "141")]:
            extra = self.index["figures"][label]["additional_locations"]
            self.assertIn((page, printed), [(v["pdf_page"], v["printed_page"]) for v in extra])
        self.assertEqual(self.index["tables"]["7.2"]["additional_locations"][0]["printed_page"], "134")
        expected = {
            "7.20": ("152", "left", "Feuchtegehalts"),
            "7.21": ("153", "right", "Sekundärluft"),
        }
        for label, (printed, side, phrase) in expected.items():
            entry = self.index["figures"][label]
            self.assertEqual(entry["status"], "indexed")
            self.assertEqual((entry["pdf_page"], entry["printed_page"], entry["side"]),
                             (114, printed, side))
            self.assertEqual(entry["section"], "7.2")
            self.assertIn(phrase, entry["caption_de"])
            self.assertEqual(entry["caption_source"], "page_image")
        self.assertEqual(self.index["page_map"]["114"], {"left": "152", "right": "153"})
        self.assertEqual(self.index["printed_page_map"]["152"][0]["pdf_page"], 114)
        self.assertEqual(self.index["printed_page_map"]["153"][0]["side"], "right")
        self.assertEqual(self.index["coverage"]["7"]["located_figures"], 21)
        self.assertEqual(self.index["coverage"]["7"]["missing_figures"], 0)
        kinds = [item["kind"] for item in self.index["source_issues"]]
        self.assertIn("supplement_scan", kinds)
        self.assertNotIn("missing_pages", kinds)

    def test_duplicate_and_rotated_scan_page_map(self):
        self.assertEqual(self.index["page_map"]["84"], self.index["page_map"]["85"])
        self.assertEqual(len(self.index["printed_page_map"]["148"]), 2)
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

    def test_skip_render_does_not_call_renderer_or_touch_png(self):
        called = []
        original = CORPUS.render_supplement_png
        CORPUS.render_supplement_png = lambda *args: called.append(args)
        try:
            missing = Path(tempfile.gettempdir()) / "hamel-skip-missing-114.png"
            if missing.exists():
                missing.unlink()
            note = CORPUS.place_supplement_png(Path("unused.pdf"), missing, True)
            self.assertFalse(missing.exists())
            self.assertEqual(called, [])
            self.assertIn("没有", note)
            kept = Path(tempfile.gettempdir()) / "hamel-skip-kept-114.png"
            kept.write_bytes(b"original-png")
            note = CORPUS.place_supplement_png(Path("unused.pdf"), kept, True)
            self.assertEqual(kept.read_bytes(), b"original-png")
            self.assertEqual(called, [])
            self.assertIn("沿用", note)
            CORPUS.place_supplement_png(Path("unused.pdf"), kept, False)
            self.assertEqual(len(called), 1)
        finally:
            CORPUS.render_supplement_png = original
            kept = Path(tempfile.gettempdir()) / "hamel-skip-kept-114.png"
            if kept.exists():
                kept.unlink()

    def test_cli_skip_render_does_not_rewrite_png(self):
        pdf = ROOT / "docs" / "hamel_supplement_S152-153.pdf"
        if not pdf.exists():
            self.skipTest("本地没有补页 PDF")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            pages = output / "pages"
            pages.mkdir()
            png = pages / "114.png"
            png.write_bytes(b"do-not-rewrite")
            shutil.copy(ROOT / "refs" / "hamel" / "index.json", output / "index.json")
            command = [sys.executable, str(ROOT / "tools" / "build_hamel_corpus.py"),
                       "--supplement-only", "--skip-render", "--supplement", str(pdf),
                       "--output", str(output)]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(png.read_bytes(), b"do-not-rewrite")
            index = json.loads((output / "index.json").read_text(encoding="utf-8"))
            self.assertIn("未渲染，沿用已有", index["_meta"]["supplement"]["render"])
            png.unlink()
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(png.exists())
            index = json.loads((output / "index.json").read_text(encoding="utf-8"))
            self.assertIn("本地没有", index["_meta"]["supplement"]["render"])


if __name__ == "__main__":
    unittest.main()
