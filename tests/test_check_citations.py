# -*- coding: utf-8 -*-
"""引用校验的构造样例。不往仓库的 hamel/ 里放函数。"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "check_citations", ROOT / "tools" / "check_citations.py"
)
check_citations = importlib.util.module_from_spec(SPEC)
sys.modules["check_citations"] = check_citations
SPEC.loader.exec_module(check_citations)

INDEX = ROOT / "refs" / "hamel" / "index.json"

CORRECT = '''\
def rate_r6(temperature):
    """R6 与表 5.2–5.4。

    Source: Hamel Eq.5.39, S. 91
    Source: Hamel Tab.5.2, S. 88
    Source: Hamel Tab. 5.3, S. 91
    Source: Hamel Tab.5.4, S. 97
    """
    return temperature


def repeated_numbers(value):
    """重复编号写明页码后，主位置和 additional_locations 都可以。

    Source: Hamel Eq.3.49, S. 35
    Source: Hamel Eq.3.50, S. 36
    Source: Hamel Tab.5.1, S. 74
    Source: Hamel Tab.5.1, S. 79
    Source: Hamel Bild 5.1, S. 73
    """
    return value


def continuation_pages(value):
    """续页记在 additional_locations，没有 ambiguous_number 也是已登记页。

    Source: Hamel Tab.7.2, S. 134
    Source: Hamel Bild 7.10, S. 137
    Source: Hamel Bild 7.13, S. 141
    """
    return value


def _private_without_citation():
    return None


class Rates:
    def documented(self):
        """Source: Hamel Eq.5.40, S. 91"""
        return 1

    def _hidden(self):
        return 0
'''


def write_module(root: Path, relative: str, source: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(source, encoding="utf-8")


def codes(findings) -> list[str]:
    return [item.code for item in findings]


class CitationTests(unittest.TestCase):
    def check(self, source: str, relative: str = "kinetics/sample.py"):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_module(root, relative, source)
            return check_citations.check(root, INDEX)

    def test_correct_sample_passes(self):
        self.assertEqual(self.check(CORRECT), [])

    def test_repository_hamel_passes(self):
        self.assertEqual(check_citations.check(ROOT / "hamel", INDEX), [])

    def test_missing_citation(self):
        findings = self.check('def rate_r6(temperature):\n    """没有来源。"""\n    return temperature\n')
        self.assertIn("missing_citation", codes(findings))
        self.assertIn("rate_r6", findings[0].message)

    def test_unknown_equation(self):
        findings = self.check(
            'def rate_r6(temperature):\n    """Source: Hamel Eq.9.99, S. 1"""\n    return temperature\n'
        )
        self.assertIn("unknown_label", codes(findings))
        self.assertTrue(any("Eq.9.99" in item.message for item in findings))

    def test_page_mismatch(self):
        findings = self.check(
            'def rate_r6(temperature):\n    """Source: Hamel Eq.5.39, S. 90"""\n    return temperature\n'
        )
        self.assertIn("page_mismatch", codes(findings))
        self.assertTrue(any("S. 90" in item.message and "91" in item.message for item in findings))

    def test_ambiguous_equation_requires_page(self):
        findings = self.check(
            'def phi(value):\n    """Source: Hamel Eq.3.49"""\n    return value\n'
        )
        self.assertIn("ambiguous_page_required", codes(findings))
        self.assertTrue(any("Eq.3.49" in item.message for item in findings))

    def test_ambiguous_table_requires_page_and_accepts_either_location(self):
        missing = self.check(
            'def parameters():\n    """Source: Hamel Tab.5.1"""\n    return None\n'
        )
        self.assertIn("ambiguous_page_required", codes(missing))
        wrong = self.check(
            'def parameters():\n    """Source: Hamel Tab.5.1, S. 88"""\n    return None\n'
        )
        self.assertIn("page_mismatch", codes(wrong))

    def test_unknown_table_and_figure(self):
        findings = self.check(
            'def parameters():\n'
            '    """Source: Hamel Tab.9.99, S. 1\\n'
            '    Source: Hamel Bild 9.99, S. 1"""\n'
            '    return None\n'
        )
        messages = " ".join(item.message for item in findings)
        self.assertIn("Tab.9.99", messages)
        self.assertIn("Bild 9.99", messages)
        self.assertEqual(codes(findings).count("unknown_label"), 2)

    def test_registered_continuation_page_passes_and_unregistered_page_fails(self):
        accepted = self.check(
            'def figure():\n    """Source: Hamel Bild 7.10, S. 137"""\n    return None\n'
        )
        self.assertEqual(accepted, [])
        findings = self.check(
            'def figure():\n    """Source: Hamel Bild 7.10, S. 138"""\n    return None\n'
        )
        self.assertIn("page_mismatch", codes(findings))

    def test_missing_page_is_rejected(self):
        findings = self.check(
            'def rate_r6():\n    """Source: Hamel Eq.5.39"""\n    return None\n'
        )
        self.assertEqual(codes(findings), ["page_required"])
        self.assertTrue(any("Eq.5.39" in item.message for item in findings))
        # 下一行另一条引用的页码不能被这条缺页引用拿走。
        adjacent = self.check(
            'def rate_r6():\n'
            '    """Source: Hamel Eq.5.39\\n'
            '    Source: Hamel Eq.5.40, S. 91"""\n'
            '    return None\n'
        )
        self.assertEqual(codes(adjacent), ["page_required"])

    def test_broken_page_is_a_format_error(self):
        for citation in ("Source: Hamel Eq.5.39, S. -90", "Source: Hamel Eq.5.39, S. ninety"):
            findings = self.check(f'def rate_r6():\n    """{citation}"""\n    return None\n')
            self.assertEqual(codes(findings), ["page_format"], citation)
            self.assertTrue(any("Eq.5.39" in item.message for item in findings), citation)

    def test_public_method_without_citation(self):
        findings = self.check(
            "class Rates:\n"
            "    def rate_r6(self):\n"
            "        '''没有来源。'''\n"
            "        return 1\n"
        )
        self.assertIn("missing_citation", codes(findings))
        self.assertIn("Rates.rate_r6", findings[0].message)

    def test_other_modules_do_not_require_citations(self):
        findings = self.check(
            "def assemble():\n    return 1\n",
            relative="solver/sample.py",
        )
        self.assertEqual(findings, [])
        bad = self.check(
            'def assemble():\n    """Source: Hamel Eq.9.99, S. 1"""\n    return 1\n',
            relative="cases/sample.py",
        )
        self.assertEqual(codes(bad), ["unknown_label"])

    def test_cli_reports_missing_citation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_module(
                root,
                "freeboard/sample.py",
                "def throw_off():\n    return 1\n",
            )
            proc = subprocess.run(
                [sys.executable, str(ROOT / "tools" / "check_citations.py"), "--hamel", str(root), "--index", str(INDEX)],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(proc.returncode, 1)
        self.assertIn("missing_citation", proc.stdout)
        self.assertIn("throw_off", proc.stdout)


if __name__ == "__main__":
    unittest.main()
