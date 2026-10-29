# -*- coding: utf-8 -*-
"""偏离检查的违规样例：阶段后缀、未登记配置字段，以及章程 §9 的其余信号。"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("drift_check", ROOT / "tools" / "drift_check.py")
drift_check = importlib.util.module_from_spec(SPEC)
sys.modules["drift_check"] = drift_check
SPEC.loader.exec_module(drift_check)

TODAY = date(2026, 10, 13)
STATUS_TEMPLATE = """# 项目状态（STATUS）

## 工况指标

上次草稿写过 2099-01-01，以表内更新日期为准。

| 工况 | 收敛 | T | CO | CO₂ | H₂ | H₂O | CH₄ | X_C | 通过 | 更新日期 |
|---|---|---|---|---|---|---|---|---|---|---|
| WSV400 A | — | — | — | — | — | — | — | — | 否 | {day} |
| WSV400 C | — | — | — | — | — | — | — | — | 否 | 2020-01-01 |

## 其他

| 更新日期 |
|---|
| 2099-01-01 |
"""


def py_lines(count: int) -> str:
    return "".join(f"v{index} = {index}\n" for index in range(count))


def signal_codes(report) -> list[str]:
    return [signal.code for signal in report.signals]


class ParseTests(unittest.TestCase):
    def test_dict_class_and_private_fields(self):
        source = """
PHYSICAL = {"bed_diameter_m": 0.4, "_internal": 1}
PHYSICAL = dict(bed_diameter_m=0.4, fuel_hhv=1)
SOLVER = {"damping": 0.5}

class PhysicalConfig:
    particle_density: float
    _hidden = 1

    def ignored(self):
        self.not_a_field = 1

class SolverConfig:
    max_iterations = 20
"""
        physical, solver = drift_check.extract_config_fields(source)
        self.assertEqual(physical, {"bed_diameter_m", "fuel_hhv", "particle_density"})
        self.assertEqual(solver, {"damping", "max_iterations"})

    def test_non_literal_config_is_unreadable(self):
        with self.assertRaises(drift_check.ConfigError):
            drift_check.extract_config_fields("PHYSICAL = load_config()\n")

    def test_later_mutation_is_unreadable(self):
        samples = (
            'PHYSICAL = {}\nPHYSICAL.update({"invented_knob": 1})\n',
            'PHYSICAL = {}\nPHYSICAL["invented_knob"] = 1\n',
            'SOLVER = {"damping": 0.5}\nSOLVER |= {"new_tol": 1e-8}\n',
        )
        for source in samples:
            with self.subTest(source=source):
                with self.assertRaises(drift_check.ConfigError):
                    drift_check.extract_config_fields(source)

    def test_stage_suffix_matches_path_components(self):
        self.assertTrue(drift_check.has_stage_suffix("hamel/kinetics/rate_v7.py"))
        self.assertTrue(drift_check.has_stage_suffix("notes_r12_r5.md"))
        self.assertTrue(drift_check.has_stage_suffix("hamel/freeboard_stage3.py"))
        self.assertTrue(drift_check.has_stage_suffix("hamel/kinetics_v7/rate.py"))
        self.assertFalse(drift_check.has_stage_suffix("hamel/hydrodynamics.py"))
        self.assertFalse(drift_check.has_stage_suffix("hamel/validation.py"))
        self.assertFalse(drift_check.has_stage_suffix("hamel/freeboard_stage.py"))

    def test_adr_match_is_a_whole_identifier(self):
        self.assertTrue(drift_check.field_in_adr("bed_diameter_m", "采用 bed_diameter_m。"))
        self.assertFalse(drift_check.field_in_adr("damping", "见 newton_damping 与 the"))
        self.assertFalse(drift_check.field_in_adr("t", "the temperature"))

    def test_metrics_date_ignores_decoys_outside_the_table(self):
        found = drift_check.parse_metrics_date(STATUS_TEMPLATE.format(day="2026-10-13"))
        self.assertEqual(found, date(2026, 10, 13))

    def test_repository_status_table(self):
        text = (ROOT / "docs" / "STATUS.md").read_text(encoding="utf-8")
        found = drift_check.parse_metrics_date(text)
        self.assertIsNotNone(found)
        section = text.split("## 工况指标", 1)[1].split("\n## ", 1)[0]
        self.assertIn(found.isoformat(), section)


class DriftRepoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.env = os.environ.copy()
        self.env.update(
            {
                "GIT_AUTHOR_NAME": "drift-check",
                "GIT_AUTHOR_EMAIL": "drift-check@example.com",
                "GIT_COMMITTER_NAME": "drift-check",
                "GIT_COMMITTER_EMAIL": "drift-check@example.com",
                "GIT_CONFIG_NOSYSTEM": "1",
                "GIT_CONFIG_GLOBAL": os.devnull,
            }
        )
        self.git(["init", "-b", "main"])
        self.write_status("2026-10-13")
        self.commit("init")

    def git(self, args: list[str]) -> subprocess.CompletedProcess:
        result = subprocess.run(
            ["git", *args],
            cwd=self.root,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            raise AssertionError(f"git {args} 失败：{result.stderr}")
        return result

    def write(self, relative: str, text: str) -> None:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def write_status(self, day: str) -> None:
        self.write("docs/STATUS.md", STATUS_TEMPLATE.format(day=day))

    def commit(self, message: str) -> None:
        self.git(["add", "-A"])
        self.git(["commit", "-m", message])

    def branch(self) -> None:
        self.git(["checkout", "-b", "task/example"])

    def check(self, today: date = TODAY):
        return drift_check.check(self.root, "main", today)

    def add_case_field(self, name: str, field: dict, *, solver: bool = False):
        self.branch()
        group = "SOLVER" if solver else "PHYSICAL"
        self.write("hamel/config.py", f'{group} = {{"{name}": 1}}\n')
        self.write("data/cases/wsv400.json", json.dumps({name: field}, ensure_ascii=False))
        self.commit("case field")
        return self.check()

    def test_clean_tree_has_no_signals(self):
        report = self.check()
        self.assertEqual(signal_codes(report), [])
        text = drift_check.format_text(report)
        self.assertIn("结论：未触发偏离信号。", text)
        self.assertIn("物理参数 0 → 0", text)
        self.assertIn("工况指标最近变化 2026-10-13，距检查日 0 天", text)

    def test_dirty_worktree(self):
        self.write("README.md", "改动\n")
        report = self.check()
        self.assertEqual(signal_codes(report), ["dirty"])
        self.assertIn("工作树不干净", drift_check.format_text(report))

    def test_unstaged_edits_count_toward_net_lines(self):
        self.write("hamel/a.py", py_lines(2))
        self.write("hamel/b.py", py_lines(3))
        self.commit("two files")
        self.write("hamel/a.py", py_lines(12))
        report = self.check()
        self.assertEqual(report.net_lines, 10)
        self.assertEqual(report.hamel_lines, 15)
        self.assertEqual(signal_codes(report), ["dirty"])

    def test_deleted_code_reduces_net_lines(self):
        self.write("hamel/core.py", py_lines(10))
        self.commit("core")
        self.branch()
        (self.root / "hamel" / "core.py").unlink()
        self.commit("drop")
        report = self.check()
        self.assertEqual(report.net_lines, -10)
        self.assertEqual(report.hamel_lines, 0)
        self.assertEqual(signal_codes(report), [])

    def test_net_line_limit_is_exclusive_at_800(self):
        self.branch()
        self.write("hamel/core.py", py_lines(800))
        self.commit("boundary")
        self.assertEqual(signal_codes(self.check()), [])
        self.write("hamel/core.py", py_lines(801))
        self.commit("over")
        report = self.check()
        self.assertEqual(report.net_lines, 801)
        self.assertEqual(report.hamel_lines, 801)
        self.assertEqual(signal_codes(report), ["net_lines"])
        self.assertIn("超过 800 行", drift_check.format_text(report))

    def test_status_approval_raises_the_net_limit_only_on_that_branch(self):
        self.git(["checkout", "-b", "task/T25a-assembly"])
        self.write(
            "docs/STATUS.md",
            STATUS_TEMPLATE.format(day="2026-10-13")
            + "\nT25a 净增代码上限 1000 行，目标所有者 2026-10-20 派发前批准\n",
        )
        self.write("hamel/core.py", py_lines(1000))
        self.commit("approved ceiling")
        report = self.check()
        self.assertEqual(report.net_lines, 1000)
        self.assertEqual(report.approved_net_limit, 1000)
        self.assertEqual(signal_codes(report), [])
        self.assertIn("批准上限 1000", drift_check.format_text(report))
        self.write("hamel/core.py", py_lines(1001))
        self.commit("over approval")
        report = self.check()
        self.assertEqual(signal_codes(report), ["net_lines"])
        self.assertIn("超过批准上限 1000 行", drift_check.format_text(report))
        self.git(["checkout", "-b", "task/other"])
        report = self.check()
        self.assertIsNone(report.approved_net_limit)
        self.assertEqual(signal_codes(report), ["net_lines"])
        self.assertIn("超过 800 行", drift_check.format_text(report))

    def test_hamel_total_limit_is_exclusive_at_11000(self):
        self.write("hamel/core.py", py_lines(11000))
        self.commit("bulk")
        self.assertEqual(self.check().hamel_lines, 11000)
        self.assertEqual(signal_codes(self.check()), [])
        self.branch()
        self.write("hamel/extra.py", py_lines(1))
        self.commit("one")
        report = self.check()
        self.assertEqual(report.net_lines, 1)
        self.assertEqual(report.hamel_lines, 11001)
        self.assertEqual(signal_codes(report), ["hamel_lines"])

    def test_untracked_hamel_code_counts(self):
        self.write("hamel/sneak.py", py_lines(801))
        report = self.check()
        self.assertEqual(report.net_lines, 801)
        self.assertEqual(signal_codes(report), ["dirty", "net_lines"])

    def test_tests_under_hamel_do_not_count_as_code(self):
        self.branch()
        self.write("hamel/tests/huge_v7.py", py_lines(900))
        self.write("hamel/test_hidden.py", py_lines(900))
        self.commit("tests")
        report = self.check()
        self.assertEqual(report.net_lines, 0)
        self.assertEqual(report.hamel_lines, 0)
        self.assertEqual(signal_codes(report), ["stage_suffix"])
        self.assertEqual(report.stage_files, ["hamel/tests/huge_v7.py"])

    def test_stage_suffix_filenames_are_reported(self):
        self.branch()
        self.write("hamel/kinetics/rate_v7.py", "x = 1\n")
        self.write("hamel/freeboard_stage3.py", "x = 1\n")
        self.write("notes_r12_r5.md", "记录\n")
        self.write("hamel/hydrodynamics.py", "x = 1\n")
        self.commit("names")
        report = self.check()
        self.assertEqual(signal_codes(report), ["stage_suffix"])
        self.assertEqual(
            report.stage_files,
            ["hamel/freeboard_stage3.py", "hamel/kinetics/rate_v7.py", "notes_r12_r5.md"],
        )
        text = drift_check.format_text(report)
        for name in report.stage_files:
            self.assertIn(name, text)

    def test_existing_stage_suffix_on_main_is_not_new(self):
        self.write("legacy_v7.py", "x = 1\n")
        self.commit("legacy")
        report = self.check()
        self.assertEqual(report.stage_files, [])
        self.assertEqual(signal_codes(report), [])

    def test_renamed_file_uses_the_new_name(self):
        self.write("hamel/keep.py", "x = 1\n")
        self.commit("keep")
        self.branch()
        self.git(["mv", "hamel/keep.py", "hamel/keep_v7.py"])
        self.commit("rename")
        report = self.check()
        self.assertEqual(report.stage_files, ["hamel/keep_v7.py"])
        self.assertEqual(signal_codes(report), ["stage_suffix"])

    def test_new_config_field_without_adr(self):
        self.branch()
        self.write(
            "hamel/config.py",
            'PHYSICAL = {"bed_diameter_m": 0.4, "invented_knob": 1}\nSOLVER = {"damping": 0.5}\n',
        )
        self.write("docs/adr/0001-known.md", "采用 bed_diameter_m 与 damping。\n")
        self.write("docs/adr/CANDIDATES.md", "invented_knob 仍待核实，不能当作 ADR。\n")
        self.commit("config")
        report = self.check()
        self.assertEqual(report.unregistered_fields, ["invented_knob"])
        self.assertEqual(report.physical_base, 0)
        self.assertEqual(report.physical_head, 2)
        self.assertEqual(report.solver_head, 1)
        self.assertEqual(signal_codes(report), ["config_without_adr"])
        text = drift_check.format_text(report)
        self.assertIn("没有对应 ADR：invented_knob", text)
        self.assertIn("物理参数 0 → 2", text)

    def test_new_field_with_formal_adr_is_accepted(self):
        self.branch()
        self.write("hamel/config.py", 'PHYSICAL = {"invented_knob": 1}\n')
        self.write("docs/adr/0002-knob.md", "解释性选择：invented_knob。\n")
        self.commit("adr")
        report = self.check()
        self.assertEqual(signal_codes(report), [])
        self.assertIn("新增 invented_knob", drift_check.format_text(report))

    def test_field_already_on_main_is_not_treated_as_new(self):
        self.write("hamel/config.py", 'PHYSICAL = {"orphan_knob": 1}\n')
        self.commit("old knob")
        report = self.check()
        self.assertEqual(report.added_physical, [])
        self.assertEqual(signal_codes(report), [])

    def test_config_files_are_unioned(self):
        self.write("hamel/config.py", 'PHYSICAL = {"bed_diameter_m": 0.4}\n')
        self.write("docs/adr/0001-diameter.md", "bed_diameter_m\n")
        self.commit("base config")
        self.branch()
        self.write("hamel/cases/config.py", 'PHYSICAL = {"bed_diameter_m": 0.4, "fuel_hhv": 1}\n')
        self.write("docs/adr/0002-hhv.md", "fuel_hhv\n")
        self.commit("second file")
        report = self.check()
        self.assertEqual(report.physical_head, 2)
        self.assertEqual(report.added_physical, ["fuel_hhv"])
        self.assertEqual(signal_codes(report), [])

    def test_physical_field_limit(self):
        names = [f"field_{index:02d}" for index in range(31)]
        body = "PHYSICAL = {" + ", ".join(f'"{name}": {index}' for index, name in enumerate(names)) + "}\n"
        self.branch()
        self.write("hamel/config.py", body)
        self.write("docs/adr/0001-fields.md", " ".join(names) + "\n")
        self.commit("thirty one")
        report = self.check()
        self.assertEqual(report.physical_head, 31)
        self.assertEqual(signal_codes(report), ["physical_config_limit"])
        self.assertIn("超过 30 项", drift_check.format_text(report))

        self.write(
            "hamel/config.py",
            "PHYSICAL = {" + ", ".join(f'"{name}": 1' for name in names[:30]) + "}\n",
        )
        self.commit("thirty")
        report = self.check()
        self.assertEqual(report.physical_head, 30)
        self.assertEqual(signal_codes(report), [])

    def test_solver_fields_do_not_use_the_physical_limit(self):
        names = [f"solver_{index:02d}" for index in range(31)]
        body = "SOLVER = {" + ", ".join(f'"{name}": {index}' for index, name in enumerate(names)) + "}\n"
        self.branch()
        self.write("hamel/config.py", body)
        self.write("docs/adr/0001-solver.md", " ".join(names) + "\n")
        self.commit("solver")
        report = self.check()
        self.assertEqual(report.solver_head, 31)
        self.assertEqual(report.physical_head, 0)
        self.assertEqual(signal_codes(report), [])

    def test_unreadable_config_is_reported(self):
        self.branch()
        self.write("hamel/config.py", "PHYSICAL = {\n")
        self.commit("broken")
        report = self.check()
        self.assertEqual(signal_codes(report), ["config_unreadable"])
        self.assertIn("hamel/config.py", drift_check.format_text(report))

    def test_config_update_after_assignment_is_reported(self):
        self.branch()
        self.write(
            "hamel/config.py",
            'PHYSICAL = {}\nPHYSICAL.update({"invented_knob": 1})\n',
        )
        self.commit("update")
        report = self.check()
        self.assertEqual(report.physical_head, 0)
        self.assertEqual(report.unregistered_fields, [])
        self.assertEqual(signal_codes(report), ["config_unreadable"])
        self.assertIn("invented_knob", self.root.joinpath("hamel/config.py").read_text(encoding="utf-8"))
        self.assertIn("无法静态读取", drift_check.format_text(report))

    def test_main_advancing_does_not_offset_this_card(self):
        self.branch()
        self.write("hamel/task_card.py", py_lines(899))
        self.write("hamel/config.py", 'PHYSICAL = {"invented_knob": 1}\n')
        self.write("card_v7.py", "x = 1\n")
        self.commit("task card")
        self.git(["checkout", "main"])
        self.write("hamel/other_card.py", py_lines(600))
        self.write("hamel/config.py", 'PHYSICAL = {"from_main": 1}\n')
        self.write("main_stage3.py", "x = 1\n")
        self.commit("other card")
        self.git(["checkout", "task/example"])
        report = self.check()
        fork = self.git(["merge-base", "main", "HEAD"]).stdout.strip()
        main_tip = self.git(["rev-parse", "main"]).stdout.strip()
        self.assertEqual(report.compared, fork)
        self.assertNotEqual(fork, main_tip)
        self.assertEqual(report.net_lines, 900)
        self.assertEqual(report.hamel_lines, 900)
        self.assertEqual(report.unregistered_fields, ["invented_knob"])
        self.assertEqual(report.stage_files, ["card_v7.py"])
        self.assertEqual(signal_codes(report), ["net_lines", "config_without_adr", "stage_suffix"])
        self.assertIn(fork[:12], drift_check.format_text(report))

    def test_default_base_prefers_origin_main_and_notes_divergence(self):
        self.write("README.md", "remote\n")
        self.commit("remote")
        self.git(["update-ref", "refs/remotes/origin/main", "HEAD"])
        self.git(["reset", "--hard", "HEAD~1"])
        self.assertEqual(drift_check.default_base(self.root), "origin/main")
        report = drift_check.check(self.root, "origin/main", TODAY)
        self.assertEqual(signal_codes(report), [])
        self.assertTrue(any("落后 origin/main 1 个提交" in note for note in report.base_notes))
        text = drift_check.format_text(report)
        self.assertIn("不会自动 fetch", text)
        self.assertIn("落后 origin/main 1 个提交", text)

    def test_metrics_stale_after_14_days(self):
        self.write_status("2026-09-29")
        self.commit("old")
        report = self.check()
        self.assertEqual(report.metrics_age_days, 14)
        self.assertEqual(signal_codes(report), ["metrics_stale"])
        self.write_status("2026-09-30")
        self.commit("newer")
        report = self.check()
        self.assertEqual(report.metrics_age_days, 13)
        self.assertEqual(signal_codes(report), [])

    def test_sourced_physical_field_does_not_need_adr(self):
        report = self.add_case_field("bed_diameter_m", {
            "value": 0.4, "unit": "m", "source": "表 7.4, S. 143", "status": "source",
        })
        self.assertEqual(signal_codes(report), [])

    def test_source_without_page_still_needs_adr(self):
        report = self.add_case_field("bed_diameter_m", {
            "value": 0.4, "unit": "m", "source": "表 7.4", "status": "source",
        })
        self.assertEqual(report.unregistered_fields, ["bed_diameter_m"])

    def test_placeholder_source_still_needs_adr(self):
        report = self.add_case_field("bed_diameter_m", {
            "value": None, "source": "待补", "status": "source",
        })
        self.assertEqual(report.unregistered_fields, ["bed_diameter_m"])

    def test_derived_without_formula_still_needs_adr(self):
        report = self.add_case_field("fuel_hhv", {
            "value": 1, "source": "表 7.4, S. 143", "status": "derived",
        })
        self.assertEqual(report.unregistered_fields, ["fuel_hhv"])

    def test_blank_formula_still_needs_adr(self):
        report = self.add_case_field("fuel_hhv", {
            "value": 1, "source": "表 7.4, S. 143", "status": "derived", "formula": "  ",
        })
        self.assertEqual(report.unregistered_fields, ["fuel_hhv"])

    def test_derived_with_page_and_formula_is_exempt(self):
        report = self.add_case_field("fuel_hhv", {
            "value": 1, "source": "表 7.4, S. 143", "status": "derived", "formula": "Q = V / 22.4",
        })
        self.assertEqual(signal_codes(report), [])

    def test_missing_status_still_needs_adr(self):
        report = self.add_case_field("bed_height_m", {
            "value": None, "source": "查过 S. 142", "status": "missing",
        })
        self.assertEqual(report.unregistered_fields, ["bed_height_m"])

    def test_mixed_case_copies_do_not_exempt(self):
        self.branch()
        self.write("hamel/config.py", 'PHYSICAL = {"bed_diameter_m": 0.4}\n')
        self.write("data/cases/wsv400.json", json.dumps({
            "shared": {"bed_diameter_m": {"value": 0.4, "source": "表 7.4, S. 143", "status": "source"}},
            "A": {"bed_diameter_m": {"value": None, "source": "待补", "status": "missing"}},
        }, ensure_ascii=False))
        self.commit("mixed")
        report = self.check()
        self.assertEqual(report.unregistered_fields, ["bed_diameter_m"])

    def test_solver_field_is_not_exempted_by_case_file(self):
        report = self.add_case_field("bed_diameter_m", {
            "value": 0.4, "source": "表 7.4, S. 143", "status": "source",
        }, solver=True)
        self.assertEqual(report.unregistered_fields, ["bed_diameter_m"])
        self.assertEqual(signal_codes(report), ["config_without_adr"])

    def test_identical_skill_copies_are_quiet(self):
        self.write(".claude/skills/hamel-source-check/SKILL.md", "同样\n")
        self.write(".cursor/skills/hamel-source-check/SKILL.md", "同样\n")
        self.commit("skills")
        self.assertEqual(signal_codes(self.check()), [])

    def test_skill_content_mismatch(self):
        self.write(".claude/skills/hamel-source-check/SKILL.md", "甲\n")
        self.write(".cursor/skills/hamel-source-check/SKILL.md", "乙\n")
        self.commit("mismatch")
        report = self.check()
        self.assertEqual(signal_codes(report), ["skill_mismatch"])
        self.assertIn("hamel-source-check：内容不一致", drift_check.format_text(report))

    def test_skill_present_on_only_one_side(self):
        self.write(".claude/skills/hamel-source-check/SKILL.md", "甲\n")
        self.commit("claude only")
        report = self.check()
        self.assertEqual(report.skill_mismatches, ["hamel-source-check：仅有 .claude/skills"])
        self.write(".cursor/skills/other/SKILL.md", "乙\n")
        (self.root / ".claude/skills/hamel-source-check/SKILL.md").unlink()
        self.commit("cursor only")
        report = self.check()
        self.assertEqual(report.skill_mismatches, ["other：仅有 .cursor/skills"])

    def test_missing_metrics_table(self):
        self.branch()
        (self.root / "docs" / "STATUS.md").unlink()
        self.commit("drop status")
        report = self.check()
        self.assertEqual(signal_codes(report), ["metrics_missing"])

    def test_cli_reports_stage_suffix_and_unregistered_field(self):
        self.branch()
        self.write("model_stage3.py", "x = 1\n")
        self.write("hamel/config.py", 'PHYSICAL = {"invented_knob": 1}\n')
        self.commit("violations")
        result = self.run_cli("--root", str(self.root), "--base", "main", "--today", "2026-10-13")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("model_stage3.py", result.stdout)
        self.assertIn("invented_knob", result.stdout)
        self.assertIn("须总司令审查", result.stdout)
        payload = json.loads(self.run_cli(
            "--root", str(self.root), "--base", "main", "--today", "2026-10-13", "--json"
        ).stdout)
        self.assertFalse(payload["ok"])
        self.assertEqual(
            [item["code"] for item in payload["signals"]],
            ["config_without_adr", "stage_suffix"],
        )

    def test_cli_clean_exit_and_usage_errors(self):
        result = self.run_cli("--root", str(self.root), "--base", "main", "--today", "2026-10-13")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("未触发偏离信号", result.stdout)
        missing = self.run_cli("--root", str(self.root / "missing"))
        self.assertEqual(missing.returncode, 2)
        self.assertIn("找不到目录", missing.stderr)
        bad_day = self.run_cli("--root", str(self.root), "--today", "yesterday")
        self.assertEqual(bad_day.returncode, 2)
        self.assertIn("无法解析日期", bad_day.stderr)

    def run_cli(self, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(ROOT / "tools" / "drift_check.py"), *args],
            capture_output=True,
            text=True,
            check=False,
        )


if __name__ == "__main__":
    unittest.main()
