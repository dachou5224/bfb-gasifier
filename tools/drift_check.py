#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""章程 §9 偏离检查。

对照基准分支（默认 origin/main，否则 main）检查当前工作树。净增行数、新文件名和新增配置字段以该基准与 HEAD 的共同祖先为起点，避免 main 在分叉后前进把本卡的新增抵扣掉。任一信号触发则退出码为 1。
本脚本不修改仓库，也不自动 fetch，也不检查文献引用。

信号与上限：

- 工作树不干净。
- ``hamel/`` 下非测试 Python 相对共同祖先的净增行数超过 800（CHARTER §9）。
  未提交改动计入，避免超限代码只留在工作树里。目标所有者在派发前书面批准、
  并且批准记在 ``docs/STATUS.md`` 的任务卡除外：同一行要有任务号、
  「净增代码上限 N 行」和「派发前批准」，当前分支名含该任务号，净增不超过 N。
- ``hamel/`` 非测试 Python 的当前总行数超过 11000（CHARTER v0.14 §6.2）。
- 新增配置字段没有出处。字段只从 ``hamel/**/config.py`` 读取：模块级字面量字典
  ``PHYSICAL`` / ``PHYSICAL_CONFIG``、``SOLVER`` / ``SOLVER_CONFIG``，以及类
  ``PhysicalConfig``、``SolverConfig`` 的公开字段。求解器字段必须在正式 ADR
  （``docs/adr/NNNN-*.md``，不含 ``CANDIDATES.md``）里按完整标识符出现。物理字段
  同样可以靠 ADR；也可以靠 ``data/cases/*.json`` 里的同名键豁免，但该键的
  ``status`` 须为 ``source`` 或 ``derived``，``source`` 含 ``S. <数字>``，
  ``derived`` 还要有非空 ``formula``（CHARTER v0.6 §9）。
- 相对共同祖先新出现的路径组件匹配 ``_v<数字>``、``_r<数字>_r<数字>``、``_stage<数字>``。
- ``docs/STATUS.md`` 工况指标表的最近更新日期距检查日已满 14 天。
- ``.claude/skills/<名称>/SKILL.md`` 与 ``.cursor/skills/<名称>/SKILL.md`` 内容不一致，
  或只存在于其中一边。两边都没有则不检查。

退出码：0 未触发，1 触发，2 无法检查（不是 git 工作树、基准不存在、日期无效）。
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path


# CHARTER §9：单张任务卡净增代码超过 800 行。派发前批准的上限见 _approved_net_limit。
NET_CODE_LIMIT = 800
_NET_APPROVAL = re.compile(
    r"(T\d+[a-z]?)\D{0,40}净增代码上限\s*(\d+)\s*行[^|\n]{0,80}派发前批准"
)
# CHARTER v0.14 §6.2：第 1 关通过前，hamel/ 代码总量上限 11000 行（不含测试）。
HAMEL_LINE_LIMIT = 11000
# CHARTER §6.2：物理参数配置不超过 30 项。求解器参数单独成组，不计入此项。
PHYSICAL_FIELD_LIMIT = 30
# CHARTER §9：连续 2 周没有任何工况指标变化。
STALE_DAYS = 14

# 阶段后缀写在文件名或目录名里，不要求位于字符串末尾。
_STAGE_SUFFIX = re.compile(r"_v\d+|_r\d+_r\d+|_stage\d+")
_DATE = re.compile(r"(\d{4})-(\d{2})-(\d{2})")
# T7：印刷页码写成 “S. 143”；允许 S. 与数字之间有空白。
_PRINTED_PAGE = re.compile(r"S\.\s*\d+")
_GROUP_NAMES = {
    "PHYSICAL": "physical",
    "PHYSICAL_CONFIG": "physical",
    "SOLVER": "solver",
    "SOLVER_CONFIG": "solver",
}
_CLASS_GROUPS = {"PhysicalConfig": "physical", "SolverConfig": "solver"}
_CONFIG_SIGNAL_CODES = ("config_without_adr", "physical_config_limit", "config_unreadable")


class DriftError(Exception):
    """仓库状态不足以完成检查。"""


class ConfigError(Exception):
    """config.py 里有配置绑定，但读不出字段名。"""


@dataclass(frozen=True)
class Signal:
    code: str
    message: str


@dataclass
class Report:
    root: str
    base: str
    compared: str
    base_notes: list[str]
    today: str
    dirty: bool
    net_lines: int
    hamel_lines: int
    physical_base: int
    physical_head: int
    solver_base: int
    solver_head: int
    added_physical: list[str]
    removed_physical: list[str]
    added_solver: list[str]
    removed_solver: list[str]
    unregistered_fields: list[str]
    stage_files: list[str]
    skill_mismatches: list[str]
    config_errors: list[str]
    metrics_date: str | None
    metrics_age_days: int | None
    signals: list[Signal]
    approved_net_limit: int | None = None

    def as_dict(self) -> dict:
        return {
            "ok": not self.signals,
            "root": self.root,
            "base": self.base,
            "compared": self.compared,
            "base_notes": self.base_notes,
            "today": self.today,
            "dirty": self.dirty,
            "net_lines": self.net_lines,
            "hamel_lines": self.hamel_lines,
            "physical_fields": {"base": self.physical_base, "head": self.physical_head},
            "solver_fields": {"base": self.solver_base, "head": self.solver_head},
            "added_physical": self.added_physical,
            "removed_physical": self.removed_physical,
            "added_solver": self.added_solver,
            "removed_solver": self.removed_solver,
            "unregistered_fields": self.unregistered_fields,
            "stage_files": self.stage_files,
            "skill_mismatches": self.skill_mismatches,
            "config_errors": self.config_errors,
            "metrics_date": self.metrics_date,
            "metrics_age_days": self.metrics_age_days,
            "signals": [{"code": signal.code, "message": signal.message} for signal in self.signals],
            "limits": {
                "net_lines": self.approved_net_limit or NET_CODE_LIMIT,
                "hamel_lines": HAMEL_LINE_LIMIT,
                "physical_fields": PHYSICAL_FIELD_LIMIT,
                "stale_days": STALE_DAYS,
            },
        }


def resolve_root(explicit: Path | None) -> Path:
    start = explicit.resolve() if explicit is not None else Path.cwd()
    if explicit is not None and not start.exists():
        raise DriftError(f"找不到目录：{start}")
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=start,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise DriftError("找不到 git") from exc
    if proc.returncode != 0:
        raise DriftError(f"不是 git 工作树：{start}")
    return Path(proc.stdout.strip()).resolve()


def _ref_exists(root: Path, ref: str) -> bool:
    proc = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", ref],
        cwd=root,
        capture_output=True,
        check=False,
    )
    return proc.returncode == 0


def default_base(root: Path) -> str:
    if _ref_exists(root, "origin/main"):
        return "origin/main"
    if _ref_exists(root, "main"):
        return "main"
    raise DriftError("找不到基准分支 origin/main 或 main，请用 --base 指定")


def _comparison_commit(root: Path, base: str) -> str:
    return _git(root, ["merge-base", base, "HEAD"]).decode().strip()


def _rev_list_count(root: Path, revision_range: str) -> int:
    try:
        blob = _git(root, ["rev-list", "--count", revision_range])
    except DriftError:
        return 0
    return int(blob.decode().strip() or "0")


def _reflog_time(root: Path, ref: str) -> str | None:
    proc = subprocess.run(
        ["git", "reflog", "show", "-1", "--date=iso", "--format=%gd", ref],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    text = proc.stdout.strip()
    if proc.returncode != 0 or "@{" not in text or not text.endswith("}"):
        return None
    return text.split("@{", 1)[1][:-1]


def _base_notes(root: Path) -> list[str]:
    notes: list[str] = []
    if not _ref_exists(root, "origin/main"):
        return notes
    fetched = _reflog_time(root, "origin/main")
    if fetched:
        notes.append(f"origin/main 上次更新于 {fetched}。脚本不会自动 fetch")
    else:
        notes.append("origin/main 没有 reflog，无法判断上次 fetch。脚本不会自动 fetch")
    if _ref_exists(root, "main"):
        behind = _rev_list_count(root, "main..origin/main")
        ahead = _rev_list_count(root, "origin/main..main")
        if behind:
            notes.append(f"本地 main 落后 origin/main {behind} 个提交")
        if ahead:
            notes.append(f"本地 main 比 origin/main 超前 {ahead} 个提交")
    return notes


def _git(root: Path, args: list[str]) -> bytes:
    try:
        proc = subprocess.run(
            ["git", "-c", "core.quotepath=false", *args],
            cwd=root,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise DriftError("找不到 git") from exc
    if proc.returncode != 0:
        detail = proc.stderr.decode("utf-8", "replace").strip()
        raise DriftError(detail or f"git {' '.join(args)} 失败")
    return proc.stdout


def _nul_paths(blob: bytes) -> list[str]:
    return [part.decode("utf-8") for part in blob.split(b"\0") if part]


def _ensure_base(root: Path, base: str) -> None:
    proc = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", base],
        cwd=root,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise DriftError(f"找不到基准 {base}")


def _is_hamel_code(path: str) -> bool:
    parts = Path(path).parts
    if len(parts) < 2 or parts[0] != "hamel" or not path.endswith(".py"):
        return False
    if "__pycache__" in parts or "tests" in parts:
        return False
    return not Path(path).name.startswith("test_")


def _is_config_path(path: str) -> bool:
    parts = Path(path).parts
    return len(parts) >= 2 and parts[0] == "hamel" and parts[-1] == "config.py" and "__pycache__" not in parts


def _count_lines(blob: bytes) -> int:
    if not blob:
        return 0
    return len(blob.splitlines())


def _worktree_code_lines(root: Path) -> dict[str, int]:
    hamel = root / "hamel"
    if not hamel.is_dir():
        return {}
    counts: dict[str, int] = {}
    for path in hamel.rglob("*.py"):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if _is_hamel_code(relative):
            counts[relative] = _count_lines(path.read_bytes())
    return counts


def _base_hamel_paths(root: Path, base: str) -> list[str]:
    blob = _git(root, ["ls-tree", "-r", "-z", "--name-only", base, "--", "hamel"])
    return _nul_paths(blob)


def _base_line_counts(root: Path, base: str, paths: list[str]) -> dict[str, int]:
    if not paths:
        return {}
    specs = [f"{base}:{path}" for path in paths]
    proc = subprocess.run(
        ["git", "cat-file", "--batch"],
        cwd=root,
        input=b"".join(spec.encode("utf-8") + b"\n" for spec in specs),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        detail = proc.stderr.decode("utf-8", "replace").strip()
        raise DriftError(detail or "git cat-file 失败")
    counts: dict[str, int] = {}
    cursor = 0
    data = proc.stdout
    for spec, path in zip(specs, paths):
        newline = data.find(b"\n", cursor)
        if newline < 0:
            raise DriftError("git cat-file 输出截断")
        header = data[cursor:newline]
        cursor = newline + 1
        if header.endswith(b" missing"):
            counts[path] = 0
            continue
        parts = header.split()
        if len(parts) < 3 or parts[1] != b"blob":
            raise DriftError(f"无法解析 git cat-file 头：{header.decode('utf-8', 'replace')}")
        size = int(parts[2])
        blob = data[cursor:cursor + size]
        cursor += size
        if cursor < len(data) and data[cursor:cursor + 1] == b"\n":
            cursor += 1
        counts[path] = _count_lines(blob)
    return counts


def _line_delta(root: Path, base: str) -> tuple[int, int]:
    base_paths = [path for path in _base_hamel_paths(root, base) if _is_hamel_code(path)]
    base_counts = _base_line_counts(root, base, base_paths)
    work_counts = _worktree_code_lines(root)
    net = sum(work_counts.get(path, 0) - base_counts.get(path, 0) for path in set(base_counts) | set(work_counts))
    return net, sum(work_counts.values())


def _const_str(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    str_node = getattr(ast, "Str", None)
    if str_node is not None and isinstance(node, str_node):
        return node.s
    return None


def _dict_keys(value: ast.AST) -> set[str] | None:
    if isinstance(value, ast.Dict):
        keys: set[str] = set()
        for key in value.keys:
            if key is None:
                return None
            text = _const_str(key)
            if text is None:
                return None
            if not text.startswith("_"):
                keys.add(text)
        return keys
    if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict":
        if value.args or any(keyword.arg is None for keyword in value.keywords):
            return None
        return {keyword.arg for keyword in value.keywords if keyword.arg and not keyword.arg.startswith("_")}
    return None


def _target_names(target: ast.AST) -> list[str]:
    if isinstance(target, ast.Name):
        return [target.id]
    if isinstance(target, (ast.Tuple, ast.List)):
        names: list[str] = []
        for element in target.elts:
            names.extend(_target_names(element))
        return names
    return []


def _class_fields(stmt: ast.ClassDef) -> set[str]:
    fields: set[str] = set()
    for item in stmt.body:
        names: list[str] = []
        if isinstance(item, ast.AnnAssign):
            names = _target_names(item.target)
        elif isinstance(item, ast.Assign):
            for target in item.targets:
                names.extend(_target_names(target))
        fields.update(name for name in names if not name.startswith("_"))
    return fields


def _config_root(node: ast.AST) -> str | None:
    current = node
    while isinstance(current, (ast.Attribute, ast.Subscript)):
        current = current.value
    if isinstance(current, ast.Name) and current.id in _GROUP_NAMES:
        return current.id
    return None


def _config_mutation(tree: ast.AST) -> str | None:
    """字面量赋值之后的 update、下标写入等会漏掉字段，直接判为不可读。"""
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            name = _config_root(node.func.value)
            if name:
                return name
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                name = _config_root(target) if isinstance(target, ast.Subscript) else None
                if name:
                    return name
        elif isinstance(node, ast.AugAssign):
            if isinstance(node.target, ast.Name) and node.target.id in _GROUP_NAMES:
                return node.target.id
            name = _config_root(node.target) if isinstance(node.target, ast.Subscript) else None
            if name:
                return name
        elif isinstance(node, ast.Delete):
            for target in node.targets:
                name = _config_root(target) if isinstance(target, ast.Subscript) else None
                if name:
                    return name
    return None


def extract_config_fields(source: str) -> tuple[set[str], set[str]]:
    """返回 (物理字段, 求解器字段)。非字面量绑定或事后修改视为不可读。"""
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        raise ConfigError(f"语法错误：{exc.msg}") from exc
    mutated = _config_mutation(tree)
    if mutated:
        raise ConfigError(f"{mutated} 在字面量之外被修改，无法静态读取全部字段")
    dicts: dict[str, set[str]] = {}
    class_fields = {"physical": set(), "solver": set()}
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign):
            targets, value = stmt.targets, stmt.value
        elif isinstance(stmt, ast.AnnAssign) and stmt.value is not None:
            targets, value = [stmt.target], stmt.value
        elif isinstance(stmt, ast.ClassDef) and stmt.name in _CLASS_GROUPS:
            class_fields[_CLASS_GROUPS[stmt.name]] |= _class_fields(stmt)
            continue
        else:
            continue
        names = [name for target in targets for name in _target_names(target)]
        interesting = [name for name in names if name in _GROUP_NAMES]
        if not interesting:
            continue
        if len(names) != 1:
            raise ConfigError(f"{interesting[0]} 必须单独赋值为字面量字典")
        keys = _dict_keys(value)
        if keys is None:
            raise ConfigError(f"{interesting[0]} 必须是字面量字典")
        dicts[interesting[0]] = keys
    physical = set(class_fields["physical"])
    solver = set(class_fields["solver"])
    for name, keys in dicts.items():
        if _GROUP_NAMES[name] == "physical":
            physical |= keys
        else:
            solver |= keys
    return physical, solver


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise ConfigError("不是 UTF-8 文本") from exc
    except OSError as exc:
        raise ConfigError(str(exc)) from exc


def _show_text(root: Path, spec: str) -> str | None:
    proc = subprocess.run(
        ["git", "-c", "core.quotepath=false", "show", spec],
        cwd=root,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    try:
        return proc.stdout.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ConfigError("不是 UTF-8 文本") from exc


def _work_config_paths(root: Path) -> list[str]:
    hamel = root / "hamel"
    if not hamel.is_dir():
        return []
    found = []
    for path in hamel.rglob("config.py"):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if _is_config_path(relative):
            found.append(relative)
    return sorted(found)


def _load_fields(items: list[tuple[str, str]]) -> tuple[set[str], set[str], list[str]]:
    physical: set[str] = set()
    solver: set[str] = set()
    errors: list[str] = []
    for path, source in items:
        try:
            found_physical, found_solver = extract_config_fields(source)
        except ConfigError as exc:
            errors.append(f"{path}：{exc}")
            continue
        physical |= found_physical
        solver |= found_solver
    return physical, solver, errors


def _config_sides(root: Path, base: str) -> tuple[set[str], set[str], set[str], set[str], list[str]]:
    errors: list[str] = []
    base_items: list[tuple[str, str]] = []
    for path in _base_hamel_paths(root, base):
        if not _is_config_path(path):
            continue
        try:
            source = _show_text(root, f"{base}:{path}")
        except ConfigError as exc:
            errors.append(f"{path}：{exc}")
            continue
        if source is None:
            errors.append(f"{path}：基准中无法读取")
            continue
        base_items.append((path, source))
    work_items: list[tuple[str, str]] = []
    for path in _work_config_paths(root):
        file = root / path
        try:
            work_items.append((path, _read_text(file)))
        except ConfigError as exc:
            errors.append(f"{path}：{exc}")
    base_physical, base_solver, base_errors = _load_fields(base_items)
    work_physical, work_solver, work_errors = _load_fields(work_items)
    return base_physical, work_physical, base_solver, work_solver, errors + base_errors + work_errors


def case_field_qualifies(obj: dict) -> bool:
    """CHARTER v0.6 §9：带来源页码的物理输入可以不另写 ADR。"""
    status = obj.get("status")
    if status not in {"source", "derived"}:
        return False
    source = obj.get("source")
    if not isinstance(source, str) or _PRINTED_PAGE.search(source) is None:
        return False
    if status == "derived":
        formula = obj.get("formula")
        if not isinstance(formula, str) or not formula.strip():
            return False
    return True


def _walk_case_fields(node: object, found: dict[str, list[dict]]) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            if isinstance(value, dict) and "status" in value:
                found.setdefault(str(key), []).append(value)
            _walk_case_fields(value, found)
    elif isinstance(node, list):
        for item in node:
            _walk_case_fields(item, found)


def case_exemptions(root: Path) -> tuple[set[str], list[str]]:
    """返回可豁免的物理字段名，以及读不了的工况文件。同名键必须全部合格。"""
    folder = root / "data" / "cases"
    if not folder.is_dir():
        return set(), []
    found: dict[str, list[dict]] = {}
    errors: list[str] = []
    for path in sorted(folder.glob("*.json")):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            errors.append(f"{relative}：无法读取工况字段（{exc.__class__.__name__}）")
            continue
        _walk_case_fields(payload, found)
    exempt = {
        name
        for name, objects in found.items()
        if objects and all(case_field_qualifies(obj) for obj in objects)
    }
    return exempt, errors


def field_in_adr(name: str, corpus: str) -> bool:
    """字段名按标识符整体匹配，避免 damping 命中 newton_damping。"""
    return re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])", corpus) is not None


def _adr_corpus(root: Path) -> str:
    folder = root / "docs" / "adr"
    if not folder.is_dir():
        return ""
    chunks: list[str] = []
    for path in sorted(folder.glob("[0-9][0-9][0-9][0-9]-*.md")):
        chunks.append(path.name)
        chunks.append(path.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(chunks)


def skill_mismatches(root: Path) -> list[str]:
    """同名 SKILL.md 须在 .claude 与 .cursor 两边内容一致。两边都没有则不检查。"""
    roots = {
        ".claude/skills": root / ".claude" / "skills",
        ".cursor/skills": root / ".cursor" / "skills",
    }
    names: set[str] = set()
    for folder in roots.values():
        if not folder.is_dir():
            continue
        for child in folder.iterdir():
            if (child / "SKILL.md").is_file():
                names.add(child.name)
    problems: list[str] = []
    for name in sorted(names):
        left = roots[".claude/skills"] / name / "SKILL.md"
        right = roots[".cursor/skills"] / name / "SKILL.md"
        has_left = left.is_file()
        has_right = right.is_file()
        if has_left and not has_right:
            problems.append(f"{name}：仅有 .claude/skills")
        elif has_right and not has_left:
            problems.append(f"{name}：仅有 .cursor/skills")
        elif left.read_bytes() != right.read_bytes():
            problems.append(f"{name}：内容不一致")
    return problems


def has_stage_suffix(path: str) -> bool:
    return any(_STAGE_SUFFIX.search(part) for part in Path(path).parts)


def _added_paths(blob: bytes) -> list[str]:
    parts = [part for part in blob.split(b"\0") if part]
    added: list[str] = []
    index = 0
    while index < len(parts):
        status = parts[index].decode("utf-8", "replace")
        index += 1
        if not status:
            continue
        code = status[0]
        if code in {"R", "C"}:
            if index + 1 >= len(parts):
                raise DriftError("git diff --name-status 输出不完整")
            added.append(parts[index + 1].decode("utf-8"))
            index += 2
            continue
        if index >= len(parts):
            raise DriftError("git diff --name-status 输出不完整")
        path = parts[index].decode("utf-8")
        index += 1
        if code == "A":
            added.append(path)
    return added


def _new_paths(root: Path, base: str) -> list[str]:
    diff = _git(root, ["diff", "--name-status", "-z", "--find-renames", base])
    found = set(_added_paths(diff))
    found.update(_nul_paths(_git(root, ["ls-files", "-z", "--others", "--exclude-standard"])))
    return sorted(found)


def parse_metrics_date(text: str) -> date | None:
    """读取「工况指标」表里最新的更新日期。表外日期忽略。"""
    lines = text.splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.startswith("## 工况指标"):
            start = index + 1
            break
    if start is None:
        return None
    column = None
    found: list[date] = []
    for line in lines[start:]:
        if line.startswith("## "):
            break
        stripped = line.strip()
        if not stripped.startswith("|"):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if column is None:
            if "更新日期" not in cells:
                continue
            column = cells.index("更新日期")
            continue
        if column >= len(cells):
            continue
        match = _DATE.fullmatch(cells[column])
        if not match:
            continue
        try:
            found.append(date(int(match.group(1)), int(match.group(2)), int(match.group(3))))
        except ValueError:
            continue
    if not found:
        return None
    return max(found)


def _metrics(root: Path, today: date) -> tuple[date | None, int | None]:
    path = root / "docs" / "STATUS.md"
    if not path.is_file():
        return None, None
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None, None
    found = parse_metrics_date(text)
    if found is None:
        return None, None
    return found, (today - found).days


def _field_delta(base: int, head: int, added: list[str], removed: list[str]) -> str:
    text = f"{base} → {head}"
    notes = []
    if added:
        notes.append("新增 " + "、".join(added))
    if removed:
        notes.append("移除 " + "、".join(removed))
    if notes:
        text += "（" + "；".join(notes) + "）"
    return text


def _config_fact(report: Report) -> str:
    physical = _field_delta(
        report.physical_base, report.physical_head, report.added_physical, report.removed_physical
    )
    solver = _field_delta(report.solver_base, report.solver_head, report.added_solver, report.removed_solver)
    return f"配置字段：物理参数 {physical}；求解器 {solver}"


def _worktree_text(report: Report) -> str:
    return "工作树不干净" if report.dirty else "工作树干净"


def _net_text(report: Report) -> str:
    limit = report.approved_net_limit or NET_CODE_LIMIT
    if report.net_lines > limit:
        if report.approved_net_limit:
            return f"相对共同祖先净增代码 {report.net_lines} 行，超过批准上限 {limit} 行"
        return f"相对共同祖先净增代码 {report.net_lines} 行，超过 {NET_CODE_LIMIT} 行"
    if report.approved_net_limit and report.approved_net_limit > NET_CODE_LIMIT:
        return (
            f"相对共同祖先净增代码 {report.net_lines} 行"
            f"（章程 §9 批准上限 {report.approved_net_limit}，默认上限 {NET_CODE_LIMIT}）"
        )
    return f"相对共同祖先净增代码 {report.net_lines} 行（上限 {NET_CODE_LIMIT}）"


def _task_in_branch(task: str, branch: str) -> bool:
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(task)}(?![A-Za-z0-9])", branch) is not None


def _approved_net_limit(root: Path) -> int | None:
    """STATUS 里记了派发前批准、且当前分支就是该任务卡时，返回批准上限。"""
    branch = _git(root, ["rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
    if branch in ("", "HEAD"):
        return None
    path = root / "docs" / "STATUS.md"
    if not path.is_file():
        return None
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None
    approved: int | None = None
    for match in _NET_APPROVAL.finditer(text):
        task, limit = match.group(1), int(match.group(2))
        if _task_in_branch(task, branch):
            approved = limit if approved is None else max(approved, limit)
    return approved


def _hamel_text(report: Report) -> str:
    if report.hamel_lines > HAMEL_LINE_LIMIT:
        return f"hamel/ 代码 {report.hamel_lines} 行，超过 {HAMEL_LINE_LIMIT} 行（不含测试）"
    return f"hamel/ 代码 {report.hamel_lines} 行（上限 {HAMEL_LINE_LIMIT}，不含测试）"


def _stage_text(report: Report) -> str:
    if not report.stage_files:
        return "新文件名无阶段后缀"
    return "新文件名呈阶段后缀：" + "、".join(report.stage_files)


def _skill_text(report: Report) -> str:
    if not report.skill_mismatches:
        return "skill 副本一致，或两边都没有"
    return "skill 副本不一致：" + "、".join(report.skill_mismatches)


def _metrics_text(report: Report) -> str:
    if report.metrics_date is None:
        return "工况指标表缺少可读取的更新日期"
    if report.metrics_age_days is not None and report.metrics_age_days < 0:
        return f"工况指标最近变化 {report.metrics_date}，晚于检查日 {report.today}"
    if report.metrics_age_days is not None and report.metrics_age_days >= STALE_DAYS:
        return (
            f"工况指标最近变化 {report.metrics_date}，"
            f"距检查日已 {report.metrics_age_days} 天，连续 2 周没有变化"
        )
    return f"工况指标最近变化 {report.metrics_date}，距检查日 {report.metrics_age_days} 天"


def _mark(text: str, triggered: bool) -> str:
    return f"{'[触发]' if triggered else '[通过]'} {text}"


def format_text(report: Report) -> str:
    triggered = {signal.code for signal in report.signals}
    config_hits = [signal for signal in report.signals if signal.code in _CONFIG_SIGNAL_CODES]
    if config_hits:
        config_lines = [f"[触发] {signal.message}" for signal in config_hits]
        config_lines.append(f"       {_config_fact(report)}")
    else:
        config_lines = [_mark(_config_fact(report), False)]
    metric_triggered = "metrics_stale" in triggered or "metrics_missing" in triggered
    lines = [
        "偏离检查（CHARTER §9）",
        f"根目录：{report.root}",
        f"基准：{report.base}",
        f"比较提交：{report.compared[:12]}",
        *[f"基准说明：{note}" for note in report.base_notes],
        f"检查日：{report.today}",
        "",
        _mark(_worktree_text(report), "dirty" in triggered),
        _mark(_net_text(report), "net_lines" in triggered),
        _mark(_hamel_text(report), "hamel_lines" in triggered),
        *config_lines,
        _mark(_stage_text(report), "stage_suffix" in triggered),
        _mark(_skill_text(report), "skill_mismatch" in triggered),
        _mark(_metrics_text(report), metric_triggered),
        "",
        (
            f"结论：触发 {len(report.signals)} 项，须总司令审查。"
            if report.signals
            else "结论：未触发偏离信号。"
        ),
    ]
    return "\n".join(lines) + "\n"


def check(root: Path, base: str, today: date) -> Report:
    root = resolve_root(root)
    _ensure_base(root, base)
    compared = _comparison_commit(root, base)
    notes = _base_notes(root)
    dirty = bool(
        _git(root, ["-c", "status.showUntrackedFiles=all", "status", "--porcelain=v1", "-z"]).strip(b"\0")
    )
    net_lines, hamel_lines = _line_delta(root, compared)
    base_physical, work_physical, base_solver, work_solver, config_errors = _config_sides(root, compared)
    added_physical = sorted(work_physical - base_physical)
    removed_physical = sorted(base_physical - work_physical)
    added_solver = sorted(work_solver - base_solver)
    removed_solver = sorted(base_solver - work_solver)
    corpus = _adr_corpus(root)
    exempt, case_errors = case_exemptions(root)
    config_errors = [*config_errors, *case_errors]
    unregistered = sorted(
        name
        for name in added_physical
        if not field_in_adr(name, corpus) and name not in exempt
    ) + sorted(name for name in added_solver if not field_in_adr(name, corpus))
    unregistered = sorted(set(unregistered))
    stage_files = sorted(path for path in _new_paths(root, compared) if has_stage_suffix(path))
    mismatched_skills = skill_mismatches(root)
    metrics_date, metrics_age = _metrics(root, today)

    signals: list[Signal] = []
    if dirty:
        signals.append(Signal("dirty", "工作树不干净"))
    approved_net = _approved_net_limit(root)
    net_limit = approved_net if approved_net is not None and approved_net > NET_CODE_LIMIT else NET_CODE_LIMIT
    if net_lines > net_limit:
        if approved_net is not None and approved_net > NET_CODE_LIMIT:
            signals.append(Signal("net_lines", f"相对共同祖先净增代码 {net_lines} 行，超过批准上限 {net_limit} 行"))
        else:
            signals.append(Signal("net_lines", f"相对共同祖先净增代码 {net_lines} 行，超过 {NET_CODE_LIMIT} 行"))
    if hamel_lines > HAMEL_LINE_LIMIT:
        signals.append(
            Signal("hamel_lines", f"hamel/ 代码 {hamel_lines} 行，超过 {HAMEL_LINE_LIMIT} 行（不含测试）")
        )
    if unregistered:
        signals.append(Signal("config_without_adr", "新增配置项没有对应 ADR：" + "、".join(unregistered)))
    if len(work_physical) > PHYSICAL_FIELD_LIMIT:
        signals.append(
            Signal(
                "physical_config_limit",
                f"物理参数配置 {len(work_physical)} 项，超过 {PHYSICAL_FIELD_LIMIT} 项",
            )
        )
    if config_errors:
        signals.append(Signal("config_unreadable", "无法读取配置字段：" + "；".join(config_errors)))
    if stage_files:
        signals.append(Signal("stage_suffix", "新文件名呈阶段后缀：" + "、".join(stage_files)))
    if mismatched_skills:
        signals.append(Signal("skill_mismatch", "skill 副本不一致：" + "、".join(mismatched_skills)))
    if metrics_date is None:
        signals.append(Signal("metrics_missing", "工况指标表缺少可读取的更新日期"))
    elif metrics_age is not None and metrics_age >= STALE_DAYS:
        signals.append(
            Signal(
                "metrics_stale",
                f"工况指标最近变化 {metrics_date.isoformat()}，距检查日已 {metrics_age} 天，连续 2 周没有变化",
            )
        )

    return Report(
        root=str(root),
        base=base,
        compared=compared,
        base_notes=notes,
        today=today.isoformat(),
        dirty=dirty,
        net_lines=net_lines,
        hamel_lines=hamel_lines,
        physical_base=len(base_physical),
        physical_head=len(work_physical),
        solver_base=len(base_solver),
        solver_head=len(work_solver),
        added_physical=added_physical,
        removed_physical=removed_physical,
        added_solver=added_solver,
        removed_solver=removed_solver,
        unregistered_fields=unregistered,
        stage_files=stage_files,
        skill_mismatches=mismatched_skills,
        config_errors=config_errors,
        metrics_date=None if metrics_date is None else metrics_date.isoformat(),
        metrics_age_days=metrics_age,
        approved_net_limit=approved_net if approved_net is not None and approved_net > NET_CODE_LIMIT else None,
        signals=signals,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="检查章程 §9 偏离信号。")
    parser.add_argument("--root", type=Path, default=None, help="仓库根目录，默认当前工作树")
    parser.add_argument("--base", default=None, help="对比基准，默认 origin/main，否则 main")
    parser.add_argument("--today", default=None, help="检查日 YYYY-MM-DD，默认今天")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出")
    args = parser.parse_args(argv)
    try:
        today = date.fromisoformat(args.today) if args.today else date.today()
    except ValueError:
        print(f"无法解析日期：{args.today}", file=sys.stderr)
        return 2
    try:
        root = resolve_root(args.root)
        base = args.base or default_base(root)
        report = check(root, base, today)
    except DriftError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if args.json:
        json.dump(report.as_dict(), sys.stdout, ensure_ascii=False, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(format_text(report))
    return 1 if report.signals else 0


if __name__ == "__main__":
    sys.exit(main())
