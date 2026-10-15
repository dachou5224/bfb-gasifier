#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""检查 hamel/ 里的 Hamel 引用是否对得上 refs/hamel/index.json。

扫描全部 ``Source: Hamel Eq./Tab./Bild …``。式号、表号、图号必须出现在
``equations``、``tables``、``figures`` 里。写了印刷页时，页码须与该条目的
``printed_page`` 一致。

``ambiguous_number: true`` 的重复编号（如式 3.49、3.50，表 5.1）还可以对上
``additional_locations`` 里的某一页。这种编号不写页码就报错，不默认采用主位置。
没有该标记的条目只认主位置的 ``printed_page``；续页写在 ``additional_locations``
里、但没有 ``ambiguous_number`` 的，不拿来放行页码。

``balances/``、``kinetics/``、``hydrodynamics/``、``pyrolysis/``、``freeboard/``
里的公开函数（模块级，以及公开类的公开方法；名字不以 ``_`` 开头）的 docstring
必须至少有一处上述引用。其它模块里的函数不因缺引用报错，但文件里出现的引用仍要核对。

退出码：0 无问题，1 有问题，2 无法检查。
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path


CORE_MODULES = (
    "balances",
    "kinetics",
    "hydrodynamics",
    "pyrolysis",
    "freeboard",
)
KIND_TO_SECTION = {"eq": "equations", "tab": "tables", "bild": "figures"}

# Source: Hamel Eq.5.39, S. 91 / Tab. 5.2 / Bild 5.1 / Tab.A.2.1.2
_CITATION = re.compile(
    r"Source:\s*Hamel\s+"
    r"(?P<kind>Eq|Tab|Bild)\.?\s*"
    r"(?P<label>(?:[A-Za-z]\.)?\d+(?:\.\d+)+)"
    r"(?:\s*,\s*S\.\s*(?P<page>\d+))?",
    re.IGNORECASE,
)


class CheckError(Exception):
    """索引或目录不足以完成检查。"""


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    code: str
    message: str

    def format(self) -> str:
        return f"{self.path}:{self.line}: {self.code}: {self.message}"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_index(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckError(f"无法读取索引 {path}：{exc}") from exc
    missing = [name for name in ("equations", "tables", "figures") if not isinstance(data.get(name), dict)]
    if missing:
        raise CheckError(f"索引缺少可用的 {', '.join(missing)}")
    return data


def _pages_match(cited: str, indexed: object) -> bool:
    indexed_text = str(indexed).strip()
    if cited == indexed_text:
        return True
    return cited.isdigit() and indexed_text.isdigit() and int(cited) == int(indexed_text)


def _page_list(entry: dict, *, include_additional: bool) -> list[str]:
    pages: list[str] = []
    primary = entry.get("printed_page")
    if primary is not None:
        pages.append(str(primary).strip())
    if include_additional:
        for location in entry.get("additional_locations") or []:
            printed = location.get("printed_page")
            if printed is not None:
                pages.append(str(printed).strip())
    return pages


def _format_pages(pages: list[str]) -> str:
    return "、".join(f"S. {page}" for page in pages)


def check_citation(kind: str, label: str, page: str | None, index: dict) -> str | None:
    """返回错误说明；通过时返回 None。"""
    section = KIND_TO_SECTION[kind.lower()]
    entry = index[section].get(label)
    kind_text = {"eq": "Eq.", "tab": "Tab.", "bild": "Bild "}[kind.lower()]
    shown = f"{kind_text}{label}"
    if entry is None:
        return f"{shown} 不在 index.json 的 {section} 中"
    ambiguous = bool(entry.get("ambiguous_number"))
    allowed = _page_list(entry, include_additional=ambiguous)
    if page is None:
        if ambiguous:
            return f"{shown} 是重复编号，引用必须带页码（{_format_pages(allowed)}）"
        return None
    if any(_pages_match(page, item) for item in allowed):
        return None
    if ambiguous:
        return f"{shown} 的 S. {page} 不在主位置或 additional_locations（{_format_pages(allowed)}）中"
    primary = _format_pages(allowed) or "（索引未给 printed_page）"
    return f"{shown} 的 S. {page} 与 printed_page {primary} 不一致"


def _line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _citation_findings(text: str, rel: str, index: dict) -> list[Finding]:
    findings: list[Finding] = []
    for match in _CITATION.finditer(text):
        message = check_citation(match.group("kind"), match.group("label"), match.group("page"), index)
        if message is None:
            continue
        page = match.group("page")
        code = "unknown_label" if "不在 index.json" in message else (
            "ambiguous_page_required" if page is None else "page_mismatch"
        )
        findings.append(Finding(rel, _line_number(text, match.start()), code, message))
    return findings


def _public_functions(tree: ast.AST) -> list[ast.AST]:
    found: list[ast.AST] = []
    for node in getattr(tree, "body", []):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
            found.append(node)
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and not item.name.startswith("_"):
                    found.append(item)
    return found


def _missing_citation_findings(text: str, rel: str) -> list[Finding]:
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        return [Finding(rel, exc.lineno or 1, "syntax_error", f"无法解析：{exc.msg}")]
    findings: list[Finding] = []
    for node in _public_functions(tree):
        docstring = ast.get_docstring(node, clean=False) or ""
        if _CITATION.search(docstring):
            continue
        findings.append(
            Finding(
                rel,
                node.lineno,
                "missing_citation",
                f"公开函数 {_qualified_name(tree, node)} 的 docstring 没有 Source 引用",
            )
        )
    return findings


def _qualified_name(tree: ast.AST, node: ast.AST) -> str:
    for parent in getattr(tree, "body", []):
        if isinstance(parent, ast.ClassDef) and any(item is node for item in parent.body):
            return f"{parent.name}.{node.name}"
    return node.name


def _is_core(rel: Path) -> bool:
    return len(rel.parts) >= 2 and rel.parts[0] in CORE_MODULES


def check(hamel_root: Path, index_path: Path) -> list[Finding]:
    if not hamel_root.is_dir():
        raise CheckError(f"找不到 hamel 目录：{hamel_root}")
    index = load_index(index_path)
    findings: list[Finding] = []
    for path in sorted(hamel_root.rglob("*.py")):
        rel = path.relative_to(hamel_root).as_posix()
        text = path.read_text(encoding="utf-8")
        findings.extend(_citation_findings(text, rel, index))
        if _is_core(path.relative_to(hamel_root)):
            findings.extend(_missing_citation_findings(text, rel))
    findings.sort(key=lambda item: (item.path, item.line, item.code, item.message))
    return findings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="检查 hamel/ 的 Hamel 式号、表号、图号引用。")
    parser.add_argument("--hamel", type=Path, default=None, help="hamel 包目录，默认仓库内的 hamel/")
    parser.add_argument("--index", type=Path, default=None, help="index.json 路径")
    args = parser.parse_args(argv)
    root = repo_root()
    hamel_root = args.hamel or (root / "hamel")
    index_path = args.index or (root / "refs" / "hamel" / "index.json")
    try:
        findings = check(hamel_root, index_path)
    except CheckError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    for finding in findings:
        print(finding.format())
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
