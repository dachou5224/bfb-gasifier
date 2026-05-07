#!/usr/bin/env python3
"""流程图一致性只读审计：校验 ``describe_hamel_flowchart_mapping`` 中路径可导入。

用法：在项目根目录执行 ``python3 scripts/audit_flowchart_mapping.py``。
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

# 保证可导入 src.*
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def _resolve_dotted(qualname: str) -> object:
    """解析 ``module.submodule.Class.attr`` 式限定名。"""
    parts = qualname.split(".")
    if len(parts) < 2:
        raise ImportError(f"invalid qualname: {qualname!r}")
    for k in range(len(parts), 0, -1):
        modname = ".".join(parts[:k])
        try:
            mod = importlib.import_module(modname)
        except ImportError:
            continue
        obj: object = mod
        for name in parts[k:]:
            obj = getattr(obj, name)
        return obj
    raise ImportError(f"could not import {qualname!r}")


def main() -> int:
    from src.workflow.simulation_runner import (
        describe_hamel_flowchart_mapping,
        hamel_global_nr_pipeline_order,
    )

    mapping = describe_hamel_flowchart_mapping()
    errors: list[str] = []
    for label, dotted in mapping.items():
        try:
            _resolve_dotted(dotted)
        except Exception as exc:  # noqa: BLE001 — 审计脚本需汇总全部失败
            errors.append(f"{label}: {dotted} -> {exc!r}")

    for dotted in hamel_global_nr_pipeline_order():
        try:
            _resolve_dotted(dotted)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"pipeline: {dotted} -> {exc!r}")

    mmd = _ROOT / "docs" / "gasifier_model_flowchart_trilingual.mmd"
    if not mmd.is_file():
        errors.append(f"missing mmd: {mmd}")
    else:
        text = mmd.read_text(encoding="utf-8")
        for token in ("PRECALC", "Newton-Raphson", "Check1", "Check2", "Vorabrechnung"):
            if token not in text:
                errors.append(f"mmd missing keyword: {token}")

    if errors:
        print("AUDIT FAILED", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1
    print(f"OK: {len(mapping)} mapping entries + {len(hamel_global_nr_pipeline_order())} pipeline hops + mmd keywords")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
