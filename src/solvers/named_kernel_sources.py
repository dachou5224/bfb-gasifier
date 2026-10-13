"""具名核检查点的源码指纹。

覆盖实际执行到的 P2、自由空间轨迹、热解和案例数据。
Git HEAD 和少数入口文件不能代替这份指纹。
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_SEEDS = (
    "src/solvers/wirsum_1997_named_f.py",
    "src/solvers/wirsum_1997_named_f_e8.py",
    "src/solvers/wirsum_1997_named_f_control.py",
    "src/solvers/wirsum_1997_named_f_p2.py",
    "src/solvers/wesseling_two_zone.py",
    "src/solvers/wesseling_case_contract.py",
    "src/thermal/devolatilization.py",
    "data/validation_cases.json",
)


def _local_modules(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: list[str] = []
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("src."):
            names.append(node.module)
        elif isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names if alias.name.startswith("src."))
        for name in names:
            module = name.replace(".", "/")
            for rel in (f"{module}.py", f"{module}/__init__.py"):
                if (_ROOT / rel).is_file():
                    found.append(rel)
    return found


def named_kernel_source_files() -> tuple[str, ...]:
    """入口文件及其 ``src`` 导入闭包，再加上案例数据。"""
    seen = set(_SEEDS)
    stack = [rel for rel in _SEEDS if rel.endswith(".py")]
    while stack:
        rel = stack.pop()
        for child in _local_modules(_ROOT / rel):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return tuple(sorted(seen))


def named_kernel_source_sha256() -> str:
    """按路径排序后的内容哈希。路径本身也进入摘要，避免同内容换名。"""
    digest = hashlib.sha256()
    for rel in named_kernel_source_files():
        digest.update(rel.encode())
        digest.update((_ROOT / rel).read_bytes())
    return digest.hexdigest()
