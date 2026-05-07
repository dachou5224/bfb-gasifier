"""审计当前反应核中的活跃反应网络。

目标：
1. 列出从 kinetics 模块导入的 `rate_R*` 函数；
2. 找出 `Cell.calc_reactions()` 中真正调用的反应；
3. 标出“已导入但未激活”的反应；
4. 按 thesis 主骨架与 implementation 扩展路径分别列出。

用法：
    python scripts/audit_reaction_activation.py
"""

from __future__ import annotations

import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
KINETICS_FILE = REPO_ROOT / "src" / "core" / "cell_kinetics.py"

THESIS_CORE_IMPLEMENTATION_RATES = {
    "rate_R1",
    "rate_R2",
    "rate_R3",
    "rate_R4_effective",
    "rate_R5_bubble",
    "rate_R5_suspension",
    "rate_R6",
    "rate_R7",
    "rate_R8",
    "rate_R10",
    "rate_R11_bubble",
    "rate_R11_suspension",
    "rate_R12",
}

IMPLEMENTATION_EXTENSION_RATES = {
    "rate_R9",  # H2S oxidation placeholder; not thesis R1-R11 main skeleton
}


class ReactionActivationAudit(ast.NodeVisitor):
    def __init__(self) -> None:
        self.imported_rates: set[str] = set()
        self.kernel_called_rates: set[str] = set()
        self._inside_kernel = False

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        module = node.module or ""
        if module.startswith("src.kinetics"):
            for alias in node.names:
                if alias.name.startswith("rate_R"):
                    self.imported_rates.add(alias.asname or alias.name)
        self.generic_visit(node)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        prev = self._inside_kernel
        if node.name == "build_reaction_sources":
            self._inside_kernel = True
            self.generic_visit(node)
            self._inside_kernel = prev
            return
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if self._inside_kernel and isinstance(node.func, ast.Name):
            if node.func.id.startswith("rate_R"):
                self.kernel_called_rates.add(node.func.id)
        self.generic_visit(node)


def main() -> int:
    source = KINETICS_FILE.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(KINETICS_FILE))
    audit = ReactionActivationAudit()
    audit.visit(tree)

    imported = sorted(audit.imported_rates)
    called = sorted(audit.kernel_called_rates)
    imported_not_used = sorted(audit.imported_rates - audit.kernel_called_rates)
    active_extensions = sorted(audit.kernel_called_rates & IMPLEMENTATION_EXTENSION_RATES)
    unexpected_active = sorted(
        audit.kernel_called_rates
        - THESIS_CORE_IMPLEMENTATION_RATES
        - IMPLEMENTATION_EXTENSION_RATES
    )
    missing_thesis_core = sorted(THESIS_CORE_IMPLEMENTATION_RATES - audit.kernel_called_rates)

    print("=" * 72)
    print("Reaction activation audit: src/core/cell_kinetics.py::build_reaction_sources")
    print("=" * 72)
    print(f"Kinetics kernel file: {KINETICS_FILE}")
    print()

    print("Imported rate functions:")
    for name in imported:
        print(f"  - {name}")
    print()

    print("Active rate calls in build_reaction_sources():")
    for name in called:
        print(f"  - {name}")
    print()

    print("Imported but not used in build_reaction_sources():")
    if imported_not_used:
        for name in imported_not_used:
            print(f"  - {name}")
    else:
        print("  (none)")
    print()

    print("Active implementation extensions outside thesis R1-R11 core:")
    if active_extensions:
        for name in active_extensions:
            print(f"  - {name}")
    else:
        print("  (none)")
    print()

    print("Unexpected active rate calls:")
    if unexpected_active:
        for name in unexpected_active:
            print(f"  - {name}")
    else:
        print("  (none)")
    print()

    print("Missing thesis-core implementation rates in build_reaction_sources():")
    if missing_thesis_core:
        for name in missing_thesis_core:
            print(f"  - {name}")
    else:
        print("  (none)")
    print()

    print("Thesis mapping note:")
    print("  - rate_R3  -> thesis R4 (hydrogasification)")
    print("  - rate_R4_effective -> thesis R3 (Boudouard)")
    print("  - rate_R6  -> thesis R7 (CH4 oxidation)")
    print("  - rate_R7  -> thesis R9 (methane reforming)")
    print("  - rate_R12 -> thesis R6 (H2 oxidation)")
    print("  - rate_R9  -> implementation sulfur placeholder")
    print()

    if imported_not_used or unexpected_active or missing_thesis_core:
        print("Result: ATTENTION NEEDED")
        return 1

    print("Result: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
