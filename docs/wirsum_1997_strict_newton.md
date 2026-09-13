# Wirsum 1997 strict Newton

Wirsum (1997) §2.3.3：解 ``J Δx = −F``（Eq. 2.21），再取最小 ``j≥0`` 使

``||F(x + Δx / 2^j)||₂ < ||F(x)||₂``（Eq. 2.23）。

实现：`src/solvers/wirsum_1997_strict.py`。不含 fence、fastox、syngas floor、split/energy merit、refresh-aware accept。

`ReactorConfig.nr_solver_mode='wirsum_1997_strict'` 会关掉上述可由 config 开关的工程桥，并强制 `holdup_k0_anchor='none'`。生产 `solve_global_nr` 仍用缩放 RMS 接受步；原文 ``||F||₂`` 只保证在 `wirsum_1997_strict_newton` 上。

对照：`docs/wirsum_1997_nr_solver_notes.md`；固相 `K=0` 奇异见 `docs/solid_balance_formulation_audit.md`。
