# Wirsum 1997 strict Newton

Wirsum (1997) §2.3.3：解 ``J Δx = −F``（Eq. 2.21），再取最小 ``j≥0`` 使

``||F(x + Δx / 2^j)||₂ < ||F(x)||₂``（Eq. 2.23）。

实现：`src/solvers/wirsum_1997_strict.py`。不含 fence、fastox、syngas floor、split/energy merit、refresh-aware accept。

`ReactorConfig.nr_solver_mode='wirsum_1997_strict'` 会关掉上述可由 config 开关的工程桥，并强制 `holdup_k0_anchor='none'`。生产 `solve_global_nr` 仍用缩放 RMS 接受步；原文 ``||F||₂`` 只保证在 `wirsum_1997_strict_newton` 上。

对照：`docs/wirsum_1997_nr_solver_notes.md`；固相 `K=0` 奇异见 `docs/solid_balance_formulation_audit.md`。

## P2.1 四组 A/B（三格 holdup，`x0 = 0.05 kg`）

机器表：`data/wirsum_1997_strict_ab.json`。

| 问题 | A raw+fresh | B scaled+fresh | C raw+frozen | D scaled+frozen |
|---|---|---|---|---|
| 线性 `K>0` | 1 步，α=1 | 同 A | 同 A（J 为常数） | 同 A |
| 非线性 MM 汇 | 5 步，α=1 | 同 A | 40 步未到 `1e-10`（`max_iter`） | 同 C |
| MM 行加权 `(1, 10³, 10⁻³)` | 8 步，α 含 1/8 | 5 步，α=1 | 第 2 步线搜索失败 | 40 步爬行 |
| `K=0` 无出流本构 | `singular_jacobian`，rank 0，3 零行 | — | — | — |

结论：线性闭合问题不需要缩放或刷新 J。非线性但尺度均匀时，难度在 **frozen Jacobian**，不在 raw vs scaled。残差行尺度相差 10⁶ 时，**raw `||F||₂` 会强迫回溯**，scaled RMS 仍满步；frozen J 在 raw merit 下直接线搜索失败。`K=0` 的零行是质量平衡没闭合，不是 Newton 实现问题。
