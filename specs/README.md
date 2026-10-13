# `specs/` 规格目录

Hamel 算法规格入口。**方程全文**在 `docs/hamel_submodels/00–07`；本目录保留已核验摘要。

## 权威文档集（阅读顺序）

1. **主文献**：`docs/_2001_VDI-Dis._…Brennstoffe.pdf`
2. **方程 SoT**：`docs/hamel_submodels/00–07`
3. **摘要入口**：`hamel_checked_algorithm_spec.md`（本目录，本地-only）
4. **流程 / B 档 / 求解语义**：`docs/gasifier_model_flowchart_trilingual.mmd`、`docs/hamel_b_tier_thesis_extraction.md`、`docs/vorab_init_abgleich_vs_outer_preinner.md`
5. **工程补缝**：`docs/hamel_engineering_bridge_notes.md`
6. **同构台账**：`docs/hamel_isomorphic_gap_implementation_checklist.md`（含证据矩阵）

完整索引见 [`docs/README.md`](../docs/README.md)。

## 本目录文件

| 文件 | 角色 |
|------|------|
| `hamel_checked_algorithm_spec.md` | 已核验方程 + 代码映射 + 开放差异（需自行放置，见 `.gitignore`） |

## 已清理（勿再引用）

- 旧规格：`00_implementation_scope.md` … `04_kinetics.md`、`species.md`
- 重复 meta 层：`hamel_original_algorithm_reconstruction.md`、`hamel_evidence_matrix.md`、`hamel_doc_trace_index_latest_to_earliest.md`

源码注释应指向 `docs/hamel_submodels/` 或 `hamel_checked_algorithm_spec.md`，而非已删除的 `specs/0x_*.md`。
