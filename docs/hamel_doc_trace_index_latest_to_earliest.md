# Hamel 文档溯源索引（清理后）

目的：给出清理后目录中“可直接用于原版算法核验”的最小文档集。

## 当前保留文档（按用途）

| 文档 | 角色 | 证据口径 |
|---|---|---|
| `docs/hamel_original_algorithm_reconstruction.md` | 原版算法主说明 | 论文公式/章节 + 源码对照 |
| `docs/hamel_evidence_matrix.md` | 模块-公式-实现矩阵 | 仅保留论文 PDF 与源码锚点 |
| `docs/hamel_isomorphic_gap_implementation_checklist.md` | 同构落地清单 | 仅保留实现差距与验收项 |
| `specs/hamel_checked_algorithm_spec.md` | 实现细节规格入口 | 逐条列出已核验方程与代码映射 |

## 主源与边界

- 主源：`docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf`
- 本索引不再依赖已清理历史文档。
- 任何条目若标注“部分一致/不一致”，均不得用于“完全复现”声明。
