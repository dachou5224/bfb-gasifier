# T29：把已裁定的第 1 关输入写入工况文件

> 任务状态只记在 `docs/STATUS.md` 的任务表中（CHARTER v0.10 §8），本卡不设状态栏。

| 字段 | 内容 |
|---|---|
| 执行者 | cursor（新会话） |
| 审查方 | 总司令（数据录入） |
| 前置任务 | ADR 0014、0019（均须已批准）；T18（KURKELA Fig. 1 读数） |
| 分支 | `task/T29-vtt1-case-data` |
| 时间上限 | 1 个工作日 |

## 目标

ADR 0014、0019 已经裁定的第 1 关输入，在 `data/cases/vtt_pressurized.json` 中仍是 `missing`（M1 审计）。本卡只做录入，使工况文件与已批准的裁定一致，T25 装配时不再遇到这些缺口。

## 范围内

逐项录入，每个字段写 `value`、`unit`、`status`（`derived` 或 `cited`）、`source`（ADR 编号与原文页码）、`formula`（若换算）：

1. ADR 0014：`bed_height` = 1,2 m；`cell_heights` = 8 × 0,15 m + 8 × 0,2875 m；进气组成（干空气，按 ADR 0014 的写法）；固体与气体加入温度 298,15 K；床料组成（焦炭与灰）；床料质量 M_bed（由床层压降，A3 的式子）。
2. 燃料粒径分级质量分数：由 `data/digitized/kurkela1992_fig1.json` 的粒级分数录入（ADR 0014 D2）；粒级的代表粒径按 ADR 0014 所写的换算；ADR 0014 若没有写明 <0,1 mm 一级的代表粒径或等效换算，停下上报，不自定。
3. ADR 0019：`distributor_orifice_count` = 139,5（`derived`，n_0 = φ·(D/d_h)²，KURKELA_1991 印刷页 8；note 写明“年代一致性未确认，目标所有者 2026-10-18 批准作跨时期外推假设”）；`distributor_outlet_height` = 0（`derived`，ADR 0019 第 1 项）。第二块板的 95,6 只写在 note 中，供敏感性运行。
4. ρ_P 不录主值：它是外层循环中的量（ADR 0014 A4），在字段 note 中写明由 A4 计算。
5. 同步 `data/cases/vtt_pressurized_gaps.md` 的缺口表。

## 范围外（明确禁止）

- 不录入任何未经 ADR 裁定的值（例如 ε_mf、φ_S，见 X6）。
- 不改 ADR；发现 ADR 与工况文件的字段定义不一致时停下上报。
- 偏离检查的例外条件（CHARTER §9）：`derived` 字段须在 `source` 中给出页码，`formula` 写明换算。

## 交付物

- `data/cases/vtt_pressurized.json`、`data/cases/vtt_pressurized_gaps.md` 的修改
- 若装配工具（`hamel/cases/assemble.py`）有对应检查，`pytest` 通过
- `docs/STATUS.md` 中 T29 一行

## 验收证据

- `python3 tools/drift_check.py` 通过（新字段满足 §9 的例外条件）。
- 逐字段的 ADR 引用可在原 ADR 中找到。

## 停止条件

到达时间上限，或需要越出范围时，停下并在 PR / STATUS 中写明。
