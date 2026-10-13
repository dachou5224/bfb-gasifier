# BFB Gasifier（重建中）

用 Python 复现 Hamel（1999/2001）学位论文中的鼓泡流化床气化模型与算法，并以论文第 7 章的工况，对照 Hamel 本人的模拟结果进行验证。

本仓库于 2026-10-13 按 [`docs/CHARTER.md`](docs/CHARTER.md) 重置主线，新实现将位于 `hamel/` 包中。当前进度见 [`docs/STATUS.md`](docs/STATUS.md)。

## 旧版本

重置前的全部代码、测试、脚本、审计报告与交接记录保存在：

- tag：`archive-2026-10-12`
- 分支：`archive/2026-10-12-snapshot`

```bash
git checkout archive-2026-10-12
```

## 数据

- `data/validation_cases.json`：15 个 Hamel 验证工况（WSV400、VTT、HTW）的输入与实验值。
- `data/hamel_fig7_parity_digitized.json`：第 7 章 parity 图的旧数字化读数（只读，待 T3 复核）。

## 环境

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```
