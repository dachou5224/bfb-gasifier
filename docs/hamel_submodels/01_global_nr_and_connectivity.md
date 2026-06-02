# 全局 Newton-Raphson 与连接矩阵（最小注释版）

## 论文编号锚点

- Kapitel 2.3
- `Bild 2.2`
- `(2.9)`

## 原文直引（OCR）

- `[可直接采信] ... Newton-Raphson-Algorithmus ...`
- `[可直接采信] ... blocktridiagonale Struktur der Funktionalmatrix ...`
- `[可直接采信] ... Nebenelementen ... Feststoffrezirkulation ...`
- `[可直接采信] ... gilt ... (2.9)`
- `[可直接采信] ... „gedämpfte“ Variante des Newton-Verfahrens ...`

## 算法与机理（按论文原文重述）

1. 由第 2.2 节各类守恒方程组装得到**非线性方程组**（Kapitel 2.3）。
2. 由于单元沿主流向一维排列，且每个单元主要依赖“本单元 + 相邻单元”变量，函数矩阵（Jacobian）呈**块三对角结构**（blocktridiagonale Struktur）。
3. 若存在来自非相邻单元的回流（论文举例：旋风分离器返回床层的固体再循环），则产生**边带非零项**（Nebenelemente）。
4. 论文采用 Wozny (1983) 的 Newton-Raphson 解法，并采用 Wirsum (1998) 的**阻尼 Newton**改善远离解初值下的收敛。

## 公式锚点 (2.9)

按 `PAGE 20` OCR 可直接确认的更新关系为：

$$x^{v+1} = x^{v} + \Delta x^{v} \qquad (2.9)$$

> 注：`(2.9)` 在提取文本中与函数矩阵线性化式存在版面断裂，但更新式与“Newton + 阻尼变体”的求解框架可直接对应原文段落。

## 参数/结构来源口径

- 结构来源：`PAGE 20` 对 block-tridiagonal + side-elements（recirculation）的文字定义。
- 算法来源：同页明确给出 Newton-Raphson 与 damped Newton 的文献归属（Wozny 1983; Wirsum 1998）。
- 本文档不引入论文外的矩阵重排、预条件或伪代码。

## 一致性约束

- 仅使用 `(2.9)` 编号，不改写为自定义公式编号。
- 不给出论文外矩阵重排与伪代码。

## PDF 页锚点

- `===== PAGE 20 =====`

## 修改标注（2026-05-07）

- `[M01-1]` 增补“算法与机理”四步链路，明确 block-tridiagonal 与 side-elements 的物理来源。
- `[M01-2]` 增补 `(2.9)` 的可确认更新式与 OCR 断裂说明，避免论文外推导。
- `[M01-3]` 增补“参数/结构来源口径”段，明确来源页与文献归属。
