# X3：输入来源

本地论文是 `docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf`。表 7.1 在印刷页 128，图 7.4 在印刷页 129，扫描见同目录 `table71.png`。

表 7.1 有三列：Wesseling 空气/蒸汽、Wesseling 氧/蒸汽、Berrenrath 氧/蒸汽。列出的是燃料、压力、给料、气化剂类型、\(\lambda_{O_2}\)、是否循环，以及水分、灰、挥发分和干基 C/H/O/N/S。表上没有蒸汽 kg/h、蒸汽温度、热值、筛分质量、二次流 Nm³/h 和喷口坐标。

表下脚注是：气化剂最多从 10 个不同层面进入反应器。这是进气层面，不是 10 个给料粒级。同页正文写自由空间加入空气/蒸汽混合物后局部温度升高。图 7.4 上的温度峰没有在这一轮重新数字化。旧的 \(\xi=0.6\) 不是喷口测量。

表题引用 Hamel、Krumm、Abraham、Mittelstädt、Karkowski（1999），第四届煤科学国际会议，北京，第 411–420 页。这篇会议论文的全文不在本工作区。表 7.1 实测数据的直接来源是 Rheinbraun 与 Krupp Uhde 1997 年的内部报告，不公开。Adlhoch 等（1992）和 Schumacher 等（1993）是装置背景，本工作区的书目记录里没有把缺的蒸汽、热值、筛分或二次流流量补上。

因此这些缺项仍空着。数值目标里的蒸汽 0、壁面热损 0 和示意生成焓 \(-3.96\times10^{6}\) J/kg 是假设，不是测量。加载接口不把 30.1 或 37.5 Nm³/h 从总风里扣出。

实验比较还不能做。顺序是：先在这组条件下把 Wesseling 1 的数值方程求到 \(s=1\)，再接入真实蒸汽、热值、筛分和分级进料，并单独做多 cell 映射，然后才和出口测量比。Wesseling 2 用同一模型和同一求解策略作同装置对照。两区核过门不等于生产 `Reactor` 已认证。当前路径停在联合通过的 \(s=0.32\)。
