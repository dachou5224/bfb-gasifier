# WSV400 图上读数清单

印刷数字已经写入 `wsv400.json`。下面两项要由目标所有者用数字化工具读取。不要用目测代替。

## Bild 7.14（PDF 81 左 / S. 142）

图题：Wirbelschichtvergaser WSV 400，模拟所用的 cell 划分。图内文字框已经印有 Betriebsdruck = 0,1 MPa、Reaktordurchmesser = 0,4 m / 0,6 m、Reaktorhöhe = 3 m、Brennstoffleistung = 100 kW、Vergasungsmittel = Luft。这些与 §7.1.5 正文一致，不在本清单里重读。

待读：

1. 从 Düsenboden 向上，cell 的个数。锥形过渡段如果单独成格，单独计。
2. 每一格的高度，单位 m。若图上没有尺寸线，按图比例读出并记下标定。
3. 床面记号（Wirbelschicht 旁的液面符号）到底部分配器的高度。这是表 2.2 的 Wirbelbetthöhe。
4. 燃料箭头、气化剂箭头相对分配器的高度。正文已写燃料在 Düsenboden 以上 200 mm；若图上位置与这句话不一致，两个数都保留。
5. 除底部气化剂外，是否还有第二处进气箭头及其高度。

## Bild 7.15（PDF 82 左 / S. 144）

图题写明：木片样品折成等表面积圆球的粒径分布，图内文字框为 Holzhackschnitzelprobe zum Versuch A。只用于工况 A。

横坐标是印刷的当量直径分级，已记入 `wsv400.json` 的 `fuel_size_class_edges`：bis 2、2–5、5–8、8–11、11–14、14–17、17–20、20–23、23–26 mm。

待读：这 9 根柱的质量分数，单位 kg/kg（图上为 Gew.-%）。纵轴标题是 Massenanteile in Gew.-%。按柱顶读取，并记下网格线像素。
