# P5 续：源项 H/O 开口分解

实现：`src/solvers/wirsum_1997_elemental_decomp.py`。

reduced 3-cell LU 的 Vorabrechnung **x0** 上，`combined_element_molar_rates` 给出 C ~ `10⁻¹²`，H/O ~ `10⁻³`–`10⁻²`。

两件分开的事：

1. **质量残差 10¹³**：裸 thesis R6 / 仓库 R12（Hamel Eq. 5.38）在 2.5 MPa、~1300 K 下 `ξ ~ 10¹³` mol/s。`H2 + 0.5 O2 → H2O` **元素中性**，不制造 H/O 净口。
2. **H/O 净口**：
   - bed0：气相原子流 ~ 35 mol-atom H/s，固相 VM+水分约 −35，差 0.007（相对固相 H ~ `2×10⁻⁴`）。`_vm_*_cache` 在这条 init 上是空的，活热解写在 `R_solid`。
   - bed1/2：无固相源项；净口来自 CO/CO₂/CH₄ 等均相反应在大数对消后的残差。

不是 sanity 第 10 项的 `An=b` 不可行。工程 fastox / limiter 是 mixing-control 桥，须单独命名。本核不加 clip。
