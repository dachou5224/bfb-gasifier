# Wirsum 1997 对齐 P5：LU 重诊

实现：`src/solvers/wirsum_1997_lu_reentry.py`。
默认只看 reduced 3-cell Vorabrechnung **x0**，不跑 `reactor.solve()`，不追 150-iter 收敛。

## 清单（固定顺序）

1. elemental closure — 气+固源项 C/H/O/N/S [mol-atom/s]
2. mass closure — 气相 / 固相 cell residual
3. energy closure — 能量 residual [W]
4. inner NR — packed residual（可选 Jacobian）
5. outer — 单帧 `ε_b, K_bd, K_s, u_b, T`（P4 核；一次快照不算 ρ）
6. freeboard — 生产名仍解析到 `engineering_freeboard`
7. composition observation — 对照 **P0 baseline**，不是实验
8. temperature observation — 对照 **P0 轴向 T**，不是测点

第 7、8 项的 `role` 是 `observation_vs_p0_baseline`。禁止用实验 CO 批准 clip、renormalize 或 rate scale。

## x0 已看到的开口

- C 源项近机器精度；H/O 开口，与 sanity「Pyrolysis elemental allocator」一致
- 气相 residual 与 P0 packed `residual_l2 ~ 1e13` 同类
- 能量 residual 在 MW 量级
- `analytical_wirsum` 仍是工程别名；`wirsum_1997_exact` 未接入全炉
