---
name: bfb-diagnosing-bugs
description: BFB 气化炉数值与物理 bug 诊断循环。用于 NR 不收敛、物种预算偏差、单格/全炉异常、重构回归。先建 tight 反馈环再假设；调参前须排除 Hamel parity 问题。
---

# bfb-diagnosing-bugs（BFB 诊断循环）

Matt Pocock `diagnosing-bugs` 的 BFB 适配版：**先建可重复的 tight 反馈环**，再最小化、假设、插桩、修复。禁止无环读代码猜原因。

## 与兄弟技能的分工

| 技能 | 何时用 |
|------|--------|
| **本技能** | 有具体症状（不收敛、CO/CO₂ 偏、守恒 gap、回归）时的**调试纪律** |
| `bfb-numerical-guard` | 已定位到求解器/Jacobian/缩放/`np.exp` 时的**检查清单** |
| `bfb-hamel-parity-audit` | 怀疑**算法与 Hamel 不一致**（非纯数值）时的方程级审计 |
| `bfb-scientific-audit` | 修复完成、合并前的**门禁** |

**调参门禁**（`docs/CLAUDE.md` §一致性优先）：若症状可能是**方程/单位/closure 错误**，先跑 `bfb-hamel-parity-audit` 或 `audit_hamel_consistency.py`，**不得**用 scale 掩盖 parity 问题。

## 何时启用

- global NR / fsolve **不收敛**、**假收敛**（`dT` 小但 `rms_scaled` 大）
- **物种/能量预算**异常（CO/CO₂、炭、O₂ 分配、bed-top gap）
- **单格可解、全炉不可解**（或反之）
- 重构/改开关后的**回归**（`benchmark_phase2_*` JSON 变差）
- 用户说「诊断 / debug / 排查 / 为什么 LU 不收敛」

探索代码前，读 **`docs/vorab_init_abgleich_vs_outer_preinner.md`**（外层 Init / Pre-Inner / Check2 分工）与 **`docs/module_dependency_overview.md`**（模块依赖）。

---

## Phase 1 — 建反馈环（本技能核心）

没有 **tight、red-capable** 的单条命令，禁止进入 Phase 2。

### BFB 建环顺序（由快到慢）

1. **Sanity 基线**（排除低级破坏）  
   `cd bfb-gasifier && python3 tests/sanity_checks.py`

2. **Seam 级 pytest**（单方程/单模块）  
   见 [`SYMPTOM-MAP.md`](SYMPTOM-MAP.md) §A；例：  
   `pytest tests/test_mass_transfer.py tests/test_cell_kinetics.py -q --tb=short`

3. **专项 audit 脚本**（断言具体物理量）  
   按主症状从 [`SYMPTOM-MAP.md`](SYMPTOM-MAP.md) 选一条；例：  
   `python3 scripts/audit_bed_co_co2_budget_lu.py`

4. **Throwaway harness**（最小复现）  
   单格 `solve_cell`、单速率 `rate_R5_*`、固定 `T,P,C_*` 的 10 行脚本；勿一上来跑全炉。

5. **Differential loop**（重构回归）  
   同一 LU case + 同一 ladder level，diff 两次 `data/phase2_*` JSON 的 `rms_scaled_final`、出口组分、bed-top char gap。

6. **Convergence ladder**（全炉，慢）  
   `python3 scripts/benchmark_phase2_convergence_ladder.py --level p2 --max-global-iter 20`  
   **硬超时 420s**；超时记 `TIMEOUT`，勿无限等待。

7. **Git bisect**（已知好/坏提交）  
   `git bisect run` + 上表中最窄 audit/pytest 作为 exit code 探针。

8. **模块审计入口**（症状模糊）  
   `python3 scripts/run_module_audits.py --module <physics|kinetics|thermal|solver|...>`

**不要**用 curl/Playwright；**不要**把「能跑完 app.py」当作 green。

### Tighten the loop

- **更快**：`--max-global-iter 1`、单 bed 审计、`pytest -k` 子集
- **更尖**：断言用户**原话症状**（如 `bed-top char gap < 1e-4 kg/s`），非泛化 `exit code == 0`
- **更稳**：固定 `LADDER_P2` / `_lu_harness` 配置；记录 `data/*.json` 路径

### 非确定性

NR 线搜索偶发失败：连跑 N 次统计 `line_search_failures`；或 pin `numpy` 随机种子（若路径含随机初值）。

### 无法建环时

明确列出已试命令；向用户索要：**(a)**  repro 用的 ladder level / `--max-global-iter`，**(b)** 已有 `data/phase2_*.json` 或 log，**(c)** 是否允许临时 `[DEBUG-bfb-…]` 插桩。**无环不进 Phase 2。**

### Phase 1 完成标准

已运行至少一次，且满足：

- [ ] **Red-capable** — 断言**用户的具体症状**（残差阈值、物种分数、守恒 gap）
- [ ] **Deterministic** — 同配置同 verdict（或 flake 率已量化）
- [ ] **Fast** — 内层假设验证用 **seam/audit**（秒~十秒）；全炉仅作确认
- [ ] **Agent-runnable** — 一条完整 shell 命令，输出已粘贴

---

## Phase 2 — 复现 + 最小化

- 确认 loop 失败模式与用户描述**一致**（勿 fix 错 bug）
- **最小化**：减 `--max-global-iter`、关自由板/回收、单 zone audit；每砍一步重跑 loop
- 完成当**每个剩余要素都 load-bearing**

### BFB 常见最小化轴

| 轴 | 做法 |
|----|------|
| 空间 | 单 cell → 单 bed → 全炉 |
| 化学 | 关 R10 / 关 freeboard reactions 对照 |
| 求解 | init only（`max-global-iter 1`）vs 多 outer |
| 配置 | `LADDER_LEVELS["baseline"]` vs `p2` |

---

## Phase 3 — 假设（先列再测）

生成 **3–5 条可 falsify 假设**，展示给用户后再测。BFB 高频假设模板：

- 「若 **Stale freeboard closure** 是原因，则 refresh bed-top char 后 gap→0 但 RMS 不变/变」
- 「若 **Pre-Inner 与 Init Abgleich 职责重叠** 是原因，则只开/关 `run_outer_preinner_start_value_bridge` 会移动 init split RMS」
- 「若 **R5 suspension 缺 C_H2O^0.5** 是原因，则 `audit_hamel_consistency` R5d 项 FAIL」
- 「若 **Jacobian 漏 side-element 耦合** 是原因，则 dense vs structured `nnz` diff ≠ 0（freeze 上下文）」
- 「若 **假收敛** 是原因，则 `--strict` 下 `converged_inner_nr=False` 但 outer 已停」

**展示排序后再 instrument。** 若假设指向方程/文献，转 `bfb-hamel-parity-audit`，勿继续调 scale。

---

## Phase 4 — 插桩

- **一次只改一个变量**；每条 log 对应 Phase 3 一条预测
- 日志前缀 **`[DEBUG-bfb-xxxx]`**，收尾 `rg '\[DEBUG-bfb-'` 删净
- 优先：**断点/REPL** > 边界 log（bed-top、Check1/Check2 前后 `rms_scaled`）> 禁止 log 全残差向量
- **性能回归**：用 `perf_counter` 包 `global_nr` 单步；勿靠肉眼

### 求解器专项

插桩 NR 时同步查阅 **`bfb-numerical-guard`**（缩放、`np.exp`、Jacobian pattern、假收敛）。

### 外层语义

Check1/Check2/Pre-Inner 混淆时，**先读** `docs/vorab_init_abgleich_vs_outer_preinner.md`，再插桩 `outer_loop.py` / `reactor._solve_global_nr`。

---

## Phase 5 — 修复 + 回归

1. 若有 **correct seam**（pytest 或 audit 脚本），**先写 failing test/断言**，再改代码
2. seam 太浅（全炉才现、单格不现）→ 记录为架构发现，回归用 **audit JSON 快照** 或 ladder 子集
3. 修复后依次跑：最小 loop → 原 loop → `sanity_checks.py` → 相关 module audit
4. 若动方程/单位，补 **`bfb-hamel-parity-audit`** 或更新 gap checklist

**禁止**以「出口更接近 Table 2」 alone 证明 fix；须 loop green + 科学门禁。

---

## Phase 6 — 清理 + 事后

- [ ] 原 repro loop green
- [ ] 回归 test/audit 通过（或 documented 无 seam）
- [ ] 无 `[DEBUG-bfb-…]` 残留
- [ ] 临时 harness 已删或移入 `scripts/_deprecated/`
- [ ] 结论写入 commit/PR：**根因 + 哪条假设成立**

**事后一问**：能否用更窄的 audit 永久锁定？若缺 seam，建议在 `tests/` 或 `scripts/audit_*.py` 增一条（不强制本轮实现）。

---

## 输出模板

```markdown
## BFB 诊断 — {症状摘要} — {日期}

**Phase 1 loop**: `{命令}` → {RED/GREEN/TIMEOUT}  
**最小 repro**: …  
**根因**: …  
**假设命中**: #N — …  
**修复**: {文件}  
**回归**: {pytest/audit 命令} → PASS  
**parity/调参**: {已审计 / 未闭合，禁止调参}  
**开放**: …
```

---

## 相关文件

- [`SYMPTOM-MAP.md`](SYMPTOM-MAP.md) — 症状 → 脚本/pytest 索引
- `scripts/_lu_harness.py` — Phase2 LU 配置 SSOT（`LADDER_LEVELS`）
- `scripts/run_module_audits.py` — 分模块 pytest + audit
- `docs/vorab_init_abgleich_vs_outer_preinner.md` — 外层求解语义
- `docs/optimization_handoff_summary.md` §8–9 — 已知陷阱与建议顺序
- `.cursor/skills/bfb-numerical-guard/SKILL.md`
- `.cursor/skills/bfb-hamel-parity-audit/SKILL.md`
