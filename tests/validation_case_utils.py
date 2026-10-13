"""HTW Wesseling Table 2 LU（``CASE_HTW_WESSELING_1``）验证数据加载与 Phase 1 基线配置。

单一事实来源供 ``test_table2_LU.py``、``scripts/audit_phase1_htw_lu.py`` 等复用，
避免 JSON 展平与 ``ReactorConfig`` 字段在多处漂移。

**Phase 1 基线**（与 ``test_table2_LU`` 中 slow 用例一致）：
- ``n_cells=10``, ``H_bed=6.0`` m（JSON ``bed_height_m``，非全炉 ``height_m``）, ``H_freeboard=0``（自由板未集成）
- ``heat_loss_frac=0.1``（``ReactorConfig`` 默认，中试炉热损失文献带 5–10% 见 validation_cases 注释）
- ``u0_target=None``（表观气速走内部默认/闭环）
- 求解：``Reactor.solve(..., solver=\"global_nr\", ...)``（NR-only；历史 GS 路径已移除，见 ``scripts/_deprecated/nr/``）

**Phase 1 global NR 开发口径**：
- 与同一 LU shared config 复用几何/进料/动力学窗口
- 求解：``Reactor.solve(..., solver=\"global_nr\", max_global_iter=20, tol_global=1.0, nr_init_strategy=\"vorabrechnung\", nr_jacobian_strategy=\"block_tridiag_structured\")``
- 该口径当前用于 NR 审计与开发验证，尚不等同于 GS 基线稳定门控

**Phase 2 freeboard-aware 开发口径**：
- 在 shared ``global_nr`` 配置基础上恢复 ``H_freeboard = reactor_height - H_bed``
- 将床层顶部之后的 freeboard 作为独立 gas-only 轴向段求解
- Jacobian 通过 ``nr_jacobian_strategy=None`` 委托 ``Reactor`` 的拓扑默认选择
  （等价于开发口径里的“auto”，side-block 路径避免误用纯 block-tridiag）
- 用于区分 ``bed exit`` 与 ``reactor exit`` 的验证/审计

Source: ``docs/hamel_submodels/00_readme_and_citation_rules.md``；
``docs/hamel_submodels/03_hydrodynamics_core_chain.md``；
``docs/hamel_submodels/04_freeboard_particle_trajectory_analytic.md``；
``data/validation_cases.json``

实现已拆分至 ``tests/validation/``；本模块仅 re-export 以保持 import 路径不变。
"""

from __future__ import annotations

from tests.validation.case_lu_loader import (
    CASE_LU_VALIDATION_KEY,
    MULTI_CASE_MATRIX_KEYS,
    air_nm3_h_to_species_mol_s,
    estimate_gas_feeds,
    extract_case_inlet_topology_hints,
    flatten_validation_case,
    json_numeric_or_none,
    load_case_LU,
    load_validation_case_flat,
    load_validation_case_node,
    load_validation_json_root,
    nm3_h_to_mol_s_stp,
    o2_nm3_h_to_species_mol_s,
    strict_validation_gate,
    validation_numeric_tolerances,
)
from tests.validation.config_mutators import (
    apply_case_secondary_stream,
    apply_secondary_air_repartition,
)
from tests.validation.config_phase1 import (
    build_phase1_htw_lu_global_nr_reactor_config,
    build_phase1_htw_lu_refined_config,
    build_phase1_htw_lu_reactor_config,
)
from tests.validation.config_phase2_freeboard import (
    VTT_MINIMAL_LOOP_SOLVE_KWARGS,
    apply_enable_r12_temperature_seed_policy,
    build_phase2_htw_lu_freeboard_reactor_config,
    build_phase2_htw_lu_local_refined_damped_config,
    build_vtt_minimal_loop_nr_acceptance_config,
)
from tests.validation.solve_kwargs import (
    PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS,
    PHASE1_HTW_LU_SOLVE_KWARGS,
    PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS,
)

__all__ = [
    "CASE_LU_VALIDATION_KEY",
    "MULTI_CASE_MATRIX_KEYS",
    "PHASE1_HTW_LU_GLOBAL_NR_SOLVE_KWARGS",
    "PHASE1_HTW_LU_SOLVE_KWARGS",
    "PHASE2_HTW_LU_FREEBOARD_SOLVE_KWARGS",
    "VTT_MINIMAL_LOOP_SOLVE_KWARGS",
    "air_nm3_h_to_species_mol_s",
    "apply_case_secondary_stream",
    "apply_enable_r12_temperature_seed_policy",
    "apply_secondary_air_repartition",
    "build_phase1_htw_lu_global_nr_reactor_config",
    "build_phase1_htw_lu_refined_config",
    "build_phase1_htw_lu_reactor_config",
    "build_phase2_htw_lu_freeboard_reactor_config",
    "build_phase2_htw_lu_local_refined_damped_config",
    "build_vtt_minimal_loop_nr_acceptance_config",
    "estimate_gas_feeds",
    "extract_case_inlet_topology_hints",
    "flatten_validation_case",
    "json_numeric_or_none",
    "load_case_LU",
    "load_validation_case_flat",
    "load_validation_case_node",
    "load_validation_json_root",
    "nm3_h_to_mol_s_stp",
    "o2_nm3_h_to_species_mol_s",
    "strict_validation_gate",
    "validation_numeric_tolerances",
]
