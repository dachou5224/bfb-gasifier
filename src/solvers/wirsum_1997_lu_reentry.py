"""P5：在干净架构上重诊 LU，按 Wirsum 对齐清单，不追收敛。

顺序固定：

1. elemental closure
2. mass closure
3. energy closure
4. inner NR
5. outer Vorabrechnung
6. freeboard closure
7. composition observation
8. temperature observation

第 7、8 项只对照 P0 baseline，**不用实验 CO/T 给 numerical bridge 背书**。
默认只看 Vorabrechnung x0，不跑 ``reactor.solve()``。

Source: Wirsum (1997) §2.2–2.3；handoff Phase P5
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from src.core.cell import S_CHAR, S_MOISTURE, S_VM
from src.core.composition import wet_to_dry_mole_fractions
from src.core.elemental_ledger import (
    TRACKED_ELEMENTS,
    combined_element_molar_rates,
    max_absolute_element_rate,
)
from src.core.species import GAS_SPECIES
from src.physics.wirsum_1997_freeboard import resolve_freeboard_model
from src.solvers.wirsum_1997_outer_diagnostics import snapshot_vorab_fields

CHECKLIST_ORDER: tuple[str, ...] = (
    "elemental_closure",
    "mass_closure",
    "energy_closure",
    "inner_nr",
    "outer_convergence",
    "freeboard_closure",
    "composition_observation",
    "temperature_observation",
)

# 这些键不得出现在 composition / temperature 的判据里
FORBIDDEN_BRIDGE_KEYS: tuple[str, ...] = (
    "experiment_CO",
    "measured_CO",
    "y_CO_target",
    "calibration_objective",
    "lu_plant_parity",
)

EXIT_SPECIES: tuple[str, ...] = ("H2", "CO", "CO2", "CH4")


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _elemental_rates_from_cell(cell: Any) -> dict[str, float]:
    gas_rates = np.asarray(cell.R_gas_b, dtype=np.float64) + np.asarray(cell.R_gas_d, dtype=np.float64)
    return combined_element_molar_rates(
        gas_rates=gas_rates,
        solid_rates=cell.R_solid,
        fuel_type=cell.fuel_type,
        ash_dry_wt=cell.solid.ash_dry_wt,
        C_dry=cell.solid.C_dry,
        H_dry=cell.solid.H_dry,
        O_dry=cell.solid.O_dry,
        nitrogen_fraction=cell.solid.nitrogen_fraction,
        sulfur_fraction=cell.solid.sulfur_fraction,
        sulfur_volatile_frac=cell.solid.sulfur_volatile_frac,
        char_index=S_CHAR,
        vm_index=S_VM,
        moisture_index=S_MOISTURE,
    )


def _refresh_cell_residuals(cell: Any) -> None:
    cell.residuals()


def audit_elemental_closure(cells: list[Any]) -> dict[str, Any]:
    """气+固源项的 C/H/O/N/S 净原子流 [mol-atom/s]。"""
    per_cell: list[dict[str, Any]] = []
    max_abs = 0.0
    max_elem = ""
    for i, cell in enumerate(cells):
        rates = _elemental_rates_from_cell(cell)
        local = max_absolute_element_rate(rates)
        if local >= max_abs:
            max_abs = float(local)
            max_elem = max(rates, key=lambda k: abs(float(rates[k])))
        per_cell.append(
            {
                "cell": i,
                "rates_mol_atom_s": {el: _finite(rates[el]) for el in TRACKED_ELEMENTS},
                "max_abs_mol_atom_s": _finite(local),
            }
        )
    return {
        "item": "elemental_closure",
        "units": "mol-atom/s",
        "per_cell": per_cell,
        "max_abs_mol_atom_s": _finite(max_abs),
        "dominant_element": max_elem,
        "open": bool(max_abs > 1e-8),
        "note": (
            "C 近机器精度而 H/O 开口，对应 sanity「Pyrolysis elemental allocator」。"
            "这是源项闭合，不是实验 CO 偏差。"
        ),
    }


def audit_mass_closure(cells: list[Any]) -> dict[str, Any]:
    """气相 / 固相 cell residual 的欧氏范数。"""
    per_cell: list[dict[str, Any]] = []
    gas_max = 0.0
    solid_max = 0.0
    for i, cell in enumerate(cells):
        F_gas = np.asarray(cell.residuals(groups=("gas",)), dtype=np.float64)
        F_solid = np.asarray(cell.residuals(groups=("solid",)), dtype=np.float64)
        gas_l2 = float(np.linalg.norm(F_gas))
        solid_l2 = float(np.linalg.norm(F_solid))
        gas_max = max(gas_max, gas_l2)
        solid_max = max(solid_max, solid_l2)
        per_cell.append(
            {
                "cell": i,
                "gas_residual_l2": _finite(gas_l2),
                "solid_residual_l2": _finite(solid_l2),
            }
        )
    return {
        "item": "mass_closure",
        "units": {"gas": "mol/s residual L2", "solid": "kg/s residual L2"},
        "per_cell": per_cell,
        "max_gas_residual_l2": _finite(gas_max),
        "max_solid_residual_l2": _finite(solid_max),
        "open": bool(gas_max > 1.0 or solid_max > 1e-3),
        "note": "x0 气相残差量级与 P0 packed residual 同类；不据此调 rate multiplier。",
    }


def audit_energy_closure(cells: list[Any]) -> dict[str, Any]:
    """能量 residual [W]。"""
    per_cell: list[dict[str, Any]] = []
    max_abs = 0.0
    for i, cell in enumerate(cells):
        F_e = np.asarray(cell.residuals(groups=("energy",)), dtype=np.float64)
        val = float(F_e.reshape(-1)[0])
        max_abs = max(max_abs, abs(val))
        per_cell.append({"cell": i, "energy_residual_W": _finite(val)})
    return {
        "item": "energy_closure",
        "units": "W",
        "per_cell": per_cell,
        "max_abs_W": _finite(max_abs),
        "open": bool(max_abs > 1.0),
        "note": "能量口未闭合时先看 Eq.2.7 残差，不先加 T clip。",
    }


def audit_inner_nr(reactor: Any, *, include_jacobian: bool = False) -> dict[str, Any]:
    """x0 上的 packed residual；默认不算 Jacobian。"""
    from src.solvers.global_nr_solver import build_equation_scales, global_residual, pack_reactor
    from src.solvers.nr_indexing import n_var

    cfg = reactor.config
    cells = reactor._solver_cells_for_nr()
    x = pack_reactor(cells)
    F = global_residual(x, cells, reactor._apply_all_bc_for_nr)
    eq_scale = build_equation_scales(
        cells,
        ref_gas_mol_s=float(cfg.O2_feed + cfg.H2O_feed + cfg.N2_feed),
        ref_solid_kg_s=float(cfg.fuel_feed),
        ref_energy_W=float(cfg.fuel_feed) * 20e6,
    )
    F_hat = F / eq_scale
    out: dict[str, Any] = {
        "item": "inner_nr",
        "solved": False,
        "n_var": int(sum(n_var(c) for c in cells)),
        "residual_l2": _finite(np.linalg.norm(F, ord=2)),
        "residual_linf": _finite(np.linalg.norm(F, ord=np.inf)),
        "rms_scaled": _finite(np.sqrt(np.mean(np.square(F_hat)))),
        "note": "只诊断 x0 residual，不跑 Newton，不追收敛。",
    }
    if include_jacobian:
        from src.solvers.global_nr_solver import build_jacobian_fd

        _J, meta = build_jacobian_fd(
            x,
            F,
            cells,
            reactor._apply_all_bc_for_nr,
            eq_scale,
            strategy="block_tridiag_structured",
        )
        out["jacobian_zero_rows"] = int(meta.get("zero_rows", 0))
        out["jacobian_zero_cols"] = int(meta.get("zero_cols", 0))
        out["jacobian_nnz"] = int(meta.get("nnz", 0))
    return out


def audit_outer_convergence(cells: list[Any]) -> dict[str, Any]:
    """单帧 Vorabrechnung 场。一次快照不能判定 ρ。"""
    snap = snapshot_vorab_fields(cells)
    return {
        "item": "outer_convergence",
        "status": "single_frame_only",
        "n_cells": snap.n_cells,
        "snapshot": snap.as_lists(),
        "contractivity": None,
        "note": "P4 核已就绪；全炉 outer 历史尚未接入生产 loop。",
    }


def audit_freeboard_closure(cfg: Any, *, n_freeboard_cells: int) -> dict[str, Any]:
    configured = str(getattr(cfg, "freeboard_trajectory_model", "analytical_wirsum"))
    spec = resolve_freeboard_model(configured)
    return {
        "item": "freeboard_closure",
        "configured": configured,
        "resolved_name": spec.name,
        "u_p0_over_u_b": spec.u_p0_over_u_b,
        "ejection": spec.ejection,
        "n_freeboard_cells": int(n_freeboard_cells),
        "wirsum_1997_exact_wired": False,
        "note": "生产仍是 engineering 别名；exact 核未接入 simulate_freeboard。",
    }


def _exit_observation(cells: list[Any]) -> dict[str, Any]:
    top = cells[-1]
    y_wet = top._mole_fractions("combined")
    wet = {sp: float(y_wet[j]) for j, sp in enumerate(GAS_SPECIES)}
    dry = wet_to_dry_mole_fractions(wet)
    return {
        "exit_T_K": _finite(top.T),
        "bed_T_K": [_finite(c.T) for c in cells],
        "exit_gas_dry": {sp: _finite(dry.get(sp)) for sp in EXIT_SPECIES},
        "exit_gas_wet": {sp: _finite(wet.get(sp)) for sp in EXIT_SPECIES},
    }


def _delta_vs_p0(current: dict[str, Any], baseline: dict[str, Any] | None) -> dict[str, Any] | None:
    if not baseline:
        return None
    out: dict[str, Any] = {}
    if "exit_T_K" in current and "exit_T_K" in baseline:
        out["d_exit_T_K"] = _finite(current["exit_T_K"] - baseline["exit_T_K"])
    cur_dry = current.get("exit_gas_dry") or {}
    base_dry = baseline.get("exit_gas_dry") or {}
    out["d_exit_gas_dry"] = {
        sp: _finite(cur_dry.get(sp) - base_dry[sp]) if base_dry.get(sp) is not None else None
        for sp in EXIT_SPECIES
    }
    return out


def audit_composition_observation(
    cells: list[Any],
    *,
    p0_exit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    obs = _exit_observation(cells)
    return {
        "item": "composition_observation",
        "role": "observation_vs_p0_baseline",
        "forbidden": list(FORBIDDEN_BRIDGE_KEYS),
        "exit_gas_dry": obs["exit_gas_dry"],
        "exit_gas_wet": obs["exit_gas_wet"],
        "delta_vs_p0": _delta_vs_p0(obs, p0_exit),
        "note": "对照 P0 x0，不是实验组成，不能用来批准 clip / renormalize / rate scale。",
    }


def audit_temperature_observation(
    cells: list[Any],
    *,
    p0_exit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    obs = _exit_observation(cells)
    return {
        "item": "temperature_observation",
        "role": "observation_vs_p0_baseline",
        "forbidden": list(FORBIDDEN_BRIDGE_KEYS),
        "exit_T_K": obs["exit_T_K"],
        "bed_T_K": obs["bed_T_K"],
        "delta_vs_p0": _delta_vs_p0(obs, p0_exit),
        "note": "对照 P0 轴向 T，不是工业测点，不据此加 fence。",
    }


def load_p0_reduced_exit(path: str | Path | None = None) -> dict[str, Any] | None:
    root = Path(__file__).resolve().parents[2]
    json_path = Path(path) if path is not None else root / "data" / "wirsum_alignment_p0_baseline.json"
    if not json_path.is_file():
        return None
    data = json.loads(json_path.read_text(encoding="utf-8"))
    case = (data.get("cases") or {}).get("reduced_3cell_lu_init")
    if not isinstance(case, dict):
        return None
    return case


def audit_lu_reentry(
    reactor: Any,
    *,
    p0_exit: dict[str, Any] | None = None,
    include_jacobian: bool = False,
) -> dict[str, Any]:
    """按固定清单审计已初始化的 LU reactor。不调用 ``solve()``。"""
    if p0_exit is None:
        p0_exit = load_p0_reduced_exit()
    reactor._apply_all_bc_for_nr()
    cells = list(reactor.cells)
    for cell in cells:
        _refresh_cell_residuals(cell)

    items = {
        "elemental_closure": audit_elemental_closure(cells),
        "mass_closure": audit_mass_closure(cells),
        "energy_closure": audit_energy_closure(cells),
        "inner_nr": audit_inner_nr(reactor, include_jacobian=include_jacobian),
        "outer_convergence": audit_outer_convergence(cells),
        "freeboard_closure": audit_freeboard_closure(
            reactor.config,
            n_freeboard_cells=int(getattr(reactor.config, "n_freeboard_cells", 0) or 0),
        ),
        "composition_observation": audit_composition_observation(cells, p0_exit=p0_exit),
        "temperature_observation": audit_temperature_observation(cells, p0_exit=p0_exit),
    }
    assert tuple(items) == CHECKLIST_ORDER
    open_items = [
        name
        for name, payload in items.items()
        if payload.get("open") is True
    ]
    return {
        "source": "Wirsum alignment P5 LU re-entry",
        "case": "reduced_3cell_lu_init",
        "solved": False,
        "checklist_order": list(CHECKLIST_ORDER),
        "open_items": open_items,
        "items": items,
        "rule": (
            "Do not treat CO vs plant data as proof that a numerical bridge is correct. "
            "Close elemental/mass/energy and diagnose inner/outer/freeboard first."
        ),
    }


def build_p5_lu_reentry_payload(*, include_jacobian: bool = False) -> dict[str, Any]:
    from tests.global_nr_fixtures import build_initialized_lu_reactor

    reactor = build_initialized_lu_reactor(n_cells=3)
    return audit_lu_reentry(reactor, include_jacobian=include_jacobian)
