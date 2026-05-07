"""单 cell 守恒方程求解器（scipy.optimize.fsolve / root / least_squares）。

给定上游条件（入口摩尔流率、进料），求解当前 cell 的出口状态
（N_b, N_d, m_solid, T），使得所有守恒方程残差为零。

Source: docs/CLAUDE.md Phase 5.2
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from src.core.cell import Cell, N_SOLID_COMP
from src.core.species import N_GAS, GAS_SPECIES_INDEX


def _build_residual_scales(cell: Cell) -> np.ndarray:
    """构造单 cell 残差缩放尺度。"""
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    scales = np.ones(2 * N_GAS + nv_sol + 1, dtype=np.float64)

    gas_ref_d = float(np.sum(np.maximum(cell.N_zu_d + cell.N_rez_d + cell.N_d_in, 0.0)))
    gas_ref_b = float(np.sum(np.maximum(cell.N_zu_b + cell.N_rez_b + cell.N_b_in, 0.0)))
    solid_ref_state = cell.m_solid_zu + cell.m_solid_rez + cell.m_solid_in
    if str(getattr(cell, "solid_state_model", "legacy_stream")) in {"holdup_transport", "freeboard_closure"}:
        solid_ref_state = solid_ref_state + cell.m_solid_auf_in + cell.m_solid_ab_in
    solid_ref = float(np.sum(np.maximum(solid_ref_state, 0.0)))

    H_in = (
        cell._calc_gas_enthalpy_flow(cell.N_b_in + cell.N_d_in, cell.T_in_gas)
        + cell._calc_gas_enthalpy_flow(cell.N_zu_b + cell.N_zu_d, cell.T_zu_gas)
        + cell._calc_gas_enthalpy_flow(cell.N_rez_b + cell.N_rez_d, cell.T_rez_gas)
        + cell._calc_solid_enthalpy_flow(cell.m_solid_rez, cell.T_rez_solid)
        + cell._calc_solid_enthalpy_flow(cell.m_solid_in, cell.T_in_solid)
        + cell._calc_solid_enthalpy_flow(cell.m_solid_zu, cell.T_zu_solid)
    )

    scales[:N_GAS] = max(gas_ref_d, 1e-6)
    scales[N_GAS:2 * N_GAS] = max(gas_ref_b, 1e-6)
    scales[2 * N_GAS:2 * N_GAS + nv_sol] = max(solid_ref, 1e-8)
    scales[-1] = max(abs(H_in), 1e5)
    return scales

def _pack_state(cell: Cell) -> np.ndarray:
    """将 cell 的自由变量打包为一维向量。"""
    return np.concatenate([
        cell.N_d,
        cell.N_b,
        cell.m_solid.flatten(),
        np.array([cell.T]),
    ])

def _unpack_state(x: np.ndarray, cell: Cell) -> None:
    """将一维向量解包回 cell 状态。"""
    nk = cell.solid.n_size_classes
    nv_sol = nk * N_SOLID_COMP
    # 数值防护：强制非负并处理 NaN
    x_safe = np.nan_to_num(x, nan=0.0, posinf=1e10, neginf=0.0)
    
    np.maximum(x_safe[:N_GAS], 0.0, out=cell.N_d)
    np.maximum(x_safe[N_GAS:2 * N_GAS], 0.0, out=cell.N_b)
    m_flat = np.maximum(x_safe[2 * N_GAS:2 * N_GAS + nv_sol], 0.0)
    cell.m_solid[:] = m_flat.reshape((nk, N_SOLID_COMP))
    cell.T = float(np.clip(x_safe[2 * N_GAS + nv_sol], 300.0, 3000.0))

def evaluate_cell_state(
    cell: Cell,
    scales: np.ndarray | None = None,
    rate_multiplier: float = 1.0,
) -> dict[str, np.ndarray | float]:
    """Evaluate residual metrics for the cell's current primitive state."""
    scales_local = _build_residual_scales(cell) if scales is None else scales
    residual_vector = cell.residuals(rate_multiplier=rate_multiplier)
    max_res = float(np.max(np.abs(residual_vector)))
    rms_scaled = float(np.sqrt(np.mean((residual_vector / scales_local) ** 2)))
    return {
        "residual_vector": residual_vector,
        "residual": max_res,
        "rms_scaled": rms_scaled,
    }


def _residual_wrapper(x: np.ndarray, cell: Cell, rate_multiplier: float = 1.0) -> np.ndarray:
    """求解器的目标函数。"""
    _unpack_state(x, cell)
    return cell.residuals(rate_multiplier=rate_multiplier)


def _residual_wrapper_scaled(x: np.ndarray, cell: Cell, scales: np.ndarray) -> np.ndarray:
    """按方程尺度归一后的残差。"""
    _unpack_state(x, cell)
    return cell.residuals() / scales


def _prefer_candidate_metrics(
    *,
    incumbent_rms: float,
    incumbent_res: float,
    candidate_rms: float,
    candidate_res: float,
    residual_tol_scaled: float,
    incumbent_success: bool = True,
    candidate_success: bool = True,
    allow_near_tie: bool = False,
    near_tie_margin: float = 1.05,
) -> bool:
    """Accept new branches only if scaled and physical metrics are aligned."""
    inc_rms = float(incumbent_rms) if np.isfinite(incumbent_rms) else np.inf
    inc_res = float(abs(incumbent_res)) if np.isfinite(incumbent_res) else np.inf
    cand_rms = float(candidate_rms) if np.isfinite(candidate_rms) else np.inf
    cand_res = float(abs(candidate_res)) if np.isfinite(candidate_res) else np.inf
    if not np.isfinite(cand_rms) or not np.isfinite(cand_res):
        return False

    clear_rms_improve = cand_rms < 0.98 * inc_rms
    clear_res_improve = cand_res < 0.98 * inc_res
    res_not_worse = cand_res <= 1.05 * inc_res

    if cand_rms <= residual_tol_scaled:
        return res_not_worse or clear_res_improve
    if clear_rms_improve and (res_not_worse or clear_res_improve):
        return True
    if (not incumbent_success) and candidate_success and clear_res_improve:
        return True
    if allow_near_tie and cand_rms <= near_tie_margin * inc_rms and clear_res_improve:
        return True
    return False

def _homotopy_warmup(cell: Cell, n_steps: int = 7) -> None:
    """Source-term homotopy warm-up for stiff bottom cells.

    Gradually ramps the reaction rate multiplier from a small fraction to 1.0
    via ``n_steps`` explicit balance passes.  Each pass uses a mini
    ``least_squares`` solve at the current multiplier level to find a
    consistent state, then hands that state to the next level.  This prevents
    least_squares from jumping to the cold/zero-combustion branch that is
    often the nearest local minimum from a cold initial guess.

    The cell state is updated in-place; the caller then runs the full primary
    solve from this warm-started state.
    """
    import numpy as np  # already imported at top of module

    lower_bounds = np.zeros(len(_pack_state(cell)))
    lower_bounds[-1] = 300.0
    upper_bounds = np.full(len(_pack_state(cell)), np.inf)
    upper_bounds[-1] = 3000.0

    # Bias more points at low-rate end to reduce first-step shock in bottom cells.
    if n_steps >= 7:
        multipliers = np.array([0.05, 0.10, 0.20, 0.35, 0.55, 0.75, 1.00], dtype=float)
    else:
        multipliers = np.linspace(1.0 / n_steps, 1.0, n_steps)
    for rm in multipliers:
        x0 = _pack_state(cell)

        def _res_homotopy(x, _cell=cell, _rm=rm):
            return _residual_wrapper(x, _cell, rate_multiplier=float(_rm))

        try:
            res = least_squares(
                _res_homotopy,
                x0,
                bounds=(lower_bounds, upper_bounds),
                jac="2-point",
                loss="soft_l1",
                f_scale=1.0,
                ftol=1e-4,
                xtol=1e-4,
                gtol=1e-4,
                max_nfev=120,
            )
            x_out = np.clip(res.x, lower_bounds, upper_bounds)
            _unpack_state(x_out, cell)
        except Exception:
            # If a homotopy step fails, keep current state and continue.
            pass
    # Restore full-rate reactions so the primary solve uses correct residuals.
    cell.calc_reactions(rate_multiplier=1.0)


def solve_cell(
    cell: Cell,
    max_iter: int = 50,
    tol: float = 1e-6,
    residual_tol_scaled: float = 0.1,
    stiff_stabilization: bool = False,
    skip_homotopy: bool = False,
    skip_multistart: bool = False,
    temperature_cap_K: float | None = None,
    verbose: bool = False,
) -> dict:
    """稳健地求解单个 cell 的守恒方程（有界 least_squares）。

    针对强非线性/刚性状态，可开启 ``stiff_stabilization=True`` 的两段式策略：
    1) 先解缩放残差（平衡气/固/能量方程量级）；
    2) 再解原始残差（收敛到真实守恒方程）。

    Parameters
    ----------
    temperature_cap_K:
        可选温度上界 [K]。用于 GS 底部 stiff continuation，将 full solve
        约束在 Vorabrechnung 参考窗口附近，避免从同一 warm-start 直接跳到
        明显过热的局部分支。
    """
    x0 = _pack_state(cell)
    x0_raw = x0.copy()
    n_vars = len(x0)

    lower_bounds = np.zeros(n_vars)
    lower_bounds[-1] = 300.0
    upper_bounds = np.full(n_vars, np.inf)
    upper_bounds[-1] = 3000.0
    if temperature_cap_K is not None:
        upper_bounds[-1] = float(
            np.clip(float(temperature_cap_K), lower_bounds[-1] + 1e-6, upper_bounds[-1])
        )

    # Minor-species guardrails:
    # H2S/NH3 are weakly constrained in the current reduced mechanism and can
    # numerically blow up, contaminating dry-basis major-gas fractions.
    # Limit each phase by local available inflow + modest cushion.
    idx_h2s = GAS_SPECIES_INDEX["H2S"]
    idx_nh3 = GAS_SPECIES_INDEX["NH3"]
    h2s_in = float(
        max(cell.N_d_in[idx_h2s], 0.0)
        + max(cell.N_b_in[idx_h2s], 0.0)
        + max(cell.N_zu_d[idx_h2s], 0.0)
        + max(cell.N_zu_b[idx_h2s], 0.0)
        + max(cell.N_rez_d[idx_h2s], 0.0)
        + max(cell.N_rez_b[idx_h2s], 0.0)
    )
    nh3_in = float(
        max(cell.N_d_in[idx_nh3], 0.0)
        + max(cell.N_b_in[idx_nh3], 0.0)
        + max(cell.N_zu_d[idx_nh3], 0.0)
        + max(cell.N_zu_b[idx_nh3], 0.0)
        + max(cell.N_rez_d[idx_nh3], 0.0)
        + max(cell.N_rez_b[idx_nh3], 0.0)
    )
    minor_pad = 0.5  # [mol/s] per cell, allows release but prevents runaway
    h2s_cap = max(
        h2s_in + minor_pad,
        float(max(cell.N_d[idx_h2s], 0.0) + max(cell.N_b[idx_h2s], 0.0)),
        1e-6,
    )
    nh3_cap = max(
        nh3_in + minor_pad,
        float(max(cell.N_d[idx_nh3], 0.0) + max(cell.N_b[idx_nh3], 0.0)),
        1e-6,
    )
    upper_bounds[idx_h2s] = h2s_cap
    upper_bounds[N_GAS + idx_h2s] = h2s_cap
    upper_bounds[idx_nh3] = nh3_cap
    upper_bounds[N_GAS + idx_nh3] = nh3_cap

    # Ensure initial guess is feasible under current bounds (especially
    # after minor-species caps are applied), otherwise least_squares raises
    # "x0 is infeasible".
    # Also enforce bounds consistency so lb <= ub always holds.
    upper_bounds = np.maximum(upper_bounds, lower_bounds + 1e-12)
    x0 = np.clip(x0, lower_bounds, upper_bounds)
    x0_raw = np.clip(x0_raw, lower_bounds, upper_bounds)  # probe fallback safety

    diag = np.ones(n_vars)
    diag[:2 * N_GAS] = 10.0
    diag[-1] = 1000.0
    diag = np.maximum(diag, 1e-6)

    scales = _build_residual_scales(cell)
    max_nfev = max(80, min(max_iter * 8, 240))

    # Source homotopy warm-up: ramp reaction rates 0→1 to avoid cold-branch local minima.
    # Activated via stiff_stabilization=True (same flag, same intent — solver robustness).
    # skip_homotopy=True disables the warmup in GS context to save ~840 nfev per cell.
    if stiff_stabilization and not skip_homotopy:
        _homotopy_warmup(cell, n_steps=7)
        # Re-clip after homotopy: the warmup uses loose [0, inf] bounds and may
        # push minor species (H2S, NH3) above the caps set above.
        x0 = np.clip(_pack_state(cell), lower_bounds, upper_bounds)

    # Primary solve: keep legacy behavior as default trajectory.
    res_primary = least_squares(
        _residual_wrapper,
        x0,
        args=(cell,),
        bounds=(lower_bounds, upper_bounds),
        x_scale=diag,
        jac="2-point",
        loss="linear",
        ftol=tol,
        xtol=tol,
        gtol=tol,
        max_nfev=max_nfev,
        verbose=2 if verbose else 0,
    )
    x_best = np.clip(res_primary.x, lower_bounds, upper_bounds)
    _unpack_state(x_best, cell)
    final_res_best = cell.residuals()
    max_res_best = float(np.max(np.abs(final_res_best)))
    rms_scaled_best = float(np.sqrt(np.mean((final_res_best / scales) ** 2)))
    res_best = res_primary

    # Stabilization fallback:
    # For stiff-mode diagnostics we should not wait for catastrophic failure.
    # As long as primary solve is not physically converged, attempt scaled-stage
    # correction from the current best state.
    need_stabilize = stiff_stabilization and (
        (not bool(res_primary.success))
        or (not np.isfinite(rms_scaled_best))
        or (rms_scaled_best >= residual_tol_scaled)
    )
    stiff_attempted = False
    stiff_accepted = False
    if need_stabilize:
        stiff_attempted = True
        res_scaled = least_squares(
            _residual_wrapper_scaled,
            x_best,
            args=(cell, scales),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="soft_l1",
            f_scale=1.0,
            ftol=tol,
            xtol=tol,
            gtol=tol,
            max_nfev=max(100, max_nfev),
            verbose=2 if verbose else 0,
        )

        x_stab0 = np.clip(res_scaled.x, lower_bounds, upper_bounds)
        res_stab = least_squares(
            _residual_wrapper,
            x_stab0,
            args=(cell,),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="linear",
            ftol=tol,
            xtol=tol,
            gtol=tol,
            max_nfev=max(120, max_nfev),
            verbose=2 if verbose else 0,
        )
        x_stab = np.clip(res_stab.x, lower_bounds, upper_bounds)
        _unpack_state(x_stab, cell)
        final_res_stab = cell.residuals()
        max_res_stab = float(np.max(np.abs(final_res_stab)))
        rms_scaled_stab = float(np.sqrt(np.mean((final_res_stab / scales) ** 2)))

        # Optional polishing pass: tighten around stabilized point to reduce
        # scaled RMS for physical convergence gate.
        x_polish0 = x_stab
        res_polish = least_squares(
            _residual_wrapper,
            x_polish0,
            args=(cell,),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="linear",
            ftol=max(tol * 0.2, 1e-9),
            xtol=max(tol * 0.2, 1e-9),
            gtol=max(tol * 0.2, 1e-9),
            max_nfev=max(80, max_nfev // 2),
            verbose=2 if verbose else 0,
        )
        x_polish = np.clip(res_polish.x, lower_bounds, upper_bounds)
        _unpack_state(x_polish, cell)
        final_res_polish = cell.residuals()
        max_res_polish = float(np.max(np.abs(final_res_polish)))
        rms_scaled_polish = float(np.sqrt(np.mean((final_res_polish / scales) ** 2)))

        if np.isfinite(rms_scaled_polish) and (
            (rms_scaled_polish < rms_scaled_stab)
            or (max_res_polish < max_res_stab)
        ):
            x_stab = x_polish
            final_res_stab = final_res_polish
            max_res_stab = max_res_polish
            rms_scaled_stab = rms_scaled_polish
            res_stab = res_polish

        use_stab = _prefer_candidate_metrics(
            incumbent_rms=rms_scaled_best,
            incumbent_res=max_res_best,
            candidate_rms=rms_scaled_stab,
            candidate_res=max_res_stab,
            residual_tol_scaled=residual_tol_scaled,
            incumbent_success=bool(res_primary.success),
            candidate_success=bool(res_stab.success),
        )
        if use_stab:
            x_best = x_stab
            final_res_best = final_res_stab
            max_res_best = max_res_stab
            rms_scaled_best = rms_scaled_stab
            res_best = res_stab
            stiff_accepted = True

    # Probe fallback: if stiff path still fails physical threshold, retry from
    # the original pre-homotopy initial guess and keep the better solution.
    if stiff_stabilization and (
        (not np.isfinite(rms_scaled_best))
        or (rms_scaled_best >= residual_tol_scaled)
    ):
        res_probe = least_squares(
            _residual_wrapper,
            x0_raw,
            args=(cell,),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="linear",
            ftol=tol,
            xtol=tol,
            gtol=tol,
            max_nfev=max(80, max_nfev),
            verbose=2 if verbose else 0,
        )
        x_probe = np.clip(res_probe.x, lower_bounds, upper_bounds)
        _unpack_state(x_probe, cell)
        final_res_probe = cell.residuals()
        max_res_probe = float(np.max(np.abs(final_res_probe)))
        rms_scaled_probe = float(np.sqrt(np.mean((final_res_probe / scales) ** 2)))

        if _prefer_candidate_metrics(
            incumbent_rms=rms_scaled_best,
            incumbent_res=max_res_best,
            candidate_rms=rms_scaled_probe,
            candidate_res=max_res_probe,
            residual_tol_scaled=residual_tol_scaled,
            incumbent_success=bool(res_best.success),
            candidate_success=bool(res_probe.success),
        ):
            x_best = x_probe
            final_res_best = final_res_probe
            max_res_best = max_res_probe
            rms_scaled_best = rms_scaled_probe
            res_best = res_probe

    # Final micro-refinement for near-threshold cases.
    if stiff_stabilization and np.isfinite(rms_scaled_best) and rms_scaled_best >= residual_tol_scaled:
        res_refine_scaled = least_squares(
            _residual_wrapper_scaled,
            x_best,
            args=(cell, scales),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="soft_l1",
            f_scale=0.5,
            ftol=max(tol * 0.5, 1e-9),
            xtol=max(tol * 0.5, 1e-9),
            gtol=max(tol * 0.5, 1e-9),
            max_nfev=80,
            verbose=2 if verbose else 0,
        )
        x_ref0 = np.clip(res_refine_scaled.x, lower_bounds, upper_bounds)
        res_refine = least_squares(
            _residual_wrapper,
            x_ref0,
            args=(cell,),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="linear",
            ftol=max(tol * 0.5, 1e-9),
            xtol=max(tol * 0.5, 1e-9),
            gtol=max(tol * 0.5, 1e-9),
            max_nfev=80,
            verbose=2 if verbose else 0,
        )
        x_ref = np.clip(res_refine.x, lower_bounds, upper_bounds)
        _unpack_state(x_ref, cell)
        final_res_ref = cell.residuals()
        max_res_ref = float(np.max(np.abs(final_res_ref)))
        rms_scaled_ref = float(np.sqrt(np.mean((final_res_ref / scales) ** 2)))
        if _prefer_candidate_metrics(
            incumbent_rms=rms_scaled_best,
            incumbent_res=max_res_best,
            candidate_rms=rms_scaled_ref,
            candidate_res=max_res_ref,
            residual_tol_scaled=residual_tol_scaled,
            incumbent_success=bool(res_best.success),
            candidate_success=bool(res_refine.success),
        ):
            x_best = x_ref
            final_res_best = final_res_ref
            max_res_best = max_res_ref
            rms_scaled_best = rms_scaled_ref
            res_best = res_refine

    # Temperature-seeded multi-start (branch recovery):
    # Bottom stiff cells can be sensitive to initial temperature and get trapped
    # in a high-residual branch. Seeds must start from a chemically populated
    # state (`x_best`), not the raw inlet state (`x0_raw`), because the raw
    # state often has many zero species and an ill-conditioned Jacobian.
    # skip_multistart=True disables this section to save time in GS context.
    if (
        stiff_stabilization
        and not skip_multistart
        and np.isfinite(rms_scaled_best)
        and rms_scaled_best >= residual_tol_scaled
    ):
        t_seed_best = float(np.clip(x_best[-1], 300.0, 3000.0))
        h_center = getattr(cell.geo, "h_center", 1.0)
        dh = max(float(getattr(cell.geo, "dh", 0.1)), 1e-6)
        in_bottom_band = h_center <= 1.5 * dh
        in_lower_bed_band = h_center <= 3.5 * dh
        T_MIN_PHYS = 950.0
        if in_bottom_band:
            T_floor = 800.0
            T_cap = 1800.0
        elif in_lower_bed_band:
            T_floor = 900.0
            T_cap = 2200.0
        else:
            T_floor = T_MIN_PHYS
            T_cap = 2600.0
        T_cap = min(T_cap, float(upper_bounds[-1]))

        # Reduced seed candidates to limit per-cell cost in GS context.
        # Use 2 relative seeds (±15%) + 2 absolute seeds for stiff bottom cells.
        seed_candidates = [
            (x_best, 0.85 * t_seed_best),
            (x_best, 1.15 * t_seed_best),
        ]
        if in_bottom_band:
            seed_candidates.extend([
                (x_best, 950.0),
                (x_best, 1000.0),
                (x_best, 1100.0),
                (x_best, 1200.0),
            ])
        elif in_lower_bed_band:
            seed_candidates.extend([
                (x_best, 950.0),
                (x_best, 1000.0),
                (x_best, 1100.0),
                (x_best, 1200.0),
            ])

        t_candidates_with_base = [
            (base, t) for base, t in seed_candidates if T_floor <= t <= upper_bounds[-1]
        ]
        for base_vec, t0 in t_candidates_with_base:
            x_try0 = np.array(base_vec, dtype=float)
            x_try0[-1] = float(np.clip(t0, lower_bounds[-1], upper_bounds[-1]))

            # First pass: scaled robust solve with automatic Jacobian-based
            # variable scaling to handle ill-conditioned bottom-cell states.
            res_try_scaled = least_squares(
                _residual_wrapper_scaled,
                x_try0,
                args=(cell, scales),
                bounds=(lower_bounds, upper_bounds),
                x_scale="jac",
                jac="2-point",
                loss="soft_l1",
                f_scale=0.5,
                ftol=max(tol * 0.5, 1e-9),
                xtol=max(tol * 0.5, 1e-9),
                gtol=max(tol * 0.5, 1e-9),
                max_nfev=80,  # reduced budget per seed to keep GS iteration cost bounded
                verbose=2 if verbose else 0,
            )
            x_try1 = np.clip(res_try_scaled.x, lower_bounds, upper_bounds)

            # Second pass: unscaled polish onto the true residual system.
            res_try = least_squares(
                _residual_wrapper,
                x_try1,
                args=(cell,),
                bounds=(lower_bounds, upper_bounds),
                x_scale="jac",
                jac="2-point",
                loss="linear",
                ftol=max(tol * 0.2, 1e-10),
                xtol=max(tol * 0.2, 1e-10),
                gtol=max(tol * 0.2, 1e-10),
                max_nfev=120,  # reduced budget per seed (second pass)
                verbose=2 if verbose else 0,
            )
            x_try = np.clip(res_try.x, lower_bounds, upper_bounds)
            _unpack_state(x_try, cell)
            final_res_try = cell.residuals()
            max_res_try = float(np.max(np.abs(final_res_try)))
            rms_scaled_try = float(np.sqrt(np.mean((final_res_try / scales) ** 2)))
            t_final_try = float(x_try[-1])
            temp_phys_ok = (T_floor <= t_final_try <= T_cap)
            t_best_current = float(x_best[-1])
            best_temp_bad = t_best_current < T_floor

            if temp_phys_ok and _prefer_candidate_metrics(
                incumbent_rms=rms_scaled_best,
                incumbent_res=max_res_best,
                candidate_rms=rms_scaled_try,
                candidate_res=max_res_try,
                residual_tol_scaled=residual_tol_scaled,
                incumbent_success=bool(res_best.success),
                candidate_success=bool(res_try.success),
                allow_near_tie=best_temp_bad,
            ):
                x_best = x_try
                final_res_best = final_res_try
                max_res_best = max_res_try
                rms_scaled_best = rms_scaled_try
                res_best = res_try

    # Final polish: if still not converged after all seeds, do one last
    # high-budget solve from the overall best point found.
    # Only active when skip_multistart=False (i.e., not in the fast GS context).
    if stiff_stabilization and not skip_multistart and np.isfinite(rms_scaled_best) and rms_scaled_best >= residual_tol_scaled:
        res_polish = least_squares(
            _residual_wrapper,
            x_best,
            args=(cell,),
            bounds=(lower_bounds, upper_bounds),
            x_scale=diag,
            jac="2-point",
            loss="linear",
            ftol=max(tol * 0.2, 1e-10),
            xtol=max(tol * 0.2, 1e-10),
            gtol=max(tol * 0.2, 1e-10),
            max_nfev=300,
            verbose=2 if verbose else 0,
        )
        x_pol = np.clip(res_polish.x, lower_bounds, upper_bounds)
        _unpack_state(x_pol, cell)
        final_res_pol = cell.residuals()
        rms_pol = float(np.sqrt(np.mean((final_res_pol / scales) ** 2)))
        if _prefer_candidate_metrics(
            incumbent_rms=rms_scaled_best,
            incumbent_res=max_res_best,
            candidate_rms=rms_pol,
            candidate_res=float(np.max(np.abs(final_res_pol))),
            residual_tol_scaled=residual_tol_scaled,
            incumbent_success=bool(res_best.success),
            candidate_success=bool(res_polish.success),
        ):
            x_best = x_pol
            final_res_best = final_res_pol
            max_res_best = float(np.max(np.abs(final_res_pol)))
            rms_scaled_best = rms_pol
            res_best = res_polish

    _unpack_state(x_best, cell)
    optimizer_success = bool(res_best.success)
    physically_converged = (
        optimizer_success
        and np.isfinite(rms_scaled_best)
        and rms_scaled_best <= residual_tol_scaled
    )

    return {
        "converged": optimizer_success,
        "optimizer_success": optimizer_success,
        "physically_converged": physically_converged,
        "residual": max_res_best,
        "rms_scaled": rms_scaled_best,
        "stiff_attempted": stiff_attempted,
        "stiff_accepted": stiff_accepted,
        "message": str(res_best.message),
    }
