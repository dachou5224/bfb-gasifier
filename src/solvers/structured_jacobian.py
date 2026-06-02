"""Structured Jacobian helpers for thesis/freeboard global NR.

This module keeps the solver-graph view and the structured linear algebra
separate from ``global_nr_solver.py`` so we can evolve the implementation
without disturbing the Newton control flow.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Callable

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from src.core.cell import Cell
from src.solvers.nr_indexing import build_var_index_to_cell_map, cell_offsets, unpack_reactor

_BLOCK_NNZ_TOL = 1e-14


@dataclass(frozen=True)
class SolverGraphBlock:
    block_id: str
    cell_index: int
    chain_pos: int
    kind: str


@dataclass(frozen=True)
class SolverGraph:
    blocks: tuple[SolverGraphBlock, ...]
    main_chain: tuple[int, ...]
    n_bed: int
    n_freeboard: int
    has_cyclone: bool
    has_return_leg: bool
    tail_start: int | None
    head_stop: int | None


@dataclass(frozen=True)
class JacobianStructure:
    offsets: tuple[int, ...]
    band_block_pairs: frozenset[tuple[int, int]]
    side_block_pairs: frozenset[tuple[int, int]]
    main_chain: tuple[int, ...]
    side_tail_blocks: tuple[int, ...]
    side_head_blocks: tuple[int, ...]


@dataclass(frozen=True)
class StructuredJacobian:
    matrix: sp.csr_matrix
    band_matrix: sp.csr_matrix
    structure: JacobianStructure
    side_blocks: tuple[tuple[int, int, np.ndarray], ...]
    unexpected_blocks: tuple[tuple[int, int], ...]
    residual_cell_calls: int
    zero_cols: int
    zero_rows: int
    row_support_by_col: tuple[np.ndarray, ...] | None = None

    @property
    def structure_validation_ok(self) -> bool:
        return len(self.unexpected_blocks) == 0

    def to_uv(self) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(U, V)`` so that ``matrix = band_matrix + U @ V.T`` on known side blocks."""
        n_total = self.matrix.shape[0]
        if not self.side_blocks:
            return np.zeros((n_total, 0), dtype=np.float64), np.zeros((n_total, 0), dtype=np.float64)

        offsets = self.structure.offsets
        u_parts: list[np.ndarray] = []
        v_parts: list[np.ndarray] = []
        for row_block, col_block, block in self.side_blocks:
            row_start, row_end = offsets[row_block], offsets[row_block + 1]
            col_start, col_end = offsets[col_block], offsets[col_block + 1]
            if block.size == 0 or not np.any(np.abs(block) > _BLOCK_NNZ_TOL):
                continue
            n_cols = int(col_end - col_start)
            u_local = np.zeros((n_total, n_cols), dtype=np.float64)
            v_local = np.zeros((n_total, n_cols), dtype=np.float64)
            u_local[row_start:row_end, :] = block
            v_local[col_start:col_end, :] = np.eye(n_cols, dtype=np.float64)
            u_parts.append(u_local)
            v_parts.append(v_local)
        if not u_parts:
            return np.zeros((n_total, 0), dtype=np.float64), np.zeros((n_total, 0), dtype=np.float64)
        return np.concatenate(u_parts, axis=1), np.concatenate(v_parts, axis=1)


def build_solver_graph(cells: list[Cell]) -> SolverGraph:
    """Infer the linear solver graph from the ordered NR cell list."""
    blocks: list[SolverGraphBlock] = []
    n_bed = 0
    n_freeboard = 0
    has_cyclone = False
    has_return_leg = False

    for idx, cell in enumerate(cells):
        kind = str(getattr(cell, "cell_type", "bed")).strip().lower()
        if kind == "bed":
            block_id = f"bed_cell_{n_bed}"
            n_bed += 1
        elif kind == "freeboard":
            block_id = f"freeboard_cell_{n_freeboard}"
            n_freeboard += 1
        elif kind == "cyclone":
            block_id = "cyclone_block"
            has_cyclone = True
        elif kind == "return_leg":
            block_id = "return_leg_block"
            has_return_leg = True
        else:
            block_id = f"{kind}_cell_{idx}"
        blocks.append(
            SolverGraphBlock(
                block_id=block_id,
                cell_index=idx,
                chain_pos=idx,
                kind=kind,
            )
        )

    tail_start: int | None = None
    head_stop: int | None = None
    if has_return_leg and n_bed > 0:
        head_stop = n_bed
        tail_start = max(n_bed - 1, 0)

    return SolverGraph(
        blocks=tuple(blocks),
        main_chain=tuple(range(len(blocks))),
        n_bed=n_bed,
        n_freeboard=n_freeboard,
        has_cyclone=has_cyclone,
        has_return_leg=has_return_leg,
        tail_start=tail_start,
        head_stop=head_stop,
    )


def build_jacobian_structure(graph: SolverGraph, cells: list[Cell]) -> JacobianStructure:
    """Build the block sparsity pattern from the inferred solver graph."""
    offsets = tuple(cell_offsets(cells))
    n_blocks = len(graph.blocks)

    band_pairs: set[tuple[int, int]] = set()
    for i in range(n_blocks):
        for j in (i - 1, i, i + 1):
            if 0 <= j < n_blocks:
                band_pairs.add((i, j))

    side_pairs: set[tuple[int, int]] = set()
    side_tail_blocks: tuple[int, ...] = tuple()
    side_head_blocks: tuple[int, ...] = tuple()
    if graph.n_freeboard == 0 and graph.n_bed > 1:
        # Bed-only thesis/validation cases can still have a top-cell -> bottom-cell
        # recycle boundary, and bottom-cell transport can enter upstream energy
        # rows through axial solid enthalpy.  Keep those nonlocal closures as
        # explicit side blocks so the default Jacobian need not fall back to dense
        # finite differences.
        top_bed_idx = graph.n_bed - 1
        if abs(top_bed_idx - 0) > 1:
            side_pairs.add((0, top_bed_idx))
            for row_idx in range(graph.n_bed):
                if abs(row_idx - 0) > 1:
                    side_pairs.add((row_idx, 0))
            side_head_blocks = (0,)
            side_tail_blocks = (top_bed_idx,)
    if graph.tail_start is not None and graph.head_stop is not None:
        source_tail_idx = graph.n_bed + graph.n_freeboard - 1 if graph.n_freeboard > 0 else graph.n_bed - 1
        cyclone_idx = next((b.cell_index for b in graph.blocks if b.kind == "cyclone"), None)
        return_leg_idx = next((b.cell_index for b in graph.blocks if b.kind == "return_leg"), None)

        if graph.n_freeboard == 0:
            # Hamel side-element core path for bed + cyclone + return-leg graph:
            # recycle nonlocal couplings are limited to the physical closure chain
            # instead of a dense head×tail closure.
            side_head_blocks = (0,) if graph.n_bed > 0 else tuple()
            explicit_tail: list[int] = [source_tail_idx]
            if cyclone_idx is not None:
                explicit_tail.append(int(cyclone_idx))
            if return_leg_idx is not None:
                explicit_tail.append(int(return_leg_idx))
            side_tail_blocks = tuple(
                dict.fromkeys(int(idx) for idx in explicit_tail if 0 <= int(idx) < n_blocks)
            )

            if source_tail_idx >= 0 and side_head_blocks and abs(side_head_blocks[0] - source_tail_idx) > 1:
                side_pairs.add((side_head_blocks[0], source_tail_idx))
            if side_head_blocks and cyclone_idx is not None and abs(side_head_blocks[0] - cyclone_idx) > 1:
                side_pairs.add((side_head_blocks[0], int(cyclone_idx)))
            if (
                source_tail_idx >= 0
                and return_leg_idx is not None
                and abs(int(return_leg_idx) - source_tail_idx) > 1
            ):
                side_pairs.add((int(return_leg_idx), source_tail_idx))
        else:
            # Explicit freeboard-in-graph path still keeps the conservative closure.
            side_head_blocks = tuple(range(graph.head_stop))
            explicit_tail = list(range(graph.tail_start, n_blocks))
            side_tail_blocks = tuple(dict.fromkeys(int(idx) for idx in explicit_tail))

            # Explicit recycle side-chain couplings:
            # top source -> cyclone -> return leg -> bed_0 -> bed_1
            if source_tail_idx >= 0 and cyclone_idx is not None and abs(cyclone_idx - source_tail_idx) > 1:
                side_pairs.add((cyclone_idx, source_tail_idx))
            if source_tail_idx >= 0 and return_leg_idx is not None and abs(return_leg_idx - source_tail_idx) > 1:
                side_pairs.add((return_leg_idx, source_tail_idx))
            if cyclone_idx is not None and return_leg_idx is not None and abs(return_leg_idx - cyclone_idx) > 1:
                side_pairs.add((return_leg_idx, cyclone_idx))

            top_bed_idx = max(graph.n_bed - 1, 0)
            first_freeboard_idx = graph.n_bed if graph.n_freeboard > 0 else None
            last_freeboard_idx = graph.n_bed + graph.n_freeboard - 1 if graph.n_freeboard > 0 else None
            if graph.n_bed > 1:
                side_pairs.add((0, top_bed_idx))
            if last_freeboard_idx is not None and abs(top_bed_idx - last_freeboard_idx) > 1:
                side_pairs.add((top_bed_idx, last_freeboard_idx))
            for col_idx in range(n_blocks):
                affected_rows: set[int] = {col_idx}
                if col_idx < graph.n_bed:
                    if col_idx > 0:
                        affected_rows.add(col_idx - 1)
                    if col_idx + 1 < graph.n_bed:
                        affected_rows.add(col_idx + 1)
                    if col_idx == 0:
                        affected_rows.update(range(graph.n_bed))
                    if col_idx == top_bed_idx and first_freeboard_idx is not None:
                        affected_rows.add(first_freeboard_idx)
                elif first_freeboard_idx is not None and col_idx < graph.n_bed + graph.n_freeboard:
                    fb_idx = col_idx - graph.n_bed
                    if fb_idx > 0:
                        affected_rows.add(col_idx - 1)
                    if fb_idx + 1 < graph.n_freeboard:
                        affected_rows.add(col_idx + 1)
                    if fb_idx == 0:
                        affected_rows.update(range(graph.n_bed))
                    if col_idx == last_freeboard_idx and cyclone_idx is not None:
                        affected_rows.add(int(cyclone_idx))
                        if graph.n_bed > 0:
                            affected_rows.add(0)
                        if graph.n_bed > 1:
                            affected_rows.add(1)
                elif cyclone_idx is not None and col_idx == int(cyclone_idx):
                    if last_freeboard_idx is not None:
                        affected_rows.add(last_freeboard_idx)
                    if return_leg_idx is not None:
                        affected_rows.add(int(return_leg_idx))
                    if graph.n_bed > 0:
                        affected_rows.add(0)
                        affected_rows.add(top_bed_idx)
                    if graph.n_bed > 1:
                        affected_rows.add(1)
                elif return_leg_idx is not None and col_idx == int(return_leg_idx):
                    if graph.n_bed > 0:
                        affected_rows.add(0)
                        affected_rows.add(top_bed_idx)
                    if graph.n_bed > 1:
                        affected_rows.add(1)

                for row_idx in affected_rows:
                    if 0 <= row_idx < n_blocks and abs(row_idx - col_idx) > 1:
                        side_pairs.add((row_idx, col_idx))

    return JacobianStructure(
        offsets=offsets,
        band_block_pairs=frozenset(band_pairs),
        side_block_pairs=frozenset(side_pairs),
        main_chain=graph.main_chain,
        side_tail_blocks=side_tail_blocks,
        side_head_blocks=side_head_blocks,
    )


def structured_strategy_name_from_cells(cells: list[Cell]) -> str:
    graph = build_solver_graph(cells)
    if graph.has_return_leg:
        return "band_plus_side_elements_structured"
    return "block_tridiag_structured"


def _row_support_by_col_from_matrix(matrix: sp.csr_matrix) -> tuple[np.ndarray, ...]:
    csc = matrix.tocsc()
    return tuple(
        csc.indices[csc.indptr[j] : csc.indptr[j + 1]].astype(np.intp, copy=True)
        for j in range(matrix.shape[1])
    )


def assemble_jacobian_fd_structured(
    x: np.ndarray,
    F0: np.ndarray,
    cells: list[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
    *,
    structure: JacobianStructure,
    apply_local_bc_fn: Callable[[int], None] | None = None,
    affected_residual_cells_fn: Callable[[int], tuple[int, ...]] | None = None,
    verbose: bool = False,
) -> StructuredJacobian:
    """Assemble a structured FD Jacobian directly on the declared block graph."""
    from src.solvers.global_nr_solver import _project_cell_residual_to_main_nr

    n_total = len(x)
    offsets = structure.offsets
    var_to_cell = build_var_index_to_cell_map(list(offsets))
    n_blocks = len(cells)
    affected_row_blocks_by_col: dict[int, tuple[int, ...]] = {}
    for col_block in range(n_blocks):
        allowed_rows = set(
            row_block
            for row_block in range(n_blocks)
            if (row_block, col_block) in structure.band_block_pairs
            or (row_block, col_block) in structure.side_block_pairs
        )
        if affected_residual_cells_fn is not None:
            local_rows = set(int(i) for i in affected_residual_cells_fn(int(col_block)))
            restricted_rows = sorted(allowed_rows.intersection(local_rows))
            # Keep a safe fallback for unexpected callback output.
            affected_row_blocks_by_col[col_block] = tuple(restricted_rows if restricted_rows else sorted(allowed_rows))
        else:
            affected_row_blocks_by_col[col_block] = tuple(sorted(allowed_rows))

    rows_all: list[int] = []
    cols_all: list[int] = []
    vals_all: list[float] = []
    band_rows: list[int] = []
    band_cols: list[int] = []
    band_vals: list[float] = []
    side_block_arrays: dict[tuple[int, int], np.ndarray] = {}
    residual_cell_calls = 0
    unexpected_blocks: set[tuple[int, int]] = set()

    unpack_reactor(x, cells)
    apply_bc_fn()
    F0 = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
    F0s = F0 / eq_scale
    if verbose and len(x) > 100:
        print(f"DEBUG: Building structured Jacobian for {n_total} variables...")

    x_p = x.copy()
    for j in range(n_total):
        if verbose and j > 0 and j % 100 == 0:
            print(f"DEBUG:   Structured Jacobian progress: {j}/{n_total}")

        col_block = int(var_to_cell[j])
        xj = float(x[j])
        h = max(1e-6 * max(abs(xj), 1.0), 1e-8)
        x_p[j] = xj
        x_p[j] += h

        unpack_reactor(x_p, cells)
        if apply_local_bc_fn is not None:
            apply_local_bc_fn(col_block)
        else:
            apply_bc_fn()
        affected_row_blocks = affected_row_blocks_by_col[col_block]
        if len(affected_row_blocks) == n_blocks:
            F_p_all = np.concatenate([_project_cell_residual_to_main_nr(c) for c in cells])
            residual_cell_calls += len(cells)
        else:
            F_p_all = None

        for row_block in affected_row_blocks:
            cell = cells[row_block]
            row_start = offsets[row_block]
            row_end = offsets[row_block + 1]
            if F_p_all is None:
                F_local = _project_cell_residual_to_main_nr(cell)
                residual_cell_calls += 1
            else:
                F_local = F_p_all[row_start:row_end]
            dF_local = (F_local / eq_scale[row_start:row_end] - F0s[row_start:row_end]) / h
            nz = np.flatnonzero(np.abs(dF_local) > _BLOCK_NNZ_TOL)
            if nz.size == 0:
                continue
            rows_all.extend((row_start + nz).tolist())
            cols_all.extend([j] * int(nz.size))
            vals_all.extend(dF_local[nz].astype(float).tolist())
            if (row_block, col_block) in structure.band_block_pairs:
                band_rows.extend((row_start + nz).tolist())
                band_cols.extend([j] * int(nz.size))
                band_vals.extend(dF_local[nz].astype(float).tolist())
                continue
            if (row_block, col_block) in structure.side_block_pairs:
                block_key = (row_block, col_block)
                block_arr = side_block_arrays.get(block_key)
                if block_arr is None:
                    block_arr = np.zeros((row_end - row_start, offsets[col_block + 1] - offsets[col_block]), dtype=np.float64)
                    side_block_arrays[block_key] = block_arr
                block_arr[nz, j - offsets[col_block]] = dF_local[nz]
                continue
            unexpected_blocks.add((row_block, col_block))

        x_p[j] = xj

    unpack_reactor(x, cells)
    apply_bc_fn()

    J = sp.csr_matrix((vals_all, (rows_all, cols_all)), shape=(n_total, n_total))
    J_band = sp.csr_matrix((band_vals, (band_rows, band_cols)), shape=(n_total, n_total))
    side_blocks = tuple(
        (row_block, col_block, block)
        for (row_block, col_block), block in sorted(side_block_arrays.items())
        if np.any(np.abs(block) > _BLOCK_NNZ_TOL)
    )
    zero_cols = int(np.count_nonzero(np.diff(J.tocsc().indptr) == 0))
    zero_rows = int(np.count_nonzero(np.diff(J.indptr) == 0))

    return StructuredJacobian(
        matrix=J,
        band_matrix=J_band,
        structure=structure,
        side_blocks=side_blocks,
        unexpected_blocks=tuple(sorted(unexpected_blocks)),
        residual_cell_calls=int(residual_cell_calls),
        zero_cols=zero_cols,
        zero_rows=zero_rows,
        row_support_by_col=_row_support_by_col_from_matrix(J),
    )


def solve_structured_linear_step(
    structured: StructuredJacobian,
    rhs: np.ndarray,
    *,
    regularize: bool = False,
) -> tuple[np.ndarray, str, bool, dict[str, float | int | bool]]:
    """Solve ``(J_band + side) dx = rhs`` using Woodbury-style side correction."""
    timing = {
        "band_lu_s": 0.0,
        "side_update_s": 0.0,
        "schur_size": 0,
        "fallback_used": False,
    }
    if not structured.structure_validation_ok:
        dx, method, regularized = solve_sparse_direct_linear_step(
            structured.matrix,
            rhs,
            regularize=regularize,
        )
        timing["fallback_used"] = True
        return dx, "sparse_direct_fallback", regularized, timing

    t0 = perf_counter()
    A = structured.band_matrix
    if regularize:
        A = A + sp.eye(A.shape[0], format="csr") * 1e-8
    try:
        lu = splu(A.tocsc())
    except Exception:
        dx, method, regularized = solve_sparse_direct_linear_step(
            structured.matrix,
            rhs,
            regularize=True,
        )
        timing["fallback_used"] = True
        return dx, "sparse_direct_fallback", True or regularized, timing
    timing["band_lu_s"] = perf_counter() - t0

    U, V = structured.to_uv()
    if U.shape[1] == 0:
        return lu.solve(rhs), "structured_direct", regularize, timing

    t1 = perf_counter()
    y = lu.solve(rhs)
    Y = lu.solve(U)
    schur = np.eye(U.shape[1], dtype=np.float64) + V.T @ Y
    timing["schur_size"] = int(schur.shape[0])
    try:
        z = np.linalg.solve(schur, V.T @ y)
    except Exception:
        z, _, _, _ = np.linalg.lstsq(schur, V.T @ y, rcond=None)
    timing["side_update_s"] = perf_counter() - t1
    return y - Y @ z, "structured_direct", regularize, timing


def _tikhonov_svd_damped(
    J_dense: np.ndarray,
    rhs: np.ndarray,
    t_dof_indices: np.ndarray,
    t_max_step_K: float,
) -> tuple[np.ndarray, float]:
    """Return (dx, mu) from scalar Tikhonov regularization chosen so max|dx_T| ≤ t_max_step_K.

    Uses a single SVD of J_dense and binary-searches for the regularisation
    parameter mu.  This is purely arithmetic (no additional matrix factorisation)
    and adds negligible overhead compared with the FD Jacobian build.

    The scalar Tikhonov filter is:  dx = Vh^T * (σ/(σ²+μ)) * U^T rhs
    For μ → 0  : recovers minimum-norm lstsq.
    For μ → ∞  : dx → 0 (fully suppressed).
    """
    U, s, Vh = np.linalg.svd(J_dense, full_matrices=False)
    UTF = U.T @ rhs

    def _dx_for_mu(mu: float) -> np.ndarray:
        s_filt = s / (s**2 + mu)
        return Vh.T @ (s_filt * UTF)

    # Quick check: is the unregularised solution already within limit?
    dx0 = _dx_for_mu(0.0)
    if float(np.max(np.abs(dx0[t_dof_indices]))) <= t_max_step_K:
        return dx0, 0.0

    # Binary search: find smallest mu such that max|dx_T| ≤ t_max_step_K.
    lo, hi = 0.0, float(np.max(s) ** 2 * 1e8)  # upper bound kills all modes
    for _ in range(60):
        mid = (lo + hi) * 0.5
        max_T = float(np.max(np.abs(_dx_for_mu(mid)[t_dof_indices])))
        if max_T > t_max_step_K:
            lo = mid
        else:
            hi = mid
        if hi - lo < lo * 1e-6:
            break
    dx_reg = _dx_for_mu(hi)
    return dx_reg, hi


def solve_sparse_direct_linear_step(
    J: sp.csr_matrix,
    rhs: np.ndarray,
    *,
    regularize: bool = False,
    large_step_threshold: float = 5000.0,
    t_dof_indices: np.ndarray | None = None,
    t_max_step_K: float = 600.0,
) -> tuple[np.ndarray, str, bool]:
    """Sparse direct solve with minimum-norm lstsq fallback for near-singular Jacobians.

    When the Jacobian is nearly singular (holdup_transport + frozen Vorabrechnung creates a
    near-degenerate temperature mode), splu amplifies the near-null direction and produces
    non-physical Newton steps (|ΔT| >> 1e5 K, |Δgas| >> 1e3 mol/s).  The fallback uses
    minimum-norm least-squares (rcond=1e-12) to find the minimum-norm solution.

    If t_dof_indices is provided and the resulting max|ΔT| still exceeds t_max_step_K,
    a scalar Tikhonov regularisation is applied via SVD binary search.  This damps the
    near-null temperature modes so that the accepted Newton step stays within a physically
    plausible temperature range (default 600 K per cell), while barely affecting the
    well-conditioned gas/solid modes (σ ≫ √μ).

    large_step_threshold
        If max(|dx|) after splu exceeds this value, minimum-norm lstsq is used instead.
        Default 5000 separates the holdup_transport blow-up (~2e5 K) from the legacy_stream
        path (~886 K max ΔT, where splu is fine).
    t_dof_indices
        Indices of temperature DOFs in the solution vector.  Required for Tikhonov fallback.
    t_max_step_K
        Maximum allowed raw |ΔT| per cell before Tikhonov regularisation kicks in.
    """
    J_solve = J
    if regularize:
        J_solve = J + sp.eye(J.shape[0], format="csr") * 1e-8
    J_dense: np.ndarray | None = None
    try:
        lu = splu(J_solve.tocsc())
        dx = lu.solve(rhs)
        if float(np.max(np.abs(dx))) > large_step_threshold:
            # Near-singular Jacobian: use minimum-norm lstsq.
            J_dense = J_solve.toarray()
            dx, _, _, _ = np.linalg.lstsq(J_dense, rhs, rcond=1e-12)
            # Tikhonov fallback: if T steps are still too large, apply SVD damping.
            if (
                t_dof_indices is not None
                and len(t_dof_indices) > 0
                and float(np.max(np.abs(dx[t_dof_indices]))) > t_max_step_K
            ):
                dx, mu = _tikhonov_svd_damped(J_dense, rhs, t_dof_indices, t_max_step_K)
                return dx, f"lstsq_tikhonov(mu={mu:.2e})", regularize
            return dx, "lstsq", regularize
        return dx, "splu", regularize
    except Exception:
        if J_dense is None:
            J_dense = J_solve.toarray()
        dx, _, _, _ = np.linalg.lstsq(J_dense, rhs, rcond=None)
        return dx, "lstsq", regularize
