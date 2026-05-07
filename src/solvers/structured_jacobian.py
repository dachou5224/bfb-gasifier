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
    if graph.tail_start is not None and graph.head_stop is not None:
        # Wraparound recycle couplings do not stop at bed_0: once the tail-side
        # blocks (bed top / freeboard / cyclone / return leg) feed back to the
        # lower bed, the downstream bed chain induces a small block-closure over
        # all bed cells. We therefore keep the physical side chain explicit, then
        # close it over the bed head/tail partition instead of using a generic
        # dense Jacobian.
        side_head_blocks = tuple(range(graph.head_stop))
        source_tail_idx = graph.n_bed + graph.n_freeboard - 1 if graph.n_freeboard > 0 else graph.n_bed - 1
        cyclone_idx = next((b.cell_index for b in graph.blocks if b.kind == "cyclone"), None)
        return_leg_idx = next((b.cell_index for b in graph.blocks if b.kind == "return_leg"), None)

        explicit_tail: list[int] = list(range(graph.tail_start, n_blocks))
        side_tail_blocks = tuple(dict.fromkeys(int(idx) for idx in explicit_tail))

        # Explicit recycle side-chain couplings:
        # top source -> cyclone -> return leg -> bed_0 -> bed_1
        if source_tail_idx >= 0 and cyclone_idx is not None and abs(cyclone_idx - source_tail_idx) > 1:
            side_pairs.add((cyclone_idx, source_tail_idx))
        if source_tail_idx >= 0 and return_leg_idx is not None and abs(return_leg_idx - source_tail_idx) > 1:
            side_pairs.add((return_leg_idx, source_tail_idx))
        if cyclone_idx is not None and return_leg_idx is not None and abs(return_leg_idx - cyclone_idx) > 1:
            side_pairs.add((return_leg_idx, cyclone_idx))

        for head_idx in side_head_blocks:
            for tail_idx in side_tail_blocks:
                if abs(head_idx - tail_idx) <= 1:
                    continue
                side_pairs.add((head_idx, tail_idx))
                side_pairs.add((tail_idx, head_idx))

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


def assemble_jacobian_fd_structured(
    x: np.ndarray,
    F0: np.ndarray,
    cells: list[Cell],
    apply_bc_fn: Callable[[], None],
    eq_scale: np.ndarray,
    *,
    structure: JacobianStructure,
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
        allowed_rows = sorted(
            row_block
            for row_block in range(n_blocks)
            if (row_block, col_block) in structure.band_block_pairs
            or (row_block, col_block) in structure.side_block_pairs
        )
        affected_row_blocks_by_col[col_block] = tuple(allowed_rows)

    rows_all: list[int] = []
    cols_all: list[int] = []
    vals_all: list[float] = []
    band_rows: list[int] = []
    band_cols: list[int] = []
    band_vals: list[float] = []
    side_block_arrays: dict[tuple[int, int], np.ndarray] = {}
    residual_cell_calls = 0
    unexpected_blocks: set[tuple[int, int]] = set()

    F0s = F0 / eq_scale
    if verbose and len(x) > 100:
        print(f"DEBUG: Building structured Jacobian for {n_total} variables...")

    x_p = x.copy()
    for j in range(n_total):
        if verbose and j > 0 and j % 100 == 0:
            print(f"DEBUG:   Structured Jacobian progress: {j}/{n_total}")

        xj = float(x[j])
        h = max(1e-6 * max(abs(xj), 1.0), 1e-8)
        x_p[j] = xj
        x_p[j] += h

        unpack_reactor(x_p, cells)
        apply_bc_fn()
        col_block = int(var_to_cell[j])
        affected_row_blocks = affected_row_blocks_by_col[col_block]
        residual_cell_calls += len(affected_row_blocks)

        for row_block in affected_row_blocks:
            cell = cells[row_block]
            row_start = offsets[row_block]
            row_end = offsets[row_block + 1]
            F_local = _project_cell_residual_to_main_nr(cell)
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


def solve_sparse_direct_linear_step(
    J: sp.csr_matrix,
    rhs: np.ndarray,
    *,
    regularize: bool = False,
) -> tuple[np.ndarray, str, bool]:
    """Sparse direct solve with least-squares fallback."""
    J_solve = J
    if regularize:
        J_solve = J + sp.eye(J.shape[0], format="csr") * 1e-8
    try:
        lu = splu(J_solve.tocsc())
        return lu.solve(rhs), "splu", regularize
    except Exception:
        dx, _, _, _ = np.linalg.lstsq(J_solve.toarray(), rhs, rcond=None)
        return dx, "lstsq", regularize
