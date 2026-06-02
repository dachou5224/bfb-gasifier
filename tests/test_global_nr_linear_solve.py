from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from src.solvers.global_nr_solver import _solve_linear_step
from src.solvers.structured_jacobian import (
    JacobianStructure,
    StructuredJacobian,
    solve_structured_linear_step,
)

def test_solve_linear_step_regularizes_rank_deficient_jacobian():
    J = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 0.0]], dtype=float))
    rhs = np.array([-1.0, 0.0], dtype=float)
    dx, method, regularized, _ = _solve_linear_step(J, rhs, zero_rows=1, zero_cols=0)
    assert method in {"splu", "lstsq"}
    assert regularized is True
    assert np.isfinite(dx).all()
    assert abs(float(dx[0]) + 1.0) < 1e-6


def test_solve_linear_step_identity_without_regularization():
    J = sp.csr_matrix(np.eye(3, dtype=float))
    rhs = np.array([1.0, -2.0, 0.5], dtype=float)
    dx, method, regularized, diag = _solve_linear_step(J, rhs, zero_rows=0, zero_cols=0)
    assert method == "splu"
    assert regularized is False
    assert np.allclose(dx, rhs)
    assert diag["schur_size"] == 0


def test_structured_linear_step_matches_sparse_direct_for_band_plus_one_side_block():
    band = sp.csr_matrix(
        np.array(
            [
                [4.0, 1.0, 0.0, 0.0],
                [1.0, 3.5, 1.0, 0.0],
                [0.0, 1.0, 3.0, 1.0],
                [0.0, 0.0, 1.0, 2.5],
            ],
            dtype=float,
        )
    )
    side_block = np.array([[0.25]], dtype=float)
    full = band.toarray()
    full[0:1, 3:4] += side_block
    structured = StructuredJacobian(
        matrix=sp.csr_matrix(full),
        band_matrix=band,
        structure=JacobianStructure(
            offsets=(0, 1, 2, 3, 4),
            band_block_pairs=frozenset({(0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 2), (3, 3)}),
            side_block_pairs=frozenset({(0, 3)}),
            main_chain=(0, 1, 2, 3),
            side_tail_blocks=(3,),
            side_head_blocks=(0,),
        ),
        side_blocks=((0, 3, side_block),),
        unexpected_blocks=(),
        residual_cell_calls=0,
        zero_cols=0,
        zero_rows=0,
    )
    rhs = np.array([1.0, -0.5, 0.75, 2.0], dtype=float)
    dx_struct, method, regularized, diag = solve_structured_linear_step(structured, rhs, regularize=False)
    dx_ref = np.linalg.solve(full, rhs)
    assert method == "structured_direct"
    assert regularized is False
    assert diag["schur_size"] == 1
    assert np.allclose(dx_struct, dx_ref)
