from __future__ import annotations

from src.core.connectivity_graph import affected_nr_residual_cells


def test_affected_cells_simple_bed_chain():
    assert affected_nr_residual_cells(
        changed_cell_idx=2,
        n_bed=4,
        n_freeboard=0,
        explicit_freeboard_graph=False,
        side_blocks_in_boundary_path=False,
        has_side_blocks=False,
    ) == (1, 2, 3)


def test_affected_cells_bed_with_side_blocks():
    assert affected_nr_residual_cells(
        changed_cell_idx=3,
        n_bed=4,
        n_freeboard=0,
        explicit_freeboard_graph=False,
        side_blocks_in_boundary_path=True,
        has_side_blocks=True,
    ) == (2, 3, 4)


def test_affected_cells_explicit_freeboard_graph_with_recycle():
    assert affected_nr_residual_cells(
        changed_cell_idx=5,
        n_bed=3,
        n_freeboard=1,
        explicit_freeboard_graph=True,
        side_blocks_in_boundary_path=True,
        has_side_blocks=True,
    ) == (0, 1, 2, 4, 5)
