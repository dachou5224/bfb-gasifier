"""NR residual graph neighborhood helpers.

Keep adjacency/topology logic separate from ``Reactor`` orchestration so the
Hamel-style explicit graph can be audited independently.
"""

from __future__ import annotations


def affected_nr_residual_cells(
    *,
    changed_cell_idx: int,
    n_bed: int,
    n_freeboard: int,
    explicit_freeboard_graph: bool,
    side_blocks_in_boundary_path: bool,
    has_side_blocks: bool,
) -> tuple[int, ...]:
    """Return affected residual cell indices for one changed solver cell."""
    idx = int(changed_cell_idx)

    if explicit_freeboard_graph:
        n_total = int(n_bed) + int(n_freeboard) + (2 if has_side_blocks else 0)
        if n_total <= 0 or idx < 0 or idx >= n_total:
            return tuple()

        cyclone_idx = (n_bed + n_freeboard) if has_side_blocks else None
        return_leg_idx = (n_bed + n_freeboard + 1) if has_side_blocks else None
        top_bed_idx = max(int(n_bed) - 1, 0)
        first_freeboard_idx = int(n_bed) if int(n_freeboard) > 0 else None
        last_freeboard_idx = (int(n_bed) + int(n_freeboard) - 1) if int(n_freeboard) > 0 else None

        neighbors: set[int] = {idx}
        if idx < int(n_bed):
            if idx > 0:
                neighbors.add(idx - 1)
            if idx + 1 < int(n_bed):
                neighbors.add(idx + 1)
            if idx == 0:
                neighbors.update(range(int(n_bed)))
            if idx == top_bed_idx:
                if first_freeboard_idx is not None:
                    neighbors.add(first_freeboard_idx)
                elif cyclone_idx is not None:
                    neighbors.add(cyclone_idx)
            return tuple(sorted(neighbors))

        if first_freeboard_idx is not None and idx < int(n_bed) + int(n_freeboard):
            fb_idx = idx - int(n_bed)
            if fb_idx > 0:
                neighbors.add(idx - 1)
            if fb_idx + 1 < int(n_freeboard):
                neighbors.add(idx + 1)
            if fb_idx == 0:
                neighbors.update(range(int(n_bed)))
            if idx == last_freeboard_idx and cyclone_idx is not None:
                neighbors.add(cyclone_idx)
                if int(n_bed) > 0:
                    neighbors.add(0)
                if int(n_bed) > 1:
                    neighbors.add(1)
            return tuple(sorted(neighbors))

        if cyclone_idx is not None and idx == cyclone_idx:
            if last_freeboard_idx is not None:
                neighbors.add(last_freeboard_idx)
            elif int(n_bed) > 0:
                neighbors.add(top_bed_idx)
            if return_leg_idx is not None:
                neighbors.add(return_leg_idx)
            if int(n_bed) > 0:
                neighbors.add(0)
                neighbors.add(top_bed_idx)
            if int(n_bed) > 1:
                neighbors.add(1)
            return tuple(sorted(neighbors))

        if return_leg_idx is not None and idx == return_leg_idx:
            if int(n_bed) > 0:
                neighbors.add(0)
                neighbors.add(top_bed_idx)
            if int(n_bed) > 1:
                neighbors.add(1)

        return tuple(sorted(neighbors))

    if not side_blocks_in_boundary_path:
        if n_bed <= 0:
            return tuple()
        affected: list[int] = []
        if idx == 0:
            affected.extend(range(n_bed))
        if idx > 0:
            affected.append(idx - 1)
        if 0 <= idx < n_bed:
            affected.append(idx)
        if idx + 1 < n_bed:
            affected.append(idx + 1)
        if idx == n_bed - 1 and 0 not in affected:
            affected.append(0)
        return tuple(sorted(set(affected)))

    cyclone_idx = n_bed
    return_leg_idx = n_bed + 1
    affected: set[int] = {idx}
    if idx < n_bed:
        if idx > 0:
            affected.add(idx - 1)
        if idx + 1 < n_bed:
            affected.add(idx + 1)
        if idx == n_bed - 1:
            affected.add(cyclone_idx)
    elif idx == cyclone_idx:
        affected.add(return_leg_idx)
    elif idx == return_leg_idx:
        affected.add(0)
    return tuple(sorted(affected))
