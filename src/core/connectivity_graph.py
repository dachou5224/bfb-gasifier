"""Hamel-style Verbindungsmatrix schema and graph helpers.

Source: Hamel (2001) Kap 2.1–2.3; ``docs/hamel_b_tier_thesis_extraction.md``.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class ConnectivityBlock:
    """One bilanzierten Bilanzraum or auxiliary feed/exit node."""

    block_id: str
    kind: str
    axial_index: int | None = None
    xi_reactor: float | None = None
    solver_coupling: str = "explicit"
    solver_cell_index: int | None = None


@dataclass(frozen=True)
class ConnectivityEdge:
    """Directed stream link between blocks (Verbindungsmatrix row semantics)."""

    from_block: str
    to_block: str
    kind: str
    direction: str
    solver_coupling: str
    phase_scope: str | None = None
    mechanism: str | None = None
    stream_class: str | None = None  # zu | rez | auf | ab | boundary | postprocess


@dataclass
class ConnectivityGraph:
    """Frozen connectivity topology for thesis / B-tier NR assembly."""

    mode: str
    blocks: tuple[ConnectivityBlock, ...]
    edges: tuple[ConnectivityEdge, ...]
    freeboard_active: bool = False
    recycle_active: bool = False
    secondary_injection_active: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    def block_ids(self) -> frozenset[str]:
        return frozenset(block.block_id for block in self.blocks)

    def block_by_id(self) -> dict[str, ConnectivityBlock]:
        return {block.block_id: block for block in self.blocks}

    def inbound_edges(self, block_id: str) -> tuple[ConnectivityEdge, ...]:
        return tuple(edge for edge in self.edges if edge.to_block == block_id)

    def outbound_edges(self, block_id: str) -> tuple[ConnectivityEdge, ...]:
        return tuple(edge for edge in self.edges if edge.from_block == block_id)

    def cross_link_edges(self) -> tuple[ConnectivityEdge, ...]:
        """Edges with Hamel ``stream_class == rez`` or inferred non-adjacent coupling."""
        out: list[ConnectivityEdge] = []
        for edge in self.edges:
            sc = edge.stream_class or infer_stream_class(edge.kind, edge.direction)
            if sc == "rez":
                out.append(edge)
        return tuple(out)

    @property
    def n_cross_links(self) -> int:
        """``NCL``: count of non-adjacent / recirculation couplings."""
        return len(self.cross_link_edges())

    def solver_cell_block_ids(self) -> tuple[str, ...]:
        """Block ids that map to NR solver cells (explicit_state coupling)."""
        return tuple(
            block.block_id
            for block in self.blocks
            if block.solver_cell_index is not None
        )

    def validate(self) -> None:
        ids = self.block_ids()
        if len(ids) != len(self.blocks):
            raise ValueError("duplicate block_id in connectivity graph")
        for edge in self.edges:
            if edge.from_block not in ids:
                raise ValueError(f"edge from unknown block: {edge.from_block!r}")
            if edge.to_block not in ids:
                raise ValueError(f"edge to unknown block: {edge.to_block!r}")

    def to_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "blocks": [asdict(block) for block in self.blocks],
            "edges": [asdict(edge) for edge in self.edges],
            "counts": {
                "n_blocks": len(self.blocks),
                "n_edges": len(self.edges),
                "n_bed_blocks": sum(1 for b in self.blocks if b.kind == "bed_cell"),
                "n_freeboard_blocks": sum(1 for b in self.blocks if b.kind == "freeboard_cell"),
                "n_cross_links": self.n_cross_links,
            },
            "freeboard_active": bool(self.freeboard_active),
            "recycle_active": bool(self.recycle_active),
            "secondary_injection_active": bool(self.secondary_injection_active),
            "metadata": dict(self.metadata),
        }


def infer_stream_class(kind: str, direction: str) -> str:
    """Map edge kind to Hamel zu/rez/auf/ab semantics."""
    k = str(kind).strip().lower()
    if k in {"side_feed", "primary_gas_feed", "fuel_feed"} or direction == "into_reactor":
        return "zu"
    if k in {"external_zirkulation", "solid_separation"} and direction in {
        "downward_to_bottom",
        "downward",
    }:
        return "rez"
    if k in {"gas_exit"} or direction == "out_of_system":
        return "postprocess"
    if k in {"internal_zirkulation"} or direction == "downward":
        return "ab"
    if direction in {"upward", "upward_to_top"}:
        return "auf"
    if direction in {"into_reactor", "downward_to_bottom"}:
        return "boundary"
    return "auf"


def _edge_from_dict(raw: Mapping[str, Any]) -> ConnectivityEdge:
    kind = str(raw["kind"])
    direction = str(raw["direction"])
    stream_class = raw.get("stream_class")
    if stream_class is None:
        stream_class = infer_stream_class(kind, direction)
    return ConnectivityEdge(
        from_block=str(raw["from_block"]),
        to_block=str(raw["to_block"]),
        kind=kind,
        direction=direction,
        solver_coupling=str(raw["solver_coupling"]),
        phase_scope=raw.get("phase_scope"),
        mechanism=raw.get("mechanism"),
        stream_class=str(stream_class) if stream_class is not None else None,
    )


def _block_from_dict(raw: Mapping[str, Any]) -> ConnectivityBlock:
    return ConnectivityBlock(
        block_id=str(raw["block_id"]),
        kind=str(raw["kind"]),
        axial_index=raw.get("axial_index"),
        xi_reactor=raw.get("xi_reactor"),
        solver_coupling=str(raw.get("solver_coupling", "explicit")),
        solver_cell_index=raw.get("solver_cell_index"),
    )


def connectivity_graph_from_dict(payload: Mapping[str, Any]) -> ConnectivityGraph:
    """Parse a connectivity graph from JSON/YAML-compatible dict."""
    blocks = tuple(_block_from_dict(item) for item in payload.get("blocks", ()))
    edges = tuple(_edge_from_dict(item) for item in payload.get("edges", ()))
    graph = ConnectivityGraph(
        mode=str(payload.get("mode", "hamel_connectivity")),
        blocks=blocks,
        edges=edges,
        freeboard_active=bool(payload.get("freeboard_active", False)),
        recycle_active=bool(payload.get("recycle_active", False)),
        secondary_injection_active=bool(payload.get("secondary_injection_active", False)),
        metadata=dict(payload.get("metadata", {})),
    )
    graph.validate()
    return graph


def load_connectivity_graph_json(path: str | Path) -> ConnectivityGraph:
    """Load ``ConnectivityGraph`` from a JSON file."""
    text = Path(path).read_text(encoding="utf-8")
    return connectivity_graph_from_dict(json.loads(text))


def side_block_id(kind: str, index: int, *, count: int) -> str:
    """Return canonical side-block id (legacy singular name when ``count == 1``)."""
    if count <= 1:
        return "cyclone_block" if kind == "cyclone" else "return_leg_block"
    return f"{'cyclone_block' if kind == 'cyclone' else 'return_leg_block'}_{int(index)}"


def ordered_solver_block_ids(graph: ConnectivityGraph, kind: str) -> tuple[str, ...]:
    """Return solver-coupled block ids of one kind sorted by ``solver_cell_index``."""
    blocks = [
        block
        for block in graph.blocks
        if block.kind == kind and block.solver_cell_index is not None
    ]
    blocks.sort(key=lambda block: int(block.solver_cell_index))
    return tuple(block.block_id for block in blocks)


def build_connectivity_graph_from_reactor_params(
    *,
    n_bed_cells: int,
    n_freeboard_cells: int,
    H_bed: float,
    H_freeboard: float,
    recycle_gas: bool,
    recirculation_frac: float,
    cyclone_capture_char_frac: float,
    cyclone_capture_ash_frac: float,
    secondary_injection_applied: bool,
    secondary_injection_segment: int | None,
    explicit_side_blocks: bool = False,
    preset: str = "hamel_minimal",
    include_connection_pipe: bool = False,
    n_cyclone_stages: int = 1,
    n_return_legs: int = 1,
    return_leg_bed_targets: tuple[int, ...] | None = None,
) -> ConnectivityGraph:
    """Build the default Hamel-style graph from reactor segment counts.

    Ref: Hamel (2001) Bild 2.2 Konfiguration; Kap 7 VTT cyclone loop.
    """
    total_height = max(float(H_bed) + max(float(H_freeboard), 0.0), 1e-12)
    freeboard_active = int(max(n_freeboard_cells, 0)) > 0 and float(H_freeboard) > 0.0
    solver_idx = 0

    blocks: list[ConnectivityBlock] = [
        ConnectivityBlock("primary_gas_feed", "side_feed", solver_coupling="boundary_condition"),
        ConnectivityBlock("fuel_feed", "side_feed", solver_coupling="boundary_condition"),
    ]
    edges: list[ConnectivityEdge] = [
        ConnectivityEdge(
            "primary_gas_feed",
            "bed_cell_0",
            "side_feed",
            "into_reactor",
            "boundary_condition",
            stream_class="zu",
        ),
        ConnectivityEdge(
            "fuel_feed",
            "bed_cell_0",
            "side_feed",
            "into_reactor",
            "boundary_condition",
            stream_class="zu",
        ),
    ]

    for i in range(int(max(n_bed_cells, 0))):
        xi = ((i + 0.5) * float(H_bed) / max(float(n_bed_cells), 1.0)) / total_height
        blocks.append(
            ConnectivityBlock(
                block_id=f"bed_cell_{i}",
                kind="bed_cell",
                axial_index=i,
                xi_reactor=float(xi),
                solver_coupling="explicit_state",
                solver_cell_index=solver_idx,
            )
        )
        solver_idx += 1
        if i > 0:
            edges.extend(
                [
                    ConnectivityEdge(
                        f"bed_cell_{i - 1}",
                        f"bed_cell_{i}",
                        "axial_upflow",
                        "upward",
                        "explicit_state",
                        stream_class="auf",
                    ),
                    ConnectivityEdge(
                        f"bed_cell_{i - 1}",
                        f"bed_cell_{i}",
                        "wake_solid_upflow",
                        "upward",
                        "explicit_state",
                        phase_scope="suspension_phase_only",
                        mechanism="bubble_wake_carrying",
                        stream_class="auf",
                    ),
                    ConnectivityEdge(
                        f"bed_cell_{i}",
                        f"bed_cell_{i - 1}",
                        "internal_zirkulation",
                        "downward",
                        "explicit_state",
                        phase_scope="suspension_phase_only",
                        mechanism="mass_continuity_backmixing",
                        stream_class="ab",
                    ),
                ]
            )

    top_source_block = f"bed_cell_{max(int(n_bed_cells) - 1, 0)}"
    if freeboard_active:
        interface_xi = float(H_bed / total_height)
        blocks.append(
            ConnectivityBlock(
                "bed_freeboard_interface",
                "interface_block",
                axial_index=int(max(n_bed_cells, 0)),
                xi_reactor=interface_xi,
                solver_coupling="explicit_handoff",
            )
        )
        edges.append(
            ConnectivityEdge(
                top_source_block,
                "bed_freeboard_interface",
                "bubble_rupture_handoff",
                "upward",
                "explicit_handoff",
                stream_class="auf",
            )
        )
        for i in range(int(max(n_freeboard_cells, 0))):
            xi = (
                float(H_bed) + (i + 0.5) * float(H_freeboard) / max(float(n_freeboard_cells), 1.0)
            ) / total_height
            blocks.append(
                ConnectivityBlock(
                    block_id=f"freeboard_cell_{i}",
                    kind="freeboard_cell",
                    axial_index=i,
                    xi_reactor=float(xi),
                    solver_coupling="explicit_state",
                    solver_cell_index=solver_idx,
                )
            )
            solver_idx += 1
            if i == 0:
                edges.append(
                    ConnectivityEdge(
                        "bed_freeboard_interface",
                        "freeboard_cell_0",
                        "freeboard_entry",
                        "upward",
                        "explicit_state",
                        stream_class="auf",
                    )
                )
            else:
                edges.append(
                    ConnectivityEdge(
                        f"freeboard_cell_{i - 1}",
                        f"freeboard_cell_{i}",
                        "axial_upflow",
                        "upward",
                        "explicit_state",
                        stream_class="auf",
                    )
                )
        top_source_block = f"freeboard_cell_{max(int(n_freeboard_cells) - 1, 0)}"

    if secondary_injection_applied and secondary_injection_segment is not None:
        target_block = (
            f"freeboard_cell_{int(secondary_injection_segment)}"
            if freeboard_active
            else top_source_block
        )
        blocks.append(
            ConnectivityBlock(
                "secondary_inlet_block",
                "side_feed",
                axial_index=int(secondary_injection_segment),
                solver_coupling="explicit_side_feed",
            )
        )
        edges.append(
            ConnectivityEdge(
                "secondary_inlet_block",
                target_block,
                "side_feed",
                "into_reactor",
                "explicit_side_feed",
                stream_class="zu",
            )
        )

    cyclone_coupling = "explicit_state" if explicit_side_blocks else "post_freeboard_separation"
    return_coupling = "explicit_state" if explicit_side_blocks else "boundary_condition"
    n_cyc = max(1, int(n_cyclone_stages))
    n_ret = max(1, int(n_return_legs))
    cyclone_ids = tuple(side_block_id("cyclone", i, count=n_cyc) for i in range(n_cyc))
    return_ids = tuple(side_block_id("return_leg", i, count=n_ret) for i in range(n_ret))
    if return_leg_bed_targets is None:
        bed_targets = tuple(0 for _ in range(n_ret))
    else:
        bed_targets = tuple(int(v) for v in return_leg_bed_targets)
        if len(bed_targets) != n_ret:
            raise ValueError(
                f"return_leg_bed_targets length ({len(bed_targets)}) must equal n_return_legs ({n_ret})"
            )

    blocks.append(ConnectivityBlock("system_exit", "system_exit", solver_coupling="postprocess"))
    for i, cyclone_id in enumerate(cyclone_ids):
        blocks.append(
            ConnectivityBlock(
                cyclone_id,
                "cyclone",
                axial_index=i,
                solver_coupling=cyclone_coupling,
                solver_cell_index=solver_idx if explicit_side_blocks else None,
            )
        )
        if explicit_side_blocks:
            solver_idx += 1
    for i, return_id in enumerate(return_ids):
        blocks.append(
            ConnectivityBlock(
                return_id,
                "return_leg",
                axial_index=i,
                solver_coupling=return_coupling,
                solver_cell_index=solver_idx if explicit_side_blocks else None,
            )
        )
        if explicit_side_blocks:
            solver_idx += 1

    edges.append(
        ConnectivityEdge(
            top_source_block,
            cyclone_ids[0],
            "cyclone_inlet",
            "upward",
            cyclone_coupling,
            stream_class="auf",
        )
    )
    for prev_id, next_id in zip(cyclone_ids[:-1], cyclone_ids[1:]):
        edges.append(
            ConnectivityEdge(
                prev_id,
                next_id,
                "cyclone_stage",
                "downward",
                cyclone_coupling,
                stream_class="auf",
            )
        )
    edges.append(
        ConnectivityEdge(
            cyclone_ids[-1],
            "system_exit",
            "gas_exit",
            "out_of_system",
            "postprocess",
            stream_class="postprocess",
        )
    )

    recycle_active = (
        bool(recycle_gas)
        or float(recirculation_frac) > 0.0
        or float(cyclone_capture_char_frac) > 0.0
        or float(cyclone_capture_ash_frac) > 0.0
    )
    if recycle_active:
        for return_id in return_ids:
            edges.append(
                ConnectivityEdge(
                    cyclone_ids[-1],
                    return_id,
                    "solid_separation",
                    "downward",
                    cyclone_coupling,
                    stream_class="rez",
                )
            )
        for return_id, bed_idx in zip(return_ids, bed_targets):
            bed_idx = int(np.clip(bed_idx, 0, max(int(n_bed_cells) - 1, 0)))
            edges.append(
                ConnectivityEdge(
                    return_id,
                    f"bed_cell_{bed_idx}",
                    "external_zirkulation",
                    "downward_to_bottom",
                    return_coupling,
                    stream_class="rez",
                )
            )

    graph = ConnectivityGraph(
        mode=str(preset),
        blocks=tuple(blocks),
        edges=tuple(edges),
        freeboard_active=bool(freeboard_active),
        recycle_active=bool(recycle_active),
        secondary_injection_active=bool(
            secondary_injection_applied and secondary_injection_segment is not None
        ),
        metadata={
            "preset": preset,
            "explicit_side_blocks": bool(explicit_side_blocks),
            "n_cyclone_stages": n_cyc,
            "n_return_legs": n_ret,
            "return_leg_bed_targets": bed_targets,
        },
    )
    graph.validate()
    if include_connection_pipe:
        graph = add_connection_pipe_passthrough(graph)
    return graph


def add_connection_pipe_passthrough(
    graph: ConnectivityGraph,
    *,
    pipe_block_id: str = "connection_pipe_0",
) -> ConnectivityGraph:
    """Insert a passthrough Verbindungsleitung block before the cyclone inlet edge."""
    cyclone_edges = [edge for edge in graph.edges if edge.kind == "cyclone_inlet"]
    if not cyclone_edges:
        return graph
    inlet_edge = cyclone_edges[0]
    blocks = list(graph.blocks)
    blocks.append(
        ConnectivityBlock(
            pipe_block_id,
            "connection_pipe",
            solver_coupling="passthrough",
        )
    )
    edges = [edge for edge in graph.edges if edge is not inlet_edge]
    edges.extend(
        [
            ConnectivityEdge(
                inlet_edge.from_block,
                pipe_block_id,
                "connection_duct",
                "upward",
                "passthrough",
                stream_class="auf",
            ),
            ConnectivityEdge(
                pipe_block_id,
                inlet_edge.to_block,
                "cyclone_inlet",
                "upward",
                inlet_edge.solver_coupling,
                stream_class="auf",
            ),
        ]
    )
    out = ConnectivityGraph(
        mode=graph.mode,
        blocks=tuple(blocks),
        edges=tuple(edges),
        freeboard_active=graph.freeboard_active,
        recycle_active=graph.recycle_active,
        secondary_injection_active=graph.secondary_injection_active,
        metadata={**dict(graph.metadata), "connection_pipe": pipe_block_id},
    )
    out.validate()
    return out


def build_connectivity_from_reactor_config(cfg: Any, *, reactor: Any | None = None) -> ConnectivityGraph:
    """Konfiguration-stage graph builder for ``ReactorConfig``."""
    path = getattr(cfg, "connectivity_graph_path", None)
    if path:
        return load_connectivity_graph_json(path)

    preset = str(getattr(cfg, "connectivity_preset", "from_config")).strip().lower()
    if preset in {"vtt_minimal_loop", "hamel_minimal", "from_config"}:
        n_fb = int(getattr(cfg, "n_freeboard_cells", 0))
        H_fb = float(getattr(cfg, "H_freeboard", 0.0))
        freeboard_active = bool(getattr(cfg, "thesis_mode", False)) and H_fb > 0.0 and n_fb > 0
        explicit_side = False
        if reactor is not None and hasattr(reactor, "_use_explicit_side_block_cells"):
            explicit_side = bool(reactor._use_explicit_side_block_cells())
        elif bool(getattr(cfg, "thesis_mode", False)):
            explicit_side = bool(getattr(cfg, "explicit_side_blocks_enabled", True))

        sec_mode = str(getattr(cfg, "freeboard_secondary_injection_mode", "lumped")).strip().lower()
        sec_applied = freeboard_active and (
            float(getattr(cfg, "freeboard_secondary_O2_mol_s", 0.0)) > 0.0
            or float(getattr(cfg, "freeboard_secondary_N2_mol_s", 0.0)) > 0.0
            or float(getattr(cfg, "freeboard_secondary_H2O_mol_s", 0.0)) > 0.0
        ) and sec_mode != "lumped"
        sec_segment = int(getattr(cfg, "freeboard_secondary_local_refine", 1)) - 1 if sec_applied else None

        mode_name = "hamel_minimal_connectivity" if preset == "from_config" else preset
        n_cyc = int(max(getattr(cfg, "n_cyclone_stages_thesis", 1), 1))
        n_ret = int(max(getattr(cfg, "n_return_legs_thesis", 1), 1))
        raw_targets = getattr(cfg, "return_leg_bed_targets_thesis", None)
        bed_targets = tuple(int(v) for v in raw_targets) if raw_targets else None
        return build_connectivity_graph_from_reactor_params(
            n_bed_cells=int(getattr(cfg, "n_cells", 0)),
            n_freeboard_cells=n_fb if freeboard_active else 0,
            H_bed=float(getattr(cfg, "H_bed", 0.0)),
            H_freeboard=H_fb if freeboard_active else 0.0,
            recycle_gas=bool(getattr(cfg, "recycle_gas", False)),
            recirculation_frac=float(getattr(cfg, "recirculation_frac", 0.0)),
            cyclone_capture_char_frac=float(getattr(cfg, "freeboard_cyclone_capture_char_frac", 0.0)),
            cyclone_capture_ash_frac=float(getattr(cfg, "freeboard_cyclone_capture_ash_frac", 0.0)),
            secondary_injection_applied=bool(sec_applied),
            secondary_injection_segment=sec_segment,
            explicit_side_blocks=explicit_side,
            preset=mode_name,
            n_cyclone_stages=n_cyc,
            n_return_legs=n_ret,
            return_leg_bed_targets=bed_targets,
        )

    raise ValueError(f"unsupported connectivity_preset: {preset!r}")


@dataclass(frozen=True)
class SolverTopologyIndices:
    """NR solver indices derived from a frozen ``ConnectivityGraph``."""

    n_bed: int
    n_freeboard: int
    n_solver: int
    cyclone_idx: int | None
    return_leg_idx: int | None
    cyclone_indices: tuple[int, ...]
    return_leg_indices: tuple[int, ...]
    top_bed_idx: int
    first_freeboard_idx: int | None
    last_freeboard_idx: int | None
    source_tail_idx: int
    explicit_freeboard_graph: bool
    has_side_blocks: bool


def solver_topology_from_graph(graph: ConnectivityGraph) -> SolverTopologyIndices:
    """Map connectivity blocks to NR solver-cell indices."""
    n_bed = sum(1 for block in graph.blocks if block.kind == "bed_cell")
    n_freeboard = sum(1 for block in graph.blocks if block.kind == "freeboard_cell")
    cyclone_indices = tuple(
        sorted(
            int(block.solver_cell_index)
            for block in graph.blocks
            if block.kind == "cyclone" and block.solver_cell_index is not None
        )
    )
    return_leg_indices = tuple(
        sorted(
            int(block.solver_cell_index)
            for block in graph.blocks
            if block.kind == "return_leg" and block.solver_cell_index is not None
        )
    )
    cyclone_idx = cyclone_indices[0] if cyclone_indices else None
    return_leg_idx = return_leg_indices[0] if return_leg_indices else None
    solver_indices = [
        int(block.solver_cell_index)
        for block in graph.blocks
        if block.solver_cell_index is not None
    ]
    n_solver = (max(solver_indices) + 1) if solver_indices else 0
    top_bed_idx = max(n_bed - 1, 0)
    first_freeboard_idx = n_bed if n_freeboard > 0 else None
    last_freeboard_idx = (n_bed + n_freeboard - 1) if n_freeboard > 0 else None
    source_tail_idx = last_freeboard_idx if last_freeboard_idx is not None else top_bed_idx
    has_side_blocks = bool(cyclone_indices) and bool(return_leg_indices)
    return SolverTopologyIndices(
        n_bed=n_bed,
        n_freeboard=n_freeboard,
        n_solver=n_solver,
        cyclone_idx=cyclone_idx,
        return_leg_idx=return_leg_idx,
        cyclone_indices=cyclone_indices,
        return_leg_indices=return_leg_indices,
        top_bed_idx=top_bed_idx,
        first_freeboard_idx=first_freeboard_idx,
        last_freeboard_idx=last_freeboard_idx,
        source_tail_idx=source_tail_idx,
        explicit_freeboard_graph=n_freeboard > 0,
        has_side_blocks=has_side_blocks,
    )


def build_solver_cell_registry(reactor: Any) -> dict[str, Any]:
    """Build block_id -> Cell registry for the active NR solver graph."""
    registry: dict[str, Any] = {}
    for i, cell in enumerate(getattr(reactor, "cells", ())):
        registry[f"bed_cell_{i}"] = cell
    for i, cell in enumerate(getattr(reactor, "freeboard_cells", ())):
        registry[f"freeboard_cell_{i}"] = cell
    by_id = getattr(reactor, "side_block_cells_by_id", None)
    if isinstance(by_id, dict):
        registry.update(by_id)
    else:
        cyclone = getattr(reactor, "cyclone_cell", None)
        if cyclone is not None:
            registry["cyclone_block"] = cyclone
        return_leg = getattr(reactor, "return_leg_cell", None)
        if return_leg is not None:
            registry["return_leg_block"] = return_leg
        for i, cell in enumerate(getattr(reactor, "cyclone_cells", ())):
            registry[side_block_id("cyclone", i, count=max(len(reactor.cyclone_cells), 1))] = cell
        for i, cell in enumerate(getattr(reactor, "return_leg_cells", ())):
            registry[side_block_id("return_leg", i, count=max(len(reactor.return_leg_cells), 1))] = cell
    return registry


def return_leg_upstream_block_id(graph: ConnectivityGraph, return_block_id: str) -> str | None:
    """Return the upstream block feeding one return-leg inlet edge."""
    for edge in graph.inbound_edges(return_block_id):
        stream_class = edge.stream_class or infer_stream_class(edge.kind, edge.direction)
        if stream_class == "rez" or edge.kind in {"solid_separation", "cyclone_stage"}:
            return edge.from_block
    return None


def cyclone_inlet_source_block_id(graph: ConnectivityGraph) -> str | None:
    """Return the block id feeding the cyclone inlet edge (skip passthrough segments)."""
    for edge in graph.edges:
        if edge.kind == "cyclone_inlet":
            return resolve_physical_block_id(graph, edge.from_block, direction="upstream")
    return None


_AUF_EDGE_KINDS = frozenset(
    {"axial_upflow", "freeboard_entry", "bubble_rupture_handoff", "cyclone_inlet"}
)
_AUF_SOLID_EDGE_KINDS = frozenset({"wake_solid_upflow", "axial_upflow"})
_AB_SOLID_EDGE_KINDS = frozenset({"internal_zirkulation"})
_PASSTHROUGH_BLOCK_KINDS = frozenset({"connection_pipe", "interface_block"})


def block_id_bed_index(block_id: str) -> int | None:
    prefix = "bed_cell_"
    if block_id.startswith(prefix):
        return int(block_id[len(prefix) :])
    return None


def block_id_freeboard_index(block_id: str) -> int | None:
    prefix = "freeboard_cell_"
    if block_id.startswith(prefix):
        return int(block_id[len(prefix) :])
    return None


def axial_upstream_block_id(graph: ConnectivityGraph, block_id: str) -> str | None:
    """Return the upstream ``auf`` source block for a bed/freeboard cell."""
    for edge in graph.inbound_edges(block_id):
        stream_class = edge.stream_class or infer_stream_class(edge.kind, edge.direction)
        if stream_class == "auf" or edge.kind in _AUF_EDGE_KINDS:
            return edge.from_block
    return None


def axial_downstream_bed_block_id(graph: ConnectivityGraph, block_id: str) -> str | None:
    """Return the downstream bed neighbor connected by ``axial_upflow``."""
    for edge in graph.outbound_edges(block_id):
        if edge.kind != "axial_upflow":
            continue
        if block_id_bed_index(edge.to_block) is not None:
            return edge.to_block
    return None


def resolve_bed_neighbor_block_ids(graph: ConnectivityGraph, bed_index: int) -> tuple[str | None, str | None]:
    """Return ``(upstream_block_id, downstream_bed_block_id)`` for one bed cell."""
    block_id = f"bed_cell_{int(bed_index)}"
    upstream_id = axial_upstream_block_id(graph, block_id)
    downstream_id = axial_downstream_bed_block_id(graph, block_id)
    return upstream_id, downstream_id


def ordered_freeboard_block_ids(graph: ConnectivityGraph) -> tuple[str, ...]:
    ids = [block.block_id for block in graph.blocks if block.kind == "freeboard_cell"]
    return tuple(sorted(ids, key=lambda bid: block_id_freeboard_index(bid) or 0))


def freeboard_chain_inlet_source_block_id(graph: ConnectivityGraph) -> str | None:
    """Resolve the physical source cell for ``freeboard_cell_0`` through interface hops."""
    upstream = axial_upstream_block_id(graph, "freeboard_cell_0")
    if upstream is None:
        return None
    if upstream == "bed_freeboard_interface":
        return axial_upstream_block_id(graph, upstream)
    return upstream


def solid_auf_source_block_id(graph: ConnectivityGraph, block_id: str) -> str | None:
    """Return the upstream solid ``K_auf m`` / upflow source block (Eq. 2.6 ``auf``)."""
    for edge in graph.inbound_edges(block_id):
        if edge.kind in _AUF_SOLID_EDGE_KINDS:
            return edge.from_block
    return None


def solid_ab_source_block_id(graph: ConnectivityGraph, block_id: str) -> str | None:
    """Return the upstream solid ``K_ab m`` source block (Eq. 2.6 ``ab`` backmix)."""
    for edge in graph.inbound_edges(block_id):
        if edge.kind in _AB_SOLID_EDGE_KINDS:
            return edge.from_block
        stream_class = edge.stream_class or infer_stream_class(edge.kind, edge.direction)
        if stream_class == "ab":
            return edge.from_block
    return None


def solid_ab_freeboard_downstream_block_id(graph: ConnectivityGraph, block_id: str) -> str | None:
    """Return the freeboard block above ``block_id`` for ``_solid_downflow_rates`` coupling."""
    for edge in graph.outbound_edges(block_id):
        if edge.kind != "axial_upflow":
            continue
        if block_id_freeboard_index(edge.to_block) is not None:
            return edge.to_block
    return None


def bed_top_freeboard_ab_source_block_id(graph: ConnectivityGraph) -> str | None:
    """Freeboard cell feeding the top bed ``m_solid_ab_in`` when no higher bed exists."""
    n_bed = sum(1 for block in graph.blocks if block.kind == "bed_cell")
    if n_bed <= 0:
        return None
    top_bed_id = f"bed_cell_{n_bed - 1}"
    if axial_downstream_bed_block_id(graph, top_bed_id) is not None:
        return None
    fb_ids = ordered_freeboard_block_ids(graph)
    return fb_ids[0] if fb_ids else None


def resolve_physical_block_id(
    graph: ConnectivityGraph,
    block_id: str,
    *,
    direction: str = "upstream",
) -> str | None:
    """Skip passthrough blocks (connection pipe / interface) to a physical cell id."""
    current = str(block_id)
    visited: set[str] = set()
    block_map = graph.block_by_id()
    while current in block_map and block_map[current].kind in _PASSTHROUGH_BLOCK_KINDS:
        if current in visited:
            return None
        visited.add(current)
        if direction == "upstream":
            nxt = axial_upstream_block_id(graph, current)
        else:
            nxt = None
            for edge in graph.outbound_edges(current):
                if edge.kind in _AUF_EDGE_KINDS or edge.kind == "axial_upflow":
                    nxt = edge.to_block
                    break
        if nxt is None:
            return None
        current = nxt
    return current


def resolve_registry_cell(
    graph: ConnectivityGraph,
    cell_registry: dict[str, Any],
    block_id: str | None,
    *,
    direction: str = "upstream",
) -> Any | None:
    """Map a graph block id to a ``Cell``, following passthrough segments."""
    if block_id is None:
        return None
    physical_id = resolve_physical_block_id(graph, block_id, direction=direction)
    if physical_id is None:
        return None
    return cell_registry.get(physical_id)


def jacobian_side_block_pairs_from_connectivity_graph(
    graph: ConnectivityGraph,
) -> frozenset[tuple[int, int]]:
    """Derive Nebenelemente block pairs from BC-propagation closure on the graph."""
    topo = solver_topology_from_graph(graph)
    if topo.n_solver <= 0:
        return frozenset()

    side_pairs: set[tuple[int, int]] = set()
    if topo.n_freeboard == 0 and topo.n_bed > 1 and not topo.has_side_blocks:
        if abs(topo.top_bed_idx - 0) > 1:
            side_pairs.add((0, topo.top_bed_idx))
            for row_idx in range(topo.n_bed):
                if abs(row_idx - 0) > 1:
                    side_pairs.add((row_idx, 0))
        return frozenset(side_pairs)

    if topo.has_side_blocks:
        source_tail_idx = int(topo.source_tail_idx)
        for cyclone_idx in topo.cyclone_indices:
            if abs(int(cyclone_idx) - source_tail_idx) > 1:
                side_pairs.add((int(cyclone_idx), source_tail_idx))
        for return_leg_idx in topo.return_leg_indices:
            if abs(int(return_leg_idx) - source_tail_idx) > 1:
                side_pairs.add((int(return_leg_idx), source_tail_idx))
        for cyclone_idx in topo.cyclone_indices:
            for return_leg_idx in topo.return_leg_indices:
                if abs(int(return_leg_idx) - int(cyclone_idx)) > 1:
                    side_pairs.add((int(return_leg_idx), int(cyclone_idx)))
        if topo.n_bed > 1:
            side_pairs.add((0, topo.top_bed_idx))
        if topo.last_freeboard_idx is not None and abs(topo.top_bed_idx - int(topo.last_freeboard_idx)) > 1:
            side_pairs.add((topo.top_bed_idx, int(topo.last_freeboard_idx)))

    for col_idx in range(topo.n_solver):
        affected = _affected_nr_with_topology(
            changed_cell_idx=col_idx,
            topo=topo,
        )
        for row_idx in affected:
            if abs(row_idx - col_idx) > 1:
                side_pairs.add((row_idx, col_idx))
    return frozenset(side_pairs)


def _affected_nr_with_topology(
    *,
    changed_cell_idx: int,
    topo: SolverTopologyIndices,
) -> tuple[int, ...]:
    """BC-propagation neighborhood using graph-derived solver topology."""
    idx = int(changed_cell_idx)
    n_bed = int(topo.n_bed)
    n_freeboard = int(topo.n_freeboard)
    cyclone_idx = topo.cyclone_idx
    return_leg_idx = topo.return_leg_idx
    cyclone_indices = topo.cyclone_indices
    return_leg_indices = topo.return_leg_indices
    top_bed_idx = int(topo.top_bed_idx)
    first_freeboard_idx = topo.first_freeboard_idx
    last_freeboard_idx = topo.last_freeboard_idx

    if topo.n_solver <= 0 or idx < 0 or idx >= topo.n_solver:
        return tuple()

    if topo.explicit_freeboard_graph:
        neighbors: set[int] = {idx}
        if idx < n_bed:
            if idx > 0:
                neighbors.add(idx - 1)
            if idx + 1 < n_bed:
                neighbors.add(idx + 1)
            if idx == 0:
                neighbors.update(range(n_bed))
            if idx == top_bed_idx:
                if first_freeboard_idx is not None:
                    neighbors.add(first_freeboard_idx)
                elif cyclone_idx is not None:
                    neighbors.add(cyclone_idx)
            return tuple(sorted(neighbors))

        if first_freeboard_idx is not None and idx < n_bed + n_freeboard:
            fb_idx = idx - n_bed
            if fb_idx > 0:
                neighbors.add(idx - 1)
            if fb_idx + 1 < n_freeboard:
                neighbors.add(idx + 1)
            if fb_idx == 0:
                neighbors.update(range(n_bed))
            if idx == last_freeboard_idx and cyclone_indices:
                neighbors.update(cyclone_indices)
                if n_bed > 0:
                    neighbors.add(0)
                if n_bed > 1:
                    neighbors.add(1)
            return tuple(sorted(neighbors))

        if idx in cyclone_indices:
            if last_freeboard_idx is not None:
                neighbors.add(last_freeboard_idx)
            elif n_bed > 0:
                neighbors.add(top_bed_idx)
            neighbors.update(return_leg_indices)
            if n_bed > 0:
                neighbors.add(0)
                neighbors.add(top_bed_idx)
            if n_bed > 1:
                neighbors.add(1)
            return tuple(sorted(neighbors))

        if idx in return_leg_indices:
            if n_bed > 0:
                neighbors.add(0)
                neighbors.add(top_bed_idx)
            if n_bed > 1:
                neighbors.add(1)
        return tuple(sorted(neighbors))

    if not topo.has_side_blocks:
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

    affected_side: set[int] = {idx}
    if idx < n_bed:
        if idx > 0:
            affected_side.add(idx - 1)
        if idx + 1 < n_bed:
            affected_side.add(idx + 1)
        if idx == n_bed - 1 and cyclone_indices:
            affected_side.update(cyclone_indices)
    elif cyclone_indices and idx in cyclone_indices:
        affected_side.update(return_leg_indices)
    elif return_leg_indices and idx in return_leg_indices:
        affected_side.add(0)
    return tuple(sorted(affected_side))


def affected_nr_residual_cells(
    *,
    changed_cell_idx: int,
    n_bed: int,
    n_freeboard: int,
    explicit_freeboard_graph: bool,
    side_blocks_in_boundary_path: bool,
    has_side_blocks: bool,
    connectivity_graph: ConnectivityGraph | None = None,
) -> tuple[int, ...]:
    """Return affected residual cell indices for one changed solver cell."""
    if connectivity_graph is not None:
        topo = solver_topology_from_graph(connectivity_graph)
        return _affected_nr_with_topology(changed_cell_idx=int(changed_cell_idx), topo=topo)

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

        neighbors_legacy: set[int] = {idx}
        if idx < int(n_bed):
            if idx > 0:
                neighbors_legacy.add(idx - 1)
            if idx + 1 < int(n_bed):
                neighbors_legacy.add(idx + 1)
            if idx == 0:
                neighbors_legacy.update(range(int(n_bed)))
            if idx == top_bed_idx:
                if first_freeboard_idx is not None:
                    neighbors_legacy.add(first_freeboard_idx)
                elif cyclone_idx is not None:
                    neighbors_legacy.add(cyclone_idx)
            return tuple(sorted(neighbors_legacy))

        if first_freeboard_idx is not None and idx < int(n_bed) + int(n_freeboard):
            fb_idx = idx - int(n_bed)
            if fb_idx > 0:
                neighbors_legacy.add(idx - 1)
            if fb_idx + 1 < int(n_freeboard):
                neighbors_legacy.add(idx + 1)
            if fb_idx == 0:
                neighbors_legacy.update(range(int(n_bed)))
            if idx == last_freeboard_idx and cyclone_idx is not None:
                neighbors_legacy.add(cyclone_idx)
                if int(n_bed) > 0:
                    neighbors_legacy.add(0)
                if int(n_bed) > 1:
                    neighbors_legacy.add(1)
            return tuple(sorted(neighbors_legacy))

        if cyclone_idx is not None and idx == cyclone_idx:
            if last_freeboard_idx is not None:
                neighbors_legacy.add(last_freeboard_idx)
            elif int(n_bed) > 0:
                neighbors_legacy.add(top_bed_idx)
            if return_leg_idx is not None:
                neighbors_legacy.add(return_leg_idx)
            if int(n_bed) > 0:
                neighbors_legacy.add(0)
                neighbors_legacy.add(top_bed_idx)
            if int(n_bed) > 1:
                neighbors_legacy.add(1)
            return tuple(sorted(neighbors_legacy))

        if return_leg_idx is not None and idx == return_leg_idx:
            if int(n_bed) > 0:
                neighbors_legacy.add(0)
                neighbors_legacy.add(top_bed_idx)
            if int(n_bed) > 1:
                neighbors_legacy.add(1)

        return tuple(sorted(neighbors_legacy))

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
    affected_side: set[int] = {idx}
    if idx < n_bed:
        if idx > 0:
            affected_side.add(idx - 1)
        if idx + 1 < n_bed:
            affected_side.add(idx + 1)
        if idx == n_bed - 1:
            affected_side.add(cyclone_idx)
    elif idx == cyclone_idx:
        affected_side.add(return_leg_idx)
    elif idx == return_leg_idx:
        affected_side.add(0)
    return tuple(sorted(affected_side))
