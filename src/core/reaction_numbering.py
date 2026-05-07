"""Reaction numbering mapping utilities.

Authority:
- Thesis numbering follows Hamel (1999).
- Implementation numbering follows current Python legacy labels.
"""

from __future__ import annotations

from typing import Dict

# impl id -> thesis id (None means not in thesis R1-R11 backbone)
IMPL_TO_THESIS: dict[str, str | None] = {
    "R1": "R1",
    "R2": "R2",
    "R3": "R4",   # impl R3 = hydrogenating gasification
    "R4": "R3",   # impl R4 = Boudouard
    "R5": "R5",
    "R6": "R7",   # impl R6 = CH4 oxidation
    "R7": "R9",   # impl R7 = methane steam reforming
    "R8": "R8",
    "R9": None,   # impl sulfur placeholder path
    "R10": "R10",
    "R11": "R11",
    "R11b": "R11",
    "R11d": "R11",
    "R12": "R6",  # impl R12 = H2 oxidation
}

THESIS_TO_IMPL_PRIMARY: dict[str, str] = {
    "R1": "R1",
    "R2": "R2",
    "R3": "R4",
    "R4": "R3",
    "R5": "R5",
    "R6": "R12",
    "R7": "R6",
    "R8": "R8",
    "R9": "R7",
    "R10": "R10",
    "R11": "R11",
}


def reaction_numbering_metadata() -> dict[str, object]:
    """Return a stable metadata dict for reports/diagnostics."""
    return {
        "authority": "hamel_1999_thesis",
        "impl_to_thesis": dict(IMPL_TO_THESIS),
        "thesis_to_impl_primary": dict(THESIS_TO_IMPL_PRIMARY),
    }


def map_impl_diag_to_thesis(diag_impl: Dict[str, float]) -> dict[str, float]:
    """Aggregate implementation-diag entries into thesis labels.

    For aliases like impl R11b/R11d, values are accumulated into thesis R11.
    Entries without thesis mapping (impl R9) are kept under `_impl_only:R9`.
    """
    out: dict[str, float] = {}
    for impl_id, val in diag_impl.items():
        thesis_id = IMPL_TO_THESIS.get(str(impl_id))
        if thesis_id is None:
            out[f"_impl_only:{impl_id}"] = float(out.get(f"_impl_only:{impl_id}", 0.0) + float(val))
            continue
        out[thesis_id] = float(out.get(thesis_id, 0.0) + float(val))
    return out


def map_impl_profile_diag_to_thesis(profile_diag_impl: dict[str, list[float]]) -> dict[str, list[float]]:
    """Map per-axial reaction profile diag from impl ids to thesis ids."""
    out: dict[str, list[float]] = {}
    for impl_id, seq in profile_diag_impl.items():
        thesis_id = IMPL_TO_THESIS.get(str(impl_id))
        key = thesis_id if thesis_id is not None else f"_impl_only:{impl_id}"
        if key not in out:
            out[key] = [0.0 for _ in seq]
        out[key] = [float(a) + float(b) for a, b in zip(out[key], seq)]
    return out
