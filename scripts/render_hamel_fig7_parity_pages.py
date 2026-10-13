#!/usr/bin/env python3
"""Render Hamel PDF parity figure pages for digitization (Fig 7.10, 7.17)."""

from __future__ import annotations

import argparse
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PDF = ROOT / (
    "docs/_2001_VDI-Dis._Mathematische Modellierung und experimentelle "
    "Untersuchung der Vergasung verschiedener fester Brennstoffe.pdf"
)
OUT_DIR = ROOT / "data" / "hamel_fig7_parity_crops"

# 1-based PDF page numbers (from fitz scan)
FIG_PAGES = {
    "7.10": 78,
    "7.13": 80,
    "7.17": 83,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--scale", type=float, default=3.0)
    parser.add_argument("--figures", nargs="*", default=["7.10", "7.17"])
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(args.pdf)
    mat = fitz.Matrix(args.scale, args.scale)

    for fig in args.figures:
        pno = FIG_PAGES[fig]
        page = doc[pno - 1]
        pix = page.get_pixmap(matrix=mat, alpha=False)
        out = OUT_DIR / f"page_{pno}_fig{fig.replace('.', '')}_parity.png"
        pix.save(out)
        print(f"Fig {fig} → {out} ({pix.width}×{pix.height})")


if __name__ == "__main__":
    main()
