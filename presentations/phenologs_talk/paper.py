"""Panels cut straight from the preprint PDF, so result plots on slides are the published ones.

Each crop is a page number and a box in fractions of the page (left, top, right, bottom). If the
preprint is revised, drop the new PDF into assets/ and adjust the boxes here; `python -m
phenologs_talk.build crops` writes them to out/crops/ for a quick look.

Needs `pdftoppm` (poppler: `brew install poppler`).
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "assets" / "phenologs_preprint_biorxiv.pdf"
PAGE_W_IN, PAGE_H_IN = 8.5, 11.0

# name: (pdf page, (left, top, right, bottom)), what it is
CROPS = {
    "fig2a": (7, (0.096, 0.062, 0.500, 0.316)),     # disease network, real vs randomized, human-mouse
    "fig2b": (7, (0.508, 0.062, 0.912, 0.316)),     # phenotype network, same
    "fig3": (10, (0.086, 0.078, 0.920, 0.352)),     # temporal validation, panels A and B
    "fig3a": (10, (0.086, 0.078, 0.502, 0.352)),
    "fig3b": (10, (0.505, 0.078, 0.920, 0.352)),
    "fig4abc": (13, (0.118, 0.190, 0.880, 0.368)),  # species contributions
    "fig4d": (13, (0.118, 0.388, 0.845, 0.818)),    # Carney complex additive decomposition
    "fig5": (15, (0.118, 0.262, 0.880, 0.860)),     # Bardet-Biedl additive decomposition
    "figS3": (27, (0.090, 0.078, 0.905, 0.482)),    # FDR cutoff histogram
    "figS7ab": (35, (0.086, 0.066, 0.935, 0.468)),  # leave-one-out heatmaps
    "figS7cd": (35, (0.086, 0.482, 0.935, 0.800)),  # recall and AUPRC, hg vs nb (complete model, no genes withheld)
    "figS7c": (35, (0.086, 0.482, 0.492, 0.800)),   # recall of known genes within the top N
    "figS7d": (35, (0.500, 0.482, 0.935, 0.800)),   # per-disease AUPRC scatter, nb against hg
}


def crop(name: str, out_dir: Path, dpi: int = 300) -> Path:
    """Render one crop to out_dir/<name>.png (cached on the PDF's mtime) and return its path."""
    page, (l, t, r, b) = CROPS[name]
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{name}.png"
    if out.exists() and out.stat().st_mtime >= max(PDF.stat().st_mtime, Path(__file__).stat().st_mtime):
        return out
    if not shutil.which("pdftoppm"):
        raise RuntimeError("pdftoppm not found. Install poppler (brew install poppler) to cut panels from the preprint.")
    x, y = round(l * PAGE_W_IN * dpi), round(t * PAGE_H_IN * dpi)
    w, h = round((r - l) * PAGE_W_IN * dpi), round((b - t) * PAGE_H_IN * dpi)
    subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-f", str(page), "-l", str(page), "-x", str(x), "-y", str(y),
                    "-W", str(w), "-H", str(h), "-singlefile", str(PDF), str(out.with_suffix(""))], check=True)
    return out


def crop_all(out_dir: Path, dpi: int = 300) -> list[Path]:
    return [crop(name, out_dir, dpi) for name in CROPS]
