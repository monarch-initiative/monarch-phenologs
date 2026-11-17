"""Render every registered figure to out/figures/. `python -m phenologs_talk.render_figures [name ...]`"""
import sys
from pathlib import Path

from . import figures as F
from . import figures_results  # noqa: F401  (registers its figures)

OUT = Path(__file__).resolve().parents[1] / "out" / "figures"
BUILDS = {"phenolog_cartoon": [1, 2, 3], "phenolog_idea": [1, 2], "kg_components": [1, 2, 3], "phenolog_concept": [1, 2],
          "marbles": [1, 2, 3], "eval_results": [1, 2]}     # figures drawn in steps: one image per step


def render(names=None, data=None):
    D = data or F.load_data()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, (fn, _alt) in F.REGISTRY.items():
        if names and name not in names:
            continue
        for step in BUILDS.get(name, [None]):
            key = f"{name}-{step}" if step else name
            paths[key] = fn(D, step=step).save(OUT / f"{key}.png")
    return paths


if __name__ == "__main__":
    for key, path in render(sys.argv[1:] or None).items():
        print(f"{key:24s} {path}")
