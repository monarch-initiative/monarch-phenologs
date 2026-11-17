"""The talk's diagrams, one function per figure.

Each function takes the data (data.yaml as a Dot) and returns a Canvas. Nothing here holds a
number: change data.yaml and rebuild. Register a new figure with @figure("name", "alt text")
and refer to it from deck.yaml as `figure: name`.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import yaml
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path as MplPath

from . import theme as T
from .draw import Canvas, sci

HERE = Path(__file__).resolve().parent
REGISTRY: dict[str, tuple] = {}


class Dot(dict):
    """dict with attribute access, so deck text can say {totals.disease_phenologs:,}."""
    __getattr__ = dict.__getitem__


def dot(o):
    if isinstance(o, dict):
        return Dot({k: dot(v) for k, v in o.items()})
    if isinstance(o, list):
        return [dot(v) for v in o]
    return o


def load_data(path: Path | None = None) -> Dot:
    D = dot(yaml.safe_load((path or HERE / "data.yaml").read_text()))
    assert [s.key for s in D.species] == T.SPECIES_ORDER, "species in data.yaml must follow theme.SPECIES_ORDER"
    D["sp"] = Dot({s.key: s for s in D.species})
    return D


def figure(name: str, alt: str):
    def deco(fn):
        REGISTRY[name] = (fn, alt)
        return fn
    return deco


def n(x) -> str:
    return f"{x:,}"


def short_binomial(b: str) -> str:
    genus, rest = b.split(" ", 1)
    return f"{genus[0]}. {rest}"


# =============================================================================================
# 1. The problem
# =============================================================================================
@figure("diagnostic_gap", "Ten patient figures, five to seven of them shaded as still undiagnosed after exome or "
                          "genome sequencing, beside four reasons a patient stays undiagnosed.")
def fig_diagnostic_gap(D, step=None):
    c = Canvas(h=T.FIG_H_TALL)                       # this slide has no takeaway line, so the figure is taller
    lo, hi = 100 - D.problem.yield_high, 100 - D.problem.yield_low        # undiagnosed = 100 - diagnostic yield
    big = f"{lo}–{hi}%"
    c.text(0.1, 1.05, big, size=60, weight="bold")
    bw = c.measure(big, 60, "bold")[0]
    cap = c.wrap("of rare-disease patients with sequencing data (exome or genome) still go undiagnosed",
                 6.7 - (0.1 + bw + 0.28), 14.5)
    c.text(0.1 + bw + 0.28, 1.05, cap, size=14.5, color=T.MUTED)

    count, x0, dx, top, h = 10, 0.5, 0.67, 2.4, 1.75
    n_lo, n_hi = round(lo / 10), round(hi / 10)
    for i in range(count):
        col = T.NAVY if i < n_lo else (T.CAROLINA if i < n_hi else "#CDD5DD")
        c.person(x0 + i * dx, top, h, color=col)

    def bracket(i0, i1, label, color):
        xa, xb, y = x0 + i0 * dx - 0.27, x0 + i1 * dx + 0.27, top + h + 0.22
        c.line(xa, y, xb, y, color=color, lw=2.2)
        c.text((xa + xb) / 2, y + 0.28, label, size=12.5, color=color, ha="center", weight="bold")

    bracket(0, n_lo - 1, "undiagnosed", T.NAVY)
    bracket(n_lo, n_hi - 1, "varies by cohort", T.CAROLINA_DARK)
    bracket(n_hi, count - 1, "diagnosed", T.MUTED)
    c.text(0.1, c.h - 0.22, D.problem.yield_source, size=10, color=T.MUTED)

    px, ph = 7.0, c.h - 0.2
    c.panel(px, 0.1, 5.0, ph)
    c.text(px + 0.3, 0.5, "Why no diagnosis?", size=16, weight="bold")
    cards = D.problem.reasons
    y0 = 0.9
    pitch = (0.1 + ph - 0.26 - y0 + 0.14) / len(cards)
    assert cards[-1].get("focus") and sum(bool(r.get("focus")) for r in cards) == 1, "the talk's reason goes last"
    for i, r in enumerate(cards):
        y = y0 + i * pitch
        hot = bool(r.get("focus"))
        c.box(px + 0.3, y, 4.4, pitch - 0.14, fill=T.ACCENT_TINT if hot else T.WHITE,
              edge=T.ACCENT if hot else T.RULE, lw=1.8 if hot else 1.0, r=0.12)
        mid = y + (pitch - 0.14) / 2
        subs = [r.sub] if isinstance(r.sub, str) else list(r.sub)
        size = min(c.fit(s, 3.9, 12) for s in subs)
        top = mid - 0.18 - 0.12 * (len(subs) - 1)
        c.text(px + 0.55, top, r.head, size=c.fit(r.head, 3.9, 15, "bold"), weight="bold",
               color=T.ACCENT_DARK if hot else T.INK)
        for k, s in enumerate(subs):
            c.text(px + 0.55, top + 0.37 + 0.24 * k, s, size=size, color=T.INK if hot else T.MUTED)
    return c


@figure("bbs_overview", "A patient figure with six features of Bardet-Biedl syndrome, beside three boxes: the number "
                        "of known genes, the number of known phenotypes, and that about 1 in 5 clinically diagnosed "
                        "patients has no identifiable pathogenic variant.")
def fig_bbs_overview(D, step=None):
    c = Canvas()
    b = D.bbs
    cx, top, h = 3.6, 0.5, 3.75
    c.person(cx, top, h, color=T.NAVY)
    head_y = top + h * 0.15
    body_l, body_r = cx - h * 0.18, cx + h * 0.18
    f = b.features
    left = [(f[0], 0.95, (cx - h * 0.15, head_y)), (f[1], 2.6, (body_l, 2.7)), (f[2], 3.8, (body_l, 3.6))]
    right = [(f[5], 0.95, (cx + h * 0.15, head_y)), (f[3], 2.6, (body_r, 2.7)), (f[4], 3.8, (body_r, 3.7))]
    lx, rx = 2.5, 4.7
    for label, y, (tx, ty) in left:
        w = c.measure(label, 13.5)[0] + 0.5
        c.pill(lx - w / 2, y, label, w=w, size=13.5)
        c.scurve(lx, y, tx, ty, color=T.CAROLINA, lw=1.5)
        c.circle(tx, ty, 0.055, fill=T.WHITE, edge=T.CAROLINA, lw=1.5, zorder=6)
    for label, y, (tx, ty) in right:
        w = c.measure(label, 13.5)[0] + 0.5
        c.pill(rx + w / 2, y, label, w=w, size=13.5)
        c.scurve(rx, y, tx, ty, color=T.CAROLINA, lw=1.5)
        c.circle(tx, ty, 0.055, fill=T.WHITE, edge=T.CAROLINA, lw=1.5, zorder=6)

    px, pw, bh, gap = 7.6, 4.4, 1.25, 0.17

    def stat(y, big, label, sub):
        c.panel(px, y, pw, bh)
        c.text(px + 0.3, y + bh / 2 + 0.02, big, size=54, weight="bold")
        tx = px + 0.3 + c.measure(big, 54, "bold")[0] + 0.25
        c.text(tx, y + 0.37, label, size=17, weight="bold")
        c.text(tx, y + 0.83, sub, size=11.5, color=T.MUTED)

    stat(0.1, n(b.known_genes), "genes known", "BBSome · chaperonin complex ·\nintraflagellar transport")
    ex = [e[0].lower() + e[1:] for e in b.phenotype_examples] + ["…"]
    tw = pw - 0.6 - c.measure(n(b.phenotypes_known), 54, "bold")[0] - 0.25
    lines = [ex[0]]
    for e in ex[1:]:                                  # break between terms, never inside one
        if c.measure(f"{lines[-1]} · {e}", 11.5)[0] <= tw:
            lines[-1] += f" · {e}"
        else:
            lines[-1] += " ·"
            lines.append(e)
    stat(0.1 + bh + gap, n(b.phenotypes_known), "phenotypes known", "\n".join(lines))

    y3 = 0.1 + 2 * (bh + gap)
    c.panel(px, y3, pw, bh)
    hot = {"color": T.ACCENT_DARK, "weight": "bold"}
    c.rich(px + 0.3, y3 + bh / 2, [("About", None), ("1 in 5 clinically diagnosed patients", hot),
                                   ("has no identifiable pathogenic variant.", None)], 3.5, size=15)

    c.text(px + 0.05, 4.52, "Full list of causal and modifier genes is still incomplete", size=12.5, style="italic",
           color=T.MUTED)
    src = "Forsythe & Beales 2013 · Niederlova 2019 · Perea-Romero 2022 · Monarch KG"
    c.text(px + 0.05, 4.85, src, size=c.fit(src, pw - 0.1, 9.5, min_size=7.5), color=T.FAINT)
    return c


# =============================================================================================
# 2. The dataset
# =============================================================================================
def _cluster(c, cx, cy, R, n_genes, n_phen, color, rng, node_r=0.055, term_fill=T.WHITE):
    """Cartoon gene-phenotype network: phenotype squares inside, gene circles around.
    term_fill: the squares' fill; pass the network colour to draw them solid (diseases)."""
    P = []
    for i in range(n_phen):
        a = 2 * math.pi * i / n_phen + rng.uniform(-0.35, 0.35)
        rr = R * 0.40 * rng.uniform(0.55, 1.1) if n_phen > 1 else 0.0
        P.append((cx + rr * math.cos(a), cy + rr * math.sin(a)))
    G = []
    for i in range(n_genes):
        a = 2 * math.pi * i / n_genes + rng.uniform(-0.10, 0.10)
        rr = R * rng.uniform(0.80, 1.0)
        G.append((cx + rr * math.cos(a), cy + rr * math.sin(a)))
    for gx, gy in G:
        order = sorted(range(n_phen), key=lambda j: (P[j][0] - gx) ** 2 + (P[j][1] - gy) ** 2)
        links = order[:1] + (order[1:2] if n_phen > 1 and rng.random() < 0.45 else [])
        for j in links:
            c.line(gx, gy, P[j][0], P[j][1], color=color, lw=0.8, alpha=0.45, zorder=1)
    for gx, gy in G:
        c.circle(gx, gy, node_r, fill=color, edge=T.WHITE, lw=0.6, zorder=3)
    for px, py in P:
        c.square(px, py, node_r * 2.7, fill=term_fill, edge=color, lw=1.5, zorder=4, r=0.025)
    return G, P


def _kg_networks(D, orthology):
    """Human and the seven model species as separate cartoon networks, in order of divergence.

    orthology=True adds the dashed arcs that join every model network to the human one (and the
    room they need above the clusters); without them the clusters are a little larger.
    """
    c = Canvas()
    rng = np.random.default_rng(11)
    H = D.human
    k = 1.0 if orthology else 1.1                     # cluster scale
    cy, hx, hR = (2.1 if orthology else 1.66), 1.0, 0.72 * k
    xs = [2.80, 4.25, 5.58, 6.91, 8.39, 9.88, 11.24]  # spaced so neighbouring clusters and captions both clear

    # legend
    c.circle(0.2, 0.22, 0.06, fill=T.MUTED, edge=None)
    c.text(0.34, 0.22, "gene", size=11, color=T.MUTED)
    c.square(1.0, 0.22, 0.15, fill=T.WHITE, edge=T.MUTED, lw=1.5)
    c.text(1.16, 0.22, "phenotype or disease", size=11, color=T.MUTED)
    c.line(2.95, 0.22, 3.4, 0.22, color=T.MUTED, lw=1.0)
    c.text(3.5, 0.22, "curated edge", size=11, color=T.MUTED)
    if orthology:
        c.line(4.75, 0.22, 5.2, 0.22, color=T.MUTED, lw=1.2, ls=(0, (3, 3)))
        c.text(5.3, 0.22, "orthology (PANTHER family)", size=11, color=T.MUTED)

    _cluster(c, hx, cy, hR, 26, 6, T.NAVY, rng)
    lo, hi = math.log10(600), math.log10(260000)
    for i, (sp, x) in enumerate(zip(D.species, xs)):
        frac = (math.log10(sp.annotations) - lo) / (hi - lo)
        R = k * (0.27 + 0.43 * frac)
        n_genes = int(round(5 + 23 * frac))
        _cluster(c, x, cy, R, n_genes, max(2, n_genes // 4), sp.color, rng)
        if orthology:
            c.arc(hx, cy - hR - 0.08, x, cy - R - 0.08, lift=0.2 + 0.085 * i, color=T.FAINT, lw=1.1, ls=(0, (3, 3)))

    by = cy + hR + 0.32                               # the row of species badges
    cap = "gene–phenotype edges"
    cap_size = c.fit(cap, min(b - a for a, b in zip(xs, xs[1:])) - 0.14, 9.5, min_size=8.0)

    def label_col(x, code, color, name, binom, lines):
        c.badge(x, by, code, color, r=0.17)
        c.text(x, by + 0.35, name, size=12.5, weight="bold", ha="center")
        c.text(x, by + 0.58, short_binomial(binom), size=9.5, style="italic", color=T.MUTED, ha="center")
        y = by + 0.87
        for num, what in lines:
            c.text(x, y, num, size=13.5, weight="bold", color=color, ha="center")
            c.text(x, y + 0.2, what, size=cap_size, color=T.MUTED, ha="center")
            y += 0.45

    label_col(hx, H.code, T.NAVY, H.common, H.binomial,
              [(n(H.phenotype_annotations), cap), (n(H.gene_disease_pairs), "gene–disease edges")])
    for sp, x in zip(D.species, xs):
        label_col(x, sp.code, sp.color, sp.common, sp.binomial, [(n(sp.annotations), cap)])

    # divergence axis
    ay = by + 1.52
    c.arrow(xs[0] - 0.55, ay, xs[-1] + 0.62, ay, color=T.FAINT, lw=1.2, head=0.1)
    for sp, x in zip(D.species, xs):
        c.line(x, ay - 0.05, x, ay + 0.05, color=T.FAINT, lw=1.2)
        c.text(x, ay + 0.2, n(sp.mya), size=10, color=T.MUTED, ha="center")
    c.text((xs[0] + xs[-1]) / 2, ay + 0.43, "million years since the last common ancestor with human (TimeTree, approximate)",
           size=9.5, color=T.FAINT, ha="center")
    return c


@figure("kg_networks", "Eight separate cartoon gene-phenotype networks in a row, human then the seven model species "
                       "ordered by evolutionary distance from mouse to slime mold, each coloured by species, sized by "
                       "its number of gene-phenotype edges and labelled with that number.")
def fig_kg_networks(D, step=None):
    return _kg_networks(D, orthology=False)


@figure("kg_networks_orthology", "The eight gene-phenotype networks, human then the seven model species, with every "
                                 "model network joined to the human one by a dashed orthology arc.")
def fig_kg_networks_orthology(D, step=None):
    return _kg_networks(D, orthology=True)


@figure("fill_gaps", "The gene-phenotype networks of seven model organisms on the left, with a dashed link through a "
                     "question mark to human diseases and human phenotypes on the right. Below, a bar chart: 463,460 "
                     "non-human gene-phenotype edges against 131,495 human gene-phenotype edges and 12,426 human "
                     "gene-disease edges, with the question: how can we leverage this information?")
def fig_fill_gaps(D, step=None):
    c = Canvas()
    rng = np.random.default_rng(5)
    H = D.human

    # --- left: one cartoon network per model species (circles are genes, squares are phenotypes)
    px, pw, py, ph = 0.1, 5.3, 0.05, 3.42
    c.panel(px, py, pw, ph, title="Model organism gene–phenotype networks", title_color=T.INK)
    step_x = pw / 4
    rows = [(1.22, D.species[:4]), (2.52, D.species[4:])]
    lx = px + 3.5 * step_x - 0.4                      # legend, in the free cell at the end of the second row
    c.circle(lx, rows[1][0] - 0.13, 0.05, fill=T.MUTED, edge=None)
    c.text(lx + 0.14, rows[1][0] - 0.13, "gene", size=10, color=T.MUTED)
    c.square(lx, rows[1][0] + 0.13, 0.12, fill=T.WHITE, edge=T.MUTED, lw=1.3)
    c.text(lx + 0.14, rows[1][0] + 0.13, "phenotype", size=10, color=T.MUTED)
    lo, hi = math.log10(600), math.log10(260000)
    for cy, group in rows:
        for j, sp in enumerate(group):
            frac = (math.log10(sp.annotations) - lo) / (hi - lo)
            R = 0.21 + 0.26 * frac
            n_genes = int(round(5 + 15 * frac))
            x = px + step_x / 2 + j * step_x
            _cluster(c, x, cy, R, n_genes, max(2, n_genes // 4), sp.color, rng, node_r=0.044)
            c.text(x, cy + 0.65, sp.common, size=10.5, weight="bold", ha="center")

    # --- right: the human side, rolled up to its two kinds of term; the patient stays in the background
    qx, qy = 6.25, (rows[0][0] + rows[1][0]) / 2
    tl, tw, th = 7.05, 3.0, 0.6
    y_d, y_p = qy - 0.52, qy + 0.52
    c.pill(tl + tw / 2, y_d, "Human diseases", w=tw, h=th, size=16, weight="bold", fill=T.NAVY, edge=None, color=T.WHITE)
    c.pill(tl + tw / 2, y_p, "Human phenotypes", w=tw, h=th, size=16, weight="bold", lw=1.8)
    c.person(10.95, 0.47, 2.75, color="#DCE3EA")

    # --- the open question: can the networks on the left say something about the terms on the right?
    dash = (0, (3, 2.5))
    c.line(px + pw + 0.06, qy, qx, qy, color=T.CAROLINA, lw=1.6, ls=dash)
    for y in (y_d, y_p):
        c.scurve(qx, qy, tl - 0.08, y, color=T.CAROLINA, lw=1.6, ls=dash, head=0.1)
    c.circle(qx, qy, 0.27, fill=T.ACCENT, edge=T.WHITE, lw=2.0, zorder=7)
    c.text(qx, qy + 0.01, "?", size=19, weight="bold", color=T.WHITE, ha="center", zorder=8)

    # --- how many edges of each kind the graph holds (one scale, labelled directly)
    labels = ["Non-human gene–phenotype edges", "Human gene–phenotype edges", "Human gene–disease edges"]
    x0 = 0.1 + max(c.measure(s, 12.5, "bold")[0] for s in labels) + 0.15
    wmax, bh, by = 6.5, 0.28, [4.03, 4.45, 4.87]
    total = D.totals.model_annotations
    for lab, y in zip(labels, by):
        c.text(x0 - 0.15, y, lab, size=12.5, weight="bold", ha="right")
    x = x0
    for sp in D.species:                              # the non-human bar is split by species, same colours as above
        w = wmax * sp.annotations / total
        c.rect(x, by[0] - bh / 2, max(w - 0.02, 0.012), bh, fill=sp.color)
        if w > 0.2:
            c.text(x + w / 2, by[0] - bh / 2 - 0.14, sp.common.lower(), size=9.5, color=T.MUTED, ha="center")
        x += w
    nx = x0 + wmax + 0.12
    c.text(nx, by[0], n(total), size=13.5, weight="bold")
    for val, y in zip([H.phenotype_annotations, H.gene_disease_pairs], by[1:]):
        w = wmax * val / total
        c.box(x0, y - bh / 2, w, bh, fill=T.NAVY, edge=None, r=0.04)
        c.text(x0 + w + 0.12, y, n(val), size=13.5, weight="bold")

    # --- the question the talk answers, pointed at the total (all species, not one segment of the bar)
    q = "How can we leverage this information?"
    ax_ = nx + c.measure(n(total), 13.5, "bold")[0] / 2
    ty = by[2]
    c.text(ax_ - 0.42, ty, q, size=15, weight="bold", color=T.ACCENT_DARK, ha="right")
    c.line(ax_ - 0.32, ty, ax_, ty, color=T.ACCENT_DARK, lw=1.8)
    c.arrow(ax_, ty + 0.01, ax_, by[0] + 0.17, color=T.ACCENT_DARK, lw=1.8, head=0.12)
    return c


@figure("bbs_fill_gaps", "Two halves. Left: Bardet-Biedl syndrome as an example complex disease, a small patient figure "
                         "with six clinical features, 26 known genes and 95 known phenotypes, and the gap: about 1 in "
                         "5 clinically diagnosed patients has no identifiable pathogenic variant. Right: the "
                         "gene-phenotype networks of seven model organisms, linked back to the disease through a "
                         "question mark, above a bar chart: 463,460 non-human gene-phenotype edges against 131,495 "
                         "human gene-phenotype edges and 12,426 human gene-disease edges, with the question: how can "
                         "we leverage this information?")
def fig_bbs_fill_gaps(D, step=None):
    """The disease and its gaps (bbs_overview) and the model organism data (fill_gaps) on one slide."""
    c = Canvas()
    rng = np.random.default_rng(5)
    H, b = D.human, D.bbs

    # ---------------------------------------------------------------- left: the disease and what is missing
    lw_ = 5.55
    c.text(0.05, 0.27, f"Example complex disease: {b.label}", size=14, weight="bold")
    c.text(0.05, 0.56, f"autosomal recessive ciliopathy, about {b.prevalence} people", size=11, color=T.MUTED)
    cx, top, h = 2.62, 0.84, 1.9
    c.person(cx, top, h, color=T.NAVY)
    head_y, body = top + h * 0.15, h * 0.18
    f = b.features
    spots = [(f[0], -1, top + 0.28, head_y), (f[1], -1, top + 1.05, top + 1.05), (f[2], -1, top + 1.7, top + 1.55),
             (f[5], 1, top + 0.28, head_y), (f[3], 1, top + 1.05, top + 1.05), (f[4], 1, top + 1.7, top + 1.62)]
    for label, side, y, ty in spots:
        w = c.measure(label, 10.5)[0] + 0.36
        px_ = cx + side * 0.78
        tx = cx + side * (h * 0.15 if ty == head_y else body)
        c.pill(px_ + side * w / 2, y, label, w=w, h=0.34, size=10.5, lw=1.3)
        c.scurve(px_, y, tx, ty, color=T.CAROLINA, lw=1.3)
        c.circle(tx, ty, 0.045, fill=T.WHITE, edge=T.CAROLINA, lw=1.3, zorder=6)

    ty0, th, tg = 2.98, 0.68, 0.12
    tw_ = (lw_ - tg) / 2
    for i, (big, label) in enumerate([(n(b.known_genes), "genes known"), (n(b.phenotypes_known), "phenotypes known")]):
        x = i * (tw_ + tg)
        c.panel(x, ty0, tw_, th)
        c.text(x + 0.22, ty0 + th / 2 + 0.01, big, size=28, weight="bold")
        c.text(x + 0.22 + c.measure(big, 28, "bold")[0] + 0.16, ty0 + th / 2, label, size=13, weight="bold")
    gy = ty0 + th + tg
    c.box(0, gy, lw_, th, fill=T.ACCENT_TINT, edge=T.ACCENT, lw=1.4, r=0.12)
    hot = {"color": T.ACCENT_DARK, "weight": "bold"}
    c.rich(0.22, gy + th / 2, [("About", None), ("1 in 5 clinically diagnosed patients", hot),
                               ("has no identifiable pathogenic variant.", None)], lw_ - 0.44, size=12.5)
    c.text(0.05, gy + th + 0.25, "Full list of causal and modifier genes is still incomplete", size=11.5, style="italic",
           color=T.MUTED)
    src = "Forsythe & Beales 2013 · Niederlova 2019 · Perea-Romero 2022 · Monarch KG"
    c.text(0.05, c.h - 0.1, src, size=c.fit(src, lw_ - 0.1, 9, min_size=7.5), color=T.FAINT)

    # ---------------------------------------------------------------- right: what the model organisms hold
    px, pw, py, ph = 6.3, c.w - 6.3, 0.05, 3.22
    c.panel(px, py, pw, ph, title="Model organism gene–phenotype networks", title_color=T.INK)
    step_x = pw / 4
    rows = [(1.14, D.species[:4]), (2.36, D.species[4:])]
    lx = px + 3.5 * step_x - 0.38                     # legend, in the free cell at the end of the second row
    c.circle(lx, rows[1][0] - 0.13, 0.05, fill=T.MUTED, edge=None)
    c.text(lx + 0.14, rows[1][0] - 0.13, "gene", size=10, color=T.MUTED)
    c.square(lx, rows[1][0] + 0.13, 0.12, fill=T.WHITE, edge=T.MUTED, lw=1.3)
    c.text(lx + 0.14, rows[1][0] + 0.13, "phenotype", size=10, color=T.MUTED)
    lo, hi = math.log10(600), math.log10(260000)
    for cy, group in rows:
        for j, sp in enumerate(group):
            frac = (math.log10(sp.annotations) - lo) / (hi - lo)
            R = 0.2 + 0.25 * frac
            n_genes = int(round(5 + 15 * frac))
            x = px + step_x / 2 + j * step_x
            _cluster(c, x, cy, R, n_genes, max(2, n_genes // 4), sp.color, rng, node_r=0.042)
            c.text(x, cy + 0.62, sp.common, size=10.5, weight="bold", ha="center")

    # the open question: can the networks on the right say something about the disease on the left?
    qx, qy = (lw_ + px) / 2 + 0.04, top + 1.05
    c.line(px - 0.04, qy, qx, qy, color=T.CAROLINA, lw=1.6, ls=(0, (3, 2.5)))
    c.arrow(qx, qy, lw_ - 0.08, qy, color=T.CAROLINA, lw=1.6, head=0.1, ls=(0, (3, 2.5)))
    c.circle(qx, qy, 0.22, fill=T.ACCENT, edge=T.WHITE, lw=2.0, zorder=7)
    c.text(qx, qy + 0.01, "?", size=16, weight="bold", color=T.WHITE, ha="center", zorder=8)

    # how many edges of each kind the graph holds (one scale, labelled directly)
    labels = ["Non-human gene–phenotype edges", "Human gene–phenotype edges", "Human gene–disease edges"]
    x0 = px + max(c.measure(s_, 11, "bold")[0] for s_ in labels) + 0.12
    total = D.totals.model_annotations
    wmax = c.w - x0 - c.measure(n(total), 12.5, "bold")[0] - 0.14
    bh, by = 0.24, [py + ph + 0.36, py + ph + 0.72, py + ph + 1.08]
    for lab, y in zip(labels, by):
        c.text(x0 - 0.12, y, lab, size=11, weight="bold", ha="right")
    x = x0
    for sp in D.species:                              # the non-human bar is split by species, same colours as above
        w = wmax * sp.annotations / total
        c.rect(x, by[0] - bh / 2, max(w - 0.02, 0.012), bh, fill=sp.color)
        x += w
    c.text(x0 + wmax + 0.1, by[0], n(total), size=12.5, weight="bold")
    for val, y in zip([H.phenotype_annotations, H.gene_disease_pairs], by[1:]):
        w = wmax * val / total
        c.box(x0, y - bh / 2, w, bh, fill=T.NAVY, edge=None, r=0.04)
        c.text(x0 + w + 0.1, y, n(val), size=12.5, weight="bold")
    c.text(c.w, c.h - 0.2, "How can we leverage this information?", size=14, weight="bold", color=T.ACCENT_DARK, ha="right")
    return c


# =============================================================================================
# 3. What a phenolog is
# =============================================================================================
@figure("mcgary2010", "The 2010 PNAS paper by McGary and colleagues that defined phenologs as orthologous phenotypes, "
                      "what data it had, and three of its non-obvious disease models: yeast for angiogenesis, a plant "
                      "for Waardenburg syndrome and a worm for breast and ovarian cancer.")
def fig_mcgary2010(D, step=None):
    c = Canvas()
    m = D.history.mcgary2010
    y0, h = 0.1, 2.35

    # the paper
    c.box(0, y0, 3.7, h, fill=T.WHITE, edge=T.RULE, lw=1.2, r=0.12)
    c.text(0.25, y0 + 0.32, f"{m.journal} · 2010", size=11.5, weight="bold", color=T.CAROLINA_DARK)
    tt = c.wrap(m.title, 3.2, 13, "bold")
    c.text(0.25, y0 + 0.56, tt, size=13, weight="bold", va="top", linespacing=1.2)
    ay = y0 + 0.56 + 0.235 * (tt.count("\n") + 1) + 0.14
    c.text(0.25, ay, c.wrap(m.authors, 3.2, 10.5), size=10.5, color=T.MUTED, va="top", linespacing=1.2)
    c.text(0.25, y0 + h - 0.26, f"PMID {m.pmid}", size=10, color=T.FAINT)

    # the definition
    dx, dw = 3.95, 4.35
    c.panel(dx, y0, dw, h)
    c.text(dx + 0.28, y0 + 0.5, "phenolog", size=32, weight="bold")
    c.text(dx + 0.3, y0 + 1.0, "an orthologous phenotype", size=15, weight="bold", color=T.CAROLINA_DARK)
    body = c.wrap("Two phenotypes, in different species, whose genes are orthologs of each other more often than "
                  "chance allows.", dw - 0.6, 13)
    c.text(dx + 0.3, y0 + 1.3, body, size=13, va="top", linespacing=1.3)

    # what they had
    tx, tw, th, g = 8.55, 1.65, 1.1, 0.15
    tiles = [(str(m.species), "species, human to plant"), (n(m.gene_phenotype_associations), "gene–phenotype associations"),
             (m.human_diseases_label, "human diseases"), (f"> {n(m.phenologs_total)}", "phenologs at 5% FDR")]
    for i, (big, cap) in enumerate(tiles):
        r, k = divmod(i, 2)
        x, y = tx + k * (tw + g), y0 + r * (th + g)
        c.box(x, y, tw, th, fill=T.TINT, edge=None, r=0.1)
        c.text(x + 0.18, y + 0.38, big, size=c.fit(big, tw - 0.36, 21, "bold"), weight="bold")
        c.text(x + 0.18, y + 0.66, c.wrap(cap, tw - 0.32, 10.5), size=10.5, color=T.MUTED, va="top", linespacing=1.15)

    # the models nobody would have picked
    c.text(0.05, 2.8, "The non-obvious models that made the paper", size=12.5, weight="bold", color=T.MUTED)
    cw, gap, cy, ch = 3.86, 0.21, 3.02, 2.0
    grey = "#7A8794"
    for i, cl in enumerate(m.classics):
        x = i * (cw + gap)
        c.box(x, cy, cw, ch, fill=T.WHITE, edge=T.RULE, lw=1.2, r=0.12)
        rows = [(cl.a, cy + 0.4), (cl.b, cy + 0.93)]
        c.line(x + 0.38, rows[0][1], x + 0.38, rows[1][1], color=T.FAINT, lw=1.4)
        for side, y in rows:
            key = side.get("key")
            col = T.NAVY if key == "human" else (D.sp[key].color if key else grey)
            c.badge(x + 0.38, y, side.code, col, r=0.17)
            pb = c.pill_left(x + 0.68, y, side.label, size=c.fit(side.label, 2.2, 12, "bold"), weight="bold", h=0.4,
                             pad=0.18, edge=col, lw=1.5)
            c.text(pb.right + 0.12, y, side.species, size=10, color=T.MUTED, style="italic")
        c.text(x + 0.25, cy + 1.3, c.wrap(cl.outcome, cw - 0.5, 11), size=11, va="top", linespacing=1.25)
    return c


@figure("phenolog_papers", "The three papers that introduced and developed phenologs: McGary 2010 in PNAS, Cha 2012 in "
                           "PLoS Biology and Woods 2013 in BMC Bioinformatics. Below them the definition: two phenotypes "
                           "in different species that share significantly more orthologous genes than expected by chance.")
def fig_phenolog_papers(D, step=None):
    c = Canvas()
    h = D.history
    cw, gap, y0, ch = 3.86, 0.21, 0.06, 2.3
    for i, (p, year) in enumerate([(h.mcgary2010, 2010), (h.cha2012, 2012), (h.woods2013, 2013)]):
        x = i * (cw + gap)
        c.box(x, y0, cw, ch, fill=T.WHITE, edge=T.RULE, lw=1.2, r=0.12)
        c.text(x + 0.25, y0 + 0.32, f"{p.journal} · {year}", size=11.5, weight="bold", color=T.CAROLINA_DARK)
        c.text(x + 0.25, y0 + 0.58, c.wrap(p.title, cw - 0.5, 13, "bold"), size=13, weight="bold", va="top", linespacing=1.2)
        c.text(x + 0.25, y0 + ch - 0.56, p.authors, size=c.fit(p.authors, cw - 0.5, 10.5), color=T.MUTED)
        c.text(x + 0.25, y0 + ch - 0.28, f"PMID {p.pmid}", size=10, color=T.FAINT)

    # the definition, across the slide
    dy = y0 + ch + 0.3
    dh = c.h - dy - 0.06
    c.panel(0, dy, c.w, dh)
    c.text(0.42, dy + dh / 2 - 0.22, "phenolog", size=38, weight="bold")
    c.text(0.46, dy + dh / 2 + 0.4, "an orthologous phenotype", size=15, weight="bold", color=T.CAROLINA_DARK)
    tx = 3.95
    tw = c.w - tx - 0.35
    c.line(tx - 0.4, dy + 0.38, tx - 0.4, dy + dh - 0.38, color=T.RULE, lw=1.4)
    key = {"weight": "bold"}
    c.rich(tx, dy + dh / 2 - 0.3, [("Two phenotypes in", None), ("different species", key), ("that share significantly more", None),
                                  ("orthologous genes", key), ("than expected by chance.", None)], tw, size=20)
    who = f"{h.mcgary2010.cite}:"
    qy = dy + dh - 0.5
    c.text(tx, qy, who, size=11.5, weight="bold", color=T.MUTED)
    quote = f"“{h.mcgary2010.definition}”"
    qx = tx + c.measure(who, 11.5, "bold")[0] + 0.1
    c.text(qx, qy, quote, size=c.fit(quote, tx + tw - qx, 11.5, style="italic"), color=T.MUTED, style="italic")
    return c


@figure("phenolog_idea", "Left: the idea of the 2010 study as a cartoon. Four model organism gene-phenotype networks "
                         "(mouse 74,250 associations, worm 27,065, budding yeast 86,383, Arabidopsis 22,921) surround a "
                         "much smaller human gene-disease network (1,923 associations, about 300 diseases). Dashed "
                         "orthology lines join genes of each model network to human genes, and each link carries the "
                         "number of phenologs found with human disease: 3,755, 9, 154 and 48. Right: one of the 48 "
                         "human-plant phenologs. Waardenburg syndrome (3 known genes) and negative gravitropism defects "
                         "in Arabidopsis (5 genes) share 2 genes; the 3 plant-only genes give candidates, and knocking "
                         "down one of them, sec23ip, in frog embryos stops neural crest cells from migrating.")
def fig_phenolog_idea(D, step=None):
    c = Canvas()
    rng = np.random.default_rng(3)
    m = D.history.mcgary2010
    W = m.waardenburg
    dash = (0, (3, 2.2))
    cw = 7.4                                          # the cartoon's width; the example panel takes the rest
    hx, hy = cw / 2, 2.52

    def size(assoc):
        """Radius and gene count of a cartoon network, from its number of associations (log scale)."""
        f = (math.log10(assoc) - math.log10(1500)) / (math.log10(100000) - math.log10(1500))
        return 0.47 + 0.23 * f, int(round(9 + 15 * f))

    def label(x, y, code, col, name, lines, ha):
        """A network's caption: badge and name on one row, then (bold, plain) text rows. ha: where x is."""
        rows = [(0.28 + 0.08 + c.measure(name, 12, "bold")[0], None)]
        rows += [(c.measure(b, 10.5, "bold")[0] + (0.06 if b and t else 0) + c.measure(t, 9.5)[0], (b, t)) for b, t in lines]
        for i, (w, row) in enumerate(rows):
            x0 = {"left": x, "right": x - w, "center": x - w / 2}[ha]
            yy = y + i * 0.24
            if row is None:
                c.badge(x0 + 0.14, yy, code, col, r=0.14)
                c.text(x0 + 0.36, yy, name, size=12, weight="bold")
            else:
                b, t = row
                c.text(x0, yy + 0.04, b, size=10.5, weight="bold", va="baseline")
                c.text(x0 + c.measure(b, 10.5, "bold")[0] + (0.06 if b and t else 0), yy + 0.04, t, size=9.5, color=T.MUTED,
                       va="baseline")

    # ---- the human gene-disease network, small, in the middle
    hR, hn = size(m.associations.human)
    HG, _ = _cluster(c, hx, hy, hR, hn, 3, T.NAVY, rng, term_fill=T.NAVY)
    label(hx, hy + hR + 0.3, D.human.code, T.NAVY, D.human.common,
          [(n(m.associations.human), "gene–disease associations"), ("", f"{m.human_diseases_label} diseases · {m.human_source}")],
          "center")

    # ---- the four model organism networks around it, each joined to human genes by orthology
    top, bottom, cap_top, cap_bottom = 1.44, 3.72, 0.2, 4.68
    spots = [(0.95, top, "left", cap_top), (0.95, bottom, "left", cap_bottom), (cw - 0.95, bottom, "right", cap_bottom),
             (cw - 0.95, top, "right", cap_top)]
    for mo, (x, y, ha, cap_y) in zip(m.models, spots):
        col = D.sp[mo.key].color if mo.get("key") else mo.color
        hot = mo.name == W.model
        R, ng = size(m.associations[mo.name])
        G, _ = _cluster(c, x, y, R, ng, max(3, ng // 4), col, rng)
        near = sorted(G, key=lambda g: math.hypot(g[0] - hx, g[1] - hy))[:3]
        used = []
        for gx, gy in near:                           # the genes facing human, each to the human gene nearest to it
            t = min((h for h in HG if h not in used), key=lambda h: math.hypot(h[0] - gx, h[1] - gy))
            used.append(t)
            c.line(gx, gy, t[0], t[1], color=T.ACCENT_DARK if hot else T.MUTED, lw=1.3, ls=dash, zorder=2)
        d = math.hypot(hx - x, hy - y)
        ux, uy = (hx - x) / d, (hy - y) / d
        mx, my = (x + ux * R + hx - ux * hR) / 2, (y + uy * R + hy - uy * hR) / 2
        c.pill(mx, my, n(mo.phenologs_with_human), h=0.3, size=10.5, weight="bold", pad=0.12, lw=1.3, zorder=5,
               fill=T.ACCENT if hot else T.WHITE, edge=T.ACCENT if hot else T.MUTED, color=T.WHITE if hot else T.INK)
        label(0.05 if ha == "left" else cw - 0.05, cap_y, mo.code, col, mo.common,
              [(n(m.associations[mo.name]), "gene–phenotype associations")], ha)
        if hot:                                       # this comparison is the one worked through on the right
            c.arrow(x + R + 0.1, y, cw + 0.14, y, color=T.ACCENT_DARK, lw=1.6, head=0.11)

    # ---- legend, above the human network
    lx, ly = hx - 1.22, 0.2
    c.circle(lx + 0.07, ly, 0.055, fill=T.MUTED, edge=None)
    c.text(lx + 0.2, ly, "gene", size=9.5, color=T.MUTED)
    c.square(lx + 0.74, ly, 0.13, fill=T.WHITE, edge=T.MUTED, lw=1.3)
    c.text(lx + 0.87, ly, "phenotype", size=9.5, color=T.MUTED)
    c.square(lx + 1.72, ly, 0.13, fill=T.MUTED, edge=T.MUTED, lw=1.3)
    c.text(lx + 1.85, ly, "disease", size=9.5, color=T.MUTED)
    c.line(lx, ly + 0.24, lx + 0.34, ly + 0.24, color=T.MUTED, lw=1.0)
    c.text(lx + 0.44, ly + 0.24, "curated association", size=9.5, color=T.MUTED)
    c.line(lx, ly + 0.48, lx + 0.34, ly + 0.48, color=T.MUTED, lw=1.3, ls=dash)
    c.text(lx + 0.44, ly + 0.48, "orthologous genes", size=9.5, color=T.MUTED)
    c.pill(lx + 0.17, ly + 0.74, "n", h=0.2, w=0.34, size=8.5, weight="bold", lw=1.1, edge=T.MUTED)
    c.text(lx + 0.44, ly + 0.74, "phenologs with human disease", size=9.5, color=T.MUTED)

    # ---- how much more there was to borrow
    models_total = sum(m.associations[mo.name] for mo in m.models)
    ratio = f"≈ {round(models_total / m.associations.human, -1):.0f}×"
    c.text(hx, 4.42, ratio, size=20, weight="bold", ha="center")
    c.text(hx, 4.8, "more associations in the four\nmodel organisms than in human", size=9.5, color=T.MUTED, ha="center",
           linespacing=1.2)

    # ---- right: one comparison worked through (the paper's Fig. 4)
    px, pw, ph = cw + 0.2, c.w - cw - 0.2, c.h - 0.02
    plant = next(mo for mo in m.models if mo.name == W.model)
    c.panel(px, 0, pw, ph)
    c.text(px + 0.22, 0.31, f"One of the {plant.phenologs_with_human} human–plant phenologs", size=12.5, weight="bold",
           color=T.MUTED)

    # the two gene sets as a Venn: human above, plant below, the plant-only genes marked as candidates
    vx, vr, y1, y2 = px + 0.82, 0.5, 1.18, 1.86
    c.circle(vx, y2, vr, fill=T.ACCENT_TINT, edge=None, zorder=2)
    c.circle(vx, y1, vr, fill=T.WHITE, edge=T.NAVY, lw=1.8, zorder=3)
    c.circle(vx, y2, vr, fill=None, edge=plant.color, lw=1.8, zorder=4)
    c.text(vx, y1 - 0.2, str(W.human_genes - W.shared), size=13, weight="bold", ha="center")
    c.text(vx, (y1 + y2) / 2, str(W.shared), size=13, weight="bold", ha="center")
    c.text(vx, y2 + 0.22, str(W.candidates), size=13, weight="bold", color=T.ACCENT_DARK, ha="center")

    tx = px + 1.62
    tw = px + pw - 0.2 - tx

    def term(y, code, col, name, sub):
        c.badge(tx + 0.13, y, code, col, r=0.13)
        c.text(tx + 0.35, y, name, size=c.fit(name, tw - 0.35, 12.5, "bold"), weight="bold")
        c.text(tx, y + 0.25, sub, size=c.fit(sub, tw, 10), color=T.MUTED)

    term(y1 - 0.3, D.human.code, T.NAVY, W.disease, f"{W.human_genes} known genes · {W.disease_gloss}")
    c.text(tx, (y1 + y2) / 2, f"{W.shared} orthologous genes shared", size=10, style="italic", color=T.MUTED)
    c.line(vx + 0.16, (y1 + y2) / 2, tx - 0.1, (y1 + y2) / 2, color=T.FAINT, lw=1.0, zorder=5)
    term(y2 + 0.08, plant.code, plant.color, W.phenotype, f"{W.plant_genes} genes · {W.phenotype_gloss}")

    # the plant-only genes become candidates for the disease
    cy = 2.84
    c.arrow(vx, y2 + 0.4, vx, cy - 0.2, color=T.ACCENT_DARK, lw=1.6, head=0.1)
    pb = c.pill(vx + 0.12, cy, W.tested_gene, h=0.32, size=11, weight="bold", pad=0.14, fill=T.ACCENT, edge=None, color=T.WHITE)
    c.text(pb.right + 0.12, cy, f"one of {W.candidates} new candidate genes", size=10.5, color=T.ACCENT_DARK, weight="bold")

    # the test in frog embryos
    bx, by, bw = px + 0.2, 3.32, pw - 0.4
    bh = ph - by - 0.18
    c.arrow(vx, cy + 0.18, vx, by - 0.02, color=T.ACCENT_DARK, lw=1.6, head=0.1)
    c.box(bx, by, bw, bh, fill=T.WHITE, edge=T.RULE, lw=1.0, r=0.12, zorder=1)
    frog = D.sp.frog.color
    c.badge(bx + 0.33, by + 0.32, W.frog_code, frog, r=0.13)
    c.text(bx + 0.55, by + 0.32, "Tested in frog embryos", size=12, weight="bold")
    ty = by + 0.86                                    # the embryo counts are in the notes, not on the slide (round 1)
    c.text(bx + 0.2, ty, W.frog_gene, size=11.5, style="italic")
    c.text(bx + 0.2 + c.measure(W.frog_gene, 11.5, style="italic")[0] + 0.05, ty, "knockdown:", size=11.5)
    c.text(bx + 0.2, ty + 0.27, "neural crest cells fail to migrate", size=11.5)

    # two sides of one embryo, schematic (after the paper's Fig. 4E): a tailbud embryo in side view, head to the
    # left. Streams of neural crest cells leave the back and move down; on the knockdown side they barely start.
    def embryo(x0, yc, w, h, reach):
        pts = [(0, 0), (0, -0.62), (0.34, -0.56), (0.55, -0.40), (0.76, -0.24), (0.90, -0.20), (1.0, -0.02),
               (0.90, 0.26), (0.70, 0.50), (0.45, 0.52), (0.20, 0.55), (0, 0.50), (0, 0)]
        verts = [(x0 + a * w, yc + b * h) for a, b in pts]
        codes = [MplPath.MOVETO] + [MplPath.CURVE4] * 12
        c.ax.add_patch(PathPatch(MplPath(verts, codes), facecolor="#F3EEDD", edgecolor="#B9AE8C", lw=1.2, zorder=2,
                                 joinstyle="round"))
        c.circle(x0 + 0.1 * w, yc - 0.08 * h, 0.022, fill="#B9AE8C", edge=None, zorder=3)
        for k in range(3):
            sx = x0 + (0.22 + 0.11 * k) * w
            y0 = yc - 0.36 * h                         # just inside the back
            c.line(sx, y0, sx - 0.02, y0 + reach * h, color="#5A4A9C", lw=3.0, zorder=3)

    ew, eh = 1.12, 0.42
    ex = bx + bw - ew - 0.2
    for i, (cap, reach) in enumerate([("control side", 0.62), ("knockdown side", 0.12)]):
        yc = by + 0.38 + i * 0.7
        embryo(ex, yc, ew, eh, reach)
        c.text(ex + ew / 2, yc + 0.35, cap, size=8.5, color=T.MUTED, ha="center")
    return c


MONO = "Menlo"          # Biolink CURIEs are set in a monospace face (figures are placed as images, so any local font works)


@figure("kg_components", "The three parts of the Monarch knowledge graph that a phenolog calculation uses, in three "
                         "columns. Human: a gene-disease network (12,426 edges, 8,040 Mondo diseases) and a "
                         "gene-phenotype network (131,495 edges, 8,750 HPO phenotypes), genes as HGNC ids. Model "
                         "organisms: one gene-phenotype network for each of seven species, with its phenotype "
                         "ontology, gene id source, term count and edge count. Gene equivalence: PANTHER orthology "
                         "edges, 19,159 human genes in 7,489 protein families, and the families each species shares "
                         "with human. Under each column, the Biolink node and edge categories that select it.")
def fig_kg_components(D, step=None):
    c = Canvas()
    H, K, B = D.human, D.kg, D.kg.biolink
    MODEL = T.CAROLINA_DARK                       # a model organism in general (the species keep their own colours)
    ph = c.h - 0.02
    x1, w1, x2, w2, x3, w3 = 0.0, 3.8, 3.95, 4.65, 8.75, 3.25
    ys, by = 0.9, 3.3                             # the row of schematic nodes; the top of the Biolink row
    bh = ph - by - 0.12

    def gene(x, y, color, r=0.11):
        c.circle(x, y, r, fill=color, edge=None, zorder=4)

    def term(x, y, color, solid, s=0.23):
        c.square(x, y, s, fill=color if solid else T.WHITE, edge=color, lw=1.6, zorder=4)

    def edge(xa, xb, color, label, ls="-", lift=0.19):
        c.line(xa, ys, xb, ys, color=color, lw=1.6, ls=ls, zorder=2)
        c.text((xa + xb) / 2, ys - lift, label, size=9.5, color=T.MUTED, ha="center")

    def under(x, s):
        c.text(x, ys + 0.3, s, size=10, color=T.MUTED, ha="center")

    def band(x, w, rows):
        """The Biolink row of a column: a white box of (glyph, label, CURIE) lines. A glyph is a function that
        draws the schematic's own symbol at (x, y), so the row doubles as the schematic's legend. The heading
        says Biolink, so the lines drop the "biolink:" prefix that data.yaml keeps."""
        bx = x + 0.12
        c.box(bx, by, w - 0.24, bh, fill=T.WHITE, edge=T.RULE, lw=1.0, r=0.1, zorder=1)
        c.text(bx + 0.18, by + 0.25, "Biolink categories", size=10.5, weight="bold", color=T.CAROLINA_DARK)
        y = by + 0.52
        for glyph, label, curie in rows:
            tx = bx + 0.5
            if glyph:
                glyph(bx + 0.28, y)
            if label:
                c.text(tx, y, label, size=9, color=T.MUTED, style="italic")
                tx += c.measure(label, 9, style="italic")[0] + 0.07
            c.text(tx, y + 0.005, curie.removeprefix("biolink:"), size=9.5, family=MONO, color=T.INK if glyph else T.MUTED)
            y += 0.186

    def g_node(color):
        return lambda x, y: gene(x, y, color, r=0.06)

    def t_node(color, solid):
        return lambda x, y: term(x, y, color, solid, s=0.125)

    def e_line(color, ls="-"):
        return lambda x, y: c.line(x - 0.1, y, x + 0.1, y, color=color, lw=1.6, ls=ls, zorder=4)

    # ---- human: two networks around one set of genes
    c.panel(x1, 0, w1, ph)
    c.badge(x1 + 0.36, 0.31, H.code, T.NAVY, r=0.14)
    c.text(x1 + 0.6, 0.31, "Human networks", size=12.5, weight="bold", color=T.MUTED)
    dx, gx, px = x1 + 0.8, x1 + w1 / 2, x1 + w1 - 0.8
    edge(dx, gx, T.NAVY, "gene–disease")
    edge(gx, px, T.NAVY, "gene–phenotype")
    term(dx, ys, T.NAVY, True)
    gene(gx, ys, T.NAVY)
    term(px, ys, T.NAVY, False)
    under(dx, "disease")
    under(gx, "gene")
    under(px, "phenotype")
    c.line(gx, 1.46, gx, 3.04, color=T.RULE, lw=1.2)
    for sx, head, edges, what, count, onto in [
            (x1 + 0.28, "Gene–disease", H.gene_disease_pairs, "diseases", H.diseases_with_gene, K.human.disease_ontology),
            (gx + 0.22, "Gene–phenotype", H.phenotype_annotations, "phenotypes", H.phenotype_terms, K.human.phenotype_ontology)]:
        c.text(sx, 1.58, head, size=12, weight="bold")
        num = n(edges)
        c.text(sx, 1.95, num, size=17, weight="bold", va="baseline")
        c.text(sx + c.measure(num, 17, "bold")[0] + 0.07, 1.95, "edges", size=10.5, color=T.MUTED, va="baseline")
        for i, (k, v) in enumerate([(what, n(count)), ("ontology", onto), ("gene ids", K.human.gene_ids)]):
            y = 2.33 + i * 0.26
            c.text(sx, y, k, size=10.5, color=T.MUTED)
            c.text(sx + 0.92, y, v, size=10.5, weight="bold")
    band(x1, w1, [(g_node(T.NAVY), None, B.gene), (t_node(T.NAVY, True), None, B.disease),
                  (t_node(T.NAVY, False), None, B.phenotype)] +
                 [(e_line(T.NAVY), None, cat) for cat in list(B.gene_disease) + [B.gene_phenotype]])

    # ---- model organisms: one gene-phenotype network each
    c.panel(x2, 0, w2, ph)
    c.text(x2 + 0.2, 0.31, "Model organism networks", size=12.5, weight="bold", color=T.MUTED)
    mg, mp = x2 + 0.75, x2 + 2.05
    edge(mg, mp, MODEL, "gene–phenotype")
    gene(mg, ys, MODEL)
    term(mp, ys, MODEL, False)
    under(mg, "gene")
    under(mp, "phenotype")
    tot = n(D.totals.model_annotations)
    c.text(x2 + 2.85, ys + 0.07, tot, size=17, weight="bold", va="baseline")
    c.text(x2 + 2.85 + c.measure(tot, 17, "bold")[0] + 0.07, ys + 0.07, "edges", size=10.5, color=T.MUTED, va="baseline")
    c.text(x2 + 2.85, ys + 0.3, f"in {D.totals.model_species} species", size=10, color=T.MUTED)
    cols = [("ontology", x2 + 1.58, "left"), ("gene ids", x2 + 2.64, "left"), ("terms", x2 + 3.76, "right"), ("edges", x2 + w2 - 0.22, "right")]
    hy, ry, rh = 1.5, 1.77, 0.225
    for head, cx, ha in cols:
        c.text(cx, hy, head, size=10, color=T.MUTED, ha=ha)
    c.line(x2 + 0.2, hy + 0.13, x2 + w2 - 0.2, hy + 0.13, color=T.RULE, lw=1.0)
    for i, sp in enumerate(D.species):
        y = ry + i * rh
        c.badge(x2 + 0.31, y, sp.code, sp.color, r=0.1)
        c.text(x2 + 0.48, y, sp.common, size=11, weight="bold")
        for (head, cx, ha), v in zip(cols, [sp.ontology, K.gene_ids[sp.key], n(sp.phenotype_terms), n(sp.annotations)]):
            c.text(cx, y, v, size=10.5, ha=ha)
    band(x2, w2, [(g_node(MODEL), None, B.gene), (t_node(MODEL, False), None, B.phenotype),
                  (e_line(MODEL), None, B.gene_phenotype), (None, "predicate", B.gene_phenotype_predicate)])

    # ---- gene equivalence: orthology edges that carry a PANTHER family
    c.panel(x3, 0, w3, ph)
    c.text(x3 + 0.2, 0.31, "Gene equivalence", size=12.5, weight="bold", color=T.MUTED)
    ha_, hb_ = x3 + 0.62, x3 + w3 - 0.62
    dash = (0, (3, 2.2))
    edge(ha_, hb_, T.MUTED, "gene–gene orthology", ls=dash, lift=0.31)
    gene(ha_, ys, T.NAVY)
    gene(hb_, ys, MODEL)
    c.pill((ha_ + hb_) / 2, ys, f"{K.orthology.source} family id", h=0.27, size=9.5, edge=T.MUTED, lw=1.1, pad=0.1, zorder=3)
    under(ha_, "human gene")
    under(hb_, "model gene")
    nx, fx = x3 + 0.25, x3 + 1.9
    c.text(nx, 1.7, n(H.genes_with_ortholog), size=17, weight="bold", va="baseline")
    c.text(nx, 1.81, "human genes with\nan ortholog", size=10, color=T.MUTED, va="top", linespacing=1.2)
    c.arrow(x3 + 1.42, 1.62, x3 + 1.74, 1.62, color=T.MUTED, lw=1.5, head=0.1)
    c.text(fx, 1.7, n(H.families), size=17, weight="bold", va="baseline")
    c.text(fx, 1.81, f"{K.orthology.source} protein\nfamilies", size=10, color=T.MUTED, va="top", linespacing=1.2)
    c.text(nx, 2.42, "Families shared with human", size=10, color=T.MUTED)
    c.line(x3 + 0.2, 2.55, x3 + w3 - 0.2, 2.55, color=T.RULE, lw=1.0)
    cell = (w3 - 0.4) / 4
    for i, sp in enumerate(D.species):
        cx_, y = x3 + 0.2 + (i % 4) * cell, 2.77 + (i // 4) * 0.29
        c.badge(cx_ + 0.11, y, sp.code, sp.color, r=0.1)
        c.text(cx_ + 0.27, y, n(sp.families_shared), size=10.5)
    band(x3, w3, [(g_node(T.NAVY), None, B.gene), (e_line(T.MUTED, dash), None, B.orthology),
                  (None, "predicate", B.orthology_predicate), (None, "evidence", B.orthology_evidence)])
    return c


@figure("phenolog_concept", "Top: the human gene-phenotype network times the mouse gene-phenotype network gives a grid "
                            "with a row for each of 8,750 human phenotypes and a column for each of 10,998 mouse "
                            "phenotypes, about 96 million pairs to test; one cell, rod-cone dystrophy against retina "
                            "degeneration, is ringed. Bottom: that test. The 152 genes of rod-cone dystrophy and the 115 "
                            "genes of mouse retina degeneration are collapsed to PANTHER families; 40 are shared, "
                            "P = 3.8e-45, so the two are phenologs, and the 57 families seen only in mouse, 58 genes, "
                            "become candidates for the human phenotype by guilt by association.")
def fig_phenolog_concept(D, step=None):
    c = Canvas()
    rng = np.random.default_rng(4)
    H, mouse, cm = D.human, D.sp.mouse, D.comparisons

    # ------------------------------------------------------------------ a) one network against the other
    d, rc = D.rod_cone.derivation, D.rod_cone
    ah = 2.3
    c.panel(0, 0, c.w, ah, tag="a)", title="Test every human phenotype term against every mouse phenotype term")
    cy, R = 1.3, 0.5

    def net(x, code, color, name, n_genes, count):
        _cluster(c, x, cy, R, n_genes, max(3, n_genes // 4), color, rng, node_r=0.048)
        tx = x + R + 0.24
        c.badge(tx + 0.15, cy - 0.31, code, color, r=0.15)
        c.text(tx + 0.4, cy - 0.31, name, size=13.5, weight="bold")
        c.text(tx, cy - 0.01, "gene–phenotype network", size=10.5, color=T.MUTED)
        c.text(tx, cy + 0.25, f"{n(count)} phenotypes", size=11.5, weight="bold")

    net(0.8, H.code, T.NAVY, H.common, 20, H.phenotype_terms)
    c.text(3.52, cy, "×", size=26, weight="bold", color=T.CAROLINA_DARK, ha="center")
    net(4.4, mouse.code, mouse.color, mouse.common, 23, cm.mouse.phenotypes)
    c.text(7.12, cy, "=", size=26, weight="bold", color=T.CAROLINA_DARK, ha="center")

    # the grid: a row per human phenotype, a column per mouse phenotype, a cell per pair. A few cells are
    # significant; the ringed one is the pair walked through in panel b.
    gx0, gyt, pitch, cell, ncol, nrow = 8.82, 0.98, 0.132, 0.1, 13, 7
    ri, ci = 4, 8
    hot = {(0, 2), (1, 10), (2, 5), (3, 12), (5, 0), (6, 6)}
    for r_ in range(nrow):
        for k in range(ncol):
            fill = T.ACCENT if (r_, k) in hot or (r_, k) == (ri, ci) else "#A9B8C8" if r_ == ri or k == ci else "#CBD6E0"
            c.rect(gx0 + k * pitch, gyt + r_ * pitch, cell, cell, fill=fill)
    c.box(gx0 + ci * pitch - 0.028, gyt + ri * pitch - 0.028, cell + 0.056, cell + 0.056, fill=None, edge=T.NAVY, lw=1.6,
          r=0.02, zorder=5)
    gw, ghh = ncol * pitch - (pitch - cell), nrow * pitch - (pitch - cell)
    colx, rowy = gx0 + ci * pitch + cell / 2, gyt + ri * pitch + cell / 2
    c.text(gx0 + gw / 2, gyt - 0.43, f"{n(cm.mouse.phenotypes)} mouse phenotypes", size=10.5, ha="center")
    c.text(colx, gyt - 0.19, d.phenotype, size=9.5, color=T.MUTED, style="italic", ha="center")
    c.line(colx, gyt - 0.1, colx, gyt - 0.03, color=T.MUTED, lw=1.0)
    c.text(gx0 - 0.13, gyt + 0.05, n(H.phenotype_terms), size=10.5, ha="right")
    c.text(gx0 - 0.13, gyt + 0.23, "human phenotypes", size=10.5, ha="right")
    c.text(gx0 - 0.13, rowy, rc.label, size=9.5, color=T.MUTED, style="italic", ha="right")
    c.line(gx0 - 0.1, rowy, gx0 - 0.03, rowy, color=T.MUTED, lw=1.0)
    ly = gyt + ghh + 0.2
    c.rect(gx0, ly - 0.05, cell, cell, fill=T.ACCENT)
    c.text(gx0 + 0.17, ly, "overlap beyond chance: a phenolog", size=9.5, color=T.MUTED)
    # human count: HPO terms on a gene with an ortholog in any of the seven species, so the product is a slight
    # upper bound for mouse alone. Hence the "≈".
    pairs = H.phenotype_terms * cm.mouse.phenotypes
    nx = gx0 + gw + 0.16
    big = f"≈ {pairs / 1e6:.0f} million"
    c.text(nx, gyt + 0.3, big, size=c.fit(big, c.w - nx - 0.1, 15, "bold"), weight="bold")
    c.text(nx, gyt + 0.55, "pairs to test", size=10, color=T.MUTED)

    # ------------------------------------------------------------------ b) one of those tests
    b0 = ah + 0.14
    c.panel(0, b0, c.w, c.h - b0, tag="b)", title="Example comparison")
    gx, fx, tx0, tw, px, kx = 2.78, 4.95, 6.2, 1.95, 9.28, 11.05
    for x, head in [(gx, "Genes"), (fx, "PANTHER protein families"), (tx0 + tw / 2, "Hypergeometric test"),
                    (px, "Phenolog"), (kx, "Guilt by association")]:
        c.text(x, b0 + 0.28, head, size=11.5, weight="bold", color=T.MUTED, ha="center")

    y1, y2 = b0 + 0.98, b0 + 1.86
    ym = (y1 + y2) / 2
    c.badge(0.3, y1, H.code, T.NAVY, r=0.17)
    c.text(0.56, y1, rc.label.replace(" ", "\n", 1), size=12.5, weight="bold")
    c.badge(0.3, y2, mouse.code, mouse.color, r=0.17)
    c.text(0.56, y2, d.phenotype.replace(" ", "\n", 1), size=12.5, weight="bold")
    fam = []
    for y, genes, eg, fams in [(y1, d.genes, d.gene_examples, d.families),
                               (y2, d.phenotype_genes, d.phenotype_gene_examples, d.phenotype_families)]:
        g = c.pill(gx, y, f"{n(genes)} genes", sub=f"e.g. {eg}", size=13, weight="bold", sub_size=10.5, sub_color=T.MUTED,
                   sub_style="italic", w=1.95, edge=T.RULE, lw=1.3)
        f = c.pill(fx, y, f"{n(fams)} families", sub=f"e.g. {d.family_example}", size=13, weight="bold", sub_size=10.5,
                   sub_color="#C9DDF0", w=1.55, fill=T.NAVY, edge=None, color=T.WHITE)
        c.arrow(g.right + 0.07, y, f.x - 0.07, y, color=T.FAINT, lw=1.5)
        fam.append(f)
    bx = fam[0].right + 0.1
    for y in (y1, y2):
        c.line(bx, y, bx + 0.16, y, color=T.NAVY, lw=1.6)
    c.line(bx + 0.16, y1, bx + 0.16, y2, color=T.NAVY, lw=1.6)
    tb = c.box(tx0, y1 - 0.42, tw, y2 - y1 + 0.78, fill=T.WHITE, edge=T.RULE, lw=1.3, r=0.12)
    c.arrow(bx + 0.16, ym, tb.x - 0.05, ym, color=T.FAINT, lw=1.5)

    # the overlap: shared in blue, the mouse-only part in the accent colour (it is where the new candidates come from)
    vy, r, off = y1 + 0.2, 0.5, 0.28
    hc = Circle((tb.cx - off, vy), r, facecolor=T.TINT_2, edgecolor="none", zorder=3)
    mc = Circle((tb.cx + off, vy), r, facecolor=T.ACCENT_TINT, edgecolor="none", zorder=3)
    both = Circle((tb.cx + off, vy), r, facecolor="#A9CCE8", edgecolor="none", zorder=3.5)
    for patch in (hc, mc, both):
        c.ax.add_patch(patch)
    both.set_clip_path(hc)
    c.circle(tb.cx - off, vy, r, fill=None, edge=T.NAVY, lw=1.5, zorder=4)
    c.circle(tb.cx + off, vy, r, fill=None, edge=T.ACCENT, lw=1.5, zorder=4)
    c.badge(tb.x + 0.2, tb.y + 0.2, H.code, T.NAVY, r=0.11)
    c.badge(tb.right - 0.2, tb.y + 0.2, mouse.code, mouse.color, r=0.11)
    only_m = d.phenotype_families - d.shared
    c.text(tb.cx - 0.53, vy, n(d.families - d.shared), size=12, color=T.MUTED, ha="center")
    c.text(tb.cx, vy, n(d.shared), size=16, weight="bold", ha="center")
    c.text(tb.cx + 0.53, vy, n(only_m), size=13, weight="bold", color=T.ACCENT_DARK, ha="center")
    c.text(tb.cx, tb.bottom - 0.27, f"of {n(d.N)} families\nin both species", size=10, color=T.MUTED, ha="center")
    c.text(tb.cx, tb.bottom + 0.23, f"P = {sci(d.p, bold=True)}", size=14, weight="bold", ha="center")

    # the phenolog, and what it buys: the genes of one term become candidates for the other
    pw = 1.6
    c.arrow(tb.right + 0.06, ym, px - pw / 2 - 0.08, ym, color=T.FAINT, lw=1.5)
    top = c.pill(px, y1, rc.label, size=11, w=pw, h=0.4, lw=1.8)
    bot = c.pill(px, y2, d.phenotype, size=11, w=pw, h=0.4, lw=1.8, edge=mouse.color)
    c.line(px, top.bottom, px, ym - 0.13, color=T.NAVY, lw=1.8)
    c.line(px, ym + 0.13, px, bot.y, color=T.NAVY, lw=1.8)
    c.text(px, ym, "phenolog", size=12.5, weight="bold", ha="center")
    c.ccurve(bot.right + 0.05, y2, top.right + 0.07, y1, bulge=0.62, color=T.ACCENT, lw=2.2, head=0.13, zorder=4)
    kx0 = top.right + 0.62
    c.text(kx0, ym - 0.36, f"{n(only_m)} families\n→ {n(d.mouse_only_genes)} genes", size=12.5, weight="bold",
           color=T.ACCENT_DARK, linespacing=1.2)
    c.text(kx0, ym + 0.22, "seen only in mouse\nbecome candidates\nfor the human term", size=10, va="center", linespacing=1.25)
    c.text(top.x, tb.bottom + 0.23, f"The reverse holds for the {n(d.families - d.shared)} human-only families.",
           size=10, color=T.MUTED)
    return c


@figure("gene_equivalence", "Four cases of how a human gene and a model-organism gene can correspond: one to one, one to "
                            "many, many to one and many to many. In every case all the genes belong to one PANTHER "
                            "protein family, so each side is counted as a single family.")
def fig_gene_equivalence(D, step=None):
    c = Canvas()
    H = D.human
    cases = [("One to one", "the simple case", "mouse", 1, 1),
             ("One to many", "duplicated in the model lineage", "zebrafish", 1, 2),
             ("Many to one", "duplicated in the human lineage", "worm", 2, 1),
             ("Many to many", "paralogs on both sides", "pombe", 2, 2)]
    pw, gap, py, ph = 2.82, 0.24, 0.05, 3.38
    for i, (title, sub, key, nh, nm) in enumerate(cases):
        sp = D.sp[key]
        px = i * (pw + gap)
        c.panel(px, py, pw, ph)
        c.text(px + 0.22, py + 0.3, title, size=13.5, weight="bold")
        c.text(px + 0.22, py + 0.56, sub, size=10.5, color=T.MUTED)
        hx, mx, gy = px + 0.78, px + pw - 0.78, 1.6
        c.text(hx, 0.92, "human", size=10.5, weight="bold", color=T.MUTED, ha="center")
        c.text(mx, 0.92, sp.common.lower(), size=10.5, weight="bold", color=T.MUTED, ha="center")
        ys_h = [gy] if nh == 1 else [gy - 0.27, gy + 0.27]
        ys_m = [gy] if nm == 1 else [gy - 0.27, gy + 0.27]
        c.box(px + 0.3, gy - 0.52, pw - 0.6, 1.04, fill=T.WHITE, edge=T.FAINT, lw=1.1, r=0.2, ls=(0, (3, 2.5)), zorder=1)
        for yh in ys_h:
            for ym in ys_m:
                c.line(hx + 0.19, yh, mx - 0.19, ym, color=T.CAROLINA, lw=1.4, ls=(0, (3, 2.5)), zorder=2)
        for j, yh in enumerate(ys_h):
            c.circle(hx, yh, 0.19, fill=T.NAVY, edge=T.WHITE, lw=1.5, zorder=4)
            c.text(hx, yh + 0.005, "A" if j == 0 else "A′", size=10.5, weight="bold", color=T.WHITE, ha="center", zorder=5)
        for j, ym in enumerate(ys_m):
            c.circle(mx, ym, 0.19, fill=sp.color, edge=T.WHITE, lw=1.5, zorder=4)
            c.text(mx, ym + 0.005, "a" if j == 0 else "a′", size=10.5, weight="bold", color=T.WHITE, ha="center", zorder=5)
        c.text(px + pw / 2, gy + 0.68, "one PANTHER family", size=9.5, color=T.MUTED, ha="center")
        c.arrow(px + pw / 2, gy + 0.84, px + pw / 2, gy + 1.08, color=T.FAINT, lw=1.5, head=0.1)
        ty, ts = gy + 1.36, 0.36
        c.line(px + pw / 2 - 0.45, ty, px + pw / 2 + 0.45, ty, color=T.NAVY, lw=1.8, zorder=2)
        for tx_, col in [(px + pw / 2 - 0.45, T.NAVY), (px + pw / 2 + 0.45, sp.color)]:
            c.square(tx_, ty, ts, fill=col, edge=T.WHITE, lw=1.5, r=0.07, zorder=4)
            c.text(tx_, ty + 0.005, "F", size=11, weight="bold", color=T.WHITE, ha="center", zorder=5)

    # what the roll-up buys
    by = 3.66
    c.box(0, by, 3.3, 1.36, fill=T.NAVY, edge=None, r=0.12)
    c.text(0.25, by + 0.46, f"{n(H.genes_with_ortholog)}  →  {n(H.families)}", size=22, weight="bold", color=T.WHITE)
    c.text(0.25, by + 0.8, "human genes with an ortholog,\nrolled up to PANTHER families", size=11, color="#C9DDF0",
           va="top", linespacing=1.25)
    notes = [("One unit of comparison", "Every gene is swapped for its PANTHER family, so a phenotype's genes become "
                                        "a set of families in every species."),
             ("One gene, one family", "PANTHER builds a gene tree and an HMM per family and places each gene in "
                                      "exactly one family and subfamily."),
             ("Gene detail comes back", "Families are only for finding phenologs. Candidates are still ranked gene "
                                        "by gene.")]
    nx, nw = 3.62, 2.82
    for k, (head, body) in enumerate(notes):
        x = nx + k * nw
        c.circle(x + 0.07, by + 0.24, 0.06, fill=T.CAROLINA, edge=None)
        c.text(x + 0.24, by + 0.22, head, size=12.5, weight="bold")
        c.text(x + 0.24, by + 0.45, c.wrap(body, nw - 0.5, 11), size=11, color=T.MUTED, va="top", linespacing=1.25)
    return c


@figure("phenolog_cartoon", "Two small networks, a human phenotype and a mouse phenotype each with their genes. Dashed "
                            "lines join the orthologous genes they share; the overlap is tested with a hypergeometric "
                            "test, and the genes not shared become candidates.")
def fig_phenolog_cartoon(D, step=3):
    step = step or 3
    c = Canvas()
    rc, ex, mouse = D.rod_cone, D.rod_cone.examples[0], D.sp.mouse
    L, R, ring, nr = (2.6, 2.05), (8.9, 2.05), 1.4, 0.21
    shared_l, only_l = [-35, 0, 35], [112, 180, 248]
    shared_r, only_r = [215, 180, 145], [68, 0, -68]

    def pos(center, deg):
        a = math.radians(deg)
        return center[0] + ring * math.cos(a), center[1] + ring * math.sin(a)

    c.badge(L[0] - 1.72, 0.22, "Hs", T.NAVY, r=0.17)
    c.text(L[0] - 1.45, 0.22, "Human gene–phenotype network", size=13, weight="bold")
    c.badge(R[0] - 1.72, 0.22, mouse.code, mouse.color, r=0.17)
    c.text(R[0] - 1.45, 0.22, "Mouse gene–phenotype network", size=13, weight="bold")

    hot = step >= 2
    if hot:
        for a, b in zip(shared_l, shared_r):
            (x0, y0), (x1, y1) = pos(L, a), pos(R, b)
            c.line(x0 + nr + 0.04, y0, x1 - nr - 0.04, y1, color=T.CAROLINA, lw=2.4, ls=(0, (4, 3)), zorder=1)
        c.text(6.0 - 0.25, 0.86, "orthologs", size=12.5, weight="bold", color=T.CAROLINA_DARK, ha="center")
        c.text(6.0 - 0.25, 1.08, "same PANTHER family", size=10.5, color=T.MUTED, ha="center")
    else:                       # step 1: nothing joins the two networks yet. The same "?" as on the fill_gaps slide.
        c.circle(6.0 - 0.25, L[1], 0.27, fill=T.ACCENT, edge=T.WHITE, lw=2.0, zorder=7)
        c.text(6.0 - 0.25, L[1] + 0.01, "?", size=19, weight="bold", color=T.WHITE, ha="center", zorder=8)

    def node(center, deg, label, color, shared):
        x, y = pos(center, deg)
        c.line(center[0], center[1], x, y, color=color, lw=1.3, alpha=0.55, zorder=1)
        if not hot or shared:
            c.circle(x, y, nr, fill=color, edge=None, zorder=4)
            c.text(x, y + 0.01, label, size=12, weight="bold", color=T.WHITE, ha="center", zorder=7)
        elif step >= 3:
            c.circle(x, y, nr + 0.02, fill=T.ACCENT_TINT, edge=T.ACCENT, lw=2.4, zorder=4)
            c.text(x, y + 0.01, label, size=12, weight="bold", color=T.ACCENT, ha="center", zorder=7)
        else:
            c.circle(x, y, nr, fill=T.WHITE, edge=color, lw=1.6, zorder=4)
            c.text(x, y + 0.01, label, size=12, weight="bold", color=color, ha="center", zorder=7)

    for deg, lab in zip(shared_l, "ABC"):
        node(L, deg, lab, T.NAVY, True)
    for deg, lab in zip(only_l, "DEF"):
        node(L, deg, lab, T.NAVY, False)
    for deg, lab in zip(shared_r, "ABC"):
        node(R, deg, lab, mouse.color, True)
    for deg, lab in zip(only_r, "XYZ"):
        node(R, deg, lab, mouse.color, False)

    c.pill(L[0], L[1], rc.label, sub=f"human phenotype · {rc.id}", size=13, weight="bold", lw=2.0, sub_size=9,
           sub_color=T.MUTED, w=2.28, zorder=5)
    c.pill(R[0], R[1], ex.label, sub=f"mouse phenotype · {ex.id}", size=13, weight="bold", lw=2.0, edge=mouse.color,
           sub_size=9, sub_color=T.MUTED, w=2.28, zorder=5)

    if step >= 3:
        x, y = pos(R, 0)
        c.text(x + nr + 0.18, y, "candidates for\nthe human\nphenotype", size=11.5, weight="bold", color=T.ACCENT)

    # bottom band
    if step == 1:
        c.text(6.0, 4.3, "Two species, two vocabularies, two lists of genes. Which genes are the same gene?",
               size=15, color=T.MUTED, ha="center")
        return c
    chips = [(n(ex.m), "human families"), (n(ex.n), "mouse families"), (n(ex.c), "shared"), (n(ex.N), "in both species")]
    cw, gap, arrow_w, res_w = 1.62, 0.14, 0.55, 2.5
    total = len(chips) * cw + (len(chips) - 1) * gap + arrow_w + res_w
    x = (c.w - total) / 2
    y0 = 3.86
    for i, (num, cap) in enumerate(chips):
        is_shared = cap == "shared"
        c.box(x, y0, cw, 0.8, fill=T.TINT_2 if is_shared else T.TINT, edge=T.CAROLINA if is_shared else None, lw=1.6, r=0.1)
        c.text(x + cw / 2, y0 + 0.29, num, size=19, weight="bold", ha="center")
        c.text(x + cw / 2, y0 + 0.61, cap, size=10.5, color=T.MUTED, ha="center")
        x += cw + gap
    c.arrow(x - gap + 0.1, y0 + 0.4, x - gap + arrow_w - 0.06, y0 + 0.4, color=T.MUTED, lw=1.8, head=0.13)
    x += arrow_w - gap
    c.box(x, y0, res_w, 0.8, fill=T.NAVY, edge=None, r=0.1)
    c.text(x + res_w / 2, y0 + 0.29, f"P = {sci(ex.p, bold=True)}", size=17, weight="bold", color=T.WHITE, ha="center")
    c.text(x + res_w / 2, y0 + 0.61, "hypergeometric: phenologs", size=10.5, color="#C9DDF0", ha="center")
    if step >= 3:
        c.text(6.0, 4.9, f"The {ex.n - ex.c} mouse families not yet tied to the human phenotype become candidates for it, "
                         "and the reverse.", size=12.5, color=T.ACCENT, weight="bold", ha="center")
    return c


def _nonobvious(D, ex_block):
    """One human term on the left, four of its phenologs on the right: the expected one, then the surprises."""
    c = Canvas()
    rc = ex_block
    hx, hy, hw = 1.6, 2.2, 2.9
    lab = c.wrap(rc.label, hw - 0.5, 14, "bold")
    nl = lab.count("\n") + 1
    hh = 0.56 + 0.27 * nl
    c.badge(hx, hy - hh / 2 - 0.34, "Hs", T.NAVY, r=0.19)
    c.box(hx - hw / 2, hy - hh / 2, hw, hh, fill=T.WHITE, edge=T.NAVY, lw=2.0, r=0.3)
    c.text(hx, hy - 0.13, lab, size=14, weight="bold", ha="center", linespacing=1.15)
    c.text(hx, hy + hh / 2 - 0.2, f"{n(rc.known_genes)} genes · {rc.id}", size=10.5, color=T.MUTED, ha="center")
    c.text(hx, hy + hh / 2 + 0.28, f"{n(rc.phenologs)} phenologs in {len(rc.by_species)} species", size=11.5,
           color=T.MUTED, ha="center")

    ys = [0.55, 1.62, 2.69, 3.76]
    max_w = 3.1
    for ex, y in zip(rc.examples, ys):
        sp = D.sp[ex.species]
        c.scurve(hx + hw / 2, hy, 3.55, y, color=sp.color, lw=1.8, bend=0.55)
        c.badge(3.78, y, sp.code, sp.color, r=0.2)
        edge, fill = (T.NAVY, T.WHITE) if ex.obvious else (T.ACCENT, T.ACCENT_TINT)
        if c.measure(ex.label, 13.5, "bold")[0] <= max_w:
            pb = c.pill_left(4.12, y, ex.label, size=13.5, edge=edge, fill=fill, lw=1.8, weight="bold")
            cap_y = y + 0.42
        else:                                          # a long ontology label: two lines, never shortened
            wrapped = c.wrap(ex.label, max_w, 12.5, "bold")
            w = max(c.measure(line, 12.5, "bold")[0] for line in wrapped.split("\n")) + 0.48
            pb = c.pill(4.12 + w / 2, y, wrapped, w=w, h=0.64, size=12.5, edge=edge, fill=fill, lw=1.8, weight="bold",
                        r=0.24)
            cap_y = y + 0.5
        tag, tcol = ("as expected", T.MUTED) if ex.obvious else ("non-obvious", T.ACCENT_DARK)
        c.text(pb.right + 0.18, y, tag, size=11.5, color=tcol, weight="bold", style="italic")
        c.text(4.18, cap_y, f"{sp.common} · {ex.c} shared families · P = {sci(ex.p)}", size=11, color=T.MUTED)

    kx, ky, kw, kh = 9.2, 0.55, 2.78, 3.3
    c.box(kx, ky, kw, kh, fill=T.NAVY, edge=None, r=0.14)
    c.text(kx + 0.28, ky + 0.62, "Same genes,\ndifferent readout", size=17, weight="bold", color=T.WHITE, linespacing=1.15)
    body = c.wrap("Molecular function is conserved across species. The organism-level phenotype is not: it is an "
                  "emergent property of the system those genes sit in.", kw - 0.56, 12.5)
    c.text(kx + 0.28, ky + 1.22, body, size=12.5, color="#DCE8F5", va="top", linespacing=1.3)
    c.text(kx + 0.28, ky + kh - 0.3, "after McGary et al. 2010", size=10, color=T.TINT_2, style="italic")

    c.text(0.15, 4.72, rc.punchline, size=c.fit(rc.punchline, 11.7, 13.5, "bold"), color=T.INK, weight="bold")
    return c


@figure("nonobvious", "Rod-cone dystrophy in human linked to four model phenotypes: mouse retina degeneration as "
                      "expected, then zebrafish kinked notochord and two worm phenotypes that are not obvious, with the "
                      "idea that the same genes give a different organism-level readout in each species.")
def fig_nonobvious(D, step=None):
    return _nonobvious(D, D.rod_cone)


@figure("nonobvious_hboc", "Hereditary breast ovarian cancer syndrome linked to four model phenotypes: mouse mammary "
                           "adenocarcinoma as expected, then zebrafish with no females, a worm germ-cell response to UV "
                           "and radiation-sensitive fission yeast, which are not obvious.")
def fig_nonobvious_hboc(D, step=None):
    return _nonobvious(D, D.hboc)


@figure("history", "Timeline of phenolog work: McGary 2010 in PNAS, Cha 2012 in PLoS Biology, Woods 2013 in BMC "
                   "Bioinformatics, and this work in 2026, with the two flavours of phenolog underneath.")
def fig_history(D, step=None):
    c = Canvas()
    h = D.history
    m, ch, w = h.mcgary2010, h.cha2012, h.woods2013
    cards = [
        ("2010", f"{m.cite} · {m.journal}", m.title, f"PMID {m.pmid}",
         [f"{m.species} species, {n(m.gene_phenotype_associations)} gene–phenotype associations",
          f"more than {n(m.phenologs_total)} phenologs at 5% FDR",
          "plant gravitropism defects ↔ Waardenburg syndrome: SEC23IP"], False),
        ("2012", f"{ch.cite} · {ch.journal}", ch.title, f"PMID {ch.pmid}",
         ["yeast lovastatin sensitivity ↔ vertebrate angiogenesis defects",
          "led to thiabendazole, an antifungal, as a vascular disrupting agent"], False),
        ("2013", f"{w.cite} · {w.journal}", w.title, f"PMID {w.pmid}",
         ["ranks candidate genes from a phenotype's k nearest phenologs",
          "naive Bayes and additive scores", "leave-one-out cross validation"], False),
        ("2026", "This work · bioRxiv", "Phenologs at knowledge-graph scale", f"doi {D.meta.doi}",
         [f"8 species from one knowledge graph", f"{n(D.totals.disease_phenologs)} disease phenologs",
          "gene-level hypergeometric score, time-split validation"], True),
    ]
    cw, gap, x, ty = 2.88, 0.16, 0.0, 0.42
    c.line(0.2, ty, c.w - 0.2, ty, color=T.RULE, lw=2)
    for year, head, title, ref, facts, ours in cards:
        fg, sub, fill, edge = (T.WHITE, "#C9DDF0", T.NAVY, None) if ours else (T.INK, T.MUTED, T.WHITE, T.RULE)
        c.circle(x + 0.5, ty, 0.1, fill=T.NAVY if ours else T.CAROLINA, edge=T.WHITE, lw=2, zorder=5)
        c.text(x + 0.72, ty - 0.2, year, size=15, weight="bold", color=T.NAVY)
        top = 0.75
        c.box(x, top, cw, 3.5, fill=fill, edge=edge, lw=1.2, r=0.12)
        c.text(x + 0.2, top + 0.3, head, size=c.fit(head, cw - 0.4, 11.5, "bold"), weight="bold",
               color=T.WHITE if ours else T.CAROLINA_DARK)
        tt = c.wrap(title, cw - 0.4, 11, "bold")
        c.text(x + 0.2, top + 0.52, tt, size=11, weight="bold", color=fg, va="top", linespacing=1.2)
        ty2 = top + 0.52 + 0.195 * (tt.count("\n") + 1) + 0.12
        c.text(x + 0.2, ty2, ref, size=9.5, color=sub, va="top")
        fy = max(ty2 + 0.42, top + 1.62)
        for fact in facts:
            ft = c.wrap(fact, cw - 0.58, 11)
            c.circle(x + 0.27, fy + 0.09, 0.035, fill=T.CAROLINA if not ours else T.WHITE, edge=None)
            c.text(x + 0.4, fy, ft, size=11, color=fg, va="top", linespacing=1.2)
            fy += 0.195 * (ft.count("\n") + 1) + 0.13
        x += cw + gap

    by = 4.4
    c.text(0.05, by + 0.36, "Two flavours", size=12.5, weight="bold", color=T.MUTED)
    a = c.pill_left(1.45, by + 0.36, "classical:  phenotype ↔ phenotype", size=12.5, h=0.5, weight="bold")
    b2 = c.pill_left(a.right + 0.3, by + 0.36, "disease-framed:  disease ↔ phenotype", size=12.5, h=0.5, weight="bold",
                     fill=T.TINT_2, edge=T.CAROLINA)
    c.text(b2.right + 0.2, by + 0.36, "a disease is just another gene set", size=12, color=T.MUTED, style="italic")
    return c


@figure("derivation", "How one phenolog is derived: Bardet-Biedl syndrome's 26 genes and the 115 genes of mouse retina "
                      "degeneration are collapsed to PANTHER families, 13 of which are shared out of 7,394, giving a "
                      "hypergeometric P of 1.3e-18.")
def fig_derivation(D, step=None):
    c = Canvas()
    d, b, mouse = D.bbs.derivation, D.bbs, D.sp.mouse
    y1, y2 = 1.55, 3.0
    gx, fx = 3.55, 6.05                      # centres of the gene and family columns
    c.text(gx, 0.42, "Genes", size=14, weight="bold", color=T.MUTED, ha="center")
    c.text(fx, 0.42, "PANTHER protein families", size=14, weight="bold", color=T.MUTED, ha="center")
    c.text(8.72, 0.42, "Hypergeometric test", size=14, weight="bold", color=T.MUTED, ha="center")
    c.text(11.1, 0.42, "Phenolog", size=14, weight="bold", color=T.MUTED, ha="center")

    c.badge(0.3, y1, "Hs", T.NAVY, r=0.2)
    c.text(0.62, y1, "Bardet–Biedl\nsyndrome", size=14, weight="bold")
    c.badge(0.3, y2, mouse.code, mouse.color, r=0.2)
    c.text(0.62, y2, d.phenotype.replace(" ", "\n", 1), size=14, weight="bold")

    rows = [(y1, d.disease_genes, d.disease_gene_examples, d.disease_families, T.NAVY),
            (y2, d.phenotype_genes, d.phenotype_gene_examples, d.phenotype_families, mouse.color)]
    fam_boxes = []
    for y, genes, eg, fams, col in rows:
        g = c.pill(gx, y, f"{n(genes)} genes", sub=f"e.g. {eg}", size=14.5, weight="bold", sub_size=11.5,
                   sub_color=T.MUTED, sub_style="italic", w=2.5, edge=T.RULE, lw=1.3, h=0.86)
        f = c.pill(fx, y, f"{n(fams)} families", sub=f"e.g. {d.family_example}", size=14.5, weight="bold", sub_size=11.5,
                   sub_color="#C9DDF0", w=2.05, fill=T.NAVY, edge=None, color=T.WHITE, h=0.86)
        c.arrow(g.right + 0.08, y, f.x - 0.08, y, color=T.FAINT, lw=1.5)
        fam_boxes.append(f)

    bx = fam_boxes[0].right + 0.12
    c.line(bx, y1, bx + 0.18, y1, color=T.NAVY, lw=1.6)
    c.line(bx, y2, bx + 0.18, y2, color=T.NAVY, lw=1.6)
    c.line(bx + 0.18, y1, bx + 0.18, y2, color=T.NAVY, lw=1.6)
    tb = c.box(7.6, 0.95, 2.25, 2.72, fill=T.WHITE, edge=T.RULE, lw=1.3, r=0.12)
    c.arrow(bx + 0.18, (y1 + y2) / 2, tb.x - 0.05, (y1 + y2) / 2, color=T.FAINT, lw=1.5)

    vy, r, off = 2.05, 0.66, 0.36
    c.circle(tb.cx - off, vy, r, fill=T.CAROLINA, edge=None, alpha=0.30, zorder=3)
    c.circle(tb.cx + off, vy, r, fill=mouse.color, edge=None, alpha=0.22, zorder=3)
    c.circle(tb.cx - off, vy, r, fill=None, edge=T.CAROLINA_DARK, lw=1.6, zorder=4)
    c.circle(tb.cx + off, vy, r, fill=None, edge=mouse.color, lw=1.6, zorder=4)
    c.text(tb.cx - off - 0.36, vy, n(d.disease_families - d.shared), size=13, color=T.MUTED, ha="center")
    c.text(tb.cx + off + 0.36, vy, n(d.phenotype_families - d.shared), size=13, color=T.MUTED, ha="center")
    c.text(tb.cx, vy, n(d.shared), size=18, weight="bold", ha="center")
    c.text(tb.cx, vy - r - 0.19, "shared", size=10.5, color=T.MUTED, ha="center")
    c.text(tb.cx, 3.2, f"of {n(d.N)} families\nin both species", size=12, color=T.MUTED, ha="center")
    c.text(tb.cx, 4.02, f"P = {sci(d.p, bold=True)}", size=17, weight="bold", ha="center")

    px = 11.1
    c.arrow(tb.right + 0.08, (y1 + y2) / 2, px - 0.62, (y1 + y2) / 2, color=T.FAINT, lw=1.5)
    top = c.pill(px, 1.52, "Bardet–Biedl\nsyndrome", size=12, w=1.62, h=0.72, lw=1.8, r=0.3)
    bot = c.pill(px, 3.02, d.phenotype.replace(" ", "\n", 1), size=12, w=1.62, h=0.72, lw=1.8, edge=mouse.color, r=0.3)
    c.line(px, top.bottom, px, 2.1, color=T.NAVY, lw=1.8)
    c.line(px, 2.44, px, bot.y, color=T.NAVY, lw=1.8)
    c.text(px + 0.05, 2.27, "phenolog", size=13.5, weight="bold", ha="center")

    c.text(0.1, 4.75, f"Kept at {D.totals.fdr_pct}% FDR: the human–mouse cutoff from {n(D.totals.randomized_trials)} randomized "
                      f"networks is P ≤ {sci(mouse.p_cutoff)}.", size=12.5, color=T.MUTED)
    return c


def _phenolog_tree(c, D, block, flag, ry, top=0.28, row_h=0.362, quiet=True):
    """One human term on the left, a branch per species, a leaf per phenolog with its P value.
    block: bbs or rod_cone from data.yaml; flag: the row key that marks a leaf as the expected kind (navy outline)
    rather than the non-obvious kind (accent). Returns the y just below the last thing drawn."""
    rx = 1.12
    c.badge(rx, ry - 0.72, "Hs", T.NAVY, r=0.19)
    root = c.pill(rx, ry, block.label.replace(" ", "\n", 1), size=14, weight="bold", w=2.0, h=0.82, lw=2.2, r=0.26)
    c.text(rx, ry + 0.66, f"{n(block.known_genes)} genes", size=12.5, color=T.MUTED, ha="center")

    y = top
    groups = []
    for sp in D.species:
        if sp.key not in block.tree:
            continue
        rows = block.tree[sp.key]
        ys = [y + i * row_h for i in range(len(rows))]
        groups.append((sp, rows, ys))
        y = ys[-1] + row_h + 0.15
    for sp, rows, ys in groups:
        gy = sum(ys) / len(ys)
        c.scurve(root.right, ry, 2.95, gy, color=T.NAVY, lw=1.8, bend=0.6)
        c.badge(3.2, gy, sp.code, sp.color, r=0.19)
        c.text(3.5, gy - 0.12, sp.common, size=13.5, weight="bold")
        c.text(3.5, gy + 0.14, f"{n(block.by_species[sp.key])} phenologs", size=11, color=T.MUTED)
        for row, yy in zip(rows, ys):
            c.scurve(4.72, gy, 5.2, yy, color=sp.color, lw=1.3, bend=0.55)
            c.text(6.42, yy, sci(row.p), size=11.5, color=T.MUTED, ha="right")
            edge, fill = (T.NAVY, T.WHITE) if row[flag] else (T.ACCENT, T.ACCENT_TINT)
            c.line(6.48, yy, 6.6, yy, color=edge, lw=1.3)
            c.pill_left(6.6, yy, row.label, size=12, h=0.31, edge=edge, fill=fill, lw=1.5, pad=0.18)
    if not quiet:
        return y - 0.15 - row_h + 0.155
    names = ", ".join(s.common for s in D.species if s.key not in block.by_species)
    qy = y + 0.02
    c.scurve(root.right, ry, 2.95, qy, color=T.FAINT, lw=1.3, ls=(0, (3, 3)), bend=0.6)
    c.text(3.05, qy, names, size=12, weight="bold", color=T.MUTED)
    c.text(3.05 + c.measure(names, 12, "bold")[0] + 0.2, qy, f"no phenolog at {D.totals.fdr_pct}% FDR", size=11.5,
           color=T.FAINT)
    return qy + 0.12


@figure("bbs_tree", "Bardet-Biedl syndrome connected to its phenologs in mouse, zebrafish and worm with P values. "
                    "Mouse phenotypes read like the clinical features; the zebrafish and worm ones do not.")
def fig_bbs_tree(D, step=None):
    c = Canvas()
    _phenolog_tree(c, D, D.bbs, "clinic", ry=2.55)
    lx = 10.05
    c.text(lx, 0.34, "model phenotype…", size=11.5, color=T.MUTED, weight="bold")
    c.pill_left(lx, 0.74, "reads like the clinic", size=11.5, h=0.31, edge=T.NAVY, lw=1.5, pad=0.16)
    c.pill_left(lx, 1.14, "is non-obvious", size=11.5, h=0.31, edge=T.ACCENT, fill=T.ACCENT_TINT, lw=1.5, pad=0.16)
    c.text(lx, 1.62, "P: hypergeometric,\nfamily overlap", size=10.5, color=T.FAINT)
    return c


@figure("rodcone_tree", "Rod-cone dystrophy, a human phenotype, connected to ten of its phenologs with P values. The three "
                        "mouse ones and two of the zebrafish ones are retina phenotypes. Zebrafish kinked notochord and "
                        "cystic pronephros, and the three worm phenotypes (osmotic avoidance defective, cilia morphology "
                        "variant, ivermectin resistant) are not, and are explained only by the genes they share.")
def fig_rodcone_tree(D, step=None):
    c = Canvas()
    rc = D.rod_cone
    bottom = _phenolog_tree(c, D, rc, "obvious", ry=1.83, top=0.18, row_h=0.325, quiet=False)
    c.text(1.12, 1.83 + 0.92, f"human phenotype\n{rc.id}", size=10, color=T.FAINT, ha="center", linespacing=1.25)

    # the two kinds of counterpart
    by, bh, lw_, gap = bottom + 0.2, c.h - bottom - 0.24, 5.55, 0.2
    rx, rw = lw_ + gap, c.w - lw_ - gap
    left, right = c.wrap(rc.similar_text, lw_ - 0.56, 11.5), c.wrap(rc.nonobvious_text, rw - 0.56, 11.5)
    assert max(left.count("\n"), right.count("\n")) <= 2, "the two cards under the tree hold three lines each"
    c.box(0, by, lw_, bh, fill=T.WHITE, edge=T.NAVY, lw=1.5, r=0.12)
    c.text(0.28, by + 0.28, "Semantically similar", size=13.5, weight="bold")
    c.text(0.28, by + 0.52, left, size=11.5, color=T.MUTED, va="top", linespacing=1.28)
    c.box(rx, by, rw, bh, fill=T.ACCENT_TINT, edge=T.ACCENT, lw=1.5, r=0.12)
    c.text(rx + 0.28, by + 0.28, "Non-obvious", size=13.5, weight="bold", color=T.ACCENT_DARK)
    c.text(rx + 0.28, by + 0.52, right, size=11.5, va="top", linespacing=1.28)
    return c
