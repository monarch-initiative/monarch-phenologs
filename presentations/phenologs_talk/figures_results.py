"""Figures for the second half of the talk: scale, scoring, evaluation, resource, backup.

Result plots come from the preprint PDF (paper.py); everything else is drawn from data.yaml.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from matplotlib.patches import PathPatch, Polygon
from matplotlib.path import Path as MPath

from . import paper
from . import theme as T
from .draw import Canvas, sci
from .figures import figure, n

ROOT = Path(__file__).resolve().parents[1]
CROP_DIR = ROOT / "out" / "crops"
LIGHT = "#C9DDF0"          # text on navy


def _tile(c, x, y, w, h, big, cap, sub=None, fill=T.TINT, fg=T.INK, sub_fg=T.MUTED, big_size=30, edge=None):
    """Stat tile: a big number, a caption, an optional smaller line."""
    c.box(x, y, w, h, fill=fill, edge=edge, lw=1.4, r=0.12)
    c.text(x + 0.25, y + 0.2 + big_size / 144, big, size=big_size, weight="bold", color=fg)
    ty = y + 0.2 + big_size / 72 + 0.18
    cap_w = c.wrap(cap, w - 0.5, 12.5, "bold")
    c.text(x + 0.25, ty, cap_w, size=12.5, weight="bold", color=fg, va="top", linespacing=1.2)
    ty += 0.21 * (cap_w.count("\n") + 1) + 0.08
    if sub:
        c.text(x + 0.25, ty, c.wrap(sub, w - 0.5, 11), size=11, color=sub_fg, va="top", linespacing=1.25)


def _fit_pill(c, x, y, label, max_w, edge, fill, size=12, floor=10.5):
    """Pill with its left edge at x. Text shrinks to `floor`, then wraps to two lines."""
    if c.measure(label, floor)[0] <= max_w:
        return c.pill_left(x, y, label, size=c.fit(label, max_w, size, min_size=floor), h=0.42, lw=1.5, pad=0.18,
                           edge=edge, fill=fill)
    wrapped = c.wrap(label, max_w, floor)
    w = max(c.measure(line, floor)[0] for line in wrapped.split("\n")) + 0.36
    return c.pill(x + w / 2, y, wrapped, w=w, size=floor, h=0.52, lw=1.5, edge=edge, fill=fill, r=0.2)


def _note(c, x, y, w, head, body, dot=T.CAROLINA, head_fg=T.INK, body_fg=T.MUTED, size=13):
    """Heading with a dot and wrapped body text. Returns the y below it."""
    c.circle(x + 0.07, y + 0.02, 0.06, fill=dot, edge=None)
    c.text(x + 0.26, y, head, size=size, weight="bold", color=head_fg)
    bt = c.wrap(body, w - 0.3, size - 1.5)
    c.text(x + 0.26, y + 0.2, bt, size=size - 1.5, color=body_fg, va="top", linespacing=1.25)
    return y + 0.2 + 0.2 * (bt.count("\n") + 1) + 0.3


# =============================================================================================
# 4. Phenologs at knowledge-graph scale
# =============================================================================================
@figure("then_now", "Paired bars comparing the data behind the 2010 phenolog study with this work: species, human "
                    "diseases, human gene-disease annotations, human phenotypes, human gene-phenotype annotations and "
                    "model organism gene-phenotype annotations. The 2010 study had no human phenotype annotations. "
                    "Beside it, how the annotations and gene orthology were assembled then and now.")
def fig_then_now(D, step=None):
    c = Canvas()
    m, t, H = D.history.mcgary2010, D.totals, D.human
    then_c, now_c = "#B4C0CC", T.NAVY
    model_then = sum(v for k, v in m.associations.items() if k != "human")
    absent = "none: human data were OMIM diseases only"
    # label, 2010 value, 2026 value, 2010 text, 2026 text, fold-change prefix (None: no fold shown)
    rows = [("Species", m.species, 1 + t.model_species, f"{m.species}: {m.species_names}",
             f"{1 + t.model_species}: human + {t.model_species} model organisms", None),
            ("Human diseases", m.human_diseases_approx, H.diseases_with_gene, m.human_diseases_label,
             n(H.diseases_with_gene), "≈ "),
            ("Human gene–disease\nannotations", m.associations.human, H.gene_disease_pairs, n(m.associations.human),
             n(H.gene_disease_pairs), ""),
            ("Human phenotypes", None, H.phenotype_terms, absent, n(H.phenotype_terms), None),
            ("Human gene–phenotype\nannotations", None, H.phenotype_annotations, absent, n(H.phenotype_annotations), None),
            ("Model organism gene–\nphenotype annotations", model_then, t.model_annotations, n(model_then),
             n(t.model_annotations), "")]
    bx, bmax, fx = 2.62, 2.35, 7.5
    c.rect(bx, 0.1, 0.22, 0.16, fill=then_c)
    c.text(bx + 0.32, 0.18, m.cite, size=11.5, color=T.MUTED)
    c.rect(bx + 2.0, 0.1, 0.22, 0.16, fill=now_c)
    c.text(bx + 2.32, 0.18, "this work, 2026", size=11.5, color=T.MUTED)
    y, rh = 0.46, 0.765
    for label, old, new, old_s, new_s, fold in rows:
        c.text(0.05, y + 0.34, label, size=12.5, weight="bold", linespacing=1.15)
        if old is None:
            c.text(bx, y + 0.17, old_s, size=11, color=T.FAINT, style="italic")
        else:
            w_old = max(0.06, bmax * old / new)
            c.box(bx, y + 0.06, w_old, 0.22, fill=then_c, edge=None, r=0.05)
            c.text(bx + w_old + 0.12, y + 0.17, old_s, size=11, color=T.MUTED)
        c.box(bx, y + 0.37, bmax, 0.22, fill=now_c, edge=None, r=0.05)
        c.text(bx + bmax + 0.12, y + 0.48, new_s, size=c.fit(new_s, fx - bx - bmax - 0.14, 11.5, "bold"), weight="bold")
        if old is None:
            c.text(fx, y + 0.17, "new", size=15, weight="bold", color=T.CAROLINA_DARK, ha="right")
        elif fold is not None:
            f = new / old
            c.text(fx, y + 0.17, fold + (f"{f:.0f}×" if f >= 10 else f"{f:.1f}×"), size=17, weight="bold",
                   color=T.CAROLINA_DARK, ha="right")
        if label != rows[-1][0]:
            c.line(0.05, y + rh - 0.06, fx, y + rh - 0.06, color=T.RULE, lw=0.8)
        y += rh

    px, pw = 7.95, 4.05
    notes = [("2010 and 2013", then_c, "Gene–phenotype annotations were assembled by hand from each source database, and "
                                       "gene homology was computed for each study with the InParanoid algorithm."),
             ("2026", now_c, "The Monarch KG and the Biolink model harmonize the annotation sets, and PANTHER protein "
                             "families resolve gene equivalence across species.")]
    c.panel(px, 0.1, pw, 3.62)
    c.text(px + 0.3, 0.5, "Methodological differences", size=15, weight="bold")
    yy = 1.02
    for head, dot, body in notes:
        yy = _note(c, px + 0.3, yy, pw - 0.6, head, body, dot=dot, size=14) + 0.08
    return c


@figure("pipeline_fdr", "Five pipeline steps from the Monarch KG subgraph to 61,243 phenologs, beside the published "
                        "plot of real against randomized p-value distributions for human diseases and mouse phenotypes.")
def fig_pipeline_fdr(D, step=None):
    c = Canvas()
    t, H, mouse = D.totals, D.human, D.sp.mouse
    steps = [("Extract the subgraph", "8 species: genes, phenotypes, diseases, orthology"),
             ("Collapse genes to families", f"{n(H.genes_with_ortholog)} human genes become {n(H.families)} PANTHER families"),
             ("Test every pair", f"{t.pairwise_tests_disease / 1e6:.0f} million disease × phenotype overlaps, hypergeometric"),
             ("Shuffle the networks", f"{n(t.randomized_trials)} randomized networks per species pair give the null"),
             ("Keep the real signal", f"{t.fdr_pct}% FDR cutoff per species: {n(t.disease_phenologs)} phenologs")]
    h, gap, y, w = 0.84, 0.13, 0.12, 5.75
    for i, (head, body) in enumerate(steps):
        last = i == len(steps) - 1
        c.box(0.0, y, w, h, fill=T.NAVY if last else T.TINT, edge=None, r=0.12)
        c.circle(0.42, y + h / 2, 0.2, fill=T.WHITE if last else T.NAVY, edge=None)
        c.text(0.42, y + h / 2 + 0.01, str(i + 1), size=13, weight="bold", color=T.NAVY if last else T.WHITE, ha="center")
        c.text(0.8, y + 0.28, head, size=13.5, weight="bold", color=T.WHITE if last else T.INK)
        c.text(0.8, y + 0.58, body, size=c.fit(body, w - 1.0, 11.5), color=LIGHT if last else T.MUTED)
        y += h + gap

    img = paper.crop("fig2a", CROP_DIR)
    b = c.image(img, 6.35, 0.05, h=4.25)
    cap = c.wrap(f"Human diseases × mouse phenotypes. P ≤ {sci(mouse.p_cutoff)} keeps {n(mouse.disease_phenologs)} pairs "
                 f"at a realised FDR of {mouse.realised_fdr}%.", 5.5, 11.5)
    c.text(6.4, b.bottom + 0.14, cap, size=11.5, color=T.MUTED, va="top", linespacing=1.25)
    return c


@figure("fdr_published", "Published histogram of hypergeometric P values for human against mouse phenotype pairs: "
                         "real data in blue and permuted networks in grey, with a dashed line at the P-value cutoff "
                         "where permuted pairs are 5% of real ones (preprint Fig. S3). Beside it, how the cutoff is "
                         "computed.")
def fig_fdr_published(D, step=None):
    c = Canvas()
    t, mouse = D.totals, D.sp.mouse
    b = c.image(paper.crop("figS3", CROP_DIR), 0.0, 0.05, w=7.6)
    px = b.right + 0.4
    pw = c.w - px
    c.text(px, 0.3, "How the cutoff is computed", size=15, weight="bold")
    yy = 0.85
    for head, body in [(f"{n(t.randomized_trials)} permuted networks",
                        "per species pair: the same tests are run again on permuted gene–phenotype links"),
                       ("Real against permuted",
                        "at each P-value cutoff, count the real pairs and the permuted pairs (averaged per trial) that pass"),
                       (f"{t.fdr_pct}% FDR cutoff",
                        f"the largest P value where permuted pairs are at most {t.fdr_pct}% of real pairs"),
                       (f"Human vs. mouse: P ≤ {sci(mouse.p_cutoff)}",
                        f"{n(mouse.disease_phenologs)} disease phenologs kept, realised FDR {mouse.realised_fdr}%")]:
        yy = _note(c, px, yy, pw, head, body)
    return c


@figure("counts", "Table of significant phenologs per model species for the disease and phenotype networks, with "
                  "log-scale bars and totals of 61,243 and 2,187,518, above the 3,966 human disease phenologs of "
                  "McGary et al. 2010. Beside it, four counts: human disease terms and human phenotype terms with at "
                  "least one phenolog, and the model organism phenotypes linked to each.")
def fig_counts(D, step=None):
    c = Canvas()
    t, m = D.totals, D.history.mcgary2010
    sp_sorted = sorted(D.species, key=lambda s: -s.disease_phenologs)
    cx1, cx2, bw, tx1 = 2.2, 4.75, 1.5, 7.2
    c.text(cx1, 0.22, "Disease phenologs", size=13.5, weight="bold")
    c.text(cx1, 0.48, "human disease ↔ model phenotype", size=10.5, color=T.MUTED)
    c.text(cx2, 0.22, "Phenotype phenologs", size=13.5, weight="bold")
    c.text(cx2, 0.48, "human phenotype ↔ model phenotype", size=10.5, color=T.MUTED)

    def bar(x, y, v, lo, hi, color):
        frac = (math.log10(v) - lo) / (hi - lo)
        c.box(x, y - 0.13, max(0.08, bw * frac), 0.26, fill=color, edge=None, r=0.05)
        c.text(x + max(0.08, bw * frac) + 0.1, y, n(v), size=12.5, weight="bold")

    y, pitch = 0.95, 0.44
    for s in sp_sorted:
        c.badge(0.24, y, s.code, s.color, r=0.17)
        c.text(0.54, y, s.common, size=13, weight="bold")
        bar(cx1, y, s.disease_phenologs, 1.0, 5.0, s.color)
        bar(cx2, y, s.phenotype_phenologs, 2.0, 6.6, s.color)
        y += pitch
    c.line(0.05, y - 0.16, tx1, y - 0.16, color=T.NAVY, lw=1.4)
    c.text(0.1, y + 0.14, "All seven species", size=13.5, weight="bold")
    c.text(cx1, y + 0.14, n(t.disease_phenologs), size=19, weight="bold")
    c.text(cx2, y + 0.14, n(t.phenotype_phenologs), size=19, weight="bold")

    # the earlier study, for scale: its phenologs between a model organism and human disease (it had no human phenotypes)
    then = sum(mo.phenologs_with_human for mo in m.models)
    assert m.phenologs_human_mouse == next(mo.phenologs_with_human for mo in m.models if mo.name == "mouse")
    y2 = y + 0.62
    c.line(0.05, y2 - 0.26, tx1, y2 - 0.26, color=T.RULE, lw=0.8)
    c.text(0.1, y2, m.cite, size=13.5, weight="bold", color=T.MUTED)
    c.text(cx1, y2, n(then), size=19, weight="bold", color=T.MUTED)
    c.text(cx1 + c.measure(n(then), 19, "bold")[0] + 0.14, y2 + 0.02, f"{len(m.models)} model species", size=10.5,
           color=T.MUTED)
    c.text(0.1, c.h - 0.16, f"{t.fdr_pct}% FDR, estimated per species from {n(t.randomized_trials)} randomized networks. "
                            "Bars are on a log scale.", size=10.5, color=T.FAINT)

    # how many terms on each side take part: the disease network (navy), then the same two counts for the phenotype network
    px = 7.5
    pw, th, gap = c.w - px, 1.1, 0.13
    tiles = [(t.diseases_with_phenologs, "human disease terms with\n≥ 1 phenolog association", True),
             (t.model_phenotypes_with_diseases, "model organism phenotypes\nlinked to ≥ 1 human disease", True),
             (t.hp_terms_with_phenologs, "human phenotype terms with\n≥ 1 phenolog association", False),
             (t.model_phenotypes_with_hp, "model organism phenotypes\nlinked to ≥ 1 human phenotype", False)]
    num_w = max(c.measure(n(v), 28, "bold")[0] for v, _, _ in tiles)
    cap_x = px + 0.28 + num_w + 0.28
    size = min(c.fit(line, px + pw - 0.2 - cap_x, 12.5, "bold") for _, cap, _ in tiles for line in cap.split("\n"))
    for i, (v, cap, disease) in enumerate(tiles):
        ty = 0.1 + i * (th + gap)
        c.box(px, ty, pw, th, fill=T.NAVY if disease else T.TINT, edge=None, r=0.12)
        fg = T.WHITE if disease else T.INK
        c.text(px + 0.28, ty + th / 2 + 0.01, n(v), size=28, weight="bold", color=fg)
        c.text(cap_x, ty + th / 2, cap, size=size, weight="bold", color=fg, linespacing=1.25)
    return c


@figure("top_species", "One top-ranked disease phenolog for each of the seven model species, from mouse to slime mold, "
                       "with shared family counts and P values.")
def fig_top_species(D, step=None):
    c = Canvas()
    c.text(1.95, 0.2, "Human disease", size=12, weight="bold", color=T.MUTED)
    c.text(5.95, 0.2, "Model-organism phenotype", size=12, weight="bold", color=T.MUTED)
    c.text(10.35, 0.2, "shared", size=12, weight="bold", color=T.MUTED, ha="center")
    c.text(11.95, 0.2, "P", size=12, weight="bold", color=T.MUTED, ha="right")
    y = 0.74
    for row in D.top_per_species:
        sp = D.sp[row.species]
        c.badge(0.24, y, sp.code, sp.color, r=0.19)
        c.text(0.54, y, sp.common, size=12.5, weight="bold")
        d = _fit_pill(c, 1.95, y, row.disease, 3.3, T.NAVY, T.WHITE)
        edge, fill = (sp.color, T.WHITE) if row.obvious else (T.ACCENT, T.ACCENT_TINT)
        p = _fit_pill(c, 5.95, y, row.phenotype, 3.3, edge, fill)
        c.line(d.right + 0.06, y, p.x - 0.06, y, color=T.FAINT, lw=1.4, ls=(0, (2, 2)))
        c.text(10.35, y, f"{row.c} of {row.m}", size=12, ha="center")
        c.text(11.95, y, sci(row.p, bold=True), size=12, ha="right", weight="bold")
        y += 0.57
    c.pill_left(1.95, y + 0.05, "phenotype reads like the disease", size=10.5, h=0.3, edge=T.MUTED, lw=1.2, pad=0.14,
                color=T.MUTED)
    c.pill_left(4.6, y + 0.05, "non-obvious", size=10.5, h=0.3, edge=T.ACCENT, fill=T.ACCENT_TINT, lw=1.2, pad=0.14)
    c.text(5.85, y + 0.05, "shared: PANTHER families in common, of the disease's families", size=10.5, color=T.FAINT)
    return c


# =============================================================================================
# 5. From phenologs to candidate genes
# =============================================================================================
@figure("pooling", "Three panels and a strip. One: Bardet-Biedl syndrome's 109 significant phenologs, ranked by P. Two: each "
                   "gene counted by how many of the 109 phenologs carry it, BBS1 in 50 of 109. Three: the output, a ranked "
                   "list of candidate genes whose score column is left as question marks. The strip below names three "
                   "ways to compute that score: naive Bayes, additive and hypergeometric.")
def fig_pooling(D, step=None):
    c = Canvas()
    b = D.bbs
    y0, ph = 0.05, 3.75                       # the three panels
    sy, sh = 3.94, 1.08                       # the strip under them
    x1, w1 = 0.0, 4.1
    x2, w2 = 4.38, 3.3
    x3, w3 = 7.96, 4.04
    head_y, row_y = y0 + 0.61, y0 + 0.97

    def vdots(x, y):
        for dy in (-0.06, 0.0, 0.06):
            c.circle(x, y + dy, 0.012, fill=T.FAINT, edge=None)

    # panel 1: the phenologs
    c.panel(x1, y0, w1, ph, tag="1", title=f"{n(b.phenologs)} significant phenologs")
    c.text(0.62, head_y, "model phenotype", size=9.5, color=T.FAINT)
    c.text(2.85, head_y, "shared", size=9.5, color=T.FAINT, ha="center")
    c.text(3.92, head_y, "P", size=9.5, color=T.FAINT, ha="right")
    y = row_y
    for row in b.top_phenologs:
        sp = D.sp[row.species]
        c.badge(0.33, y, sp.code, sp.color, r=0.15)
        c.text(0.62, y, row.label, size=c.fit(row.label, 1.9, 11.5, "bold"), weight="bold")
        c.text(2.85, y, f"{row.c}/{row.n}", size=11, color=T.MUTED, ha="center")
        c.text(3.92, y, sci(row.p), size=10.5, color=T.MUTED, ha="right")
        y += 0.44
    c.text(0.62, y - 0.04, f"… {n(b.phenologs - len(b.top_phenologs))} more", size=11, color=T.FAINT, style="italic")
    chips_y, x = y0 + ph - 0.3, 0.3
    for key, cnt in b.by_species.items():
        sp = D.sp[key]
        c.badge(x + 0.12, chips_y, sp.code, sp.color, r=0.14)
        c.text(x + 0.33, chips_y, f"{cnt}", size=11.5, weight="bold")
        x += 0.95

    # panel 2: recurrence, as "50 of 109"
    c.panel(x2, y0, w2, ph, tag="2", title="Count each gene's recurrence")
    c.text(x2 + 0.2, head_y, "phenologs carrying an ortholog of the gene", size=9.5, color=T.FAINT)
    ranked = sorted(b.candidates, key=lambda g: g.rank)
    new = [g for g in ranked if not g.known][:2]
    rec = lambda g: sum(g.support.values())
    shown = sorted([g for g in ranked if g.known], key=lambda g: -rec(g))[:6] + new
    top, bx, bw = max(rec(g) for g in shown), x2 + 1.04, 1.25
    y = row_y
    for g in shown:
        c.text(x2 + 0.2, y, g.gene, size=11.5, weight="bold")
        c.box(bx, y - 0.11, bw * rec(g) / top, 0.22, fill=T.NAVY, edge=None, r=0.04)
        c.text(bx + bw * rec(g) / top + 0.09, y, f"{rec(g)} of {n(b.phenologs)}", size=10.5, color=T.MUTED)
        y += 0.36

    # panel 3: the output we are after. The ranks are the ones the hypergeometric score gives; the score itself is
    # left open, because the next slides are about how to compute it.
    c.panel(x3, y0, w3, ph, tag="3", title="Rank the candidate genes")
    cr, cg, cs, cn, ct = x3 + 0.2, x3 + 0.62, x3 + 1.74, x3 + 2.42, x3 + 2.84
    rows, y = [], row_y
    for g in ranked[:5]:
        rows.append((y, g))
        y += 0.32
    for g in new:
        vdots(cr + 0.05, y - 0.095)
        y += 0.13
        rows.append((y, g))
        y += 0.32
    last = rows[-1][0]
    c.box(cs - 0.31, head_y - 0.17, 0.62, last - head_y + 0.34, fill=T.ACCENT_TINT, edge=T.ACCENT, lw=1.2, r=0.1,
          zorder=1)
    c.text(cr, head_y, "rank", size=9.5, color=T.FAINT)
    c.text(cg, head_y, "gene", size=9.5, color=T.FAINT)
    c.text(cs, head_y, "score", size=10, weight="bold", color=T.ACCENT_DARK, ha="center")
    c.text(cn, head_y, "recurs in", size=9.5, color=T.FAINT, ha="center")
    c.text(ct, head_y, "status", size=9.5, color=T.FAINT)
    for y, g in rows:
        c.text(cr, y, f"{g.rank}", size=11, color=T.MUTED)
        c.text(cg, y, g.gene, size=11.5, weight="bold")
        c.text(cs, y, "?", size=13, weight="bold", color=T.ACCENT_DARK, ha="center")
        c.text(cn, y, f"{rec(g)}", size=11, color=T.MUTED, ha="center")
        if g.known:
            c.text(ct, y, "known", size=10.5, color=T.MUTED)
        else:
            c.text(ct, y, "new candidate", size=10.5, weight="bold")
    vdots(cr + 0.05, last + 0.225)
    c.text(cr, last + 0.43, f"{n(b.genes_ranked)} genes ranked", size=11, color=T.FAINT, style="italic")

    ay = y0 + ph / 2
    c.arrow(x1 + w1 + 0.04, ay, x2 - 0.04, ay, color=T.MUTED, lw=1.8, head=0.12)
    c.arrow(x2 + w2 + 0.04, ay, x3 - 0.04, ay, color=T.MUTED, lw=1.8, head=0.12)

    # the strip: the question the "?" column asks, and the three ways to answer it
    c.box(0.0, sy, T.FIG_W, sh, fill=T.NAVY, edge=None, r=0.14, zorder=0)
    c.text(0.22, sy + 0.27, "Gene scoring / ranking problem", size=12.5, weight="bold", color=LIGHT)
    big = n(b.genes_ranked)
    c.text(0.22, sy + 0.73, big, size=26, weight="bold", color=T.WHITE)
    tx = 0.22 + c.measure(big, 26, "bold")[0] + 0.18
    c.text(tx, sy + 0.62, "genes pooled for this disease", size=12.5, weight="bold", color=T.WHITE)
    c.text(tx, sy + 0.85, "Which are most likely?", size=12.5, color=LIGHT)
    pw, gap, pill_h, cy = 2.25, 0.12, 0.64, sy + sh / 2
    px = T.FIG_W - 0.22 - 3 * pw - 2 * gap
    c.arrow(px - 0.62, cy, px - 0.14, cy, color=LIGHT, lw=1.6, head=0.11)
    for i, (name, who, ours) in enumerate([("naive Bayes", "McGary 2010, Woods 2013", False),
                                           ("additive", "Woods 2013", False), ("hypergeometric", "this work", True)]):
        x = px + i * (pw + gap)
        c.box(x, cy - pill_h / 2, pw, pill_h, fill=T.ACCENT if ours else "#24406B", edge=None, r=0.1)
        c.text(x + 0.17, cy - 0.11, name, size=12.5, weight="bold", color=T.WHITE)
        c.text(x + 0.17, cy + 0.14, who, size=9.5, color=T.WHITE if ours else LIGHT)
    return c


@figure("naive_bayes", "A candidate gene carried by three phenologs. Each phenolog contributes a weight; naive Bayes "
                       "combines them as one minus the product of the complements, the additive score sums them.")
def fig_naive_bayes(D, step=None):
    c = Canvas()
    b = D.bbs
    wit = [b.top_phenologs[0], b.top_phenologs[3], b.top_phenologs[4]]
    assert b.known_genes == b.families and all(w.m == b.families for w in wit), "the cards say genes for families"
    short = b.label.removesuffix(" syndrome")
    weights = [(w.c / w.m) * (1 - w.p) for w in wit]
    nb = 1 - float(np.prod([1 - x for x in weights]))
    add = sum(weights)

    gx, gy = 1.05, 2.5
    g = c.pill(gx, gy, "candidate\ngene X", size=14, weight="bold", w=1.7, h=0.86, edge=T.ACCENT, fill=T.ACCENT_TINT,
               lw=2.2, r=0.26)
    c.text(gx, gy + 0.72, c.wrap(f"An example arbitrary gene candidate for {b.label}", 1.95, 10.5), size=10.5,
           color=T.MUTED, ha="center", va="top", linespacing=1.2)

    cx, cw, ch = 2.55, 4.85, 1.3
    ys = [0.25, 1.85, 3.45]
    nbx, bus, rw, rh = 8.05, cx + cw + 0.3, 3.95, 2.05
    r1, r2 = ys[0], ys[-1] + ch - rh                   # the two result boxes: top with the first card, bottom with the last
    for w, wt, y in zip(wit, weights, ys):
        sp = D.sp[w.species]
        c.scurve(g.right, gy, cx, y + ch / 2, color=T.ACCENT, lw=1.6)
        c.box(cx, y, cw, ch, fill=T.WHITE, edge=T.RULE, lw=1.3, r=0.12)
        c.badge(cx + 0.36, y + 0.45, sp.code, sp.color, r=0.19)
        c.text(cx + 0.68, y + 0.45, w.label, size=14.5, weight="bold")
        # "disease genes": this disease's known genes sit one to a protein family, so shared families are shared genes
        val = f"{wt:.2f}"
        line = f"{w.c}/{w.m} {short} disease genes shared in common  ×  (1 − P)"
        c.text(cx + 0.18, y + 0.95, line, size=c.fit(line, cw - 0.4, 12, min_size=10.5), color=T.MUTED)
        c.text(cx + cw - 0.22, y + 0.45, val, size=24, weight="bold", ha="right", color=T.CAROLINA_DARK)   # on the title's row
        c.line(cx + cw, y + ch / 2, bus, y + ch / 2, color=T.FAINT, lw=1.5)
    c.line(bus, ys[0] + ch / 2, bus, ys[-1] + ch / 2, color=T.FAINT, lw=1.5)
    c.arrow(bus, r1 + rh / 2, nbx - 0.05, r1 + rh / 2, color=T.FAINT, lw=1.5)
    c.arrow(bus, r2 + rh / 2, nbx - 0.05, r2 + rh / 2, color=T.FAINT, lw=1.5)

    c.box(nbx, r1, rw, rh, fill=T.NAVY, edge=None, r=0.12)
    c.text(nbx + 0.28, r1 + 0.38, "naive Bayes", size=15, weight="bold", color=T.WHITE)
    c.text(nbx + rw - 0.28, r1 + 0.38, "McGary 2010 · Woods 2013", size=10, color=LIGHT, ha="right")
    expr = "1 − " + "".join(f"(1 − {x:.2f})" for x in weights)
    c.text(nbx + 0.28, r1 + 1.0, expr, size=c.fit(expr, rw - 0.56, 14), color=LIGHT)
    c.text(nbx + 0.28, r1 + 1.58, f"= {nb:.2f}", size=26, weight="bold", color=T.WHITE)

    c.box(nbx, r2, rw, rh, fill=T.WHITE, edge=T.NAVY, lw=1.5, r=0.12)
    c.text(nbx + 0.28, r2 + 0.38, "additive", size=15, weight="bold")
    c.text(nbx + rw - 0.28, r2 + 0.38, "Woods 2013", size=10, color=T.MUTED, ha="right")
    c.text(nbx + 0.28, r2 + 1.0, " + ".join(f"{x:.2f}" for x in weights), size=14, color=T.MUTED)
    total = f"= {add:.2f}"
    c.text(nbx + 0.28, r2 + 1.58, total, size=26, weight="bold")
    c.text(nbx + 0.28 + c.measure(total, 26, "bold")[0] + 0.3, r2 + 1.6, "visualization friendly", size=11.5,
           color=T.MUTED)
    return c


@figure("hypergeom", "An urn of all gene-phenotype annotations. A gene's own annotations are marked; the disease's "
                     "phenologs draw a subset. The score asks whether the gene appears in the draw more often than "
                     "its background frequency predicts, which separates a specific gene from a hub gene.")
def fig_hypergeom(D, step=None):
    c = Canvas()
    ill = D.hg_illustration
    rng = np.random.default_rng(5)
    ux, uy, uw, uh = 0.1, 0.62, 5.5, 2.9
    c.text(ux, 0.28, "N", size=14, weight="bold")
    c.text(ux + 0.25, 0.28, "every gene–phenotype annotation in the cross-species graph", size=11.5, color=T.MUTED)
    c.box(ux, uy, uw, uh, fill=T.TINT, edge=None, r=0.14)
    cols, rows = 24, 11
    dx, dy = (uw - 0.5) / (cols - 1), (uh - 0.5) / (rows - 1)
    box_c, box_r = (14, 22), (2, 8)                  # the draw: columns and rows, inclusive
    cells = [(i, j) for i in range(cols) for j in range(rows)]
    inside = [p for p in cells if box_c[0] <= p[0] <= box_c[1] and box_r[0] <= p[1] <= box_r[1]]
    outside = [p for p in cells if p not in inside]
    marked = [inside[k] for k in rng.choice(len(inside), 7, replace=False)] + \
             [outside[k] for k in rng.choice(len(outside), 5, replace=False)]
    for (i, j) in cells:
        x, y = ux + 0.25 + i * dx, uy + 0.25 + j * dy
        if (i, j) in marked:
            c.circle(x, y, 0.075, fill=T.ACCENT, edge=T.WHITE, lw=0.8, zorder=4)
        else:
            c.circle(x, y, 0.045, fill="#B9C9DA", edge=None, zorder=3)
    bx0, by0 = ux + 0.25 + box_c[0] * dx - 0.14, uy + 0.25 + box_r[0] * dy - 0.14
    bw, bh = (box_c[1] - box_c[0]) * dx + 0.28, (box_r[1] - box_r[0]) * dy + 0.28
    c.box(bx0, by0, bw, bh, fill=None, edge=T.NAVY, lw=1.8, r=0.1, ls=(0, (4, 3)), zorder=5)

    ly = uy + uh + 0.3
    c.circle(ux + 0.1, ly, 0.075, fill=T.ACCENT, edge=None)
    c.text(ux + 0.28, ly, "m", size=13, weight="bold", color=T.ACCENT)
    c.text(ux + 0.52, ly, "this gene's annotations, all species", size=11.5, color=T.MUTED)
    c.line(ux + 0.02, ly + 0.36, ux + 0.2, ly + 0.36, color=T.NAVY, lw=1.8, ls=(0, (3, 2)))
    c.text(ux + 0.28, ly + 0.36, "n", size=13, weight="bold")
    c.text(ux + 0.52, ly + 0.36, "annotations pulled in by the disease's phenologs", size=11.5, color=T.MUTED)
    c.text(ux + 0.28, ly + 0.72, "c", size=13, weight="bold", color=T.ACCENT)
    c.text(ux + 0.52, ly + 0.72, "of the gene's annotations landed in that draw", size=11.5, color=T.MUTED)
    c.text(ux, ly + 1.1, "score = P(X ≥ c)", size=14, weight="bold")
    c.text(ux + 1.85, ly + 1.1, "hypergeometric upper tail, given N, m, n", size=11.5, color=T.MUTED)

    # two genes
    gx, gw = 6.1, 2.85
    c.text(gx, 0.28, "Same recurrence, different verdict", size=13.5, weight="bold")
    c.text(11.95, 0.28, "illustrative numbers", size=10, color=T.FAINT, ha="right", style="italic")
    for k, (g, hot) in enumerate([(ill.specific, True), (ill.hub, False)]):
        x = gx + k * (gw + 0.2)
        expected = ill.n * g.m / ill.N
        c.box(x, 0.62, gw, 2.9, fill=T.ACCENT_TINT if hot else T.WHITE, edge=T.ACCENT if hot else T.RULE, lw=1.6, r=0.12)
        c.text(x + 0.22, 0.95, g.label, size=13.5, weight="bold", color=T.ACCENT if hot else T.INK)
        c.text(x + 0.22, 1.27, f"annotated to {n(g.m)} phenotypes", size=11, color=T.MUTED)
        top = max(g.c, expected)
        for r_, (lab, val, col) in enumerate([("expected", expected, "#B4C0CC"), ("seen", g.c, T.ACCENT if hot else T.NAVY)]):
            yy = 1.85 + r_ * 0.5
            c.text(x + 0.22, yy, lab, size=11, color=T.MUTED)
            c.box(x + 1.0, yy - 0.13, max(0.05, 1.2 * val / top), 0.26, fill=col, edge=None, r=0.04)
            c.text(x + 1.0 + max(0.05, 1.2 * val / top) + 0.08, yy, f"{val:.1f}" if val < 10 else f"{val:.0f}", size=11.5,
                   weight="bold")
        c.text(x + 0.22, 3.08, f"{g.c / expected:.0f}× expected" if g.c / expected >= 10 else f"{g.c / expected:.1f}× expected",
               size=15, weight="bold", color=T.ACCENT if hot else T.MUTED)

    c.box(gx, 3.78, 5.9, 1.22, fill=T.NAVY, edge=None, r=0.12)
    c.text(gx + 0.28, 4.14, "Every gene is judged against the same background,", size=12.5, color=LIGHT)
    c.text(gx + 0.28, 4.42, "so scores compare across diseases:", size=12.5, color=LIGHT)
    c.text(gx + 0.28, 4.74, f"one ranked list of {n(D.totals.ranked_pairs)} gene–disease pairs", size=13.5, weight="bold",
           color=T.WHITE)
    return c


# --- the hypergeometric score as a bag of marbles ------------------------------------------------
# `marbles` works one real gene through one real disease (bbs.score_example in data.yaml): its bag is schematic,
# its numbers are the release's. `marbles_two_genes` is drawn to made-up counts (hg_illustration.marbles),
# so the P values printed on it are exact for the picture.
OTHER_MARBLES = ["#B5C6D8", "#C6D4E2", "#A7BACE"]      # every other gene's marbles
SHARE_WORDS = {2: "half", 3: "a third", 4: "a quarter", 5: "a fifth", 6: "a sixth", 8: "an eighth", 10: "a tenth"}


def _hg_tail(c_, N, m, n_):
    """Exact hypergeometric upper tail P(X >= c_): N marbles, m of one colour, n_ drawn."""
    return sum(math.comb(m, k) * math.comb(N - m, n_ - k) for k in range(c_, min(m, n_) + 1)) / math.comb(N, n_)


def _p_text(p, bold=False):
    return sci(p, bold=bold) if p < 0.01 else f"{p:.2f}"


def _share(part, whole):
    """20 of 120 -> 'a sixth'."""
    return SHARE_WORDS.get(whole // part, f"{100 * part / whole:.0f}%") if whole % part == 0 else f"{100 * part / whole:.0f}%"


def _marble(c, x, y, r, color, zorder=4, hollow=False):
    if hollow:
        c.circle(x, y, r * 0.88, fill=T.WHITE, edge=color, lw=1.1, zorder=zorder)
        return
    c.circle(x, y, r, fill=color, edge=None, zorder=zorder)
    c.circle(x - r * 0.3, y - r * 0.32, r * 0.27, fill=T.WHITE, edge=None, zorder=zorder + 1, alpha=0.5)


def _bag(c, x, y, w, h, fill=T.WHITE, edge=T.NAVY, lw=2.0):
    """A drawstring bag with its top left at (x, y). Returns the outline of its body, for packing marbles into."""
    P = lambda u, v: (x + u * w, y + v * h)
    ny = 0.21                                                    # height of the tie
    ruffle, codes = [P(0.38, ny), P(0.27, 0.06)], [MPath.MOVETO, MPath.LINETO]
    xs = [0.27, 0.385, 0.5, 0.615, 0.73]
    for a, b in zip(xs, xs[1:]):
        ruffle += [P((a + b) / 2, -0.02), P(b, 0.06)]
        codes += [MPath.CURVE3, MPath.CURVE3]
    ruffle += [P(0.62, ny), P(0.38, ny)]
    codes += [MPath.LINETO, MPath.CLOSEPOLY]
    c.ax.add_patch(PathPatch(MPath(ruffle, codes), facecolor=fill, edgecolor=edge, lw=lw, zorder=1, joinstyle="round"))
    body = MPath([P(0.38, ny), P(0.14, 0.27), P(0.0, 0.47), P(0.0, 0.70), P(0.0, 0.95), P(0.12, 1.0), P(0.5, 1.0),
                  P(0.88, 1.0), P(1.0, 0.95), P(1.0, 0.70), P(1.0, 0.47), P(0.86, 0.27), P(0.62, ny), P(0.38, ny)],
                 [MPath.MOVETO] + [MPath.CURVE4] * 12 + [MPath.CLOSEPOLY])
    c.ax.add_patch(PathPatch(body, facecolor=fill, edgecolor=edge, lw=lw, zorder=1, joinstyle="round"))
    for end in ((0.80, ny + 0.10), (0.73, ny + 0.13)):           # the loose ends of the string
        string = MPath([P(0.64, ny), P(end[0], ny - 0.01), P(*end)], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
        c.ax.add_patch(PathPatch(string, facecolor="none", edgecolor=edge, lw=lw * 0.8, zorder=2, capstyle="round"))
    c.box(x + 0.34 * w, y + (ny - 0.026) * h, 0.32 * w, 0.052 * h, fill=edge, edge=None, r=0.05, zorder=3)
    return body


def _pack(path, count, r_max=0.16, top=0.17, margin=0.04, gap=0.016):
    """Hex-pack `count` marbles into `path` from the bottom up, as large as will fit below `top` (a fraction of the
    path's height, so the bag is not full to the neck). Returns (centres, radius)."""
    (x0, y0), (x1, y1) = path.get_extents().get_points()
    cx, y_lo = (x0 + x1) / 2, y0 + top * (y1 - y0)
    ring = np.array([(math.cos(a), math.sin(a)) for a in np.linspace(0, 2 * math.pi, 16, endpoint=False)])
    r = r_max
    while True:
        pitch = 2 * r + gap
        ys = np.arange(y1 - margin - r, y_lo + r, -pitch * 0.866)
        cen = np.array([(cx + (pitch / 2 if i % 2 else 0.0) + k * pitch, y) for i, y in enumerate(ys) for k in range(-30, 31)])
        if len(cen):
            probe = (cen[:, None, :] + (r + margin) * ring[None, :, :]).reshape(-1, 2)
            ok = path.contains_points(probe).reshape(len(cen), len(ring)).all(axis=1)
            cen = cen[ok]
            while len(cen) and (np.isclose(cen[:, 1], cen[:, 1].max())).sum() < 4:      # no lone marble under the pile
                cen = cen[~np.isclose(cen[:, 1], cen[:, 1].max())]
        if len(cen) >= count or r < 0.03:
            break
        r -= 0.0005
    pts = sorted(map(tuple, cen), key=lambda p: (-round(p[1], 4), abs(p[0] - cx)))
    return pts[:count], r


def _marble_colours(pts, groups, rng, spread=0.0):
    """One colour per marble. groups is [(colour, count), ...]; the first group's marbles are kept `spread` apart so
    a rare colour does not clump. Everything else gets a shade of OTHER_MARBLES."""
    cols = [OTHER_MARBLES[k] for k in rng.integers(0, len(OTHER_MARBLES), len(pts))]
    free = list(rng.permutation(len(pts)))
    for gi, (colour, count) in enumerate(groups):
        chosen = []
        for i in list(free):
            if len(chosen) == count:
                break
            if gi == 0 and any(math.dist(pts[i], pts[j]) < spread for j in chosen):
                continue
            chosen.append(i)
            free.remove(i)
        assert len(chosen) == count, "could not place the marbles; lower `spread`"
        for i in chosen:
            cols[i] = colour
    return cols


def _tray(c, cx, cy, rows, r, cols, gap=0.03, fill=T.WHITE, edge=T.NAVY):
    """Marbles lying in a tray, seen from above. rows: marbles per row, e.g. [7, 6, 7]."""
    pitch = 2 * r + gap
    w = max(rows) * pitch + 0.26
    h = (len(rows) - 1) * pitch * 0.866 + 2 * r + 0.3
    b = c.box(cx - w / 2, cy - h / 2, w, h, fill=fill, edge=edge, lw=1.8, r=0.18, zorder=1)
    k = 0
    for i, cnt in enumerate(rows):
        y = cy + (i - (len(rows) - 1) / 2) * pitch * 0.866
        for j in range(cnt):
            _marble(c, cx + (j - (cnt - 1) / 2) * pitch, y, r, cols[k])
            k += 1
    return b


def _tray_rows(count):
    """A compact three-row layout for the handful: 20 -> [7, 6, 7]."""
    mid = (count - 2) // 3
    side = (count - mid) // 2
    assert 2 * side + mid == count and side - mid == 1, "choose a handful of 3k + 2 marbles (17, 20, 23, ...)"
    return [side, mid, side]


def _marble_row(c, x, y, text, colour, r=0.075, size=11, text_color=T.MUTED, weight="normal", center=False):
    """A marble and a caption beside it. center: x is the middle of the pair."""
    tw = c.measure(text, size, weight)[0]
    if center:
        x -= (2 * r + 0.1 + tw) / 2
    _marble(c, x + r, y, r, colour)
    c.text(x + 2 * r + 0.1, y, text, size=size, color=text_color, weight=weight)
    return x + 2 * r + 0.1 + tw


def _panel_head(c, x, y, w, tag, title, size=12.5, floor=10.5):
    """A numbered panel header that may run to two lines, split where the two lines come out most even.
    Returns the number of lines used."""
    c.text(x + 0.2, y + 0.28, tag, size=12.5, weight="bold", color=T.NAVY)
    tx = x + 0.2 + c.measure(tag, 12.5, "bold")[0] + 0.12
    room = x + w - 0.18 - tx
    words = title.split(" ")
    while True:
        if c.measure(title, size, "bold")[0] <= room:
            lines = [title]
            break
        splits = [(" ".join(words[:k]), " ".join(words[k:])) for k in range(1, len(words))]
        width = lambda pair: max(c.measure(s, size, "bold")[0] for s in pair)
        best = min(splits, key=width)
        if width(best) <= room or size <= floor:
            lines = list(best)
            break
        size -= 0.5
    for k, line in enumerate(lines):
        c.text(tx, y + 0.28 + 0.215 * k, line, size=size, weight="bold", color=T.MUTED)
    return len(lines)


@figure("marbles", "The hypergeometric score as a bag of marbles, worked through for one gene and Bardet-Biedl syndrome "
                   "in three panels. One: model organism gene-phenotype annotations drop into a bag, one marble each, "
                   "401,464 in all; the gene ARL2BP owns 36 of them. Two: the disease's 109 significant phenologs "
                   "as a ranked table of model phenotypes with their P values and the marbles each one draws "
                   "(gastrulation disrupted 47, notochord kinked 51, male infertility 1,027, osmotic avoidance "
                   "defective 29, and 105 more); all 109 pool into 6,436 marbles spread across 3,075 candidate "
                   "genes, 14 of them ARL2BP's. Three: 401,464 marbles in the bag, 36 of them ARL2BP's; 6,436 "
                   "drawn, 14 of them ARL2BP's; the chance of drawing 14 or more at random, 2.0 x 10^-16, is the "
                   "score.")
def fig_marbles(D, step=None):
    """Three builds, one panel more each time."""
    step = 3 if step is None else step
    c = Canvas()
    b, ex = D.bbs, D.bbs.score_example
    cand = next(g for g in b.candidates if g.gene == ex.gene)
    rng = np.random.default_rng(3)
    y0, ph = 0.05, T.FIG_H - 0.1
    x1, w1 = 0.0, 3.8
    x2, w2 = 4.08, 4.34
    x3, w3 = 8.7, 3.3
    p = _hg_tail(ex.drawn, ex.bag, ex.own, ex.pulled)
    assert abs(p / cand.score - 1) < 0.05, "bbs.score_example no longer reproduces the score in bbs.candidates"
    assert sum(cand.support.values()) == ex.drawn == sum(len(v) for v in ex.carrying.values())
    assert [r.rank for r in ex.phenolog_rows] == sorted(r.rank for r in ex.phenolog_rows)

    # 1: what goes into the bag
    c.panel(x1, y0, w1, ph)
    dy = 0.2 * (_panel_head(c, x1, y0, w1, "1", f"Pool all {len(D.species)} model organism gene–phenotype annotations together") - 1)
    cx1, mx = x1 + w1 / 2, x1 + w1 - 0.36
    y = y0 + 0.7 + dy
    for k, a in enumerate(ex.annotations):
        sp = D.sp[a.species]
        c.badge(x1 + 0.34, y, sp.code, sp.color, r=0.125)
        c.text(x1 + 0.56, y, a.gene, size=10.5, weight="bold")
        gx = x1 + 0.56 + c.measure(a.gene, 10.5, "bold")[0] + 0.1
        size = c.fit(a.phenotype, mx - 0.28 - gx, 10.5)
        c.text(gx, y, a.phenotype, size=size, color=T.MUTED)
        c.line(gx + c.measure(a.phenotype, size)[0] + 0.1, y, mx - 0.16, y, color=T.FAINT, lw=1.0, ls=":")
        _marble(c, mx, y, 0.078, T.ACCENT if a.get("ours") else OTHER_MARBLES[k % len(OTHER_MARBLES)])
        y += 0.255
    c.text(x1 + 0.56, y - 0.01, f"… {n(ex.bag)} annotations, {len(D.species)} model species", size=10, color=T.FAINT,
           style="italic")
    bag_y = y + 0.52
    for dx, dy, col in ((-0.06, 0.2, OTHER_MARBLES[0]), (0.08, 0.33, T.ACCENT), (-0.03, 0.46, OTHER_MARBLES[2])):
        _marble(c, cx1 + dx, y + dy, 0.07, col)                    # marbles dropping into the bag
    cap_y = y0 + ph - 0.66
    body = _bag(c, cx1 - 1.3, bag_y, 2.6, cap_y - 0.2 - bag_y)
    pts, r = _pack(body, 330, r_max=0.075, top=0.13, margin=0.03, gap=0.012)
    for (px, py), col in zip(pts, _marble_colours(pts, [(T.ACCENT, 4)], rng, spread=9 * r)):
        _marble(c, px, py, r, col)
    cap = "one marble = one gene–phenotype annotation"
    c.text(cx1, cap_y, cap, size=c.fit(cap, w1 - 0.3, 11.5, "bold"), weight="bold", ha="center")
    own = c.wrap(f"{ex.gene} is associated with {ex.own} model organism phenotypes: {ex.own} marbles.", w1 - 0.78, 10.5)
    c.text(x1 + 0.5, cap_y + 0.19, own, size=10.5, color=T.MUTED, va="top", linespacing=1.25)
    _marble(c, x1 + 0.33, cap_y + 0.28, 0.075, T.ACCENT)

    ay = y0 + ph / 2                                               # the arrows from one panel to the next
    if step < 2:
        return c
    c.arrow(x1 + w1 + 0.04, ay, x2 - 0.04, ay, color=T.MUTED, lw=1.8, head=0.12)

    # 2: the disease's phenologs as the paper's Fig. 1C lists them, and what they pool into
    c.panel(x2, y0, w2, ph)
    dy = 0.2 * (_panel_head(c, x2, y0, w2, "2", "Sample from this bag of annotations via phenologs") - 1)
    tx2, rx2 = x2 + 0.22, x2 + w2 - 0.22
    c.text(tx2, y0 + 0.72 + dy, b.label, size=12, weight="bold")
    chip_x = rx2
    for key, cnt in reversed(list(b.by_species.items())):
        sp = D.sp[key]
        c.text(chip_x, y0 + 0.72 + dy, f"{cnt}", size=10.5, weight="bold", ha="right")
        chip_x -= c.measure(f"{cnt}", 10.5, "bold")[0] + 0.19
        c.badge(chip_x, y0 + 0.72 + dy, sp.code, sp.color, r=0.115)
        chip_x -= 0.115 + 0.16
    hy, pcol, mcx = y0 + 1.04 + dy, rx2 - 0.9, rx2 - 0.64
    c.text(tx2, hy, f"{n(b.phenologs)} significant phenologs, ranked by P", size=9.5, color=T.FAINT)
    c.text(pcol - 0.36, hy, "P", size=9.5, color=T.FAINT, ha="center")
    c.text(rx2, hy, "marbles drawn", size=9.5, color=T.FAINT, ha="right")

    def vdots(x, y):
        for dy in (-0.06, 0.0, 0.06):
            c.circle(x, y + dy, 0.012, fill=T.FAINT, edge=None)

    lab_size = min(c.fit(r.label, pcol - 0.86 - (tx2 + 0.32), 10.5, "bold") for r in ex.phenolog_rows)
    y, prev = y0 + 1.34 + dy, None
    for row in ex.phenolog_rows:
        if prev is not None and row.rank != prev + 1:              # a gap in the ranked list
            vdots(tx2 + 0.12, y - 0.07)
            y += 0.16
        sp = D.sp[row.species]
        c.badge(tx2 + 0.12, y, sp.code, sp.color, r=0.115)
        c.text(tx2 + 0.32, y, row.label, size=lab_size, weight="bold")
        c.text(pcol, y, sci(row.p), size=10, color=T.MUTED, ha="right")
        for k in (-1, 0, 1):
            _marble(c, mcx + k * 0.09, y, 0.042, T.ACCENT if row.get("ours") and k == 0 else OTHER_MARBLES[k + 1])
        c.text(rx2, y, n(row.marbles), size=11, weight="bold", ha="right")
        prev = row.rank
        y += 0.3 - dy / 8
    vdots(tx2 + 0.12, y - 0.07)
    c.text(tx2 + 0.32, y - 0.07, f"… {n(b.phenologs - len(ex.phenolog_rows))} more", size=10.5, color=T.FAINT,
           style="italic")
    c.text(tx2, y + 0.34, f"All {n(b.phenologs)} phenologs pool into", size=10.5, color=T.MUTED)
    tray_n = 23
    tray = _tray(c, rx2 - 0.86, y + 1.04, _tray_rows(tray_n), 0.075,
                 _marble_colours([(k, 0) for k in range(tray_n)], [(T.ACCENT, 2)], rng, spread=6))
    c.arrow(mcx, y + 0.02, mcx, tray.y - 0.06, color=T.MUTED, lw=1.6, head=0.11)
    c.text(tx2, y + 0.76, f"{n(ex.pulled)} marbles", size=19, weight="bold")
    c.text(tx2, y + 1.09, f"spread across {n(b.genes_ranked)}", size=10.5, color=T.MUTED)
    c.text(tx2, y + 1.3, "candidate genes", size=10.5, color=T.MUTED)
    _marble_row(c, tx2, y + 1.68, f"{ex.drawn} link to {ex.gene}", T.ACCENT, size=11.5, text_color=T.INK, weight="bold")

    if step < 3:
        return c
    c.arrow(x2 + w2 + 0.04, ay, x3 - 0.04, ay, color=T.MUTED, lw=1.8, head=0.12)

    # 3: the score, from the four numbers: the bag and the gene's share of it, the draw and the gene's share of that
    c.panel(x3, y0, w3, ph)
    dy = 0.2 * (_panel_head(c, x3, y0, w3, "3", "Test the significance of how often each gene shows up") - 1)
    tx, rx = x3 + 0.22, x3 + w3 - 0.22
    y = y0 + 1.0 + dy
    for lab, total, own in (("marbles in the bag", ex.bag, ex.own), ("marbles drawn", ex.pulled, ex.drawn)):
        c.text(tx, y, lab, size=11.5, color=T.MUTED)
        c.text(rx, y, n(total), size=13, weight="bold", ha="right")
        _marble(c, tx + 0.07, y + 0.36, 0.07, T.ACCENT)
        c.text(tx + 0.23, y + 0.36, ex.gene, size=11.5)
        c.text(rx, y + 0.36, n(own), size=13, weight="bold", ha="right")
        y += 1.1
    c.line(tx, y - 0.24, rx, y - 0.24, color=T.RULE, lw=1.2)
    head = f"Chance of drawing {ex.drawn} or more at random"
    c.text(tx, y + 0.2, head, size=c.fit(head, w3 - 0.44, 11.5, "bold"), weight="bold")
    c.text(tx, y + 0.72, f"P = {_p_text(p, bold=True)}", size=20, weight="bold", color=T.ACCENT_DARK)
    c.text(tx, y + 1.12, "hypergeometric test", size=12, weight="bold", color=T.ACCENT_DARK)
    return c


@figure("marbles_two_genes", "One bag of 120 marbles and one handful of 20, with two genes marked. A specific gene has "
                             "six marbles in the bag and a hub gene has thirty; both have five in the handful. Two "
                             "cards give the verdicts: five of six is far more than chance predicts, five of thirty "
                             "is exactly what chance predicts. A box below says every gene is judged against the same "
                             "bag, giving one ranked list of 1,485,624 gene-disease pairs.")
def fig_marbles_two_genes(D, step=None):
    c = Canvas()
    M = D.hg_illustration.marbles
    rng = np.random.default_rng(8)
    genes = [(M.specific, "A specific gene", "a specific gene", T.ACCENT, T.ACCENT_DARK),
             (M.hub, "A hub gene", "a hub gene", T.NAVY, T.NAVY)]
    share = _share(M.handful, M.bag)

    # left: the bag and the handful, both genes marked
    lw_ = 5.95
    c.panel(0.0, 0.05, lw_, 4.95, title="One bag, one handful")
    body = _bag(c, 0.25, 0.55, 3.0, 3.15)
    pts, r = _pack(body, M.bag)
    for (x, y), col in zip(pts, _marble_colours(pts, [(col_, g.own) for g, _, _, col_, _ in genes], rng,
                                                spread=5 * r)):
        _marble(c, x, y, r, col)
    tcx, tcy = 4.66, 2.5
    tcols = _marble_colours([(k, 0) for k in range(M.handful)], [(col_, g.drawn) for g, _, _, col_, _ in genes], rng,
                            spread=1.5)
    tray = _tray(c, tcx, tcy, _tray_rows(M.handful), 0.12, tcols)
    c.arrow(3.31, tcy, tray.x - 0.07, tcy, color=T.MUTED, lw=1.8, head=0.12)
    c.text(tcx, tray.y - 0.52, "the disease's phenologs", size=11.5, weight="bold", ha="center")
    c.text(tcx, tray.y - 0.28, "draw a handful", size=11.5, weight="bold", ha="center")
    c.text(tcx, tray.bottom + 0.26, f"{M.handful} marbles: {share} of the bag", size=11, color=T.MUTED, ha="center")
    c.text(1.75, 3.95, "every gene–phenotype annotation", size=11.5, weight="bold", ha="center")
    c.text(1.75, 4.19, f"in the seven model species: {n(M.bag)} marbles", size=11, color=T.MUTED, ha="center")
    lx = 0.25
    for g, _, low, col_, _ in genes:
        lx = _marble_row(c, lx, 4.68, f"{low}: {g.own} marbles", col_, size=10.5) + 0.24
    _marble_row(c, lx, 4.68, "all other genes", OTHER_MARBLES[0], size=10.5)

    # right: the two verdicts
    gx = 6.25
    right = T.FIG_W - 0.03
    gw = (right - gx - 0.2) / 2
    c.text(gx, 0.28, "Same count in the handful, different verdict", size=13.5, weight="bold")
    c.text(right, 0.28, "made-up numbers", size=10, color=T.FAINT, ha="right", style="italic")
    cy0, ch = 0.62, 3.05
    for k, (g, label, _, col_, dark) in enumerate(genes):
        x = gx + k * (gw + 0.2)
        hot = k == 0
        expected = M.handful * g.own / M.bag
        p = _hg_tail(g.drawn, M.bag, g.own, M.handful)
        c.box(x, cy0, gw, ch, fill=T.ACCENT_TINT if hot else T.WHITE, edge=T.ACCENT if hot else T.RULE, lw=1.6, r=0.12)
        c.text(x + 0.22, cy0 + 0.33, label, size=13.5, weight="bold", color=dark)
        c.text(x + 0.22, cy0 + 0.61, f"{g.own} marbles in the bag", size=11, color=T.MUTED)
        mr, per = 0.085, 10
        pitch = (gw - 0.44 - 2 * mr) / (per - 1)
        for j in range(g.own):                                    # its marbles: filled came out, hollow stayed
            _marble(c, x + 0.22 + mr + (j % per) * pitch, cy0 + 0.98 + (j // per) * 0.215, mr, col_, hollow=j >= g.drawn)
        c.text(x + 0.22, cy0 + 1.8, f"{g.drawn} of its {g.own} came out", size=11.5, weight="bold")
        c.text(x + 0.22, cy0 + 2.05, f"chance predicts {expected:g}", size=11, color=T.MUTED)
        c.text(x + 0.22, cy0 + 2.47, f"P = {_p_text(p, bold=True)}", size=17, weight="bold", color=dark if hot else T.MUTED)
        c.text(x + 0.22, cy0 + 2.8, "a strong candidate" if hot else "just what chance predicts", size=11.5,
               weight="bold" if hot else "normal", color=dark if hot else T.MUTED)

    by = cy0 + ch + 0.18
    c.box(gx, by, right - gx, 5.0 - by, fill=T.NAVY, edge=None, r=0.12)
    c.text(gx + 0.28, by + 0.3, "Every gene is judged against the same bag,", size=12.5, color=LIGHT)
    c.text(gx + 0.28, by + 0.56, "so scores compare across diseases:", size=12.5, color=LIGHT)
    c.text(gx + 0.28, by + 0.86, f"one ranked list of {n(D.totals.ranked_pairs)} gene–disease pairs", size=13.5,
           weight="bold", color=T.WHITE)
    return c


# =============================================================================================
# 6. Evaluation
# =============================================================================================
@figure("time_split_simple", "How the scores were evaluated. Two versions of the Monarch KG one year apart: on the 2025 "
                             "version phenologs are computed and candidate genes scored per disease; the 2026 version "
                             "adds new gene-disease pairs. Below: 73 newly added gene-disease pairs across 43 "
                             "diseases can be evaluated; for each, the rank of the gene in that disease's predictions "
                             "from one year ago. Scoring methods evaluated: naive Bayes, hypergeometric, and a random "
                             "gene for comparison.")
def fig_time_split_simple(D, step=None):
    c = Canvas()
    e, m = D.evaluation, D.meta
    bw, y0, hh = 5.0, 0.42, 1.5
    lx, rx = 0.0, c.w - bw

    # the two versions of the graph
    c.box(lx, y0, bw, hh, fill=T.TINT, edge=None, r=0.14)
    c.text(lx + 0.32, y0 + 0.52, f"Monarch KG {m.kg_prior}", size=18, weight="bold")
    a = "compute phenologs, score candidate genes per disease"
    c.text(lx + 0.32, y0 + 1.02, a, size=c.fit(a, bw - 0.64, 13), color=T.MUTED)
    c.box(rx, y0, bw, hh, fill=T.NAVY, edge=None, r=0.14)
    c.text(rx + 0.32, y0 + 0.52, f"Monarch KG {m.kg_release}", size=18, weight="bold", color=T.WHITE)
    c.text(rx + 0.32, y0 + 1.02, "newly added gene–disease pairs", size=13, color=LIGHT)
    c.arrow(lx + bw + 0.2, y0 + hh / 2 + 0.12, rx - 0.2, y0 + hh / 2 + 0.12, color=T.MUTED, lw=2.2, head=0.15)
    c.text(c.w / 2, y0 + hh / 2 - 0.16, "one year", size=13.5, weight="bold", color=T.MUTED, ha="center")

    # what is compared, and the scoring methods it is done for
    sy, sh = y0 + hh + 0.55, 2.35
    c.box(0.0, sy, c.w, sh, fill=T.NAVY, edge=None, r=0.14, zorder=0)
    pw, gap, pill_h = 1.85, 0.1, 0.72
    px = c.w - 0.3 - 3 * pw - 2 * gap
    big = n(e.genuine)
    c.text(0.4, sy + sh / 2, big, size=44, weight="bold", color=T.WHITE)
    tx = 0.4 + c.measure(big, 44, "bold")[0] + 0.3
    tw = px - 0.45 - tx
    head = c.wrap(f"newly added gene–disease pairs that can be evaluated, across {n(e.genuine_diseases)} diseases", tw, 13.5, "bold")
    sub = c.wrap("For each pair, find the rank of the gene in that disease's predictions from one year ago.", tw, 12)
    hl, sl = head.count("\n") + 1, sub.count("\n") + 1
    ty = sy + sh / 2 - (0.25 * hl + 0.12 + 0.22 * sl) / 2
    c.text(tx, ty, head, size=13.5, weight="bold", color=T.WHITE, va="top", linespacing=1.25)
    c.text(tx, ty + 0.25 * hl + 0.12, sub, size=12, color=LIGHT, va="top", linespacing=1.25)
    c.line(px - 0.25, sy + 0.35, px - 0.25, sy + sh - 0.35, color="#3B5580", lw=1.2, zorder=3)
    my = sy + sh / 2 + 0.2
    c.text(px, my - pill_h / 2 - 0.28, "Scoring methods evaluated", size=12.5, weight="bold", color=T.WHITE)
    for i, (name, who, ours) in enumerate([("naive Bayes", "McGary 2010, Woods 2013", False),
                                           ("hypergeometric", "this work", True), ("random gene", "for comparison", False)]):
        x = px + i * (pw + gap)
        c.box(x, my - pill_h / 2, pw, pill_h, fill=T.ACCENT if ours else "#24406B", edge=None, r=0.1)
        c.text(x + 0.15, my - 0.13, name, size=12.5, weight="bold", color=T.WHITE)
        c.text(x + 0.15, my + 0.15, who, size=9.5, color=T.WHITE if ours else LIGHT)
    return c


@figure("time_split", "Two snapshots of the Monarch KG a year apart. Of 524 gene-disease pairs added, 402 are not "
                      "testable, 49 are redundant within Mondo, and 73 are genuine new pairs used for evaluation.")
def fig_time_split(D, step=None):
    c = Canvas()
    e, m = D.evaluation, D.meta
    bw, bh = 4.3, 1.3
    c.box(0.3, 0.2, bw, bh, fill=T.TINT, edge=None, r=0.12)
    c.text(0.55, 0.52, f"Monarch KG {m.kg_prior}", size=15, weight="bold")
    c.text(0.55, 0.9, "build phenologs, rank candidate genes", size=12, color=T.MUTED)
    c.text(0.55, 1.2, "the prediction", size=11, color=T.CAROLINA_DARK, weight="bold")
    c.box(7.4, 0.2, bw, bh, fill=T.NAVY, edge=None, r=0.12)
    c.text(7.65, 0.52, f"Monarch KG {m.kg_release}", size=15, weight="bold", color=T.WHITE)
    c.text(7.65, 0.9, f"{n(e.new_pairs)} gene–disease pairs added since", size=12, color=LIGHT)
    c.text(7.65, 1.2, "the answer key", size=11, color=T.WHITE, weight="bold")
    c.arrow(0.3 + bw + 0.15, 0.85, 7.25, 0.85, color=T.MUTED, lw=2.0, head=0.14)
    c.text((0.3 + bw + 7.4) / 2, 0.62, "one year", size=12.5, weight="bold", color=T.MUTED, ha="center")

    x0, y0, W, H = 0.3, 2.45, 11.4, 0.72
    c.text(x0, 2.1, f"Which of the {n(e.new_pairs)} new pairs can fairly test a prediction?", size=13.5, weight="bold")
    segs = [(e.not_testable, "#C4CED8", T.INK), (e.redundant, T.CAROLINA, T.WHITE), (e.genuine, T.ACCENT, T.WHITE)]
    x = x0
    edges = []
    for v, col, fg in segs:
        w = W * v / e.new_pairs
        c.rect(x, y0, w, H, fill=col)
        c.text(x + w / 2, y0 + H / 2, n(v), size=17, weight="bold", color=fg, ha="center")
        edges.append((x, x + w))
        x += w
    c.text(edges[0][0], y0 + H + 0.3, "not testable", size=13, weight="bold", color=T.MUTED)
    c.text(edges[0][0], y0 + H + 0.58, "the disease had no gene a year earlier, so no phenolog could exist", size=11.5,
           color=T.MUTED)
    mid = (edges[1][0] + edges[1][1]) / 2
    c.line(mid, y0 + H, mid, y0 + H + 0.95, color=T.CAROLINA_DARK, lw=1.3)
    c.text(mid - 0.12, y0 + H + 1.05, "redundant within Mondo", size=13, weight="bold", color=T.CAROLINA_DARK, ha="right")
    c.text(mid - 0.12, y0 + H + 1.33, "gene already on a parent, child or sibling term  (MLH1: Lynch syndrome → Lynch syndrome 1)",
           size=11.5, color=T.MUTED, ha="right")
    mid2 = (edges[2][0] + edges[2][1]) / 2
    c.line(mid2, y0 + H, mid2, y0 + H + 1.72, color=T.ACCENT, lw=1.3)
    c.text(mid2 - 0.12, y0 + H + 1.82, f"{n(e.genuine)} genuine new pairs across {n(e.genuine_diseases)} diseases: the test set",
           size=14, weight="bold", color=T.ACCENT, ha="right")
    return c


@figure("eval_results", "The published recall curves for the temporal validation, each under its headline. Left, "
                        "looking at genes within each of the 43 diseases: the top 100 candidates per disease recover "
                        "9 of 73 new pairs, 23.5 times a random gene; the hypergeometric score and naive Bayes "
                        "perform the same. Right, pooled across all gene-disease pairs: the top 10,000 predictions, "
                        "0.7% of the list, recover 4 new pairs, 14 times a random ordering; the hypergeometric score "
                        "can rank each gene-disease association across the entire set of predictions, AUC 0.58 "
                        "against 0.45 for naive Bayes.")
def fig_eval_results(D, step=None):
    """Two builds: 1 the per-disease panel, 2 the pooled panel beside it. The slide has no lede: each panel's
    header carries its own headline (round 2 edits)."""
    from PIL import Image
    step = 2 if step is None else step
    c = Canvas(h=T.FIG_H_TALL)
    e = D.evaluation
    top = next(r for r in e.recall if r.n == 100)
    bold = {"weight": "bold"}
    cw, gap, hh, ih = 5.85, 0.3, 1.22, 3.42
    cols = [("fig3a", f"Looking at genes within each of the {n(e.genuine_diseases)} diseases",
             [(f"The top {top.n} candidates per disease recover", None), (f"{top.hits} of {n(e.genuine)} new pairs:", bold),
              (f"{top.fold:g} times a random gene.", None)],
             f"hg and naive Bayes perform the same (P = {e.hg_vs_nb_p})"),
            ("fig3b", "Pooled across all gene–disease pairs",
             [(f"The top {n(e.pooled_top)} predictions ({e.pooled_top_pct}% of the list) recover", None),
              (f"{e.pooled_hits} new pairs:", bold), (f"{e.pooled_fold:g} times a random ordering.", None)],
             "hg can rank each gene–disease association across the entire\n"
             f"set of predictions (AUC {e.pooled_auc_hg} vs {e.pooled_auc_nb} for naive Bayes)")]
    for k, (crop, head, headline, result) in enumerate(cols):
        if k + 1 > step:
            continue
        x = k * (cw + gap)
        c.box(x, 0.05, cw, hh, fill=T.TINT, edge=None, r=0.12)
        c.text(x + 0.25, 0.36, head, size=c.fit(head, cw - 0.5, 14, "bold"), weight="bold")
        c.rich(x + 0.25, 0.86, headline, cw - 0.5, size=12.5)
        img = paper.crop(crop, CROP_DIR)
        iw_, ih_ = Image.open(img).size
        b = c.image(img, x + (cw - ih * iw_ / ih_) / 2, 0.05 + hh + 0.1, h=ih)
        c.text(x + cw / 2, b.bottom + 0.12, result, size=12, weight="bold", ha="center", va="top", linespacing=1.25)
    c.text(0.05, c.h - 0.1, f"hg: hypergeometric score. k = {n(D.totals.k_phenologs)} phenologs per disease, "
                            f"{n(e.genuine)} new pairs.", size=10, color=T.FAINT)
    return c


@figure("eval_results_tiles", "Published recall curves for the temporal validation, per disease and pooled, beside the "
                              "headline numbers: 9 of 73 new pairs in the top 100, 23.5 times random, and FAAP100.")
def fig_eval_results_tiles(D, step=None):
    c = Canvas()
    e = D.evaluation
    r100 = next(r for r in e.recall if r.n == 100)
    b = c.image(paper.crop("fig3", CROP_DIR), 0.0, 0.25, w=7.55)
    cap = c.wrap(f"hg score, k = {n(D.totals.k_phenologs)} phenologs per disease. Left: recall of the {e.genuine} new pairs "
                 "within each disease's top N candidates. Right: every disease pooled into one ranked list.", 7.3, 11)
    c.text(0.15, b.bottom + 0.2, cap, size=11, color=T.MUTED, va="top", linespacing=1.25)

    px, pw = 7.95, 4.05
    c.box(px, 0.1, pw, 1.72, fill=T.NAVY, edge=None, r=0.12)
    c.text(px + 0.25, 0.52, f"{r100.hits} of {e.genuine}", size=30, weight="bold", color=T.WHITE)
    c.text(px + 0.25, 1.02, "new pairs sit in their disease's top 100", size=12.5, weight="bold", color=T.WHITE)
    c.text(px + 0.25, 1.4, f"{r100.fold}× a random gene", size=12.5, weight="bold", color=T.WHITE)
    c.text(px + 0.25, 1.63, f"95% CI {e.top100_ci} · P = {sci(e.top100_p)}", size=10.5, color=LIGHT)

    c.box(px, 1.97, pw, 1.5, fill=T.TINT, edge=None, r=0.12)
    c.text(px + 0.25, 2.25, "Within a disease", size=12.5, weight="bold")
    c.text(px + 0.25, 2.5, f"hg and naive Bayes tie (P = {e.hg_vs_nb_p})", size=11.5, color=T.MUTED)
    c.text(px + 0.25, 2.88, "Pooled across diseases", size=12.5, weight="bold")
    c.text(px + 0.25, 3.13, f"only hg beats random: AUC {e.pooled_auc_hg} vs {e.pooled_auc_nb}", size=11.5, color=T.MUTED)

    p = e.prospective
    c.box(px, 3.62, pw, 1.38, fill=T.ACCENT_TINT, edge=T.ACCENT, lw=1.5, r=0.12)
    c.text(px + 0.25, 3.92, f"{p.gene} → {p.disease}", size=13.5, weight="bold", color=T.ACCENT)
    c.text(px + 0.25, 4.2, c.wrap(f"Absent from the 2025 graph and ranked on mouse phenologs alone. "
                                  f"The human link was reported in {p.reported}.", pw - 0.5, 11),
           size=11, color=T.INK, va="top", linespacing=1.25)
    return c


@figure("recover_known", "Two tests on known disease genes. Left: leave-one-out cross validation in three steps. One "
                         "protein family's genes are removed from the human gene-disease network, every phenolog and "
                         "every candidate ranking is recomputed, and the rank of the hidden genes is read from the "
                         "disease's candidate list; repeated for each of 2,001 protein families. Right: for the 414 "
                         "diseases with at least four known genes, with nothing hidden, a ranked candidate list is "
                         "turned into a precision-recall curve and the area under it is taken, once per score.")
def fig_recover_known(D, step=None):
    c = Canvas()
    L = D.loo
    hide, cut, grey = T.ACCENT, T.ACCENT_DARK, "#C4CED8"
    dash = (0, (2.2, 1.8))

    # ---- a) leave-one-out: three steps and a loop
    aw = 7.0
    c.panel(0, 0, aw, c.h - 0.02, tag="a)", title="Leave-one-out: hide a protein family, then look for its genes")
    steps = [(0.25, "Hide one family", "human gene–disease network"),
             (2.78, "Recompute", "phenologs, then rankings"),
             (4.92, "Read their rank", "a disease's candidates")]
    for i, (x, head, sub) in enumerate(steps):
        c.circle(x + 0.16, 0.9, 0.16, fill=T.NAVY, edge=None)
        c.text(x + 0.16, 0.91, str(i + 1), size=11, weight="bold", color=T.WHITE, ha="center")
        c.text(x + 0.42, 0.8, head, size=13, weight="bold")
        c.text(x + 0.42, 1.03, sub, size=9.5, color=T.MUTED)

    # 1: one disease and its genes. Two of them are paralogs in one protein family, and both are taken out.
    dx, dy, R = 1.1, 2.3, 0.74

    def at(deg):
        return dx + R * math.cos(math.radians(deg)), dy + R * math.sin(math.radians(deg))

    for deg in (132, 168, 204, 240):
        gx, gy = at(deg)
        c.line(dx, dy, gx, gy, color=T.NAVY, lw=1.4, zorder=2)
        c.circle(gx, gy, 0.09, fill=T.NAVY, edge=T.WHITE, lw=0.8, zorder=3)
    fam = [at(-31), at(29)]
    fx0, fx1 = min(g[0] for g in fam) - 0.2, max(g[0] for g in fam) + 0.2
    fy0, fy1 = fam[0][1] - 0.22, fam[1][1] + 0.22
    c.box(fx0, fy0, fx1 - fx0, fy1 - fy0, fill=T.ACCENT_TINT, edge=hide, lw=1.5, r=0.2, ls=(0, (4, 2.5)), zorder=1)
    c.text((fx0 + fx1) / 2 + 0.12, fy0 - 0.16, "one protein family", size=9.5, weight="bold", color=cut, ha="center")
    for gx, gy in fam:
        ex, ey = fx0 - 0.03, dy + (fx0 - 0.03 - dx) * (gy - dy) / (gx - dx)       # where the edge meets the family's box
        c.line(dx, dy, ex, ey, color=T.FAINT, lw=1.3, ls=dash, zorder=2)
        mx, my, s = (dx + ex) / 2 + 0.06, (dy + ey) / 2 + 0.06 * (gy - dy) / (gx - dx), 0.055
        c.line(mx - s, my - s, mx + s, my + s, color=cut, lw=2.0, zorder=5)
        c.line(mx - s, my + s, mx + s, my - s, color=cut, lw=2.0, zorder=5)
        c.circle(gx, gy, 0.09, fill=T.WHITE, edge=hide, lw=1.6, zorder=6)
    c.square(dx, dy, 0.3, fill=T.NAVY, edge=T.NAVY, zorder=4)
    ky = 3.2
    c.square(0.36, ky, 0.14, fill=T.NAVY, edge=T.NAVY)
    c.text(0.49, ky, "disease", size=9.5, color=T.MUTED)
    c.circle(1.14, ky, 0.065, fill=T.NAVY, edge=None)
    c.text(1.26, ky, "gene", size=9.5, color=T.MUTED)
    c.circle(1.77, ky, 0.065, fill=T.WHITE, edge=hide, lw=1.4)
    c.text(1.89, ky, "hidden", size=9.5, color=T.MUTED)
    c.text(0.27, ky + 0.24, c.wrap("The whole family goes, so a paralog cannot stand in for a hidden gene.", 2.25, 9.5),
           size=9.5, color=T.MUTED, va="top", linespacing=1.25)

    # 2: everything downstream is recomputed
    cx2 = 3.68
    for j, key in enumerate(T.SPECIES_ORDER):
        c.badge(cx2 + (j - 3) * 0.26, 1.5, D.sp[key].code, D.sp[key].color, r=0.12)
    c.text(cx2, 1.76, "against seven model species", size=9.5, color=T.MUTED, ha="center")
    cols, rows, pitch, cell = 7, 4, 0.19, 0.15
    gx0, gy0 = cx2 - (cols * pitch - (pitch - cell)) / 2, 2.0
    hits = {(1, 0), (4, 1), (2, 3), (6, 2)}
    for i in range(cols):
        for j in range(rows):
            c.rect(gx0 + i * pitch, gy0 + j * pitch, cell, cell, fill=T.NAVY if (i, j) in hits else "#C9D5E2", zorder=2)
    c.text(cx2, gy0 + rows * pitch + 0.13, "every phenolog, from scratch", size=9.5, color=T.MUTED, ha="center")
    py = 3.42
    c.arrow(cx2, gy0 + rows * pitch + 0.28, cx2, py - 0.2, color=T.MUTED, lw=1.5, head=0.09)
    for k, lab in enumerate(["nb", "hg"]):
        c.pill(cx2 + (k - 0.5) * 0.74, py, lab, w=0.62, h=0.3, size=10.5, weight="bold", lw=1.3)
    c.text(cx2, py + 0.27, "every disease re-ranked,\nwith both scores", size=9.5, color=T.MUTED, ha="center", va="top",
           linespacing=1.25)

    # 3: where the hidden genes come back in the disease's candidate list
    lx, y0, rh = 4.97, 1.46, 0.25
    lens = [1.2, 1.08, 0.98, 0.88, 0.8, 0.72, 0.65, 0.58]
    for i, ln in enumerate(lens):
        y, hot = y0 + i * rh, i in (2, 5)
        c.text(lx + 0.2, y, str(i + 1), size=9.5, color=cut if hot else T.FAINT, ha="right", weight="bold" if hot else "normal")
        c.box(lx + 0.3, y - 0.085, ln, 0.17, fill=hide if hot else grey, edge=None, r=0.04)
        if hot:
            c.text(lx + 0.3 + ln + 0.08, y, "hidden", size=9.5, weight="bold", color=cut)
    c.text(lx + 0.3, y0 + len(lens) * rh - 0.04, "…", size=12, color=T.FAINT)
    c.text(lx + 0.02, y0 + len(lens) * rh + 0.22, "Did the hidden genes come\nback, and how high?", size=9.5,
           color=T.MUTED, va="top", linespacing=1.25)
    for x in (2.52, 4.66):
        c.arrow(x, 2.3, x + 0.2, 2.3, color=T.MUTED, lw=1.8, head=0.12)

    # the loop
    ly, lx0, lx1 = 4.42, dx, 5.85
    c.line(lx1, 4.12, lx1, ly, color=T.MUTED, lw=1.6)
    c.line(lx1, ly, lx0, ly, color=T.MUTED, lw=1.6)
    c.arrow(lx0, ly, lx0, 4.12, color=T.MUTED, lw=1.6, head=0.12)
    c.pill((lx0 + lx1) / 2, ly, f"repeated for each of {n(L.families_hidden)} protein families", h=0.4, size=12.5,
           weight="bold", fill=T.NAVY, edge=None, color=T.WHITE, pad=0.26)
    c.text((lx0 + lx1) / 2, ly + 0.38, "every family that carries a human disease gene  ·  after Woods et al. 2013",
           size=9.5, color=T.MUTED, ha="center")

    # ---- b) precision of a ranking: the complete model, nothing hidden
    px, pw = aw + 0.2, c.w - aw - 0.2
    c.panel(px, 0, pw, c.h - 0.02, tag="b)", title="Precision: are known genes near the top?")
    c.pill_left(px + 0.25, 0.86, f"{L.compared_diseases} diseases with ≥ {L.min_known_genes} known genes", size=11.5,
                weight="bold", h=0.36, pad=0.16, lw=1.4)
    c.text(px + 0.27, 1.17, "the complete model: no genes are hidden here", size=9.5, color=T.MUTED, style="italic")

    ex = L.pr_example
    known = list(ex.known_at)
    ry0, rdy = 1.82, 0.245
    c.text(px + 0.45, ry0 - 0.3, "rank", size=9.5, color=T.FAINT, ha="right")
    c.text(px + 0.86, ry0 - 0.3, "precision", size=9.5, color=T.FAINT)
    pts, found = [], 0
    for i in range(ex.shown):
        y, rank = ry0 + i * rdy, i + 1
        hot = rank in known
        c.text(px + 0.45, y, str(rank), size=9.5, color=T.FAINT, ha="right")
        c.circle(px + 0.64, y, 0.075, fill=T.NAVY if hot else T.WHITE, edge=T.NAVY if hot else T.FAINT, lw=1.2)
        if hot:
            found += 1
            c.text(px + 0.86, y, f"{found} of {rank}", size=10, weight="bold")
            pts.append((found / len(known), found / rank, y))
    ky2 = ry0 + ex.shown * rdy + 0.12
    c.circle(px + 0.32, ky2, 0.065, fill=T.NAVY, edge=None)
    c.text(px + 0.44, ky2, "known gene", size=9.5, color=T.MUTED)
    c.circle(px + 1.3, ky2, 0.065, fill=T.WHITE, edge=T.FAINT, lw=1.2)
    c.text(px + 1.42, ky2, "other", size=9.5, color=T.MUTED)

    # the curve those four points draw: precision holds between known genes, then steps down
    ox, oy, gw, gh = px + 2.62, 3.5, 1.86, 1.72
    c.box(ox - 0.08, oy - gh - 0.14, gw + 0.2, gh + 0.2, fill=T.WHITE, edge=None, r=0.08, zorder=1)

    def X(r): return ox + r * gw
    def Y(p): return oy - p * gh

    curve = [(0.0, pts[0][1])]
    for r, p, _ in pts:
        curve += [(curve[-1][0], p), (r, p)]
    xy = [(X(r), Y(p)) for r, p in curve]
    c.ax.add_patch(Polygon(xy + [(X(1), Y(0)), (X(0), Y(0))], closed=True, facecolor=T.TINT_2, edgecolor="none", zorder=2))
    c.ax.plot([q[0] for q in xy], [q[1] for q in xy], color=T.NAVY, lw=2.2, zorder=4, solid_joinstyle="miter")
    for r, p, _ in pts:
        c.circle(X(r), Y(p), 0.055, fill=T.NAVY, edge=T.WHITE, lw=0.8, zorder=5)
    c.line(ox, oy, X(1) + 0.06, oy, color=T.NAVY, lw=1.3, zorder=3)
    c.line(ox, oy, ox, Y(1) - 0.08, color=T.NAVY, lw=1.3, zorder=3)
    c.text(ox - 0.08, Y(1), "1", size=9, color=T.MUTED, ha="right")
    c.text(X(1), oy + 0.15, "1", size=9, color=T.MUTED, ha="center")
    c.text(X(0), oy + 0.15, "0", size=9, color=T.MUTED, ha="center")
    c.text(ox - 0.04, Y(1) - 0.24, "precision", size=10, color=T.MUTED)
    c.text(X(0.5), oy + 0.17, "recall: known genes found", size=10, color=T.MUTED, ha="center")
    c.text(X(0.4), Y(0.33), "AUPRC", size=13.5, weight="bold", ha="center")
    c.text(X(0.4), Y(0.33) + 0.22, "area under the curve", size=9.5, color=T.MUTED, ha="center")
    third = pts[2]
    c.scurve(px + 1.42, third[2], X(third[0]) - 0.07, Y(third[1]), color=T.MUTED, lw=1.1, ls=dash, zorder=3)

    c.text(px + 0.25, 4.42, "One AUPRC per disease, for each score", size=12.5, weight="bold")
    c.text(px + 0.25, 4.72, "1 would put every known gene at the very top of the list", size=10, color=T.MUTED)
    return c


@figure("precision", "Published plots comparing the two gene scores across the 414 diseases with at least four known "
                     "genes, with no genes hidden: mean recall of known genes within the top N candidates, and a "
                     "scatter of per-disease area under the precision-recall curve, naive Bayes against "
                     "hypergeometric. Beside them: hg is more precise in 67% of diseases, finds 45% against 39% of "
                     "known genes in the top 20, and in the leave-one-out runs places the true disease in the top 20 "
                     "for 165 diseases against 12.")
def fig_precision(D, step=None):
    c = Canvas()
    L = D.loo
    b = c.image(paper.crop("figS7cd", CROP_DIR), 0.0, 0.12, w=7.75)
    cap = c.wrap(f"The complete model, no genes hidden: k = {L.precision_k} phenologs per disease, {L.compared_diseases} diseases with at "
                 f"least {L.min_known_genes} known genes. Left: mean share of a disease's known genes within its top N candidates. "
                 "Right: area under the precision–recall curve (AUPRC), one point per disease.", 7.5, 11)
    c.text(0.15, b.bottom + 0.2, cap, size=11, color=T.MUTED, va="top", linespacing=1.25)

    px, pw = 8.0, 4.0
    _tile(c, px, 0.08, pw, 1.42, f"{L.more_precise_pct}%", "of diseases: hg is more precise",
          f"{L.more_precise_diseases} of {L.compared_diseases}; median AUPRC gain {L.auprc_median_gain}", fill=T.NAVY,
          fg=T.WHITE, sub_fg=LIGHT, big_size=28)
    _tile(c, px, 1.63, pw, 1.42, f"{L.recall_top20_hg}% vs {L.recall_top20_nb}%", "of known genes are in the top 20",
          f"hg vs naive Bayes. {L.recall_ceiling_pct}% are in the list at any rank.", big_size=28)
    _tile(c, px, 3.18, pw, 1.88, f"{L.true_disease_top20_hg} vs {L.true_disease_top20_nb}",
          "leave-one-out: diseases where the hidden gene's true disease ranks in its top 20",
          f"hg vs naive Bayes, k = {L.k}. Within a disease the two tie.", big_size=28)
    return c


@figure("species_contrib", "Bars of genuine new gene-disease pairs reached by each species alone: mouse 38, fission "
                           "yeast 7, worm 3, zebrafish 3, each with at least one pair no other species reaches.")
def fig_species_contrib(D, step=None):
    c = Canvas()
    e = D.evaluation
    c.text(0.1, 0.25, f"Genuine new pairs each species' own ranking reaches (of {e.genuine})", size=13.5, weight="bold")
    x0, scale = 2.3, 3.0 / e.ceiling
    y = 0.95
    for row in e.species_recall:
        sp = D.sp[row.species]
        c.badge(0.3, y, sp.code, sp.color, r=0.19)
        c.text(0.62, y, sp.common, size=13, weight="bold")
        c.box(x0, y - 0.19, row.reached * scale, 0.38, fill=sp.color, edge=None, r=0.05, alpha=0.38)
        c.box(x0, y - 0.19, row.unique * scale, 0.38, fill=sp.color, edge=None, r=0.05)
        c.text(x0 + row.reached * scale + 0.12, y, f"{row.reached}", size=14, weight="bold")
        c.text(x0 + row.reached * scale + 0.12 + c.measure(f"{row.reached}", 14, "bold")[0] + 0.1, y,
               f"{row.unique} reached by this species alone", size=11, color=T.MUTED)
        y += 0.66
    c.text(x0, y - 0.02, f"{e.ceiling} of {e.genuine} are reached by at least one species", size=12, weight="bold")
    c.rect(x0, y + 0.42, 0.3, 0.16, fill=T.MUTED)
    c.text(x0 + 0.4, y + 0.5, "only this species", size=10.5, color=T.MUTED)
    c.rect(x0 + 1.9, y + 0.42, 0.3, 0.16, fill=T.MUTED, alpha=0.38)
    c.text(x0 + 2.3, y + 0.5, "also reached by another", size=10.5, color=T.MUTED)
    c.text(0.1, 4.75, e.silent_species + ".", size=11.5, color=T.FAINT, style="italic")

    px, pw = 8.0, 4.0
    _tile(c, px, 0.1, pw, 1.55, f"+{e.extra_over_mouse}", "pairs reached only when all seven species are pooled",
          "at the same ranking quality as mouse alone", fill=T.NAVY, fg=T.WHITE, sub_fg=LIGHT)
    zf = next(r for r in e.species_recall if r.species == "zebrafish")
    _tile(c, px, 1.8, pw, 1.45, f"AUC {zf.auc}", "zebrafish: few hits, ranked sharply",
          "its three pairs sit in the top 6.5% of its list")
    _tile(c, px, 3.4, pw, 1.6, "ATM", "fission yeast: hereditary breast ovarian cancer syndrome",
          "top 1.2% of its list: DNA repair, conserved to yeast", big_size=26)
    return c


@figure("bbs_ranking", "The top 25 ranked candidate genes for Bardet-Biedl syndrome. The first eleven are all known "
                       "genes; new candidates such as IFT140 and ARL2BP rank among them.")
def fig_bbs_ranking(D, step=None):
    c = Canvas()
    b = D.bbs
    cands = b.candidates
    # the slide's lede quotes the known genes' ranks: keep them in step with the list drawn here
    ranks = sorted(b.known_ranks.values())
    assert len(ranks) == b.known_genes == b.known_recovered and ranks[-2:] == [b.known_rank_all_but_one, b.known_rank_last]
    assert {g.gene: g.rank for g in cands if g.known} == {g: r for g, r in b.known_ranks.items() if r <= len(cands)}
    top = max(-math.log10(g.score) for g in cands)
    per_col = math.ceil(len(cands) / 2)
    col_w = 6.0
    c.rect(0.1, 0.13, 0.18, 0.16, fill=T.NAVY)
    c.text(0.36, 0.21, "known Bardet–Biedl gene", size=11.5, color=T.MUTED)
    c.rect(2.75, 0.13, 0.18, 0.16, fill=T.ACCENT)
    c.text(3.01, 0.21, "new candidate", size=11.5, color=T.MUTED)
    c.text(4.5, 0.21, "bar: −log10 hg score     numbers: phenologs in each species that carry an ortholog of the gene", size=10.5,
           color=T.FAINT)
    for col in range(2):
        for k, key in enumerate(["mouse", "zebrafish", "worm"]):
            c.badge(0.1 + col * col_w + 4.1 + k * 0.56, 0.62, D.sp[key].code, D.sp[key].color, r=0.13)
    for i, g in enumerate(cands):
        col, row = divmod(i, per_col)
        x, y = 0.1 + col * col_w, 0.98 + row * 0.315
        fg = T.NAVY if g.known else T.ACCENT
        c.text(x + 0.3, y, f"{g.rank}", size=10.5, color=T.FAINT, ha="right")
        c.text(x + 0.42, y, g.gene, size=12.5, weight="bold", color=fg)
        w = 2.45 * (-math.log10(g.score)) / top
        c.box(x + 1.45, y - 0.11, w, 0.22, fill=fg, edge=None, r=0.04)
        sx = x + 4.05
        for key in ["mouse", "zebrafish", "worm"]:
            if key in g.support:
                sp = D.sp[key]
                c.text(sx + 0.05, y, f"{g.support[key]}", size=11.5, color=sp.color, weight="bold", ha="center")
            sx += 0.56
    return c


@figure("bbs_decomposition", "Two parts. Left, in the style of the published figure: for Bardet-Biedl syndrome, "
                             "sixteen protein families as stacked bars on an additive score axis, one block per "
                             "top-10 phenolog that carries the family, orange if non-obvious and navy if "
                             "semantically similar; the eleven highest bars are all known "
                             "Bardet-Biedl genes (BBS7, BBS10, BBS1, BBS4 and others). Right: the ten phenologs in two "
                             "lists. Non-obvious, zebrafish embryo: gastrulation disrupted, notochord kinked, somite "
                             "increased width, body axis decreased length, notochord increased width. Semantically "
                             "similar, mouse: retina degeneration, obese, photoreceptor outer and inner segment "
                             "degeneration, absent sperm flagellum.")
def fig_bbs_decomposition(D, step=None):
    c = Canvas()
    b, dec = D.bbs, D.bbs.decomposition
    phs = {p.rank: p for p in dec.phenologs}
    assert sorted(phs) == list(range(1, len(phs) + 1))
    weight = {r: p.shared / b.families * (1 - p.p) for r, p in phs.items()}      # the additive score's term
    colour = {r: T.NAVY if p.obvious else T.ACCENT for r, p in phs.items()}     # two colours: easier to read (round 3)
    order = [p.rank for p in dec.phenologs if not p.obvious] + [p.rank for p in dec.phenologs if p.obvious]
    fams = dec.families
    totals = [sum(weight[r] for r in f.ranks) for f in fams]
    assert all(t0 >= t1 - 1e-9 for t0, t1 in zip(totals, totals[1:])), "bbs.decomposition.families is not in score order"
    y0, ph = 0.05, T.FIG_H - 0.1

    # a) the decomposition, laid out like the preprint's Fig. 5 (round 2): white ground, the family on the left and
    #    its genes at the end of the bar, a score axis; fewer phenologs and families than published. The blocks are
    #    in two colours by kind of phenolog, as before round 2, not one colour per phenolog (round 3).
    xa, wa = 0.0, 6.42
    c.box(xa + 0.01, y0, wa - 0.02, ph, fill=T.WHITE, edge=T.RULE, lw=1.2, r=0.12, zorder=0)
    c.text(xa + 0.22, y0 + 0.3, f"Additive decomposition (top {len(phs)} phenologs, top {len(fams)} protein families)",
           size=12.5, weight="bold")
    lx = xa + 0.22
    for lab, col in (("non-obvious", T.ACCENT), ("semantically similar", T.NAVY)):
        c.box(lx, y0 + 0.49, 0.22, 0.14, fill=col, edge=None, r=0.03)
        c.text(lx + 0.29, y0 + 0.56, lab, size=10.5, color=T.MUTED)
        lx += 0.29 + c.measure(lab, 10.5)[0] + 0.22
    c.text(xa + wa - 0.22, y0 + 0.56, "number: phenolog rank   ·   bold: known gene", size=10, color=T.FAINT, ha="right")

    def names(f):
        """A family's genes as (known, the rest): the known ones are set in bold."""
        known, other = list(f.get("known", [])), list(f.get("other", []))
        rest = f.get("n_other", len(other)) - len(other)
        tail = ", ".join(other) + (f" +{rest}" if rest else "")
        return ", ".join(known), (", " + tail if known and tail else tail)

    gs = 10.5
    id_x, bar_x = xa + 1.1, xa + 1.2
    avail = xa + wa - 0.18 - bar_x
    scale = min((avail - 0.1 - c.measure(k, gs, "bold")[0] - c.measure(t_, gs)[0]) / tot
                for (k, t_), tot in zip(map(names, fams), totals))
    row0, pitch, bar_h, gap = y0 + 0.88, 0.208, 0.158, 0.018
    axis_y = row0 + (len(fams) - 1) * pitch + pitch / 2 + 0.05
    for tick in range(int(avail / scale) + 1):                       # the score axis, with the paper's dashed grid
        x = bar_x + tick * scale
        if tick:
            c.line(x, row0 - pitch / 2, x, axis_y, color=T.RULE, lw=0.9, ls=(0, (3, 2.5)), zorder=1)
        c.line(x, axis_y, x, axis_y + 0.05, color=T.INK, lw=0.9)
        c.text(x, axis_y + 0.17, str(tick), size=9.5, ha="center")
    c.line(bar_x, row0 - pitch / 2, bar_x, axis_y, color=T.INK, lw=0.9, zorder=3)
    c.line(bar_x, axis_y, bar_x + avail, axis_y, color=T.INK, lw=0.9, zorder=3)
    c.text(bar_x + avail / 2, axis_y + 0.38, "Additive prediction score", size=10.5, ha="center")
    for i, f in enumerate(fams):
        y = row0 + i * pitch
        c.text(id_x, y, f.family, size=9, color=T.MUTED, ha="right")
        x = bar_x
        for r in order:
            if r in f.ranks:
                w = weight[r] * scale
                c.rect(x, y - bar_h / 2, w - gap, bar_h, fill=colour[r])
                c.text(x + (w - gap) / 2, y + 0.004, str(r), size=8.5, weight="bold", color=T.WHITE, ha="center")
                x += w
        known, tail = names(f)
        c.text(x + 0.06, y, known, size=gs, weight="bold")
        c.text(x + 0.06 + c.measure(known, gs, "bold")[0], y, tail, size=gs, color=T.MUTED)
    with_known = sum(1 for f in fams if f.get("known"))
    c.text(xa + 0.22, y0 + ph - 0.24, f"{with_known} of the top {len(fams)} families hold a known {b.label} gene.", size=12,
           weight="bold")

    # b) what the ten phenotypes are: the two kinds of counterpart, as on the rod-cone slide
    xb = xa + wa + 0.2
    wb, bgap = c.w - xb - 0.02, 0.14
    bh = (ph - bgap) / 2
    kinds = ((True, "Semantically similar", T.WHITE, T.NAVY, T.INK, T.MUTED, dec.why_similar, ""),
             (False, "Non-obvious", T.ACCENT_TINT, T.ACCENT, T.ACCENT_DARK, T.INK, dec.why_nonobvious, "embryo"))
    lab_w = wb - 0.7 - 0.25 - 0.45
    lab_size = min(c.fit(p.label, lab_w, 12, min_size=9.0) for p in dec.phenologs)
    lines = [c.wrap(kd[6], wb - 0.5, 11).count("\n") + 1 for kd in kinds]
    shift = 0.095 * (lines[1] - lines[0])               # the box with the longer sentence gets the room
    tops, heights = (y0, y0 + bh - shift + bgap), (bh - shift, bh + shift)
    for k, (obvious, head, fill, edge, head_col, why_col, why, what) in enumerate(kinds):
        by, bh = tops[k], heights[k]
        rows = [p for p in dec.phenologs if p.obvious == obvious]
        assert len({p.species for p in rows}) == 1, "the two lists are drawn as one species each"
        sp = D.sp[rows[0].species]
        c.box(xb, by, wb, bh, fill=fill, edge=edge, lw=1.5, r=0.12)
        c.text(xb + 0.25, by + 0.3, head, size=13.5, weight="bold", color=head_col)
        hx = xb + 0.25 + c.measure(head, 13.5, "bold")[0] + 0.3
        c.badge(hx, by + 0.3, sp.code, sp.color, r=0.115)
        c.text(hx + 0.19, by + 0.3, f"{sp.common.lower()} {what}".strip(), size=10.5, color=T.MUTED)
        c.text(xb + wb - 0.25, by + 0.3, f"known families shared, of {b.families}", size=9.5, color=T.FAINT, ha="right")
        for j, p in enumerate(rows):
            y = by + 0.68 + j * 0.265
            c.box(xb + 0.25, y - 0.095, 0.32, 0.19, fill=colour[p.rank], edge=None, r=0.03)
            c.text(xb + 0.41, y + 0.004, str(p.rank), size=9, weight="bold", color=T.WHITE, ha="center")
            c.text(xb + 0.7, y, p.label, size=lab_size)
            if p.get("hpo"):                            # the disease's own HPO term for the same thing
                c.text(xb + 0.7 + c.measure(p.label, lab_size)[0] + 0.16, y, f"↔ HPO: {p.hpo}", size=10,
                       weight="bold", color=T.CAROLINA_DARK)
            c.text(xb + wb - 0.25, y, f"{p.shared}", size=12, weight="bold", ha="right")
        c.text(xb + 0.25, by + 0.68 + (len(rows) - 1) * 0.265 + 0.2, c.wrap(why, wb - 0.5, 11), size=11, color=why_col,
               va="top", linespacing=1.25)
    return c


# =============================================================================================
# 7. Resource
# =============================================================================================
@figure("dashboard", "Screenshot of the phenologs dashboard page for Bardet-Biedl syndrome in a browser frame, beside "
                     "what the release contains.")
def fig_dashboard(D, step=None):
    c = Canvas()
    shot = ROOT / "assets" / "dashboard_bbs.png"
    fx, fy, fw = 0.0, 0.05, 7.6
    bar_h = 0.34
    fh = 4.55
    if shot.exists():
        from PIL import Image
        img = Image.open(shot)
        keep = min(img.height, round(img.width * fh / (fw - 0.06)))
        shot = ROOT / "out" / "dashboard_cropped.png"
        shot.parent.mkdir(exist_ok=True)
        img.crop((0, 0, img.width, keep)).save(shot)
    c.box(fx, fy, fw, bar_h + fh + 0.06, fill="#DDE3EA", edge=None, r=0.1)
    for k, col in enumerate(["#F26B43", "#FBAE40", "#7BC67B"]):
        c.circle(fx + 0.22 + k * 0.2, fy + bar_h / 2, 0.055, fill=col, edge=None)
    c.box(fx + 0.95, fy + 0.07, fw - 1.2, bar_h - 0.14, fill=T.WHITE, edge=None, r=0.08)
    c.text(fx + 1.1, fy + bar_h / 2, D.meta.dashboard_url, size=9.5, color=T.MUTED)
    if shot.exists():
        c.image(shot, fx + 0.03, fy + bar_h, w=fw - 0.06, h=fh)
    else:
        c.box(fx + 0.03, fy + bar_h, fw - 0.06, fh, fill="#0D1117", edge=None, r=0.04)
        c.text(fx + fw / 2, fy + bar_h + fh / 2, "dashboard screenshot\n(run tools/screenshot_dashboard.mjs)", size=14,
               color=T.WHITE, ha="center")

    t = D.totals
    px, pw = 8.0, 4.0
    c.text(px, 0.3, "Release content", size=15, weight="bold")
    yy = 0.8
    for head, body in [(f"{n(t.disease_phenologs)} disease phenologs", "every significant pair and its p-value"),
                       (f"{n(t.ranked_pairs)} ranked gene–disease pairs",
                        f"across {n(t.diseases_with_phenologs)} diseases, known and novel labeled"),
                       (f"{n(t.phenotype_phenologs)} phenotype phenologs", "top 100 candidate genes per HPO term"),
                       ("Parquet tables", "one release per Monarch KG version")]:
        yy = _note(c, px, yy, pw, head, body)
    c.box(px, 4.05, pw, 0.9, fill=T.NAVY, edge=None, r=0.12)
    c.text(px + 0.22, 4.34, D.meta.dashboard_url, size=c.fit(D.meta.dashboard_url, pw - 0.44, 12.5, "bold"), weight="bold",
           color=T.WHITE)
    c.text(px + 0.22, 4.66, D.meta.repo_url, size=c.fit(D.meta.repo_url, pw - 0.44, 11), color=LIGHT)
    return c


# =============================================================================================
# Backup
# =============================================================================================
@figure("loo", "Published leave-one-out heatmaps for naive Bayes and hypergeometric scores, with the headline "
               "comparison beside them.")
def fig_loo(D, step=None):
    c = Canvas()
    L = D.loo
    c.image(paper.crop("figS7ab", CROP_DIR), 0.0, 0.1, w=7.7)
    px, pw = 8.0, 4.0
    _tile(c, px, 0.1, pw, 1.35, n(L.families_hidden), "protein families hidden, one at a time",
          "phenologs and rankings recomputed each time", fill=T.NAVY, fg=T.WHITE, sub_fg=LIGHT, big_size=28)
    _tile(c, px, 1.6, pw, 1.9, f"{L.true_disease_top20_hg} vs {L.true_disease_top20_nb}",
          "diseases where the hidden gene's true disease ranks in its top 20",
          f"hg vs naive Bayes, k = {L.k}. Within a disease the two tie.", big_size=28)
    _tile(c, px, 3.65, pw, 1.35, f"{L.more_precise_pct}%", "of diseases: hg is more precise",
          f"top-20 recall {L.recall_top20_hg}% vs {L.recall_top20_nb}%", big_size=28)
    return c


@figure("additive_bbs", "Published additive decomposition for Bardet-Biedl syndrome: the top 20 protein families as "
                        "stacked bars, each segment one supporting phenolog.")
def fig_additive_bbs(D, step=None):
    c = Canvas()
    b = c.image(paper.crop("fig5", CROP_DIR), 0.0, 0.0, h=5.1)
    px = b.right + 0.5
    pw = 12.0 - px
    yy = 0.5
    c.text(px, 0.2, "Reading the bars", size=15, weight="bold")
    for head, body in [("One bar, one PANTHER family", "Its length is the family's additive score."),
                       ("One segment, one phenolog", "Each of the disease's top 20 phenologs that carries a gene of "
                                                     "the family adds (shared ÷ disease families) × (1 − P)."),
                       ("Bold gene names are known genes", "16 of the top 20 families hold an established "
                                                          "Bardet–Biedl gene."),
                       ("Mouse and zebrafish drive it jointly", "Retinal, renal, sperm and obesity phenotypes from mouse; "
                                                                "gastrulation, notochord and pronephros from zebrafish.")]:
        yy = _note(c, px, yy + 0.2, pw, head, body)
    return c


@figure("carney", "Published additive decomposition for Carney complex: five phenologs from mouse, slime mold and worm; "
                  "only the two causal families are supported by all five.")
def fig_carney(D, step=None):
    c = Canvas()
    b = c.image(paper.crop("fig4d", CROP_DIR), 0.0, 0.05, h=5.0)
    px = b.right + 0.4
    pw = 12.0 - px
    c.text(px, 0.25, "A slime mold breaks the tie", size=15, weight="bold")
    yy = 0.75
    for head, body in [("Carney complex, caused by PRKAR1A", "An endocrine tumour syndrome. Five significant phenologs."),
                       ("Three mouse adiposity phenologs", "They support the causal families, and eight others just "
                                                           "as strongly."),
                       ("Slime mold and worm add two", "Aberrant chemotaxis to cAMP and a quiescence variant back only "
                                                       "the phosphodiesterase and PKA regulatory families."),
                       ("Only those two have all five", "An organism with no endocrine system ranks the gene of an "
                                                        "endocrine syndrome first.")]:
        yy = _note(c, px, yy, pw, head, body, dot=T.ACCENT if "Slime" in head else T.CAROLINA)
    return c
