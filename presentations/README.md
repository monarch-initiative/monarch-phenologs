# Presentations

Talks about this project, built from text files so they can be remade in one command.

## The phenologs talk (20 to 25 minutes)

```
make setup      # once: a small environment for slide tooling, separate from the pipeline's uv.lock
make deck       # figures + deck  ->  out/phenologs_talk.pptx
make pdf        # also a PDF and one PNG per slide, exported through PowerPoint (macOS)
```

`out/` is not tracked; everything in it is rebuilt.

### Two files to add first

The UNC template and the preprint PDF are not kept in git. Put them in `assets/` before building:

- `assets/UNC_NEW_Template_2026.pptx`: the UNC PowerPoint template.
- `assets/phenologs_preprint_biorxiv.pdf`: the preprint PDF from bioRxiv (its doi is `meta.doi` in `data.yaml`), the
  version the crop boxes in `paper.py` were set on.

The build stops and names the file if either is missing.

### Where things live

| file | holds | edit it to |
|---|---|---|
| `phenologs_talk/deck.yaml` | the slides, in order: titles, one-sentence ledes, speaker notes | reword, reorder, drop (`skip: true`) or add slides |
| `phenologs_talk/data.yaml` | every number on the slides, with its source | move to a new release |
| `phenologs_talk/figures.py`, `figures_results.py` | one function per diagram | change a picture |
| `phenologs_talk/paper.py` | crop boxes for panels cut from the preprint PDF | re-cut after a preprint revision |
| `phenologs_talk/theme.py`, `draw.py` | colours, slide geometry, drawing helpers | restyle everything at once |
| `phenologs_talk/build.py` | puts it on the UNC template | add a slide kind |
| `assets/` | the dashboard screenshot, plus the UNC template and the preprint PDF (those two are not in git) | swap a template or source |

Text in `deck.yaml` can quote `data.yaml`: `{totals.disease_phenologs:,}` prints 61,243 and
`{bbs.best_p!p}` prints 9.6 × 10⁻²⁸. No number is typed into a slide or a figure.

After a wording change only, `make text` rebuilds the deck in a second without redrawing figures.
A title may break onto a second line with `\n`; such a slide has no room for a lede.

### How the slides are made

- **Template.** `assets/UNC_NEW_Template_2026.pptx`. The build keeps its master and layouts and drops its
  example slides. Titles, footers and slide numbers are the template's own placeholders.
- **Text** (titles, ledes, cards, notes) is real text, so it stays editable.
- **Diagrams** are drawn with matplotlib at the size they appear on the slide (a 14 pt label is 14 pt on
  the slide) and placed as pictures. A figure can be drawn in steps, one slide each
  (`phenolog_cartoon`, on slides that are parked for now).
- **Result plots** (FDR curves, temporal validation, leave-one-out, additive decompositions) are cut from
  the preprint PDF at 300 dpi, so they are the published figures. Needs `pdftoppm` (`brew install poppler`).
- **Dashboard screenshot.** `make screenshot` drives headless Chrome and waits for the page's data
  before capturing (Node 22+).

### Google Slides

Upload `out/phenologs_talk.pptx` to Drive and open it with Google Slides (or File > Import slides).
The deck uses Arial only, pictures for diagrams and plain shapes for cards, all of which import as they are.
Speaker notes come across. Edits made in Google Slides are not carried back here: change `deck.yaml` and
rebuild, or treat the Google copy as the final hand-edited version.

### Editing the deck in PowerPoint

Reordering or rewording slides in `out/phenologs_talk.pptx` is fine, but the next build writes that file again.
The build notices a deck that changed since it last wrote it and keeps a copy in `hand_edited/` (not tracked)
first. To make hand edits stick, carry them into `deck.yaml`: the slide order is the order of the list, and each
slide's words are under its `id`. A slide deleted in PowerPoint is parked with `skip: true`, not removed from the file.
Layout changed by hand has its own keys on the slide: `title_size`, `title_box`, `lede_size`, `lede_box` and
`figure_dy` (described at the top of `deck.yaml`).

PowerPoint shrinks text that does not fit its placeholder when it opens a deck (the cover's subtitle to 70%, the
footer's line spacing by 10%). The build writes the same autofit values, so the built file and a copy saved from
PowerPoint look alike.

### Moving to a new data release

1. Update `data.yaml`. The species block, the totals and the Bardet-Biedl and rod-cone dystrophy examples
   were read from the dashboard (Release, Species, and the two term pages); the evaluation block and the
   `comparisons` block (Table S4) from the preprint. Each block names its source.
2. If the preprint changed, replace `assets/phenologs_preprint_biorxiv.pdf` and run `make crops` to check
   the boxes in `paper.py`.
3. `make screenshot`, then `make deck`.

### Open items in the current deck

- The problem slide shows 50 to 70% undiagnosed, the complement of the 30 to 50% diagnostic yield in the outline;
  that yield still needs a citation.
- `future` ("Can we leverage this large dataset in meaningful ways?") holds four placeholder directions; the ncRNA
  figure is not in yet.
- Which phenologs are tagged "non-obvious" is a judgement call recorded as `clinic:` / `obvious:` in `data.yaml`.
- The order, titles and ledes follow the deck as edited in PowerPoint on 2026-11-16 (kept in `hand_edited/`): the
  problem, Bardet-Biedl, the model organism data that could fill the gap (`gaps`), the three papers (`mcgary`), the
  idea as a cartoon (`idea`); Part 1 (`kg_components`, `then_now`, `concept`, `rodcone_tree`, `counts`); Part 2
  (`pooling`, `naive_bayes`, `hypergeom_marbles`, `time_split`, `eval_results`, `bbs_model`, `resource`); Part 3
  (`future`), `thanks`. 35 slides, 23 before the backup divider.
- The lede on `bbs_model` was typed with two placeholders ("ranks below XYZ except for one at rank XYZ"). They are
  filled from `bbs.known_ranks` in `data.yaml`: 25 of the 26 known genes rank 233 or better, SCAPER is at 1,228.
- The footer names a slide's section. It follows where the slide now sits: `gaps` says "The problem", and
  `recover_known` and `precision` say "Backup".
- `recover_known` (figure `recover_known`, in backup since 2026-11-16) has two halves that are two different analyses in the preprint
  (Supplementary Note 7). Left: leave-one-out, one protein family hidden at a time, 2,001 families. Right: the
  precision-recall curve per disease, which the preprint computes on the complete model with no genes withheld,
  for the 414 diseases with at least four known genes. If the precision analysis was in fact run on the
  leave-one-out results, the right half, the slide's notes and the caption on `precision` need changing.
  The list and curve in the right half are a made-up example (`loo.pr_example` in `data.yaml`).
- `precision` shows Fig. S7C and D as published (crop `figS7cd`); `figS7c` and `figS7d` are also defined in
  `paper.py` if only one panel is wanted. The leave-one-out heatmaps (Fig. S7A and B) are the backup slide `loo`.
- `bbs_model` (figure `bbs_decomposition`) is the slide on how well the method models one disease, and why. Left:
  the dashboard's additive decomposition for Bardet-Biedl syndrome at "top phenologs = 10", 16 protein families as
  stacked bars. Right: the ten phenologs in two lists, non-obvious (zebrafish embryo) and semantically similar
  (mouse), the same split as on `rodcone_tree`. Everything drawn is `bbs.decomposition` in `data.yaml`, read from the
  release's `phenologs` and `additive` tables through the dashboard's database; the block widths and bar lengths are
  computed from it. Which phenotypes are non-obvious (`obvious:`) and the two sentences in the boxes are judgement
  calls, and reading the five fish phenotypes as convergent extension is an interpretation (see the slide's notes).
  The gene-by-gene top 25 it replaced (`bbs_ranking`) is in backup.
- Backup holds what left the main flow: `recover_known`, `precision`, `loo`, `bbs_ranking`, `additive`, `carney`, `top_species`,
  `species`, `networks`, and a second copy of the `gaps` slide (`gaps_backup`), as in the hand-edited deck.
- On `idea`: the network sizes follow the 2010 association counts on a log scale; "≈ 110×" is the four model
  organisms' 210,619 associations over the 1,923 human ones; the Venn (1 | 2 | 3) and the candidate genes are from
  the paper's Fig. 4B; the embryo drawings are schematic, not the paper's photographs.
- On `kg_components` the model organism table counts every annotated term (11,009 for mouse), while `concept` shows
  the terms testable against human (10,998); the human disease count is the dashboard's 8,040, where the preprint's
  Fig. 1A says 7,935.
- Parked in `deck.yaml` with `skip: true`, at the end of the file. Taken out of the deck by hand on 2026-11-16:
  the two-genes follow-up to the score slide (`hypergeom_two_genes`), the take-home numbers (`takehome`), the
  dot-grid score explainer (`hypergeom`), how significance is set (`fdr`) and two earlier section dividers
  (`part1_old`, `part2_old`). Earlier versions of slides: the earlier 2010 slide with its data tiles and
  non-obvious models (`mcgary_models`), the two-network cartoon in three builds (`classical`, `classical_test`),
  gene equivalence (`equivalence`), the breast-ovarian cancer and rod-cone "nobody would guess" slides
  (`nonobvious`, `nonobvious_rodcone`), the timeline (`history`), the Bardet-Biedl walk-through and tree
  (`derivation`, `tree`) and the networks with orthology arcs (`dataset`).
- On `pooling`, the third panel is the output table of the preprint's Fig. 1C with its score column left as question
  marks on purpose (the next two slides compute it); the ranks are the hypergeometric ones at k = 1,000. The three
  scoring methods sit in the strip under the panels.
- On `rodcone_tree`, the ten phenologs shown are a pick from 688, and the cilium explanation of the non-obvious
  ones is an interpretation, see its notes.
- "Then and now" compares with McGary et al. 2010 only: Woods et al. 2013 does not state its totals. The 2010 paper
  had no human phenotype annotations, so those two rows have no "then" bar.
- The significance slide (`fdr`) is parked, so the deck states the 5% FDR without showing how it is set.
- The 2010 paper's facts (`history.mcgary2010` in `data.yaml`) were read from the main text on PMC, not the supplement.
- The score explainer is `hypergeom_marbles` (figure `marbles`): a bag of marbles in three panels, worked through
  for one real gene (ARL2BP) and Bardet-Biedl syndrome. Panel 2 is the phenolog table of the preprint's Fig. 1C with
  a column for the marbles each phenolog draws, pooling into 6,436; panel 3 is the four numbers and the score. The
  numbers are `bbs.score_example` in `data.yaml`: the bag (N = 401,464 annotations), the marbles drawn (n = 6,436,
  every gene on the 109 phenologs, spread across the 3,075 candidates), the gene's own (m = 36) and how many of
  those were drawn (c = 14), plus the per-phenolog counts in `phenolog_rows`. They were read from the release tables
  through the dashboard's database, and they reproduce the published score for seven genes, so the build
  recomputes the score from them and checks it. The bag and the tray are schematic. On a new release these numbers
  have to be read again.
- `hypergeom_two_genes` (figure `marbles_two_genes`, parked) is a follow-up to it: a specific and a hub gene in one
  bag, drawn to made-up counts (`hg_illustration.marbles`, 120 marbles, 20 drawn) and labelled as such. The earlier
  dot-grid explainer (`hypergeom`) is parked too.
- The cover's date line is `meta.date` in `deck.yaml`.
