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

- The problem slide shows 50 to 70% undiagnosed, the complement of the 30 to 50% diagnostic yield in the outline.
  It now cites Chung et al. 2023 (Genet Med, PMID 37191093) and Clark et al. 2018 (npj Genom Med, PMID 30002876);
  their pooled yields (34 to 41%) put the undiagnosed share at about 59 to 66%, inside the range shown.
- `future` ("Can we leverage this dataset in meaningful ways?") holds four directions as worded in round 1; the
  ncRNA figure is not in yet.
- Round 1 of edits (Notion, 2026-11-16) is applied. The old slides 3 and 4 (`bbs`, `gaps`) are one slide,
  `bbs_gaps` (figure `bbs_fill_gaps`). `time_split` is a simpler drawing (`time_split_simple`) that only mentions
  the 73 evaluable pairs; the version with all 524 pairs is the first backup slide (`time_split_full`).
  The earlier one-slide layout of `eval_results` is the figure `eval_results_tiles`, parked as `eval_results_old`.
  The FAAP100 box is no longer on the slide; it is in the notes.
- Round 3 (Notion, 2026-11-17): the wording of two cards on the problem slide, and the colours on `bbs_model`.
  "Complex diseases" on the problem slide follows the preprint, which also calls Bardet-Biedl syndrome a complex
  disease (`bbs_gaps` does the same); strictly it is Mendelian and genetically heterogeneous, see that slide's notes.
- Round 2 of edits (Notion, 2026-11-17) is applied. Five slides are shown in steps, each slide adding to the last
  (`steps:` in `deck.yaml`, `BUILDS` in `render_figures.py`): `idea` (networks, then the worked example),
  `kg_components` (one panel per step), `concept` (panel a, then b), `hypergeom_marbles` (one panel per step) and
  `eval_results` (left plot, then right plot). `time_split` no longer has the ranked-list example in the middle.
  `eval_results` has no lede: each panel's header carries its headline. Its left panel is titled "within each of
  the 43 diseases" where the edit list said 73: the 73 new pairs fall on 43 diseases. The note about leave-one-out
  cross validation is off the slide. `fdr_plot` (preprint Fig. S3) is a new backup slide.
- On `bbs_gaps`, "3.5× more model organism annotations" is the non-human gene-phenotype edges over the human
  gene-phenotype edges (the first two bars), computed in the figure.
- On `bbs_model`, "↔ HPO" marks a mouse phenotype whose meaning matches a term annotated to the disease in HPO
  (obese and Obesity, retina degeneration and Retinal dystrophy). The match is by meaning against the disease's
  95 HPO terms from the Monarch API, not a formal MP-to-HP mapping.
- On `naive_bayes`, the lede says the second factor is the fraction of the disease's known gene families the
  phenolog shares, which is what `python/phenologs_utils.py` computes; it is not a per-gene probability.
- On `counts`, the McGary et al. 2010 row (3,966) is the sum of that paper's phenologs with human disease over its
  four model organisms, from `history.mcgary2010` in `data.yaml`.
- Which phenologs are tagged "non-obvious" is a judgement call recorded as `clinic:` / `obvious:` in `data.yaml`.
- The order, titles and ledes follow the deck as edited in PowerPoint on 2026-11-16 (kept in `hand_edited/`): the
  problem, Bardet-Biedl with the model organism data that could fill the gap (`bbs_gaps`), the three papers (`mcgary`), the
  idea as a cartoon (`idea`); Part 1 (`kg_components`, `then_now`, `concept`, `rodcone_tree`, `counts`); Part 2
  (`pooling`, `naive_bayes`, `hypergeom_marbles`, `time_split`, `eval_results`, `bbs_model`, `resource`); Part 3
  (`future`), `thanks`. 43 slides, 29 before the backup divider.
- The lede on `bbs_model` was typed with two placeholders ("ranks below XYZ except for one at rank XYZ"). They are
  filled from `bbs.known_ranks` in `data.yaml`: 25 of the 26 known genes rank 233 or better, SCAPER is at 1,228.
- The footer names a slide's section. It follows where the slide now sits: `bbs_gaps` says "The problem", and
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
  stacked bars, laid out since round 2 like the preprint's Fig. 5 (white ground, family ids on the left, genes at
  the bar ends, score axis). The blocks are in two colours by kind of phenolog; round 2's one colour per phenolog
  was dropped in round 3 as harder to read. Right: the ten phenologs in two lists, semantically similar (mouse) on top and non-obvious
  (zebrafish embryo) below, the same split as on `rodcone_tree`. Everything drawn is `bbs.decomposition` in `data.yaml`, read from the
  release's `phenologs` and `additive` tables through the dashboard's database; the block widths and bar lengths are
  computed from it. Which phenotypes are non-obvious (`obvious:`) and the two sentences in the boxes are judgement
  calls, and reading the five fish phenotypes as convergent extension is an interpretation (see the slide's notes).
  The gene-by-gene top 25 it replaced (`bbs_ranking`) is in backup.
- Backup holds what left the main flow: `fdr_plot`, `time_split_full`, `recover_known`, `precision`, `loo`, `bbs_ranking`, `additive`,
  `carney`, `top_species`, `species`, `networks`, and the old slide 4 (`gaps_backup`), as in the hand-edited deck.
- On `idea`: the network sizes follow the 2010 association counts on a log scale; "≈ 110×" is the four model
  organisms' 210,619 associations over the 1,923 human ones; the Venn (1 | 2 | 3) and the candidate genes are from
  the paper's Fig. 4B; the embryo drawings are schematic, not the paper's photographs.
- On `kg_components` the model organism table counts every annotated term (11,009 for mouse), while `concept` shows
  the terms testable against human (10,998); the human disease count is the dashboard's 8,040, where the preprint's
  Fig. 1A says 7,935.
- Parked in `deck.yaml` with `skip: true`, at the end of the file. Replaced in round 1: `bbs`, `gaps`,
  `eval_results_old`. Taken out of the deck by hand on 2026-11-16:
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
- The significance slide with the pipeline steps (`fdr`) is parked; the backup slide `fdr_plot` shows how the 5% FDR
  cutoff is set, with the preprint's Fig. S3.
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
