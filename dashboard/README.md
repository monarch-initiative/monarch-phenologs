# Phenologs dashboard

A static site over one pipeline run: the release at a glance, and every disease, gene, phenotype,
protein family and species as a page you can search for. There is no server side. A build step turns the
pipeline's TSVs into parquet tables, and the browser queries those with DuckDB-WASM.

```
python dashboard/build_release.py -p path/to/project_dir --phenio_version 2025-11-14   # ~10 s, ~30 MB
python dashboard/serve.py                                                             # http://localhost:8765
```

`build_release.py` finds everything under the project directory the pipeline wrote (`species_data/`,
the pooled phenologs, the FDR table, the candidate rankings) for the FDR, k and rank metric you ask
for (defaults 0.05, 1000, hg). The graph version is taken from a date in the project path unless given.

## Layout

```
dashboard/
  build_release.py       CLI: project dir -> site/data/*.parquet + overview.json
  release/
    sources.py           where a run keeps things; DuckDB views over its files; species display details and the tree
    tables.py            one SQL statement per entity table (phenologs, candidates, genes, phenotypes, families, diseases, species, search)
    additive.py          the Woods 2013 additive decomposition, per disease (what the case study figure shows)
    overview.py          overview.json: tiles, species stats, FDR table, tree, top lists, graph inputs, known pair funnel
  serve.py               local static server with HTTP range support
  site/
    index.html, style.css
    js/app.js            hash router, search box
    js/db.js             DuckDB-WASM: tables are fetched on first use
    js/format.js         formatting, species colours, outbound links (the DisMech hook lives here)
    js/charts/           radar, species tree, additive stacked bars, bar lists, tables (in memory, or paged by query for the 1.5M candidates)
    js/views/            one module per page: overview, candidates, disease, gene, phenotype, family, species, browse, search
    data/                built output (ignored by git)
```

## Data tables

| table | one row per | notes |
|---|---|---|
| `phenologs` | significant phenolog | ranked within its disease by p |
| `candidates` | (disease, gene) the ranking scored | rank within the disease and global rank across all diseases, ties sharing a place; score, known flag, and support: per species, the disease's phenologs that are phenotypes of an ortholog of the gene, of how many, through how many orthologs, and the (phenolog, ortholog) hits the score counts (the pipeline's `XSpecies_occurence`). The pipeline also appends a disease's known genes its phenologs never reached, without a p-value; those are left out |
| `human_orthologs` | (model gene, human gene) | ortholog edges from the run's species to the base species |
| `additive` | (disease, family, phenolog) | top 20 phenologs x top 40 families per disease; the page re-slices by k, species and order |
| `diseases`, `genes`, `phenotypes`, `families`, `species` | entity | counts the pages and browse lists show. `genes.in_reach`: a base species gene one of whose orthologs carries a phenotype the ranking can use |
| `known_associations` | curated gene-disease pair | from the base species' gene_to_disease table |
| `gene_phenotypes` | gene-phenotype annotation | every species |
| `search` | searchable thing | kind, id, label, detail, weight |

Counts are 32 bit so they arrive in the browser as plain numbers. The tables are also a reasonable
distribution format on their own: `duckdb` can read them from a bucket URL.

## Adding a page or a chart

A view exports `render(el, ...params)` and is wired in `ROUTES` in `app.js`. It asks `db.js` for the tables it
needs with `ensure(...)`, queries with `query(sql)`, and builds HTML with the helpers in `views/common.js`
and the charts in `js/charts/`. Anything that needs a new column goes in `tables.py` first.

## Hosting

Any static host works (a bucket, GitHub Pages). The parquet files are fetched whole on first use (the
largest is 17 MB), so range support is not required; `serve.py` has it anyway for when the data layer
moves to range reads. Scripts come from jsdelivr (DuckDB-WASM) and cdnjs (D3); fonts from Google Fonts.

## Later

- `links.dismech()` in `format.js` returns null for now. Once DisMech's disorder ids are known it can return
  a URL and the disease page will show the link.
- Leave-one-out validation and the temporal comparison between releases are not in the build yet.
