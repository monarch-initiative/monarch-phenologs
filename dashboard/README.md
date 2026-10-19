# Phenologs dashboard

A static site over the pipeline's runs: a release at a glance, and every disease, gene, phenotype,
protein family and species as a page you can search for. There is no server side. A build step turns a
run's TSVs into a release folder of parquet tables, and the browser queries those with DuckDB-WASM.

```
python dashboard/build_release.py -p path/to/project_dir --phenio_version 2025-11-14   # ~20 s, ~150 MB
python dashboard/serve.py                                                             # http://localhost:8765
```

`build_release.py` finds everything under the project directory the pipeline wrote (`species_data/`,
the phenologs, the FDR table, the candidate rankings). The site serves the FDR, k and rank metric you ask
for (defaults 0.05, 1000, hg). The graph version is taken from a date in the project path unless given; with
neither, the build stops.
The release is written to `dashboard/releases/<graph_version>/<prediction_network>/` (`-o` for another
releases folder) and added to `releases.json` there. A rebuild replaces the folder. `--threads` and
`--memory_limit` cap DuckDB, for a container.

## Layout

```
dashboard/
  build_release.py       CLI: project dir -> releases/<graph_version>/<prediction_network>/
  release/
    sources.py           where a run keeps things; DuckDB views over its files; species display details and the tree
    tables.py            one SQL statement per entity table (phenologs, candidates, genes, phenotypes, families, diseases, species, search)
    additive.py          the Woods 2013 additive decomposition, per disease (what the case study figure shows)
    overview.py          overview.json: tiles, species stats, FDR table, tree, top lists, graph inputs, known pair funnel
    downloads.py         the run's raw results as single tables, and a README on their columns
    manifest.py          manifest.json (every file of a release) and releases.json (every release)
  serve.py               local static server: site/ at /, releases/ at /releases/, with HTTP range support
  site/                  HTML, JS and CSS only
    index.html, style.css
    js/app.js            hash router, search box
    js/db.js             DuckDB-WASM: opens a release, fetches its tables on first use
    js/format.js         formatting, species colours, outbound links (the DisMech hook lives here)
    js/charts/           radar, species tree, additive stacked bars, bar lists, tables (in memory, or paged by query for the 1.5M candidates)
    js/views/            one module per page: overview, candidates, disease, gene, phenotype, family, species, browse, search
  releases/              built releases (ignored by git), laid out as in the bucket
    releases.json        every release's release block and folder, newest first
    <graph_version>/
      disease/           overview.json, manifest.json, the parquet tables below, downloads/
      phenotype/         (later) the same for the phenotype network
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

`downloads/` holds the run's raw results, for people who want everything. Its `README.md` describes every column.

| download | one row per | notes |
|---|---|---|
| `candidate_rankings.parquet` | (ranking, disease, gene) | every ranking the run made at its FDR (k x rank metric), as the pipeline wrote it, known genes with no score included. 7.6M rows, ~107 MB |
| `phenologs_all_tested.parquet` | (disease, species, phenotype) | every pair tested, with `significant` at the run's FDR. 2.9M rows, ~13 MB |
| `phenologs_significant.tsv.gz` | significant phenolog | the `phenologs` table as TSV, for spreadsheets |
| `fdr_table.tsv` | (species, FDR level) | the p-value cutoffs |

## Adding a page or a chart

A view exports `render(el, ...params)` and is wired in `ROUTES` in `app.js`. It asks `db.js` for the tables it
needs with `ensure(...)`, queries with `query(sql)`, and builds HTML with the helpers in `views/common.js`
and the charts in `js/charts/`. Anything that needs a new column goes in `tables.py` first.

## Hosting

The site and the releases are hosted apart: `site/` on GitHub Pages, the releases in the Monarch data bucket,
`gs://data-public-monarchinitiative/phenologs/alpha/`, which is public at
`https://data.monarchinitiative.org/phenologs/alpha/` (CORS open, range requests work, no listing). The site
reads `releases.json` there, opens the newest disease release (or `?release=<graph_version>`,
`?network=phenotype`), and fetches that release's files. Served locally it reads `releases/` instead, and
`?data=<url>` points it at any other releases folder. A release folder never changes once uploaded. When a
table's layout changes, `SCHEMA` in `overview.py` goes up and the older releases are rebuilt.

Publishing a release. The bucket's `releases.json` is fetched first, so the build adds the new release to the
published list rather than replacing it (skip the fetch for the very first release):

```
curl -fsS https://data.monarchinitiative.org/phenologs/alpha/releases.json -o dashboard/releases/releases.json
python dashboard/build_release.py -p path/to/project_dir --graph_version 2026-03-06 --phenio_version 2025-11-14
gsutil -m rsync -r dashboard/releases/2026-03-06/disease gs://data-public-monarchinitiative/phenologs/alpha/2026-03-06/disease
gsutil -h "Cache-Control:no-cache" cp dashboard/releases/releases.json gs://data-public-monarchinitiative/phenologs/alpha/releases.json
```

Upload the release folder before `releases.json`, so the site never lists a release that is not there yet.
`no-cache` makes a new release show up at once; the bucket otherwise caches files for an hour.

The parquet files are fetched whole on first use (the largest is 17 MB), so range support is not required;
`serve.py` has it anyway for when the data layer moves to range reads. Scripts come from jsdelivr
(DuckDB-WASM) and cdnjs (D3); fonts from Google Fonts.

## Later

- `links.dismech()` in `format.js` returns null for now. Once DisMech's disorder ids are known it can return
  a URL and the disease page will show the link.
- Leave-one-out validation and the temporal comparison between releases are not in the build yet.
