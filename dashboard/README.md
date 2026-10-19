# Phenologs dashboard

A static site over the pipeline's releases: a release at a glance, and every disease, gene, phenotype, protein
family and species as a page you can search for. There is no server side. The browser reads a release's parquet
tables with DuckDB-WASM. Releases are made by the pipeline's last step, `python/07_build_release.py`.

```
uv run python python/07_build_release.py -p path/to/project_dir -gv 2026-03-06 -pv 2025-11-14   # ~20 s
python dashboard/serve.py -p path/to/project_dir                                                  # http://localhost:8765
```

## Layout

```
dashboard/
  index.html, style.css
  js/
    app.js           hash router, search box
    db.js            DuckDB-WASM: opens a release, fetches its tables on first use
    format.js        formatting, species colours, outbound links (the DisMech hook lives here)
    common.js        page pieces the pages share
    overview.js  candidates.js  disease.js  gene.js  phenotype.js  family.js  species.js  browse.js  search.js
                     one module per page
    bars.js  radar.js  stacked.js  table.js  tooltip.js  tree.js
                     charts (tables are in memory, or paged by query for the 1.5M candidates)
  serve.py           local server: the dashboard, with a run's release/ at /releases/
```

## The release

`07_build_release.py` reads the project directory a run wrote (`species_data/`, the phenologs, the FDR table,
the candidate rankings) and writes the release next to them, laid out as in the bucket:

```
<project_dir>/release/
  releases.json          every release's release block and folder, newest first
  2026-03-06/
    disease/             the release, the folder that goes to the bucket
      overview.json      everything the release page shows
      manifest.json      every file with its size, rows, sha256 and what it holds
      *.parquet          the tables below
      downloads/         the run's raw results as single tables
```

- **Settings.** The site shows the FDR, k and rank metric given (defaults 0.05, 1000, hg); the downloads hold every ranking the run made.
- **Graph version.** Taken from a date in the project path unless `-gv` is given; with neither, the build stops.
- **Rebuilds.** A rebuild replaces the release folder.
- **Other options.** `-o` writes the releases somewhere else. `-c` and `-m` cap DuckDB's cores and memory, for a container.

The code is `python/phenologs_release_utils.py`. A change in the pipeline's files is a change to the `Run` views there.

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

`downloads/` holds the run's raw results, for people who want everything. Their `README.md` describes every column.

| download | one row per | notes |
|---|---|---|
| `candidate_rankings.parquet` | (ranking, disease, gene) | every ranking the run made at its FDR (k x rank metric), as the pipeline wrote it, known genes with no score included. 7.6M rows, ~107 MB |
| `phenologs_all_tested.parquet` | (disease, species, phenotype) | every pair tested, with `significant` at the run's FDR. 2.9M rows, ~13 MB |
| `phenologs_significant.tsv.gz` | significant phenolog | the `phenologs` table as TSV, for spreadsheets |
| `fdr_table.tsv` | (species, FDR level) | the p-value cutoffs |

## Adding a page or a chart

A page module exports `render(el, ...params)` and is wired in `ROUTES` in `app.js`. It asks `db.js` for the tables
it needs with `ensure(...)`, queries with `query(sql)`, and builds HTML with the helpers in `common.js` and the
charts. Anything that needs a new column goes in the table SQL of `phenologs_release_utils.py` first.

## Hosting

The dashboard and the releases are hosted apart: this folder on GitHub Pages, the releases in the Monarch data
bucket, `gs://data-public-monarchinitiative/phenologs/alpha/`. That bucket is public at
`https://data.monarchinitiative.org/phenologs/alpha/`: CORS is open and range requests work, but it can't be listed.
- **Which release.** The site reads `releases.json` there and opens the newest disease release, or `?release=<graph_version>` and `?network=phenotype`.
- **Local data.** Served by `serve.py`, it reads the run's `release/` folder instead. `?data=<url>` points it at any other releases folder.
- **Release folders never change once uploaded.** When a table's layout changes, `SCHEMA` in `phenologs_release_utils.py` goes up and the older releases are rebuilt.

**Publishing a release.** Fetch the bucket's `releases.json` first, so the build adds the new release to the published
list rather than replacing it. Skip the fetch for the very first release.

```
curl -fsS https://data.monarchinitiative.org/phenologs/alpha/releases.json -o path/to/project_dir/release/releases.json
uv run python python/07_build_release.py -p path/to/project_dir -gv 2026-03-06 -pv 2025-11-14
gsutil -m rsync -r path/to/project_dir/release/2026-03-06/disease gs://data-public-monarchinitiative/phenologs/alpha/2026-03-06/disease
gsutil -h "Cache-Control:no-cache" cp path/to/project_dir/release/releases.json gs://data-public-monarchinitiative/phenologs/alpha/releases.json
```

- **Order.** Upload the release folder before `releases.json`, so the site never lists a release that is not there yet.
- **Caching.** `no-cache` makes a new release show up at once; the bucket otherwise caches files for an hour.

**Fetching.** The parquet files are fetched whole on first use, so range support is not required; the largest is 17 MB. `serve.py` has range support anyway, for when the data layer moves to range reads.

**Outside sources.** Scripts come from jsdelivr (DuckDB-WASM) and cdnjs (D3), and fonts from Google Fonts.

## Later

- `links.dismech()` in `format.js` returns null for now. Once DisMech's disorder ids are known it can return
  a URL and the disease page will show the link.
- Leave-one-out validation and the temporal comparison between releases are not in the build yet.
