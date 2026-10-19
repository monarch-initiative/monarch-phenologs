# Phenologs dashboard

A static site over the pipeline's releases: a release at a glance, and every disease, gene, phenotype, protein
family and species as a page you can search for. There is no server side. The browser reads a release's parquet
tables with DuckDB-WASM. Releases are made by the pipeline's step `python/07_build_release.py`.

```
uv run python python/07_build_release.py -p path/to/project_dir     # ~20 s
python dashboard/serve.py -p path/to/project_dir                     # http://localhost:8765
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
  serve.py           local server: the dashboard, with a run's release/ at /releases/ and a releases.json listing it
```

## The release

`07_build_release.py` reads the project directory a run wrote and writes the release next to it, laid out as in
the bucket. It reads `monarch_kg/versions.json`, `species_data/`, the phenologs, the FDR table and the candidate
ranking.

```
<project_dir>/release/
  2026-03-06/
    disease/             the release, the folder that goes to the bucket (~30 MB)
      overview.json      everything the release page shows; its release block has the versions and settings
      README.md          what every table and column holds
      manifest.json      every file with its size, rows, sha256 and what it holds
      *.parquet          the tables below, which are also the downloads
  releases.json          the published list with this run's releases in it, written by 08_update_release_list.py
```

- **Settings.** The release holds one ranking: the FDR, k and rank metric given (defaults 0.05, 1000, hg).
- **Versions.** The KG and PHENIO versions come from `monarch_kg/versions.json`, which step 01 writes. The code
  version comes from the image, or from git.
- **Rebuilds.** A rebuild replaces the release folder.
- **Other options.** `-c` and `-m` cap DuckDB's cores and memory, for a container.

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

- **Counts.** Counts are 32 bit, so they arrive in the browser as plain numbers.
- **Downloads.** The tables are also the downloads. `README.md` in each release describes every column (a test
  checks that none is missing), and `duckdb` can read the tables straight from a bucket URL.

## Adding a page or a chart

A page module exports `render(el, ...params)` and is wired in `ROUTES` in `app.js`. It asks `db.js` for the tables
it needs with `ensure(...)`, queries with `query(sql)`, and builds HTML with the helpers in `common.js` and the
charts. Anything that needs a new column goes in the table SQL of `phenologs_release_utils.py` first, and in the
release README there.

## Hosting

The dashboard and the releases are hosted apart. This folder goes on GitHub Pages. The releases go in the Monarch
data bucket, `gs://data-public-monarchinitiative/phenologs/alpha/`, which is public at
`https://data.monarchinitiative.org/phenologs/alpha/`. CORS is open and range requests work, but the bucket can't be
listed, so `releases.json` lists the releases.
- **Which release.** The site reads `releases.json` there and opens the newest disease release, or
  `?release=<graph_version>` and `?network=phenotype`. It asks for the list past any cache, so a new release shows
  at once.
- **Local data.** Served by `serve.py`, the site reads the run's `release/` folder instead, and `serve.py` lists the
  releases found there. `?data=<url>` points the site at any other releases folder.
- **Publishing.** See `RELEASING.md`.
- **Layout only grows.** A release folder never changes once uploaded. Later builds may add tables and columns, but
  never rename or remove one. Each release block's `schema` (`SCHEMA` in `phenologs_release_utils.py`) tells the site
  what that release has.

**Fetching.** The parquet files are fetched whole on first use, so range support is not required; the largest is 17 MB. `serve.py` has range support anyway, for when the data layer moves to range reads.

**Outside sources.** Scripts come from jsdelivr (DuckDB-WASM) and cdnjs (D3), and fonts from Google Fonts.

## Later

- `links.dismech()` in `format.js` returns null for now. Once DisMech's disorder ids are known it can return
  a URL and the disease page will show the link.
- Leave-one-out validation and the temporal comparison between releases are not in the build yet.
