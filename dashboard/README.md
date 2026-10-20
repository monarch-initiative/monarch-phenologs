# Phenologs dashboard

A static site over the pipeline's releases. A Monarch KG version has two: disease phenologs (human diseases against
the model species' phenotypes) and phenotype phenologs (human phenotypes against them). The landing page shows both
side by side; each has its own release page, and every disease or human phenotype, gene, phenotype, protein family
and species is a page you can search for. There is no server side. The browser reads a release's parquet tables with
DuckDB-WASM. Releases are made by the pipeline's step `python/07_build_release.py`.

```
uv run python python/07_build_release.py -p path/to/project_dir -prd disease      # ~10 s
uv run python python/07_build_release.py -p path/to/project_dir -prd phenotype    # a few minutes
python dashboard/serve.py -p path/to/project_dir                                  # http://localhost:8765
```

## Layout

```
dashboard/
  index.html, style.css
  js/
    app.js           hash router, top bar (a row of pages per network), search box
    db.js            DuckDB-WASM: opens a release, fetches its tables on first use
    format.js        formatting, the open network's words, species colours, outbound links (the DisMech hook lives here)
    common.js        page pieces the pages share
    landing.js  overview.js  candidates.js  term.js  gene.js  phenotype.js  family.js  species.js  browse.js
    search.js  downloads.js
                     one module per page (term.js is the disease or human phenotype page)
    bars.js  radar.js  stacked.js  table.js  tooltip.js  tree.js
                     charts (tables are in memory, or paged by query for the 1.5M candidates)
  serve.py           local server: the dashboard, with a run's release/ at /releases/ and a releases.json listing it
```

## The release

`07_build_release.py -prd <network>` reads the project directory a run wrote and writes that network's release next
to it, laid out as in the bucket. It reads `monarch_kg/versions.json`, `species_data/`, the phenologs, the FDR table
and the candidate ranking.

```
<project_dir>/release/
  2026-03-06/
    disease/             one network's release, the folder that goes to the bucket (~30 MB)
      overview.json      everything the release page shows; its release block has the versions and settings
      README.md          what every table and column holds
      manifest.json      every file with its size, rows, sha256 and what it holds
      *.parquet          the tables below. The two the Downloads page lists are named for the network:
                         disease-phenologs.parquet and gene-disease-candidates.parquet here
    phenotype/           the same files for the phenotype network
  releases.json          the published list with this run's releases in it, written by 08_update_release_list.py
```

- **Terms.** Both networks have the same tables and columns. The human side of a phenolog, a disease or a human
  phenotype, is a term in both: `term_id`, `terms.parquet`, `n_known_terms`. Known means a curated edge of the base
  species' gene_to_disease or gene_to_phenotype table.
- **Settings.** The release holds one ranking: the FDR, k and rank metric given (defaults 0.05, 1000, hg).
- **Kept candidates.** A human phenotype's ranking scores thousands of genes, so the phenotype release keeps each
  one's best 100 (`-top_candidates`, recorded in the release block). Counts of scored and recovered genes, and a known
  gene's rank, are still over the whole ranking.
- **Versions.** The KG and PHENIO versions come from `monarch_kg/versions.json`, which step 01 writes. The code
  version comes from the image, or from git.
- **Rebuilds.** A rebuild replaces the release folder.
- **Other options.** `-c` and `-m` cap DuckDB's cores and memory, for a container.

The code is `python/phenologs_release_utils.py`. A change in the pipeline's files is a change to the `Run` views there.

| table | one row per | notes |
|---|---|---|
| `phenologs` (`<network>-phenologs.parquet`) | significant phenolog | ranked within its term by p |
| `candidates` (`gene-<network>-candidates.parquet`) | (term, gene) the release keeps | rank within the term and global rank across all terms, ties sharing a place; score, known flag, and support: per species, the term's phenologs that are phenotypes of an ortholog of the gene, of how many, through how many orthologs, and the (phenolog, ortholog) hits the score counts (the pipeline's `XSpecies_occurence`). The pipeline also appends a term's known genes its phenologs never reached, without a p-value; those are left out |
| `human_orthologs` | (model gene, human gene) | ortholog edges from the run's species to the base species |
| `additive` | (term, family, phenolog) | top 20 phenologs x top 40 families per term; the page re-slices by k, species and order |
| `terms`, `genes`, `phenotypes`, `families`, `species` | entity | counts the pages and browse lists show. `genes.in_reach`: a base species gene one of whose orthologs carries a phenotype the ranking can use |
| `known_associations` | curated gene-term pair | from the base species' gene_to_disease or gene_to_phenotype table, with the gene's rank and score in the term's whole ranking |
| `gene_phenotypes` | gene-phenotype annotation | every species |
| `search` | searchable thing | kind, id, label, detail, weight |

- **Counts.** Counts are 32 bit, so they arrive in the browser as plain numbers.
- **Downloads.** The site's Downloads page lists each network's phenologs and candidates and its `README.md`, with a
  DuckDB command per file that downloads it as TSV. `README.md` describes every column of every table in the release
  (a test checks that none is missing). Tables on the pages have a TSV button for the rows shown, as filtered and
  sorted.

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
- **Which release.** The site reads `releases.json` there and opens the newest KG version, or
  `?release=<graph_version>`. Its landing page and Downloads cover both networks; the rest of the site shows one,
  the disease network or `?network=phenotype`. The top bar has a row of pages for each. It asks for the list past any cache, so a new release
  shows at once, and says so when there is none yet.
- **Local data.** Served by `serve.py`, the site reads the run's `release/` folder instead, and `serve.py` lists the
  releases found there. `?data=<url>` points the site at any other releases folder.
- **Publishing.** See `RELEASING.md`.
- **Layout only grows.** A release folder never changes once uploaded. Later builds may add tables and columns, but
  never rename or remove one. Each release block's `schema` (`SCHEMA` in `phenologs_release_utils.py`) tells the site
  what that release has, and the site says when a release has a layout it doesn't read.
- **Deploying.** `.github/workflows/pages.yml` publishes this folder on every push to main that changes it.

**Fetching.** The parquet files are fetched whole on first use, so range support is not required. `serve.py` has range support anyway, for when the data layer moves to range reads.

**Outside sources.** Scripts come from jsdelivr (DuckDB-WASM) and cdnjs (D3), and fonts from Google Fonts.

## Later

- `links.dismech()` in `format.js` returns null for now. Once DisMech's disorder ids are known it can return
  a URL and the disease page will show the link.
- Leave-one-out validation and the temporal comparison between releases are not in the build yet.
