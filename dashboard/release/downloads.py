"""
The run's raw results as single tables, for download: every candidate ranking it made, every disease-phenotype
pair it tested, the significant phenologs and the FDR table, with a README on their columns. COPY streams each
straight from the pipeline's files, so nothing big is held in memory. Depends on the phenologs table.
"""
import os

DOWNLOADS = {
    # Grouped by ranking, then by disease. Each disease's genes keep the pipeline's order
    "candidate_rankings.parquet": "SELECT * FROM rankings",

    # A pair is a phenolog when its p-value is at or below its species' cutoff at the run's FDR (as the pipeline pools them)
    "phenologs_all_tested.parquet": """
        SELECT a.*, a.pval <= f.pvalue_cutoff AS significant
        FROM all_phenologs a LEFT JOIN fdr_table f ON f.species = a.species AND f.fdr = {fdr}
        ORDER BY a.species, a.pval, a.disease_id, a.phenotype_id""",

    "phenologs_significant.tsv.gz": "SELECT * FROM phenologs ORDER BY pval, disease_id, species, phenotype_id",
    "fdr_table.tsv": "SELECT * FROM fdr_table ORDER BY species, fdr",
}


def build(con, run, outdir):
    """ Write every download and the README to <outdir>/downloads """

    folder = os.path.join(outdir, "downloads")
    os.makedirs(folder, exist_ok=True)
    for name, query in DOWNLOADS.items():
        options = "FORMAT PARQUET, COMPRESSION ZSTD" if name.endswith(".parquet") else "FORMAT CSV, DELIMITER '\\t', HEADER"
        con.execute("COPY ({}) TO '{}' ({})".format(query.format(fdr=run.fdr), os.path.join(folder, name), options))

    with open(os.path.join(folder, "README.md"), "w") as fh:
        fh.write(README.format(graph_version=run.graph_version, fdr=run.fdr, kneighbs=run.kneighbs, rank_metric=run.rank_metric,
                               phenio=" and phenio {}".format(run.phenio_version) if run.phenio_version else ""))


README = """# Phenologs downloads, Monarch KG {graph_version}

The raw results of one phenologs pipeline run, on Monarch KG {graph_version}{phenio}. Human diseases are compared
with the phenotypes of each model species through the PANTHER protein families of their genes. A disease and a
phenotype whose families overlap more than chance at FDR {fdr} make a phenolog. Each disease's phenologs then rank
the human genes they point to as candidates for that disease.

| file | one row per |
|---|---|
| `candidate_rankings.parquet` | (ranking, disease, gene): every candidate gene ranking the run made |
| `phenologs_all_tested.parquet` | (disease, species, phenotype): every pair the run tested |
| `phenologs_significant.tsv.gz` | phenolog: the significant pairs, as a tab separated table |
| `fdr_table.tsv` | (species, FDR level): the p-value cutoffs |

`../manifest.json` lists every file of the release with its size, row count and sha256.

The parquet files open in DuckDB, Python and R:

    duckdb -c "SELECT * FROM 'candidate_rankings.parquet' WHERE kneighbs = {kneighbs} AND rank_metric = '{rank_metric}' LIMIT 10"
    python -c "import pandas as pd; print(pd.read_parquet('phenologs_all_tested.parquet').head())"
    Rscript -e "print(head(arrow::read_parquet('phenologs_all_tested.parquet')))"

## candidate_rankings.parquet

For each disease with a phenolog, the human genes its k best phenologs (lowest p-value first) reach. A phenolog
reaches a human gene through a gene of the model species that carries the phenotype and is an ortholog of it: one
hit. Every ranking the run made is here, each k with each of two scores. The site shows k = {kneighbs}, {rank_metric}.
Rows are grouped by ranking, then by disease. Within a disease they keep the pipeline's order: the scored genes,
strongest first, then the disease's known genes no phenolog reached.

| column | |
|---|---|
| `kneighbs` | k, the number of phenologs the ranking draws on |
| `rank_metric` | `hg` (hypergeometric) or `nb` (naive Bayes) |
| `disease_id` | Mondo disease |
| `gene_id` | human gene (HGNC) |
| `symbol` | gene symbol |
| `known` | the gene is curated for the disease |
| `species_hits` | the gene's hits per model species, as `Mus-musculus:3,Danio-rerio:1` |
| `score` | smaller is stronger. `hg`: hypergeometric p-value of the gene's hits, out of every hit the ranking drew. `nb`: the product, over the gene's hits, of 1 - (the phenolog's overlap / the disease's known genes) x (1 - the phenolog's p-value). Empty for known genes no phenolog reached |
| `phenotypes` | names of the model phenotypes behind the hits, `;` separated |
| `family` | the gene's PANTHER family |

## phenologs_all_tested.parquet

Every disease-phenotype pair the run tested, for every model species, by species and then p-value.

| column | |
|---|---|
| `disease_id`, `disease_name` | Mondo disease |
| `species` | model species, as `Genus-species` |
| `phenotype_id`, `phenotype_name` | the model species' phenotype |
| `overlap` | protein families carrying both the disease and the phenotype |
| `n_disease_orthogroups` | protein families of the disease's genes, counting only families the two species share |
| `n_phenotype_orthogroups` | protein families of the phenotype's genes, counting only families the two species share |
| `n_common` | protein families present in both species: the background of the test |
| `pval` | hypergeometric p-value of the overlap |
| `significant` | `pval` is at or below the species' cutoff at FDR {fdr} (`fdr_table.tsv`): the pair is a phenolog |

## phenologs_significant.tsv.gz

The phenologs: the rows of `phenologs_all_tested.parquet` with `significant` true, without that column, ordered by
p-value. One more column, `rank`: the phenolog's place among its disease's phenologs, by p-value.

## fdr_table.tsv

| column | |
|---|---|
| `species` | model species |
| `fdr` | false discovery rate level |
| `pvalue_cutoff` | largest p-value still called a phenolog at this FDR |
| `phenolog_count` | real disease-phenotype comparisons at or below the cutoff |
| `expected_false_positives` | randomised comparisons at or below the cutoff, per trial (mean over the trials) |
| `realized_fdr` | expected false positives / phenologs |
| `n_random_trials` | randomised trials the cutoffs come from |
"""
