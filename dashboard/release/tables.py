"""
The entity tables the dashboard searches and browses. Each is one SQL statement over the views in
sources.py (and over tables built before it), materialised in DuckDB and then written out as parquet.
Counts are 32 bit so they arrive in the browser as plain numbers.
"""
import os

# Built in this order. Later tables may refer to earlier ones
TABLES = ["phenologs", "candidates", "genes", "known_associations", "gene_phenotypes",
          "phenotypes", "families", "diseases", "species", "search"]


def sql(name, run):
    base = run.base
    return {
        # Every significant phenolog of the run, ranked within its disease by p-value
        "phenologs": """
            SELECT p.*, row_number() OVER (PARTITION BY disease_id ORDER BY pval, species, phenotype_id)::INTEGER AS rank
            FROM pooled p""",

        # Every candidate gene the ranking wrote, ranked within its disease. One row per (disease, gene)
        "candidates": """
            SELECT disease_id, gene_id, symbol, family, score, known, species_occurrence,
                   CASE WHEN species_occurrence = '' THEN 0
                        ELSE length(species_occurrence) - length(replace(species_occurrence, ',', '')) + 1 END::INTEGER AS n_species,
                   row_number() OVER (PARTITION BY disease_id ORDER BY score, symbol)::INTEGER AS rank
            FROM cand_raw""",

        # Every gene of every species in the run that has a protein family or a phenotype
        "genes": """
            WITH base AS (
                SELECT species, gene_id, max(family) AS family
                FROM (SELECT species, gene_id, family FROM g2o UNION ALL SELECT species, gene_id, NULL FROM g2p)
                GROUP BY species, gene_id),
            np AS (SELECT gene_id, count(DISTINCT phenotype_id) AS n FROM g2p GROUP BY gene_id),
            nk AS (SELECT gene_id, count(DISTINCT disease_id) AS n FROM g2d GROUP BY gene_id),
            pr AS (SELECT gene_id, count(DISTINCT disease_id) AS n, count(DISTINCT disease_id) FILTER (WHERE NOT known) AS novel, min(rank) AS best_rank
                   FROM candidates GROUP BY gene_id)
            SELECT b.species, b.gene_id, coalesce(n.symbol, b.gene_id) AS symbol, b.family,
                   coalesce(np.n, 0)::INTEGER AS n_phenotypes,
                   coalesce(nk.n, 0)::INTEGER AS n_known_diseases,
                   coalesce(pr.n, 0)::INTEGER AS n_predicted_diseases,
                   coalesce(pr.novel, 0)::INTEGER AS n_novel_diseases,
                   pr.best_rank::INTEGER AS best_rank
            FROM base b
            LEFT JOIN gene_names n USING (gene_id)
            LEFT JOIN np USING (gene_id) LEFT JOIN nk USING (gene_id) LEFT JOIN pr USING (gene_id)""",

        # The base species' curated gene to disease associations the predictions are measured against
        "known_associations": """
            SELECT DISTINCT gene_id, symbol, disease_id, disease_name, predicate, family FROM g2d""",

        # Gene to phenotype annotations of every species (what a phenotype's gene list and a gene's phenotype list show)
        "gene_phenotypes": """
            SELECT DISTINCT species, gene_id, phenotype_id, phenotype_name FROM g2p""",

        # Every annotated phenotype of every species, with how often it turns up as a phenolog
        "phenotypes": """
            WITH ann AS (
                SELECT species, phenotype_id, any_value(phenotype_name) AS name,
                       count(DISTINCT gene_id) AS n_genes, count(DISTINCT family) AS n_orthogroups
                FROM g2p GROUP BY species, phenotype_id),
            ph AS (
                SELECT species, phenotype_id, count(*) AS n_phenologs, count(DISTINCT disease_id) AS n_diseases, min(pval) AS best_pval
                FROM phenologs GROUP BY species, phenotype_id)
            SELECT a.species, a.phenotype_id, a.name, a.n_genes::INTEGER AS n_genes, a.n_orthogroups::INTEGER AS n_orthogroups,
                   coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_diseases, 0)::INTEGER AS n_diseases, ph.best_pval
            FROM ann a LEFT JOIN ph USING (species, phenotype_id)""",

        # PANTHER families: who is in them, and what they predict
        "families": """
            WITH fam AS (SELECT family, species, gene_id FROM genes WHERE family IS NOT NULL),
            hum AS (SELECT family, string_agg(symbol, ', ' ORDER BY symbol) AS human_symbols, count(*) AS n
                    FROM genes WHERE species = '{base}' AND family IS NOT NULL GROUP BY family),
            pr AS (SELECT family, count(DISTINCT disease_id) AS n, count(DISTINCT disease_id) FILTER (WHERE NOT known) AS novel, min(score) AS best_score
                   FROM candidates GROUP BY family),
            kn AS (SELECT family, count(DISTINCT disease_id) AS n FROM g2d WHERE family IS NOT NULL GROUP BY family)
            SELECT f.family, count(*)::INTEGER AS n_genes, count(DISTINCT f.species)::INTEGER AS n_species,
                   coalesce(any_value(hum.n), 0)::INTEGER AS n_base_genes, any_value(hum.human_symbols) AS base_symbols,
                   coalesce(any_value(pr.n), 0)::INTEGER AS n_predicted_diseases, coalesce(any_value(pr.novel), 0)::INTEGER AS n_novel_diseases,
                   any_value(pr.best_score) AS best_score, coalesce(any_value(kn.n), 0)::INTEGER AS n_known_diseases
            FROM fam f LEFT JOIN hum USING (family) LEFT JOIN pr USING (family) LEFT JOIN kn USING (family)
            GROUP BY f.family""",

        # The base species' diseases that have at least one phenolog
        "diseases": """
            WITH bysp AS (SELECT disease_id, species, count(*) AS n FROM phenologs GROUP BY disease_id, species),
            bysp_json AS (
                SELECT disease_id, '{{' || string_agg(format('"{{}}":{{}}', species, n), ',' ORDER BY n DESC) || '}}' AS by_species,
                       arg_max(species, n) AS top_species
                FROM bysp GROUP BY disease_id),
            ph AS (SELECT disease_id, any_value(disease_name) AS name, count(*) AS n_phenologs, count(DISTINCT species) AS n_species, min(pval) AS best_pval
                   FROM phenologs GROUP BY disease_id),
            kn AS (SELECT disease_id, count(DISTINCT gene_id) AS n FROM g2d GROUP BY disease_id),
            cd AS (SELECT disease_id, count(*) AS n, count(*) FILTER (WHERE NOT known) AS novel, count(*) FILTER (WHERE known) AS recovered
                   FROM candidates GROUP BY disease_id)
            SELECT ph.disease_id, ph.name, coalesce(kn.n, 0)::INTEGER AS n_known_genes, ph.n_phenologs::INTEGER AS n_phenologs,
                   ph.n_species::INTEGER AS n_species, ph.best_pval,
                   coalesce(cd.n, 0)::INTEGER AS n_candidates, coalesce(cd.novel, 0)::INTEGER AS n_novel_candidates,
                   coalesce(cd.recovered, 0)::INTEGER AS n_known_recovered, bj.top_species, bj.by_species
            FROM ph LEFT JOIN kn USING (disease_id) LEFT JOIN cd USING (disease_id) LEFT JOIN bysp_json bj USING (disease_id)""",

        # The species of the run with everything the release page and the species pages show
        "species": """
            WITH ph AS (SELECT species, count(*) AS n_phenologs, count(DISTINCT disease_id) AS n_diseases, count(DISTINCT phenotype_id) AS n_phenotypes_matched
                        FROM phenologs GROUP BY species),
            occ AS (SELECT known, unnest(string_split(species_occurrence, ',')) AS item FROM candidates WHERE species_occurrence <> ''),
            cs AS (SELECT split_part(item, ':', 1) AS species, count(*) AS n_candidates, count(*) FILTER (WHERE NOT known) AS n_novel
                   FROM occ GROUP BY 1),
            cc AS (SELECT species_b AS species, n_common FROM common_counts WHERE species_a = '{base}'),
            fd AS (SELECT species, pvalue_cutoff, phenolog_count FROM fdr_table WHERE fdr = {fdr})
            SELECT i.*, coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_diseases, 0)::INTEGER AS n_diseases,
                   coalesce(ph.n_phenotypes_matched, 0)::INTEGER AS n_phenotypes_matched,
                   coalesce(cs.n_candidates, 0)::INTEGER AS n_candidates, coalesce(cs.n_novel, 0)::INTEGER AS n_novel_candidates,
                   cc.n_common::INTEGER AS n_common_with_base, fd.pvalue_cutoff AS fdr_cutoff
            FROM species_info i
            LEFT JOIN ph USING (species) LEFT JOIN cs USING (species) LEFT JOIN cc USING (species) LEFT JOIN fd USING (species)""",

        # One row per thing the search box can find
        "search": """
            SELECT 'disease' AS kind, disease_id AS id, name AS label,
                   format('{{}} phenologs · {{}} known genes', n_phenologs, n_known_genes) AS detail, n_phenologs::INTEGER AS weight
            FROM diseases
            UNION ALL
            SELECT 'gene', gene_id, symbol, species || coalesce(' · ' || family, ''),
                   (n_predicted_diseases + 10 * n_known_diseases + n_phenotypes)::INTEGER
            FROM genes WHERE n_phenotypes > 0 OR n_known_diseases > 0 OR n_predicted_diseases > 0
            UNION ALL
            SELECT 'phenotype', phenotype_id, name, format('{{}} · {{}} phenologs', species, n_phenologs), n_phenologs::INTEGER
            FROM phenotypes
            UNION ALL
            SELECT 'family', family, family, coalesce(base_symbols, ''), n_predicted_diseases::INTEGER
            FROM families
            UNION ALL
            SELECT 'species', species, binomial, common, 1000000 FROM species""",
    }[name].format(base=base, fdr=run.fdr)


def build(con, run, outdir):
    """ Materialise every table in order and write each as parquet. Returns {name: row count} """

    counts = {}
    for name in TABLES:
        con.execute("CREATE TABLE {} AS {}".format(name, sql(name, run)))
        counts[name] = write(con, name, outdir)
    return counts


def write(con, name, outdir):
    path = os.path.join(outdir, name + ".parquet")
    con.execute("COPY {} TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)".format(name, path))
    return con.execute("SELECT count(*) FROM {}".format(name)).fetchone()[0]
