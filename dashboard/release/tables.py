"""
The entity tables the dashboard searches and browses. Each is one SQL statement over the views in
sources.py (and over tables built before it), materialised in DuckDB and then written out as parquet.
Counts are 32 bit so they arrive in the browser as plain numbers.
"""
import os

# Built in this order. Later tables may refer to earlier ones. Internal ones are steps on the way, not written out
TABLES = ["phenologs", "human_orthologs", "support", "candidates", "genes", "known_associations", "gene_phenotypes",
          "phenotypes", "families", "diseases", "species", "search"]
INTERNAL = {"support"}


def sql(name, run):
    base = run.base
    return {
        # Every significant phenolog of the run, ranked within its disease by p-value
        "phenologs": """
            SELECT p.*, row_number() OVER (PARTITION BY disease_id ORDER BY pval, species, phenotype_id)::INTEGER AS rank
            FROM pooled p""",

        # How a disease's k nearest phenologs reach each base species gene, per model species: the phenologs whose
        # phenotype is annotated to an ortholog of the gene, how many orthologs that takes, and the (phenolog,
        # ortholog) hits the hg score counts (n_hits is the pipeline's XSpecies_occurence). n_used is how many of
        # the disease's phenologs with that species the ranking used
        "support": """
            WITH topk AS (SELECT disease_id, species, phenotype_id FROM phenologs WHERE rank <= {k}),
            used AS (SELECT disease_id, species, count(*) AS n FROM topk GROUP BY disease_id, species),
            hits AS (SELECT t.disease_id, t.species, t.phenotype_id, g.gene_id AS ortholog_id, o.human_gene_id AS gene_id
                     FROM topk t
                     JOIN ranking_g2p g ON g.species = t.species AND g.phenotype_id = t.phenotype_id
                     JOIN human_orthologs o ON o.species = g.species AND o.gene_id = g.gene_id)
            SELECT h.disease_id, h.gene_id, h.species,
                   count(DISTINCT h.phenotype_id)::INTEGER AS n_phenologs, any_value(u.n)::INTEGER AS n_used,
                   count(DISTINCT h.ortholog_id)::INTEGER AS n_orthologs, count(*)::INTEGER AS n_hits
            FROM hits h JOIN used u USING (disease_id, species)
            GROUP BY h.disease_id, h.gene_id, h.species""",

        # Every candidate gene the ranking scored, ranked within its disease by score (ties share a place).
        # One row per (disease, gene). The pipeline also appends the disease's known genes its phenologs never
        # reached, with no p-value. Those are not candidates and are left out here (known_associations still has them).
        # support is "species:phenologs:of the disease's:orthologs:hits" per species, most phenologs first.
        # global_rank places every pair, known or novel, in one list across all diseases by score (ties share a rank)
        "candidates": """
            WITH c AS (
                SELECT disease_id, gene_id, symbol, family, score, known,
                       rank() OVER (PARTITION BY disease_id ORDER BY score)::INTEGER AS rank
                FROM cand_raw WHERE score IS NOT NULL),
            s AS (
                SELECT disease_id, gene_id,
                       string_agg(format('{{}}:{{}}:{{}}:{{}}:{{}}', species, n_phenologs, n_used, n_orthologs, n_hits), ','
                                  ORDER BY n_phenologs DESC, species) AS support,
                       count(*) AS n_species, sum(n_phenologs) AS n_phenologs
                FROM support GROUP BY disease_id, gene_id),
            names AS (SELECT disease_id, any_value(disease_name) AS disease_name FROM pooled GROUP BY disease_id)
            SELECT c.*, coalesce(s.support, '') AS support, coalesce(s.n_species, 0)::INTEGER AS n_species,
                   coalesce(s.n_phenologs, 0)::INTEGER AS n_phenologs, n.disease_name,
                   rank() OVER (ORDER BY c.score)::INTEGER AS global_rank
            FROM c LEFT JOIN s USING (disease_id, gene_id) LEFT JOIN names n USING (disease_id)
            ORDER BY c.score, c.disease_id, c.symbol""",

        # Every gene of every species in the run that has a protein family, a phenotype or a curated disease.
        # n_ortholog_species counts the other species of the run the gene has an ortholog edge to. in_reach (base
        # species genes only): one of its orthologs carries a phenotype the ranking can use, so a disease's
        # ranking can reach it at all
        "genes": """
            WITH base AS (
                SELECT species, gene_id, max(family) AS family
                FROM (SELECT species, gene_id, family FROM g2o UNION ALL SELECT species, gene_id, NULL FROM g2p
                      UNION ALL SELECT '{base}', gene_id, family FROM g2d)
                GROUP BY species, gene_id),
            reach AS (SELECT DISTINCT o.human_gene_id AS gene_id
                      FROM ranking_g2p g JOIN human_orthologs o ON o.species = g.species AND o.gene_id = g.gene_id),
            os AS (SELECT o.gene_id, count(DISTINCT s.species) AS n
                   FROM (SELECT gene_a AS gene_id, taxon_b AS taxon FROM edges UNION SELECT gene_b, taxon_a FROM edges) o
                   JOIN species_info s ON s.taxon_id = o.taxon GROUP BY o.gene_id),
            np AS (SELECT gene_id, count(DISTINCT phenotype_id) AS n FROM g2p GROUP BY gene_id),
            nk AS (SELECT gene_id, count(DISTINCT disease_id) AS n FROM g2d GROUP BY gene_id),
            pr AS (SELECT gene_id, count(DISTINCT disease_id) AS n, count(DISTINCT disease_id) FILTER (WHERE NOT known) AS novel, min(rank) AS best_rank
                   FROM candidates GROUP BY gene_id)
            SELECT b.species, b.gene_id, coalesce(n.symbol, b.gene_id) AS symbol, b.family,
                   coalesce(os.n, 0)::INTEGER AS n_ortholog_species,
                   coalesce(np.n, 0)::INTEGER AS n_phenotypes,
                   coalesce(nk.n, 0)::INTEGER AS n_known_diseases,
                   coalesce(pr.n, 0)::INTEGER AS n_predicted_diseases,
                   coalesce(pr.novel, 0)::INTEGER AS n_novel_diseases,
                   pr.best_rank::INTEGER AS best_rank,
                   CASE WHEN b.species = '{base}' THEN b.gene_id IN (SELECT gene_id FROM reach) END AS in_reach
            FROM base b
            LEFT JOIN gene_names n USING (gene_id) LEFT JOIN os USING (gene_id)
            LEFT JOIN np USING (gene_id) LEFT JOIN nk USING (gene_id) LEFT JOIN pr USING (gene_id)""",

        # Every gene of the run's other species with an ortholog edge to a base species gene, one row per pair
        "human_orthologs": """
            WITH base AS (SELECT taxon_id FROM species_info WHERE species = '{base}'),
            pairs AS (SELECT taxon_a AS taxon, gene_a AS gene_id, gene_b AS human_gene_id FROM edges WHERE taxon_b = (SELECT taxon_id FROM base)
                      UNION SELECT taxon_b, gene_b, gene_a FROM edges WHERE taxon_a = (SELECT taxon_id FROM base))
            SELECT s.species, p.gene_id, p.human_gene_id, coalesce(n.symbol, p.human_gene_id) AS human_symbol
            FROM pairs p JOIN species_info s ON s.taxon_id = p.taxon LEFT JOIN gene_names n ON n.gene_id = p.human_gene_id
            ORDER BY s.species, p.gene_id, human_symbol""",

        # The base species' curated gene to disease associations the predictions are measured against
        "known_associations": """
            SELECT DISTINCT gene_id, symbol, disease_id, disease_name, predicate, family FROM g2d""",

        # Gene to phenotype annotations of every species (what a phenotype's gene list and a gene's phenotype list show)
        "gene_phenotypes": """
            SELECT DISTINCT species, gene_id, phenotype_id, phenotype_name FROM g2p""",

        # Every annotated phenotype of every species, with how often it turns up as a phenolog.
        # n_known_diseases counts the curated diseases of the genes annotated with the phenotype. Curated gene to
        # disease associations are the base species' only, so model organism phenotypes have none
        "phenotypes": """
            WITH ann AS (
                SELECT g.species, g.phenotype_id, any_value(g.phenotype_name) AS name,
                       count(DISTINCT g.gene_id) AS n_genes, count(DISTINCT g.family) AS n_orthogroups,
                       CASE WHEN g.species <> '{base}' THEN count(DISTINCT c.family) END AS n_shared_orthogroups
                FROM g2p g LEFT JOIN common_families c ON c.species = g.species AND c.family = g.family
                GROUP BY g.species, g.phenotype_id),
            kd AS (
                SELECT p.species, p.phenotype_id, count(DISTINCT d.disease_id) AS n
                FROM g2p p JOIN g2d d USING (gene_id) GROUP BY p.species, p.phenotype_id),
            ph AS (
                SELECT species, phenotype_id, count(*) AS n_phenologs, count(DISTINCT disease_id) AS n_diseases, min(pval) AS best_pval
                FROM phenologs GROUP BY species, phenotype_id)
            SELECT a.species, a.phenotype_id, a.name, a.n_genes::INTEGER AS n_genes, a.n_orthogroups::INTEGER AS n_orthogroups,
                   a.n_shared_orthogroups::INTEGER AS n_shared_orthogroups,
                   coalesce(kd.n, 0)::INTEGER AS n_known_diseases,
                   coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_diseases, 0)::INTEGER AS n_diseases, ph.best_pval
            FROM ann a LEFT JOIN kd USING (species, phenotype_id) LEFT JOIN ph USING (species, phenotype_id)""",

        # PANTHER families: who is in them, and what they predict. A disease the family is ranked for is novel
        # for the family when none of its base species genes is curated for that disease
        "families": """
            WITH fam AS (SELECT family, species, gene_id FROM genes WHERE family IS NOT NULL),
            hum AS (SELECT family, string_agg(symbol, ', ' ORDER BY symbol) AS human_symbols, count(*) AS n
                    FROM genes WHERE species = '{base}' AND family IS NOT NULL GROUP BY family),
            kpairs AS (SELECT DISTINCT family, disease_id FROM g2d WHERE family IS NOT NULL),
            pr AS (SELECT c.family, count(DISTINCT c.disease_id) AS n, min(c.score) AS best_score,
                          count(DISTINCT c.disease_id) FILTER (WHERE k.disease_id IS NULL) AS novel
                   FROM candidates c LEFT JOIN kpairs k ON k.family = c.family AND k.disease_id = c.disease_id
                   GROUP BY c.family),
            kn AS (SELECT family, count(*) AS n FROM kpairs GROUP BY family)
            SELECT f.family, count(*)::INTEGER AS n_genes, count(DISTINCT f.species)::INTEGER AS n_species,
                   coalesce(any_value(hum.n), 0)::INTEGER AS n_base_genes, any_value(hum.human_symbols) AS base_symbols,
                   coalesce(any_value(pr.n), 0)::INTEGER AS n_predicted_diseases, coalesce(any_value(pr.novel), 0)::INTEGER AS n_novel_diseases,
                   any_value(pr.best_score) AS best_score, coalesce(any_value(kn.n), 0)::INTEGER AS n_known_diseases
            FROM fam f LEFT JOIN hum USING (family) LEFT JOIN pr USING (family) LEFT JOIN kn USING (family)
            GROUP BY f.family""",

        # The base species' diseases that have at least one phenolog
        "diseases": """
            WITH bysp AS (SELECT disease_id, species, count(*) AS n FROM phenologs GROUP BY disease_id, species),
            bysp_json AS (              -- ties go to the species first in plot order
                SELECT b.disease_id, '{{' || string_agg(format('"{{}}":{{}}', b.species, b.n), ',' ORDER BY b.n DESC, s.plot_order) || '}}' AS by_species,
                       first(b.species ORDER BY b.n DESC, s.plot_order) AS top_species
                FROM bysp b JOIN species_info s USING (species) GROUP BY b.disease_id),
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

        # The species of the run with everything the release page and the species pages show. n_phenotypes_annotated
        # counts every phenotype the species' genes carry (the species table's "Unique phenotypes" only counts those
        # on genes with an ortholog)
        "species": """
            WITH ph AS (SELECT species, count(*) AS n_phenologs, count(DISTINCT disease_id) AS n_diseases, count(DISTINCT phenotype_id) AS n_phenotypes_matched
                        FROM phenologs GROUP BY species),
            cs AS (SELECT s.species, count(*) AS n_candidates, count(*) FILTER (WHERE NOT c.known) AS n_novel
                   FROM support s JOIN candidates c USING (disease_id, gene_id) GROUP BY s.species),
            cc AS (SELECT species_b AS species, n_common FROM common_counts WHERE species_a = '{base}'),
            ho AS (SELECT species, count(DISTINCT gene_id) AS n FROM human_orthologs GROUP BY species),
            fd AS (SELECT species, pvalue_cutoff, phenolog_count FROM fdr_table WHERE fdr = {fdr}),
            pa AS (SELECT species, count(DISTINCT phenotype_id) AS n FROM g2p GROUP BY species)
            SELECT i.*, coalesce(pa.n, 0)::INTEGER AS n_phenotypes_annotated,
                   coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_diseases, 0)::INTEGER AS n_diseases,
                   coalesce(ph.n_phenotypes_matched, 0)::INTEGER AS n_phenotypes_matched,
                   coalesce(cs.n_candidates, 0)::INTEGER AS n_candidates, coalesce(cs.n_novel, 0)::INTEGER AS n_novel_candidates,
                   cc.n_common::INTEGER AS n_common_with_base, ho.n::INTEGER AS n_human_ortholog_genes, fd.pvalue_cutoff AS fdr_cutoff
            FROM species_info i
            LEFT JOIN ph USING (species) LEFT JOIN cs USING (species) LEFT JOIN cc USING (species)
            LEFT JOIN ho USING (species) LEFT JOIN fd USING (species) LEFT JOIN pa USING (species)""",

        # One row per thing the search box can find
        "search": """
            SELECT 'disease' AS kind, disease_id AS id, name AS label, NULL AS species,
                   format('{{:,}} phenologs · {{:,}} known genes', n_phenologs, n_known_genes) AS detail, n_phenologs::INTEGER AS weight
            FROM diseases
            UNION ALL
            SELECT 'gene', gene_id, symbol, species, replace(species, '-', ' ') || coalesce(' · ' || family, ''),
                   (n_predicted_diseases + 10 * n_known_diseases + n_phenotypes)::INTEGER
            FROM genes WHERE n_phenotypes > 0 OR n_known_diseases > 0 OR n_predicted_diseases > 0
            UNION ALL
            SELECT 'phenotype', phenotype_id, name, species,
                   CASE WHEN species = '{base}' THEN format('{{}} · {{:,}} genes', replace(species, '-', ' '), n_genes)
                        ELSE format('{{}} · {{:,}} phenologs', replace(species, '-', ' '), n_phenologs) END,
                   n_phenologs::INTEGER
            FROM phenotypes
            UNION ALL
            SELECT 'family', family, family, NULL, coalesce(base_symbols, ''), n_predicted_diseases::INTEGER
            FROM families
            UNION ALL
            SELECT 'species', species, binomial, species, common, 1000000 FROM species""",
    }[name].format(base=base, fdr=run.fdr, k=run.kneighbs)


def build(con, run, outdir):
    """ Materialise every table in order and write each as parquet. Returns {name: row count} """

    counts = {}
    for name in TABLES:
        con.execute("CREATE TABLE {} AS {}".format(name, sql(name, run)))
        if name not in INTERNAL:
            counts[name] = write(con, name, outdir)
    return counts


def support_mismatches(con):
    """
    (pairs, diseases) whose support hits differ from the pipeline's XSpecies_occurence. Runs from before
    phenologs were sorted deterministically can break p-value ties at the k boundary another way
    """

    return con.execute("""
        WITH occ AS (SELECT disease_id, gene_id, split_part(item, ':', 1) AS species, split_part(item, ':', 2)::INTEGER AS n_hits
                     FROM (SELECT disease_id, gene_id, unnest(string_split(species_hits, ',')) AS item
                           FROM cand_raw WHERE score IS NOT NULL AND species_hits <> ''))
        SELECT count(DISTINCT (disease_id, gene_id)), count(DISTINCT disease_id)
        FROM occ FULL JOIN support s USING (disease_id, gene_id, species)
        WHERE occ.n_hits IS DISTINCT FROM s.n_hits""").fetchone()


def write(con, name, outdir):
    path = os.path.join(outdir, name + ".parquet")
    con.execute("COPY {} TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)".format(name, path))
    return con.execute("SELECT count(*) FROM {}".format(name)).fetchone()[0]
