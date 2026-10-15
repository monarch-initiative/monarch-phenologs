"""
The additive decomposition behind a disease's candidates (Woods et al. 2013, as drawn in the paper's
case study figure). For each disease: its top phenologs by p-value; for each phenolog, the protein
families carrying that model phenotype that the species shares with the base species (the families the
phenolog test counts); each (family, phenolog) pair contributes (overlap / disease orthogroups) * (1 - p).
Families are ranked by their summed contribution.

One row per (disease, family, phenolog) is kept so the page can re-slice by species or by a smaller
k without going back to the pipeline. Depends on the phenologs, genes and known_associations tables.
"""
from .tables import write


def build(con, run, outdir, top_phenologs=20, top_families=40):
    con.execute("""
        CREATE TABLE additive AS
        WITH top_ph AS (
            SELECT disease_id, phenotype_id, phenotype_name, species, overlap, n_disease_orthogroups, pval, rank
            FROM phenologs WHERE rank <= {top_phenologs}),
        phen_fam AS (                       -- the protein families carrying each model phenotype, shared with the base species
            SELECT DISTINCT p.species, p.phenotype_id, o.family
            FROM g2p p JOIN g2o o ON o.species = p.species AND o.gene_id = p.gene_id
            JOIN common_families c ON c.species = o.species AND c.family = o.family),
        contrib AS (
            SELECT t.disease_id, f.family, t.rank AS phenolog_rank, t.phenotype_id, t.phenotype_name, t.species,
                   (t.overlap::DOUBLE / t.n_disease_orthogroups) * (1 - t.pval) AS score
            FROM top_ph t JOIN phen_fam f ON f.species = t.species AND f.phenotype_id = t.phenotype_id),
        fam_tot AS (                        -- summed as decimals, so the totals and the order of tied families are the same every build
            SELECT disease_id, family, sum(score::DECIMAL(38, 18))::DOUBLE AS family_total,
                   row_number() OVER (PARTITION BY disease_id ORDER BY sum(score::DECIMAL(38, 18)) DESC, family)::INTEGER AS family_rank
            FROM contrib GROUP BY disease_id, family),
        base_genes AS (                     -- the base species' genes of each family, split by whether the disease already knows them
            SELECT ft.disease_id, ft.family,
                   string_agg(g.symbol, ', ' ORDER BY g.symbol) FILTER (WHERE k.gene_id IS NOT NULL) AS known_symbols,
                   string_agg(g.symbol, ', ' ORDER BY g.symbol) FILTER (WHERE k.gene_id IS NULL)     AS other_symbols
            FROM fam_tot ft
            JOIN genes g ON g.family = ft.family AND g.species = '{base}'
            LEFT JOIN known_associations k ON k.gene_id = g.gene_id AND k.disease_id = ft.disease_id
            WHERE ft.family_rank <= {top_families}
            GROUP BY ft.disease_id, ft.family)
        SELECT c.disease_id, c.family, ft.family_rank, ft.family_total, c.phenolog_rank::INTEGER AS phenolog_rank,
               c.phenotype_id, c.phenotype_name, c.species, c.score,
               bg.known_symbols, bg.other_symbols
        FROM contrib c
        JOIN fam_tot ft USING (disease_id, family)
        LEFT JOIN base_genes bg USING (disease_id, family)
        WHERE ft.family_rank <= {top_families}
        ORDER BY c.disease_id, ft.family_rank, c.phenolog_rank
    """.format(top_phenologs=top_phenologs, top_families=top_families, base=run.base))
    return write(con, "additive", outdir)
