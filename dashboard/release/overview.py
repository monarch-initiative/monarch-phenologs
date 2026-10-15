"""
overview.json: everything the release page shows before any parquet is opened. Headline tiles,
per species stats (the metric cards and the tree), the FDR table, common ortholog counts between
species, the top lists, and what the graph gave the pipeline (every taxon's extraction counts and how
far the curated disease genes get). Depends on the tables built by tables.py.
"""
import json
import os
import datetime


def rows(con, query):
    cur = con.execute(query)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


def build(con, run, outdir, table_counts):
    base = run.base
    species = rows(con, "SELECT * FROM species ORDER BY plot_order")
    for s in species:
        s["is_base"] = bool(s["is_base"])

    funnel = known_funnel(con)
    tiles = rows(con, """
        SELECT (SELECT count(*) FROM diseases)                                   AS diseases,
               (SELECT count(*) FROM phenologs)                                  AS phenologs,
               (SELECT count(DISTINCT (species, phenotype_id)) FROM phenologs)   AS model_phenotypes,
               (SELECT count(*) FROM candidates)                                 AS candidate_pairs,
               (SELECT count(*) FROM candidates WHERE NOT known)                 AS novel_pairs,
               (SELECT count(DISTINCT gene_id) FROM candidates WHERE NOT known)  AS candidate_genes,
               (SELECT count(DISTINCT family) FROM candidates)                   AS families_implicated,
               (SELECT count(*) FROM species WHERE NOT is_base)                  AS model_species""")[0]
    tiles.update(known_recovered=funnel["recovered"], known_pairs=funnel["curated"], known_in_reach=funnel["in_reach"])

    top = {
        "diseases": rows(con, """SELECT disease_id AS id, name AS label, n_phenologs AS value, n_species, n_novel_candidates, top_species
                                 FROM diseases ORDER BY n_phenologs DESC, name LIMIT 15"""),
        # Phenotypes carried by a single gene match every disease that gene is tied to, so the list is
        # kept to phenotypes with a few protein families behind them
        "phenotypes": rows(con, """SELECT phenotype_id AS id, name AS label, species, n_phenologs AS value, n_diseases, n_orthogroups
                                   FROM phenotypes WHERE n_phenologs > 0 AND n_orthogroups >= 5 ORDER BY n_phenologs DESC, name LIMIT 15"""),
        "genes": rows(con, """SELECT gene_id AS id, symbol AS label, family, n_predicted_diseases AS value, n_novel_diseases, n_known_diseases
                              FROM genes WHERE species = '{}' ORDER BY n_predicted_diseases DESC, symbol LIMIT 15""".format(base)),
        "families": rows(con, """SELECT family AS id, family AS label, base_symbols, n_predicted_diseases AS value, n_known_diseases, n_species
                                 FROM families ORDER BY n_predicted_diseases DESC, family LIMIT 15"""),
    }

    out = {
        "release": {"graph_version": run.graph_version,
                    "phenio_version": run.phenio_version,
                    "built": datetime.date.today().isoformat(),
                    "base_species": base,
                    "prediction_network": run.prediction_network,
                    "fdr": run.fdr, "kneighbs": run.kneighbs, "rank_metric": run.rank_metric,
                    "project_dir": os.path.basename(run.project_dir.rstrip("/")),
                    "table_rows": table_counts},
        "tiles": tiles,
        "species": species,
        "tree": run.species_tree(species),
        "fdr": rows(con, "SELECT * FROM fdr_table ORDER BY species, fdr"),
        "common_orthologs": rows(con, "SELECT species_a, species_b, n_common FROM common_counts"),
        "top": top,
        "inputs": run.graph_inputs(),
        "known_funnel": funnel,
    }
    path = os.path.join(outdir, "overview.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=_json_default)
    return path


def known_funnel(con):
    """
    How far the base species' curated gene-disease pairs (and their diseases) get. Each step keeps the
    pairs of the step before: the gene has an ortholog in a model species of the run, one of those orthologs
    carries a phenotype the ranking can use (only then can any disease's ranking reach the gene), the disease
    has a significant phenolog (the pair is within reach: it can be recovered), and the gene is in the disease's
    candidate ranking (recovered, in-sample)
    """

    return rows(con, """
        WITH pairs AS (SELECT DISTINCT gene_id, disease_id FROM known_associations),
        f AS (SELECT p.disease_id,
                     p.gene_id IN (SELECT human_gene_id FROM human_orthologs) AS has_ortholog,
                     p.gene_id IN (SELECT gene_id FROM genes WHERE in_reach) AS gene_in_reach,
                     p.disease_id IN (SELECT disease_id FROM diseases) AS has_phenolog,
                     (p.disease_id, p.gene_id) IN (SELECT (disease_id, gene_id) FROM candidates WHERE known) AS recovered
              FROM pairs p)
        SELECT count(*) AS curated,
               count(*) FILTER (WHERE has_ortholog) AS has_ortholog,
               count(*) FILTER (WHERE gene_in_reach) AS gene_in_reach,
               count(*) FILTER (WHERE gene_in_reach AND has_phenolog) AS in_reach,
               count(*) FILTER (WHERE recovered) AS recovered,
               count(DISTINCT disease_id) AS diseases,
               count(DISTINCT disease_id) FILTER (WHERE has_ortholog) AS diseases_ortholog,
               count(DISTINCT disease_id) FILTER (WHERE gene_in_reach) AS diseases_gene_in_reach,
               count(DISTINCT disease_id) FILTER (WHERE gene_in_reach AND has_phenolog) AS diseases_in_reach,
               count(DISTINCT disease_id) FILTER (WHERE recovered) AS diseases_recovered
        FROM f""")[0]


def _json_default(value):
    """ numpy scalars from DuckDB/pandas --> plain json """
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(type(value))
