"""
overview.json: everything the release page shows before any parquet is opened. Headline tiles,
per species stats (the player cards and the tree), the FDR table, common ortholog counts between
species, and the top lists. Depends on the tables built by tables.py.
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

    tiles = rows(con, """
        SELECT (SELECT count(*) FROM diseases)                                   AS diseases,
               (SELECT count(*) FROM phenologs)                                  AS phenologs,
               (SELECT count(DISTINCT (species, phenotype_id)) FROM phenologs)   AS model_phenotypes,
               (SELECT count(*) FROM candidates)                                 AS candidate_pairs,
               (SELECT count(*) FROM candidates WHERE NOT known)                 AS novel_pairs,
               (SELECT count(*) FROM candidates WHERE known)                     AS known_recovered,
               (SELECT count(DISTINCT gene_id) FROM candidates)                  AS candidate_genes,
               (SELECT count(DISTINCT family) FROM candidates)                   AS families_implicated,
               (SELECT count(*) FROM known_associations)                         AS known_pairs,
               (SELECT count(*) FROM species WHERE NOT is_base)                  AS model_species""")[0]

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

    # How strongly each species' phenologs score: -log10 p binned per species, for the overview strip
    pval_bins = rows(con, """
        SELECT species, least(floor(-log10(greatest(pval, 1e-300))), 50)::INTEGER AS neglog10, count(*) AS n
        FROM phenologs GROUP BY species, neglog10 ORDER BY species, neglog10""")

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
        "pval_bins": pval_bins,
    }
    path = os.path.join(outdir, "overview.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=_json_default)
    return path


def _json_default(value):
    """ numpy scalars from DuckDB/pandas --> plain json """
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(type(value))
