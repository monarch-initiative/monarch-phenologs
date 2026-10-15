"""
Where a pipeline run keeps the things the dashboard needs, and DuckDB views over them.
Everything downstream (tables.py, additive.py, overview.py) reads these views only, so a change
in the pipeline's file layout is a change here and nowhere else.
"""
import os
import re
import sys

import duckdb
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "python"))
from phenologs_utils import SPECIES_NAMES, species_plot_order  # noqa: E402  (the pipeline's species names)


# Display details for species we know. Colors are the paper's Okabe-Ito palette (human lightened for a
# dark page). A species the run has beyond this list gets its binomial, a neutral color and no tree position
SPECIES_DISPLAY = {"NCBITaxon:9606":   dict(common="human",               ontology="HP",          color="#F2F2F2"),
                   "NCBITaxon:10090":  dict(common="mouse",               ontology="MP",          color="#0072B2"),
                   "NCBITaxon:7955":   dict(common="zebrafish",           ontology="ZP",          color="#009E73"),
                   "NCBITaxon:6239":   dict(common="nematode",            ontology="WBPhenotype", color="#E69F00"),
                   "NCBITaxon:4896":   dict(common="fission yeast",       ontology="FYPO",        color="#CC79A7"),
                   "NCBITaxon:10116":  dict(common="rat",                 ontology="MP",          color="#56B4E9"),
                   "NCBITaxon:8364":   dict(common="western clawed frog", ontology="XPO",         color="#D55E00"),
                   "NCBITaxon:44689":  dict(common="slime mold",          ontology="DDPHENO",     color="#B0B0B0"),
                   "NCBITaxon:7227":   dict(common="fruit fly",           ontology="FBcv",        color="#8C564B"),
                   "NCBITaxon:559292": dict(common="budding yeast",       ontology="APO",         color="#E377C2"),
                   "NCBITaxon:8355":   dict(common="African clawed frog", ontology="XPO",         color="#FF7F0E"),
                   "NCBITaxon:9031":   dict(common="chicken",             ontology="",            color="#BCBD22"),
                   "NCBITaxon:9615":   dict(common="dog",                 ontology="",            color="#7F7F7F"),
                   "NCBITaxon:9823":   dict(common="pig",                 ontology="",            color="#7F7F7F"),
                   "NCBITaxon:9913":   dict(common="cattle",              ontology="",            color="#7F7F7F"),
                   "NCBITaxon:227321": dict(common="Aspergillus",         ontology="",            color="#7F7F7F")}

# Approximate divergence from human in millions of years (TimeTree medians), for the species tree.
# Species sharing a value are joined to each other first at the split given in SISTER_SPLITS
DIVERGENCE_MYA = {"NCBITaxon:9606": 0, "NCBITaxon:10090": 87, "NCBITaxon:10116": 87, "NCBITaxon:9615": 94,
                  "NCBITaxon:9823": 94, "NCBITaxon:9913": 94, "NCBITaxon:9031": 319, "NCBITaxon:8364": 352,
                  "NCBITaxon:8355": 352, "NCBITaxon:7955": 429, "NCBITaxon:6239": 686, "NCBITaxon:7227": 686,
                  "NCBITaxon:4896": 1105, "NCBITaxon:559292": 1105, "NCBITaxon:227321": 1105, "NCBITaxon:44689": 1300}
SISTER_SPLITS = {("NCBITaxon:10090", "NCBITaxon:10116"): 20,
                 ("NCBITaxon:8364", "NCBITaxon:8355"): 48,
                 ("NCBITaxon:6239", "NCBITaxon:7227"): 600,
                 ("NCBITaxon:4896", "NCBITaxon:559292"): 580,
                 ("NCBITaxon:9823", "NCBITaxon:9913"): 62}


def read_tsv_sql(path, **options):
    """ read_csv(...) call text for a tab separated file with a header """

    opts = ", ".join(["delim='\\t'", "header=true"] + ["{}={}".format(k, v) for k, v in options.items()])
    return "read_csv('{}', {})".format(path, opts)


class Run:
    """ One pipeline run (a project directory) and the parameter set the dashboard is built from """

    def __init__(self, project_dir, base="Homo-sapiens", prediction_network="disease", fdr=0.05,
                 kneighbs=1000, rank_metric="hg", graph_version=None, phenio_version=None):
        self.project_dir = os.path.abspath(project_dir)
        self.base = base
        self.prediction_network = prediction_network
        self.fdr = fdr
        self.kneighbs = kneighbs
        self.rank_metric = rank_metric
        self.graph_version = graph_version or self._version_from_path()
        self.phenio_version = phenio_version or ""

        self.species_data = os.path.join(self.project_dir, "species_data")
        self.results_dir = os.path.join(self.project_dir, "phenologs_results",
                                        "{}_{}_results".format(base.replace("-", "_"), prediction_network))
        self.pooled = os.path.join(self.results_dir, "{}_pooled_phenologs_fdr{}.tsv".format(base, fdr))
        self.fdr_table = os.path.join(self.results_dir, "{}_fdr_table.tsv".format(base))
        self.candidates_dir = os.path.join(self.results_dir, "{}_gene_to_{}_candidates_fdr{}_{}kneighbs_{}".format(
            base.replace("-", "_"), prediction_network, fdr, kneighbs, rank_metric))
        self.species_table = os.path.join(self.species_data, "species_information_table.tsv")

        for path in (self.species_table, self.pooled, self.fdr_table, self.candidates_dir):
            if not os.path.exists(path):
                raise FileNotFoundError("Expected pipeline output missing: {}".format(path))

    def _version_from_path(self):
        m = re.search(r"(\d{4}-\d{2}-\d{2})", self.project_dir)
        return m.group(1) if m else "unknown"

    # -- species -----------------------------------------------------------------------------------
    def species_rows(self):
        """
        The species this run has, in the pipeline's plot order, with display details and the counts
        extraction tallied. Works for tables with and without the Included column
        """

        df = pd.read_csv(self.species_table, sep="\t")
        if "Included" in df.columns:
            df = df[df["Included"] == True]
        else:
            df = df[df["Genes with >= 1 phenotype edge"] > 0]

        rows = []
        for _, r in df.iterrows():
            taxon, token = r["Taxon ID"], r["Taxon label"]
            disp = SPECIES_DISPLAY.get(taxon, dict(common="", ontology=str(r.get("Phenotype prefix", "")), color="#7F7F7F"))
            rows.append(dict(taxon_id=taxon,
                             species=token,
                             binomial=token.replace("-", " "),
                             common=disp["common"],
                             ontology=disp["ontology"] or str(r.get("Phenotype prefix", "")),
                             color=disp["color"],
                             is_base=(token == self.base),
                             divergence_mya=DIVERGENCE_MYA.get(taxon),
                             n_genes=int(r["Taxon gene count"]),
                             n_phenotyped_genes=int(r["Genes with >= 1 phenotype edge"]),
                             n_phenotypes_annotated=int(r["Unique phenotypes"]),
                             n_ortholog_genes=int(r["Genes with >= 1 ortholog edge"]),
                             n_orthogroups=int(r["Unique ortholog ids (orthogroups)"])))
        order = species_plot_order([r["species"] for r in rows])
        for r in rows:
            r["plot_order"] = order.index(r["species"])
        return sorted(rows, key=lambda r: r["plot_order"])

    def graph_inputs(self):
        """
        Every taxon of the species table, taken part or not, with what extraction found for it: what the graph
        gave the pipeline. Included species first in plot order. Older tables have no Included, Notes or
        phenotype-and-ortholog columns
        """

        df = pd.read_csv(self.species_table, sep="\t", keep_default_na=False)
        rows = []
        for _, r in df.iterrows():
            taxon, token = r["Taxon ID"], r["Taxon label"]
            disp = SPECIES_DISPLAY.get(taxon, dict(common="", color="#7F7F7F"))
            both = r.get("Genes with >= 1 phenotype and ortholog edge", "")
            rows.append(dict(taxon_id=taxon,
                             species=token,
                             binomial=token.replace("-", " "),
                             common=disp["common"],
                             color=disp["color"],
                             included=str(r.get("Included", True)) == "True",
                             notes=str(r.get("Notes", "")),
                             n_genes=int(r["Taxon gene count"]),
                             n_phenotyped_genes=int(r["Genes with >= 1 phenotype edge"]),
                             n_ortholog_genes=int(r["Genes with >= 1 ortholog edge"]),
                             n_both=int(both) if both != "" else None,
                             n_phenotypes=int(r["Unique phenotypes"]),
                             n_orthogroups=int(r["Unique ortholog ids (orthogroups)"]),
                             n_g2p_edges=int(r["Total gene2phenotype edges"])))
        order = species_plot_order([r["species"] for r in rows])
        return sorted(rows, key=lambda r: (not r["included"], order.index(r["species"])))

    def species_tree(self, rows):
        """
        Nested {name, taxon_id, time, children} joining species by divergence from the base species.
        Species without a known divergence are left off (the overview lists them separately)
        """

        placed = {r["taxon_id"]: r for r in rows if r["divergence_mya"] is not None}
        base = next((r for r in rows if r["is_base"]), None)
        if base is None or base["taxon_id"] not in placed:
            return None

        def leaf(r):
            return dict(name=r["species"], taxon_id=r["taxon_id"], time=0)

        root = leaf(base)
        leaves = {base["taxon_id"]: root}
        others = sorted((r for r in placed.values() if not r["is_base"]), key=lambda r: r["divergence_mya"])
        for r in others:
            sister = next((b if a == r["taxon_id"] else a for (a, b), t in SISTER_SPLITS.items()
                           if r["taxon_id"] in (a, b) and (b if a == r["taxon_id"] else a) in leaves), None)
            if sister is not None:
                split = SISTER_SPLITS.get((sister, r["taxon_id"]), SISTER_SPLITS.get((r["taxon_id"], sister)))
                node = leaves[sister]
                # The sister leaf becomes an inner node holding both
                inner = dict(name="", time=split, children=[dict(node), leaf(r)])
                node.clear()
                node.update(inner)
                leaves[r["taxon_id"]] = inner["children"][1]
            else:
                new_leaf = leaf(r)
                root = dict(name="", time=r["divergence_mya"], children=[root, new_leaf])
                leaves[r["taxon_id"]] = new_leaf
        return root

    # -- duckdb ------------------------------------------------------------------------------------
    def connect(self):
        """ A DuckDB connection with a view over every pipeline file the dashboard reads """

        con = duckdb.connect()
        con.register("species_info", pd.DataFrame(self.species_rows()))
        sd = self.species_data

        con.execute("CREATE VIEW species_table AS SELECT * FROM {}".format(read_tsv_sql(self.species_table)))

        con.execute("""
            CREATE VIEW pooled AS
            SELECT "Species A Phenotype ID"              AS disease_id,
                   "Species A Phenotype Name"            AS disease_name,
                   split_part("X Species Comparison", ',', 2) AS species,
                   "Species B Phenotype ID"              AS phenotype_id,
                   "Species B Phenotype Name"            AS phenotype_name,
                   "Overlap Count"::INTEGER              AS overlap,
                   "Ortholog Count A"::INTEGER           AS n_disease_orthogroups,
                   "Ortholog Count B"::INTEGER           AS n_phenotype_orthogroups,
                   "Common Ortholog Count"::INTEGER      AS n_common,
                   hg_pval::DOUBLE                       AS pval
            FROM {}""".format(read_tsv_sql(self.pooled)))

        con.execute("""
            CREATE VIEW fdr_table AS
            SELECT split_part(species_comparison, ',', 2) AS species, fdr::DOUBLE AS fdr, pvalue_cutoff::DOUBLE AS pvalue_cutoff,
                   phenolog_count::INTEGER AS phenolog_count, expected_false_positives::DOUBLE AS expected_false_positives,
                   realized_fdr::DOUBLE AS realized_fdr, n_random_trials::INTEGER AS n_random_trials
            FROM {}""".format(read_tsv_sql(self.fdr_table)))

        cand_types = "{'Gene':'VARCHAR','Gene_name':'VARCHAR','Association_exists':'BOOLEAN','XSpecies_occurence':'VARCHAR'," \
                     "'p-value':'DOUBLE','XSpecies_phenotypes':'VARCHAR','Protein_family_id':'VARCHAR'}"
        con.execute("""
            CREATE VIEW cand_raw AS
            SELECT 'MONDO:' || regexp_extract(filename, 'MONDO-(\\d+)_', 1) AS disease_id,
                   Gene AS gene_id, Gene_name AS symbol, Association_exists AS known,
                   XSpecies_occurence AS species_occurrence, "p-value" AS score, Protein_family_id AS family
            FROM {}""".format(read_tsv_sql(os.path.join(self.candidates_dir, "*.tsv"), filename="true", types=cand_types)))

        con.execute("""
            CREATE VIEW g2p AS
            SELECT regexp_extract(filename, '([^/]+)_gene_to_phenotype\\.tsv$', 1) AS species,
                   gene AS gene_id, gene_name AS symbol, phenotype AS phenotype_id, phenotype_name,
                   nullif(ortholog_id, '') AS family
            FROM {}""".format(read_tsv_sql(os.path.join(sd, "*_gene_to_phenotype.tsv"), filename="true", all_varchar="true")))

        con.execute("""
            CREATE VIEW g2o AS
            SELECT regexp_extract(filename, '([^/]+)_gene_to_ortholog\\.tsv$', 1) AS species, gene AS gene_id, ortholog_id AS family
            FROM {}""".format(read_tsv_sql(os.path.join(sd, "*_gene_to_ortholog.tsv"), filename="true", all_varchar="true")))

        con.execute("""
            CREATE VIEW g2d AS
            SELECT gene AS gene_id, gene_name AS symbol, predicate, disease AS disease_id, disease_name, nullif(ortholog_id, '') AS family
            FROM {}""".format(read_tsv_sql(os.path.join(sd, "{}_gene_to_disease.tsv".format(self.base)), all_varchar="true")))

        con.execute("CREATE VIEW edges AS SELECT * FROM {}".format(
            read_tsv_sql(os.path.join(sd, "ortholog_edges.tsv"), all_varchar="true")))

        con.execute("""
            CREATE VIEW common_counts AS
            SELECT Species_A AS species_a, Species_B AS species_b, Common_Orthologs::INTEGER AS n_common
            FROM {}""".format(read_tsv_sql(os.path.join(sd, "species_common_orthologs_counts.tsv"))))

        # One symbol per gene id, from wherever the run named it. Empty names (older graphs) lose
        con.execute("""
            CREATE VIEW gene_names AS
            SELECT gene_id, any_value(symbol) AS symbol
            FROM (SELECT gene_a AS gene_id, gene_a_name AS symbol FROM edges
                  UNION ALL SELECT gene_b, gene_b_name FROM edges
                  UNION ALL SELECT gene_id, symbol FROM g2p
                  UNION ALL SELECT gene_id, symbol FROM g2d)
            WHERE symbol IS NOT NULL AND symbol <> ''
            GROUP BY gene_id""")
        return con
