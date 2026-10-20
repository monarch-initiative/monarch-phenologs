# Builds a release from a finished pipeline run, for 07_build_release.py. The sections run in the order the
# build calls them: the run's files as DuckDB views, the tables the dashboard reads, the additive decomposition,
# overview.json, the README on every table and column, and manifest.json listing every file

# General imports
import os
import json
import hashlib
import datetime
import subprocess
import duckdb
import pandas as pd

# Custom imports
from phenologs_utils import species_plot_order

# The layout of a release folder. Bump it when a table or file changes, and rebuild the older releases with it
SCHEMA = 2

# A release holds one prediction network, named for the human side of its phenologs: diseases or human phenotypes.
# Tables and columns call them terms in both (term_id, terms.parquet, n_known_terms); these are the words for the README
TERM_WORDS = {"disease": dict(term="disease", terms="diseases", humans="Human diseases", pair="gene-disease", ontology="Mondo"),
              "phenotype": dict(term="human phenotype", terms="human phenotypes", humans="Human phenotypes", pair="gene-phenotype", ontology="HPO")}

# Candidates a release keeps per term, best first (0 keeps every one). A human phenotype's ranking scores thousands
# of genes, so a phenotype release keeps the top of each
TOP_CANDIDATES = {"disease": 0, "phenotype": 100}


########################################################################
### Where a pipeline run keeps its files, and DuckDB views over them ###

# Run knows where a run keeps its files and puts DuckDB views over them. Everything after this section reads those
# views only, so a change in the pipeline's file layout is a change here and nowhere else

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
                 kneighbs=1000, rank_metric="hg", top_candidates=0):
        self.project_dir = os.path.abspath(project_dir)
        self.base = base
        self.prediction_network = prediction_network
        self.words = TERM_WORDS[prediction_network]
        self.fdr = fdr
        self.kneighbs = kneighbs
        self.rank_metric = rank_metric
        self.top_candidates = top_candidates or None        # None keeps every candidate

        # The kg and phenio versions 01_download_and_setup.py settled on
        versions_path = os.path.join(self.project_dir, "monarch_kg", "versions.json")
        if not os.path.isfile(versions_path):
            raise FileNotFoundError("No {}. 01_download_and_setup.py writes it".format(versions_path))
        with open(versions_path) as fh:
            versions = json.load(fh)
        self.graph_version = versions["graph_version"]
        self.phenio_version = versions["phenio_version"]
        self.graph_url = versions["graph_url"]

        self.species_data = os.path.join(self.project_dir, "species_data")
        self.results_dir = os.path.join(self.project_dir, "phenologs_results",
                                        "{}_{}_results".format(base.replace("-", "_"), prediction_network))
        self.pooled = os.path.join(self.results_dir, "{}_pooled_phenologs_fdr{}.tsv".format(base, fdr))
        self.fdr_table = os.path.join(self.results_dir, "{}_fdr_table.tsv".format(base))
        self.candidates_dir = os.path.join(self.results_dir, "{}_gene_to_{}_candidates_fdr{}_{}kneighbs_{}".format(
            base.replace("-", "_"), prediction_network, fdr, kneighbs, rank_metric))
        self.species_table = os.path.join(self.species_data, "species_information_table.tsv")
        self.known_table = os.path.join(self.species_data, "{}_gene_to_{}.tsv".format(base, prediction_network))

        for path in (self.species_table, self.pooled, self.fdr_table, self.candidates_dir, self.known_table):
            if not os.path.exists(path):
                raise FileNotFoundError("Expected pipeline output missing: {}".format(path))

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
            SELECT "Species A Phenotype ID"              AS term_id,
                   "Species A Phenotype Name"            AS term_name,
                   split_part("X Species Comparison", ',', 2) AS species,
                   "Species B Phenotype ID"              AS phenotype_id,
                   "Species B Phenotype Name"            AS phenotype_name,
                   "Overlap Count"::INTEGER              AS overlap,
                   "Ortholog Count A"::INTEGER           AS n_term_orthogroups,
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

        # The candidate ranking the release serves, one file per term, named for its id (MONDO-0000107_<name>,
        # HP-0001249_<name>). After the genes it scored, the pipeline lists the term's known genes no phenolog
        # reached, with no score
        cand_types = "{'Gene':'VARCHAR','Gene_name':'VARCHAR','Association_exists':'BOOLEAN','XSpecies_occurence':'VARCHAR'," \
                     "'p-value':'DOUBLE','XSpecies_phenotypes':'VARCHAR','Protein_family_id':'VARCHAR'}"
        con.execute("""
            CREATE VIEW cand_raw AS
            SELECT replace(regexp_extract(filename, '([A-Za-z]+-[0-9]+)_[^/]*$', 1), '-', ':') AS term_id,
                   Gene AS gene_id, Gene_name AS symbol, Association_exists AS known, XSpecies_occurence AS species_hits,
                   "p-value" AS score, Protein_family_id AS family
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

        # The families each model species shares with the base species, and the phenotype annotations the
        # candidate ranking can use: only genes of those families count, as in the pipeline
        con.execute("""
            CREATE VIEW common_families AS
            SELECT regexp_extract(filename, '_vs_([^/]+)\\.tsv$', 1) AS species, ortholog_id AS family
            FROM {}""".format(read_tsv_sql(os.path.join(sd, "common_orthologs_{}_vs_*.tsv".format(self.base)),
                                           filename="true", all_varchar="true")))
        con.execute("""
            CREATE VIEW ranking_g2p AS
            SELECT DISTINCT g.species, g.gene_id, g.phenotype_id
            FROM g2p g JOIN common_families c ON c.species = g.species AND c.family = g.family""")

        # The base species' curated gene to term edges, its gene_to_disease or gene_to_phenotype table: what known means
        con.execute("""
            CREATE VIEW known AS
            SELECT gene AS gene_id, gene_name AS symbol, predicate, {prd} AS term_id, {prd}_name AS term_name,
                   nullif(ortholog_id, '') AS family
            FROM {table}""".format(prd=self.prediction_network, table=read_tsv_sql(self.known_table, all_varchar="true")))

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
                  UNION ALL SELECT gene_id, symbol FROM known)
            WHERE symbol IS NOT NULL AND symbol <> ''
            GROUP BY gene_id""")
        return con


#####################################################
### The tables the dashboard searches and browses ###

# Each is one SQL statement over the views above (and over tables built before it), materialised in DuckDB and
# then written out as parquet. Counts are 32 bit so they arrive in the browser as plain numbers

# Built in this order. Later tables may refer to earlier ones. Internal ones are steps on the way, not written out
TABLES = ["phenologs", "human_orthologs", "scored", "support", "candidates", "genes", "known_associations", "gene_phenotypes",
          "phenotypes", "families", "terms", "species", "search"]
INTERNAL = {"scored", "support"}

# The two tables meant for download carry their network in the file name, so the four files of a KG version keep
# apart once downloaded: disease-phenologs.parquet, gene-disease-candidates.parquet and the phenotype pair. The site
# reads them by these names (FILE_NAMES in dashboard/js/db.js)
FILE_NAMES = {"phenologs": "{prd}-phenologs", "candidates": "gene-{prd}-candidates"}


def file_name(name, prd):
    """ The file of table name in a release of the prd network """
    return FILE_NAMES.get(name, name).format(prd=prd) + ".parquet"


def table_sql(name, run):
    base = run.base
    kept = "rank <= {}".format(run.top_candidates) if run.top_candidates else "true"
    return {
        # Every significant phenolog of the run, ranked within its term by p-value
        "phenologs": """
            SELECT p.*, row_number() OVER (PARTITION BY term_id ORDER BY pval, species, phenotype_id)::INTEGER AS rank
            FROM pooled p""",

        # Every (term, gene) the ranking scored, ranked within its term by score (ties share a place). The pipeline
        # also appends the term's known genes its phenologs never reached, with no p-value. Those are not candidates
        # and are left out (known_associations still has them)
        "scored": """
            SELECT term_id, gene_id, symbol, family, score, known, species_hits,
                   rank() OVER (PARTITION BY term_id ORDER BY score)::INTEGER AS rank
            FROM cand_raw WHERE score IS NOT NULL""",

        # How a term's k nearest phenologs reach each candidate the release keeps, per model species: the phenologs
        # whose phenotype is annotated to an ortholog of the gene, how many orthologs that takes, and the (phenolog,
        # ortholog) hits the hg score counts (n_hits is the pipeline's XSpecies_occurence). n_used is how many of the
        # term's phenologs with that species the ranking used
        "support": """
            WITH kept AS (SELECT term_id, gene_id FROM scored WHERE {kept}),
            topk AS (SELECT term_id, species, phenotype_id FROM phenologs WHERE rank <= {k}),
            used AS (SELECT term_id, species, count(*) AS n FROM topk GROUP BY term_id, species),
            hits AS (SELECT t.term_id, t.species, t.phenotype_id, g.gene_id AS ortholog_id, k.gene_id
                     FROM kept k
                     JOIN human_orthologs o ON o.human_gene_id = k.gene_id
                     JOIN ranking_g2p g ON g.species = o.species AND g.gene_id = o.gene_id
                     JOIN topk t ON t.term_id = k.term_id AND t.species = g.species AND t.phenotype_id = g.phenotype_id)
            SELECT h.term_id, h.gene_id, h.species,
                   count(DISTINCT h.phenotype_id)::INTEGER AS n_phenologs, any_value(u.n)::INTEGER AS n_used,
                   count(DISTINCT h.ortholog_id)::INTEGER AS n_orthologs, count(*)::INTEGER AS n_hits
            FROM hits h JOIN used u USING (term_id, species)
            GROUP BY h.term_id, h.gene_id, h.species""",

        # The candidates the release keeps: every scored gene of a disease, a human phenotype's best (run.top_candidates,
        # with ties at the last place). One row per (term, gene). support is "species:phenologs:of the term's:orthologs:hits"
        # per species, most phenologs first. global_rank places every kept pair, known or novel, in one list across all
        # terms by score (ties share a rank)
        "candidates": """
            WITH c AS (SELECT term_id, gene_id, symbol, family, score, known, rank FROM scored WHERE {kept}),
            s AS (
                SELECT term_id, gene_id,
                       string_agg(format('{{}}:{{}}:{{}}:{{}}:{{}}', species, n_phenologs, n_used, n_orthologs, n_hits), ','
                                  ORDER BY n_phenologs DESC, species) AS support,
                       count(*) AS n_species, sum(n_phenologs) AS n_phenologs
                FROM support GROUP BY term_id, gene_id),
            names AS (SELECT term_id, any_value(term_name) AS term_name FROM pooled GROUP BY term_id)
            SELECT c.*, coalesce(s.support, '') AS support, coalesce(s.n_species, 0)::INTEGER AS n_species,
                   coalesce(s.n_phenologs, 0)::INTEGER AS n_phenologs, n.term_name,
                   rank() OVER (ORDER BY c.score)::INTEGER AS global_rank
            FROM c LEFT JOIN s USING (term_id, gene_id) LEFT JOIN names n USING (term_id)
            ORDER BY c.score, c.term_id, c.symbol""",

        # Every gene of every species in the run that has a protein family, a phenotype or a curated term.
        # n_ortholog_species counts the other species of the run the gene has an ortholog edge to. in_reach (base
        # species genes only): one of its orthologs carries a phenotype the ranking can use, so a term's
        # ranking can reach it at all
        "genes": """
            WITH base AS (
                SELECT species, gene_id, max(family) AS family
                FROM (SELECT species, gene_id, family FROM g2o UNION ALL SELECT species, gene_id, NULL FROM g2p
                      UNION ALL SELECT '{base}', gene_id, family FROM known)
                GROUP BY species, gene_id),
            reach AS (SELECT DISTINCT o.human_gene_id AS gene_id
                      FROM ranking_g2p g JOIN human_orthologs o ON o.species = g.species AND o.gene_id = g.gene_id),
            os AS (SELECT o.gene_id, count(DISTINCT s.species) AS n
                   FROM (SELECT gene_a AS gene_id, taxon_b AS taxon FROM edges UNION SELECT gene_b, taxon_a FROM edges) o
                   JOIN species_info s ON s.taxon_id = o.taxon GROUP BY o.gene_id),
            np AS (SELECT gene_id, count(DISTINCT phenotype_id) AS n FROM g2p GROUP BY gene_id),
            nk AS (SELECT gene_id, count(DISTINCT term_id) AS n FROM known GROUP BY gene_id),
            pr AS (SELECT gene_id, count(DISTINCT term_id) AS n, count(DISTINCT term_id) FILTER (WHERE NOT known) AS novel, min(rank) AS best_rank
                   FROM candidates GROUP BY gene_id)
            SELECT b.species, b.gene_id, coalesce(n.symbol, b.gene_id) AS symbol, b.family,
                   coalesce(os.n, 0)::INTEGER AS n_ortholog_species,
                   coalesce(np.n, 0)::INTEGER AS n_phenotypes,
                   coalesce(nk.n, 0)::INTEGER AS n_known_terms,
                   coalesce(pr.n, 0)::INTEGER AS n_predicted_terms,
                   coalesce(pr.novel, 0)::INTEGER AS n_novel_terms,
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

        # The base species' curated gene to term associations the predictions are measured against, with the gene's
        # place and score in the term's ranking when the ranking scored it, kept or not
        "known_associations": """
            SELECT DISTINCT k.gene_id, k.symbol, k.term_id, k.term_name, k.predicate, k.family, s.rank, s.score
            FROM known k LEFT JOIN scored s ON s.term_id = k.term_id AND s.gene_id = k.gene_id""",

        # Gene to phenotype annotations of every species (what a phenotype's gene list and a gene's phenotype list show)
        "gene_phenotypes": """
            SELECT DISTINCT species, gene_id, phenotype_id, phenotype_name FROM g2p""",

        # Every annotated phenotype of every species, with how often it turns up as a phenolog.
        # n_known_terms counts the curated terms of the genes annotated with the phenotype. Curated gene to
        # term associations are the base species' only, so model organism phenotypes have none
        "phenotypes": """
            WITH ann AS (
                SELECT g.species, g.phenotype_id, any_value(g.phenotype_name) AS name,
                       count(DISTINCT g.gene_id) AS n_genes, count(DISTINCT g.family) AS n_orthogroups,
                       CASE WHEN g.species <> '{base}' THEN count(DISTINCT c.family) END AS n_shared_orthogroups
                FROM g2p g LEFT JOIN common_families c ON c.species = g.species AND c.family = g.family
                GROUP BY g.species, g.phenotype_id),
            kd AS (
                SELECT p.species, p.phenotype_id, count(DISTINCT d.term_id) AS n
                FROM g2p p JOIN known d USING (gene_id) GROUP BY p.species, p.phenotype_id),
            ph AS (
                SELECT species, phenotype_id, count(*) AS n_phenologs, count(DISTINCT term_id) AS n_terms, min(pval) AS best_pval
                FROM phenologs GROUP BY species, phenotype_id)
            SELECT a.species, a.phenotype_id, a.name, a.n_genes::INTEGER AS n_genes, a.n_orthogroups::INTEGER AS n_orthogroups,
                   a.n_shared_orthogroups::INTEGER AS n_shared_orthogroups,
                   coalesce(kd.n, 0)::INTEGER AS n_known_terms,
                   coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_terms, 0)::INTEGER AS n_terms, ph.best_pval
            FROM ann a LEFT JOIN kd USING (species, phenotype_id) LEFT JOIN ph USING (species, phenotype_id)""",

        # PANTHER families: who is in them, and what they predict. A term the family is ranked for is novel
        # for the family when none of its base species genes is curated for that term
        "families": """
            WITH fam AS (SELECT family, species, gene_id FROM genes WHERE family IS NOT NULL),
            hum AS (SELECT family, string_agg(symbol, ', ' ORDER BY symbol) AS human_symbols, count(*) AS n
                    FROM genes WHERE species = '{base}' AND family IS NOT NULL GROUP BY family),
            kpairs AS (SELECT DISTINCT family, term_id FROM known WHERE family IS NOT NULL),
            pr AS (SELECT c.family, count(DISTINCT c.term_id) AS n, min(c.score) AS best_score,
                          count(DISTINCT c.term_id) FILTER (WHERE k.term_id IS NULL) AS novel
                   FROM candidates c LEFT JOIN kpairs k ON k.family = c.family AND k.term_id = c.term_id
                   GROUP BY c.family),
            kn AS (SELECT family, count(*) AS n FROM kpairs GROUP BY family)
            SELECT f.family, count(*)::INTEGER AS n_genes, count(DISTINCT f.species)::INTEGER AS n_species,
                   coalesce(any_value(hum.n), 0)::INTEGER AS n_base_genes, any_value(hum.human_symbols) AS base_symbols,
                   coalesce(any_value(pr.n), 0)::INTEGER AS n_predicted_terms, coalesce(any_value(pr.novel), 0)::INTEGER AS n_novel_terms,
                   any_value(pr.best_score) AS best_score, coalesce(any_value(kn.n), 0)::INTEGER AS n_known_terms
            FROM fam f LEFT JOIN hum USING (family) LEFT JOIN pr USING (family) LEFT JOIN kn USING (family)
            GROUP BY f.family""",

        # The base species' terms that have at least one phenolog. The candidate counts are over the whole
        # ranking, also where the release keeps only the top of it
        "terms": """
            WITH bysp AS (SELECT term_id, species, count(*) AS n FROM phenologs GROUP BY term_id, species),
            bysp_json AS (              -- ties go to the species first in plot order
                SELECT b.term_id, '{{' || string_agg(format('"{{}}":{{}}', b.species, b.n), ',' ORDER BY b.n DESC, s.plot_order) || '}}' AS by_species,
                       first(b.species ORDER BY b.n DESC, s.plot_order) AS top_species
                FROM bysp b JOIN species_info s USING (species) GROUP BY b.term_id),
            ph AS (SELECT term_id, any_value(term_name) AS name, count(*) AS n_phenologs, count(DISTINCT species) AS n_species, min(pval) AS best_pval
                   FROM phenologs GROUP BY term_id),
            kn AS (SELECT term_id, count(DISTINCT gene_id) AS n FROM known GROUP BY term_id),
            cd AS (SELECT term_id, count(*) AS n, count(*) FILTER (WHERE NOT known) AS novel, count(*) FILTER (WHERE known) AS recovered
                   FROM scored GROUP BY term_id)
            SELECT ph.term_id, ph.name, coalesce(kn.n, 0)::INTEGER AS n_known_genes, ph.n_phenologs::INTEGER AS n_phenologs,
                   ph.n_species::INTEGER AS n_species, ph.best_pval,
                   coalesce(cd.n, 0)::INTEGER AS n_candidates, coalesce(cd.novel, 0)::INTEGER AS n_novel_candidates,
                   coalesce(cd.recovered, 0)::INTEGER AS n_known_recovered, bj.top_species, bj.by_species
            FROM ph LEFT JOIN kn USING (term_id) LEFT JOIN cd USING (term_id) LEFT JOIN bysp_json bj USING (term_id)""",

        # The species of the run with everything the release page and the species pages show. n_phenotypes_annotated
        # counts every phenotype the species' genes carry (the species table's "Unique phenotypes" only counts those
        # on genes with an ortholog)
        "species": """
            WITH ph AS (SELECT species, count(*) AS n_phenologs, count(DISTINCT term_id) AS n_terms, count(DISTINCT phenotype_id) AS n_phenotypes_matched
                        FROM phenologs GROUP BY species),
            cs AS (SELECT s.species, count(*) AS n_candidates, count(*) FILTER (WHERE NOT c.known) AS n_novel
                   FROM support s JOIN candidates c USING (term_id, gene_id) GROUP BY s.species),
            cc AS (SELECT species_b AS species, n_common FROM common_counts WHERE species_a = '{base}'),
            ho AS (SELECT species, count(DISTINCT gene_id) AS n FROM human_orthologs GROUP BY species),
            fd AS (SELECT species, pvalue_cutoff, phenolog_count FROM fdr_table WHERE fdr = {fdr}),
            pa AS (SELECT species, count(DISTINCT phenotype_id) AS n FROM g2p GROUP BY species)
            SELECT i.*, coalesce(pa.n, 0)::INTEGER AS n_phenotypes_annotated,
                   coalesce(ph.n_phenologs, 0)::INTEGER AS n_phenologs, coalesce(ph.n_terms, 0)::INTEGER AS n_terms,
                   coalesce(ph.n_phenotypes_matched, 0)::INTEGER AS n_phenotypes_matched,
                   coalesce(cs.n_candidates, 0)::INTEGER AS n_candidates, coalesce(cs.n_novel, 0)::INTEGER AS n_novel_candidates,
                   cc.n_common::INTEGER AS n_common_with_base, ho.n::INTEGER AS n_human_ortholog_genes, fd.pvalue_cutoff AS fdr_cutoff
            FROM species_info i
            LEFT JOIN ph USING (species) LEFT JOIN cs USING (species) LEFT JOIN cc USING (species)
            LEFT JOIN ho USING (species) LEFT JOIN fd USING (species) LEFT JOIN pa USING (species)""",

        # One row per thing the search box can find. A human phenotype that is a term of the release is found as the term
        "search": """
            SELECT 'term' AS kind, term_id AS id, name AS label, NULL AS species,
                   format('{{:,}} phenologs · {{:,}} known genes', n_phenologs, n_known_genes) AS detail, n_phenologs::INTEGER AS weight
            FROM terms
            UNION ALL
            SELECT 'gene', gene_id, symbol, species, replace(species, '-', ' ') || coalesce(' · ' || family, ''),
                   (n_predicted_terms + 10 * n_known_terms + n_phenotypes)::INTEGER
            FROM genes WHERE n_phenotypes > 0 OR n_known_terms > 0 OR n_predicted_terms > 0
            UNION ALL
            SELECT 'phenotype', phenotype_id, name, species,
                   CASE WHEN species = '{base}' THEN format('{{}} · {{:,}} genes', replace(species, '-', ' '), n_genes)
                        ELSE format('{{}} · {{:,}} phenologs', replace(species, '-', ' '), n_phenologs) END,
                   n_phenologs::INTEGER
            FROM phenotypes WHERE NOT (species = '{base}' AND phenotype_id IN (SELECT term_id FROM terms))
            UNION ALL
            SELECT 'family', family, family, NULL, coalesce(base_symbols, ''), n_predicted_terms::INTEGER
            FROM families
            UNION ALL
            SELECT 'species', species, binomial, species, common, 1000000 FROM species""",
    }[name].format(base=base, fdr=run.fdr, k=run.kneighbs, kept=kept)


def build_tables(con, run, outdir):
    """ Materialise every table in order and write each as parquet. Returns {name: row count} """

    counts = {}
    for name in TABLES:
        con.execute("CREATE TABLE {} AS {}".format(name, table_sql(name, run)))
        if name not in INTERNAL:
            counts[name] = write_parquet(con, run, name, outdir)
    return counts


def support_mismatches(con):
    """
    (pairs, terms) among the kept candidates whose support hits differ from the pipeline's XSpecies_occurence. Runs
    from before phenologs were sorted deterministically can break p-value ties at the k boundary another way
    """

    return con.execute("""
        WITH occ AS (SELECT term_id, gene_id, split_part(item, ':', 1) AS species, split_part(item, ':', 2)::INTEGER AS n_hits
                     FROM (SELECT term_id, gene_id, unnest(string_split(species_hits, ',')) AS item
                           FROM scored SEMI JOIN candidates USING (term_id, gene_id) WHERE species_hits <> ''))
        SELECT count(DISTINCT (term_id, gene_id)), count(DISTINCT term_id)
        FROM occ FULL JOIN support s USING (term_id, gene_id, species)
        WHERE occ.n_hits IS DISTINCT FROM s.n_hits""").fetchone()


def write_parquet(con, run, name, outdir):
    path = os.path.join(outdir, file_name(name, run.prediction_network))
    con.execute("COPY {} TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD)".format(name, path))
    return con.execute("SELECT count(*) FROM {}".format(name)).fetchone()[0]


#############################################################
### The additive decomposition behind a term's candidates ###

# For each term: its top phenologs by p-value; for each phenolog, the protein families carrying that model
# phenotype that the species shares with the base species (the families the phenolog test counts). Each (family,
# phenolog) pair contributes (overlap / term orthogroups) * (1 - p), and families are ranked by their sum. One row
# per (term, family, phenolog) is kept, so the page can re-slice by species or a smaller k. Needs the phenologs,
# genes and known_associations tables

def build_additive(con, run, outdir, top_phenologs=20, top_families=40):
    con.execute("""
        CREATE TABLE additive AS
        WITH top_ph AS (
            SELECT term_id, phenotype_id, phenotype_name, species, overlap, n_term_orthogroups, pval, rank
            FROM phenologs WHERE rank <= {top_phenologs}),
        phen_fam AS (                       -- the protein families carrying each model phenotype, shared with the base species
            SELECT DISTINCT p.species, p.phenotype_id, o.family
            FROM g2p p JOIN g2o o ON o.species = p.species AND o.gene_id = p.gene_id
            JOIN common_families c ON c.species = o.species AND c.family = o.family),
        contrib AS (
            SELECT t.term_id, f.family, t.rank AS phenolog_rank, t.phenotype_id, t.phenotype_name, t.species,
                   (t.overlap::DOUBLE / t.n_term_orthogroups) * (1 - t.pval) AS score
            FROM top_ph t JOIN phen_fam f ON f.species = t.species AND f.phenotype_id = t.phenotype_id),
        fam_tot AS (                        -- summed as decimals, so the totals and the order of tied families are the same every build
            SELECT term_id, family, sum(score::DECIMAL(38, 18))::DOUBLE AS family_total,
                   row_number() OVER (PARTITION BY term_id ORDER BY sum(score::DECIMAL(38, 18)) DESC, family)::INTEGER AS family_rank
            FROM contrib GROUP BY term_id, family),
        base_genes AS (                     -- the base species' genes of each family, split by whether the term already knows them
            SELECT ft.term_id, ft.family,
                   string_agg(g.symbol, ', ' ORDER BY g.symbol) FILTER (WHERE k.gene_id IS NOT NULL) AS known_symbols,
                   string_agg(g.symbol, ', ' ORDER BY g.symbol) FILTER (WHERE k.gene_id IS NULL)     AS other_symbols
            FROM fam_tot ft
            JOIN genes g ON g.family = ft.family AND g.species = '{base}'
            LEFT JOIN (SELECT DISTINCT gene_id, term_id FROM known_associations) k ON k.gene_id = g.gene_id AND k.term_id = ft.term_id
            WHERE ft.family_rank <= {top_families}
            GROUP BY ft.term_id, ft.family)
        SELECT c.term_id, c.family, ft.family_rank, ft.family_total, c.phenolog_rank::INTEGER AS phenolog_rank,
               c.phenotype_id, c.phenotype_name, c.species, c.score,
               bg.known_symbols, bg.other_symbols
        FROM contrib c
        JOIN fam_tot ft USING (term_id, family)
        LEFT JOIN base_genes bg USING (term_id, family)
        WHERE ft.family_rank <= {top_families}
        ORDER BY c.term_id, ft.family_rank, c.phenolog_rank
    """.format(top_phenologs=top_phenologs, top_families=top_families, base=run.base))
    return write_parquet(con, run, "additive", outdir)


#####################################################################################
### overview.json, everything the release page shows before any parquet is opened ###

# Headline tiles, per species stats (the metric cards and the tree), the FDR table, common ortholog counts between
# species, the top lists, and what the graph gave the pipeline (every taxon's extraction counts and how far the
# curated term genes get). Needs the tables above

def rows(con, query):
    cur = con.execute(query)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


def build_overview(con, run, outdir, table_counts):
    base = run.base
    species = rows(con, "SELECT * FROM species ORDER BY plot_order")
    for s in species:
        s["is_base"] = bool(s["is_base"])

    funnel = known_funnel(con)
    tiles = rows(con, """
        SELECT (SELECT count(*) FROM terms)                                      AS terms,
               (SELECT count(*) FROM phenologs)                                  AS phenologs,
               (SELECT count(DISTINCT (species, phenotype_id)) FROM phenologs)   AS model_phenotypes,
               (SELECT count(*) FROM scored)                                     AS scored_pairs,
               (SELECT count(*) FROM candidates)                                 AS candidate_pairs,
               (SELECT count(*) FROM candidates WHERE NOT known)                 AS novel_pairs,
               (SELECT count(DISTINCT gene_id) FROM candidates WHERE NOT known)  AS candidate_genes,
               (SELECT count(DISTINCT family) FROM candidates)                   AS families_implicated,
               (SELECT count(*) FROM species WHERE NOT is_base)                  AS model_species""")[0]
    tiles.update(known_recovered=funnel["recovered"], known_pairs=funnel["curated"], known_in_reach=funnel["in_reach"])

    top = {
        "terms": rows(con, """SELECT term_id AS id, name AS label, n_phenologs AS value, n_species, n_novel_candidates, top_species
                              FROM terms ORDER BY n_phenologs DESC, name LIMIT 15"""),
        # Phenotypes carried by a single gene match every term that gene is tied to, so the list is
        # kept to phenotypes with a few protein families behind them
        "phenotypes": rows(con, """SELECT phenotype_id AS id, name AS label, species, n_phenologs AS value, n_terms, n_orthogroups
                                   FROM phenotypes WHERE n_phenologs > 0 AND n_orthogroups >= 5 ORDER BY n_phenologs DESC, name LIMIT 15"""),
        "genes": rows(con, """SELECT gene_id AS id, symbol AS label, family, n_predicted_terms AS value, n_novel_terms, n_known_terms
                              FROM genes WHERE species = '{}' ORDER BY n_predicted_terms DESC, symbol LIMIT 15""".format(base)),
        "families": rows(con, """SELECT family AS id, family AS label, base_symbols, n_predicted_terms AS value, n_known_terms, n_species
                                 FROM families ORDER BY n_predicted_terms DESC, family LIMIT 15"""),
    }

    out = {
        "release": {"schema": SCHEMA,
                    "graph_version": run.graph_version,
                    "phenio_version": run.phenio_version,
                    "graph_url": run.graph_url,
                    "code_version": code_version(),
                    "built": datetime.date.today().isoformat(),
                    "base_species": base,
                    "prediction_network": run.prediction_network,
                    "fdr": run.fdr, "kneighbs": run.kneighbs, "rank_metric": run.rank_metric,
                    "top_candidates": run.top_candidates,
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
    How far the base species' curated gene-term pairs (and their terms) get. Each step keeps the pairs of the
    step before: the gene has an ortholog in a model species of the run, one of those orthologs carries a
    phenotype the ranking can use (only then can any term's ranking reach the gene), the term has a significant
    phenolog (the pair is within reach: it can be recovered), and the term's ranking scored the gene (recovered,
    in-sample; anywhere in the ranking, also where the release keeps only the top of it)
    """

    return rows(con, """
        WITH pairs AS (SELECT gene_id, term_id, bool_or(rank IS NOT NULL) AS recovered FROM known_associations GROUP BY gene_id, term_id),
        f AS (SELECT p.term_id, p.recovered,
                     p.gene_id IN (SELECT human_gene_id FROM human_orthologs) AS has_ortholog,
                     p.gene_id IN (SELECT gene_id FROM genes WHERE in_reach) AS gene_in_reach,
                     p.term_id IN (SELECT term_id FROM terms) AS has_phenolog
              FROM pairs p)
        SELECT count(*) AS curated,
               count(*) FILTER (WHERE has_ortholog) AS has_ortholog,
               count(*) FILTER (WHERE gene_in_reach) AS gene_in_reach,
               count(*) FILTER (WHERE gene_in_reach AND has_phenolog) AS in_reach,
               count(*) FILTER (WHERE recovered) AS recovered,
               count(DISTINCT term_id) AS terms,
               count(DISTINCT term_id) FILTER (WHERE has_ortholog) AS terms_ortholog,
               count(DISTINCT term_id) FILTER (WHERE gene_in_reach) AS terms_gene_in_reach,
               count(DISTINCT term_id) FILTER (WHERE gene_in_reach AND has_phenolog) AS terms_in_reach,
               count(DISTINCT term_id) FILTER (WHERE recovered) AS terms_recovered
        FROM f""")[0]


def code_version():
    """
    The code a release is built with: the image's date tag and commit, e.g. 2026-10-19-2315 (1a2b3c4), set when the
    image is built. Outside an image, the commit git reports for this checkout
    """

    if os.environ.get("PHENOLOGS_IMAGE"):
        return "{} ({})".format(os.environ["PHENOLOGS_IMAGE"], os.environ.get("PHENOLOGS_COMMIT") or "unknown")
    try:
        return subprocess.run(["git", "describe", "--always", "--dirty", "--exclude=*"], cwd=os.path.dirname(os.path.abspath(__file__)),
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _json_default(value):
    """ numpy scalars from DuckDB/pandas --> plain json """
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(type(value))


#################################################################
### README.md, what every table and column of a release holds ###

# Written next to the tables, so whoever downloads a release has its documentation with it. Every table's columns
# are listed (tests/test_pipeline.py checks that none is missing)

def build_readme(run, outdir, top_phenologs=20, top_families=40):
    """ Write <outdir>/README.md """

    words = dict(run.words, Term=run.words["term"].capitalize(), Terms=run.words["terms"].capitalize())
    kept = KEPT.format(top=run.top_candidates, **words) if run.top_candidates else ""
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write(RELEASE_README.format(graph_version=run.graph_version, phenio_version=run.phenio_version, code_version=code_version(),
                                       fdr=run.fdr, kneighbs=run.kneighbs, metric=METRICS[run.rank_metric][0],
                                       score=METRICS[run.rank_metric][1], top_phenologs=top_phenologs, top_families=top_families,
                                       kept=kept, prd=run.prediction_network, **words))


# Each rank metric's name, and what its score is
METRICS = {"hg": ("hypergeometric", "the hypergeometric p-value of the gene's hits, out of every hit the ranking drew"),
           "nb": ("naive Bayes", "the product, over the gene's hits, of 1 - (the phenolog's overlap / the term's known genes) x "
                                 "(1 - the phenolog's p-value)")}

# Said where a release keeps only the top of each ranking
KEPT = """ This release keeps each {term}'s best {top} candidates, with any ties at the last place; `n_candidates` in
`terms.parquet` counts every gene its ranking scored."""

RELEASE_README = """# {Term} phenologs from Monarch KG {graph_version}

Phenologs and gene candidates computed from Monarch KG {graph_version} (PHENIO {phenio_version}) by monarch-phenologs
{code_version} (https://github.com/monarch-initiative/monarch-phenologs).

{humans} are compared with the phenotypes of each model species through the PANTHER protein families of their
genes. A {term} and a model phenotype whose families overlap more than chance make a phenolog: each species' p-value
cutoff is set from randomised trials at an FDR of {fdr} (`fdr` in `overview.json`). Each {term}'s {kneighbs:,} most
significant phenologs then rank the human genes they point to as candidates for the {term}, by the {metric} score.

The tables call the {terms} terms (`term_id`, `terms.parquet`, `n_known_terms`), the same names a disease release and a
human phenotype release use. Only the phenologs and candidates files are named for the network, so the four files of a
KG version keep apart once downloaded.

| file | one row per |
|---|---|
| `gene-{prd}-candidates.parquet` | ({term}, gene) the ranking scored: the gene candidates |
| `{prd}-phenologs.parquet` | significant phenolog |
| `known_associations.parquet` | curated human {pair} association |
| `gene_phenotypes.parquet` | gene-phenotype annotation, every species |
| `human_orthologs.parquet` | (model species gene, human gene) ortholog pair |
| `additive.parquet` | ({term}, protein family, phenolog), for each {term}'s top {top_phenologs} phenologs and top {top_families} families |
| `terms.parquet`, `genes.parquet`, `phenotypes.parquet`, `families.parquet`, `species.parquet` | {term}, gene, phenotype, protein family, species, with their counts |
| `search.parquet` | thing the site's search box finds |

`overview.json` holds what the site's release page shows; its `release` block has the versions and settings.
`manifest.json` lists every file with its size, row count and sha256.

The tables open in DuckDB, Python and R:

    duckdb -c "SELECT * FROM 'gene-{prd}-candidates.parquet' WHERE NOT known ORDER BY global_rank LIMIT 10"
    python -c "import pandas as pd; print(pd.read_parquet('gene-{prd}-candidates.parquet').head())"
    Rscript -e "print(head(arrow::read_parquet('gene-{prd}-candidates.parquet')))"

{Terms} are {ontology} ids, human genes HGNC ids, protein families PANTHER ids, and species are written `Genus-species`.

## gene-{prd}-candidates.parquet

The gene candidates: every ({term}, gene) the ranking scored, best first.{kept} The pipeline also lists each {term}'s
curated genes its phenologs never reached, without a score. Those are left out here (`known_associations.parquet` has
every curated pair).

| column | |
|---|---|
| `term_id`, `term_name` | the {term} |
| `gene_id`, `symbol` | the human gene |
| `family` | the gene's PANTHER family |
| `score` | smaller is stronger: {score} |
| `known` | the gene is curated for the {term} |
| `rank` | the gene's place in the {term}'s ranking, by score. Ties share a place |
| `global_rank` | the pair's place among every pair in this table, by score. Ties share a place |
| `support` | per model species, most phenologs first, as `species:phenologs:of:orthologs:hits`: the {term}'s phenologs with that species whose phenotype is annotated to an ortholog of the gene, of how many of the {term}'s phenologs with that species the ranking used, the orthologs that takes, and the (phenolog, ortholog) hits the score counts |
| `n_species` | model species in `support` |
| `n_phenologs` | phenologs in `support`, summed over species |

## {prd}-phenologs.parquet

Every significant phenolog, ranked within its {term} by p-value.

| column | |
|---|---|
| `term_id`, `term_name` | the {term} |
| `species` | the model species |
| `phenotype_id`, `phenotype_name` | the model species' phenotype |
| `overlap` | protein families carrying both the {term} and the model phenotype |
| `n_term_orthogroups` | protein families of the {term}'s genes, counting only families the two species share |
| `n_phenotype_orthogroups` | protein families of the model phenotype's genes, counting only families the two species share |
| `n_common` | protein families present in both species: the background of the test |
| `pval` | hypergeometric p-value of the overlap, at or below the species' cutoff |
| `rank` | the phenolog's place among its {term}'s phenologs, by p-value |

## known_associations.parquet

The curated human {pair} associations the candidates are measured against.

| column | |
|---|---|
| `gene_id`, `symbol` | the human gene |
| `term_id`, `term_name` | the {term} |
| `predicate` | the association's Biolink predicate |
| `family` | the gene's PANTHER family |
| `rank`, `score` | the gene's place and score in the {term}'s ranking, empty when the ranking did not score it |

## gene_phenotypes.parquet

| column | |
|---|---|
| `species` | the gene's species |
| `gene_id` | the gene |
| `phenotype_id`, `phenotype_name` | a phenotype annotated to it |

## human_orthologs.parquet

| column | |
|---|---|
| `species` | the model species |
| `gene_id` | the model species gene |
| `human_gene_id`, `human_symbol` | its human ortholog |

## additive.parquet

What each protein family adds to a {term}'s candidates. For each of the {term}'s top {top_phenologs} phenologs: the
families carrying that model phenotype that the species shares with humans (the families the phenolog test counts).
Each (family, phenolog) pair adds (overlap / {term} orthogroups) x (1 - p-value), and families are ranked by their sum.
Each {term} keeps its top {top_families} families.

| column | |
|---|---|
| `term_id` | the {term} |
| `family` | the PANTHER family |
| `family_rank` | the family's place for the {term}, by `family_total` |
| `family_total` | the family's sum over the {term}'s top phenologs |
| `phenolog_rank` | the phenolog's place among the {term}'s phenologs (`rank` in `{prd}-phenologs.parquet`) |
| `phenotype_id`, `phenotype_name`, `species` | the phenolog's model phenotype |
| `score` | what the phenolog adds to the family |
| `known_symbols` | the family's human genes curated for the {term} |
| `other_symbols` | the family's other human genes |

## terms.parquet

The {terms} with at least one phenolog.

| column | |
|---|---|
| `term_id`, `name` | the {term} |
| `n_known_genes` | its curated genes |
| `n_phenologs` | its phenologs |
| `n_species` | model species with one of its phenologs |
| `best_pval` | its lowest phenolog p-value |
| `n_candidates` | genes its ranking scored |
| `n_novel_candidates` | of those, genes not curated for it |
| `n_known_recovered` | of those, genes curated for it |
| `top_species` | the model species with most of its phenologs |
| `by_species` | its phenologs per model species, as JSON: `{{"Mus-musculus":12,"Danio-rerio":3}}` |

## genes.parquet

Every gene of every species with a protein family, a phenotype or a curated {term}.

| column | |
|---|---|
| `species`, `gene_id`, `symbol` | the gene |
| `family` | its PANTHER family |
| `n_ortholog_species` | other species of the release the gene has an ortholog in |
| `n_phenotypes` | phenotypes annotated to it |
| `n_known_terms` | {terms} it is curated for (human genes only) |
| `n_predicted_terms` | {terms} whose candidates here include it |
| `n_novel_terms` | of those, {terms} it is not curated for |
| `best_rank` | its best place in a {term}'s ranking |
| `in_reach` | human genes only: one of its orthologs carries a phenotype the ranking can use, so a {term}'s ranking can reach it |

## phenotypes.parquet

Every phenotype annotated to a gene, in every species.

| column | |
|---|---|
| `species`, `phenotype_id`, `name` | the phenotype |
| `n_genes` | genes annotated with it |
| `n_orthogroups` | protein families of those genes |
| `n_shared_orthogroups` | model species only: of those, the families shared with humans (the ones the phenolog test counts) |
| `n_known_terms` | curated {terms} of the genes annotated with it (so human phenotypes only) |
| `n_phenologs` | phenologs with it |
| `n_terms` | {terms} with such a phenolog |
| `best_pval` | the lowest p-value among them |

## families.parquet

The PANTHER protein families.

| column | |
|---|---|
| `family` | the family |
| `n_genes`, `n_species` | its genes, and the species they are in |
| `n_base_genes`, `base_symbols` | its human genes |
| `n_predicted_terms` | {terms} whose candidates here include one of its genes |
| `n_novel_terms` | of those, {terms} none of its human genes is curated for |
| `best_score` | the best score of one of its genes |
| `n_known_terms` | {terms} one of its human genes is curated for |

## species.parquet

| column | |
|---|---|
| `taxon_id`, `species`, `binomial`, `common` | the species |
| `ontology` | its phenotype ontology |
| `is_base` | the base species (human) |
| `color`, `plot_order`, `divergence_mya` | for the site: its colour, its place in lists, and when its lineage split from the human one, in million years (TimeTree) |
| `n_genes` | its genes in the graph |
| `n_phenotyped_genes` | of those, genes with a phenotype |
| `n_ortholog_genes` | of those, genes with an ortholog |
| `n_orthogroups` | protein families of its genes |
| `n_phenotypes_annotated` | phenotypes annotated to its genes |
| `n_phenologs`, `n_terms`, `n_phenotypes_matched` | its phenologs, the {terms} they are with, and its phenotypes in them |
| `n_candidates`, `n_novel_candidates` | candidate pairs here it supports, and of those the ones not curated |
| `n_common_with_base` | protein families it shares with humans |
| `n_human_ortholog_genes` | its genes with a human ortholog |
| `fdr_cutoff` | its p-value cutoff at FDR {fdr} |

## search.parquet

| column | |
|---|---|
| `kind` | term, gene, phenotype, family or species |
| `id`, `label`, `species` | what it is |
| `detail` | the line the search box shows under it |
| `weight` | orders the results |
"""


#####################
### manifest.json ###

# Every file of a release folder with its size, row count, sha256 and what it holds (DESCRIPTIONS)
DESCRIPTIONS = {
    "overview.json":                 "Everything the release page shows: headline numbers, species, FDR table, top lists, graph inputs",
    "{prd}-phenologs.parquet":       "Every significant phenolog, ranked within its {term} by p-value",
    "gene-{prd}-candidates.parquet": "The {best}gene candidates of each {term}: rank, score, known flag and support per species",
    "additive.parquet":              "Per {term}, each protein family's score from each of its top phenologs",
    "terms.parquet":                 "The {terms} with a phenolog, with their counts",
    "genes.parquet":                 "Every gene with a protein family, a phenotype or a known {term}, with its counts",
    "phenotypes.parquet":            "Every annotated phenotype of every species, with its counts",
    "families.parquet":              "The PANTHER protein families: their genes and the {terms} they are ranked for",
    "species.parquet":               "The species of the run, with their counts",
    "known_associations.parquet":    "The curated human {pair} associations",
    "gene_phenotypes.parquet":       "Gene to phenotype annotations of every species",
    "human_orthologs.parquet":       "Ortholog pairs between model species genes and human genes",
    "search.parquet":                "What the search box finds",
    "README.md":                     "What every table and column holds",
}


def build_manifest(con, run, outdir):
    """ Write <outdir>/manifest.json. Returns its list of files """

    descriptions = {path.format(prd=run.prediction_network): text for path, text in DESCRIPTIONS.items()}
    best = "best {:,} ".format(run.top_candidates) if run.top_candidates else ""
    files = []
    for folder, _, names in sorted(os.walk(outdir)):
        for name in sorted(names):
            if name.startswith(".") or name == "manifest.json":
                continue
            path = os.path.join(folder, name)
            rel = os.path.relpath(path, outdir)
            files.append(dict(path=rel, bytes=os.path.getsize(path), rows=count_rows(con, path), sha256=sha256(path),
                              description=descriptions.get(rel, "").format(best=best, **run.words)))

    with open(os.path.join(outdir, "manifest.json"), "w") as fh:
        json.dump(files, fh, indent=1)
    return files


def count_rows(con, path):
    if path.endswith(".parquet"):
        return con.execute("SELECT count(*) FROM read_parquet('{}')".format(path)).fetchone()[0]
    if path.endswith((".tsv", ".tsv.gz")):
        return con.execute("SELECT count(*) FROM {}".format(read_tsv_sql(path, all_varchar="true"))).fetchone()[0]
    return None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
