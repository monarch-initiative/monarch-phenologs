# Species identity for a phenologs run.
#
# Three things live here and nowhere else:
#   1. The registry: a hand curated list of taxa the monarch kg carries gene nodes for, with the one
#      token each goes by in filenames, data keys and figure labels. The graph fills in any taxon the
#      registry doesn't know, and every gap is written down rather than guessed at silently.
#   2. The naming conventions: how species tokens are combined into filenames, directory names and
#      the comparison keys stored inside result tables, plus the inverse parsers.
#   3. The manifest: species_information_table.tsv as written by 02_kg_data_extraction.py is the
#      record of which species a run has, under what name, and whether each one carries enough data to
#      be used. Everything after extraction asks the manifest, never the graph or the registry.

# General imports
import os
import re
import sys
import pandas as pd


###################################################################################################
### Registry ###

TAXON_PREFIX = "NCBITaxon:"

# Registry order is plot order (the base species of a run is skipped when plotting its comparisons).
# token: letters, digits and hyphens only, as "_" and "," are the separators our filenames and columns use
# phenotype_prefixes: the ontologies we expect that species' phenotype nodes under. Only used to flag
# what extraction actually saw against what was expected, never to select nodes
SPECIES_REGISTRY = [
    {"taxon_id":"NCBITaxon:9606",   "token":"Homo-sapiens",              "binomial":"Homo sapiens",                   "common":"human",          "phenotype_prefixes":["HP"]},
    {"taxon_id":"NCBITaxon:10090",  "token":"Mus-musculus",              "binomial":"Mus musculus",                   "common":"mouse",          "phenotype_prefixes":["MP"]},
    {"taxon_id":"NCBITaxon:7955",   "token":"Danio-rerio",               "binomial":"Danio rerio",                    "common":"zebrafish",      "phenotype_prefixes":["ZP"]},
    {"taxon_id":"NCBITaxon:6239",   "token":"Caenorhabditis-elegans",    "binomial":"Caenorhabditis elegans",         "common":"worm",           "phenotype_prefixes":["WBPhenotype"]},
    {"taxon_id":"NCBITaxon:4896",   "token":"Schizosaccharomyces-pombe", "binomial":"Schizosaccharomyces pombe",      "common":"fission yeast",  "phenotype_prefixes":["FYPO"]},
    {"taxon_id":"NCBITaxon:10116",  "token":"Rattus-norvegicus",         "binomial":"Rattus norvegicus",              "common":"rat",            "phenotype_prefixes":["MP"]},
    {"taxon_id":"NCBITaxon:8364",   "token":"Xenopus-tropicalis",        "binomial":"Xenopus tropicalis",             "common":"frog",           "phenotype_prefixes":["XPO"]},
    {"taxon_id":"NCBITaxon:44689",  "token":"Dictyostelium-discoideum",  "binomial":"Dictyostelium discoideum",       "common":"slime mold",     "phenotype_prefixes":["DDPHENO"]},

    # Carried by the graph with orthology but no usable phenotype annotations so far
    {"taxon_id":"NCBITaxon:7227",   "token":"Drosophila-melanogaster",   "binomial":"Drosophila melanogaster",        "common":"fruit fly",      "phenotype_prefixes":["FBcv"]},
    {"taxon_id":"NCBITaxon:559292", "token":"Saccharomyces-cerevisiae",  "binomial":"Saccharomyces cerevisiae S288C", "common":"budding yeast",  "phenotype_prefixes":["APO"]},
    {"taxon_id":"NCBITaxon:8355",   "token":"Xenopus-laevis",            "binomial":"Xenopus laevis",                 "common":"clawed frog",    "phenotype_prefixes":["XPO"]},
    {"taxon_id":"NCBITaxon:9031",   "token":"Gallus-gallus",             "binomial":"Gallus gallus",                  "common":"chicken",        "phenotype_prefixes":[]},
    {"taxon_id":"NCBITaxon:9615",   "token":"Canis-lupus-familiaris",    "binomial":"Canis lupus familiaris",         "common":"dog",            "phenotype_prefixes":[]},
    {"taxon_id":"NCBITaxon:9823",   "token":"Sus-scrofa",                "binomial":"Sus scrofa",                     "common":"pig",            "phenotype_prefixes":[]},
    {"taxon_id":"NCBITaxon:9913",   "token":"Bos-taurus",                "binomial":"Bos taurus",                     "common":"cow",            "phenotype_prefixes":[]},
    {"taxon_id":"NCBITaxon:227321", "token":"Aspergillus-nidulans",      "binomial":"Aspergillus nidulans FGSC A4",   "common":"",               "phenotype_prefixes":[]},
]

REGISTRY_BY_TAXON = {entry["taxon_id"]:entry for entry in SPECIES_REGISTRY}
REGISTRY_ORDER = {entry["token"]:i for i, entry in enumerate(SPECIES_REGISTRY)}

# What a token is allowed to look like. Enforced at every point a token is minted or read back in
TOKEN_PATTERN = re.compile(r"^[A-Za-z0-9-]+$")


def taxon_key(taxon_id):
    """ 9606, "9606" or "NCBITaxon:9606" --> "NCBITaxon:9606" """

    key = str(taxon_id).strip()
    if not key.startswith(TAXON_PREFIX):
        key = TAXON_PREFIX + key
    return key


def taxon_number(taxon_id):
    """ "NCBITaxon:9606" or 9606 --> "9606" """

    return taxon_key(taxon_id)[len(TAXON_PREFIX):]


def make_token(label):
    """ "Mus musculus" --> "Mus-musculus". Letters, digits and hyphens only. None when nothing usable is left """

    if label is None or pd.isna(label):
        return None
    token = re.sub(r"[^A-Za-z0-9]+", "-", str(label).strip()).strip("-")
    return token or None


def validate_token(token):
    """ A token that can't round trip through our filenames and comparison keys is an error, not a warning """

    if not isinstance(token, str) or not TOKEN_PATTERN.match(token):
        print("- ERROR, species token {!r} must be letters, digits and hyphens only. Exiting...".format(token))
        sys.exit(1)
    return token


def resolve_species(taxon_id, graph_label=""):
    """
    One species record for a taxon. The registry wins, the graph's own label for the taxon fills in a
    taxon the registry doesn't know, and the taxon id itself is the last resort. Every path hands back
    the same fields, and name_source says which one won, so a run can never be half named without
    that being written down in the manifest.
    """

    key = taxon_key(taxon_id)
    if key in REGISTRY_BY_TAXON:
        return dict(REGISTRY_BY_TAXON[key], name_source="registry")

    # A graph label that is just the taxon id again (what in_taxon_label holds in some releases) is no label
    label = "" if (graph_label is None or pd.isna(graph_label)) else str(graph_label).strip()
    token = make_token(label) if (label and not label.startswith(TAXON_PREFIX)) else None
    if token:
        return {"taxon_id":key, "token":token, "binomial":label, "common":"", "phenotype_prefixes":[], "name_source":"graph"}

    return {"taxon_id":key, "token":key.replace(":", "-"), "binomial":"", "common":"", "phenotype_prefixes":[], "name_source":"taxon_id"}


def plot_order(tokens):
    """ Tokens in registry order, then anything the registry doesn't know in the order given """

    return sorted(tokens, key=lambda t: (REGISTRY_ORDER.get(t, len(REGISTRY_ORDER)), tokens.index(t)))


###################################################################################################
### Naming conventions ###

VS = "_vs_"

# Per species tables written by extraction into species_data/, as {token}_{suffix}
SPECIES_FILE_SUFFIXES = ["gene_to_phenotype.tsv",
                         "gene_to_disease.tsv",
                         "gene_to_ortholog.tsv",
                         "phenotype_to_ortholog.pkl",
                         "disease_to_ortholog.pkl"]


def species_filename(token, suffix):
    return "{}_{}".format(token, suffix)


def comparison_stem(species_a, species_b):
    """ Homo-sapiens, Mus-musculus --> Homo-sapiens_vs_Mus-musculus """

    return "{}{}{}".format(species_a, VS, species_b)


def split_comparison_stem(stem):
    """ Homo-sapiens_vs_Mus-musculus --> (Homo-sapiens, Mus-musculus). Tokens never contain "_" so the split is exact """

    parts = stem.split(VS)
    if len(parts) != 2:
        raise ValueError("- {!r} is not a species_a_vs_species_b stem".format(stem))
    return parts[0], parts[1]


def common_orthologs_filename(species_a, species_b):
    return "common_orthologs_{}.tsv".format(comparison_stem(species_a, species_b))


def parse_common_orthologs_filename(fname):
    """ common_orthologs_A_vs_B.tsv --> (A, B), or None for any other file """

    if not (fname.startswith("common_orthologs_") and fname.endswith(".tsv")):
        return None
    return split_comparison_stem(fname[len("common_orthologs_"):-len(".tsv")])


def all_phenologs_filename(species_a, species_b):
    return "{}_all_phenologs.tsv".format(comparison_stem(species_a, species_b))


def parse_all_phenologs_filename(fname):
    """ A_vs_B_all_phenologs.tsv --> (A, B), or None for any other file """

    if not fname.endswith("_all_phenologs.tsv"):
        return None
    return split_comparison_stem(fname[:-len("_all_phenologs.tsv")])


def trial_filename(species_a, species_b, trial_num):
    return "{}_{}.tsv.gz".format(comparison_stem(species_a, species_b), trial_num)


def parse_trial_filename(fname):
    """ A_vs_B_17.tsv.gz (or .tsv) --> (A, B, 17), or None for any other file """

    stem = fname
    for ext in (".tsv.gz", ".tsv"):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
            break
    else:
        return None

    if VS not in stem or "_" not in stem.split(VS)[1]:
        return None
    stem, trial = stem.rsplit("_", 1)
    if not trial.isdigit():
        return None
    species_a, species_b = split_comparison_stem(stem)
    return species_a, species_b, int(trial)


def token_dirname(token):
    """ Homo-sapiens --> Homo_sapiens. Directory names carry the token with underscores, filenames with hyphens """

    return token.replace("-", "_")


def results_dir_name(token, prediction_network):
    """ Homo-sapiens, disease --> Homo_sapiens_disease_results (the phenologs_results/ subdirectory of a base species) """

    return "{}_{}_results".format(token_dirname(token), prediction_network)


def trials_dir_name(token, prediction_network):
    """ Homo-sapiens, disease --> Homo_sapiens_disease_trials (the random_trials/ subdirectory of a base species) """

    return "{}_{}_trials".format(token_dirname(token), prediction_network)


def candidates_dir_name(token, rank_type, prediction_network, fdr, kneighbs, rank_metric):
    """ Homo-sapiens, gene, disease, 0.05, 1000, hg --> Homo_sapiens_gene_to_disease_candidates_fdr0.05_1000kneighbs_hg """

    return "{}_{}_to_{}_candidates_fdr{}_{}kneighbs_{}".format(token_dirname(token), rank_type, prediction_network, fdr, kneighbs, rank_metric)


def xspecies_key(species_a, species_b):
    """ The "X Species Comparison" column value: Homo-sapiens,Mus-musculus """

    return "{},{}".format(species_a, species_b)


def split_xspecies_key(key):
    species_a, species_b = key.split(",")
    return species_a, species_b


def occurrence_string(counts):
    """ {"Mus-musculus":30, "Danio-rerio":5} --> "Mus-musculus:30,Danio-rerio:5" (the XSpecies_occurence column) """

    return ",".join("{}:{}".format(k, v) for k, v in counts.items())


def parse_occurrence_string(value):
    """ "Mus-musculus:30,Danio-rerio:5" --> {"Mus-musculus":30, "Danio-rerio":5}. Empty / missing --> {} """

    if value is None or pd.isna(value) or str(value).strip() == "":
        return {}
    counts = {}
    for piece in str(value).split(","):
        token, count = piece.rsplit(":", 1)
        counts[token] = int(count)
    return counts


###################################################################################################
### Manifest ###

SPECIES_TABLE_FILENAME = "species_information_table.tsv"

# Column names of the manifest. "Taxon label" holds the token and keeps its historical name so
# existing runs and anything reading them keep working
TAXON_COL = "Taxon ID"
TOKEN_COL = "Taxon label"
BINOMIAL_COL = "Binomial"
COMMON_COL = "Common name"
SOURCE_COL = "Name source"
INCLUDED_COL = "Included"
REASON_COL = "Exclusion reason"
PREFIXES_COL = "Phenotype prefixes"
FLAGS_COL = "Flags"


def species_table_path(project_dir):
    return os.path.join(project_dir, "species_data", SPECIES_TABLE_FILENAME)


class RunSpecies:
    """
    One species of a run: its identity from the manifest row, and the species_data/ files it owns.
    species_name is the token, kept under that attribute name as the calculation classes use it
    """

    def __init__(self, row, project_dir):
        self.taxon_id = row[TAXON_COL]
        self.taxon_number = taxon_number(self.taxon_id)
        self.token = validate_token(row[TOKEN_COL])
        self.species_name = self.token
        self.binomial = row.get(BINOMIAL_COL, "")
        self.common = row.get(COMMON_COL, "")
        self.name_source = row.get(SOURCE_COL, "")
        self.included = str(row.get(INCLUDED_COL, True)).lower() == "true"
        self.exclusion_reason = row.get(REASON_COL, "")
        self.phenotype_prefixes = [p for p in str(row.get(PREFIXES_COL, "")).split(";") if p]
        self.flags = [f for f in str(row.get(FLAGS_COL, "")).split(";") if f]
        self.counts = {k:v for k, v in row.items() if k not in (TAXON_COL, TOKEN_COL, BINOMIAL_COL, COMMON_COL,
                                                                 SOURCE_COL, INCLUDED_COL, REASON_COL, PREFIXES_COL, FLAGS_COL)}
        self.project_dir = project_dir

        self.disease_to_ortholog_path = None
        self.gene_to_disease_path = None
        self.gene_to_ortholog_path = None
        self.gene_to_phenotype_path = None
        self.phenotype_to_ortholog_path = None
        self.common_orthologs_paths = []

    def initiate_filepaths(self):
        """ Resolve this species' tables. Disease tables are optional (only the base species tends to have them) """

        species_data_dir = os.path.join(self.project_dir, "species_data")
        for suffix in SPECIES_FILE_SUFFIXES:
            fpath = os.path.join(species_data_dir, species_filename(self.token, suffix))
            if os.path.isfile(fpath):
                setattr(self, "{}_path".format(suffix.split(".")[0]), fpath)
            elif "disease" not in suffix:
                print("- ERROR, Species {} is missing file: {}".format(self.token, fpath))
                sys.exit(1)

        for fname in sorted(os.listdir(species_data_dir)):
            parsed = parse_common_orthologs_filename(fname)
            if parsed and parsed[0] == self.token:
                self.common_orthologs_paths.append(os.path.join(species_data_dir, fname))

    def __repr__(self):
        return "RunSpecies({} {} included={})".format(self.taxon_id, self.token, self.included)


class SpeciesManifest:
    """
    The species of a run, read from species_information_table.tsv. Lookups accept a taxon id with or
    without the NCBITaxon: prefix, a token, or a binomial. Only included species are loaded unless
    asked otherwise, so an excluded species is unreachable downstream by any name
    """

    def __init__(self, project_dir, included_only=True, resolve_files=True):
        self.project_dir = project_dir
        self.path = species_table_path(project_dir)
        if not os.path.isfile(self.path):
            print("- ERROR, Project species information table doesn't seem to exist at {}. Exiting...".format(self.path))
            sys.exit(1)

        df = pd.read_csv(self.path, sep='\t', dtype=str, keep_default_na=False)
        if INCLUDED_COL not in df.columns:
            # A table from before the manifest columns existed: every row it holds was usable
            df[INCLUDED_COL] = "True"

        self.species = []
        self._by_key = {}
        for row in df.to_dict(orient="records"):
            sp = RunSpecies(row, project_dir)
            if included_only and not sp.included:
                continue
            if resolve_files and sp.included:
                sp.initiate_filepaths()
            self.species.append(sp)
            for key in (sp.taxon_id, sp.taxon_number, sp.token, sp.binomial):
                if key:
                    self._by_key[key] = sp

    def __iter__(self):
        return iter(self.species)

    def __len__(self):
        return len(self.species)

    def __contains__(self, key):
        return str(key).strip() in self._by_key

    def get(self, key):
        """ A RunSpecies by taxon id, number, token or binomial. Unknown is an error that names what is known """

        key = str(key).strip()
        if key not in self._by_key:
            print("- ERROR, {!r} is not a species of this run. Known: {}. Exiting...".format(key, ", ".join(self.tokens())))
            sys.exit(1)
        return self._by_key[key]

    def tokens(self):
        return [sp.token for sp in self.species]

    def taxon_ids(self):
        return [sp.taxon_id for sp in self.species]

    def plot_order(self, exclude=()):
        """ Tokens in registry order, minus any given (typically the base species of a comparison) """

        return [t for t in plot_order(self.tokens()) if t not in set(exclude)]

    def resolve_many(self, value):
        """ "9606,10090" or "Homo-sapiens,Mus-musculus" (or a list) --> [RunSpecies, ...]. Empty --> [] """

        if value is None or (isinstance(value, str) and value.strip() == ""):
            return []
        keys = value if isinstance(value, (list, tuple)) else str(value).split(",")
        return [self.get(k) for k in keys if str(k).strip()]
