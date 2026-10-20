# End to end checks on the kg slices under test_data/. Each slice runs through the python pipeline once (session
# scoped), then the tests read what it wrote. The current slice runs every step of both networks; the legacy one, whose
# differences are in the kg, runs the disease network's steps to its release.
#
#   python -m pytest tests/ -x -q
#
# Small: ten random trials and five leave one out families.

import os
import sys
import csv
import json
import gzip
import shutil
import importlib
import subprocess
import duckdb
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
CODE = os.path.join(REPO, "python")
sys.path.insert(0, CODE)

from phenologs_utils import build_species_info

CURRENT = os.path.join(REPO, "test_data", "kg_2026-03-06")
LEGACY = os.path.join(REPO, "test_data", "kg_2024-03-07")


def run_pipeline(slice_dir, project_dir, *extra, quick=False):
    env = dict(os.environ, PY=sys.executable, QUICK="1" if quick else "0")
    cmd = [os.path.join(HERE, "run_pipeline.sh"), slice_dir, project_dir] + list(extra)
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        raise AssertionError("pipeline failed:\n" + res.stdout[-4000:] + "\n" + res.stderr[-4000:])
    return res.stdout


@pytest.fixture(scope="session")
def current_run(tmp_path_factory):
    project = str(tmp_path_factory.mktemp("current"))
    log = run_pipeline(CURRENT, project)
    return project, log


@pytest.fixture(scope="session")
def legacy_run(tmp_path_factory):
    project = str(tmp_path_factory.mktemp("legacy"))
    log = run_pipeline(LEGACY, project, quick=True)
    return project, log


def stage_slice(slice_dir, project):
    """ Unpack a slice's kg files into project/monarch_kg/ the way run_pipeline.sh does """

    kg = os.path.join(project, "monarch_kg")
    os.makedirs(kg, exist_ok=True)
    for f in ("monarch-kg_nodes.tsv", "monarch-kg_edges.tsv"):
        with gzip.open(os.path.join(slice_dir, "monarch_kg", f + ".gz"), "rb") as src, open(os.path.join(kg, f), "wb") as dst:
            shutil.copyfileobj(src, dst)
    shutil.copy(os.path.join(slice_dir, "monarch_kg", "phenio-relation-graph.gz"), kg)


def extract(project, *extra):
    """ Step 02 alone, on an already staged project """

    subprocess.run([sys.executable, "02_kg_data_extraction.py", "-p", project] + list(extra), cwd=CODE, check=True, capture_output=True)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_table(project):
    """ {species name: row} """

    return {r["Taxon label"]:r for r in read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))}


# ---- the current graph -------------------------------------------------------------------------

def test_species_table_decides_who_takes_part(current_run):
    project, _ = current_run
    rows = species_table(project)
    assert rows["Homo-sapiens"]["Included"] == "True"
    assert all(r["Notes"] for r in rows.values() if r["Included"] == "False")

    # Every later step reads the same table and sees only the included species, by id, number or name
    org_taxon_ids, species = build_species_info(project)
    assert set(org_taxon_ids.values()) == {name for name, r in rows.items() if r["Included"] == "True"}
    assert species["9606"] is species["NCBITaxon:9606"] is species["Homo-sapiens"]
    assert species["9606"].gene_to_disease_path
    assert len(species["9606"].common_orthologs_paths) == len(org_taxon_ids) - 1


def test_outputs_land_where_the_next_step_looks(current_run):
    project, _ = current_run
    names = set(build_species_info(project)[0].values())
    res = os.path.join(project, "phenologs_results", "Homo_sapiens_disease_results")

    assert os.path.isfile(os.path.join(res, "Homo-sapiens_pooled_phenologs_fdr0.05.tsv"))
    pairs = [f.replace("_all_phenologs.tsv", "").split("_vs_") for f in os.listdir(res) if f.endswith("_all_phenologs.tsv")]
    assert pairs and all(a == "Homo-sapiens" and b in names for a, b in pairs)

    trials = os.listdir(os.path.join(project, "random_trials", "Homo_sapiens_disease_trials"))
    assert trials and all(f.startswith("Homo-sapiens_vs_") and f.endswith(".tsv.gz") for f in trials)

    for row in read_tsv(os.path.join(res, "Homo-sapiens_fdr_table.tsv")):
        a, b = row["species_comparison"].split(",")
        assert a == "Homo-sapiens" and b in names

    # Species ablation (protein family ranking restricted to mouse) writes a species suffixed matrix
    assert os.path.isfile(os.path.join(res, "Homo-sapiens_ortholog_to_disease_Mus-musculus_5kneighbs_hg_0.05.pkl"))


def test_candidate_files_name_the_species_behind_each_gene(current_run):
    project, _ = current_run
    names = set(build_species_info(project)[0].values())
    cdir = os.path.join(project, "phenologs_results", "Homo_sapiens_disease_results",
                        "Homo_sapiens_gene_to_disease_candidates_fdr0.05_5kneighbs_hg")
    files = os.listdir(cdir)
    assert len(files) > 0
    seen = set()
    for fname in files[:50]:
        for row in read_tsv(os.path.join(cdir, fname)):
            for item in filter(None, row["XSpecies_occurence"].split(",")):
                species, count = item.split(":")
                assert species in names and int(count) > 0
                seen.add(species)
    assert "Mus-musculus" in seen


def test_release_holds_what_the_run_made(current_run):
    project, _ = current_run
    with open(os.path.join(project, "monarch_kg", "versions.json")) as fh:
        assert json.load(fh) == {"graph_version": "2000-01-01", "phenio_version": "2000-01-01", "graph_url": None}

    # 08 started a new list, so it holds this run's two releases alone, with the versions that made them
    with open(os.path.join(project, "release", "releases.json")) as fh:
        releases = json.load(fh)
    assert [r["path"] for r in releases] == ["2000-01-01/disease/", "2000-01-01/phenotype/"]
    assert all(r["phenio_version"] == "2000-01-01" and r["code_version"] for r in releases)

    # One flat folder per network: the tables, overview.json, the README and the manifest
    for network in ("disease", "phenotype"):
        release = os.path.join(project, "release", "2000-01-01", network)
        with open(os.path.join(release, "manifest.json")) as fh:
            files = {f["path"]: f["rows"] for f in json.load(fh)}
        pooled = read_tsv(os.path.join(project, "phenologs_results", "Homo_sapiens_{}_results".format(network), "Homo-sapiens_pooled_phenologs_fdr0.05.tsv"))
        assert files["{}-phenologs.parquet".format(network)] == len(pooled) > 0
        assert files["gene-{}-candidates.parquet".format(network)] > 0
        assert "README.md" in files and "overview.json" in files and not any("/" in f for f in files)


def test_release_readme_lists_every_column(current_run):
    project, _ = current_run
    for network in ("disease", "phenotype"):
        release = os.path.join(project, "release", "2000-01-01", network)
        with open(os.path.join(release, "README.md")) as fh:
            sections = {s.split("\n", 1)[0]: s for s in fh.read().split("\n## ")[1:]}
        for name in sorted(f for f in os.listdir(release) if f.endswith(".parquet")):
            columns = [c[0] for c in duckdb.sql("DESCRIBE SELECT * FROM '{}'".format(os.path.join(release, name))).fetchall()]
            missing = [c for c in columns if "`{}`".format(c) not in sections.get(name, "")]
            assert not missing, "{} README.md does not describe {} of {}".format(network, missing, name)


def test_phenotype_release_keeps_the_top_of_each_ranking(current_run):
    """ run_pipeline.sh keeps 3 candidates per human phenotype. Recovery still counts the whole ranking """

    project, _ = current_run
    release = os.path.join(project, "release", "2000-01-01", "phenotype")
    sql = lambda q: duckdb.sql(q.replace("@", release + os.sep)).fetchone()
    with open(os.path.join(release, "overview.json")) as fh:
        overview = json.load(fh)
    assert overview["release"]["top_candidates"] == 3

    # Human phenotypes, named for their HP ids, each with its best 3 (more only when tied at the last place)
    assert sql("SELECT count(*) FROM '@terms.parquet' WHERE NOT starts_with(term_id, 'HP:')") == (0,)
    assert sql("SELECT max(rank) FROM '@gene-phenotype-candidates.parquet'") == (3,)
    kept, scored = sql("SELECT count(*), (SELECT sum(n_candidates) FROM '@terms.parquet') FROM '@gene-phenotype-candidates.parquet'")
    assert 0 < kept < scored == overview["tiles"]["scored_pairs"]

    # A known gene keeps its place in the whole ranking, kept or not, and recovery counts all of them
    recovered = sql("SELECT count(DISTINCT (term_id, gene_id)) FROM '@known_associations.parquet' WHERE rank IS NOT NULL")[0]
    assert recovered == overview["known_funnel"]["recovered"] == sql("SELECT sum(n_known_recovered) FROM '@terms.parquet'")[0]
    assert sql("SELECT count(*) FROM '@known_associations.parquet' WHERE rank > 3")[0] > 0


def test_release_list_keeps_the_published_releases(current_run, tmp_path):
    project, _ = current_run
    shutil.copytree(os.path.join(project, "release"), tmp_path / "release")
    published = tmp_path / "published.json"
    published.write_text(json.dumps([dict(path="2001-01-01/disease/", graph_version="2001-01-01", prediction_network="disease"),
                                     dict(path="2000-01-01/disease/", graph_version="2000-01-01", prediction_network="disease", built="earlier"),
                                     dict(path="1999-01-01/disease/", graph_version="1999-01-01", prediction_network="disease")]))
    subprocess.run([sys.executable, "08_update_release_list.py", "-p", str(tmp_path), "-r", str(published)], cwd=CODE, check=True, capture_output=True)

    # The run's builds take the place of the published one, the others stay, newest first
    with open(tmp_path / "release" / "releases.json") as fh:
        releases = json.load(fh)
    assert [r["path"] for r in releases] == ["2001-01-01/disease/", "2000-01-01/disease/", "2000-01-01/phenotype/", "1999-01-01/disease/"]
    assert releases[1]["built"] != "earlier"


def test_kg_versions_read_both_metadata_layouts():
    setup = importlib.import_module("01_download_and_setup")

    # Since 2026-05-09 phenio sits under kg-phenio (whose own version can differ), before that under data
    current = ("id: monarch-kg\nversion: '2026-10-07'\nsources:\n- id: kg-phenio\n  version: v2026-10-07\n  sources:\n"
               "  - id: phenio\n    name: PHENIO\n    urls:\n    - https://github.com/monarch-initiative/phenio\n    version: v2026-10-06\n")
    earlier = "kg-version: 2026-03-06\n\npackages:\n  biolink: 4.3.6\n\ndata:\n  phenio: v2025-11-14\n  alliance: 8.3.0\n"
    assert setup.kg_versions(current) == ("2026-10-07", "2026-10-06")
    assert setup.kg_versions(earlier) == ("2026-03-06", "2025-11-14")


# ---- the 2024 graph: no type column, empty gene names, no DDPHENO, a nearly empty pombe ---------

def test_legacy_graph_runs_and_gates_the_thin_species(legacy_run):
    project, log = legacy_run
    rows = species_table(project)
    assert "Dictyostelium-discoideum" not in rows          # absent from the graph entirely
    pombe = rows["Schizosaccharomyces-pombe"]
    assert pombe["Included"] == "False" and "< 100" in pombe["Notes"]
    assert "genes with empty name" in pombe["Notes"]
    assert "Warning, DDPHENO parent term" in log
    assert "Schizosaccharomyces-pombe" not in build_species_info(project)[0].values()


def test_gate_can_be_switched_off(tmp_path):
    project = str(tmp_path / "nogate")
    stage_slice(LEGACY, project)
    extract(project, "-min_genes", "0")
    assert species_table(project)["Schizosaccharomyces-pombe"]["Included"] == "True"


# ---- a taxon SPECIES_NAMES does not know -------------------------------------------------------

def test_unregistered_taxon_is_named_from_the_graph(tmp_path):
    """
    Clone the slice, relabel a copy of the mouse genes as a made up taxon the graph names "Felis catus",
    and check extraction names it from that label and says so, rather than stopping
    """

    project = str(tmp_path / "felis")
    os.makedirs(os.path.join(project, "monarch_kg"))
    src = os.path.join(CURRENT, "monarch_kg")
    shutil.copy(os.path.join(src, "phenio-relation-graph.gz"), os.path.join(project, "monarch_kg"))

    with gzip.open(os.path.join(src, "monarch-kg_nodes.tsv.gz"), "rt") as fh, open(os.path.join(project, "monarch_kg", "monarch-kg_nodes.tsv"), "w") as out:
        header = fh.readline()
        cols = header.rstrip("\n").split("\t")
        out.write(header)
        for line in fh:
            out.write(line)
            row = dict(zip(cols, line.rstrip("\n").split("\t")))
            if row["in_taxon"] == "NCBITaxon:10090":
                row["id"] = "FELIS:" + row["id"].split(":")[1]
                row["in_taxon"] = "NCBITaxon:99999"
                row["in_taxon_label"] = "NCBITaxon:99999"      # the gene level label is no help, as in the 2026 graph
                out.write("\t".join(row[c] for c in cols) + "\n")
        taxon_row = {c:"" for c in cols}
        taxon_row.update({"id":"NCBITaxon:99999", "category":"biolink:OrganismTaxon", "name":"Felis catus"})
        out.write("\t".join(taxon_row[c] for c in cols) + "\n")

    with gzip.open(os.path.join(src, "monarch-kg_edges.tsv.gz"), "rt") as fh, open(os.path.join(project, "monarch_kg", "monarch-kg_edges.tsv"), "w") as out:
        header = fh.readline()
        cols = header.rstrip("\n").split("\t")
        out.write(header)
        for line in fh:
            out.write(line)
            row = dict(zip(cols, line.rstrip("\n").split("\t")))
            if row["subject"].startswith("MGI:") or row["object"].startswith("MGI:"):
                row["subject"] = row["subject"].replace("MGI:", "FELIS:")
                row["object"] = row["object"].replace("MGI:", "FELIS:")
                row["id"] = row["id"] + "-felis"
                out.write("\t".join(row[c] for c in cols) + "\n")

    extract(project)
    rows = {r["Taxon ID"]:r for r in species_table(project).values()}
    felis = rows["NCBITaxon:99999"]
    assert felis["Taxon label"] == "Felis-catus" and felis["Included"] == "True"
    assert "not in SPECIES_NAMES" in felis["Notes"]
    assert os.path.isfile(os.path.join(project, "species_data", "Felis-catus_gene_to_phenotype.tsv"))
    assert "Felis-catus" in build_species_info(project)[1]


# ---- random trials with a seed are reproducible ------------------------------------------------

def test_random_trials_repeat_with_a_seed(current_run, tmp_path):
    """ Two runs with the same seed write the same trial files; another seed writes different ones """

    project, _ = current_run

    def trials(seed, where):
        where = str(tmp_path / where)
        shutil.copytree(os.path.join(project, "species_data"), os.path.join(where, "species_data"))
        subprocess.run([sys.executable, "03_compute_phenologs_randomized_trials.py", "-p", where, "-n", "2", "-c", "2",
                        "-taxon_id", "9606", "-prd", "disease", "-seed", str(seed)], cwd=CODE, check=True, capture_output=True)
        trial_dir = os.path.join(where, "random_trials", "Homo_sapiens_disease_trials")
        out = {}
        for f in sorted(os.listdir(trial_dir)):
            with gzip.open(os.path.join(trial_dir, f), "rb") as fh:
                out[f] = fh.read()
        return out

    a, b, c = trials(7, "a"), trials(7, "b"), trials(8, "c")
    assert len(a) > 0 and a.keys() == b.keys() == c.keys()
    assert a == b
    assert any(a[f] != c[f] for f in a)
