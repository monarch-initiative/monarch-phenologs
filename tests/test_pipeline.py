# End to end checks on the kg slices under test_data/. Each slice runs through the whole python
# pipeline once (session scoped), then the tests read what it wrote.
#
#   python -m pytest tests/ -x -q
#
# Slow but small: ten random trials and five leave one out families per slice.

import os
import sys
import csv
import gzip
import shutil
import subprocess
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
CODE = os.path.join(REPO, "python")
sys.path.insert(0, CODE)

from phenologs_utils import build_species_info

CURRENT = os.path.join(REPO, "test_data", "kg_2026-03-06")
LEGACY = os.path.join(REPO, "test_data", "kg_2024-03-07")


def run_pipeline(slice_dir, project_dir, *extra):
    env = dict(os.environ, PY=sys.executable)
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
    log = run_pipeline(LEGACY, project)
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
