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
sys.path.insert(0, os.path.join(REPO, "python"))

from species import (SpeciesManifest, INCLUDED_COL, REASON_COL, FLAGS_COL, SOURCE_COL, TOKEN_COL, TAXON_COL,
                     results_dir_name, candidates_dir_name, parse_occurrence_string, split_xspecies_key,
                     parse_trial_filename, parse_all_phenologs_filename)

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


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---- the current graph -------------------------------------------------------------------------

def test_manifest_lists_every_taxon_with_a_verdict(current_run):
    project, _ = current_run
    rows = read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))
    assert {INCLUDED_COL, REASON_COL, FLAGS_COL, SOURCE_COL} <= set(rows[0])
    included = [r for r in rows if r[INCLUDED_COL] == "True"]
    excluded = [r for r in rows if r[INCLUDED_COL] == "False"]
    assert "Homo-sapiens" in {r[TOKEN_COL] for r in included}
    assert all(r[REASON_COL] for r in excluded)
    assert all(not r[REASON_COL] for r in included)
    # Every registered species is named from the registry, never from the graph
    assert all(r[SOURCE_COL] == "registry" for r in rows if r[TOKEN_COL] in ("Homo-sapiens", "Mus-musculus"))


def test_manifest_loader_sees_only_included_species(current_run):
    project, _ = current_run
    manifest = SpeciesManifest(project)
    rows = read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))
    excluded = [r[TOKEN_COL] for r in rows if r[INCLUDED_COL] == "False"]
    assert all(t not in manifest for t in excluded)
    human = manifest.get("9606")
    assert manifest.get("NCBITaxon:9606") is human and manifest.get("Homo-sapiens") is human and manifest.get("Homo sapiens") is human
    assert human.gene_to_disease_path and len(human.common_orthologs_paths) == len(manifest) - 1


def test_outputs_follow_the_naming_conventions(current_run):
    project, _ = current_run
    manifest = SpeciesManifest(project)
    res = os.path.join(project, "phenologs_results", results_dir_name("Homo-sapiens", "disease"))
    assert os.path.isfile(os.path.join(res, "Homo-sapiens_pooled_phenologs_fdr0.05.tsv"))
    for fname in os.listdir(res):
        if fname.endswith("_all_phenologs.tsv"):
            a, b = parse_all_phenologs_filename(fname)
            assert a == "Homo-sapiens" and b in manifest
    trials = os.path.join(project, "random_trials", "Homo_sapiens_disease_trials")
    parsed = [parse_trial_filename(f) for f in os.listdir(trials)]
    assert all(p and p[0] == "Homo-sapiens" and p[1] in manifest for p in parsed)
    for row in read_tsv(os.path.join(res, "Homo-sapiens_fdr_table.tsv")):
        a, b = split_xspecies_key(row["species_comparison"])
        assert a == "Homo-sapiens" and b in manifest


def test_candidate_files_carry_parseable_species_occurrences(current_run):
    project, _ = current_run
    cdir = os.path.join(project, "phenologs_results", results_dir_name("Homo-sapiens", "disease"),
                        candidates_dir_name("Homo-sapiens", "gene", "disease", 0.05, 5, "hg"))
    files = os.listdir(cdir)
    assert len(files) > 0
    manifest = SpeciesManifest(project)
    seen = set()
    for fname in files[:50]:
        for row in read_tsv(os.path.join(cdir, fname)):
            for token in parse_occurrence_string(row["XSpecies_occurence"]):
                assert token in manifest
                seen.add(token)
    assert "Mus-musculus" in seen


def test_species_ablation_writes_a_species_suffixed_matrix(current_run):
    project, _ = current_run
    res = os.path.join(project, "phenologs_results", results_dir_name("Homo-sapiens", "disease"))
    assert os.path.isfile(os.path.join(res, "Homo-sapiens_ortholog_to_disease_Mus-musculus_5kneighbs_hg_0.05.pkl"))


# ---- the 2024 graph: no type column, empty gene names, no DDPHENO, a nearly empty pombe ---------

def test_legacy_graph_runs_and_gates_the_thin_species(legacy_run):
    project, log = legacy_run
    rows = {r[TOKEN_COL]:r for r in read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))}
    assert "Dictyostelium-discoideum" not in rows          # absent from the graph entirely
    pombe = rows["Schizosaccharomyces-pombe"]
    assert pombe[INCLUDED_COL] == "False" and "<" in pombe[REASON_COL]
    assert "genes with empty name" in pombe[FLAGS_COL]
    assert "Warning, DDPHENO parent term" in log
    assert "Schizosaccharomyces-pombe" not in SpeciesManifest(project)


def test_gate_can_be_switched_off(tmp_path):
    project = str(tmp_path / "nogate")
    env = dict(os.environ, PY=sys.executable)
    # Only extraction is needed to see the verdicts. The slice goes in first so 01 skips its downloads
    stage_slice(LEGACY, project)
    subprocess.run([sys.executable, os.path.join(REPO, "python", "01_download_and_setup.py"), "-p", project, "-gv", "2000-01-01", "-pv", "2000-01-01"], check=True, capture_output=True)
    subprocess.run([sys.executable, "02_kg_data_extraction.py", "-p", project, "-min_phenotyped_genes", "0", "-min_ortholog_genes", "0", "-min_phenotypes", "0"],
                   cwd=os.path.join(REPO, "python"), check=True, capture_output=True)
    rows = {r[TOKEN_COL]:r for r in read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))}
    assert rows["Schizosaccharomyces-pombe"][INCLUDED_COL] == "True"


# ---- a taxon the registry does not know --------------------------------------------------------

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

    subprocess.run([sys.executable, "02_kg_data_extraction.py", "-p", project], cwd=os.path.join(REPO, "python"), check=True, capture_output=True)
    rows = {r[TAXON_COL]:r for r in read_tsv(os.path.join(project, "species_data", "species_information_table.tsv"))}
    felis = rows["NCBITaxon:99999"]
    assert felis[TOKEN_COL] == "Felis-catus" and felis[SOURCE_COL] == "graph" and "name from graph" in felis[FLAGS_COL]
    assert felis[INCLUDED_COL] == "True"
    assert os.path.isfile(os.path.join(project, "species_data", "Felis-catus_gene_to_phenotype.tsv"))
    assert SpeciesManifest(project).get("Felis catus").token == "Felis-catus"
