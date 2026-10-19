"""
manifest.json: every file of a release folder, with its size, row count, sha256 and what it holds.
releases.json: the release block of every published release, with its folder, newest first. The site opens
a release from it.
"""
import os
import json
import hashlib

from .sources import read_tsv_sql

# What each file of a release folder holds
DESCRIPTIONS = {
    "overview.json":              "Everything the release page shows: headline numbers, species, FDR table, top lists, graph inputs",
    "phenologs.parquet":          "Every significant phenolog, ranked within its disease by p-value",
    "candidates.parquet":         "Every (disease, gene) the site's ranking scored: rank, score, known flag and support per species",
    "additive.parquet":           "Per disease, each protein family's score from each of its top phenologs",
    "diseases.parquet":           "The human diseases with a phenolog, with their counts",
    "genes.parquet":              "Every gene with a protein family, a phenotype or a known disease, with its counts",
    "phenotypes.parquet":         "Every annotated phenotype of every species, with its counts",
    "families.parquet":           "The PANTHER protein families: their genes and the diseases they are ranked for",
    "species.parquet":            "The species of the run, with their counts",
    "known_associations.parquet": "The curated human gene-disease associations",
    "gene_phenotypes.parquet":    "Gene to phenotype annotations of every species",
    "human_orthologs.parquet":    "Ortholog pairs between model species genes and human genes",
    "search.parquet":             "What the search box finds",
    "downloads/README.md":                    "What every download and column holds",
    "downloads/candidate_rankings.parquet":   "Every candidate gene ranking the run made (k x rank metric)",
    "downloads/phenologs_all_tested.parquet": "Every disease-phenotype pair tested, with its p-value and whether it is significant",
    "downloads/phenologs_significant.tsv.gz": "The significant phenologs, as a tab separated table",
    "downloads/fdr_table.tsv":                "Each species' p-value cutoff at each FDR level",
}


def build(con, outdir):
    """ Write <outdir>/manifest.json. Returns its list of files """

    files = []
    for folder, _, names in sorted(os.walk(outdir)):
        for name in sorted(names):
            if name.startswith(".") or name == "manifest.json":
                continue
            path = os.path.join(folder, name)
            rel = os.path.relpath(path, outdir)
            files.append(dict(path=rel, bytes=os.path.getsize(path), rows=count_rows(con, path), sha256=sha256(path),
                              description=DESCRIPTIONS.get(rel, "")))

    with open(os.path.join(outdir, "manifest.json"), "w") as fh:
        json.dump(files, fh, indent=1)
    return files


def index(root, folder):
    """
    Add the release in folder to <root>/releases.json, in place of an older build of it. To add it to the
    published list, fetch the bucket's releases.json into root before building. Returns the list of releases
    """

    with open(os.path.join(folder, "overview.json")) as fh:
        release = dict(json.load(fh)["release"], path=os.path.relpath(folder, root) + "/")

    path = os.path.join(root, "releases.json")
    releases = []
    if os.path.isfile(path):
        with open(path) as fh:
            releases = [r for r in json.load(fh) if r["path"] != release["path"]]
    releases = sorted(releases + [release], key=lambda r: r["graph_version"], reverse=True)

    with open(path, "w") as fh:
        json.dump(releases, fh, indent=1)
    return releases


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
