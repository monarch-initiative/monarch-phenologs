# Build a small monarch kg slice for testing the pipeline end to end.
#
# Picks PANTHER families that human and other species carry phenotype annotated genes in, takes a capped
# number of genes per family per species, and writes nodes / edges / phenio files holding only what
# those genes touch. Everything is read straight from a full kg download, so the slice carries that
# release's own quirks (missing columns, empty names, terms typed under another category).
#
#   python tests/build_test_kg.py -p path/to/project_dir -o test_data/kg_2026-03-06
#
# project_dir is a pipeline project directory whose monarch_kg/ holds monarch-kg_nodes.tsv,
# monarch-kg_edges.tsv and phenio-relation-graph.gz (what 01_download_and_setup.py produces).
# The slice's nodes and edges are written gzipped (tests/run_pipeline.sh unpacks them into a project dir).

# General imports
import os
import sys
import csv
import gzip
import json
import argparse
from collections import defaultdict, Counter

import pandas as pd

csv.field_size_limit(sys.maxsize)

# Edge categories the extractor consumes, plus the ones we write so the slice is self contained
GENE_PHEN = "biolink:GeneToPhenotypicFeatureAssociation"
GENE_DISEASE = {"biolink:CausalGeneToDiseaseAssociation", "biolink:CorrelatedGeneToDiseaseAssociation"}
ORTHOLOGY = "biolink:GeneToGeneHomologyAssociation"
EDGE_COLS = ["id", "subject", "predicate", "object", "category", "has_evidence", "primary_knowledge_source"]
NODE_COLS = ["id", "category", "name", "in_taxon", "in_taxon_label", "type"]

# Ontologies the extractor filters through phenio (same list as 02_kg_data_extraction.py)
PHENIO_ONTOLOGIES = {"MONDO", "HP", "FYPO", "MP", "DDPHENO"}


def panther_id(has_evidence):
    import re
    m = re.search(r"PANTHER\.FAMILY:(PTHR\d+)", str(has_evidence))
    return m.group(1) if m else ""


def scan_edges(edges_path):
    """
    One pass over the edges file. Returns per gene: its PANTHER family (from its orthology edges), its
    phenotype neighbors and its disease neighbors. Orthology edges are not kept (that is the second pass)
    """

    gene_family = {}
    gene_phens = defaultdict(set)
    gene_diseases = defaultdict(set)
    with open(edges_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        h = {c:i for i, c in enumerate(header)}
        for i, cols in enumerate(reader):
            cat = cols[h["category"]]
            if cat == GENE_PHEN:
                gene_phens[cols[h["subject"]]].add(cols[h["object"]])
            elif cat in GENE_DISEASE:
                gene_diseases[cols[h["subject"]]].add(cols[h["object"]])
            elif cat == ORTHOLOGY:
                fam = panther_id(cols[h["has_evidence"]])
                if fam:
                    gene_family.setdefault(cols[h["subject"]], fam)
                    gene_family.setdefault(cols[h["object"]], fam)
            if i % 5_000_000 == 0 and i > 0:
                print("- {} edges scanned...".format(format(i, ",")))
    return gene_family, gene_phens, gene_diseases


def select_genes(gene_taxon, gene_family, gene_phens, gene_diseases, n_families, genes_per_family, min_per_species):
    """
    Families are ranked by how many species carry a phenotype annotated member, human disease
    annotated members first. Families are added until n_families is reached and every species has at
    least min_per_species phenotype annotated genes in the slice (or the graph runs out)
    """

    human = "NCBITaxon:9606"
    fam_species = defaultdict(lambda: defaultdict(list))   # family -> taxon -> [genes with phenotypes]
    fam_human_disease = set()
    for gene, fam in gene_family.items():
        tx = gene_taxon.get(gene)
        if tx is None:
            continue
        if gene in gene_phens:
            fam_species[fam][tx].append(gene)
        if tx == human and gene in gene_diseases:
            fam_human_disease.add(fam)

    def rank(fam):
        return (fam in fam_human_disease, len(fam_species[fam]), fam)

    families = sorted(fam_species, key=rank, reverse=True)
    taxa = sorted(set(gene_taxon.values()))
    per_species = Counter()
    chosen_fams = []
    for fam in families:
        need_more = len(chosen_fams) < n_families or any(per_species[tx] < min_per_species for tx in taxa
                                                          if any(fam_species[f].get(tx) for f in families))
        if not need_more:
            break
        helps = len(chosen_fams) < n_families or any(per_species[tx] < min_per_species and fam_species[fam].get(tx) for tx in taxa)
        if not helps:
            continue
        chosen_fams.append(fam)
        for tx, genes in fam_species[fam].items():
            per_species[tx] += min(len(genes), genes_per_family)

    genes = set()
    for fam in chosen_fams:
        for tx, members in fam_species[fam].items():
            # Most annotated members first, ties by id, so the pick is deterministic
            members = sorted(members, key=lambda g: (-len(gene_phens[g]), g))
            genes.update(members[:genes_per_family])
        # Human disease genes of the family that carry no phenotype edge still matter for the disease network
        for gene in gene_diseases:
            if gene_family.get(gene) == fam and gene_taxon.get(gene) == human:
                genes.add(gene)

    print("- {} families chosen, {} genes".format(len(chosen_fams), len(genes)))
    for tx, n in sorted(per_species.items()):
        print("  - {} genes with phenotypes: {}".format(tx, n))
    return genes


def write_edges(edges_path, out_path, genes):
    """ Second pass: every consumed edge touching the chosen genes (orthology only between two chosen genes) """

    n_written = 0
    with open(edges_path, newline="") as fh, gzip.open(out_path, "wt", newline="") as out:
        reader = csv.reader(fh, delimiter="\t")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        header = next(reader)
        h = {c:i for i, c in enumerate(header)}
        present = [c for c in EDGE_COLS if c in h]
        writer.writerow(present)
        for cols in reader:
            cat = cols[h["category"]]
            s, o = cols[h["subject"]], cols[h["object"]]
            if cat == ORTHOLOGY:
                keep = s in genes and o in genes
            elif cat == GENE_PHEN or cat in GENE_DISEASE:
                keep = s in genes
            else:
                keep = False
            if keep:
                writer.writerow([cols[h[c]] for c in present])
                n_written += 1
    print("- {} edges written".format(format(n_written, ",")))


def write_nodes(nodes_df, out_path, keep_ids):
    present = [c for c in NODE_COLS if c in nodes_df.columns]
    sub = nodes_df[nodes_df["id"].isin(keep_ids) | nodes_df["id"].str.startswith("NCBITaxon:")][present]
    sub.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print("- {} nodes written ({})".format(format(len(sub), ","), ", ".join(present)))
    return sub


def write_phenio(phenio_path, out_path, terms):
    """
    The rdfs:subClassOf closure above every chosen term, for the ontologies the extractor filters on.
    Enough for the include / exclude roots to resolve exactly as they do on the full file
    """

    parents = defaultdict(set)
    with gzip.open(phenio_path, "rt") as fh:
        for line in fh:
            e1, pred, e2 = line.rstrip("\n").split("\t")
            if pred != "rdfs:subClassOf":
                continue
            if e1.split(":")[0] in PHENIO_ONTOLOGIES and e1.split(":")[0] == e2.split(":")[0]:
                parents[e1].add(e2)

    keep = set()
    stack = [t for t in terms if t.split(":")[0] in PHENIO_ONTOLOGIES]
    while stack:
        t = stack.pop()
        for p in parents.get(t, ()):
            if (t, p) not in keep:
                keep.add((t, p))
                stack.append(p)

    with gzip.open(out_path, "wt") as out:
        for t, p in sorted(keep):
            out.write("{}\trdfs:subClassOf\t{}\n".format(t, p))
    print("- {} phenio subClassOf edges written".format(format(len(keep), ",")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-p", "--project_dir", required=True, help="pipeline project dir holding monarch_kg/")
    ap.add_argument("-o", "--out_dir", required=True, help="where to write the slice (gets a monarch_kg/ subdirectory)")
    ap.add_argument("--families", type=int, default=250, help="PANTHER families to take (more are added until every species reaches --min_per_species)")
    ap.add_argument("--genes_per_family", type=int, default=3, help="phenotype annotated genes kept per family per species")
    ap.add_argument("--min_per_species", type=int, default=150, help="phenotype annotated genes wanted per species, where the graph has them")
    args = ap.parse_args()

    kg = os.path.join(args.project_dir, "monarch_kg")
    nodes_path, edges_path, phenio_path = [os.path.join(kg, f) for f in ("monarch-kg_nodes.tsv", "monarch-kg_edges.tsv", "phenio-relation-graph.gz")]
    out_kg = os.path.join(args.out_dir, "monarch_kg")
    os.makedirs(out_kg, exist_ok=True)

    print("- Reading nodes...")
    nodes_df = pd.read_csv(nodes_path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    gene_taxon = dict(zip(nodes_df.loc[nodes_df["category"] == "biolink:Gene", "id"],
                          nodes_df.loc[nodes_df["category"] == "biolink:Gene", "in_taxon"]))

    print("- Scanning edges (pass 1)...")
    gene_family, gene_phens, gene_diseases = scan_edges(edges_path)

    genes = select_genes(gene_taxon, gene_family, gene_phens, gene_diseases, args.families, args.genes_per_family, args.min_per_species)
    terms = set()
    for g in genes:
        terms |= gene_phens.get(g, set()) | gene_diseases.get(g, set())

    print("- Writing edges (pass 2)...")
    write_edges(edges_path, os.path.join(out_kg, "monarch-kg_edges.tsv.gz"), genes)
    write_nodes(nodes_df, os.path.join(out_kg, "monarch-kg_nodes.tsv.gz"), genes | terms)
    write_phenio(phenio_path, os.path.join(out_kg, "phenio-relation-graph.gz"), terms)

    with open(os.path.join(args.out_dir, "provenance.json"), "w") as fh:
        json.dump({"source_project_dir":os.path.abspath(args.project_dir),
                   "families":args.families, "genes_per_family":args.genes_per_family, "min_per_species":args.min_per_species,
                   "genes":len(genes), "terms":len(terms)}, fh, indent=2)
    print("- Slice written to {}".format(args.out_dir))


if __name__ == "__main__":
    main()
