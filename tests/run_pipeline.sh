#!/usr/bin/env bash
# Run the python pipeline end to end on a kg slice with tiny parameters (10 random trials, a handful of
# leave one out families). For checking that changes still run, not for producing results.
#
#   tests/run_pipeline.sh <slice_dir> <project_dir> [extra args for 02_kg_data_extraction.py]
#
# slice_dir holds a monarch_kg/ directory (see tests/build_test_kg.py). Outputs land in project_dir.
# PY selects the interpreter (default: python on PATH).
set -euo pipefail

SLICE="$(cd "$1" && pwd)"
PROJECT="$2"
shift 2
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
CODE="$HERE/../python"

TAXON=9606
PRD=disease
N_TRIALS="${N_TRIALS:-10}"
SEED="${SEED:-1}"
CORES="${CORES:-2}"
FDR=0.05
K="${K:-5}"
LOO_MAX="${LOO_MAX:-5}"

mkdir -p "$PROJECT/monarch_kg"
for f in monarch-kg_nodes.tsv monarch-kg_edges.tsv; do
    [[ -e "$PROJECT/monarch_kg/$f" ]] || gunzip -c "$SLICE/monarch_kg/$f.gz" > "$PROJECT/monarch_kg/$f"
done
[[ -e "$PROJECT/monarch_kg/phenio-relation-graph.gz" ]] || ln -s "$SLICE/monarch_kg/phenio-relation-graph.gz" "$PROJECT/monarch_kg/phenio-relation-graph.gz"
PROJECT="$(cd "$PROJECT" && pwd)"

cd "$CODE"
echo "== [01] project layout (downloads skipped, the slice is already in place) =="
$PY 01_download_and_setup.py -p "$PROJECT" -gv 2000-01-01 -pv 2000-01-01

echo "== [02] extraction =="
$PY 02_kg_data_extraction.py -p "$PROJECT" "$@"

echo "== [03] $N_TRIALS random trials =="
$PY 03_compute_phenologs_randomized_trials.py -p "$PROJECT" -n "$N_TRIALS" -c "$CORES" -taxon_id "$TAXON" -prd "$PRD" -seed "$SEED"

echo "== [04] real phenologs =="
$PY 04_compute_phenologs.py -p "$PROJECT" -c "$CORES" -taxon_id "$TAXON" -prd "$PRD"

echo "== [05] fdr table =="
$PY 05_compute_phenologs_fdr_data.py -p "$PROJECT" -taxon_id "$TAXON" -prd "$PRD"

for M in nb hg; do
    echo "== [06] gene rankings k=$K $M =="
    $PY 06_phenologs_compute_ortholog_to_phenotype_rankings.py -p "$PROJECT" -taxon_id "$TAXON" -prd "$PRD" -fdr "$FDR" -kneighbs "$K" -rank_metric "$M" -rank_type gene
    echo "== candidate rankings k=$K $M =="
    $PY phenologs_disease_gene_candidate_rankings.py -p "$PROJECT" -taxon_id "$TAXON" -prd "$PRD" -fdr "$FDR" -kneighbs "$K" -rank_metric "$M"
done

echo "== [06] protein family ranking, mouse only (species ablation path) =="
$PY 06_phenologs_compute_ortholog_to_phenotype_rankings.py -p "$PROJECT" -taxon_id "$TAXON" -prd "$PRD" -fdr "$FDR" -kneighbs "$K" -rank_metric hg -rank_type protein_family -xtaxon_ids 10090

echo "== leave one out ($LOO_MAX families) =="
$PY phenologs_leave_one_out_xvalidation.py -p "$PROJECT" -c "$CORES" -taxon_id "$TAXON" -prd "$PRD" -fdr "$FDR" -max_orthologs "$LOO_MAX"
$PY phenologs_xvalidation_ortholog_to_phenotype_ranking.py -p "$PROJECT" -c "$CORES" -taxon_id "$TAXON" -prd "$PRD" -fdr "$FDR" -kneighbs "$K" -rank_metric hg -rank_type gene

echo "== [07] release =="
$PY 07_build_release.py -p "$PROJECT" -prd "$PRD" -fdr "$FDR" -kneighbs "$K" -c "$CORES"

echo "== [08] release list (a new one, nothing is fetched) =="
$PY 08_update_release_list.py -p "$PROJECT" -r none

echo "== DONE -> $PROJECT =="
