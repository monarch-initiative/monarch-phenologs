# Test slices of the Monarch KG

Small cuts of real Monarch KG releases, enough to run the whole python pipeline in minutes.
Each directory is what `01_download_and_setup.py` would leave under `monarch_kg/`, with the
nodes and edges gzipped (`tests/run_pipeline.sh` unpacks them into a project directory).

| Slice | Source release | Why it is here |
|---|---|---|
| `kg_2026-03-06` | monarch-kg 2026-03-06, phenio 2025-11-14 | the paper's graph: 8 species, `type` column present, Dictyostelium present |
| `kg_2024-03-07` | monarch-kg 2024-03-07 | an older schema: no `type` column, genes with empty names, no DDPHENO in phenio, and a fission yeast with almost no annotations (exercises the species gate) |

Each slice keeps every PANTHER family chosen by `tests/build_test_kg.py` (families with human disease
genes and phenotype annotated members in several species, a few genes per family per species, topped
up until every species has at least ~150 phenotype annotated genes where the release has them), all
the gene to phenotype, gene to disease and orthology edges among those genes, the nodes they touch, and
the `rdfs:subClassOf` closure above every phenotype and disease term for the ontologies the extractor
filters on. `provenance.json` records the parameters and counts.

Rebuild a slice from a full download with:

```
python tests/build_test_kg.py -p path/to/project_dir -o test_data/kg_<date>
```

Run the pipeline over a slice with:

```
PY=python tests/run_pipeline.sh test_data/kg_2026-03-06 /tmp/phenologs_test
python -m pytest tests/ -x -q      # both slices, ten random trials, five leave one out families
```
