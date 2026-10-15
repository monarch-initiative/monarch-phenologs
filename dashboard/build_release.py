"""
Builds everything the phenologs dashboard serves from one pipeline run: a parquet table per searchable
thing (species, diseases, genes, phenotypes, protein families, phenologs, candidates, the additive
decomposition, a search index) and overview.json for the release page.

    python dashboard/build_release.py -p path/to/project_dir -o dashboard/site/data \\
        --graph_version 2026-03-06 --phenio_version 2025-11-14

Then serve the site:  python dashboard/serve.py
"""
import os
import sys
import time
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from release.sources import Run
from release import tables, additive, overview


if __name__ == "__main__":
    def parse_input_command():
        parser = argparse.ArgumentParser(description="Build the dashboard data release from a phenologs project directory")
        parser.add_argument("-p", "--project_dir", help="Top most project directory of a finished pipeline run", required=True, type=str)
        parser.add_argument("-o", "--outdir", help="Where to write the parquet tables and overview.json", required=False, type=str,
                            default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "site", "data"))
        parser.add_argument("--graph_version", help="Monarch KG version of the run (defaults to a date found in the project path)", required=False, type=str, default=None)
        parser.add_argument("--phenio_version", help="Phenio version of the run", required=False, type=str, default=None)
        parser.add_argument("--base", help="Base species name", required=False, type=str, default="Homo-sapiens")
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease", required=False, type=str, default="disease")
        parser.add_argument("-fdr", help="FDR level of the pooled phenologs and candidates to serve", required=False, type=float, default=0.05)
        parser.add_argument("-kneighbs", help="k of the candidate ranking to serve", required=False, type=int, default=1000)
        parser.add_argument("-rank_metric", help="hg or nb", required=False, type=str, default="hg")
        parser.add_argument("--top_phenologs", help="Phenologs per disease kept in the additive decomposition", required=False, type=int, default=20)
        parser.add_argument("--top_families", help="Protein families per disease kept in the additive decomposition", required=False, type=int, default=40)
        return parser.parse_args()

    args = parse_input_command()

    run = Run(args.project_dir, base=args.base, prediction_network=args.prediction_network, fdr=args.fdr,
              kneighbs=args.kneighbs, rank_metric=args.rank_metric,
              graph_version=args.graph_version, phenio_version=args.phenio_version)
    os.makedirs(args.outdir, exist_ok=True)
    print("- Building release for graph {} from {}".format(run.graph_version, run.project_dir))

    t0 = time.time()
    con = run.connect()
    counts = tables.build(con, run, args.outdir)
    pairs, diseases = tables.support_mismatches(con)
    if pairs:
        print("- Note, supporting phenolog counts differ from the pipeline's XSpecies_occurence for {:,} pairs in {:,} "
              "disease(s) (p-value ties at the k boundary, broken differently before phenologs were sorted deterministically)".format(pairs, diseases))
    counts["additive"] = additive.build(con, run, args.outdir, args.top_phenologs, args.top_families)
    for name, n in counts.items():
        print("  {:<20} {:>10,} rows  {:>7.1f} MB".format(name, n, os.path.getsize(os.path.join(args.outdir, name + ".parquet")) / 1e6))

    path = overview.build(con, run, args.outdir, counts)
    print("- Wrote {} ({:.1f} s total)".format(path, time.time() - t0))
