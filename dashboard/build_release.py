"""
Builds one release of the phenologs dashboard from a pipeline run, in its own folder: a parquet table per
searchable thing (species, diseases, genes, phenotypes, protein families, phenologs, candidates, the additive
decomposition, a search index), overview.json for the release page, the run's raw results as single tables
under downloads/, and manifest.json listing every file. releases.json, in the releases folder, lists them all.

    python dashboard/build_release.py -p path/to/project_dir --graph_version 2026-03-06 --phenio_version 2025-11-14

writes dashboard/releases/2026-03-06/disease/ (a rebuild replaces it) and adds it to dashboard/releases/releases.json.
Then serve the site:  python dashboard/serve.py
"""
import os
import sys
import time
import shutil
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from release.sources import Run
from release import tables, additive, overview, downloads, manifest


if __name__ == "__main__":
    def parse_input_command():
        parser = argparse.ArgumentParser(description="Build a dashboard release from a phenologs project directory")
        parser.add_argument("-p", "--project_dir", help="Top most project directory of a finished pipeline run", required=True, type=str)
        parser.add_argument("-o", "--outdir", help="Releases folder. The release is written to <outdir>/<graph_version>/<prediction_network>", required=False, type=str,
                            default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "releases"))
        parser.add_argument("--graph_version", help="Monarch KG version of the run (defaults to a date found in the project path)", required=False, type=str, default=None)
        parser.add_argument("--phenio_version", help="Phenio version of the run", required=False, type=str, default=None)
        parser.add_argument("--base", help="Base species name", required=False, type=str, default="Homo-sapiens")
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease", required=False, type=str, default="disease")
        parser.add_argument("-fdr", help="FDR level of the pooled phenologs and candidates to serve", required=False, type=float, default=0.05)
        parser.add_argument("-kneighbs", help="k of the candidate ranking to serve", required=False, type=int, default=1000)
        parser.add_argument("-rank_metric", help="hg or nb", required=False, type=str, default="hg")
        parser.add_argument("--top_phenologs", help="Phenologs per disease kept in the additive decomposition", required=False, type=int, default=20)
        parser.add_argument("--top_families", help="Protein families per disease kept in the additive decomposition", required=False, type=int, default=40)
        parser.add_argument("--threads", help="Threads DuckDB may use (defaults to every core)", required=False, type=int, default=None)
        parser.add_argument("--memory_limit", help="Memory DuckDB may use, e.g. 6GB (defaults to 80%% of the machine's)", required=False, type=str, default=None)
        return parser.parse_args()

    args = parse_input_command()

    run = Run(args.project_dir, base=args.base, prediction_network=args.prediction_network, fdr=args.fdr,
              kneighbs=args.kneighbs, rank_metric=args.rank_metric,
              graph_version=args.graph_version, phenio_version=args.phenio_version)
    outdir = os.path.join(args.outdir, run.graph_version, run.prediction_network)
    shutil.rmtree(outdir, ignore_errors=True)   # nothing of an older build of this release is left behind
    os.makedirs(outdir)
    print("- Building the {} {} release from {}".format(run.graph_version, run.prediction_network, run.project_dir))

    t0 = time.time()
    con = run.connect()
    if args.threads:
        con.execute("SET threads = {}".format(args.threads))
    if args.memory_limit:
        con.execute("SET memory_limit = '{}'".format(args.memory_limit))
    counts = tables.build(con, run, outdir)
    pairs, diseases = tables.support_mismatches(con)
    if pairs:
        print("- Note, supporting phenolog counts differ from the pipeline's XSpecies_occurence for {:,} pairs in {:,} "
              "disease(s) (p-value ties at the k boundary, broken differently before phenologs were sorted deterministically)".format(pairs, diseases))
    counts["additive"] = additive.build(con, run, outdir, args.top_phenologs, args.top_families)
    overview.build(con, run, outdir, counts)
    downloads.build(con, run, outdir)

    for f in manifest.build(con, outdir):
        rows = "{:,} rows".format(f["rows"]) if f["rows"] is not None else ""
        print("  {:<40} {:>16} {:>8.1f} MB".format(f["path"], rows, f["bytes"] / 1e6))
    releases = manifest.index(args.outdir, outdir)
    print("- Wrote {} ({:.1f} s total). releases.json lists {} release(s)".format(outdir, time.time() - t0, len(releases)))
