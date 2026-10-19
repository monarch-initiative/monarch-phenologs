# General imports
import os
import time
import shutil
import argparse

# Custom imports
from phenologs_release_utils import (Run,
                                     build_tables,
                                     support_mismatches,
                                     build_additive,
                                     build_overview,
                                     build_downloads,
                                     build_manifest,
                                     add_release)




if __name__ == '__main__':
    ################
    ## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Builds the release of a finished run: the tables the dashboard reads, the run\'s raw \
                                                      results as single tables for download, and manifest.json listing every file. \
                                                      Written to <project_dir>/release/<graph_version>/<prediction_network>/ (laid out \
                                                      as in the bucket) and added to releases.json there. A rebuild replaces it.')

        parser.add_argument("-p", "--project_dir", help="Top most project directory of a finished run", required=True, type=str)
        parser.add_argument("-o", "--outdir", help="Releases folder to write to (default: <project_dir>/release)", required=False, type=str, default=None)
        parser.add_argument("-gv", "--graph_version", help="Monarch KG version of the run (default: a date found in the project path)", required=False, type=str, default=None)
        parser.add_argument("-pv", "--phenio_version", help="Phenio version of the run", required=False, type=str, default=None)
        parser.add_argument("-base", help="Base species name", required=False, type=str, default="Homo-sapiens")
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease", required=False, type=str, default="disease")
        parser.add_argument("-fdr", help="FDR level of the phenologs and rankings the dashboard shows", required=False, type=float, default=0.05)
        parser.add_argument("-kneighbs", help="k of the ranking the dashboard shows", required=False, type=int, default=1000)
        parser.add_argument("-rank_metric", help="hg or nb, the ranking the dashboard shows", required=False, type=str, default="hg")
        parser.add_argument("-top_phenologs", help="Phenologs per disease kept in the additive decomposition", required=False, type=int, default=20)
        parser.add_argument("-top_families", help="Protein families per disease kept in the additive decomposition", required=False, type=int, default=40)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores DuckDB may use (default: all of them)", required=False, type=int, default=None)
        parser.add_argument("-m", "--memory", help="Memory DuckDB may use, e.g. 6GB (default: 80%% of the machine's)", required=False, type=str, default=None)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    run = Run(args.project_dir, base=args.base, prediction_network=args.prediction_network, fdr=args.fdr,
              kneighbs=args.kneighbs, rank_metric=args.rank_metric,
              graph_version=args.graph_version, phenio_version=args.phenio_version)

    # A rebuild replaces the release, so nothing of an older build is left behind
    releases_dir = args.outdir or os.path.join(run.project_dir, "release")
    outdir = os.path.join(releases_dir, run.graph_version, run.prediction_network)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    print("- Building the {} {} release from {}".format(run.graph_version, run.prediction_network, run.project_dir))

    t0 = time.time()
    con = run.connect()
    if args.cpu_cores:
        con.execute("SET threads = {}".format(args.cpu_cores))
    if args.memory:
        con.execute("SET memory_limit = '{}'".format(args.memory))

    counts = build_tables(con, run, outdir)
    pairs, diseases = support_mismatches(con)
    if pairs:
        print("- Note, supporting phenolog counts differ from the pipeline's XSpecies_occurence for {:,} pairs in {:,} "
              "disease(s) (p-value ties at the k boundary, broken differently before phenologs were sorted deterministically)".format(pairs, diseases))
    counts["additive"] = build_additive(con, run, outdir, args.top_phenologs, args.top_families)
    build_overview(con, run, outdir, counts)
    build_downloads(con, run, outdir)

    for f in build_manifest(con, outdir):
        rows = "{:,} rows".format(f["rows"]) if f["rows"] is not None else ""
        print("  {:<40} {:>16} {:>8.1f} MB".format(f["path"], rows, f["bytes"] / 1e6))
    releases = add_release(releases_dir, outdir)
    print("- Wrote {} ({:.1f} s total). releases.json lists {} release(s)".format(outdir, time.time() - t0, len(releases)))
