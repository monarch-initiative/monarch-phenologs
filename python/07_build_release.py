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
                                     build_readme,
                                     build_manifest)




if __name__ == '__main__':
    ################
    ## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Builds the release of a finished run: the tables the dashboard reads (which are also \
                                                      the downloads), a README on every table and column, and manifest.json listing \
                                                      every file. Written to <project_dir>/release/<graph_version>/<prediction_network>/, \
                                                      laid out as in the bucket. The versions come from monarch_kg/versions.json. A \
                                                      rebuild replaces it. 08_update_release_list.py adds it to the published list.')

        parser.add_argument("-p", "--project_dir", help="Top most project directory of a finished run", required=True, type=str)
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
              kneighbs=args.kneighbs, rank_metric=args.rank_metric)

    # A rebuild replaces the release, so nothing of an older build is left behind
    outdir = os.path.join(run.project_dir, "release", run.graph_version, run.prediction_network)
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
    build_readme(run, outdir, args.top_phenologs, args.top_families)

    for f in build_manifest(con, outdir):
        rows = "{:,} rows".format(f["rows"]) if f["rows"] is not None else ""
        print("  {:<30} {:>16} {:>8.1f} MB".format(f["path"], rows, f["bytes"] / 1e6))
    print("- Wrote {} ({:.1f} s total)".format(outdir, time.time() - t0))
