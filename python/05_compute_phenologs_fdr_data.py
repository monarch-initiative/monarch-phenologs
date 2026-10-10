# General imports
import os
import sys
import argparse
import pickle
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from typing import Dict, List

# Custom imports
from phenologs_utils import (FDR_LEVELS,
                             DEFAULT_FDR,
                             TRIAL_PARAM_COLS,
                             PHENOLOG_PARAM_COLS,
                             bulk_compute_hyper_geom)
from species import (SpeciesManifest,
                     parse_trial_filename,
                     all_phenologs_filename,
                     results_dir_name,
                     trials_dir_name,
                     xspecies_key)




def gather_trial_data(input_dir):
    """
    Grabs all files within the input directory where
    - "_vs_" is contained in the filename,
    - the filename ends with ".tsv" or ".tsv.gz"
    - the last charcter of the filename before ".tsv" must pertain to a trial number (i.e. an integer)

    {(species_a, species_b):{trial_number:path/to/trial/results/file}} is the format we store all our our data
    """

    # Gather all filepaths to trial data
    trials_by_species = {}
    for fname in sorted(os.listdir(input_dir)):

        # Anything that isn't a species_a_vs_species_b_<trial>.tsv(.gz) file is skipped
        parsed = parse_trial_filename(fname)
        if parsed == None:
            continue
        specA, specB, trial_num = parsed

        # Keep track of trial information by species-species comparison
        key = tuple((specA, specB))
        if trials_by_species.get(key, None) == None:
            trials_by_species.update({key:{}})

        if trials_by_species[key].get(trial_num, None) == None:
            trials_by_species[key].update({trial_num:os.path.join(input_dir, fname)})

    print("- Total filepaths gathered {}".format(format(sum([len(v) for v in trials_by_species.values()]), ',')))
    print("- Total pairwise species comparison types found {}".format(len(trials_by_species)))
    for k, v in trials_by_species.items():
        print("- {} -- Trials found {}".format(k, format(len(v), ',')))

    return trials_by_species


def collapse_cross_species_trials(trial_paths: Dict, common_ortholog_count: int):
    """
    Keeps track of how many times each (c,N,m,n) hypergeometric parameter set is found across all of
    our trials. 10,000 Human phenotypes x 10,000 mouse phenotypes = 100,000,000 trial phenolog
    calculations. We do 1,000 trials so we would have 100 Billion data points. Instead, we collapse
    our data by unique parameter set and count the occurence of them across all trials. By doing this,
    we end up with ~56,000 key, value pairs that we are able to work with much more effectively.

    Parameter sets are the key here rather than the hg_pval column on purpose. Pandas reads floats
    back with a parser that isn't correctly rounded, so a pvalue written by step 03 and the same
    pvalue written by step 04 can come back differing in the last bit. That splits a single discrete
    pvalue in two and leaves the random and the real distribution unable to be lined up against
    each other. Recomputing from the integer parameters sidesteps the round trip entirely.
    """

    param_counts = Counter()
    for trial_num, trial_path in sorted(trial_paths.items()):

        df = pd.read_csv(trial_path, sep='\t') # Compression is inferred from the file extension

        # N is a property of the species pair, so a trial file disagreeing with the results file is
        # either from a different run, or predates the trial file parameter column labelling fix
        trial_world = set(df[TRIAL_PARAM_COLS[1]].unique().tolist())
        if trial_world != {common_ortholog_count}:
            print("- ERROR, {} holds common ortholog count(s) {} where {} was expected.".format(os.path.basename(trial_path),
                                                                                                trial_world,
                                                                                                common_ortholog_count))
            print("- Random trials written before the trial file column labelling fix must be regenerated (step 03). Exiting...")
            sys.exit(1)

        # Tally up parameter set occurrences across this trial
        for params, occ in zip(map(tuple, df[TRIAL_PARAM_COLS].to_numpy().tolist()),
                               df["occurrence"].to_numpy().tolist()):
            param_counts[params] += occ

    print("- Random trial comparisons {} collapsed to {} unique hg parameter sets".format(format(sum(param_counts.values()), ','),
                                                                                         format(len(param_counts), ',')))
    return param_counts


def collapse_cross_species_phenologs(input_dir, comp_key):
    """
    The same parameter set tally that collapse_cross_species_trials produces, for the one real
    (non permuted) comparison. Also hands back the common ortholog count N, which the random trials
    for this species pair are validated against.
    """

    # Read in data to memory
    phenologs_filepath = os.path.join(input_dir, all_phenologs_filename(comp_key[0], comp_key[1]))
    if not os.path.isfile(phenologs_filepath):
        print("- ERROR, No real phenologs results found at {}. Exiting...".format(phenologs_filepath))
        sys.exit(1)

    phen_df = pd.read_csv(phenologs_filepath, sep='\t')
    world_counts = set(phen_df[PHENOLOG_PARAM_COLS[1]].unique().tolist())
    if len(world_counts) != 1:
        print("- ERROR, {} holds more than one common ortholog count {}. Exiting...".format(os.path.basename(phenologs_filepath),
                                                                                           world_counts))
        sys.exit(1)

    # Grouped rather than tallied row by row, as this file runs to millions of comparisons
    param_counts = {tuple(int(p) for p in params):int(occ)
                    for params, occ in phen_df.groupby(PHENOLOG_PARAM_COLS).size().items()}

    # Every phenolog has to end up in exactly one parameter set, otherwise D(t) is understated and
    # every fdr computed from it is optimistic
    if sum(param_counts.values()) != len(phen_df):
        print("- ERROR, {} comparisons tallied from the {} rows of {}. Exiting...".format(format(sum(param_counts.values()), ','),
                                                                                         format(len(phen_df), ','),
                                                                                         os.path.basename(phenologs_filepath)))
        sys.exit(1)

    print("- Real phenolog comparisons {} collapsed to {} unique hg parameter sets".format(format(sum(param_counts.values()), ','),
                                                                                          format(len(param_counts), ',')))
    return param_counts, world_counts.pop()


def param_counts_to_pval_counts(param_counts: Dict, hg_pvals: Dict):
    """
    Rolls a {(c,N,m,n):occurrence} tally up into the {pvalue:occurrence} distribution that the fdr
    calculation and the random vs real figure both work from. Different parameter sets routinely land
    on the same pvalue, so this is a many to one collapse.
    """

    pval_counts = Counter()
    for params, occ in param_counts.items():
        pval_counts[hg_pvals[params]] += occ

    return pval_counts


# For visualization
def collapsed_data_to_cumulitive(pval_counts: Dict):
    """
    [[pvalue, occurrence], ...] ordered from least to most significant, alongside the fraction of the
    distribution that is more significant than each pvalue (what the random vs real figure plots)
    """

    total = float(sum(pval_counts.values()))
    sorted_data = np.asarray(sorted([[k,v] for k,v in pval_counts.items()], key=lambda x: x[0], reverse=True))
    data_frac = 1. - (np.cumsum(sorted_data.T[1]) / total)

    return sorted_data, data_frac, total


def get_fdr_cuttoffs(random_pval_counts: Dict, real_pval_counts: Dict, n_trials: int, fdr_levels: List[float]):
    """
    False discovery rate as laid out in McGary et. al. 2010. For a candidate pvalue threshold t,
        P(t)   = permuted comparisons at or below t, pooled over all n_trials trials
        D(t)   = real comparisons at or below t (the phenologs we would be calling significant)
        fdr(t) = (P(t) / n_trials) / D(t)

    Only the pvalues the real data actually achieves are candidate thresholds, because a threshold
    placed between two of them takes on false positives without adding a single discovery. fdr(t) is
    not monotonic in t (both counts grow, and D grows in discrete jumps), so every candidate is
    scored and the most permissive threshold still within each requested level is taken. That is the
    threshold giving the most phenologs at that level.

    Both distributions leave out the zero overlap comparisons, which is consistent and doesn't skew
    the ratio. A zero overlap has a pvalue of 1.0, so it can never be a discovery nor a false
    positive at any threshold we would consider.
    """

    if len(real_pval_counts) == 0 or len(random_pval_counts) == 0:
        print("- ERROR, An empty real or random pvalue distribution can't be assigned an fdr. Exiting...")
        sys.exit(1)

    # Candidate thresholds least to most permissive, and the discovery count at each
    candidates, real_counts = zip(*sorted(real_pval_counts.items()))
    candidates, discoveries = np.asarray(candidates), np.cumsum(real_counts)

    # Expected false positives per trial at each candidate threshold
    random_pvals, random_counts = zip(*sorted(random_pval_counts.items()))
    random_cumul = np.cumsum(random_counts)
    inds = np.searchsorted(np.asarray(random_pvals), candidates, side='right')
    false_positives = np.where(inds > 0, random_cumul[np.maximum(inds - 1, 0)], 0.) / float(n_trials)

    fdr_curve = false_positives / discoveries

    # Most permissive threshold that still holds for each requested level
    fdr_cutoffs = {}
    for level in fdr_levels:

        # No threshold reaches this level. A 0.0 cutoff lets nothing through, which is the only safe
        # reading of "unreachable" (the old code defaulted to 1.0, which passes everything)
        passing = np.flatnonzero(fdr_curve <= level)
        if len(passing) == 0:
            fdr_cutoffs[level] = {"pvalue_cutoff":0.,
                                  "phenolog_count":0,
                                  "expected_false_positives":0.,
                                  "realized_fdr":0.}
            continue

        # Candidates run least to most permissive, so the last one passing is the one yielding the
        # most phenologs at this level
        ind = passing[-1]
        fdr_cutoffs[level] = {"pvalue_cutoff":float(candidates[ind]),
                              "phenolog_count":int(discoveries[ind]),
                              "expected_false_positives":float(false_positives[ind]),
                              "realized_fdr":float(fdr_curve[ind])}

    return fdr_cutoffs, {"pvalue":candidates,
                         "phenolog_count":discoveries,
                         "expected_false_positives":false_positives,
                         "fdr":fdr_curve}


def get_trial_fdr_data(trials_by_species: Dict, results_dir: str, fdr_levels: List[float]):

    species_trial_info = {}
    for k in trials_by_species:

        print("- Collapsing {} vs. {}...".format(k[0], k[1]))

        # Collapse real results and random trials to {(c,N,m,n):occurrence} tallies
        # NOTE - Collapsing the random trials is the step that should be parallelized!
        real_params, common_ortholog_count = collapse_cross_species_phenologs(results_dir, k)
        random_params = collapse_cross_species_trials(trials_by_species[k], common_ortholog_count)

        # Everything from here down is fine single core...
        # One hg calculation per unique parameter set across both distributions, so that a given
        # discrete pvalue is the exact same float on the real and on the random side of the comparison
        hg_pvals = bulk_compute_hyper_geom(set(real_params) | set(random_params))

        random_pval_counts = param_counts_to_pval_counts(random_params, hg_pvals)
        real_pval_counts = param_counts_to_pval_counts(real_params, hg_pvals)

        random_bulk_info = collapsed_data_to_cumulitive(random_pval_counts)
        phen_bulk_info = collapsed_data_to_cumulitive(real_pval_counts)

        # Generate false discovery rate pvalue cutoffs by weighing the real phenolog counts at each
        # candidate threshold against the false positives the random trials produce there.
        # Trial count is taken per species comparison, as it is the divisor turning the pooled
        # permuted hits into an expected count of false positives for a single trial
        n_trials = len(trials_by_species[k])
        fdr_cutoffs, fdr_curve = get_fdr_cuttoffs(random_pval_counts, real_pval_counts, n_trials, fdr_levels)

        species_trial_info.update({k:{"fdr_cutoffs":fdr_cutoffs,
                                      "fdr_curve":fdr_curve,
                                      "n_trials":n_trials,
                                      "common_ortholog_count":common_ortholog_count,
                                      "random_data":random_bulk_info,
                                      "real_data":phen_bulk_info}})

        for level in fdr_levels:
            cutoff = fdr_cutoffs[level]
            print("- fdr {} -- pvalue <= {:.4g} yields {} phenologs ({:.1f} expected false positives)".format(level,
                                                                                                              cutoff["pvalue_cutoff"],
                                                                                                              format(cutoff["phenolog_count"], ','),
                                                                                                              cutoff["expected_false_positives"]))

    return species_trial_info


def format_pval(p, sigfigs=2):
    if p <= 0:
        return r"$p = 0$"

    exp = int(np.floor(np.log10(abs(p))))
    base = p / (10 ** exp)
    decimals = max(sigfigs - 1, 0)

    if exp == 0:
        return rf"$p = {p:.{decimals}f}$"
    else:
        return rf"$p = {base:.{decimals}f} \times 10^{{{exp}}}$"


def plot_fdr_cumul_results(trial_info: Dict,
                           plot_order: List,
                           fdr_default=DEFAULT_FDR,
                           savefig1=False,
                           savefig2=False,
                           display_figs=False):

    # Lay comparisons out in plot_order (species tokens), then append anything present but unlisted so
    # that a newly added species can't be silently left out of the figures
    plot_keys = [k for name in plot_order for k in sorted(trial_info) if k[1] == name]
    plot_keys += [k for k in sorted(trial_info) if k not in plot_keys]

    # nrows x 2 columns formatted multiplot figure
    nrows, ncols = int(math.ceil(len(plot_keys)/2)), 2
    fig = plt.figure()
    fig.set_size_inches(ncols*5, nrows*3) # Width x height
    for i, k in enumerate(plot_keys):

        ax = plt.subplot2grid((nrows, ncols), (i // ncols, i % ncols))

        # Random, real, fdr data
        random_data, random_data_frac, random_total = trial_info[k]["random_data"]
        phen_data, phen_data_frac, phen_total = trial_info[k]["real_data"]
        pval_cutoff = trial_info[k]["fdr_cutoffs"][fdr_default]["pvalue_cutoff"]

        ax.plot(random_data.T[0], random_data_frac, label="Random")
        ax.plot(phen_data.T[0], phen_data_frac, label="Real")

        # FDR data. A cutoff of 0.0 means no threshold reaches this fdr level, so there is no line to draw
        if pval_cutoff > 0.:
            ax.axvline(pval_cutoff,
                       linestyle='--',
                       color='red',
                       alpha=.5,
                       label="FDR = {}{}{}".format(fdr_default, '\n', format_pval(pval_cutoff)))
        else:
            ax.plot([], [], linestyle='--', color='red', alpha=.5,
                    label="FDR = {}{}not reached".format(fdr_default, '\n'))

        ax.set_xlim(min(random_data.T[0]) / 2., 1.)

        ax.invert_xaxis()
        ###ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel("Hypergeometric probability of observing{} gene intersection by chance".format('\n'))
        ax.set_ylabel("Fraction of phenologs")

        ax.legend()
        ax.set_title("{} vs. {}".format(k[0], k[1]))

    # Save and display options
    plt.tight_layout()
    if savefig1 != False:
        plt.savefig(savefig1)
        plt.savefig(savefig1.replace(".pdf", ".png"), dpi=600)
    if display_figs != False:
        plt.show()
    else:
        plt.close()

    ############################################################################
    ### Plot the false discovery rate against phenolog yield for all species ###
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))

    for k in plot_keys:
        fdr_curve = trial_info[k]["fdr_curve"]
        ax.plot(fdr_curve["phenolog_count"], fdr_curve["fdr"], label=k[1])

    ax.axhline(fdr_default, linestyle='--', color='red', alpha=.5, label="FDR = {}".format(fdr_default))

    ax.set_xlabel("Number of phenologs at or below pvalue threshold")
    ax.set_ylabel("False discovery rate")
    ax.legend(loc='lower right', fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Save and display options
    if savefig2 != False:
        plt.savefig(savefig2, bbox_inches='tight')
        plt.savefig(savefig2.replace(".pdf", ".png"), bbox_inches='tight', dpi=600)
    if display_figs != False:
        plt.show()
    else:
        plt.close()


# Leveraged in phenologs_randomized_fdr_pvals
def initiate_phenologs_fdr_fpaths(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through.
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Ensures part1 & part2 of pipeline have been completed
    check_outdir = os.path.join(input_args.project_dir, "random_trials")
    if not os.path.isdir(check_outdir):
        print("- ERROR, Project random_trials directory doesn't seem to exist. Exiting...")
        sys.exit(1)

    # Which species this run has, from the manifest. "all" means every included species in turn
    manifest = SpeciesManifest(input_args.project_dir)
    if input_args.taxon_id == "all":
        base_species = list(manifest)
    else:
        base_species = [manifest.get(input_args.taxon_id)]

    # Now we need to filter our comparison species list by relevant / available prediction networks
    # "phenotype" or "disease" networks are available for use
    res_paths, rand_paths, sp_names = [], [], []

    for sp in base_species:
        rand_dir = os.path.join(input_args.project_dir, "random_trials", trials_dir_name(sp.token, input_args.prediction_network))
        res_dir = os.path.join(input_args.project_dir, "phenologs_results", results_dir_name(sp.token, input_args.prediction_network))

        if os.path.isdir(res_dir) and os.path.isdir(rand_dir):
            res_paths.append(res_dir)
            rand_paths.append(rand_dir)
            sp_names.append(sp.token)
        else:
            print("- Warning, random_trials and or phenologs_results appear to be missing for {}".format(sp.token))

    if len(res_paths) == 0:
        print('- ERROR, No relevant data found for taxon_id {} of type {}. Exiting...'.format(input_args.taxon_id,
                                                                                              input_args.prediction_network))
        sys.exit(1)

    return res_paths, rand_paths, sp_names, manifest




if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='-p argument OR -i & -o arguments are required. Data will  be pulled \
                                                      from an assumed directory / project structure. This script will gather \
                                                      all randomized trial results from speciesA-speciesB comparisons, and finds \
                                                      the pvalues corresponding to different false discovery rates. \
                                                      SpeciesA, SpeciesB, fdr, pval... for all data is written to _fdr_table.tsv')

        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-display_figs", help='Option to walkthrough analysis figure by figure in separate windows', required=False, action='store_true')
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    fpaths, rpaths, snames, manifest = initiate_phenologs_fdr_fpaths(args)
    for res_path, rand_path, sname in zip(fpaths, rpaths, snames):

        # Gather random trial data and compute fdr cuttoff values from hg pvalues
        trial_data = gather_trial_data(rand_path)
        trial_info = get_trial_fdr_data(trial_data, res_path, FDR_LEVELS)

        # Generate our outfile paths for data and plots
        fdr_table_outpath = os.path.join(res_path, "{}_fdr_table.tsv".format(sname))
        plot_data_outpath = os.path.join(res_path, "{}_fdr_plot_data.pkl".format(sname))
        plot1_data_outpath = os.path.join(res_path, "{}_xspecies_random_vs_real.pdf".format(sname))
        plot2_data_outpath = os.path.join(res_path, "{}_xspecies_phenolog_counts.pdf".format(sname))

        # Generate and write fdr table
        fdr_table = [{"species_comparison":xspecies_key(*k),
                      "fdr":level,
                      "pvalue_cutoff":cutoff["pvalue_cutoff"],
                      "phenolog_count":cutoff["phenolog_count"],
                      "expected_false_positives":cutoff["expected_false_positives"],
                      "realized_fdr":cutoff["realized_fdr"],
                      "n_random_trials":trial_info[k]["n_trials"]}
                     for k in trial_info for level, cutoff in sorted(trial_info[k]["fdr_cutoffs"].items())]

        pd.DataFrame(fdr_table).to_csv(fdr_table_outpath, sep='\t', index=False)

        # Write out data used for making plots
        pickle.dump(trial_info, open(plot_data_outpath, 'wb'))

        # Plot our results and save figures
        plot_fdr_cumul_results(trial_info,
                               plot_order=manifest.plot_order(exclude=[sname]),
                               fdr_default=DEFAULT_FDR,
                               savefig1=plot1_data_outpath,
                               savefig2=plot2_data_outpath,
                               display_figs=args.display_figs)

        print("- FDR table and figure(s) generated... Done!")
