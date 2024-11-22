import os
import pickle
import random
import copy
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from collections import Counter
from IPython.display import display
from pydantic import BaseModel
from typing import Union, Literal, Dict, List, Any, Optional


def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


def initiate_random_species_comparison_configs(input_args, fdr_path=None):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    """

    # Ensures part1 & part2 of pipeline have been completed
    check_file = os.path.join(input_args.project_dir, "species_data", "species_information_table.tsv")
    check_outdir = os.path.join(input_args.project_dir, "random_trials")
    
    if not os.path.isfile(check_file):
        print("- ERROR, Project species information table doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isdir(check_outdir):
        print("- ERROR, Project random_trials directory doesn't seem to exist. Exiting...")
        sys.exit()
    
    # Figure out which species ids / names we have and which ones are relevant
    species_df = pd.read_csv(check_file, sep='\t')
    species_df = species_df[species_df['Total Phenotype Edges'] > 0] # Only want species whith non zero phenotype information
    

    # Pull out taxon_id information (twice as two separate variables)
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    org_taxon_ids = copy.copy(ids_to_name)
    
    
    # Format, and add additional keys to make more friendly to input arguments
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    ids_to_name.update({sp_id.split(":")[1]:v for sp_id,v in ids_to_name.items()}) # Add keys without NCBITaxon: prefix
    display(species_df)

    # Simple table mapping species name to the total number of genes (with at least one phenotype term associated with them)
    species_g_count = {t:int(v) for t,v in zip(list(species_df["Taxon Label"]), list(species_df["Genes >= 1 Phenotype"]))}

    # Figure out which species ids we are tasked with comparing to one another
    valid_species_ids = set(ids_to_name.keys())
    if not input_args.all:
        if not input_args.taxon_ids:
            print('- ERROR, -all or -taxon_ids argument must be specified. Exiting...')
            sys.exit()
        
        t_ids = input_args.taxon_ids.split(",")
        tot_species = len(t_ids)
        print(valid_species_ids)

        # Ensure input is compatible with data in graph
        if tot_species < 2:
            print('- ERROR, Total number of input taxon ids must be greater than 1. For example 9606,8355... Exiting...')
            sys.exit()
        
        if len(set(t_ids) & valid_species_ids) != len(t_ids):
            print('- ERROR, One or more incompatible taxon_ids input. Exiting...')
            sys.exit()
        
    else:
        t_ids = list(org_taxon_ids.keys())
    
    tot_species = len(t_ids)
    tot_comps = (tot_species*(tot_species-1))/2
    print("- Species found for comparison {}, Total comparisons to make {}".format(tot_species, tot_comps))

    # Generate our input configurations list
    # This ensures that all filepaths exist for all comparisons before we start making any computations

    # The following double for loop performs all pairwise comparisons with out repeats (i.e. a->b but not b->a)
    # TO DO: Do we want to double this? (i.e. a->b, and b->a)
    configs = []
    for i, species_a in enumerate(t_ids[:-1]):
        a_name = ids_to_name[species_a]
        for species_b in t_ids[i+1:]:
            b_name = ids_to_name[species_b]

            # Base level config
            cpath = os.path.join(input_args.project_dir, "species_data", "common_orthologs_{}_vs_{}.tsv".format(a_name, b_name))
            apath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(a_name))
            bpath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(b_name))

            # Ensure all files exist
            if (not os.path.isfile(cpath)) or (not os.path.isfile(apath)) or (not os.path.isfile(bpath)):
                print("- ERROR, Species {} vs {} is missing common orthologs and or phenotype files...".format(a_name,
                                                                                                               b_name))
                sys.exit()
            
            # Generate two configs for comparison (Species A-->B, and B-->A)
            config_a = {"species_a":a_name,
                        "species_b":b_name,
                        "species_a_phenotype_path":apath,
                        "species_b_phenotype_path":bpath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir,
                        "fdr_path":fdr_path}
            
            # Swap a & b values to get comparison in other direction
            config_b = {"species_a":b_name,
                        "species_b":a_name,
                        "species_a_phenotype_path":bpath,
                        "species_b_phenotype_path":apath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir,
                        "fdr_path":fdr_path}
            
            # Just do one config for now
            configs.append(config_a)
            #configs.append(RandomSpeciesComparison.parse_obj(config_a))
    
    print("- {} Pairwise comparisons configurations created...".format(len(configs)))
    return t_ids, configs


def initiate_phenologs_species_comparison_configs(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    """

    # Ensures parts 1, 2, and 3 have been completed (species info table, random trial results and fdr table)
    check_file1 = os.path.join(input_args.project_dir, "species_data", "species_information_table.tsv")
    check_file2 = os.path.join(input_args.project_dir, "random_trials_fdr", "fdr_table.tsv")
    check_inputdir = os.path.join(input_args.project_dir, "random_trials")
    check_outdir = os.path.join(input_args.project_dir, "phenologs_results")
    
    if not os.path.isfile(check_file1):
        print("- ERROR, Project species information table doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isfile(check_file2):
        print("- ERROR, fdr table file doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isdir(check_inputdir):
        print("- ERROR, Project random_trials directory doesn't seem to exist. Exiting...")
        sys.exit()
    
    elif len([fname for fname in os.listdir(check_inputdir) if "_vs_" in fname]) == 0:
        print("- ERROR, Project random_trials directory has no relevant file types ('_vs_' in the name). Exiting...")
    
    if not os.path.isdir(check_outdir):
        print("- Making output directory at {}".format(check_outdir))
        os.makedirs(check_outdir)
    
    # Figure out which species ids / names we have and which ones are relevant
    species_df = pd.read_csv(check_file1, sep='\t')
    species_df = species_df[species_df['Total Phenotype Edges'] > 0] # Only want species whith non zero phenotype information
    
    # Pull out taxon_id information (twice as two separate variables)
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    org_taxon_ids = copy.copy(ids_to_name)
    
    # Format, and add additional keys to make more friendly to input arguments
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    ids_to_name.update({sp_id.split(":")[1]:v for sp_id,v in ids_to_name.items()}) # Add keys without NCBITaxon: prefix
    display(species_df)

    # Simple table mapping species name to the total number of genes (with at least one phenotype term associated with them)
    species_g_count = {t:int(v) for t,v in zip(list(species_df["Taxon Label"]), list(species_df["Genes >= 1 Phenotype"]))}

    # Figure out which species ids we are tasked with comparing to one another
    valid_species_ids = set(ids_to_name.keys())
    if not input_args.all:
        if not input_args.taxon_ids:
            print('- ERROR, -all or -taxon_ids argument must be specified. Exiting...')
            sys.exit()
        
        t_ids = input_args.taxon_ids.split(",")
        tot_species = len(t_ids)
        print(valid_species_ids)

        # Ensure input is compatible with data in graph
        if tot_species < 2:
            print('- ERROR, Total number of input taxon ids must be greater than 1. For example 9606,8355... Exiting...')
            sys.exit()
        
        if len(set(t_ids) & valid_species_ids) != len(t_ids):
            print('- ERROR, One or more incompatible taxon_ids input. Exiting...')
            sys.exit()
        
    else:
        t_ids = list(org_taxon_ids.keys())
    
    tot_species = len(t_ids)
    tot_comps = (tot_species*(tot_species-1))/2
    print("- Species found for comparison {}, Total comparisons to make {}".format(tot_species, tot_comps))

    # Generate our input configurations list
    # This ensures that all filepaths exist for all comparisons before we start making any computations

    # The following double for loop performs all pairwise comparisons with out repeats (i.e. a->b but not b->a)
    # TO DO: Do we want to double this? (i.e. a->b, and b->a)
    configs = []
    for i, species_a in enumerate(t_ids[:-1]):
        a_name = ids_to_name[species_a]
        for species_b in t_ids[i+1:]:
            b_name = ids_to_name[species_b]

            # Base level config
            cpath = os.path.join(input_args.project_dir, "species_data", "common_orthologs_{}_vs_{}.tsv".format(a_name, b_name))
            apath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(a_name))
            bpath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(b_name))

            # Ensure all files exist
            if (not os.path.isfile(cpath)) or (not os.path.isfile(apath)) or (not os.path.isfile(bpath)):
                print("- ERROR, Species {} vs {} is missing common orthologs and or phenotype files...".format(a_name,
                                                                                                               b_name))
                sys.exit()

            # Generate two configs for comparison (Species A-->B, and B-->A)
            config_a = {"species_a":a_name,
                        "species_b":b_name,
                        "species_a_phenotype_path":apath,
                        "species_b_phenotype_path":bpath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir,
                        "fdr_path":check_file2} # This is the fdr table file computed from the randomized trials
            
            # Swap a & b values to get comparison in other direction
            config_b = {"species_a":b_name,
                        "species_b":a_name,
                        "species_a_phenotype_path":bpath,
                        "species_b_phenotype_path":apath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir,
                        "fdr_path":check_file2}
            
            # Just do one config for now
            configs.append(config_a)
            #configs.append(RandomSpeciesComparison.parse_obj(config_a))
    
    print("- {} Pairwise comparisons configurations created...".format(len(configs)))
    return t_ids, configs


def bulk_compute_hyper_geom(params):
    return {pr:float(hypergeom.pmf(*pr)) for pr in params}


class SpeciesComparison(BaseModel):
    species_a: str
    species_b: str
    species_a_phenotype_path: str
    species_b_phenotype_path: str
    common_orthologs_path: str

    base_a_p2g: Optional[dict] = None
    base_b_p2g: Optional[dict] = None
    base_pool: Optional[list] = None
    common_ortholog_count: Optional[int] = None
    output_directory: Optional[str] = None
    
    fdr_path: Optional[str] = None
    
    
    def load_species_ortholog_information(self):
        """
        Initiates necessary attributes for downstream computations.
        - Number of common orthologs
        - The common ortholog pool to sample from for randomized data
        - Base phenotype-->ortholog dictionaries for each species
        """
        
        # Load species dict and common orthologs and set remaining filepaths
        base_a_p2g = pickle.load(open(self.species_a_phenotype_path, 'rb'))
        base_b_p2g = pickle.load(open(self.species_b_phenotype_path, 'rb'))
        corth_df = pd.read_csv(self.common_orthologs_path, sep='\t', header=0, low_memory=False)
        common_orthologs = list(corth_df["ortholog_id"])
        common_ortholog_count = len(common_orthologs)
        
        # Set our reusable attributes
        self.base_a_p2g = base_a_p2g
        self.base_b_p2g = base_b_p2g
        self.base_pool = common_orthologs
        self.common_ortholog_count = common_ortholog_count

        print("- Species {} vs. {} comparison data loaded... {} common orthologs found".format(self.species_a,
                                                                                               self.species_b,
                                                                                               format(self.common_ortholog_count, ',')))
     
    
    def subset_species_to_common_orthologs(self):

        common_orths = set(self.base_pool)
        a_sub, b_sub = {}, {}
        for k,v in self.base_a_p2g.items():
            orth_list = list(set([vv for vv in v if vv in common_orths]))
            if len(orth_list) > 0:
                a_sub.update({k:orth_list})

        for k,v in self.base_b_p2g.items():
            orth_list = list(set([vv for vv in v if vv in common_orths]))
            if len(orth_list) > 0:
                b_sub.update({k:orth_list})

        return common_orths, a_sub, b_sub


class RandomSpeciesComparison(SpeciesComparison):
    super(SpeciesComparison)
    
    def run_randomized_comparison_pipeline(self, seed: Optional[int] = None):
        """
        Generates a randomized phenotype-->ortholog dataset for each species.
        Pairwise comparisons between inter species phenotypes are made
        by computing the number of orthologs that overlap between them. Note that
        ths does NOT compute pvalues, but rather gathers the set(s) of parameters used
        to calculate pvalues and full dataset of non zero overlap count parameters.
        """
        
        # Optional random seed
        if seed:
            np.random.seed(seed)
        
        # Load initial data
        self.load_species_ortholog_information()
        
        # NEW / Current way - Subset out only the common orthologs for each species
        common_orths, a_sub, b_sub = self.subset_species_to_common_orthologs()

        #print("- {}, {}".format(len(a_sub), len(b_sub)))
        #print("- {}".format(self.common_ortholog_count))
        #print("- {}, {}".format(max([len(v) for v in a_sub.values()]), max([len(v) for v in b_sub.values()])))

        # Generate randomized dataset
        randomized_a = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in a_sub.items()}
        randomized_b = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in b_sub.items()}

        # OLD WAY (Will leave until confirmation of proper hg params is found)
        # Generate randomized dataset
        ##randomized_a = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in self.base_a_p2g.items()}
        ##randomized_b = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in self.base_b_p2g.items()}
        ##print("- Randomized data generation complete")
        
        # Convert to compact list data structure (instead of big numpy array)
        b_matrix = [v for v in randomized_b.values()]
        b_matrix_lengths = {i:len(v) for i, v in enumerate(b_matrix)}
        hg_a_count_params = {k:len(v) for k,v in randomized_a.items()}
        hg_b_count_params = {k:len(v) for k,v in randomized_b.items()}

        # Map each unique value (ortholog id) in our data to the rows it belongs to (pre processing step)
        orth_to_coords = {}
        for i, row in enumerate(b_matrix):
            for v in row:
                if v not in orth_to_coords:
                    orth_to_coords.update({v:[]})
                orth_to_coords[v].append(i)
        
        #print("- Total rows {}".format(format(len(orth_to_coords), ',')))
        #print("- Total row associations {}".format(format(sum([len(v) for v in orth_to_coords.values()]), ',')))
        
        # Precompute b matrix length for initiation of zeros array within for loop
        b_phenotype_count = len(b_matrix)
        processed = 0
        t_comps = 0
        
        # For keeping track of our hg params
        hg_a_counts = []
        hg_b_counts = []
        hg_overlap_counts = []
        hg_world_counts = []
        
        # Loop through all phenotype-ortholog associatoins for species a, and compute overlap with species b 
        for a_phen, a_orthologs in randomized_a.items():
            
            # Initialize datastructure to count the number of overlaps found for each of species b phen-orthologs
            counts = np.zeros(b_phenotype_count)
            a_ortho_count = hg_a_count_params[a_phen]
            
            # Compute commonality between each ortholog and species b orthologs
            for orth in a_orthologs:
                if orth in orth_to_coords:
                    counts[orth_to_coords[orth]] += 1
                t_comps += 1
            
            # This gives back an oddly shaped data array that we can flatten to get all non-zero row indices
            inds = np.argwhere(counts > 0).flatten()
            ind_count = len(inds)
            

            # TO DO: Confirm proper hg parameters are pulled here. The old way (previous code base) would break
            # because in the most literal sense, if you physically assigned our data to marbles and bags, you would 
            # have actually been drawing from different bags, which the test is not designed for. 
            
            # World size parameter is set to the number of common orthologs between species a and b.
            # Therefore, our phenotype<-->gene networks must ONLY consist of genes that are orthologous between a and b.
            # In reality, species a phenotypeXYZ can be associated with thousands of genes (but only a handful of them
            # might be common orthologs). If we were to take all genes from a and b (not just common orthologs) then
            # we could wind up in a situation where the number of genes associated with XYZ phenotype (in species a or b)
            # is MORE than the number of orthologs those species share in common. By doing it that way, you not only
            # run into an edge case that will break that hyper geometric test, but it doesn't make sense to include
            # all genes, given our initial world size parameter is derived from only the common orothologs, and not
            # the totality of species a or species b genes. Should only be pulling from the common pool of orthologs
            # between them and the hyper geometric test should take care of the rest

            # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
            # c = number of common orthologs between phenotypes (ortholog matches)
            # N = total number of orthologs shared between species
            # m = number of orthologs in species B phenotype
            # n = number of orthologs in species A phenotype
            
            hg_overlap_counts += list(counts[inds])
            hg_world_counts += [self.common_ortholog_count] * ind_count
            hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            hg_a_counts += [a_ortho_count] * ind_count

            #processed += 1
            #if processed % 1000 == 0:
                #print("- Processed {}/10,000".format(processed))

        
        # Combine our data into a list of hg params
        tt = np.asarray([hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]).T.astype(int).tolist()
        #print("- Total non zero hg params computed {}".format(format(len(tt), ',')))
        
        # Create unique set of tuples
        hg_formatted = list(map(tuple, tt))
        hg_params = set(hg_formatted)
        hg_param_counts = Counter(hg_formatted)
        #print("- Unique hg params computed {}".format(format(len(hg_params), ',')))

        return hg_params, hg_param_counts
        

    def run_randomized_comparison_trials(self, n_trials: List[int]):
        
        # Keep track of unique sets of hyper geometric tests to compute so we don't recalculate
        hg_params = set()
        data = {}
        for i, n_trial in enumerate(n_trials):

            # Generate randomized data trials
            trial_hg_params, trial_data = self.run_randomized_comparison_pipeline()
            
            # Update our hire level output datastructures
            hg_params = hg_params | trial_hg_params
            data.update({n_trial:trial_data})
            print("- Total hg_params computed {} for {}/{} trials".format(format(len(hg_params), ','), 
                                                                          i+1, 
                                                                          format(len(n_trials), ',')))
        
        print("- {} vs. {} pairwise phenotype-->ortholog networks overlaps computed for {} trials".format(self.species_a,
                                                                                                          self.species_b,
                                                                                                          len(n_trials)))
        
        # Now compute our hyper geometric tests in bulk
        hg_pvals = bulk_compute_hyper_geom(hg_params)
        print("- Hypergeometric tests computed {}...".format(format(len(hg_pvals), ',')))
        
        # Map to pvalues, and write data files
        for trial_num, trial_data in data.items():
            
            # Define our outfile path and map params to pvalues 
            outfile_path = os.path.join(self.output_directory, 
                                        "{}_vs_{}_{}.tsv.gz".format(self.species_a, self.species_b, trial_num))
                                        
            # Map data back to computed pvalues, and format our trial data for easy writing to file
            pvals = [hg_pvals[tuple(k)] for k in trial_data]
            #param_data = np.asarray(trial_data).T
            
            # Reduced file size version
            param_data = np.asarray(list(trial_data.keys())).T
            occurrence = list(trial_data.values())
            
            # Write data using pandas
            pd.DataFrame({"a_ortholog_count":param_data[0],
                          "b_ortholog_count":param_data[1],
                          "overlap_count":param_data[2],
                          "common_ortholog_count":param_data[3],
                          "hg_pval":pvals, ## Full file size version).to_csv(outfile_path, sep='\t', index=False)

                         # Compact / reduced file size version (10 fold line count reduction for human mouse comparison)
                         "occurrence":occurrence}).to_csv(outfile_path, sep='\t', index=False, compression='gzip')

            print("- Data written to file for trial {}...".format(trial_num))
        
        ##return hg_pvals, data
        return
    

    # EXPERIMENTAL (Need a separate function to write data for run_comparisons_parallel_v2)
    def write_trial_data(self, trial_list, hg_pvals):
        
        for data in trial_list:

            # Map to pvalues, and write data files
            for trial_num, trial_data in data.items():
                
                # Define our outfile path and map params to pvalues 
                outfile_path = os.path.join(self.output_directory, 
                                            "{}_vs_{}_{}.tsv.gz".format(self.species_a, self.species_b, trial_num))
                                            
                # Map data back to computed pvalues, and format our trial data for easy writing to file
                pvals = [hg_pvals[tuple(k)] for k in trial_data]
                #param_data = np.asarray(trial_data).T
                
                # Reduced file size version
                param_data = np.asarray(list(trial_data.keys())).T
                occurrence = list(trial_data.values())
                
                # Write data using pandas
                pd.DataFrame({"a_ortholog_count":param_data[0],
                            "b_ortholog_count":param_data[1],
                            "overlap_count":param_data[2],
                            "common_ortholog_count":param_data[3],
                            "hg_pval":pvals, ## Full file size version).to_csv(outfile_path, sep='\t', index=False)

                            # Compact / reduced file size version (10 fold line count reduction for human mouse comparison)
                            "occurrence":occurrence}).to_csv(outfile_path, sep='\t', index=False, compression='gzip')

                print("- Data written to file for trial {}...".format(trial_num))
        
        ##return hg_pvals, data
        return


    # EXPERIMENTAL This version is more piece mail and will not compute hg pvals or write data
    def run_randomized_comparison_trials_v2(self, n_trials: List[int]):
        
        # Keep track of unique sets of hyper geometric tests to compute so we don't recalculate
        hg_params = set()
        data = {}
        for i, n_trial in enumerate(n_trials):
            # Generate randomized data
            trial_hg_params, trial_data = self.run_randomized_comparison_pipeline()
            
            # Update our hire level output datastructures
            hg_params = hg_params | trial_hg_params
            data.update({n_trial:trial_data})
            print("- Total hg_params computed {} for {}/{} trials".format(format(len(hg_params), ','), 
                                                                          i+1, 
                                                                          format(len(n_trials), ',')))
        
        print("- {} vs. {} pairwise phenotype-->ortholog networks overlaps computed for {} trials".format(self.species_a,
                                                                                                          self.species_b,
                                                                                                          len(n_trials)))
        
        return data, hg_params
    

class PhenologsSpeciesComparison(SpeciesComparison):
    super(SpeciesComparison)
    
    
    def get_fdr_info_from_fdr_file(self):
        """
        Assumes fdr table file path is being passed in where the header looks like
        fdr table file should obey the following format -- Sinlge header line, with a single row of information

        For example:
        fdr:0.05    fdr:0.01    fdr:0.001    fdr:0.0001
        .005       .0005        .00005       .000001
        """
        
        fdrs = {}
        fdr_df = pd.read_csv(self.fdr_path, sep='\t')
        for fval in fdr_df.columns:
            if len(fdr_df[fval]) > 1:
                print("- ERROR, Likely wrong fdr table file type/path provided...")
                sys.exit()
            
            pval_cutoff = float(fdr_df[fval][0])
            fdrs.update({pval_cutoff:fval})
            
        return fdrs
    
    
    def get_phenolog_params(self):
        """
        Pairwise comparisons between inter species phenotypes are made
        by computing the number of orthologs that overlap between them. Note that
        ths does NOT compute pvalues, but rather gathers the set(s) of parameters used
        to calculate pvalues and full dataset of non zero overlap count parameters.
        """

        # Load initial data
        self.load_species_ortholog_information()
        common_orths, a_sub, b_sub = self.subset_species_to_common_orthologs()
        a_phens = list(a_sub.keys())
        b_phens = list(b_sub.keys())


        # Convert to compact list data structure 
        # (instead of big numpy array or sparse array which introduce lots of uncessary overhead here)
        b_matrix = [v for v in b_sub.values()]
        b_matrix_lengths = {i:len(v) for i, v in enumerate(b_matrix)}

        hg_a_count_params = {k:len(v) for k,v in a_sub.items()}
        hg_b_count_params = {k:len(v) for k,v in b_sub.items()}

        # Map each unique value (ortholog id) in our data to the rows it belongs to (pre processing step)
        orth_to_coords = {}
        for i, row in enumerate(b_matrix):
            for v in row:
                if v not in orth_to_coords:
                    orth_to_coords.update({v:[]})
                orth_to_coords[v].append(i)

        #print("- Total rows {}".format(format(len(orth_to_coords), ',')))
        #print("- Total row associations {}".format(format(sum([len(v) for v in orth_to_coords.values()]), ',')))

        # Precompute b matrix length for initiation of zeros array within for loop
        b_phenotype_count = len(b_matrix)
        processed = 0
        t_comps = 0

        # For keeping track of our hg params
        hg_a_counts = []
        hg_b_counts = []
        hg_overlap_counts = []
        hg_world_counts = []
        

        # Loop through all phenotype-ortholog associatoins for species a, and compute overlap with species b 
        species_a_phens, species_b_phens = [], []
        for a_phen, a_orthologs in a_sub.items():

            # Initialize datastructure to count the number of overlaps found for each of species b phen-orthologs
            counts = np.zeros(b_phenotype_count)
            ind_count = len(counts)
            a_ortho_count = hg_a_count_params[a_phen]

            # Compute commonality between each ortholog and species b orthologs
            for orth in a_orthologs:
                if orth in orth_to_coords:
                    counts[orth_to_coords[orth]] += 1
                t_comps += 1
            
            # This gives back an oddly shaped data array that we can flatten to get all non-zero row indices
            inds = np.argwhere(counts > 0).flatten()
            ind_count = len(inds)

            # If we want non zero overlap terms
            hg_overlap_counts += list(counts[inds].astype(int))
            hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            hg_a_counts += [a_ortho_count] * ind_count
            hg_world_counts += [self.common_ortholog_count] * ind_count
            species_a_phens += [a_phen] * ind_count
            species_b_phens += [b_phens[ind] for ind in inds]

            # This is if we want ALL data 
            # TO DO: Confirm proper hg params are pulled here
            #hg_overlap_counts += list(counts.astype(int))
            #hg_b_counts += [b_matrix_lengths[ind] for ind in range(0, b_phenotype_count)]
            #hg_a_counts += [a_ortho_count] * b_phenotype_count
            #hg_world_counts += [self.common_ortholog_count] * b_phenotype_count
            #species_a_phens += [a_phen] * b_phenotype_count
            #species_b_phens += b_phens

            #processed += 1
            #if processed % 1000 == 0:
                #print("- Processed {}/10,000".format(processed)

        # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
        # c = number of common orthologs between phenotypes (ortholog matches)
        # N = total number of orthologs shared between species
        # m = number of orthologs in species B phenotype
        # n = number of orthologs in species A phenotype
        
        # Combine our data into a list of hg params
        tt = np.asarray([hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]).T.astype(int).tolist()
        
        # Create unique set of tuples
        hg_formatted = list(map(tuple, tt))
        hg_params = set(hg_formatted)

        # Collate data into dictionary (easy format to work with downstream)
        return_data = {"Species A Phenotype ID":species_a_phens,
                       "Species B Phenotype ID":species_b_phens,
                       "Ortholog Count A":hg_a_counts,
                       "Ortholog Count B":hg_b_counts,
                       "Overlap Count":hg_overlap_counts,
                       "Common Ortholog Count":hg_world_counts}

        return hg_params, return_data
    
    
    def compute_cross_species_phenologs(self, out_directory):
        
        # Read in fdr table to get pvalue cutoffs
        fdr_info = self.get_fdr_info_from_fdr_file()
        
        # Compute parameters for hyper geometric tests
        hg_params, comparison_data = self.get_phenolog_params()
        print("- {} vs. {} -- hg parameters computed {}...".format(self.species_a,
                                                                   self.species_b,
                                                                   format(len(hg_params), ',')))
        
        # Compute hyper geometric tests
        hg_pvals = bulk_compute_hyper_geom(hg_params)
        
        # Map each row's set of hg params back to the repsective pvalue
        pval_col, sig_col = [], []
        for c_param, N_param, m_param, n_param in zip(comparison_data["Overlap Count"], 
                                                      comparison_data["Common Ortholog Count"],
                                                      comparison_data["Ortholog Count B"],
                                                      comparison_data["Ortholog Count A"]):
            
            # Note that ordering here matters
            key = tuple((c_param, N_param, m_param, n_param))
            pval = hg_pvals.get(key, 1.0)
            pval_col.append(pval)
            
        # Add pvalues and convert to data frame
        comparison_data.update({"P-Value":pval_col})
        comparison_data = pd.DataFrame(comparison_data)
        
        # Write select pvalue / fdr cuttoff phenologs
        for pval_cutoff,f_ext in fdr_info.items():

            # Format output filenames
            fname = "{}_vs_{}_phenologs_{}.tsv".format(self.species_a, self.species_b, f_ext)
            outpath_name = os.path.join(out_directory, fname)

            # Subset df and write data
            comparison_data[comparison_data["P-Value"] <= pval_cutoff].to_csv(outpath_name, sep='\t', index=False)
            print("- {} written...".format(outpath_name))
        
        # TO DO: Do we need this? It's a bit excessive in terms of filesize...
        # Write all data
        ##outpath_name = os.path.join(out_directory, "{}_vs_{}_fulldata.tsv".format(self.species_a, self.species_b))                                                                   
        ##comparison_data.to_csv(outpath_name)