# General imports
import os
import sys
import argparse
import pandas as pd
import copy
from collections import Counter
from IPython.display import display

# Custom imports
from phenologs_utils import (Species,
                             build_species_info,
                             hyper_geom)


def get_gene_phenotype_networks(species_info_dict, taxon_id):
    
    # Read our common orthologs into memory for our species of interest (should be human)
    base_species_name = species_info_dict[str(taxon_id)].species_name
    common_orths = {}
    global_orths = set()
    for fpath in species_info[str(taxon_id)].common_orthologs_paths:
        species_name = fpath.split('/')[-1].split("_vs_")[-1].replace(".tsv", "")
        orths = set(pd.read_csv(fpath)["ortholog_id"])
        common_orths.update({species_name:orths})
        global_orths = global_orths | orths
    
    print("- Using input arg taxon_id {} --> {} as base species...".format(taxon_id, base_species_name))
    print("- Total protein families with >= 1 cross species orthology {}...".format(format(len(global_orths), ',')))
    for k,v in common_orths.items():
        print("- {}, Common protein families found with {}".format(k, format(len(v), ',')))
    
    # Read gene-phenotype networks into memory for each species, but only include the genes
    # that can be mapped to our species of interest (i.e. only common protein families)
    
    # Need to figure out how many times a gene shows up (and sum across species)
    xspecies_g2p = {}
    xspecies_p2g = {}
    for species_name in common_orths.keys():
        g2p_df = pd.read_csv(species_info[species_name].gene_to_phenotype_path, sep='\t')
        for g_id, p_id, orth_id in zip(g2p_df["gene"],  
                                       g2p_df["phenotype"], 
                                       g2p_df["ortholog_id"]):
            
            gkey = (species_name, g_id)
            pkey = (species_name, p_id)
            
            # We only want genes that are orthologous between our base species and comparison species
            if orth_id in common_orths[species_name]:
                
                # Gene to phenotype network 
                if gkey not in xspecies_g2p:
                    xspecies_g2p.update({gkey:{}})
                xspecies_g2p[gkey].update({p_id:''})
                
                # Phenotype to gene network
                if pkey not in xspecies_p2g:
                    xspecies_p2g.update({pkey:{}})
                xspecies_p2g[pkey].update({g_id:''})
    
    return xspecies_g2p, xspecies_p2g


def read_ortholog_edges_to_gene_maps(orths_edges_file, base_taxon_id: str, taxon_id_to_name: dict):
    
    # Read data into memory as df
    orth_edges_df = pd.read_csv(orths_edges_file, sep='\t')
    
    base_taxon = base_taxon_id
    if "NCBITaxon:" not in base_taxon_id:
        base_taxon = "NCBITaxon:{}".format(base_taxon_id)
    base_taxon_name = taxon_id_to_name[base_taxon]
        

    base_taxon_map = {}
    species_gene_names = {}
    for tx_a, g_a, g_an, tx_b, g_b, g_bn in zip(orth_edges_df["taxon_a"],
                                                orth_edges_df["gene_a"],
                                                orth_edges_df["gene_a_name"],
                                                orth_edges_df["taxon_b"],
                                                orth_edges_df["gene_b"],
                                                orth_edges_df["gene_b_name"]):
        

        # Pull in data for select species
        if (tx_a not in taxon_id_to_name) or (tx_b not in taxon_id_to_name):
            continue
        
        # Format our taxon_id (NCBITaxon:...--> human readable name)
        tx_a = taxon_id_to_name[tx_a]
        tx_b = taxon_id_to_name[tx_b]
        
        # Update our gene_id --> gene_names
        if g_a not in species_gene_names:
            species_gene_names.update({g_a:g_an})
        elif g_an != species_gene_names[g_a]:
            raise ValueError("- Different name found for gene id {}-->{}, {}".format(g_a, 
                                                                                     species_gene_names[g_a],
                                                                                     g_an))
        if g_b not in species_gene_names:
            species_gene_names.update({g_b:g_bn})
        elif g_bn != species_gene_names[g_b]:
            raise ValueError("- Different name found for gene id {}-->{}, {}".format(g_b, 
                                                                                     species_gene_names[g_b],
                                                                                     g_bn))
        
        # Only deal with edges that include a base_taxon gene node
        if (base_taxon_name != tx_a) and (base_taxon_name != tx_b):
            continue

        if base_taxon_name == tx_a:
            map_to = g_a
            map_from = (tx_b, g_b)
        else:
            map_to = g_b
            map_from = (tx_a, g_a)

        if map_from not in base_taxon_map:
            base_taxon_map.update({map_from:{}})

        base_taxon_map[map_from].update({map_to:''})
    
    
    base_genes_mapped = len(set([vv for v in base_taxon_map.values() for vv in v]))
    print("- Cross species orthologs genes mapped back to taxon {} = {}".format(base_taxon,
                                                                                format(len(base_taxon_map), ',')))
    print("- Taxon {} genes with mapping from >= 1 species = {}".format(base_taxon, 
                                                                        format(base_genes_mapped, ',')))
    return base_taxon_map, species_gene_names


def compute_background_base_g2p(xspecies_g2p, base_taxon_gmap):
    
    base_taxon_gene_background = Counter()
    for k,v in xspecies_g2p.items():

        if k not in base_taxon_gmap:
            continue

        base_genes = base_taxon_gmap[k]
        for gg in base_genes:
            base_taxon_gene_background[gg] += len(v)
    
    tot = len(base_taxon_gene_background)
    print("- {} genes found >= 1 phenotype association across all species...".format(format(tot, ',')))
    return base_taxon_gene_background
    

def group_sig_phenologs_by_phen(sig_phens_path):

    if sig_phens_path.endswith(".gz"):
        sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t', compression="gzip")
    else:
        sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t')

    # Build "k-nearest neighbor" data structure for each "phenotype" we want to assign gene rankings to
    # First, map row indices to each phenotype
    ind = 0
    p2_phens = {}

    xspecies_filtered = 0
    #Species A Phenotype ID  Species B Phenotype ID  Ortholog Count A        Ortholog Count B        Overlap Count   Common Ortholog Count   hg_pval Species A Phenotype Name        Species B Phenotype Name             X Species Comparison
    for phen_id, phen_dist_val, xcomp_name in zip(list(sig_phenologs_df["Species A Phenotype ID"]), 
                                                  list(sig_phenologs_df["hg_pval"]),
                                                  list(sig_phenologs_df["X Species Comparison"])):

        if phen_id not in p2_phens:
            p2_phens.update({phen_id:[]})
        p2_phens[phen_id].append(ind)
        ind += 1

    # Now make mini dataframes for each "phenotype" identifier and sort by best "distance" 
    # for whichever metric is chosen
    zero_sig_phens = {}
    for k in p2_phens:
        p2_phens[k] = copy.copy(sig_phenologs_df.iloc[p2_phens[k]])
        p2_phens[k] = p2_phens[k].sort_values(by="hg_pval")

        # Ensure data is of proper type for calculations
        p2_phens[k]["Overlap Count"] = p2_phens[k]["Overlap Count"].astype(float)
        p2_phens[k]["Ortholog Count B"] = p2_phens[k]["Ortholog Count B"].astype(float)
    
    return p2_phens


def rank_disease_top_phenolog_genes(species_phenotype2gene, 
                                    phenotype2phenotype_results_df, 
                                    map_to_base, 
                                    gene_name_map,
                                    gene_disease_map,
                                    disease_name_map,
                                    outdirectory=False):
    
    # If we want to write results files or not we need to make our output directory
    if outdirectory != False:
        os.makedirs(outdirectory, exist_ok=True)

    hg_N = len(species_phenotype2gene) # Global phenotype count across all species
    count = 0
    for k,v in phenotype2phenotype_results_df.items():
        
        # "Draw size" parameter for hg test
        hg_n = len(v) # How many significant phenologs found for this phenotype

        # Phenotype and species names
        phen_ids = list(v["Species B Phenotype ID"])
        phen_names = list(v["Species B Phenotype Name"])
        xspecies = [vv.split(",")[1] for vv in list(v["X Species Comparison"])]

        # Now loop through each phenolog (phenotype) and pull out genes associated with it
        base_mapped_genes = {}
        for p_id, p_name, xspe in zip(phen_ids, phen_names, xspecies):
            phenolog_genes = species_phenotype2gene[(xspe, p_id)]

            #if xspe != "Dictyostelium-discoideum":
            #    continue

            # Loop through genes associated with said phenolog/phenotype
            for pgene in phenolog_genes:
                key = (xspe, pgene)

                if key not in map_to_base:
                    continue

                base_genes = map_to_base[key]

                # Tally how many species show up for this gene, and how many times it shows up
                for bgene in base_genes:
                    if bgene not in base_mapped_genes:
                        base_mapped_genes.update({bgene:{"xspecies":Counter(), "occurence":0, "phenotypes":Counter()}})

                    base_mapped_genes[bgene]["xspecies"][xspe] += 1
                    base_mapped_genes[bgene]["occurence"] += 1
                    base_mapped_genes[bgene]["phenotypes"][p_name] += 1
                    

        #print("- {} gene candidate ranking info | Significant phenolog count {}".format(k, len(v)))
        #print("- Total genes {}".format(len(base_mapped_genes)))
        
        out_df = {"Gene":[], 
                  "Gene_name":[], 
                  "Association_exists":[], 
                  "XSpecies_occurence":[], 
                  "p-value":[], 
                  "XSpecies_phenotypes":[]}
        
        for bgene, gdata in base_mapped_genes.items():
            hg_m = base_taxon_g2p_occurence[bgene]
            hg_c = gdata["occurence"]
            pval = hyper_geom(hg_c, hg_N, hg_m, hg_n)

            out_df["Gene"].append(bgene)
            out_df["Gene_name"].append(gene_name_map[bgene])
            
            # Tells us weather this gene-->disease association has been found before or not
            if bgene in gene_disease_map:
                out_df["Association_exists"].append("True")
            else:
                out_df["Association_exists"].append("False")
            
            out_df["XSpecies_occurence"].append(','.join(["{}:{}".format(kk,vv) for kk,vv in gdata["xspecies"].items()]))
            out_df["p-value"].append(pval)
            out_df["XSpecies_phenotypes"].append(';'.join(list(gdata["phenotypes"].keys())))
            #print("  - Gene={} | Species={} | p-value={}".format(bgene, gdata["xspecies"], pval))

        out_df = pd.DataFrame(out_df)
        out_df_sorted = out_df.sort_values(by="p-value")
        
        # Write data
        if outdirectory != False:
            dname_formatted = disease_name_map[k].replace(" ", "-").replace(",", "-").replace("/", "-")
            outname = "{}_{}_gene_candidates.tsv".format(k, dname_formatted)
            outfile_path = os.path.join(outdirectory, outname)
            out_df_sorted.to_csv(outfile_path, sep='\t', index=False)

        
        # Progress statement
        count += 1
        if count % 1_000 == 0:
            print('- Processed {}/{}'.format(count, len(phenotype2phenotype_results_df)))
    
    print("- Creation of disease candidate gene rankings complete!")
    return


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes randomized comparison trials between species a vs species b')
        parser.add_argument("-p", "--project_dir", help="Top most project directory", required=True, type=str)
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=False, type=float, default=.95)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
###############
### PROGRAM ###


# Currently this is hard coded to human, but could be applied to non-human diseases like horse when they are
# eventualy inorporated into the graph
taxon_id = "9606"
specieces_data_dir = os.path.join(args.project_dir, "species_data")
disease_results_dir = os.path.join(args.project_dir, "phenologs_results", "Homo_sapiens_disease_results")
sig_phens_file = os.path.join(disease_results_dir, "Homo-sapiens_pooled_phenologs_fdr{}.tsv".format(args.fdr))
orth_edges_file = os.path.join(specieces_data_dir, "ortholog_edges.tsv")
output_directory = os.path.join(disease_results_dir, "Homo_sapiens_disease_gene_candidates_fdr{}".format(args.fdr))

# Make sure our relevant files exist beforehand
if not os.path.isfile(sig_phens_file):
    raise FileNotFoundError("- pooled phenologs results file not found {}...".format(sig_phens_file))

if not os.path.isfile(sig_phens_file):
    raise FileNotFoundError("- orthologs edges file not found {}...".format(orth_edges_file))


# Initiate our species level information (so we can easily pull relevant filepaths)
tx_ids, species_info = build_species_info(args.project_dir)
species_g2p, species_p2g = get_gene_phenotype_networks(species_info, taxon_id)
map_to_base, species_gnames = read_ortholog_edges_to_gene_maps(orth_edges_file, base_taxon_id=str(taxon_id), taxon_id_to_name=tx_ids)
base_taxon_g2p_occurence = compute_background_base_g2p(xspecies_g2p=species_g2p, base_taxon_gmap=map_to_base)
p2p_dict_df = group_sig_phenologs_by_phen(sig_phens_file)

# Grab human disease name info
g2d = pd.read_csv(species_info[str(taxon_id)].gene_to_disease_path, sep='\t')
disease_name_map = {d_id:dname for d_id,dname in zip(g2d["disease"], g2d["disease_name"])}
gene2disease_map = {g_id:dname for g_id,dname in zip(g2d["gene"], g2d["disease_name"])}

# Now pull everything together and compute our disease candidate gene rankings
print("- Computing disease candidate gene rankings for taxon id {}...".format(taxon_id))
rank_disease_top_phenolog_genes(species_phenotype2gene=species_p2g, 
                                phenotype2phenotype_results_df=p2p_dict_df,
                                map_to_base=map_to_base,
                                gene_name_map=species_gnames,
                                gene_disease_map=gene2disease_map,
                                disease_name_map=disease_name_map,
                                outdirectory=output_directory)