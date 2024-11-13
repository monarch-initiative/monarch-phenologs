import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

from collections import Counter
from pydantic import BaseModel
from typing import Optional, Any, List, Dict, Union

from IPython.display import display


class Taxa(BaseModel):
    label: Optional[str]
    ID: Optional[str]
    name: Optional[str]
    gene_prefix: Optional[str]
    phenotype_prefix: Optional[str]
    gene_to_phenotype: Optional[Dict] = {}
    phenotype_ortholog: Optional[Dict] = {}
        

def pull_taxon_gene_information(nodes_filepath, write_taxon_info=False):

    # Our taxon information comes from gene nodes
    valid_nodes = {"biolink:Gene":'', "biolink:PhenotypicFeature":''}
    taxa_objs = {}
    genes_to_taxa = {}
    genes_to_name = {}
    phenotype_to_name = {}
    
    # Read data into memory via pandas and fill missing values with an empty string (instead of NaN)
    df = pd.read_csv(nodes_filepath, sep='\t', low_memory=False).fillna('')
    for n_id, n_cat, n_name, n_tax, n_tax_lab in zip(list(df["id"]), 
                                                     list(df["category"]),
                                                     list(df["name"]),
                                                     list(df["in_taxon"]), 
                                                     list(df["in_taxon_label"])):
        
        
        # Irrelivant node_types
        if n_cat not in valid_nodes:
            continue
        
        # We cannot link phenotypic feature nodes to taxa (until we read in edges)
        if n_cat == "biolink:PhenotypicFeature":
            phenotype_to_name.update({n_id:n_name})
            continue
            
        # First encounter of xyz species. Need to update data structures
        if (n_tax not in taxa_objs) and (n_tax != ''):
            
            # Initiate Taxa obj config and add new taxa to our data structure
            g_prefix = n_id.split(":")[0]
            tconfig = {"ID":n_tax, 
                       "label":n_tax_lab,
                       "gene_prefix":g_prefix}
            
            taxa_objs.update({n_tax:Taxa.parse_obj(tconfig)})
        
        # Link our gene_id to its respective taxon
        genes_to_taxa.update({n_id:taxa_objs[n_tax]})
        genes_to_name.update({n_id:n_name})
    
    # TO DO: Sanity check to ensure multiple different gene_prefixes are not associated with the same taxon
    # (Old code to do this no longer is compatible with this so I got rid of it)
    
    # Nicely format some general taxon information about what we just read into memory
    taxon_table = {"Taxon ID":list(taxa_objs.keys()),
                   "Taxon Label":[v.label for k, v in taxa_objs.items()],
                   "Taxon Gene Prefix":[v.gene_prefix for k, v in taxa_objs.items()],
                   "Taxon Gene Count":[len([1 for vv in genes_to_taxa.values() if v.ID == vv.ID]) 
                                       for k, v in taxa_objs.items()]}
    
    # Write data if desired
    if write_taxon_info != False:
        pd.DataFrame(taxon_table).to_csv(write_taxon_info, sep='\t')
    
    return taxa_objs, genes_to_taxa, genes_to_name, phenotype_to_name, taxon_table


def fill_in_gene_to_phenotype_networks(edges_filepath, genes_to_taxa, taxa_objs, taxa_table):
    
    # Figure out which node prefixes (i.e. gene nodes) we need to select for
    valid_node_prefixes = {g.split(":")[0]:'' for g in genes_to_taxa.keys()}
    #valid_node_prefixes.update({"MONDO":''})
    #valid_node_prefixes.update({"ORPHA":''})
    #valid_node_prefixes.update({"HGNC":''})
    
    
    # Open monarch kg edges file, and pull out has_phenotype predicates.
    # These types of relationships occur between gene and phenotype nodes
    valid_preds = {"biolink:has_phenptype":''}
    with open(edges_filepath, 'r') as infile:
        
        # Read header line and creat data structure to map column names back to column indices
        h = {v:i for i,v in enumerate(infile.readline().strip('\r').strip('\n').split('\t'))}
        
        # Read through file and fill in taxon specific gene-->phenotype networks
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            pred = cols[h["predicate"]]

            # Gene has phenotype edges
            if pred == 'biolink:has_phenotype':
                
                # Grab gene and phenotype ids
                g_id = cols[h["subject"]] 
                p_id = cols[h["object"]]
                node_pref = g_id.split(":")[0]
                
                # Select for gene nodes only
                if node_pref in valid_node_prefixes:
                
                    # Link gene back to taxa obj and update its network (taxa-->gene-->phenetype)
                    if g_id not in taxa_objs[genes_to_taxa[g_id].ID].gene_to_phenotype:
                        taxa_objs[genes_to_taxa[g_id].ID].gene_to_phenotype.update({g_id:{}})
                    taxa_objs[genes_to_taxa[g_id].ID].gene_to_phenotype[g_id].update({p_id:''})
    
    # Tally up how many relationships we made and add them to our table
    g_with_p, tot_g_to_p, uniq_p = [], [], []
    for t_id in taxa_table["Taxon ID"]:
        tg2p = sum([len(v) for v in taxa_objs[t_id].gene_to_phenotype.values()])
        u_p = len(set([k for v in taxa_objs[t_id].gene_to_phenotype.values() for k in v]))
        g_with_p.append(len(taxa_objs[t_id].gene_to_phenotype))
        tot_g_to_p.append(tg2p)
        uniq_p.append(u_p)
    
    taxa_table.update({"Genes >= 1 Phenotype":g_with_p})
    taxa_table.update({"Total Phenotype Edges":tot_g_to_p})
    taxa_table.update({"Unique Phenotypes":uniq_p})
    
    return taxa_objs, taxa_table


def generate_orthologs_mapping_file(edges_filepath: str,
                                    outfile_path: str,
                                    genes_to_taxa: dict, 
                                    phen_to_name: dict, 
                                    taxa_objs: dict, 
                                    taxa_table: dict):
    
    # Figure out which node prefixes (i.e. gene nodes) we need to select for
    valid_node_prefixes = {g.split(":")[0]:'' for g in genes_to_taxa.keys()}
    
    # Open monarch kg edges file, and pull out has_phenotype predicates.
    # These types of relationships occur between gene and phenotype nodes
    valid_preds = {"biolink:has_phenptype":''}
    ortho_data = []
    with open(edges_filepath, 'r') as infile:
        
        # Read header line and creat data structure to map column names back to column indices
        h = {v:i for i,v in enumerate(infile.readline().strip('\r').strip('\n').split('\t'))}
        
        # Read through file and fill in taxon specific gene-->phenotype networks
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            pred = cols[h["predicate"]]

            # Gene has phenotype edges
            if pred == 'biolink:orthologous_to':
                
                # Grab gene and phenotype ids
                subj = cols[h["subject"]] 
                obj = cols[h["object"]]
                ortho_id = cols[h["has_evidence"]].split(":")[1]
                
                # Panther ortholog data is one directional. (i.e. SpeciesA-->SpeciesB)
                # We need both directions, so we swap subject and object fields
                ortho_data.append([subj, pred, obj, ortho_id])
                ortho_data.append([obj, pred, subj, ortho_id])
                
    print("- Ortho relationships created {} ({} if you divide by 2)".format(format(len(ortho_data), ','),
                                                                            format(len(ortho_data)/2, ',')))
    
    ortho_data = np.asarray(ortho_data).T
    pd.DataFrame({"geneA":ortho_data[0],
                  "predicate":ortho_data[1],
                  "geneB":ortho_data[2],
                  "ortholog_id":ortho_data[3]}).to_csv(outfile_path, sep='\t', index=False)
    
    print("- Orthologs converted to array of shape {}".format(ortho_data.shape))
    print("- Orthologs file written to... {}".format(outfile_path))


def generate_common_orthologs_files(orthologs_filepath, outpath_prefix, taxa_objs):

    ortho_df = pd.read_csv(orthologs_filepath, sep='\t')
    gA, gB, o_ids = list(ortho_df["geneA"]), list(ortho_df["geneB"]), list(ortho_df["ortholog_id"])
    taxa_ids = list(taxa_objs.keys())
    
    cc = 0
    for t_id in taxa_ids:
        pref_a = taxa_objs[t_id].gene_prefix
        for t_id2 in taxa_ids:
            if t_id == t_id2:
                continue

            # format our outfile path
            pref_b = taxa_objs[t_id2].gene_prefix
            outpath = "{}_{}_vs_{}.tsv".format(outpath_prefix,
                                               '-'.join(taxa_objs[t_id].label.split()), 
                                               '-'.join(taxa_objs[t_id2].label.split()))

            # Pull common orths and write
            common_orths = list(set([orth for a,b,orth in zip(gA, gB, o_ids) if (a.startswith(pref_a) and  
                                                                                 b.startswith(pref_b)) or 
                                                                                (a.startswith(pref_b) and  
                                                                                 b.startswith(pref_a))]))
            pd.DataFrame({"ortholog_id":common_orths}).to_csv(outpath, index=False)
            cc += 1
    
    print("- Common orthologs files written (A-->B and B-->A) {}".format(format(cc, ',')))

    
def generate_gene_and_phenotype_ortholog_files(orthologs_filepath, 
                                               outpath_dir, 
                                               gene_to_name, 
                                               phen_to_name,
                                               taxa_objs):
    
    # Read in ortholog mapping and create phenotype-->orthologs map
    ortho_df = pd.read_csv(orthologs_filepath, sep='\t')
    gA, gB, o_ids = list(ortho_df["geneA"]), list(ortho_df["geneB"]), list(ortho_df["ortholog_id"])
    taxa_ids = list(taxa_objs.keys())
    ortho_map = {g:o_id for g,o_id in zip(gA, o_ids)}
    
    # Loop through taxa and write gene_to_phenotype and phenotype_to_ortholog files
    cc = 0
    for t_id in taxa_ids:
        pref_a = taxa_objs[t_id].gene_prefix
        
        # Identify relevant gene_ids here and map them to ortholog_id
        rel_genes = set([ortho_map[v] for v,o_id in zip(gA, o_ids) if v.startswith(pref_a)])
        
        # Build phenotype-->gene/ortholog network here
        phen_to_orth = {}
        for k,v in taxa_objs[t_id].gene_to_phenotype.items():
            
            ortho_id = ortho_map.get(k, None)
            if not ortho_id:
                continue
            
            # Gene is not found in other speces 
            if ortho_id not in rel_genes:
                continue
            
            # Add this ortholog id to each phenotype entry
            for p_id in v:
                if p_id not in phen_to_orth:
                    phen_to_orth.update({p_id:[]})
                phen_to_orth[p_id].append(ortho_id)
        
        
        # Format our outfile path and write to file
        outpath_p2g = os.path.join(outpath_dir, 
                                   "{}_phenotype_to_ortholog.pkl".format("-".join(taxa_objs[t_id].label.split())))
        pickle.dump(phen_to_orth, open(outpath_p2g, "wb"))
        
        # Build data in table for and write via pandas
        outpath_g2p = os.path.join(outpath_dir, 
                                   "{}_gene_to_phenotype.tsv".format("-".join(taxa_objs[t_id].label.split())))
        
        taxa_pcount = sum([len(v) for k,v in taxa_objs[t_id].gene_to_phenotype.items()])
        genes = [k for k,v in taxa_objs[t_id].gene_to_phenotype.items() for kk in v]
        gnames = [gene_to_name[k] for k,v in taxa_objs[t_id].gene_to_phenotype.items() for kk in v]
        preds = ["biolink:has_phenotype" for i in range(0, taxa_pcount)]
        phens = [kk for k,v in taxa_objs[t_id].gene_to_phenotype.items() for kk in v]
        pnames = [phen_to_name[kk] for k,v in taxa_objs[t_id].gene_to_phenotype.items() for kk in v]
        
        pd.DataFrame({"gene":genes,
                      "gene_name":gnames,
                      "predicate":preds,
                      "phenotype":phens,
                      "phenotype_name":pnames}).to_csv(outpath_g2p, sep='\t', index=False)
        
        print("- Data written to {}".format(outpath_p2g))
        print("- Data written to {}".format(outpath_g2p))



if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Generates necessary files for all pairwaise species \
                                                      phenologs comparisons within the monarch kg. Gene to phenotype \
                                                      and phenotype to ortholog files are created along with a file \
                                                      that maps all gene ids to an ortholog id. A file containing the \
                                                      set of common orthologs (between any two species) is also \
                                                      written for each species pairwise comparison.')

        parser.add_argument("-o", "--out_dir", help="Directory to write files to", required=True, type=str)
        parser.add_argument("-n", "--nodes_file", help="File path to monarch kg nodes file", required=True, type=str, default=None)
        parser.add_argument("-e", "--edges_file", help="File path to monarch kg edges file", required=True, type=str, default=None)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Create initial taxa objects and "mapping tables" from monarch nodes file
    taxa_objs, genes_to_taxa, genes_to_name, phen_to_name, taxa_table = pull_taxon_gene_information(args.nodes_file)

    # Fill in taxa object gene_to_phenotype dictionary datastructure with relevant edges from monarch kg edge file
    taxa_objs, taxa_table = fill_in_gene_to_phenotype_networks(args.edges_file, genes_to_taxa, taxa_objs, taxa_table)

    # Write out orthologos mappings of genes across species. A-->B, AND B-->A
    ortho_map_file = os.path.join(args.out_dir, "panther_orthologs.tsv")
    generate_orthologs_mapping_file(edges_filepath=args.edges_file,
                                    outfile_path=ortho_map_file,
                                    genes_to_taxa=genes_to_taxa, 
                                    phen_to_name=phen_to_name, 
                                    taxa_objs=taxa_objs, 
                                    taxa_table=taxa_table)

    # Generate common orthologs files between all pairwise species comparisons
    common_ortho_prefix = os.path.join(args.out_dir, "common_orthologs")
    generate_common_orthologs_files(ortho_map_file, common_ortho_prefix, taxa_objs)

    # Generate gene-->phenotype, and phenotype-->ortholog_ids files (tsv and pk respectivly)
    generate_gene_and_phenotype_ortholog_files(ortho_map_file, 
                                               args.out_dir, 
                                               genes_to_name, 
                                               phen_to_name,
                                               taxa_objs)

    display(pd.DataFrame(taxa_table))