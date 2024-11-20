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
    genes: Optional[Dict] = {}
        

def pull_taxon_gene_information(nodes_filepath):

    # Our taxon information comes from gene nodes
    valid_nodes = {"biolink:Gene":'',
                   "biolink:Genotype":'',
                   "biolink:PhenotypicFeature":''}
    
    taxa_node_prefix = "NCBITaxon:"
    
    taxa_objs = {}
    genes_to_taxa = {}
    genes_to_name = {}
    phenotype_to_name = {}
    
    # Read data into memory via pandas and fill missing values with an empty string (instead of NaN)
    df = pd.read_csv(nodes_filepath, sep='\t', low_memory=False).fillna('')
    print("- Nodes file read into memory...")
    
    # Create NCBITaxon --> name mapping dict 
    # (Not all nodes will have a corresponding taxa name so we can store that information here)
    taxon_to_name = {n_id:n_name for n_id, n_name in zip(list(df["id"]), 
                                                         list(df["name"])) if n_id.startswith(taxa_node_prefix)}
    print("- NCBITaxon id nodes --> name table created of size {}...".format(format(len(taxon_to_name), ',')))
    print("- Parsing node information...")
    # Filter out irrelivant node types, and fill in each taxa information that we have gene_ids for
    for n_id, n_cat, n_name, n_tax, n_tax_lab in zip(list(df["id"]), 
                                                     list(df["category"]),
                                                     list(df["name"]),
                                                     list(df["in_taxon"]), 
                                                     list(df["in_taxon_label"])):
        
        
        # Irrelivant node_types
        if n_cat not in valid_nodes:
            continue
        
        # We cannot link phenotypic feature nodes to taxa (until we read in edges)
        # But we can build a table that maps the id to a name 
        if n_cat == "biolink:PhenotypicFeature":
            phenotype_to_name.update({n_id:n_name})
            continue
        
        # First encounter of xyz species. Need to update data structures
        if (n_tax not in taxa_objs) and (n_tax != ''):
            
            # Leverage lookup table we made earlier. This deals with rows that don't have taxon labels
            lk_lab = taxon_to_name.get(n_tax, n_tax_lab)
            if lk_lab != "":
                n_tax_lab = lk_lab
            
            # Replace a non-existant label with the taxon_id instead
            n_tax_lab = n_tax if n_tax_lab == "" else n_tax_lab
            
            # Initiate Taxa obj config and add new taxa to our data structure
            g_prefix = n_id.split(":")[0]
            tconfig = {"ID":n_tax, 
                       "label":n_tax_lab,
                       "gene_prefix":g_prefix}
            
            taxa_objs.update({n_tax:Taxa.parse_obj(tconfig)})
        
        # Link our gene_id to its respective taxon (NOTE - We only link genes/genotypes that have a taxon id)
        if n_tax != '':
            taxa_objs[n_tax].genes.update({n_id:''})
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
    
    print("- Initial dataset table created from nodes file...")
    return taxa_objs, genes_to_taxa, genes_to_name, phenotype_to_name, taxon_table


def fill_in_gene_to_phenotype_networks(edges_filepath, genes_to_taxa, taxa_objs, taxa_table):
    
    # Figure out which node prefixes (i.e. gene nodes) we need to select for
    valid_node_prefixes = {g.split(":")[0]:'' for g in genes_to_taxa.keys()}
    
    # TO DO: Need to figure out how to deal with these
    ##valid_node_prefixes.update({"MONDO":''})
    ##valid_node_prefixes.update({"ORPHA":''})
    ##valid_node_prefixes.update({"HGNC":''})
    
    
    # Open monarch kg edges file, and pull out has_phenotype predicates.
    # These types of relationships occur between gene and phenotype nodes
    valid_preds = {"biolink:has_phenptype":''}
    uncertain_nodes_count = 0
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
    phen_prefixes = []
    for t_id in taxa_table["Taxon ID"]:
        tg2p = sum([len(v) for v in taxa_objs[t_id].gene_to_phenotype.values()])
        u_p = len(set([k for v in taxa_objs[t_id].gene_to_phenotype.values() for k in v]))
        g_with_p.append(len(taxa_objs[t_id].gene_to_phenotype))
        tot_g_to_p.append(tg2p)
        uniq_p.append(u_p)
        
        # Phenotype prefixes per species (some can have multiple)
        phen_prefs = list({vv.split(':')[0]:'' for v in taxa_objs[t_id].gene_to_phenotype.values() for vv in v}.keys())
        phen_prefixes.append(",".join(phen_prefs))
        
    taxa_table.update({"Phenotype Prefix":phen_prefixes})
    taxa_table.update({"Genes >= 1 Phenotype":g_with_p})
    taxa_table.update({"Total Phenotype Edges":tot_g_to_p})
    taxa_table.update({"Unique Phenotypes":uniq_p})
    
    return taxa_objs, taxa_table


def generate_orthologs_mapping_file(edges_filepath: str,
                                    outfile_path: str,
                                    genes_to_taxa: dict, 
                                    phen_to_name: dict, 
                                    taxa_objs: dict,
                                    taxa_table):
    
    # Build mapping tables to compute ortholog counts per taxa
    genes_to_taxa_label = {g:v.label for g, v in genes_to_taxa.items()}
    taxa_orth_counts = {taxa_obj.label:0 for taxa_obj in taxa_objs.values()}
    
    # Figure out which node prefixes (i.e. gene nodes) we need to select for
    valid_node_prefixes = {g.split(":")[0]:'' for g in genes_to_taxa.keys()}
    
    # Open monarch kg edges file, and pull out has_phenotype predicates.
    # These types of relationships occur between gene and phenotype nodes
    valid_preds = {"biolink:has_phenptype":''}
    ortho_data = []
    with open(edges_filepath, 'r') as infile:
        
        # Read header line and create data structure to map column names back to column indices
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

                # Tally orth taxa counts
                taxa_orth_counts[genes_to_taxa_label[subj]] += 1
                taxa_orth_counts[genes_to_taxa_label[obj]] += 1
                
    print("- Ortho relationships created {} ({} if you divide by 2)".format(format(len(ortho_data), ','),
                                                                            format(len(ortho_data)/2, ',')))
    
    # Add ortholog mapping counts to our dataframe (on a per taxa basis)
    txa_counts = [taxa_orth_counts[t_label] for t_label in list(taxa_table["Taxon Label"])]
    taxa_table["Ortholog Mapping Count"] = txa_counts
    
    ortho_data = np.asarray(ortho_data).T
    pd.DataFrame({"geneA":ortho_data[0],
                  "predicate":ortho_data[1],
                  "geneB":ortho_data[2],
                  "ortholog_id":ortho_data[3]}).to_csv(outfile_path, sep='\t', index=False)
    
    print("- Orthologs converted to array of shape {}".format(ortho_data.shape))
    print("- Orthologs file written to... {}".format(outfile_path))
    return taxa_table


def generate_common_orthologs_files(orthologs_filepath, outpath_prefix, taxa_objs, ortho_counts_outpath):
    """
    In order for a species to be included in our dataset, there must exist some degree of orthology between
    other species within the knowledge graph. This orthology is derived from panther
    """

    ortho_df = pd.read_csv(orthologs_filepath, sep='\t')
    gA, gB, o_ids = list(ortho_df["geneA"]), list(ortho_df["geneB"]), list(ortho_df["ortholog_id"])
    taxa_ids = list(taxa_objs.keys())
    
    cc = 0
    comps_count_table = []
    for t_id in taxa_ids:
        pref_a = taxa_objs[t_id].gene_prefix
        genes_a = taxa_objs[t_id].genes
        
        for t_id2 in taxa_ids:
            if t_id == t_id2:
                continue
            
            # format our outfile path
            pref_b = taxa_objs[t_id2].gene_prefix
            genes_b = taxa_objs[t_id2].genes
            outpath = "{}_{}_vs_{}.tsv".format(outpath_prefix,
                                               '-'.join(taxa_objs[t_id].label.split()), 
                                               '-'.join(taxa_objs[t_id2].label.split()))

            # Pull common orths and write 
            # First deal with edge case where both species use the same gene id space (Xenopus for example)
            # TO DO: Is this right? Or should we be comparing globally instead of just genes with phenotypes
            if pref_a == pref_b:
                common_orths = []###list(gene_keys_a & gene_keys_b)

            else:
                # Gene A-->B or vice versa
                common_orths = list(set([orth for a,b,orth in zip(gA, gB, o_ids) if (a in genes_a and b in genes_b) or
                                                                                    (b in genes_a and a in genes_b)]))
                
                comps_count_table.append([t_id, t_id2, str(len(common_orths))])
 
            # Covert to df and write                                                     
            pd.DataFrame({"ortholog_id":common_orths}).to_csv(outpath, index=False)
            cc += 1
    
    # Convert our orthologs count file to nparray-->dataframe and write
    comps_count_table = np.asarray(comps_count_table).T
    pd.DataFrame({"Species_A":comps_count_table[0],
                  "Species_B":comps_count_table[1],
                  "Common_Orthologs":comps_count_table[2]}).to_csv(ortho_counts_outpath, sep='\t', index=False)

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
        
        print("- Processing taxa {}, {} gene/phen/ortholog mappings".format(t_id, taxa_objs[t_id].label))
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
        
        # Write our data via pandas
        pd.DataFrame({"gene":genes,
                      "gene_name":gnames,
                      "predicate":preds,
                      "phenotype":phens,
                      "phenotype_name":pnames}).to_csv(outpath_g2p, sep='\t', index=False)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Generates necessary files for all pairwaise species \
                                                      phenologs comparisons within the monarch kg. Gene to phenotype \
                                                      and phenotype to ortholog files are created along with a file \
                                                      that maps all gene ids to an ortholog id. A file containing the \
                                                      set of common orthologs (between any two species) is \
                                                      written for each species pairwise comparison. \
                                                      And a table detailing the number of common orthologs shared between any two species is also \
                                                      written. Species within the monarch-kg that do not have any phenotype edges or \
                                                      do not have any orthologous mappings then they are excluded from the final output table.')

        parser.add_argument("-p", "--project_dir", help="Phenologs top level project directory", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Input files
    nodes_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_nodes.tsv")
    edges_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")

    # Output files / info
    species_data_dir = os.path.join(args.project_dir, "species_data")
    ortho_map_file = os.path.join(species_data_dir, "panther_orthologs.tsv")
    common_ortho_prefix = os.path.join(species_data_dir, "common_orthologs")
    species_info_table_file = os.path.join(species_data_dir, "species_information_table.tsv")
    species_info_table_file_unfiltered = os.path.join(species_data_dir, "species_information_table_unfiltered.tsv")
    common_ortho_counts_file = os.path.join(species_data_dir, "species_common_orthologs_counts.tsv")


    # Create initial taxa objects and "mapping tables" from monarch nodes file
    taxa_objs, genes_to_taxa, genes_to_name, phen_to_name, taxa_table = pull_taxon_gene_information(nodes_file)

    # Fill in taxa object gene_to_phenotype dictionary datastructure with relevant edges from monarch kg edge file
    taxa_objs, taxa_table = fill_in_gene_to_phenotype_networks(edges_file, genes_to_taxa, taxa_objs, taxa_table)

    # Write out orthologos mappings of genes across species. A-->B, AND B-->A
    taxa_table = generate_orthologs_mapping_file(edges_filepath=edges_file,
                                                 outfile_path=ortho_map_file,
                                                 genes_to_taxa=genes_to_taxa, 
                                                 phen_to_name=phen_to_name, 
                                                 taxa_objs=taxa_objs,
                                                 taxa_table=taxa_table)
    
    # Convert to df and write unfiltered table (for record keeping)
    taxa_table = pd.DataFrame(taxa_table)
    taxa_table.to_csv(species_info_table_file_unfiltered, sep='\t', index=False)

    # Filter taxa_table based on phenotype edges and orthologous gene mappings and write file (for downstream analysis)
    taxa_table = taxa_table[(taxa_table["Total Phenotype Edges"].astype(int) > 0) & 
                            (taxa_table["Ortholog Mapping Count"].astype(int) > 0)]
    taxa_table.to_csv(species_info_table_file, sep='\t', index=False)

    # Filter out taxa that do not have the necessary information, so we don't generate excess / useless files
    relevant_taxa_ids = set(taxa_table["Taxon ID"])
    relevant_taxa_objs = {t:v for t,v in taxa_objs.items() if t in relevant_taxa_ids}

    # Generate common orthologs files between all pairwise species comparisons
    generate_common_orthologs_files(ortho_map_file, 
                                    common_ortho_prefix, 
                                    relevant_taxa_objs,
                                    ortho_counts_outpath=common_ortho_counts_file)

    # Generate gene-->phenotype, and phenotype-->ortholog_ids files (tsv and pk respectivly)
    generate_gene_and_phenotype_ortholog_files(ortho_map_file, 
                                               species_data_dir, 
                                               genes_to_name, 
                                               phen_to_name,
                                               relevant_taxa_objs)
    
    print("- Species information table (used in downstream steps) written to {}...".format(species_info_table_file))
    print("- Data extraction pipeline complete...")