import os
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle

from collections import Counter
from pydantic import BaseModel
from typing import Optional, Any, List, Dict, Union

from IPython.display import display


# General function to read an .obo ontontolgy file into memory using pronto to gather all terms that do not fall under a particular parent class
def read_ontology_to_accepted_terms(ontology_obo_file, filter_term="HP:0000118", include=True):
    
    # Read ontology file into memory
    onto = Ontology(ontology_obo_file)
    exclude_terms = {}
    term_count = len(list(onto.terms()))
    
    for term in onto.terms():
        
        # Gather ancestor terms and update our filtering datastructure accordingly
        parent_terms = {ancestor.id: ancestor.name for ancestor in term.superclasses()}
        if filter_term not in parent_terms:
            exclude_terms.update({term.id:term.name})

    
    print("- Terms from ontology found that do not belong to parent class {} {}/{}".format(filter_term, 
                                                                                           format(len(exclude_terms)), 
                                                                                           format(term_count)))
    return exclude_terms


# Step 1)
def initiate_graph_with_nodes(nodes_filepath, node_colname="id", node_type_colname="category"):
    
    # Our valid node types here
    valid_node_types = {"biolink:Gene":'',
                        "biolink:PhenotypicFeature":'',
                        "biolink:Disease":''}
    
    taxa_node_prefix = "NCBITaxon:"
    taxons = {}

    # Read data into memory via pandas and fill missing values with an empty string (instead of NaN)
    nodes_df = pd.read_csv(nodes_filepath, sep='\t', low_memory=False).fillna('')
    print("- Nodes file read into memory...")
    
    taxon_to_name = {n_id:n_name for n_id, n_name in zip(list(nodes_df["id"]), 
                                                         list(nodes_df["name"])) if n_id.startswith(taxa_node_prefix)}
    
    cols = list(nodes_df.columns)
    col_inds = {c:i for i,c in enumerate(cols)}
    
    
    # Initiate nx graph
    graph = nx.Graph()
    
    attr_count = 0
    repeat_nodes = {}
    for cc in zip(*[list(nodes_df[c]) for c in cols]):
        
        # Add node based on column index that maps to node_colname
        # and keep track of repeate nodes
        node_id = cc[col_inds[node_colname]]
        node_cat = cc[col_inds[node_type_colname]]
        
        # Filter for relevant node types
        if node_cat not in valid_node_types:
            continue
            
        # Check for repeate nodes
        if graph.has_node(node_id):
            
            if node_id not in repeat_nodes:
                repeat_nodes.update({node_id:0})
            repeat_nodes[node_id] += 1
        
        else:
            graph.add_node(node_id)
            
            # Add attributes
            attributes = {node_id:{cname:cval for cname, cval in zip(cols, cc) if not pd.isna(cval)}}
            nx.set_node_attributes(graph, attributes)
            attr_count += len(attributes[node_id])
    
    print("- Initial graph created with {} nodes of types:".format(format(len(list(graph.nodes())), ',')))
    for k,v in Counter([graph.nodes[n]['category'] for n in graph.nodes()]).items():
        print("  - {} {}".format(k, format(v, ',')))
    
    print("- Repeat nodes found {}...".format(format(len(repeat_nodes), ',')))
    return graph, taxon_to_name


# Step 2)
def fill_in_graph_with_edges(edges_filepath, graph):
    
    uncertain_nodes_count = 0
    attr_count = 0
    repeat_edges = {}
    
    valid_edges = set(["biolink:DiseaseToPhenotypicFeatureAssociation",
                       "biolink:GeneToPhenotypicFeatureAssociation",
                       "biolink:CausalGeneToDiseaseAssociation",
                       "biolink:CorrelatedGeneToDiseaseAssociation",
                       "biolink:GeneToGeneHomologyAssociation"])
    
    with open(edges_filepath, 'r') as infile:
        
        # Read header line and create data structure to map column names back to column indices
        h = {v:i for i,v in enumerate(infile.readline().strip('\r').strip('\n').split('\t'))}
        
        # Read through file and fill in taxon specific gene-->phenotype networks
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            pred = cols[h["predicate"]]
            cat = cols[h['category']]
            
            # We only care about certain edge types "Association" types from biolink
            if cat not in valid_edges:
                continue
        
            # Ensure both nodes of edge exist in graph
            e1, e2 = (cols[h["subject"]], cols[h["object"]])
            eid = (e1, e2)
            if (not graph.has_node(e1)) or (not graph.has_node(e2)):
                continue

            # Check for potential repeat edge
            elif graph.has_edge(e1, e2):
                if (e1, e2) not in repeat_edges and (e2, e1):
                    repeat_edges.update({eid:0})
                repeat_edges[eid] += 1

            # Panther is the ONLY source of orthology we want to use for orthologous edges
            elif (pred == 'biolink:orthologous_to') and ("PANTHER.FAMILY:" not in cols[h["has_evidence"]]):
                continue
            
            # Add edge along with attributes
            else:
                
                graph.add_edge(e1, e2)
                attributes = {eid:{cname:cols[cindex] for cname, cindex in h.items() if not pd.isna(cols[cindex])}}
                nx.set_edge_attributes(graph, attributes)
                attr_count += len(attributes[eid])
    
    
    print("- Successfully added {} edges to graph of predicate types:".format(format(graph.number_of_edges(), ',')))
    for k,v in Counter([graph.edges[e]["predicate"] for e in graph.edges()]).items():
        print("  - {} {}".format(k, format(v, ',')))
    return graph


# Step 3 repeat call function
def compute_taxon_info(input_graph, taxon_id: str, taxon_nodes: list, outpath_dir: str):
    
    # Select for gene-->disease or gene-->phenotype edges for specific taxon
    gene_to_phens_lines = []
    gene_to_disease_lines = []
    g2p_edge_count, g2g_edge_count, g2d_edge_count = 0, 0, 0
    
    # Relevant mapping tables for relationships we want
    g2p = {n:[] for n in taxon_nodes} # gene-->phenotype
    g2d = {n:[] for n in taxon_nodes} # gene-->disease
    
    g2o = {} # gene-->ortholog_id
    p2o = {} # phenotype-->ortholog_id
    d2o = {} # disease-->ortholog_id
    gtype = "biolink:Gene"
    
    # Grab taxon name / label / gene prefix
    tx_id = input_graph.nodes[taxon_nodes[0]]["in_taxon"]
    tx_name = input_graph.nodes[taxon_nodes[0]]["in_taxon_label"]
    tx_gprefix = taxon_nodes[0].split(":")[0] # Gene namespace prefix
    tx_pprefix = "" # Phenotype ontology namespace prefix (Needs to be generated from a phenotype neighbor node)
    
    if tx_name == "": # Deal with empty name
        tx_name = tx_id
    
    # Gather up relevant information pertaining to this taxa
    for node_id in taxon_nodes:
        
        # Pull up neighbor nodes
        neighbs = list(input_graph.neighbors(node_id))
        
        # First, Grab panther protein family id (if possible)
        orths = [neib for neib in neighbs if input_graph.nodes[neib]["category"] == "biolink:Gene"]
        if len(orths) > 0:
            orth_id = input_graph.edges[(node_id, orths[0])]["has_evidence"].replace("PANTHER.FAMILY:", "")
        else:
            orth_id = ""
        
        # Update our ortholog mappings
        if orth_id != "":
            g2o.update({node_id:orth_id})
            
        # Loop through neighboring nodes and deal with the different associations we want to tally
        for neib in neighbs:
            
            nn = input_graph.nodes[neib]
            nn_cat = nn["category"]
            
            # Gene-->Phenotype relationships
            if nn_cat == "biolink:PhenotypicFeature":
                
                # Grab FIRST taxon phenotype ontology prefix we find for this species
                ppref = neib.split(":")[0]
                if tx_pprefix == "":
                    tx_pprefix = ppref
                elif tx_pprefix != ppref:
                    print("- Warning, taxon {} has multiple phenotype prefixes... {}, {}".format(tx_id, 
                                                                                                 tx_pprefix,
                                                                                                 ppref))
                
                # Output format should be the following
                #gene, gene_name, predicate, phenotype, phenotype_name, ortholog_id
                g2p_rels = [node_id, 
                            input_graph.nodes[node_id]["name"], 
                            input_graph.edges[(node_id, neib)]["predicate"],
                            neib,
                            nn["name"],
                            orth_id]
                
                # Update phenotype-->ortholog
                if orth_id != "":
                    if neib not in p2o:
                        p2o.update({neib:[]})
                    p2o[neib].append(orth_id)
                
                # Add data to relevant datastuctures
                g2p[node_id].append(neib)
                gene_to_phens_lines.append(g2p_rels)
                g2p_edge_count += 1
                
            
            # Gene-->Disease relationships
            elif nn_cat == "biolink:Disease":
                
                # Output format should be the following
                #gene, gene_name, predicate, disease, disease_name, ortholog_id
                g2d_rels = [node_id, 
                            input_graph.nodes[node_id]["name"], 
                            input_graph.edges[(node_id, neib)]["predicate"],
                            neib,
                            nn["name"],
                            orth_id]
                
                # Update phenotype-->ortholog
                if orth_id != "":
                    if neib not in d2o:
                        d2o.update({neib:[]})
                    d2o[neib].append(orth_id)
                
                # Add data to relevant datastuctures
                g2d[node_id].append(neib)
                gene_to_disease_lines.append(g2d_rels)
                g2d_edge_count += 1
            
            # Tally up orthology edges
            elif nn_cat == "biolink:Gene":
                g2g_edge_count += 1
            
            # Nothing else on the radar
            else:
                continue
    
    # Don't do anything if no phenotype data is present
    if (len(gene_to_phens_lines) == 0) or (len(g2o) == 0):
        print("- Skipping taxon id/name {} {}... Missing phenotype or orthology data".format(tx_id, tx_name))
        return {}

    # Generate outputfile paths
    tx_formatted = "-".join(tx_name.split())
    outpath_g2p = os.path.join(outpath_dir, "{}_gene_to_phenotype.tsv".format(tx_formatted))
    outpath_g2d = os.path.join(outpath_dir, "{}_gene_to_disease.tsv".format(tx_formatted))
    outpath_g2o = os.path.join(outpath_dir, "{}_gene_to_ortholog.tsv".format(tx_formatted))
    outpath_p2o = os.path.join(outpath_dir, "{}_phenotype_to_ortholog.pkl".format(tx_formatted))
    outpath_d2o = os.path.join(outpath_dir, "{}_disease_to_ortholog.pkl".format(tx_formatted))
    
    # Write gene-->phenotype file
    gene_to_phens_lines = np.asarray(gene_to_phens_lines).T
    pd.DataFrame({"gene":gene_to_phens_lines[0],
                  "gene_name":gene_to_phens_lines[1],
                  "predicate":gene_to_phens_lines[2],
                  "phenotype":gene_to_phens_lines[3],
                  "phenotype_name":gene_to_phens_lines[4],
                  "ortholog_id":gene_to_phens_lines[5]}).to_csv(outpath_g2p, sep='\t', index=False)
    
    # Not all species have gene-->disease associations
    if len(gene_to_disease_lines) != 0:
        
        # Write out gene-->disease
        gene_to_disease_lines = np.asarray(gene_to_disease_lines).T
        pd.DataFrame({"gene":gene_to_disease_lines[0],
                      "gene_name":gene_to_disease_lines[1],
                      "predicate":gene_to_disease_lines[2],
                      "disease":gene_to_disease_lines[3],
                      "disease_name":gene_to_disease_lines[4],
                      "ortholog_id":gene_to_disease_lines[5]}).to_csv(outpath_g2d, sep='\t', index=False)
    
    # Write phenotype-->ortholog .pkl file ({phenotype_id:[panther_id1, panther_id1, panther_id2, etc...]})
    pickle.dump(p2o, open(outpath_p2o, "wb"))

    # Write disease-->ortholog .pkl file ({phenotype_id:[panther_id1, panther_id1, panther_id2, etc...]})
    if len(d2o) > 0:
        pickle.dump(d2o, open(outpath_d2o, "wb"))

    # Write out gene-->ortholog_id
    genes_to_orths = np.asarray([[k,v] for k,v in g2o.items()]).T
    pd.DataFrame({"gene":genes_to_orths[0],
                  "ortholog_id":genes_to_orths[1]}).to_csv(outpath_g2o, sep='\t', index=False)
    
    # Tally up stats we just computed
    output_stats = {"Taxon ID":tx_id,
                    "Taxon label":tx_name,
                    "Taxon gene prefix":tx_gprefix,
                    "Taxon gene count":len(taxon_nodes),
                    
                    "Phenotype prefix":tx_pprefix,
                    "Genes with >= 1 phenotype edge":len([1 for v in g2p.values() if len(v) > 0]),
                    "Genes with >= 1 disease edge":len([1 for v in g2d.values() if len(v) > 0]),
                    "Total gene2phenotype edges":g2p_edge_count,
                    "Total gene2disease edges":g2d_edge_count,
                    "Unique phenotypes":len(p2o),
                    "Unique diseases":len(d2o),
                    "Ortholog mapping count":g2g_edge_count,
                    "Genes with >= 1 ortholog edge":len(g2o),
                    "Unique ortholog ids (orthogroups)":len({v:'' for v in g2o.values()})}
    
    return output_stats


# Step 3)
def compute_association_tables(input_graph, outpath_dir: str):
    
    # Output datastructure
    stats_table = {}
    
    # Subset out gene nodes. Apart from taxon and genotype nodes, 
    # gene nodes are the only nodes that have species level information (i.e. taxon)
    taxon_gene_nodes = {}
    for n in input_graph.nodes():
        if "biolink:Gene" == input_graph.nodes[n]["category"]:
            tx = input_graph.nodes[n]["in_taxon"]
            if tx not in taxon_gene_nodes:
                taxon_gene_nodes.update({tx:[]})
            taxon_gene_nodes[tx].append(n)
    
    # Process each taxon we identified from the gene nodes within the subgraph that we made
    for tx_id, tx_gene_nodes in taxon_gene_nodes.items():
        print("- Processing taxon {}".format(tx_id))
        tx_stats = compute_taxon_info(input_graph, tx_id, tx_gene_nodes, outpath_dir)
        
        if len(tx_stats) > 0:
            for k,v in tx_stats.items():
                if k not in stats_table:
                    stats_table.update({k:[]})
                stats_table[k].append(v)
    
    return stats_table


# Step 4)
def generate_common_orthologs_files(tables_dir, taxon_stats):
    """
    All input/output files will be read in from and written to the input argument tables_dir.
    
    Each specific gene-->ortholog mapping file is read into memory
    All pairwise comparisons are made between species to find the common set of orthologs shared between the two.
    
    Two files are made for each species to species comparison that are identical in content but only differ in name.
    For example...
    common_orthologs_Homo-sapiens_vs_Danio-rerio.tsv
    common_orthologs_Danio-rerio_vs_Homo-sapiens.tsv
    
    A single file is made that tallies the common ortholog counts across each species comparison.
    species_common_orthologs_counts.tsv
    """
    
    # Read ortholog files into memory
    tx_orths = {}
    for fname in os.listdir(tables_dir):
        if fname.endswith("_gene_to_ortholog.tsv"):
            tx = fname.replace("_gene_to_ortholog.tsv", "")
            fpath = os.path.join(tables_dir, fname)
            orth_ids = set(list(pd.read_csv(fpath, sep='\t')["ortholog_id"]))
            tx_orths.update({tx:orth_ids})
    
    # Write out all common orthologs between any two species (twice for naming directionality)
    cc = 0
    common_orths_data = []
    tx_ids = list(tx_orths.keys())
    for tx1 in tx_ids:
        for tx2 in tx_ids:
            if tx1 == tx2:
                continue
            
            # Write out species specific file
            comm_orths = list(tx_orths[tx1] & tx_orths[tx2])
            comm_orth_count = len(comm_orths)
            
            outpath = os.path.join(tables_dir, "common_orthologs_{}_vs_{}.tsv".format(tx1, tx2))
            pd.DataFrame({"ortholog_id":comm_orths}).to_csv(outpath, sep='\t', index=False)
            
            # Tally common ortholog counts between species a and species b
            common_orths_data.append([tx1, tx2, comm_orth_count])
            cc += 1
    
    # Lastly, write out common ortholog counts file
    outpath = os.path.join(tables_dir, "species_common_orthologs_counts.tsv")
    common_orths_data = np.asarray(common_orths_data).T
    pd.DataFrame({"Species_A":common_orths_data[0],
                  "Species_B":common_orths_data[1],
                  "Common_Orthologs":common_orths_data[2]}).to_csv(outpath, sep='\t', index=False)
    
    print("- Common orthologs files written (A-->B and B-->A) {}".format(format(cc, ',')))


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
    species_info_table_file = os.path.join(species_data_dir, "species_information_table.tsv")

    # Run data extraction pipeline (create nx graph by selecting for specific node and edge types)
    graph, taxon_names = initiate_graph_with_nodes(nodes_file, node_colname="id", node_type_colname="category")
    graph = fill_in_graph_with_edges(edges_file, graph)
    taxon_stats = compute_association_tables(graph, species_data_dir)
    generate_common_orthologs_files(tables_dir=species_data_dir, taxon_stats=taxon_stats)
    
    # Write taxon stats table to file
    pd.DataFrame(taxon_stats).to_csv(species_info_table_file, sep='\t', index=False)

    print("- Species information table (used in downstream steps) written to {}...".format(species_info_table_file))
    print("- Data extraction pipeline complete...")