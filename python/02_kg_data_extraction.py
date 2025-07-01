# General imports
import os
import argparse
import gzip
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from collections import Counter


# Functions to read phenio-relation-graph.gz file and build directed graphs for each ontology
def read_phenio_to_ontologies(phenio_filepath, relevant_ontologies):
    """
    Reads through phenio-relation-graph.gz file line by line,
    and builds a directed graph for each ontology within the relevant_ontolgies argument.
    Each ontology graph only contains nodes specific to said ontology, 
    and only edges of the 'rdfs:subClassOf'. Brings in many edges for each ontology, but 
    because we only care about the descendents of particular terms, this is not an issue as those relationships
    will still hold true.
    """
    
    # Initiate an empty directed graph for each ontology
    onto_graphs = {k:nx.DiGraph() for k in relevant_ontologies}

    with gzip.open(phenio_filepath, 'rt') as infile: # rt is for read text
        
        # No header line for this file
        
        # Read through file and fill in taxon specific gene-->phenotype networks
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            e1, pred, e2 = cols # Expand into respective variables
            
            # We only want rdfs:subClassOf
            if pred != "rdfs:subClassOf":
                continue
        
            # Ensure both nodes of edge exist in graph
            eid = (e1, e2)
            
            onto1 = e1.split(":")[0]
            onto2 = e2.split(":")[0]
            
            if (onto1 == onto2) and (onto1 in onto_graphs):
                onto_graphs[onto1].add_edge(e2, e1)
    
    
    for k,v in onto_graphs.items():
        print("- Ontology {} pulled from phenio with {} nodes, and {} edges...".format(k, 
                                                                                          v.number_of_nodes(), 
                                                                                          v.number_of_edges()))
              
    return onto_graphs


def merge_ontology_nodes_to_exclusion_terms(ontology_graphs, ontology_terms):
    """
    The goal of this function is to create a set of nodes (consiting of ontology terms) that we need to 
    exclude from our analysis. The reason being is we are only interested in "abnormal" phenotypes,
    and there are wild type, and other "normal" phenotypes and or mode of inheritence phenotypes that do not
    make sense to include given the goal of the overarching analysis.
    In the original phenologs paper (https://www.pnas.org/doi/10.1073/pnas.0910200107) they only 
    use "mutational phenotypes between different organisms". We can then take the non-abnormal phenotype nodes 
    and avoid them during the construction of the subgraph of relevant data in the following functions.
    
    Depending on the ontology, we will need to filter one of two ways to get the nodes that we need to exclude.
    We either just pull all ancestors + orignal term, or we take all other nodes in the ontology not belonging to
    the parent term and filter those out. 
    
    ontology_graphs = {ontology_name:nx.DiGraph()}
    ontology_terms = {ontology_name:[parent_term, "Filter"/"Include"], ...}
    """
    
    nodes_to_exclude = set()
    for onto_name, onto_graph in ontology_graphs.items():
        parent_term = ontology_terms[onto_name][0]
        include = ontology_terms[onto_name][1]
        
        # Means all terms under parent term (and parent term) need to be avoided in the future when building graph.
        # So we grab all descendants of this node
        if include == "Filter":
            select_nodes = set(nx.descendants(onto_graph, parent_term)) | set([parent_term])
        
        # Means we need to include these nodes in the future. 
        # So we need to grab all other nodes in the ontology so that we know to avoid them in the future.
        elif include == "Include":
            select_nodes = set(onto_graph.nodes()) - set(nx.descendants(onto_graph, parent_term))
            select_nodes = select_nodes | set([parent_term])
        
        else:
            raise ValueError("- ERROR, only 'Filter' or 'Include' allowed... {} found".format(include))
        
        
        # Edge case where there is only one term we need to exclude. In this case, we simply use the input term
        if len(select_nodes) == 0:
            select_nodes = set([parent_term])
        
        # Add the nodes that we need to avoid in the future to our output datastructure
        nodes_to_exclude = nodes_to_exclude | select_nodes
        print("- {} nodes added to filtering {}/{} ".format(onto_name, len(select_nodes), len(onto_graph)))
        
    print("- Total nodes found from ontologies to avoid in the future {}".format(format(len(nodes_to_exclude), ',')))
    return nodes_to_exclude
 

# Step 1)
def initiate_graph_with_nodes(nodes_filepath, 
                              node_colname="id", 
                              node_cat_colname="category", 
                              node_type_colname="type",
                              filter_nodes={}):
    
    # Our valid node types here
    valid_node_types = {"biolink:Gene":'',
                        "biolink:PhenotypicFeature":'',
                        "biolink:Disease":''}
    
    # Protein coding from http://www.sequenceontology.org/browser/current_svn/term/SO:0001217
    ###protein_coding_term = "SO:0001217" # NOTE, this does not work for all species within the graph.
    ### But I am leaving this here for future reference if we are able to get this information incorporated for all species
    ### within the graph 
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
        node_cat = cc[col_inds[node_cat_colname]]
        node_type = cc[col_inds[node_type_colname]]
        
        # Filter for relevant node types
        if node_cat not in valid_node_types:
            continue

        # Check if our node_id is part of an input set that we DO NOT want to include
        if node_id in filter_nodes:
            continue
        
        # Check if node is a protein coding gene.
        # Note, this only works for a handful of species within the graph.
        # If / when the SO term is incorportated across all genes for all species, we can leverage this methodology
        ###if node_cat == "biolink:Gene":
        ###    if node_type != protein_coding_term:
        ###        continue
            
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


# Step 5)
def generate_ortho_edges_file(input_graph, outfile):
    
    # Output datastructure
    ortholog_edges = {"taxon_a":[],
                      "gene_a":[],
                      "gene_a_name":[],
                      "taxon_b":[],
                      "gene_b":[],
                      "gene_b_name":[],
                      "panther_id":[]}

    rel_edges = [e for e in input_graph.edges() if input_graph.edges[e]["predicate"] == "biolink:orthologous_to"]
    for e in rel_edges:
        # Gene ids are node ids, and gene names are nodes attribute "name"
        n1, n2 = e[0], e[1]
        taxon_a = input_graph.nodes[n1]["in_taxon"]
        taxon_b = input_graph.nodes[n2]["in_taxon"]
        gene_a_name = input_graph.nodes[n1]["name"]
        gene_b_name = input_graph.nodes[n2]["name"]
        panth_id = input_graph.edges[e]["has_evidence"].replace("PANTHER.FAMILY:", "")

        ortholog_edges["taxon_a"].append(taxon_a)
        ortholog_edges["gene_a"].append(n1)
        ortholog_edges["gene_a_name"].append(gene_a_name)
        ortholog_edges["taxon_b"].append(taxon_b)
        ortholog_edges["gene_b"].append(n2)
        ortholog_edges["gene_b_name"].append(gene_b_name)
        ortholog_edges["panther_id"].append(panth_id)
    
    # Write out ortholog edges file
    pd.DataFrame(ortholog_edges).to_csv(outfile, sep='\t', index=False)
    print("- Orthologs edges file written with {} edges...".format(format(len(ortholog_edges["taxon_a"]), ',')))


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

    # List out our ontology specific filtering criteria in the form of parent nodes/classes that each node from an ontology
    onto_parent_classes = {'MONDO':['MONDO:0700096', 'Include'],      ### Human disease (we WANT these terms)
                           'HP':['HP:0000118', 'Include'],            ### Phenotypic abnormality (we WANT these terms)
                           'FYPO':['FYPO:0001985', 'Include'],        ### Abnormal phenotype (we WANT these terms)
                           'MP':['MP:0002873', 'Filter'],             ### Normal phenotype (we do NOT want these terms)
                           'DDPHENO':['DDPHENO:0000142', 'Filter']}   ### Wild type (we do NOT want these terms)

    # Input files
    nodes_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_nodes.tsv")
    edges_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    phenio_file = os.path.join(args.project_dir, "monarch_kg", "phenio-relation-graph.gz")

    # Output files / info
    species_data_dir = os.path.join(args.project_dir, "species_data")
    species_info_table_file = os.path.join(species_data_dir, "species_information_table.tsv")
    ortholog_edges_file = os.path.join(species_data_dir, "ortholog_edges.tsv")

    # Read through phenio-relation-graph to pull in relevant ontologies so we can filter out nodes downstream
    phenio_onto_graphs = read_phenio_to_ontologies(phenio_file, relevant_ontologies=onto_parent_classes)
    nodes_to_exclude = merge_ontology_nodes_to_exclusion_terms(phenio_onto_graphs, onto_parent_classes)

    # Run data extraction pipeline (create nx graph by selecting for specific node and edge types)
    graph, taxon_names = initiate_graph_with_nodes(nodes_filepath=nodes_file, 
                                                   node_colname="id", 
                                                   node_cat_colname="category", 
                                                   node_type_colname="type",
                                                   filter_nodes=nodes_to_exclude)

    graph = fill_in_graph_with_edges(edges_file, graph)
    taxon_stats = compute_association_tables(graph, species_data_dir)
    generate_common_orthologs_files(tables_dir=species_data_dir, taxon_stats=taxon_stats)
    
    # Write taxon stats table to file
    pd.DataFrame(taxon_stats).to_csv(species_info_table_file, sep='\t', index=False)

    # Write orthologs edges file
    generate_ortho_edges_file(graph, outfile=ortholog_edges_file)

    print("- Species information table (used in downstream steps) written to {}...".format(species_info_table_file))
    print("- Data extraction pipeline complete...")