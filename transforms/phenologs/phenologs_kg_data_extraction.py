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


class Taxa(BaseModel):
    label: Optional[str] = None
    ID: Optional[str] = None
    name: Optional[str] = None
    gene_prefix: Optional[str] = None
    phenotype_prefix: Optional[str] = None
    gene_to_ortholog_family: Optional[Dict] = {}
    gene_to_phenotype: Optional[Dict] = {}
    phenotype_ortholog: Optional[Dict] = {}
    genes: Optional[Dict] = {}
        

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
def pull_taxon_gene_information(nodes_filepath):

    # Our taxon information comes from gene nodes
    # TO DO: Add in "biolink:Genotype":'', 
    # We need / want to add in Genotype nodes which can connect to disease nodes. 
    # Currently we cannot connect Genotype nodes to gene nodes though, which is an issue that needs resolving still
    valid_nodes = {"biolink:Gene":'',
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
            
            taxa_objs.update({n_tax:Taxa.model_validate(tconfig)})
        
        # Link our gene_id to its respective taxon (NOTE - We only link genes/genotypes that have a taxon id)
        if (n_tax != '') and ("biolink:Gene" == n_cat):
            taxa_objs[n_tax].genes.update({n_id:''})
            genes_to_taxa.update({n_id:taxa_objs[n_tax]})
            genes_to_name.update({n_id:n_name})
            
            g_prefix = n_id.split(":")[0]
            if g_prefix != taxa_objs[n_tax].gene_prefix:
                print("- WARNING -- Multiple gene prefixes found for {}, {} and {}...".format(n_tax, 
                                                                                              taxa_objs[n_tax].gene_prefix))
            
    
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


# Step 2)
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
        
        # Read header line and create data structure to map column names back to column indices
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
                    txk = genes_to_taxa[g_id].ID
                    if g_id not in taxa_objs[txk].gene_to_phenotype:
                        taxa_objs[txk].gene_to_phenotype.update({g_id:{}})
                    taxa_objs[txk].gene_to_phenotype[g_id].update({p_id:''})
    
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


# Step 3)
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
    ortho_data = []
    non_panther_ortho_relationships = 0 # Keep track of other sources of orthology encountered beside panther...
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
                
                # Hard code for now, but we only want orthologous relationships that have a panther id
                panther_id = cols[h["has_evidence"]]
                if "PANTHER.FAMILY:" not in panther_id:
                    non_panther_ortho_relationships += 1
                    continue
                    
                ortho_id = panther_id.split(":")[1]
                
                # Panther ortholog data is one directional. (i.e. SpeciesA-->SpeciesB)
                # We need both directions, so we swap subject and object fields
                ortho_data.append([subj, pred, obj, ortho_id])
                ortho_data.append([obj, pred, subj, ortho_id])

                # Tally orth taxa counts
                taxa_orth_counts[genes_to_taxa_label[subj]] += 1
                taxa_orth_counts[genes_to_taxa_label[obj]] += 1

                # Update Taxa object gene-->ortholog mapping table
                # Link gene back to taxa obj and update its network (taxa-->gene-->phenetype)
                txk_subj = genes_to_taxa[subj].ID
                txk_obj = genes_to_taxa[obj].ID

                if subj not in taxa_objs[txk_subj].gene_to_ortholog_family:
                    taxa_objs[txk_subj].gene_to_ortholog_family.update({subj:ortho_id})
                
                if obj not in taxa_objs[txk_obj].gene_to_ortholog_family:
                    taxa_objs[txk_obj].gene_to_ortholog_family.update({obj:ortho_id})
    
    print("- Non panther orthology relationships found {} (and removed)".format(
                                                                         format(non_panther_ortho_relationships, ',')))
    print("- Panther ortho relationships created {} ({} if you divide by 2)".format(format(len(ortho_data), ','),
                                                                            format(len(ortho_data)/2, ',')))
    
    # Add ortholog mapping counts to our dataframe (on a per taxa basis)
    txa_counts = [taxa_orth_counts[t_label] for t_label in list(taxa_table["Taxon Label"])]
    taxa_table["Ortholog Mapping Count"] = txa_counts
    taxa_table["Genes with >= 1 ortholog edge"] = [len(taxa_objs[t_id].gene_to_ortholog_family) for t_id in list(taxa_table["Taxon ID"])]
    taxa_table["Unique ortholog ids (orthogroups)"] = [len(set(taxa_objs[t_id].gene_to_ortholog_family.values())) for t_id in list(taxa_table["Taxon ID"])]

    ortho_data = np.asarray(ortho_data).T
    pd.DataFrame({"geneA":ortho_data[0],
                  "predicate":ortho_data[1],
                  "geneB":ortho_data[2],
                  "ortholog_id":ortho_data[3]}).to_csv(outfile_path, sep='\t', index=False)
    
    print("- Orthologs converted to array of shape {}".format(ortho_data.shape))
    print("- Orthologs file written to... {}".format(outfile_path))
    return taxa_objs, taxa_table


# Step 4)
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
            # This currently will write out all orthologs between the species
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


# Step 5)  
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
            
            # Only want genes that have an ortholog edges
            ortho_id = ortho_map.get(k, None)
            if not ortho_id:
                continue
            
            # Gene is not found in other species 
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


#################
### Version 2 ###

# Step 1)?
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
        print("- {} {}".format(k, format(v, ',')))
    
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
        print("- {} {}".format(k, format(v, ',')))
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


# Step 3
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
def generate_common_orthologs_files_v2(tables_dir, taxon_stats):
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
    
    # ##################
    # ### PROGRAM V1 ###

    # # Input files
    # nodes_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_nodes.tsv")
    # edges_file = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")

    # # Output files / info
    # species_data_dir = os.path.join(args.project_dir, "species_data")
    # ortho_map_file = os.path.join(species_data_dir, "panther_orthologs.tsv")
    # common_ortho_prefix = os.path.join(species_data_dir, "common_orthologs")
    # species_info_table_file = os.path.join(species_data_dir, "species_information_table.tsv")
    # species_info_table_file_unfiltered = os.path.join(species_data_dir, "species_information_table_unfiltered.tsv")
    # common_ortho_counts_file = os.path.join(species_data_dir, "species_common_orthologs_counts.tsv")


    # # Create initial taxa objects and "mapping tables" from monarch nodes file
    # taxa_objs, genes_to_taxa, genes_to_name, phen_to_name, taxa_table = pull_taxon_gene_information(nodes_file)

    # # Fill in taxa object gene_to_phenotype dictionary datastructure with relevant edges from monarch kg edge file
    # taxa_objs, taxa_table = fill_in_gene_to_phenotype_networks(edges_file, genes_to_taxa, taxa_objs, taxa_table)

    # # Write out orthologos mappings of genes across species. A-->B, AND B-->A
    # taxa_objs, taxa_table = generate_orthologs_mapping_file(edges_filepath=edges_file,
    #                                                         outfile_path=ortho_map_file,
    #                                                         genes_to_taxa=genes_to_taxa, 
    #                                                         phen_to_name=phen_to_name, 
    #                                                         taxa_objs=taxa_objs,
    #                                                         taxa_table=taxa_table)
    
    # # Convert to df and write unfiltered table (for record keeping)
    # taxa_table = pd.DataFrame(taxa_table)
    # taxa_table.to_csv(species_info_table_file_unfiltered, sep='\t', index=False)

    # # Filter taxa_table based on phenotype edges and orthologous gene mappings and write file (for downstream analysis)
    # taxa_table = taxa_table[(taxa_table["Total Phenotype Edges"].astype(int) > 0) & 
    #                         (taxa_table["Ortholog Mapping Count"].astype(int) > 0)]
    # taxa_table.to_csv(species_info_table_file, sep='\t', index=False)

    # # Filter out taxa that do not have the necessary information, so we don't generate excess / useless files
    # relevant_taxa_ids = set(taxa_table["Taxon ID"])
    # relevant_taxa_objs = {t:v for t,v in taxa_objs.items() if t in relevant_taxa_ids}

    # # Generate common orthologs files between all pairwise species comparisons
    # generate_common_orthologs_files(ortho_map_file, 
    #                                 common_ortho_prefix, 
    #                                 relevant_taxa_objs,
    #                                 ortho_counts_outpath=common_ortho_counts_file)

    # # Generate gene-->phenotype, and phenotype-->ortholog_ids files (tsv and pk respectivly)
    # generate_gene_and_phenotype_ortholog_files(ortho_map_file, 
    #                                            species_data_dir, 
    #                                            genes_to_name, 
    #                                            phen_to_name,
    #                                            relevant_taxa_objs)
    
    # print("- Species information table (used in downstream steps) written to {}...".format(species_info_table_file))
    # print("- Data extraction pipeline complete...")


    ##################
    ### PROGRAM V2 ###

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
    generate_common_orthologs_files_v2(tables_dir=species_data_dir, taxon_stats=taxon_stats)
    
    # Write taxon stats table to file
    pd.DataFrame(taxon_stats).to_csv(species_info_table_file, sep='\t', index=False)

    print("- Species information table (used in downstream steps) written to {}...".format(species_info_table_file))
    print("- Data extraction pipeline complete...")