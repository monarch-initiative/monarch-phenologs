'''
Using the randomized datasets, we need to calculate a false discovery rate that will be used to identify the
"significant" phenologs. Because this involves processing 1000s of randomized datasets, will need to
utilize python's multiprocessing capabilities or seek out other methods for improving performance.

Caveats:
Need to check and see if there's any issues with p-values being so small that they are effectively 0.


"Significant phenologs were identified at a FDR of 0.05 by ranking real and permuted phenologs on the basis of
the associated hypergeometric probabilities and selecting a threshold of probability where the proportion of
permuted phenologs above the cutoff accounted for 5% of the phenologs."


Steps:
-For each species pair, iterate through the randomized datasets.
-Calculate p-values, and assemble into list? dataframe?.
-Once completed with all phenolog calculations, will need to assemble
-Identify which p-value represents the 5% significance cutoff.


NOTE: Phenolog calculation script: Could that be generalizable enough to create a utility script to be called from multiple processing scripts?

'''


panther_filepath = "../datasets/intermediate/panther/panther_orthologs.tsv"

human_gene_prefix = 'HGNC:'
human_gene_to_phenotype_filepath = "../datasets/intermediate/human/human_gene_to_phenotype.tsv"
human_random_hvm_filepath = "../datasets/intermediate/random/human/human_hvm_"
human_random_hvr_filepath = "../datasets/intermediate/random/human/human_hvr_"
human_random_hvw_filepath = "../datasets/intermediate/random/human/human_hvw_"
human_random_hvz_filepath = "../datasets/intermediate/random/human/human_hvz_"

mouse_gene_prefix = 'MGI:'
mouse_gene_to_phenotype_filepath = "../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
mouse_random_mvh_filepath = "../datasets/intermediate/random/mouse/mouse_mvh_"
mouse_random_mvr_filepath = "../datasets/intermediate/random/mouse/mouse_mvr_"
mouse_random_mvw_filepath = "../datasets/intermediate/random/mouse/mouse_mvw_"
mouse_random_mvz_filepath = "../datasets/intermediate/random/mouse/mouse_mvz_"

rat_gene_prefix = 'RGD:'
rat_gene_to_phenotype_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
rat_random_rvh_filepath = "../datasets/intermediate/random/rat/rat_rvh_"
rat_random_rvm_filepath = "../datasets/intermediate/random/rat/rat_rvm_"
rat_random_rvw_filepath = "../datasets/intermediate/random/rat/rat_rvw_"
rat_random_rvz_filepath = "../datasets/intermediate/random/rat/rat_rvz_"

worm_gene_prefix = 'WB:'
worm_gene_to_phenotype_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
worm_random_wvh_filepath = "../datasets/intermediate/random/worm/worm_wvh_"
worm_random_wvm_filepath = "../datasets/intermediate/random/worm/worm_wvm_"
worm_random_wvr_filepath = "../datasets/intermediate/random/worm/worm_wvr_"
worm_random_wvz_filepath = "../datasets/intermediate/random/worm/worm_wvz_"

zebrafish_gene_prefix = 'ZFIN:'
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
zebrafish_random_zvh_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvh_"
zebrafish_random_zvm_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvm_"
zebrafish_random_zvr_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvr_"
zebrafish_random_zvw_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvw_"

organism_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']

for i in range(0, 10):








