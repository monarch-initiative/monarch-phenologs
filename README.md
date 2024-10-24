# monarch-phenologs


Phenologs extends the concept of orthologous genes to orthologous phenotypes. 

By comparing the genes associated with phenotypes across species, utilizing  looking for enrichment
Additionally, gene candidate predictions can be made using a  k-nearest neighbors approach 
to compare sets of gene-phenotype relationships across species using gene orthology connections.

Inclusion of model organisms into the calculations to generate phenologs is dependent on the availability of 
gene-phenotype associations as well as gene orthology annotations. For ease of data inclusion, 
we will attempt to utilize data already ingested and formatted within the Monarch KG.

Included organisms:
- Humans
- Mouse
- Rats
- Zebrafish
- C. Elegans
- Xenopus

Source datasets (parsed from the Monarch KG):
- Mondo
- Human Phenotype Ontology
- Mouse Genome Database
- Rat Genome Database
- Zebrafish Information Network
- worm?
- Xenopus?





## Methods

### Data Sources:

This implementation of phenologs uses the Monarch Knowledge Graph as the data source, 
with the intention of integrating the calculated phenologs and phenolog gene candidate predictions 
back to that same instance/version of the Monarch KG after calculations are complete.



### Calculation of phenologs

The calculation of phenologs begins by acquiring the latest version of the Monarch KG and extracting data from the Monarch KG.

Data is parsed into phenotype-to-gene tables, and these are then converted to phenotype to ortholog files, 
collapsing genes with the same mapped ortholog.

(Still to implement: Collapsing phenotypes with 80-90% gene overlaps into a single phenotype for 
the phenolog calculation, to avoid redundant comparisons. See supplemental materials of the McGary et. al. paper for details.)

Before we can begin calculating the phenologs from the phenotype-ortholog files, we must first calculate a False Discovery Rate (FDR).
This requires generation of 1000 randomized datasets, calculating the p-values for every phenolog comparison in those randomized datasets,
and determining the p-value at which 5% of these randomized phenolog calculations are significant.

Once the FDR has been calculated, the final phenolog calculations can be performed, 
and the significant phenologs can be identified simply as those with a p-value less than the calculated FDR.

### Calculation of gene candidate predictions by k-nearest neighbors

Calculation of false discovery rate

#### Distance matrix
Distance matrix is calculated by...

#### Weight matrix
The weight matrix is calculated by...

#### k-nearest neighbors





References

[McGary KL, Park TJ, Woods JO, Cha HJ, Wallingford JB, Marcotte EM. Systematic discovery of nonobvious human disease models through orthologous phenotypes. Proc Natl Acad Sci U S A. 2010 Apr 6;107(14):6544-9. doi: 10.1073/pnas.0910200107. Epub 2010 Mar 22. PMID: 20308572; PMCID: PMC2851946.](https://www.pnas.org/doi/10.1073/pnas.0910200107)

[Woods JO, Singh-Blom UM, Laurent JM, McGary KL, Marcotte EM. Prediction of gene-phenotype associations in humans, mice, and plants using phenologs. BMC Bioinformatics. 2013 Jun 21;14:203. doi: 10.1186/1471-2105-14-203. PMID: 23800157; PMCID: PMC3704650.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-203)



