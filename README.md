# monarch-phenologs



This repo will contain code to calculate phenologs and gene candidate predictions.



Phenologs 
Phenologs extends the concept of orthologous genes to orthologous phenotypes. By using a  k-nearest neighbors approach to compare sets of gene-phenotype relationships across species using gene orthology connections.

Inclusion of model organisms into the calculations to generate phenologs is dependent on the availability of gene-phenotype associations as well as gene orthology annotations. For ease of data inclusion, we will attempt to utilize data already ingested and formatted within the Monarch KG.

Source datasets:
- Mondo
- Human Phenotype Ontology
- Mouse Genome Database
- Zebrafish Information Network

Rat
Mouse
Zebrafish
Worm
Chicken
Fission yeast?



Methods

Identification of phenologs through gene enrichment analysis


Gene candidate predictions by k-nearest neighbors
Calculation of false discovery rate
Distance matrix
Weight matrix






References

[McGary KL, Park TJ, Woods JO, Cha HJ, Wallingford JB, Marcotte EM. Systematic discovery of nonobvious human disease models through orthologous phenotypes. Proc Natl Acad Sci U S A. 2010 Apr 6;107(14):6544-9. doi: 10.1073/pnas.0910200107. Epub 2010 Mar 22. PMID: 20308572; PMCID: PMC2851946.](https://www.pnas.org/doi/10.1073/pnas.0910200107)

[Woods JO, Singh-Blom UM, Laurent JM, McGary KL, Marcotte EM. Prediction of gene-phenotype associations in humans, mice, and plants using phenologs. BMC Bioinformatics. 2013 Jun 21;14:203. doi: 10.1186/1471-2105-14-203. PMID: 23800157; PMCID: PMC3704650.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-203)



