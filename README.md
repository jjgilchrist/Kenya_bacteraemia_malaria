# Kenya_bacteraemia_malaria
**BIRC6 modifies risk of invasive bacterial infection in Kenyan children.** (eLife, 2022)

https://elifesciences.org/articles/77461

Code and source data for GWAS of invasive bacterial disease in Kenyan children.

**Dataset**
5,400 Kenyan children: 1,445 bacteraemia cases, 1,143 severe malaria, 2,812 controls.

**Genome-wide association analysis analysis**
In this analysis we utilised probabilistic diagnostic models to identify children with a high probability of invasive bacterial disease among critically unwell Kenyan children with P. falciparum parasitaemia. We construct a joint dataset including 1,445 bacteraemia cases and 1,143 severe malaria cases, and population controls, among critically unwell Kenyan children that have previously been genotyped for human genetic variation. Using these data we perform a cross-trait genome-wide association study of invasive bacterial infection, weighting cases according to their probability of bacterial disease.

**Overview of repository**
* Figure folders: contains R script and source data to reproduce each main and supplementary Figure from the manuscript.
* Weighted_GWAS: contains R script and example data to run weighted logistic regression association analysis across a 1000 SNP chunk on chromosome 1.

**Data availability**
Patient level genotype and phenotype data are available are available via the European Genome-Phenome Archive, with accession codes EGAD00010000950 (WTCCC2: bacteraemia cases and controls) and EGAD00010000904 (MalariaGEN Consortium: severe malaria cases and controls). Full GWAS summary statistics will be made available on publication via the GWAS Catalog with accession code GCST90094632.
