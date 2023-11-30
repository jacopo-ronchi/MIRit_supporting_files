# MIRit Supporting Files <img src="logo.png?raw=true" alt="MIRit_logo" height="139" align="right"/>

This repository contains the analysis files that were used in the paper to showcase the performance of MIRit for integrative miRNA-mRNA analyses. For details regarding the functioning of MIRit, please refer to the original publication. The source code of the R package is available on GitHub in the main repository of [MIRit](https://github.com/jacopo-ronchi/MIRit).

## Authors

__Dr. Jacopo Ronchi__ <a itemprop="sameAs" content="https://orcid.org/0000-0001-5520-4631" href="https://orcid.org/0000-0001-5520-4631" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup> (author and maintainer)

__Dr. Maria Foti__ <a itemprop="sameAs" content="https://orcid.org/0000-0002-4481-1900" href="https://orcid.org/0000-0002-4481-1900" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup>

<sup>1</sup>School of Medicine and Surgery, University of Milano-Bicocca, Italy

## Reference

The corresponding publication where we describe the implementation of MIRit and its usage is available on [bioRxiv](doi:10.1101/2023.11.24.568528):

> Ronchi J and Foti M. ‘MIRit: an integrative R framework for the
> identification of impaired miRNA-mRNA regulatory networks in complex
> diseases’. bioRxiv (2023). <doi:10.1101/2023.11.24.568528>

## Datasets description

MIRit's behavior has been tested with two different datasets:

- sequencing data obtained from Riesco-Eizaguirre et al. (2015) who profiled gene and miRNA expression in 8 papillary thyroid carcinoma tumors and corresponding contralateral tissue from the same patients;
- a small RNA-Seq experiment that examines post-mortem BA9 tissues from 6 patients with AD and 7 healthy controls, in addition to transcriptomic data where the BA9 region has been profiled in 9 patients with AD and 9 healthy donors through microarray technology.

All the datasets are obtained from Gene Expression Omnibus (GEO), with accession numbers: [GSE63511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63511), [GSE63501](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63501), and [GSE150696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150696). Before the analyses, data preprocessing was carried out as described in the manuscript. Only preprocessed expression data are included in this repository, whereas raw data are available on GEO.

## File organization

The analysis files for each example are stored in the relative folders: "Thyroid cancer" and "Alzheimer's disease". Specifically, each example contains:

- The R script that was used to generate the results, `thyroid_cancer_analysis.R` for the first example and `alzheimer_analysis.R` for the second one,
- A `data` folder with processed expression matrices and sample metadata,
- A `results` folder containing the raw results of the analysis;
- A `plots` folder that stores the generated figures used in the manuscript.
