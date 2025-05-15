# MIRit Supporting Files <img src="logo.png?raw=true" alt="MIRit_logo" height="139" align="right"/>

This repository contains the analysis files used in the manuscript to demonstrate the effectiveness of MIRit in integrative miRNA-mRNA analyses. For more information on how MIRit works, please refer to the original paper. The source code of MIRit is accessible on [GitHub](https://github.com/jacopo-ronchi/MIRit).

## Authors

__Dr. Jacopo Ronchi__ <a itemprop="sameAs" content="https://orcid.org/0000-0001-5520-4631" href="https://orcid.org/0000-0001-5520-4631" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup> (author and maintainer)

__Dr. Maria Foti__ <a itemprop="sameAs" content="https://orcid.org/0000-0002-4481-1900" href="https://orcid.org/0000-0002-4481-1900" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a><sup>1</sup>

<sup>1</sup>School of Medicine and Surgery, University of Milano-Bicocca, Italy

## Reference

The publication that details the implementation of MIRit and its usage is currently available on [bioRxiv](doi:10.1101/2023.11.24.568528):

> Ronchi J and Foti M. ‘MIRit: an integrative R framework for the
> identification of impaired miRNA-mRNA regulatory networks in complex
> diseases’. bioRxiv (2023). <doi:10.1101/2023.11.24.568528>

## Datasets description

The performance of MIRit has been tested with three different datasets:

1. A sample-matched experiment in which miRNA and mRNA expression were evaluated in patients with dilated cardiomyopathy using miRNA-Seq and RNA-Seq, respectively.
2. A sample-matched experiment in which miRNA and mRNA expression were profiled in clear cell renal cell carcinoma (ccRCC) samples and normal adjacent renal tissue.
3. A dataset without matched samples, in which miRNA expression was evaluated in postmortem BA9 tissue from a cohort of patients with Alzheimer's disease (AD) and gene expression was assessed using microarray technology in the same brain region, but in a different cohort of patients.

The datasets used in this study are publicly available on GEO under the accession numbers [GSE243406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243406), [GSE16441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi), [GSE63501](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi), and [GSE150696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi).

## File organization

The R scripts used to generate the results presented in the manuscript can be found in the "Dilated Cardiomyopathy," "Clear Cell Renal Cell Carcinoma," and "Alzheimer's Disease" folders. R scripts used for data preprocessing are also included when necessary.
