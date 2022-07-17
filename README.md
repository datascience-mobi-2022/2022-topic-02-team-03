# 2022-topic-02-team-03: lung adenocarcinoma (LUAD)

**Team:** Paul Brunner, Marie Kleinert, Felipe Stünkel, Chloé Weiler
<br/> **Supervisor:** PD. Dr. Carl Herrmann
<br/> **Tutor:** Wangjun Hu

# Abstract 
## Biological Background and Objectives
Lung adenocarcinoma (LUAD) is the most common type of non-small cell lung cancer. It distinguishes itself through its particularly low 5-year overall survival rate of merely 18%. Knowing those facts it becomes evident as to why LUAD is an intriguing research topic and so over the course of our project, we set out to analyse gene expression profiles of LUAD patients before and after developing the disease as well as compare LUAD to 32 other cancer types. 
<br/>By doing so we managed to gain a better understanding over which specific gain-of-function and loss-of-function gene mutations drive the development of LUAD through gene expression deregulation and we were able to create a logistic regression model that can predict whether a person will eventually develop LUAD based on their current gene expression profiles.
## Our Data
Apart from the analysis of cancer hallmark pathways which was done on a dataset conatining a list of gene sets, all our analyses were performed on one of two datasets: **tcga_exp_log2TPM** and **tcga_tumor_normal**.
- **tcga_exp_log2TPM** contains RNA-seq data from almost 10,000 TCGA cancer patients for 33 different tumor types, retrieved from The Cancer Genome Atlas (TCGA) and normalized through a TPM ('transcripts per million') normalization and a  log<sub>2</sub> transformation.
- **tcga_tumor_normal** contains TCGA expression data of tumor tissue and of the corresponding healthy tissue for five different cancer types.

# Folder structure

- **[data cleaning](/data_cleaning):** pre-cleaning of tcga_exp and tcga_tumor_normal data

- **[cancer hallmark pathways](/cancer_hallmark_pathways):** comparison of similarity between provided cancer hallmarks gene sets and the metabolic gene sets selected from MSigDB database

- **[pan cancer](/pan_cancer):** All the code used to perform a pan cancer analysis on tcga_exp cleaned data

- **[normal vs tumor](/normal_vs_tumor):** All the code used to perform a focused analysis on tcga_tumor_normal cleaned data

- **[output](/output):** all of our plots

- **[report](/report):** our final report and .bib files of our citations

- **data:** stores raw data locally
