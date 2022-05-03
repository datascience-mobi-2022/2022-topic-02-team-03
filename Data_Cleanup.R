# Please set working directory to folder containing the following folders with the datasets
tcga_exp <- readRDS("./RNA-SeqAllTumors/tcga_tumor_log2TPM.RDS")
tcga_anno <- readRDS("./RNA-SeqAllTumors/tcga_tumor_annotation.RDS")
tumor_vs_norm <- readRDS("./tumorvsnormal/tcga_tumor_normal_datascience_proj_2022.RDS")
genesets <- readRDS("./hallmarkGenes/hallmarks_genesets.rds")

tcga_exp_copy <- tcga_exp

IDs <- c()
for (element in rownames(tcga_exp)) {
  IDs <- c(IDs, strsplit(element, "|", fixed = TRUE)[[1]][2])
}

ensemble <- c()
for (element in rownames(tcga_exp)) {
  ensemble <- c(ensemble, strsplit(element, "|", fixed = TRUE)[[1]][1])
}

rownames(tcga_exp_copy) <- make.names(IDs, unique = TRUE)
