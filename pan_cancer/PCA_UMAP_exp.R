library(parallel)
library(FactoMineR)
library(Seurat)
library(ggplot2)
library(uwot)
exp_highvar <- t(readRDS("./data/tcga_exp_small.RDS"))


#tcga_exp_pca <- PCA(exp_highvar, ncp = 2)

tcga_exp_pca <- RunPCA(exp_highvar)
tcga_exp_pca <- as.data.frame(tcga_exp_pca@cell.embeddings)

ggplot(tcga_exp_pca, aes(x = PC_1, y = PC_2, color = tcga_anno$cancer_type_abbreviation)) +
  geom_point()


tcga_exp_umap <- as.data.frame(umap(tcga_exp_pca, n_threads = 8, metric = "euclidean"))
ggplot(tcga_exp_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

