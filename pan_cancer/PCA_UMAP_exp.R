library(parallel)
library(FactoMineR)
library(Seurat)
library(ggplot2)
library(uwot)
library(ggpubr)
library(babyplots)

exp_highvar <- t(readRDS("./data/tcga_exp_small.RDS"))
tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")


#tcga_exp_pca <- PCA(exp_highvar, ncp = 2)

tcga_exp_pca <- RunPCA(exp_highvar)
tcga_exp_pca <- as.data.frame(tcga_exp_pca@cell.embeddings)

ggplot(tcga_exp_pca, aes(x = PC_1, y = PC_2, color = tcga_anno$cancer_type_abbreviation)) +
  geom_point()


tcga_exp_umap <- as.data.frame(umap(tcga_exp_pca, n_threads = 8, metric = "cosine", n_components = 3))
ggplot(tcga_exp_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

pointCloud(as.matrix(tcga_exp_umap), colorBy = "categories", colorVar = tcga_anno$cancer_type_abbreviation, turntable = TRUE, rotationRate = 0.001, xScale = 0.25, yScale = 0.25, zScale = 0.25 )

###################################
# save different tumor types in one list and run separate PCAs

tumor_type_dfs <- list()
tumor_types <- unique(tcga_anno$cancer_type_abbreviation)


anno_tumortypes <- list()
for(types in tumor_types){
  anno_tumortypes[[types]] <- as.data.frame(tcga_anno[which(tcga_anno$cancer_type_abbreviation == types),])
}
anno_tumortypes <- anno_tumortypes[names(tumor_type_umap)]

for(types in tumor_types){
  tumor_type_dfs[[types]] <- as.data.frame(exp_highvar[,which(tcga_anno$cancer_type_abbreviation == types)])
}


tumor_type_pca <- list()
tumor_type_pca <- lapply(tumor_type_dfs, FUN=RunPCA)     

RunPCA(exp_highvar[,which(tcga_anno$cancer_type_abbreviation=="LUAD")])

tumor_type_umap <- list()
for(types in tumor_types){
  tumor_type_umap[[types]] <- as.data.frame(umap(RunPCA(exp_highvar[,which(tcga_anno$cancer_type_abbreviation== types)])@cell.embeddings))
}

tumor_type_umap <- list()
for(types in tumor_types){
  tumor_type_umap[[types]] <- as.data.frame(umap(exp_highvar[,which(tcga_anno$cancer_type_abbreviation== types)]))
}


ggplots <- list()
for(plot in 1:length(tumor_type_umap)){
  print(ggplot(tumor_type_umap[[plot]], aes(x = V1, y=V2, color=anno_tumortypes[[plot]]$gender))+
                      geom_point()+
                      ggtitle(names(tumor_type_umap)[plot]))
}

ggarrange(plotlist = ggplots, ncol=5, nrow = 6)



ggplot(tumor_type_umap[[1]], aes(x = V1, y=V2, color=anno_tumortypes[[1]]$gender))+
  geom_point()+
  ggtitle(names(tumor_type_umap)[1])
