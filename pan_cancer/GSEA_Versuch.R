library(parallel)
library(FactoMineR)
library(Seurat)
library(ggplot2)
library(uwot)
library(ggpubr)
library(babyplots)
library(cluster)
library(fgsea)
library(GSVA)
library(pheatmap)
library(RColorBrewer)


######################## enrichment test using ssGSVA
gsva_list <- list()
gsva_list <- mclapply(tumor_type_dfs, function(x){gsva(as.matrix(x), geneset, method = "zscore")}, mc.cores = 6)

# quality control

qc_umap_gsva <- list()
for(types in tumor_types){
  qc_umap_gsva[[types]] <- as.data.frame(umap(RunPCA(as.matrix(gsva_list[[types]]), npcs = 25)@cell.embeddings, metric = "cosine"))
}


ggplot(as.data.frame(qc_umap_gsva$LUAD), aes(x=V1, y=V2))+
  geom_point()


###############################


percent_metabol <- sapply(metabolism_list_gs, function(genesets){len.geneset = length(genesets); intersection = length(sum(which(genesets %in% rownames(exp_highvar)))); return(intersection/len.geneset)})
percent_hallmark <- sapply(geneset, function(genesets){len.geneset = length(genesets); intersection = length(sum(which(genesets %in% rownames(exp_highvar)))); return(intersection/len.geneset)})


write.csv(c(percent_hallmark, percent_metabol),file="output.csv",row.names=TRUE)



#################################

BRCA_and_LUAD <- saveRDS(list("BRCA" = tumor_type_dfs[["BRCA"]], "LUAD" =tumor_type_dfs[["LUAD"]]), "LUADandBRCA.RDS")


test <- readRDS("LUADandBRCA.RDS")

all_genesets <- append(geneset, metabolism_list_gs)

saveRDS(all_genesets, "allGenesets.RDS")
tes <- readRDS("allGenesets.RDS")


######################## Test with all genes / uncleaned dataset
big_tumor_dfs <- list()
tumor_types <- unique(tcga_anno$cancer_type_abbreviation)

tcga_exp_backup <- tcga_exp
#remove accession numbers
IDs <- c()
for (element in rownames(tcga_exp)) {
  IDs <- c(IDs, strsplit(element, "|", fixed = TRUE)[[1]][2])
}
rownames(tcga_exp) <- make.names(IDs, unique = TRUE)


for(types in tumor_types){
  big_tumor_dfs[[types]] <- as.data.frame(tcga_exp[,which(tcga_anno$cancer_type_abbreviation == types)])
}

big_umap <- list()
for(types in tumor_types){
  big_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(big_tumor_dfs[[types]]), npcs = 25)@cell.embeddings, metric = "cosine"))
}

big_ggplots <- list()
for(plot in 1:length(big_umap)){
  big_ggplots[[plot]] <- ggplot(big_umap[[plot]], aes(x = V1, y=V2))+
          geom_point(size = 0.8)+
          ggtitle(names(big_umap)[plot])
}

ggarrange(plotlist = big_ggplots)



#gsva and umap


big_gsva_list <- list()
big_gsva_list <- mclapply(big_tumor_dfs, function(x){gsva(as.matrix(x), all_genesets, method = "zscore")}, mc.cores = 6)

big_gsva_umap <- list()
for(types in tumor_types){
  big_gsva_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(big_gsva_list[[types]]))@cell.embeddings, metric = "cosine"))
}

big_ggplots_gsva <- list()
for(plot in names(big_gsva_umap)){
  big_ggplots_gsva[[plot]] <- ggplot(big_gsva_umap[[plot]], aes(x = V1, y=V2))+
    geom_point(size = 0.8)+
    ggtitle(plot)
}

ggarrange(plotlist = big_ggplots_gsva)


################################################## pathway activity heatmap



pathway_enrichment_means <- matrix(nrow = 87, ncol = 33, byrow = FALSE)
for( type in 1:length(tumor_types)){
  pathway_enrichment_means[,type] <- apply(as.data.frame(big_gsva_list[[type]]), 1, mean)
}
pathway_enrichment_means <- as.data.frame(pathway_enrichment_means)
colnames(pathway_enrichment_means) <- tumor_types
rownames(pathway_enrichment_means) <- rownames(big_gsva_list[["LUAD"]])

heatmap <- pheatmap(as.matrix(pathway_enrichment_means), color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 18, fontsize = 18 )
png(filename = "./output/enrichmentHeatmap.png", height = 1400, width = 2000, units = "px")
heatmap
dev.off()
