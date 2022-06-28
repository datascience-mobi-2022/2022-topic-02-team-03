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
library(msigdbr)
library(ComplexHeatmap)
library(grid)



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


#################################### GSVA for all Tumor types
all_umap <- as.data.frame(umap(RunPCA(as.matrix(tcga_exp))@cell.embeddings, metric = "cosine"))

ggplot(all_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

all_gsva <- gsva(as.matrix(tcga_exp),all_genesets, method = "zscore", min.sz = 5)



#################################### Intersections for Troubleshooting

gene_intersection <- function(names.dataset = rownames(tcga_exp), pathway){
  return(c(sum(pathway %in% names.dataset)/(length(pathway))))
}

pathway_intersections <- lapply(all_genesets, function(x){gene_intersection(pathway = x)})


all_genes_in_dataset<- sum(rownames(tcga_exp) %in% unlist(unique(all_genesets)))/length(rownames(tcga_exp))
all_genes_in_dataset



################################# get more pathways for more coverage
GO_pathways <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

GO_list = list() #this creates an empty list to store the for loop in
for (gene.name in unique(GO_pathways$gs_name)){
  pathway = c(GO_pathways[GO_pathways$gs_name == gene.name, 7])
  names(pathway) = substring(gene.name, 6) #to extract the "GOBP_" from each name, substring is used
  GO_list[substring(gene.name, 6)] = pathway}

pathway_intersections_GO <- lapply(GO_list, function(x){gene_intersection(pathway = x)})
intersections_GO <- unlist(pathway_intersections_GO)
number_of_genes <- unlist(lapply(GO_list, function(x){return(length(x))}))
GO_intersect_overview <- data.frame("intersect" = intersections_GO, "number_of_genes" = number_of_genes)


boxplot(intersections_GO[which(number_of_genes > 20)])
GO_intersect_overview <- GO_intersect_overview[which(number_of_genes > 20),]

# 25% quantile = 0.95
GO_intersect_overview <- GO_intersect_overview[which(GO_intersect_overview$intersect > 0.95),]
boxplot(GO_intersect_overview$intersect)


# combine pathway
all_genesets_copy <- append(all_genesets, GO_list[which(names(GO_list) %in% rownames(GO_intersect_overview))])


############################# delete dataset entris for genes not in pathways
tcga_exp_short <- tcga_exp[which(rownames(tcga_exp) %in% unlist(unique(all_genesets_copy))),]


############################## new GSVA with short dataset
short_types_dfs <- list()
for(types in tumor_types){
  short_types_dfs[[types]] <- as.data.frame(tcga_exp_short[,which(tcga_anno$cancer_type_abbreviation == types)])
}

short_gsva_list <- list()
short_gsva_list <- mclapply(short_types_dfs, function(x){gsva(as.matrix(x), all_genesets_copy, method = "zscore")}, mc.cores = 6)


pathway_enrichment_means <- matrix(nrow = 2941, ncol = 33, byrow = FALSE)
for( type in 1:length(tumor_types)){
  pathway_enrichment_means[,type] <- apply(as.data.frame(short_gsva_list[[type]]), 1, mean)
}
pathway_enrichment_means <- as.data.frame(pathway_enrichment_means)
colnames(pathway_enrichment_means) <- tumor_types
rownames(pathway_enrichment_means) <- rownames(short_gsva_list[["LUAD"]])

heatmap <- pheatmap(as.matrix(pathway_enrichment_means), color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 18, fontsize = 18 )
png(filename = "./output/enrichmentHeatmap.png", height = 6000, width = 5000, units = "px")
heatmap
dev.off()

######################### GSVA for all tumor types at once

all_gsva <- gsva(as.matrix(tcga_exp_short), all_genesets_copy, method = "zscore")

patients_heatmap <- pheatmap(all_gsva[, 1:200], color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 18, fontsize = 18 )
png(filename = "./output/patientsenrichmentHeatmap.png", height = 30000, width = 5000, units = "px")
patients_heatmap
dev.off()

# subset the gsva df into the different tumor types
all_gsva_devided <- list()
for(types in tumor_types){
  all_gsva_devided[[types]] <- as.data.frame(all_gsva[,which(tcga_anno$cancer_type_abbreviation == types)])
}

# calculate means for every pathway
all_gsva_means <- lapply(all_gsva_devided, function(x){apply(x, 1, mean)})

# plot heatmap
pathway_enrichment_means <- matrix(nrow = 2941, ncol = 33, byrow = FALSE)
for( type in 1:length(tumor_types)){
  pathway_enrichment_means[,type] <- as.vector(all_gsva_means[[type]])
}
pathway_enrichment_means <- as.data.frame(pathway_enrichment_means)
colnames(pathway_enrichment_means) <- tumor_types
rownames(pathway_enrichment_means) <- rownames(short_gsva_list[["LUAD"]])


#kick out pathways that have a small sd
SDs <- apply(pathway_enrichment_means, 1, sd)
pathway_enrichment_means_highSD <- pathway_enrichment_means[which(SDs > 3),]

heatmap_try <-  pheatmap(as.matrix(pathway_enrichment_means_highSD[1:200,]), color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 12, fontsize = 12 )
png(filename = "./output/enrichmentHeatmapTry.png", height = 4000, width = 2000, units = "px")
heatmap_try
dev.off()


heatmap_try2 <- Heatmap(as.matrix(pathway_enrichment_means_highSD[1:200,]), column_km = 3,column_names_gp =  grid::gpar(fontsize=8), row_names_gp = grid::gpar(fontsize = 4), height = nrow(pathway_enrichment_means_highSD)*unit(0.25, "mm"), width = ncol(pathway_enrichment_means_highSD)*unit(5, "mm"))
png(filename = "./output/enrichmentHeatmapTry2.png", height = 5200, width = 5500, units = "px", res = 500)
heatmap_try2
dev.off()
