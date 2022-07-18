########################################
# This script was used to perform GSVA on the big tcga exp dataset using different geneset lists to evaluate their information density and compare LUAD to other cancer types 
########################################


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
library(gplots)
library(gtools)
library(EnhancedVolcano)
library(BiocParallel)
library(colorRamps)

tcga_exp_short <- readRDS("./data/tcga_exp_small.RDS")
all_genesets_c5 <- readRDS("./data/genesetlist_whole_C5.RDS")

all_genesets <- append(geneset, metabolism_list_gs)

#################################### GSVA for all Tumor types
all_umap <- as.data.frame(umap(RunPCA(as.matrix(tcga_exp_short))@cell.embeddings, metric = "cosine"))

ggplot(all_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

####### get C5 pathway --> old metabolism genesets were replaced with C5 genesets
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

############################# delete dataset entries for genes not in pathways
tcga_exp_short <- tcga_exp[which(rownames(tcga_exp) %in% unlist(unique(all_genesets_copy))),]

######################### GSVA for all tumor types at once

all_gsva <- as.data.frame(gsva(as.matrix(tcga_exp_short), all_genesets_copy, method = "gsva")) #calculated on cluster

whole_gsva_umap <- as.data.frame(umap(t(all_gsva), metric = "cosine"))
ggplot(whole_gsva_umap, aes(x=V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()+
  

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

heatmap_try <-  pheatmap(as.matrix(pathway_enrichment_means_highSD[1:200,]), color = gplots::bluered(10), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 12, fontsize = 12 )
png(filename = "./output/enrichmentHeatmapTryTest.png", height = 4000, width = 2000, units = "px")
heatmap_try
dev.off()


heatmap_try2 <- Heatmap(as.matrix(pathway_enrichment_means_highSD[1:200,]), column_km = 3,column_names_gp =  grid::gpar(fontsize=8), row_names_gp = grid::gpar(fontsize = 4), height = nrow(pathway_enrichment_means_highSD)*unit(0.25, "mm"), width = ncol(pathway_enrichment_means_highSD)*unit(5, "mm"))
png(filename = "./output/enrichmentHeatmapTry2.png", height = 5200, width = 5500, units = "px", res = 500)
heatmap_try2
dev.off()


# extract LUADS most over and under-expressed pathways
LUAD_underexpressed <- rownames(pathway_enrichment_means_highSD[1:46,])[which(pathway_enrichment_means_highSD$LUAD[1:46] %in% sort(pathway_enrichment_means_highSD$LUAD[1:46])[1:20])]
LUAD_overexpressed <- rownames(pathway_enrichment_means_highSD[1:46,])[which(pathway_enrichment_means_highSD$LUAD[1:46] %in% sort(pathway_enrichment_means_highSD$LUAD[1:46], decreasing = TRUE)[1:20])]


# compare top 20 under and overexpressed genesets with other tumor types to find specific pathways
rest_means <- apply(as.data.frame(pathway_enrichment_means_C5)[,-4], 1, mean)

rest_underexpressed <- rownames(pathway_enrichment_means_highSD[1:46,])[which(rest_means[1:46] %in% sort(rest_means[1:46])[1:20])]
specific_underexpressed_LUAD <- LUAD_underexpressed[which(!(LUAD_underexpressed %in% rest_underexpressed))]

rest_overexpressed <- rownames(pathway_enrichment_means_highSD[1:46,])[which(rest_means[1:46] %in% sort(rest_means[1:46], decreasing = TRUE)[1:20])]
specific_overexpressed_LUAD <- LUAD_overexpressed[which(!(LUAD_overexpressed %in% rest_overexpressed))]

# plot results and compare plots
LUAD_expression_df <- data.frame("expressions" = c(sort(pathway_enrichment_means_highSD$LUAD[1:46])[1:20], sort(pathway_enrichment_means_highSD$LUAD[1:46], decreasing = TRUE)[1:20]))
LUAD_expression_df$genes <- c(LUAD_underexpressed, LUAD_overexpressed)

ggplot(LUAD_expression_df, aes( y = expressions, x = reorder(genes, expressions)))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45))

# calculate fold change between LUAD and all other cancer types and perform a vulcano plot
foldchange_LUADvsRest <- as.data.frame(pathway_enrichment_means_C5)$LUAD - rest_means
foldchange_LUADvsRest <- foldchange(as.data.frame(pathway_enrichment_means_C5)$LUAD, rest_means )
foldchange_LUADvsRest <- sapply(foldchange_LUADvsRest, function(x){log2(x)})

all_gsva_rest <- as.data.frame(all_gsva_C5_highSD[,which(tcga_anno$cancer_type_abbreviation != "LUAD")])

pvalues_LUADvsRest <- vector()
for(geneset in rownames(C5_devided$LUAD)){
  pvalues_LUADvsRest <- append(pvalues_LUADvsRest, wilcox.test(as.numeric(C5_devided$LUAD[geneset,]), as.numeric(all_gsva_rest[geneset,]), alternative = "two.sided")$p.value)
}

volcano_df <- data.frame("genset" = rownames(pathway_enrichment_means_C5), "pvalues" = pvalues_LUADvsRest, "foldchange" = foldchange_LUADvsRest)
volcano_df <- volcano_df[-17,]
EnhancedVolcano(volcano_df, lab = volcano_df$genset, x = "foldchange", y="pvalues", pointSize = 3.5, labSize = 4)

#################################### C5 geneset
all_gsva_C5 <- readRDS("./data/all_gsva_c5.RDS") #was calculated on a cluster
all_gsva_C5 <- as.data.frame(all_gsva_C5)

SDs <- apply(all_gsva_C5, 1, sd)

all_gsva_C5_highSD <- all_gsva_C5[which(SDs > 0.28),]

C5_PCAUMAP <- as.data.frame(umap(RunPCA(as.matrix(all_gsva_C5))@cell.embeddings, metric = "cosine"))

ggplot(C5_PCAUMAP, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()+
  guides(color=guide_legend(title="cancer type"))+
  theme_classic()+
  theme(legend.text=element_text(size=14))

C5_devided <- list()
for(types in tumor_types){
  C5_devided[[types]] <- as.data.frame(all_gsva_C5_highSD[,which(tcga_anno$cancer_type_abbreviation == types)])
}

C5_gsva_means <- lapply(C5_devided, function(x){apply(x, 1, mean)})

pathway_enrichment_means_C5 <- matrix(nrow = 90, ncol = 33, byrow = FALSE)

for( type in 1:length(tumor_types)){
  pathway_enrichment_means_C5[,type] <- as.vector(C5_gsva_means[[type]])
}


colnames(pathway_enrichment_means_C5) <- tumor_types
rownames(pathway_enrichment_means_C5) <- rownames(all_gsva_C5_highSD)

Heatmap(pathway_enrichment_means_C5, column_km = 3, row_names_gp =  grid::gpar(fontsize=7), row_names_side = "left", show_row_dend = FALSE, heatmap_legend_param = list(title = "relative expression"), row_names_max_width = unit(14, "cm"), height = unit(19, "cm"))


#################################### create second geneset list with C2 genesets from msigdb
C2_pathways <- msigdbr(species = "Homo sapiens", category = "C2")

C2_list = list() #this creates an empty list to store the for loop in
for (gene.name in unique(C2_pathways$gs_name)){
  pathway = c(C2_pathways[C2_pathways$gs_name == gene.name, 7])
  names(pathway) = paste(c(strsplit(gene.name, "_", fixed = TRUE)[[1]][-1]), collapse = "_") #to extract the "GOBP_" from each name, substring is used
  C2_list[paste(c(strsplit(gene.name, "_", fixed = TRUE)[[1]][-1]), collapse = "_")] = pathway}

#append C2 sets to rest of the genesets
combined_genesets_C2 <- append(all_genesets, C2_list)


#clean up c2 genesets
combined_genesets_C2 <- combined_genesets_C2[[which(sum(duplicated(combined_genesets_C2))<1)]]

pathway_intersections_C2 <- lapply(combined_genesets_C2, function(x){gene_intersection(pathway = x)})

pathway_intersections_C2 <- pathway_intersections_C2[which(as.numeric(pathway_intersections_C2)>0.97)]

combined_genesets_C2_short <- combined_genesets_C2[which(names(combined_genesets_C2) %in% names(pathway_intersections_C2))]

saveRDS(combined_genesets_C2_short, "./data/pathwaysC2.RDS")

#heatmap for c2 genesets
all_gsva_c2 <- readRDS("./data/all_gsva_c2.RDS")
all_gsva_c2 <- as.data.frame(all_gsva_c2)

SDs <- apply(all_gsva_c2, 1, sd)
all_gsva_C2_highSD <- all_gsva_c2[which(SDs > 0.295),]

C2_PCAUMAP <- as.data.frame(umap(RunPCA(as.matrix(all_gsva_c2))@cell.embeddings, metric = "cosine"))

ggplot(C2_PCAUMAP, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()+
  guides(color=guide_legend(title="cancer type"))+
  theme_classic()+
  theme(legend.text=element_text(size=14))

C2_devided <- list()
for(types in tumor_types){
  C2_devided[[types]] <- as.data.frame(all_gsva_C2_highSD[,which(tcga_anno$cancer_type_abbreviation == types)])
}

C2_gsva_means <- lapply(C2_devided, function(x){apply(x, 1, mean)})

pathway_enrichment_means_C2 <- matrix(nrow = 86, ncol = 33, byrow = FALSE)

for( type in 1:length(tumor_types)){
  pathway_enrichment_means_C2[,type] <- as.vector(C2_gsva_means[[type]])
}


colnames(pathway_enrichment_means_C2) <- tumor_types
rownames(pathway_enrichment_means_C2) <- rownames(all_gsva_C2_highSD)

Heatmap(pathway_enrichment_means_C2, column_km = 3, row_names_gp =  grid::gpar(fontsize=7), row_names_side = "left", show_row_dend = FALSE, heatmap_legend_param = list(title = "relative expression"), row_names_max_width = unit(14, "cm"))

#################################### Intersections for Troubleshooting, was discarded later

#gene_intersection <- function(names.dataset = rownames(tcga_exp_short), pathway){
#  return(c(sum(pathway %in% names.dataset)/(length(pathway))))
#}

#pathway_intersections <- lapply(all_genesets, function(x){gene_intersection(pathway = x)})


#all_genes_in_dataset<- sum(rownames(tcga_exp) %in% unlist(unique(all_genesets)))/length(rownames(tcga_exp))
#all_genes_in_dataset




