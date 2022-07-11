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

tcga_exp_short <- readRDS("./data/tcga_exp_small.RDS")
all_genesets_c5 <- readRDS("./data/genesetlist_whole_C5.RDS")
######################## enrichment test using ssGSVA
# gsva_list <- list()
# gsva_list <- mclapply(tumor_type_dfs, function(x){gsva(as.matrix(x), geneset, method = "zscore")}, mc.cores = 6)

# quality control

# qc_umap_gsva <- list()
# for(types in tumor_types){
#   qc_umap_gsva[[types]] <- as.data.frame(umap(RunPCA(as.matrix(gsva_list[[types]]), npcs = 25)@cell.embeddings, metric = "cosine"))
# }
# 
# 
# ggplot(as.data.frame(qc_umap_gsva$LUAD), aes(x=V1, y=V2))+
#   geom_point()


#################################

all_genesets <- append(geneset, metabolism_list_gs)


######################## Test with all genes / uncleaned dataset
# big_tumor_dfs <- list()
# tumor_types <- unique(tcga_anno$cancer_type_abbreviation)
# 
# tcga_exp_backup <- tcga_exp
# #remove accession numbers
# IDs <- c()
# for (element in rownames(tcga_exp)) {
#   IDs <- c(IDs, strsplit(element, "|", fixed = TRUE)[[1]][2])
# }
# rownames(tcga_exp) <- make.names(IDs, unique = TRUE)
# 
# tcga_exp_clean <- tcga_exp[,]
# 
# for(types in tumor_types){
#   big_tumor_dfs[[types]] <- as.data.frame(tcga_exp[,which(tcga_anno$cancer_type_abbreviation == types)])
# }
# 
# big_umap <- list()
# for(types in tumor_types){
#   big_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(big_tumor_dfs[[types]]), npcs = 25)@cell.embeddings, metric = "cosine"))
# }
# 
# big_ggplots <- list()
# for(plot in 1:length(big_umap)){
#   big_ggplots[[plot]] <- ggplot(big_umap[[plot]], aes(x = V1, y=V2))+
#           geom_point(size = 0.8)+
#           ggtitle(names(big_umap)[plot])
# }
# 
# ggarrange(plotlist = big_ggplots)
# 
# 
# 
# #gsva and umap
# 
# 
# big_gsva_list <- list()
# big_gsva_list <- mclapply(big_tumor_dfs, function(x){gsva(as.matrix(x), all_genesets, method = "zscore")}, mc.cores = 6)
# 
# all_gsva_umap <- list()
# for(types in names(all_gsva_devided)){
#   all_gsva_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(all_gsva_devided[[types]]), npcs = 25, assay = "RNA-Seq")@cell.embeddings, metric = "cosine"))
# }
# 
# all_ggplots_gsva <- list()
# for(plot in names(all_gsva_umap)){
#   all_ggplots_gsva[[plot]] <- ggplot(all_gsva_umap[[plot]], aes(x = V1, y=V2))+
#     geom_point(size = 0.8)+
#     ggtitle(plot)
# }
# 
# ggarrange(plotlist = all_ggplots_gsva)
# 

################################################## pathway activity heatmap



# pathway_enrichment_means <- matrix(nrow = 87, ncol = 33, byrow = FALSE)
# for( type in 1:length(tumor_types)){
#   pathway_enrichment_means[,type] <- apply(as.data.frame(big_gsva_list[[type]]), 1, mean)
# }
# pathway_enrichment_means <- as.data.frame(pathway_enrichment_means)
# colnames(pathway_enrichment_means) <- tumor_types
# rownames(pathway_enrichment_means) <- rownames(big_gsva_list[["LUAD"]])
# 
# heatmap <- pheatmap(as.matrix(pathway_enrichment_means), color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 18, fontsize = 18 )
# png(filename = "./output/enrichmentHeatmap.png", height = 1400, width = 2000, units = "px")
# heatmap
# dev.off()
#################################### create second geneset list with C2 genesets from msigdb
C2_pathways <- msigdbr(species = "Homo sapiens", category = "C2")

C2_list = list() #this creates an empty list to store the for loop in
for (gene.name in unique(C2_pathways$gs_name)){
  pathway = c(C2_pathways[C2_pathways$gs_name == gene.name, 7])
  names(pathway) = paste(c(strsplit(gene.name, "_", fixed = TRUE)[[1]][-1]), collapse = "_") #to extract the "GOBP_" from each name, substring is used
  C2_list[paste(c(strsplit(gene.name, "_", fixed = TRUE)[[1]][-1]), collapse = "_")] = pathway}

#append C2 sets to rest of the genesets
combined_genesets_C2 <- append(all_genesets, C2_list)

#################################### GSVA for all Tumor types
all_umap <- as.data.frame(umap(RunPCA(as.matrix(tcga_exp_short))@cell.embeddings, metric = "cosine"))

ggplot(all_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

all_gsva_C5 <- gsva(as.matrix(tcga_exp_short),all_genesets_c5, method = "gsva", min.sz = 20, parallel.sz = 4)

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


############################# delete dataset entries for genes not in pathways
tcga_exp_short <- tcga_exp[which(rownames(tcga_exp) %in% unlist(unique(all_genesets_copy))),]


############################## new GSVA with short dataset
# short_types_dfs <- list()
# for(types in tumor_types){
#   short_types_dfs[[types]] <- as.data.frame(tcga_exp_short[,which(tcga_anno$cancer_type_abbreviation == types)])
# }
# 
# short_gsva_list <- list()
# short_gsva_list <- mclapply(short_types_dfs, function(x){gsva(as.matrix(x), all_genesets_copy, method = "zscore")}, mc.cores = 6)
# 
# 
# pathway_enrichment_means <- matrix(nrow = 2941, ncol = 33, byrow = FALSE)
# for( type in 1:length(tumor_types)){
#   pathway_enrichment_means[,type] <- apply(as.data.frame(short_gsva_list[[type]]), 1, mean)
# }
# pathway_enrichment_means <- as.data.frame(pathway_enrichment_means)
# colnames(pathway_enrichment_means) <- tumor_types
# rownames(pathway_enrichment_means) <- rownames(short_gsva_list[["LUAD"]])
# 
# heatmap <- pheatmap(as.matrix(pathway_enrichment_means), color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), angle_col = "45", fontsize_row = 12, cellheight = 13, fontsize_col = 18, fontsize = 18 )
# png(filename = "./output/enrichmentHeatmap.png", height = 6000, width = 5000, units = "px")
# heatmap
# dev.off()

######################### GSVA for all tumor types at once

all_gsva <- as.data.frame(gsva(as.matrix(tcga_exp_short), all_genesets_copy, method = "zscore"))

whole_gsva_umap <- as.data.frame(umap(t(all_gsva), metric = "cosine"))
ggplot(whole_gsva_umap, aes(x=V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()

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
rest_means <- apply(pathway_enrichment_means_highSD[,-4], 1, mean)

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

foldchange_LUADvsRest <- foldchange(rest_means, pathway_enrichment_means_highSD$LUAD)
foldchange_LUADvsRest <- sapply(foldchange_LUADvsRest, function(x){log2(x)})

all_gsva_rest <- as.data.frame(all_gsva[,which(tcga_anno$cancer_type_abbreviation != "LUAD")])

pvalues_LUADvsRest <- vector()
for(geneset in rownames(pathway_enrichment_means_highSD)){
  pvalues_LUADvsRest <- append(pvalues_LUADvsRest, wilcox.test(as.numeric(all_gsva_devided$LUAD[geneset,]), as.numeric(all_gsva_rest[geneset,]), alternative = "two.sided")$p.value)
}

volcano_df <- data.frame("genset" = rownames(pathway_enrichment_means_highSD), "pvalues" = pvalues_LUADvsRest, "foldchange" = foldchange_LUADvsRest)
volcano_df_small <- volcano_df[1:50,]
EnhancedVolcano(volcano_df_small, lab = volcano_df_small$genset, x = "foldchange", y="pvalues", pointSize = 1.5, labSize = 4)
