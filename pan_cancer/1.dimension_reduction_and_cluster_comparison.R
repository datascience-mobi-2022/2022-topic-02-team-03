#################################################
# This script was used to perform dimension reduction on the original expression data in order to find clusters and similarities between tumour types, as well as comparing intra-cancer clusters to each other
################################################

library(parallel)
library(FactoMineR)
library(Seurat)
library(ggplot2)
library(uwot)
library(ggpubr)
library(babyplots)
library(cluster)
library(EnhancedVolcano)


tcga_anno <- readRDS("./data/tcga_tumor_annotation.RDS")
tcga_exp_short <- readRDS("./data/tcga_exp_small.RDS")


#tcga_exp_pca <- PCA(exp_highvar, ncp = 2)

tcga_exp_pca <- RunPCA(as.matrix(tcga_exp_short), assay = "RNA-Seq", npcs = 20)
tcga_exp_pca <- as.data.frame(tcga_exp_pca@cell.embeddings)

ggplot(tcga_exp_pca, aes(x = PC_1, y = PC_2, color = tcga_anno$cancer_type_abbreviation)) +
  geom_point()


tcga_exp_umap <- as.data.frame(umap(tcga_exp_pca, n_threads = 8, metric = "cosine", n_components = 3))
ggplot(tcga_exp_umap, aes(x = V1, y = V2, color = tcga_anno$cancer_type_abbreviation))+
  geom_point()+
  guides(color=guide_legend(title="cancer type"))+
  theme_classic()+
  theme(legend.text=element_text(size=14))
  

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
  tumor_type_dfs[[types]] <- as.data.frame(tcga_exp_short[,which(tcga_anno$cancer_type_abbreviation == types)])
}


tumor_type_umap <- list()
for(types in tumor_types){
  tumor_type_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(tcga_exp_short[,which(tcga_anno$cancer_type_abbreviation== types)]), npcs = 20)@cell.embeddings, metric = "cosine"))
}

#tumor_type_umap <- list()
#for(types in tumor_types){
#  tumor_type_umap[[types]] <- as.data.frame(umap(exp_highvar[,which(tcga_anno$cancer_type_abbreviation== types)]))
#}

ggplot(tumor_type_umap$LUAD, aes(x = V1, y = V2, color = anno_tumortypes$LUAD$gender))+
  geom_point(size = 3)+
  guides(color=guide_legend(title="gender"))+
  theme(legend.text=element_text(size=20))

ggplots <- list()
for(plot in 1:length(tumor_type_umap)){
  ggplots[[plot]] <- ggplot(tumor_type_umap[[plot]], aes(x = V1, y=V2))+
                      geom_point(size = 1)+
                      ggtitle(names(tumor_type_umap)[plot])
}

ggarrange(plotlist = ggplots, ncol=6, nrow = 6)



ggplot(tumor_type_umap[[1]], aes(x = V1, y=V2, color=anno_tumortypes[[1]]$gender))+
  geom_point()+
  ggtitle(names(tumor_type_umap)[1])


##############################################
# k-means on UMAP for comparison of clusters within cancer types

tumor_type_kmeans <- list()
kmeans_type <- c("LUAD","BRCA", "ESCA", "SARC")

temp_means <- vector()

for (type in kmeans_type) {
  temp <- matrix(nrow= 4, ncol = nrow(tumor_type_umap[[type]]), byrow = TRUE)
  temp_means <- vector()
  for (i in 2:5){
    dist_temp <- dist(tumor_type_umap[[type]])
    temp[i-1,] <- silhouette(kmeans(as.matrix(tumor_type_umap[[type]]), centers = i, nstart = 100)$cluster, dist = dist_temp)[,3]
    temp_means <- append(temp_means, mean(temp[i-1,]))
  }
  max_temp <- order(temp_means, decreasing = TRUE)[1]+1
  tumor_type_kmeans[[type]] <- kmeans(as.matrix(tumor_type_umap[[type]]), centers = max_temp, nstart = 100)
}


for (type in kmeans_type){
  png(file = paste("/Users/paulbrunner/Pictures/BioInfo Zwischenspeicher/top4clusters/umap_clustering", type), width= 297,
      height    = 210,
      units     = "mm",
      res       = 1200)
  print(ggplot(tumor_type_umap[[type]], aes(x = V1, y = V2, color = as.character(tumor_type_kmeans[[type]]$cluster)))+
    geom_point()+
    scale_color_discrete(name="Cluster")+
    ggtitle(type))
  
  dev.off()
}

##########################################
#extract patients from exp_highvar that are in one cluster for one cancer and calculating foldchange and p-values

vulcano_list <- list()
for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    vulcano_list[[paste("Cluster",j)]][[type]] <- tumor_type_dfs[[type]][, which(tumor_type_kmeans[[type]]$cluster==j)]}
}

for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    vulcano_list[[paste("Cluster",j,"Means")]][[type]] <- apply(vulcano_list[[j]][[type]], 1, mean)}
}

for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
      vulcano_list[[paste("foldchange",j,"toRest")]][[type]] <- vulcano_list[[j+3]][[type]] - apply(tumor_type_dfs[[type]][, -which(tumor_type_kmeans[[type]]$cluster==j)], 1, mean)}
}


for(type in names(tumor_type_kmeans)){
  for(j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    for(gene in rownames(vulcano_list[[j]][[type]]))
    vulcano_list[[paste("wilcox",j,"toRest")]][[type]][gene] <- wilcox.test(as.vector(t(vulcano_list[[j]][[type]][gene,])), as.vector(t(tumor_type_dfs[[type]][, -which(tumor_type_kmeans[[type]]$cluster==j)][gene,])), paired = FALSE, alternative = "two.sided")$p.value}
}

volcano_LUAD <- data.frame("gene" = rownames(vulcano_list$`Cluster 1`$LUAD), "pval" = vulcano_list$`wilcox 1 toRest`$LUAD, "foldchange" = vulcano_list$`foldchange 1 toRest`$LUAD)
EnhancedVolcano(volcano_LUAD, lab = volcano_LUAD$gene, x = "foldchange", y="pval", pointSize = 1.5, labSize = 4, title = "Comparison of the two LUAD clusters", selectLab = c("FGB","KCNU1", "INSL4", "FGA", "CALCA", "SLC16A14", "S100P", "FGL1", "AKR1C2", "CX3CL1", "ADGRF1", "KRT32", "GLB1L3", "IVL", "SLC6A20", "CST6", "HPGDS", "IL87"))


#######################

for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    alpha <- -log10(0.005/nrow(vulcano_list[[j]][[type]]))
    
    png(file = paste("/Users/paulbrunner/Pictures/BioInfo Zwischenspeicher/Vulcanoplots_cluster/cluster_new",j, type), width= 297,
        height    = 210,
        units     = "mm",
        res       = 800)
    print(ggplot(vulcano_list[[paste("plot_prep",j)]][[type]], aes(x=foldchange, y = -log10(logpval), color = expr))+
      geom_point()+
      geom_hline(yintercept = alpha, linetype = "dashed")+
      scale_color_manual(values=c("blue", "gray", "red"))+
      ggtitle(paste(type,"comparison between cluster",j,"and the other clusters"))+
      theme(legend.position = "none"))
    
    dev.off()
  }
}

for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    distance_up <- c()
    names_up <- c()
    for (gene in 1:nrow(exp_highvar)){
    if(vulcano_list[[paste("expression",j)]][[type]][gene]=="up"){
      distance_up <- append(distance_up, sqrt((vulcano_list[[paste("plot_prep",j)]][[type]]$logpval[gene])^2+(vulcano_list[[paste("plot_prep",j)]][[type]]$foldchange[gene])^2))
      names_up <- append(names_up, rownames(exp_highvar)[gene])
      }
    }
    names(distance_up) <- names_up
    vulcano_list[[paste("upregulated_genes",j)]][[type]] <- sort(distance_up)[1:10]
   }
}

for(type in names(tumor_type_kmeans)){
  for (j in 1:length(unique(tumor_type_kmeans[[type]]$cluster))){
    distance_down <- c()
    names_down <- c()
    for (gene in 1:nrow(exp_highvar)){
      if(vulcano_list[[paste("expression",j)]][[type]][gene]=="down"){
        distance_down <- append(distance_down, sqrt((vulcano_list[[paste("plot_prep",j)]][[type]]$logpval[gene])^2+(vulcano_list[[paste("plot_prep",j)]][[type]]$foldchange[gene])^2))
        names_down <- append(names_down, rownames(exp_highvar)[gene])
      }
    }
    names(distance_down) <- names_down
    vulcano_list[[paste("downregulated_genes",j)]][[type]] <- sort(distance_down)[1:10]
  }
}
