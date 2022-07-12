#----------------------------------------------------
# Boxplots for TCGA_exp

library(ggplot2)
library(reshape2)
library(gtools)
library(Seurat)


exp_highvar <- readRDS("./data/tcga_exp_small.RDS")
genesets <- readRDS("./data/hallmarks_genesets.rds")


# create subset for LUAD patients only
LUAD_patients <- tcga_exp_copy[which(tcga_anno$cancer_type_abbreviation == "LUAD"),]

gene_means_LUAD <- apply(LUAD_patients, 2, mean)


# function to create a vector of mean gene expressions of pathways

pathway_mean_activity <- function(pathway) {
  
  means <- gene_means_LUAD[which(names(gene_means_LUAD) %in% genesets$genesets[[pathway]])]
  return(means)
  
}
pathway_activities <- list()
for(entry in names(genesets[[1]])) {
  pathway_activities <- append(pathway_activities, data.frame(entry = pathway_mean_activity(entry)))
}
names(pathway_activities) <- names(genesets$genesets)

melted_pathways <- melt(pathway_activities, id = "Pathway")
colnames(melted_pathways) <- c("value", "pathway")

ggplot(melted_pathways, aes(x = value, y = pathway )) +
  geom_boxplot(color = "dodgerblue3") +
  geom_vline(xintercept = 0)+
  theme(legend.position="none")+
  ggtitle("pathway activities by mean expression in LUAD patients")


means <- apply(tcga_exp_copy, 2, mean)
variances <- apply(tcga_exp_copy, 2, var)
mv_quotient <- variances / means
mean_variance_df <- data.frame("names" = colnames(exp_highvar), "means" = means, "variances" = variances, "quotient" = mv_quotient)

ggplot(mean_variance_df, aes(x = means, y = variances))+
  geom_point()+
  ggtitle("mean vs variance plot of all genes")+
  geom_smooth()


LUAD_df <- data.frame("means" = gene_means_LUAD, "names" = names(gene_means_LUAD))
ggplot(LUAD_df, aes(x = means, y = names))+
  geom_violin()
  


seurat_tcga <- CreateSeuratObject(tcga_exp_copy, project = "data exploration", min.cells = 3, min.features = 200)
seurat_tcga <- FindVariableFeatures(tcga_exp, selection.method = "vst", nfeatures = 500)
seurat_var_plot <- VariableFeaturePlot(seurat_tcga)
var_plot_labels <- LabelPoints(plot = seurat_var_plot, points = head(VariableFeatures(seurat_tcga), repel = TRUE))
