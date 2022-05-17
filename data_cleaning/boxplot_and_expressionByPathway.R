#----------------------------------------------------
# Boxplots for TCGA_exp

library(ggplot2)
library(reshape2)
library(gtools)


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
  geom_point(aes(colour = cut(quotient, c(-Inf, -1, 1, Inf))))+
  scale_color_manual(name = "quotient",
                     values = c("(-Inf,-1]" = "blue",
                                "(-1,1]" = "grey",
                                "(1, Inf]" = "red"))+
  ggtitle("mean vs variance plot of all genes")
  
