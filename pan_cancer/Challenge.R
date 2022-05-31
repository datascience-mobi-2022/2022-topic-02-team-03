library(ggplot2)
library(reshape2)
library(gtools)
library(Seurat)


dfs <- list(tcga_exp_copy[which(tcga_anno$cancer_type_abbreviation=="LUAD"),], tcga_exp_copy[which(tcga_anno$cancer_type_abbreviation=="BRCA"),])
geneset <- genesets$genesets

create_pvalue_dfs <- function(df_list, geneset_list){
  p_vector <- vector()
  pvalue_df_list <- list()
  for (df in 1:length(df_list)){
    temp_df <- df_list[[df]]
    pvalue_df <- data.frame()
    for (patient in 1: nrow(temp_df)){
      for (geneset in 1:length(geneset_list)){
        temp_geneset <- geneset_list[[geneset]]
        genes_in <- as.vector(temp_df[patient, which(colnames(temp_df) %in% temp_geneset)])
        genes_out <- as.vector(temp_df[patient, -which(colnames(temp_df) %in% temp_geneset)])
        p_vector <- append(p_vector, wilcox.test(genes_in, genes_out, paired = FALSE)$statistic$p.value)
      }
      pvalue_df <- rbind(pvalue_df, p_vector)
      p_vector <- vector()
    }
    pvalue_df_list[df_list] <- pvalue_df 
  }
  return(pvalue_df_list)
}



test_challenge <- create_pvalue_dfs(df_list = dfs, geneset_list = geneset)


