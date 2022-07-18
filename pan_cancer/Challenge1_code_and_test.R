library(parallel)
library(tidyverse)
library(ggplot2)
library(uwot)
library(Seurat)

pvaluelist_combined <- readRDS("./data/pvaluelist_pan_exp.RDS")

#enrichment funktion
enrichment = function(expressiondata, genesets = genesets_ids){
  ESmatrix = sapply(genesets, FUN = function(x){
    ins = na.omit(which(rownames(expressiondata) %in% x == TRUE))#indices der gene im aktuellen set
    outs = -ins#indices der gene nicht im aktuellen set
    #gibt einen vektor der für jeden patienten den pval für das aktuelle gene enthält
    res = NULL
    for (i in 1:ncol(expressiondata)){#testet für jeden patienten
      res[i] = wilcox.test(expressiondata[ins,i],expressiondata[outs,i],'two.sided')$p.value
    }
    return(res)
  })
  row.names(ESmatrix) = colnames(expressiondata); return(ESmatrix)
}
pvalueslist_metabol = mclapply(tumor_type_dfs, FUN = function(x){enrichment(x,metabolism_list_gs)}, mc.cores = 6)


pvaluelist_combined <- list()
for(type in names(pvalueslist)){
  pvaluelist_combined[[type]] <- as.data.frame(cbind(pvalueslist[[type]], pvalueslist_metabol[[type]]))
}

pvaluelist_combined_copy <- pvaluelist_combined

####################### Umwandlung in Score über -log10
scaled_pvalues <- pvaluelist_combined

for(type in names(pvaluelist_combined)){
  for(i in 1:ncol(pvaluelist_combined[[type]])){
    for(j in 1:nrow(pvaluelist_combined[[type]])){
      scaled_pvalues[[type]][j,i] <- -log10(pvaluelist_combined[[type]][j,i])
    }
  }
}


######################## example umap plot with Angio_AACR
ggplot(tumor_type_umap$ESCA, aes(x = V1, y = V2, color = scaled_pvalues$ESCA$Growth_Tumor_Supp))+
  geom_point()

######################### Quality control - rerun pca and umap on new scores

qc_umap <- list()
for(types in tumor_types){
  qc_umap[[types]] <- as.data.frame(umap(RunPCA(as.matrix(t(pvaluelist_combined[[types]])))@cell.embeddings, metric = "cosine"))
}



qc_plots <- list()
for(type in names(qc_umap)){
  print(ggplot(qc_umap[[type]], aes(x = V1, y=V2))+
          geom_point()+
          ggtitle(type))
}
ggarrange(plotlist = qc_plots)

